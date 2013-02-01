module class_io
! -----------------------------------------------------------------------
! Name :
! class io 
! -----------------------------------------------------------------------
  use netcdf
  use class_field
  use class_mesh
  implicit none
  
  !-> Declare everything private by default
  private

  !-> Declare exported procedure
  public :: write_field1,write_grid
  public :: read_grid
  public :: write_var3d

contains

  subroutine write_field1(name,x)
! -----------------------------------------------------------------------
! io : create/open and add data in netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    type(field),intent(in) :: x
    character(len=*),intent(in) :: name
    logical :: file_exists
    
    inquire(file=name,exist=file_exists)

    if (file_exists) then
       call write_field_add(name,x)
    else
       call write_field_init(name,x)
    endif

  end subroutine write_field1

  subroutine write_var3d(file_name,ndim,dim_name,dim_len,var_name,var)
! -----------------------------------------------------------------------
! io : write 3d variable in a netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    character(len=*),intent(in) :: file_name,var_name
    integer(ik),intent(in) :: ndim
    character(len=*),intent(in) :: dim_name(ndim)

    real(rk) :: var(:,:,:)
    logical :: file_exist
    integer(ik) :: varid(1),i
    integer(ik) :: ncid,dim_len(ndim),dimid(ndim),dim_len_check
    
    !-> check if file exist
    inquire(file=file_name,exist=file_exist)

    !-> open/create file
    if (file_exist) then
       call io_check(nf90_open(path=file_name,mode=nf90_write,ncid=ncid))
       call io_check(nf90_redef(ncid))
    else
       call io_check(nf90_create(path=file_name,cmode=nf90_clobber,ncid=ncid))
    endif

    !-> create/add dimensions
    do i=1,ndim
       if (nf90_inq_dimid(ncid,dim_name(i),dimid(i))/=nf90_noerr) then 
          call io_check(nf90_def_dim(ncid,dim_name(i),dim_len(i),dimid(i)))
       else
          call io_check(nf90_inquire_dimension(ncid,dimid(i),len=dim_len_check))
          if (dim_len_check/=dim_len(i)) call error_stop("NETCDF Error : wrong dimensions")
       endif
    enddo

    !-> if variable exist         : get variable id
    !-> if variable doesn't exist : define variable
    if (nf90_inq_varid(ncid,var_name,varid(1))/=nf90_noerr) then 
       call io_check(nf90_def_var(ncid,var_name,nf90_real,dimid,varid(1)))
    endif

    !-> end of definition
    call io_check(nf90_enddef(ncid))

    !-> write field variable
    call io_check(nf90_put_var(ncid,varid(1),var,start=(/1,1,1/),count=get_dim_size(var)))

    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine write_var3d

  subroutine write_field_init(name,x)
! -----------------------------------------------------------------------
! io : create a netcdf file and write one variable 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    type(field),intent(in) :: x
    character(len=*),intent(in) :: name
    integer :: ncid, status
    integer :: dimid(4),varid(1)

    !-> create file
!    call io_check(nf90_create(path=name,cmode=nf90_hdf5,ncid=ncid))
    call io_check(nf90_create(path=name,cmode=nf90_clobber,ncid=ncid))

    !-> begin of definition

    !-> define dimensions
    call io_check(nf90_def_dim(ncid,'resolution_x',x%nx,dimid(1)))
    call io_check(nf90_def_dim(ncid,'resolution_y',x%ny,dimid(2)))
    call io_check(nf90_def_dim(ncid,'resolution_z',x%nz,dimid(3)))
!    call io_check(nf90_def_dim(ncid,'resolution_z',nf90_unlimited,dimid(3)))
    call io_check(nf90_def_dim(ncid,'dim_one',1,dimid(4)))

    !-> define field variable
    call io_check(nf90_def_var(ncid,x%name,nf90_real,put_id(x%f,dimid),varid(1)))
    
    !-> end of definition
    call io_check(nf90_enddef(ncid))

    !-> write field variable
    call io_check(nf90_put_var(ncid,varid(1),x%f,start=(/1,1,1/),count=get_dim_size(x%f)))
    
    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine write_field_init

  subroutine write_grid(name,grid)
! -----------------------------------------------------------------------
! io : create a netcdf file and write grid
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    type(mesh_grid),intent(in) :: grid
    character(len=*),intent(in) :: name
    integer :: ncid, status
    integer :: dimid(4),varid(3)

    !-> create file
!    call io_check(nf90_create(path=name,cmode=nf90_hdf5,ncid=ncid))
    call io_check(nf90_create(path=name,cmode=nf90_clobber,ncid=ncid))

    !-> begin of definition

    !-> define dimensions
    call io_check(nf90_def_dim(ncid,'resolution_x',grid%nx,dimid(1)))
    call io_check(nf90_def_dim(ncid,'resolution_y',grid%ny,dimid(2)))
    call io_check(nf90_def_dim(ncid,'resolution_z',grid%nz,dimid(3)))
    call io_check(nf90_def_dim(ncid,'dim_one',1,dimid(4)))

    !-> define grid variables    
    call io_check(nf90_def_var(ncid,grid%gridx3dn,nf90_real,put_id(grid%gridx3d,dimid),varid(1)))
    call io_check(nf90_def_var(ncid,grid%gridy3dn,nf90_real,put_id(grid%gridy3d,dimid),varid(2)))
    call io_check(nf90_def_var(ncid,grid%gridz3dn,nf90_real,put_id(grid%gridz3d,dimid),varid(3)))

    !-> end of definition
    call io_check(nf90_enddef(ncid))

    !-> write grid variables  
    call io_check(nf90_put_var(ncid,varid(1),grid%gridx3d,start=(/1,1,1/), &
         count=get_dim_size(grid%gridx3d)))
    call io_check(nf90_put_var(ncid,varid(2),grid%gridy3d,start=(/1,1,1/), &
         count=get_dim_size(grid%gridy3d)))
    call io_check(nf90_put_var(ncid,varid(3),grid%gridz3d,start=(/1,1,1/), &
         count=get_dim_size(grid%gridz3d)))
    
    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine write_grid

  subroutine write_field_add(name,x)
! -----------------------------------------------------------------------
! io : overwrite or add variable in existing netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    type(field),intent(in) :: x
    character(len=*),intent(in) :: name
    integer :: ncid, status
    integer :: dimid(4),varid(1)

    !-> open file
    call io_check(nf90_open(path=name,mode=nf90_write,ncid=ncid))

    !-> get dimensions id
    call io_check(nf90_inq_dimid(ncid,'resolution_x',dimid(1)))
    call io_check(nf90_inq_dimid(ncid,'resolution_y',dimid(2)))
    call io_check(nf90_inq_dimid(ncid,'resolution_z',dimid(3)))
    call io_check(nf90_inq_dimid(ncid,'dim_one',dimid(4)))

    !-> if variable exist         : get variable id
    !-> if variable doesn't exist : define variable
    if (nf90_inq_varid(ncid,x%name,varid(1))/=nf90_noerr) then 

       !-> begin of redefinition
       call io_check(nf90_redef(ncid))
       
       !-> define field variable
       call io_check(nf90_def_var(ncid,x%name,nf90_real,put_id(x%f,dimid),varid(1)))
    
       !-> end of definition
       call io_check(nf90_enddef(ncid))

    endif

    !-> write field variable
    call io_check(nf90_put_var(ncid,varid(1),x%f,start=(/1,1,1/),count=get_dim_size(x%f)))
    
    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine write_field_add

  subroutine read_grid(name,grid)
! -----------------------------------------------------------------------
! io : create a netcdf file and write grid
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    type(mesh_grid),intent(out) :: grid
    character(len=*),intent(in) :: name
    integer :: ncid, status
    integer :: dimid(4),varid(3)

    !-> open file
    call io_check(nf90_open(path=name,mode=nf90_nowrite,ncid=ncid))

    !-> get dimensions id
    call io_check(nf90_inq_dimid(ncid,'resolution_x',dimid(1)))
    call io_check(nf90_inq_dimid(ncid,'resolution_y',dimid(2)))
    call io_check(nf90_inq_dimid(ncid,'resolution_z',dimid(3)))
    call io_check(nf90_inq_dimid(ncid,'dim_one',dimid(4)))
    
    !-> get dimensions length

    !-> get variables id
    call io_check(nf90_inq_varid(ncid,grid%gridx3dn,varid(1)))
    call io_check(nf90_inq_varid(ncid,grid%gridy3dn,varid(2)))
    call io_check(nf90_inq_varid(ncid,grid%gridz3dn,varid(3)))

    !-> get values
    print*,grid%gridx3d
    call io_check(nf90_get_var(ncid,varid(1),grid%gridx3d))
    call io_check(nf90_get_var(ncid,varid(2),grid%gridy3d))
    call io_check(nf90_get_var(ncid,varid(3),grid%gridz3d))
    print*,grid%gridx3d

    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine read_grid

  function get_dim_size(x)
! -----------------------------------------------------------------------
! io : return the size of each dimension for 3d field
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer :: get_dim_size(3)
    real(rk),intent(in) :: x(:,:,:)
    integer(ik) :: i
    do i=1,3
       get_dim_size(i)=size(x,dim=i)
    enddo
  end function get_dim_size

  function put_id(x,id)
! -----------------------------------------------------------------------
! io : return correct id for each dimension
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer :: put_id(3)
    real(rk),intent(in) :: x(:,:,:)
    integer :: n(3)
    integer :: id(4)
    
    n=get_dim_size(x)
    !n1=size(x,dim=1) ; n2=size(x,dim=2) ; n3=size(x,dim=3) 
    put_id=id(4)
    if (n(1)/=1) put_id(1)=id(1)
    if (n(2)/=1) put_id(2)=id(2)
    if (n(3)/=1) put_id(3)=id(3)

  end function put_id

  subroutine io_check(status)
! -----------------------------------------------------------------------
! io : check netcdf error output and stop code if needed 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    integer,intent(in) :: status
     
    if(status /= nf90_noerr) then
       print*,trim(nf90_strerror(status))
       print'(a)',"Netcdf Error : aborting"
       stop
    end if
  end subroutine io_check

end module class_io
