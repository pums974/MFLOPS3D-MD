module class_io
! -----------------------------------------------------------------------
! Name :
! class io 
! -----------------------------------------------------------------------
  use netcdf
!  use class_field
!  use class_mesh
  use precision
  implicit none
  
  !-> Declare everything private by default
  private

  !-> Declare exported procedure
  public :: write_var3d,read_var3d
  public :: get_dim_size
  public :: get_var3d_info

contains

  subroutine write_var3d(file_name,dim_name,dim_len,var_name,var,mode,mpid,&
       dbl,inter)
! -----------------------------------------------------------------------
! io : write 3d variable in a netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    use mpi
    use class_md
    implicit none
    character(len=*),intent(in) :: file_name,var_name
    integer(ik),parameter :: ndim=3
    character(len=*),intent(in) :: dim_name(ndim)
    character(len=*),intent(in),optional :: mode
    type(mpi_data),optional :: mpid
    character(*),optional :: dbl,inter

    real(rk) :: var(:,:,:)
    logical :: file_exist
    integer(ik) :: varid(1),i
    integer(ik) :: ncid,dim_len(ndim),dimid(ndim),dim_len_check
    integer(ik) :: dimt(3),coord(3,2)
    integer(ik) :: startv(3),countv(3)
    
    !-> check if file exist
    inquire(file=file_name,exist=file_exist)

    !-> open/create file
    if (file_exist) then
       if (present(mpid)) then
          call io_check(nf90_open(path=file_name,&
!               mode=IOR(NF90_WRITE,NF90_MPIPOSIX),ncid=ncid,&
               mode=IOR(NF90_WRITE,NF90_MPIIO),ncid=ncid,&
               comm=mpid%comm,info=MPI_INFO_NULL))
       else
          call io_check(nf90_open(path=file_name,mode=nf90_write,ncid=ncid))
       endif
       call io_check(nf90_redef(ncid))
    else
       if (present(mpid)) then
          call io_check(nf90_create(path=file_name,&
!               cmode=IOR(NF90_NETCDF4,NF90_MPIPOSIX),ncid=ncid,&
               cmode=IOR(NF90_NETCDF4,NF90_MPIIO),ncid=ncid,&
               comm=mpid%comm,info=MPI_INFO_NULL))
       else
          call io_check(nf90_create(path=file_name,cmode=nf90_clobber,ncid=ncid))
       endif
    endif

    !-> recompute dimensions, start and count if mpi
    if (present(mpid)) then
       if (present(inter)) then
          call md_mpi_global_coord(mpid,dimt,coord,inter=inter)
       else
          call md_mpi_global_coord(mpid,dimt,coord)
       endif
       startv=(/1,1,1/)
       countv=(/1,1,1/)
       do i=1,3
          if (dim_len(i)>1) then 
             dim_len(i)=dimt(i)
             startv(i)=coord(i,1)
             countv(i)=coord(i,2)
          endif
       enddo
    else
       startv=(/1,1,1/)
       countv=get_dim_size(var)
    endif

    !-> create/add dimensions
    do i=1,ndim
       if (nf90_inq_dimid(ncid,dim_name(i),dimid(i))/=nf90_noerr) then 
          call io_check(nf90_def_dim(ncid,dim_name(i),dim_len(i),dimid(i)))
       else
          call io_check(nf90_inquire_dimension(ncid,dimid(i),len=dim_len_check))
          if (dim_len_check/=dim_len(i)) &
               call error_stop("NETCDF Error : wrong dimensions")
       endif
    enddo

    !-> if variable exist         : get variable id
    !-> if variable doesn't exist : define variable
    if (nf90_inq_varid(ncid,var_name,varid(1))/=nf90_noerr) then
       if (present(dbl)) then
          call io_check(nf90_def_var(ncid,var_name,nf90_double,dimid,varid(1)))
       else
          call io_check(nf90_def_var(ncid,var_name,nf90_real,dimid,varid(1)))
       endif
    endif

    !-> end of definition
    call io_check(nf90_enddef(ncid))

    !-> write field variable
    call io_check(nf90_put_var(ncid,varid(1),var,start=startv,count=countv))

    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine write_var3d

  subroutine read_var3d(file_name,var_name,var,mpid,inter)
! -----------------------------------------------------------------------
! io : read 3d variable in a netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    use mpi
    use class_md
    implicit none
    character(len=*),intent(in) :: file_name,var_name
    real(rk) :: var(:,:,:)
    integer(ik) :: varid(1),i
    integer(ik) :: ncid

    integer(ik),parameter :: ndim=3
    type(mpi_data),optional :: mpid
    integer(ik) :: dim_len(ndim),dimid(ndim),dim_len_check
    integer(ik) :: dimt(3),coord(3,2)
    integer(ik) :: startv(3),countv(3)
    character(*),optional :: inter
    
    !-> open file
    if (present(mpid)) then
       call io_check(nf90_open(path=file_name,&
!               mode=IOR(NF90_WRITE,NF90_MPIPOSIX),ncid=ncid,&
            mode=IOR(NF90_WRITE,NF90_MPIIO),ncid=ncid,&
            comm=mpid%comm,info=MPI_INFO_NULL))
    else
       call io_check(nf90_open(path=file_name,mode=nf90_write,ncid=ncid))
    endif
    
    !-> get variable id
    call io_check(nf90_inq_varid(ncid,var_name,varid(1)))

    !-> recompute dimensions, start and count if mpi
    if (present(mpid)) then
       if (present(inter)) then
          call md_mpi_global_coord(mpid,dimt,coord,inter=inter)
       else
          call md_mpi_global_coord(mpid,dimt,coord)
       endif
       startv=(/1,1,1/)
       countv=(/1,1,1/)
       do i=1,3
          !if (dim_len(i)>1) then 
             dim_len(i)=dimt(i)
             startv(i)=coord(i,1)
             countv(i)=coord(i,2)
          !endif
       enddo
    else
       startv=(/1,1,1/)
       countv=get_dim_size(var)
    endif

    !-> read field variable
!    call io_check(nf90_get_var(ncid,varid(1),var))
    call io_check(nf90_get_var(ncid,varid(1),var,start=startv,count=countv))

    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine read_var3d

  subroutine get_var3d_info(file_name,var_name,dim_name,dim_len)
! -----------------------------------------------------------------------
! io : get 3d variable dimension length in a netcdf file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    character(len=*),intent(in) :: file_name,var_name
    integer(ik),parameter :: ndim=3
    character(len=*),intent(inout) :: dim_name(ndim)

    integer(ik) :: varid(1)
    integer(ik) :: i,ncid,dim_len(ndim),dimid(ndim)
    
    !-> open file
    call io_check(nf90_open(path=file_name,mode=nf90_write,ncid=ncid))
    
    !-> get variable id
    call io_check(nf90_inq_varid(ncid,var_name,varid(1)))

    !-> get variable dimensions id
    call io_check(nf90_inquire_variable(ncid,varid(1),dimids=dimid))

    !-> get dimensions length
    do i=1,ndim
       call io_check(nf90_inquire_dimension(ncid,dimid(i),name=dim_name(i),len=dim_len(i)))
    enddo

    !-> close file
    call io_check(nf90_close(ncid))

  end subroutine get_var3d_info

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
