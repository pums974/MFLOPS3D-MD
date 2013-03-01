module class_mesh
  use precision
  implicit none
  
  type mesh_grid
     !-> dimensions grid1d
     integer(ik) :: n
     !-> dimensions grid3d
     integer(ik) :: nx,ny,nz
     !-> dimensions names
     character(len=512) :: nxn="resolution_x",nyn="resolution_y",nzn="resolution_z"
     !-> mean spatial step
     real(rk) :: pas
     !-> 1d grid array
     real(rk),allocatable :: grid1d(:)
     !-> 3d grid array
     real(rk),allocatable :: grid3d(:,:,:)
     !-> grid names
     character(len=512) :: gridn='noname'
  end type mesh_grid

contains

! =======================================================================
! =======================================================================
! mesh : initialization and destruction methods
! =======================================================================
! =======================================================================

  subroutine grid1d_allocate(n,grid)
! -----------------------------------------------------------------------
! mesh : allocate or reallocate 1d grid in type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer,intent(in) :: n
    real(rk),allocatable :: grid(:)

    if (.not.allocated(grid)) then
       allocate(grid(n))
    elseif(allocated(grid)) then
       deallocate(grid) ; allocate(grid(n))
    endif

  end subroutine grid1d_allocate

  subroutine grid3d_allocate(n1,n2,n3,grid)
! -----------------------------------------------------------------------
! mesh : allocate or reallocate 3d grid in type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer,intent(in) :: n1,n2,n3
    real(rk),allocatable :: grid(:,:,:)

    if (.not.allocated(grid)) then
       allocate(grid(n1,n2,n3))
    elseif(allocated(grid)) then
       deallocate(grid)
       allocate(grid(n1,n2,n3))
    endif

  end subroutine grid3d_allocate

  subroutine mesh_init(grid,gn,choice,nx,ny,nz,n1n,n2n,n3n)
! -----------------------------------------------------------------------
! mesh : initialize type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer(ik),intent(in) :: nx,ny,nz
    type(mesh_grid),intent(inout) :: grid
    character(len=*),intent(in) :: choice
    character(len=*),optional,intent(in) :: n1n,n2n,n3n
    character(len=*),intent(in) :: gn
    real(rk) :: xa,xb,ya,yb,za,zb
    real(rk) :: xi,yi,zi,alpha,beta,delta,pi
    integer(ik) :: i,j,k
    
    !-> dimensions length
    if (choice=="x") grid%n=nx
    if (choice=="y") grid%n=ny
    if (choice=="z") grid%n=nz
    grid%nx=nx ; grid%ny=ny ; grid%nz=nz 

    !-> dimensions name
    if (present(n1n)) grid%nxn=n1n
    if (present(n2n)) grid%nyn=n2n
    if (present(n3n)) grid%nzn=n3n

    !-> grid name
     grid%gridn=gn

    !-> allocate grid1d
    call grid1d_allocate(grid%n,grid%grid1d)

    !-> allocate grid3d
    call grid3d_allocate(grid%nx,grid%ny,grid%nz,grid%grid3d)

  end subroutine mesh_init

  subroutine mesh_grid_init(grid,choice,nx,ny,nz,mpid)
! -----------------------------------------------------------------------
! mesh : initialize type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    use class_md
    implicit none
    integer(ik),intent(in) :: nx,ny,nz
    type(mesh_grid),intent(inout) :: grid
    type(mpi_data),optional :: mpid
    character(len=*),intent(in) :: choice
    integer(ik) :: dom_coord(3),nd(3)
    real(rk) :: xa,xb,ya,yb,za,zb
    real(rk) :: xi,yi,zi,alpha,beta,delta,pi,gam,ksi
    integer(ik) :: i,j,k
    real(rk) :: size_dom,num_dom,num_dom_tot,xat,xbt
    real(rk) :: r(grid%n),aux 
    
    !-> parameters
    pi=4._rk*atan(1._rk)
    if (choice=="x") then
       xa=-aux ; xb=aux
    endif
    if (choice=="y") then
       xa=-aux ; xb=aux
    endif
    if (choice=="z") then
       xa=-aux ; xb=aux
    endif
!    aux=0.5_rk
!    xa=-aux ; xb=aux
!    ya=-aux ; yb=aux
!    za=-aux ; zb=aux
!    xa=0._rk ; xb=aux
!    ya=0._rk ; yb=aux
!    za=0._rk ; zb=aux

    if (present(mpid)) then
       call md_mpi_getcoord(mpid,dom_coord)
       call md_mpi_getnumberdomains(mpid,nd)
       if (choice=="x") num_dom_tot=nd(1)
       if (choice=="y") num_dom_tot=nd(2)
       if (choice=="z") num_dom_tot=nd(3)
       if (choice=="x") num_dom=dom_coord(1)
       if (choice=="y") num_dom=dom_coord(2)
       if (choice=="z") num_dom=dom_coord(3)
       size_dom=(xb-xa)/num_dom_tot
       xa=xa+num_dom*size_dom
       xb=xa+size_dom
!       print*,choice,dom_coord,xa,xb
    endif

    !-> gridx1d
    grid%pas=(xb-xa)/(real(grid%n,rk)-1._rk)

    beta=0.1_rk
    alpha=0.99_rk
    !alpha=1._rk
    delta=1.1
    gam=3._rk
    ksi=0.5_rk*(xb-xa)
    do i=1,grid%n
       xi=xa+real(i-1,rk)*grid%pas

!       if (choice=="x") grid%grid1d(i)=xi   
!       if (choice=="y") grid%grid1d(i)=xi
!       if (choice=="z") grid%grid1d(i)=xi

!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+(xb-xa)*&
!            (sinh((2._rk*(xi-xa)/(xb-xa)-1._rk)*(1._rk-beta))&
!            /sinh(1._rk-beta)+1._rk)*0.5_rk

       if (choice=="x".or.choice=="y".or.choice=="z") grid%grid1d(i)=xi

!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+(xb-xa)*(asin(-alpha*cos(pi*(xi-xa)/(xb-xa))) &
!            /asin(alpha)+1._rk)/2._rk 

!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+(xb-xa)*(1._rk+tanh(delta*(2*(xi-xa)/(xb-xa)-1._rk)) &
!       /tanh(delta))/2._rk
      
!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+ksi*(1._rk-tanh(gam*(ksi-(xi-xa)))/tanh(gam*ksi))
    enddo


    !-> piecewise geometric stretching (only work with even number of points)

    goto 100    
    if (choice=="x".or.choice=="y") then
    r=1.1
    r(1:5)=1.1_rk ; r(6:)=1._rk
    grid%grid1d(1)=xa
    grid%grid1d(2)=xa+grid%pas
    do i=3,grid%n/2+1
       grid%grid1d(i)=grid%grid1d(i-1)+r(i)*(grid%grid1d(i-1)-grid%grid1d(i-2))
    enddo
    grid%grid1d(2:grid%n/2+1)=grid%grid1d(1)+ &
         (grid%grid1d(2:grid%n/2+1)-grid%grid1d(1))*0.5_rk*(xb-xa) &
         /(grid%grid1d(grid%n/2+1)-grid%grid1d(1))
    grid%grid1d(grid%n/2+2:grid%n)=grid%grid1d(1)+xb-grid%grid1d(grid%n/2:1:-1)
!    print*,grid%grid1d
    endif
100 continue

    !-> grid3d

    if (choice=="x") grid%grid3d(:,1,1)=grid%grid1d(:)
    if (choice=="y") grid%grid3d(1,:,1)=grid%grid1d(:)
!    if (choice=="y") then
!       do i=1,grid%nx
!          grid%grid3d(i,:,1)=grid%grid1d(:)
!       enddo
!    endif
    if (choice=="z") grid%grid3d(1,1,:)=grid%grid1d(:)

  end subroutine mesh_grid_init

! =======================================================================
! =======================================================================
! mesh : derivatives methods
! =======================================================================
! =======================================================================

  subroutine derivatives_coefficients_init(grid,dc,n,name,solver)
! -----------------------------------------------------------------------
! mesh : Initialisation of derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    use class_derivatives
    implicit none
    type(derivatives_coefficients),intent(out) :: dc
    type(mesh_grid),intent(in) :: grid
    integer(ik),intent(in) :: n
    character(len=*),optional :: name
    character(len=500) :: error_msg
    character(len=*),optional :: solver
    logical :: solve

    !-> set preference for solver type
    if (present(solver)) then 
       if (solver=='yes') then
          solve=.true.
       else
          solve=.false.
       endif
    else
       solve=.false.
    endif

    !-> compute derivatives coefficients
    if (present(name)) then
       call der_coeffs_init(grid%grid1d,dc,grid%n,solve,name)
    else
       call der_coeffs_init(grid%grid1d,dc,grid%n,solve)
    endif

  end subroutine derivatives_coefficients_init

! =======================================================================
! =======================================================================
! mesh : solver methods
! =======================================================================
! =======================================================================

  subroutine solver_coefficients_init_poisson(grid,dc,n,name)
! -----------------------------------------------------------------------
! mesh : Initialisation of solver coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    use class_solver_1d
    implicit none
    type(solver_coefficients),intent(out) :: dc
    type(mesh_grid),intent(in) :: grid
    integer(ik),intent(in) :: n
    character(len=*),optional :: name
    character(len=500) :: error_msg

    if (present(name)) then
       call solver_coeffs_init_poisson(grid%grid1d,dc,grid%n,name)
    else
       call solver_coeffs_init_poisson(grid%grid1d,dc,grid%n)
    endif

  end subroutine solver_coefficients_init_poisson

  subroutine solver_init_3d(grid1,grid2,grid3,sc,bct,name)
! -----------------------------------------------------------------------
! mesh : Initialisation of solver coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    use class_solver_3d
    implicit none
    type(solver_coeffs_3d),intent(out) :: sc
    type(mesh_grid),intent(in) :: grid1,grid2,grid3
    integer(ik),intent(in) :: bct(6)
    character(len=*),optional :: name
    character(len=500) :: error_msg

    if (present(name)) then
       call solver_coeffs_init_3d(grid1%grid1d,grid2%grid1d,grid3%grid1d,&
            sc,grid1%n,grid2%n,grid3%n,bct,name)
    else
       call solver_coeffs_init_3d(grid1%grid1d,grid2%grid1d,grid3%grid1d,&
            sc,grid1%n,grid2%n,grid3%n,bct)
    endif

  end subroutine solver_init_3d

! =======================================================================
! =======================================================================
! mesh : input and output methods
! =======================================================================
! =======================================================================

  subroutine write_mesh(file_name,grid,mpid)
! -----------------------------------------------------------------------
! mesh : create/open and add grid variable in output file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    use class_md
    use class_io
    implicit none
    type(mesh_grid),intent(in) :: grid
    character(len=*),intent(in) :: file_name
    type(mpi_data),intent(in),optional :: mpid

    !-> write grid 
    if (present(mpid)) then
       call write_var3d(file_name//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),grid%gridn,grid%grid3d,mpid=mpid)
    else
       call write_var3d(file_name//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),grid%gridn,grid%grid3d)
    endif

  end subroutine write_mesh

  subroutine read_mesh(file_name,choice,grid,grid_name)
! -----------------------------------------------------------------------
! mesh : read grid variable in input file
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    use class_io
    implicit none
    type(mesh_grid),intent(inout) :: grid
    character(len=*),intent(in) :: file_name,grid_name
    character(len=*),intent(in) :: choice
    integer(ik) :: dim_len(3)
    character(len=512) :: dim_name(3)
    
    !-> get grid information from file
    call get_var3d_info(file_name//".nc",grid_name,dim_name(1),dim_len(1))

    !-> initialize type mesh
    call mesh_init(grid,grid_name,choice,dim_len(1),dim_len(2),dim_len(3),&
         n1n=dim_name(1),n2n=dim_name(2),n3n=dim_name(3))

    !-> read grid3d
    call read_var3d(file_name//".nc",grid%gridn,grid%grid3d)
    
    !-> copy grid3d in grid1d
    if (choice=='x') grid%grid1d(:)=grid%grid3d(:,1,1)
    if (choice=='y') grid%grid1d(:)=grid%grid3d(1,:,1)
    if (choice=='z') grid%grid1d(:)=grid%grid3d(1,1,:)


  end subroutine read_mesh

end module class_mesh
