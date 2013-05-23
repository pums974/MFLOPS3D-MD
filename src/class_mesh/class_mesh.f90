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
     real(rk),allocatable :: dgrid1d(:)
     !-> 3d grid array
     real(rk),allocatable :: grid3d(:,:,:)
     real(rk),allocatable :: dgrid3d(:,:,:)
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
    call grid1d_allocate(grid%n,grid%dgrid1d)

    !-> allocate grid3d
    call grid3d_allocate(grid%nx,grid%ny,grid%nz,grid%grid3d)
    call grid3d_allocate(grid%nx,grid%ny,grid%nz,grid%dgrid3d)

  end subroutine mesh_init

  subroutine mesh_grid_init(grid,choice,nx,ny,nz,mpid,stretch1,taille1)
! -----------------------------------------------------------------------
! mesh : initialize type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    use class_md
    use mpi
    implicit none
    integer(ik),intent(in) :: nx,ny,nz
    type(mesh_grid),intent(inout) :: grid
    type(mpi_data),optional :: mpid
    character(len=*),intent(in) :: choice
    integer(ik) :: dom_coord(3),nd(3)
    real(rk) :: xa,xb,ya,yb,za,zb,x1,x2,dx1,dx2,xc,xd,dx0,xm,x0
    real(rk) :: xi,yi,zi,alpha,alpha1,alpha2,beta,delta,pi,gam,ksi
    integer(ik) :: i,j,k,n_deb,n_tot,num_dom,num_dom_tot
    real(rk) :: size_dom,xat,xbt
    real(rk) :: r(grid%n),aux 
    real(rk),dimension(:),allocatable :: x_tot
    real(rk),optional :: taille1
    logical,optional :: stretch1

    
    !-> parameters
    pi=4._rk*atan(1._rk)
    if (choice=="x") then
!       xa=0._rk ; xb=5._rk*pi
!       xa=0._rk ; xb=4._rk*pi
       xa=0._rk ; xb=3._rk*pi
!       xa=0._rk ; xb=2._rk*pi
    endif
    if (choice=="y") then
!       xa=-1._rk ; xb=1._rk
       xa=0._rk ; xb=2._rk
    endif
    if (choice=="z") then
!       xa=0._rk ; xb=3._rk*pi/2._rk
       xa=0._rk ; xb=1._rk*pi
    endif

    xc=xa ; xd=xb ! dimensions globales (pour stretch global)
    n_tot=grid%n
    n_deb=1

    if (choice=="y") then
!    beta=2.77_rk ! ratio max = 2 ratio bord = 1.7
!    beta=8.77_rk ! ratio max = 3 ratio bord = 2.2
!    beta=10.7_rk ! ratio max = 4 ratio bord = 2.5
    beta=33.3_rk ! ratio max = 5 ratio bord = 2.6
    else
!    beta=2.77_rk ! ratio max = 2 ratio bord = 1.7
!    beta=8.77_rk ! ratio max = 3 ratio bord = 2.2
!    beta=10.7_rk ! ratio max = 4 ratio bord = 2.5
    beta=33.3_rk ! ratio max = 5 ratio bord = 2.6
    endif
    alpha1=1._rk-beta/real((grid%n**2),rk)
    alpha2=4.0_rk

    if (present(mpid)) then
       if (mpid%dims.ne.0) then
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
         n_tot=num_dom_tot*grid%n-num_dom_tot+1
         n_deb=1+num_dom*(grid%n-1)
      endif
    endif

    if (choice=="y") then
      allocate(x_tot(n_tot))
      !-> gridx1d
      grid%pas=(xd-xc)/real(n_tot-1,rk)

      do i=1,n_tot
         xi=xc+real(i-1,rk)*grid%pas ! grille [xc,xd] stretch global
       if(.true.) then
         x_tot(i)=(1._rk+tanh(alpha2*((xi-xc)/(xd-xc)-0.5_rk))/tanh(alpha2))*0.5_rk
       elseif(.false.) then
         x_tot(i)=sign(1._rk-sinh(3.25_rk*(1._rk-abs(xi)))/ sinh(3.25_rk),xi)
       elseif(.false.) then

       if(i==1) then
          x_tot(i)=0._rk
       elseif(i.le.9.or.i.ge.n_tot-8) then 
          x_tot(i)=x_tot(i-1)+5e-3
       elseif(i.le.18) then
          x_tot(i)=x_tot(i-1)+5e-3+1.75e-2*(i-9)/8
       elseif(i.ge.n_tot-17) then
          x_tot(i)=x_tot(i-1)+5e-3+1.75e-2*(n_tot-9-i)/8
       else
          x_tot(i)=x_tot(i-1)+3e-2
       endif
       endif
      enddo
      x_tot=xc + (xd-xc)* (x_tot-x_tot(1))/(x_tot(n_tot)-x_tot(1))

      xa=x_tot(n_deb) 
      xb=x_tot(n_deb+grid%n-1) 
    endif
    !-> gridx1d
    grid%pas=(xb-xa)/real(grid%n-1,rk)

    grid%grid1d(1)=0._rk
    xi=xa
    do i=2,grid%n
       xm=xi 
       xi=xa+real(i-1,rk)*grid%pas ! grille [xa,xb] avant stretching
       if (choice=="y")  xi=x_tot(n_deb+i-1)

       x0=(xi-xa)/(xb-xa)
       dx0=(xi-xm)/(xb-xa)

!grille [0,1] stretche par bloc
       x1=(asin(-alpha1*cos(pi*x0))/asin(alpha1)+1._rk)*0.5_rk

!grille [0,1] stretche dans domaine
       x2=(1._rk+tanh(alpha2*(x0-0.5_rk))/tanh(alpha2))*0.5_rk


       xm=(xm+xi)*0.5_rk
!pour cummuler les stretch, derive les stretch precedent
       dx1=alpha1*pi*sin(pi*(xm-xa)/(xb-xa))/&
           (asin(alpha1)*(xb-xa)*sqrt(1._rk-alpha1**2*cos(pi*(xm-xa)/(xb-xa))**2))*0.5_rk
         dx2=alpha2*sech(alpha2*((xm-xc)/(xd-xc)-0.5_rk))**2/(tanh(alpha2)*(xd-xc))


!choix de la grille [0,1]
       if (choice=="x".or.choice=="z") then
         grid%grid1d(i)=x1
       else
!les ajoute pour faire une grille [0,?] doublement stretche
!         grid%grid1d(i)=grid%grid1d(i-1)+ dx0 * (dx1/dx0)*(dx2/dx0)
!           grid%grid1d(i)=x_tot(n_deb+i-2)-xa + dx1 -2*dx0
           grid%grid1d(i)=grid%grid1d(i-1)+(x_tot(n_deb+i-1)-x_tot(n_deb+i-2))*dx1/dx0

!         grid%grid1d(i)=x0
       endif

!        grid%grid1d(i)=xa+2._rk*xi-(xb-xa)*(sinh((2._rk*(xi-xa)/(xb-xa)-1._rk)*(1._rk-alpha))&
!            /sinh(1._rk-alpha)+1._rk)*0.5_rk
!        grid%grid1d(i)=xa+(xb-xa)*(1._rk+tanh(alpha*(2*(xi-xa)/(xb-xa)-1._rk))/tanh(alpha))/2._rk
        
!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+(xb-xa)*(asin(-alpha*cos(pi*(xi-xa)/(xb-xa))) &
!            /asin(alpha)+1._rk)/2._rk 

!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+(xb-xa)*(1._rk+tanh(delta*(2*(xi-xa)/(xb-xa)-1._rk)) &
!       /tanh(delta))/2._rk
      
!       if (choice=="x".or.choice=="y".or.choice=="z") &
!            grid%grid1d(i)=xa+ksi*(1._rk-tanh(gam*(ksi-(xi-xa)))/tanh(gam*ksi))
    enddo
!transforme grille [0,?] en grille [xa,xb]
    grid%grid1d=xa + (xb-xa)*(grid%grid1d-grid%grid1d(1))/(grid%grid1d(grid%n)-grid%grid1d(1))

    !-> piecewise geometric stretching (only work with even number of points)

    goto 100  

!->alpha + delta -> k
!    alpha=cmd%stretch_value2 ! ratio bord
!    delta=cmd%stretch_value3  ! ratio max
!    j=cmd%stretch_value1
!    delta=3.5_rk  ! ratio max
!    k=max(0,grid%n/2-1)    ! nb mailles max de stretch


!-> alpha +k -> delta
!    alpha=cmd%stretch_value2 ! ratio bord
!    k=cmd%stretch_value1    ! nb mailles max de stretch
!    delta=0  ! ratio max
!    j=0


 !->k + delta -> alpha
    delta=0._rk !1.1_rk  ! ratio max
    k=1    ! nb mailles max de stretch
    j=0
    alpha=1.01_rk !delta**(1._rk/k) ! ratio bord

    k=min(k,grid%n/2-1)    ! nb mailles max de stretch
200 r=1._rk
    beta=1._rk    ! ratio max = prod r(i)
    do i=1,k
      r(i)= 1._rk + (alpha-1._rk)*((real(k-1,rk)-i)/real(k-2,rk))**j
      beta=beta*r(i)
      if(delta.gt.0._rk) then
      if (beta.gt.(delta+1e-8)) then
        k=k-1
        goto 200
      endif
      endif
    enddo
    grid%grid1d(1)=0._rk
    grid%grid1d(2)=grid%pas
    do i=3,(grid%n)/2+1
       grid%grid1d(i)=grid%grid1d(i-1)+r(i-2)*(grid%grid1d(i-1)-grid%grid1d(i-2))
    enddo
    do i=(grid%n)/2+2,grid%n
      grid%grid1d((grid%n)/2+2:grid%n)=grid%grid1d((grid%n)/2+1)+grid%grid1d((grid%n+1)/2)-grid%grid1d(grid%n-(grid%n)/2-1:1:-1)
    enddo
    grid%grid1d(:)=xa+ grid%grid1d(:)*(xb-xa)  /grid%grid1d(grid%n)
!    print*,grid%grid1d
100 continue

    !-> dgrid1d
    do i=1,grid%n
       if(i.eq.1) then
         grid%dgrid1d(i)= grid%grid1d(i+1)- grid%grid1d(i)
       elseif(i.eq.grid%n) then
         grid%dgrid1d(i)= grid%grid1d(i  )- grid%grid1d(i-1)
       else
         grid%dgrid1d(i)=(grid%grid1d(i+1)- grid%grid1d(i-1))*0.5_rk
       endif
    enddo

    !-> grid3d
    if (choice=="x") then
       do k=1,grid%nz
          do j=1,grid%ny
             grid%grid3d(:,j,k)=grid%grid1d(:)
             grid%dgrid3d(:,j,k)=grid%dgrid1d(:)
          enddo
       enddo
    endif

    if (choice=="y") then
       do k=1,grid%nz
          do i=1,grid%nx
             grid%grid3d(i,:,k)=grid%grid1d(:)
             grid%dgrid3d(i,:,k)=grid%dgrid1d(:)
          enddo
       enddo
    endif

    if (choice=="z") then
       do j=1,grid%ny
          do i=1,grid%nx
             grid%grid3d(i,j,:)=grid%grid1d(:)
             grid%dgrid3d(i,j,:)=grid%dgrid1d(:)
          enddo
       enddo
    endif

xi=maxval(grid%dgrid1d(1:grid%n))
x1=xi
    if (present(mpid)) then
    if (mpid%dims.ne.0) then
    call mpi_reduce(xi,x1,1,mpi_double_precision,mpi_max,0,&
         mpi_comm_world,mpid%code)
    endif
    endif

xi=minval(grid%dgrid1d(1:grid%n))
x2=xi
    if (present(mpid)) then
    if (mpid%dims.ne.0)  then
    call mpi_reduce(xi,x2,1,mpi_double_precision,mpi_min,0,&
         mpi_comm_world,mpid%code)
   endif
   endif

    if (present(mpid)) then
if (mpid%rank==0) then
!     write(*,'(2A)') trim(grid%gridn),"    rmax              rbord"
!     write(*,'(2es17.8)') &
!         grid%dgrid1d(grid%n/2)/grid%dgrid1d(grid%n), &
!         grid%dgrid1d(grid%n-1)/grid%dgrid1d(grid%n)
     write(*,'(2A)') trim(grid%gridn),"    dxmax              dxmin"
     write(*,'(2es17.8)') x1, x2 
endif
   endif

  end subroutine mesh_grid_init

! =======================================================================
! =======================================================================
! mesh : derivatives methods
! =======================================================================
! =======================================================================

  subroutine derivatives_coefficients_init(grid,dc,n,so,name,solver)
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
    integer(ik),intent(in) :: n,so
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
       call der_coeffs_init(grid%grid1d,dc,grid%n,solve,so,name)
    else
       call der_coeffs_init(grid%grid1d,dc,grid%n,solve,so)
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

  subroutine solver_init_3d(grid1,grid2,grid3,sc,bct,so,name)
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
    integer(ik),intent(in) :: bct(6),so
    character(len=*),optional :: name
    character(len=500) :: error_msg

    if (present(name)) then
       call solver_coeffs_init_3d(grid1%grid1d,grid2%grid1d,grid3%grid1d,&
            sc,grid1%n,grid2%n,grid3%n,bct,so,name)
    else
       call solver_coeffs_init_3d(grid1%grid1d,grid2%grid1d,grid3%grid1d,&
            sc,grid1%n,grid2%n,grid3%n,bct,so)
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
      if(mpid%dims.ne.0) then
       call write_var3d(trim(file_name)//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),grid%gridn,grid%grid3d,mpid=mpid)
       call write_var3d(trim(file_name)//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),'d'//grid%gridn,grid%dgrid3d,mpid=mpid)
    else
       call write_var3d(trim(file_name)//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),grid%gridn,grid%grid3d)
       call write_var3d(trim(file_name)//".nc",&
            [character(len=512) :: grid%nxn,grid%nyn,grid%nzn],&
            get_dim_size(grid%grid3d),'d'//grid%gridn,grid%dgrid3d)
    endif
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
    integer(ik) :: dim_len(3),i
    character(len=512) :: dim_name(3)
    
    !-> get grid information from file
    call get_var3d_info(file_name//".nc",grid_name,dim_name(1),dim_len(1))
    dim_len(1)=max(dim_len(1),1)
    dim_len(2)=max(dim_len(2),1)
    dim_len(3)=max(dim_len(3),1)

    !-> initialize type mesh
    call mesh_init(grid,grid_name,choice,dim_len(1),dim_len(2),dim_len(3),&
         n1n=dim_name(1),n2n=dim_name(2),n3n=dim_name(3))

    !-> read grid3d
    call read_var3d(file_name//".nc",grid%gridn,grid%grid3d)
    
    !-> copy grid3d in grid1d
    if (choice=='x') grid%grid1d(:)=grid%grid3d(:,1,1)
    if (choice=='y') grid%grid1d(:)=grid%grid3d(1,:,1)
    if (choice=='z') grid%grid1d(:)=grid%grid3d(1,1,:)

    !-> dgrid1d
    do i=1,grid%n
       if(i.eq.1) then
         grid%dgrid1d(i)= grid%grid1d(i+1)- grid%grid1d(i)
       elseif(i.eq.grid%n) then
         grid%dgrid1d(i)= grid%grid1d(i  )- grid%grid1d(i-1)
       else
         grid%dgrid1d(i)=(grid%grid1d(i+1)- grid%grid1d(i-1))*0.5_rk
       endif
    enddo

    !-> dgrid3d
    if (choice=="x") grid%dgrid3d(:,1,1)=grid%dgrid1d(:)
    if (choice=="y") grid%dgrid3d(1,:,1)=grid%dgrid1d(:)
    if (choice=="z") grid%dgrid3d(1,1,:)=grid%dgrid1d(:)

  end subroutine read_mesh

function sech(x)
real(rk)::x,sech
sech=2._rk*exp(-x)/(1._rk+exp(-2._rk*x))
end function sech
end module class_mesh
