module class_navier_3D
! -----------------------------------------------------------------------
! Name :
! class navier_3D
! -----------------------------------------------------------------------
! Object :
! Resolution of incompressible Navier-Stokes equations in 3D
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
  use class_field
  use class_mesh
  use class_derivatives
  use class_solver_3d
  use class_md
  use color_print
!  use command_line
  implicit none

  !-> time scheme order
  integer(ik),parameter :: nt=4
  integer(ik) :: nullv=0

  type navier3d
     !-> number of time steps of the time scheme
     integer(ik) :: nt=nt
     !-> order of spatial discretization
     integer(ik) :: so
     !-> dimensions
     integer(ik) :: nx,ny,nz
     !-> velocity fields
     type(field) :: u(nt),v(nt),w(nt),sub_u,sub_v,sub_w,sub_p
     !-> pressure
     type(field) :: p(nt),phi(nt)
     !-> velocity rhs
     type(field) :: fu(nt),fv(nt),fw(nt),rhs_u(2),rhs_v(2),rhs_w(2),rhs_px,rhs_py,rhs_pz
     !-> pressure rhs
     type(field) :: fp(nt),fphi
     !-> turbulent viscosity
     type(field):: les_nu,du(3,3),les_f(nt)
     !-> auxiliary field
     type(field) :: aux,aux1,aux2
     !-> velocity boundary conditions
     type(boundary_condition) :: bcu(nt),bcv(nt),bcw(nt)
     !-> pressure boundary conditions
     type(boundary_condition) :: bcp(nt),bcphi(nt)
     !-> mesh
     type(mesh_grid) :: gridx,gridy,gridz
     !-> derivatives coefficients
     type(derivatives_coefficients) :: dcx,dcy,dcz
     !-> sigma
     real(rk) :: sigmau,sigmap
     !-> solver coefs
     type(solver_coeffs_3d) :: scu,scv,scw,scp
     !-> influence matrixes
!     type(mpi_inf_mat) :: infu,infv,infw,infp
     type(mpi_inf_mat),pointer :: infu,infv,infw,infp
     !-> guess sol
     type(mpi_inf_sol) :: infsolu,infsolv,infsolw,infsolphi
     !-> time, time-step, and number of time steps
     real(rk) :: time,ts,fac(nt)
     integer(ik) :: it(nt),ntime
     !-> reynolds number
     real(rk) :: rey
     !-> constant**2 for turbulent viscosity
     real(rk) :: les_c
     !-> nonlinear type, projection type, time order
     integer(ik) :: nlt,pt,tou,top
     integer(ik) :: subite,nsubite,les_type
  end type navier3d

contains

  subroutine restart_write(mpid,nav)
! -----------------------------------------------------------------------
! navier : restart write
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: it(nav%nt),nt,i
    character(1) :: stepn
    character(100) :: dirname
    character(100) :: timen,proc
    logical :: exist
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    
    !-> create dirname and test existence
    dirname='restart_'//trim(tn(nav))
    inquire(file=dirname,exist=exist)
    if (.not.exist) then
       if (mpid%rank==0) then
          call system('mkdir '//trim(dirname))
       endif
    endif

    if (mpid%rank==0) then
       open(10,file=trim(dirname)//'/param.dat',access='stream')
       write(10)nav%time
       write(10)nav%it(:)
       close(10)
    endif

    do i=1,nt
       write(stepn,'(i0)')i
       call write_field(trim(dirname)//'/vel_u_'//stepn,nav%u(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/vel_v_'//stepn,nav%v(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/vel_w_'//stepn,nav%w(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/vel_p_'//stepn,nav%p(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/vel_phi_'//stepn,nav%phi(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/rhs_u_'//stepn,nav%fu(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/rhs_v_'//stepn,nav%fv(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/rhs_w_'//stepn,nav%fw(it(i)),&
            mpid,dbl='y',inter='y')
       call write_field(trim(dirname)//'/rhs_p_'//stepn,nav%fp(it(i)),&
            mpid,dbl='y',inter='y')
    enddo
    call write_field(trim(dirname)//'/rhs_phi_'//stepn,nav%fphi,&
         mpid,dbl='y',inter='y')

    if (mpid%dims.ne.0) then
      call md_influence_guess_write(mpid,nav%infsolu,trim(dirname)//'/guess_u.dat')
      call md_influence_guess_write(mpid,nav%infsolv,trim(dirname)//'/guess_v.dat')
      call md_influence_guess_write(mpid,nav%infsolw,trim(dirname)//'/guess_w.dat')
      call md_influence_guess_write(mpid,nav%infsolphi,trim(dirname)//'/guess_phi.dat')
    endif

    call write_field(trim(dirname)//'/vel_u',nav%u(it(nt)),mpid)
    call write_field(trim(dirname)//'/vel_v',nav%v(it(nt)),mpid)
    call write_field(trim(dirname)//'/vel_w',nav%w(it(nt)),mpid)
    call write_field(trim(dirname)//'/vel_p',nav%w(it(nt)),mpid)

  end subroutine restart_write

  subroutine restart_read(mpid,nav,exist)
! -----------------------------------------------------------------------
! navier : restart read
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: it(nav%nt),nt,i
    character(1) :: stepn
    character(100) :: dirname
    character(100) :: timen,proc
    logical :: exist
    integer(ik) :: l,m,c(3),inter(3,2)

    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infu,c,inter)
    else
      inter=0
    endif

    !-> create dirname and test existence
    !write(timen,'(i8.8,f0.8)')floor(nav%time),nav%time-floor(nav%time)
    dirname='restart'
    inquire(file=trim(dirname)//'/param.dat',exist=exist)
    if (.not.exist) then
       return
    endif

    !-> read parameters
    open(10,file=trim(dirname)//'/param.dat',access='stream')
    read(10)nav%time
    read(10)nav%it(:)
    close(10)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    
    do i=1,nt
       write(stepn,'(i0)')i
       call read_field(trim(dirname)//'/vel_u_'//stepn,nav%u(it(i)),&
            nav%u(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/vel_v_'//stepn,nav%v(it(i)),&
            nav%v(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/vel_w_'//stepn,nav%w(it(i)),&
            nav%w(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/vel_p_'//stepn,nav%p(it(i)),&
            nav%p(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/vel_phi_'//stepn,nav%phi(it(i)),&
            nav%phi(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/rhs_u_'//stepn,nav%fu(it(i)),&
            nav%fu(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/rhs_v_'//stepn,nav%fv(it(i)),&
            nav%fv(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/rhs_w_'//stepn,nav%fw(it(i)),&
            nav%fw(it(i))%name,mpid,inter='y')
       call read_field(trim(dirname)//'/rhs_p_'//stepn,nav%fp(it(i)),&
            nav%fp(it(i))%name,mpid,inter='y')
    enddo
    call read_field(trim(dirname)//'/rhs_phi_'//stepn,nav%fphi,&
         nav%fphi%name,mpid,inter='y')

    do i=1,nt
       call boundary_put_field(nav%u(it(i)),nav%bcu(it(i)),inter)
       call boundary_put_field(nav%v(it(i)),nav%bcv(it(i)),inter)
       call boundary_put_field(nav%w(it(i)),nav%bcw(it(i)),inter)
       call boundary_put_field(nav%p(it(i)),nav%bcp(it(i)),inter)
       call boundary_put_field(nav%phi(it(i)),nav%bcphi(it(i)),inter)
    enddo
    if (mpid%dims.ne.0) then
      call md_influence_guess_read(mpid,nav%infsolu,trim(dirname)//'/guess_u.dat')
      call md_influence_guess_read(mpid,nav%infsolv,trim(dirname)//'/guess_v.dat')
      call md_influence_guess_read(mpid,nav%infsolw,trim(dirname)//'/guess_w.dat')
      call md_influence_guess_read(mpid,nav%infsolphi,trim(dirname)//'/guess_phi.dat')
    endif
  end subroutine restart_read

 subroutine navier_write_fields(mpid,nav,inc,ite)
! -----------------------------------------------------------------------
! navier : restart write
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: inc,ite


    !-> write fields
!    if (mod(ite,inc)==0.and.ite>0) then
    if (mod(ite,inc)==0) then
       if (mpid%rank==0) print*,'Write field'
       !-> write fields
       call write_field('vel_u_'//trim(tn(nav)),nav%u(nav%it(nav%nt)),mpid)
       call write_field('vel_v_'//trim(tn(nav)),nav%v(nav%it(nav%nt)),mpid)
       call write_field('vel_w_'//trim(tn(nav)),nav%w(nav%it(nav%nt)),mpid)
       call write_field('vel_p_'//trim(tn(nav)),nav%p(nav%it(nav%nt)),mpid)
       call write_field('vel_phi_'//trim(tn(nav)),nav%phi(nav%it(nav%nt)),mpid)
       call write_field('les_nu_'//trim(tn(nav)),nav%les_nu,mpid)

!       call write_mesh('grid_x_'//trim(tn(nav)),nav%gridx,mpid)
!       call write_mesh('grid_y_'//trim(tn(nav)),nav%gridy,mpid)
!       call write_mesh('grid_z_'//trim(tn(nav)),nav%gridz,mpid)
    endif

  end subroutine navier_write_fields

  function tn(nav)
! -----------------------------------------------------------------------
! navier : put time in character
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2013
!
    implicit none
    type(navier3d) :: nav
    character(100) :: tn

    write(tn,'(i8.8,a,i8.8)')floor(nav%time),'.',&
         int((nav%time-floor(nav%time))*1.d8)

  end function tn

  subroutine navier_bc_pressure(mpid,nav)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Alexandre Poux
! 01/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: i,l,m,c(3),inter(3,2)
    integer(ik) :: it(nav%nt),nt
    
    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infp,c,inter)
    else
      inter=0
    endif

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> m : left-right ; l : directions (x,y,z)
    do m=1,2
      do l=1,3
        if (inter(l,m)<=0) then
          select case (l)
            case (1)
              nav%bcphi(it(1))%bcx(:,:,m)=nav%rhs_px%f((m-1)*(nav%nx-1)+1,2:nav%ny-1,2:nav%nz-1)&
                                         -nav%fac(1)*nav%bcu(it(1))%bcx(:,:,m)
            case (2)
              nav%bcphi(it(1))%bcy(:,:,m)=nav%rhs_py%f(2:nav%nx-1,(m-1)*(nav%ny-1)+1,2:nav%nz-1)&
                                         -nav%fac(1)*nav%bcv(it(1))%bcy(:,:,m)
            case (3)
              nav%bcphi(it(1))%bcz(:,:,m)=nav%rhs_pz%f(2:nav%nx-1,2:nav%ny-1,(m-1)*(nav%nz-1)+1)&
                                         -nav%fac(1)*nav%bcw(it(1))%bcz(:,:,m)
          end select
        endif
      enddo
    enddo

   if(nav%pt==3)      call add_boundary_rotrot(mpid,nav)

!nav%bcphi(it(1))%bcx(:,:,:)=0._rk
!nav%bcphi(it(1))%bcy(:,:,:)=0._rk
!nav%bcphi(it(1))%bcz(:,:,:)=0._rk

  call erase_boundary_inter(inter,nav%bcphi(nav%it(1)))

  end subroutine navier_bc_pressure

  subroutine navier_bc_velocity(mpid,nav)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 11/2012
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: l,m,c(3),inter(3,2)
    
    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infu,c,inter)
    else
      inter=0
    endif

    call navier_bc_velocity_utils(inter,nav%bcu(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'u',&
         nav%nx,nav%ny,nav%nz,nav%rey)

    call navier_bc_velocity_utils(inter,nav%bcv(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'v',&
         nav%nx,nav%ny,nav%nz,nav%rey)

    call navier_bc_velocity_utils(inter,nav%bcw(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'w',&
         nav%nx,nav%ny,nav%nz,nav%rey)

    call erase_boundary_inter(inter,nav%bcu(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcv(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcw(nav%it(1)))

  end subroutine navier_bc_velocity

  subroutine erase_boundary_inter(inter,bc)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 11/2012
!
    implicit none
    integer(ik) :: l,m,inter(3,2)
    type(boundary_condition) :: bc
    
    !-> m : left-right ; l : directions (x,y,z)
    do m=1,2
       do l=1,3
          if (inter(l,m)>0) then
             if (l==1) bc%bcx(:,:,m)=0._rk
             if (l==2) bc%bcy(:,:,m)=0._rk
             if (l==3) bc%bcz(:,:,m)=0._rk
          endif
       enddo
    enddo

  end subroutine erase_boundary_inter

  subroutine navier_nullify_boundary(mpid,nav,x,i)
! -----------------------------------------------------------------------
! navier :
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    type(field) :: x
    integer(ik) :: l,m,c(3),inter(3,2),i
    integer(ik) :: it(nav%nt),nt


    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infp,c,inter)
    else
      inter=0
    endif

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

!    !-> m : left-right ; l : directions (x,y,z)
!    do m=1,2
!      do l=1,3
!        if ((i.ge.0.and.inter(l,m)>0).or.(i.le.0.and.inter(l,m)<=0)) then
!          select case (l)
!            case (1)
!              x%f((m-1)*(nav%nx-1)+1,:,:)=0._rk
!            case (2)
!              x%f(:,(m-1)*(nav%ny-1)+1,:)=0._rk
!            case (3)
!              x%f(:,:,(m-1)*(nav%nz-1)+1)=0._rk
!          end select
!        endif
!      enddo
!    enddo

  end subroutine navier_nullify_boundary

  subroutine add_boundary_gradient(mpid,nav)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt
    real(rk) :: fac1,fac2,fac3
    integer(ik) :: l,m,c(3),inter(3,2)
    
    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infu,c,inter)
    else
      inter=0
    endif

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    

    !-> bcx
!    nav%aux%f=0._rk
!    nav%aux=derx(nav%dcy,navier_extrapol(nav,nav%phi,type='p'))
!    nav%bcu(it(1))%bcx(:,:,1)=nav%bcu(it(1))%bcx(:,:,1)&
!                             +nav%aux%f(1,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)
!    nav%bcu(it(1))%bcx(:,:,2)=nav%bcu(it(1))%bcx(:,:,2)&
!                             +nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcv(it(1))%bcx(:,:,1)=nav%bcv(it(1))%bcx(:,:,1)&
                             +nav%aux%f(1,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)
    nav%bcv(it(1))%bcx(:,:,2)=nav%bcv(it(1))%bcx(:,:,2)&
                             +nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derz(nav%dcz,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcw(it(1))%bcx(:,:,1)=nav%bcw(it(1))%bcx(:,:,1)&
                             +nav%aux%f(1,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)
    nav%bcw(it(1))%bcx(:,:,2)=nav%bcw(it(1))%bcx(:,:,2)&
                             +nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)

    !-> bcy
!    nav%aux%f=0._rk
!    nav%aux=dery(nav%dcx,navier_extrapol(nav,nav%phi,type='p'))
!    nav%bcv(it(1))%bcy(:,:,1)=nav%bcv(it(1))%bcy(:,:,1)&
!                             +nav%aux%f(2:nav%nx-1,1,2:nav%nz-1)/nav%fac(1)
!    nav%bcv(it(1))%bcy(:,:,2)=nav%bcv(it(1))%bcy(:,:,2)&
!                             +nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcu(it(1))%bcy(:,:,1)=nav%bcu(it(1))%bcy(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,1,2:nav%nz-1)/nav%fac(1)
    nav%bcu(it(1))%bcy(:,:,2)=nav%bcu(it(1))%bcy(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derz(nav%dcz,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcw(it(1))%bcy(:,:,1)=nav%bcw(it(1))%bcy(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,1,2:nav%nz-1)/nav%fac(1)
    nav%bcw(it(1))%bcy(:,:,2)=nav%bcw(it(1))%bcy(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%fac(1)

    !-> bcz
!    nav%aux%f=0._rk
!    nav%aux=derz(nav%dcx,navier_extrapol(nav,nav%phi,type='p'))
!    nav%bcw(it(1))%bcz(:,:,1)=nav%bcw(it(1))%bcz(:,:,1)&
!                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,1)/nav%fac(1)
!    nav%bcw(it(1))%bcz(:,:,2)=nav%bcw(it(1))%bcz(:,:,2)&
!                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcu(it(1))%bcz(:,:,1)=nav%bcu(it(1))%bcz(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,1)/nav%fac(1)
    nav%bcu(it(1))%bcz(:,:,2)=nav%bcu(it(1))%bcz(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%phi,type='p'))
    nav%bcv(it(1))%bcz(:,:,1)=nav%bcv(it(1))%bcz(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,1)/nav%fac(1)
    nav%bcv(it(1))%bcz(:,:,2)=nav%bcv(it(1))%bcz(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%fac(1)

    call erase_boundary_inter(inter,nav%bcu(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcv(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcw(nav%it(1)))

  end subroutine add_boundary_gradient



  subroutine add_boundary_rotrot(mpid,nav)
! -----------------------------------------------------------------------
! navier :
! -----------------------------------------------------------------------
! Alexandre Poux
! 01/2013
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt,i,j,k
    real(rk) :: fac1,fac2,fac3,x,y,z,t
    integer(ik) :: l,m,c(3),inter(3,2)

    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infu,c,inter)
    else
      inter=0
    endif
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> bcx
    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%v,type='p'))&
           +derz(nav%dcz,navier_extrapol(nav,nav%w,type='p'))
    nav%aux=derx(nav%dcx,nav%aux)-ddery(nav%dcy,navier_extrapol(nav,nav%u,type='p'))&
                                 -dderz(nav%dcz,navier_extrapol(nav,nav%u,type='p'))
    nav%bcphi(it(1))%bcx(:,:,1)=nav%bcphi(it(1))%bcx(:,:,1)-nav%aux%f(1     ,2:nav%ny-1,2:nav%nz-1)/nav%rey
    nav%bcphi(it(1))%bcx(:,:,2)=nav%bcphi(it(1))%bcx(:,:,2)-nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%rey


    !-> bcy
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%u,type='p'))&
           +derz(nav%dcz,navier_extrapol(nav,nav%w,type='p'))
    nav%aux=dery(nav%dcy,nav%aux)-dderx(nav%dcx,navier_extrapol(nav,nav%v,type='p'))&
                                 -dderz(nav%dcz,navier_extrapol(nav,nav%v,type='p'))
    nav%bcphi(it(1))%bcy(:,:,1)=nav%bcphi(it(1))%bcy(:,:,1)-nav%aux%f(2:nav%nx-1,1     ,2:nav%nz-1)/nav%rey
    nav%bcphi(it(1))%bcy(:,:,2)=nav%bcphi(it(1))%bcy(:,:,2)-nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%rey

    !-> bcz
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%u,type='p'))&
           +dery(nav%dcy,navier_extrapol(nav,nav%v,type='p'))
    nav%aux=derz(nav%dcz,nav%aux)-dderx(nav%dcx,navier_extrapol(nav,nav%w,type='p'))&
                                 -ddery(nav%dcy,navier_extrapol(nav,nav%w,type='p'))
    nav%bcphi(it(1))%bcz(:,:,1)=nav%bcphi(it(1))%bcz(:,:,1)-nav%aux%f(2:nav%nx-1,2:nav%ny-1,1     )/nav%rey
    nav%bcphi(it(1))%bcz(:,:,2)=nav%bcphi(it(1))%bcz(:,:,2)-nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%rey

  end subroutine add_boundary_rotrot
  
  subroutine navier_bc_velocity_utils(inter,bc,gridx,gridy,gridz,t,&
       var,nx,ny,nz,rey)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    integer(ik) :: l,m,inter(3,2)
    type(boundary_condition) :: bc
    type(mesh_grid) :: gridx,gridy,gridz
    real(rk) :: x,y,z,t,rey
    integer(ik) :: i,j,k,nx,ny,nz
    character(*) :: var

    !-> boundary condition
    !-> x-direction
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
    do k=2,nz-1
       do j=2,ny-1
          y=gridy%grid1d(j)
          z=gridz%grid1d(k)
          
          x=gridx%grid1d(1)
          bc%bcx(j-1,k-1,1)=sol(x,y,z,t,var,rey)
!          bc%bcx(j-1,k-1,1)=0._rk
          x=gridx%grid1d(nx)
          bc%bcx(j-1,k-1,2)=sol(x,y,z,t,var,rey)
!          bc%bcx(j-1,k-1,2)=0._rk
       enddo
    enddo
!$OMP END PARALLEL DO
    !print*,x,y,z,t
    !-> y-direction
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
    do k=2,nz-1
       do i=2,nx-1
          x=gridx%grid1d(i)
          z=gridz%grid1d(k)

          y=gridy%grid1d(1)
          bc%bcy(i-1,k-1,1)=sol(x,y,z,t,var,rey)
!          bc%bcy(i-1,k-1,1)=0._rk
          y=gridy%grid1d(ny)
          bc%bcy(i-1,k-1,2)=sol(x,y,z,t,var,rey)
!          if (var=='u') then
!             bc%bcy(i-1,k-1,2)=(1._rk-x**2)
!          else
!             bc%bcy(i-1,k-1,2)=0._rk
!          endif
       enddo
    enddo
!$OMP END PARALLEL DO
    !-> z-direction
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
    do j=2,ny-1
       do i=2,nx-1
          x=gridx%grid1d(i)
          y=gridy%grid1d(j)
          
          z=gridz%grid1d(1)
          bc%bcz(i-1,j-1,1)=sol(x,y,z,t,var,rey)
!          bc%bcz(i-1,j-1,1)=0._rk
          z=gridz%grid1d(nz)
          bc%bcz(i-1,j-1,2)=sol(x,y,z,t,var,rey)
!          bc%bcz(i-1,j-1,2)=0._rk
       enddo
    enddo
!$OMP END PARALLEL DO

  end subroutine navier_bc_velocity_utils

  function sol(x,y,z,t,type,rey)
! -----------------------------------------------------------------------
! exact solution : function, derivatives and rhs
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    real(rk) :: sol,rey
    integer(ik) ::i,n,j,k
    PARAMETER(n=5)     !espace vitesse
!    PARAMETER(i=0)     !espace pression
!    PARAMETER(j=0)     !temps vitesse
!    PARAMETER(k=0)     !temps pression
    real(rk) :: x,y,z,t,pi,wx(n),wt(n),wy(n),wp(n),a(n),b(n)
    character(*) :: type

    sol=0._rk

!    if (type=="u") then
!       if(n>0) sol= (n*(x**n)+n*n*x*y**(n-1))*t**j
!    endif
!    if (type=="v") then
!       if(n>0) sol=-(n*(y**n)+n*n*y*x**(n-1))*t**j
!    endif
!    if (type=="w") then
!       sol=0._rk
!    endif
!    if (type=="p") then
!       if(i>0) sol=(i*(x**i+y**i+x*y))*cos(k*t)
!    endif
!    if (type=="dxp") then
!       if(i>0) sol=i*cos(k*t)*(y+i*x**(i-1))
!    endif
!    if (type=="dyp") then
!       if(i>0) sol=i*cos(k*t)*(i*y**(i-1)+x)
!    endif
!    if (type=="dzp") then
!       sol=0._rk
!    endif

!    if (type=="rhsu") then
!       if(n>0) then
!         if(j>0) sol=sol+n*j*(n*x*(y**(n-1))+x**n)*(t**(j-1))    !temps
!         if(n>2) sol=sol-(n-1)*(n**2)*(t**j)*(n-2)*x*(y**(n-3))/rey  !lap
!         if(n>1) sol=sol-(n-1)*(n**2)*(t**j)*(x**(n-2))/rey  !lap
!       endif
!       if(i>0) sol=sol+i*cos(k*t)*(y+i*(x**(i-1)))                 !grad
!       if(n>0) sol=sol+n**3*t**(2*j)*x*y**(2*n-2)&
!               -n**5*t**(2*j)*x**n*y**(n-1)&
!               +2*n**4*t**(2*j)*x**n*y**(n-1)&
!               +n**3*t**(2*j)*x**n*y**(n-1)&
!               +n**3*t**(2*j)*x**(2*n-1) !nonlinear
!    endif
!    if (type=="rhsv") then
!       if(n>0) then
!         if(j>0) sol=sol-n*j*((y**n)+n*(x**(n-1))*y)*(t**(j-1))    !temps
!         if(n>2) sol=sol+(n-2)*(n-1)*(n**2)*(t**j)*y*(x**(n-3))/rey  !lap
!         if(n>1) sol=sol+(n-1)*(n**2)*(t**j)*(y**(n-2))/rey  !lap
!       endif
!       if(i>0) sol=sol+i*cos(k*t)*(i*(y**(i-1))+x)                 !grad
!       if(n>0) sol=sol+n**3*t**(2*j)*y**(2*n-1)&
!               -n**5*t**(2*j)*x**(n-1)*y**n&
!               +2*n**4*t**(2*j)*x**(n-1)*y**n&
!               +n**3*t**(2*j)*x**(n-1)*y**n&
!               +n**3*t**(2*j)*x**(2*n-2)*y !nonlinear
!    endif
!    if (type=="rhsw") then
!       sol=0._rk
!    endif
!    if (type=="rhsp") then
!       sol=0._rk
!    endif

    pi=4._rk*atan(1._rk)
    do i=1,n
      a(i)=1._rk ; b(i)=1._rk
      wx(i)=1._rk*pi*i ; wy(i)=1._rk*pi*i ; wp(i)=1._rk*pi*i ; wt(i)=0._rk*i*pi

    select case (type)
    case ("u")
       sol=sol+a(i)*sin(wx(i)*x)*cos(wy(i)*y)*cos(wt(i)*t)/wx(i)
    case ("v")
       sol=sol-a(i)*cos(wx(i)*x)*sin(wy(i)*y)*cos(wt(i)*t)/wy(i)
    case ("w")
       sol=sol+0._rk
    case ("p")
       sol=sol+b(i)*sin(wp(i)*(x-y))*cos(wt(i)*t)
    case ("dxp")
       sol=sol+b(i)*wp(i)*cos(t*wt(i))*cos(wp(i)*(x-y))
    case ("dyp") 
       sol=sol-b(i)*wp(i)*cos(t*wt(i))*cos(wp(i)*(x-y))
    case ("dzp")
       sol=sol+0._rk
    case ("rhsu")
sol=sol+b(i)*wp(i)*cos(wt(i)*t)*cos(wp(i)*(y-x)) & !grap p
+a(i)**2*cos(wt(i)*t)**2*cos(wx(i)*x)*sin(wx(i)*x)*(sin(wy(i)*y)**2+cos(wy(i)*y)**2)/wx(i) & !nonlinear
+a(i)*cos(wt(i)*t)*sin(wx(i)*x)*cos(wy(i)*y)*(wy(i)**2/wx(i)+wx(i))/rey & !lap
-a(i)*wt(i)*sin(wt(i)*t)*sin(wx(i)*x)*cos(wy(i)*y)/wx(i)  !temps
    case ("rhsv")
sol=sol-b(i)*wp(i)*cos(wt(i)*t)*cos(wp(i)*(y-x)) & !grap p
+a(i)**2*cos(wt(i)*t)**2*cos(wy(i)*y)*sin(wy(i)*y)*(sin(wx(i)*x)**2+cos(wx(i)*x)**2)/wy(i)& !nonlinear
-a(i)*cos(wt(i)*t)*cos(wx(i)*x)*sin(wy(i)*y)*(wy(i)+wx(i)**2/wy(i))/rey& !lap
+a(i)*wt(i)*sin(wt(i)*t)*cos(wx(i)*x)*sin(wy(i)*y)/wy(i)  !temps
    case ("rhsw")
       sol=sol+0._rk
    case ("rhsp")
       sol=sol+0._rk
    case default
       sol=sol+0._rk
    end select
    end do


  end function sol

  function f(nav,mpid,var)
! -----------------------------------------------------------------------
! field : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
!$ use OMP_LIB
    implicit none
    type(field) :: f
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: i,j,k,nt,it(nav%nt)
    real(rk) :: x,y,z,t
    character(*) :: var
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    call field_init(f,"F",nav%nx,nav%ny,nav%nz)
    f%f=0._rk

    t=nav%time
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
    do k=1,nav%nz
       do j=1,nav%ny
          do i=1,nav%nx
             x=nav%gridx%grid1d(i)
             y=nav%gridy%grid1d(j)
             z=nav%gridz%grid1d(k)
             f%f(i,j,k)=sol(x,y,z,t,var,nav%rey)
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

if(nav%les_type.ne.0) then
    nav%les_f(it(nt))%f=0._rk    

    if(var=='rhsu') i=1
    if(var=='rhsv') i=2
    if(var=='rhsw') i=3
  
    if(.true.) then 
    do j=1,3
      !-> tau = nu * 2*S
      nav%aux =nav%les_nu*(nav%du(i,j)+nav%du(j,i))

      !->f=-div(tau)
      if(j==1) nav%les_f(it(nt))=nav%les_f(it(nt))-derx(nav%dcx,nav%aux)
      if(j==2) nav%les_f(it(nt))=nav%les_f(it(nt))-dery(nav%dcy,nav%aux)
      if(j==3) nav%les_f(it(nt))=nav%les_f(it(nt))-derz(nav%dcz,nav%aux)

    enddo

else  !->div(tau)=2S*grad(nu) + 2nu*lap(u) 

!    nav%les_f(it(nt))=nav%les_f(it(nt))-derx(nav%dcx,nav%les_nu)*(nav%du(i,1)+nav%du(1,i))
!    nav%les_f(it(nt))=nav%les_f(it(nt))-dery(nav%dcy,nav%les_nu)*(nav%du(i,2)+nav%du(2,i))
!    nav%les_f(it(nt))=nav%les_f(it(nt))-derz(nav%dcz,nav%les_nu)*(nav%du(i,3)+nav%du(3,i))

if(var=='rhsu') then
      nav%aux=nav%u(it(nt))
!      nav%aux=navier_extrapol(nav,nav%u,type='v')
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderx(nav%dcx,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*ddery(nav%dcy,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderz(nav%dcz,nav%aux)
    endif
    if(var=='rhsv') then
      nav%aux=nav%v(it(nt))
!      nav%aux=navier_extrapol(nav,nav%v,type='v')
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderx(nav%dcx,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*ddery(nav%dcy,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderz(nav%dcz,nav%aux)
    endif
    if(var=='rhsw') then
      nav%aux=nav%w(it(nt))
!      nav%aux=navier_extrapol(nav,nav%w,type='v')
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderx(nav%dcx,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*ddery(nav%dcy,nav%aux)
      nav%les_f(it(nt))=nav%les_f(it(nt))-nav%les_nu*2._rk*dderz(nav%dcz,nav%aux)
    endif

endif

!    f=-navier_extrapol(nav,nav%les_f,type='v')
    f=f-nav%les_f(it(nt))
endif

  end function f

 subroutine navier_nonlinear(mpid,nav,x,f)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 11/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt,i
    type(field) :: x(nav%nt),f
    type(field),dimension(:),allocatable,save ::tmp

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    if(.not.allocated(tmp)) then
      allocate(tmp(nav%nt))
      do i=1,nav%nt
         call field_init(tmp(i),"NLT",nav%nx,nav%ny,nav%nz)
      enddo
    endif

    !-> nonlinear terms
    if (nav%nlt==1) then
      do i=1,nav%nt
          tmp(i)=(nav%u(i)*derx(nav%dcx,x(i))+&
                  nav%v(i)*dery(nav%dcy,x(i))+&
                  nav%w(i)*derz(nav%dcz,x(i)))
      enddo
    elseif (nav%nlt==2) then
      do i=1,nav%nt
          tmp(i)=(nav%u(i)*derx(nav%dcx,x(i))+&
                  nav%v(i)*dery(nav%dcy,x(i))+&
                  nav%w(i)*derz(nav%dcz,x(i)))*0.5_rk

          nav%aux=x(i)*nav%u(i)
          tmp(i)=tmp(i)+derx(nav%dcx,nav%aux)*0.5_rk
          nav%aux=x(i)*nav%v(i)
          tmp(i)=tmp(i)+dery(nav%dcy,nav%aux)*0.5_rk
          nav%aux=x(i)*nav%w(i)
          tmp(i)=tmp(i)+derz(nav%dcz,nav%aux)*0.5_rk

      enddo
    else
      do i=1,nav%nt
        tmp(i)%f=0._rk
      enddo
    endif
    f=f+navier_extrapol(nav,tmp,type='v')

  end subroutine navier_nonlinear


  subroutine navier_LES(mpid,nav)
! -----------------------------------------------------------------------
! navier : compute les_nu 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    real(rk)       :: dx,dy,dz,g2(3,3)
    real(rk)       :: sigma1,sigma2,sigma3
    integer(ik)    :: i,j,k,l,m
    integer(ik) :: it(nav%nt),nt

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)


nav%aux=nav%u(it(nt))
!nav%aux=navier_extrapol(nav,nav%u,type='v')
nav%du(1,1)=derx(nav%dcx,nav%aux)
nav%du(1,2)=dery(nav%dcy,nav%aux)
nav%du(1,3)=derz(nav%dcz,nav%aux)

nav%aux=nav%v(it(nt))
!nav%aux=navier_extrapol(nav,nav%v,type='v')
nav%du(2,1)=derx(nav%dcx,nav%aux)
nav%du(2,2)=dery(nav%dcy,nav%aux)
nav%du(2,3)=derz(nav%dcz,nav%aux)

nav%aux=nav%w(it(nt))
!nav%aux=navier_extrapol(nav,nav%w,type='v')
nav%du(3,1)=derx(nav%dcx,nav%aux)
nav%du(3,2)=dery(nav%dcy,nav%aux)
nav%du(3,3)=derz(nav%dcz,nav%aux)

do i=1,3
do j=1,3
call field_zero_edges(nav%du(i,j))
enddo
enddo

nav%les_nu%f=0._rk

if(nav%les_type.ne.0)then
if(nav%les_type.ne.3)then
!-> S:S
nav%aux%f=0._rk
do i=1,3
   nav%aux=nav%aux+nav%du(i,i)**2
do j=i+1,3
   nav%aux=nav%aux+0.5_rk*(nav%du(i,j)+nav%du(j,i))**2
enddo
enddo

!nav%les_c=0.165_rk**2
if(nav%les_type.eq.1) then    !-> Smagorinski
    nav%les_nu%f=nav%les_c*sqrt(2._rk*nav%aux%f)
!print*,'smag',mpid%rank,minval(nav%les_nu%f),maxval(nav%les_nu%f),sum(nav%les_nu%f)


elseif(nav%les_type.eq.2) then    !-> WALE
!!-> Omega:Omega
!nav%les_nu%f=0._rk
!do i=1,3
!do j=i+1,3
!   nav%les_nu=nav%les_nu+0.5_rk*(nav%du(i,j)-nav%du(j,i))**2
!    enddo
!enddo
!
!!-> Sd:Sd
!nav%les_nu=(1._rk/6._rk) * (nav%aux*nav%aux + nav%les_nu*nav%les_nu) &
!          +(2._rk/3._rk) * (nav%aux*nav%les_nu)
!
!do i=1,3
!do j=1,3
!do k=1,3
!do l=1,3
!   nav%les_nu=nav%les_nu+0.125_rk*( (nav%du(i,k)+nav%du(k,i))  &
!                                *(nav%du(k,j)+nav%du(j,k))  &
!                                *(nav%du(j,l)-nav%du(l,j))  &
!                                *(nav%du(l,i)-nav%du(i,l)))/2._rk 
!enddo
!enddo
!enddo
!enddo
!nav%les_nu%f=max(nav%les_nu%f,0._rk) ! round-off
!print*,'1',mpid%rank,minval(nav%les_nu%f),maxval(nav%les_nu%f)
!
!!-> Omega:Omega
!nav%les_nu%f=0._rk  
!do i=1,3
!do j=i+1,3
!   nav%les_nu=nav%les_nu+0.5_rk*(nav%du(i,j)-nav%du(j,i))**2
!enddo
!enddo
!
!nav%aux2%f=0._rk
!do i=1,3
!do j=1,3
!
!  nav%aux1%f=0._rk
!  do k=1,3
!   nav%aux1=nav%aux1+0.25_rk*( (nav%du(i,k)+nav%du(k,i))  &
!                              *(nav%du(k,j)+nav%du(j,k))  )
!
!   nav%aux1=nav%aux1+0.25_rk*( (nav%du(i,k)-nav%du(k,i))  &
!                              *(nav%du(k,j)-nav%du(j,k)) ) 
!  enddo
!
!  if(i==j) nav%aux1=nav%aux1 - (nav%aux - nav%les_nu)/3._rk
!  nav%aux2= nav%aux2 + nav%aux1**2
!
!enddo
!enddo
!
!  nav%les_nu= nav%aux2
!
!
!
!print*,'2',mpid%rank,minval(nav%les_nu%f),maxval(nav%les_nu%f)

nav%les_nu%f=0._rk
do i=1,3
do j=1,3

  nav%aux1%f=0._rk
  do k=1,3
    nav%aux1=nav%aux1+0.5_rk*(nav%du(i,k)*nav%du(k,j) + nav%du(j,k)*nav%du(k,i) )
    if(i==j)   nav%aux1=nav%aux1  - (nav%du(1,k)*nav%du(k,1)+ &
                                     nav%du(2,k)*nav%du(k,2)+ & 
                                     nav%du(3,k)*nav%du(k,3))/3._rk
  enddo
  nav%les_nu= nav%les_nu + nav%aux1**2

enddo
enddo

!nav%les_c=0.5_rk**2
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k) &
!$OMP SCHEDULE(RUNTIME)
    do k=1,nav%nz
    do j=1,nav%ny
    do i=1,nav%nx   
    if (nav%les_nu%f(i,j,k).lt.1e-10) then
       nav%les_nu%f(i,j,k)=0._rk
    else
    !-> WALE
    nav%les_nu%f(i,j,k)=          nav%les_c*nav%les_nu%f(i,j,k)**(3._rk/2._rk)          & 
                  /(nav%aux%f(i,j,k)**(5._rk/2._rk) + nav%les_nu%f(i,j,k)**(5._rk/4._rk))
    endif
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
!print*,'wale',mpid%rank,minval(nav%les_nu%f),maxval(nav%les_nu%f),sum(nav%les_nu%f)

endif
elseif(nav%les_type.eq.3)then ! sigma

!nav%les_c=1.35_rk**2
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,g2,sigma1,sigma2,sigma3) &
!$OMP SCHEDULE(RUNTIME)
    do k=1,nav%nz
       do j=1,nav%ny
          do i=1,nav%nx

            do l=1,3
            do m=1,3
              g2(l,m)=nav%du(l,m)%f(i,j,k)
            enddo 
            enddo
            call les_sigma(g2,sigma1,sigma2,sigma3)
            if(sigma1.lt.1e-10) then
                nav%les_nu%f(i,j,k)=0._rk
            else
                nav%les_nu%f(i,j,k)=nav%les_c*sigma3*(sigma1-sigma2)*(sigma2-sigma3)/(sigma1**2)
            endif
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO
!print*,'sigm',mpid%rank,minval(nav%les_nu%f),maxval(nav%les_nu%f),sum(nav%les_nu%f)

endif

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP SCHEDULE(RUNTIME)
    do k=1,nav%nz
       if(k.eq.1) then
         dz= nav%gridz%grid1d(k+1)- nav%gridz%grid1d(k)
       elseif(k.eq.nav%nz) then
         dz= nav%gridz%grid1d(k  )- nav%gridz%grid1d(k-1)
       else
         dz=(nav%gridz%grid1d(k+1)- nav%gridz%grid1d(k-1))*0.5_rk
       endif
       do j=1,nav%ny
          if(j.eq.1) then
            dy= nav%gridy%grid1d(j+1)- nav%gridy%grid1d(j)
          elseif(j.eq.nav%ny) then
            dy= nav%gridy%grid1d(j  )- nav%gridy%grid1d(j-1)
          else
            dy=(nav%gridy%grid1d(j+1)- nav%gridy%grid1d(j-1))*0.5_rk
          endif
          do i=1,nav%nx
            if(i.eq.1) then
              dx= nav%gridx%grid1d(i+1)- nav%gridx%grid1d(i)
            elseif(i.eq.nav%nx) then
              dx= nav%gridx%grid1d(i  )- nav%gridx%grid1d(i-1)
            else
              dx=(nav%gridx%grid1d(i+1)- nav%gridx%grid1d(i-1))*0.5_rk
            endif

            !-> delta ** 2
!            nav%les_nu%f(i,j,k)=nav%les_nu%f(i,j,k) * ((dx*dy*dz)**(2._rk/3._rk))
            nav%les_nu%f(i,j,k)=nav%les_nu%f(i,j,k) * (dx**2 + dy**2 + dz**2)/3._rk

          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

call field_zero_edges(nav%les_nu)

endif
  end subroutine navier_LES


  subroutine les_sigma(g2,sigma1,sigma2,sigma3)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    real(rk)       :: g(3,3),g2(3,3)
    real(rk)       :: sigma1,sigma2,sigma3
    real(rk)       :: alpha1,alpha2,alpha3
    real(rk)       :: i1,i2,i3,pi
    integer(ik)    :: i,j
    
    pi=4._rk*atan(1._rk)

g=matmul(g2,transpose(g2))

    i1=g(1,1)+g(2,2)+g(3,3)    ! trace

    i2=g(1,1)*g(2,2) + g(1,1)*g(3,3) + g(2,2)*g(3,3) &    !( tr(g)^2 -tr(g^2) )/ 2&
      -g(1,2)*g(2,1) - g(1,3)*g(3,1) - g(2,3)*g(3,2)

    i3=g(1,1)*g(2,2)*g(3,3) + g(2,1)*g(3,2)*g(1,3) &  ! determinant
      +g(3,1)*g(1,2)*g(2,3) - g(1,3)*g(2,2)*g(3,1) &
      -g(2,3)*g(3,2)*g(1,1) - g(3,3)*g(1,2)*g(2,1) 

    alpha1=(i1**2)/9._rk - i2/3._rk
    alpha2=(i1**3)/27._rk - i1*i2/6._rk + i3/2._rk
    if(alpha1.lt.1e-10) then ! alors tout est nul
        alpha3=0._rk
    else
        alpha3=acos(alpha2/(alpha1**(3._rk/2._rk)))/3._rk
    endif

    sigma1=(i1/3._rk + 2._rk*sqrt(alpha1)*cos(alpha3))
    sigma2=(i1/3._rk - 2._rk*sqrt(alpha1)*cos(alpha3+pi/3._rk))
    sigma3=(i1-sigma1-sigma2) !(i1/3._rk - 2._rk*sqrt(alpha1)*cos(alpha3-pi/3._rk))

    sigma1=sqrt(max(sigma1,0._rk)) ! round-off
    sigma2=sqrt(max(sigma2,0._rk)) ! round-off
    sigma3=sqrt(max(sigma3,0._rk)) ! round-off

  end subroutine les_sigma
  subroutine navier_presolve_u(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid

    !--------------------------------------------------------------------
    !-> compute rhs

    nav%rhs_u(2)=nav%rhs_u(1)
    nav%rhs_u(1)%f=0._rk
    call navier_nonlinear(mpid,nav,nav%u,nav%rhs_u(1))
    nav%rhs_u(1)=f(nav,mpid,'rhsu') - nav%rhs_u(1)

    call field_zero_edges(nav%rhs_u(1))

  end subroutine navier_presolve_u

  subroutine navier_presolve_v(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid

    !--------------------------------------------------------------------
    !-> compute rhs

    nav%rhs_v(2)=nav%rhs_v(1)
    nav%rhs_v(1)%f=0._rk
    call navier_nonlinear(mpid,nav,nav%v,nav%rhs_v(1))
    nav%rhs_v(1)=f(nav,mpid,'rhsv') - nav%rhs_v(1)

    call field_zero_edges(nav%rhs_v(1))

  end subroutine navier_presolve_v

  subroutine navier_presolve_w(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid

    !--------------------------------------------------------------------
    !-> compute rhs

    nav%rhs_w(2)=nav%rhs_w(1)
    nav%rhs_w(1)%f=0._rk
    call navier_nonlinear(mpid,nav,nav%w,nav%rhs_w(1))
    nav%rhs_w(1)=f(nav,mpid,'rhsw') -  nav%rhs_w(1)

    call field_zero_edges(nav%rhs_w(1))

  end subroutine navier_presolve_w

  subroutine navier_solve_u(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fu(it(nt))%f=0._rk
    
    !-> time
    nav%fu(it(nt))=nav%fu(it(nt))+nav%fac(2)*nav%u(nav%it(nav%nt  )) &
                                 +nav%fac(3)*nav%u(nav%it(nav%nt-1))
    if(nav%tou==3) nav%fu(it(nt))=nav%fu(it(nt))+nav%fac(4)*nav%u(nav%it(nav%nt-2))
    
    !-> pressure 
    if (nav%pt>=2) then
       if (nav%pt==2) then
          nav%fu(it(nt))=nav%fu(it(nt))+derx(nav%dcx,navier_extrapol(nav,nav%p,type='p'))
       else
          nav%fu(it(nt))=nav%fu(it(nt))+derx(nav%dcx,nav%p(it(1)))
       endif
    endif
    
    !-> nonlinear terms and function
    nav%fu(it(nt))=nav%fu(it(nt))-nav%rhs_u(1)

    !-> reynolds number multiplication
    nav%fu(it(nt))=nav%rey*nav%fu(it(nt))

    
    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infu,nt,it,nav%bcu,nav%u)
    if (mpid%dims.ne.0) then
      call multidomain_solve(mpid,nav%infu,nav%scu,nav%bcu(it(1)),nav%u(it(1)),&
          nav%fu(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolu)
    else
      call solver_3d(nav%scu,nav%fu(it(nt)),nav%u(it(1)),nav%bcu(it(1)),nav%sigmau)
    endif


  end subroutine navier_solve_u

  subroutine navier_solve_v(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve v helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fv(it(nt))%f=0._rk
    
    !-> time
    nav%fv(it(nt))=nav%fv(it(nt))+nav%fac(2)*nav%v(nav%it(nav%nt  )) &
                                 +nav%fac(3)*nav%v(nav%it(nav%nt-1))
    if(nav%tou==3) nav%fv(it(nt))=nav%fv(it(nt))+nav%fac(4)*nav%v(nav%it(nav%nt-2))
    
    !-> pressure 
    if (nav%pt>=2) then
       if (nav%pt==2) then
          nav%fv(it(nt))=nav%fv(it(nt))+dery(nav%dcy,navier_extrapol(nav,nav%p,type='p'))
       else
          nav%fv(it(nt))=nav%fv(it(nt))+dery(nav%dcy,nav%p(it(1)))
       endif
    endif
    
    !-> nonlinear terms and function
    nav%fv(it(nt))=nav%fv(it(nt))-nav%rhs_v(1)

    !-> reynolds number multiplication
    nav%fv(it(nt))=nav%rey*nav%fv(it(nt))

    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infv,nt,it,nav%bcv,nav%v)
    if (mpid%dims.ne.0) then
      call multidomain_solve(mpid,nav%infv,nav%scv,nav%bcv(it(1)),nav%v(it(1)),&
          nav%fv(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolv)
    else
      call solver_3d(nav%scv,nav%fv(it(nt)),nav%v(it(1)),nav%bcv(it(1)),nav%sigmau)
    endif

  end subroutine navier_solve_v

  subroutine navier_solve_w(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve w helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fw(it(nt))%f=0._rk
    
    !-> time
    nav%fw(it(nt))=nav%fw(it(nt))+nav%fac(2)*nav%w(nav%it(nav%nt  )) &
                                 +nav%fac(3)*nav%w(nav%it(nav%nt-1))
    if(nav%tou==3) nav%fw(it(nt))=nav%fw(it(nt))+nav%fac(4)*nav%w(nav%it(nav%nt-2))
    
    !-> pressure 
    if (nav%pt>=2) then
       if (nav%pt==2) then
          nav%fw(it(nt))=nav%fw(it(nt))+derz(nav%dcz,navier_extrapol(nav,nav%p,type='p'))
       else
          nav%fw(it(nt))=nav%fw(it(nt))+derz(nav%dcz,nav%p(it(1)))
       endif
    endif
    
    !-> nonlinear terms and function
    nav%fw(it(nt))=nav%fw(it(nt))-nav%rhs_w(1)

    !-> reynolds number multiplication
    nav%fw(it(nt))=nav%rey*nav%fw(it(nt))

    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infw,nt,it,nav%bcw,nav%w)
    if (mpid%dims.ne.0) then
      call multidomain_solve(mpid,nav%infw,nav%scw,nav%bcw(it(1)),nav%w(it(1)),&
          nav%fw(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolw)
    else
      call solver_3d(nav%scw,nav%fw(it(nt)),nav%w(it(1)),nav%bcw(it(1)),nav%sigmau)
    endif

  end subroutine navier_solve_w


  subroutine navier_presolve_phi(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Alexandre Poux
! 01/2013
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid

    !--------------------------------------------------------------------
    !-> compute rhs
  nav%rhs_px=navier_phi_rhs(mpid,nav,nav%u,nav%rhs_u)
  nav%rhs_py=navier_phi_rhs(mpid,nav,nav%v,nav%rhs_v)
  nav%rhs_pz=navier_phi_rhs(mpid,nav,nav%w,nav%rhs_w)

  call field_zero_edges(nav%rhs_px)
  call field_zero_edges(nav%rhs_py)
  call field_zero_edges(nav%rhs_pz)

  end subroutine navier_presolve_phi

  subroutine navier_solve_phi(mpid,nav)
! -----------------------------------------------------------------------
! navier : solve phi helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt,i,j,k
    real(rk) :: fac
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fphi=nav%sigmap*navier_extrapol(nav,nav%phi,type='p')
!    nav%phi(it(1))%f=0._rk


    !-> 
    nav%fphi=nav%fphi+ derx(nav%dcx,nav%rhs_px)&
                     + dery(nav%dcy,nav%rhs_py)&
                     + derz(nav%dcz,nav%rhs_pz)

    goto 101
    nav%fphi=(1.5_rk/nav%ts)*(&
         derx(nav%dcx,nav%u(it(1)))+&
         dery(nav%dcy,nav%v(it(1)))+&
         derz(nav%dcz,nav%w(it(1))))

    nav%fphi=nav%fphi+(-2._rk/nav%ts)*(&
         derx(nav%dcx,nav%u(it(nt)))+&
         dery(nav%dcy,nav%v(it(nt)))+&
         derz(nav%dcz,nav%w(it(nt))))
    
    nav%fphi=nav%fphi+(0.5_rk/nav%ts)*(&
         derx(nav%dcx,nav%u(it(nt-1)))+&
         dery(nav%dcy,nav%v(it(nt-1)))+&
         derz(nav%dcz,nav%w(it(nt-1))))
101 continue

    !-> function
!    nav%fphi=nav%fphi+f(nav,'rhsp')
    
    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infp,nt,it,nav%bcphi,nav%phi)
    if (mpid%dims.ne.0) then
     call multidomain_solve(mpid,nav%infp,nav%scp,nav%bcphi(it(1)),nav%phi(it(1)),&
          nav%fphi,nav%aux,nav%sigmap,nav%dcx,nav%dcy,nav%dcz,null=nullv,&
          inf_sol=nav%infsolphi,var='p')
    else
      call solver_3d(nav%scp,nav%fphi,nav%phi(it(1)),nav%bcphi(it(1)),nav%sigmap)
    endif

  end subroutine navier_solve_phi

  subroutine navier_projection(mpid,nav)
! -----------------------------------------------------------------------
! navier : compute divergence null u,v,w,p
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    real(rk) :: fac
    integer(ik) :: iaux,i,ex(3,2)
    integer(ik) :: it(nav%nt),nt
    integer(ik) :: l,m,c(3),inter(3,2)
    
    !-> get interface type
    if (mpid%dims.ne.0) then
      call md_mpi_getcoord(mpid,c)
      call md_get_interfaces_number(nav%infu,c,inter)
    else
      inter=0
    endif

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> compute extrema
    ex(1,1)=2 ; ex(2,1)=2 ; ex(3,1)=2 
    ex(1,2)=nav%nx-1 ; ex(2,2)=nav%ny-1 ; ex(3,2)=nav%nz-1

    call field_zero_edges(nav%phi(it(1)))
    !-> pressure
    if(nav%pt==2.or.nav%pt==4) nav%aux=navier_extrapol(nav,nav%p,type='p')
    nav%p(it(1))=nav%phi(it(1))
    if(nav%pt==2.or.nav%pt==4) nav%p(it(1))=nav%p(it(1)) + nav%aux

    call field_zero_edges(nav%p(it(1)))

    !-> rotationnal
    if(nav%pt<=2) then
      nav%aux=(derx(nav%dcx,nav%u(it(1)))+&
               dery(nav%dcy,nav%v(it(1)))+&
               derz(nav%dcz,nav%w(it(1))))/nav%rey
      nav%p(it(1))=nav%p(it(1))-nav%aux
    elseif(nav%pt==4) then
      nav%aux=(derx(nav%dcx,navier_extrapol(nav,nav%u,type='p'))+&
               dery(nav%dcy,navier_extrapol(nav,nav%v,type='p'))+&
               derz(nav%dcz,navier_extrapol(nav,nav%w,type='p')))/nav%rey
      nav%p(it(1))=nav%p(it(1))-nav%aux
    endif
!    nav%p(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
!         nav%p(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
!         -nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))

 
    !-> Brown
!    nav%aux=fac*(dderx(nav%dcx,nav%phi(it(1)))+&
!         ddery(nav%dcy,nav%phi(it(1)))+&
!         dderz(nav%dcz,nav%phi(it(1))))/nav%rey
!    nav%p(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
!         nav%p(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
!         -nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))
!    nav%p(it(1))=nav%p(it(1))-nav%aux

    call field_zero_edges(nav%p(it(1)))
   

!    goto 102
    do m=1,2
       do l=1,3
          if (inter(l,m)>0) then
             if (l==1.and.m==1) ex(l,m)=1
             if (l==1.and.m==2) ex(l,m)=nav%nx
             if (l==2.and.m==1) ex(l,m)=1
             if (l==2.and.m==2) ex(l,m)=nav%ny
             if (l==3.and.m==1) ex(l,m)=1
             if (l==3.and.m==2) ex(l,m)=nav%nz
          endif
       enddo
    enddo
!102 continue


!    call navier_bc_velocity(mpid,nav)
!    call field_put_boundary(nav%u(it(1)),nav%bcu(it(1)),inter)
!    call field_put_boundary(nav%v(it(1)),nav%bcv(it(1)),inter)
!    call field_put_boundary(nav%w(it(1)),nav%bcw(it(1)),inter)

    !-> velocity
    goto 101
    nav%aux=derx(nav%dcx,nav%phi(it(1)))
    nav%u(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%u(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))

    nav%aux=dery(nav%dcy,nav%phi(it(1)))
    nav%v(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%v(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))

    nav%aux=derz(nav%dcz,nav%phi(it(1)))
    nav%w(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%w(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))
101 continue

    nav%u(it(1))=(nav%rhs_px-derx(nav%dcx,nav%phi(it(1))))/nav%fac(1)
    nav%v(it(1))=(nav%rhs_py-dery(nav%dcy,nav%phi(it(1))))/nav%fac(1)
    nav%w(it(1))=(nav%rhs_pz-derz(nav%dcz,nav%phi(it(1))))/nav%fac(1)

    call field_zero_edges(nav%u(it(1)))
    call field_zero_edges(nav%v(it(1)))
    call field_zero_edges(nav%w(it(1)))


  end subroutine navier_projection

  function navier_phi_rhs(mpid,nav,var,fvar)
! -----------------------------------------------------------------------
! navier : compute rhs for phi equation with velocity correction method
! -----------------------------------------------------------------------
! Alexandre Poux
! 01/2013
!
    implicit none
    type(navier3d),intent(in) :: nav
    type(mpi_data),intent(in) :: mpid
    type(field),intent(in)    :: var(nav%nt),fvar(2)
    type(field)    :: navier_phi_rhs
    integer(ik) :: it(nav%nt),nt

    call field_init(navier_phi_rhs,"phi_rhs",nav%nx,nav%ny,nav%nz)


    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    navier_phi_rhs%f=0._rk
    if(nav%pt<=2) then
       navier_phi_rhs=nav%fac(1)*var(it(1))
    elseif(nav%pt==3) then
       navier_phi_rhs=fvar(1) &
            -nav%fac(2)*var(nav%it(nav%nt  )) &
            -nav%fac(3)*var(nav%it(nav%nt-1)) &
            -nav%fac(4)*var(nav%it(nav%nt-2)) 
    elseif(nav%pt==4) then
      if(nav%subite==1) then
       navier_phi_rhs=fvar(1) &
            -nav%fac(2)*var(nav%it(nav%nt  )) &
            -nav%fac(3)*var(nav%it(nav%nt-1)) &
            -nav%fac(4)*var(nav%it(nav%nt-2)) 

       navier_phi_rhs=navier_phi_rhs - fvar(2) &
            +nav%fac(1)*var(nav%it(nav%nt  )) &
            +nav%fac(2)*var(nav%it(nav%nt-1)) &
            +nav%fac(3)*var(nav%it(nav%nt-2)) &
            +nav%fac(4)*var(nav%it(nav%nt-3)) 
      else
        navier_phi_rhs=fvar(1)- fvar(2) + nav%fac(1)*var(nav%it(1))
      endif
    endif
       call field_zero_edges(navier_phi_rhs)

  end function navier_phi_rhs


  function navier_extrapol(nav,var,type,ordre)
    implicit none
    type(navier3d),intent(in) :: nav
    type(field),intent(in)    :: var(nav%nt)
    type(field)               :: navier_extrapol
    character(*),optional     :: type
    integer,optional          :: ordre
    integer                   :: ordre1

    call field_init(navier_extrapol,"extrapol",nav%nx,nav%ny,nav%nz)

    ordre1=0
    if(present(ordre)) then
      ordre1=ordre
    elseif(present(type)) then
      if(type=='v') ordre1=nav%tou
      if(type=='p') ordre1=nav%top-1
    endif
       
    navier_extrapol%f=0._rk
    if(nav%subite==1) then
      if (ordre1==-1) then
         navier_extrapol=1._rk*var(nav%it(1))
			elseif (ordre1==1) then
			   navier_extrapol=1._rk*var(nav%it(nav%nt))
			elseif (ordre1==2) then
			   navier_extrapol=2._rk*var(nav%it(nav%nt))&
                        -1._rk*var(nav%it(nav%nt-1))
      elseif(ordre1==3) then
			   navier_extrapol=3._rk*var(nav%it(nav%nt))&
			                  -3._rk*var(nav%it(nav%nt-1))&
                        +1._rk*var(nav%it(nav%nt-2))
			endif
    else
  	    navier_extrapol%f=var(nav%it(1))%f
    endif

    call field_zero_edges(navier_extrapol)


  end function navier_extrapol


  subroutine navier_time(nav)
! -----------------------------------------------------------------------
! navier : update time variable
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(navier3d) :: nav
    
    nav%time=nav%time+nav%ts

  end subroutine navier_time

  subroutine navier_initialization(cmd,mpid,nav)
! -----------------------------------------------------------------------
! navier : initialize navier type
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    use command_line
    implicit none
    type(cmd_line) :: cmd
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: nx,ny,nz,i,j
    integer(ik) :: bctu(6),bctv(6),bctw(6),bctp(6)

    if(cmd%ndx*cmd%ndy*cmd%ndz.ne.1.or.maxval(cmd%periods).ne.0) then
    !--------------------------------------------------------------------
    !-> initialize mpi
      call md_mpi_init(mpid,cmd)
    !-> initialize petsc
      call md_petsc_initialize()
    else
      mpid%rank=0
      mpid%dims=0
    endif
    if (mpid%rank==0) then
       call color(ired);print'(a)','Precomputation : ';call color(color_off)
    endif

    !-> pressure singular method
    if (cmd%psm==0) nullv=0
    if (cmd%psm==1) nullv=1
    if (cmd%psm==2) nullv=2

    !-> put dimensions in variables for ease of use
    nav%nx=cmd%nx ; nav%ny=cmd%ny ; nav%nz=cmd%nz
    nx=cmd%nx ; ny=cmd%ny ; nz=cmd%nz

    !-> time, time step, nlt
    nav%time=0._rk
    nav%ts=cmd%ts
    nav%ntime=int(cmd%ntime)
    nav%it(:)=(/(i,i=1,nt)/)
    nav%nlt=cmd%nlt
    nav%pt=cmd%pt
    nav%tou=cmd%tou
    nav%top=cmd%top
    nav%nsubite=cmd%nsubite
    nav%les_type=cmd%les_type
    nav%les_c=cmd%les_c**2
    nav%so=cmd%so

    !-> reynolds number 
    nav%rey=cmd%reynolds

    if (nav%tou==1) then
       nav%fac(1)= 1.0_rk/nav%ts
       nav%fac(2)=-1.0_rk/nav%ts
       nav%fac(3)= 0.0_rk
       nav%fac(4)= 0.0_rk
    elseif (nav%tou==2) then
       nav%fac(1)= 1.5_rk/nav%ts
       nav%fac(2)=-2.0_rk/nav%ts
       nav%fac(3)= 0.5_rk/nav%ts
       nav%fac(4)= 0.0_rk
    elseif(nav%tou==3) then
       nav%fac(1)= 11._rk/(6._rk*nav%ts)
       nav%fac(2)=-18._rk/(6._rk*nav%ts)
       nav%fac(3)=  9._rk/(6._rk*nav%ts)
       nav%fac(4)=- 2._rk/(6._rk*nav%ts)
    endif

    !-> compute sigma
    nav%sigmau=-nav%rey*nav%fac(1)
    nav%sigmap=0.0_rk

    bctu=(/1,1,1,1,1,1/)
    bctv=(/1,1,1,1,1,1/)
    bctw=(/1,1,1,1,1,1/)
    bctp=(/2,2,2,2,2,2/)

    !--------------------------------------------------------------------
    !-> initialize mesh
    call mesh_init(nav%gridx,'gridx','x',nx,1,1)
    call mesh_init(nav%gridy,'gridy','y',1,ny,1)
    call mesh_init(nav%gridz,'gridz','z',1,1,nz)

    !-> initialize grid
    call mesh_grid_init(nav%gridx,'x',nx,1,1,mpid)
    call mesh_grid_init(nav%gridy,'y',1,ny,1,mpid)
    call mesh_grid_init(nav%gridz,'z',1,1,nz,mpid)

    call write_mesh('grid_x',nav%gridx,mpid)
    call write_mesh('grid_y',nav%gridy,mpid)
    call write_mesh('grid_z',nav%gridz,mpid)

    if (mpid%dims.ne.0) then
    !    allocate(nav%infu,nav%infv,nav%infw,nav%infp)
        allocate(nav%infu,nav%infp)
        nav%infv=>nav%infu
        nav%infw=>nav%infu

        !--------------------------------------------------------------------
        !-> start initialization of u influence matrix
        call influence_matrix_init_start(mpid,nav%infu,nav%scu,nav%bcu(1),&
             nav%u(1),nav%fu(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'u')

        !-> start initialization of v influence matrix
    !    call influence_matrix_init_start(mpid,nav%infv,nav%scv,nav%bcv(1),&
    !         nav%v(1),nav%fv(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'v')

        !-> start initialization of w influence matrix
    !    call influence_matrix_init_start(mpid,nav%infw,nav%scw,nav%bcw(1),&
    !         nav%w(1),nav%fw(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'w')

        !-> start initialization of pressure influence matrix
        call influence_matrix_init_start(mpid,nav%infp,nav%scp,nav%bcp(1),&
             nav%p(1),nav%fp(1),nav%sigmap,nav%dcx,nav%dcy,nav%dcz,'p')
    endif

    !--------------------------------------------------------------------
    !-> initialize md poisson solvers coefficients
    if (mpid%dims.ne.0) then
      call md_boundary_condition_init(mpid,nav%infu,bctu)
      call md_boundary_condition_init(mpid,nav%infv,bctv)
      call md_boundary_condition_init(mpid,nav%infw,bctw)
      call md_boundary_condition_init(mpid,nav%infp,bctp)
    endif

    !--------------------------------------------------------------------
    !-> initialize poisson solvers coefficients
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scu,bctu,nav%so)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scv,bctv,nav%so)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scw,bctw,nav%so)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scp,bctp,nav%so)

    !--------------------------------------------------------------------
    !-> initialize type field
    do i=1,nav%nt
       call field_init(nav%u(i),"U",nx,ny,nz)
       call field_init(nav%v(i),"V",nx,ny,nz)
       call field_init(nav%w(i),"W",nx,ny,nz)
       call field_init(nav%p(i),"P",nx,ny,nz)
       call field_init(nav%fu(i),"RHS_U",nx,ny,nz)
       call field_init(nav%fv(i),"RHS_V",nx,ny,nz)
       call field_init(nav%fw(i),"RHS_W",nx,ny,nz)
       call field_init(nav%fp(i),"RHS_P",nx,ny,nz)
       call field_init(nav%phi(i),"P",nx,ny,nz)
       call field_init(nav%les_f(i),"F_LES",nav%nx,nav%ny,nav%nz)
    enddo
    do i=1,2
       call field_init(nav%rhs_u(i),"RHS2_U",nx,ny,nz)
       call field_init(nav%rhs_v(i),"RHS2_V",nx,ny,nz)
       call field_init(nav%rhs_w(i),"RHS2_W",nx,ny,nz)
    enddo
    do i=1,3
      do j=1,3
       call field_init(nav%du(i,j),"GRAD_U",nx,ny,nz)
      enddo
    enddo
    call field_init(nav%fphi,"RHS_PHI",nx,ny,nz)
    call field_init(nav%aux,"AUX",nx,ny,nz)
    call field_init(nav%aux1,"AUX",nx,ny,nz)
    call field_init(nav%aux2,"AUX",nx,ny,nz)
    call field_init(nav%rhs_px,"RHS_PX",nx,ny,nz)
    call field_init(nav%rhs_py,"RHS_PY",nx,ny,nz)
    call field_init(nav%rhs_pz,"RHS_PZ",nx,ny,nz)
    call field_init(nav%sub_u,"SUB_U",nx,ny,nz)
    call field_init(nav%sub_v,"SUB_V",nx,ny,nz)
    call field_init(nav%sub_w,"SUB_W",nx,ny,nz)
    call field_init(nav%sub_p,"SUB_P",nx,ny,nz)
    call field_init(nav%les_nu,"les_nu",nx,ny,nz)

    !--------------------------------------------------------------------
    !-> initialize type boundary_conditions
    do i=1,nav%nt
       call boundary_condition_init(nav%bcu(i),nx,ny,nz)
       call boundary_condition_init(nav%bcv(i),nx,ny,nz)
       call boundary_condition_init(nav%bcw(i),nx,ny,nz)
       call boundary_condition_init(nav%bcp(i),nx,ny,nz)
       call boundary_condition_init(nav%bcphi(i),nx,ny,nz)
    enddo

    !--------------------------------------------------------------------
    !-> initialisation of derivatives coefficients
    call derivatives_coefficients_init(nav%gridx,nav%dcx,nx,nav%so,solver='yes')
    call derivatives_coefficients_init(nav%gridy,nav%dcy,ny,nav%so,solver='yes')
    call derivatives_coefficients_init(nav%gridz,nav%dcz,nz,nav%so,solver='yes')

    if (mpid%dims.ne.0) then
        !--------------------------------------------------------------------
        !-> end initialize u influence matrix
        call influence_matrix_init_end(mpid,nav%infu,nav%scu,nav%bcu(1),&
             nav%u(1),nav%fu(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'u')

        !-> end initialize v influence matrix
    !    call influence_matrix_init_end(mpid,nav%infv,nav%scv,nav%bcv(1),&
    !         nav%v(1),nav%fv(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'v')

        !-> end initialize w influence matrix
    !    call influence_matrix_init_end(mpid,nav%infw,nav%scw,nav%bcw(1),&
    !         nav%w(1),nav%fw(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'w')

        !-> end initialize velocity influence matrix
        call influence_matrix_init_end(mpid,nav%infp,nav%scp,nav%bcp(1),&
             nav%p(1),nav%fp(1),nav%sigmap,nav%dcx,nav%dcy,nav%dcz,'p',&
             null=nullv)

        !--------------------------------------------------------------------
        !-> initialize u multidomain solver
    !    call md_solve_init(mpid,nav%infu)
        call md_solve_init(mpid,nav%infu,kspn='u_')

        !-> initialize v multidomain solver
    !    call md_solve_init(mpid,nav%infv)
     
        !-> initialize w multidomain solver
    !    call md_solve_init(mpid,nav%infw)
     
        !-> initialize pressuremultidomain solver
    !    call md_solve_init(mpid,nav%infp,null=nullv)
        call md_solve_init(mpid,nav%infp,null=nullv,kspn='p_')

        !--------------------------------------------------------------------
        !-> initialize guess sol
        call md_guess_init(mpid,nav%infu,nav%infsolu)
        call md_guess_init(mpid,nav%infv,nav%infsolv)
        call md_guess_init(mpid,nav%infw,nav%infsolw)
        call md_guess_init(mpid,nav%infp,nav%infsolphi)
    endif

    !--------------------------------------------------------------------
    !-> initialize fields

    do i=1,nav%nt
       nav%u(i)%f=0._rk ; nav%v(i)%f=0._rk ; nav%w(i)%f=0._rk ; nav%p(i)%f=0._rk
       nav%fu(i)%f=0._rk ; nav%fv(i)%f=0._rk ; nav%fw(i)%f=0._rk ; nav%fp(i)%f=0._rk
       nav%phi(i)%f=0._rk ; nav%les_f(i)%f=0._rk
       nav%bcu(i)%bcx=0._rk ; nav%bcu(i)%bcy=0._rk ; nav%bcu(i)%bcz=0._rk
       nav%bcv(i)%bcx=0._rk ; nav%bcv(i)%bcy=0._rk ; nav%bcv(i)%bcz=0._rk 
       nav%bcw(i)%bcx=0._rk ; nav%bcw(i)%bcy=0._rk ; nav%bcw(i)%bcz=0._rk 
       nav%bcp(i)%bcx=0._rk ; nav%bcp(i)%bcy=0._rk ; nav%bcp(i)%bcz=0._rk  
       nav%bcphi(i)%bcx=0._rk ; nav%bcphi(i)%bcy=0._rk ; nav%bcphi(i)%bcz=0._rk
    enddo
    do i=1,2
       nav%rhs_u(i)%f=0._rk ; nav%rhs_v(i)%f=0._rk ; nav%rhs_w(i)%f=0._rk
    enddo
    nav%fphi%f=0._rk

  end subroutine navier_initialization

  subroutine navier_finalization(cmd,mpid,nav)
! -----------------------------------------------------------------------
! navier : finalize navier type
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    use command_line
    implicit none
    type(cmd_line) :: cmd
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: i,j
    if (mpid%dims.ne.0) then
        !--------------------------------------------------------------------
        !-> deallocate velocity influence matrix
        call md_influence_matrix_destroy(mpid,nav%infu)

        !-> deallocate velocity influence matrix
        call md_influence_matrix_destroy(mpid,nav%infv)

        !-> deallocate velocity influence matrix
        call md_influence_matrix_destroy(mpid,nav%infw)

        !-> deallocate pressure influence matrix
        call md_influence_matrix_destroy(mpid,nav%infp)
    endif
    !--------------------------------------------------------------------
    !-> destroy type field
    do i=1,nav%nt
       call field_destroy(nav%u(i))
       call field_destroy(nav%v(i))
       call field_destroy(nav%w(i))
       call field_destroy(nav%p(i))
       call field_destroy(nav%phi(i))
       call field_destroy(nav%fu(i))
       call field_destroy(nav%fv(i))
       call field_destroy(nav%fw(i))
       call field_destroy(nav%fp(i))
       call field_destroy(nav%les_f(i))
    enddo
    do i=1,2
       call field_destroy(nav%rhs_u(i))
       call field_destroy(nav%rhs_v(i))
       call field_destroy(nav%rhs_w(i))
    enddo
    do i=1,3
      do j=1,3
       call field_destroy(nav%du(i,j))
      enddo
    enddo
    call field_destroy(nav%fphi)
    call field_destroy(nav%aux)
    call field_destroy(nav%aux1)
    call field_destroy(nav%aux2)
    call field_destroy(nav%rhs_px)
    call field_destroy(nav%rhs_py)
    call field_destroy(nav%rhs_pz)
    call field_destroy(nav%sub_u)
    call field_destroy(nav%sub_v)
    call field_destroy(nav%sub_w)
    call field_destroy(nav%sub_p)
    call field_destroy(nav%les_nu)
    if (mpid%dims.ne.0) then
        !--------------------------------------------------------------------
        !-> finalize petsc
        call md_petsc_finalize()
        !-> finalize mpi
        call md_mpi_finalize(mpid)
    endif
  end subroutine navier_finalization

function integrale(mpid,nav,x)
  use class_field
  use class_md
  use precision
!$ use OMP_LIB
  implicit none
  type(field),intent(in) ::x
  type(navier3d),intent(in) :: nav
  type(mpi_data) :: mpid
  real(rk) :: integrale1,integrale,dx,dy,dz
  integer(ik) :: i,j,k

  integrale1=0._rk

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:integrale1)
  do k=2,x%nz-1
     do j=2,x%ny-1
        do i=2,x%nx-1
           dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
           dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
           dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
           integrale1=integrale1+x%f(i,j,k)*dx*dy*dz
        enddo
     enddo
  enddo
!$OMP  END PARALLEL DO

  !-> x boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:integrale1)
   do j=2,x%ny-1
      do k=2,x%nz-1
         dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
         dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
         dx=(nav%gridx%grid1d(2)-nav%gridx%grid1d(1))*0.5_rk
         integrale1=integrale1+x%f(1,j,k)*dx*dy*dz
         dx=(nav%gridx%grid1d(x%nx)-nav%gridx%grid1d(x%nx-1))*0.5_rk
         integrale1=integrale1+x%f(x%nx,j,k)*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  !-> y boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:integrale1)
   do i=2,x%nx-1
      do k=2,x%nz-1
         dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
         dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
         dy=(nav%gridy%grid1d(2)-nav%gridy%grid1d(1))*0.5_rk
         integrale1=integrale1+x%f(i,1,k)*dx*dy*dz
         dy=(nav%gridy%grid1d(x%ny)-nav%gridy%grid1d(x%ny-1))*0.5_rk
         integrale1=integrale1+x%f(i,x%ny,k)*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  !-> z boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:integrale1)
   do i=2,x%nx-1
      do j=2,x%ny-1
         dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
         dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
         dz=(nav%gridz%grid1d(2)-nav%gridz%grid1d(1))*0.5_rk
         integrale1=integrale1+x%f(i,j,1)*dx*dy*dz
         dz=(nav%gridz%grid1d(x%nz)-nav%gridz%grid1d(x%nz-1))*0.5_rk
         integrale1=integrale1+x%f(i,j,x%nz)*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  if (mpid%dims.ne.0) then
    call md_mpi_reduce_double(mpid,integrale1,integrale)
    call md_mpi_bcast_double(mpid,integrale,0)
  else
    integrale=integrale1
  endif

end function integrale


function norme2(mpid,nav,x)
  use class_field
  use class_md
  use precision
!$ use OMP_LIB
  implicit none
  type(field),intent(in) ::x
  type(mpi_data) :: mpid
  type(navier3d),intent(in) :: nav
  real(rk) :: norme2,som1,dx,dy,dz
  integer(ik) :: i,j,k

!dx=maxval(abs(x%f(2:x%nx-1,2:x%ny-1,:)))
!dy=maxval(abs(x%f(2:x%nx-1,:,2:x%nz-1)))
!dz=maxval(abs(x%f(:,2:x%ny-1,2:x%nz-1)))
!norme2=max(dx,dy,dz)
!norme2=maxval(abs(x%f(2:x%nx-1,2:x%ny-1,2:x%nz-1)))
!return

  som1=0._rk

!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:som1)
  do k=2,x%nz-1
     do j=2,x%ny-1
        do i=2,x%nx-1
           dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
           dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
           dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
           som1=som1+x%f(i,j,k)**2*dx*dy*dz
        enddo
     enddo
  enddo
!$OMP  END PARALLEL DO

  !-> x boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:som1)
   do j=2,x%ny-1
      do k=2,x%nz-1
         dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
         dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
         dx=(nav%gridx%grid1d(2)-nav%gridx%grid1d(1))*0.5_rk
         som1=som1+x%f(1,j,k)**2*dx*dy*dz
         dx=(nav%gridx%grid1d(x%nx)-nav%gridx%grid1d(x%nx-1))*0.5_rk
         som1=som1+x%f(x%nx,j,k)**2*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  !-> y boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:som1)
   do i=2,x%nx-1
      do k=2,x%nz-1
         dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
         dz=(nav%gridz%grid1d(k+1)-nav%gridz%grid1d(k-1))*0.5_rk
         dy=(nav%gridy%grid1d(2)-nav%gridy%grid1d(1))*0.5_rk
         som1=som1+x%f(i,1,k)**2*dx*dy*dz
         dy=(nav%gridy%grid1d(x%ny)-nav%gridy%grid1d(x%ny-1))*0.5_rk
         som1=som1+x%f(i,x%ny,k)**2*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  !-> z boundary
!$OMP  PARALLEL DO &
!$OMP  DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP  SCHEDULE(RUNTIME) &
!$OMP  REDUCTION(+:som1)
   do i=2,x%nx-1
      do j=2,x%ny-1
         dx=(nav%gridx%grid1d(i+1)-nav%gridx%grid1d(i-1))*0.5_rk
         dy=(nav%gridy%grid1d(j+1)-nav%gridy%grid1d(j-1))*0.5_rk
         dz=(nav%gridz%grid1d(2)-nav%gridz%grid1d(1))*0.5_rk
         som1=som1+x%f(i,j,1)**2*dx*dy*dz
         dz=(nav%gridz%grid1d(x%nz)-nav%gridz%grid1d(x%nz-1))*0.5_rk
         som1=som1+x%f(i,j,x%nz)**2*dx*dy*dz
      enddo
  enddo
!$OMP  END PARALLEL DO

  if (mpid%dims.ne.0) then
    call md_mpi_reduce_double(mpid,som1,norme2)
    call md_mpi_bcast_double(mpid,norme2,0)
  else
    norme2=som1
  endif

  norme2=sqrt(norme2)

end function norme2

subroutine error_stop2(error_mesg)
!  use mpi_utils, only : code,rang
  implicit none
  character(*) :: error_mesg
!  call MPI_FINALIZE(code)
!  if (rang==0) then
     print'(a)',error_mesg(1:len_trim(error_mesg))
!  endif
!  stop
  
end subroutine error_stop2

end module class_navier_3D
