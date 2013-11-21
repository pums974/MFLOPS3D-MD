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
  use class_mapping
  use class_solver_3d
  use class_md
  use color_print
!  use command_line
  implicit none

  !-> time scheme order
  integer(ik),parameter :: nt=4
!  integer(ik),parameter :: nullv=0
  integer(ik) :: nullv=0

  type navier3d
     !-> number of time steps of the time scheme
     integer(ik) :: nt=nt
     !-> dimensions
     integer(ik) :: nx,ny,nz
     !-> velocity fields
     type(field) :: u(nt),v(nt),w(nt)
     !-> pressure
     type(field) :: p(nt),phi(nt)
     !-> velocity rhs
     type(field) :: fu(nt),fv(nt),fw(nt)
     !-> pressure rhs
     type(field) :: fp(nt),fphi
     !-> auxiliary field
     type(field) :: aux,aux1
     !-> velocity boundary conditions
!     type(boundary_condition) :: bcu(0:nt),bcv(0:nt),bcw(0:nt)
     type(boundary_condition) :: bcu(nt),bcv(nt),bcw(nt)
     !-> pressure boundary conditions
!     type(boundary_condition) :: bcp(0:nt),bcphi
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
     real(rk) :: time,ts
     integer(ik) :: it(nt),ntime
     !-> reynolds number
     real(rk) :: rey
     !-> nonlinear type, projection type, time order
     integer(ik) :: nlt,pt,tou,top
     !-> mapping : 1 -> yes, 0-> no
     type(mapping) :: dcm
     !-> scheme order : 1 : solver, 2 : derivatives
     integer(ik) :: so(2)
     !-> iteration number for mapping
     integer(ik) :: iterm
  end type navier3d

contains

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
    if (mod(ite,inc)==0.and.ite>0) then
       call write_field('vel_u_'//trim(tn(nav)),nav%u(nav%it(nav%nt)),mpid)
       call write_field('vel_v_'//trim(tn(nav)),nav%v(nav%it(nav%nt)),mpid)
       call write_field('vel_w_'//trim(tn(nav)),nav%w(nav%it(nav%nt)),mpid)
       call write_field('vel_p_'//trim(tn(nav)),nav%p(nav%it(nav%nt)),mpid)
!       call write_field('vel_phi_'//trim(tn(nav)),nav%phi(nav%it(nav%nt)),mpid)
    endif

    !-> write mesh 
!    if (mod(ite,inc)==0.and.ite>0) then
!       call write_mesh('grid_x_'//trim(tn(nav)),nav%gridx,mpid)
!       call write_mesh('grid_y_'//trim(tn(nav)),nav%gridy,mpid)
!       call write_mesh('grid_z_'//trim(tn(nav)),nav%gridz,mpid)
!    endif


  end subroutine navier_write_fields

  subroutine navier_residual(mpid,nav,inc,ite)
! -----------------------------------------------------------------------
! navier : residuals computation
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: inc,ite
    real(rk) :: conv(5),convt(5),stopcode

    !-> compute max of U velocity
    conv(1)=maxval(abs(nav%u(nav%it(nav%nt))%f))

    !-> compute residuals
    conv(2)=maxval(abs(nav%u(nav%it(nav%nt))%f)&
         -abs(nav%u(nav%it(nav%nt-1))%f))
    conv(3)=maxval(abs(nav%v(nav%it(nav%nt))%f)&
         -abs(nav%v(nav%it(nav%nt-1))%f))
    conv(4)=maxval(abs(nav%w(nav%it(nav%nt))%f)&
         -abs(nav%w(nav%it(nav%nt-1))%f))
    
    !-> compute divergence
    nav%aux=derxm(nav%dcm,nav%u(nav%it(nav%nt)))+&
         derym(nav%dcm,nav%v(nav%it(nav%nt)))+&
         derzm(nav%dcm,nav%w(nav%it(nav%nt)))
    conv(5)=maxval(abs(nav%aux%f(2:nav%nx-1,2:nav%ny-1,2:nav%nz-1)))

    !-> mpi max
    call md_mpi_reduce_double_max(mpid,conv,convt,5)
    convt(2)=convt(2)/(convt(1)*nav%ts)
    convt(3)=convt(3)/(convt(1)*nav%ts)
    convt(4)=convt(4)/(convt(1)*nav%ts)
    convt(5)=convt(5)/(convt(1))

    !-> print 
    if (mpid%rank==0) then
       print'(i8,5es17.8,5i8)',ite,nav%time,convt(2),convt(3),convt(4),convt(5),&
            nav%infsolu%iter,nav%infsolv%iter,nav%infsolw%iter,nav%infsolphi%iter,&
            nav%iterm
       write(450,'(i8,5es17.8,5i8)')ite,nav%time,convt(2),convt(3),convt(4),convt(5),&
            nav%infsolu%iter,nav%infsolv%iter,nav%infsolw%iter,nav%infsolphi%iter,&
            nav%iterm
    endif

    stopcode=1._rk/nav%ts-1.e-7_rk
    if (ite>10.and.mpid%rank==0) then
       if (convt(2)>stopcode) then
          print*,'Residuals reach limit : stopping code'
          !-> finalize mpi
          call md_mpi_abort(mpid)
          stop
       endif
    endif

  end subroutine navier_residual

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
!    write(timen,'(i8.8,a,i8.8)')floor(nav%time),'.',&
!         int((nav%time-floor(nav%time))*1.d8)
!    dirname='restart_'//trim(timen)
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
!       call write_field(trim(dirname)//'/rhs_u_'//stepn,nav%fu(it(i)),&
!            mpid,dbl='y',inter='y')
!       call write_field(trim(dirname)//'/rhs_v_'//stepn,nav%fv(it(i)),&
!            mpid,dbl='y',inter='y')
!       call write_field(trim(dirname)//'/rhs_w_'//stepn,nav%fw(it(i)),&
!            mpid,dbl='y',inter='y')
!       call write_field(trim(dirname)//'/rhs_p_'//stepn,nav%fp(it(i)),&
!            mpid,dbl='y',inter='y')
    enddo
!    call write_field(trim(dirname)//'/rhs_phi_'//stepn,nav%fphi,&
!         mpid,dbl='y',inter='y')

    call md_influence_guess_write(mpid,nav%infsolu,trim(dirname)//'/guess_u.dat')
    call md_influence_guess_write(mpid,nav%infsolv,trim(dirname)//'/guess_v.dat')
    call md_influence_guess_write(mpid,nav%infsolw,trim(dirname)//'/guess_w.dat')
    call md_influence_guess_write(mpid,nav%infsolphi,trim(dirname)//'/guess_phi.dat')

    write(proc,'(i0)')mpid%rank
!    call write_solver_coefficient(nav%scu,trim(dirname)//'/scu'//proc)
!    call write_solver_coefficient(nav%scv,trim(dirname)//'/scv'//proc)
!    call write_solver_coefficient(nav%scw,trim(dirname)//'/scw'//proc)
!    call write_solver_coefficient(nav%scp,trim(dirname)//'/scp'//proc)
!    call write_derivatives_coefficient(nav%dcx,trim(dirname)//'/dcx'//proc)
!    call write_derivatives_coefficient(nav%dcy,trim(dirname)//'/dcy'//proc)
!    call write_derivatives_coefficient(nav%dcz,trim(dirname)//'/dcz'//proc)

  end subroutine restart_write

  subroutine restart_read(mpid,nav)
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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> create dirname and test existence
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
!       call read_field(trim(dirname)//'/rhs_u_'//stepn,nav%fu(it(i)),&
!            nav%fu(it(i))%name,mpid,inter='y')
!       call read_field(trim(dirname)//'/rhs_v_'//stepn,nav%fv(it(i)),&
!            nav%fv(it(i))%name,mpid,inter='y')
!       call read_field(trim(dirname)//'/rhs_w_'//stepn,nav%fw(it(i)),&
!            nav%fw(it(i))%name,mpid,inter='y')
!       call read_field(trim(dirname)//'/rhs_p_'//stepn,nav%fp(it(i)),&
!            nav%fp(it(i))%name,mpid,inter='y')
    enddo
!    call read_field(trim(dirname)//'/rhs_phi_'//stepn,nav%fphi,&
!         nav%fphi%name,mpid,inter='y')

    do i=1,nt
       call boundary_put_field(nav%u(it(i)),nav%bcu(it(i)),inter)
       call boundary_put_field(nav%v(it(i)),nav%bcv(it(i)),inter)
       call boundary_put_field(nav%w(it(i)),nav%bcw(it(i)),inter)
       call boundary_put_field(nav%p(it(i)),nav%bcp(it(i)),inter)
       call boundary_put_field(nav%phi(it(i)),nav%bcphi(it(i)),inter)
    enddo

    call md_influence_guess_read(mpid,nav%infsolu,trim(dirname)//'/guess_u.dat')
    call md_influence_guess_read(mpid,nav%infsolv,trim(dirname)//'/guess_v.dat')
    call md_influence_guess_read(mpid,nav%infsolw,trim(dirname)//'/guess_w.dat')
    call md_influence_guess_read(mpid,nav%infsolphi,trim(dirname)//'/guess_phi.dat')

!    write(proc,'(i0)')mpid%rank
!    call read_solver_coefficient(nav%scu,trim(dirname)//'/scu'//proc)
!    call read_solver_coefficient(nav%scv,trim(dirname)//'/scv'//proc)
!    call read_solver_coefficient(nav%scw,trim(dirname)//'/scw'//proc)
!    call read_solver_coefficient(nav%scp,trim(dirname)//'/scp'//proc)
!    call read_derivatives_coefficient(nav%dcx,trim(dirname)//'/dcx'//proc)
!    call read_derivatives_coefficient(nav%dcy,trim(dirname)//'/dcy'//proc)
!    call read_derivatives_coefficient(nav%dcz,trim(dirname)//'/dcz'//proc)

  end subroutine restart_read

  subroutine navier_bc_pressure(mpid,nav)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: i,l,m,c(3),inter(3,2)
    integer(ik) :: it(nav%nt),nt
    
    !-> get interface type
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

!    do i=1,nav%nt
!       nav%bcp(i)%bcx=0._rk
!       nav%bcp(i)%bcy=0._rk
!       nav%bcp(i)%bcz=0._rk
!    enddo

    if (nav%dcm%mapt==0) then
       nav%bcphi(it(1))%bcx=0._rk
       nav%bcphi(it(1))%bcy=0._rk
       nav%bcphi(it(1))%bcz=0._rk
    elseif (nav%dcm%mapt==1) then
       call mapping_bcphi(nav%dcm,nav%aux,nav%phi(it(nt)),nav%bcphi(it(1)))
!       nav%bcphi(it(1))%bcx=0._rk
!       nav%bcphi(it(1))%bcy=0._rk
!       nav%bcphi(it(1))%bcz=0._rk
    endif

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    call navier_bc_velocity_utils(inter,nav%bcu(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'u',&
         nav%nx,nav%ny,nav%nz,nav%rey,mpid)

    call navier_bc_velocity_utils(inter,nav%bcv(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'v',&
         nav%nx,nav%ny,nav%nz,nav%rey,mpid)

    call navier_bc_velocity_utils(inter,nav%bcw(nav%it(1)),&
         nav%gridx,nav%gridy,nav%gridz,nav%time,'w',&
         nav%nx,nav%ny,nav%nz,nav%rey,mpid)

    call navier_advection(mpid,nav)

    call erase_boundary_inter(inter,nav%bcu(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcv(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcw(nav%it(1)))

  end subroutine navier_bc_velocity

  subroutine navier_mean(mpid,nav,x,meant)
! -----------------------------------------------------------------------
! navier :  compute mean value in plan xz
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 08/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    type(field) :: x
    real(rk) :: mean,meant,loc
    
    !-> compute mean value
    mean=(sum(x%f))/(x%nx*x%ny*x%nz)

    !-> mpi max
    call md_mpi_reduce_double_sum(mpid,mean,meant)
    meant=meant/(mpid%nd(1)*mpid%nd(2)*mpid%nd(3))
    call md_mpi_bcast_double(mpid,meant,0)

  end subroutine navier_mean

  subroutine navier_testmax(mpid,nav,x1,x2,test)
! -----------------------------------------------------------------------
! navier :  compute mean value in plan xz
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 08/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    type(field) :: x1,x2
    real(rk) :: max(2),maxt(2),test
    
    !-> compute max value
    max(1)=maxval(abs(x1%f(2:nav%nx,2:nav%ny,2:nav%nz)))
    max(2)=maxval(abs(x1%f(2:nav%nx,2:nav%ny,2:nav%nz)&
         -x2%f(2:nav%nx,2:nav%ny,2:nav%nz)))
    call md_mpi_reduce_double_max(mpid,max,maxt,2)
    test=maxt(2)/maxt(1)
    call md_mpi_bcast_double(mpid,test,0)

  end subroutine navier_testmax

  subroutine navier_mean_yz(mpid,nav,x,meant)
! -----------------------------------------------------------------------
! navier :  compute mean value in plan xz
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 08/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    type(field) :: x
    real(rk) :: mean,meant,loc
    
    !-> compute max of U velocity
    if (mpid%coord(1)==mpid%nd(1)-1) then
       mean=abs((sum(x%f(2,1:x%ny,1:x%nz)))/(x%ny*x%nz))
       mean=maxval(x%f(2,1:x%ny,1:x%nz))
    else
       mean=0._rk
    endif

    !-> mpi max
    call md_mpi_reduce_double_sum(mpid,mean,meant)
    meant=meant/(mpid%nd(2)*mpid%nd(3))
    call md_mpi_bcast_double(mpid,meant,0)

  end subroutine navier_mean_yz

  subroutine navier_advection(mpid,nav)
! -----------------------------------------------------------------------
! navier :  advection velocity
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2013
!
    implicit none
    type(mpi_data) :: mpid
    type(navier3d) :: nav
    integer(ik) :: it(nav%nt),nt,ie
    real(rk) :: uc(3),mean(3)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    call navier_mean_yz(mpid,nav,nav%u(it(nt)),mean(1))
    call navier_mean_yz(mpid,nav,nav%v(it(nt)),mean(2))
    call navier_mean_yz(mpid,nav,nav%w(it(nt)),mean(3))
!    if (mpid%rank==0) print*,mpid%coord,'mean :',mean(:)

    uc(:)=0.8_rk
!    uc(:)=mean(1)

    if (mpid%coord(1)==mpid%nd(1)-1) then
       !-> x-direction
       nav%aux=2._rk*derxm(nav%dcm,nav%u(it(nt)))-derxm(nav%dcm,nav%u(it(nt-1)))
       ie=nav%nx
       nav%bcu(it(1))%bcx(:,:,2)=(4._rk*nav%u(it(nt))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -nav%u(it(nt-1))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -2._rk*uc(1)*nav%ts*nav%aux%f(ie,2:nav%ny-1,2:nav%nz-1))/3._rk

       !-> y-direction
       nav%aux=2._rk*derxm(nav%dcm,nav%v(it(nt)))-derxm(nav%dcm,nav%v(it(nt-1)))
       ie=nav%nx
       nav%bcv(it(1))%bcx(:,:,2)=(4._rk*nav%v(it(nt))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -nav%v(it(nt-1))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -2._rk*uc(2)*nav%ts*nav%aux%f(ie,2:nav%ny-1,2:nav%nz-1))/3._rk

       !-> z-direction
       nav%aux=2._rk*derxm(nav%dcm,nav%w(it(nt)))-derxm(nav%dcm,nav%w(it(nt-1)))
       ie=nav%nx
       nav%bcw(it(1))%bcx(:,:,2)=(4._rk*nav%w(it(nt))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -nav%w(it(nt-1))%f(ie,2:nav%ny-1,2:nav%nz-1)&
            -2._rk*uc(3)*nav%ts*nav%aux%f(ie,2:nav%ny-1,2:nav%nz-1))/3._rk
    endif

  end subroutine navier_advection

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
             if (l==1.and.m==1) bc%bcx(:,:,1)=0._rk
             if (l==1.and.m==2) bc%bcx(:,:,2)=0._rk
             if (l==2.and.m==1) bc%bcy(:,:,1)=0._rk
             if (l==2.and.m==2) bc%bcy(:,:,2)=0._rk
             if (l==3.and.m==1) bc%bcz(:,:,1)=0._rk
             if (l==3.and.m==2) bc%bcz(:,:,2)=0._rk
          endif
       enddo
    enddo

  end subroutine erase_boundary_inter

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    
    fac1=1._rk
    fac2=1._rk/sqrt(2._rk)

    if (nav%tou==2) then
       fac3=2._rk*nav%ts/3._rk
    elseif(nav%tou==3) then
       fac3=6._rk*nav%ts/11._rk
    endif
    !fac3=2._rk*nav%ts/3._rk
    !fac3=1._rk

    if (nav%pt==1) then
    !-> bcx
    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derym(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derym(nav%dcm,nav%phi(it(nt)))-derym(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derym(nav%dcm,nav%phi(it(nt)))-3._rk*derym(nav%dcm,nav%phi(it(nt-1)))&
            +derym(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcv(it(1))%bcx(:,:,1)=fac1*(nav%bcv(it(1))%bcx(:,:,1)&
         +fac3*nav%aux%f(1,2:nav%ny-1,2:nav%nz-1))
    nav%bcv(it(1))%bcx(:,:,2)=fac1*(nav%bcv(it(1))%bcx(:,:,2)&
         +fac3*nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1))

    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derzm(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derzm(nav%dcm,nav%phi(it(nt)))-derzm(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derzm(nav%dcm,nav%phi(it(nt)))-3._rk*derzm(nav%dcm,nav%phi(it(nt-1)))&
            +derzm(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcw(it(1))%bcx(:,:,1)=fac1*(nav%bcw(it(1))%bcx(:,:,1)&
         +fac3*nav%aux%f(1,2:nav%ny-1,2:nav%nz-1))
    nav%bcw(it(1))%bcx(:,:,2)=fac1*(nav%bcw(it(1))%bcx(:,:,2)&
         +fac3*nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1))

    !-> bcy
    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derxm(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derxm(nav%dcm,nav%phi(it(nt)))-derxm(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derxm(nav%dcm,nav%phi(it(nt)))-3._rk*derxm(nav%dcm,nav%phi(it(nt-1)))&
            +derxm(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcu(it(1))%bcy(:,:,1)=fac1*(nav%bcu(it(1))%bcy(:,:,1)&
         +fac3*nav%aux%f(2:nav%nx-1,1,2:nav%nz-1))
    nav%bcu(it(1))%bcy(:,:,2)=fac1*(nav%bcu(it(1))%bcy(:,:,2)&
         +fac3*nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1))
    
    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derzm(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derzm(nav%dcm,nav%phi(it(nt)))-derzm(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derzm(nav%dcm,nav%phi(it(nt)))-3._rk*derzm(nav%dcm,nav%phi(it(nt-1)))&
            +derzm(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcw(it(1))%bcy(:,:,1)=fac1*(nav%bcw(it(1))%bcy(:,:,1)&
         +fac3*nav%aux%f(2:nav%nx-1,1,2:nav%nz-1))
    nav%bcw(it(1))%bcy(:,:,2)=fac1*(nav%bcw(it(1))%bcy(:,:,2)&
         +fac3*nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1))

    !-> bcz
    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derxm(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derxm(nav%dcm,nav%phi(it(nt)))-derxm(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derxm(nav%dcm,nav%phi(it(nt)))-3._rk*derxm(nav%dcm,nav%phi(it(nt-1)))&
            +derxm(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcu(it(1))%bcz(:,:,1)=fac1*(nav%bcu(it(1))%bcz(:,:,1)&
         +fac3*nav%aux%f(2:nav%nx-1,2:nav%ny-1,1))
    nav%bcu(it(1))%bcz(:,:,2)=fac1*(nav%bcu(it(1))%bcz(:,:,2)&
         +fac3*nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz))

    nav%aux%f=0._rk
    if (nav%top==2) then
       nav%aux=derym(nav%dcm,nav%phi(it(nt)))
    elseif (nav%top==3) then
       nav%aux=2._rk*derym(nav%dcm,nav%phi(it(nt)))-derym(nav%dcm,nav%phi(it(nt-1)))
    elseif (nav%top==4) then
       nav%aux=3._rk*derym(nav%dcm,nav%phi(it(nt)))-3._rk*derym(nav%dcm,nav%phi(it(nt-1)))&
            +derym(nav%dcm,nav%phi(it(nt-2)))
    endif
    nav%bcv(it(1))%bcz(:,:,1)=fac1*(nav%bcv(it(1))%bcz(:,:,1)&
         +fac3*nav%aux%f(2:nav%nx-1,2:nav%ny-1,1))
    nav%bcv(it(1))%bcz(:,:,2)=fac1*(nav%bcv(it(1))%bcz(:,:,2)&
         +fac3*nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz))

    call erase_boundary_inter(inter,nav%bcu(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcv(nav%it(1)))
    call erase_boundary_inter(inter,nav%bcw(nav%it(1)))

    endif

  end subroutine add_boundary_gradient

  subroutine navier_bc_velocity_utils(inter,bc,gridx,gridy,gridz,t,&
       var,nx,ny,nz,rey,mpid)
! -----------------------------------------------------------------------
! navier : 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(mpi_data) :: mpid
    integer(ik) :: l,m,inter(3,2)
    type(boundary_condition) :: bc
    type(mesh_grid) :: gridx,gridy,gridz
    real(rk) :: x,y,z,t,rey,hy,hz,umax
    integer(ik) :: i,j,k,nx,ny,nz
    character(*) :: var

    hy=3._rk ; hz=2._rk
!    hy=0.2_rk ; hz=0.2_rk
    umax=1._rk
!    umax=0.1_rk

    !-> boundary condition
    !-> x-direction
    do k=2,nz-1
       do j=2,ny-1
          y=gridy%grid1d(j)
          z=gridz%grid1d(k)
          
          x=gridx%grid1d(1)
!          bc%bcx(j-1,k-1,1)=sol(x,y,z,t,var,rey)
          if (var=='u'.and.mpid%coord(1)==0) then
             bc%bcx(j-1,k-1,1)=umax*(1._rk-((y-hy/2._rk)/(hy/2._rk))**2)!&
!                  *(1._rk-((z-hz/2._rk)/(hz/2._rk))**2)
             !bc%bcx(j-1,k-1,1)=1._rk*(1._rk-(y)**2)&
             !     *(1._rk-(z)**2)
          else
             bc%bcx(j-1,k-1,1)=0._rk
          endif

          x=gridx%grid1d(nx)
!          bc%bcx(j-1,k-1,2)=sol(x,y,z,t,var,rey)
!          if (var=='v') &
!          bc%bcx(j-1,k-1,2)=sol(x,y,z,t,'dx'//var,rey)
          bc%bcx(j-1,k-1,2)=0._rk
       enddo
    enddo
    
    !print*,x,y,z,t
    !-> y-direction
    do k=2,nz-1
       do i=2,nx-1
          x=gridx%grid1d(i)
          z=gridz%grid1d(k)

          y=gridy%grid1d(1)
!          bc%bcy(i-1,k-1,1)=sol(x,y,z,t,var,rey)
          bc%bcy(i-1,k-1,1)=0._rk
          y=gridy%grid1d(ny)
!          bc%bcy(i-1,k-1,2)=sol(x,y,z,t,var,rey)
!          if (var=='u') then
!             bc%bcy(i-1,k-1,2)=1._rk
!          else
             bc%bcy(i-1,k-1,2)=0._rk
!          endif
       enddo
    enddo
    !-> z-direction
    do j=2,ny-1
       do i=2,nx-1
          x=gridx%grid1d(i)
          y=gridy%grid1d(j)
          
          z=gridz%grid1d(1)
!          bc%bcz(i-1,j-1,1)=sol(x,y,z,t,var,rey)
          bc%bcz(i-1,j-1,1)=0._rk
          z=gridz%grid1d(nz)
!          bc%bcz(i-1,j-1,2)=sol(x,y,z,t,var,rey)
          bc%bcz(i-1,j-1,2)=0._rk
       enddo
    enddo

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
    real(rk) :: x,y,z,t,a,g,pi
    character(*) :: type
    
    pi=4._rk*atan(1._rk)
    a=1._rk*pi ; g=2._rk*pi

    if (type=="u") then
       sol=0._rk
    endif
    if (type=="v") then
       sol=0._rk
    endif
    if (type=="w") then
       sol=0._rk
    endif
    if (type=="p") then
       sol=0._rk
    endif
    if (type=="dxp") then
       sol=0._rk
    endif
    if (type=="dyp") then
       sol=0._rk
    endif
    if (type=="dzp") then
       sol=0._rk
    endif

    if (type=="rhsu") then
       sol=0._rk
    endif
    if (type=="rhsv") then
       sol=0._rk
    endif
    if (type=="rhsw") then
       sol=0._rk
    endif
    if (type=="rhsp") then
       sol=0._rk
    endif

  end function sol

  function f(nav,var)
! -----------------------------------------------------------------------
! field : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    type(field) :: f
    type(navier3d),intent(in) :: nav
    integer(ik) :: i,j,k
    real(rk) :: x,y,z,t
    character(*) :: var
    call field_init(f,"F",nav%nx,nav%ny,nav%nz)

    t=nav%time
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

  end function f

  subroutine navier_nonlinear(mpid,nav,x,f,choice)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 11/2012
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt
    type(field) :: x(nav%nt),f
    character(*) :: choice
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> nonlinear terms
    if (nav%nlt==1) then
       if (nav%tou==2) then
          f=f+2._rk*(&
               nav%u(it(nt))*derxm(nav%dcm,x(it(nt)))+&
               nav%v(it(nt))*derym(nav%dcm,x(it(nt)))+&
               nav%w(it(nt))*derzm(nav%dcm,x(it(nt))))

          f=f-1._rk*(&
               nav%u(it(nt-1))*derxm(nav%dcm,x(it(nt-1)))+&
               nav%v(it(nt-1))*derym(nav%dcm,x(it(nt-1)))+&
               nav%w(it(nt-1))*derzm(nav%dcm,x(it(nt-1))))
       elseif(nav%tou==3) then
          f=f+3._rk*(&
               nav%u(it(nt))*derxm(nav%dcm,x(it(nt)))+&
               nav%v(it(nt))*derym(nav%dcm,x(it(nt)))+&
               nav%w(it(nt))*derzm(nav%dcm,x(it(nt))))

          f=f-3._rk*(&
               nav%u(it(nt-1))*derxm(nav%dcm,x(it(nt-1)))+&
               nav%v(it(nt-1))*derym(nav%dcm,x(it(nt-1)))+&
               nav%w(it(nt-1))*derzm(nav%dcm,x(it(nt-1))))

          f=f+1._rk*(&
               nav%u(it(nt-2))*derxm(nav%dcm,x(it(nt-2)))+&
               nav%v(it(nt-2))*derym(nav%dcm,x(it(nt-2)))+&
               nav%w(it(nt-2))*derzm(nav%dcm,x(it(nt-2))))
       endif
    elseif (nav%nlt==2) then
       if (nav%tou==2) then
          f=f+1._rk*(&
               nav%u(it(nt))*derxm(nav%dcm,x(it(nt)))+&
               nav%v(it(nt))*derym(nav%dcm,x(it(nt)))+&
               nav%w(it(nt))*derzm(nav%dcm,x(it(nt))))

          nav%aux=1._rk*x(it(nt))*nav%u(it(nt))
          f=f+derxm(nav%dcm,nav%aux)
          nav%aux=1._rk*x(it(nt))*nav%v(it(nt))
          f=f+derym(nav%dcm,nav%aux)
          nav%aux=1._rk*x(it(nt))*nav%w(it(nt))
          f=f+derzm(nav%dcm,nav%aux)
          
          f=f-0.5_rk*(&
               nav%u(it(nt-1))*derxm(nav%dcm,x(it(nt-1)))+&
               nav%v(it(nt-1))*derym(nav%dcm,x(it(nt-1)))+&
               nav%w(it(nt-1))*derzm(nav%dcm,x(it(nt-1))))
          
          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%u(it(nt-1))
          f=f+derxm(nav%dcm,nav%aux)
          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%v(it(nt-1))
          f=f+derym(nav%dcm,nav%aux)
          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%w(it(nt-1))
          f=f+derzm(nav%dcm,nav%aux)
       elseif(nav%tou==3) then
          f=f+1.5_rk*(&
               nav%u(it(nt))*derxm(nav%dcm,x(it(nt)))+&
               nav%v(it(nt))*derym(nav%dcm,x(it(nt)))+&
               nav%w(it(nt))*derzm(nav%dcm,x(it(nt))))

          nav%aux=1.5_rk*x(it(nt))*nav%u(it(nt))
          f=f+derxm(nav%dcm,nav%aux)
          nav%aux=1.5_rk*x(it(nt))*nav%v(it(nt))
          f=f+derym(nav%dcm,nav%aux)
          nav%aux=1.5_rk*x(it(nt))*nav%w(it(nt))
          f=f+derzm(nav%dcm,nav%aux)

          f=f-1.5_rk*(&
               nav%u(it(nt-1))*derxm(nav%dcm,x(it(nt-1)))+&
               nav%v(it(nt-1))*derym(nav%dcm,x(it(nt-1)))+&
               nav%w(it(nt-1))*derzm(nav%dcm,x(it(nt-1))))
          
          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%u(it(nt-1))
          f=f+derxm(nav%dcm,nav%aux)
          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%v(it(nt-1))
          f=f+derym(nav%dcm,nav%aux)
          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%w(it(nt-1))
          f=f+derzm(nav%dcm,nav%aux)

          f=f+0.5_rk*(&
               nav%u(it(nt-2))*derxm(nav%dcm,x(it(nt-2)))+&
               nav%v(it(nt-2))*derym(nav%dcm,x(it(nt-2)))+&
               nav%w(it(nt-2))*derzm(nav%dcm,x(it(nt-2))))

          nav%aux=(0.5_rk)*x(it(nt-2))*nav%u(it(nt-2))
          f=f+derxm(nav%dcm,nav%aux)
          nav%aux=(0.5_rk)*x(it(nt-2))*nav%v(it(nt-2))
          f=f+derym(nav%dcm,nav%aux)
          nav%aux=(0.5_rk)*x(it(nt-2))*nav%w(it(nt-2))
          f=f+derzm(nav%dcm,nav%aux)
       endif
    endif

  end subroutine navier_nonlinear

  subroutine navier_lap_nc(mpid,nav,x,f)
! -----------------------------------------------------------------------
! navier : solve u helmholtz problem
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 08/2013
!
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: it(nav%nt),nt
    type(field) :: x(nav%nt),f
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> non-cartesian part of laplacian
    if (nav%tou==2) then
       f=f+(-2._rk*(dderxm_nc(nav%dcm,x(it(nt)))+dderym_nc(nav%dcm,x(it(nt))) &
            +dderzm_nc(nav%dcm,x(it(nt))))+ &
            (dderxm_nc(nav%dcm,x(it(nt-1)))+dderym_nc(nav%dcm,x(it(nt-1))) &
            +dderzm_nc(nav%dcm,x(it(nt-1)))))/nav%rey
    elseif(nav%tou==3) then
       f=f+(-3._rk*(dderxm_nc(nav%dcm,x(it(nt)))+dderym_nc(nav%dcm,x(it(nt))) &
            +dderzm_nc(nav%dcm,x(it(nt))))+ &
            3._rk*(dderxm_nc(nav%dcm,x(it(nt-1)))+dderym_nc(nav%dcm,x(it(nt-1))) &
            +dderzm_nc(nav%dcm,x(it(nt-1))))- &
            (dderxm_nc(nav%dcm,x(it(nt-2)))+dderym_nc(nav%dcm,x(it(nt-2))) &
            +dderzm_nc(nav%dcm,x(it(nt-2)))))/nav%rey
    endif

  end subroutine navier_lap_nc

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
!    nav%fu(it(nt))=-nav%u(it(nt))/nav%ts
    if (nav%tou==2) then
       nav%fu(it(nt))=0.5_rk*(-4._rk*nav%u(it(nt))+nav%u(it(nt-1)))/nav%ts
    elseif(nav%tou==3) then
       nav%fu(it(nt))=(-18._rk*nav%u(it(nt))+9._rk*nav%u(it(nt-1))&
            -2._rk*nav%u(it(nt-2)))/(6._rk*nav%ts)
    endif
    nav%fu(it(nt))=nav%fu(it(nt))+dertm_nc(nav%dcm,nav%u(it(nt)))

    !-> pressure 
    if (nav%pt==2) then
       if (nav%top==2) then
          nav%fu(it(nt))=nav%fu(it(nt))+derxm(nav%dcm,nav%p(it(nt)))
       elseif(nav%top==3) then
          nav%fu(it(nt))=nav%fu(it(nt))+2._rk*derxm(nav%dcm,nav%p(it(nt)))&
               -derxm(nav%dcm,nav%p(it(nt-1)))
       endif
    endif
    
    !-> nonlinear terms
    call navier_nonlinear(mpid,nav,nav%u,nav%fu(it(nt)),'u')

    !-> non-cartesian laplacian
    call navier_lap_nc(mpid,nav,nav%u,nav%fu(it(nt)))

    !-> function
    nav%fu(it(nt))=nav%fu(it(nt))-f(nav,'rhsu')

    !-> reynolds number multiplication
    nav%fu(it(nt))=nav%rey*nav%fu(it(nt))
    
    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infu,nt,it,nav%bcu,nav%u)
    call multidomain_solve(mpid,nav%infu,nav%scu,nav%bcu(it(1)),nav%u(it(1)),&
          nav%fu(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolu)

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
!    nav%fv(it(nt))=-nav%v(it(nt))/nav%ts
    if (nav%tou==2) then
       nav%fv(it(nt))=0.5_rk*(-4._rk*nav%v(it(nt))+nav%v(it(nt-1)))/nav%ts
    elseif(nav%tou==3) then
       nav%fv(it(nt))=(-18._rk*nav%v(it(nt))+9._rk*nav%v(it(nt-1))&
            -2._rk*nav%v(it(nt-2)))/(6._rk*nav%ts)
    endif
    nav%fv(it(nt))=nav%fv(it(nt))+dertm_nc(nav%dcm,nav%v(it(nt)))
    
    !-> pressure 
    if (nav%pt==2) then
       if (nav%top==2) then
          nav%fv(it(nt))=nav%fv(it(nt))+derym(nav%dcm,nav%p(it(nt)))
       elseif(nav%top==3) then
          nav%fv(it(nt))=nav%fv(it(nt))+2._rk*derym(nav%dcm,nav%p(it(nt)))&
               -derym(nav%dcm,nav%p(it(nt-1)))
       endif
    endif
    
    !-> nonlinear terms
    call navier_nonlinear(mpid,nav,nav%v,nav%fv(it(nt)),'v')

    !-> non-cartesian laplacian
    call navier_lap_nc(mpid,nav,nav%v,nav%fv(it(nt)))

    !-> function
    nav%fv(it(nt))=nav%fv(it(nt))-f(nav,'rhsv')

    !-> reynolds number multiplication
    nav%fv(it(nt))=nav%rey*nav%fv(it(nt))

    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infv,nt,it,nav%bcv,nav%v)
    call multidomain_solve(mpid,nav%infv,nav%scv,nav%bcv(it(1)),nav%v(it(1)),&
          nav%fv(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolv)

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
!    nav%fw(it(nt))=-nav%w(it(nt))/nav%ts
    if (nav%tou==2) then
       nav%fw(it(nt))=0.5_rk*(-4._rk*nav%w(it(nt))+nav%w(it(nt-1)))/nav%ts
    elseif(nav%tou==3) then
       nav%fw(it(nt))=(-18._rk*nav%w(it(nt))+9._rk*nav%w(it(nt-1))&
            -2._rk*nav%w(it(nt-2)))/(6._rk*nav%ts)
    endif
    nav%fw(it(nt))=nav%fw(it(nt))+dertm_nc(nav%dcm,nav%w(it(nt)))
    
    !-> pressure 
    if (nav%pt==2) then
       if (nav%top==2) then
          nav%fw(it(nt))=nav%fw(it(nt))+derzm(nav%dcm,nav%p(it(nt))) 
       elseif(nav%top==3) then
          nav%fw(it(nt))=nav%fw(it(nt))+2._rk*derzm(nav%dcm,nav%p(it(nt)))&
               -derzm(nav%dcm,nav%p(it(nt-1)))
       endif
    endif
    
    !-> nonlinear terms
    call navier_nonlinear(mpid,nav,nav%w,nav%fw(it(nt)),'w')

    !-> non-cartesian laplacian
    call navier_lap_nc(mpid,nav,nav%w,nav%fw(it(nt)))

    !-> function
    nav%fw(it(nt))=nav%fw(it(nt))-f(nav,'rhsw')

    !-> reynolds number multiplication
    nav%fw(it(nt))=nav%rey*nav%fw(it(nt))

    !--------------------------------------------------------------------
    !-> solve
!    call md_set_guess(mpid,nav%infw,nt,it,nav%bcw,nav%w)
    call multidomain_solve(mpid,nav%infw,nav%scw,nav%bcw(it(1)),nav%w(it(1)),&
          nav%fw(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz,&
          inf_sol=nav%infsolw)

  end subroutine navier_solve_w

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
    integer(ik) :: it(nav%nt),nt,k,kmax
    real(rk) :: fac,max(2),maxt(2),meant,test

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    nav%aux1=nav%phi(it(nt))
    kmax=1
    if (nav%dcm%mapt==1) kmax=50
    do k=1,kmax
    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fphi%f=0._rk

    if (nav%tou==2) then
       fac=1.5_rk/nav%ts
    elseif(nav%tou==3) then
       fac=11._rk/(6._rk*nav%ts)
    endif

    !-> 
    nav%fphi=fac*(&
         derxm(nav%dcm,nav%u(it(1)))+&
         derym(nav%dcm,nav%v(it(1)))+&
         derzm(nav%dcm,nav%w(it(1))))
    
    !-> non-cartesian laplacian
       nav%fphi=nav%fphi-(&
            dderxm_nc(nav%dcm,nav%aux1)&
            +dderym_nc(nav%dcm,nav%aux1) &
            +dderzm_nc(nav%dcm,nav%aux1))

    !--------------------------------------------------------------------
    !-> solve
    call multidomain_solve(mpid,nav%infp,nav%scp,nav%bcphi(it(1)),nav%phi(it(1)),&
          nav%fphi,nav%aux,nav%sigmap,nav%dcx,nav%dcy,nav%dcz,null=nullv,&
          inf_sol=nav%infsolphi,var='p')

    call navier_bc_pressure(mpid,nav)

!    call navier_mean(mpid,nav,nav%phi(it(1)),meant)
!    nav%phi(it(1))=nav%phi(it(1))-meant

    call navier_testmax(mpid,nav,nav%phi(it(1)),nav%aux1,test)
!    if (mpid%rank==1) print*,test
    nav%iterm=k
    if (test<nav%ts) exit
    
    nav%aux1=nav%phi(it(1))
 enddo

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> compute extrema
    ex(1,1)=2 ; ex(2,1)=2 ; ex(3,1)=2 
    ex(1,2)=nav%nx-1 ; ex(2,2)=nav%ny-1 ; ex(3,2)=nav%nz-1

    !-> coefficient
!    fac=nav%ts
    if (nav%tou==2) then
       fac=2._rk*nav%ts/3._rk
    elseif(nav%tou==3) then
       fac=6._rk*nav%ts/11._rk
    endif

    call field_zero_edges(nav%phi(it(1)))
    !-> pressure
    if (nav%pt==1) then
       nav%p(it(1))=nav%phi(it(1))
    elseif(nav%pt==2) then
       if (nav%top==2) then
          nav%p(it(1))=nav%phi(it(1))+nav%p(nav%it(nav%nt))
       elseif(nav%top==3) then
          nav%p(it(1))=nav%phi(it(1))+2._rk*nav%p(nav%it(nav%nt))-nav%p(nav%it(nav%nt-1))
       endif
    endif

    !-> rotationnal
    nav%aux=(derxm(nav%dcm,nav%u(it(1)))+&
         derym(nav%dcm,nav%v(it(1)))+&
         derzm(nav%dcm,nav%w(it(1))))/nav%rey
    nav%p(it(1))=nav%p(it(1))-nav%aux

    !-> Brown
!    nav%aux=fac*(dderxm(nav%dcm,nav%phi(it(1)))+&
!         dderym(nav%dcm,nav%phi(it(1)))+&
!         dderzm(nav%dcm,nav%phi(it(1))))/nav%rey
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

!    goto 101
!    call navier_bc_velocity(mpid,nav)
!    call field_put_boundary(nav%u(it(1)),nav%bcu(it(1)),inter)
!    call field_put_boundary(nav%v(it(1)),nav%bcv(it(1)),inter)
!    call field_put_boundary(nav%w(it(1)),nav%bcw(it(1)),inter)

    !-> velocity
    goto 101
    nav%aux=derxm(nav%dcm,nav%phi(it(1)))
    nav%u(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%u(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))

    nav%aux=derym(nav%dcm,nav%phi(it(1)))
    nav%v(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%v(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))

    nav%aux=derzm(nav%dcm,nav%phi(it(1)))
    nav%w(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))=&
         nav%w(it(1))%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))&
         -fac*nav%aux%f(ex(1,1):ex(1,2),ex(2,1):ex(2,2),ex(3,1):ex(3,2))
101 continue

    nav%u(it(1))=nav%u(it(1))-fac*derxm(nav%dcm,nav%phi(it(1)))
    nav%v(it(1))=nav%v(it(1))-fac*derym(nav%dcm,nav%phi(it(1)))
    nav%w(it(1))=nav%w(it(1))-fac*derzm(nav%dcm,nav%phi(it(1)))

    
    call field_zero_edges(nav%u(it(1)))
    call field_zero_edges(nav%v(it(1)))
    call field_zero_edges(nav%w(it(1)))

    !-> switch it    
    iaux=nav%it(1)
    do i=1,nav%nt-1
       nav%it(i)=nav%it(i+1)
    enddo
    nav%it(nt)=iaux

  end subroutine navier_projection


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

  subroutine navier_mapping(mpid,nav)
! -----------------------------------------------------------------------
! navier : compute mapping
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 08/2013
!
    use command_line
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid

    !-> initialize mapping
    call mapping_init(mpid,nav%gridx,nav%gridy,nav%dcx,nav%dcy,nav%dcz,&
         nav%aux,nav%dcm,nav%time)

  end subroutine navier_mapping

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
    integer(ik) :: nx,ny,nz,i
    integer(ik) :: bctu(6),bctv(6),bctw(6),bctp(6)

    !--------------------------------------------------------------------
    !-> initialize mpi
    call md_mpi_init(mpid,cmd)
    !-> initialize petsc
    call md_petsc_initialize()

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
    nav%ntime=cmd%ntime
    nav%it(:)=(/(i,i=1,nt)/)
    nav%nlt=cmd%nlt
    nav%pt=cmd%pt
    nav%tou=cmd%tou
    nav%top=cmd%top

    !-> mapping
    nav%dcm%mapt=cmd%mapt

    !-> scheme order : default 6
    nav%so=cmd%so
    if (nav%so(1)==0) nav%so(1)=6
    if (nav%so(2)==0) nav%so(2)=6

    !-> reynolds number 
    nav%rey=cmd%reynolds

    !--------------------------------------------------------------------
    !-> initialize mesh
!    call mesh_init(nav%gridx,'gridx','x',nx,1,1)
!    call mesh_init(nav%gridy,'gridy','y',1,ny,1)
!    call mesh_init(nav%gridz,'gridz','z',1,1,nz)

    call mesh_init(nav%gridx,'gridx','x',nx,ny,nz)
    call mesh_init(nav%gridy,'gridy','y',nx,ny,nz)
    call mesh_init(nav%gridz,'gridz','z',nx,ny,nz)

    !-> initialize grid
    call mesh_grid_init(nav%gridx,'x',nx,1,1,mpid)
    call mesh_grid_init(nav%gridy,'y',nx,ny,1,mpid)
    call mesh_grid_init(nav%gridz,'z',1,1,nz,mpid)

    !-> compute sigma
    if (nav%tou==2) then
!       nav%sigmau=-nav%rey/nav%ts
       nav%sigmau=-nav%rey*1.5_rk/nav%ts
       nav%sigmap=0._rk
    endif
    if (nav%tou==3) then
       nav%sigmau=-nav%rey*11._rk/(6._rk*nav%ts)
       nav%sigmap=0._rk
    endif

!    allocate(nav%infu,nav%infv,nav%infw,nav%infp)
    allocate(nav%infu,nav%infp)
!    allocate(nav%infu,nav%infv,nav%infp)
    nav%infv=>nav%infu
    nav%infw=>nav%infu

    !--------------------------------------------------------------------
    !-> start initialization of u influence matrix
    call influence_matrix_init_start(mpid,nav%infu,nav%scu,nav%bcu(1),&
         nav%u(1),nav%fu(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'u')

    !-> start initialization of u influence matrix
!    call influence_matrix_init_start(mpid,nav%infv,nav%scv,nav%bcv(1),&
!         nav%v(1),nav%fv(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'v')

    !-> start initialization of u influence matrix
!    call influence_matrix_init_start(mpid,nav%infw,nav%scw,nav%bcw(1),&
!         nav%w(1),nav%fw(1),nav%sigmau,nav%dcx,nav%dcy,nav%dcz,'w')

    !-> start initialization of pressure influence matrix
    call influence_matrix_init_start(mpid,nav%infp,nav%scp,nav%bcp(1),&
         nav%p(1),nav%fp(1),nav%sigmap,nav%dcx,nav%dcy,nav%dcz,'p')

    !--------------------------------------------------------------------
    !-> initialize poisson solver coefficient for u
    bctu=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infu,bctu)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scu,bctu,nav%so(1))

    !-> initialize poisson solver coefficient for v
    bctv=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infv,bctv)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scv,bctv,nav%so(1))

    !-> initialize poisson solver coefficient for w
    bctw=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infw,bctw)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scw,bctw,nav%so(1))

    !-> initialize poisson solver coefficient for pressure
    bctp=(/2,2,2,2,2,2/)
    call md_boundary_condition_init(mpid,nav%infp,bctp)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scp,bctp,nav%so(1))

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
    enddo
    do i=1,nav%nt
       call field_init(nav%phi(i),"P",nx,ny,nz)
    enddo
    call field_init(nav%fphi,"RHS_P",nx,ny,nz)
    call field_init(nav%aux,"AUX",nx,ny,nz)
    call field_init(nav%aux1,"AUX",nx,ny,nz)

    !--------------------------------------------------------------------
    !-> initialize type boundary_condition for velocity
    do i=1,nav%nt
       call boundary_condition_init(nav%bcu(i),nx,ny,nz)
       call boundary_condition_init(nav%bcv(i),nx,ny,nz)
       call boundary_condition_init(nav%bcw(i),nx,ny,nz)
    enddo

    !-> initialize type boundary_condition for pressure
    do i=1,nav%nt
       call boundary_condition_init(nav%bcp(i),nx,ny,nz)
    call boundary_condition_init(nav%bcphi(i),nx,ny,nz)
    enddo

    !--------------------------------------------------------------------
    !-> initialisation of derivatives coefficients
    call derivatives_coefficients_init(nav%gridx,nav%dcx,nx,nav%so(2),&
         solver='yes')
    call derivatives_coefficients_init(nav%gridy,nav%dcy,ny,nav%so(2),&
         solver='yes')
    call derivatives_coefficients_init(nav%gridz,nav%dcz,nz,nav%so(2),&
         solver='yes')

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

!    call md_influence_matrix_view(mpid,nav%infp)

    !--------------------------------------------------------------------
    !-> initialize u multidomain solver
!    call md_solve_init(mpid,nav%infu)
    call md_solve_init(mpid,nav%infu,kspn='u_')
 
    !-> initialize v multidomain solver
!    call md_solve_init(mpid,nav%infv,kspn='u_')
 
    !-> initialize w multidomain solver
!    call md_solve_init(mpid,nav%infw,kspn='u_')
 
    !-> initialize pressuremultidomain solver
!    call md_solve_init(mpid,nav%infp,null=nullv)
    call md_solve_init(mpid,nav%infp,null=nullv,kspn='p_')

    !-> initialize mapping
    call navier_mapping(mpid,nav)

    !--------------------------------------------------------------------
    !-> initialize fields

    do i=1,nav%nt
       nav%u(i)%f=0._rk ; nav%v(i)%f=0._rk ; nav%w(i)%f=0._rk ; nav%p(i)%f=0._rk
       nav%fu(i)%f=0._rk ; nav%fv(i)%f=0._rk ; nav%fw(i)%f=0._rk ; nav%fp(i)%f=0._rk
    enddo
    do i=1,nav%nt
       nav%phi(i)%f=0._rk 
    enddo
    nav%fphi%f=0._rk
    do i=1,nav%nt
!    do i=0,nav%nt
       nav%bcu(i)%bcx=0._rk ; nav%bcu(i)%bcy=0._rk ; nav%bcu(i)%bcz=0._rk
       nav%bcv(i)%bcx=0._rk ; nav%bcv(i)%bcy=0._rk ; nav%bcv(i)%bcz=0._rk 
       nav%bcw(i)%bcx=0._rk ; nav%bcw(i)%bcy=0._rk ; nav%bcw(i)%bcz=0._rk 
       nav%bcp(i)%bcx=0._rk ; nav%bcp(i)%bcy=0._rk ; nav%bcp(i)%bcz=0._rk  
       nav%bcphi(i)%bcx=0._rk ; nav%bcphi(i)%bcy=0._rk ; nav%bcphi(i)%bcz=0._rk
    enddo

    !--------------------------------------------------------------------
    !-> initialize guess sol
    call md_guess_init(mpid,nav%infu,nav%infsolu)
    call md_guess_init(mpid,nav%infv,nav%infsolv)
    call md_guess_init(mpid,nav%infw,nav%infsolw)
    call md_guess_init(mpid,nav%infp,nav%infsolphi)

  end subroutine navier_initialization

  subroutine navier_finalization(mpid,nav)
! -----------------------------------------------------------------------
! navier : finalize navier type
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    use command_line
    implicit none
    type(navier3d) :: nav
    type(mpi_data) :: mpid
    integer(ik) :: i

    !--------------------------------------------------------------------
    !-> deallocate velocity influence matrix
    call md_influence_matrix_destroy(mpid,nav%infu)

    !-> deallocate velocity influence matrix
    call md_influence_matrix_destroy(mpid,nav%infv)

    !-> deallocate velocity influence matrix
    call md_influence_matrix_destroy(mpid,nav%infw)

    !-> deallocate pressure influence matrix
    call md_influence_matrix_destroy(mpid,nav%infp)

    !--------------------------------------------------------------------
    !-> destroy type field
    do i=1,nav%nt
       call field_destroy(nav%u(i))
       call field_destroy(nav%v(i))
       call field_destroy(nav%w(i))
       call field_destroy(nav%p(i))
       call field_destroy(nav%fu(i))
       call field_destroy(nav%fv(i))
       call field_destroy(nav%fw(i))
       call field_destroy(nav%fp(i))
    enddo

    !--------------------------------------------------------------------
    !-> finalize petsc
    call md_petsc_finalize()
    !-> finalize mpi
    call md_mpi_finalize(mpid)

  end subroutine navier_finalization

end module class_navier_3D

