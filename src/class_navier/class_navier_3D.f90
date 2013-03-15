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
     !-> auxiliary field
     type(field) :: aux
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
     !-> nonlinear type, projection type, time order
     integer(ik) :: nlt,pt,tou,top
     integer(ik) :: subite,nsubite
  end type navier3d

contains

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infp,c,inter)

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infp,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> m : left-right ; l : directions (x,y,z)
    do m=1,2
      do l=1,3
        if ((i.ge.0.and.inter(l,m)>0).or.(i.le.0.and.inter(l,m)<=0)) then
          select case (l)
            case (1)
              x%f((m-1)*(nav%nx-1)+1,:,:)=0._rk
            case (2)
              x%f(:,(m-1)*(nav%ny-1)+1,:)=0._rk
            case (3)
              x%f(:,:,(m-1)*(nav%nz-1)+1)=0._rk
          end select
        endif
      enddo
    enddo

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
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)
    

    !-> bcx
    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%phi,'p'))
    nav%bcv(it(1))%bcx(:,:,1)=nav%bcv(it(1))%bcx(:,:,1)&
                             +nav%aux%f(1,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)
    nav%bcv(it(1))%bcx(:,:,2)=nav%bcv(it(1))%bcx(:,:,2)&
                             +nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derz(nav%dcz,navier_extrapol(nav,nav%phi,'p'))
    nav%bcw(it(1))%bcx(:,:,1)=nav%bcw(it(1))%bcx(:,:,1)&
                             +nav%aux%f(1,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)
    nav%bcw(it(1))%bcx(:,:,2)=nav%bcw(it(1))%bcx(:,:,2)&
                             +nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%fac(1)

    !-> bcy
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%phi,'p'))
    nav%bcu(it(1))%bcy(:,:,1)=nav%bcu(it(1))%bcy(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,1,2:nav%nz-1)/nav%fac(1)
    nav%bcu(it(1))%bcy(:,:,2)=nav%bcu(it(1))%bcy(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=derz(nav%dcz,navier_extrapol(nav,nav%phi,'p'))
    nav%bcw(it(1))%bcy(:,:,1)=nav%bcw(it(1))%bcy(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,1,2:nav%nz-1)/nav%fac(1)
    nav%bcw(it(1))%bcy(:,:,2)=nav%bcw(it(1))%bcy(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%fac(1)

    !-> bcz
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%phi,'p'))
    nav%bcu(it(1))%bcz(:,:,1)=nav%bcu(it(1))%bcz(:,:,1)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,1)/nav%fac(1)
    nav%bcu(it(1))%bcz(:,:,2)=nav%bcu(it(1))%bcz(:,:,2)&
                             +nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%fac(1)

    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%phi,'p'))
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
    integer(ik) :: it(nav%nt),nt
    real(rk) :: fac1,fac2,fac3
    integer(ik) :: l,m,c(3),inter(3,2)

    !-> get interface type
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(nav%infu,c,inter)

    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !-> bcx
    nav%aux%f=0._rk
    nav%aux=dery(nav%dcy,navier_extrapol(nav,nav%v,'p'))&
           +derz(nav%dcz,navier_extrapol(nav,nav%w,'p'))
    nav%aux=derx(nav%dcx,nav%aux)-ddery(nav%dcy,navier_extrapol(nav,nav%u,'p'))&
                                 -dderz(nav%dcz,navier_extrapol(nav,nav%u,'p'))
    nav%bcphi(it(1))%bcx(:,:,1)=nav%bcphi(it(1))%bcx(:,:,1)-nav%aux%f(1     ,2:nav%ny-1,2:nav%nz-1)/nav%rey
    nav%bcphi(it(1))%bcx(:,:,2)=nav%bcphi(it(1))%bcx(:,:,2)-nav%aux%f(nav%nx,2:nav%ny-1,2:nav%nz-1)/nav%rey

    !-> bcy
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%u,'p'))&
           +derz(nav%dcz,navier_extrapol(nav,nav%w,'p'))
    nav%aux=dery(nav%dcy,nav%aux)-dderx(nav%dcx,navier_extrapol(nav,nav%v,'p'))&
                                 -dderz(nav%dcz,navier_extrapol(nav,nav%v,'p'))
    nav%bcphi(it(1))%bcy(:,:,1)=nav%bcphi(it(1))%bcy(:,:,1)-nav%aux%f(2:nav%nx-1,1     ,2:nav%nz-1)/nav%rey
    nav%bcphi(it(1))%bcy(:,:,2)=nav%bcphi(it(1))%bcy(:,:,2)-nav%aux%f(2:nav%nx-1,nav%ny,2:nav%nz-1)/nav%rey

    !-> bcz
    nav%aux%f=0._rk
    nav%aux=derx(nav%dcx,navier_extrapol(nav,nav%u,'p'))&
           +dery(nav%dcy,navier_extrapol(nav,nav%v,'p'))
    nav%aux=derz(nav%dcz,nav%aux)-dderx(nav%dcx,navier_extrapol(nav,nav%w,'p'))&
                                 -ddery(nav%dcy,navier_extrapol(nav,nav%w,'p'))
    nav%bcphi(it(1))%bcz(:,:,1)=nav%bcphi(it(1))%bcz(:,:,1)-nav%aux%f(2:nav%nx-1,2:nav%ny-1,1     )/nav%rey
    nav%bcphi(it(1))%bcz(:,:,2)=nav%bcphi(it(1))%bcz(:,:,2)-nav%aux%f(2:nav%nx-1,2:nav%ny-1,nav%nz)/nav%rey

    call erase_boundary_inter(inter,nav%bcphi(nav%it(1)))


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
    !print*,x,y,z,t
    !-> y-direction
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
!             bc%bcy(i-1,k-1,2)=1._rk
!          else
!             bc%bcy(i-1,k-1,2)=0._rk
!          endif
       enddo
    enddo
    !-> z-direction
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
    integer(ik) ::i,n
    PARAMETER(n=1)
    real(rk) :: x,y,z,t,pi,wx(n),wt(n),wy(n),wp(n),a(n),b(n)
    character(*) :: type

    pi=4._rk*atan(1._rk)
    sol=0._rk
    do i=1,n
      a(i)=1._rk*i ; b(i)=2._rk*i
      wx(i)=1._rk*pi*i ; wy(i)=2._rk*pi*i ; wt(i)=1._rk*i*pi ; wp(i)=1.5_rk*pi*i

    if (type=="u") then
       sol=sol+a(i)*sin(wx(i)*x)*cos(wy(i)*y)*cos(wt(i)*t)/wx(i)
    endif
    if (type=="v") then
       sol=sol-a(i)*cos(wx(i)*x)*sin(wy(i)*y)*cos(wt(i)*t)/wy(i)
    endif
    if (type=="w") then
       sol=sol+0._rk
    endif
    if (type=="p") then
       sol=sol+b(i)*sin(wp(i)*(x-y))*cos(wt(i)*t)
    endif
    if (type=="dxp") then
       sol=sol+b(i)*wp(i)*cos(t*wt(i))*cos(wp(i)*(x-y))
    endif
    if (type=="dyp") then
       sol=sol-b(i)*wp(i)*cos(t*wt(i))*cos(wp(i)*(x-y))
    endif
    if (type=="dzp") then
       sol=sol+0._rk
    endif

    if (type=="rhsu") then
       sol=sol+(a(i)*cos(t*wt(i))*wy(i)**2*sin(wx(i)*x)*cos(wy(i)*y)/wx(i)+a(i)*cos(t*wt(i))*wx(i)*sin(wx(i)*x&
            )*cos(wy(i)*y))/rey-a(i)*wt(i)*sin(t*wt(i))*sin(wx(i)*x)*cos(wy(i)*y)/wx(i)+b(i)*wp(i)*cos(t&
            *wt(i))*cos(wp(i)*(x-y))
    endif
    if (type=="rhsv") then
	sol=sol+(-a(i)*cos(t*wt(i))*wy(i)*cos(wx(i)*x)*sin(wy(i)*y)-a(i)*cos(t*wt(i))*wx(i)**2*cos(wx(i)*x)*&
            sin(wy(i)*y)/wy(i))/rey+a(i)*wt(i)*sin(t*wt(i))*cos(wx(i)*x)*sin(wy(i)*y)/wy(i)-b(i)*wp(i)*cos(&
            t*wt(i))*cos(wp(i)*(x-y))
    endif
    if (type=="rhsw") then
       sol=sol+0._rk
    endif
    if (type=="rhsp") then
       sol=sol-2*(a(i)**2*cos(t*wt(i))**2*sin(wx(i)*x)**2*sin(wy(i)*y)**2-a(i)**2*cos(t*wt(i))**2&
           *cos(wx(i)*x)**2*cos(wy(i)*y)**2)-2*b(i)*wp(i)**2*cos(t*wt(i))*sin(wp(i)*(x-y))
    endif
    end do
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
    type(field) :: x(nav%nt),f,tmp(nav%nt)
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    do i=1,nav%nt
      call field_init(tmp(i),"NLT",nav%nx,nav%ny,nav%nz)
    enddo

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

          nav%aux=x(i)*nav%u(i)*0.5_rk
          tmp(i)=tmp(i)+derx(nav%dcx,nav%aux)
          nav%aux=x(i)*nav%v(i)*0.5_rk
          tmp(i)=tmp(i)+dery(nav%dcy,nav%aux)
          nav%aux=x(i)*nav%w(i)*0.5_rk
          tmp(i)=tmp(i)+derz(nav%dcz,nav%aux)

      enddo
    endif
    f=f+navier_extrapol(nav,tmp,'v')

    do i=1,nav%nt
       call field_destroy(tmp(i))
    enddo


!    if (nav%nlt==1) then
!       if (nav%tou==2) then
!          f=f+2._rk*(&
!               nav%u(it(nt))*derx(nav%dcx,x(it(nt)))+&
!               nav%v(it(nt))*dery(nav%dcy,x(it(nt)))+&
!               nav%w(it(nt))*derz(nav%dcz,x(it(nt))))

!          f=f-1._rk*(&
!               nav%u(it(nt-1))*derx(nav%dcx,x(it(nt-1)))+&
!               nav%v(it(nt-1))*dery(nav%dcy,x(it(nt-1)))+&
!               nav%w(it(nt-1))*derz(nav%dcz,x(it(nt-1))))
!       elseif(nav%tou==3) then
!          f=f+3._rk*(&
!               nav%u(it(nt))*derx(nav%dcx,x(it(nt)))+&
!               nav%v(it(nt))*dery(nav%dcy,x(it(nt)))+&
!               nav%w(it(nt))*derz(nav%dcz,x(it(nt))))

!          f=f-3._rk*(&
!               nav%u(it(nt-1))*derx(nav%dcx,x(it(nt-1)))+&
!               nav%v(it(nt-1))*dery(nav%dcy,x(it(nt-1)))+&
!               nav%w(it(nt-1))*derz(nav%dcz,x(it(nt-1))))

!          f=f+1._rk*(&
!               nav%u(it(nt-2))*derx(nav%dcx,x(it(nt-2)))+&
!               nav%v(it(nt-2))*dery(nav%dcy,x(it(nt-2)))+&
!               nav%w(it(nt-2))*derz(nav%dcz,x(it(nt-2))))
!       endif
!    elseif (nav%nlt==2) then
!       if (nav%tou==2) then
!          f=f+1._rk*(&
!               nav%u(it(nt))*derx(nav%dcx,x(it(nt)))+&
!               nav%v(it(nt))*dery(nav%dcy,x(it(nt)))+&
!               nav%w(it(nt))*derz(nav%dcz,x(it(nt))))

!          nav%aux=1._rk*x(it(nt))*nav%u(it(nt))
!          f=f+derx(nav%dcx,nav%aux)
!          nav%aux=1._rk*x(it(nt))*nav%v(it(nt))
!          f=f+dery(nav%dcy,nav%aux)
!          nav%aux=1._rk*x(it(nt))*nav%w(it(nt))
!          f=f+derz(nav%dcz,nav%aux)
!          
!          f=f-0.5_rk*(&
!               nav%u(it(nt-1))*derx(nav%dcx,x(it(nt-1)))+&
!               nav%v(it(nt-1))*dery(nav%dcy,x(it(nt-1)))+&
!               nav%w(it(nt-1))*derz(nav%dcz,x(it(nt-1))))
!          
!          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%u(it(nt-1))
!          f=f+derx(nav%dcx,nav%aux)
!          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%v(it(nt-1))
!          f=f+dery(nav%dcy,nav%aux)
!          nav%aux=(-0.5_rk)*x(it(nt-1))*nav%w(it(nt-1))
!          f=f+derz(nav%dcz,nav%aux)
!       elseif(nav%tou==3) then
!          f=f+1.5_rk*(&
!               nav%u(it(nt))*derx(nav%dcx,x(it(nt)))+&
!               nav%v(it(nt))*dery(nav%dcy,x(it(nt)))+&
!               nav%w(it(nt))*derz(nav%dcz,x(it(nt))))

!          nav%aux=1.5_rk*x(it(nt))*nav%u(it(nt))
!          f=f+derx(nav%dcx,nav%aux)
!          nav%aux=1.5_rk*x(it(nt))*nav%v(it(nt))
!          f=f+dery(nav%dcy,nav%aux)
!          nav%aux=1.5_rk*x(it(nt))*nav%w(it(nt))
!          f=f+derz(nav%dcz,nav%aux)

!          f=f-1.5_rk*(&
!               nav%u(it(nt-1))*derx(nav%dcx,x(it(nt-1)))+&
!               nav%v(it(nt-1))*dery(nav%dcy,x(it(nt-1)))+&
!               nav%w(it(nt-1))*derz(nav%dcz,x(it(nt-1))))
!          
!          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%u(it(nt-1))
!          f=f+derx(nav%dcx,nav%aux)
!          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%v(it(nt-1))
!          f=f+dery(nav%dcy,nav%aux)
!          nav%aux=(-1.5_rk)*x(it(nt-1))*nav%w(it(nt-1))
!          f=f+derz(nav%dcz,nav%aux)

!          f=f+0.5_rk*(&
!               nav%u(it(nt-2))*derx(nav%dcx,x(it(nt-2)))+&
!               nav%v(it(nt-2))*dery(nav%dcy,x(it(nt-2)))+&
!               nav%w(it(nt-2))*derz(nav%dcz,x(it(nt-2))))

!          nav%aux=(0.5_rk)*x(it(nt-2))*nav%u(it(nt-2))
!          f=f+derx(nav%dcx,nav%aux)
!          nav%aux=(0.5_rk)*x(it(nt-2))*nav%v(it(nt-2))
!          f=f+dery(nav%dcy,nav%aux)
!          nav%aux=(0.5_rk)*x(it(nt-2))*nav%w(it(nt-2))
!          f=f+derz(nav%dcz,nav%aux)
!       endif
!    endif

  end subroutine navier_nonlinear

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
    nav%rhs_u(1)=f(nav,'rhsu') - nav%rhs_u(1)

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
    nav%rhs_v(1)=f(nav,'rhsv') - nav%rhs_v(1)

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
    nav%rhs_w(1)=f(nav,'rhsw') -  nav%rhs_w(1)

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
          nav%fu(it(nt))=nav%fu(it(nt))+derx(nav%dcx,navier_extrapol(nav,nav%p,'p'))
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
    call multidomain_solve(mpid,nav%infu,nav%scu,nav%bcu(it(1)),nav%u(it(1)),&
          nav%fu(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz ,&
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
    nav%fv(it(nt))=nav%fv(it(nt))+nav%fac(2)*nav%v(nav%it(nav%nt  )) &
                                 +nav%fac(3)*nav%v(nav%it(nav%nt-1))
    if(nav%tou==3) nav%fv(it(nt))=nav%fv(it(nt))+nav%fac(4)*nav%v(nav%it(nav%nt-2))
    
    !-> pressure 
    if (nav%pt>=2) then
       if (nav%pt==2) then
          nav%fv(it(nt))=nav%fv(it(nt))+dery(nav%dcy,navier_extrapol(nav,nav%p,'p'))
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
    call multidomain_solve(mpid,nav%infv,nav%scv,nav%bcv(it(1)),nav%v(it(1)),&
          nav%fv(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz ,&
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
    nav%fw(it(nt))=nav%fw(it(nt))+nav%fac(2)*nav%w(nav%it(nav%nt  )) &
                                 +nav%fac(3)*nav%w(nav%it(nav%nt-1))
    if(nav%tou==3) nav%fw(it(nt))=nav%fw(it(nt))+nav%fac(4)*nav%w(nav%it(nav%nt-2))
    
    !-> pressure 
    if (nav%pt>=2) then
       if (nav%pt==2) then
          nav%fw(it(nt))=nav%fw(it(nt))+derz(nav%dcz,navier_extrapol(nav,nav%p,'p'))
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
    call multidomain_solve(mpid,nav%infw,nav%scw,nav%bcw(it(1)),nav%w(it(1)),&
          nav%fw(it(nt)),nav%aux,nav%sigmau,nav%dcx,nav%dcy,nav%dcz ,&
          inf_sol=nav%infsolw)

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
    integer(ik) :: it(nav%nt),nt
    real(rk) :: fac
    
    !-> put nav%nt in nt for ease of use
    nt=nav%nt
    it(:)=nav%it(:)

    !--------------------------------------------------------------------
    !-> compute rhs
    nav%fphi=nav%sigmap*navier_extrapol(nav,nav%phi,'p')
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
    call multidomain_solve(mpid,nav%infp,nav%scp,nav%bcphi(it(1)),nav%phi(it(1)),&
          nav%fphi,nav%aux,nav%sigmap,nav%dcx,nav%dcy,nav%dcz,null=nullv,&
          inf_sol=nav%infsolphi,var='p')

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

    call field_zero_edges(nav%phi(it(1)))
    !-> pressure
    if(nav%pt==2.or.nav%pt==4) nav%aux=navier_extrapol(nav,nav%p,'p')
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
      nav%aux=(derx(nav%dcx,navier_extrapol(nav,nav%u,'p'))+&
               dery(nav%dcy,navier_extrapol(nav,nav%v,'p'))+&
               derz(nav%dcz,navier_extrapol(nav,nav%w,'p')))/nav%rey
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


  function navier_extrapol(nav,var,type)
    implicit none
    type(navier3d),intent(in) :: nav
    type(field),intent(in) :: var(nav%nt)
    type(field)           ::navier_extrapol
    character(*)          :: type

    call field_init(navier_extrapol,"extrapol",nav%nx,nav%ny,nav%nz)

    navier_extrapol%f=0._rk
    if(nav%subite==1) then
      if(type=='p')then
				if (nav%top==2) then
				   navier_extrapol=1._rk*var(nav%it(nav%nt))
				elseif (nav%top==3) then
				   navier_extrapol=2._rk*var(nav%it(nav%nt))&
                          -1._rk*var(nav%it(nav%nt-1))
				endif
      elseif(type=='v')then
				if (nav%tou==2) then
				   navier_extrapol=2._rk*var(nav%it(nav%nt))&
                          -1._rk*var(nav%it(nav%nt-1))
        elseif(nav%tou==3) then
				   navier_extrapol=3._rk*var(nav%it(nav%nt))&
				                  -3._rk*var(nav%it(nav%nt-1))&
                          +1._rk*var(nav%it(nav%nt-2))
				endif
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
    nav%nsubite=cmd%nsubite

    !-> reynolds number 
    nav%rey=cmd%reynolds

    !--------------------------------------------------------------------
    !-> initialize mesh
    call mesh_init(nav%gridx,'gridx','x',nx,1,1)
    call mesh_init(nav%gridy,'gridy','y',1,ny,1)
    call mesh_init(nav%gridz,'gridz','z',1,1,nz)

    !-> initialize grid
    call mesh_grid_init(nav%gridx,'x',nx,1,1,mpid)
    call mesh_grid_init(nav%gridy,'y',nx,ny,1,mpid)
    call mesh_grid_init(nav%gridz,'z',1,1,nz,mpid)

    if (nav%tou==2) then
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

    !--------------------------------------------------------------------
    !-> initialize poisson solver coefficient for u
    bctu=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infu,bctu)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scu,bctu)

    !-> initialize poisson solver coefficient for v
    bctv=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infv,bctv)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scv,bctv)

    !-> initialize poisson solver coefficient for w
    bctw=(/1,1,1,1,1,1/)
    call md_boundary_condition_init(mpid,nav%infw,bctw)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scw,bctw)

    !-> initialize poisson solver coefficient for pressure
    bctp=(/2,2,2,2,2,2/)
    call md_boundary_condition_init(mpid,nav%infp,bctp)
    call solver_init_3d(nav%gridx,nav%gridy,nav%gridz,nav%scp,bctp)

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
       call field_init(nav%phi(i),"PHI",nx,ny,nz)
    enddo
    do i=1,2
       call field_init(nav%rhs_u(i),"RHS2_U",nx,ny,nz)
       call field_init(nav%rhs_v(i),"RHS2_V",nx,ny,nz)
       call field_init(nav%rhs_w(i),"RHS2_W",nx,ny,nz)
    enddo
    call field_init(nav%fphi,"RHS_PHI",nx,ny,nz)
    call field_init(nav%aux,"AUX",nx,ny,nz)
    call field_init(nav%rhs_px,"RHS_PX",nx,ny,nz)
    call field_init(nav%rhs_py,"RHS_PY",nx,ny,nz)
    call field_init(nav%rhs_pz,"RHS_PZ",nx,ny,nz)
    call field_init(nav%sub_u,"SUB_U",nx,ny,nz)
    call field_init(nav%sub_v,"SUB_V",nx,ny,nz)
    call field_init(nav%sub_w,"SUB_W",nx,ny,nz)
    call field_init(nav%sub_p,"SUB_P",nx,ny,nz)

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
    call derivatives_coefficients_init(nav%gridx,nav%dcx,nx,solver='yes')
    call derivatives_coefficients_init(nav%gridy,nav%dcy,ny,solver='yes')
    call derivatives_coefficients_init(nav%gridz,nav%dcz,nz,solver='yes')

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
    !-> initialize fields

    do i=1,nav%nt
       nav%u(i)%f=0._rk ; nav%v(i)%f=0._rk ; nav%w(i)%f=0._rk ; nav%p(i)%f=0._rk
       nav%fu(i)%f=0._rk ; nav%fv(i)%f=0._rk ; nav%fw(i)%f=0._rk ; nav%fp(i)%f=0._rk
       nav%phi(i)%f=0._rk 
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

    !--------------------------------------------------------------------
    !-> initialize guess sol
    call md_guess_init(mpid,nav%infu,nav%infsolu)
    call md_guess_init(mpid,nav%infv,nav%infsolv)
    call md_guess_init(mpid,nav%infw,nav%infsolw)
    call md_guess_init(mpid,nav%infp,nav%infsolphi)

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
       call field_destroy(nav%phi(i))
       call field_destroy(nav%fu(i))
       call field_destroy(nav%fv(i))
       call field_destroy(nav%fw(i))
       call field_destroy(nav%fp(i))
    enddo
    do i=1,2
       call field_destroy(nav%rhs_u(i))
       call field_destroy(nav%rhs_v(i))
       call field_destroy(nav%rhs_w(i))
    enddo
    call field_destroy(nav%fphi)
    call field_destroy(nav%aux)
    call field_destroy(nav%rhs_px)
    call field_destroy(nav%rhs_py)
    call field_destroy(nav%rhs_pz)
    call field_destroy(nav%sub_u)
    call field_destroy(nav%sub_v)
    call field_destroy(nav%sub_w)
    call field_destroy(nav%sub_p)

    !--------------------------------------------------------------------
    !-> finalize petsc
    call md_petsc_finalize()
    !-> finalize mpi
    call md_mpi_finalize(mpid)

  end subroutine navier_finalization

end module class_navier_3D
