module class_derivatives
! -----------------------------------------------------------------------
! Name :
! class derivatives 
! -----------------------------------------------------------------------
! Object : 
! Computation of first and second derivatives with an 8-order compact 
! finite difference scheme of maximum stencil 5-7 on irregular grid
! -----------------------------------------------------------------------
! Files :
! class_derivatives.f90
! class_derivatives_coefficient.f90
! -----------------------------------------------------------------------
! Public type :
! derivatives_coefficients -> coefficients for derivatives
! -----------------------------------------------------------------------
! Public Procedure :
! derivatives_coefficients_init -> initialize derivatives coefficients 
! derx,dery,derz -> compute first derivatives in x,y,z directions
! dderx,ddery,dderz -> compute second derivatives in x,y,z directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
!  use precision
  use class_derivatives_coefficient
  implicit none

  !-> Define type containing all derivatives coefficients
  type,public :: derivatives_coefficients
     private
     !-> dimension
     integer(ik) :: n
     !-> coefficient
     real(rk),allocatable :: dl(:,:),ddl(:,:),dr(:,:),ddr(:,:)
     !-> solver derivatives or not
     logical :: solver
     !-> dimension for solver type
     integer(ik) :: n_s
     !-> coefficient for solver type
     real(rk),allocatable :: dl_s(:,:),ddl_s(:,:),dr_s(:,:),ddr_s(:,:)
!     real(rk) :: dl(n,5),ddl(n,5),dr(n,7),ddr(nx,7)
  end type derivatives_coefficients
 
  !-> Declare everything private by default
  private
  
  !-> Declare exported procedure
  public :: der,dder
  public :: der_s,dder_s
  public :: der_coeffs_init
  public :: der_solver_init_type,dertype

contains

  subroutine print(n,dl,dr)
! -----------------------------------------------------------------------
! Derivatives : print derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2010
!
    implicit none
    integer,intent(in) :: n
    real(rk),intent(in) :: dl(n,5),dr(n,7)
    integer(ik) :: i,j

    print*,'dl'
    do i=1,n
       print'(5es17.8)',(dl(i,j),j=1,5)
    enddo
    print*,'dr'
    do i=1,n
       print'(7es17.8)',(dr(i,j),j=1,7)
    enddo

  end subroutine print

  subroutine der_solver_init_type(dc,solver_type)
! -----------------------------------------------------------------------
! Derivatives : initialize solver type derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
    implicit none
    type(derivatives_coefficients),intent(out) :: dc
    logical :: solver_type
    
    dc%solver=solver_type

  end subroutine der_solver_init_type

  function dertype(dc)
! -----------------------------------------------------------------------
! Derivatives : initialize solver type derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
    implicit none
    type(derivatives_coefficients),intent(in) :: dc
    logical :: dertype
    
    dertype=dc%solver
    
  end function dertype

  subroutine der_coeffs_init(grid,dc,n,solve,out_name)
! -----------------------------------------------------------------------
! Derivatives : Initialisation of derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    implicit none
    type(derivatives_coefficients),intent(out) :: dc
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: grid(n)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    logical :: solve

    dc%solver=solve
    !-> initialize and allocate dc

    !-> coefficients for n
    if (.not.allocated(dc%dl)) then
       dc%n=n 
       allocate(dc%dl(n,5),dc%ddl(n,5),dc%dr(n,7),dc%ddr(n,7))
    elseif (allocated(dc%dl).and.n/=dc%n) then
       deallocate(dc%dl,dc%ddl,dc%dr,dc%ddr)
       dc%n=n 
       allocate(dc%dl(n,5),dc%ddl(n,5),dc%dr(n,7),dc%ddr(n,7))
    endif

    !-> coefficients for solver type : n_s
    if (dc%solver) then
       if (.not.allocated(dc%dl_s)) then
          dc%n_s=n-2 
          allocate(dc%dl_s(n-2,5),dc%ddl_s(n-2,5),dc%dr_s(n-2,7),dc%ddr_s(n-2,7))
       elseif (allocated(dc%dl_s).and.dc%n_s/=n-2) then
          deallocate(dc%dl_s,dc%ddl_s,dc%dr_s,dc%ddr_s)
          dc%n_s=n-2 
          allocate(dc%dl_s(n-2,5),dc%ddl_s(n-2,5),dc%dr_s(n-2,7),dc%ddr_s(n-2,7))
       endif
    endif

    !-> compute coefficients for first derivatives
    if (present(out_name)) open(333,file='d_'//out_name)
    dc%dl=0._rk ; dc%dr=0._rk
    call der_coeffs(dc%n,dc%dl,dc%dr,grid)
    call der_lu_factor(dc%n,dc%dl)
    if (present(out_name)) close(333)

    !-> compute coefficients for first derivatives solver type
    if (dc%solver) then
       if (present(out_name)) open(333,file='ds_'//out_name)
       dc%dl_s=0._rk ; dc%dr_s=0._rk
       call der_coeffs(dc%n_s,dc%dl_s,dc%dr_s,grid(2))
       call der_lu_factor(dc%n_s,dc%dl_s)
       if (present(out_name)) close(333)
    endif

    !-> compute coefficients for second derivatives
    if (present(out_name)) open(333,file='dd_'//out_name)
    dc%ddl=0._rk ; dc%ddr=0._rk
    call dder_coeffs(dc%n,dc%ddl,dc%ddr,grid)
    call der_lu_factor(dc%n,dc%ddl)
    if (present(out_name)) close(333)

    !-> compute coefficients for second derivatives solver type
    if (dc%solver) then
       if (present(out_name)) open(333,file='dds_'//out_name)
       dc%ddl_s=0._rk ; dc%ddr_s=0._rk
       call dder_coeffs(dc%n_s,dc%ddl_s,dc%ddr_s,grid(2))
       call der_lu_factor(dc%n_s,dc%ddl_s)
       if (present(out_name)) close(333)
    endif

  end subroutine der_coeffs_init

  subroutine der(dc,x,dx)
! -----------------------------------------------------------------------
! Derivatives : compute first derivative 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    implicit none
    type(derivatives_coefficients),intent(in) :: dc
    real(rk),intent(out) :: dx(:)
    real(rk),intent(in) :: x(:)
    integer(ik) :: i
    real(rk) :: rhs(dc%n)

    call derx_rhs(dc%n,dc%dr,x,rhs)
    call der_lu_solve(dc%n,dc%dl,dc%dr,dx,rhs)

  end subroutine der

  subroutine der_s(dc,x,dx)
! -----------------------------------------------------------------------
! Derivatives : compute first derivative 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
    implicit none
    type(derivatives_coefficients),intent(in) :: dc
    real(rk),intent(out) :: dx(dc%n)
    real(rk),intent(in) :: x(dc%n)
    integer(ik) :: i
    real(rk) :: rhs(dc%n)

    call derx_rhs(dc%n_s,dc%dr_s,x(2),rhs(2))
    call der_lu_solve(dc%n_s,dc%dl_s,dc%dr_s,dx(2),rhs(2))

  end subroutine der_s

  subroutine dder(dc,x,ddx)
! -----------------------------------------------------------------------
! Derivatives : compute second derivative 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    implicit none
    type(derivatives_coefficients),intent(in) :: dc
    real(rk),intent(out) :: ddx(:)
    real(rk),intent(in) :: x(:)
    integer(ik) :: i
    real(rk) :: rhs(dc%n)

    call derx_rhs(dc%n,dc%ddr,x,rhs)
    call der_lu_solve(dc%n,dc%ddl,dc%ddr,ddx,rhs)

  end subroutine dder

  subroutine dder_s(dc,x,ddx)
! -----------------------------------------------------------------------
! Derivatives : compute second derivative 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    implicit none
    type(derivatives_coefficients),intent(in) :: dc
    real(rk),intent(out) :: ddx(dc%n)
    real(rk),intent(in) :: x(dc%n)
    integer(ik) :: i
    real(rk) :: rhs(dc%n)

    call derx_rhs(dc%n_s,dc%ddr_s,x(2),rhs(2))
    call der_lu_solve(dc%n_s,dc%ddl_s,dc%ddr_s,ddx(2),rhs(2))

  end subroutine dder_s

  subroutine derx_rhs(n1,dr,f,rhs)
! -----------------------------------------------------------------------
! Derivatives : computation of rhs of compact finite difference scheme
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
! TODO : possible optimization using intrinsic dot_product
!
    implicit none
    integer(ik),intent(in) :: n1
    integer(ik) :: i
    real(rk),intent(in) :: dr(n1,7),f(n1)
    real(rk),intent(out) :: rhs(n1)

    rhs(1)=dr(1,1)*f(1)+dr(1,2)*f(2)+dr(1,3)*f(3)+dr(1,4)*f(4)&
         +dr(1,5)*f(5)+dr(1,6)*f(6)+dr(1,7)*f(7)

    rhs(2)=dr(2,1)*f(1)+dr(2,2)*f(2)+dr(2,3)*f(3)+dr(2,4)*f(4) &
         +dr(2,5)*f(5)+dr(2,6)*f(6)+dr(2,7)*f(7)

    do i=3,n1-2
        rhs(i)=dot_product(dr(i,1:5),f(i-2:i+2))
!       rhs(i)=dr(i,1)*f(i-2)+dr(i,2)*f(i-1)+dr(i,3)*f(i)+dr(i,4)*f(i+1) &
!            +dr(i,5)*f(i+2)
    enddo

    rhs(n1-1)=dr(n1-1,1)*f(n1-6)+dr(n1-1,2)*f(n1-5)+dr(n1-1,3)*f(n1-4) &
         +dr(n1-1,4)*f(n1-3)+dr(n1-1,5)*f(n1-2)+dr(n1-1,6)*f(n1-1)+dr(n1-1,7)*f(n1)

    rhs(n1)=dr(n1,1)*f(n1-6)+dr(n1,2)*f(n1-5)+dr(n1,3)*f(n1-4) &
         +dr(n1,4)*f(n1-3)+dr(n1,5)*f(n1-2)+dr(n1,6)*f(n1-1)+dr(n1,7)*f(n1)

  end subroutine derx_rhs

  subroutine der_lu_factor(n,dl)
! -----------------------------------------------------------------------
! Derivatives : LU factorization
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer(ik),intent(in) :: n
    real(rk),intent(inout) :: dl(n,5)
    
    call pentalu(n,dl(1,1),dl(1,2),dl(1,3),dl(1,4),dl(1,5),dl(1,1),&
         dl(n,5))

  end subroutine der_lu_factor

  subroutine der_lu_solve(n,dl,dr,x,b)
! -----------------------------------------------------------------------
! Derivatives : LU solve
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: dl(n,5),dr(n,5),b(n)    
    real(rk),intent(out) :: x(n)

    call pentares(n,dl(1,1),dl(1,2),dl(1,3),dl(1,4),dl(1,5),dl(1,1),&
         dl(n,5),x,b)

  end subroutine der_lu_solve
 
  subroutine der_coeffs(n,dl,dr,grid)
! -----------------------------------------------------------------------
! Derivatives : computation of coefficients for first derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: dl(n,5),dr(n,7)
    real(rk),intent(in) :: grid(n)

    ! alpha,beta,gam,del : central coefficients lhs
    ! a,b,c,d,e : central coefficients rhs
    real(rk) :: alp,bet,gam,del
    real(rk) :: a,b,c,d,e
    ! alp1,bet1,gam1,del1 : 1st points coefficients lhs
    ! alp2,bet2,gam2,del2 : 2st points coefficients lhs
    real(rk) :: alp1,bet1,gam1,del1
    real(rk) :: alp2,bet2,gam2,del2
    ! b1,c1,d1,e1,f1,g1 : 1st points coefficients rhs
    ! b2,c2,d2,e2,f1,g1 : 2st points coefficients lhs
    real(rk) :: a1,b1,c1,d1,e1,f1,g1
    real(rk) :: a2,b2,c2,d2,e2,f2,g2
    ! h1,h2,h3,h4,h5,h6 : grid intervals
    real(rk) :: h1,h2,h3,h4,h5,h6

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call der_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6)
    dl(1,1)=del1 ; dl(1,2)=0._rk ; dl(1,3)=alp1 ; dl(1,4)=bet1 ; dl(1,5)=gam1
    dr(1,1)=a1 ; dr(1,2)=b1 ; dr(1,3)=c1 ; dr(1,4)=d1 ; dr(1,5)=e1 
    dr(1,6)=f1 ; dr(1,7)=g1 ; 

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call der_coeffs_i2(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6)
    dl(2,1)=0._rk ; dl(2,2)=alp2 ; dl(2,3)=bet2 ; dl(2,4)=gam2 ; dl(2,5)=del2
    dr(2,1)=a2 ; dr(2,2)=b2 ; dr(2,3)=c2 ; dr(2,4)=d2 ; dr(2,5)=e2 
    dr(2,6)=f2 ; dr(2,7)=g2 ; 

    do i=3,n-2
       h1=grid(i-1)-grid(i-2)
       h2=grid(i)-grid(i-1)
       h3=grid(i+1)-grid(i)
       h4=grid(i+2)-grid(i+1)
       call der_coeffs_c(a,b,c,d,e,alp,bet,gam,del,h1,h2,h3,h4)
       dl(i,1)=alp ; dl(i,2)=bet ; dl(i,3)=1._rk ; dl(i,4)=gam ; dl(i,5)=del
       dr(i,1)=a ; dr(i,2)=b ; dr(i,3)=c ; dr(i,4)=d ; dr(i,5)=e
    enddo

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call der_coeffs_i2(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6)
    dl(n-1,1)=del2 ; dl(n-1,2)=gam2 ; dl(n-1,3)=bet2 ; dl(n-1,4)=alp2 ; dl(n-1,5)=0._rk
    dr(n-1,1)=-g2 ; dr(n-1,2)=-f2
    dr(n-1,3)=-e2 ; dr(n-1,4)=-d2 ; dr(n-1,5)=-c2 ; dr(n-1,6)=-b2 ; dr(n-1,7)=-a2
 
    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call der_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6)
    dl(n,1)=gam1 ; dl(n,2)=bet1 ; dl(n,3)=alp1 ; dl(n,4)=0._rk ; dl(n,5)=del1
    dr(n,1)=-g1 ; dr(n,2)=-f1
    dr(n,3)=-e1 ; dr(n,4)=-d1 ; dr(n,5)=-c1 ; dr(n,6)=-b1 ; dr(n,7)=-a1

  end subroutine der_coeffs

  subroutine dder_coeffs(n,ddl,ddr,grid)
! -----------------------------------------------------------------------
! Derivatives : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2011
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,7)
    real(rk),intent(in) :: grid(n)

    ! alpha,beta : central coefficients lhs
    ! a,b : central coefficients rhs
    real(rk) :: alp,bet,gam,del
    real(rk) :: a,b,c,d,e
    ! alp1,bet1,gam1,del1 : 1st points coefficients lhs
    ! alp2,bet2,gam2,del2 : 2st points coefficients lhs
    real(rk) :: alp1,bet1,gam1,del1
    real(rk) :: alp2,bet2,gam2,del2
    ! b1,c1,d1,e1,f1,g1 : 1st points coefficients rhs
    ! b2,c2,d2,e2,f2,g2 : 2st points coefficients lhs
    real(rk) :: a1,b1,c1,d1,e1,f1,g1
    real(rk) :: a2,b2,c2,d2,e2,f2,g2
    ! h1,h2,h3,h4,h5,h6 : grid intervals
    real(rk) :: h1,h2,h3,h4,h5,h6

    f1=0._rk ; g1=0._rk ; f2=0._rk ; g2=0._rk

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call dder_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6)
    alp1=1._rk
    ddl(1,1)=del1 ; ddl(1,2)=0._rk ; ddl(1,3)=alp1 ; ddl(1,4)=bet1 ; ddl(1,5)=gam1
    ddr(1,1)=a1 ; ddr(1,2)=b1 ; ddr(1,3)=c1 ; ddr(1,4)=d1 ; ddr(1,5)=e1
    ddr(1,6)=f1 ; ddr(1,7)=g1

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call dder_coeffs_i2(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6)
    bet2=1._rk
    ddl(2,1)=0._rk ; ddl(2,2)=alp2 ; ddl(2,3)=bet2 ; ddl(2,4)=gam2 ; ddl(2,5)=del2
    ddr(2,1)=a2 ; ddr(2,2)=b2 ; ddr(2,3)=c2 ; ddr(2,4)=d2 ; ddr(2,5)=e2
    ddr(2,6)=f2 ; ddr(2,7)=g2

    do i=3,n-2
       h1=grid(i-1)-grid(i-2)
       h2=grid(i)-grid(i-1)
       h3=grid(i+1)-grid(i)
       h4=grid(i+2)-grid(i+1)
       call dder_coeffs_c(a,b,c,d,e,alp,bet,gam,del,h1,h2,h3,h4)
       ddl(i,1)=alp ; ddl(i,2)=bet ; ddl(i,3)=1._rk ; ddl(i,4)=gam ; ddl(i,5)=del
       ddr(i,1)=a ; ddr(i,2)=b ; ddr(i,3)=c ; ddr(i,4)=d ; ddr(i,5)=e
    enddo

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call dder_coeffs_i2(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6)
    bet2=1._rk
    ddl(n-1,1)=del2 ; ddl(n-1,2)=gam2 ; ddl(n-1,3)=bet2 ; ddl(n-1,4)=alp2 ; ddl(n-1,5)=0._rk
    ddr(n-1,1)=g2 ; ddr(n-1,2)=f2
    ddr(n-1,3)=e2 ; ddr(n-1,4)=d2 ; ddr(n-1,5)=c2 ; ddr(n-1,6)=b2 ; ddr(n-1,7)=a2

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call dder_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6)
    alp1=1._rk
    ddl(n,1)=gam1 ; ddl(n,2)=bet1 ; ddl(n,3)=alp1 ; ddl(n,4)=0._rk ; ddl(n,5)=del1
    ddr(n,1)=g1 ; ddr(n,2)=f1
    ddr(n,3)=e1 ; ddr(n,4)=d1 ; ddr(n,5)=c1 ; ddr(n,6)=b1 ; ddr(n,7)=a1

  end subroutine dder_coeffs

  subroutine pentalu(n,p,f,c,d,e,alpha,beta)
! -----------------------------------------------------------------------
! Derivatives : Resolution of [a](x)=(b) with [a] quasi-pentadiagonal
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2009
!
    implicit double precision (a-h,o-z)
    integer(ik) :: i
    integer(ik),intent(in) :: n
    real(rk),intent(inout) :: p(n),f(n),c(n),d(n),e(n)
    real(rk),intent(inout) :: alpha,beta

    !-> step 1

    f(2)=f(2)/c(1)
    c(2)=c(2)-d(1)*f(2)
    d(2)=d(2)-e(1)*f(2)
    e(2)=e(2)-alpha*f(2)

    !-> step 2

    i=3
    p(i)=p(i)/c(i-2)
    f(i)=(f(i)-p(i)*d(i-2))/c(i-1)
    c(i)=c(i)-d(i-1)*f(i)-e(i-2)*p(i)
    d(i)=d(i)-e(i-1)*f(i)-alpha*p(i)

    do i=4,n-1
       p(i)=p(i)/c(i-2)
       f(i)=(f(i)-p(i)*d(i-2))/c(i-1)
       c(i)=c(i)-e(i-2)*p(i)-d(i-1)*f(i)
       d(i)=d(i)-e(i-1)*f(i)
    enddo

    !-> step 3

    beta=beta/c(n-3)
    p(n)=p(n)-beta*d(n-3)
    f(n)=f(n)-beta*e(n-3)
    p(n)=p(n)/c(n-2)
    f(n)=(f(n)-p(n)*d(n-2))/c(n-1)
    c(n)=c(n)-e(n-2)*p(n)-d(n-1)*f(n)

  end subroutine pentalu

  subroutine pentares(n,p,f,c,d,e,alpha,beta,x,b)
! -----------------------------------------------------------------------
! Derivatives : Resolution of [a](x)=(b) with [a] quasi-pentadiagonal
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 09/2009
!
    implicit none
    integer(ik) :: i
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: p(n),f(n),c(n),d(n),e(n),alpha,beta,b(n)   
    real(rk),intent(out) :: x(n)

    !-> step 4

    x(1)=b(1)
    x(2)=b(2)-x(1)*f(2)

    !-> step 5

    do  i=3,n-1
       x(i)=b(i)-p(i)*x(i-2)-f(i)*x(i-1)
    enddo

    i=n
    x(i)=b(i)-p(i)*x(i-2)-f(i)*x(i-1)-beta*x(i-3)

    !-> step 6

    x(n)=x(n)/c(n)
    x(n-1)=(x(n-1)-d(n-1)*x(n))/c(n-1)

    !-> step 7

    do i=n-2,2,-1
       x(i)=(x(i)-e(i)*x(i+2)-d(i)*x(i+1))/c(i)
    enddo

    i=1
    x(i)=(x(i)-e(i)*x(i+2)-d(i)*x(i+1)-alpha*x(i+3))/c(i)

end subroutine pentares

end module class_derivatives
