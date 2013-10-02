module class_derivatives
! -----------------------------------------------------------------------
! Name :
! class derivatives 
! -----------------------------------------------------------------------
! Object : 
! Computation of first and second derivatives with an so-order compact 
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
     !-> order of spatial discretization
     integer(ik) :: so
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

  subroutine der_coeffs_init(grid,dc,n,solve,so,out_name)
! -----------------------------------------------------------------------
! Derivatives : Initialisation of derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2011
!
    implicit none
    type(derivatives_coefficients),intent(out) :: dc
    integer(ik),intent(in) :: n,so
    real(rk),intent(in) :: grid(n)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    logical :: solve

    dc%solver=solve
    dc%so=so
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
    call der_coeffs_generic_gobal(dc%n,dc%dl,dc%dr,grid,1,dc%so)
    call der_lu_factor(dc%n,dc%dl)
    if (present(out_name)) close(333)

    !-> compute coefficients for first derivatives solver type
    if (dc%solver) then
       if (present(out_name)) open(333,file='ds_'//out_name)
       dc%dl_s=0._rk ; dc%dr_s=0._rk
       call der_coeffs_generic_gobal(dc%n_s,dc%dl_s,dc%dr_s,grid(2),1,dc%so)
       call der_lu_factor(dc%n_s,dc%dl_s)
       if (present(out_name)) close(333)
    endif

    !-> compute coefficients for second derivatives
    if (present(out_name)) open(333,file='dd_'//out_name)
    dc%ddl=0._rk ; dc%ddr=0._rk
    call der_coeffs_generic_gobal(dc%n,dc%ddl,dc%ddr,grid,1,dc%so)
    call der_lu_factor(dc%n,dc%ddl)
    if (present(out_name)) close(333)

    !-> compute coefficients for second derivatives solver type
    if (dc%solver) then
       if (present(out_name)) open(333,file='dds_'//out_name)
       dc%ddl_s=0._rk ; dc%ddr_s=0._rk
       call der_coeffs_generic_gobal(dc%n_s,dc%ddl_s,dc%ddr_s,grid(2),1,dc%so)
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
 
  subroutine der_coeffs_generic_gobal(n,dl,dr,grid,der,so)
! -----------------------------------------------------------------------
! Derivatives : computation of coefficients for first derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer(ik),intent(in) :: n,der,so
!  der : order of derivation
    integer(ik) :: i,j,l
    real(rk),intent(out) :: dl(n,5),dr(n,7)
    real(rk),intent(in) :: grid(n)
    real(rk) :: h(7),coef_exp(7),coef_imp(5)
    dl=0._rk
    dr=0._rk

    h=grid(1)-grid(1:7) ! point on the edge
    call der_coeffs_generic(coef_imp,coef_exp,h,1,der,so)
    dl(1,:)=CSHIFT( coef_imp, shift=-2 )
    dr(1,:)=-coef_exp(:)

    h=grid(2)-grid(1:7) ! first point after the edge
    call der_coeffs_generic(coef_imp,coef_exp,h,2,der,so)
    dl(2,:)=CSHIFT( coef_imp, shift=-1 )
    dr(2,:)=-coef_exp(:)

    h=0._rk
    do i=3,n-2
      h(1:5)=grid(i)-grid(i-2:i+2) ! every other points
      call der_coeffs_generic(coef_imp,coef_exp,h,3,der,so)
      dl(i,:)=coef_imp
      dr(i,:)=-coef_exp
    enddo

    h=grid(n-1)-grid(n:n-6:-1) ! last point before the edge
    call der_coeffs_generic(coef_imp,coef_exp,h,2,der,so)
    dl(n-1,:)=CSHIFT( coef_imp(5:1:-1), shift=1 )
    dr(n-1,:)=-coef_exp(7:1:-1)

    h=grid(n)-grid(n:n-6:-1) ! point on the edge
    call der_coeffs_generic(coef_imp,coef_exp,h,1,der,so)
    dl(n,:)=CSHIFT( coef_imp(5:1:-1), shift=2 )
    dr(n,:)=-coef_exp(7:1:-1)

    if (der==2)  dr=-dr

  end subroutine der_coeffs_generic_gobal

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
