module class_solver_1d
! -----------------------------------------------------------------------
! Name :
! class solver 1d
! -----------------------------------------------------------------------
! Object : 
! Solve poisson and helmholtz equation with a 5-5 compact finite 
! differences scheme on irregular grid
! -----------------------------------------------------------------------
! Files :
! class_solver_1d.f90
! class_solver_coefficient_1d.f90
! -----------------------------------------------------------------------
! Public type :
! solver_coefficients -> coefficients for solver
! -----------------------------------------------------------------------
! Public Procedure :
! solver_coefficients_init -> initialize solver coefficients 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
!  use precision
  use class_solver_coefficient_1d
  implicit none

  !-> Define type containing all derivatives coefficients
  type,public :: solver_coefficients
     private
     integer(ik) :: nall,n
     real(rk),allocatable :: ddl(:,:),ddr(:,:)
     real(rk) :: bc(4)
  end type solver_coefficients
 
  !-> Declare everything private by default
  private
  
  !-> Declare exported procedure
  public :: solver_coeffs_init_poisson
  public :: solve_poisson_1d

contains

  subroutine print(n,ddl,ddr)
! -----------------------------------------------------------------------
! solver : print derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer,intent(in) :: n
    real(rk),intent(in) :: ddl(n,5),ddr(n,5)
    integer(ik) :: i,j

    print*,'ddl'
    do i=1,n
       print'(5es17.8)',(ddl(i,j),j=1,5)
    enddo
    print*,'ddr'
    do i=1,n
       print'(5es17.8)',(ddr(i,j),j=1,5)
    enddo

  end subroutine print

  subroutine solver_coeffs_init_poisson(grid,dc,n,out_name)
! -----------------------------------------------------------------------
! solver : Initialisation of derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    type(solver_coefficients),intent(out) :: dc
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: grid(n)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg

    !-> initialize and allocate dc

    if (.not.allocated(dc%ddl)) then
       dc%nall=n
       dc%n=n-2 
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,5))
    elseif (allocated(dc%ddl).and.n/=dc%nall) then
       deallocate(dc%ddl,dc%ddr)
       dc%nall=n 
       dc%n=n-2        
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,5))
    endif

    !-> compute coefficients for second derivatives

    if (present(out_name)) open(333,file='sd1d_'//out_name)
    dc%ddl=0._rk ; dc%ddr=0._rk
    call dder_coeffs(dc%n,dc%ddl,dc%ddr,grid)
    call print(dc%n,dc%ddl,dc%ddr)

    call dder_coeffs_poisson(dc%n,dc%ddl,dc%ddr)
    
    !-> boundary coefficients
    dc%bc(1)=dc%ddl(1,2) ; dc%bc(2)=dc%ddl(2,1)
    dc%bc(3)=dc%ddl(dc%n-1,5) ; dc%bc(4)=dc%ddl(dc%n,4) 
    
    !-> factorize linear system
!    call print(dc%n,dc%ddl,dc%ddr)
    call der_lu_factor(dc%n,dc%ddl)
!    call print(dc%n,dc%ddl,dc%ddr)
    if (present(out_name)) close(333)

  end subroutine solver_coeffs_init_poisson

  subroutine solve_poisson_1d(dc,x,sol)
! -----------------------------------------------------------------------
! solver : compute second derivative 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    type(solver_coefficients),intent(in) :: dc
    real(rk),intent(inout) :: sol(dc%n+2)
    real(rk),intent(in) :: x(dc%n+2)
    integer(ik) :: i
    real(rk) :: rhs(dc%n)

    !-> compute rhs
    call solve_poisson_rhs(dc%n,dc%ddr,x,rhs)

    !-> boundary condition
    rhs(1)=rhs(1)-dc%bc(1)*sol(1)
    rhs(2)=rhs(2)-dc%bc(2)*sol(1)
    rhs(dc%n-1)=rhs(dc%n-1)-dc%bc(3)*sol(dc%n+2)
    rhs(dc%n)=rhs(dc%n)-dc%bc(4)*sol(dc%n+2)

    !-> solve linear system
    call der_lu_solve(dc%n,dc%ddl,dc%ddr,sol(2),rhs)

  end subroutine solve_poisson_1d

  subroutine solve_poisson_rhs(n,dr,f,rhs)
! -----------------------------------------------------------------------
! solver : computation of rhs of compact finite difference scheme
!          for poisson equation
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
! TODO : possible optimization using intrinsic dot_product
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(in) :: dr(n,5),f(0:n+1)
    real(rk),intent(out) :: rhs(n)

    rhs(1)=dr(1,1)*f(1)+dr(1,2)*f(2)+dr(1,3)*f(3)+dr(1,4)*f(4)&
         +dr(1,5)*f(5)

    rhs(2)=dr(2,1)*f(1)+dr(2,2)*f(2)+dr(2,3)*f(3)+dr(2,4)*f(4) &
         +dr(2,5)*f(5)

    do i=3,n-2
       rhs(i)=dr(i,1)*f(i-2)+dr(i,2)*f(i-1)+dr(i,3)*f(i)+dr(i,4)*f(i+1) &
            +dr(i,5)*f(i+2)
    enddo

    rhs(n-1)=dr(n-1,1)*f(n-4)+dr(n-1,2)*f(n-3)+dr(n-1,3)*f(n-2) &
         +dr(n-1,4)*f(n-1)+dr(n-1,5)*f(n)

    rhs(n)=dr(n,1)*f(n-4)+dr(n,2)*f(n-3)+dr(n,3)*f(n-2) &
         +dr(n,4)*f(n-1)+dr(n,5)*f(n)

  end subroutine solve_poisson_rhs

  subroutine der_lu_factor(n,dl)
! -----------------------------------------------------------------------
! solver : LU factorization
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer(ik),intent(in) :: n
    real(rk),intent(inout) :: dl(n,5)
    
    call pentalu(n,dl(1,1),dl(1,2),dl(1,3),dl(1,4),dl(1,5), &
         dl(1,1),dl(n,5))

  end subroutine der_lu_factor

  subroutine der_lu_solve(n,dl,dr,x,b)
! -----------------------------------------------------------------------
! solver : LU solve
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: dl(n,5),dr(n,5),b(n)    
    real(rk),intent(out) :: x(n)

    call pentares(n,dl(1,1),dl(1,2),dl(1,3),dl(1,4),dl(1,5),dl(1,1),&
         dl(n,5),x,b)

  end subroutine der_lu_solve
 
  subroutine dder_coeffs_poisson(n,ddl,ddr)
! -----------------------------------------------------------------------
! solver : rearrange ddl and ddr coefficients for poisson solver
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,5)
    real(rk) :: tmp(n,5)
   
    !-> save ddl coeffs in temporary variable 
    tmp=ddl

    !-> create lhs of poisson equation from ddr coeffs

    ddl=ddr
    
    ddl(1,1)=ddr(1,5) ; ddl(1,2)=ddr(1,1) ; ddl(1,3)=ddr(1,2) 
    ddl(1,4)=ddr(1,3) ; ddl(1,5)=ddr(1,4)

!    ddl(2,1)=ddr(2,5) ; ddl(2,2)=ddr(2,1) ; ddl(2,3)=ddr(2,2) 
!    ddl(2,4)=ddr(2,3) ; ddl(2,5)=ddr(2,4)

!    ddl(3:n-2,:)=ddr(3:n-2,:)
    ddl(2:n-1,:)=ddr(2:n-1,:)

!    ddl(n-1,5)=ddr(n-1,1) ; ddl(n-1,4)=ddr(n-1,5) ; ddl(n-1,3)=ddr(n-1,4) 
!    ddl(n-1,2)=ddr(n-1,3) ; ddl(n-1,1)=ddr(n-1,2)

    ddl(n,5)=ddr(n,1) ; ddl(n,4)=ddr(n,5) ; ddl(n,3)=ddr(n,4) 
    ddl(n,2)=ddr(n,3) ; ddl(n,1)=ddr(n,2)

!    ddl=ddr

    !-> create rhs of poisson equation from ddl coeffs

    ddr(1,1)=tmp(1,3) ; ddr(1,2)=tmp(1,4) ; ddr(1,3)=tmp(1,5)
    ddr(1,4)=tmp(1,1) ; ddr(1,5)=tmp(1,2)

    ddr(2,1)=tmp(2,2) ; ddr(2,2)=tmp(2,3) ; ddr(2,3)=tmp(2,4)
    ddr(2,4)=tmp(2,5) ; ddr(2,5)=tmp(2,1)

    ddr(3:n-2,:)=tmp(3:n-2,:)

    ddr(n-1,1)=tmp(n-1,5) ; ddr(n-1,2)=tmp(n-1,1) ; ddr(n-1,3)=tmp(n-1,2)
    ddr(n-1,4)=tmp(n-1,3) ; ddr(n-1,5)=tmp(n-1,4) 

    ddr(n,1)=tmp(n,4) ; ddr(n,2)=tmp(n,5) ; ddr(n,3)=tmp(n,1)
    ddr(n,4)=tmp(n,2) ; ddr(n,5)=tmp(n,3) 
    
!    ddr=tmp

  end subroutine dder_coeffs_poisson

  subroutine dder_coeffs(n,ddl,ddr,grid)
! -----------------------------------------------------------------------
! solver : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,5)
    real(rk),intent(in) :: grid(0:n+1)

    ! alpha,beta : central coefficients lhs
    ! a,b : central coefficients rhs
    real(rk) :: alp,bet,gam,del
    real(rk) :: a,b,c,d,e
    ! alp1,bet1,gam1,del1 : 1st points coefficients lhs
    ! alp2,bet2,gam2,del2 : 2st points coefficients lhs
    real(rk) :: alp1,bet1,gam1,del1
    real(rk) :: alp2,bet2,gam2,del2
    ! b1,c1,d1,e1 : 1st points coefficients rhs
    ! b2,c2,d2,e2 : 2st points coefficients lhs
    real(rk) :: a1,b1,c1,d1,e1
    real(rk) :: a2,b2,c2,d2,e2
    ! h1,h2,h3,h4 : grid intervals
    real(rk) :: h1,h2,h3,h4

    ddl=0._rk ; ddr=0._rk

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,alp1,bet1,gam1,del1,h1,h2,h3,h4)
    ddl(1,1)=del1 ; ddl(1,2)=0._rk ; ddl(1,3)=alp1 ; ddl(1,4)=bet1 ; ddl(1,5)=gam1
    ddr(1,1)=a1 ; ddr(1,2)=b1 ; ddr(1,3)=c1 ; ddr(1,4)=d1 ; ddr(1,5)=e1

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    call dder_coeffs_i3(a2,b2,c2,d2,e2,alp2,bet2,gam2,del2,h1,h2,h3,h4)
    ddl(2,1)=0._rk ; ddl(2,2)=alp2 ; ddl(2,3)=bet2 ; ddl(2,4)=gam2 ; ddl(2,5)=del2
    ddr(2,1)=a2 ; ddr(2,2)=b2 ; ddr(2,3)=c2 ; ddr(2,4)=d2 ; ddr(2,5)=e2

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
    call dder_coeffs_i3(a2,b2,c2,d2,e2,alp2,bet2,gam2,del2,h1,h2,h3,h4)
    ddl(n-1,1)=del2 ; ddl(n-1,2)=gam2 ; ddl(n-1,3)=bet2 ; ddl(n-1,4)=alp2 ; ddl(n-1,5)=0._rk
    ddr(n-1,1)=e2 ; ddr(n-1,2)=d2 ; ddr(n-1,3)=c2 ; ddr(n-1,4)=b2 ; ddr(n-1,5)=a2

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,alp1,bet1,gam1,del1,h1,h2,h3,h4)
    ddl(n,1)=gam1 ; ddl(n,2)=bet1 ; ddl(n,3)=alp1 ; ddl(n,4)=0._rk ; ddl(n,5)=del1
    ddr(n,1)=e1 ; ddr(n,2)=d1 ; ddr(n,3)=c1 ; ddr(n,4)=b1 ; ddr(n,5)=a1

  end subroutine dder_coeffs

  subroutine pentalu(n,p,f,c,d,e,alpha,beta)
! -----------------------------------------------------------------------
! solver : Resolution of [a](x)=(b) with [a] quasi-pentadiagonal
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
! solver : Resolution of [a](x)=(b) with [a] quasi-pentadiagonal
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

end module class_solver_1d
