module class_solver_coefficient_1d
! -----------------------------------------------------------------------
! Name :
! class solver 
! -----------------------------------------------------------------------
! Object : 
! Solve poisson and helmholtz equation with a 5-5 compact finite 
! differences scheme on irregular grid
! -----------------------------------------------------------------------
! Files :
! class_solver.f90
! class_solver_coefficient.f90
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
  use precision
  implicit none
!  private
  real(rk),private,save :: h(7)
  
contains
  
subroutine dder_coeffs_c(a1,b1,c1,d1,e1,alpha,beta,gamma,delta,h1,h2,h3,h4)
! -----------------------------------------------------------------------
! solver : Coefficients for second derivatives : non-uniform grid,centered 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
  implicit none
  real(rk),intent(out) :: a1,b1,c1,d1,e1,alpha,beta,gamma,delta
  real(rk),intent(in) :: h1,h2,h3,h4
  !
  ! n : size of linear system
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=9
  character(len=1) :: fact,trans,equed
  integer(ik) :: nrhs,lda,ldb,ldaf,ldx
  real(rk) :: a(n,n),af(n,n),ipiv(n),b(n,1),r(n),c(n),x(n),work(4*n)
  real(rk) :: rcond,ferr,berr,dum
  integer(ik) :: iwork(n),info

  !-> space interval
  h(1)=h1 ; h(2)=h2 ; h(3)=h3 ; h(4)=h4

  !-> define matrix
  a(1,1)=0;a(1,2)=0;a(1,3)=0;a(1,4)=0;a(1,5)=1;a(1,6)=1;a(1,7)=1;a(1,8)=1;a(1,9)=1
  
  a(2,1)=0;a(2,2)=0;a(2,3)=0;a(2,4)=0;a(2,5)=-h2-h1;a(2,6)=-h2;a(2,7)=0;a(2,8)=h3;a(2,9)=h4+h3
  
  a(3,1)=-1;a(3,2)=-1;a(3,3)=-1;a(3,4)=-1;a(3,5)=(h2+h1)**2/2.0;a(3,6)=h2**2/2.0
  a(3,7)=0;a(3,8)=h3**2/2.0;a(3,9)=(-h4-h3)**2/2.0
  
  a(4,1)=h2+h1;a(4,2)=h2;a(4,3)=-h3;a(4,4)=-h4-h3;a(4,5)=-(h2+h1)**3/6.0;a(4,6)=-h2**3/6.0
  a(4,7)=0;a(4,8)=h3**3/6.0;a(4,9)=-(-h4-h3)**3/6.0
  
  a(5,1)=-(h2+h1)**2/2.0;a(5,2)=-h2**2/2.0;a(5,3)=-h3**2/2.0;a(5,4)=-(-h4-h3)**2/2.0
  a(5,5)=(h2+h1)**4/24.0;a(5,6)=h2**4/24.0;a(5,7)=0;a(5,8)=h3**4/24.0;a(5,9)=(-h4-h3)**4/24.0
  
  a(6,1)=(h2+h1)**3/6.0;a(6,2)=h2**3/6.0;a(6,3)=-h3**3/6.0;a(6,4)=(-h4-h3)**3/6.0
  a(6,5)=-(h2+h1)**5/120.0;a(6,6)=-h2**5/120.0;a(6,7)=0;a(6,8)=h3**5/120.0;a(6,9)=-(-h4-h3)**5/120.0
  
  a(7,1)=-(h2+h1)**4/24.0;a(7,2)=-h2**4/24.0;a(7,3)=-h3**4/24.0;a(7,4)=-(-h4-h3)**4/24.0
  a(7,5)=(h2+h1)**6/720.0;a(7,6)=h2**6/720.0;a(7,7)=0;a(7,8)=h3**6/720.0;a(7,9)=(-h4-h3)**6/720.0
  
  a(8,1)=(h2+h1)**5/120.0;a(8,2)=h2**5/120.0;a(8,3)=-h3**5/120.0;a(8,4)=(-h4-h3)**5/120.0
  a(8,5)=-(h2+h1)**7/5040.0;a(8,6)=-h2**7/5040.0;a(8,7)=0;a(8,8)=h3**7/5040.0;a(8,9)=-(-h4-h3)**7/5040.0
  
  a(9,1)=-(h2+h1)**6/720.0;a(9,2)=-h2**6/720.0;a(9,3)=-h3**6/720.0;a(9,4)=-(-h4-h3)**6/720.0
  a(9,5)=(h2+h1)**8/40320.0;a(9,6)=h2**8/40320.0;a(9,7)=0;a(9,8)=h3**8/40320.0
  a(9,9)=(-h4-h3)**8/40320.0; 

  !-> define rhs
  b(1,1) = 0 ; b(2,1) = 0 ; b(3,1) = 1 ; b(4,1) = 0 ; b(5,1) = 0 ; b(6,1) = 0 ; b(7,1) = 0
  b(8,1) = 0 ; b(9,1) = 0
 
  fact='N' ; trans='N' ; nrhs=1 ; lda=n ; ldb=n ; ldaf=n ; equed='N' ; ldx=n
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,equed, r, c, b, ldb, x, &
       ldx, rcond, ferr, berr, work, iwork, info )

  !-> output coefficient
  alpha=x(1) ; beta=x(2) ; gamma=x(3) ; delta=x(4) 
  a1=x(5) ; b1=x(6) ; c1=x(7) ; d1=x(8) ; e1=x(9) 

  dum=0._rk
  write(333,'(11es17.8,i4,4es17.8)')alpha,beta,1._rk,gamma,delta,a1,b1,c1,d1,e1,dum,&
       info,rcond,ferr,berr,work(1)

end subroutine dder_coeffs_c

subroutine dder_coeffs_i2(a1,b1,c1,d1,e1,alpha,beta,gamma,delta,h1,h2,h3,h4)
! -----------------------------------------------------------------------
! Derivatives : Coefficients for second derivatives : non-uniform grid, second point
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
  implicit none
  real(rk),intent(out) :: a1,b1,c1,d1,e1,alpha,beta,gamma,delta
  real(rk),intent(in) :: h1,h2,h3,h4
  !
  ! n : size of linear system
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=8
  character(len=1) :: fact,trans,equed
  integer(ik) :: nrhs,lda,ldb,ldaf,ldx
  real(rk) :: a(n,n),af(n,n),ipiv(n),b(n,1),r(n),c(n),x(n),work(4*n)
  real(rk) :: rcond,ferr,berr
  integer(ik) :: iwork(n),info

  !-> space interval
  h=0._rk
  h(1)=h1 ; h(2)=h2 ; h(3)=h3 ; h(4)=h4 

  !-> define matrix

  A(1,1)=0 ; A(1,2)=0 ; A(1,3)=0 ; A(1,4)=1 ; A(1,5)=1 ; A(1,6)=1 ; A(1,7)=1 ; A(1,8)=1

  A(2,1)=0 ; A(2,2)=0 ; A(2,3)=0 ; A(2,4)=-h1 ; A(2,5)=0 ; A(2,6)=h2 ; A(2,7)=h3+h2
  A(2,8)=h4+h3+h2

  A(3,1)=-1 ; A(3,2)=-1 ; A(3,3)=-1 ; A(3,4)=h1**2/2.0 ; A(3,5)=0 ; A(3,6)=h2**2/2.0
  A(3,7)=(-h3-h2)**2/2.0 ; A(3,8)=(-h4-h3-h2)**2/2.0

  A(4,1)=-h2 ; A(4,2)=-h3-h2 ; A(4,3)=-h4-h3-h2 ; A(4,4)=-h1**3/6.0 ; A(4,5)=0
  A(4,6)=h2**3/6.0 ; A(4,7)=-(-h3-h2)**3/6.0 ; A(4,8)=-(-h4-h3-h2)**3/6.0

  A(5,1)=-h2**2/2.0 ; A(5,2)=-(-h3-h2)**2/2.0 ; A(5,3)=-(-h4-h3-h2)**2/2.0 ; A(5,4)=h1**4/24.0
  A(5,5)=0 ; A(5,6)=h2**4/24.0 ; A(5,7)=(-h3-h2)**4/24.0 ; A(5,8)=(-h4-h3-h2)**4/24.0

  A(6,1)=-h2**3/6.0 ; A(6,2)=(-h3-h2)**3/6.0 ; A(6,3)=(-h4-h3-h2)**3/6.0 ; A(6,4)=-h1**5/120.0
  A(6,5)=0 ; A(6,6)=h2**5/120.0 ; A(6,7)=-(-h3-h2)**5/120.0 ; A(6,8)=-(-h4-h3-h2)**5/120.0

  A(7,1)=-h2**4/24.0 ; A(7,2)=-(-h3-h2)**4/24.0 ; A(7,3)=-(-h4-h3-h2)**4/24.0
  A(7,4)=h1**6/720.0 ; A(7,5)=0 ; A(7,6)=h2**6/720.0 ; A(7,7)=(-h3-h2)**6/720.0
  A(7,8)=(-h4-h3-h2)**6/720.0

  A(8,1)=-h2**5/120.0 ; A(8,2)=(-h3-h2)**5/120.0 ; A(8,3)=(-h4-h3-h2)**5/120.0
  A(8,4)=-h1**7/5040.0 ; A(8,5)=0 ; A(8,6)=h2**7/5040.0 ; A(8,7)=-(-h3-h2)**7/5040.0
  A(8,8)=-(-h4-h3-h2)**7/5040.0

  !-> define rhs

  B(1,1)=0 ; B(2,1)=0 ; B(3,1)=1 ; B(4,1)=0 ; B(5,1)=0 ; B(6,1)=0 ; B(7,1)=0
  B(8,1)=0

  fact='N' ; trans='N' ; nrhs=1 ; lda=n ; ldb=n ; ldaf=n ; equed='N' ; ldx=n
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,equed, r, c, b, ldb, x, &
       ldx, rcond, ferr, berr, work, iwork, info )

  !-> output coefficient
  alpha=x(1) ; beta=1._rk ; gamma=x(2) 
  a1=x(3) ; b1=x(4) ; c1=x(5) ; d1=x(6) ; e1=x(7)

  alpha=1._rk ; beta=x(1) ; gamma=x(2) ; delta=x(3)
  a1=x(4) ; b1=x(5) ; c1=x(6) ; d1=x(7) ; e1=x(8)

  write(333,'(8es17.8,i4,4es17.8)')alpha,beta,gamma,a1,b1,c1,d1,e1,&
       info,rcond,ferr,berr,work(1)

end subroutine dder_coeffs_i2

subroutine dder_coeffs_i3(a1,b1,c1,d1,e1,alpha,beta,gamma,delta,h1,h2,h3,h4)
! -----------------------------------------------------------------------
! solver : Coefficients for second derivatives : non-uniform grid, first point
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
  implicit none
  real(rk),intent(out) :: a1,b1,c1,d1,e1,alpha,beta,gamma,delta
  real(rk),intent(in) :: h1,h2,h3,h4
  !
  ! n : size of linear system
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=6
  character(len=1) :: fact,trans,equed
  integer(ik) :: nrhs,lda,ldb,ldaf,ldx
  real(rk) :: a(n,n),af(n,n),ipiv(n),b(n,1),r(n),c(n),x(n),work(4*n)
  real(rk) :: rcond,ferr,berr
  integer(ik) :: iwork(n),info

  !-> space interval
  h(1)=h1 ; h(2)=h2 ; h(3)=h3 ; h(4)=h4 

  !-> define matrix

  A(1,1)=0 ; A(1,2)=1 ; A(1,3)=1 ; A(1,4)=1 ; A(1,5)=1 ; A(1,6)=1 

  A(2,1)=0 ; A(2,2)=-h2-h1 ; A(2,3)=-h2 ; A(2,4)=0 ; A(2,5)=h3 ; A(2,6)=h4+h3

  A(3,1)=-1 ; A(3,2)=(h2+h1)**2/2.0 ; A(3,3)=h2**2/2.0 ; A(3,4)=0 ; A(3,5)=h3**2/2.0
  A(3,6)=(-h4-h3)**2/2.0

  A(4,1)=h2 ; A(4,2)=-(h2+h1)**3/6.0 ; A(4,3)=-h2**3/6.0 ; A(4,4)=0 
  A(4,5)=h3**3/6.0 ; A(4,6)=-(-h4-h3)**3/6.0

  A(5,1)=-h2**2/2.0 ; A(5,2)=(h2+h1)**4/24.0 ; A(5,3)=h2**4/24.0 ; A(5,4)=0
  A(5,5)=h3**4/24.0 ; A(5,6)=(-h4-h3)**4/24.0

  A(6,1)=h2**3/6.0 ; A(6,2)=-(h2+h1)**5/120.0 ; A(6,3)=-h2**5/120.0
  A(6,4)=0 ; A(6,5)=h3**5/120.0 ; A(6,6)=-(-h4-h3)**5/120.0

  !-> define rhs

  B(1,1)=0 ; B(2,1)=0 ; B(3,1)=1 ; B(4,1)=0 ; B(5,1)=0 ; B(6,1)=0 !; B(7,1)=0
 
  fact='N' ; trans='N' ; nrhs=1 ; lda=n ; ldb=n ; ldaf=n ; equed='N' ; ldx=n
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,equed, r, c, b, ldb, x, &
       ldx, rcond, ferr, berr, work, iwork, info )

  !-> output coefficient
  alpha=x(1) ; beta=1._rk ; gamma=0._rk ; delta=0._rk
  a1=x(2) ; b1=x(3) ; c1=x(4) ; d1=x(5) ; e1=x(6) 

  write(333,'(8es17.8,i4,4es17.8)')alpha,beta,gamma,a1,b1,c1,d1,e1,&
       info,rcond,ferr,berr,work(1)
  
end subroutine dder_coeffs_i3

subroutine dder_coeffs_i3b(a1,b1,c1,d1,e1,alpha,beta,gamma,delta,h1,h2,h3,h4)
! -----------------------------------------------------------------------
! solver : Coefficients for second derivatives : non-uniform grid, first point
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
  implicit none
  real(rk),intent(out) :: a1,b1,c1,d1,e1,alpha,beta,gamma,delta
  real(rk),intent(in) :: h1,h2,h3,h4
  !
  ! n : size of linear system
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=5
  character(len=1) :: fact,trans,equed
  integer(ik) :: nrhs,lda,ldb,ldaf,ldx
  real(rk) :: a(n,n),af(n,n),ipiv(n),b(n,1),r(n),c(n),x(n),work(4*n)
  real(rk) :: rcond,ferr,berr
  integer(ik) :: iwork(n),info

  !-> space interval
  h(1)=h1 ; h(2)=h2 ; h(3)=h3 ; h(4)=h4 

  !-> define matrix

A(1,1) = 0
A(1,2) = 1
A(1,3) = 1
A(1,4) = 1
A(1,5) = 1
A(2,1) = 0
A(2,2) = -h2-h1
A(2,3) = -h2
A(2,4) = 0
A(2,5) = h3
A(3,1) = -1
A(3,2) = (h2+h1)**2/2.0
A(3,3) = h2**2/2.0
A(3,4) = 0
A(3,5) = h3**2/2.0
A(4,1) = h2
A(4,2) = -(h2+h1)**3/6.0
A(4,3) = -h2**3/6.0
A(4,4) = 0
A(4,5) = h3**3/6.0
A(5,1) = -h2**2/2.0
A(5,2) = (h2+h1)**4/24.0
A(5,3) = h2**4/24.0
A(5,4) = 0
A(5,5) = h3**4/24.0

  !-> define rhs

  B(1,1)=0 ; B(2,1)=0 ; B(3,1)=1 ; B(4,1)=0 ; B(5,1)=0 ! ; B(6,1)=0 !; B(7,1)=0
 
  fact='N' ; trans='N' ; nrhs=1 ; lda=n ; ldb=n ; ldaf=n ; equed='N' ; ldx=n
  call dgesvx( fact, trans, n, nrhs, a, lda, af, ldaf, ipiv,equed, r, c, b, ldb, x, &
       ldx, rcond, ferr, berr, work, iwork, info )

  !-> output coefficient
  alpha=x(1) ; beta=1._rk ; gamma=0._rk ; delta=0._rk
  a1=x(2) ; b1=x(3) ; c1=x(4) ; d1=x(5) ; e1=0._rk !x(6) 

  write(333,'(8es17.8,i4,4es17.8)')alpha,beta,gamma,a1,b1,c1,d1,e1,&
       info,rcond,ferr,berr,work(1)
  
end subroutine dder_coeffs_i3b

end module class_solver_coefficient_1d
