module class_solver_coefficient_3d
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
  real(rk),private,save :: h(10)
  
contains

subroutine der_coeffs_generic(coef_imp,coef_exp,h0,cl,der,bct,so)
! -----------------------------------------------------------------------
! Derivatives : Coefficients for first derivatives : non-uniform grid,centered 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2011
!
  implicit none
  real(rk),intent(out) :: coef_imp(5),coef_exp(8)
  real(rk),intent(in) :: h0(7)
  ! cl : position of the point
  ! der : order of derivative
  ! bct : type of boundary condition
  integer(ik),intent(in) :: cl,der,bct,so
  !
  ! n : size of linear system
  ! os : order of discretization
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=10
  integer(ik) :: info,n1,n2,n3,i,j,l,k
  real(rk),allocatable :: a(:,:), r(:)
  integer(ik),allocatable :: iwork(:)

h=0._rk

!->  ----- Stencil and space step -------

if (cl==1) then ! point on the edge

  if(so==8) then
    n1=3 ! implicite
    n2=6 ! explicite
!    n1=1 ! implicite
!    n2=6 ! explicite
  elseif(so==6) then
     n1=1 ! implicite
     n2=6 ! explicite
  elseif(so==4) then
     n1=1 ! implicite
     n2=4 ! explicite
  elseif(so==2) then
    n1=0 ! implicite
    n2=3 ! explicite
  endif
  n3=n1+n2

  h(1:n1)=h0(3:n1+2)
  h(n1+1:n3)=h0(1:n2)

elseif (cl==4) then ! extrapolation for point on the edge
  !Stencil
  if(so==8) then
    n1=0 ! implicite
    n2=7 ! explicite
  elseif(so==6) then
    n1=0 ! implicite
    n2=7 ! explicite
  elseif(so==4) then
    n1=0 ! implicite
    n2=5 ! explicite
  elseif(so==2) then
    n1=0 ! implicite
    n2=3 ! explicite
  endif
  n3=n1+n2

  h(1:n3)=h0(1:n3)

elseif (cl==2) then ! first point after the edge

  if(so==8) then
!Stencil
    if (bct==1) then
!       n1=3 ! implicite
!       n2=7 ! explicite
       n1=1 ! implicite
       n2=5 ! explicite
    elseif (bct==2)then
       n1=1 ! implicite
       n2=5 ! explicite
    elseif (bct==3) then
       n1=2 ! implicite
       n2=7 ! explicite
    endif
  elseif(so==6) then
    if (bct==1) then
        n1=1 ! implicite
        n2=5 ! explicite
    elseif (bct==2)then
        n1=1 ! implicite
        n2=5 ! explicite
    elseif (bct==3) then
        n1=1 ! implicite
        n2=3 ! explicite
    endif
  elseif(so==4) then
     if (bct==1) then
         n1=1 ! implicite
         n2=5 ! explicite
     elseif (bct==2)then
         n1=1 ! implicite
         n2=5 ! explicite
     elseif (bct==3) then
         n1=1 ! implicite
         n2=3 ! explicite
     endif 
  elseif(so==2) then
    n1=0 ! implicite
    n2=3 ! explicite
  endif
  n3=n1+n2

  h(1)=h0(2) ; h(2:n1)=h0(4:n1+2)
  i=max(0,(5-n2)/2) ! shift
  h(n1+1:n3)=h0(1+i:n2+i)

elseif (cl==3) then ! every other points
  !Stencil
  if(so==8) then
     n1=4 ! implicite
     n2=5 ! explicite
  elseif(so==6) then
     n1=2 ! implicite
     n2=5 ! explicite
  elseif(so==4) then
     n1=2 ! implicite
     n2=3 ! explicite
  elseif(so==2) then
    n1=0 ! implicite
    n2=3 ! explicite
  endif
  n3=n1+n2

  i=(4-n1)/2 ! shift
  h(1:2)=h0(1+i:2+i) ; h(3-i:n1)=h0(4:n1+1+i)
  i=(5-n2)/2 ! shift
  h(n1+1:n3)=h0(1+i:n2+i)

endif

!->  ----- matrix -------

allocate(a(n3,n3))
allocate(r(n3))
allocate(iwork(n3))

a=0._rk

!column for implicit points
a(der+1,1:n1)=-1._rk 
do i=der+2,n3
do j=1,n1
  a(i,j)=a(i-1,j)*h(j)/(i-der-1)
enddo
enddo

!column for explicit points
a(1,n1+1:n3)=1._rk 
do i=2,n3
do j=n1+1,n3
  a(i,j)=a(i-1,j)*h(j)/(i-1)
enddo
enddo


!boundary condition
if((cl==2.and.n2>4).or.cl==1) then
  a(1:bct-1,n1+1)=0._rk
  a(bct,n1+1)=1._rk
  if(bct==3) a(bct,n1+1)=-1._rk
  do i=bct+1,n3
    a(i,n1+1)=a(i-1,n1+1)*h(n1+1)/(i-bct)
  enddo
endif


r=0._rk
r(1+der)=1._rk ! der

!->  ----- solve -------

iwork=0
call dgesv( n3, 1, a, n3, iwork, r, n3, info )
if(info/=0) print*,'ERROR in discretization'

coef_imp=0._rk
coef_exp=0._rk

!-> ----- export ----

j=0
l=n1
k=0
if(cl==3) k=(4-n1)/2  ! shift
if(n1+k.ge.cl) l=n1+1 ! quasi-always true
do i=1,l
     if(i/=cl-k) then
    j=j+1
    coef_imp(i+k)=r(j) 
  endif
enddo
coef_imp(cl)=1._rk

i=0
if(cl==3) i=(5-n2)/2 ! shift
  if(cl==2) i=max(0,(5-n2)/2) ! shift
coef_exp(1+i:n2+i)=r(n1+1:n3)


!-> ---- end ----

deallocate(a)
deallocate(r)
deallocate(iwork)

end subroutine der_coeffs_generic

end module class_solver_coefficient_3d
