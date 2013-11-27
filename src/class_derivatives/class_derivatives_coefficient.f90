module class_derivatives_coefficient
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
  use precision
  implicit none
!  private
  real(rk),private,save :: h(10),h0(10)
  
contains

subroutine der_coeffs_generic(coef_imp,coef_exp,h0,cl,der,so)
! -----------------------------------------------------------------------
! Derivatives : Coefficients for first derivatives : non-uniform grid,centered 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2011
!
  implicit none
  real(rk),intent(out) :: coef_imp(5),coef_exp(7)
  real(rk),intent(in) :: h0(7)
  ! cl : position of the point
  ! der : order of derivative
  integer(ik),intent(in) :: cl,der,so
  !
  ! n : size of linear system
  ! so : order of discretization
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=10
  integer(ik) :: info,n1,n2,nt,i,j,l,k
  real(rk),allocatable :: a(:,:), r(:)
  integer(ik),allocatable :: iwork(:)

h=0._rk

  !->  ----- Stencil and space step -------
  !->  ----- n1+n2-n3 MUST be ge 3 (dder) -------

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
  nt=n1+n2

  h(1:n1)=h0(2:n1+1)
  h(n1+1:nt)=h0(1:n2)

elseif (cl==2) then ! first point after the edge

  if(so==8) then
!     n1=3 ! implicite
!     n2=7 ! explicite
     n1=1 ! implicite
     n2=5 ! explicite
  elseif(so==6) then
     n1=2 ! implicite
     n2=6 ! explicite
!     n1=1 ! implicite
!     n2=5 ! explicite
  elseif(so==4) then
     n1=1 ! implicite
     n2=5 ! explicite
  elseif(so==2) then
    n1=0 ! implicite
    n2=3 ! explicite
  endif
  nt=n1+n2

  h(1)=h0(1) ; h(2:n1)=h0(3:n1+1)
  h(n1+1:nt)=h0(1:n2)

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
  nt=n1+n2

  i=(4-n1)/2 ! shift
  h(1:2)=h0(1+i:2+i) ; h(3-i:n1)=h0(4:n1+1+i)
  i=(5-n2)/2 ! shift
  h(n1+1:nt)=h0(1+i:n2+i)

endif

!->  ----- matrix -------

  allocate(a(nt,nt))
  allocate(r(nt))
  allocate(iwork(nt))

a=0._rk

!column for implicit points
a(der+1,1:n1)=-1._rk 
  do i=der+2,nt
do j=1,n1
  a(i,j)=a(i-1,j)*h(j)/(i-der-1)
enddo
enddo

!column for explicit points
  a(1,n1+1:nt)=1._rk
  do i=2,nt
     do j=n1+1,nt
  a(i,j)=a(i-1,j)*h(j)/(i-1)
enddo
enddo

r=0._rk
r(1+der)=1._rk ! der

!->  ----- solve -------

iwork=0
  call dgesv( nt, 1, a, nt, iwork, r, nt, info )
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
  coef_exp(1+i:n2+i)=r(n1+1:nt)


!-> ---- end ----

deallocate(a)
deallocate(r)
deallocate(iwork)

end subroutine der_coeffs_generic


subroutine der_coeffs_mac_1(coefu,coefp,hu,hp,cl,der)
! -----------------------------------------------------------------------
! Derivatives : Coefficients for first derivatives : non-uniform grid,centered 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2011
!
  implicit none
  real(rk),intent(out) :: coefu(5),coefp(7)
  real(rk),intent(in) :: hu(5),hp(7)
  integer(ik),intent(in) :: cl,der
  !
  ! n : size of linear system
  ! info : output of solver (zero is ok) 
  integer,parameter :: n=10
  real(rk) :: h0(5),x(n)
  integer(ik) :: info,os,n1,n2,n3,i,j,l
  real(rk),allocatable :: a(:,:), r(:)
  integer(ik),allocatable :: iwork(:)
os=2
h=0._rk

if (cl==1) then

!Stencil
if(os==8) then
n1=2 ! implicite
n2=6 ! explicite
elseif(os==2) then
n1=0 ! implicite
n2=3 ! explicite
endif
n3=n1+n2

h(1:n1)=hu(2:n1+1)
h(n1+1:n3)=hp(1:n2)

elseif (cl==2) then

!Stencil
if(os==8) then
n1=3 ! implicite
n2=7 ! explicite
elseif(os==2) then
n1=0 ! implicite
n2=3 ! explicite
endif
n3=n1+n2

h(1)=hu(1) ; h(2:n1)=hu(3:n1+1)
h(n1+1:n3)=hp(1:n2)

elseif (cl==3) then
!Stencil
if(os==8) then
n1=4 ! implicite
n2=4 ! explicite
elseif(os==2) then
n1=0 ! implicite
n2=2 ! explicite
endif
n3=n1+n2

i=(5-n2)/2
h(1:2)=hu(1:2) ; h(3:n1)=hu(4:n1+1)
h(n1+1:n3)=hp(1+i:n2+i)

endif

allocate(a(n3,n3))
allocate(r(n3))
allocate(iwork(n3))

a=0._rk

a(2,1:n1)=-1._rk 
do i=3,n3
do j=1,n1
  a(i,j)=a(i-1,j)*h(j)/(i-2)
enddo
enddo

a(1,n1+1:n3)=1._rk 
do i=2,n3
do j=n1+1,n3
  a(i,j)=a(i-1,j)*h(j)/(i-1)
enddo
enddo

r=0._rk
r(1+der)=1._rk ! der

!print*,cl
!do i=1,n3
!write(*,'(I4,30es17.8)') i,a(i,:),r(i)
!enddo

iwork=0
  call dgesv( n3, 1, a, n3, iwork, r, n3, info )
if(info/=0) print*,'ERROR'

!do i=1,n3
!write(*,'(I4,30es17.8)') i,r(i)
!enddo
!if(cl==1) stop

x=0._rk
x(1:n3)=r(1:n3)
coefu=0._rk
coefp=0._rk

  !-> output coefficient
j=0
l=n1
if(n1.ge.cl) l=n1+1
do i=1,l
  if(i/=cl) then
    j=j+1
    coefu(i)=x(j) 
  endif
enddo
coefu(cl)=1._rk

if(cl==3) then
  i=(4-n2)/2
  coefp(1+i:n2+i)=x(n1+1:n3)
else
  coefp(1:n2)=x(n1+1:n3)
endif


deallocate(a)
deallocate(r)
deallocate(iwork)

end subroutine der_coeffs_mac_1

end module class_derivatives_coefficient
