! =====================================================================
! =                  TEST OF LEAST SQUARES SPLINE FILTER              = 
! =====================================================================
!goal    : to filter (ni-1) data in 1D (because the original function  *
!          is periodic,and a period is taken to research here, so the  *
!          first and the last points are coincident completely,so one  *
!          point of these two is enough)                               *
!   k    : the order of spline, which can be 1,2,3,4,5, but k=3 is the *
!          best by experience                                          *
!yi(ni+1): the size of the array yi is required to be the points number*
!          plus 2 in order to do FFT,but here the points number used   *
!          to do FFT is (ni-1), so the size of fi should be ni+1.      *
!          this rule applies to yo(no+1) too.                          *
!fftfax  : to initialize the values in order to simplify the fft in the* 
!          following computation.so the initializing points number is  *
!          (ni-1) and (no-1), in order to coincide with above-mentioned* 
!          research goal                                               *
!xi(ni)   : from -1 to 1, equidistantly,the same to xo(no)             *
!xn(nn+1) : from -1 to 1, equidistantly, in order to keep nn intervals,* 
!          the size of xn should be nn+1                               *
!    j    : the scope of j is from 1 to ni/2, because the largest      *
!           wavenumber after FFT is ni/2, and it's better to be ni/2   *
!           in order to see the overall view of all the effective k    *
!=======================================================================
program lsspline  
!  use spline
  implicit none

                                                       

  integer(4), parameter :: ni=4097,no=4097,nn=(ni-1)/32,k=5      ! 2 <= nn <= ni-k (if not periodic)
                                                               ! 2 <= nn <= ni-1 (if periodic)


  real(8) :: xi(ni),yi(ni+1),wi(ni+1),spi(0:(ni-1)/2)
  real(8) :: xo(no),yo(no+1),ypo(no+1),yppo(no+1),spo(0:(no-1)/2)
  real(8) :: xn(nn+1)

  real(8)    :: wki(4*(ni+1)),wko(4*(no+1))
  complex(8) :: exi(3*ni/2),exo(3*no/2)
  integer(4) :: ifaxi(100),ifaxo(100)

  real(8)    :: pi,ang
  integer(4) :: i,j,l,ier,npt

  real(8)    ::  c1(4,ni),ivec(no)

  pi = acos(-1.)

  spo(:) = 0.
  spi(:) = 0.

!  call randini

!  call fftfax(ni-1,ifaxi,exi)
!  call fftfax(no-1,ifaxo,exo)

  do i=1,ni
   xi(i) = 1.*float(2*i-ni-1)/float(ni-1)
  enddo

  do i=1,no
   xo(i) = 1.*float(2*i-no-1)/float(no-1)
  enddo

  do i=1,nn+1
   xn(i) = 1.*float(2*(i-1)-nn)/float(nn)
  enddo


 do l=1,100 

 yi(:) = 0.
  do j=1,ni/2
!    call randa(ang)
    do i=1,ni
      yi(i) = yi(i)+ sqrt(2.0)*cos(pi*xi(i)*float(j)+2.*pi*ang)  ! *(abs(xi(i))-1.)**8
    enddo
  enddo

  wi(:)  = 1.

!  npt = (ni-1)/nn
!  do i=1,nn+1
!  do j=0,npt/2
!   wi((i-1)*npt+j+1) = (j+1)**2
!   wi((i  )*npt-j+1) = (j+1)**2
!  enddo
!  enddo


!  call cubic_spline(yi,xi,yo,ypo,yppo,xo,c1,ivec,ni,no)


  call LSQUnivariateSpline(xi,yi,xn,xo,yo,ypo,wi,k,ni,nn+1,no,ier)


  OPEN(10,file='functions.dat')   
    do i=1,no
       WRITE(10,'(i6,3(3x,e12.4))')i,yi(i),yo(i), yi(i)-yo(i) !  ypo(i)/10.
    end do
  CLOSE(10)


 enddo

  end
