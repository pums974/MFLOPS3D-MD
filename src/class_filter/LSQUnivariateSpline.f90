

  subroutine LSQUnivariateSpline(xi,yi,xn,xo,yo,ypo,wi,k,ni,nn,no,ier)

!---------------------------------------------------------------------
!
! INPUT :
!
!   xi : coordinate of the function to filter
!   yi : function to filter
!   xn : coordinate of the nodes for Spline interpolation
!   wi : weight factor for each points of the function to filter
!   k : order of polynomial
!   ni : number of points of the function to filter
!   nn : number of nodes for Spline interpolation
!   no : number of points of the  filter function
!   xo : coordinate of the filtered function
!
! OUTPUT :
!
!   yo : filtered function
!   ier : error code (0 for normal)
!
! REQUIREMENT : 
!
! xi[1]  = xn[1]  = xb 
! xi[ni] = xn[nn] = xe
!
!---------------------------------------------------------------------

  implicit none

  integer(4) :: no,ni,nn,k,ier
  real(8) :: xi(ni),yi(ni),wi(ni)   
  real(8) :: xo(no),yo(no),ypo(no)              
  real(8) :: xn(nn)    

  integer(4) :: nest,m
  integer(4) :: iwrk(ni+k+1),lwrk
  real(8) :: xb,xe,fp,wrk(ni*(k+1)+(ni+k+1)*(7+3*k)),wrk2(no)
  real(8) :: t(ni+k+1),c(ni+k+1)

  nest = ni+k+1
  lwrk = ni*(k+1)+nest*(7+3*k)
  xb = xi(1)        
  xe = xi(ni)        
  
  m = nn+2*k

  t(k+2:m-k-1) = xn(2:nn-1)

  t(1:k+1) = xb
  t(m-k:) = xe

  call curfit(-1,ni,xi,yi,wi,xb,xe,k,0.,nest,m,t,c,fp,wrk,lwrk,iwrk,ier)

  print *,'fp=',fp

  if (ier .ne. 0) then
    print *,'Error with curfit:',ier
    stop
  endif

  call splev(t,m,c,k,xo,yo,no,ier)

  call splder(t,m,c,k,1,xo,ypo,no,wrk2,ier)


end subroutine LSQUnivariateSpline
