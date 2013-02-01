module class_filter
! -----------------------------------------------------------------------
! Name :
! class filter 
! -----------------------------------------------------------------------
  use precision
  implicit none

  type,public :: filter
     !-> input dimensions
     integer(ik) :: nx_in,ny_in,nz_in
     !-> original grid
     real(rk),allocatable :: grid_in(:,:,:)
     !-> filtered dimensions
     integer(ik) :: nx_out,ny_out,nz_out
     !-> filtered grid
     real(rk),allocatable :: grid_out(:,:,:)
     !-> filter width dimensions
     integer(ik) :: nx_nodes,ny_nodes,nz_nodes
     !-> filter width grid
     real(rk),allocatable :: grid_nodes(:,:,:)
     !-> filter order
     integer(ik) :: order
  end type filter

  !-> Declare everything private by default
  private
  
  !-> Declare exported procedure
  public :: flt_init,LSQUnivariateSpline
  public :: print_filter
contains

  subroutine flt_init(filt_coeffs,na1,na2,na3,nb1,nb2,nb3,nc1,nc2,nc3,&
       grid1,grid2,grid3,order)
    implicit none
    type(filter),intent(out) :: filt_coeffs
    integer(ik),intent(in) :: na1,na2,na3,nb1,nb2,nb3,nc1,nc2,nc3
    real(rk),intent(in) :: grid1(:,:,:),grid2(:,:,:),grid3(:,:,:)
    integer(ik) :: order

    !-> input dimensions
    filt_coeffs%nx_in=na1
    filt_coeffs%ny_in=na2
    filt_coeffs%nz_in=na3

    !-> filtered dimensions
    filt_coeffs%nx_out=nb1
    filt_coeffs%ny_out=nb2
    filt_coeffs%nz_out=nb3

    !-> filter width dimensions
    filt_coeffs%nx_nodes=nc1
    filt_coeffs%ny_nodes=nc2
    filt_coeffs%nz_nodes=nc3

    !-> allocate grids
    call grid3d_allocate(filt_coeffs%nx_in,filt_coeffs%ny_in,filt_coeffs%nz_in,&
         filt_coeffs%grid_in)
    call grid3d_allocate(filt_coeffs%nx_out,filt_coeffs%ny_out,filt_coeffs%nz_out,&
         filt_coeffs%grid_out)
    call grid3d_allocate(filt_coeffs%nx_nodes,filt_coeffs%ny_nodes,filt_coeffs%nz_nodes,&
         filt_coeffs%grid_nodes)

    !-> grids
    filt_coeffs%grid_in=grid1
    filt_coeffs%grid_out=grid2
    filt_coeffs%grid_nodes=grid3
    
    !-> order
    filt_coeffs%order=order

  end subroutine flt_init

  subroutine print_filter(filt)
    implicit none
    type(filter) :: filt

    print*,filt%grid_in,filt%grid_out,filt%grid_nodes

  end subroutine print_filter


  subroutine grid3d_allocate(n1,n2,n3,grid)
! -----------------------------------------------------------------------
! filter : allocate or reallocate 3d grid in type mesh
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2011
!
    implicit none
    integer,intent(in) :: n1,n2,n3
    real(rk),allocatable :: grid(:,:,:)

    if (.not.allocated(grid)) then
       allocate(grid(n1,n2,n3))
    elseif(allocated(grid)) then
       deallocate(grid)
       allocate(grid(n1,n2,n3))
    endif

  end subroutine grid3d_allocate

  subroutine LSQUnivariateSpline(ni,xi,yi,wi,nn,xn,no,xo,yo,ypo,k)

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

    integer(ik) :: no,ni,nn,k,ier
    real(rk) :: xi(ni),yi(ni),wi(ni)   
    real(rk) :: xo(no),yo(no),ypo(no)              
    real(rk) :: xn(nn)    

    integer(ik) :: nest,m
    integer(ik) :: iwrk(ni+k+1),lwrk
    real(rk) :: xb,xe,fp,wrk(ni*(k+1)+(ni+k+1)*(7+3*k)),wrk2(no)
    real(rk) :: t(ni+k+1),c(ni+k+1)

    nest = ni+k+1
    lwrk = ni*(k+1)+nest*(7+3*k)
    xb = xi(1)        
    xe = xi(ni)        
  
    m = nn+2*k

    t(k+2:m-k-1) = xn(2:nn-1)
    
    t(1:k+1) = xb
    t(m-k:) = xe

    !-> compute filter coefficients
    call curfit(-1,ni,xi,yi,wi,xb,xe,k,0.,nest,m,t,c,fp,wrk,lwrk,iwrk,ier)


!    print *,'fp=',fp
    
    if (ier .ne. 0) then
       print *,'FILTER : Error with curfit:',ier
       stop
    endif

    !-> give filtered function yo on coordinates xo
    call splev(t,m,c,k,xo,yo,no,ier)
    
    !-> give derivatives of filtered function ypo on coordinates xo
    call splder(t,m,c,k,1,xo,ypo,no,wrk2,ier)


end subroutine LSQUnivariateSpline


end module class_filter
