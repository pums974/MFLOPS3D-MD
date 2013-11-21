module class_mapping
!------------------------------------------------------------------------
! Name :
! class mapping
! -----------------------------------------------------------------------
! Object :
! Mapping of Navier-stokes equation
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
  use precision
  use class_derivatives
  use class_field
  use class_mesh
  implicit none
  
  type mapping
     !-> dimension
     integer(ik) :: n
     !-> mapping : 1 -> yes, 0-> no
     integer(ik) :: mapt
     !-> dimensions name
     character(len=512) :: nxn="resolution"
     !-> mapping functions : 1,: -> function, :,2 -> lower and upper
     real(rk),allocatable :: etha(:,:)
     !-> derivatives of mapping function
     real(rk),allocatable :: detha(:,:),ddetha(:,:),dtetha(:,:)
     !-> mapping coefficients
     type(field) :: ydx,ydy,yddx,ydx2,ydy2,ydt
     !-> tangent and normal vectors
     type(boundary_condition) :: tan(3),nor(3)
     !-> derivatives coefficients
     type(derivatives_coefficients),pointer :: dcx,dcy,dcz
  end type mapping

contains

  function dertm_nc(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dertm_nc
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dertm_nc,"DXU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dertm_nc=0._rk
    elseif (dc%mapt==1) then
       dertm_nc=dc%ydt*dery(dc%dcy,x)
    endif

  end function dertm_nc

  function derxm(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: derxm
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(derxm,"DXU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       derxm=derx(dc%dcx,x)
    elseif (dc%mapt==1) then
       derxm=derx(dc%dcx,x)+dc%ydx*dery(dc%dcy,x)
    endif

  end function derxm

  function dderxm(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderxm
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderxm,"DDXU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dderxm=dderx(dc%dcx,x)
    elseif (dc%mapt==1) then
       dderxm=dderx(dc%dcx,x)+dc%ydx*dc%ydx*ddery(dc%dcy,x) &
            +2._rk*dc%ydx*derx(dc%dcx,dery(dc%dcy,x)) &
            +dc%yddx*dery(dc%dcy,x)
    endif

  end function dderxm

  function dderxm_nc(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderxm_nc
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderxm_nc,"DDXU",x%nx,x%ny,x%nz)

    if (dc%mapt==0) then
       dderxm_nc=0._rk
    elseif (dc%mapt==1) then
       dderxm_nc=dc%ydx2*ddery(dc%dcy,x) &
            +2._rk*dc%ydx*derx(dc%dcx,dery(dc%dcy,x)) &
            +dc%yddx*dery(dc%dcy,x)
    endif

  end function dderxm_nc

  function derym(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: derym
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(derym,"DYU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       derym=dery(dc%dcy,x)
    elseif (dc%mapt==1) then
       derym=dc%ydy*dery(dc%dcy,x)
    endif

  end function derym

  function dderym(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderym
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderym,"DDYU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dderym=ddery(dc%dcy,x)
    elseif (dc%mapt==1) then
       dderym=dc%ydy*dc%ydy*ddery(dc%dcy,x)
    endif

  end function dderym

  function dderym_nc(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderym_nc
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderym_nc,"DDYU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dderym_nc=0._rk
    elseif (dc%mapt==1) then
       dderym_nc=dc%ydy2*ddery(dc%dcy,x)
    endif

  end function dderym_nc

  function derzm(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: derzm
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(derzm,"DZU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       derzm=derz(dc%dcz,x)
    elseif (dc%mapt==1) then
       derzm=derz(dc%dcz,x)
    endif

  end function derzm

  function dderzm(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderzm
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderzm,"DZU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dderzm=dderz(dc%dcz,x)
    elseif (dc%mapt==1) then
       dderzm=dderz(dc%dcz,x)
    endif

  end function dderzm

  function dderzm_nc(dc,x)
! -----------------------------------------------------------------------
! mapping : compute first derivative in x direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 06/2013
!
    use class_derivatives
    implicit none
    type(field) :: dderzm_nc
    type(mapping),intent(in) :: dc
    type(field),intent(in) :: x
    call field_init(dderzm_nc,"DZU",x%nx,x%ny,x%nz)
  
    if (dc%mapt==0) then
       dderzm_nc=0._rk
    elseif (dc%mapt==1) then
       dderzm_nc=0._rk
    endif

  end function dderzm_nc

  subroutine mapping_bcphi(map,aux,x,bc)
! -----------------------------------------------------------------------
! mapping : mapping allocate or reallocate 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    use class_derivatives
    use class_md
    implicit none
    type(field),intent(in) :: x
    type(field),intent(inout) :: aux
    type(mapping),intent(in) :: map
    type(boundary_condition),intent(inout) :: bc

    !-> x direction
    aux=dery(map%dcy,x)
    bc%bcx(:,:,1)=-map%ydx%f(1,2:x%ny-1,2:x%nz-1)*&
         aux%f(1,2:x%ny-1,2:x%nz-1)
    bc%bcx(:,:,2)=-map%ydx%f(x%nx,2:x%ny-1,2:x%nz-1)*&
         aux%f(x%nx,2:x%ny-1,2:x%nz-1)
    
    !-> y direction
    aux=derx(map%dcx,x)
    bc%bcy(:,:,1)=map%ydx%f(2:x%nx-1,1,2:x%nz-1)* &
         aux%f(2:x%nx-1,1,2:x%nz-1)/(map%ydy%f(2:x%nx-1,1,2:x%nz-1)-&
         map%ydx%f(2:x%nx-1,1,2:x%nz-1)**2)
    
    bc%bcy(:,:,2)=map%ydx%f(2:x%nx-1,x%ny,2:x%nz-1)* &
         aux%f(2:x%nx-1,x%ny,2:x%nz-1)/(map%ydy%f(2:x%nx-1,x%ny,2:x%nz-1)-&
         map%ydx%f(2:x%nx-1,x%ny,2:x%nz-1)**2)
    
    !-> z direction
    bc%bcz=0._rk

  end subroutine mapping_bcphi

! =======================================================================
! =======================================================================
! mapping : initialization and destruction methods
! =======================================================================
! =======================================================================

  subroutine mapping_init(mpid,gridx,gridy,dcx,dcy,dcz,aux,map,time)
! -----------------------------------------------------------------------
! mapping : mapping allocate or reallocate 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    use class_derivatives
    use class_md
    implicit none
    type(mpi_data) :: mpid
    type(field) :: aux
    type(derivatives_coefficients),intent(in),target :: dcx,dcy,dcz
    type(mapping),intent(inout) :: map
    type(mesh_grid),intent(inout) :: gridx,gridy
    integer(ik) :: i,j,k,c(3)
    real(rk) :: amp,sig,xo,x,vnorm,ya,yb,a,ha,hb,pi,time
    
    pi=4._rk*atan(1._rk)

    !--------------------------------------------------------------------
    !-> allocate type mapping
    call mapping_allocate(gridx%n,map)

    !-> allocate type derivatives 
    allocate(map%dcx,map%dcy,map%dcz)
    map%dcx=>dcx ; map%dcy=>dcy ; map%dcz=>dcz 

    !--------------------------------------------------------------------
    c(1)=mpid%coord(1) ; c(2)=mpid%coord(2) ; c(3)=mpid%coord(3)
    ya=gridy%grid1d(1)
    yb=gridy%grid1d(gridy%n)

    !--------------------------------------------------------------------
    !-> define lower mapping
    if (c(2)==0) amp=0.6_rk 
    if (c(2)==1) amp=0._rk 
    sig=1.5_rk ; xo=5._rk 
    do i=1,map%n
       x=gridx%grid1d(i)
       map%etha(i,1)=gridy%grid1d(1)+&
            amp*exp(-(x-xo)**2/(2._rk*sig**2)) !*sin(pi*time)
       map%dtetha(i,1)=amp*exp(-(x-xo)**2/(2._rk*sig**2))*pi*cos(pi*time)
       map%dtetha(i,1)=0._rk

       goto 101
       if (c(1)==0) map%etha(i,1)=0._rk
       ha=0._rk
       hb=0.2_rk
       a=(hb-ha)/&
            (gridx%grid1d(gridx%n)-gridx%grid1d(1))
       if (c(1)==1) map%etha(i,1)=a*x+hb-&
            gridx%grid1d(gridx%n)*a
       ha=0.2_rk
       hb=0._rk
       a=(hb-ha)/&
            (gridx%grid1d(gridx%n)-gridx%grid1d(1))
       if (c(1)==2) map%etha(i,1)=a*x+hb-&
            gridx%grid1d(gridx%n)*a
       if (c(1)>2) map%etha(i,1)=hb
101    continue

    enddo

    !--------------------------------------------------------------------
    !-> define upper mapping
    if (c(2)==0) amp=0._rk
    if (c(2)==1) amp=0._rk 
    sig=1.5_rk ; xo=5._rk 
    do i=1,map%n
       x=gridx%grid1d(i)
       map%etha(i,2)=gridy%grid1d(gridy%n)+&
            amp*exp(-(x-xo)**2/(2._rk*sig**2))
       map%dtetha(i,2)=0._rk
    enddo
    
    if (map%mapt==0) then
       map%etha(:,1)=gridy%grid1d(1)
       map%dtetha(:,1)=0._rk
       map%etha(:,2)=gridy%grid1d(gridy%n)
       map%dtetha(:,2)=0._rk
    endif

    !--------------------------------------------------------------------
    !-> compute lower mapping first derivative
    do k=1,aux%nz
       do j=1,aux%ny
          aux%f(:,j,k)=map%etha(:,1)
       enddo
    enddo
    aux=derx(map%dcx,aux)
    map%detha(:,1)=aux%f(:,2,2)

    !-> compute upper mapping first derivative
    do k=1,aux%nz
       do j=1,aux%ny
          aux%f(:,j,k)=map%etha(:,2)
       enddo
    enddo
    aux=derx(map%dcx,aux)
    map%detha(:,2)=aux%f(:,2,2)

    !-> compute lower mapping second derivative
    do k=1,aux%nz
       do j=1,aux%ny
          aux%f(:,j,k)=map%detha(:,1)
       enddo
    enddo
    aux=derx(map%dcx,aux)
    map%ddetha(:,1)=aux%f(:,2,2)
    
    !-> compute upper mapping second derivative
    do k=1,aux%nz
       do j=1,aux%ny
          aux%f(:,j,k)=map%detha(:,2)
       enddo
    enddo
    aux=derx(map%dcx,aux)
    map%ddetha(:,2)=aux%f(:,2,2)
    
    !--------------------------------------------------------------------
    !-> allocate mapping coefficients
    call field_init(map%ydx,"ydx",aux%nx,aux%ny,aux%nz)
    call field_init(map%ydy,"ydy",aux%nx,aux%ny,aux%nz)
    call field_init(map%yddx,"yddx",aux%nx,aux%ny,aux%nz)
    call field_init(map%ydx2,"ydx2",aux%nx,aux%ny,aux%nz)
    call field_init(map%ydy2,"ydy2",aux%nx,aux%ny,aux%nz)
    call field_init(map%ydt,"ydt",aux%nx,aux%ny,aux%nz)

    !--------------------------------------------------------------------
    !-> compute mapping coefficients
    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%ydx%f(i,j,k)=-(yb-ya)/(map%etha(i,2)-map%etha(i,1))**2* &
             (map%detha(i,2)*(gridy%grid1d(j)-map%etha(i,1))- &
             map%detha(i,1)*(gridy%grid1d(j)-map%etha(i,2)))
          enddo
       enddo
    enddo
    
    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%ydy%f(i,j,k)=(yb-ya)/(map%etha(i,2)-map%etha(i,1))
          enddo
       enddo
    enddo
    
    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%yddx%f(i,j,k)=&
                  (map%ddetha(i,2)*ya-map%ddetha(i,1)*yb)/ &
                  (map%etha(i,2)-map%etha(i,1))- &
                  (map%ddetha(i,2)-map%ddetha(i,1))* &
                  (gridy%grid1d(j)*(yb-ya)-map%etha(i,1)*yb+map%etha(i,2)*ya)/ &
                  (map%etha(i,2)-map%etha(i,1))**2- &
                  2._rk*(map%detha(i,2)-map%detha(i,1))* &
                  (map%detha(i,2)*ya-map%detha(i,1)*yb)/ &
                  (map%etha(i,2)-map%etha(i,1))**2+ &
                  2._rk*(gridy%grid1d(j)*(yb-ya)- &
                  map%etha(i,1)*yb+map%etha(i,2)*ya)* &
                  (map%detha(i,2)-map%detha(i,1))**2/ &
                  (map%etha(i,2)-map%etha(i,1))**3
          enddo
       enddo
    enddo
    
   do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%ydx2%f(i,j,k)=map%ydx%f(i,j,k)*map%ydx%f(i,j,k)
!             map%ydx2%f(i,j,k)=((yb-ya)**2*&
!                  (map%detha(i,2)*(gridy%grid1d(j)-map%etha(i,1))-&
!                  map%detha(i,1)*(gridy%grid1d(j)-map%etha(i,2)))**2-&
!                  (map%etha(i,2)-map%etha(i,1))**4)/&
!                  (map%etha(i,2)-map%etha(i,1))**4
          enddo
       enddo
    enddo
 
    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%ydy2%f(i,j,k)=((yb-ya)**2-&
                  (map%etha(i,2)-map%etha(i,1))**2)/&
                  (map%etha(i,2)-map%etha(i,1))**2
          enddo
       enddo
    enddo

    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             map%ydt%f(i,j,k)=-(yb-ya)/(map%etha(i,2)-map%etha(i,1))**2* &
             (map%dtetha(i,2)*(gridy%grid1d(j)-map%etha(i,1))- &
             map%dtetha(i,1)*(gridy%grid1d(j)-map%etha(i,2)))
          enddo
       enddo
    enddo

    !-> change gridy
    do k=1,aux%nz
       do j=1,aux%ny
          do i=1,aux%nx
             gridy%grid3d(i,j,k)=((map%etha(i,2)-map%etha(i,1))*&
                  gridy%grid1d(j)+map%etha(i,1)*yb-map%etha(i,2)*ya)/&
                  (yb-ya)
             enddo
          enddo
       enddo
    
    goto 100
    !--------------------------------------------------------------------
    !-> allocate vectors WE WILL SEE IF NEEDED
    do i=1,3
       call boundary_condition_init(map%tan(3),aux%nx,aux%ny,aux%nz)
       call boundary_condition_init(map%nor(3),aux%nx,aux%ny,aux%nz)
    enddo

    !-> compute tangent and normal vectors 
    do k=1,aux%nz-2
       do j=1,aux%ny-2
          map%tan(1)%bcx(j,k,1)=1._rk
          map%tan(2)%bcx(j,k,1)=0._rk
          map%tan(3)%bcx(j,k,1)=0._rk
          map%nor(1)%bcx(j,k,1)=0._rk
          map%nor(2)%bcx(j,k,1)=1._rk
          map%nor(3)%bcx(j,k,1)=0._rk
       enddo
    enddo

    do k=1,aux%nz-2
       do i=1,aux%nx-2
          vnorm=sqrt(1._rk/(1._rk+map%ydx%f(i+1,1,k+1)*map%ydx%f(i+1,1,k+1)))
          map%tan(1)%bcy(i,k,1)=1._rk/vnorm
          map%tan(2)%bcy(i,k,1)=map%ydx%f(i+1,1,k+1)/vnorm
          map%tan(3)%bcy(i,k,1)=0._rk

          vnorm=sqrt(1._rk/(1._rk+map%ydx%f(i+1,1,k+1)*map%ydx%f(i+1,1,k+1)))
          map%nor(1)%bcy(i,k,1)=map%ydx%f(i+1,1,k+1)/vnorm
          map%nor(2)%bcy(i,k,1)=1._rk/vnorm
          map%nor(3)%bcy(i,k,1)=0._rk             
       enddo
    enddo
100 continue
 

  end subroutine mapping_init

  subroutine mapping_allocate(n,map)
! -----------------------------------------------------------------------
! mapping : mapping allocate or reallocate 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
    implicit none
    integer,intent(in) :: n
    type(mapping),intent(inout) :: map
    logical :: realloc
    
    !-> number of points
    map%n=n

    !-> memory allocation
    if (.not.allocated(map%etha)) then
       allocate(map%etha(n,2))
    elseif(allocated(map%etha)) then
       deallocate(map%etha) ; allocate(map%etha(n,2))
    endif

    if (.not.allocated(map%detha)) then
       allocate(map%detha(n,2))
    elseif(allocated(map%detha)) then
       deallocate(map%detha) ; allocate(map%detha(n,2))
    endif

    if (.not.allocated(map%ddetha)) then
       allocate(map%ddetha(n,2))
    elseif(allocated(map%ddetha)) then
       deallocate(map%ddetha) ; allocate(map%ddetha(n,2))
    endif

    if (.not.allocated(map%dtetha)) then
       allocate(map%dtetha(n,2))
    elseif(allocated(map%dtetha)) then
       deallocate(map%dtetha) ; allocate(map%dtetha(n,2))
    endif
  end subroutine mapping_allocate

end module class_mapping
