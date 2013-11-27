program test_derivatives
  !  use parameters
  use class_field
  use class_mesh
  use class_derivatives
  use color_print
  implicit none

  integer(ik),parameter :: nt=1
  integer(ik) :: nx,ny,nz
  type(derivatives_coefficients) :: dc(5,5,3)
  type(mesh_grid) :: grid(5,3)
  type(field) :: u,v,w,p,r,p1,p2,u1,dxu

  integer(ik) :: i,j,k,l,m
  real(rk) :: x,y,z
  real(rk),allocatable :: tmpx(:),tmpy(:),tmpz(:)
  real(rk),allocatable :: uex(:,:,:,:),uey(:,:,:,:),uez(:,:,:,:)
  real(rk) :: xa,xb,pi,alp,bet,derror,dderror,test
  real(rk) :: sumderr,sumdderr,derror1,dderror1,errcoef=2._rk
  character(5) :: out 
  integer(ik) :: mindim=11,maxdim=4050,nerr

  ! define level of output files : LOW or FULL
  out='FULL'

  !-> initial dimensions
  nx=11 ; ny=mindim ; nz=mindim

        !-> initialize grid
        do i =1,5 !var
          j=0 ; k=0 ; l=0
          if(i==2) j=1
          if(i==3) k=1
          if(i==4) l=1
          if(i==5) j=-1
          if(i==5) k=-1
          if(i==5) l=-1
          call mesh_init(grid(i,1),'gridx','x',nx+j,1,1)
          call mesh_init(grid(i,2),'gridy','y',1,ny+k,1)
          call mesh_init(grid(i,3),'gridz','z',1,1,nz+l)
        enddo

        !-> define grid
        call mesh_grid_init(grid(1,1),'x',nx,1,1)
        call mesh_grid_init(grid(1,2),'y',1,ny,1)
        call mesh_grid_init(grid(1,3),'z',1,1,nz)

        grid(2,2:3)=grid(1,2:3)
        grid(3,1:3:2)=grid(1,1:3:2)
        grid(4,2:3)=grid(1,2:3)

        call mesh_grid_mac(grid(2,1),grid(1,1),1,'x')
        call mesh_grid_mac(grid(3,2),grid(1,2),1,'y')
        call mesh_grid_mac(grid(4,3),grid(1,3),1,'z')
        call mesh_grid_mac(grid(5,1),grid(1,1),-1,'x')
        call mesh_grid_mac(grid(5,2),grid(1,2),-1,'y')
        call mesh_grid_mac(grid(5,3),grid(1,3),-1,'z')

        !-> initialisation of derivatives coefficients
        do l=1,4
          call derivatives_coefficients_init(grid(l,1),dc(l,l,1),nx,6)
          call derivatives_coefficients_init(grid(l,2),dc(l,l,2),ny,6)
          call derivatives_coefficients_init(grid(l,3),dc(l,l,3),nz,6)
          if(l/=1) then
            call derivatives_coefficients_init_mac(grid(l,l-1),grid(1,l-1),dc(1,l,l-1))
            call derivatives_coefficients_init_mac(grid(1,l-1),grid(l,l-1),dc(l,1,l-1))
          endif
        enddo

      !-> initialize type field
        call field_init(p,"P",nx,ny,nz)
        call field_init(dxu,"dxu",nx,ny,nz)
        call field_init(p1,"P",nx,ny,nz)
        call field_init(p2,"P",nx,ny,nz)
        call field_init(u,"U",nx+1,ny,nz)
        call field_init(u1,"U",nx+1,ny,nz)
        call field_init(v,"V",nx,ny+1,nz)
        call field_init(w,"W",nx,ny,nz+1)
        call field_init(r,"R",nx-1,ny-1,nz-1)

        !-> define exact solution and derivatives
        pi=4._rk*atan(1._rk)
        do l=1,5
        do k=1,grid(l,3)%nz
           do j=1,grid(l,2)%ny
              do i=1,grid(l,1)%nx
                 !-> grid
                 x=grid(l,1)%grid1d(i) ; y=grid(l,2)%grid1d(j) ; z=grid(l,3)%grid1d(k)
                 !-> test function
                 select case(l)
                 case(1) 
                        p%f(i,j,k)= 5._rk * x**2
                        dxu%f(i,j,k)= pi*cos(2.3_rk*pi*x)*cos(1.3_rk*pi*y)
                 case(2) 
                        u%f(i,j,k)=  (1._rk/2.3_rk)*sin(2.3_rk*pi*x)*cos(1.3_rk*pi*y)!*2._rk*cos(0.5_rk*pi*z)
                 case(3)
                        v%f(i,j,k)= -(1._rk/1.3_rk)*cos(2.3_rk*pi*x)*sin(1.3_rk*pi*y)!*cos(0.5_rk*pi*z)
                 case(4)
                        w%f(i,j,k)= 0._rk !-(1._rk/0.5_rk)*cos(2.5_rk*pi*x)*cos(1.5_rk*pi*y)*sin(0.5_rk*pi*z)
                 case(5)
                        r%f(i,j,k)= 0._rk
                 end select
              enddo
           enddo
        enddo
        enddo

       call der_mac(dc(2,1,1),'x',u,p1)
       call write_field('grad_u',p1)  ;        print*,p1%f(p1%nx-2:p1%nx:2,5,5)
       p1=dxu-p1
       call write_field('grad_u',p1)  ;        print*,p1%f(p1%nx-2:p1%nx:2,5,5)

       call der_mac(dc(1,2,1),'x',p,u1)
       call write_field('grad_p',u1)  ;        print*,u1%f(u1%nx-2:u1%nx:2,5,5)
       call der_mac(dc(2,1,1),'x',u1,p1)
       call write_field('lap_p',p1)   ;        print*,p1%f(p1%nx-2:p1%nx:2,5,5)


       call der_mac(dc(2,1,1),'x',u,p1)
       call der_mac(dc(3,1,2),'y',v,p2)
       p=p1+p2
       call write_field('div_u',p)    ;        print*,maxval(abs(p%f)) 


end program test_derivatives 








!  print'(a)','Testing divergence : '
!  print'(a)','Testing gradient : '
!  print'(a)','Testing rotationnal : '



!  do m=1,3
!call color(ired)
!     if (m==1) print'(a)','Testing first direction : '
!     if (m==2) print'(a)','Testing second direction : '
!     if (m==3) print'(a)','Testing third direction : '
!     call color(color_off)
!     sumderr=0._rk ; sumdderr=0._rk ; nerr=0
!     do l=1,1000
!        !-> initialize work array
!        allocate(x(nx),y(ny),z(nz))
!        allocate(tmpx(nx),tmpy(ny),tmpz(nz))
!        allocate(uex(nx,ny,nz,2),uey(nx,ny,nz,2),uez(nx,ny,nz,2))

!        !-> initialize type field
!        call field_init(u(1),"U",nx,ny,nz)
!        call field_init(dxu,"DXU",nx,ny,nz)
!        call field_init(ddxu,"DDXU",nx,ny,nz)
!        call field_init(dyu,"DYU",nx,ny,nz)
!        call field_init(ddyu,"DDYU",nx,ny,nz)
!        call field_init(dzu,"DZU",nx,ny,nz)
!        call field_init(ddzu,"DDZU",nx,ny,nz)

!        !-> initialize grid
!        call mesh_init(gridx,'gridx','x',nx,1,1)
!        call mesh_init(gridy,'gridy','y',1,ny,1)
!        call mesh_init(gridz,'gridz','z',1,1,nz)
!        !-> define grid
!        call mesh_grid_init(gridx,'x',nx,1,1)
!        call mesh_grid_init(gridy,'y',1,ny,1)
!        call mesh_grid_init(gridz,'z',1,1,nz)

!        !-> define exact solution and derivatives
!        alp=1._rk ; bet=2._rk
!        pi=4._rk*atan(1._rk)

!        do k=1,u(1)%nz
!           do j=1,u(1)%ny
!              do i=1,u(1)%nx
!                 !-> grid
!                 x(i)=gridx%grid1d(i) ; y(j)=gridy%grid1d(j) ; z(k)=gridz%grid1d(k)
!                 !-> test function
!                 u(1)%f(i,j,k)=&
!                      sin(2._rk*pi*alp*x(i))+cos(2._rk*pi*bet*x(i))+ &
!                      sin(2._rk*pi*alp*y(j))+cos(2._rk*pi*bet*y(j))+ &
!                      sin(2._rk*pi*alp*z(k))+cos(2._rk*pi*bet*z(k))                
!                 !-> first derivative of test function
!                 uex(i,j,k,1)=2._rk*pi*(alp*cos(2._rk*pi*alp*x(i)) &
!                      -bet*sin(2._rk*pi*bet*x(i)))
!                 uey(i,j,k,1)=2._rk*pi*(alp*cos(2._rk*pi*alp*y(j)) &
!                      -bet*sin(2._rk*pi*bet*y(j)))
!                 uez(i,j,k,1)=2._rk*pi*(alp*cos(2._rk*pi*alp*z(k)) &
!                      -bet*sin(2._rk*pi*bet*z(k)))
!                 !-> second derivative of test function
!                 uex(i,j,k,2)=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*x(i)) &
!                      -bet*bet*cos(2._rk*pi*bet*x(i)))
!                 uey(i,j,k,2)=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*y(j)) &
!                      -bet*bet*cos(2._rk*pi*bet*y(j)))
!                 uez(i,j,k,2)=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*z(k)) &
!                      -bet*bet*cos(2._rk*pi*bet*z(k)))
!              enddo
!           enddo
!        enddo

!        !-> initialisation of derivatives coefficients
!        call derivatives_coefficients_init(gridx,dcx,nx)
!        call derivatives_coefficients_init(gridy,dcy,ny)
!        call derivatives_coefficients_init(gridz,dcz,nz)

!        !-> computation of derivatives 
!        dxu=derx(dcx,u(1)) ; ddxu=dderx(dcx,u(1))
!        dyu=dery(dcy,u(1)) ; ddyu=ddery(dcy,u(1))
!        dzu=derz(dcz,u(1)) ; ddzu=dderz(dcz,u(1))
!!        dxu=u(1)%derx(dcx) ; ddxu=u(1)%dderx(dcx)

!        !-> write results
!        if (out=='FULL') then
!           do i=1,nx
!              write(100,'(14e17.8)')x(i),u(1)%f(i,1,1),&
!                   dxu%f(i,1,1),uex(i,1,1,1),abs(uex(i,1,1,1)-dxu%f(i,1,1)), &
!                   ddxu%f(i,1,1),uex(i,1,1,2),abs(uex(i,1,1,2)-ddxu%f(i,1,1))
!           enddo
!           do j=1,ny
!              write(101,'(14e17.8)')y(j),u(1)%f(1,j,1),&
!                   dyu%f(1,j,1),uey(1,j,1,1),abs(uey(1,j,1,1)-dyu%f(1,j,1)), &
!                   ddyu%f(1,j,1),uey(1,j,1,2),abs(uey(1,j,1,2)-ddyu%f(1,j,1))
!           enddo
!           do k=1,nz
!              write(102,'(14e17.8)')z(k),u(1)%f(1,1,k),&
!                   dzu%f(1,1,k),uez(1,1,k,1),abs(uez(1,1,k,1)-dzu%f(1,1,k)), &
!                   ddzu%f(1,1,k),uez(1,1,k,2),abs(uez(1,1,k,2)-ddzu%f(1,1,k))
!           enddo
!           do i=1,nx-1
!              write(200,'(10e24.13)')x(i),x(i+1)-x(i) 
!           enddo
!        endif

!        !-> compute errors

!        if (m==1) then
!           if (l>1) derror1=derror
!           call err(dxu%f,uex(1,1,1,1),nx,ny,nz,derror)
!           call err(ddxu%f,uex(1,1,1,2),nx,ny,nz,dderror)
!           write(*,'(a,f8.2,a,3(a,i6),3(a,es9.2))')'Memory usage: ', &
!                nx*ny*nz*8._rk*7._rk/(1024._rk*1024._rk),' Mb, ', &
!                'Dimensions : ',nx,' x ',ny,' x ',nz, &
!                ', pasx=',gridx%pas,', dx error=',derror,', ddx error=',dderror
!           write(300,'(15es17.8)')real(nx),gridx%pas,derror,dderror
!           if (l>1.and.derror<derror1) then
!              nerr=nerr+1 ; sumderr=sumderr+log(derror1/derror)/log(errcoef) 
!           endif
!        endif

!        if (m==2) then
!           if (l>1) derror1=derror
!           call err(dyu%f,uey(1,1,1,1),nx,ny,nz,derror)
!           call err(ddyu%f,uey(1,1,1,2),nx,ny,nz,dderror)
!           write(*,'(a,f8.2,a,3(a,i6),3(a,es9.2))')'Memory usage: ', &
!                nx*ny*nz*8._rk*7._rk/(1024._rk*1024._rk),' Mb, ', &
!                'Dimensions : ',nx,' x ',ny,' x ',nz, & 
!                ', pasy=',gridy%pas,', dy error=',derror,', ddy error=',dderror
!           write(301,'(15es17.8)')real(ny),gridy%pas,derror,dderror
!           if (l>1.and.derror<derror1) then
!              nerr=nerr+1 ; sumderr=sumderr+log(derror1/derror)/log(errcoef) 
!           endif
!        endif

!        if (m==3) then
!           if (l>1) derror1=derror
!           call err(dzu%f,uez(1,1,1,1),nx,ny,nz,derror)
!           call err(ddzu%f,uez(1,1,1,2),nx,ny,nz,dderror)
!           write(*,'(a,f8.2,a,3(a,i6),3(a,es9.2))')'Memory usage: ', &
!                nx*ny*nz*8._rk*7._rk/(1024._rk*1024._rk),' Mb, ', &
!                'Dimensions : ',nx,' x ',ny,' x ',nz, & 
!                ', pasz=',gridz%pas,', dz error=',derror,', ddz error=',dderror
!           write(302,'(15es17.8)')real(nz),gridz%pas,derror,dderror
!           if (l>1.and.derror<derror1) then
!              nerr=nerr+1 ; sumderr=sumderr+log(derror1/derror)/log(errcoef) 
!           endif
!        endif

!        !-> deallocate work arrays
!        deallocate(x,y,z)
!        deallocate(tmpx,tmpy,tmpz)
!        deallocate(uex,uey,uez)

!        !-> deallocate tyep field
!        call field_destroy(u(1))
!        call field_destroy(dxu)
!        call field_destroy(ddxu)
!        call field_destroy(dyu)
!        call field_destroy(ddyu)
!        call field_destroy(dzu)
!        call field_destroy(ddzu)

!        !-> define dimensions
!        if (m==1.and.nx<maxdim) then 
!           nx=nx*int(errcoef)-1 ; ny=mindim ; nz=mindim 
!        elseif (m==1.and.nx>maxdim) then
!           nx=mindim ; print'(a,f5.2)','Approximate scheme order : ',sumderr/real(nerr,rk) ; exit
!        endif
!        if (m==2.and.ny<maxdim) then 
!           ny=ny*int(errcoef)-1 ; nx=mindim ; nz=mindim
!        elseif (m==2.and.ny>maxdim) then
!           ny=mindim ; print'(a,f5.2)','Approximate scheme order : ',sumderr/real(nerr,rk) ; exit
!        endif
!        if (m==3.and.nz<maxdim) then  
!           nz=nz*int(errcoef)-1 ; nx=mindim ; ny=mindim 
!        elseif (m==3.and.nz>maxdim) then
!           nz=mindim ; print'(a,f5.2)','Approximate scheme order : ',sumderr/real(nerr,rk) ; exit
!        endif
!     enddo
!  enddo

!end program test_derivatives

subroutine error_stop(error_mesg)
!  use mpi_utils, only : code,rang
  implicit none
  character(*) :: error_mesg
!  call MPI_FINALIZE(code)
!  if (rang==0) then
     print'(a)',error_mesg(1:len_trim(error_mesg))
!  endif
!  stop
  
end subroutine error_stop

subroutine err(u,uex,nx,ny,nz,error)
  use precision
  implicit none
  integer(ik) :: nx,ny,nz,i,j,k
  real(rk) :: u(nx,ny,nz),uex(nx,ny,nz)
  real(rk) :: som1,som2,error

  som1=0._rk
  som2=0._rk
  do k=1,nz
     do j=1,ny
        do i=1,nx
           som1=som1+(u(i,j,k)-uex(i,j,k))**2
           som2=som2+uex(i,j,k)**2
        enddo
     enddo
  enddo
  error=sqrt(som1)/sqrt(som2)
  !error=sqrt(som1/(nx*ny*nz))

end subroutine err
