program test_solver
  !$ use OMP_LIB
  !  use parameters
  use class_field
  use class_mesh
  use class_solver_3d
  use color_print
  use command_line
  implicit none

  integer(ik) :: n,nx,ny,nz,i,j,k
  type(solver_coeffs_3d) :: scx
  type(mesh_grid) :: gridx,gridy,gridz
  type(field) :: u,g
  type(boundary_condition) :: bc
  type(cmd_line) :: cmd
  
  character(5) :: out 

  real(rk),allocatable :: uex(:,:,:)
  real(rk) :: x,y,z,error,op,time,sol,sigma
  integer(8) :: t1,t2,irate
  integer(ik) :: bctype(6)

  !-> get command line informations
  call commandline(cmd)

  ! define level of output files : LOW or FULL
  out='LOW'

  !-> define sigma : equal to zero -> poisson , different of zero -> helmholtz 
  sigma=2._rk

  !-> initial dimensions
  nx=cmd%nx ; ny=cmd%ny ; nz=cmd%nz

  allocate(uex(nx,ny,nz))

  call color(ired) ; print'(a)','Initialization : ' ; call color(color_off)

  !-> initialize mesh
  call mesh_init(gridx,'gridx','x',nx,1,1)
  call mesh_init(gridy,'gridy','y',nx,ny,1)
  call mesh_init(gridz,'gridz','z',1,1,nz)

  !-> initialize grid
  call mesh_grid_init(gridx,'x',nx,1,1)
  call mesh_grid_init(gridy,'y',nx,ny,1)
  call mesh_grid_init(gridz,'z',1,1,nz)

  !-> initialize poisson solver coefficient
  bctype=(/1,1,1,1,1,1/)
  call solver_init_3d(gridx,gridy,gridz,scx,bctype)

  !-> initialize type field
  call field_init(u,"U",nx,ny,nz)
  call field_init(g,"G",nx,ny,nz)

  !-> initialize type boundary_condition
  call boundary_condition_init(bc,nx,ny,nz)

  !-> exact solution and rhs
  call solution_init(gridx,gridy,gridz,u,g,uex,bc,bctype,sigma)

  ! solve poisson equation  
  call color(ired) ; print'(a)','Solver : ' ; call color(color_off)

  call system_clock(t1,irate)
  call solver_3d(scx,g,u,bc,sigma)
  call system_clock(t2,irate)

  !-> rescale solution for all neumann bc
  if (all(bctype==2)) then
     u%f(2:u%nx-1,2:u%ny-1,2:u%nz-1)=u%f(2:u%nx-1,2:u%ny-1,2:u%nz-1) &
          +(uex(u%nx/2,u%ny/2,u%nz/2)-u%f(u%nx/2,u%ny/2,u%nz/2))
       !+sum(uex(2:u%nx-1,2:u%ny-1,2:u%nz-1)-u%f(2:u%nx-1,2:u%ny-1,2:u%nz-1)) &
       !     /((u%nx-2)*(u%ny-2)*(u%nz-2))
  endif

  do i=1,u%nx
     write(10,'(5es17.8)')gridx%grid1d(i),uex(i,2,2),u%f(i,2,2)
  enddo
  do j=1,u%ny
     write(11,'(5es17.8)')gridy%grid1d(j),uex(2,j,2),u%f(2,j,2)
  enddo

  do j=1,u%ny
     do i=1,u%nx
        write(100,'(10es17.8)')gridx%grid1d(i),gridy%grid1d(j),uex(i,j,2),u%f(i,j,2),&
             uex(i,j,2)-u%f(i,j,2)
     enddo
     write(100,*)
  enddo

  do i=1,u%nx-1
     write(20,'(3es17.8)')gridx%grid1d(i+1)-gridx%grid1d(i)
  enddo
  do i=1,u%ny-1
     write(21,'(3es17.8)')gridy%grid1d(i+1)-gridy%grid1d(i)
  enddo
  do i=1,u%nz-1
     write(22,'(3es17.8)')gridz%grid1d(i+1)-gridz%grid1d(i)
  enddo

  call err(u%f,uex,nx,ny,nz,error) 
  op=2*real(nx)*real(ny)*real(nz)*(real(nx)+real(ny)+real(nz)) &
       +3*real(nx)*real(ny)*real(nz) &
       +12._rk*(real(ny)*real(nz)+real(nx)*real(nz)+real(nx)*real(ny))
  time=real(t2-t1)/real(irate)
  print'(3es17.8,f10.4,f10.4)',gridx%pas,error,op,time,(op/(1000._rk*1000._rk))/time

end program test_solver

subroutine solution_init(gridx,gridy,gridz,u,g,uex,bc,bctype,sigma)
  use class_field
  use class_mesh
  implicit none
  type(mesh_grid) :: gridx,gridy,gridz
  type(field) :: u,g
  type(boundary_condition) :: bc
  real(rk) :: uex(u%nx,u%ny,u%nz)
  real(rk) :: x,y,z,sol,sigma
  integer(ik) :: bctype(6),i,j,k
  
!$OMP PARALLEL PRIVATE(i,j,k,x,y,z)
!$OMP DO SCHEDULE(RUNTIME)
  do k=1,u%nz
     do j=1,u%ny
        do i=1,u%nx
           x=gridx%grid1d(i)
           y=gridy%grid1d(j)
           z=gridz%grid1d(k)
           !-> exact function
           uex(i,j,k)=sol(x,y,z,'func')
        enddo
     enddo
  enddo
!$OMP END DO
!$OMP DO SCHEDULE(RUNTIME)
  do k=2,u%nz-1
     do j=2,u%ny-1
        do i=2,u%nx-1
           x=gridx%grid1d(i)
           y=gridy%grid1d(j)
           z=gridz%grid1d(k)

           !-> exact function
           !uex(i,j,k)=sol(x,y,z,'func')

           !-> solution
           u%f(i,j,k)=sol(x,y,z,'func')

           !-> rhs
           g%f(i,j,k)=sol(x,y,z,'rhs')+sigma*sol(x,y,z,'func')

        enddo
     enddo
  enddo
!$OMP END DO

  !-> boundary condition
  !-> x-direction
!$OMP DO SCHEDULE(RUNTIME)
  do k=2,u%nz-1
     do j=2,u%ny-1
        y=gridy%grid1d(j)
        z=gridz%grid1d(k)

        x=gridx%grid1d(1)
        if (bctype(1)==1) bc%bcx(j-1,k-1,1)=sol(x,y,z,'func')
        if (bctype(1)==2) bc%bcx(j-1,k-1,1)=sol(x,y,z,'derx')

        x=gridx%grid1d(u%nx)
        if (bctype(2)==1) bc%bcx(j-1,k-1,2)=sol(x,y,z,'func')
        if (bctype(2)==2) bc%bcx(j-1,k-1,2)=sol(x,y,z,'derx')
     enddo
  enddo
!$OMP END DO
  !-> y-direction
!$OMP DO SCHEDULE(RUNTIME)
  do k=2,u%nz-1
     do i=2,u%nx-1
        x=gridx%grid1d(i)
        z=gridz%grid1d(k)

        y=gridy%grid1d(1)
        if (bctype(3)==1) bc%bcy(i-1,k-1,1)=sol(x,y,z,'func')
        if (bctype(3)==2) bc%bcy(i-1,k-1,1)=sol(x,y,z,'dery')

        y=gridy%grid1d(u%ny)
        if (bctype(4)==1) bc%bcy(i-1,k-1,2)=sol(x,y,z,'func')
        if (bctype(4)==2) bc%bcy(i-1,k-1,2)=sol(x,y,z,'dery')
     enddo
  enddo
!$OMP END DO
  !-> z-direction
!$OMP DO SCHEDULE(RUNTIME)
  do j=2,u%ny-1
     do i=2,u%nx-1
        x=gridx%grid1d(i)
        y=gridy%grid1d(j)

        z=gridz%grid1d(1)
        if (bctype(5)==1) bc%bcz(i-1,j-1,1)=sol(x,y,z,'func')
        if (bctype(5)==2) bc%bcz(i-1,j-1,1)=sol(x,y,z,'derz')

        z=gridz%grid1d(u%nz)
        if (bctype(6)==1) bc%bcz(i-1,j-1,2)=sol(x,y,z,'func')
        if (bctype(6)==2) bc%bcz(i-1,j-1,2)=sol(x,y,z,'derz')
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

end subroutine solution_init

function sol(x,y,z,type)
! -----------------------------------------------------------------------
! exact solution : function, derivatives and rhs
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2012
!
  use parameters
  implicit none
  real(rk) :: sol
  real(rk) :: x,y,z,alp,bet,gam,pi
  character(*) :: type

  alp=2._rk ; bet=2._rk ; gam=2._rk
  pi=4._rk*atan(1._rk)
  
  if (type=="func") then
     sol=sin(2._rk*pi*alp*x)+sin(2._rk*pi*bet*y)+sin(2._rk*pi*gam*z) &
          + 2._rk*x**3+3._rk*y**4+4._rk*z**5
  endif
  if (type=="derx") then
     sol=2._rk*pi*alp*cos(2._rk*pi*alp*x)+6._rk*x**2
  endif
  if (type=="dery") then
     sol=2._rk*pi*bet*cos(2._rk*pi*bet*y)+12._rk*y**3
  endif
  if (type=="derz") then
     sol=2._rk*pi*gam*cos(2._rk*pi*gam*z)+20._rk*z**4
  endif
  if (type=="rhs") then
     sol=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*x))+ &
          4._rk*pi*pi*(-bet*bet*sin(2._rk*pi*bet*y))+ &
          4._rk*pi*pi*(-gam*gam*sin(2._rk*pi*gam*z))+ &
          12._rk*x+80._rk*z**3+36._rk*y**2
  endif

end function sol


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
  do k=2,nz-1
     do j=2,ny-1
        do i=2,nx-1
           som1=som1+(u(i,j,k)-uex(i,j,k))**2
           som2=som2+uex(i,j,k)**2
        enddo
     enddo
  enddo

  !-> x boundary
  som1=som1+sum((u(1,2:ny-1,2:nz-1)-uex(1,2:ny-1,2:nz-1))**2)
  som2=som2+sum(uex(1,2:ny-1,2:nz-1)**2)
  som1=som1+sum((u(nx,2:ny-1,2:nz-1)-uex(nx,2:ny-1,2:nz-1))**2)
  som2=som2+sum(uex(nx,2:ny-1,2:nz-1)**2)

  !-> y boundary
  som1=som1+sum((u(2:nx-1,1,2:nz-1)-uex(2:nx-1,1,2:nz-1))**2)
  som2=som2+sum(uex(2:nx-1,1,2:nz-1)**2)
  som1=som1+sum((u(2:nx-1,ny,2:nz-1)-uex(2:nx-1,ny,2:nz-1))**2)
  som2=som2+sum(uex(2:nx-1,ny,2:nz-1)**2)

  !-> z boundary
  som1=som1+sum((u(2:nx-1,2:ny-1,1)-uex(2:nx-1,2:ny-1,1))**2)
  som2=som2+sum(uex(2:nx-1,2:nz-1,1)**2)
  som1=som1+sum((u(2:nx-1,2:ny-1,nz)-uex(2:nx-1,2:ny-1,nz))**2)
  som2=som2+sum(uex(2:nx-1,2:ny-1,nz)**2)

  error=sqrt(som1)/sqrt(som2)

!  som1=0._rk
!  som2=0._rk
!  do k=1,nz
!     do j=1,ny
!        do i=1,nx
!           som1=som1+(u(i,j,k)-uex(i,j,k))**2
!           som2=som2+uex(i,j,k)**2
!        enddo
!     enddo
!  enddo
!  error=sqrt(som1)/sqrt(som2)
  !error=sqrt(som1/(nx*ny*nz))

end subroutine err
