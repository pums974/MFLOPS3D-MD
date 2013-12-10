program test_md
  use class_field
  use class_mesh
  use class_derivatives
  use class_solver_3d
  use color_print
  use class_md
  use command_line
  implicit none
  type(mpi_data) :: mpid
  type(cmd_line) :: cmd
  type(mpi_inf_mat) :: inf
  integer(ik) :: i,j,k,nx,ny,nz
  type(mesh_grid) :: gridx,gridy,gridz
  type(solver_coeffs_3d) :: sch
  type(boundary_condition) :: bc
  type(field) :: u,g,du
  type(derivatives_coefficients) :: dcx,dcy,dcz
  integer(ik) :: bctype(6),dom_coord(3),neum
  real(rk) :: sigma,error,t,aux
  real(rk),allocatable :: uex(:,:,:)

  !-> get command line informations
  call commandline(cmd)

  !------------------------------------------------------------------------ 

  !-> initialize mpi
  call md_mpi_init(mpid,cmd)

  !-> initialize petsc
  call md_petsc_initialize()

  !-> put dimensions in variables for ease of use
  nx=mpid%nx ; ny=mpid%ny ; nz=mpid%nz

  !------------------------------------------------------------------------
  if (mpid%rank==0) then
     call color(ired) ; print'(a)','Initialize grid, field, solver : ' ; call color(color_off)
  endif
  !-> initialize mesh
  call mesh_init(gridx,'gridx','x',nx,ny,nz)
  call mesh_init(gridy,'gridy','y',nx,ny,nz)
  call mesh_init(gridz,'gridz','z',nx,ny,nz)

  !-> initialize grid
  call mesh_grid_init(gridx,'x',nx,1,1,mpid)
  call mesh_grid_init(gridy,'y',nx,ny,1,mpid)
  call mesh_grid_init(gridz,'z',1,1,nz,mpid)

  !-> define sigma : equal to zero -> poisson , different of zero -> helmholtz 
  sigma=0._rk

  !------------------------------------------------------------------------
  if (mpid%rank==0) then
     call color(ired) ; print'(a)','Start Initialize influence matrix : ' ; call color(color_off)
  endif
  !-> start initialization of influence matrix
  call influence_matrix_init_start(mpid,inf,sch,bc,u,g,sigma,dcx,dcy,dcz,'u')

  !-> initialize poisson solver coefficient
  neum=0
  bctype=(/1,1,1,1,1,1/)
!  bctype=(/2,2,2,2,2,2/)
  if (all(bctype==2)) neum=1
  call md_boundary_condition_init(mpid,inf,bctype)

  call solver_init_3d(gridx,gridy,gridz,sch,bctype,cmd%so(1))

  !-> initialize type field
  call field_init(u,"U",nx,ny,nz)
  call field_init(g,"G",nx,ny,nz)
  call field_init(du,"G",nx,ny,nz)

  !-> initialize type boundary_condition
  call boundary_condition_init(bc,nx,ny,nz)

  !-> initialisation of derivatives coefficients
  !call derivatives_coefficients_init(gridx,dcx,nx)
  !call derivatives_coefficients_init(gridy,dcy,ny)
  !call derivatives_coefficients_init(gridz,dcz,nz)
  call derivatives_coefficients_init(gridx,dcx,nx,solver='yes',so=cmd%so(2))
  call derivatives_coefficients_init(gridy,dcy,ny,solver='yes',so=cmd%so(2))
  call derivatives_coefficients_init(gridz,dcz,nz,solver='yes',so=cmd%so(2))

  !------------------------------------------------------------------------
  if (mpid%rank==0) then
     call color(ired) ; print'(a)','End Initialize influence matrix : ' ; call color(color_off)
  endif
  !-> end initialize influence matrix
  call influence_matrix_init_end(mpid,inf,sch,bc,u,g,sigma,dcx,dcy,dcz,'u')
!  call influence_matrix_init(mpid,inf,sch,bc,u,g,sigma,dcx,dcy,dcz)
  
  u%f=0._rk
  g%f=0._rk

  !-> initialize solver multidomain
  if (neum==1) then
     call md_solve_init(mpid,inf,null=1)
  else
     call md_solve_init(mpid,inf)
  endif
 
 !------------------------------------------------------------------------
  if (mpid%rank==0) then
     call color(ired) ; print'(a)','Initialize solution : ' ; call color(color_off)
  endif
  allocate(uex(nx,ny,nz))

  t=0._rk
  do i=1,3
     
     call solution_init(mpid,inf,gridx,gridy,gridz,u,g,uex,bc,bctype,sigma,t)
     t=t+0.01_rk
  !------------------------------------------------------------------------
     if (mpid%rank==0) then
        call color(ired) ; print'(a)','Solve influence matrix : ' ; call color(color_off)
     endif
     if (neum==1) then
        call multidomain_solve(mpid,inf,sch,bc,u,g,du,sigma,dcx,dcy,dcz,null=1)
     else
        call multidomain_solve(mpid,inf,sch,bc,u,g,du,sigma,dcx,dcy,dcz)
     endif

     if (i==2) call md_solve_guess_nonzero(mpid,inf)

     if (neum==1) then
        if (mpid%rank==0) then
           aux=(u%f(u%nx-1,u%ny-1,u%nz-1)-uex(u%nx-1,u%ny-1,u%nz-1))
        endif
        call md_mpi_bcast_double(mpid,aux,0)
        u%f=u%f-aux
     endif

     call err(u%f,uex,nx,ny,nz,error) 
     print*,'rank : ',mpid%rank,', error : ',error 

  enddo
!  call md_check(mpid,inf)

  !if (mpid%rank==0) then
     do j=2,u%ny-1
        do i=2,u%nx-1
        !   write(10+mpid%rank,*)gridx%grid1d(i),gridy%grid1d(j),uex(i,j,5),&
        !        u%f(i,j,5)
        enddo
        !write(10+mpid%rank,*)
     enddo
     j=9
     do k=1,u%nz
        do i=1,u%nx
           write(10+mpid%rank,*)gridx%grid1d(i),gridz%grid1d(k),uex(i,j,k),&
                u%f(i,j,k),abs(u%f(i,j,k)-uex(i,j,k))
        enddo
        write(10+mpid%rank,*)
     enddo
  !endif

  !------------------------------------------------------------------------
  !-> print various informations
  if (mpid%rank==0) then

     print*,'Size of influence matrix',inf%ninf     
     print*,'Number of interfaces (x,y,z,total)',inf%ix,inf%iy,inf%iz,inf%it
!     print*,'Number of interfaces points (x,y,z,total)',&
          !inf%ix*inf%ny*inf%nz,inf%iy*inf%nx*inf%nz,inf%iz*inf%nx*inf%ny,&
          !inf%ix*inf%ny*inf%nz+inf%iy*inf%nx*inf%nz+inf%iz*inf%nx*inf%ny
!     do i=1,inf%it
!        print*,'inf type,size',i,inf%dm_int_size(i,1),inf%dm_int_size(i,2),&
!             inf%dm_int_size(i,3)
!     enddo
  endif
!  print*,mpid%rank,'MATRIX SIZE',mpid%coord,inf%rows,inf%d_nz,'<->',inf%o_nnz
!  print*,'number of rows',mpid%coord+1,inf%rows
!  print*,'number of off diagonal',mpid%coord+1,inf%rows%dm_int
!  print*,mpid%rank,'<->',mpid%coord

  !------------------------------------------------------------------------
  if (mpid%rank==0) then
     call color(ired) ; print'(a)','Destroy influence matrix : ' ; call color(color_off)
  endif
  !-> deallocate influence matrix
  call md_influence_matrix_destroy(mpid,inf)


  call field_destroy(u)
  call field_destroy(g)
  call field_destroy(du)

  !------------------------------------------------------------------------
  !-> finalize petsc
  call md_petsc_finalize()
  !-> finalize mpi
  call md_mpi_finalize(mpid)

end program test_md

subroutine solution_init(mpid,inf,gridx,gridy,gridz,u,g,uex,bc,bctype,sigma,t)
  use class_field
  use class_mesh
  use class_md
  implicit none
  type(mpi_inf_mat) :: inf
  type(mpi_data) :: mpid
  type(mesh_grid) :: gridx,gridy,gridz
  type(field) :: u,g
  type(boundary_condition) :: bc
  real(rk) :: uex(u%nx,u%ny,u%nz)
  real(rk) :: x,y,z,sol,sigma,t
  integer(ik) :: bctype(6),i,j,k
  integer(ik) :: l,m,c(3),inter(3,2)

  u%f=0._rk
  g%f=0._rk
  bc%bcx=0._rk
  bc%bcy=0._rk
  bc%bcz=0._rk

!$OMP PARALLEL PRIVATE(i,j,k,x,y,z)
!$OMP DO SCHEDULE(RUNTIME)
  do k=1,u%nz
     do j=1,u%ny
        do i=1,u%nx
           x=gridx%grid1d(i)
           y=gridy%grid1d(j)
           z=gridz%grid1d(k)

           !-> exact function
           uex(i,j,k)=sol(x,y,z,t,'func')
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
           !uex(i,j,k)=sol(x,y,z,t,'func')

           !-> solution
           u%f(i,j,k)=sol(x,y,z,t,'func')

           !-> rhs
           g%f(i,j,k)=sol(x,y,z,t,'rhs')+sigma*sol(x,y,z,t,'func')

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
        if (bctype(1)==1) bc%bcx(j-1,k-1,1)=sol(x,y,z,t,'func')
        if (bctype(1)==2) bc%bcx(j-1,k-1,1)=sol(x,y,z,t,'derx')

        x=gridx%grid1d(u%nx)
        if (bctype(2)==1) bc%bcx(j-1,k-1,2)=sol(x,y,z,t,'func')
        if (bctype(2)==2) bc%bcx(j-1,k-1,2)=sol(x,y,z,t,'derx')
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
        if (bctype(3)==1) bc%bcy(i-1,k-1,1)=sol(x,y,z,t,'func')
        if (bctype(3)==2) bc%bcy(i-1,k-1,1)=sol(x,y,z,t,'dery')

        y=gridy%grid1d(u%ny)
        if (bctype(4)==1) bc%bcy(i-1,k-1,2)=sol(x,y,z,t,'func')
        if (bctype(4)==2) bc%bcy(i-1,k-1,2)=sol(x,y,z,t,'dery')
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
        if (bctype(5)==1) bc%bcz(i-1,j-1,1)=sol(x,y,z,t,'func')
        if (bctype(5)==2) bc%bcz(i-1,j-1,1)=sol(x,y,z,t,'derz')

        z=gridz%grid1d(u%nz)
        if (bctype(6)==1) bc%bcz(i-1,j-1,2)=sol(x,y,z,t,'func')
        if (bctype(6)==2) bc%bcz(i-1,j-1,2)=sol(x,y,z,t,'derz')
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  !-> set interface boundary condition to zero
  call md_mpi_getcoord(mpid,c)
  call md_get_interfaces_number(inf,c,inter)

  
  !-> m : left-right ; l : directions (x,y,z)
  do m=1,2
     do l=1,3
        if (inter(l,m)>0) then
           if (l==1.and.m==1) bc%bcx(:,:,1)=0._rk
           if (l==1.and.m==2) bc%bcx(:,:,2)=0._rk
           if (l==2.and.m==1) bc%bcy(:,:,1)=0._rk
           if (l==2.and.m==2) bc%bcy(:,:,2)=0._rk
           if (l==3.and.m==1) bc%bcz(:,:,1)=0._rk
           if (l==3.and.m==2) bc%bcz(:,:,2)=0._rk
        endif
     enddo
  enddo

!  u%f(1,2:u%ny-1,2:u%nz-1)=bc%bcx(:,:,1)
!  u%f(u%nx,2:u%ny-1,2:u%nz-1)=bc%bcx(:,:,2)
!  u%f(2:u%nx-1,1,2:u%nz-1)=bc%bcy(:,:,1)
!  u%f(2:u%nx-1,1,2:u%nz-1)=bc%bcy(:,:,2)
!  u%f(2:u%nx-1,2:u%ny-1,1)=bc%bcz(:,:,1)
!  u%f(2:u%nx-1,2:u%ny-1,1)=bc%bcz(:,:,2)

!  u%f(1,:,:)=0._rk
!  u%f(u%nx,:,:)=0._rk
!  u%f(:,1,:)=0._rk
!  u%f(:,u%ny,:)=0._rk
!  u%f(:,:,1)=0._rk
!  u%f(:,:,u%nz)=0._rk

end subroutine solution_init

function sol1(x,y,z,type)
! -----------------------------------------------------------------------
! exact solution : function, derivatives and rhs
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2012
!
  use parameters
  implicit none
  real(rk) :: sol,sol1
  real(rk) :: x,y,z,alp,bet,gam,pi
  character(*) :: type

  alp=1._rk ; bet=1._rk ; gam=2._rk
  pi=4._rk*atan(1._rk)
  
  if (type=="func") then
     sol=sin(2._rk*pi*gam*z)+5.
  endif
  if (type=="derx") then
     sol=0.
  endif
  if (type=="dery") then
     sol=0.
  endif
  if (type=="derz") then
     sol=2._rk*pi*gam*cos(2._rk*pi*gam*z)
  endif
  if (type=="rhs") then
     sol=&
          4._rk*pi*pi*(-gam*gam*sin(2._rk*pi*gam*z)) !+ &
          !12._rk*x+80._rk*z**3+36._rk*y**2
  endif

end function sol1

function sol(x,y,z,t,type)
! -----------------------------------------------------------------------
! exact solution : function, derivatives and rhs
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 05/2012
!
  use parameters
  implicit none
  real(rk) :: sol1,sol
  real(rk) :: x,y,z,t,alp,bet,gam,delt,pi
  character(*) :: type

  alp=2._rk ; bet=1._rk ; gam=1._rk ; delt=1._rk
  pi=4._rk*atan(1._rk)
  
  if (type=="func") then
     sol=sin(2._rk*pi*alp*x)+sin(2._rk*pi*bet*y)+sin(2._rk*pi*gam*z) &
          + cos(2._rk*pi*delt*t) &
          + 2._rk*x**3+3._rk*y**4+4._rk*z**5+5._rk
  endif
  if (type=="derx") then
     sol=2._rk*pi*alp*cos(2._rk*pi*alp*x) +6._rk*x**2
  endif
  if (type=="dery") then
     sol=2._rk*pi*bet*cos(2._rk*pi*bet*y) +12._rk*y**3
  endif
  if (type=="derz") then
     sol=2._rk*pi*gam*cos(2._rk*pi*gam*z) +20._rk*z**4
  endif
  if (type=="rhs") then
     sol=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*x))+ &
          4._rk*pi*pi*(-bet*bet*sin(2._rk*pi*bet*y))+ &
          4._rk*pi*pi*(-gam*gam*sin(2._rk*pi*gam*z)) + &
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
  som2=som2+sum(uex(2:nx-1,2:ny-1,1)**2)
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
