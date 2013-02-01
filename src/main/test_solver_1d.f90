program test_solver
  !  use parameters
  use class_field
  use class_mesh
  use class_solver_1d
  use color_print
  implicit none

  integer(ik) :: nx,ny,nz,i,j,k
  type(solver_coefficients) :: scx
  type(mesh_grid) :: gridx
  type(field) :: u,g
  
  character(5) :: out 

  real(rk),allocatable :: uex(:,:,:)
  real(rk) :: x,alp,pi,error

  ! define level of output files : LOW or FULL
  out='LOW'

  !-> initial dimensions
  nx=12 ; ny=3 ; nz=3

  allocate(uex(nx,ny,nz))

  call color(ired) ; print'(a)','Initialization : ' ; call color(color_off)

  !-> initialize mesh
  call mesh_init(gridx,'gridx','x',nx,1,1)
  call mesh_grid_init(gridx,'x',nx,1,1)

  !-> initialize poisson solver coefficient
  call solver_coefficients_init_poisson(gridx,scx,nx)

  !-> initialize type field
  call field_init(u,"U",nx,ny,nz)
  call field_init(g,"G",nx,ny,nz)

  !-> exact solution
  
  alp=2._rk
  pi=4._rk*atan(1._rk)

  do k=1,u%nz
     do j=1,u%ny
        do i=1,u%nx
           x=gridx%grid1d(i)

           !-> exact function
           uex(i,j,k)=2._rk*x**3+sin(2._rk*pi*alp*x)+2._rk
           !uex(i,j,k)=sin(2._rk*pi*alp*x)

           !-> solution
           u%f(i,j,k)=2._rk*x**3+sin(2._rk*pi*alp*x)+2._rk
           !u%f(i,j,k)=sin(2._rk*pi*alp*x)

           !-> rhs
           g%f(i,j,k)=12._rk*x+4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*x))
           !g%f(i,j,k)=4._rk*pi*pi*(-alp*alp*sin(2._rk*pi*alp*x))

        enddo
     enddo
  enddo

  call color(ired) ; print'(a)','Solver : ' ; call color(color_off)

  ! solve poisson equation
  call solve_poisson(scx,g,u)

  do i=1,u%nx
     write(10,'(3es17.8)')gridx%grid1d(i),uex(i,2,2),u%f(i,2,2)
  enddo
  do i=1,u%nx-1
     write(20,'(3es17.8)')gridx%grid1d(i+1)-gridx%grid1d(i)
  enddo

  call err(u%f,uex,nx,ny,nz,error) 
  print*,gridx%pas,error
  

end program test_solver

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
