program test_io
  use class_field
  use class_mesh
  use color_print
  implicit none

  integer(ik),parameter :: nt=1
  integer(ik) :: nx,ny,nz
  type(mesh_grid) :: gridx,gridy,gridz
  type(field) :: u(nt),v(nt),w(nt)

  integer(ik) :: i,j,k,l,m
  real(rk) :: r,p,alp,bet,pi
  real(rk) :: x,y,z

  !-> initial dimensions
  nx=11 ; ny=11 ; nz=11

  !-> initialize fields
  call field_init(u(1),"U",nx,ny,nz)
  call field_init(v(1),"V",nx,ny,nz)
  call field_init(w(1),"W",2*nx,2*ny,2*nz,n1n='dimx',n2n='dimy',n3n='dimz')
  
  !-> initialize grids
  call mesh_init(gridx,'gridx','x',nx,ny,nz)
  call mesh_init(gridy,'gridy','y',nx,ny,nz)
  call mesh_init(gridz,'gridz','z',nx,ny,nz)
  !-> define grid
  call mesh_grid_init(gridx,'x',nx,1,1)
  call mesh_grid_init(gridy,'y',nx,ny,1)
  call mesh_grid_init(gridz,'z',1,1,nz)

  !-> define exact solution and derivatives
  alp=0._rk ; bet=2._rk
  pi=4._rk*atan(1._rk)

  do k=1,u(1)%nz
     do j=1,u(1)%ny
        do i=1,u(1)%nx
           !-> grid
           x=gridx%grid1d(i) ; y=gridy%grid1d(j) ; z=gridz%grid1d(k)
           !-> test function
           u(1)%f(i,j,k)=&
                sin(2._rk*pi*alp*x)+cos(2._rk*pi*bet*x)+ &
                sin(2._rk*pi*alp*y)+cos(2._rk*pi*bet*y)+ &
                sin(2._rk*pi*alp*z)+cos(2._rk*pi*bet*z)  
           v(1)%f(i,j,k)=2.*(&
                sin(2._rk*pi*alp*x)+cos(2._rk*pi*bet*x)+ &
                sin(2._rk*pi*alp*y)+cos(2._rk*pi*bet*y)+ &
                sin(2._rk*pi*alp*z)+cos(2._rk*pi*bet*z)) 
        enddo
     enddo
  enddo


  !-> create netcdf file and write grid
  call color(bired) ; print'(a)','Write grid in netcdf file' ; call color(color_off)
  call write_mesh('grid_x',gridx)
  call write_mesh('grid_y',gridy)
  call write_mesh('grid_z',gridz)

  !-> create netcdf file and write U
  call color(bired) ; print'(a)','Write U in netcdf file' ; call color(color_off)
  call write_field('field',u(1))

  !-> append V variable in netcdf file
  call color(bired) ; print'(a)','Append V in netcdf file' ; call color(color_off)
  call write_field('field',v(1))

  !-> append W variable in netcdf file
  call color(bired) ; print'(a)','Append W in netcdf file' ; call color(color_off)
  call write_field('field',w(1))

  !-> read U variable from netcdf file
  call color(bired) ; print'(a)','Read U in netcdf file' ; call color(color_off)
  call read_field('field',u(1),'W')
  call write_field('field_new',u(1))

  !-> read grid from netcdf file
  gridx%grid1d=0._rk
  gridy%grid1d=0._rk
  gridz%grid1d=0._rk
  call color(bired) ; print'(a)','Read grid in netcdf file' ; call color(color_off)
  call read_mesh('grid_x','x',gridx,'gridx')
  call read_mesh('grid_y','y',gridy,'gridy')
  call read_mesh('grid_z','z',gridz,'gridz')
  
  print*,gridx%grid1d
  print*,gridy%grid1d
  print*,gridz%grid1d

  !-> destroy fieds
  call field_destroy(u(1))
  call field_destroy(v(1))
  call field_destroy(w(1))

end program test_io

subroutine error_stop(error_mesg)
!  use mpi_utils, only : code,rang
  implicit none
  character(*) :: error_mesg
!  call MPI_FINALIZE(code)
!  if (rang==0) then
     print'(a)',error_mesg(1:len_trim(error_mesg))
!  endif
  stop
  
end subroutine error_stop
