program test_operators
  use class_field
  use color_print
  use command_line
  implicit none
  type(cmd_line) :: cmd
  type(field) :: u1,u2,u3
  integer(ik) :: i,j,k,nx,ny,nz
  real(rk) :: mem


  !-> get command line informations
  call commandline(cmd)

  !-> put dimensions in variables for ease of use
  nx=cmd%nx ; ny=cmd%ny ; nz=cmd%nz

  mem=(real(nx)*real(ny)*real(nz)*8._rk)/1024._rk
  print*,'Memory size of 1 field : ',mem
  

  !-> initialize type field
  call field_init(u1,"U1",nx,ny,nz)
  call field_init(u2,"U2",nx,ny,nz)
  call field_init(u3,"U3",nx,ny,nz)
  
  u1%f=2._rk
  u2%f=3._rk

  u3=u1+u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=-u1-u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1*u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1/u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1**2._rk+u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1+10._rk+u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1-10._rk-u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1*3._rk-3._rk*u2
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1/3._rk
  print*,trim(u3%name),u3%f(1,1,1)

  u3=u1*3._rk-u2/2._rk+u2*2.5_rk
  print*,trim(u3%name),u3%f(1,1,1)

end program test_operators

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
