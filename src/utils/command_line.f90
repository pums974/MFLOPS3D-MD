module command_line
  use precision
  implicit none

  !-> Define type containing command line argument
  type,public :: cmd_line
     !-> monodomain dimensions
     integer(ik) :: nx,ny,nz
     !-> number of domains in each directions
     integer(ik) :: ndx,ndy,ndz
     !-> periodicity in each directions
     integer(ik) :: periods(3)
     !-> reynolds number
     real(rk) :: reynolds
     !-> time step
     real(rk) :: ts
     !-> number of time iterations
     real(rk) :: ntime
     !-> nonlinear type : 1:convective, 2:skew-symetric
     integer(ik) :: nlt
     !-> projection type : 1:moin, 2:with pressure gradient
     integer(ik) :: pt
     !-> time order u
     integer(ik) :: tou
     !-> time order p
     integer(ik) :: top
     !-> pressure singular method : 0:nothing, 1:petsc nullspace
     !   2:dirichlet at one point 
     integer(ik) :: psm
     integer(ik) :: nsubite=1
  end type cmd_line


contains
!------------------------------------------------------------------------
! get command line informations
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine commandline(cmd)
    implicit none
    type(cmd_line),intent(inout) :: cmd 
    integer,parameter :: narg=100,ncheck=12
    character(50) :: cmdl(narg),getarg
    integer :: num,len,status,nargt,check(ncheck)

    !-> get command line arguments

    nargt=0
    do num=1,narg
       call get_command_argument(num,getarg,len,status)
       if (status /= 0) exit
       cmdl(num)=getarg
       nargt=nargt+1
    enddo

    do num=1,nargt
       if (cmdl(num)=='-h') then
          call usage
          stop
       endif
    enddo

    ! -> parsing of command line arguments

    check=1
    do num=1,nargt
       if (cmdl(num)=='-dim') then
          read(cmdl(num+1),*,err=10)cmd%nx,cmd%ny,cmd%nz ; check(1)=0
       endif
       if (cmdl(num)=='-dom') then
          read(cmdl(num+1),*,err=10)cmd%ndx,cmd%ndy,cmd%ndz ; check(2)=0
       endif
       if (cmdl(num)=='-period') then
          read(cmdl(num+1),*,err=10)cmd%periods(1),cmd%periods(2), &
               cmd%periods(3) ; check(3)=0
       endif
       if (cmdl(num)=='-reynolds') then
          read(cmdl(num+1),*,err=10)cmd%reynolds ; check(4)=0
       endif
       if (cmdl(num)=='-ts') then
          read(cmdl(num+1),*,err=10)cmd%ts ; check(5)=0
       endif
       if (cmdl(num)=='-ntime') then
          read(cmdl(num+1),*,err=10)cmd%ntime ; check(6)=0
       endif
       if (cmdl(num)=='-nlt') then
          read(cmdl(num+1),*,err=10)cmd%nlt ; check(7)=0
       endif
       if (cmdl(num)=='-pt') then
          read(cmdl(num+1),*,err=10)cmd%pt ; check(8)=0
       endif
       if (cmdl(num)=='-to') then
          read(cmdl(num+1),*,err=10)cmd%tou,cmd%top ; check(9)=0
       endif
       if (cmdl(num)=='-psm') then
          read(cmdl(num+1),*,err=10)cmd%psm ; check(10)=0
       endif
       if (cmdl(num)=='-nsub') then
          read(cmdl(num+1),*,err=10)cmd%nsubite ; check(11)=0
       endif
!       if (cmdl(num)=='-file') then
!          namel=cmdl(num+1) ; check(4)=0
!       endif
    enddo


!    if (sum(check)/=0.or.nx1<0.or.nx2<0) then
!       print*,check
!       call usage
!       stop
!    endif

    return

10  call usage
    stop
  end subroutine commandline
!------------------------------------------------------------------------
!  Give usage of command line
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine usage
    implicit none
    character(40) :: form
    form='(a40,a)'
    print form,"USAGE           ",": process options file"
    print form," "
    print form,"OPTIONS (*->mandatory)"
    print form," "
    print form,"-dim int,int,int        ",": longitudinal interval (*)"
  end subroutine usage


end module command_line
