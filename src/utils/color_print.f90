module color_print
  implicit none

  ! define colors
  integer,parameter :: lenc=10

  ! Reset
  character(len=lenc),parameter :: color_off='\e[0m'       ! text reset

  ! regular colors
  character(len=lenc),parameter :: black='\e[0;30m'        ! black
  character(len=lenc),parameter :: red='\e[0;31m'          ! red
  character(len=lenc),parameter :: green='\e[0;32m'        ! green
  character(len=lenc),parameter :: yellow='\e[0;33m'       ! yellow
  character(len=lenc),parameter :: blue='\e[0;34m'         ! blue
  character(len=lenc),parameter :: purple='\e[0;35m'       ! purple
  character(len=lenc),parameter :: cyan='\e[0;36m'         ! cyan
  character(len=lenc),parameter :: white='\e[0;37m'        ! white

  ! bold
  character(len=lenc),parameter :: bblack='\e[1;30m'       ! black
  character(len=lenc),parameter :: bred='\e[1;31m'         ! red
  character(len=lenc),parameter :: bgreen='\e[1;32m'       ! green
  character(len=lenc),parameter :: byellow='\e[1;33m'      ! yellow
  character(len=lenc),parameter :: bblue='\e[1;34m'        ! blue
  character(len=lenc),parameter :: bpurple='\e[1;35m'      ! purple
  character(len=lenc),parameter :: bcyan='\e[1;36m'        ! cyan
  character(len=lenc),parameter :: bwhite='\e[1;37m'       ! white

  ! underline
  character(len=lenc),parameter :: ublack='\e[4;30m'       ! black
  character(len=lenc),parameter :: ured='\e[4;31m'         ! red
  character(len=lenc),parameter :: ugreen='\e[4;32m'       ! green
  character(len=lenc),parameter :: uyellow='\e[4;33m'      ! yellow
  character(len=lenc),parameter :: ublue='\e[4;34m'        ! blue
  character(len=lenc),parameter :: upurple='\e[4;35m'      ! purple
  character(len=lenc),parameter :: ucyan='\e[4;36m'        ! cyan
  character(len=lenc),parameter :: uwhite='\e[4;37m'       ! white

  ! background
  character(len=lenc),parameter :: on_black='\e[40m'       ! black
  character(len=lenc),parameter :: on_red='\e[41m'         ! red
  character(len=lenc),parameter :: on_green='\e[42m'       ! green
  character(len=lenc),parameter :: on_yellow='\e[43m'      ! yellow
  character(len=lenc),parameter :: on_blue='\e[44m'        ! blue
  character(len=lenc),parameter :: on_purple='\e[45m'      ! purple
  character(len=lenc),parameter :: on_cyan='\e[46m'        ! cyan
  character(len=lenc),parameter :: on_white='\e[47m'       ! white
  
  ! high intensity
  character(len=lenc),parameter :: iblack='\e[0;90m'       ! black
  character(len=lenc),parameter :: ired='\e[0;91m'         ! red
  character(len=lenc),parameter :: igreen='\e[0;92m'       ! green
  character(len=lenc),parameter :: iyellow='\e[0;93m'      ! yellow
  character(len=lenc),parameter :: iblue='\e[0;94m'        ! blue
  character(len=lenc),parameter :: ipurple='\e[0;95m'      ! purple
  character(len=lenc),parameter :: icyan='\e[0;96m'        ! cyan
  character(len=lenc),parameter :: iwhite='\e[0;97m'       ! white

  ! bold high intensity
  character(len=lenc),parameter :: biblack='\e[1;90m'      ! black
  character(len=lenc),parameter :: bired='\e[1;91m'        ! red
  character(len=lenc),parameter :: bigreen='\e[1;92m'      ! green
  character(len=lenc),parameter :: biyellow='\e[1;93m'     ! yellow
  character(len=lenc),parameter :: biblue='\e[1;94m'       ! blue
  character(len=lenc),parameter :: bipurple='\e[1;95m'     ! purple
  character(len=lenc),parameter :: bicyan='\e[1;96m'       ! cyan
  character(len=lenc),parameter :: biwhite='\e[1;97m'      ! white

  ! high intensity backgrounds
  character(len=lenc),parameter :: on_iblack='\e[0;100m'   ! black
  character(len=lenc),parameter :: on_ired='\e[0;101m'     ! red
  character(len=lenc),parameter :: on_igreen='\e[0;102m'   ! green
  character(len=lenc),parameter :: on_iyellow='\e[0;103m'  ! yellow
  character(len=lenc),parameter :: on_iblue='\e[0;104m'    ! blue
  character(len=lenc),parameter :: on_ipurple='\e[10;95m'  ! purple
  character(len=lenc),parameter :: on_icyan='\e[0;106m'    ! cyan
  character(len=lenc),parameter :: on_iwhite='\e[0;107m'   !  white

contains

  subroutine color(color_name)
    implicit none
    character(len=*) :: color_name
    
!    print*,"echo -ne '"//trim(color_name)//"'"
    call system("echo -ne '"//trim(color_name)//"'")

  end subroutine color

end module color_print
