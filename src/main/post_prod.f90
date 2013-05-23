program post_prod

  use mpi
!  use class_md
  use class_field
  use class_derivatives
  use class_mesh
!  use class_filter
  use color_print 
  use class_io
  use command_line

  implicit none

  integer(ik) :: ny
  real(rk),dimension(:),allocatable ::  moyu,moyv,moyw,rmsu,rmsv,rmsw,rmsuv,rmsuw,rmsvw,moyp,rmsp

  integer(ik) :: i,j,k,nd(3),coord(3),t
  integer(ik) :: nt
  type(cmd_line) :: cmd
!  type(mpi_data) :: mpi3
  type(mesh_grid) :: gridx, gridy, gridz
  type(field),allocatable :: u(:),v(:),w(:),p(:)
  real(rk) :: dx,dy,dz,lx,lz,utau,aux1,n,s,ref_re,ref_utau
  type(field) :: aux,deru
  type(derivatives_coefficients) :: dc
  character(20) ::t1
!  character(100) :: fich(nt)

  ref_re=3250
  ref_utau=0.57231059E-01
!  ref_re=20580
!  ref_utau=0.45390026E-01

  call commandline(cmd)
  nt=cmd%ntime
  allocate(u(nt),v(nt),w(nt),p(nt))
  !-> initialize mpi
!  call md_mpi_init(mpi3,cmd)
!  call md_mpi_getcoord(mpi3,coord)
!  call md_mpi_getnumberdomains(mpi3,nd)

!------------------- read mesh --------------------------------
  call read_mesh('results/grid_x','x',gridx,'gridx')
  call read_mesh('results/grid_y','y',gridy,'gridy')
  call read_mesh('results/grid_z','z',gridz,'gridz')

  ny=gridy%ny
  allocate(moyu(ny))
  allocate(moyv(ny))
  allocate(moyw(ny))
  allocate(moyp(ny))
  allocate(rmsp(ny))
  allocate(rmsu(ny))
  allocate(rmsv(ny))
  allocate(rmsw(ny))
  allocate(rmsuv(ny))
  allocate(rmsuw(ny))
  allocate(rmsvw(ny))

!------------------- read field --------------------------------
do t=1,nt
  write(t1,'(i20)') t
  call read_field("results/vel_u_"//trim(adjustl(t1)),u(t),'U')  
  call read_field("results/vel_v_"//trim(adjustl(t1)),v(t),'V')  
  call read_field("results/vel_w_"//trim(adjustl(t1)),w(t),'W')  
  call read_field("results/vel_p_"//trim(adjustl(t1)),p(t),'P')  
enddo

moyu=0._rk
moyv=0._rk
moyw=0._rk
moyp=0._rk
 do j=1,ny      ! calcul des moyennes
    n=0._rk
    do k=2,gridz%nz-1
      do i=2,gridx%nx-1
aux1=u(1)%f(i,j,k)
!        if(.not.(j.eq.1.or.j.eq.ny)) & 
!        aux1=min(abs(u%f(i,j,k)), &
!                abs(u%f(i+1,j,k)),abs(u%f(i-1,j,k)) , &
!                abs(u%f(i,j+1,k)),abs(u%f(i,j-1,k)) , &
!                abs(u%f(i,j,k+1)),abs(u%f(i,j,k-1)))  
        if(aux1.gt.1e-10.or.j.eq.1.or.j.eq.ny) then
        s=(gridx%grid1d(i+1)-gridx%grid1d(i-1))* &
          (gridz%grid1d(k+1)-gridz%grid1d(k-1))
do t=1,nt
          moyu(j)=moyu(j)+u(t)%f(i,j,k)*s
          moyv(j)=moyv(j)+v(t)%f(i,j,k)*s
          moyw(j)=moyw(j)+w(t)%f(i,j,k)*s
          moyp(j)=moyp(j)+p(t)%f(i,j,k)*s
          n=n+s
enddo
        endif
      enddo
    enddo
    moyu(j)=moyu(j)/n
    moyv(j)=moyv(j)/n
    moyw(j)=moyw(j)/n
    moyp(j)=moyp(j)/n
  enddo
do j=1,ny/2 ! compte tenu de la symmetricité
 moyu(j)=(moyu(j)+moyu(ny-j+1))*0.5_rk
 moyv(j)=(moyv(j)-moyv(ny-j+1))*0.5_rk
 moyw(j)=(moyw(j)+moyw(ny-j+1))*0.5_rk
 moyp(j)=(moyp(j)+moyp(ny-j+1))*0.5_rk
 moyu(ny-j+1)=moyu(j)
 moyv(ny-j+1)=-moyv(j)
 moyw(ny-j+1)=moyw(j)
 moyp(ny-j+1)=moyp(j)
enddo
  print*,'moy     : ',maxval(abs(moyv)),maxval(abs(moyw)),maxval(abs(moyp))

  ! calcul de utau=sqrt(dxu/nu)
  call derivatives_coefficients_init(gridy,dc,ny,cmd%so(2))
  call field_init(aux,"aux",ny,1,1)
  call field_init(deru,"deru",ny,1,1)
  aux%f(:,1,1)=moyu
  deru = derx(dc,aux)
  utau=sqrt((deru%f(1,1,1)-deru%f(ny,1,1))*0.5_rk/ref_re)
  print*,'utau    : ',sqrt(deru%f(1,1,1)/ref_re),sqrt(-deru%f(ny,1,1)/ref_re) ,ref_utau
  print*,'REYNOLDS: ',sqrt(deru%f(1,1,1)*ref_re),sqrt(-deru%f(ny,1,1)*ref_re) ,ref_utau*ref_re
!  utau=ref_utau

do t=1,nt
 u(t)%f = u(t)%f/utau !Normalisation
 v(t)%f = v(t)%f/utau
 w(t)%f = w(t)%f/utau
 p(t)%f = p(t)%f/(utau**2)
enddo
moyu=moyu/utau
moyv=moyv/utau
moyw=moyw/utau
moyp=moyp/(utau**2)

rmsu=0._rk ! calcul des rms
rmsv=0._rk
rmsw=0._rk
rmsp=0._rk
rmsuv=0._rk
rmsuw=0._rk
rmsvw=0._rk
  do j=1,ny
    n=0
    do k=2,gridz%nz-1
      do i=2,gridx%nx-1
aux1=u(1)%f(i,j,k)
!        if(.not.(j.eq.1.or.j.eq.ny)) & 
!        aux1=min(abs(u%f(i,j,k)), &
!                abs(u%f(i+1,j,k)),abs(u%f(i-1,j,k)) , &
!                abs(u%f(i,j+1,k)),abs(u%f(i,j-1,k)) , &
!                abs(u%f(i,j,k+1)),abs(u%f(i,j,k-1)))  
        if(aux1.gt.1e-10.or.j.eq.1.or.j.eq.ny) then
do t=1,nt
        s=(gridx%grid1d(i+1)-gridx%grid1d(i-1))* &
          (gridz%grid1d(k+1)-gridz%grid1d(k-1))
          rmsu(j)=rmsu(j)+s*((u(t)%f(i,j,k)-moyu(j)))**2
          rmsv(j)=rmsv(j)+s*((v(t)%f(i,j,k)-moyv(j)))**2
          rmsw(j)=rmsw(j)+s*((w(t)%f(i,j,k)-moyw(j)))**2
          rmsp(j)=rmsp(j)+s*((p(t)%f(i,j,k)-moyp(j)))**2
          rmsuv(j)=rmsuv(j)+((u(t)%f(i,j,k)-moyu(j)) &
                            *(v(t)%f(i,j,k)-moyv(j)))*s
          rmsuw(j)=rmsuw(j)+((u(t)%f(i,j,k)-moyu(j)) &
                            *(w(t)%f(i,j,k)-moyw(j)))*s
          rmsvw(j)=rmsvw(j)+((v(t)%f(i,j,k)-moyv(j)) &
                            *(w(t)%f(i,j,k)-moyw(j)))*s
          n=n+s
enddo
        endif
      enddo
    enddo
    rmsu(j)=sqrt(rmsu(j)/n)
    rmsv(j)=sqrt(rmsv(j)/n)
    rmsw(j)=sqrt(rmsw(j)/n)
    rmsp(j)=sqrt(rmsp(j)/n)
    rmsuv(j)=rmsuv(j)/n
    rmsuw(j)=rmsuw(j)/n
    rmsvw(j)=rmsvw(j)/n
  enddo

do j=1,ny/2 ! compte tenu de la symmetricité
 rmsu(j)=sqrt((rmsu(j)**2+rmsu(ny-j+1)**2)*0.5_rk)
 rmsv(j)=sqrt((rmsv(j)**2+rmsv(ny-j+1)**2)*0.5_rk)
 rmsw(j)=sqrt((rmsw(j)**2+rmsw(ny-j+1)**2)*0.5_rk)
 rmsp(j)=sqrt((rmsp(j)**2+rmsp(ny-j+1)**2)*0.5_rk)
 rmsuv(j)=(rmsuv(j)-rmsuv(ny-j+1))*0.5_rk
 rmsuw(j)=sqrt((rmsuw(j)**2+rmsuw(ny-j+1)**2)*0.5_rk)
 rmsvw(j)=(rmsvw(j)-rmsvw(ny-j+1))*0.5_rk
 rmsu(ny-j+1)=rmsu(j)
 rmsv(ny-j+1)=rmsv(j)
 rmsw(ny-j+1)=rmsw(j)
 rmsp(ny-j+1)=rmsp(j)
 rmsuv(ny-j+1)=-rmsuv(j)
 rmsuw(ny-j+1)=rmsuw(j)
 rmsvw(ny-j+1)=-rmsvw(j)
enddo

if(.false.) then
 print*,'rmss u  : ',maxval(abs(rmsu)),2.6590903E+00,2.6590903E+00/maxval(abs(rmsu))
 print*,'rmss v  : ',maxval(abs(rmsv)),8.4956920E-01,8.4956920E-01/maxval(abs(rmsv))
 print*,'rmss w  : ',maxval(abs(rmsw)),1.0907968E+00,1.0907968E+00/maxval(abs(rmsw))
 print*,'rmss p  : ',maxval(abs(rmsp))
 print*,'rms  uv :',maxval(abs(rmsuv)),7.3362684e-01,7.3362684e-01/maxval(abs(rmsuv))
 print*,'rms  uw :',maxval(abs(rmsuw)),7.3362684e-01,7.3362684e-01/maxval(abs(rmsuw))
 print*,'rms  vw :',maxval(abs(rmsvw)),7.3362684e-01,7.3362684e-01/maxval(abs(rmsvw))
else
 print*,'moy  u  : ',maxval(abs(moyu)) ,22.444668,22.444668/maxval(abs(moyu))
 print*,'rms  u  : ',maxval(abs(rmsu)) ,2.8271015,2.8271015/maxval(abs(rmsu))
 print*,'rms  v  : ',maxval(abs(rmsv)) ,1.0758162,1.0758162/maxval(abs(rmsv))
 print*,'rms  w  : ',maxval(abs(rmsw)) ,1.4507672,1.4507672/maxval(abs(rmsw))
 print*,'rms  p  : ',maxval(abs(rmsp)) ,2.9623945,2.9623945/maxval(abs(rmsp))
 print*,'rms  uv : ',maxval(abs(rmsuv)),0.89940619,0.89940619/maxval(abs(rmsuv))
 print*,'rms  uw : ',maxval(abs(rmsuw)),0.0086819138,0.0086819138/maxval(abs(rmsuw))
 print*,'rms  vw : ',maxval(abs(rmsvw)),0.0011614194,0.0011614194/maxval(abs(rmsvw))
endif

  do j=1,ny
    write(10,'(3e17.8)')  gridy%grid1d(j),moyu(j),rmsu(j)
    write(11,'(3e17.8)')  gridy%grid1d(j),moyv(j),rmsv(j)
    write(12,'(3e17.8)')  gridy%grid1d(j),moyw(j),rmsw(j)
    write(13,'(3e17.8)')  gridy%grid1d(j),moyv(j),rmsuv(j)
    write(14,'(3e17.8)')  gridy%grid1d(j),moyw(j),rmsuw(j)
    write(15,'(3e17.8)')  gridy%grid1d(j),moyv(j),rmsvw(j)
    write(16,'(3e17.8)')  gridy%grid1d(j),moyp(j),rmsp(j)
  enddo
!--------------------------------------------------------------

do t=1,nt
  call field_destroy(u(t))
  call field_destroy(v(t))
  call field_destroy(w(t))
  call field_destroy(p(t))
enddo
  deallocate(moyu,moyv,moyw,moyp)
  deallocate(rmsp,rmsu,rmsv,rmsw)
  deallocate(rmsuv,rmsuw,rmsvw)

end program  post_prod
