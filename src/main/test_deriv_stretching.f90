program statistics

  use mpi
  use class_md
  use class_field
  use class_derivatives
  use class_mesh
  use class_filter
  use color_print 
  use class_io
  use command_line

  implicit none

  integer(ik),parameter :: nxi=513
  integer(ik) :: nxr,nxs,nyi,nzi
  real(rk) :: fin,aux,pi,alpha,alpha1,tmp
  real(rk),allocatable :: norm_r(:),norm_s(:)

  integer(ik) :: i,j,k,nd(3),coord(3),w,f,rn
  integer(ik),allocatable :: clock(:)
  character(len=512) :: name_field='FIELD-014'
  type(cmd_line) :: cmd
  type(mpi_data) :: mpi3
  type(mesh_grid) :: gridx_i, gridx_r, gridx_s
  type(field) :: uo,ui,dur,fur,dus,fus,err_r,err_s,ur,us
  type(derivatives_coefficients) :: dcxi,dcxr,dcxs

  call commandline(cmd)
  nxr=cmd%nx
  nxs=cmd%nx
  nyi=cmd%ny!50
  nzi=cmd%nz!50
  !-> initialize mpi
  call md_mpi_init(mpi3,cmd)
  call md_mpi_getcoord(mpi3,coord)
  call md_mpi_getnumberdomains(mpi3,nd)
!------------------- read mesh --------------------------------
!  call color(bired) ; print'(a)','Read grid in netcdf file' ; call color(color_off)
!  call read_mesh('gridx_i','x',gridx_i,'gridx')

  pi=4._rk*atan(1._rk)
  call mesh_init(gridx_i,'gridx','x',nxi,1,1)
  call mesh_grid_init(gridx_i,'x',nxi,1,1,mpi3,.false.,2._rk*pi)

  fin=gridx_i%grid1d(nxi)
!------------------- read field --------------------------------
!  call color(bired) ; print'(a)','Read field in netcdf file' ; call color(color_off)
  call read_field("FIELD-014",uo,'velocity_x')  

!--------------------------------------------------------------
!  print *,'size field =',uo%nx,uo%ny,uo%nz

  call field_init(ui,"ui",nxi,uo%ny,uo%nz)
  ui%f(1:nxi-1,:,:) = uo%f(:,:,:)
  ui%f(nxi:nxi,:,:) = uo%f(1:1,:,:)
  call field_destroy(uo)


!  call color(bired) ; print'(a)','Generate initial field' ; call color(color_off)
!  if(mpi3%rank==0) then
!    CALL RANDOM_SEED(size=rn)
!    allocate(clock(rn))
!    CALL SYSTEM_CLOCK(COUNT=clock(1))
!    do i=rn,1-1
!      clock(i)=i*clock(1)
!    enddo
!    CALL RANDOM_SEED(PUT = clock)
!    deallocate(clock) 
!  endif

!  f=cmd%ny/2
!  call field_init(ui,'ui',nxi,nyi,nzi)

!  do k=1,nzi
!  do j=1,nyi
!  do w =1,f
!    if(mpi3%rank==0) call random_number(alpha) 
!    call md_mpi_bcast_double(mpi3,alpha,0)
!    ui%f(:,j,k)=ui%f(:,j,k)+cos(w*gridx_i%grid1d(:)+2._rk*pi*alpha) !*((float(w)/float(f))**(-sqrt(5._rk/3._rk)))
!  enddo
!  enddo
!  enddo

!--------------------------------------------------------------
!  call color(bired) ; print'(a)','Create other grids ans fields' ; call color(color_off)
  call mesh_init(gridx_s,'gridx','x',nxs,1,1)
  call mesh_grid_init(gridx_s,'x',nxs,1,1,mpi3,.true.,fin)
!  nxr=nint(1._rk+2._rk*(gridx_s%grid1d(nxs)-gridx_s%grid1d(1))/(gridx_s%grid1d(nxs/2+1)-gridx_s%grid1d(nxs/2-1)),ik)

  call mesh_init(gridx_r,'gridx','x',nxr,1,1)
  call mesh_grid_init(gridx_r,'x',nxr,1,1,mpi3,.false.,fin)

  call field_init(ur,"ur",nxr,nyi,nzi)
  call field_init(us,"us",nxs,nyi,nzi)
  call field_init(dur,"dur",nxr,nyi,nzi)
  call field_init(dus,"dus",nxs,nyi,nzi)
  call field_init(fus,"fus",nxs,nyi,nzi)
  call field_init(fur,"fur",nxr,nyi,nzi)

!------------------ interpol original field ------------------

!  call color(bired) ; print'(a)','Spline Interpolation' ; call color(color_off)

! do i=1,nxr
!    write(*,'(i4,2e17.8)')  mpi3%rank,gridx_r%grid1d(i),gridx_s%grid1d(i)
! enddo

  do k=1,nzi
  do j=1,nyi
!print*,'toto1'
    call interp_period_test(nxi,gridx_i%grid1d(:),ui%f(:,j,k), &
                            nxr,gridx_r%grid1d(:),ur%f(:,j,k),dur%f(:,j,k), &
                            nxs,gridx_s%grid1d(:),us%f(:,j,k),dus%f(:,j,k) )
!print*,'toto2'
  enddo
  enddo

!  call write_field("ur",ur)

!------------------ derivative by Finite Difference ------------

!  call color(bired) ; print'(a)','Derivation' ; call color(color_off)
  call derivatives_coefficients_init(gridx_r,dcxr,nxr,6)
  call derivatives_coefficients_init(gridx_s,dcxs,nxs,6)
  fur = derx(dcxr,ur)
  fus = derx(dcxs,us)

!  call mesh_grid_init(gridx_bloc_i,'x',nxi_bloc,1,1,1,tot=nxi,deb=(ibloc-1)*nxi_bloc,fin=fin)
!  call mesh_grid_init(gridx_bloc_r,'x',nxr_bloc,1,1,1,tot=nxr,deb=(ibloc-1)*nxr_bloc,fin=fin)
!  call mesh_grid_init(gridx_bloc_s,'x',nxs_bloc,1,1,2,tot=nxs,deb=(ibloc-1)*nxs_bloc,fin=fin)

! do i=1,nxs_bloc
!    write(*,'(i4,2e17.8)')  ibloc,gridx_bloc_r%grid1d(i),gridx_bloc_s%grid1d(i)
! enddo

!  gridx_r%grid1d((ibloc-1)*nxr_bloc+1:ibloc*nxr_bloc)=gridx_bloc_r%grid1d
!  gridx_s%grid1d((ibloc-1)*nxs_bloc+1:ibloc*nxs_bloc)=gridx_bloc_s%grid1d
!  gridx_r%grid3d((ibloc-1)*nxr_bloc+1:ibloc*nxr_bloc,:,:)=gridx_bloc_r%grid3d
!  gridx_s%grid3d((ibloc-1)*nxs_bloc+1:ibloc*nxs_bloc,:,:)=gridx_bloc_s%grid3d
!  dur%f((ibloc-1)*nxr_bloc+1:ibloc*nxr_bloc,:,:)=dur_bloc%f
!  dus%f((ibloc-1)*nxs_bloc+1:ibloc*nxs_bloc,:,:)=dus_bloc%f
!  fur%f((ibloc-1)*nxr_bloc+1:ibloc*nxr_bloc,:,:)=fur_bloc%f
!  fus%f((ibloc-1)*nxs_bloc+1:ibloc*nxs_bloc,:,:)=fus_bloc%f

! do i=1,nxs_bloc
!    write(*,'(i4,2e17.8)')  ibloc,dur_bloc%f(i,2,2),fur%f((ibloc-1)*nxs_bloc+i,2,2)
! enddo

!enddo


!------------------- error computation -------------------------

!  call color(bired) ; print'(a)','Error computation' ; call color(color_off)

!  call field_init(err_r,"err_r",nxr*nd(1),1,1)
!  call field_init(err_s,"err_s",nxs*nd(1),1,1)
!  allocate(norm_r(nxr*nd(1)))
!  allocate(norm_s(nxs*nd(1)))

!  err_r%f(:,:,:) = 0.0_rk
!  err_s%f(:,:,:) = 0.0_rk
!  norm_r(:) = 0.0_rk
!  norm_s(:) = 0.0_rk
!  
!  do k=1,nzi
!  do j=1,nyi
!   err_r%f(nxr*coord(1)+1:nxr*(coord(1)+1),1,1) = err_r%f(nxr*coord(1)+1:nxr*(coord(1)+1),1,1)+ (dur%f(:,j,k)-fur%f(:,j,k))**2
!   err_s%f(nxs*coord(1)+1:nxs*(coord(1)+1),1,1) = err_s%f(nxs*coord(1)+1:nxs*(coord(1)+1),1,1)+ (dus%f(:,j,k)-fus%f(:,j,k))**2
!   norm_r(nxr*coord(1)+1:nxr*(coord(1)+1)) = norm_r(nxr*coord(1)+1:nxr*(coord(1)+1))+dur%f(:,j,k)**2
!   norm_s(nxs*coord(1)+1:nxs*(coord(1)+1)) = norm_s(nxs*coord(1)+1:nxs*(coord(1)+1))+dus%f(:,j,k)**2

!  enddo
!  enddo

!  do i=1,nxr*nd(1)
!    aux=err_r%f(i,1,1)
!    call md_mpi_reduce_double(mpi3,aux,err_r%f(i,1,1))
!    call md_mpi_bcast_double(mpi3,err_r%f(i,1,1),0)
!    aux=norm_r(i)
!    call md_mpi_reduce_double(mpi3,aux,norm_r(i))
!    call md_mpi_bcast_double(mpi3,norm_r(i),0)
!  enddo
!  err_r%f(:,1,1) = sqrt(err_r%f(:,1,1))/ sqrt(norm_r(:))
!  do i=1,nxs*nd(1)
!    aux=err_s%f(i,1,1)
!    call md_mpi_reduce_double(mpi3,aux,err_s%f(i,1,1))
!    call md_mpi_bcast_double(mpi3,err_s%f(i,1,1),0)
!    aux=norm_s(i)
!    call md_mpi_reduce_double(mpi3,aux,norm_s(i))
!    call md_mpi_bcast_double(mpi3,norm_s(i),0)
!  enddo
!  err_s%f(:,1,1) = sqrt(err_s%f(:,1,1))/ sqrt(norm_s(:))


  call field_init(err_r,"err_r",nxr,1,1)
  call field_init(err_s,"err_s",nxs,1,1)
  allocate(norm_r(nxr))
  allocate(norm_s(nxs))

  err_r%f(:,:,:) = 0.0_rk
  err_s%f(:,:,:) = 0.0_rk
  norm_r(:) = 0.0_rk
  norm_s(:) = 0.0_rk
  
  do k=1,nzi
  do j=1,nyi
   err_r%f(:,1,1) = err_r%f(:,1,1)+ (dur%f(:,j,k)-fur%f(:,j,k))**2
   err_s%f(:,1,1) = err_s%f(:,1,1)+ (dus%f(:,j,k)-fus%f(:,j,k))**2
   norm_r(:) = norm_r(:)+dur%f(:,j,k)**2
   norm_s(:) = norm_s(:)+dus%f(:,j,k)**2

  enddo
  enddo
  err_r%f(:,1,1) = sqrt(err_r%f(:,1,1))/ sqrt(norm_r(:))
  err_s%f(:,1,1) = sqrt(err_s%f(:,1,1))/ sqrt(norm_s(:))

!do i=0,nbloc-1
!  if(i==mpi3%rank) 
!call read_field("FIELD-014",ui,'velocity_x',istart=(/1,1,1/),icount=(/nxi,nyi,nzi/))  
!  call MPI_BARRIER(MPI_COMM_WORLD,mpi3%code) 
!enddo

!  call color(bired) ; print'(a)','Export' ; call color(color_off)

!if(mpi3%rank==0) then 
!  call system('rm gridx_s.nc')
!  call system('rm gridx_r.nc')
!endif
!  call MPI_BARRIER(MPI_COMM_WORLD,mpi3%code) 
!  call write_mesh('gridx_r',gridx_r,istart=(/(ibloc-1)*nxr_bloc+1,1,1/),icount=(/nxr_bloc,1,1/))
!  call write_mesh('gridx_s',gridx_s,istart=(/(ibloc-1)*nxs_bloc+1,1,1/),icount=(/nxs_bloc,1,1/))
!  call write_field('gridx_r',err_r,istart=(/(ibloc-1)*nxr_bloc+1,1,1/),icount=(/nxr_bloc,1,1/))
!  call write_field('gridx_s',err_s,istart=(/(ibloc-1)*nxs_bloc+1,1,1/),icount=(/nxs_bloc,1,1/))

! do i=1,nxr
!    write(10+mpi3%rank,'(2i5,3e17.8)')  mpi3%rank,i,gridx_r%grid3d(i,1,1),err_r%f(i,1,1),norm_r(i)
! enddo
! close(10+mpi3%rank)
! do i=1,nxs
!    write(20+mpi3%rank,'(2i5,3e17.8)')  mpi3%rank,i,gridx_s%grid3d(i,1,1),err_s%f(i,1,1),norm_r(i)
! enddo
! close(20+mpi3%rank)

!do i=2,nxs
!    write(30+mpi3%rank,'(3e17.8)') & ! (ibloc-1)*nxr_bloc+1+i, &
!               (gridx_s%grid3d(i,1,1)+gridx_s%grid3d(i-1,1,1))*0.5_rk, &
!               gridx_r%grid3d(i,1,1)-gridx_r%grid3d(i-1,1,1), &
!               gridx_s%grid3d(i,1,1)-gridx_s%grid3d(i-1,1,1)
! enddo
! close(30+mpi3%rank)

if(mpi3%rank==0) then 
  i=(nxr+1)/2
  aux    = err_r%f(i,1,1)
  alpha =(gridx_r%grid3d(i+1,1,1)-gridx_r%grid3d(i-1,1,1))*0.5_rk
  fin=err_r%f(i,1,1)

  i=(nxs+1)/2
  alpha1 = err_s%f(i,1,1)
  fin=err_s%f(i,1,1) / fin


  i=nxr
  !print*,i, gridx_s%grid3d(i,1,1)-gridx_s%grid3d(i-1,1,1)          ,err_s%f(i,1,1)
  aux    = err_r%f(i,1,1) / aux
  alpha  = alpha / (gridx_r%grid3d(i,1,1)-gridx_r%grid3d(i-1,1,1))

  i=nxs
  alpha1 = err_s%f(i,1,1) / alpha1

  alpha=(gridx_s%grid3d(nxs/2+1,1,1)-gridx_s%grid3d(nxs/2-1,1,1))/(2._rk*(gridx_s%grid3d(nxs,1,1)-gridx_s%grid3d(nxs-1,1,1)))

  tmp=(gridx_s%grid3d(nxs-1,1,1)-gridx_s%grid3d(nxs-2,1,1))/(gridx_s%grid3d(nxs,1,1)-gridx_s%grid3d(nxs-1,1,1))

!do i=nxs,nxs/2,-1
!  alpha1=alpha/ (gridx_s%grid3d(i,1,1)-gridx_s%grid3d(i-1,1,1))
!  if(alpha1.lt.1.01_rk) exit
!enddo

!!  write(*,'(5e17.8)',advance='no') tmp,alpha
!  if(fin>2._rk)   call color(bired)
!  write(*,'(5e17.8)',advance='no')  fin
!  if(fin>2._rk)   call color(color_off)
!!  write(*,'(5e17.8)',advance='no') aux
!  if(alpha1<0.1_rk)   call color(bired)
!  if(alpha1>10._rk)   call color(bired)
!  write(*,'(5e17.8)',advance='no') alpha1
!  if(alpha1>10._rk)    call color(color_off)
!  if(alpha1<0.1_rk)    call color(color_off)

  write(*,'(5e17.8)',advance='no') (fin-1._rk) ,(alpha1-1._rk)
endif

!--------------------------------------------------------------


  call field_destroy(ui)

!  call field_destroy(ur)
  call field_destroy(dur)
  call field_destroy(fur)

!  call field_destroy(us)
  call field_destroy(dus)
  call field_destroy(fus)

  call MPI_BARRIER(MPI_COMM_WORLD,mpi3%code) 
  call MPI_FINALIZE(mpi3%code)

end program  statistics

subroutine error_stop(error_mesg)
  implicit none
  character(*) :: error_mesg
     print'(a)',error_mesg(1:len_trim(error_mesg))
  
end subroutine error_stop




