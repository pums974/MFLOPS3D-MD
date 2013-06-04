program testnavier3d
  use command_line
  use class_navier_3D
  use class_io
  use netcdf
  use mpi
  implicit none
  type(cmd_line)  :: cmd
  type(navier3d)  :: nav
  type(field)     :: uo,ui,us,vo,wo,po
  type(mesh_grid) :: gridxi,gridyi,gridzi,grids
  type(mpi_data)  :: mpid
  integer(ik) :: nx,ny,nz
  integer(ik) :: ite
  integer(ik) :: i,j,k,iaux,i0,j0,k0,i1,j1,k1,i2,j2,k2
  real(rk),allocatable :: uex(:,:,:,:), pex(:,:,:),vectorerror(:,:,:,:)
  real(rk) ::x,y,z,t,x1,x2,y1,y2,z1,z2
  real(rk) :: dx,dy,dz
  real(rk) ::error,aux,time,errort,ref,reft,err_vts_t,err_pre_t,init(2,49)
  integer(8) :: t1,t2,irate,subite,t3
  logical :: test

    integer(ik) :: dim_len(3)
    character(len=512) :: dim_name(3)
    integer(ik) :: ncid
    integer(ik) :: startv(3),countv(3)
    integer(ik) :: varid(1)

  !-> get command line informations
  call commandline(cmd)

  allocate(uex(cmd%nx,cmd%ny,cmd%nz,3),pex(cmd%nx,cmd%ny,cmd%nz))

  allocate(vectorerror(cmd%nx,cmd%ny,cmd%nz,5))
  !------------------------------------------------------------------------ 
  !-> pre-computation
  !------------------------------------------------------------------------ 
  call navier_initialization(cmd,mpid,nav)

  !-> read restart files if they exists
!  call restart_read(mpid,nav,test)

!initialisation
if (.false..and..not.test) then
 call initialise_navier(nav,mpid)
endif
  
  call navier_write_fields(mpid,nav,1,0)

  !------------------------------------------------------------------------ 
  !-> time loop
  !------------------------------------------------------------------------ 

  err_vts_t=0._rk
  err_pre_t=0._rk
  call system_clock(t1,irate)
temps:  do ite=1,nav%ntime

  call system_clock(t3,irate)
     if (ite==20) then
       call system_clock(t1,irate)
       err_vts_t=0._rk
       err_pre_t=0._rk
     endif

     !-> time update
     call navier_time(nav)

     if (int((ite-1)/10)*10.eq.(ite-1).and.mpid%rank==0) write(*,'(17a)', advance= 'yes') &
     "   it","    time","  subit","  it  u   res","    it  v   res","   it   w   res",&
     "   it   p   res","    div","       max","      CFLx","       CFLy","       CFLz","     T_moy","    T_inst"

subit:  do subite=1,nav%nsubite
     nav%subite=int(subite,ik)
     if (mpid%rank==0) write(*,'(i5,X,es9.2,X,i2)', advance= 'no') ite,nav%time,subite
     !---------------------------------------------------------------------
!     if (mpid%rank==0) then
!        call color(ired);print'(a)','Time : ';call color(color_off)
!     endif
     !-> define bc
     call navier_bc_velocity(mpid,nav)
     call navier_LES(mpid,nav)
     
     !---------------------------------------------------------------------
     !-> compute rhs
     call navier_presolve_u(mpid,nav)
     call navier_presolve_v(mpid,nav)
     call navier_presolve_w(mpid,nav)

     !---------------------------------------------------------------------
     !-> solve intermediate u,v,w (pressure correction)

     if(nav%pt<=2) then
       if(nav%pt==1) call add_boundary_gradient(mpid,nav)
       call navier_solve_u(mpid,nav)
       call navier_solve_v(mpid,nav)
       call navier_solve_w(mpid,nav)
     endif

     call navier_presolve_phi(mpid,nav)
     call navier_bc_pressure(mpid,nav)
     !---------------------------------------------------------------------
     !-> solve pressure increment phi

     call navier_solve_phi(mpid,nav)

     !---------------------------------------------------------------------
     !-> compute u,v,w,p
     
     call navier_projection(mpid,nav)

     !---------------------------------------------------------------------
     !-> compute final u,v,w  (velocity correction)
     if(nav%pt>=3) then
       call navier_solve_u(mpid,nav)
       call navier_solve_v(mpid,nav)
       call navier_solve_w(mpid,nav)
     endif

     test=.false.
     if (subite>1) call testconv(mpid,nav,nav%u(nav%it(1)),&
                                          nav%v(nav%it(1)),&
                                          nav%w(nav%it(1)),&
                                          nav%p(nav%it(1)),&
                                          nav%sub_u,&
                                          nav%sub_v,&
                                          nav%sub_w,&
                                          nav%sub_p,nav%aux,test,1.d-9)

     nav%sub_u=nav%u(nav%it(1))
     nav%sub_v=nav%v(nav%it(1))
     nav%sub_w=nav%w(nav%it(1))
     nav%sub_p=nav%p(nav%it(1))
 
     nav%aux%f=1._rk
     ref=norme(mpid,nav,nav%aux)

     reft=0.5_rk*(norme(mpid,nav,derx(nav%dcx,nav%u(nav%it(1))))**2 &
                 +norme(mpid,nav,dery(nav%dcy,nav%v(nav%it(1))))**2 &
                 +norme(mpid,nav,derz(nav%dcz,nav%w(nav%it(1)))))

     nav%aux%f=0._rk
     nav%aux=derx(nav%dcx,nav%u(nav%it(1)))+&
             dery(nav%dcy,nav%v(nav%it(1)))+&
             derz(nav%dcz,nav%w(nav%it(1)))

    error=norme(mpid,nav,nav%aux)
    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') error/ref
    if(error>1d10) exit temps

!    nav%aux%f=1._rk ;     reft=integrale(mpid,nav,nav%aux)
!    errort=integrale(mpid,nav,nav%u(nav%it(1)))
!    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort/reft
!    if(errort/reft>1d10) exit temps

    nav%aux=nav%u(nav%it(1))
    call navier_nullify_boundary(mpid,nav,nav%aux,0)
    error= maxval(nav%aux%f)
    call mpi_reduce(error,errort,1,mpi_double_precision,mpi_max,0,&
         mpi_comm_world,mpid%code)
    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort
    if(errort>1d10) exit temps

!    errort=0.5_rk*(norme(mpid,nav,nav%u(nav%it(1)))**2 &
!                  +norme(mpid,nav,nav%v(nav%it(1)))**2 &
!                  +norme(mpid,nav,nav%w(nav%it(1))))
!    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort/reft
!    if(errort/reft>1d10) exit temps
!
!    nav%aux=nav%p(nav%it(1))
!    call navier_nullify_boundary(mpid,nav,nav%aux,0)
!    error=integrale(mpid,nav,nav%aux)
!    nav%aux=nav%aux-error/ref
!    nav%p(nav%it(1))=nav%aux
!
!    nav%aux=nav%p(nav%it(1)) - nav%p(nav%it(nav%nt))
!    call navier_nullify_boundary(mpid,nav,nav%aux,0)
!    error=integrale(mpid,nav,nav%aux)
!    nav%aux=nav%aux- error/reft
!    error=norme(mpid,nav,nav%aux)
!    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') error
!    if(error>1d10) exit temps
!
!    errort=norme(mpid,nav,nav%u(nav%it(1))-nav%u(nav%it(nav%nt)))**2 &
!          +norme(mpid,nav,nav%v(nav%it(1))-nav%v(nav%it(nav%nt)))**2 &
!          +norme(mpid,nav,nav%w(nav%it(1))-nav%w(nav%it(nav%nt)))
!    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') sqrt(errort)
!    if(sqrt(errort)>1d10) exit temps


!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,dx,dy,dz) &
!$OMP SCHEDULE(RUNTIME)
    do k=1,nav%nz
       if(k.eq.1) then
         dz= nav%gridz%grid1d(k+1)- nav%gridz%grid1d(k)
       elseif(k.eq.nav%nz) then
         dz= nav%gridz%grid1d(k  )- nav%gridz%grid1d(k-1)
       else
         dz=(nav%gridz%grid1d(k+1)- nav%gridz%grid1d(k-1))*0.5_rk
       endif
       do j=1,nav%ny
          if(j.eq.1) then
            dy= nav%gridy%grid1d(j+1)- nav%gridy%grid1d(j)
          elseif(j.eq.nav%ny) then
            dy= nav%gridy%grid1d(j  )- nav%gridy%grid1d(j-1)
          else
            dy=(nav%gridy%grid1d(j+1)- nav%gridy%grid1d(j-1))*0.5_rk
          endif
          do i=1,nav%nx
            if(i.eq.1) then
              dx= nav%gridx%grid1d(i+1)- nav%gridx%grid1d(i)
            elseif(i.eq.nav%nx) then
              dx= nav%gridx%grid1d(i  )- nav%gridx%grid1d(i-1)
            else
              dx=(nav%gridx%grid1d(i+1)- nav%gridx%grid1d(i-1))*0.5_rk
            endif
            nav%aux%f(i,j,k) =nav%ts*nav%u(1)%f(i,j,k)/dx
            nav%aux1%f(i,j,k)=nav%ts*nav%v(1)%f(i,j,k)/dy
            nav%aux2%f(i,j,k)=nav%ts*nav%w(1)%f(i,j,k)/dz
          enddo
       enddo
    enddo
!$OMP END PARALLEL DO

    call navier_nullify_boundary(mpid,nav,nav%aux,0)
    error= maxval(nav%aux%f)
    call mpi_reduce(error,errort,1,mpi_double_precision,mpi_max,0,&
         mpi_comm_world,mpid%code)
    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort

    call navier_nullify_boundary(mpid,nav,nav%aux1,0)
    error= maxval(nav%aux1%f)
    call mpi_reduce(error,errort,1,mpi_double_precision,mpi_max,0,&
         mpi_comm_world,mpid%code)
    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort

    call navier_nullify_boundary(mpid,nav,nav%aux2,0)
    error= maxval(nav%aux2%f)
    call mpi_reduce(error,errort,1,mpi_double_precision,mpi_max,0,&
         mpi_comm_world,mpid%code)
    if (mpid%rank==0) write(*,'(X,es9.2)', advance= 'no') errort

  call system_clock(t2,irate)
  time=real(t2-t1)/real(irate)
  if (mpid%rank==0) write(*,'(2(X,es9.2))', advance= 'yes') time/(ite-19),real(t2-t3)/real(irate)

     if (test)     exit subit 
     enddo subit

     !-> switch it
     iaux=nav%it(1)
     do i=1,nav%nt-1
        nav%it(i)=nav%it(i+1)
     enddo
     nav%it(nt)=iaux

  if (int(ite/300)*300.eq.ite) then
    if (mpid%rank==0) print*,'Write fields'
    call restart_write(mpid,nav)
    call navier_write_fields(mpid,nav,ite,ite)
  endif

  enddo temps
  !------------------------------------------------------------------------ 
  !-> time loop end
  !------------------------------------------------------------------------ 
  if (mpid%rank==0) write(*,'(a)', advance= 'yes') ' STOP '
  call system_clock(t2,irate)
  time=real(t2-t1)/real(irate)
  if (mpid%rank==0) print*,'time : ',time,time/(ite-19)
stop

!  call restart_write(mpid,nav)
  call navier_write_fields(mpid,nav,ite,ite)
  !------------------------------------------------------------------------ 
  !-> time loop end
  !------------------------------------------------------------------------ 
  call restart_write(mpid,nav)

  call write_field('vel_u',nav%u(nav%it(nav%nt)),mpid)
  call write_field('vel_v',nav%v(nav%it(nav%nt)),mpid)
  call write_field('vel_w',nav%w(nav%it(nav%nt)),mpid)
  call write_field('vel_p',nav%p(nav%it(nav%nt)),mpid)
  call write_field('vel_phi',nav%phi(nav%it(nav%nt)),mpid)
  call write_field('les_nu',nav%les_nu,mpid)

  !------------------------------------------------------------------------ 
  !-> post-computation
  !------------------------------------------------------------------------ 


  call navier_finalization(cmd,mpid,nav)

contains

subroutine testconv(mpid,nav,u,v,w,p,sub_u,sub_v,sub_w,sub_p,aux,test,eps)
  use class_md
  use class_field
  use class_navier_3D
  use precision
  implicit none
  type(mpi_data) :: mpid
  type(navier3d) :: nav
  type(field) :: aux,u,    v,    w,    p
  type(field) :: sub_u,sub_v,sub_w,sub_p
  logical  :: test
  real(rk) :: ref,error1,error2,eps

!     nav%aux%f=0._rk
!     nav%aux=derx(nav%dcx,nav%u(nav%it(1)))+&
!             dery(nav%dcy,nav%v(nav%it(1)))+&
!             derz(nav%dcz,nav%w(nav%it(1)))
!    ref=norme(mpid,nav,nav%aux)
!    nav%aux=derx(nav%dcx,nav%u(nav%it(1)))+&
!            dery(nav%dcy,nav%v(nav%it(1)))+&
!            derz(nav%dcz,nav%w(nav%it(1)))-&
!            derx(nav%dcx,nav%sub_u)-&
!            dery(nav%dcy,nav%sub_v)-&
!            derz(nav%dcz,nav%sub_w)

!    error=norme(mpid,nav,nav%aux)/ref
!    if (mpid%rank==0) print*,'conv Div V       : ',error1

    aux%f=sqrt(u%f**2 + v%f**2 + w%f**2)
    ref=norme(mpid,nav,aux)

    aux%f=sqrt((sub_u%f-u%f)**2 + (sub_v%f-v%f)**2 + (sub_w%f-w%f)**2)
    error1=norme(mpid,nav,aux)/ref

!    if (mpid%rank==0) print*,'conv tot V       : ',error1

    aux%f=1._rk   ;    ref=integrale(mpid,nav,aux)
    aux=p - sub_p ; error2=integrale(mpid,nav,aux)

    sub_p=sub_p + error2/ref
    aux  =sub_p     ;    ref=norme(mpid,nav,aux)
    aux  =p - sub_p ; error2=norme(mpid,nav,aux)/ref
!    if (mpid%rank==0) print*,'conv tot P       : ',error2

    if (error1<eps.and.error2<eps)  test=.true.

end subroutine testconv

function norme(mpid,nav,x)
  use class_md

  use class_field
  use class_navier_3D
  use precision
  type(mpi_data) :: mpid
  type(navier3d) :: nav
  type(field) :: x
  real(rk) :: norme
norme=norme2(mpid,nav,x)
!norme=normeinf(mpid,nav,x)
end function norme

subroutine initialise_navier(nav,mpid)
  type(navier3d) :: nav
  type(mpi_data) :: mpid
  type(mesh_grid) :: gridxi,gridyi,gridzi
  type(field)     :: uo,vo,wo,po
  character(len=512) :: fich_grid(3),var_grid(3),fich_vel(4),var_vel(4)


  if (mpid%rank==0)     print*,'Initialisation'

if(.false.) then
if(.false.) then

  fich_grid=(/"init/gridxi.nc","init/gridyi.nc","init/gridzi.nc"/)
  var_grid=(/'grid_x','grid_y','grid_z'/)
  fich_vel=(/"init/init2","init/init2","init/init2","init/init2"/)
  var_vel=(/'velocity_x','velocity_y','velocity_z','pressure  '/)

  call read_initfiles(gridxi,gridyi,gridzi,uo,vo,wo,po,fich_grid,var_grid,fich_vel,var_vel)

elseif(.false.) then

  fich_grid=(/"init/grid_x.nc","init/grid_y.nc","init/grid_z.nc"/)
  var_grid=(/'gridx','gridy','gridz'/)
  fich_vel=(/"init/vel_u","init/vel_v","init/vel_w","init/vel_p"/)
!  var_vel=(/'U','V','W','P'/)
  var_vel=(/'U','V','W','P'/)

  call read_initfiles(gridxi,gridyi,gridzi,uo,vo,wo,po,fich_grid,var_grid,fich_vel,var_vel)
po%f=0._rk

endif
  call  interpol_initfiles(nav,gridxi,gridyi,gridzi,uo,vo,wo,po)

  call field_destroy(uo)
  call field_destroy(vo)
  call field_destroy(wo)
  call field_destroy(po)

else

     do iaux=0,nav%nt-1
     t=nav%time-(iaux+1)*nav%ts
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
     do k=1,nav%nz
        z=nav%gridz%grid1d(k)
        do j=1,nav%ny
           y=nav%gridy%grid1d(j)
           do i=1,nav%nx
             x=nav%gridx%grid1d(i)
             nav%u(nav%it(nav%nt-iaux))%f(i,j,k)=sol(x,y,z,t,'u',nav%rey)
             nav%v(nav%it(nav%nt-iaux))%f(i,j,k)=sol(x,y,z,t,'v',nav%rey)
             nav%w(nav%it(nav%nt-iaux))%f(i,j,k)=sol(x,y,z,t,'w',nav%rey)
             nav%p(nav%it(nav%nt-iaux))%f(i,j,k)=sol(x,y,z,t,'p',nav%rey)
!             if(nav%pt==1)
!             if(nav%pt==2)
             if(nav%pt==3) nav%phi(nav%it(nav%nt-iaux))%f(i,j,k)=nav%p(nav%it(nav%nt-iaux))%f(i,j,k)
!             if(nav%pt==4)
           enddo
        enddo
     enddo
!$OMP END PARALLEL DO
     enddo

endif

   do iaux=1,nav%nt 
!periodisation ! todo
!  nav%u(nav%it(iaux))%f(nav%nx,:,:)=nav%u(nav%it(iaux))%f(1,:,:)
!  nav%v(nav%it(iaux))%f(nav%nx,:,:)=nav%v(nav%it(iaux))%f(1,:,:)
!  nav%w(nav%it(iaux))%f(nav%nx,:,:)=nav%w(nav%it(iaux))%f(1,:,:)
!  nav%p(nav%it(iaux))%f(nav%nx,:,:)=nav%p(nav%it(iaux))%f(1,:,:)
!  nav%u(nav%it(iaux))%f(:,:,nav%nz)=nav%u(nav%it(iaux))%f(:,:,1)
!  nav%v(nav%it(iaux))%f(:,:,nav%nz)=nav%v(nav%it(iaux))%f(:,:,1)
!  nav%w(nav%it(iaux))%f(:,:,nav%nz)=nav%w(nav%it(iaux))%f(:,:,1)
!  nav%p(nav%it(iaux))%f(:,:,nav%nz)=nav%p(nav%it(iaux))%f(:,:,1)

  !nettoyage
     call navier_nullify_boundary(mpid,nav,nav%u(nav%it(iaux)),0)
     call navier_nullify_boundary(mpid,nav,nav%v(nav%it(iaux)),0)
     call navier_nullify_boundary(mpid,nav,nav%w(nav%it(iaux)),0)
     call navier_nullify_boundary(mpid,nav,nav%p(nav%it(iaux)),0)
 
     nav%phi(nav%it(iaux))=nav%p(nav%it(iaux))
   enddo

end subroutine initialise_navier

subroutine interpol_initfiles(nav,gridxi,gridyi,gridzi,uo,vo,wo,po)
  type(navier3d),intent(inout) :: nav
  type(mesh_grid),intent(in) :: gridxi,gridyi,gridzi
  type(field),intent(in)     :: uo,vo,wo,po

  integer(ik) :: i,j,k,iaux,i0,j0,k0,i1,j1,k1,i2,j2,k2
  real(rk) :: x,y,z,x1,x2,y1,y2,z1,z2

! simple interpolation for uo -> u and others fields
   do iaux=0,nav%nt-1
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(i,j,k,x,y,z) &
!$OMP SCHEDULE(RUNTIME)
   do i=1,nav%nx
      x=nav%gridx%grid1d(i)
!      x=nav%gridx%grid1d(i)*2._rk
!      x=nav%gridx%grid1d(i)*4._rk
      do i1=2,gridxi%nx
        i2=i1-1
        x1=gridxi%grid1d(i1)
        x2=gridxi%grid1d(i2)
        if (abs(abs(x1-x)+abs(x2-x)-abs(x1-x2))<1e-8) then
         do j=1,nav%ny
            y=nav%gridy%grid1d(j)
!            y=nav%gridy%grid1d(j)+1._rk
            do j1=2,gridyi%ny
              j2=j1-1
              y1=gridyi%grid1d(j1)
              y2=gridyi%grid1d(j2)
               if (abs(abs(y1-y)+abs(y2-y)-abs(y1-y2))<1e-8) then
               do k=1,nav%nz
                  z=nav%gridz%grid1d(k)
!                  z=nav%gridz%grid1d(k)*2._rk
!                  z=nav%gridz%grid1d(k)*4._rk
                  do k1=2,gridzi%nz
                    k2=k1-1
                    z1=gridzi%grid1d(k1)
                    z2=gridzi%grid1d(k2)
                    if (abs(abs(z1-z)+abs(z2-z)-abs(z1-z2))<1e-8) then 

!if(j.ne.1.and.j.ne.nav%ny) then !detection of null edges
!endif



      nav%u(nav%it(nav%nt-iaux))%f(i,j,k)= &
         uo%f(i1,j1,k1)*(x-x2)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i1,j1,k2)*(x-x2)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i1,j2,k1)*(x-x2)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i1,j2,k2)*(x-x2)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i2,j1,k1)*(x1-x)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i2,j1,k2)*(x1-x)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i2,j2,k1)*(x1-x)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + uo%f(i2,j2,k2)*(x1-x)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2))

      nav%v(nav%it(nav%nt-iaux))%f(i,j,k)= &
         vo%f(i1,j1,k1)*(x-x2)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i1,j1,k2)*(x-x2)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i1,j2,k1)*(x-x2)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i1,j2,k2)*(x-x2)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i2,j1,k1)*(x1-x)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i2,j1,k2)*(x1-x)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i2,j2,k1)*(x1-x)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + vo%f(i2,j2,k2)*(x1-x)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2))

      nav%w(nav%it(nav%nt-iaux))%f(i,j,k)= &
         wo%f(i1,j1,k1)*(x-x2)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i1,j1,k2)*(x-x2)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i1,j2,k1)*(x-x2)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i1,j2,k2)*(x-x2)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i2,j1,k1)*(x1-x)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i2,j1,k2)*(x1-x)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i2,j2,k1)*(x1-x)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + wo%f(i2,j2,k2)*(x1-x)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2))

      nav%p(nav%it(nav%nt-iaux))%f(i,j,k)= &
         po%f(i1,j1,k1)*(x-x2)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i1,j1,k2)*(x-x2)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i1,j2,k1)*(x-x2)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i1,j2,k2)*(x-x2)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i2,j1,k1)*(x1-x)*(y-y2)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i2,j1,k2)*(x1-x)*(y-y2)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i2,j2,k1)*(x1-x)*(y1-y)*(z-z2)/((x1-x2)*(y1-y2)*(z1-z2)) &
       + po%f(i2,j2,k2)*(x1-x)*(y1-y)*(z1-z)/((x1-x2)*(y1-y2)*(z1-z2))

                        exit
                      endif
                    enddo
                 enddo

                  exit
                endif
              enddo
           enddo

          exit
        endif
      enddo
   enddo
!$OMP END PARALLEL DO
   enddo

end subroutine interpol_initfiles

subroutine read_initfiles(gridxi,gridyi,gridzi,uo,vo,wo,po,fich_grid,var_grid,fich_vel,var_vel)
  type(mesh_grid),intent(out) :: gridxi,gridyi,gridzi
  type(field),intent(out)     :: uo,vo,wo,po
  character(len=512),intent(in) :: fich_grid(3),var_grid(3),fich_vel(4),var_vel(4)
  integer(ik) :: i
    integer(ik) :: dim_len(3)
    character(len=512) :: dim_name(3)
    integer(ik) :: ncid
    integer(ik) :: startv(3),countv(3)
    integer(ik) :: varid(1)

do i=1,3 ! read grid
    call get_var3d_info(trim(fich_grid(i)),trim(var_grid(i)),dim_name(1),dim_len(1))

    startv=1
    countv=dim_len

    if(dim_len(i).eq.1) then
       dim_len(i)=dim_len(1)
      do j=1,3
        if(j/=i) dim_len(j)=1
      enddo
    else
      do j=1,3
        if(j/=i) dim_len(j)=1
        if(j/=i) countv(j)=1
      enddo
    endif

    if(i==1) call mesh_init(gridxi,'gridx','x',dim_len(1),dim_len(2),dim_len(3))
    if(i==2) call mesh_init(gridyi,'gridy','y',dim_len(1),dim_len(2),dim_len(3))
    if(i==3) call mesh_init(gridzi,'gridz','z',dim_len(1),dim_len(2),dim_len(3))

    call io_check(nf90_open(path=trim(fich_grid(i)),mode=nf90_nowrite,ncid=ncid))

                            
    !-> get variable id
    call io_check(nf90_inq_varid(ncid,trim(var_grid(i)),varid(1)))
                                                        
    !-> read field variable
    if(i==1) call io_check(nf90_get_var(ncid,varid(1),gridxi%grid3d,start=startv,count=countv))
    if(i==2) call io_check(nf90_get_var(ncid,varid(1),gridyi%grid3d,start=startv,count=countv))
    if(i==3) call io_check(nf90_get_var(ncid,varid(1),gridzi%grid3d,start=startv,count=countv))

    !-> close file              
    call io_check(nf90_close(ncid))

    if(i==1) gridxi%grid1d(:)=gridxi%grid3d(:,1,1)
    if(i==2) gridyi%grid1d(:)=gridyi%grid3d(1,:,1)
    if(i==3) gridzi%grid1d(:)=gridzi%grid3d(1,1,:)

enddo

       call field_init(uo,"U",gridxi%n,gridyi%n,gridzi%n)
       call field_init(vo,"V",gridxi%n,gridyi%n,gridzi%n)
       call field_init(wo,"W",gridxi%n,gridyi%n,gridzi%n)
       call field_init(po,"P",gridxi%n,gridyi%n,gridzi%n)

  call read_field(trim(fich_vel(1)),uo,trim(var_vel(1)))  
  call read_field(trim(fich_vel(2)),vo,trim(var_vel(2)))  
  call read_field(trim(fich_vel(3)),wo,trim(var_vel(3)))  
  call read_field(trim(fich_vel(4)),po,trim(var_vel(4)))  

end subroutine read_initfiles

end program testnavier3d
