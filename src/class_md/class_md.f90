module class_md
! -----------------------------------------------------------------------
! Name :
! class md
! -----------------------------------------------------------------------
! Object : 
! Solve poisson and helmholtz equation with a domain decomposition
! -----------------------------------------------------------------------
! Files :
! class_md.f90
! -----------------------------------------------------------------------
! Public type :
! mpid_data -> mpi parameters
! -----------------------------------------------------------------------
! Public Procedure :
! 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  use precision
!  use mpi
  use command_line
#include <finclude/petscdef.h>
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscmat.h>
#include <finclude/petscvec.h>
#include <finclude/petscviewer.h>
#include <finclude/petscvec.h90>

  type,public :: mpi_data
     !-> mpi communicators
     integer :: code,comm,comm_3d
     !-> mpi number of processus and number of dimensions
     integer :: processus,dims
     !-> mpi rank and cartesian rank
     integer :: rank,rank_3d
     !-> number of points in x,y,z directions
     integer(ik) :: nx,ny,nz
     !-> number of domains in x,y,z directions
     integer(ik) :: nd(3)
     !-> total number of domains 
     integer(ik) :: ndom
     !-> periodicity in x,y,z directions
     logical :: periods(3)
     !-> coordinates of domain
     integer(ik) :: coord(3)
     !-> rank of neighbours (directions and left/right)
     integer(ik) :: neighbours(3,2)
     !-> rank of all process in cartesian topology
     integer(ik),allocatable :: domains(:,:,:)
  end type mpi_data

  type,public :: mpi_inf_mat
     !-> number of points in x,y,z directions
     integer(ik) :: nx,ny,nz
     !-> number of domains in x,y,z directions
     integer(ik) :: nd(3)
     !-> number of interfaces : x,y,z,total (size of influence matrix)
     integer(ik) :: ix,iy,iz,it
     !-> periodicity in x,y,z directions
     logical :: periods(3)
     !-> periodicity in x,y,z directions
     integer(ik) :: vertex(2,2,2)
     !-> interfaces rank for all domains x,y,z : 
     !   3 directions * left,right faces 
     integer,allocatable :: dm_int(:,:,:,:,:)
     !-> influence matrix dimensions
     !   (i,1) : direction (x,y,z) ; (i,2) : dimensions
     integer,allocatable :: dm_int_size(:,:)
     !-> influence matrix rows domains id
     integer,allocatable :: inf_rows_id(:,:,:)
     !-> petsc error code
     PetscErrorCode :: err
     !-> size of influence matrix
     integer :: ninf
     !-> number of rows, of non-zero diagonal and off-diagonal elements
     integer :: rows,d_nz
     integer,allocatable :: o_nnz(:),rows_size(:)
     !-> petsc influence matrix
     Mat :: inf
     !-> petsc vectors : rhs and solution
     Vec :: rhs,sol,sol_old(3)
     !-> petsc solver context
     KSP :: ksp
     character(2) :: kspn
     !-> nullspace
     MatNullSpace :: nullsp
  end type mpi_inf_mat

  type,public :: mpi_inf_sol
     !-> petsc vectors : rhs and solution
     Vec :: sol_old(3)
  end type mpi_inf_sol


  integer(ik),save,public :: md_lenght_char=100

  !-> Declare everything private by default
  private

  !-> Declare exported procedure
  public :: md_mpi_init,md_mpi_finalize,md_mpi_barrier
  public :: md_petsc_initialize,md_petsc_finalize
  public :: md_mpi_getcoord,md_get_interfaces_number
  public :: md_mpi_getnumberdomains
  public :: md_influence_matrix_init_start,md_influence_matrix_init_end
  public :: md_influence_matrix_destroy
  public :: md_influence_matrix_setvalues
  public :: md_influence_matrix_view,md_influence_matrix_infos
  public :: md_check
  public :: md_influence_matrix_start_flush,md_influence_matrix_end_flush
  public :: md_influence_matrix_write,md_influence_matrix_read
  public :: md_influence_matrix_filename
  public :: md_mpi_reduce_double,md_mpi_bcast_double
  public :: md_vector_setvalues,md_vector_getvalues
  public :: md_rhs_view,md_sol_view,md_solve,md_solve_init
  public :: md_bc_init,md_solve_guess_nonzero,md_rhs_nullspace
  public :: md_mpi_global_coord,md_vector_sol_setvalues
  public :: md_guess_init
  public :: md_add_pert,md_vector_zero_lastpoint
  public :: md_influence_guess_write,md_influence_guess_read
contains

! =======================================================================
! =======================================================================
! md : influence matrix
! =======================================================================
! =======================================================================
!------------------------------------------------------------------------
! md : check petsc solver
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_check(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    Vec :: solx
    real(rk) :: neg_one,norm,one,t(6),time(3),timet(3)
    integer(ik) :: iters(3)
    PetscLogStage :: stage

    neg_one = -1.0

    call VecDuplicate(inf_mat%sol,solx,inf_mat%err)

    one=1.0
    call VecSet(solx,one,inf_mat%err)
!    call VecSetRandom(solx,PETSC_NULL,inf_mat%err)
    call MatMult(inf_mat%inf,solx,inf_mat%rhs,inf_mat%err)

    !call VecView(inf_mat%rhs,PETSC_VIEWER_STDOUT_WORLD,inf_mat%err)

    !-> Create linear solver context
    call KSPCreate(PETSC_COMM_WORLD,inf_mat%ksp,inf_mat%err)
    
    !->  Set operators
    call KSPSetOperators(inf_mat%ksp,inf_mat%inf,inf_mat%inf,&
         SAME_PRECONDITIONER,inf_mat%err)
    
    !-> get solver type from command line
    call KSPSetFromOptions(inf_mat%ksp,inf_mat%err)


    !-> solve system
    t(1)=MPI_WTIME()
    call KSPSolve(inf_mat%ksp,inf_mat%rhs,inf_mat%sol,inf_mat%err)
    t(2)=MPI_WTIME()
    call KSPGetIterationNumber(inf_mat%ksp,iters(1),inf_mat%err)


    !-> solve system
    t(3)=MPI_WTIME()
!    call PetscLogStageRegister("Solve",stage,inf_mat%err)
!    call PetscLogStagePush(stage,inf_mat%err)
    call KSPSolve(inf_mat%ksp,inf_mat%rhs,inf_mat%sol,inf_mat%err)
!    call PetscLogStagePop(inf_mat%err)
    t(4)=MPI_WTIME()
    call KSPGetIterationNumber(inf_mat%ksp,iters(2),inf_mat%err)


    one=1.0001
    call VecSet(solx,one,inf_mat%err)
    call MatMult(inf_mat%inf,solx,inf_mat%rhs,inf_mat%err)

    !-> solve system
    call KSPSetInitialGuessNonzero(inf_mat%ksp,PETSC_TRUE,inf_mat%err)
    t(5)=MPI_WTIME()
    call KSPSolve(inf_mat%ksp,inf_mat%rhs,inf_mat%sol,inf_mat%err)
    t(6)=MPI_WTIME()
    call KSPGetIterationNumber(inf_mat%ksp,iters(3),inf_mat%err)

    !-> compare sol and solx
    call VecAXPY(solx,neg_one,inf_mat%sol,inf_mat%err)
    call VecNorm(solx,NORM_2,norm,inf_mat%err)

    !-> time output
    time(1)=t(2)-t(1) ; time(2)=t(4)-t(3) ; time(3)=t(6)-t(5)
    call md_mpi_reduce_double(mpid,time(1),timet(1))
    call md_mpi_reduce_double(mpid,time(2),timet(2))
    call md_mpi_reduce_double(mpid,time(3),timet(3))
    timet(:)=timet(:)/mpid%processus
    if (mpid%rank==0) print*,'Solution',norm,iters,timet

    call VecDestroy(solx,inf_mat%err)


    !call VecView(inf_mat%sol,PETSC_VIEWER_STDOUT_WORLD,inf_mat%err)
  end subroutine md_check

!------------------------------------------------------------------------
! md : fix one dirichlet with full neumann bc
!------------------------------------------------------------------------
! Matthieu Marquillie
! 01/2013
!
  subroutine md_add_pert(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    real(rk),allocatable :: val(:)

    if (mpid%rank==0) then
       allocate(val(inf_mat%rows+inf_mat%o_nnz(1)))
       val=0._rk
       val(1)=5._rk
!!       val=-1._rk
       call MatSetValuesRow(inf_mat%inf,0,val,inf_mat%err)
    endif

!    if (mpid%rank==mpid%processus-2) then
!       allocate(val(inf_mat%rows+inf_mat%o_nnz(inf_mat%rows)))
!       val=1._rk
!       call MatSetValuesRow(inf_mat%inf,inf_mat%ninf-1,val,inf_mat%err)
!    endif

  end subroutine md_add_pert

!------------------------------------------------------------------------
! md : fix one dirichlet with full neumann bc
!------------------------------------------------------------------------
! Matthieu Marquillie
! 01/2013
!
  subroutine md_vector_zero_lastpoint(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat

    if (mpid%rank==0) &
         call VecSetValue(inf_mat%rhs,0,0.d0,&
         INSERT_VALUES,inf_mat%err)

!    if (mpid%rank==mpid%processus-2) &
!         call VecSetValue(inf_mat%rhs,inf_mat%rows-1,0.d0,&
!         INSERT_VALUES,inf_mat%err)

    call VecAssemblyBegin(inf_mat%rhs,inf_mat%err)
    call VecAssemblyEnd(inf_mat%rhs,inf_mat%err)

  end subroutine md_vector_zero_lastpoint

!------------------------------------------------------------------------
! md : put dirichlet bc on interface
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_bc_init(mpid,inf_mat,bctype)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: bctype(2,3),l,m,c(3),inter(3,2)
    
    !-> compute indexes of where to put values
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(inf_mat,c,inter)

    do m=1,2
       do l=1,3
          if (inter(l,m)>0) then
             bctype(m,l)=1
          endif
       enddo
    enddo

  end subroutine md_bc_init
!------------------------------------------------------------------------
! md : petsc solver init
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_solve_init(mpid,inf_mat,null,kspn)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: iters
    integer(ik),optional :: null
    character(*),optional :: kspn

    !-> Create linear solver context
    call KSPCreate(PETSC_COMM_WORLD,inf_mat%ksp,inf_mat%err)
    
    !-> set option prefix name
    if (present(kspn)) then
       inf_mat%kspn=kspn
       call KSPSetOptionsPrefix(inf_mat%ksp,inf_mat%kspn,inf_mat%err)
    endif

    !->  Set operators
    call KSPSetOperators(inf_mat%ksp,inf_mat%inf,inf_mat%inf,&
         SAME_PRECONDITIONER,inf_mat%err)
    
    !-> set nullspace for singular problem
    if (present(null)) then
       if (null==1) then
          call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,&
               PETSC_NULL,inf_mat%nullsp,inf_mat%err)
          call KSPSetNullSpace(inf_mat%ksp, inf_mat%nullsp,inf_mat%err)
          call MatNullSpaceDestroy(inf_mat%nullsp,inf_mat%err)
       endif
    endif

    !-> get solver type from command line
    call KSPSetFromOptions(inf_mat%ksp,inf_mat%err)

    !-> use previous value for initial guess
    call KSPSetInitialGuessNonzero(inf_mat%ksp,PETSC_TRUE,inf_mat%err)

    !-> setup ksp
    call KSPSetUp(inf_mat%ksp,inf_mat%err)

  end subroutine md_solve_init
!------------------------------------------------------------------------
! md : petsc rhs set nullspace
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_rhs_nullspace(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat

    call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,0,inf_mat%nullsp,&
         inf_mat%err)
    call MatNullSpaceRemove(inf_mat%nullsp,inf_mat%rhs,PETSC_NULL,inf_mat%err)
    call MatNullSpaceDestroy(inf_mat%nullsp,inf_mat%err);

  end subroutine md_rhs_nullspace
!------------------------------------------------------------------------
! md : petsc solver set initial gess non zero
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_solve_guess_nonzero(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat

    !-> use previous value for initial guess
    call KSPSetInitialGuessNonzero(inf_mat%ksp,PETSC_TRUE,inf_mat%err)

  end subroutine md_solve_guess_nonzero
!------------------------------------------------------------------------
! md : petsc solver
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_solve(mpid,inf_mat,inf_sol)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    type(mpi_inf_sol),optional :: inf_sol
    integer(ik) :: iters
    Vec :: solx
    real(rk) :: neg_one,norm,fac(3)

    !--------------------------------------------------------------------
    !-> initial guess
    if (present(inf_sol)) then
       !-> first order
       call VecCopy(inf_sol%sol_old(3),inf_mat%sol,inf_mat%err)

       !-> 2nd order
       !fac(1)=-1._rk ; fac(2)=2._rk
       !call VecAXPBY(inf_mat%sol,fac(1),fac(2),inf_sol%sol_old(2),inf_mat%err)

       !-> third order
!       fac(1)=1._rk ; fac(2)=-3._rk ; fac(3)=3._rk
!       call VecAXPBYPCZ(inf_mat%sol,fac(1),fac(2),fac(3),&
!            inf_sol%sol_old(1),inf_sol%sol_old(2),inf_mat%err)
    endif

    !--------------------------------------------------------------------
    !-> solve system
    call KSPSolve(inf_mat%ksp,inf_mat%rhs,inf_mat%sol,inf_mat%err)

    !-> get numbers of iterations
    call KSPGetIterationNumber(inf_mat%ksp,iters,inf_mat%err)
!    if (mpid%rank==0) print*,'Iterations',iters

    !-> get residual norm for solution
    call KSPGetResidualNorm(inf_mat%ksp,norm,inf_mat%err)

    !--------------------------------------------------------------------
    !-> put sol in sol_old for initial guess
    if (present(inf_sol)) then
       call VecCopy(inf_sol%sol_old(2),inf_sol%sol_old(1),inf_mat%err)
       call VecCopy(inf_sol%sol_old(3),inf_sol%sol_old(2),inf_mat%err)
       call VecCopy(inf_mat%sol,inf_sol%sol_old(3),inf_mat%err)
    endif

    !-> compare sol and solx
!    call VecDuplicate(inf_mat%sol,solx,inf_mat%err)
!    call MatMult(inf_mat%inf,inf_mat%sol,solx,inf_mat%err)

!    neg_one = -1.0
!    call VecAXPY(solx,neg_one,inf_mat%rhs,inf_mat%err)
!    call VecNorm(solx,NORM_2,norm,inf_mat%err)
!    if (mpid%rank==0) print*,'Solution',norm

!    call VecDestroy(solx,inf_mat%err)

    if (mpid%rank==0) print'(a,i0,a,es17.8)',' Iterations ',iters,'  Residus solver',norm


  end subroutine md_solve
!------------------------------------------------------------------------
! md : get values in rhs vector
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_vector_getvalues(mpid,inf_mat,bcx,bcy,bcz,nx,ny,nz)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,l,m,inter(3,2),nrows,c(3)
    integer(ik) :: nx,ny,nz,c1,c2,sum
    !real(rk),intent(in) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
    real(rk),intent(inout) :: bcx(ny*nz,2),bcy(nx*nz,2),bcz(nx*ny,2)
    real(rk) :: bcxt(ny*nz),bcyt(nx*nz),bczt(nx*ny)
    integer(ik) :: idrx(ny*nz),idry(nx*nz),idrz(nx*ny)
    integer(ik) :: tag,req(6),error,it
    PetscScalar, pointer :: vec_point(:)
    integer :: status(MPI_STATUS_SIZE)
    real(rk) :: zero
    
    !bcx=0._rk ; bcy=0._rk ; bcz=0._rk

    !-> compute indexes of where to put values
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(inf_mat,c,inter)

    !-> get local vectors
    call VecGetArrayF90(inf_mat%sol,vec_point,inf_mat%err)
    !print*,mpid%rank,size(vec_point),lbound(vec_point),ubound(vec_point)


    !-> put local vectors in boundary condition array
    sum=1
    do l=1,3
       if (inter(l,2)>0) then

          !-> compute local rows locations
          if (l==1) nrows=inf_mat%ny*inf_mat%nz
          if (l==2) nrows=inf_mat%nx*inf_mat%nz
          if (l==3) nrows=inf_mat%nx*inf_mat%ny
          c1=sum
          c2=c1+nrows-1
          sum=sum+nrows

          !-> get local boundary conditions
          if (l==1) bcx(:,2)=vec_point(c1:c2)
          if (l==2) bcy(:,2)=vec_point(c1:c2)
          if (l==3) bcz(:,2)=vec_point(c1:c2)
          
       endif
    enddo

    !-> begin tranfer boundary condition to neighbours
    req=(/1,2,3,4,5,6/)
    it=0
    do m=1,2
       do l=1,3
          if (inter(l,m)>0) then
             it=it+1
             
             !-> compute local rows locations
             if (l==1) nrows=inf_mat%ny*inf_mat%nz
             if (l==2) nrows=inf_mat%nx*inf_mat%nz
             if (l==3) nrows=inf_mat%nx*inf_mat%ny

             if (l==1.and.m==1) call mpi_irecv(bcxt,nrows,mpi_double_precision, &
                  mpid%neighbours(l,1),1,mpi_comm_world,req(1),error)
             if (l==2.and.m==1) call mpi_irecv(bcyt,nrows,mpi_double_precision, &
                  mpid%neighbours(l,1),2,mpi_comm_world,req(2),error)
             if (l==3.and.m==1) call mpi_irecv(bczt,nrows,mpi_double_precision, &
                  mpid%neighbours(l,1),3,mpi_comm_world,req(3),error)

             if (l==1.and.m==2) call mpi_issend(bcx(1,2),nrows,mpi_double_precision, &
                  mpid%neighbours(l,2),1,mpi_comm_world,req(4),error)
             if (l==2.and.m==2) call mpi_issend(bcy(1,2),nrows,mpi_double_precision, &
                  mpid%neighbours(l,2),2,mpi_comm_world,req(5),error)
             if (l==3.and.m==2) call mpi_issend(bcz(1,2),nrows,mpi_double_precision, &
                  mpid%neighbours(l,2),3,mpi_comm_world,req(6),error)

          endif
       enddo
    enddo

    !-> end tranfer boundary condition to neighbours
    do m=1,2
       do l=1,3
          if (inter(l,m)>0) then
             if (l==1.and.m==2) call mpi_wait(req(4),status,error)
             if (l==2.and.m==2) call mpi_wait(req(5),status,error)
             if (l==3.and.m==2) call mpi_wait(req(6),status,error)

             if (l==1.and.m==1) call mpi_wait(req(1),status,error)
             if (l==2.and.m==1) call mpi_wait(req(2),status,error)
             if (l==3.and.m==1) call mpi_wait(req(3),status,error)
             
             if (l==1.and.m==1) bcx(:,1)=bcxt(:)
             if (l==2.and.m==1) bcy(:,1)=bcyt(:)
             if (l==3.and.m==1) bcz(:,1)=bczt(:)
          endif
       enddo
    enddo

    call VecRestoreArrayF90(inf_mat%sol,vec_point,inf_mat%err)

    zero=0._rk
    call VecSet(inf_mat%rhs,zero,inf_mat%err)

!    zero=0._rk
!    call VecSet(inf_mat%sol,zero,inf_mat%err)

  end subroutine md_vector_getvalues
!------------------------------------------------------------------------
! md : set values in rhs vector
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_vector_setvalues(mpid,inf_mat,bcx,bcy,bcz,nx,ny,nz)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,l,m,intery(3,2),nrows,c(3)
    integer(ik) :: nx,ny,nz,c1,c2
    real(rk),intent(in) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
    integer(ik) :: idrx(ny*nz),idry(nx*nz),idrz(nx*ny)

    !-> compute indexes of where to put values
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(inf_mat,c,intery)

    do m=1,2
       do l=1,3

          !-> only process effective interface
          if (intery(l,m)<0) cycle

          !-> set number of rows depending of direction of interface
          if (l==1) nrows=inf_mat%ny*inf_mat%nz
          if (l==2) nrows=inf_mat%nx*inf_mat%nz
          if (l==3) nrows=inf_mat%nx*inf_mat%ny

          !-> compute rows locations
          c1=inf_mat%dm_int_size(intery(l,m),3)-1
          c2=c1+nrows-1
          if (l==1) idrx(:)=(/ (i,i=c1,c2) /)
          if (l==2) idry(:)=(/ (i,i=c1,c2) /)
          if (l==3) idrz(:)=(/ (i,i=c1,c2) /)

          !-> put values in rhs
          if (l==1) call VecSetValues(inf_mat%rhs,nrows,idrx,bcx(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==2) call VecSetValues(inf_mat%rhs,nrows,idry,bcy(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==3) call VecSetValues(inf_mat%rhs,nrows,idrz,bcz(1,1,m),&
               ADD_VALUES,inf_mat%err)

       enddo
    enddo

    call VecAssemblyBegin(inf_mat%rhs,inf_mat%err)
    call VecAssemblyEnd(inf_mat%rhs,inf_mat%err)

end subroutine md_vector_setvalues

!------------------------------------------------------------------------
! md : set values in sol vector
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_vector_sol_setvalues(mpid,inf_mat,bcx,bcy,bcz,nx,ny,nz)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,l,m,intery(3,2),nrows,c(3)
    integer(ik) :: nx,ny,nz,c1,c2
    real(rk),intent(in) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
    integer(ik) :: idrx(ny*nz),idry(nx*nz),idrz(nx*ny)

    !-> compute indexes of where to put values
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(inf_mat,c,intery)

    do m=1,2
!    do m=2,2
       do l=1,3

          !-> only process effective interface
          if (intery(l,m)<0) cycle

          !-> set number of rows depending of direction of interface
          if (l==1) nrows=inf_mat%ny*inf_mat%nz
          if (l==2) nrows=inf_mat%nx*inf_mat%nz
          if (l==3) nrows=inf_mat%nx*inf_mat%ny

          !-> compute rows locations
          c1=inf_mat%dm_int_size(intery(l,m),3)-1
          c2=c1+nrows-1
          if (l==1) idrx(:)=(/ (i,i=c1,c2) /)
          if (l==2) idry(:)=(/ (i,i=c1,c2) /)
          if (l==3) idrz(:)=(/ (i,i=c1,c2) /)

          !-> put values in sol
          if (l==1) call VecSetValues(inf_mat%sol,nrows,idrx,bcx(1,1,m),&
               INSERT_VALUES,inf_mat%err)
          if (l==2) call VecSetValues(inf_mat%sol,nrows,idry,bcy(1,1,m),&
               INSERT_VALUES,inf_mat%err)
          if (l==3) call VecSetValues(inf_mat%sol,nrows,idrz,bcz(1,1,m),&
               INSERT_VALUES,inf_mat%err)

       enddo
    enddo

    call VecAssemblyBegin(inf_mat%sol,inf_mat%err)
    call VecAssemblyEnd(inf_mat%sol,inf_mat%err)

  end subroutine md_vector_sol_setvalues

!------------------------------------------------------------------------
! md : rhs vector view
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_rhs_view(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    
    call VecView(inf_mat%rhs,PETSC_VIEWER_STDOUT_WORLD,inf_mat%err)

  end subroutine md_rhs_view
!------------------------------------------------------------------------
! md : rhs vector view
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_sol_view(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    
    call VecView(inf_mat%sol,PETSC_VIEWER_STDOUT_WORLD,inf_mat%err)

  end subroutine md_sol_view
!------------------------------------------------------------------------
! md : PETSC initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_petsc_initialize()
    implicit none
    integer :: err

    !-> initialize petsc
    call PetscInitialize(PETSC_NULL_CHARACTER,err)

  end subroutine md_petsc_initialize
!------------------------------------------------------------------------
! md : PETSC finalization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_petsc_finalize()
    implicit none
    integer :: err

    !-> finalize petsc
    call petscfinalize(PETSC_NULL_CHARACTER,err)    

  end subroutine md_petsc_finalize
!------------------------------------------------------------------------
! md : influence matrix start initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_init_start(mpid,inf_mat,alloc)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    logical :: alloc

    !-> compute parameters of influence matrix
    call md_influence_matrix_param(mpid,inf_mat)

    !-> compute vertex informations (memory problem but not needed)
!    call md_vertex_init(mpid,inf_mat)

    !-> petsc matrix initialization 
    call md_petsc_influence_matrix_init(mpid,inf_mat,alloc)

  end subroutine md_influence_matrix_init_start
!------------------------------------------------------------------------
! md : influence matrix end initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_init_end(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: first,last,rows,columns
    real(8) :: info(MAT_INFO_SIZE)

    call MatAssemblyBegin(inf_mat%inf,MAT_FINAL_ASSEMBLY,inf_mat%err)
    call MatAssemblyEnd(inf_mat%inf,MAT_FINAL_ASSEMBLY,inf_mat%err)

  end subroutine md_influence_matrix_init_end
!------------------------------------------------------------------------
! md : influence matrix infos
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_infos(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: first,last,rows,columns
    real(8) :: info(MAT_INFO_SIZE),fullmem,scalemem
    PetscBool :: flag
    scalemem=1024._rk*1024._rk

!    call MatGetInfo(inf_mat%inf,MAT_LOCAL,info,inf_mat%err)
!    print'(a,x,a,i0,2x,10(a,f0.2,x))','Matinfo local : ',&
!         'rank  ',mpid%rank,' , mallocs  ',info(MAT_INFO_MALLOCS),&
!         ' , Nonzero allocated  ',info(MAT_INFO_NZ_ALLOCATED),&
!         ' , Nonzero used  ',info(MAT_INFO_NZ_USED),&
!         ' , Nonzero unneeded  ',info(MAT_INFO_NZ_UNNEEDED),&
!         ' , Memory used  ',info(MAT_INFO_MEMORY)/scalemem

    call MatGetInfo(inf_mat%inf,MAT_GLOBAL_SUM,info,inf_mat%err)
    fullmem=real(inf_mat%ninf,8)*real(inf_mat%ninf,8)*8._rk
    if (mpid%rank==0) then
       print'(a,x,a,i0,2x,10(a,f0.2,x))','Matinfo global : ',&
            'rank  ',mpid%rank,' , mallocs  ',info(MAT_INFO_MALLOCS),&
            ' , Nonzero allocated  ',info(MAT_INFO_NZ_ALLOCATED),&
            ' , Nonzero used  ',info(MAT_INFO_NZ_USED),&
            ' , Nonzero unneeded  ',info(MAT_INFO_NZ_UNNEEDED),&
            ' , Memory used  ',info(MAT_INFO_MEMORY)/scalemem,&
            ' , Ratio mem  ',info(MAT_INFO_MEMORY)/fullmem,&
            ' , Full inf mem  ',fullmem/scalemem
    endif

  end subroutine md_influence_matrix_infos

!------------------------------------------------------------------------
! md : influence matrix start flush
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_start_flush(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    
    call MatAssemblyBegin(inf_mat%inf,MAT_FLUSH_ASSEMBLY,inf_mat%err)

  end subroutine md_influence_matrix_start_flush
!------------------------------------------------------------------------
! md : influence matrix start flush
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_end_flush(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    
    call MatAssemblyEnd(inf_mat%inf,MAT_FLUSH_ASSEMBLY,inf_mat%err)

  end subroutine md_influence_matrix_end_flush

!------------------------------------------------------------------------
! md : influence matrix write
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_filename(mpid,filename,var)
    implicit none
    type(mpi_data) :: mpid
    character(md_lenght_char) :: filename
    character(*) :: var
    
    !-> create name of file from dimensions and number of domains
    write(filename,'(a,a,a,6(i0,a))')'matrix_',var,'_',mpid%nx, &
         '.',mpid%ny,'.',mpid%nz,'_', &
         mpid%nd(1),'.',mpid%nd(2),'.',mpid%nd(3),'.dat'

  end subroutine md_influence_matrix_filename
!------------------------------------------------------------------------
! md : influence matrix write
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_write(mpid,inf_mat,filename)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    PetscViewer :: viewer
    character(*) :: filename

    !-> open viewer
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename, &
         FILE_MODE_WRITE,viewer,inf_mat%err)
    !-> write matrix
    call MatView(inf_mat%inf,viewer,inf_mat%err)
    !-> destroy viewer
    call PetscViewerDestroy(viewer,inf_mat%err)

  end subroutine md_influence_matrix_write

!------------------------------------------------------------------------
! md : influence matrix read
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_read(mpid,inf_mat,filename)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    PetscViewer :: viewer
    character(*) :: filename

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename, &
         FILE_MODE_READ,viewer,inf_mat%err)
    call MatLoad(inf_mat%inf,viewer,inf_mat%err)
    call PetscViewerDestroy(viewer,inf_mat%err)

  end subroutine md_influence_matrix_read

!------------------------------------------------------------------------
! md : guess vector write
!------------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
  subroutine md_influence_guess_write(mpid,inf_sol,filename)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_sol) :: inf_sol
    PetscViewer :: viewer
    PetscErrorCode :: err
    character(*) :: filename

    !-> open viewer
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename, &
         FILE_MODE_WRITE,viewer,err)
    !-> write vector
    call VecView(inf_sol%sol_old(1),viewer,err)
    !-> write vector
    call VecView(inf_sol%sol_old(2),viewer,err)
    !-> write vector
    call VecView(inf_sol%sol_old(3),viewer,err)
    !-> destroy viewer
    call PetscViewerDestroy(viewer,err)

  end subroutine md_influence_guess_write

!------------------------------------------------------------------------
! md : guess vector write
!------------------------------------------------------------------------
! Matthieu Marquillie
! 04/2013
!
  subroutine md_influence_guess_read(mpid,inf_sol,filename)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_sol) :: inf_sol
    PetscViewer :: viewer
    PetscErrorCode :: err
    character(*) :: filename

    !-> open viewer
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename, &
         FILE_MODE_READ,viewer,err)
    !-> read vector
    call VecLoad(inf_sol%sol_old(1),viewer,err)
    !-> read vector
    call VecLoad(inf_sol%sol_old(2),viewer,err)
    !-> read vector
    call VecLoad(inf_sol%sol_old(3),viewer,err)
    !-> destroy viewer
    call PetscViewerDestroy(viewer,err)

  end subroutine md_influence_guess_read

!------------------------------------------------------------------------
! md : influence matrix initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_setvalues(mpid,inf_mat,c,interx,&
       bcx,bcy,bcz,nx,ny,nz,locr)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,l,m,interx,intery(3,2),ncolumns,nrows,c(3),locr
    integer(ik) :: nx,ny,nz,c1,c2
    real(rk),intent(in) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
    integer(ik) :: idrx(ny*nz),idry(nx*nz),idrz(nx*ny),idc

    !-> compute indexes of where to put values
    call md_get_interfaces_number(inf_mat,c,intery)

    do m=1,2
       do l=1,3
          !-> only process effective interface
          if (intery(l,m)<0) cycle
          !if (intery(l,m)/=interx) cycle
          
          !-> set number of rows depending of direction of interface
          if (l==1) nrows=inf_mat%ny*inf_mat%nz
          if (l==2) nrows=inf_mat%nx*inf_mat%nz
          if (l==3) nrows=inf_mat%nx*inf_mat%ny

          !-> compute column location
          idc=inf_mat%dm_int_size(interx,3)+locr-1

          !-> compute rows locations
          c1=inf_mat%dm_int_size(intery(l,m),3)-1
          c2=c1+nrows-1
          if (l==1) idrx(:)=(/ (i,i=c1,c2) /)
          if (l==2) idry(:)=(/ (i,i=c1,c2) /)
          if (l==3) idrz(:)=(/ (i,i=c1,c2) /)
          
          !-> put values in influence matrix
          if (l==1) call MatSetValues(inf_mat%inf,nrows,idrx,1,idc,bcx(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==2) call MatSetValues(inf_mat%inf,nrows,idry,1,idc,bcy(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==3) call MatSetValues(inf_mat%inf,nrows,idrz,1,idc,bcz(1,1,m),&
               ADD_VALUES,inf_mat%err)
!          print*,mpid%rank,'NEW',interx,intery(l,m),idc,c1,c2,inf_mat%err

          enddo
       enddo

  end subroutine md_influence_matrix_setvalues
!------------------------------------------------------------------------
! md : influence matrix initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_setvalues1(mpid,inf_mat,c,interx,&
       bcx,bcy,bcz,nx,ny,nz,locr)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,l,m,interx,intery(3,2),ncolumns,nrows,c(3),locr
    integer(ik) :: nx,ny,nz,c1,c2
    real(rk),intent(in) :: bcx(ny,nz,2),bcy(nx,nz,2),bcz(nx,ny,2)
    integer(ik) :: idrx(ny*nz),idry(nx*nz),idrz(nx*ny),idc

    !-> compute indexes of where to put values
    call md_get_interfaces_number(inf_mat,c,intery)

    do m=1,2
       do l=1,3
          !-> only process effective interface
          if (intery(l,m)<0) cycle
          
          !-> set number of rows depending of direction of interface
          if (l==1) nrows=inf_mat%ny*inf_mat%nz
          if (l==2) nrows=inf_mat%nx*inf_mat%nz
          if (l==3) nrows=inf_mat%nx*inf_mat%ny

          !-> compute column location
          idc=inf_mat%dm_int_size(interx,3)+locr-1

          !-> compute rows locations
          c1=inf_mat%dm_int_size(intery(l,m),3)-1
          c2=c1+nrows-1
          if (l==1) idrx(:)=(/ (i,i=c1,c2) /)
          if (l==2) idry(:)=(/ (i,i=c1,c2) /)
          if (l==3) idrz(:)=(/ (i,i=c1,c2) /)
          
          !-> put values in influence matrix
          if (l==1) call MatSetValues(inf_mat%inf,1,idc,nrows,idrx,bcx(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==2) call MatSetValues(inf_mat%inf,1,idc,nrows,idry,bcy(1,1,m),&
               ADD_VALUES,inf_mat%err)
          if (l==3) call MatSetValues(inf_mat%inf,1,idc,nrows,idrz,bcz(1,1,m),&
               ADD_VALUES,inf_mat%err)

!          print*,mpid%rank,'NEW',interx,intery(l,m),idc,c1,c2,inf_mat%err

          enddo
       enddo

     end subroutine md_influence_matrix_setvalues1
!------------------------------------------------------------------------
! md : get interface number
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_get_interfaces_number(inf,c,inter)
    implicit none
    type(mpi_inf_mat) :: inf
    integer(ik),intent(in) :: c(3)
    integer(ik),intent(out) :: inter(3,2)
    integer(ik) :: l,m

    do m=1,2
       do l=1,3
          inter(l,m)=inf%dm_int(c(1)+1,c(2)+1,c(3)+1,l,m)
       enddo
    enddo

  end subroutine md_get_interfaces_number
!------------------------------------------------------------------------
! md : petsc influence matrix initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_petsc_influence_matrix_init(mpid,inf_mat,alloc)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: sum,i,j,k,l,m,o,c(3),n,inter
    integer(ik) :: istart,iend
    logical :: alloc

    c(1)=mpid%coord(1) ; c(2)=mpid%coord(2) ; c(3)=mpid%coord(3) 

    !-> compute local number of rows
    sum=0
    if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,1,2)>0) sum=sum+inf_mat%ny*inf_mat%nz
    if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,2,2)>0) sum=sum+inf_mat%nx*inf_mat%nz
    if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,3,2)>0) sum=sum+inf_mat%nx*inf_mat%ny
    inf_mat%rows=sum
!    print*,'ROWS',c,inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,1,2),&
!         inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,2,2),&
!         inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,3,2),inf_mat%rows,sum

    !-> size of non-zeros diagonal elements
    inf_mat%d_nz=inf_mat%rows

    !-> compute size of non-zeros off-diagonal elements
    if (inf_mat%rows>0) then
       allocate(inf_mat%o_nnz(inf_mat%rows))
    else
       allocate(inf_mat%o_nnz(1))
    endif
    inf_mat%o_nnz=0
    allocate(inf_mat%rows_size(inf_mat%it))
    do m=1,inf_mat%it
       sum=0
       do o=1,11
          if (inf_mat%inf_rows_id(m,o,1)>0) then
             if (inf_mat%inf_rows_id(m,o,2)==1) sum=sum+inf_mat%ny*inf_mat%nz
             if (inf_mat%inf_rows_id(m,o,2)==2) sum=sum+inf_mat%nx*inf_mat%nz
             if (inf_mat%inf_rows_id(m,o,2)==3) sum=sum+inf_mat%nx*inf_mat%ny
          endif
       enddo
       inf_mat%rows_size(m)=sum
!       if (mpid%rank==0) print*,'size rows',inf_mat%rows_size(m)
    enddo
    iend=0
    do m=1,inf_mat%it
       istart=1
       if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,1,2)>0) then
          iend=istart+inf_mat%ny*inf_mat%nz-1
          inf_mat%o_nnz(istart:iend)=&
               inf_mat%rows_size(inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,1,2))
          istart=iend+1
       endif
       if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,2,2)>0) then
          iend=istart+inf_mat%nx*inf_mat%nz-1
          inf_mat%o_nnz(istart:iend)=&
               inf_mat%rows_size(inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,2,2))
          istart=iend+1
       endif
       if (inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,3,2)>0) then
          iend=istart+inf_mat%nx*inf_mat%ny-1
          inf_mat%o_nnz(istart:iend)=&
               inf_mat%rows_size(inf_mat%dm_int(c(1)+1,c(2)+1,c(3)+1,3,2))
          istart=iend+1
       endif
    enddo
    inf_mat%o_nnz=inf_mat%o_nnz-inf_mat%rows
    if (inf_mat%rows==0) inf_mat%o_nnz=0

!    print*,mpid%rank,'MATRIX SIZE',mpid%coord,inf_mat%rows,inf_mat%d_nz,'<->',inf_mat%o_nnz

    !-> petsc matrix initialization 
    !-> create influence matrix
    call MatCreate(PETSC_COMM_WORLD,inf_mat%inf,inf_mat%err)

    !-> set local and global dimensions of influence matrix
    call MatSetSizes(inf_mat%inf,inf_mat%rows,  inf_mat%rows,&
         inf_mat%ninf,inf_mat%ninf,inf_mat%err )

    !-> set type of matrix : parallel sparse matrix aij 
    call MatSetType(inf_mat%inf,MATMPIAIJ,inf_mat%err)
    
    if (.not.alloc) then
       !-> preallocate size of rows
       call MatMPIAIJSetPreallocation(inf_mat%inf,inf_mat%d_nz,PETSC_NULL_INTEGER,&
            PETSC_NULL_INTEGER,inf_mat%o_nnz,inf_mat%err)
    endif

    !-> prevent error of preallocation (only for debugging purposes)
    !call MatSetOption(inf_mat%inf,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,inf_mat%err)

    !-> create rhs and solution vectors
    call VecCreateMPI(PETSC_COMM_WORLD,inf_mat%rows,PETSC_DETERMINE,inf_mat%rhs,inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_mat%sol,inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_mat%sol_old(1),inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_mat%sol_old(2),inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_mat%sol_old(3),inf_mat%err)

  end subroutine md_petsc_influence_matrix_init

!------------------------------------------------------------------------
! md : petsc guess vectors initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 12/2012
!
  subroutine md_guess_init(mpid,inf_mat,inf_sol)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    type(mpi_inf_sol) :: inf_sol

    !-> create rhs and solution vectors
    call VecCreateMPI(PETSC_COMM_WORLD,inf_mat%rows,PETSC_DETERMINE,&
         inf_sol%sol_old(1),inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_sol%sol_old(2),inf_mat%err)
    call VecDuplicate(inf_mat%rhs,inf_sol%sol_old(3),inf_mat%err)

  end subroutine md_guess_init

!------------------------------------------------------------------------
! md : influence matrix finalization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_influence_matrix_destroy(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    real(8) :: mem

    !-> Free influence matrix
    call MatDestroy(inf_mat%inf,inf_mat%err)

  end subroutine md_influence_matrix_destroy

!------------------------------------------------------------------------
! md : vertex initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_vertex_init(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: i,j,k,l,m
    integer(ik) :: c(3),ct(3),inter(3,2)
    integer(ik) :: req(6)
    
    !-> get local coord
    call md_mpi_getcoord(mpid,c)
    call md_get_interfaces_number(inf_mat,c,inter)
    c=c
    !print*,mpid%rank,c,'neigh',mpid%neighbours

    inf_mat%vertex=0

    !-> eliminate vertex without neighbours
    if (mpid%neighbours(1,1)<0) inf_mat%vertex(1,:,:)=-1
    if (mpid%neighbours(2,1)<0) inf_mat%vertex(:,1,:)=-1
    if (mpid%neighbours(3,1)<0) inf_mat%vertex(:,:,1)=-1

    if (mpid%neighbours(1,2)<0) inf_mat%vertex(2,:,:)=-1
    if (mpid%neighbours(2,2)<0) inf_mat%vertex(:,2,:)=-1
    if (mpid%neighbours(3,2)<0) inf_mat%vertex(:,:,2)=-1

    !print*,mpid%rank,inf_mat%vertex

    !-> eliminate vertex with boundary conditions
    if (inter(1,1)<0) inf_mat%vertex(1,:,:)=-1
    if (inter(2,1)<0) inf_mat%vertex(:,1,:)=-1
    if (inter(3,1)<0) inf_mat%vertex(:,:,1)=-1

    if (inter(1,2)<0) inf_mat%vertex(2,:,:)=-1
    if (inter(2,2)<0) inf_mat%vertex(:,2,:)=-1
    if (inter(3,2)<0) inf_mat%vertex(:,:,2)=-1

    !-> transfer to neighbours 2 times 
    call md_vertex_info_send(mpid,inf_mat)
    call md_vertex_info_send(mpid,inf_mat)

!    do k=1,2
!       do j=1,2
!          do i=1,2
!             if (inf_mat%vertex(i,j,k)==0) then
!                print*,mpid%coord,'point',i,j,k
!             endif
!          enddo
!       enddo
!    enddo

  end subroutine md_vertex_init

!------------------------------------------------------------------------
! md : vertex initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_vertex_info_send(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer(ik) :: req(6),sendbuf(2,2,3),recvbuf(2,2,3),error,l
    integer :: status(MPI_STATUS_SIZE)

    req=(/1,2,3,4,5,6/)

    !if (mpid%neighbours(1,1)>0) print*,mpid%rank,'recvfrom',mpid%neighbours(1,1),1
    !if (mpid%neighbours(2,1)>0) print*,mpid%rank,'recvfrom',mpid%neighbours(2,1),2
    !if (mpid%neighbours(3,1)>0) print*,mpid%rank,'recvfrom',mpid%neighbours(3,1),3
    
    !if (mpid%neighbours(1,2)>0) print*,mpid%rank,'sendto  ',mpid%neighbours(1,2),1
    !if (mpid%neighbours(2,2)>0) print*,mpid%rank,'sendto  ',mpid%neighbours(2,2),2
    !if (mpid%neighbours(3,2)>0) print*,mpid%rank,'sendto  ',mpid%neighbours(3,2),3

    !print*,'NEIGH',mpid%rank,mpid%neighbours

    !-> begin tranfer boundary condition to neighbours
    if (mpid%neighbours(1,1)>=0) then
       call mpi_irecv(recvbuf(1,1,1),4,mpi_integer,mpid%neighbours(1,1),1,&
            mpi_comm_world,req(1),error)
    endif
    if (mpid%neighbours(2,1)>=0) then
       call mpi_irecv(recvbuf(1,1,2),4,mpi_integer,mpid%neighbours(2,1),2,&
            mpi_comm_world,req(2),error)
    endif
    if (mpid%neighbours(3,1)>=0) then
       call mpi_irecv(recvbuf(1,1,3),4,mpi_integer,mpid%neighbours(3,1),3,&
            mpi_comm_world,req(3),error)
    endif

    if (mpid%neighbours(1,2)>=0) then
       sendbuf(:,:,1)=inf_mat%vertex(2,:,:)
       call mpi_issend(sendbuf(1,1,1),4,mpi_integer,mpid%neighbours(1,2),1,&
            mpi_comm_world,req(4),error)
    endif
    if (mpid%neighbours(2,2)>=0) then
       sendbuf(:,:,2)=inf_mat%vertex(:,2,:)
       call mpi_issend(sendbuf(1,1,2),4,mpi_integer,mpid%neighbours(2,2),2,&
            mpi_comm_world,req(5),error)
    endif
    if (mpid%neighbours(3,2)>=0) then
       sendbuf(:,:,3)=inf_mat%vertex(:,:,2)
       call mpi_issend(sendbuf(1,1,3),4,mpi_integer,mpid%neighbours(3,2),3,&
            mpi_comm_world,req(6),error)
    endif

    !-> end tranfer boundary condition to neighbours
    if (mpid%neighbours(1,1)>=0) then 
       call mpi_wait(req(1),status,error)
       where(recvbuf(:,:,1)<0) inf_mat%vertex(1,:,:)=-1
    endif
    if (mpid%neighbours(2,1)>=0) then 
       call mpi_wait(req(2),status,error)
       where(recvbuf(:,:,2)<0) inf_mat%vertex(:,1,:)=-1
    endif
    if (mpid%neighbours(3,1)>=0) then 
       call mpi_wait(req(3),status,error)
       where(recvbuf(:,:,3)<0) inf_mat%vertex(:,:,1)=-1
    endif

    if (mpid%neighbours(1,2)>=0) call mpi_wait(req(4),error)
    if (mpid%neighbours(2,2)>=0) call mpi_wait(req(5),error)
    if (mpid%neighbours(3,2)>=0) call mpi_wait(req(6),error)

  end subroutine md_vertex_info_send

!------------------------------------------------------------------------
! md : influence matrix parameters
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_influence_matrix_param(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    integer :: i,j,k,l,n,m,o,p,sum,scale,maxgrey
    integer,allocatable :: infgrey(:,:)

    !-> number of domains
    inf_mat%nd(1)=mpid%nd(1) ; inf_mat%nd(2)=mpid%nd(2) ; inf_mat%nd(3)=mpid%nd(3) 

    !-> number of points
    inf_mat%nx=mpid%nx-2 ; inf_mat%ny=mpid%ny-2 ; inf_mat%nz=mpid%nz-2 

    !-> allocate and initialize dm_int matrix
    if (.not.allocated(inf_mat%dm_int)) then
       allocate(inf_mat%dm_int(inf_mat%nd(1),inf_mat%nd(2),inf_mat%nd(3),3,2))
    else
       deallocate(inf_mat%dm_int)
       allocate(inf_mat%dm_int(inf_mat%nd(1),inf_mat%nd(2),inf_mat%nd(3),3,2))
    endif
    inf_mat%dm_int=-1

    !-> declare inside boundary condition into dm_int
    !inf_mat%dm_int(1,1,1,1,2)=0
    !inf_mat%dm_int(1,1,1,2,2)=0
    !inf_mat%dm_int(1,1,1,3,2)=0

    !-> declare periodicity
    inf_mat%periods(1)=mpid%periods(1)
    inf_mat%periods(2)=mpid%periods(2)
    inf_mat%periods(3)=mpid%periods(3)

    !-> number of interfaces
    !-> x direction
    if (inf_mat%periods(1)) then
       inf_mat%ix=(inf_mat%nd(1))*(inf_mat%nd(2)*inf_mat%nd(3))
    else
       inf_mat%ix=(inf_mat%nd(1)-1)*(inf_mat%nd(2)*inf_mat%nd(3))
    endif
    inf_mat%ix=inf_mat%ix-count(inf_mat%dm_int(:,:,:,1,2)==0)

    !-> y direction
    if (inf_mat%periods(2)) then
       inf_mat%iy=(inf_mat%nd(2))*(inf_mat%nd(1)*inf_mat%nd(3))
    else
       inf_mat%iy=(inf_mat%nd(2)-1)*(inf_mat%nd(1)*inf_mat%nd(3))
    endif
    inf_mat%iy=inf_mat%iy-count(inf_mat%dm_int(:,:,:,2,2)==0)

    !-> z direction
    if (inf_mat%periods(3)) then
       inf_mat%iz=(inf_mat%nd(3))*(inf_mat%nd(1)*inf_mat%nd(2))
    else
       inf_mat%iz=(inf_mat%nd(3)-1)*(inf_mat%nd(1)*inf_mat%nd(2))
    endif
    inf_mat%iz=inf_mat%iz-count(inf_mat%dm_int(:,:,:,3,2)==0)

    !-> total number of interfaces 
    inf_mat%it=inf_mat%ix+inf_mat%iy+inf_mat%iz

    !-> allocate size information of influence matrix
    if (.not.allocated(inf_mat%dm_int_size)) then
       allocate(inf_mat%dm_int_size(inf_mat%it,3))
    else
       deallocate(inf_mat%dm_int_size)
       allocate(inf_mat%dm_int_size(inf_mat%it,3))
    endif

    !-> numbering of interfaces
    
    sum=0
    do i=1,inf_mat%nd(1)
       do j=1,inf_mat%nd(2)
          do k=1,inf_mat%nd(3)
             if (inf_mat%dm_int(i,j,k,1,2)/=0) then
                if (i<inf_mat%nd(1)) then 
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,1,2)=sum
                   inf_mat%dm_int_size(sum,1)=1
                elseif (inf_mat%periods(1)) then
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,1,2)=sum
                   inf_mat%dm_int_size(sum,1)=1
                endif
             endif
             if (inf_mat%dm_int(i,j,k,2,2)/=0) then
                if (j<inf_mat%nd(2)) then 
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,2,2)=sum
                   inf_mat%dm_int_size(sum,1)=2
                elseif (inf_mat%periods(2)) then
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,2,2)=sum
                   inf_mat%dm_int_size(sum,1)=2
                endif
             endif
             if (inf_mat%dm_int(i,j,k,3,2)/=0) then
                if (k<inf_mat%nd(3)) then
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,3,2)=sum
                   inf_mat%dm_int_size(sum,1)=3
                elseif (inf_mat%periods(3)) then
                   sum=sum+1
                   inf_mat%dm_int(i,j,k,3,2)=sum
                   inf_mat%dm_int_size(sum,1)=3
                endif
             endif
          enddo
       enddo
    enddo

    inf_mat%dm_int(2:inf_mat%nd(1),:,:,1,1)=inf_mat%dm_int(1:inf_mat%nd(1)-1,:,:,1,2)
    if (inf_mat%periods(1)) inf_mat%dm_int(1:1,:,:,1,1)=&
         inf_mat%dm_int(inf_mat%nd(1):inf_mat%nd(1),:,:,1,2)

    inf_mat%dm_int(:,2:inf_mat%nd(2),:,2,1)=inf_mat%dm_int(:,1:inf_mat%nd(2)-1,:,2,2)
    if (inf_mat%periods(2)) inf_mat%dm_int(:,1:1,:,2,1)=&
         inf_mat%dm_int(:,inf_mat%nd(2):inf_mat%nd(2),:,2,2)

    inf_mat%dm_int(:,:,2:inf_mat%nd(3),3,1)=inf_mat%dm_int(:,:,1:inf_mat%nd(3)-1,3,2)
    if (inf_mat%periods(3)) inf_mat%dm_int(:,:,1:1,3,1)=&
         inf_mat%dm_int(:,:,inf_mat%nd(3):inf_mat%nd(3),3,2)


    !-> replace zero values of inside interfaces to -1 
    where(inf_mat%dm_int(:,:,:,:,:)==0) inf_mat%dm_int(:,:,:,:,:)=-1

    !-> write domain interfaces informations
    if (mpid%rank==0) then
       open(4000,file='domain_interfaces.dat')
       do i=1,inf_mat%nd(1)
          do j=1,inf_mat%nd(2)
             do k=1,inf_mat%nd(3)
                write(4000,'(3i4,3x,6i4)')i,j,k,inf_mat%dm_int(i,j,k,:,:)
             enddo
          enddo
       enddo
       close(4000)
    endif

    !-> influence matrix size and id per row
    if (.not.allocated(inf_mat%inf_rows_id)) then
       allocate(inf_mat%inf_rows_id(inf_mat%it,11,2))
    else
       deallocate(inf_mat%inf_rows_id)
       allocate(inf_mat%inf_rows_id(inf_mat%it,11,2))
    endif
    inf_mat%inf_rows_id=-1

    do m=1,inf_mat%it
       sum=0
       p=1 ; inf_mat%inf_rows_id(m,p,1)=m ; inf_mat%inf_rows_id(m,p,2)=inf_mat%dm_int_size(m,1)
       do k=1,inf_mat%nd(3)
          do j=1,inf_mat%nd(2)
             do i=1,inf_mat%nd(1)
                if (any(inf_mat%dm_int(i,j,k,:,:)==m)) then
                   do n=1,3
                      do o=1,2
                         if (inf_mat%dm_int(i,j,k,n,o)>0) then
                            sum=sum+1
                            if (inf_mat%dm_int(i,j,k,n,o)/=m) then
                               p=p+1
                               inf_mat%inf_rows_id(m,p,1)=inf_mat%dm_int(i,j,k,n,o)
                               if (n==1) inf_mat%inf_rows_id(m,p,2)=1
                               if (n==2) inf_mat%inf_rows_id(m,p,2)=2
                               if (n==3) inf_mat%inf_rows_id(m,p,2)=3
                            endif
                         endif
                      enddo
                   enddo
                endif
             enddo
          enddo
       enddo
       inf_mat%dm_int_size(m,2)=sum-1
    enddo

    !-> total size of influence matrix
    sum=0
    inf_mat%dm_int_size(1,3)=1
    do i=1,inf_mat%it
       if (inf_mat%dm_int_size(i,1)==1) sum=sum+inf_mat%ny*inf_mat%nz
       if (inf_mat%dm_int_size(i,1)==2) sum=sum+inf_mat%nx*inf_mat%nz
       if (inf_mat%dm_int_size(i,1)==3) sum=sum+inf_mat%nx*inf_mat%ny
       if (i<inf_mat%it) inf_mat%dm_int_size(i+1,3)=sum+1
    enddo
    inf_mat%ninf=sum

    !-> influence matrix image
    if (mpid%rank==0) then
       allocate(infgrey(inf_mat%it,inf_mat%it))
       infgrey=256
       do m=1,inf_mat%it
          do k=1,inf_mat%nd(3)
             do j=1,inf_mat%nd(2)
                do i=1,inf_mat%nd(1)
                   if (any(inf_mat%dm_int(i,j,k,:,:)==m)) then
                      do n=1,3
                         do o=1,2
                            if (inf_mat%dm_int(i,j,k,n,o)>0.and.n==1) &
                                 infgrey(m,inf_mat%dm_int(i,j,k,n,o))=0
                            if (inf_mat%dm_int(i,j,k,n,o)>0.and.n==2) &
                                 infgrey(m,inf_mat%dm_int(i,j,k,n,o))=84
                            if (inf_mat%dm_int(i,j,k,n,o)>0.and.n==3) &
                                 infgrey(m,inf_mat%dm_int(i,j,k,n,o))=84*2
                         enddo
                      enddo
                   endif
                enddo
             enddo
          enddo
       enddo

       !-> create grey image of influence matrix
       maxgrey=maxval(infgrey)
       scale=1000/inf_mat%it
       if (scale==0) scale=1
       !scale=(inf_mat%ninf/inf_mat%it)/10
       !print*,scale,inf_mat%ninf,inf_mat%it
       open(10,file='influence_matrix.pgm')
       write(10,'(a)')'P2'
       write(10,'(a)')'#Created By M. Marquillie'
       write(10,*)inf_mat%it*scale,inf_mat%it*scale
       write(10,*)maxgrey
       write(10,'(10(i0,x))')((((infgrey(i,j),k=1,scale),j=1,inf_mat%it),&
            l=1,scale),i=1,inf_mat%it)
       close(10)
       deallocate(infgrey)
    endif

  end subroutine md_influence_matrix_param

!------------------------------------------------------------------------
! md : influence matrix view
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_influence_matrix_view(mpid,inf_mat)
    implicit none
    type(mpi_data) :: mpid
    type(mpi_inf_mat) :: inf_mat
    PetscViewer viewer
    integer(ik) :: rows,columns
    integer(ik) :: firstr,lastr,firstc,lastc
    real(rk) :: norm

    !-> get matrix rows and columns informations
!    call MatGetOwnershipRange(inf_mat%inf,firstr,lastr,inf_mat%err)
!    call MatGetOwnershipRangeColumn(inf_mat%inf,firstc,lastc,inf_mat%err)
!    call MatGetLocalSize(inf_mat%inf,rows,columns,inf_mat%err)
!    print*,'Mat local size',mpid%rank,rows,columns,firstr,lastr,lastr-firstr, &
!         firstc,lastc,lastc-firstc

    !-> get vec rows and columns informations
!    call VecGetLocalSize(inf_mat%rhs,rows,inf_mat%err)
!    print*,'Vec local size',mpid%rank,rows 


    !-> view matrix (only for small matrices < 1024 rows)
!    call PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_DENSE,inf_mat%err)
!    call MatView(inf_mat%inf,PETSC_VIEWER_STDOUT_WORLD,inf_mat%err)
!    call MatView(inf_mat%inf,PETSC_VIEWER_DRAW_WORLD,inf_mat%err)


!    call MatNorm(inf_mat%inf,NORM_INFINITY,norm,inf_mat%inf)
!    print*,'NORM',norm

    !-> write influence matrix to file
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'test.dat',viewer,inf_mat%err)
    call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_DENSE,inf_mat%err)
    call MatView(inf_mat%inf,viewer,inf_mat%err)    
    call PetscViewerDestroy(viewer,inf_mat%err)


  end subroutine md_influence_matrix_view

! =======================================================================
! =======================================================================
! md : mpi
! =======================================================================
! =======================================================================
!------------------------------------------------------------------------
! md : mpi initialization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_mpi_init(mpid,cmd)
    implicit none
    type(mpi_data) :: mpid
    type(cmd_line) :: cmd
    logical :: reorganization
    integer(ik) :: i,j,k,coord(3),err

    !-> initialize mpi
    call mpi_init(mpid%code)
    mpid%comm=MPI_COMM_WORLD
    call mpi_comm_size(mpid%comm,mpid%processus,mpid%code)
    call mpi_comm_rank(mpid%comm,mpid%rank,mpid%code)

    !-> mpi topology parameters
    !-> number of dimensions
    mpid%dims=3

    !-> number of domains
    mpid%nd(1)=cmd%ndx ; mpid%nd(2)=cmd%ndy ;mpid%nd(3)=cmd%ndz 

    !-> number of points
    mpid%nx=cmd%nx ; mpid%ny=cmd%ny ;mpid%nz=cmd%nz 

    !-> total number of domains
    mpid%ndom=mpid%nd(1)*mpid%nd(2)*mpid%nd(3)

    !-> test if number of processus = number of domains
    if (mpid%processus/=mpid%ndom) then
       if (mpid%rank==0) then
          print'(a)',"ERROR : number of processus and number of domains&
               & are different"
       endif
       call mpi_finalize(mpid%code)
       stop
    endif

    !-> declare periodicity
    mpid%periods(:)=.false.
    if (cmd%periods(1)==1) mpid%periods(1)=.true.
    if (cmd%periods(2)==1) mpid%periods(2)=.true.
    if (cmd%periods(3)==1) mpid%periods(3)=.true.

    !-> reorganization
    reorganization=.false.

    !-> create topology
    call mpi_cart_create(mpid%comm,mpid%dims,mpid%nd,mpid%periods, &
         reorganization,mpid%comm_3d,mpid%code)

    !-> get rank and coordinates in cartesian topology
    call mpi_comm_rank(mpid%comm_3d,mpid%rank_3d,mpid%code)
    call mpi_cart_coords(mpid%comm_3d,mpid%rank_3d,mpid%dims,mpid%coord,mpid%code)

    !-> get neighbours rank
    !-> x-direction
    call mpi_cart_shift(mpid%comm_3d,0,1,mpid%neighbours(1,1),mpid%neighbours(1,2),&
         mpid%code)
    !-> y-direction
    call mpi_cart_shift(mpid%comm_3d,1,1,mpid%neighbours(2,1),mpid%neighbours(2,2),&
         mpid%code)
    !-> z-direction
    call mpi_cart_shift(mpid%comm_3d,2,1,mpid%neighbours(3,1),mpid%neighbours(3,2),&
         mpid%code)

    !-> get all processus rank in cartesian topology
    if (.not.allocated(mpid%domains)) then
       allocate(mpid%domains(mpid%nd(1),mpid%nd(2),mpid%nd(3)))
    else
       deallocate(mpid%domains) 
       allocate(mpid%domains(mpid%nd(1),mpid%nd(2),mpid%nd(3)))
    endif
    do k=1,mpid%nd(3)
       do j=1,mpid%nd(2)
          do i=1,mpid%nd(1)
             coord=(/i-1,j-1,k-1/)
             call mpi_cart_rank(mpid%comm_3d,coord,mpid%domains(i,j,k),mpid%code)
          enddo
       enddo
    enddo

  end subroutine md_mpi_init
!------------------------------------------------------------------------
! md : mpi finalization
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_mpi_finalize(mpid)
    implicit none
    type(mpi_data) :: mpid

    call mpi_finalize(mpid%code)

  end subroutine md_mpi_finalize
!------------------------------------------------------------------------
! md : mpi barrier
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_mpi_barrier(mpid)
    implicit none
    type(mpi_data) :: mpid

    call mpi_barrier(mpid%comm_3d,mpid%code)

  end subroutine md_mpi_barrier

!------------------------------------------------------------------------
! md : mpi reduce single double
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_mpi_reduce_double(mpid,val,sum)
    implicit none
    type(mpi_data) :: mpid
    real(rk),intent(in) :: val
    real(rk),intent(out) :: sum

    call mpi_reduce(val,sum,1,mpi_double_precision,mpi_sum,0,&
         mpi_comm_world,mpid%code)

  end subroutine md_mpi_reduce_double
!------------------------------------------------------------------------
! md : mpi broadcast single double
!------------------------------------------------------------------------
! Matthieu Marquillie
! 07/2012
!
  subroutine md_mpi_bcast_double(mpid,val,rank)
    implicit none
    type(mpi_data) :: mpid
    integer(ik),intent(in) :: rank
    real(rk),intent(in) :: val

    call mpi_bcast(val,1,mpi_double_precision,rank,mpi_comm_world,mpid%code)

  end subroutine md_mpi_bcast_double
!------------------------------------------------------------------------
! md : mpi barrier
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_mpi_getcoord(mpid,c)
    implicit none
    type(mpi_data) :: mpid
    integer(ik),intent(out) :: c(3)

    c(1)=mpid%coord(1) ; c(2)=mpid%coord(2) ; c(3)=mpid%coord(3) 

  end subroutine md_mpi_getcoord
!------------------------------------------------------------------------
! md : mpi barrier
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_mpi_getnumberdomains(mpid,c)
    implicit none
    type(mpi_data) :: mpid
    integer(ik),intent(out) :: c(3)

    c(1)=mpid%nd(1) ; c(2)=mpid%nd(2) ; c(3)=mpid%nd(3) 

  end subroutine md_mpi_getnumberdomains
!------------------------------------------------------------------------
! md : mpi barrier
!------------------------------------------------------------------------
! Matthieu Marquillie
! 09/2012
!
  subroutine md_mpi_global_coord(mpid,dim,coord)
    implicit none
    type(mpi_data) :: mpid
    integer(ik),intent(out) :: dim(3),coord(3,2)

    !-> compute total number of points
    dim(1)=mpid%nx*mpid%nd(1)
    dim(2)=mpid%ny*mpid%nd(2)
    dim(3)=mpid%nz*mpid%nd(3)

    !-> compute global indices
    coord(1,1)=mpid%coord(1)*mpid%nx+1
    coord(2,1)=mpid%coord(2)*mpid%ny+1
    coord(3,1)=mpid%coord(3)*mpid%nz+1

    coord(1,2)=mpid%nx
    coord(2,2)=mpid%ny
    coord(3,2)=mpid%nz

  end subroutine md_mpi_global_coord

end module class_md
