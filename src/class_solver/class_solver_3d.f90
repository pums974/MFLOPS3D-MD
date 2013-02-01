module class_solver_3d
! -----------------------------------------------------------------------
! Name :
! class solver 3d
! -----------------------------------------------------------------------
! Object : 
! Solve poisson and helmholtz equation with a 5-7 compact finite 
! differences scheme on irregular grid
! -----------------------------------------------------------------------
! Files :
! class_solver_3d.f90
! class_solver_coefficient_3d.f90
! -----------------------------------------------------------------------
! Public type :
! solver_coefficients -> coefficients for solver
! -----------------------------------------------------------------------
! Public Procedure :
! solver_coefficients_init -> initialize solver coefficients 
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
!  use precision
  use class_solver_coefficient_3d
  implicit none

  !-> Define type containing solver coefficients for one direction
  type,public :: solver_coeffs
     private
     !-> nall : total size ; n : size without boundaries
     integer(ik) :: nall,n
     !-> second derivatives coefficients : lhs and rhs
     real(rk),allocatable :: ddl(:,:),ddr(:,:)
     !-> lhs : inverse matrix 
     real(rk),allocatable :: inv_lhs(:,:)
     !-> lhs : eigenvalues, eigenvectors, inverse eigenvectors
     real(rk),allocatable :: eigen(:),vect(:,:),inv_vect(:,:)
     !-> boundaries type (dirichlet or neumann)
     integer(ik) :: bc(2)
     !-> coefficients for interpolation
     real(rk),allocatable :: neum(:,:)
  end type solver_coeffs

  !-> Define type containing all solver coefficients
  type,public :: solver_coeffs_3d
     private
     type(solver_coeffs) :: cx,cy,cz
  end type solver_coeffs_3d

  !-> Declare everything private by default
  private

  !-> Declare exported procedure
  public :: solver_coeffs_init_3d,solve_3d
  
contains

  subroutine print(n,ddl,ddr)
! -----------------------------------------------------------------------
! solver3d : print derivatives coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 02/2012
!
    implicit none
    integer,intent(in) :: n
    real(rk),intent(in) :: ddl(n,5),ddr(n,8)
    integer(ik) :: i,j

    print*,'ddl'
    do i=1,n
       print'(5es17.8)',(ddl(i,j),j=1,5)
    enddo
    print*,'ddr'
    do i=1,n
       print'(8es17.8)',(ddr(i,j),j=1,8)
    enddo

  end subroutine print

  subroutine solve_3d(sc,rhs,sigma,solf,bcx,bcy,bcz,n1,n2,n3)
! -----------------------------------------------------------------------
! solver3d : solve system A*x=B with matrix diagonalization 
!            in the three directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
    !$ use OMP_LIB
    implicit none
    type(solver_coeffs_3d),intent(in) :: sc
    integer(ik),intent(in) :: n1,n2,n3
    real(rk),intent(inout) :: rhs(n1,n2,n3)
    real(rk),intent(inout) :: solf(0:n1+1,0:n2+1,0:n3+1)
    real(rk),intent(in) :: bcx(n2,n3,2),bcy(n1,n3,2),bcz(n1,n2,2)
    real(rk),intent(in) :: sigma
    real(rk) :: sol(n1,n2,n3)
    real(rk) :: aux1(n1,n2,n3),aux2(n1,n2,n3),aux3i(n2,n3),aux3f(n2,n3)
    real(rk) :: clx(n1),cly(n2),clz(n3)
    integer(ik) :: i,j,k
    real(rk) :: t1,t2

    !-> apply boundary conditions to rhs --------------------------

    call cpu_time(t1)
    !-> x-direction
    clx=0._rk
!$OMP PARALLEL PRIVATE(k,clx)
!$OMP DO SCHEDULE(RUNTIME)
    do k=1,n3
       do j=1,n2
          clx(1)=sc%cx%ddr(1,1)*bcx(j,k,1)
          clx(2)=sc%cx%ddr(2,1)*bcx(j,k,1)
          clx(n1-1)=sc%cx%ddr(n1-1,8)*bcx(j,k,2)
          clx(n1)=sc%cx%ddr(n1,7)*bcx(j,k,2)

          !rhs(:,j,k)=rhs(:,j,k)-matmul(sc%cx%inv_lhs,clx)

          rhs(:,j,k)=rhs(:,j,k)&
               -sc%cx%inv_lhs(:,1)*clx(1) &
               -sc%cx%inv_lhs(:,2)*clx(2) &
               -sc%cx%inv_lhs(:,n1-1)*clx(n1-1) &
               -sc%cx%inv_lhs(:,n1)*clx(n1) 
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

    !-> y-direction
    cly=0._rk
!$OMP PARALLEL PRIVATE(k,cly)
!$OMP DO SCHEDULE(RUNTIME)
    do k=1,n3
       do i=1,n1
          cly(1)=sc%cy%ddr(1,1)*bcy(i,k,1)
          cly(2)=sc%cy%ddr(2,1)*bcy(i,k,1)
          cly(n2-1)=sc%cy%ddr(n2-1,8)*bcy(i,k,2)
          cly(n2)=sc%cy%ddr(n2,7)*bcy(i,k,2)

          !rhs(i,:,k)=rhs(i,:,k)-matmul(sc%cy%inv_lhs,cly)

          rhs(i,:,k)=rhs(i,:,k)&
               -sc%cy%inv_lhs(:,1)*cly(1) &
               -sc%cy%inv_lhs(:,2)*cly(2) &
               -sc%cy%inv_lhs(:,n2-1)*cly(n2-1) &
               -sc%cy%inv_lhs(:,n2)*cly(n2) 
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    
    !-> z-direction
    clz=0._rk
!$OMP PARALLEL PRIVATE(j,clz)
!$OMP DO SCHEDULE(RUNTIME)
    do j=1,n2
       do i=1,n1
          clz(1)=sc%cz%ddr(1,1)*bcz(i,j,1)
          clz(2)=sc%cz%ddr(2,1)*bcz(i,j,1)
          clz(n3-1)=sc%cz%ddr(n3-1,8)*bcz(i,j,2)
          clz(n3)=sc%cz%ddr(n3,7)*bcz(i,j,2)


          !rhs(i,j,:)=rhs(i,j,:)-matmul(sc%cz%inv_lhs,clz)

          rhs(i,j,:)=rhs(i,j,:)&
               -sc%cz%inv_lhs(:,1)*clz(1) &
               -sc%cz%inv_lhs(:,2)*clz(2) &
               -sc%cz%inv_lhs(:,n3-1)*clz(n3-1) &
               -sc%cz%inv_lhs(:,n3)*clz(n3) 
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    call cpu_time(t2)
!    print*,'Boundary condition',t2-t1

    !-> compute modified rhs ------------------------------------

    !->  x-direction
    call cpu_time(t1)
    call dgx(n1,n2,n3,sc%cx%inv_vect,rhs,aux1)
    call cpu_time(t2)
!    print*,'dgemm 1',t2-t1,(2._rk*real(n1)**2*real(n2)*real(n3)/(1000.d0*1000.d0))/(t2-t1)

    !-> y-direction
    call cpu_time(t1)
    call dgy(n1,n2,n3,sc%cy%inv_vect,aux1,aux2)
    call cpu_time(t2)
!    print*,'dgemm 2',t2-t1,(2._rk*real(n1)*real(n2)**2*real(n3)/(1000.d0*1000.d0))/(t2-t1)

    !-> z-direction
    call cpu_time(t1)
    call dgz(n1,n2,n3,sc%cz%inv_vect,aux2,aux1)
    call cpu_time(t2)
!    print*,'dgemm 3',t2-t1,(2._rk*real(n1)*real(n2)*real(n3)**2/(1000.d0*1000.d0))/(t2-t1)
!    print*,'rhs matmul',t2-t1
    
    !-> compute solution with eigenvalues -------------------------
    call cpu_time(t1)
!$OMP PARALLEL PRIVATE(k)
!$OMP DO SCHEDULE(RUNTIME)
    do k=1,n3
       do j=1,n2
          do i=1,n1
             sol(i,j,k)=aux1(i,j,k)/(sc%cx%eigen(i)+sc%cy%eigen(j)+sc%cz%eigen(k)+sigma)

             !-> For 1D problems
!             sol(i,j,k)=aux1(i,j,k)/(sc%cx%eigen(i))
             !-> for neumann bc -> put zero value where all eigenvalues equals zero
             !if (i==n1-1.and.j==n2-1.and.k==n3-1) then
             !   print*,aux1(i,j,k),(sc%cx%eigen(i)+sc%cy%eigen(j)+sc%cz%eigen(k))
             !   print*,aux1(i,j,k)/(sc%cx%eigen(i)+sc%cy%eigen(j)+sc%cz%eigen(k))
             !   sol(i,j,k)=0._rk
             !endif
          enddo
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    call cpu_time(t2)
!    print*,'Eigen',t2-t1

    !-> compute solution -------------------------------------------

    !->  x-direction
    call cpu_time(t1)
    call dgx(n1,n2,n3,sc%cx%vect,sol,aux1)
    call cpu_time(t2)
!    print*,'dgemm 1',t2-t1,(2._rk*real(n1)**2*real(n2)*real(n3)/(1000.d0*1000.d0))/(t2-t1)
    
    !-> y-direction
    call cpu_time(t1)
    call dgy(n1,n2,n3,sc%cy%vect,aux1,aux2)
    call cpu_time(t2)
!    print*,'dgemm 2',t2-t1,(2._rk*real(n1)*real(n2)**2*real(n3)/(1000.d0*1000.d0))/(t2-t1)

    !-> z-direction
    call cpu_time(t1)
    call dgz(n1,n2,n3,sc%cz%vect,aux2,sol)
    call cpu_time(t2)
!    print*,'dgemm 3',t2-t1,(2._rk*real(n1)*real(n2)*real(n3)**2/(1000.d0*1000.d0))/(t2-t1)
!    print*,'U matmul',t2-t1

    call cpu_time(t1)
!$OMP PARALLEL PRIVATE(k)
!$OMP DO SCHEDULE(RUNTIME)
    do k=1,n3
       solf(1:n1,1:n2,k)=sol(:,:,k)
    enddo
!$OMP END DO
!$OMP END PARALLEL
    call cpu_time(t2)
!    print*,'copy',t2-t1

    !-> put boundary condition inside solf

    !->  x-direction
    if (sc%cx%bc(1)==1) then
       solf(0,1:n2,1:n3)=bcx(:,:,1)
    elseif (sc%cx%bc(1)==2) then
       call extrapolation_neum_x(sc%cx%neum,solf,bcx,n1,n2,n3,1)
    endif
    if (sc%cx%bc(2)==1) then
       solf(n1+1,1:n2,1:n3)=bcx(:,:,2)
    elseif (sc%cx%bc(2)==2) then
       call extrapolation_neum_x(sc%cx%neum,solf,bcx,n1,n2,n3,2)
    endif

    !-> y-direction
    if (sc%cy%bc(1)==1) then
       solf(1:n1,0,1:n3)=bcy(:,:,1)
    elseif (sc%cy%bc(1)==2) then
       call extrapolation_neum_y(sc%cy%neum,solf,bcy,n1,n2,n3,1)
    endif
    if (sc%cy%bc(2)==1) then
       solf(1:n1,n2+1,1:n3)=bcy(:,:,2)
    elseif (sc%cy%bc(2)==2) then
       call extrapolation_neum_y(sc%cy%neum,solf,bcy,n1,n2,n3,2)
    endif

    !-> z-directionn
    if (sc%cz%bc(1)==1) then
       solf(1:n1,1:n2,0)=bcz(:,:,1)
    elseif (sc%cz%bc(1)==2) then
       call extrapolation_neum_z(sc%cz%neum,solf,bcz,n1,n2,n3,1)
    endif
    if (sc%cz%bc(2)==1) then
       solf(1:n1,1:n2,n3+1)=bcz(:,:,2)
    elseif (sc%cz%bc(2)==2) then
       call extrapolation_neum_z(sc%cz%neum,solf,bcz,n1,n2,n3,2)
    endif


  end subroutine solve_3d


! -----------------------------------------------------------------------
! solver3d : solve system A*x=B with matrix diagonalization 
!            in the three directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
  subroutine dgx(n1,n2,n3,a,b,c)
    !$ use OMP_LIB
    implicit none
    integer :: n1,n2,n3
    real(rk) :: a(n1,n1),b(n1,n2,n3),c(n1,n2,n3)
    integer :: k
    
    call dgemm('n','n',n1,n2*n3,n1,1.d0,a,n1,b,n1,0.d0,c,n1)
!!$OMP PARALLEL PRIVATE(k)
!!$OMP DO SCHEDULE(RUNTIME)
!    do k=1,n3
!       call dgemm('n','n',n1,n2,n1,1.d0,a,n1,b(1,1,k),n1,0.d0,c(1,1,k),n1)
!    enddo
!!$OMP END DO
!!$OMP END PARALLEL

  end subroutine dgx

  
! -----------------------------------------------------------------------
! solver3d : solve system A*x=B with matrix diagonalization 
!            in the three directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
  subroutine dgy(n1,n2,n3,a,b,c)
    !$ use OMP_LIB
    implicit none
    integer :: n1,n2,n3
    real(rk) :: a(n2,n2),b(n1,n2,n3),c(n1,n2,n3)
    integer :: k
    
!$OMP PARALLEL PRIVATE(k)
!$OMP DO SCHEDULE(RUNTIME)
    do k=1,n3
       call dgemm('n','n',n1,n2,n2,1.d0,b(1,1,k),n1,a,n2,0.d0,c(1,1,k),n1)
!       c(:,:,k)=matmul(b(:,:,k),a)
    enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine dgy

! -----------------------------------------------------------------------
! solver3d : solve system A*x=B with matrix diagonalization 
!            in the three directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
  subroutine dgz(n1,n2,n3,a,b,c)
    !$ use OMP_LIB
    implicit none
    integer :: n1,n2,n3
    real(rk) :: a(n3,n3),b(n1,n2,n3),c(n1,n2,n3)
    real(rk) :: aux3i(n2,n3),aux3f(n2,n3)
    integer :: i

!!$OMP PARALLEL PRIVATE(i,aux3i,aux3f)
!!$OMP DO SCHEDULE(RUNTIME)
!    do i=1,n1
!       aux3i(:,:)=b(i,:,:)
!       call dgemm('n','n',n2,n3,n3,1.d0,aux3i,n2,a,n3,0.d0,aux3f,n2)
!       !call dgemm('n','n',n2,n3,n3,1.d0,b(i,:,:),n2,a,n3,0.d0,c(i,:,:),n2)
!       c(i,:,:)=aux3f(:,:)
!       !aux1(i,:,:)=matmul(aux2(i,:,:),sc%cz%inv_vect)
!    enddo
!!$OMP END DO
!!$OMP END PARALLEL

    call dgemm('n','n',n1*n2,n3,n3,1.d0,b,n1*n2,a,n3,0.d0,c,n1*n2)

  end subroutine dgz

  subroutine solver_coeffs_init_3d(grid1,grid2,grid3,sc,n1,n2,n3,bct,out_name)
! -----------------------------------------------------------------------
! solver3d : Initialisation of all 3d solver coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    type(solver_coeffs_3d),intent(out) :: sc
    integer(ik),intent(in) :: n1,n2,n3
    real(rk),intent(in) :: grid1(n1),grid2(n2),grid3(n3)
    integer(ik),intent(in) :: bct(6) ! bc type (Diri=1/Neum=2) : x1,x2,y1,y2,z1,z2
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    
!    print'(a)','Solver : x-direction -------------------------------------------'
    call solver_coeffs_init(grid1,sc%cx,n1,bct(1),out_name)

!    print'(a)','Solver : y-direction -------------------------------------------'
    call solver_coeffs_init(grid2,sc%cy,n2,bct(3),out_name)
    sc%cy%inv_vect=transpose(sc%cy%inv_vect)
    sc%cy%vect=transpose(sc%cy%vect)

!    print'(a)','Solver : z-direction -------------------------------------------'
    call solver_coeffs_init(grid3,sc%cz,n3,bct(5),out_name)
    sc%cz%inv_vect=transpose(sc%cz%inv_vect)
    sc%cz%vect=transpose(sc%cz%vect)
    


  end subroutine solver_coeffs_init_3d

  subroutine solver_coeffs_init(grid,dc,n,bct,out_name)
! -----------------------------------------------------------------------
! solver3d : Initialisation of solver coefficients in one direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    type(solver_coeffs),intent(out) :: dc
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: grid(n)
    integer(ik),intent(in) :: bct(2)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    
    !-> initialize and allocate dc

    if (.not.allocated(dc%ddl)) then
       dc%nall=n
       dc%n=n-2 
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,8))
       allocate(dc%inv_lhs(n-2,n-2))
       allocate(dc%eigen(n-2),dc%vect(n-2,n-2),dc%inv_vect(n-2,n-2))
       allocate(dc%neum(2,7))
    elseif (allocated(dc%ddl).and.n/=dc%nall) then
       deallocate(dc%ddl,dc%ddr)
       deallocate(dc%inv_lhs)
       deallocate(dc%eigen,dc%vect,dc%inv_vect)
       dc%nall=n 
       dc%n=n-2        
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,8))
       allocate(dc%inv_lhs(n-2,n-2))
       allocate(dc%eigen(n-2),dc%vect(n-2,n-2),dc%inv_vect(n-2,n-2))
       allocate(dc%neum(2,7))
    endif

    
    dc%bc(1:2)=bct(1:2)

    !-> compute coefficients for second derivatives
    
    if (present(out_name)) open(333,file='sd3d_'//out_name)
    dc%ddl=0._rk ; dc%ddr=0._rk
    call dder_coeffs(dc%n,dc%ddl,dc%ddr,grid,bct)

!    call print(dc%n,dc%ddl,dc%ddr)

    !-> compute solver matrix

    call solver_matrix_init(dc%n,dc%ddl,dc%ddr,dc%eigen,dc%vect,dc%inv_vect,&
         dc%inv_lhs,dc%bc)


    !-> compute coefficients for neumann interpolation
    call extrapolation_coeffs(dc%n,dc%neum,grid)

  end subroutine solver_coeffs_init

  subroutine solver_matrix_init(n,ddl,ddr,wr,vr,invvr,invlhs,bc)
! -----------------------------------------------------------------------
! solver3d : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    integer(ik),intent(in) :: n,bc(2)
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,8)
    real(rk) :: lhs(n,n),rhs(n,n)
    real(rk) :: invlhs(n,n),temp(n,n),temp1(n,n)
    integer :: ipiv(n,n),info ! lapack aux
    real(rk) :: work(n,n) ! lapack aux
    real(rk) :: wr(n),wi(n),vr(n,n),vl,invvr(n,n)

    lhs(:,:)=0._rk ; rhs(:,:)=0._rk

    !-> construct lhs

    lhs(1,1)=ddl(1,3) ; lhs(1,2)=ddl(1,4) ; lhs(1,3)=ddl(1,5) ; lhs(1,4)=ddl(1,1)

    lhs(2,1)=ddl(2,2) ; lhs(2,2)=ddl(2,3) ; lhs(2,3)=ddl(2,4) ; lhs(2,4)=ddl(2,5)

    do i=3,n-2
       lhs(i,i-2)=ddl(i,1) ; lhs(i,i-1)=ddl(i,2) ; lhs(i,i)=ddl(i,3) ; 
       lhs(i,i+1)=ddl(i,4) ; lhs(i,i+2)=ddl(i,5) 
    enddo

    lhs(n-1,n-3)=ddl(n-1,1) ; lhs(n-1,n-2)=ddl(n-1,2) ; lhs(n-1,n-1)=ddl(n-1,3)
    lhs(n-1,n)=ddl(n-1,4)

    lhs(n,n-3)=ddl(n,5) ; lhs(n,n-2)=ddl(n,1) ; lhs(n,n-1)=ddl(n,2)
    lhs(n,n)=ddl(n,3)

    !-> construct rhs

    rhs(1,1)=ddr(1,2) ; rhs(1,2)=ddr(1,3) ; rhs(1,3)=ddr(1,4) ; rhs(1,4)=ddr(1,5)
    rhs(1,5)=ddr(1,6) ; rhs(1,6)=ddr(1,7)

    rhs(2,1)=ddr(2,2) ; rhs(2,2)=ddr(2,3) ; rhs(2,3)=ddr(2,4) ; rhs(2,4)=ddr(2,5)
    rhs(2,5)=ddr(2,6) ; rhs(2,6)=ddr(2,7) ; rhs(2,7)=ddr(2,8) 

    do i=3,n-2
       rhs(i,i-2)=ddr(i,1) ; rhs(i,i-1)=ddr(i,2) ; rhs(i,i)=ddr(i,3) ; 
       rhs(i,i+1)=ddr(i,4) ; rhs(i,i+2)=ddr(i,5) 
    enddo

    rhs(n-1,n-6)=ddr(n-1,1)
    rhs(n-1,n-5)=ddr(n-1,2) ; rhs(n-1,n-4)=ddr(n-1,3) ; rhs(n-1,n-3)=ddr(n-1,4) 
    rhs(n-1,n-2)=ddr(n-1,5) ; rhs(n-1,n-1)=ddr(n-1,6) ; rhs(n-1,n)=ddr(n-1,7)

    rhs(n,n-5)=ddr(n,1) ; rhs(n,n-4)=ddr(n,2) ; rhs(n,n-3)=ddr(n,3)
    rhs(n,n-2)=ddr(n,4) ; rhs(n,n-1)=ddr(n,5) ; rhs(n,n)=ddr(n,6)
    
    !-> inverse lhs
    invlhs=lhs
    ! LU factorization
    call dgetrf(n,n,invlhs,n,ipiv, info)
    ! compute inverse
    call dgetri(n,invlhs,n,ipiv,work,n,info)

    !-> compute inverse-lhs*rhs

    temp=matmul(invlhs,rhs)

    !-> compute eigenvalue and eigenvector
    call dgeev('N','V',n,temp,n,wr,wi,vl,1,vr,n,work,n*n,info)

!    print*,'dgeev info=',info
!    print*,'eigenvalues : real imag'
!    do i=1,n
!       print'(100es17.8)',wr(i),wi(i)
!    enddo

    !-> For neumann bc identify position of zero valued eigenvalues
!    print*,bc
!    if (all(bc==2)) then
!       print*,'toto'
!    endif

!    print*,'eigenvectors'
!    do i=1,n
!       print'(100es17.8)',vr(i,:)
!    enddo
    

    !-> compute inverse of eigenvectors
    invvr=vr
    ! LU factorization
    call dgetrf(n,n,invvr,n,ipiv, info)
    ! compute inverse
    call dgetri(n,invvr,n,ipiv,work,n,info)

  end subroutine solver_matrix_init

  subroutine dder_coeffs(n,ddl,ddr,grid,bct)
! -----------------------------------------------------------------------
! solver3d : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    !$ use OMP_LIB
    implicit none
    integer(ik),intent(in) :: n
    integer(ik),intent(in) :: bct(2)
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,8)
    real(rk),intent(in) :: grid(0:n+1)

    ! alpha,beta : central coefficients lhs
    ! a,b : central coefficients rhs
    real(rk) :: alp,bet,gam,del
    real(rk) :: a,b,c,d,e
    ! alp1,bet1,gam1,del1 : 1st points coefficients lhs
    ! alp2,bet2,gam2,del2 : 2st points coefficients lhs
    real(rk) :: alp1,bet1,gam1,del1
    real(rk) :: alp2,bet2,gam2,del2
    ! b1,c1,d1,e1 : 1st points coefficients rhs
    ! b2,c2,d2,e2 : 2st points coefficients lhs
    real(rk) :: a1,b1,c1,d1,e1,f1,g1
    real(rk) :: a2,b2,c2,d2,e2,f2,g2,hh2
    ! h1,h2,h3,h4 : grid intervals
    real(rk) :: h1,h2,h3,h4,h5,h6

    ddl=0._rk ; ddr=0._rk

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6,bct(1))
    ddl(1,1)=del1 ; ddl(1,2)=0._rk ; ddl(1,3)=alp1 ; ddl(1,4)=bet1 ; ddl(1,5)=gam1
    ddr(1,1)=a1 ; ddr(1,2)=b1 ; ddr(1,3)=c1 ; ddr(1,4)=d1 ; ddr(1,5)=e1
    ddr(1,6)=f1 ; ddr(1,7)=g1

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    h5=grid(6)-grid(5)
    h6=grid(7)-grid(6)
    call dder_coeffs_i3(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6,bct(1))
    ddl(2,1)=0._rk ; ddl(2,2)=alp2 ; ddl(2,3)=bet2 ; ddl(2,4)=gam2 ; ddl(2,5)=del2
    ddr(2,1)=a2 ; ddr(2,2)=b2 ; ddr(2,3)=c2 ; ddr(2,4)=d2 ; ddr(2,5)=e2
    ddr(2,6)=f2 ; ddr(2,7)=g2 !; ddr(2,8)=hh2 ;

!$OMP PARALLEL PRIVATE(i,h1,h2,h3,h4,a,b,c,d,e,alp,bet,gam,del)
!$OMP DO SCHEDULE(RUNTIME)
    do i=3,n-2
       h1=grid(i-1)-grid(i-2)
       h2=grid(i)-grid(i-1)
       h3=grid(i+1)-grid(i)
       h4=grid(i+2)-grid(i+1)
       call dder_coeffs_c(a,b,c,d,e,alp,bet,gam,del,h1,h2,h3,h4)
       ddl(i,1)=alp ; ddl(i,2)=bet ; ddl(i,3)=1._rk ; ddl(i,4)=gam ; ddl(i,5)=del
       ddr(i,1)=a ; ddr(i,2)=b ; ddr(i,3)=c ; ddr(i,4)=d ; ddr(i,5)=e
    enddo
!$OMP END DO
!$OMP END PARALLEL

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call dder_coeffs_i3(a2,b2,c2,d2,e2,f2,g2,alp2,bet2,gam2,del2,h1,h2,h3,h4,h5,h6,bct(2))
    ddl(n-1,1)=del2 ; ddl(n-1,2)=gam2 ; ddl(n-1,3)=bet2 ; ddl(n-1,4)=alp2 ; ddl(n-1,5)=0._rk
    !ddr(n-1,1)=hh2 
    ddr(n-1,2)=g2 ; ddr(n-1,3)=f2 
    ddr(n-1,4)=e2 ; ddr(n-1,5)=d2 ; ddr(n-1,6)=c2 ; ddr(n-1,7)=b2 ; ddr(n-1,8)=a2
    if (bct(2)==2) ddr(n-1,8)=-a2


    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    h5=grid(n-4)-grid(n-5)
    h6=grid(n-5)-grid(n-6)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,f1,g1,alp1,bet1,gam1,del1,h1,h2,h3,h4,h5,h6,bct(2))
    ddl(n,1)=gam1 ; ddl(n,2)=bet1 ; ddl(n,3)=alp1 ; ddl(n,4)=0._rk ; ddl(n,5)=del1
    ddr(n,1)=g1 ; ddr(n,2)=f1 
    ddr(n,3)=e1 ; ddr(n,4)=d1 ; ddr(n,5)=c1 ; ddr(n,6)=b1 ; ddr(n,7)=a1
    if (bct(2)==2) ddr(n,7)=-a1

  end subroutine dder_coeffs


  subroutine extrapolation_coeffs(n,ddr,grid)
! -----------------------------------------------------------------------
! solver3d : computation of coefficients for neumann extrapolation
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddr(2,7)
    real(rk),intent(in) :: grid(0:n+1)

    ! alp1 : 1st points coefficients lhs
    real(rk) :: alp1
    ! a1,b1,c1,d1,e1,f1,g1 : 1st points coefficients rhs
    real(rk) :: a1,b1,c1,d1,e1,f1,g1
    ! h1,h2,h3,h4 : grid intervals
    real(rk) :: h1,h2,h3,h4,h5,h6

    ddr=0._rk

    h1=grid(1)-grid(0)
    h2=grid(2)-grid(1)
    h3=grid(3)-grid(2)
    h4=grid(4)-grid(3)
    h5=grid(5)-grid(4)
    h6=grid(6)-grid(5)
    call der_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,h1,h2,h3,h4,h5,h6)
    ddr(1,1)=a1 ; ddr(1,2)=b1 ; ddr(1,3)=c1 ; ddr(1,4)=d1 ; ddr(1,5)=e1
    ddr(1,6)=f1 ; ddr(1,7)=g1
    
    h1=grid(n+1)-grid(n)
    h2=grid(n)-grid(n-1)
    h3=grid(n-1)-grid(n-2)
    h4=grid(n-2)-grid(n-3)
    h5=grid(n-3)-grid(n-4)
    h6=grid(n-4)-grid(n-5)
    call der_coeffs_i1(a1,b1,c1,d1,e1,f1,g1,alp1,h1,h2,h3,h4,h5,h6)
    ddr(2,1)=-g1 ; ddr(2,2)=-f1 
    ddr(2,3)=-e1 ; ddr(2,4)=-d1 ; ddr(2,5)=-c1 ; ddr(2,6)=-b1 ; ddr(2,7)=-a1
    
  end subroutine extrapolation_coeffs

  subroutine extrapolation_neum_x(dr,sol,bc,n1,n2,n3,bcl)
! -----------------------------------------------------------------------
! solver3d : extrapolation routine for neumann boundary condition
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    integer(ik),intent(in) :: n1,n2,n3,bcl
    real(rk),intent(inout) :: sol(0:n1+1,0:n2+1,0:n3+1)
    real(rk),intent(in) :: bc(n2,n3,2)
    real(rk),intent(in) :: dr(2,7)
    real(rk) :: aux
    integer(ik) :: i,j,k

    !-> first point
    if (bcl==1) then
       aux=1._rk/dr(1,1)
       do k=1,n3
          do j=1,n2
             sol(0,j,k)=(bc(j,k,1)-dr(1,2)*sol(1,j,k)-dr(1,3)*sol(2,j,k)-dr(1,4)*sol(3,j,k)&
                  -dr(1,5)*sol(4,j,k)-dr(1,6)*sol(5,j,k)-dr(1,7)*sol(6,j,k))*aux
          enddo
       enddo
    endif

    !-> last point
    if (bcl==2) then
       aux=1._rk/dr(2,7)
       do k=1,n3
          do j=1,n2
             sol(n1+1,j,k)=(bc(j,k,2)-dr(2,6)*sol(n1,j,k)-dr(2,5)*sol(n1-1,j,k)&
                  -dr(2,4)*sol(n1-2,j,k)&
                  -dr(2,3)*sol(n1-3,j,k)-dr(2,2)*sol(n1-4,j,k)-dr(2,1)*sol(n1-5,j,k))*aux
          enddo
       enddo
    endif

  end subroutine extrapolation_neum_x

  subroutine extrapolation_neum_y(dr,sol,bc,n1,n2,n3,bcl)
! -----------------------------------------------------------------------
! solver3d : extrapolation routine for neumann boundary condition
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    integer(ik),intent(in) :: n1,n2,n3,bcl
    real(rk),intent(inout) :: sol(0:n1+1,0:n2+1,0:n3+1)
    real(rk),intent(in) :: bc(n1,n3,2)
    real(rk),intent(in) :: dr(2,7)
    real(rk) :: aux
    integer(ik) :: i,j,k

    !-> first point
    if (bcl==1) then
       aux=1._rk/dr(1,1)
       do k=1,n3
          do i=1,n1
             sol(i,0,k)=(bc(i,k,1)-dr(1,2)*sol(i,1,k)-dr(1,3)*sol(i,2,k)-dr(1,4)*sol(i,3,k)&
                  -dr(1,5)*sol(i,4,k)-dr(1,6)*sol(i,5,k)-dr(1,7)*sol(i,6,k))*aux
          enddo
       enddo
    endif

    !-> last point
    if (bcl==2) then
       aux=1._rk/dr(2,7)
       do k=1,n3
          do i=1,n1
             sol(i,n2+1,k)=(bc(i,k,2)-dr(2,6)*sol(i,n2,k)-dr(2,5)*sol(i,n2-1,k)&
                  -dr(2,4)*sol(i,n2-2,k)&
                  -dr(2,3)*sol(i,n2-3,k)-dr(2,2)*sol(i,n2-4,k)-dr(2,1)*sol(i,n2-5,k))*aux
          enddo
       enddo
    endif

  end subroutine extrapolation_neum_y

  subroutine extrapolation_neum_z(dr,sol,bc,n1,n2,n3,bcl)
! -----------------------------------------------------------------------
! solver3d : extrapolation routine for neumann boundary condition
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 10/2012
!
    implicit none
    integer(ik),intent(in) :: n1,n2,n3,bcl
    real(rk),intent(inout) :: sol(0:n1+1,0:n2+1,0:n3+1)
    real(rk),intent(in) :: bc(n1,n2,2)
    real(rk),intent(in) :: dr(2,7)
    real(rk) :: aux
    integer(ik) :: i,j,k

    !-> first point
    if (bcl==1) then
       aux=1._rk/dr(1,1)
       do j=1,n2
          do i=1,n1
             sol(i,j,0)=(bc(i,j,1)-dr(1,2)*sol(i,j,1)-dr(1,3)*sol(i,j,2)-dr(1,4)*sol(i,j,3)&
                  -dr(1,5)*sol(i,j,4)-dr(1,6)*sol(i,j,5)-dr(1,7)*sol(i,j,6))*aux
          enddo
       enddo
    endif

    !-> last point
    if (bcl==2) then
       aux=1._rk/dr(2,7)
       do j=1,n2
          do i=1,n1
             sol(i,j,n3+1)=(bc(i,j,2)-dr(2,6)*sol(i,j,n3)-dr(2,5)*sol(i,j,n3-1)&
                  -dr(2,4)*sol(i,j,n3-2)&
                  -dr(2,3)*sol(i,j,n3-3)-dr(2,2)*sol(i,j,n3-4)-dr(2,1)*sol(i,j,n3-5))*aux
          enddo
       enddo
    endif

  end subroutine extrapolation_neum_z

end module class_solver_3d
