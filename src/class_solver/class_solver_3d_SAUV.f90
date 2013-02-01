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
  use class_solver_coefficient_1d
  implicit none

  !-> Define type containing solver coefficients for one direction
  type,public :: solver_coeffs
     private
     integer(ik) :: nall,n
     !-> second derivatives coefficients : lhs and rhs
     real(rk),allocatable :: ddl(:,:),ddr(:,:)
     !-> lhs : inverse matrix 
     real(rk),allocatable :: inv_lhs(:,:)
     !-> lhs : eigenvalues, eigenvectors, inverse eigenvectors
     real(rk),allocatable :: eigen(:),vect(:,:),inv_vect(:,:)
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
    real(rk),intent(in) :: ddl(n,5),ddr(n,5)
    integer(ik) :: i,j

    print*,'ddl'
    do i=1,n
       print'(5es17.8)',(ddl(i,j),j=1,5)
    enddo
    print*,'ddr'
    do i=1,n
       print'(7es17.8)',(ddr(i,j),j=1,5)
    enddo

  end subroutine print

  subroutine solve_3d(sc,rhs,solf,n1,n2,n3)
! -----------------------------------------------------------------------
! solver3d : solve system A*x=B with matrix diagonalization 
!            in the three directions
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
    type(solver_coeffs_3d),intent(in) :: sc
    integer(ik),intent(in) :: n1,n2,n3
    real(rk),intent(inout) :: rhs(n1,n2,n3)
    real(rk),intent(inout) :: solf(0:n1+1,0:n2+1,0:n3+1)
    real(rk) :: sol(n1,n2,n3)
    real(rk) :: aux1(n1,n2,n3),aux2(n2,n2),aux3(n3,n3)
    real(rk) :: clx(n1),cly(n2),clz(n3)
    integer(ik) :: i,j,k

    !-> apply boundary conditions to rhs
    print*,solf(:,:,0)
    print*,'sumrhs',sum(rhs)
    !-> x-direction
    do k=1,n3
       do j=1,n2
          clx(1)=sc%cx%ddr(1,1)*solf(0,j,k)
          clx(2)=sc%cx%ddr(2,1)*solf(0,j,k)
          clx(4:n1-3)=0._rk
          clx(n1-1)=sc%cx%ddr(n1-1,5)*solf(n1+1,j,k)
          clx(n1)=sc%cx%ddr(n1,5)*solf(n1+1,j,k)
          rhs(:,j,k)=rhs(:,j,k)-matmul(sc%cx%inv_lhs,clx)
       enddo
    enddo

    !-> y-direction
    do k=1,n3
       do i=1,n1
          cly(1)=sc%cy%ddr(1,1)*solf(i,0,k)
          cly(2)=sc%cy%ddr(2,1)*solf(i,0,k)
          cly(4:n1-3)=0._rk
          cly(n2-1)=sc%cy%ddr(n2-1,5)*solf(i,n2+1,k)
          cly(n2)=sc%cy%ddr(n2,5)*solf(i,n2+1,k)
          rhs(i,:,k)=rhs(i,:,k)-matmul(sc%cy%inv_lhs,cly)
       enddo
    enddo
    
    !-> z-direction
    do j=1,n2
       do i=1,n1
          clz(1)=sc%cz%ddr(1,1)*solf(i,j,0)
          clz(2)=sc%cz%ddr(2,1)*solf(i,j,0)
          clz(4:n1-3)=0._rk
          clz(n3-1)=sc%cz%ddr(n3-1,5)*solf(i,j,n3+1)
          clz(n3)=sc%cz%ddr(n3,5)*solf(i,j,n3+1)
          rhs(i,j,:)=rhs(i,j,:)-matmul(sc%cz%inv_lhs,clz)
       enddo
    enddo
    print*,'sumrhs',sum(rhs)
    

    !-> compute modified rhs ------------------------------------

    !->  x-direction
!    call dgemm('n','n',n1,n2*n3,n1,1.d0,sc%cx%inv_vect,n1,rhs,n1,0.d0,aux1,n1)
!    print*,sum(aux1)/(n1*n2*n3)
    do k=1,n3
       aux1(:,:,k)=matmul(sc%cx%inv_vect,rhs(:,:,k))
    enddo
    print*,sum(aux1)/(n1*n2*n3)

    !-> y-direction
    aux2=transpose(sc%cy%inv_vect)
    do k=1,n3
       aux1(:,:,k)=matmul(aux1(:,:,k),aux2)
    enddo
    print*,sum(aux1)/(n1*n2*n3)

    !-> z-direction
    aux3=transpose(sc%cz%inv_vect)
    do i=1,n1
       aux1(i,:,:)=matmul(aux1(i,:,:),aux3)
    enddo
    print*,sum(aux1)/(n1*n2*n3)
    
    !-> compute solution with eigenvalues -------------------------
    do k=1,n3
       do j=1,n2
          do i=1,n1
             sol(i,j,k)=aux1(i,j,k)/(sc%cx%eigen(i)+sc%cy%eigen(j)+sc%cz%eigen(k))
!             sol(i,j,k)=aux1(i,j,k)/(sc%cx%eigen(i))
          enddo
       enddo
    enddo
    print*,sum(sol)/(n1*n2*n3)

    !-> compute solution -------------------------------------------
    
    !->  x-direction
    do k=1,n3
       sol(:,:,k)=matmul(sc%cx%vect,sol(:,:,k))
    enddo
    print*,sum(sol)/(n1*n2*n3)
    
    !-> y-direction
    aux2=transpose(sc%cy%vect)
    do k=1,n3
       sol(:,:,k)=matmul(sol(:,:,k),aux2)
    enddo
    print*,sum(sol)/(n1*n2*n3)

    !-> z-direction
    aux3=transpose(sc%cz%vect)
    do i=1,n1
       sol(i,:,:)=matmul(sol(i,:,:),aux3)
    enddo
    print*,sum(sol)/(n1*n2*n3)

    solf(1:n1,1:n2,1:n3)=sol(1:n1,1:n2,1:n3)
    

  end subroutine solve_3d

  subroutine solver_coeffs_init_3d(grid1,grid2,grid3,sc,n1,n2,n3,out_name)
! -----------------------------------------------------------------------
! solver3d : Initialisation of all 3d solver coefficients
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    type(solver_coeffs_3d),intent(out) :: sc
    integer(ik),intent(in) :: n1,n2,n3
    real(rk),intent(in) :: grid1(n1),grid2(n2),grid3(n3)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    
    print'(a)','Solver : x-direction -------------------------------------------'
    call solver_coeffs_init(grid1,sc%cx,n1,out_name)
    print'(a)','Solver : y-direction -------------------------------------------'
    call solver_coeffs_init(grid2,sc%cy,n2,out_name)
    print'(a)','Solver : z-direction -------------------------------------------'
    call solver_coeffs_init(grid3,sc%cz,n3,out_name)


  end subroutine solver_coeffs_init_3d

  subroutine solver_coeffs_init(grid,dc,n,out_name)
! -----------------------------------------------------------------------
! solver3d : Initialisation of solver coefficients in one direction
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    type(solver_coeffs),intent(out) :: dc
    integer(ik),intent(in) :: n
    real(rk),intent(in) :: grid(n)
    character(len=*),optional :: out_name
    character(len=500) :: error_msg
    
    !-> initialize and allocate dc

    if (.not.allocated(dc%ddl)) then
       dc%nall=n
       dc%n=n-2 
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,5))
       allocate(dc%inv_lhs(n-2,n-2))
       allocate(dc%eigen(n-2),dc%vect(n-2,n-2),dc%inv_vect(n-2,n-2))
    elseif (allocated(dc%ddl).and.n/=dc%nall) then
       deallocate(dc%ddl,dc%ddr)
       deallocate(dc%inv_lhs)
       deallocate(dc%eigen,dc%vect,dc%inv_vect)
       dc%nall=n 
       dc%n=n-2        
       allocate(dc%ddl(n-2,5),dc%ddr(n-2,5))
       allocate(dc%inv_lhs(n-2,n-2))
       allocate(dc%eigen(n-2),dc%vect(n-2,n-2),dc%inv_vect(n-2,n-2))
    endif

    !-> compute coefficients for second derivatives
    
    if (present(out_name)) open(333,file='sd3d_'//out_name)
    dc%ddl=0._rk ; dc%ddr=0._rk
    call dder_coeffs(dc%n,dc%ddl,dc%ddr,grid)

    call print(dc%n,dc%ddl,dc%ddr)

    !-> compute solver matrix

    call solver_matrix_init(dc%n,dc%ddl,dc%ddr,dc%eigen,dc%vect,dc%inv_vect,dc%inv_lhs)

  end subroutine solver_coeffs_init

  subroutine solver_matrix_init(n,ddl,ddr,wr,vr,invvr,invlhs)
! -----------------------------------------------------------------------
! solver3d : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,5)
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

    print*,'lhs'
    do i=1,n
       print'(100es17.8)',lhs(i,:)
    enddo

    !-> construct rhs

    rhs(1,1)=ddr(1,2) ; rhs(1,2)=ddr(1,3) ; rhs(1,3)=ddr(1,4) ; rhs(1,4)=ddr(1,5)

    rhs(2,1)=ddr(2,2) ; rhs(2,2)=ddr(2,3) ; rhs(2,3)=ddr(2,4) ; rhs(2,4)=ddr(2,5)

    do i=3,n-2
       rhs(i,i-2)=ddr(i,1) ; rhs(i,i-1)=ddr(i,2) ; rhs(i,i)=ddr(i,3) ; 
       rhs(i,i+1)=ddr(i,4) ; rhs(i,i+2)=ddr(i,5) 
    enddo

    rhs(n-1,n-3)=ddr(n-1,1) ; rhs(n-1,n-2)=ddr(n-1,2) ; rhs(n-1,n-1)=ddr(n-1,3)
    rhs(n-1,n)=ddr(n-1,4)

    rhs(n,n-3)=ddr(n,1) ; rhs(n,n-2)=ddr(n,2) ; rhs(n,n-1)=ddr(n,3)
    rhs(n,n)=ddr(n,4)

    print*,'rhs'
    do i=1,n
       print'(100es17.8)',rhs(i,:)
    enddo
    
    !-> inverse lhs
    invlhs=lhs
    ! LU factorization
    call dgetrf(n,n,invlhs,n,ipiv, info)
    print*,'factorization info=',info
    ! compute inverse
    call dgetri(n,invlhs,n,ipiv,work,n,info)
    print*,'inversion info=',info
    print*,'invlhs'
    do i=1,n
       print'(100es17.8)',invlhs(i,:)
    enddo

    work=matmul(lhs,invlhs)
    print*,'check'
    do i=1,n
       print'(100es17.8)',work(i,:)
    enddo

    !-> compute inverse-lhs*rhs

    temp=matmul(invlhs,rhs)
    print*,'invlhs*rhs'
    do i=1,n
       print'(100es17.8)',temp(i,:)
    enddo

    !-> compute eigenvalue and eigenvector
    temp1=temp
    call dgeev('N','V',n,temp,n,wr,wi,vl,1,vr,n,work,n*n,info)

    print*,'dgeev info=',info
    print*,'eigenvalues : real imag'
    do i=1,n
       print'(100es17.8)',wr(i),wi(i)
    enddo
    print*,'eigenvectors'
    do i=1,n
       print'(100es17.8)',vr(i,:)
    enddo
    
    temp=matmul(temp1,vr)
    print*,'check eigen'
    do i=1,n
       print'(100es17.8)',temp(i,:)
    enddo
    temp=0._rk
    do i=1,n
       temp(i,i)=wr(i)
    enddo
    temp=matmul(vr,temp)
    print*,'check eigen'
    do i=1,n
       print'(100es17.8)',temp(i,:)
    enddo
    

    !-> compute inverse of eigenvectors
    invvr=vr
    ! LU factorization
    call dgetrf(n,n,invvr,n,ipiv, info)
    print*,'factorization info=',info
    ! compute inverse
    call dgetri(n,invvr,n,ipiv,work,n,info)
    print*,'inversion info=',info
    print*,'invvr'
    do i=1,n
       print'(100es17.8)',invvr(i,:)
    enddo

    work=matmul(vr,invvr)
    print*,'check'
    do i=1,n
       print'(100es17.8)',work(i,:)
    enddo

  end subroutine solver_matrix_init

  subroutine dder_coeffs(n,ddl,ddr,grid)
! -----------------------------------------------------------------------
! solver3d : computation of coefficients for second derivatives
! -----------------------------------------------------------------------
! Matthieu Marquillie
! 03/2012
!
    implicit none
    integer(ik),intent(in) :: n
    integer(ik) :: i
    real(rk),intent(out) :: ddl(n,5),ddr(n,5)
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
    real(rk) :: a1,b1,c1,d1,e1
    real(rk) :: a2,b2,c2,d2,e2
    ! h1,h2,h3,h4 : grid intervals
    real(rk) :: h1,h2,h3,h4

    ddl=0._rk ; ddr=0._rk

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,alp1,bet1,gam1,del1,h1,h2,h3,h4)
    ddl(1,1)=del1 ; ddl(1,2)=0._rk ; ddl(1,3)=alp1 ; ddl(1,4)=bet1 ; ddl(1,5)=gam1
    ddr(1,1)=a1 ; ddr(1,2)=b1 ; ddr(1,3)=c1 ; ddr(1,4)=d1 ; ddr(1,5)=e1

    h1=grid(2)-grid(1)
    h2=grid(3)-grid(2)
    h3=grid(4)-grid(3)
    h4=grid(5)-grid(4)
    call dder_coeffs_i3(a2,b2,c2,d2,e2,alp2,bet2,gam2,del2,h1,h2,h3,h4)
    ddl(2,1)=0._rk ; ddl(2,2)=alp2 ; ddl(2,3)=bet2 ; ddl(2,4)=gam2 ; ddl(2,5)=del2
    ddr(2,1)=a2 ; ddr(2,2)=b2 ; ddr(2,3)=c2 ; ddr(2,4)=d2 ; ddr(2,5)=e2

    do i=3,n-2
       h1=grid(i-1)-grid(i-2)
       h2=grid(i)-grid(i-1)
       h3=grid(i+1)-grid(i)
       h4=grid(i+2)-grid(i+1)
       call dder_coeffs_c(a,b,c,d,e,alp,bet,gam,del,h1,h2,h3,h4)
       ddl(i,1)=alp ; ddl(i,2)=bet ; ddl(i,3)=1._rk ; ddl(i,4)=gam ; ddl(i,5)=del
       ddr(i,1)=a ; ddr(i,2)=b ; ddr(i,3)=c ; ddr(i,4)=d ; ddr(i,5)=e
    enddo

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    call dder_coeffs_i3(a2,b2,c2,d2,e2,alp2,bet2,gam2,del2,h1,h2,h3,h4)
    ddl(n-1,1)=del2 ; ddl(n-1,2)=gam2 ; ddl(n-1,3)=bet2 ; ddl(n-1,4)=alp2 ; ddl(n-1,5)=0._rk
    ddr(n-1,1)=e2 ; ddr(n-1,2)=d2 ; ddr(n-1,3)=c2 ; ddr(n-1,4)=b2 ; ddr(n-1,5)=a2

    h1=grid(n)-grid(n-1)
    h2=grid(n-1)-grid(n-2)
    h3=grid(n-2)-grid(n-3)
    h4=grid(n-3)-grid(n-4)
    call dder_coeffs_i2(a1,b1,c1,d1,e1,alp1,bet1,gam1,del1,h1,h2,h3,h4)
    ddl(n,1)=gam1 ; ddl(n,2)=bet1 ; ddl(n,3)=alp1 ; ddl(n,4)=0._rk ; ddl(n,5)=del1
    ddr(n,1)=e1 ; ddr(n,2)=d1 ; ddr(n,3)=c1 ; ddr(n,4)=b1 ; ddr(n,5)=a1

  end subroutine dder_coeffs


end module class_solver_3d
