#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module solver_mod
  use kinds, only : real_kind, int_kind
  use dimensions_mod, only : npsq, nlev
  use perf_mod, only: t_startf, t_stopf ! _EXTERNAL
  use parallel_mod, only : abortmp

  use ,intrinsic :: iso_c_binding


  implicit none
  private

  character(len=8), private, parameter :: blkjac_storage = "inverse"
  !  character(len=8), private, parameter :: blkjac_storage = "LUfactor"

  type, public :: blkjac_t
     real (kind=real_kind), dimension(npsq,npsq,nlev) :: E
     integer(kind=int_kind),     dimension(npsq,nlev) :: ipvt
  end type blkjac_t


  public  :: solver_test

#ifdef TRILINOS
  public  :: solver_test_ml
#endif


contains


#ifdef TRILINOS

  ! ================================================
  ! helm_graph:
  !
  ! ================================================
  subroutine helm_graph(N,global_index_mat,op_data_ptr) bind(C,name='homme_globalIDs')
  
    use dimensions_mod, only : nlev, np,npsq
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in)      :: N
    integer(c_int), intent(inout)      :: global_index_mat(N)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()

    integer      :: temp_global_index_mat(N)


    integer  :: nets,nete

    ! ===========
    ! Local
    ! ===========
    integer :: ie
    integer :: i,j,k
    integer gindex

    
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer


    nets=fptr%nets
    nete=fptr%nete


    gindex=1
       do ie=nets,nete
          do k=1,nlev
             do j=1,np
                do i=1,np
                  global_index_mat(gindex)=fptr%base(ie)%gdofP(i,j) 
                  gindex=gindex+1;
                enddo
             enddo
          end do
       end do
  end subroutine helm_graph





  ! ================================================
  ! helm_mat:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! return helmholtz matrix on element El
  !  form by applying to a unit vector x for each column:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !   HelmOp*x=  D QQ^t (M + a*L ) x 
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  !
  ! ================================================
  subroutine helm_mat(ie,ElMatSize,ElementMatOut,Indices,op_data_ptr) bind(C,name='helm_mat')
  
    use dimensions_mod, only : nlev, np,npsq
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in),value :: ie 
    integer(c_int), intent(in),value :: ElMatSize
    real(c_double), intent(out)      :: ElementMatOut(ElMatSize,ElMatSize)
    integer(c_int), intent(out)      :: Indices(ElMatSize)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()

    integer  :: nets,nete
!    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)

    ! ===========
    ! Local
    ! ===========
    real (kind=real_kind), allocatable :: ElMatCol(:,:)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: i,j,k
    integer ColIndex
    
    integer colindices(np,np)
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer

    ElementMatOut=0.0d0

    allocate(ElMatCol(np,np))


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Application Loop
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !do ColIndex=1,nlev*np**2
    do ColIndex=1,np**2
      do j=1,np
        do i=1,np
           x(i,j) = 0.0d0
           if(ColIndex == ((j-1)*np +i)) then 
                   x(i,j)=1 
           endif
        enddo
      enddo
             
             ! Apply x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;

          !ElMatCol(:,:)= fptr%base(ie)%rspheremp(:,:)*fptr%base(ie)%spheremp(:,:)*x(:,:) 
          !ElMatCol(:,:)= fptr%base(ie)%spheremp(:,:)*x(:,:) 
          !ElMatCol(:,:)= x(:,:) 
          ElMatCol(:,:)= fptr%base(ie)%rspheremp(:,:)*(fptr%base(ie)%spheremp(:,:)*x(:,:) +&
            alambda*laplace_sphere_wk(x,fptr%deriv,fptr%base(ie),var_coef=.false.))
           



          ElementMatOut(:,ColIndex)=reshape(ElMatCol,(/np*np/))
           
!          ElementMatOut(1+(1-k)*np**2:1+k*np**2,ColIndex)=reshape(ElMatCol,(/np*np/))
    enddo

    k=1
    do j=1,np
      do i=1,np
         Indices(k)=fptr%base(ie)%gdofP(i,j)
         !Indices(k)=k
         k=k+1
      enddo
    enddo




!       write(6,*)'vecOut=[',vecOut,']'


!    if (hybrid%masterthread) print *,'applied helmholtz'
!    call flush(6)
       deallocate(ElMatCol)
  end subroutine helm_mat


  ! ================================================
  ! helm_rhs:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! return helmholtz matrix on element El
  !  form by applying to a unit vector x for each column:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !   HelmOp*x=  D QQ^t (M + a*L ) x 
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  !
  ! ================================================
  subroutine helm_map(ie,ElMatSize,Indices,op_data_ptr) bind(C,name='helm_map')
  
 !~ elem,edge1,red,hybrid,deriv,nets,nete,vecIn)
    use dimensions_mod, only : nlev, np,npsq
!    use element_mod, only : element_t
!    use edge_mod, only : edgebuffer_t, edgevpack, edgevunpack!, edgerotate
!    use derivative_mod, only : derivative_t, laplace_sphere_wk
!    use control_mod, only : maxits, while_iter, tol, precon_method

!    use physical_constants, only : rrearth, dd_pi, rearth, omega
!    use bndry_mod, only : bndry_exchangeV
!    use linear_algebra_mod, only : matvec
!    use parallel_mod, only : haltmp
!    use hybrid_mod, only : hybrid_t
    use derived_type_mod, only :derived_type
    implicit none
    integer(c_int), intent(in),value :: ie 
    integer(c_int), intent(in),value :: ElMatSize
    !real(c_double), intent(out)      :: ElementMatOut(ElMatSize,ElMatSize)
    integer(c_int), intent(out)      :: Indices(ElMatSize)
    type(c_ptr)                      :: op_data_ptr
    type(derived_type),  pointer          :: fptr=>NULL()


    ! ===========
    ! Local
    ! ===========

    integer :: i,j,k
    
    call c_f_pointer(op_data_ptr,fptr) !convert c pointer to f pointer

     k=1
      do j=1,np
        do i=1,np
           Indices(k)=fptr%base(ie)%gdofP(i,j)
           k=k+1
        enddo
      enddo


  end subroutine helm_map




  ! ================================================
  ! solver_test_ml:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! solve for x:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D solve - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !     D QQ^t (M + a*L ) x = rhs
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  ! Note: if we solve laplace equation instead of Helmholz, we need to
  ! ensure < rhs,1>=0
  !
  ! ================================================
  subroutine solver_test_ml(elem,edge1,red,hybrid,deriv,nets,nete)
    use dimensions_mod, only : nlev, np,npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad, cg_create
    use edge_mod, only : edgevpack, edgevunpack!, edgerotate
    use edgetype_mod, only : edgebuffer_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use control_mod, only : maxits, tol, precon_method
    use physical_constants, only : rrearth, dd_pi, rearth, omega
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp
    use hybrid_mod, only : hybrid_t
    use global_norms_mod, only : linf_snorm, l2_snorm
    use derived_type_mod, only : derived_type, initialize
    use time_mod, only : TimeLevel_t


    interface

     subroutine belosfinish()  bind(C,name='belosfinish')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
     end subroutine belosfinish


     subroutine build_maps(vectorSize, nets,nete,np,nlev,datavector) &
                bind(C,name='BuildMaps')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            integer(c_int)                :: nets,nete,np,nlev
            type(c_ptr)                   :: datavector
     end subroutine build_maps


     subroutine build_matrix(vectorSize, nets,nete,np,nlev,datavector) &
                bind(C,name='BuildMatrix')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            integer(c_int)                :: nets,nete,np,nlev
            type(c_ptr)                   :: datavector
     end subroutine build_matrix



     subroutine build_rhs(nets,nete,np,nlev,rhs,datavector) &
                bind(C,name='BuildRHS')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: nets,nete,np,nlev
            real(c_double)  ,dimension(*) :: rhs
            type(c_ptr)                   :: datavector
     end subroutine build_rhs



     subroutine set_problem() &
                bind(C,name='SetProblem')
     end subroutine set_problem


     subroutine helmholtz_solve(vectorSize, lhs) &
                bind(C,name='HelmholtzSolve')
       use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
            integer(c_int)                :: vectorSize
            real(c_double)  ,dimension(*) :: lhs
     end subroutine helmholtz_solve

    end interface



    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (derivative_t)               :: deriv          ! non staggered derivative struct     (private)
    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)


    ! ===========
    ! Local
    ! ===========
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    real (kind=real_kind) :: LHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: RHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: sol(np,np,nlev,nets:nete)   ! exact solution
    real (kind=real_kind) :: solver_wts(npsq,nete-nets+1)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr
    real (kind=real_kind) :: snlat,cslat,cslon,snlon,xc,yc,zc, res, res_sol
    integer rhsVectorSize
    real(kind=real_kind):: rhsVector(np*np*nlev*(nete-nets+1))
    real(kind=real_kind):: solVector(np*np*nlev*(nete-nets+1))

    type(derived_type) ,target         :: helmholtzdata
    type(derived_type) ,pointer        :: f_ptr_helmholtzdata=>NULL()
    type(c_ptr)                        :: c_ptr_helmholtzdata

    integer linVecndx


    integer lenx
    real (kind=real_kind) pmean,dt
    type (TimeLevel_t) tl
    lenx=0
    pmean=0.0d0
    dt=0.0d0
    



!!  c_ptr_helmholtzdata should contain everything needed to evaluate the linear helmholtz operator

    call initialize(helmholtzdata, lenx, elem, pmean,edge1,edge1, edge1, &
            hybrid, deriv, dt, tl, nets, nete)


    f_ptr_helmholtzdata => helmholtzdata
    c_ptr_helmholtzdata =  c_loc(f_ptr_helmholtzdata)
!




    rhsVectorSize=np*np*nlev*(nete-nets+1)

    solVector=0.0d0



    if (hybrid%masterthread) print *,'creating manufactured solution'
    do ie=nets,nete
       iptr=1
       do j=1,np
          do i=1,np
             solver_wts(iptr,ie-nets+1) = elem(ie)%spheremp(i,j)
             iptr=iptr+1
          end do
       end do
    end do
    call cg_create(cg, npsq, nlev, nete-nets+1, hybrid, 0, solver_wts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make up an exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                snlat = SIN(elem(ie)%spherep(i,j)%lat)
                cslat = COS(elem(ie)%spherep(i,j)%lat)
                snlon = SIN(elem(ie)%spherep(i,j)%lon)
                cslon = COS(elem(ie)%spherep(i,j)%lon)
  
                xc = cslat*cslon
                yc = cslat*snlon
                zc = snlat

                ! take a couple of low-freq spherical harmonics for the solution
                sol(i,j,k,ie) = 1*xc + 2*yc + 3*zc + 4*xc*yc + 5*xc*zc + 6*yc*zc + &
                     7*(xc*yc*zc) 

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the RHS from our exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       do k=1,nlev
          RHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)
          call edgeVpack(edge1, RHS(1,1,1,ie), nlev, 0, ie)
       end do
    end do
    call bndry_exchangeV(cg%hybrid,edge1)

    linVecndx=1

    do ie=nets,nete
       ! unpack RHS
       call edgeVunpack(edge1, RHS(1,1,1,ie), nlev, 0, ie)
       do k=1,nlev
          RHS(:,:,k,ie)=RHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
       enddo
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       !
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
                rhsVector(linVecndx)=rhs(i,j,k,ie)
                ! rhsVector(linVecndx)=1.0d0
                ! rhsVector(linVecndx)=sol(i,j,k,ie) 
                linVecndx=linVecndx+1
             enddo
          enddo
       enddo
    enddo




  ! construt maps to build trilinos Matrix and rhs vectors
  call build_maps(rhsVectorSize, nets,nete,np,nlev,c_ptr_helmholtzdata) 
  ! assemble trilinos RHS vector from HOMME rhs vector - unique dof's across processors
  call build_rhs(nets,nete,np,nlev,rhsVector,c_ptr_helmholtzdata) 
  ! assemble Trilinos Matrix - unique dof's across processors
  call build_matrix(rhsVectorSize, nets,nete,np,nlev,c_ptr_helmholtzdata) 
  ! configure linear solver with ML preconditioner for solving Helmholtz
  ! with Belos iterative solver
  call set_problem() 
  ! Solve Helmholtz equation with Trilinos return solution in solVector
  call helmholtz_solve(rhsVectorSize, solVector) 


! an example of building another rhs and calling solver again
!
!    linVecndx=1
!    do linVecndx=1,rhsVectorSize
!           rhsVector(linVecndx)=1.0*rhsVector(linVecndx)
!    enddo
!
!
!  for each rhs vector send to "build_rhs" and then call "helmholtz_solve" 
!   
!  call build_rhs(nets,nete,np,nlev,rhsVector,c_ptr_helmholtzdata) 
!       call t_startf('timer_trilinosmlhelm')
!  call helmholtz_solve(rhsVectorSize, solVector) 
!       call t_stopf('timer_trilinosmlhelm')

! call belosfinish to delete memory for the trilinos Matrix maps and Helmholtz
! solver
  !end second rhs test

  call belosfinish()

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop
    ! in this version, we keep the residual C0 by DSSing the LHS
    ! and initializiont with a C0 RHS.  
    ! the update x, based on the last residual, will already be C0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cg%debug_level=1  ! 2=output every iterations
    maxits = 250
    tol=1d-10
    if (hybrid%masterthread) print *,'running solver V1 (C0 RHS) tol=',tol
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here:
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! solve x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;
             !LHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*x(:,:) + &
             !     alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)

             LHS(:,:,k,ie)=elem(ie)%rspheremp(:,:)*(elem(ie)%spheremp(:,:)*x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.))
             !LHS(:,:,k,ie)=elem(ie)%rspheremp(:,:)*(elem(ie)%spheremp(:,:)*x(:,:) )
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)

          ieptr=ie-nets+1
          do k=1,nlev
             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    tol=tol/2
!!PAL
tol=1.e-12
    if (hybrid%masterthread) print *,'running solver V2 (DG RHS) tol=',tol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! VERSION 2.  SAVE 1 DSS
    ! (important since we need to get iterations down to about 5 to be competitive)
    ! compute the RHS from our exact solution 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We should not have to apply DSS to RHS, since our equation:
    ! < PHI, lap(x) > = < PHI, DSS(RHS)>
    ! But since DSS(PHI)=PHI, and DSS is self adjoint,
    ! < PHI, DSS(RHS)>= < PHI, RHS>
    !
    ! in this version, we do not need to DSS the RHS (or LHS). 
    ! But the update x needs to be C0, so we DSS only x:
    !
    ! ONE ISSUE:  tolerence is based on <RHS,RHS>, which is not 
    ! computed correctly (and will always be larger than <DSS(RHS),DSS(RHS)>
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ! note: weak from laplace operator includes mass, so remove it
       do k=1,nlev
          RHS(:,:,k,ie)=sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)&
               / elem(ie)%spheremp(:,:)
       end do
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here: (note: r is not C0)
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! DSS x to make it C0.  use LHS for storage:
             LHS(:,:,k,ie)=x(:,:)*elem(ie)%spheremp(:,:)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             x(:,:)=LHS(:,:,k,ie) ! x() is now C0

             ! compute LHS(x) = x + laplace(x)
             LHS(:,:,k,ie)=x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)&
                  /elem(ie)%spheremp(:,:)

             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)    ! new LHS, DG
                   cg%state(ieptr)%z(iptr,k) = x(i,j)           ! z must be C0
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop
    !print *,'solver test CG iter = ',cg%iter


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================

    if (hybrid%masterthread) print *,"Trilinos Solution" , "             ", "HOMME Solution","              ","Exact Solution"
    linVecndx=1
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                if (hybrid%masterthread) print *,solVector(linVecndx) , " ", LHS(i,j,k,ie)," ",sol(i,j,k,ie)
                !set LHS to the solVector to test the l2 norm
                LHS(i,j,k,ie)=solVector(linVecndx)
                linVecndx=linVecndx+1
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo

    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    deallocate(helmholtzdata%base)

  end subroutine solver_test_ml

#endif 


  ! ================================================
  ! solver_test:
  !
  !    L(x) = laplace_sphere_wk(x) = -< grad(PHI) dot grad(x) >
  !    <   > = spheremp weighted inner-product
  !    L is self-adjoint:  < L(x),y> = < x,L(y) >
  !
  ! solve for x:
  !     <PHI,x> + a*L(x) =  < PHI, rhs >        
  !     contant a ~ 10 dx^2 (typical scaling in semi-implicit solve)
  !
  ! 2D solve - but applied to every level  k=1,nlev
  !
  ! In matrix notation, following the convention in M.T. and A.L.'s 
  ! "implicit.pdf" notes:
  !     D QQ^t (M + a*L ) x = rhs
  ! with:
  !   M    = multiply by spheremp 
  !   QQ^t = pack, exchange, unpack
  !   D    = multiply by rspheremp  D = Q V^-1 Q^-L   
  !          where V = the SEM diagonal mass matrix acting on vectors with no
  !          duplicate degrees of freedom.  V does not appear in HOMME.
  !   L is self adjoint w.r.t. M:    L(x) M y = x M L(y)
  !
  ! Note: if we solve laplace equation instead of Helmholz, we need to
  ! ensure < rhs,1>=0
  !
  ! ================================================
  subroutine solver_test(elem,edge1,red,hybrid,deriv,nets,nete)
    use dimensions_mod, only : nlev, np,npsq
    use element_mod, only : element_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use cg_mod, only : cg_t, congrad, cg_create
    use edge_mod, only : edgevpack, edgevunpack!, edgerotate
    use edgetype_mod, only : edgebuffer_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use control_mod, only : maxits, tol, precon_method
    use physical_constants, only : rrearth, dd_pi, rearth, omega
    use bndry_mod, only : bndry_exchangeV
    use linear_algebra_mod, only : matvec
    use parallel_mod, only : haltmp
    use hybrid_mod, only : hybrid_t
    use global_norms_mod, only : linf_snorm, l2_snorm

    integer, intent(in)  :: nets,nete
    type(element_t), intent(in), target :: elem(:)
    type (ReductionBuffer_ordered_1d_t)  :: red         ! CG reduction buffer   (shared memory)
    type (derivative_t), intent(in) :: deriv          ! non staggered derivative struct     (private)
    type (hybrid_t)             :: hybrid
    type (EdgeBuffer_t)               :: edge1          ! Laplacian divergence edge buffer (shared memory)


    ! ===========
    ! Local
    ! ===========
    type (cg_t)                       :: cg             ! conjugate gradient    (private)
    real (kind=real_kind) :: LHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: RHS(np,np,nlev,nets:nete)
    real (kind=real_kind) :: sol(np,np,nlev,nets:nete)   ! exact solution
    real (kind=real_kind) :: solver_wts(npsq,nete-nets+1)
    real (kind=real_kind) :: x(np,np)
    real (kind=real_kind) :: alambda = 10*250e3**2      ! low res test, dx=250km grid spacing

    integer :: ie
    integer :: i,j,k
    integer :: kptr
    integer :: iptr
    integer :: ieptr
    real (kind=real_kind) :: snlat,cslat,cslon,snlon,xc,yc,zc, res, res_sol

    if (hybrid%masterthread) print *,'creating manufactured solution'
    do ie=nets,nete
       iptr=1
       do j=1,np
          do i=1,np
             solver_wts(iptr,ie-nets+1) = elem(ie)%spheremp(i,j)
             iptr=iptr+1
          end do
       end do
    end do
    call cg_create(cg, npsq, nlev, nete-nets+1, hybrid, 0, solver_wts)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! make up an exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                snlat = SIN(elem(ie)%spherep(i,j)%lat)
                cslat = COS(elem(ie)%spherep(i,j)%lat)
                snlon = SIN(elem(ie)%spherep(i,j)%lon)
                cslon = COS(elem(ie)%spherep(i,j)%lon)
  
                xc = cslat*cslon
                yc = cslat*snlon
                zc = snlat

                ! take a couple of low-freq spherical harmonics for the solution
                sol(i,j,k,ie) = 1*xc + 2*yc + 3*zc + 4*xc*yc + 5*xc*zc + 6*yc*zc + &
                     7*(xc*yc*zc) 

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute the RHS from our exact solution
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       do k=1,nlev
          RHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)
          call edgeVpack(edge1, RHS(1,1,1,ie), nlev, 0, ie)
       end do
    end do
    call bndry_exchangeV(cg%hybrid,edge1)
    do ie=nets,nete
       ! unpack RHS
       call edgeVunpack(edge1, RHS(1,1,1,ie), nlev, 0, ie)
       do k=1,nlev
          RHS(:,:,k,ie)=RHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
       enddo
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       !
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop
    ! in this version, we keep the residual C0 by DSSing the LHS
    ! and initializiont with a C0 RHS.  
    ! the update x, based on the last residual, will already be C0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cg%debug_level=1  ! 2=output every iterations
    maxits = 250
    tol=1d-7
    if (hybrid%masterthread) print *,'running solver V1 (C0 RHS) tol=',tol
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here:
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! solve x + laplace(x)
             ! weak laplace operator already includes mass
             ! so only multiply x by mass;
             LHS(:,:,k,ie)=elem(ie)%spheremp(:,:)*x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res

    tol=tol/2
    if (hybrid%masterthread) print *,'running solver V2 (DG RHS) tol=',tol
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! VERSION 2.  SAVE 1 DSS
    ! (important since we need to get iterations down to about 5 to be competitive)
    ! compute the RHS from our exact solution 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! We should not have to apply DSS to RHS, since our equation:
    ! < PHI, lap(x) > = < PHI, DSS(RHS)>
    ! But since DSS(PHI)=PHI, and DSS is self adjoint,
    ! < PHI, DSS(RHS)>= < PHI, RHS>
    !
    ! in this version, we do not need to DSS the RHS (or LHS). 
    ! But the update x needs to be C0, so we DSS only x:
    !
    ! ONE ISSUE:  tolerence is based on <RHS,RHS>, which is not 
    ! computed correctly (and will always be larger than <DSS(RHS),DSS(RHS)>
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ie=nets,nete
       ! note: weak from laplace operator includes mass, so remove it
       do k=1,nlev
          RHS(:,:,k,ie)=sol(:,:,k,ie) + &
               alambda*laplace_sphere_wk(sol(:,:,k,ie),deriv,elem(ie),var_coef=.false.)&
               / elem(ie)%spheremp(:,:)
       end do
       !
       !  Initialize CG solver:  set %r = residual from initial guess
       !  if initial guess = 0, then we take %r=RHS
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                cg%state(ieptr)%r(iptr,k) = rhs(i,j,k,ie)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solver Loop 2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (congrad(cg,red,maxits,tol))
       do ie=nets,nete
          ieptr=ie-nets+1
          do k=1,nlev
             ! apply preconditioner here: (note: r is not C0)
             cg%state(ieptr)%z(:,k) = cg%state(ieptr)%r(:,k)
             
             !reshape(cg%state(ieptr)%z(:,k),(/np,np/))
             iptr=1
             do j=1,np
                do i=1,np
                   x(i,j) = cg%state(ieptr)%z(iptr,k)
                   iptr=iptr+1
                enddo
             enddo
             
             ! DSS x to make it C0.  use LHS for storage:
             LHS(:,:,k,ie)=x(:,:)*elem(ie)%spheremp(:,:)
             
          end do
          call edgeVpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
       end do
       call bndry_exchangeV(cg%hybrid,edge1)
       do ie=nets,nete
          ! unpack LHS
          call edgeVunpack(edge1, LHS(1,1,1,ie), nlev, 0, ie)
          do k=1,nlev
             LHS(:,:,k,ie)=LHS(:,:,k,ie)*elem(ie)%rspheremp(:,:)
          enddo

          ieptr=ie-nets+1
          do k=1,nlev
             x(:,:)=LHS(:,:,k,ie) ! x() is now C0

             ! compute LHS(x) = x + laplace(x)
             LHS(:,:,k,ie)=x(:,:) + &
                  alambda*laplace_sphere_wk(x,deriv,elem(ie),var_coef=.false.)&
                  /elem(ie)%spheremp(:,:)

             iptr=1
             do j=1,np
                do i=1,np
                   cg%state(ieptr)%s(iptr,k) = LHS(i,j,k,ie)    ! new LHS, DG
                   cg%state(ieptr)%z(iptr,k) = x(i,j)           ! z must be C0
                   iptr=iptr+1
                end do
             end do
          enddo
       enddo
       
    end do  ! CG solver while loop
    !print *,'solver test CG iter = ',cg%iter


    ! ===============================
    ! Converged! compute actual error (not residual computed in solver)
    ! ===============================
    do ie=nets,nete
       ieptr=ie-nets+1
       do k=1,nlev
          iptr=1
          do j=1,np
             do i=1,np
                LHS(i,j,k,ie) = cg%state(ieptr)%x(iptr,k)
                iptr=iptr+1
             enddo
          enddo
       enddo
    enddo
    res = l2_snorm(elem,LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized l2 error= ',res
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res


  end subroutine solver_test



end module solver_mod
