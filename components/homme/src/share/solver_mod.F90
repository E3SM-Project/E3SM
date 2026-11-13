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



contains




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
    use edge_mod, only : edgevpack, edgevunpack
    use edgetype_mod, only : edgebuffer_t
    use derivative_mod, only : derivative_t, laplace_sphere_wk
    use control_mod, only : maxits, tol
    use bndry_mod, only : bndry_exchangeV
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
    integer :: iptr
    integer :: ieptr
    real (kind=real_kind) :: snlat,cslat,cslon,snlon,xc,yc,zc, res

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
    if (res > 100*tol) call abortmp('ERROR: CG solver test failed')
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
    if (res > 100*tol) call abortmp('ERROR: DG solver test failed')
    res = linf_snorm(LHS,sol,hybrid,np,nets,nete) 
    if (hybrid%masterthread) print *,'normalized linf error= ',res


  end subroutine solver_test



end module solver_mod
