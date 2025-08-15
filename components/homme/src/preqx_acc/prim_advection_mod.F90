!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OpenACC implementations of some prim_advection_mod routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_advection_mod
  !OVERWRITING: Prim_Advec_Tracers_remap, prim_advec_init1, prim_advec_init2, Prim_Advec_Tracers_remap_rk2
  use kinds, only              : real_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, nelemd
  use physical_constants, only : rgas, Rwater_vapor, kappa, g, rearth, rrearth, cp
  use element_mod, only        : element_t
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, smooth, TimeLevel_Qdp
  use control_mod, only        : integration, test_case, hypervis_order, &
        statefreq, nu_q, nu_p, limiter_option, hypervis_subcycle_q, rsplit
  use edge_mod, only           : edgevpack, edgevunpack, initedgebuffer, initedgesbuffer, &
        edgevunpackmin, edge_g
 
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only   : abortmp
  use derivative_mod, only: derivative_t
  implicit none
  private
  real(kind=real_kind), private, allocatable :: qmin(:,:,:), qmax(:,:,:)
  real(kind=real_kind), private, allocatable :: Qtens_biharmonic(:,:,:,:,:)
  real(kind=real_kind), private, allocatable :: Qtens(:,:,:,:,:)
  real(kind=real_kind), private, allocatable :: grads_tracer(:,:,:,:,:,:)
  real(kind=real_kind), private, allocatable :: dp_star(:,:,:,:)
  type (EdgeBuffer_t), private :: edgeMinMax
  integer,parameter, private :: DSSeta = 1
  integer,parameter, private :: DSSomega = 2
  integer,parameter, private :: DSSdiv_vdp_ave = 3
  integer,parameter, private :: DSSno_var = -1
  real(kind=real_kind), allocatable, private :: data_pack(:,:,:,:), data_pack2(:,:,:,:)
  logical, private :: first_time = .true.

  public :: Prim_Advec_Tracers_remap
  public :: Prim_Advec_Tracers_observe_velocity
  public :: prim_advec_init1
  public :: prim_advec_init2

contains

  subroutine copy_qdp1_h2d( elem , tl , nets , nete )
    use element_mod, only: element_t
    use element_state, only: state_qdp
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl, nets , nete
    integer :: ie, k, j, i
    call t_startf('qdp1_pcie')
#   if USE_OPENACC
      do ie = nets , nete
        data_pack(:,:,:,ie) = state_qdp(:,:,:,1,tl,ie)
      enddo
      !$omp barrier
      !$omp master
!     do ie = 1 , nelemd
!       !$acc update device(state_qdp(:,:,:,1,tl,ie))
!     enddo
      !$acc update device(data_pack) async(1)
      !$acc parallel loop gang vector collapse(4) async(1) present(state_qdp,data_pack)
      do ie = 1 , nelemd
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              state_qdp(i,j,k,1,tl,ie) = data_pack(i,j,k,ie)
            enddo
          enddo
        enddo
      enddo
      !$acc wait(1)
      !$omp end master
      !$omp barrier
#   endif
    call t_stopf('qdp1_pcie')
  end subroutine copy_qdp1_h2d

  subroutine copy_qdp1_d2h( elem , tl , nets , nete )
    use element_mod, only: element_t
    use element_state, only: state_qdp
    use perf_mod, only: t_startf, t_stopf
    implicit none
    type(element_t), intent(in) :: elem(:)
    integer        , intent(in) :: tl, nets , nete
    integer :: ie, k, j, i
    call t_startf('qdp1_pcie')
#   if USE_OPENACC
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(4) async(1) present(state_qdp,data_pack)
      do ie = 1 , nelemd
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              data_pack(i,j,k,ie) = state_qdp(i,j,k,1,tl,ie)
            enddo
          enddo
        enddo
      enddo
      !$acc update host(data_pack) async(1)
      !$acc wait(1)
!     do ie = 1 , nelemd
!       !$acc update host(state_qdp(:,:,:,1,tl,ie))
!     enddo
      !$omp end master
      !$omp barrier
      do ie = nets , nete
        state_qdp(:,:,:,1,tl,ie) = data_pack(:,:,:,ie)
      enddo
#   endif
    call t_stopf('qdp1_pcie')
  end subroutine copy_qdp1_d2h

  subroutine Prim_Advec_Tracers_remap( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf, t_barrierf            ! _EXTERNAL
    use element_mod   , only: element_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod , only: hvcoord_t
    use hybrid_mod    , only: hybrid_t
    use time_mod      , only: TimeLevel_t, TimeLevel_Qdp
    use control_mod   , only: limiter_option, nu_p, qsplit
    use bndry_mod, only: bndry_exchangeV_timing
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete
    integer :: i
    if (first_time) then
      do i = 1 , 100
        call bndry_exchangeV_timing(hybrid,edge_g)
      enddo
      first_time = .false.
    endif
    call Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
  end subroutine Prim_Advec_Tracers_remap

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! forward-in-time 2 level vertically lagrangian step
  !  this code takes a lagrangian step in the horizontal 
  ! (complete with DSS), and then applies a vertical remap
  !
  ! This routine may use dynamics fields at timelevel np1
  ! In addition, other fields are required, which have to be 
  ! explicitly saved by the dynamics:  (in elem(ie)%derived struct)
  !
  ! Fields required from dynamics: (in 
  !    omega_p   it will be DSS'd here, for later use by CAM physics
  !              we DSS omega here because it can be done for "free"
  !    Consistent mass/tracer-mass advection (used if subcycling turned on)
  !       dp()   dp at timelevel n0
  !       vn0()  mean flux  < U dp > going from n0 to np1
  !
  ! 3 stage
  !    Euler step from t     -> t+.5
  !    Euler step from t+.5  -> t+1.0
  !    Euler step from t+1.0 -> t+1.5
  !    u(t) = u(t)/3 + u(t+2)*2/3
  !
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  subroutine Prim_Advec_Tracers_remap_rk2( elem , deriv , hvcoord , hybrid , dt , tl , nets , nete )
    use perf_mod      , only : t_startf, t_stopf, t_barrierf            ! _EXTERNAL
    use element_mod   , only: element_t
    use derivative_mod, only: derivative_t
    use hybvcoord_mod , only: hvcoord_t
    use hybrid_mod    , only: hybrid_t
    use time_mod      , only: TimeLevel_t, TimeLevel_Qdp
    use control_mod   , only: limiter_option, nu_p, qsplit
    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(inout) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    integer :: i,j,k,l,ie,q,nmin
    integer :: nfilt,rkstage,rhs_multiplier
    integer :: n0_qdp, np1_qdp

    call t_barrierf('sync_prim_advec_tracers_remap_k2', hybrid%par%comm)
    call t_startf('prim_advec_tracers_remap_rk2')
    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp) !time levels for qdp are not the same
    rkstage = 3 !   3 stage RKSSP scheme, with optimal SSP CFL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! RK2 2D advection step
    ! note: stage 3 we take the oppertunity to DSS omega
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! use these for consistent advection (preserve Q=1)
    ! derived%vdp_ave        =  mean horiz. flux:   U*dp
    ! derived%eta_dot_dpdn    =  mean vertical velocity (used for remap)
    ! derived%omega_p         =  advection code will DSS this for the physics, but otherwise 
    !                            it is not needed 
    ! Also: save a copy of div(U dp) in derived%div(:,:,:,1), which will be DSS'd 
    !       and a DSS'ed version stored in derived%div(:,:,:,2)
    call precompute_divdp( elem , hybrid , deriv , dt , nets , nete , n0_qdp )

    !rhs_multiplier is for obtaining dp_tracers at each stage:
    !dp_tracers(stage) = dp - rhs_multiplier*dt*divdp_proj

    call t_startf('euler_step_0')
    rhs_multiplier = 0
    call euler_step( np1_qdp , n0_qdp  , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSdiv_vdp_ave , rhs_multiplier )
    call t_stopf('euler_step_0')

    call t_startf('euler_step_1')
    rhs_multiplier = 1
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSeta         , rhs_multiplier )
    call t_stopf('euler_step_1')

    call t_startf('euler_step_2')
    rhs_multiplier = 2
    call euler_step( np1_qdp , np1_qdp , dt/2 , elem , hvcoord , hybrid , deriv , nets , nete , DSSomega       , rhs_multiplier )
    call t_stopf('euler_step_2')

    !to finish the 2D advection step, we need to average the t and t+2 results to get a second order estimate for t+1.  
    call qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Dissipation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ( limiter_option == 8  ) then
      ! dissipation was applied in RHS.  
    else
      call advance_hypervis_scalar(elem,hvcoord,hybrid,deriv,tl%np1,np1_qdp,nets,nete,dt)
    endif

    call t_stopf('prim_advec_tracers_remap_rk2')
  end subroutine prim_advec_tracers_remap_rk2

  subroutine prim_advec_init1(par, elem)
    use edge_mod    , only: initEdgeBuffer,initEdgeSBuffer
    use parallel_mod, only: parallel_t
    use element_mod , only: element_t
    implicit none
    type(parallel_t), intent(in) :: par
    type (element_t), intent(in) :: elem(:)
    call initEdgeSBuffer(par,edgeMinMax,elem(:),qsize*nlev*2          )

    allocate(dp_star(np,np,nlev,nelemd))
    allocate(qmin(nlev,qsize,nelemd))
    allocate(qmax(nlev,qsize,nelemd))
    allocate(Qtens_biharmonic(np,np,nlev,qsize,nelemd))
    allocate(qtens(np,np,nlev,qsize,nelemd))
    allocate(grads_tracer(np,np,2,nlev,qsize,nelemd))
    allocate(data_pack(np,np,nlev,nelemd))
    allocate(data_pack2(np,np,nlev,nelemd))
  end subroutine prim_advec_init1

  subroutine prim_advec_init2(elem,hvcoord,hybrid)
    use element_mod   , only: element_t
    use element_state , only: state_Qdp, derived_vn0, derived_divdp, derived_divdp_proj
    use hybvcoord_mod , only: hvcoord_t
    use hybrid_mod    , only: hybrid_t
    implicit none
    type(element_t)   , intent(in) :: elem(:)
    type(hvcoord_t)   , intent(in) :: hvcoord
    type(hybrid_t)    , intent(in) :: hybrid
    integer :: k, ie

    !$omp barrier
    !$omp master

    !$acc enter data pcreate(qmin,qmax,qtens_biharmonic,grads_tracer,dp_star,qtens,data_pack,data_pack2)
    !$acc enter data pcopyin(hvcoord%dp0)
    !$acc enter data pcopyin(edgeMinMax         ,edge_g)
    !$acc enter data pcopyin(edgeMinMax%desc    ,edge_g%desc)
    do ie = 1 , nelemd
    !$acc enter data pcopyin(edgeMinMax%desc(ie)    ,edge_g%desc(ie))
    enddo
    !$acc enter data pcopyin(edgeMinMax%moveLength    ,edge_g%moveLength)
    !$acc enter data pcopyin(edgeMinMax%movePtr0      ,edge_g%movePtr0)
    !$acc enter data pcopyin(edgeMinMax%srequest,edge_g%srequest)
    !$acc enter data pcopyin(edgeMinMax%rrequest,edge_g%rrequest)
    !$acc enter data pcopyin(edgeMinMax%receive ,edge_g%receive)
    !$acc enter data pcopyin(edgeMinMax%buf     ,edge_g%buf)

    !$acc enter data pcopyin(edgeMinMax%tag     ,edge_g%tag)

! this data has been removed
    !$XXX enter data pcopyin(edgeMinMax%putmap  ,edge_g%putmap)
    !$XXX enter data pcopyin(edgeMinMax%getmap  ,edge_g%getmap)
    !$XXX enter data pcopyin(edgeMinMax%reverse ,edge_g%reverse)

    !$omp end master
    !$omp barrier
  end subroutine prim_advec_init2

  subroutine Prim_Advec_Tracers_observe_velocity(elem, tl, n, nets, nete)
    type (element_t)     , intent(inout) :: elem(:)
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: n
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    ! Do nothing. Only SL transport uses this routine, and it's not supported in
    ! preqx.
  end subroutine Prim_Advec_Tracers_observe_velocity

  subroutine advance_hypervis_scalar(  elem , hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2 )
    !  hyperviscsoity operator for foward-in-time scheme
    !  take one timestep of:  
    !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    use kinds                , only: real_kind
    use dimensions_mod       , only: np, nlev
    use hybrid_mod           , only: hybrid_t
    use element_mod          , only: element_t
    use element_state        , only: derived_divdp_proj, state_qdp
    use derivative_mod       , only: derivative_t
    use perf_mod             , only: t_startf, t_stopf                          ! _EXTERNAL
    use hybvcoord_mod        , only: hvcoord_t
    use control_mod          , only: nu_q, hypervis_order, hypervis_subcycle_q, nu_p
    use viscosity_mod, only: biharmonic_wk_scalar_openacc
    use edge_mod     , only: edgeVpack_openacc, edgeVunpack_openacc
    use bndry_mod    , only: bndry_exchangeV => bndry_exchangeV_simple_overlap
    implicit none
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)     , intent(in   )         :: hvcoord
    type (hybrid_t)      , intent(in   )         :: hybrid
    type (derivative_t)  , intent(in   )         :: deriv
    integer              , intent(in   )         :: nt
    integer              , intent(in   )         :: nt_qdp
    integer              , intent(in   )         :: nets
    integer              , intent(in   )         :: nete
    real (kind=real_kind), intent(in   )         :: dt2
    ! local
    real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: min_neigh
    real (kind=real_kind), dimension(      nlev,qsize,nets:nete) :: max_neigh
    integer :: k , kptr , i , j , ie , ic , q
    real (kind=real_kind), dimension(np,np) :: lap_p
    real (kind=real_kind) :: v1,v2,dt
    real (kind=real_kind) :: dp
    integer :: density_scaling = 0
    if ( nu_q           == 0 ) return
    if ( hypervis_order /= 2 ) return
    call t_startf('advance_hypervis_scalar')
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  hyper viscosity  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dt = dt2 / hypervis_subcycle_q

    do ic = 1 , hypervis_subcycle_q
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(4) present(elem(:),derived_divdp_proj,state_qdp,hvcoord%dp0,qtens_biharmonic,qtens)
      do ie = 1 , nelemd
        ! Qtens = Q/dp   (apply hyperviscsoity to dp0 * Q, not Qdp)
        do k = 1 , nlev
          ! various options:
          !   1)  biharmonic( Qdp )
          !   2)  dp0 * biharmonic( Qdp/dp )    
          !   3)  dpave * biharmonic(Q/dp)
          ! For trace mass / mass consistenciy, we use #2 when nu_p=0
          ! and #e when nu_p>0, where dpave is the mean mass flux from the nu_p
          ! contribution from dynamics.
          do j = 1 , np
            do i = 1 , np
              dp = elem(ie)%derived%dp(i,j,k) - dt2 * derived_divdp_proj(i,j,k,ie)
              if (nu_p > 0) then
                do q = 1 , qsize
                  Qtens(i,j,k,q,ie) = elem(ie)%derived%dpdiss_ave(i,j,k)*state_Qdp(i,j,k,q,nt_qdp,ie) / dp 
                enddo
              else
                do q = 1 , qsize
                  Qtens(i,j,k,q,ie) = hvcoord%dp0(k)*state_Qdp(i,j,k,q,nt_qdp,ie) / dp
                enddo
              endif
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
      ! compute biharmonic operator. Qtens = input and output 
      call biharmonic_wk_scalar_openacc( elem , Qtens , grads_tracer , deriv , edge_g , hybrid , 1 , nelemd )
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(5) present(state_qdp,elem(:),qtens)
      do ie = 1 , nelemd
        do q = 1 , qsize
          do k = 1 , nlev
            do j = 1 , np
              do i = 1 , np
                ! advection Qdp.  For mass advection consistency:
                ! DIFF( Qdp) ~   dp0 DIFF (Q)  =  dp0 DIFF ( Qdp/dp )  
                state_Qdp(i,j,k,q,nt_qdp,ie) = state_Qdp(i,j,k,q,nt_qdp,ie) * elem(ie)%spheremp(i,j) - dt * nu_q * Qtens(i,j,k,q,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      call limiter2d_zero(state_Qdp,2,nt_qdp)
      call t_startf('ah_scalar_PEU')
      call edgeVpack_openacc(edge_g,state_qdp,qsize*nlev,0,qsize*nlev,1,nelemd,2,nt_qdp)
      !$omp end master
      !$omp barrier

      call t_startf('ah_scalar_exch')
      call bndry_exchangeV( hybrid , edge_g )
      call t_stopf('ah_scalar_exch')
      
      !$omp barrier
      !$omp master
      call edgeVunpack_openacc(edge_g,state_qdp,qsize*nlev,0,qsize*nlev,1,nelemd,2,nt_qdp)
      call t_stopf('ah_scalar_PEU')
      !$acc parallel loop gang vector collapse(5) present(state_qdp,elem(:))
      do ie = 1 , nelemd
        do q = 1 , qsize    
          ! apply inverse mass matrix
          do k = 1 , nlev
            do j = 1 , np
              do i = 1 , np
                state_Qdp(i,j,k,q,nt_qdp,ie) = elem(ie)%rspheremp(i,j) * state_Qdp(i,j,k,q,nt_qdp,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
    enddo
    call copy_qdp1_d2h( elem , nt_qdp , nets , nete )
    call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar

  subroutine qdp_time_avg( elem , rkstage , n0_qdp , np1_qdp , limiter_option , nu_p , nets , nete )
    use element_mod, only: element_t
    use element_state, only: state_qdp
    implicit none
    type(element_t)     , intent(inout) :: elem(:)
    integer             , intent(in   ) :: rkstage , n0_qdp , np1_qdp , nets , nete , limiter_option
    real(kind=real_kind), intent(in   ) :: nu_p
    integer :: ie,q,k,j,i
    !$omp barrier
    !$omp master
    !$acc parallel loop gang vector collapse(5) present(state_qdp)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          do j = 1 , np
            do i = 1 , np
              state_Qdp(i,j,k,q,np1_qdp,ie) = ( state_Qdp(i,j,k,q,n0_qdp ,ie) + (rkstage-1)*state_Qdp(i,j,k,q,np1_qdp,ie) ) / rkstage
            enddo
          enddo
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
    if (limiter_option == 8) call copy_qdp1_d2h( elem , np1_qdp , nets , nete )
    
  end subroutine qdp_time_avg

  subroutine euler_step( np1_qdp , n0_qdp , dt , elem , hvcoord , hybrid , deriv , nets , nete , DSSopt , rhs_multiplier )
  ! ===================================
  ! This routine is the basic foward
  ! euler component used to construct RK SSP methods
  !
  !           u(np1) = u(n0) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! n0 can be the same as np1.  
  !
  ! DSSopt = DSSeta or DSSomega:   also DSS eta_dot_dpdn or omega
  !
  ! ===================================
  use kinds                 , only: real_kind
  use dimensions_mod        , only: np, nlev
  use hybrid_mod            , only: hybrid_t
  use element_mod           , only: element_t
  use derivative_mod        , only: derivative_t
  use hybvcoord_mod         , only: hvcoord_t
  use control_mod           , only: limiter_option, nu_p, nu_q
  use perf_mod              , only: t_startf, t_stopf
  use element_state         , only: derived_divdp_proj, state_qdp, derived_vn0, derived_divdp
  use derivative_mod, only: divergence_sphere_openacc
  use viscosity_mod , only: biharmonic_wk_scalar_openacc, neighbor_minmax_openacc
  use edge_mod      , only: edgeVpack_openacc, edgeVunpack_openacc
  use bndry_mod     , only: bndry_exchangeV => bndry_exchangeV_simple_overlap
  implicit none
  integer              , intent(in   )         :: np1_qdp, n0_qdp
  real (kind=real_kind), intent(in   )         :: dt
  type (element_t)     , intent(inout), target :: elem(:)
  type (hvcoord_t)     , intent(in   )         :: hvcoord
  type (hybrid_t)      , intent(in   )         :: hybrid
  type (derivative_t)  , intent(in   )         :: deriv
  integer              , intent(in   )         :: nets
  integer              , intent(in   )         :: nete
  integer              , intent(in   )         :: DSSopt
  integer              , intent(in   )         :: rhs_multiplier
  ! local
  real(kind=real_kind) :: vstar(2)
  real(kind=real_kind) :: tmp, mintmp, maxtmp, dp
  integer :: ie,q,i,j,k
  integer :: rhs_viss


  do ie = nets , nete
    data_pack (:,:,:,ie) = elem(ie)%derived%dpdiss_ave
    data_pack2(:,:,:,ie) = elem(ie)%derived%dp       
  enddo
  !$omp barrier
  !$omp master
  !$acc update device(data_pack,data_pack2) async(1)
  !$acc update device(derived_divdp_proj,derived_vn0) async(2)
  !$acc parallel loop gang vector collapse(4) present(data_pack,elem) async(1)
  do ie = 1 , nelemd
    do k = 1 , nlev
      do j = 1 , np
        do i = 1 , np
          elem(ie)%derived%dpdiss_ave(i,j,k) = data_pack (i,j,k,ie)
          elem(ie)%derived%dp        (i,j,k) = data_pack2(i,j,k,ie)
        enddo
      enddo
    enddo
  enddo
  !$acc wait
  !$omp end master
  !$omp barrier

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   compute Q min/max values for lim8
  !   compute biharmonic mixing term f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  rhs_viss = 0
  if ( limiter_option == 8  ) then
    ! when running lim8, we also need to limit the biharmonic, so that term needs
    ! to be included in each euler step.  three possible algorithms here:
    ! 1) most expensive:
    !     compute biharmonic (which also computes qmin/qmax) during all 3 stages
    !     be sure to set rhs_viss=1
    !     cost:  3 biharmonic steps with 3 DSS
    !
    ! 2) cheapest:
    !     compute biharmonic (which also computes qmin/qmax) only on first stage
    !     be sure to set rhs_viss=3
    !     reuse qmin/qmax for all following stages (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps with 1 DSS
    !     main concern:  viscosity 
    !     
    ! 3)  compromise:
    !     compute biharmonic (which also computes qmin/qmax) only on last stage
    !     be sure to set rhs_viss=3
    !     compute qmin/qmax directly on first stage
    !     reuse qmin/qmax for 2nd stage stage (but update based on local qmin/qmax)
    !     cost:  1 biharmonic steps, 2 DSS
    !
    !  NOTE  when nu_p=0 (no dissipation applied in dynamics to dp equation), we should
    !        apply dissipation to Q (not Qdp) to preserve Q=1
    !        i.e.  laplace(Qdp) ~  dp0 laplace(Q)                
    !        for nu_p=nu_q>0, we need to apply dissipation to Q * diffusion_dp
    !
    ! initialize dp, and compute Q from Qdp (and store Q in Qtens_biharmonic)
    !$omp barrier
    !$omp master
    !$acc parallel loop gang vector collapse(4) private(tmp) present(elem(:),derived_divdp_proj,state_qdp,qtens_biharmonic)
    do ie = 1 , nelemd
      ! add hyperviscosity to RHS.  apply to Q at timelevel n0, Qdp(n0)/dp
      do k = 1 , nlev    !  Loop index added with implicit inversion (AAM)
        do j = 1 , np
          do i = 1 , np
            tmp = elem(ie)%derived%dp(i,j,k) - rhs_multiplier*dt*derived_divdp_proj(i,j,k,ie) 
            do q = 1 , qsize
              Qtens_biharmonic(i,j,k,q,ie) = state_Qdp(i,j,k,q,n0_qdp,ie) / tmp
            enddo
          enddo
        enddo
      enddo
    enddo
    !$acc parallel loop gang vector collapse(3) private(mintmp,maxtmp) present(qmin,qmax,qtens_biharmonic)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev    
          if (rhs_multiplier == 0 .or. rhs_multiplier == 2) then
            mintmp =  1.D20
            maxtmp = -1.D20
          else
            mintmp = qmin(k,q,ie)
            maxtmp = qmax(k,q,ie)
          endif
          do j = 1 , np
            do i = 1 , np
              mintmp = min(mintmp,qtens_biharmonic(i,j,k,q,ie))
              maxtmp = max(maxtmp,qtens_biharmonic(i,j,k,q,ie))
            enddo
          enddo
          qmin(k,q,ie) = max(mintmp,0d0)
          qmax(k,q,ie) =     maxtmp
        enddo
      enddo
    enddo
    !$omp end master
    !$omp barrier
    if ( rhs_multiplier == 0 ) call neighbor_minmax_openacc(elem,hybrid,edgeMinMax,1,nelemd,qmin,qmax)
    ! compute biharmonic mixing term
    if ( rhs_multiplier == 2 ) then
      rhs_viss = 3
      ! two scalings depending on nu_p:
      ! nu_p=0:    qtens_biharmonic *= dp0                   (apply viscsoity only to q)
      ! nu_p>0):   qtens_biharmonc *= elem()%psdiss_ave      (for consistency, if nu_p=nu_q)
      if ( nu_p > 0 ) then
        !$omp barrier
        !$omp master
        !$acc parallel loop gang vector collapse(4) present(elem(:),qtens_biharmonic,hvcoord%dp0)
        do ie = 1 , nelemd
          do k = 1 , nlev    
            do j = 1 , np
              do i = 1 , np
                tmp = elem(ie)%derived%dpdiss_ave(i,j,k) / hvcoord%dp0(k)
                do q = 1 , qsize
                  ! NOTE: divide by dp0 since we multiply by dp0 below
                  Qtens_biharmonic(i,j,k,q,ie) = Qtens_biharmonic(i,j,k,q,ie) * tmp
                enddo
              enddo
            enddo
          enddo
        enddo
        !$omp end master
        !$omp barrier
      endif
      call biharmonic_wk_scalar_openacc( elem , qtens_biharmonic , grads_tracer , deriv , edge_g , hybrid , 1 , nelemd )
      call neighbor_minmax_openacc( elem , hybrid , edgeMinMax , 1 , nelemd , qmin , qmax )
      !$omp barrier
      !$omp master
      !$acc parallel loop gang vector collapse(4) present(qtens_biharmonic,hvcoord%dp0,elem(:))
      do ie = 1 , nelemd
        do k = 1 , nlev    !  Loop inversion (AAM)
          do j = 1 , np
            do i = 1 , np
              tmp = -rhs_viss*dt*nu_q*hvcoord%dp0(k) / elem(ie)%spheremp(i,j)
              do q = 1 , qsize
                ! note: biharmonic_wk() output has mass matrix already applied. Un-apply since we apply again below:
                qtens_biharmonic(i,j,k,q,ie) = tmp*Qtens_biharmonic(i,j,k,q,ie)
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end master
      !$omp barrier
    endif
  endif  ! compute biharmonic mixing term and qmin/qmax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   2D Advection step
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute velocity used to advance Qdp 

  if ( limiter_option == 8 .and. nu_p > 0 .and. rhs_viss /= 0 ) then
    do ie = nets , nete
      data_pack(:,:,:,ie) = elem(ie)%derived%dpdiss_biharmonic
    enddo
    !$omp barrier
    !$omp master
    !$acc update device(data_pack) async(1)
    !$acc parallel loop gang vector collapse(4) present(elem,data_pack) async(1)
    do ie = 1 , nelemd
      do k = 1 , nlev
        do j = 1 , np
          do i = 1 , np
            elem(ie)%derived%dpdiss_biharmonic(i,j,k) = data_pack(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
    !$acc wait(1)
    !$omp end master
    !$omp barrier
  endif

  !$omp barrier
  !$omp master
  if (limiter_option == 8) then
    !$acc update device(derived_divdp)
  endif
  !$acc parallel loop gang vector collapse(4) present(elem(:),derived_divdp_proj,derived_vn0,dp_star,derived_divdp,grads_tracer,state_qdp) &
  !$acc& private(dp,vstar)
  do ie = 1 , nelemd
    do k = 1 , nlev    !  Loop index added (AAM)
      do j = 1 , np
        do i = 1 , np
          ! derived variable divdp_proj() (DSS'd version of divdp) will only be correct on 2nd and 3rd stage
          ! but that's ok because rhs_multiplier=0 on the first stage:
          dp = elem(ie)%derived%dp(i,j,k) - rhs_multiplier * dt * derived_divdp_proj(i,j,k,ie)
          Vstar(1) = derived_vn0(i,j,1,k,ie) / dp
          Vstar(2) = derived_vn0(i,j,2,k,ie) / dp
          if ( limiter_option == 8 ) then
            ! UN-DSS'ed dp at timelevel n0+1:  
            dp_star(i,j,k,ie) = dp - dt * derived_divdp(i,j,k,ie)  
            if ( nu_p > 0 .and. rhs_viss /= 0 ) then
              ! add contribution from UN-DSS'ed PS dissipation
              dp_star(i,j,k,ie) = dp_star(i,j,k,ie) - rhs_viss * dt * nu_q * elem(ie)%derived%dpdiss_biharmonic(i,j,k) / elem(ie)%spheremp(i,j)
            endif
          endif
          do q = 1 , qsize
            grads_tracer(i,j,1,k,q,ie) = Vstar(1) * state_Qdp(i,j,k,q,n0_qdp,ie)
            grads_tracer(i,j,2,k,q,ie) = Vstar(2) * state_Qdp(i,j,k,q,n0_qdp,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
  call divergence_sphere_openacc( grads_tracer , deriv , elem(:) , qtens , nlev*qsize , 1 , nelemd , 1 , 1 )
  !$acc parallel loop gang vector collapse(5) present(qtens,state_qdp,qtens_biharmonic)
  do ie = 1 , nelemd
    ! advance Qdp
    do q = 1 , qsize
      do k = 1 , nlev  !  dp_star used as temporary instead of divdp (AAM)
        do j = 1 , np
          do i = 1 , np
            tmp = state_Qdp(i,j,k,q,n0_qdp,ie) - dt * qtens(i,j,k,q,ie)
            ! optionally add in hyperviscosity computed above:
            if ( rhs_viss /= 0 ) tmp = tmp + Qtens_biharmonic(i,j,k,q,ie)
            Qtens(i,j,k,q,ie) = tmp
          enddo
        enddo
      enddo
    enddo
  enddo
  if ( limiter_option == 8 ) then
    call limiter_optim_iter_full( Qtens , elem(:) , qmin , qmax , dp_star )   ! apply limiter to Q = Qtens / dp_star 
  endif
  !$acc parallel loop gang vector collapse(4) present(state_Qdp,elem(:),qtens)
  do ie = 1 , nelemd
    ! apply mass matrix, overwrite np1 with solution:
    ! dont do this earlier, since we allow np1_qdp == n0_qdp 
    ! and we dont want to overwrite n0_qdp until we are done using it
    do k = 1 , nlev
      do j = 1 , np
        do i = 1 , np
          tmp = elem(ie)%spheremp(i,j)
          do q = 1 , qsize
            state_Qdp(i,j,k,q,np1_qdp,ie) =  tmp * Qtens(i,j,k,q,ie) 
          enddo
        enddo
      enddo
    enddo
  enddo
  if ( limiter_option == 4 ) then
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! sign-preserving limiter, applied after mass matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    call limiter2d_zero( state_Qdp , 2 , np1_qdp )
  endif
  ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
  ! all zero so we only have to DSS 1:nlev
  call t_startf('eus_PEU')
  call edgeVpack_openacc(edge_g , state_Qdp , nlev*qsize , 0 , nlev*qsize , 1 , nelemd , 2 , np1_qdp )
  !$omp end master
  !$omp barrier

  call t_startf('eus_exch')
  call bndry_exchangeV( hybrid , edge_g )
  call t_stopf('eus_exch')

  !$omp barrier
  !$omp master
  call edgeVunpack_openacc( edge_g , state_Qdp , nlev*qsize , 0 , nlev*qsize , 1 , nelemd , 2 , np1_qdp )
  call t_stopf('eus_PEU')
  !$acc parallel loop gang vector collapse(4) present(state_Qdp,elem(:))
  do ie = 1 , nelemd
    do k = 1 , nlev    !  Potential loop inversion (AAM)
      do j = 1 , np
        do i = 1 , np
          tmp = elem(ie)%rspheremp(i,j)
          do q = 1 , qsize
            state_Qdp(i,j,k,q,np1_qdp,ie) = tmp * state_Qdp(i,j,k,q,np1_qdp,ie)
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end master
  !$omp barrier
  end subroutine euler_step

  subroutine limiter2d_zero(Qdp,tdim,tl)
    ! mass conserving zero limiter (2D only).  to be called just before DSS
    !
    ! this routine is called inside a DSS loop, and so Q had already
    ! been multiplied by the mass matrix.  Thus dont include the mass
    ! matrix when computing the mass = integral of Q over the element
    !
    ! ps is only used when advecting Q instead of Qdp
    ! so ps should be at one timelevel behind Q
    implicit none
    integer              , intent(in   ) :: tdim
    integer              , intent(in   ) :: tl
    real (kind=real_kind), intent(inout) :: Qdp(np,np,nlev,qsize,tdim,nelemd)
    ! local
    real (kind=real_kind) :: mass,mass_new
    real (kind=real_kind) :: qtmp(np,np)
    integer i,j,k,q,ie
    !$acc parallel loop gang vector collapse(3) present(qdp) private(qtmp)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = nlev , 1 , -1
          !$acc cache(qtmp)
          qtmp = Qdp(:,:,k,q,tl,ie)
          mass = 0
          do j = 1 , np
            do i = 1 , np
              mass = mass + qtmp(i,j)
            enddo
          enddo
          ! negative mass.  so reduce all postive values to zero 
          ! then increase negative values as much as possible
          if ( mass < 0 ) qtmp = -qtmp 
          mass_new = 0
          do j = 1 , np
            do i = 1 , np
              if ( qtmp(i,j) < 0 ) then
                qtmp(i,j) = 0
              else
                mass_new = mass_new + qtmp(i,j)
              endif
            enddo
          enddo
          ! now scale the all positive values to restore mass
          if ( mass_new > 0 ) qtmp = qtmp * abs(mass) / mass_new
          if ( mass     < 0 ) qtmp = -qtmp 
          Qdp(:,:,k,q,tl,ie) = qtmp
        enddo
      enddo
    enddo
  end subroutine limiter2d_zero

  subroutine limiter_optim_iter_full(ptens,elem,minp,maxp,dpmass)
    use element_mod, only: element_t
    use kinds         , only : real_kind
    use dimensions_mod, only : np, np, nlev
    implicit none
    !THIS IS A NEW VERSION OF LIM8, POTENTIALLY FASTER BECAUSE INCORPORATES KNOWLEDGE FROM
    !PREVIOUS ITERATIONS
    !The idea here is the following: We need to find a grid field which is closest
    !to the initial field (in terms of weighted sum), but satisfies the min/max constraints.
    !So, first we find values which do not satisfy constraints and bring these values
    !to a closest constraint. This way we introduce some mass change (addmass),
    !so, we redistribute addmass in the way that l2 error is smallest. 
    !This redistribution might violate constraints thus, we do a few iterations. 
    real (kind=real_kind), intent(inout) :: ptens (np*np,nlev,qsize,nelemd)
    type(element_t)      , intent(in   ) :: elem  (:)
    real (kind=real_kind), intent(inout) :: minp  (      nlev,qsize,nelemd)
    real (kind=real_kind), intent(inout) :: maxp  (      nlev,qsize,nelemd)
    real (kind=real_kind), intent(in   ) :: dpmass(np*np,nlev      ,nelemd)
    integer :: k1, k, i, j, iter, i1, i2, q, ie
    real (kind=real_kind) :: addmass, weightssum, mass, mintmp, maxtmp, sumc
    real (kind=real_kind) :: x(np*np),c(np*np)
    integer :: maxiter = np*np-1
    real (kind=real_kind) :: tol_limiter = 5e-14

    !$acc parallel loop gang vector collapse(3) present(ptens,elem(:),minp,maxp,dpmass) private(c,x,mintmp,maxtmp,addmass,weightssum,mass,sumc) vector_length(512)
    do ie = 1 , nelemd
      do q = 1 , qsize
        do k = 1 , nlev
          !$acc cache(c,x)
          do k1 = 1 , np*np
            i = modulo(k1-1,np)+1
            j = (k1-1)/np+1
            c(k1) = elem(ie)%spheremp(i,j) * dpmass(k1,k,ie)
            x(k1) = ptens(k1,k,q,ie) / dpmass(k1,k,ie)
          enddo
          sumc = sum(c)
          mintmp = minp(k,q,ie)
          maxtmp = maxp(k,q,ie)
          if (sumc <= 0 ) CYCLE   ! this should never happen, but if it does, dont limit
          mass=sum(c*x)
          ! relax constraints to ensure limiter has a solution:
          ! This is only needed if runnign with the SSP CFL>1 or
          ! due to roundoff errors
          if( mass < mintmp*sumc ) mintmp = mass / sumc
          if( mass > maxtmp*sumc ) maxtmp = mass / sumc
          do iter=1,maxiter
            addmass=0.0d0
            do k1=1,np*np
              if((x(k1)>maxtmp)) then
                addmass=addmass+(x(k1)-maxtmp)*c(k1)
                x(k1)=maxtmp
              endif
              if((x(k1)<mintmp)) then
                addmass=addmass-(mintmp-x(k1))*c(k1)
                x(k1)=mintmp
              endif
            enddo !k1
            if(abs(addmass)<=tol_limiter*abs(mass)) exit
            weightssum=0.0d0
            if(addmass>0)then
              do k1=1,np*np
                if(x(k1)<maxtmp)then
                  weightssum=weightssum+c(k1)
                endif
              enddo !k1
              do k1=1,np*np
                if(x(k1)<maxtmp)then
                  x(k1)=x(k1)+addmass/weightssum
                endif
              enddo
            else
              do k1=1,np*np
                if(x(k1)>mintmp)then
                  weightssum=weightssum+c(k1)
                endif
              enddo
              do k1=1,np*np
                if(x(k1)>mintmp)then
                  x(k1)=x(k1)+addmass/weightssum
                endif
              enddo
            endif
          enddo!end of iteration
          minp(k,q,ie) = mintmp
          maxp(k,q,ie) = maxtmp
          ptens(:,k,q,ie) = x * dpmass(:,k,ie)
        enddo
      enddo
    enddo
  end subroutine limiter_optim_iter_full

  subroutine precompute_divdp( elem , hybrid , deriv , dt , nets , nete , n0_qdp )
    use element_mod           , only: element_t
    use element_state         , only: derived_vn0, derived_divdp, derived_divdp_proj
    use hybrid_mod            , only: hybrid_t
    use derivative_mod        , only: derivative_t
    use edge_mod              , only: edgeVpack_nlyr, edgeVunpack_nlyr
    use bndry_mod             , only: bndry_exchangeV
    use control_mod           , only: limiter_option
    use derivative_mod, only: divergence_sphere_openacc
    use openacc_utils_mod     , only: copy_ondev
    use perf_mod              , only: t_startf, t_stopf
    implicit none
    type(element_t)      , intent(inout) :: elem(:)
    type (hybrid_t)      , intent(in   ) :: hybrid
    type (derivative_t)  , intent(in   ) :: deriv
    real(kind=real_kind) , intent(in   ) :: dt
    integer              , intent(in   ) :: nets , nete , n0_qdp
    integer :: ie , k
    real(kind=real_kind), pointer, dimension(:,:,:) :: DSSvar
    call copy_qdp1_h2d(elem,n0_qdp,nets,nete)
    !$omp barrier
    !$omp master
    !$acc update device(derived_vn0)
    call divergence_sphere_openacc(derived_vn0,deriv,elem,derived_divdp,nlev,1,nelemd,1,1)
    call copy_ondev(derived_divdp_proj,derived_divdp,product(shape(derived_divdp)))
    !$acc update host(derived_divdp,derived_divdp_proj)
    !$omp end master
    !$omp barrier
    call t_startf('derived PEU')
    do ie = nets , nete
      ! note: eta_dot_dpdn is actually dimension nlev+1, but nlev+1 data is
      ! all zero so we only have to DSS 1:nlev
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k) 
        elem(ie)%derived%omega_p(:,:,k)      = elem(ie)%spheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)      
        elem(ie)%derived%divdp_proj(:,:,k)   = elem(ie)%spheremp(:,:) * elem(ie)%derived%divdp_proj(:,:,k)   
      enddo
      call edgeVpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev) , nlev , 0      , 3*nlev )
      call edgeVpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%omega_p(:,:,1:nlev)      , nlev , nlev   , 3*nlev )
      call edgeVpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%divdp_proj(:,:,1:nlev)   , nlev , nlev*2 , 3*nlev )
    enddo

    call bndry_exchangeV( hybrid , edge_g   )

    do ie = nets , nete
      call edgeVunpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%eta_dot_dpdn(:,:,1:nlev) , nlev , 0      , 3*nlev )
      call edgeVunpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%omega_p(:,:,1:nlev)      , nlev , nlev   , 3*nlev )
      call edgeVunpack_nlyr( edge_g , elem(ie)%desc, elem(ie)%derived%divdp_proj(:,:,1:nlev)   , nlev , nlev*2 , 3*nlev )
      do k = 1 , nlev
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%derived%eta_dot_dpdn(:,:,k) 
        elem(ie)%derived%omega_p(:,:,k)      = elem(ie)%rspheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)      
        elem(ie)%derived%divdp_proj(:,:,k)   = elem(ie)%rspheremp(:,:) * elem(ie)%derived%divdp_proj(:,:,k)   
      enddo
    enddo
    call t_stopf('derived PEU')
  end subroutine precompute_divdp

end module prim_advection_mod


