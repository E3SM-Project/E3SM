#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

module sl_advection
  use kinds, only              : real_kind, int_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, qsize_d
  use derivative_mod, only     : derivative_t, gradient_sphere, divergence_sphere
  use element_mod, only        : element_t
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, TimeLevel_Qdp
  use control_mod, only        : integration, test_case, hypervis_order, transport_alg, limiter_option,&
                                 vert_remap_q_alg
  use edge_mod, only           : edgevpack_nlyr, edgevunpack_nlyr, edge_g
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only       : abortmp, parallel_t
  use coordinate_systems_mod, only : cartesian3D_t
  use compose_mod

  implicit none

  private

  real(real_kind), parameter :: zero = 0.0_real_kind, fourth = 0.25_real_kind, &
       half = 0.5_real_kind, one = 1.0_real_kind, two = 2.0_real_kind, &
       eps = epsilon(1.0_real_kind)

  type (cartesian3D_t), allocatable :: dep_points_all(:,:,:,:) ! (np,np,nlev,nelemd)
  real(kind=real_kind), dimension(:,:,:,:,:), allocatable :: minq, maxq ! (np,np,nlev,qsize,nelemd)
  logical :: is_sphere

  ! For use in make_positive.
  real(kind=real_kind) :: dp_tol

  public :: prim_advec_tracers_remap_ALE, sl_init1, sl_vertically_remap_tracers, sl_unittest

  ! For testing
  public :: calc_trajectory, dep_points_all, sphere2cart

  ! For C++
  public :: sl_get_params

  logical, parameter :: barrier = .false.

contains

  !=================================================================================================!

  subroutine sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test, independent_time_steps)
    use control_mod, only: dt_remap_factor, dt_tracer_factor

    integer, intent(in) :: transport_alg
    logical, intent(out) :: slmm, cisl, qos, sl_test, independent_time_steps

    slmm = transport_alg > 1
    cisl = transport_alg == 2 .or. transport_alg == 3 .or. transport_alg >= 20
    qos  = cisl .and. (transport_alg == 3 .or. transport_alg == 39)  
    sl_test = (transport_alg >= 17 .and. transport_alg <= 19) .or. &
         transport_alg == 29 .or. transport_alg == 39
    ! Either dt_remap_factor = 0 (vertically Eulerian dynamics) or
    ! dt_remap_factor < dt_tracer_factor (vertically Lagrangian
    ! dynamics' vertical remap time step < tracer time step).
    independent_time_steps = dt_remap_factor < dt_tracer_factor
  end subroutine sl_parse_transport_alg

  subroutine sphere2cart(sphere, cart)
    use coordinate_systems_mod, only: spherical_polar_t, cartesian3D_t, change_coordinates

    type (spherical_polar_t), intent(in) :: sphere
    type (cartesian3D_t), intent(out) :: cart

    if (is_sphere) then
       cart = change_coordinates(sphere)
    else
       ! See conventions established in planar_mod::coordinates_atomic.
       cart%x = sphere%lon
       cart%y = sphere%lat
       cart%z = 0
    end if
  end subroutine sphere2cart

  subroutine sl_init1(par, elem)
    use interpolate_mod,        only : interpolate_tracers_init
    use control_mod,            only : transport_alg, semi_lagrange_cdr_alg, cubed_sphere_map, &
         nu_q, semi_lagrange_hv_q, semi_lagrange_cdr_check, geometry
    use element_state,          only : timelevels
    use coordinate_systems_mod, only : cartesian3D_t
    use perf_mod, only: t_startf, t_stopf
    use kinds, only: iulog

    type (parallel_t) :: par
    type (element_t) :: elem(:)
    type (cartesian3D_t) :: pinside
    integer :: nslots, ie, num_neighbors, need_conservation, i, j
    logical :: slmm, cisl, qos, sl_test, independent_time_steps

#ifdef HOMME_ENABLE_COMPOSE
    call t_startf('sl_init1')
    if (transport_alg > 0) then
       call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test, independent_time_steps)
       if (par%masterproc .and. nu_q > 0 .and. semi_lagrange_hv_q > 0) &
            write(iulog,*) 'COMPOSE> use HV; nu_q, all:', nu_q, semi_lagrange_hv_q
       is_sphere = trim(geometry) /= 'plane'
       nslots = nlev*qsize
       ! Technically a memory leak, but the array persists for the entire
       ! run, so not a big deal for now.
       allocate(dep_points_all(np,np,nlev,size(elem)))
       do ie = 1, size(elem)
          ! Provide a point inside the target element.
          call sphere2cart(elem(ie)%spherep(2,2), pinside)
          num_neighbors = elem(ie)%desc%actual_neigh_edges + 1
          call slmm_init_local_mesh(ie, elem(ie)%desc%neigh_corners, num_neighbors, &
               pinside, size(elem(ie)%desc%neigh_corners,2))
          if (sl_test) then
             do j = 1,np
                do i = 1,np
                   call sphere2cart(elem(ie)%spherep(i,j), pinside)
                   call slmm_check_ref2sphere(ie, pinside)
                end do
             end do
          end if
       end do
       call slmm_init_finalize()
       if (semi_lagrange_cdr_alg > 1) then
          need_conservation = 1
          call cedr_sl_init(np, nlev, qsize, qsize_d, timelevels, need_conservation)
       end if
       allocate(minq(np,np,nlev,qsize,size(elem)), maxq(np,np,nlev,qsize,size(elem)))
       dp_tol = -one
    endif
    call t_stopf('sl_init1')
#endif
  end subroutine sl_init1

  subroutine sl_get_params(nu_q_out, hv_scaling, hv_q, hv_subcycle_q, limiter_option_out, &
       cdr_check, geometry_type) bind(c)
    use control_mod, only: semi_lagrange_hv_q, hypervis_subcycle_q, semi_lagrange_cdr_check, &
         nu_q, hypervis_scaling, limiter_option, geometry
    use iso_c_binding, only: c_int, c_double

    real(c_double), intent(out) :: nu_q_out, hv_scaling
    integer(c_int), intent(out) :: hv_q, hv_subcycle_q, limiter_option_out, cdr_check, geometry_type

    nu_q_out = nu_q
    hv_scaling = hypervis_scaling
    hv_q = semi_lagrange_hv_q
    hv_subcycle_q = hypervis_subcycle_q
    limiter_option_out = limiter_option
    cdr_check = 0
    if (semi_lagrange_cdr_check) cdr_check = 1
    geometry_type = 0 ! sphere
    if (trim(geometry) == "plane") geometry_type = 1

  end subroutine sl_get_params

  subroutine prim_advec_tracers_remap_ALE(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete)
    use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
    use dimensions_mod,         only : max_neigh_edges
    use interpolate_mod,        only : interpolate_tracers, minmax_tracers
    use control_mod,            only : dt_tracer_factor, nu_q, transport_alg, semi_lagrange_hv_q, &
         semi_lagrange_cdr_alg, semi_lagrange_cdr_check
    ! For DCMIP16 supercell test case.
    use control_mod,            only : dcmip16_mu_q
    use prim_advection_base,    only : advance_physical_vis
    use compose_mod,            only : compose_h2d, compose_d2h
    use iso_c_binding,          only : c_bool

    implicit none
    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    type(cartesian3D_t)   :: dep_points  (np,np)

    integer :: i,j,k,l,n,q,ie,n0_qdp,np1_qdp
    integer :: scalar_q_bounds, info
    logical :: slmm, cisl, qos, sl_test, independent_time_steps
    logical(kind=c_bool) :: h2d, d2h

#ifdef HOMME_ENABLE_COMPOSE
    call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
    call t_startf('Prim_Advec_Tracers_remap_ALE')

    call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test, independent_time_steps)
    ! Until I get the DSS onto GPU, always need to h<->d.
    !h2d = hybrid%par%nprocs > 1 .or. semi_lagrange_cdr_check .or. & (semi_lagrange_hv_q > 0 .and. nu_q > 0)
    h2d = .true.
    d2h = compose_d2h .or. h2d
    h2d = compose_h2d .or. h2d

    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)

    call calc_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, &
         independent_time_steps, nets, nete)

    call t_startf('SLMM_csl')
    !todo Here and in the set-pointer loop for CEDR, do just in the first call.
    do ie = nets, nete
       call slmm_csl_set_elem_data(ie, elem(ie)%metdet, &
            elem(ie)%state%Qdp, n0_qdp, &
            elem(ie)%derived%dp, elem(ie)%state%Q, &
            elem(ie)%desc%actual_neigh_edges + 1, &
            h2d, d2h)
    end do
    ! edge_g buffers are shared by SLMM, CEDR, other places in HOMME, and
    ! dp_coupling in EAM. Thus, we must take care to protected threaded
    ! access. In the following, "No barrier needed" comments justify why a
    ! barrier isn't needed.
    ! No barrier needed: ale_rkdss has a horiz thread barrier at the end.
    call slmm_csl(nets, nete, dep_points_all, minq, maxq, info)
    ! No barrier needed: slmm_csl has a horiz thread barrier at the end.
    if (info /= 0) then
       call write_velocity_data(elem, nets, nete, hybrid, deriv, dt, tl)
       call abortmp('slmm_csl returned -1; see output above for more information.')
    end if
    if (barrier) call perf_barrier(hybrid)
    call t_stopf('SLMM_csl')

    if (semi_lagrange_hv_q > 0 .and. nu_q > 0) then
       n = semi_lagrange_hv_q
       call advance_hypervis_scalar(elem, hvcoord, hybrid, deriv, tl%np1, np1_qdp, nets, nete, dt, n)
       ! No barrier needed: advance_hypervis_scalar has a horiz thread barrier at the end.
    end if

    ! CEDR works with either classical SL or IR.
    if (semi_lagrange_cdr_alg > 1 .and. semi_lagrange_cdr_alg /= 42) then
       scalar_q_bounds = 0
       call cedr_sl_set_pointers_begin(nets, nete)
       do ie = nets, nete
          call cedr_sl_set_spheremp(ie, elem(ie)%spheremp)
          call cedr_sl_set_dp0(hvcoord%dp0)
          call cedr_sl_set_Qdp(ie, elem(ie)%state%Qdp, n0_qdp, np1_qdp)
          if (independent_time_steps) then
             call cedr_sl_set_dp(ie, elem(ie)%derived%divdp) ! dp_star
          else
             call cedr_sl_set_dp3d(ie, elem(ie)%state%dp3d, tl%np1)
          end if
          call cedr_sl_set_Q(ie, elem(ie)%state%Q)
       end do
       call cedr_sl_set_pointers_end(h2d, d2h)
       call t_startf('CEDR')
       ! No barrier needed: A barrier was already called.
       call cedr_sl_run_global(minq, maxq, nets, nete)
       ! No barrier needed: run_cdr has a horiz thread barrier at the end.
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR')
       call t_startf('CEDR_local')
       call cedr_sl_run_local(minq, maxq, nets, nete, scalar_q_bounds, limiter_option)
       ! Barrier needed to protect edge_g buffers use in CEDR.
#if (defined HORIZ_OPENMP)
       !$omp barrier
#endif
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR_local')
    else
       do ie = nets, nete
          do k = 1, nlev
             do q = 1, qsize
                elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Q(:,:,k,q) * &
                     elem(ie)%state%dp3d(:,:,k,tl%np1)
             enddo
          enddo
       end do
    end if
    ! Technically, dss_Qdp is needed only if semi_lagrange_cdr_alg > 1;
    ! otherwise, each Q(dp) field is already continuous. But SL transport would
    ! be run without a CDR only by a specialist doing some debugging. Meanwhile,
    ! dss_Qdp also DSSes derived%omega_p for diagnostics. It's important not to
    ! forget to do that, and we don't need to optimize for the
    ! semi_lagrange_cdr_alg <= 1 case, so just always call dss_Qdp.
    call t_startf('SL_dss')
    call dss_Qdp(elem, nets, nete, hybrid, np1_qdp)
    if (barrier) call perf_barrier(hybrid)
    call t_stopf('SL_dss')
    if (semi_lagrange_cdr_check) then
       call t_startf('CEDR_check')
       call cedr_sl_check(minq, maxq, nets, nete)
       if (barrier) call perf_barrier(hybrid)
       call t_stopf('CEDR_check')
    end if
    ! physical viscosity for supercell test case
    if (dcmip16_mu_q > 0) then
       call advance_physical_vis(elem, hvcoord, hybrid, deriv, tl%np1, np1_qdp, nets, nete, dt, dcmip16_mu_q)
    endif
    call t_stopf('Prim_Advec_Tracers_remap_ALE')
#endif
  end subroutine prim_advec_tracers_remap_ALE
  
  subroutine calc_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, &
       independent_time_steps, nets, nete)
    use vertremap_base, only : remap1

    type (element_t)     , intent(inout) :: elem(:)
    type (derivative_t)  , intent(in   ) :: deriv
    type (hvcoord_t)     , intent(in   ) :: hvcoord
    type (hybrid_t)      , intent(in   ) :: hybrid
    real(kind=real_kind) , intent(in   ) :: dt
    type (TimeLevel_t)   , intent(in   ) :: tl
    logical              , intent(in   ) :: independent_time_steps
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    real(kind=real_kind) :: wr(np,np,nlev,2)
    integer :: ie, k

    do ie = nets,nete
       elem(ie)%derived%vn0 = elem(ie)%state%v(:,:,:,:,tl%np1)
    end do
    if (independent_time_steps) then
       call t_startf('SLMM_reconstruct')
       if (dp_tol < zero) then
          ! Thread write race condition; benign b/c written value is same in all threads.
          call set_dp_tol(hvcoord, dp_tol)
       end if
       do ie = nets,nete
          ! divdp is dp_star
          call calc_vertically_lagrangian_levels(hybrid, elem(ie), ie, hvcoord, tl, dt, &
               deriv, dp_tol, elem(ie)%derived%divdp)
          wr(:,:,:,1) = elem(ie)%derived%vn0(:,:,1,:)*elem(ie)%state%dp3d(:,:,:,tl%np1)
          wr(:,:,:,2) = elem(ie)%derived%vn0(:,:,2,:)*elem(ie)%state%dp3d(:,:,:,tl%np1)
          call remap1(wr, np, 2, elem(ie)%state%dp3d(:,:,:,tl%np1), elem(ie)%derived%divdp,&
               vert_remap_q_alg)
          elem(ie)%derived%vn0(:,:,1,:) = wr(:,:,:,1)/elem(ie)%derived%divdp
          elem(ie)%derived%vn0(:,:,2,:) = wr(:,:,:,2)/elem(ie)%derived%divdp
       end do
       call t_stopf('SLMM_reconstruct')
    end if

    call ALE_RKdss(elem, nets, nete, hybrid, deriv, dt, tl, independent_time_steps)

    if (barrier) call perf_barrier(hybrid)
    call t_startf('SLMM_v2x')
    do ie = nets, nete
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k)
#endif
       do k = 1, nlev
          call ALE_departure_from_gll(dep_points_all(:,:,k,ie), &
               elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt, normalize=is_sphere)
       end do
    end do
    call t_stopf('SLMM_v2x')
  end subroutine calc_trajectory

  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE ALE_RKDSS
  ! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012
  !
  ! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid
  !
  ! CALLS:
  ! INPUT:
  !
  ! OUTPUT:
  !-----------------------------------------------------------------------------------!

  ! this will calculate the velocity at time t+1/2  along the trajectory s(t) given the velocities
  ! at the GLL points at time t and t+1 using a second order time accurate formulation.
  subroutine ALE_RKdss(elem, nets, nete, hy, deriv, dt, tl, independent_time_steps)
    use derivative_mod,  only : derivative_t, ugradv_sphere
    use edgetype_mod,    only : EdgeBuffer_t
    use bndry_mod,       only : bndry_exchangev
    use kinds,           only : real_kind
    use hybrid_mod,      only : hybrid_t
    use element_mod,     only : element_t
    use dimensions_mod,  only : np, nlev

    implicit none

    type (element_t)     , intent(inout)             :: elem(:)
    integer              , intent(in   )             :: nets
    integer              , intent(in   )             :: nete
    type (hybrid_t)      , intent(in)                :: hy ! distributed parallel structure (shared)
    type (derivative_t)  , intent(in)                :: deriv ! derivative struct
    real (kind=real_kind), intent(in)                :: dt ! timestep
    type (TimeLevel_t)   , intent(in)                :: tl
    logical              , intent(in)                :: independent_time_steps

    integer                                          :: ie, k
    real (kind=real_kind), dimension(np,np,2)        :: vtmp
    integer :: nlyr

    ! RK-SSP 2 stage 2nd order:
    !     x*(t+1) = x(t) + U(x(t),t) dt
    !     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) ) dt
    ! apply taylor series:
    !  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
    !
    ! x(t+1) = x(t) +  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) ) dt
    ! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) )
    ! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + (x*(t+1)-x(t)) gradU(x(t),t+1) )
    ! (x(t+1) - x(t)) / dt =  1/2 ( U(x(t),t+1) + U(x(t),t) + U(x(t),t) dt  gradU(x(t),t+1) )
    !
    !
    !  (x(t+1)-x(t))/dt =  1/2(U(x(t),t+1) + U(x(t),t) + dt U(x(t),t) gradU(x(t),t+1))
    !
    ! suppose dt = -ts (we go backward)
    !  (x(t-ts)-x(t))/-ts =  1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)
    !
    !  x(t-ts) = x(t)) -ts * [ 1/2( U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts) ]

    nlyr = 2*nlev
    if (independent_time_steps) nlyr = nlyr + nlev

    do ie = nets,nete
       ! vstarn0 = U(x,t)
       ! vstar   = U(x,t+1)
       do k=1,nlev
          ! vstar is v at the start of the tracer time step, and vn0
          ! is v at the end of the tracer time step.
          vtmp(:,:,:) = ugradv_sphere(elem(ie)%derived%vn0(:,:,:,k), elem(ie)%derived%vstar(:,:,:,k), &
               deriv, elem(ie))

          elem(ie)%derived%vstar(:,:,:,k) = &
               (elem(ie)%derived%vn0(:,:,:,k) + elem(ie)%derived%vstar(:,:,:,k))/2 - dt*vtmp(:,:,:)/2

          elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)* &
               elem(ie)%spheremp*elem(ie)%rspheremp
          elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)* &
               elem(ie)%spheremp*elem(ie)%rspheremp
          if (independent_time_steps) then
             ! divdp contains the reconstructed dp.
             elem(ie)%derived%divdp(:,:,k) = elem(ie)%derived%divdp(:,:,k)*elem(ie)%spheremp* &
                  elem(ie)%rspheremp
          end if
       enddo
       call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,nlyr)
       if (independent_time_steps) &
            call edgeVpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%divdp,nlev,2*nlev,nlyr)
    enddo

    call t_startf('ALE_RKdss_bexchV')
    call bndry_exchangeV(hy,edge_g)
    call t_stopf('ALE_RKdss_bexchV')

    do ie = nets,nete
       call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%vstar,2*nlev,0,nlyr)
       if (independent_time_steps) then
          call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%divdp,nlevp,2*nlev,nlyr)
       end if
    end do

#if (defined HORIZ_OPENMP)
    !$omp barrier
#endif
  end subroutine ALE_RKdss

  subroutine write_velocity_data(elem, nets, nete, hy, deriv, dt, tl)
    use derivative_mod,  only : derivative_t, ugradv_sphere
    use edgetype_mod,    only : EdgeBuffer_t
    use bndry_mod,       only : bndry_exchangev
    use kinds,           only : real_kind
    use hybrid_mod,      only : hybrid_t
    use element_mod,     only : element_t
    use dimensions_mod,  only : np, nlev
    implicit none

    type (element_t)     , intent(inout) :: elem(:)
    integer              , intent(in)    :: nets
    integer              , intent(in)    :: nete
    type (hybrid_t)      , intent(in)    :: hy
    type (derivative_t)  , intent(in)    :: deriv
    real (kind=real_kind), intent(in)    :: dt
    type (TimeLevel_t)   , intent(in)    :: tl

    integer :: ie, i, j, k, np1
    real (kind=real_kind) :: max_v, max_vstar

    np1 = tl%np1
    max_v = 0.d0
    max_vstar = 0.d0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                max_v = max(max_v, sqrt(sum(elem(ie)%state%v(i,j,:,k,np1)**2)))
                max_vstar = max(max_vstar, sqrt(sum(elem(ie)%derived%vstar(i,j,:,k)**2)))
             end do
          end do
       end do
    end do
    print *, 'max_v, max_vstar on rank', max_v, max_vstar, hy%par%rank
  end subroutine write_velocity_data

  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
  ! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
  !                                                                                   !
  ! CALLS:
  ! INPUT:
  !
  ! OUTPUT:
  !-----------------------------------------------------------------------------------!
  subroutine ALE_departure_from_gll(acart, vstar, elem, dt, normalize)
    use physical_constants,     only : scale_factor
    use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t
    use time_mod,               only : timelevel_t
    use element_mod,            only : element_t
    use kinds,                  only : real_kind
    use dimensions_mod,         only : np

    implicit none

    type (cartesian3D_t)    ,intent(out)  :: acart(np,np)
    real (kind=real_kind)   ,intent(in)   :: vstar(np,np,2)
    type (element_t)        ,intent(in)   :: elem
    real (kind=real_kind)   ,intent(in)   :: dt
    logical, intent(in) :: normalize

    integer                               :: i,j

    real (kind=real_kind)                 :: uxyz (np,np,3), norm

    ! convert velocity from lat/lon to cartesian 3D
    do i=1,3
       ! Summing along the third dimension is a sum over components for each point.
       ! (This is just a faster way of doing a dot product for each grid point,
       ! since reindexing the inputs to use the intrinsic effectively would be
       ! just asking for trouble.)
       uxyz(:,:,i)=sum( elem%vec_sphere2cart(:,:,i,:)*vstar(:,:,:) ,3)
    end do
    ! compute departure point
    ! crude, 1st order accurate approximation.  to be improved
    do i=1,np
       do j=1,np
          call sphere2cart(elem%spherep(i,j), acart(i,j))
          acart(i,j)%x = acart(i,j)%x - dt*uxyz(i,j,1)/scale_factor
          acart(i,j)%y = acart(i,j)%y - dt*uxyz(i,j,2)/scale_factor
          acart(i,j)%z = acart(i,j)%z - dt*uxyz(i,j,3)/scale_factor
          if (normalize) then
             norm = sqrt(acart(i,j)%x*acart(i,j)%x + acart(i,j)%y*acart(i,j)%y + &
                  acart(i,j)%z*acart(i,j)%z)
             acart(i,j)%x = acart(i,j)%x / norm
             acart(i,j)%y = acart(i,j)%y / norm
             acart(i,j)%z = acart(i,j)%z / norm
          end if
       enddo
    enddo

  end subroutine ALE_departure_from_gll

  subroutine dss_Qdp(elem, nets, nete, hybrid, np1_qdp)
    use edgetype_mod,    only : EdgeBuffer_t
    use bndry_mod,       only : bndry_exchangev
    use hybrid_mod,      only : hybrid_t
    use element_mod,     only : element_t
    implicit none

    type (element_t), intent(inout) :: elem(:)
    integer         , intent(in   ) :: nets, nete, np1_qdp
    type (hybrid_t) , intent(in)    :: hybrid
    integer                         :: ie, q, k

    do ie = nets, nete
       do q = 1, qsize
          do k = 1, nlev
             elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%spheremp(:,:)
          end do
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,np1_qdp), &
            qsize*nlev, 0, (qsize+1)*nlev)
       ! Also DSS omega_p for diagnostics. This has nothing to do with tracer
       ! transport; it's just a good place to do the DSS.
       do k = 1, nlev
          elem(ie)%derived%omega_p(:,:,k) = elem(ie)%spheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%omega_p(:,:,:), &
            nlev, qsize*nlev, (qsize+1)*nlev)
    enddo

    call t_startf('SLMM_bexchV')
    call bndry_exchangeV(hybrid, edge_g)
    call t_stopf('SLMM_bexchV')

    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%state%Qdp(:,:,:,:,np1_qdp), &
            qsize*nlev, 0, (qsize+1)*nlev)
       do q = 1, qsize
          do k = 1, nlev
             elem(ie)%state%Qdp(:,:,k,q,np1_qdp) = elem(ie)%state%Qdp(:,:,k,q,np1_qdp)*elem(ie)%rspheremp(:,:)
          end do
       end do
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%omega_p(:,:,:), &
            nlev, qsize*nlev, (qsize+1)*nlev)
       do k = 1, nlev
          elem(ie)%derived%omega_p(:,:,k) = elem(ie)%rspheremp(:,:) * elem(ie)%derived%omega_p(:,:,k)
       end do
    end do
  end subroutine dss_Qdp

  subroutine perf_barrier(hybrid)
    use hybrid_mod, only : hybrid_t
    implicit none
    type (hybrid_t), intent(in) :: hybrid
    integer :: ierr

#ifdef HORIZ_OPENMP
    !$OMP BARRIER
#endif
    if (hybrid%ithr == 0) call mpi_barrier(hybrid%par%comm, ierr)
#ifdef HORIZ_OPENMP
    !$OMP BARRIER
#endif
  end subroutine perf_barrier

  subroutine advance_hypervis_scalar(elem, hvcoord , hybrid , deriv , nt , nt_qdp , nets , nete , dt2, nq)
    !  hyperviscsoity operator for foward-in-time scheme
    !  take one timestep of:
    !          Q(:,:,:,np) = Q(:,:,:,np) +  dt2*nu*laplacian**order ( Q )
    !
    !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
    use kinds          , only : real_kind
    use dimensions_mod , only : np, nlev
    use hybrid_mod     , only : hybrid_t
    use element_mod    , only : element_t
    use derivative_mod , only : derivative_t
    use bndry_mod      , only : bndry_exchangev
    use perf_mod       , only : t_startf, t_stopf                          ! _EXTERNAL
    use control_mod    , only : nu_q, nu_p, hypervis_subcycle_q
    implicit none
    type (element_t)     , intent(inout), target :: elem(:)
    type (hvcoord_t)     , intent(in   )         :: hvcoord
    type (hybrid_t)      , intent(in   )         :: hybrid
    type (derivative_t)  , intent(in   )         :: deriv
    integer              , intent(in   )         :: nt
    integer              , intent(in   )         :: nt_qdp
    integer              , intent(in   )         :: nets
    integer              , intent(in   )         :: nete
    integer              , intent(in   )         :: nq
    real (kind=real_kind), intent(in   )         :: dt2

    ! local
    real (kind=real_kind), dimension(np,np,nlev,nq,nets:nete) :: Qtens
    real (kind=real_kind), dimension(np,np,nlev             ) :: dp
    real (kind=real_kind) :: dt
    integer :: k , ie , ic , q

    if ( nu_q           == 0 ) return
    if ( hypervis_order /= 2 ) return
    !   call t_barrierf('sync_advance_hypervis_scalar', hybrid%par%comm)
    call t_startf('advance_hypervis_scalar')

    dt = dt2 / hypervis_subcycle_q

    do ic = 1 , hypervis_subcycle_q
       do ie = nets , nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k) collapse(2)
#endif
          do q = 1 , nq
             do k = 1 , nlev
                Qtens(:,:,k,q,ie) = elem(ie)%state%Q(:,:,k,q)
             enddo
          enddo
       enddo ! ie loop

       ! compute biharmonic operator. Qtens = input and output
       call biharmonic_wk_scalar( elem , Qtens , deriv , edge_g , hybrid , nets , nete, nq )

       do ie = nets , nete
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k,j,i) collapse(2)
#endif
          do q = 1 , nq
             do k = 1 , nlev
                elem(ie)%state%Q(:,:,k,q) = elem(ie)%state%Q(:,:,k,q) * elem(ie)%spheremp &
                     - dt * nu_q * Qtens(:,:,k,q,ie)
             enddo
          enddo
          call edgeVpack_nlyr(edge_g , elem(ie)%desc, elem(ie)%state%Q , nq*nlev , 0 , nq*nlev )
       enddo ! ie loop

       call t_startf('ah_scalar_bexchV')
       call bndry_exchangeV( hybrid , edge_g )
       call t_stopf('ah_scalar_bexchV')

       do ie = nets , nete
          call edgeVunpack_nlyr(edge_g , elem(ie)%desc, elem(ie)%state%Q , nq*nlev , 0, nq*nlev)
#if (defined COLUMN_OPENMP)
          !$omp parallel do private(q,k) collapse(2)
#endif
          do q = 1 , nq
             ! apply inverse mass matrix
             do k = 1 , nlev
                elem(ie)%state%Q(:,:,k,q) = elem(ie)%rspheremp(:,:) * elem(ie)%state%Q(:,:,k,q)
             enddo
          enddo
       enddo ! ie loop
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
       !$OMP BARRIER
#endif
#endif
    enddo
    call t_stopf('advance_hypervis_scalar')
  end subroutine advance_hypervis_scalar

  subroutine biharmonic_wk_scalar(elem,qtens,deriv,edgeq,hybrid,nets,nete,nq)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute weak biharmonic operator
    !    input:  qtens = Q
    !    output: qtens = weak biharmonic of Q
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use control_mod, only : hypervis_scaling
    use derivative_mod, only : laplace_sphere_wk

    type (hybrid_t)      , intent(in) :: hybrid
    type (element_t)     , intent(inout), target :: elem(:)
    integer :: nets,nete,nq
    real (kind=real_kind), dimension(np,np,nlev,nq,nets:nete) :: qtens
    type (EdgeBuffer_t)  , intent(inout) :: edgeq
    type (derivative_t)  , intent(in) :: deriv

    ! local
    integer :: k,kptr,i,j,ie,ic,q
    real (kind=real_kind), dimension(np,np) :: lap_p
    logical var_coef1

    !if tensor hyperviscosity with tensor V is used, then biharmonic operator is (\grad\cdot V\grad) (\grad \cdot \grad) 
    !so tensor is only used on second call to laplace_sphere_wk
    var_coef1 = .true.
    if(hypervis_scaling > 0)    var_coef1 = .false.

    do ie=nets,nete
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k, q, lap_p)
#endif
       do q=1,nq
          do k=1,nlev    !  Potential loop inversion (AAM)
             lap_p(:,:)=qtens(:,:,k,q,ie)
             ! Original use of qtens on left and right hand sides caused OpenMP errors (AAM)
             qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=var_coef1)
          enddo
          call edgeVpack_nlyr(edgeq, elem(ie)%desc, qtens(:,:,:,q,ie),nlev,nlev*(q-1),nq*nlev)
       enddo
    enddo

    call t_startf('biwksc_bexchV')
    call bndry_exchangeV(hybrid,edgeq)
    call t_stopf('biwksc_bexchV')

    do ie=nets,nete

       ! apply inverse mass matrix, then apply laplace again
#if (defined COLUMN_OPENMP)
       !$omp parallel do private(k, q, lap_p)
#endif
       do q=1,nq      
          call edgeVunpack_nlyr(edgeq,elem(ie)%desc,qtens(:,:,:,q,ie),nlev,nlev*(q-1),nq*nlev)
          do k=1,nlev    !  Potential loop inversion (AAM)
             lap_p(:,:)=elem(ie)%rspheremp(:,:)*qtens(:,:,k,q,ie)
             qtens(:,:,k,q,ie)=laplace_sphere_wk(lap_p,deriv,elem(ie),var_coef=.true.)
          enddo
       enddo
    enddo
#ifdef DEBUGOMP
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER
#endif
#endif
  end subroutine biharmonic_wk_scalar

  subroutine calc_vertically_lagrangian_levels( &
       hybrid, elem, ie, hvcoord, tl, dt, deriv, dp_tol, dprecon)

    ! Reconstruct the vertically Lagrangian levels, thus permitting the dynamics
    ! vertical remap time step to be shorter than the tracer time step.
    !   Recall
    !     p(eta,ps) = A(eta) p0 + B(eta) ps
    !     => dp/dt = p_eta deta/dt + p_ps dps/dt
    !              = (A_eta p0 + B_eta ps) deta/dt + B(eta) dps/dt
    ! dp3d already accounts for B(eta) dps/dt, so it does not appear in what
    ! follows. We use eta as the vertical coordinate in these calculations. But
    ! here, to focus on the key ideas, consider a p = (x,z) system with
    ! velocities v = (u,w). Let 0, h, 1, suffixes denote start, middle, and end
    ! of the time step. The algorithm is as follows.
    !   The trajectory is computed forward in time, so the departure 0 points
    ! are on the grid. We need to compute z(x0,t1), the floating height at
    ! horizontal grid point x0 at time t1. We start by computing z1 = z(x1,t1),
    ! the floating height at the non-grid point x1 at time t1 or, in other
    ! words, the non-grid arrival point:
    !     z1 = z0 + dt/2 (w(p0,t0) + w(p1,t1)) + O(dt^3)
    !     w(p1,t1) = w(p0,t1) + grad w(p0,t1) (p1 - p0) + O(dt^2)
    !     x1 - x0 = dt u(p0,t0) + O(dt^2)
    !     z1 - z0 = dt w(p0,t0) + O(dt^2)
    !     z1 = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
    !                     dt (w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0))) + O(dt^3)  (*)
    ! Now we compute z(x0,t1). First, we need
    !     x0 - x1 = -dt u(p0,t0) + O(dt^2)
    ! and
    !     z_x(x1,t1) = z0_x + dt w_x(p0,t1) + O(dt^2)
    !                = dt w_x(p0,t1) + O(dt^2).
    !     z_xx(x1,t1) = z0_xx + dt w_xx(p0,t1) + O(dt^2)
    !                 = dt w_xx(p0,t1) + O(dt^2).
    ! Then we expand z(x0,t1) in Taylor series and substitute:
    !     z(x0,t1) = z(x1,t1) + z_x(x1,t1) (x0 - x1) + z_xx(x1,t1) (x0 - x1)^2
    !                + O(|x0 - x1|^3)
    !              = z(x1,t1) - dt z_x(x1,t1) u(p0,t0)
    !                + dt^2 z_xx(x1,t1) u(p0,t0)^2
    !                + O(dt^3) + O(|x0 - x1|^3)
    !              = z(x1,t1) - dt^2 w_x(p0,t1) u(p0,t0)
    !                + dt^3 w_xx(x0,t1) u(p0,t0)^2
    !                + O(dt^3) + O(|x0 - x1|^3).
    ! Gather all O(dt^3) terms, with O(|x0 - x1|^p) = O(dt^p):
    !     z(x0,t1) = z(x1,t1) - dt^2 w_x(p0,t1) u(p0,t0) + O(dt^3).
    ! Now substitute (*) for z(x1,t1):
    !     z(x0,t1) = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
    !                           dt (w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0)))
    !                - dt^2 w_x(p0,t1) u(p0,t0) + O(dt^3)
    !              = z0 + dt/2 (w(p0,t0) + w(p0,t1) +
    !                           dt (-w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0)))
    !                + O(dt^3)
    ! This is locally accurate to O(dt^3) and so globally 2nd-order
    ! accurate. Notably, compared with (*), this formula differs only in a
    ! sign. Note also that a straightforward first-order accurate formula is
    !     z(x0,t1) = z0 + dt w(p0,th) + O(dt^2).

    use control_mod, only: dt_remap_factor
    use derivative_mod, only: derivative_t, gradient_sphere
    use kinds, only: iulog

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(in) :: elem
    integer, intent(in) :: ie
    type (hvcoord_t), intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    real(kind=real_kind), intent(in) :: dt, dp_tol
    type (derivative_t), intent(in) :: deriv
    real(kind=real_kind), intent(out) :: dprecon(np,np,nlev)

    real(real_kind), dimension(np,np,nlevp) :: pref, p1r
    real(real_kind), dimension(np,np,nlevp,2) :: eta_dot_dpdn
    real(real_kind), dimension(np,np) :: dps, ptp0, v1, v2, divdp
    real(real_kind), dimension(np,np,2) :: grad, vdp
    real(real_kind) :: dp_neg_min
    integer :: i, j, k, k1, k2, d, t

#ifndef NDEBUG
    if (abs(hvcoord%hybi(1)) > 10*eps .or. hvcoord%hyai(nlevp) > 10*eps) then
       if (hybrid%masterthread) &
            print *, 'calc_vertically_lagrangian_levels: bi(1)', hvcoord%hybi(1), 'ai(nlevp)', &
            hvcoord%hyai(nlevp)
       call abortmp('hvcoord has unexpected non-0 entries at the bottom and/or top')
    end if
#endif

    ! Reconstruct an approximation to endpoint eta_dot_dpdn on
    ! Eulerian levels.
    do t = 1,2
       eta_dot_dpdn(:,:,1,t) = zero
       do k = 1,nlev
          do d = 1,2
             if (t == 1) then      
                vdp(:,:,d) = elem%derived%vstar(:,:,d,k)*elem%derived%dp(:,:,k)
             else
                vdp(:,:,d) = elem%derived%vn0(:,:,d,k)*elem%state%dp3d(:,:,k,tl%np1)
             end if
          end do
          divdp = divergence_sphere(vdp, deriv, elem)
          eta_dot_dpdn(:,:,k+1,t) = eta_dot_dpdn(:,:,k,t) + divdp
       end do
       dps = eta_dot_dpdn(:,:,nlevp,t)
       eta_dot_dpdn(:,:,nlevp,t) = zero
       do k = 2,nlev
          eta_dot_dpdn(:,:,k,t) = hvcoord%hybi(k)*dps - eta_dot_dpdn(:,:,k,t)
       end do
    end do

    ! Use p0 as the reference coordinate system. p0 differs from p1 by B(eta)
    ! (ps1 - ps0); dp3d already accounts for this term
    ! w.r.t. derived%dp. Recall
    !     eta_dot_dpdn = p_eta eta_dot = (A_eta p0 + B_eta ps) deta/dt,
    ! except that in the code eta_dot_dpdn is actually dp deta/dt rather than
    ! dp/deta deta/dt. eta_dot_dpdn is the motion of a pressure level excluding
    ! its motion due to dps/dt.
    call calc_p(hvcoord, elem%derived%dp, pref)

    do k = 2, nlev
       ! Gradient of eta_dot_dpdn = p_eta deta/dt at final time
       ! w.r.t. horizontal sphere coords.
       grad = gradient_sphere(eta_dot_dpdn(:,:,k,2), deriv, elem%Dinv)

       ! Gradient of eta_dot_dpdn = p_eta deta/dt at final time w.r.t. p at
       ! initial time.
       k1 = k-1
       k2 = k+1
       call eval_lagrange_poly_derivative(k2-k1+1, pref(:,:,k1:k2), &
            eta_dot_dpdn(:,:,k1:k2,2), &
            pref(:,:,k), ptp0)

       ! Horizontal velocity at initial time.
       k1 = k-1
       k2 = k
       v1 = half*(elem%derived%vstar(:,:,1,k1) + elem%derived%vstar(:,:,1,k2))
       v2 = half*(elem%derived%vstar(:,:,2,k1) + elem%derived%vstar(:,:,2,k2))

#define SL_ADVECTION_TRAJ_OLD
#ifdef SL_ADVECTION_TRAJ_OLD
       ! Reconstruct departure level coordinate at final time.
       p1r(:,:,k) = pref(:,:,k) + &
            half*dt*(eta_dot_dpdn(:,:,k,1) + eta_dot_dpdn(:,:,k,2) + &
            dt*(ptp0*eta_dot_dpdn(:,:,k,1) - grad(:,:,1)*v1 - grad(:,:,2)*v2))
#else
       ! Reconstruct eta_dot_dpdn over the time interval.
       eta_dot_dpdn(:,:,k,1) = &
            half*(eta_dot_dpdn(:,:,k,1) + eta_dot_dpdn(:,:,k,2) + &
            dt*(ptp0*eta_dot_dpdn(:,:,k,1) - grad(:,:,1)*v1 - grad(:,:,2)*v2))
#endif
    end do

    ! Reconstruct eta_dot_dpdn over the time interval.
#ifdef SL_ADVECTION_TRAJ_OLD
    eta_dot_dpdn(:,:,:,1) = (p1r - pref)/dt
#endif
    ! Boundary points are always 0.
    eta_dot_dpdn(:,:,1,1) = zero
    eta_dot_dpdn(:,:,nlevp,1) = zero

    dp_neg_min = reconstruct_and_limit_dp(elem%state%dp3d(:,:,:,tl%np1), &
         dt, dp_tol, eta_dot_dpdn(:,:,:,1), dprecon)
#ifndef NDEBUG
    if (dp_neg_min < dp_tol) then
       write(iulog, '(a,i7,i7,es11.4)') &
            'sl_advection: reconstruct_and_limit_dp (rank,ie) returned', &
            hybrid%par%rank, ie, dp_neg_min
    end if
#endif
  end subroutine calc_vertically_lagrangian_levels

  subroutine eval_lagrange_poly_derivative(n, xs, ys, xi, yp)
    integer, intent(in) :: n
    real(real_kind), intent(in) :: xs(np,np,n), ys(np,np,n), xi(np,np)
    real(real_kind), intent(out) :: yp(np,np)

    integer :: i, j, k
    real(real_kind) :: f(np,np), g(np,np), num(np,np)

    yp = zero
    do i = 1,n
       f = zero
       do j = 1,n
          if (j == i) cycle
          g = one
          do k = 1,n
             if (k == i) cycle
             if (k == j) then
                num = one
             else
                num = xi - xs(:,:,k)
             end if
             g = g*(num/(xs(:,:,i) - xs(:,:,k)))
          end do
          f = f + g
       end do
       yp = yp + ys(:,:,i)*f
    end do
  end subroutine eval_lagrange_poly_derivative

  subroutine calc_p(hvcoord, dp, p)
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: dp(np,np,nlev)
    real(real_kind), intent(out) :: p(np,np,nlevp)

    integer :: k

    p(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0
    do k = 1,nlev
       p(:,:,k+1) = p(:,:,k) + dp(:,:,k)
    end do
  end subroutine calc_p

  subroutine set_dp_tol(hvcoord, dp_tol)
    ! Pad by an amount ~ smallest level to keep the computed dp > 0.

    type (hvcoord_t), intent(in) :: hvcoord
    real(kind=real_kind), intent(out) :: dp_tol

    dp_tol = 10_real_kind*eps*minval(hvcoord%dp0)
  end subroutine set_dp_tol

  function reconstruct_and_limit_dp(dpref, dt, dp_tol, eta_dot_dpdn, dprecon) result(dp_neg_min)
    ! Move mass around in a column as needed to make dp nonnegative.

    real(kind=real_kind), intent(in) :: dpref(np,np,nlev), dt, dp_tol, eta_dot_dpdn(np,np,nlevp)
    real(kind=real_kind), intent(out) :: dprecon(np,np,nlev)

    integer :: k, i, j
    real(kind=real_kind) :: nmass, w(nlev), dp(nlev), dp_neg_min

    dp_neg_min = dp_tol ! < dp_tol if the limiter has to adjust dp
    do j = 1,np
       do i = 1,np
          nmass = zero
          do k = 1,nlev
             dp(k) = dpref(i,j,k) + dt*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k))
             if (dp(k) < dp_tol) then
                nmass = nmass + (dp(k) - dp_tol)
#ifndef NDEBUG
                dp_neg_min = min(dp_neg_min, dp(k))
#endif
                dp(k) = dp_tol
                w(k) = zero
             else
                w(k) = dp(k) - dp_tol
             end if
          end do
          ! Store the full update rather than reconstructing
          ! eta_dot_dpdn. See comment in sl_vertically_remap_tracers
          ! for more.
          dprecon(i,j,:) = dp
          if (nmass /= zero) dprecon(i,j,:) = dprecon(i,j,:) + nmass*(w/sum(w))
       end do
    end do
  end function reconstruct_and_limit_dp

  subroutine sl_vertically_remap_tracers(hybrid, elem, nets, nete, tl, dt_q)
    ! Remap the tracers after a tracer time step, in the case that the
    ! vertical remap time step for dynamics is shorter than the tracer
    ! time step.

    use control_mod, only: dt_tracer_factor
    use vertremap_base, only: remap1
    use parallel_mod, only: abortmp
    use kinds, only: iulog
    use perf_mod, only: t_startf, t_stopf

    type (hybrid_t), intent(in) :: hybrid
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
    real(kind=real_kind), intent(in) :: dt_q
    type (TimeLevel_t), intent(in) :: tl

    integer :: ie, i, j, k, q, n0_qdp, np1_qdp

    call t_startf('SLMM vertical remap')
    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)
    do ie = nets, nete
       ! divdp contains the reconstructed vertically Lagrangian level
       ! dp_star.
#ifndef NDEBUG
       if (any(elem(ie)%derived%divdp < zero)) then
          write(iulog,*) 'sl_vertically_remap_tracers> dp_star -ve: rank, ie, dp_tol', &
               hybrid%par%rank, ie, dp_tol
          do j = 1,np
             do i = 1,np
                if (any(elem(ie)%derived%divdp(i,j,:) < zero)) then
                   write(iulog,*) 'i,j,dp_star(i,j,:)', i, j, elem(ie)%derived%divdp(i,j,:)
                   call abortmp('sl_vertically_remap_tracers> -ve dp_star')
                end if
             end do
          end do
       end if
#endif
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp), np, qsize, elem(ie)%derived%divdp, &
            elem(ie)%state%dp3d(:,:,:,tl%np1),vert_remap_q_alg)
       do q = 1,qsize
          elem(ie)%state%Q(:,:,:,q) = elem(ie)%state%Qdp(:,:,:,q,np1_qdp)/ &
                                      elem(ie)%state%dp3d(:,:,:,tl%np1)
       enddo
    end do
    call t_stopf('SLMM vertical remap')
  end subroutine sl_vertically_remap_tracers

  function test_lagrange() result(nerr)
    use kinds, only: rt => real_kind

    real(rt), parameter :: xs(3) = (/-one, zero, half/)
    integer, parameter :: n = 3, ntrial = 10

    real(rt) :: a, b, c, x, y1, y2, alpha, ys(3), xsi(np,np,n), ysi(np,np,n), &
         xi(np,np), y2i(np,np)
    integer :: i, j, trial, nerr

    nerr = 0
    a = -half; b = 0.3_rt; c = 1.7_rt

    do i = 1,3
       x = xs(i)
       ys(i) = (a*x + b)*x + c
    end do

    do trial = 0,ntrial
       alpha = real(trial,rt)/ntrial
       x = (1-alpha)*(-2_rt) + alpha*2_rt
       y1 = 2*a*x + b
       do j = 1,np
          do i = 1,np
             xsi(i,j,:) = xs
             ysi(i,j,:) = ys
             xi(i,j) = x
          end do
       end do
       call eval_lagrange_poly_derivative(n, xsi, ysi, xi, y2i)
       if (abs(y2i(np,np) - y1) > 1d-14*abs(y1)) nerr = nerr + 1
    end do
  end function test_lagrange

  function test_reconstruct_and_limit_dp() result(nerr)
    use physical_constants, only: p0
    real(real_kind), parameter :: dt = 1800_real_kind, dp_tol = (p0/nlev)*eps, &
         tol = 1.e3_real_kind*eps

    real(real_kind) :: dpref(np,np,nlev), dpfin(np,np,nlev,2), eta_dot_dpdn(np,np,nlevp), tmp
    integer :: nerr, i, j, k

    eta_dot_dpdn(1,1,nlevp) = zero
    do k = 1,nlev
       dpref(1,1,k) = k
       eta_dot_dpdn(1,1,k) = (-one)**k*0.1_real_kind*(nlev-k)
    end do
    tmp = p0/sum(dpref(1,1,:))
    dpref(1,1,:) = dpref(1,1,:)*tmp
    eta_dot_dpdn(1,1,:) = eta_dot_dpdn(1,1,:)*tmp
    eta_dot_dpdn(1,1,1) = zero

    do j = 1,np
       do i = 1,np
          dpref(i,j,:) = dpref(1,1,:)
          eta_dot_dpdn(i,j,:) = eta_dot_dpdn(1,1,:)
       end do
    end do

    dpfin(:,:,:,1) = dpref + dt*(eta_dot_dpdn(:,:,2:) - eta_dot_dpdn(:,:,:nlev))

    tmp = reconstruct_and_limit_dp(dpref, dt, dp_tol, eta_dot_dpdn, dpfin(:,:,:,2))
    
    nerr = 0
    do j = 1,np
       do i = 1,np
          ! mass conservation
          tmp = sum(dpref(i,j,:))
          if (abs(sum(dpfin(i,j,:,1)) - tmp) > tol*tmp) nerr = nerr + 1
          if (abs(sum(dpfin(i,j,:,2)) - tmp) > tol*tmp) nerr = nerr + 1
          ! limiter needs to be active
          if (minval(dpfin(i,j,:,1)) >= zero) nerr = nerr + 1
          ! limiter succeeded
          if (minval(dpfin(i,j,:,2)) < zero) nerr = nerr + 1
       end do
    end do
  end function test_reconstruct_and_limit_dp

  subroutine sl_unittest(par)
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par

    integer :: nerr

    nerr = 0
    nerr = nerr + test_lagrange()
    nerr = nerr + test_reconstruct_and_limit_dp()

    if (nerr > 0 .and. par%masterproc) write(iulog,'(a,i3)') 'sl_unittest FAIL', nerr
  end subroutine sl_unittest

end module sl_advection
