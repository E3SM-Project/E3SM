#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

module sl_advection
  use kinds, only              : real_kind, int_kind
  use dimensions_mod, only     : nlev, nlevp, np, qsize, qsize_d
  use derivative_mod, only     : derivative_t, gradient_sphere, divergence_sphere, ugradv_sphere
  use element_mod, only        : element_t
  use hybvcoord_mod, only      : hvcoord_t
  use time_mod, only           : TimeLevel_t, TimeLevel_Qdp
  use control_mod, only        : integration, test_case, hypervis_order, transport_alg, &
       &                         limiter_option, vert_remap_q_alg, semi_lagrange_diagnostics
  use edge_mod, only           : edgevpack_nlyr, edgevunpack_nlyr, edge_g
  use edgetype_mod, only       : EdgeDescriptor_t, EdgeBuffer_t
  use hybrid_mod, only         : hybrid_t
  use bndry_mod, only          : bndry_exchangev
  use perf_mod, only           : t_startf, t_stopf, t_barrierf ! _EXTERNAL
  use parallel_mod, only       : abortmp, parallel_t
  use coordinate_systems_mod, only : cartesian3D_t
#ifdef HOMME_ENABLE_COMPOSE
  use compose_mod
#endif

  implicit none

  private

  ! Constants
  real(real_kind), parameter :: zero = 0.0_real_kind, fourth = 0.25_real_kind, &
       half = 0.5_real_kind, one = 1.0_real_kind, two = 2.0_real_kind, &
       eps = epsilon(1.0_real_kind)

  ! Configuration.
  logical :: is_sphere, enhanced_trajectory
  integer :: dep_points_ndim

  ! For use in make_positive. Set at initialization to a function of hvcoord%dp0.
  real(kind=real_kind) :: dp_tol, deta_tol

  public :: prim_advec_tracers_observe_velocity_ALE, prim_advec_tracers_remap_ALE, &
       &    sl_init1, sl_vertically_remap_tracers, sl_unittest

  ! For testing
  public :: calc_trajectory, dep_points_all, sphere2cart

  ! For C++
  public :: sl_get_params

  ! Barrier for performance analysis. Should be false in production runs.
  logical, parameter :: barrier = .false.

  ! Bounds for shape preservation.
  real(kind=real_kind), dimension(:,:,:,:,:), allocatable :: minq, maxq  ! (np,np,nlev,qsize,nelemd)

  ! Trajectory velocity data.
  real(kind=real_kind), dimension(:,:,:,:,:), allocatable :: vnode, vdep ! (ndim,np,np,nlev,nelemd)
  real(kind=real_kind), allocatable :: dep_points_all(:,:,:,:,:)         ! (ndim,np,np,nlev,nelemd)

  type :: velocity_record_t
     integer :: nvel
     ! Times to which velocity slots correspond, in reference time [0,dtf].
     real(kind=real_kind), allocatable :: t_vel(:) ! 1:nvel
     ! For n = 1:dtf, obs_slots(n,:) = [slot1, slot2], -1 if unused. These are
     ! the slots to which velocity sample n contributes. obs_slots(dtf,:) is
     ! always -1.
     integer, allocatable :: obs_slots(:,:)
     ! obs_wts(n,:) = [wt1, wt2], 0 if unused.
     real(kind=real_kind), allocatable :: obs_wts(:,:)
     ! Substep end point n in 0:nsub uses velocity slots run_step(n),
     ! run_step(n)-1.
     integer, allocatable :: run_step(:)
     ! Store state%v and state%dp3d at t_vel points.
     real(kind=real_kind), allocatable :: vel(:,:,:,:,:,:) ! (np,np,2,nlev,nss,nelemd)
     real(kind=real_kind), allocatable :: dp (:,:,  :,:,:) ! (np,np,  nlev,nss,nelemd)
  end type velocity_record_t

  type(velocity_record_t) :: vrec
  
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
         nu_q, semi_lagrange_hv_q, semi_lagrange_cdr_check, semi_lagrange_trajectory_nsubstep, &
         semi_lagrange_trajectory_nvelocity, geometry, dt_remap_factor, dt_tracer_factor, &
         semi_lagrange_halo
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
       call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test, &
            independent_time_steps)
       is_sphere = trim(geometry) /= 'plane'
       enhanced_trajectory = semi_lagrange_trajectory_nsubstep > 0
       dep_points_ndim = 3
       if (enhanced_trajectory .and. independent_time_steps) dep_points_ndim = 4
       nslots = nlev*qsize
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
       allocate(minq(np,np,nlev,qsize,size(elem)), maxq(np,np,nlev,qsize,size(elem)), &
            &   dep_points_all(dep_points_ndim,np,np,nlev,size(elem)))
       if (enhanced_trajectory) then
          allocate(vnode(dep_points_ndim,np,np,nlev,size(elem)), &
               &   vdep (dep_points_ndim,np,np,nlev,size(elem)))
       end if
       call init_velocity_record(size(elem), dt_tracer_factor, dt_remap_factor, &
            semi_lagrange_trajectory_nsubstep, semi_lagrange_trajectory_nvelocity, &
            vrec, i)
       dp_tol = -one
       deta_tol = -one
       if (par%masterproc) then
          if (nu_q > 0 .and. semi_lagrange_hv_q > 0) then
             write(iulog,'(a,es13.4,i3)') 'COMPOSE> use HV; nu_q, hv_q:', &
                  nu_q, semi_lagrange_hv_q
          end if
          if (enhanced_trajectory) then
             write(iulog,'(a,i3,i3,i3)') &
                  'COMPOSE> dt_tracer_factor, dt_remap_factor, halo:', &
                  dt_tracer_factor, dt_remap_factor, semi_lagrange_halo
             write(iulog,'(a,i3,i3)') &
                  'COMPOSE> use enhanced trajectory; nsub, nvel:', &
                  semi_lagrange_trajectory_nsubstep, vrec%nvel
          end if
       end if
    endif
    call t_stopf('sl_init1')
#endif
  end subroutine sl_init1

  subroutine sl_get_params(nu_q_out, hv_scaling, hv_q, hv_subcycle_q, limiter_option_out, &
       cdr_check, geometry_type, trajectory_nsubstep) bind(c)
    use control_mod, only: semi_lagrange_hv_q, hypervis_subcycle_q, semi_lagrange_cdr_check, &
         nu_q, hypervis_scaling, limiter_option, geometry, semi_lagrange_trajectory_nsubstep
    use iso_c_binding, only: c_int, c_double

    real(c_double), intent(out) :: nu_q_out, hv_scaling
    integer(c_int), intent(out) :: hv_q, hv_subcycle_q, limiter_option_out, cdr_check, &
         geometry_type, trajectory_nsubstep

    nu_q_out = nu_q
    hv_scaling = hypervis_scaling
    hv_q = semi_lagrange_hv_q
    hv_subcycle_q = hypervis_subcycle_q
    limiter_option_out = limiter_option
    cdr_check = 0
    if (semi_lagrange_cdr_check) cdr_check = 1
    geometry_type = 0 ! sphere
    if (trim(geometry) == "plane") geometry_type = 1
    trajectory_nsubstep = semi_lagrange_trajectory_nsubstep
  end subroutine sl_get_params

  subroutine init_velocity_record(nelemd, dtf, drf_param, nsub, nvel_param, v, error)
    integer, intent(in) :: nelemd, dtf, drf_param, nsub, nvel_param
    type (velocity_record_t), intent(out) :: v
    integer, intent(out) :: error

    real(kind=real_kind) :: t_avail(0:dtf), time
    integer :: nvel, drf, navail, n, i, iav

    error = 0
    drf = drf_param
    if (drf == 0) drf = 1 ! drf = 0 if vertically Eulerian
    nvel = nvel_param
    if (nvel == -1) nvel = 2 + ((nsub-1) / 2)
    nvel = min(nvel, nsub+1)
    navail = dtf/drf + 1
    nvel = min(nvel, navail)

    ! nsub <= 1: No substepping.
    ! nvel <= 2: Save velocity only at endpoints, as always occurs.
    if (nsub <= 1 .or. nvel <= 2) then
       v%nvel = 2
       return
    end if

    v%nvel = nvel
    allocate(v%t_vel(nvel), v%obs_slots(dtf,2), v%obs_wts(dtf,2), v%run_step(0:nsub), &
         &   v%vel(np,np,2,nlev,2:nvel-1,nelemd), v%dp(np,np,nlev,2:nvel-1,nelemd))

    ! Times at which velocity data are available.
    t_avail(0) = 0
    i = 1
    do n = 1, dtf
       if (modulo(n, drf) == 0) then
          t_avail(i) = n
          i = i + 1
       end if
    end do
    if (i /= navail) error = 1

    ! Times to which we associate velocity data.
    do n = 1, nvel
       if (modulo((n-1)*dtf, nvel-1) == 0) then
          ! Cast integer values at end of calculation.
          v%t_vel(n) = ((n-1)*dtf)/(nvel-1)
       else
          v%t_vel(n) = real((n-1)*dtf, real_kind)/(nvel-1)
       end if
    end do

    ! Build the tables mapping n in 1:dtf to velocity slots to accumulate into.
    do n = 1, dtf-1
       v%obs_slots(n,:) = -1
       v%obs_wts(n,:) = 0
       if (modulo(n, drf) /= 0) cycle
       time = n
       do i = 1, navail-1
          if (time == t_avail(i)) exit
       end do
       iav = i
       if (iav > navail-1) error = 2
       do i = 2, nvel-1
          if (t_avail(iav-1) < v%t_vel(i) .and. time > v%t_vel(i)) then
             v%obs_slots(n,1) = i
             v%obs_wts(n,1) = (v%t_vel(i) - t_avail(iav-1))/(t_avail(iav) - t_avail(iav-1))
          end if
          if (time <= v%t_vel(i) .and. t_avail(iav+1) > v%t_vel(i)) then
             v%obs_slots(n,2) = i
             v%obs_wts(n,2) = (t_avail(iav+1) - v%t_vel(i))/(t_avail(iav+1) - t_avail(iav))
          end if
       end do
    end do
    v%obs_slots(dtf,:) = -1
    v%obs_wts(dtf,:) = 0

    ! Build table mapping n to interval to use. The trajectories go backward in
    ! time, and this table reflects that.
    v%run_step(0) = nvel
    v%run_step(nsub) = 2
    do n = 1, nsub-1
       time = real((nsub-n)*dtf, real_kind)/nsub
       do i = 1, nvel-1
          if (v%t_vel(i) <= time .and. time <= v%t_vel(i+1)) exit
       end do
       if (i > nvel-1) error = 4
       v%run_step(n) = i+1
    end do
  end subroutine init_velocity_record

  subroutine prim_advec_tracers_observe_velocity_ALE(elem, tl, n, nets, nete)
    use control_mod, only: dt_remap_factor
    
    type (element_t)     , intent(inout) :: elem(:)
    type (TimeLevel_t)   , intent(in   ) :: tl
    integer              , intent(in   ) :: n       ! step in 1:dt_tracer_factor
    integer              , intent(in   ) :: nets
    integer              , intent(in   ) :: nete

    integer :: nstore, islot, slot, k, ie

    if (vrec%nvel == 2) return

    if (n == dt_remap_factor .or. (dt_remap_factor == 0 .and. n == 1)) then
       ! First observation of the tracer time step: zero accumulated quantities.
       do ie = nets, nete
          do slot = 2, vrec%nvel-1
             vrec%vel(:,:,:,:,slot,ie) = 0
             vrec%dp (:,:,  :,slot,ie) = 0
          end do
       end do
    end if

    do islot = 1, 2
       slot = vrec%obs_slots(n,islot)
       if (slot == -1) cycle
       do ie = nets, nete
          do k = 1, nlev
             vrec%vel(:,:,:,k,slot,ie) = vrec%vel(:,:,:,k,slot,ie) + &
                  vrec%obs_wts(n,islot) * elem(ie)%state%v(:,:,:,k,tl%np1)
             vrec%dp (:,:,  k,slot,ie) = vrec%dp (:,:,  k,slot,ie) + &
                  vrec%obs_wts(n,islot) * elem(ie)%state%dp3d(:,:,k,tl%np1)
          end do
       end do
    end do
  end subroutine prim_advec_tracers_observe_velocity_ALE

  subroutine prim_advec_tracers_remap_ALE(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete)
    use coordinate_systems_mod, only : cartesian3D_t, cartesian2D_t
    use dimensions_mod,         only : max_neigh_edges
    use interpolate_mod,        only : interpolate_tracers, minmax_tracers
    use control_mod,            only : dt_tracer_factor, nu_q, transport_alg, semi_lagrange_hv_q, &
         semi_lagrange_cdr_alg, semi_lagrange_cdr_check, semi_lagrange_trajectory_nsubstep
    ! For DCMIP16 supercell test case.
    use control_mod,            only : dcmip16_mu_q
    use prim_advection_base,    only : advance_physical_vis
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

    integer :: i,j,k,l,n,q,ie,n0_qdp,np1_qdp
    integer :: scalar_q_bounds, info
    logical :: slmm, cisl, qos, sl_test, independent_time_steps
    logical(kind=c_bool) :: h2d, d2h

#ifdef HOMME_ENABLE_COMPOSE
    call t_barrierf('Prim_Advec_Tracers_remap_ALE', hybrid%par%comm)
    call t_startf('Prim_Advec_Tracers_remap_ALE')

    call sl_parse_transport_alg(transport_alg, slmm, cisl, qos, sl_test, independent_time_steps)
    h2d = .true.
    d2h = compose_d2h .or. h2d
    h2d = compose_h2d .or. h2d
    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)

    if (enhanced_trajectory) then
       call calc_enhanced_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, nets, nete, &
            semi_lagrange_trajectory_nsubstep, independent_time_steps)
    else
       call calc_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, &
            independent_time_steps, nets, nete)
    end if

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
    ! dp_coupling in EAM. Thus, we must take care to protect threaded access. In
    ! the following, "No barrier needed" comments justify why a barrier isn't
    ! needed.    
    ! No barrier needed: ale_rkdss has a horiz thread barrier at the end.
    call slmm_csl(nets, nete, dep_points_all, dep_points_ndim, minq, maxq, info)
    ! No barrier needed: slmm_csl has a horiz thread barrier at the end.
    if (info /= 0) then
       call write_velocity_data(elem, nets, nete, hybrid, dt, tl)
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
       ! Barrier needed to protect edge_g buffers used in CEDR.
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
          call ALE_departure_from_gll(dep_points_all(:,:,:,k,ie), dep_points_ndim, &
               elem(ie)%derived%vstar(:,:,:,k), elem(ie), dt, normalize=is_sphere)
       end do
    end do
    call t_stopf('SLMM_v2x')
  end subroutine calc_trajectory

  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE ALE_RKDSS
  ! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012
  !
  ! DESCRIPTION: ! create a runge kutta taylor series mixture to calculate the departure grid
  !
  ! CALLS:
  ! INPUT:
  !
  ! OUTPUT:
  !-----------------------------------------------------------------------------------!

  ! this will calculate the velocity at time t+1/2  along the trajectory s(t) given the velocities
  ! at the GLL points at time t and t+1 using a second order time accurate formulation.
  subroutine ALE_RKdss(elem, nets, nete, hy, deriv, dt, tl, independent_time_steps)
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
    !     x(t+1) = x(t) + 1/2 ( U(x*(t+1),t+1) + U(x(t),t) ) dt
    ! apply taylor series:
    !     U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
    !
    ! x(t+1) = x(t) + 1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) ) dt
    ! (x(t+1) - x(t)) / dt = 1/2 ( U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1) + U(x(t),t) )
    !                      = 1/2 ( U(x(t),t+1) + U(x(t),t) + (x*(t+1)-x(t)) gradU(x(t),t+1) )
    !                      = 1/2 ( U(x(t),t+1) + U(x(t),t) + U(x(t),t) dt gradU(x(t),t+1) )
    !
    ! => (x(t+1)-x(t))/dt = 1/2 (U(x(t),t+1) + U(x(t),t) + dt U(x(t),t) gradU(x(t),t+1))
    !
    ! suppose dt = -ts (we go backward)
    !     (x(t-ts)-x(t))/-ts =  1/2 (U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)
    !     x(t-ts) = x(t)) - ts * [1/2 (U(x(t),t-ts)+U(x(t),t)) - ts 1/2 U(x(t),t) gradU(x(t),t-ts)]

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
          call edgeVunpack_nlyr(edge_g,elem(ie)%desc,elem(ie)%derived%divdp,nlev,2*nlev,nlyr)
       end if
    end do

#if (defined HORIZ_OPENMP)
    !$omp barrier
#endif
  end subroutine ALE_RKdss

  subroutine write_velocity_data(elem, nets, nete, hy, dt, tl)
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
  subroutine ALE_departure_from_gll(acart, ndim, vstar, elem, dt, normalize)
    use physical_constants,     only : scale_factor
    use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t
    use time_mod,               only : timelevel_t
    use element_mod,            only : element_t
    use kinds,                  only : real_kind
    use dimensions_mod,         only : np

    implicit none

    integer                 ,intent(in)   :: ndim
    real (kind=real_kind)   ,intent(out)  :: acart(ndim,np,np)
    real (kind=real_kind)   ,intent(in)   :: vstar(np,np,2)
    type (element_t)        ,intent(in)   :: elem
    real (kind=real_kind)   ,intent(in)   :: dt
    logical                 ,intent(in)   :: normalize

    integer                               :: i,j, d
    type (cartesian3D_t)                  :: c3d
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
          call sphere2cart(elem%spherep(i,j), c3d)
          acart(1,i,j) = c3d%x - dt*uxyz(i,j,1)/scale_factor
          acart(2,i,j) = c3d%y - dt*uxyz(i,j,2)/scale_factor
          acart(3,i,j) = c3d%z - dt*uxyz(i,j,3)/scale_factor
          if (normalize) then
             norm = sqrt(acart(1,i,j)*acart(1,i,j) + acart(2,i,j)*acart(2,i,j) + &
                  acart(3,i,j)*acart(3,i,j))
             acart(1:3,i,j) = acart(1:3,i,j)/norm
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
    !                     dt (w_x(p0,t1) u(p0,t0) + w_z(p0,t1) w(p0,t0)))
    !          + O(dt^3).                                                   (*)
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
    !                + O(dt^3).
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
    ! (ps1 - ps0); dp3d already accounts for this term w.r.t. derived%dp. Recall
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
    p1r(:,:,1) = zero
    p1r(:,:,nlevp) = zero
    eta_dot_dpdn(:,:,:,1) = (p1r - pref)/dt
#endif
    ! Boundary points are always 0.
    eta_dot_dpdn(:,:,1,1) = zero
    eta_dot_dpdn(:,:,nlevp,1) = zero

    dp_neg_min = reconstruct_and_limit_dp(elem%state%dp3d(:,:,:,tl%np1), &
         dt, dp_tol, eta_dot_dpdn(:,:,:,1), dprecon)
#if 0
    if (dp_neg_min < dp_tol) then
       write(iulog, '(a,i7,i7,es13.4)') &
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

  subroutine calc_enhanced_trajectory(elem, deriv, hvcoord, hybrid, dt, tl, &
       nets, nete, nsubstep, independent_time_steps)
    ! Top-level routine for new enhanced trajectory method. This new method
    ! permits multiple substeps, optionally using more reference-grid velocity
    ! snapshots.

    use reduction_mod, only: ParallelSum
    use kinds, only: iulog
    
    type (element_t), intent(inout) :: elem(:)
    type (derivative_t), intent(in) :: deriv
    type (hvcoord_t), intent(in) :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    real(real_kind), intent(in) :: dt
    type (TimeLevel_t), intent(in) :: tl
    integer, intent(in) :: nets, nete, nsubstep
    logical, intent(in) :: independent_time_steps

#ifdef HOMME_ENABLE_COMPOSE
    integer :: step, ie, info, limiter_active_count
    real(real_kind) :: alpha(2), dtsub

    call t_startf('SLMM_trajectory')

    call slmm_set_hvcoord(hvcoord%etai(1), hvcoord%etai(nlevp), hvcoord%etam)

    ! Set dep_points_all to level-midpoint arrival points.
    call init_dep_points_all(elem, hvcoord, nets, nete, independent_time_steps)

    limiter_active_count = 0
    dtsub = dt / nsubstep
    do step = 1, nsubstep
       ! Fill vnode.
       if (vrec%nvel == 2) then
          alpha(1) = real(nsubstep - step    , real_kind)/nsubstep
          alpha(2) = real(nsubstep - step + 1, real_kind)/nsubstep
          do ie = nets, nete
             call calc_nodal_velocities(elem(ie), deriv, hvcoord, tl, &
                  independent_time_steps, dtsub, alpha, &
                  elem(ie)%derived%vstar, elem(ie)%derived%dp, &
                  elem(ie)%state%v(:,:,:,:,tl%np1), elem(ie)%state%dp3d(:,:,:,tl%np1), &
                  vnode(:,:,:,:,ie))
          end do
       else
          call calc_nodal_velocities_using_vrec(elem, deriv, hvcoord, tl, &
               independent_time_steps, dtsub, nsubstep, step, nets, nete)
       end if

       call dss_vnode(elem, nets, nete, hybrid, vnode)

       if (step == 1) then
          call update_dep_points_all(independent_time_steps, dtsub, nets, nete, vnode)
       else
          ! Fill vdep.
          call slmm_calc_v_departure(nets, nete, step, dtsub, dep_points_all, &
               &                     dep_points_ndim, vnode, vdep, info)

          ! Using vdep, update dep_points_all to departure points.
          call update_dep_points_all(independent_time_steps, dtsub, nets, nete, vdep)
       end if
    end do

    if (independent_time_steps) then
       call interp_departure_points_to_floating_level_midpoints( &
            elem, nets, nete, tl, hvcoord, dep_points_all, limiter_active_count)
       ! Not needed in practice. Corner cases will be cleaned up by dss_Qdp.
       !call dss_divdp(elem, nets, nete, hybrid)
       if (iand(semi_lagrange_diagnostics, 1) /= 0) then
          limiter_active_count = ParallelSum(limiter_active_count, hybrid)
          if (limiter_active_count > 0 .and. hybrid%masterthread) then
             write(iulog, '(a,i11)') 'COMPOSE> limiter_active_count', &
                  limiter_active_count
          end if
       end if
    end if

    call t_stopf('SLMM_trajectory')
#endif
  end subroutine calc_enhanced_trajectory

  subroutine init_dep_points_all(elem, hvcoord, nets, nete, independent_time_steps)
    ! Initialize dep_points_all to the Eulerian coordinates.

    type (element_t), intent(inout) :: elem(:)
    type (hvcoord_t), intent(in) :: hvcoord
    integer, intent(in) :: nets, nete
    logical, intent(in) :: independent_time_steps

    type (cartesian3D_t) :: c3d
    integer :: ie, i, j, k
    
    do ie = nets, nete
       do j = 1, np
          do i = 1, np
             call sphere2cart(elem(ie)%spherep(i,j), c3d)
             dep_points_all(1,i,j,1,ie) = c3d%x
             dep_points_all(2,i,j,1,ie) = c3d%y
             dep_points_all(3,i,j,1,ie) = c3d%z
             do k = 2, nlev
                dep_points_all(1:3,i,j,k,ie) = dep_points_all(1:3,i,j,1,ie)
             end do
             if (independent_time_steps) then
                do k = 1, nlev
                   dep_points_all(4,i,j,k,ie) = hvcoord%etam(k)
                end do
             end if
          end do
       end do
    end do
  end subroutine init_dep_points_all

  subroutine calc_nodal_velocities_using_vrec(elem, deriv, hvcoord, tl, &
       independent_time_steps, dtsub, nsubstep, step, nets, nete)

    ! Wrapper to calc_nodal_velocities that orchestrates the use of the various
    ! sources of velocity data.

    type (element_t), intent(in) :: elem(:)
    type (derivative_t), intent(in) :: deriv
    type (hvcoord_t), intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    logical, intent(in) :: independent_time_steps
    real(real_kind), intent(in) :: dtsub
    integer, intent(in) :: nsubstep, step, nets, nete
    
    integer :: ie, i, k, os
    real(real_kind) :: time, alpha(2), vs(np,np,2,nlev,3), dps(np,np,nlev,3)

    do ie = nets, nete
       do i = 1, 2
          k = nsubstep - step + (i-1)
          time = (k*vrec%t_vel(vrec%nvel))/nsubstep
          os = i-1
          k = vrec%run_step(step+1-i)
          if (k == 2) then
             vs(:,:,:,:,os+1) = elem(ie)%derived%vstar
             dps(:,:,:,os+1) = elem(ie)%derived%dp
          else
             vs(:,:,:,:,os+1) = vrec%vel(:,:,:,:,k-1,ie)
             dps(:,:,:,os+1) = vrec%dp(:,:,:,k-1,ie)
          end if
          if (k == vrec%nvel) then
             vs(:,:,:,:,os+2) = elem(ie)%state%v(:,:,:,:,tl%np1)
             dps(:,:,:,os+2) = elem(ie)%state%dp3d(:,:,:,tl%np1)
          else
             vs(:,:,:,:,os+2) = vrec%vel(:,:,:,:,k,ie)
             dps(:,:,:,os+2) = vrec%dp(:,:,:,k,ie)
          end if
          alpha(1) = (vrec%t_vel(k) - time)/(vrec%t_vel(k) - vrec%t_vel(k-1))
          alpha(2) = 1 - alpha(1)
          vs(:,:,:,:,os+1) = alpha(1)*vs(:,:,:,:,os+1) + alpha(2)*vs(:,:,:,:,os+2)
          dps(:,:,:, os+1) = alpha(1)*dps(:,:,:, os+1) + alpha(2)*dps(:,:,:, os+2)
       end do
       alpha(1) = 0
       alpha(2) = 1
       call calc_nodal_velocities(elem(ie), deriv, hvcoord, tl, &
            independent_time_steps, dtsub, alpha, &
            vs(:,:,:,:,1), dps(:,:,:,1), vs(:,:,:,:,2), dps(:,:,:,2), &
            vnode(:,:,:,:,ie))
    end do    
  end subroutine calc_nodal_velocities_using_vrec

  subroutine calc_nodal_velocities(elem, deriv, hvcoord, tl, &
       independent_time_steps, dtsub, alpha, v1, dp1, v2, dp2, vnode)
    ! Evaluate a formula to provide an estimate of nodal velocities that
    ! are use to create a 2nd-order update to the trajectory. The
    ! fundamental formula for the update in position p from arrival point
    ! p1 to departure point p0 is
    !     p0 = p1 - dt/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).
    ! Here we compute the velocity estimate at the nodes:
    !     1/2 (v(p1,t0) + v(p1,t1) - dt v(p1,t1) grad v(p1,t0)).

    type (element_t), intent(in) :: elem
    type (derivative_t), intent(in) :: deriv
    type (hvcoord_t), intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    logical, intent(in) :: independent_time_steps
    real(real_kind), intent(in) :: dtsub, alpha(2)
    real(real_kind), dimension(np,np,2,nlev), intent(in) :: v1, v2
    real(real_kind), dimension(np,np,nlev), intent(in) :: dp1, dp2
    real(real_kind), intent(out) :: vnode(:,:,:,:)

    real(real_kind) :: vsph(np,np,2,nlev,2), eta_dot(np,np,nlevp,2)
    integer :: t

    if (independent_time_steps) then
       call calc_eta_dot_ref_mid(elem, deriv, tl, hvcoord, alpha, &
            &                    v1, dp1, v2, dp2, eta_dot)
    else
       eta_dot = zero
    end if

    ! Collect the horizontal nodal velocities. v1,2 are on Eulerian levels. v1
    ! is from time t1 < t2.
    do t = 1, 2
       vsph(:,:,:,:,t) = (1 - alpha(t))*v1 + alpha(t)*v2
    end do

    ! Given the vertical and horizontal nodal velocities at time
    ! endpoints, evaluate the velocity estimate formula, providing the
    ! final horizontal and vertical velocity estimates at midpoint nodes.
    call calc_vel_horiz_formula_node_ref_mid( &
         &  elem, deriv, hvcoord, dtsub, vsph, eta_dot, vnode)
    if (independent_time_steps) then
       call calc_eta_dot_formula_node_ref_mid( &
            elem, deriv, hvcoord, dtsub, vsph, eta_dot, vnode)
    end if
  end subroutine calc_nodal_velocities

  subroutine calc_eta_dot_ref_mid(elem, deriv, tl, hvcoord, alpha, v1, dp1, v2, dp2, eta_dot)
    ! Compute eta_dot at midpoint nodes at the start and end of the substep.

    type (element_t), intent(in) :: elem
    type (derivative_t), intent(in) :: deriv
    type (TimeLevel_t), intent(in) :: tl
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: alpha(2)
    real(real_kind), dimension(np,np,2,nlev), intent(in) :: v1, v2
    real(real_kind), dimension(np,np,nlev), intent(in) :: dp1, dp2
    real(real_kind), intent(out) :: eta_dot(np,np,nlevp,2)

    real(real_kind) :: vdp(np,np,2), w1(np,np)
    integer :: t, k, d

    do t = 1,2
       ! eta_dot_dpdn at interface nodes.
       eta_dot(:,:,1,t) = zero
       do k = 1,nlev
          do d = 1,2
             vdp(:,:,d) = (1 - alpha(t))*v1(:,:,d,k)*dp1(:,:,k) + &
                  &            alpha(t) *v2(:,:,d,k)*dp2(:,:,k)
          end do
          w1 = divergence_sphere(vdp, deriv, elem)
          eta_dot(:,:,k+1,t) = eta_dot(:,:,k,t) + w1
       end do
       w1 = eta_dot(:,:,nlevp,t)
       eta_dot(:,:,nlevp,t) = zero
       do k = 2,nlev
          eta_dot(:,:,k,t) = hvcoord%hybi(k)*w1 - eta_dot(:,:,k,t)
       end do
       ! Transform eta_dot_dpdn at interfaces to eta_dot at midpoints using the
       ! formula
       !     eta_dot = eta_dot_dpdn/(A_eta p0 + B_eta ps)
       !            a= eta_dot_dpdn diff(eta)/(diff(A) p0 + diff(B) ps).
       !   Compute ps.
       w1 = hvcoord%hyai(1)*hvcoord%ps0 + &
            &    (1 - alpha(t))*sum(dp1, 3) + &
            &         alpha(t) *sum(dp2, 3)
       do k = 1,nlev
          eta_dot(:,:,k,t) = half*(eta_dot(:,:,k,t) + eta_dot(:,:,k+1,t)) &
               &             * (hvcoord%etai(k+1) - hvcoord%etai(k)) &
               &             / (  (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 &
               &                + (hvcoord%hybi(k+1) - hvcoord%hybi(k))*w1)
       end do
    end do
  end subroutine calc_eta_dot_ref_mid

  subroutine calc_vel_horiz_formula_node_ref_mid( &
       elem, deriv, hvcoord, dtsub, vsph, eta_dot, vnode)

    type (element_t), intent(in) :: elem
    type (derivative_t), intent(in) :: deriv
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: dtsub, vsph(np,np,2,nlev,2), eta_dot(np,np,nlevp,2)
    real(real_kind), intent(inout) :: vnode(:,:,:,:)

    integer, parameter :: t0 = 1, t1 = 2

    real(real_kind) :: vfsph(np,np,2), w1(np,np), w2(np,np), w3(np,np,3)
    integer :: k, d, i, k1, k2

    do k = 1, nlev
       ! Horizontal terms.
       vfsph = ugradv_sphere(vsph(:,:,:,k,t1), vsph(:,:,:,k,t0), deriv, elem)
       vfsph = vsph(:,:,:,k,t0) + vsph(:,:,:,k,t1) - dtsub*vfsph
       ! Vertical term.
       do d = 1, 2 ! horiz vel dims
          if (k == 1 .or. k == nlev) then
             if (k == 1) then
                k1 = 1; k2 = 2
             else
                k1 = nlev-1; k2 = nlev
             end if
             w1 = (vsph(:,:,d,k2,t0) - vsph(:,:,d,k1,t0)) / &
                  (hvcoord%etam(k2) - hvcoord%etam(k1))
          else
             do i = 1, 3
                w3(:,:,i) = hvcoord%etam(k-2+i) ! interp support
             end do
             w2 = hvcoord%etam(k) ! derivative at this eta value
             call eval_lagrange_poly_derivative(3, w3, vsph(:,:,d,k-1:k+1,t0), w2, w1)
          end if
          vfsph(:,:,d) = vfsph(:,:,d) - dtsub*eta_dot(:,:,k,t1)*w1
       end do
       ! Finish the formula.
       vfsph = half*vfsph
       ! Transform to Cartesian.
       do d = 1, 3
          vnode(d,:,:,k) = sum(elem%vec_sphere2cart(:,:,d,:)*vfsph, 3)
       end do
    end do
  end subroutine calc_vel_horiz_formula_node_ref_mid

  subroutine calc_eta_dot_formula_node_ref_mid( &
       elem, deriv, hvcoord, dtsub, vsph, eta_dot, vnode)

    type (element_t), intent(in) :: elem
    type (derivative_t), intent(in) :: deriv
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: dtsub, vsph(np,np,2,nlev,2), eta_dot(np,np,nlevp,2)
    real(real_kind), intent(inout) :: vnode(:,:,:,:)

    integer, parameter :: t0 = 1, t1 = 2
    
    real(real_kind) :: vfsph(np,np,2), w1(np,np), w2(np,np), w3(np,np,3), w4(np,np,3)
    integer :: k, d, i, k1, k2

    do k = 1, nlev
       w2 = hvcoord%etam(k)
       if (k == 1 .or. k == nlev) then
          if (k == 1) then
             w3(:,:,1) = hvcoord%etai(1)
             w4(:,:,1) = zero
             do i = 1, 2
                w3(:,:,i+1) = hvcoord%etam(i)
                w4(:,:,i+1) = eta_dot(:,:,i,t0)
             end do
          else
             do i = 1, 2
                w3(:,:,i) = hvcoord%etam(nlev-2+i)
                w4(:,:,i) = eta_dot(:,:,nlev-2+i,t0)
             end do
             w3(:,:,3) = hvcoord%etai(nlevp)
             w4(:,:,3) = zero
          end if
          call eval_lagrange_poly_derivative(3, w3, w4, w2, w1)
       else
          k1 = k-1
          k2 = k+1
          do i = 1, 3
             w3(:,:,i) = hvcoord%etam(k1-1+i)
          end do
          call eval_lagrange_poly_derivative(k2-k1+1, w3, eta_dot(:,:,k1:k2,t0), w2, w1)
       end if
       w3(:,:,1:2) = gradient_sphere(eta_dot(:,:,k,t0), deriv, elem%Dinv)
       vnode(4,:,:,k) = &
            half*(eta_dot(:,:,k,t0) + eta_dot(:,:,k,t1) &
            &     - dtsub*(vsph(:,:,1,k,t1)*w3(:,:,1) + vsph(:,:,2,k,t1)*w3(:,:,2) &
            &              + eta_dot(:,:,k,t1)*w1))
    end do
  end subroutine calc_eta_dot_formula_node_ref_mid

  subroutine update_dep_points_all(independent_time_steps, dtsub, nets, nete, vdep)
    ! Determine the departure points corresponding to the reference grid's
    ! arrival midpoints. Reads and writes dep_points_all. Reads vdep.
    
    use physical_constants, only: scale_factor

    logical, intent(in) :: independent_time_steps
    real(real_kind), intent(in) :: dtsub
    integer, intent(in) :: nets, nete
    real(real_kind), intent(in) :: vdep(:,:,:,:,:)

    real(real_kind) :: norm, p(3)
    integer :: ie, k, j, i

    do ie = nets, nete
       do k = 1, nlev
          do j = 1, np
             do i = 1, np
                ! Update horizontal position.
                p = dep_points_all(1:3,i,j,k,ie)
                p = p - dtsub*vdep(1:3,i,j,k,ie)/scale_factor
                if (is_sphere) then
                   norm = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
                   p = p/norm
                end if
                dep_points_all(1:3,i,j,k,ie) = p
                if (independent_time_steps) then
                   ! Update vertical position.
                   dep_points_all(4,i,j,k,ie) = dep_points_all(4,i,j,k,ie) - &
                        &                       dtsub*vdep(4,i,j,k,ie)
                end if
             end do
          end do
       end do
    end do
  end subroutine update_dep_points_all

  subroutine interp_departure_points_to_floating_level_midpoints( &
       elem, nets, nete, tl, hvcoord, dep_points_all, limcnt)
    ! Determine the departure points corresponding to the vertically Lagragnian
    ! grid's arrival midpoints, where the floating levels are those that evolve
    ! over the course of the full tracer time step. Also compute
    ! elem%derived%divdp, which holds the floating levels' dp values for later
    ! use in vertical remap.

    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
    type (hvcoord_t), intent(in) :: hvcoord
    type (TimeLevel_t), intent(in) :: tl
    real(real_kind), intent(inout) :: dep_points_all(:,:,:,:,:)
    integer, intent(inout) :: limcnt

    real(real_kind) :: deta_ref(nlevp), w1(np,np), v1(np,np,nlev), &
         &             v2(np,np,nlevp), p(3)
    integer :: ie, i, j, k, d

    call set_deta_tol(hvcoord)

    deta_ref(1) = hvcoord%etam(1) - hvcoord%etai(1)
    do k = 2, nlev
       deta_ref(k) = hvcoord%etam(k) - hvcoord%etam(k-1)
    end do
    deta_ref(nlevp) = hvcoord%etai(nlevp) - hvcoord%etam(nlev)

    do ie = nets, nete
       ! Surface pressure.
       w1 = hvcoord%hyai(1)*hvcoord%ps0 + sum(elem(ie)%state%dp3d(:,:,:,tl%np1), 3)

       ! Reconstruct Lagrangian levels at t1 on arrival column:
       !     eta_arr_int = I[eta_ref_mid([0,eta_dep_mid,1])](eta_ref_int)
       call limit_etam(hvcoord, deta_ref, dep_points_all(4,:,:,:,ie), v1, limcnt)
       v2(:,:,1) = hvcoord%etai(1)
       v2(:,:,nlevp) = hvcoord%etai(nlevp)
       call eta_interp_eta(hvcoord, v1, hvcoord%etam, &
            &              nlevp-2, hvcoord%etai(2:nlev), v2(:,:,2:nlev))
       call eta_to_dp(hvcoord, w1, v2, elem(ie)%derived%divdp)

       ! Compute Lagrangian level midpoints at t1 on arrival column:
       !     eta_arr_mid = I[eta_ref_mid([0,eta_dep_mid,1])](eta_ref_mid)
       call eta_interp_eta(hvcoord, v1, hvcoord%etam, &
            &              nlev, hvcoord%etam, v2(:,:,1:nlev))
       dep_points_all(4,:,:,:,ie) = v2(:,:,1:nlev)

       ! Compute departure horizontal points corresponding to arrival
       ! Lagrangian level midpoints:
       !     p_dep_mid(eta_arr_mid) = I[p_dep_mid(eta_ref_mid)](eta_arr_mid)
       do d = 1, 3
          v1 = dep_points_all(d,:,:,:,ie)
          call eta_interp_horiz(hvcoord, hvcoord%etam, v1, &
               &                v2(:,:,1:nlev), dep_points_all(d,:,:,:,ie))
       end do
       if (is_sphere) then
          ! Normalize p = (x,y,z).
          do k = 1, nlev
             do j = 1, np
                do i = 1, np
                   p = dep_points_all(1:3,i,j,k,ie)
                   p = p/sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
                   dep_points_all(1:3,i,j,k,ie) = p
                end do
             end do
          end do
       end if
    end do
  end subroutine interp_departure_points_to_floating_level_midpoints

  subroutine set_deta_tol(hvcoord)
    type (hvcoord_t), intent(in) :: hvcoord

    real(real_kind) :: deta_ave
    integer :: k

    if (deta_tol >= 0) return

    ! Benign write race condition. A thread might see eta_tol < 0 and set it
    ! here even as another thread does the same. But because there is no read
    ! and only one value to write, the redundant writes don't matter.

    deta_ave = (hvcoord%etai(nlev+1) - hvcoord%etai(1)) / nlev
    deta_tol = 10_real_kind*eps*deta_ave
  end subroutine set_deta_tol

  subroutine limit_etam(hvcoord, deta_ref, eta, eta_lim, cnt)
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: deta_ref(nlevp), eta(np,np,nlev)
    real(real_kind), intent(out) :: eta_lim(np,np,nlev)
    integer, intent(inout) :: cnt

    real(real_kind) :: deta(nlevp)
    integer :: i, j, k
    logical :: ok

    do j = 1, np
       do i = 1, np
          ! Check nonmonotonicity in eta.
          ok = eta(i,j,1) - hvcoord%etai(1) >= deta_tol
          if (ok) then
             do k = 2, nlev
                if (eta(i,j,k) - eta(i,j,k-1) < deta_tol) then
                   ok = .false.
                   exit
                end if
             end do
             if (ok) then
                ok = hvcoord%etai(nlevp) - eta(i,j,nlev) >= deta_tol
             end if
          end if
          ! eta is monotonically increasing, so don't need to do anything
          ! further.
          if (ok) then
             eta_lim(i,j,:) = eta(i,j,:)
             cycle
          end if
          
          deta(1) = eta(i,j,1) - hvcoord%etai(1)
          do k = 2, nlev
             deta(k) = eta(i,j,k) - eta(i,j,k-1)
          end do
          deta(nlevp) = hvcoord%etai(nlevp) - eta(i,j,nlev)
          ! [0, etam(1)] and [etam(nlev),1] are half levels, but deta_tol is so
          ! small there's no reason not to use it as a lower bound for these.
          cnt = cnt + 1
          call deta_caas(nlevp, deta_ref, deta_tol, deta)
          eta_lim(i,j,1) = hvcoord%etai(1) + deta(1)
          do k = 2, nlev
             eta_lim(i,j,k) = eta_lim(i,j,k-1) + deta(k)
          end do
       end do
    end do
  end subroutine limit_etam

  subroutine deta_caas(nlp, deta_ref, lo, deta)
    integer, intent(in) :: nlp
    real(real_kind), intent(in) :: deta_ref(nlp), lo
    real(real_kind), intent(inout) :: deta(nlp)

    real(real_kind) :: nerr, w(nlp)
    integer :: k

    nerr = zero
    do k = 1, nlp
       if (deta(k) < lo) then
          nerr = nerr + (deta(k) - lo)
          deta(k) = lo
          w(k) = zero
       else
          if (deta(k) > deta_ref(k)) then
             ! Only pull mass from intervals that are larger than their
             ! reference value. This concentrates changes to intervals that, by
             ! having a lot more mass than usual, drive other levels negative.
             w(k) = deta(k) - deta_ref(k)
          else
             w(k) = zero
          end if
       end if
    end do
    if (nerr /= zero) deta = deta + nerr*(w/sum(w))
  end subroutine deta_caas

  subroutine linterp(n, x, y, ni, xi, yi, caller)
    ! Linear interpolant: yi = I[y(x)](xi).
    ! x and xi must be in ascending order.
    ! xi(1) must be >= x(1) and xi(ni) must be <= x(n).

    use kinds, only: iulog

    integer, intent(in) :: n, ni
    real(real_kind), intent(in) :: x(n), y(n), xi(ni)
    real(real_kind), intent(out) :: yi(ni)
    character(len=*), intent(in) :: caller

    integer :: k, ki
    real(real_kind) :: a

    if (xi(1) < x(1) .or. xi(ni) > x(n)) then
       write(iulog,*) 'x', x
       write(iulog,*) 'xi', xi
       call abortmp('sl_vertically_remap_tracers> linterp xi out of bounds: ' &
            // trim(caller))
    end if

    k = 2
    do ki = 1, ni
       do while (x(k) < xi(ki))
          k = k + 1
       end do
       a = (xi(ki) - x(k-1))/(x(k) - x(k-1))
       yi(ki) = (1 - a)*y(k-1) + a*y(k)
    end do
  end subroutine linterp

  subroutine eta_interp_eta(hvcoord, x, y, ni, xi, yi)
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: x(np,np,nlev), y(nlev)
    integer, intent(in) :: ni
    real(real_kind), intent(in) :: xi(ni)
    real(real_kind), intent(out) :: yi(np,np,ni)

    real(real_kind) :: x01(nlev+2), y01(nlev+2)
    integer :: i, j

    x01(1) = hvcoord%etai(1)
    x01(nlev+2) = hvcoord%etai(nlevp)
    y01(1) = hvcoord%etai(1)
    y01(2:nlev+1) = y
    y01(nlev+2) = hvcoord%etai(nlevp)
    do j = 1, np
       do i = 1, np
          x01(2:nlev+1) = x(i,j,:)
          call linterp(nlev+2, x01, y01, &
               &       ni, xi, yi(i,j,:), &
               &       'eta_interp_eta')
       end do
    end do
  end subroutine eta_interp_eta

  subroutine eta_interp_horiz(hvcoord, x, y, xi, yi)
    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: x(nlev), y(np,np,nlev), xi(np,np,nlev)
    real(real_kind), intent(out) :: yi(np,np,nlev)

    real(real_kind) :: xbdy(nlev+2), ybdy(nlev+2)
    integer :: i, j

    xbdy(1) = hvcoord%etai(1)
    xbdy(2:nlev+1) = x
    xbdy(nlev+2) = hvcoord%etai(nlevp)
    do j = 1, np
       do i = 1, np
          ! Do constant interp outside of the etam support.
          ybdy(1) = y(i,j,1)
          ybdy(2:nlev+1) = y(i,j,:)
          ybdy(nlev+2) = y(i,j,nlev)
          call linterp(nlev+2, xbdy, ybdy, &
               &       nlev, xi(i,j,:), yi(i,j,:), &
               &       'eta_interp_horiz')
       end do
    end do
  end subroutine eta_interp_horiz

  subroutine eta_to_dp(hvcoord, ps, etai, dp)
    !    e = A(e) + B(e)
    ! p(e) = A(e) p0 + B(e) ps
    !      = e p0 + B(e) (ps - p0)
    !     a= e p0 + I[Bi(eref)](e) (ps - p0)

    type (hvcoord_t), intent(in) :: hvcoord
    real(real_kind), intent(in) :: ps(np,np), etai(np,np,nlevp)
    real(real_kind), intent(out) :: dp(np,np,nlev)

    real(real_kind) :: Bi(nlevp), dps
    integer :: i, j, k

    do j = 1, np
       do i = 1, np
          call linterp(nlevp, hvcoord%etai, hvcoord%hybi, &
               &       nlevp, etai(i,j,:), Bi, &
               &       'eta_to_dp')
          dps = ps(i,j) - hvcoord%ps0
          do k = 1, nlev
             dp(i,j,k) = (etai(i,j,k+1) - etai(i,j,k))*hvcoord%ps0 + &
                  &      (Bi(k+1) - Bi(k))*dps
          end do
       end do
    end do
  end subroutine eta_to_dp

  subroutine dss_vnode(elem, nets, nete, hybrid, vnode)
    type (element_t), intent(in) :: elem(:)
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nets, nete
    real(real_kind) :: vnode(:,:,:,:,:)

    integer :: nd, nlyr, ie, k, d

    nd = size(vnode, 1)
    nlyr = nd*nlev
    
    do ie = nets, nete
       do k = 1, nlev
          do d = 1, nd
             vnode(d,:,:,k,ie) = vnode(d,:,:,k,ie)* &
                  &              elem(ie)%spheremp*elem(ie)%rspheremp
          end do
       end do
       do d = 1, nd
          call edgeVpack_nlyr(edge_g, elem(ie)%desc, vnode(d,:,:,:,ie), &
               &              nlev, nlev*(d-1), nlyr)
       end do
    end do

    call t_startf('SLMM_bexchV')
    call bndry_exchangeV(hybrid, edge_g)
    call t_stopf('SLMM_bexchV')

    do ie = nets, nete
       do d = 1, nd
          call edgeVunpack_nlyr(edge_g, elem(ie)%desc, vnode(d,:,:,:,ie), &
               &                nlev, nlev*(d-1), nlyr)
       end do
    end do

#if (defined HORIZ_OPENMP)
    !$omp barrier
#endif
  end subroutine dss_vnode

  subroutine dss_divdp(elem, nets, nete, hybrid)
    type (element_t), intent(inout) :: elem(:)
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in) :: nets, nete

    integer :: ie, k

    do ie = nets, nete
       do k = 1, nlev
          elem(ie)%derived%divdp(:,:,k) = elem(ie)%derived%divdp(:,:,k)* &
               &                          elem(ie)%spheremp*elem(ie)%rspheremp
       end do
       call edgeVpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%divdp, &
            &              nlev, 0, nlev)
    end do

    call t_startf('SLMM_bexchV')
    call bndry_exchangeV(hybrid, edge_g)
    call t_stopf('SLMM_bexchV')

    do ie = nets, nete
       call edgeVunpack_nlyr(edge_g, elem(ie)%desc, elem(ie)%derived%divdp, &
            &                nlev, 0, nlev)
    end do

#if (defined HORIZ_OPENMP)
    !$omp barrier
#endif
  end subroutine dss_divdp

  function assert(b, msg) result(nerr)
    use kinds, only: iulog

    logical, intent(in) :: b
    character(*), optional, intent(in) :: msg

    character(len=128) :: s
    integer :: nerr

    nerr = 0
    if (b) return

    s = ''
    if (present(msg)) s = msg
    write(iulog,'(a,a)') 'COMPOSE> sl_advection ASSERT: ', trim(s)
    nerr = 1
  end function assert

  function test_linterp() result (nerr)
    integer, parameter :: n = 128, ni = 111

    real(real_kind) :: x(n), y(n), xi(ni), yi(ni), yin(n), a
    integer :: k, nerr

    call random_number(x)
    do k = 2, n
       x(k) = x(k) + x(k-1)
    end do
    y = 3*x

    do k = 1, ni
       a = real(k, real_kind)/(ni+1)
       xi(k) = (1 - a)*x(1) + a*x(n)
    end do

    call linterp(n, x, y, ni, xi, yi, 'test_linterp 1')
    nerr = assert(maxval(abs( yi - 3*xi)) < 100*eps*x(n), 'linterp 1')
    
    call linterp(n, x, y, n, x, yin, 'test_linterp 2')
    nerr = nerr + assert(maxval(abs(yin - y)) < 10*eps, 'linterp 2')
  end function test_linterp

  function test_eta_to_dp(hvcoord) result(nerr)
    type (hvcoord_t), intent(in) :: hvcoord

    real(real_kind) :: ps(np,np), etai(np,np,nlevp), dp1(np,np,nlev), &
         &             dp2(np,np,nlev)
    integer :: nerr, i, j, k

    nerr = 0

    call random_number(ps)
    ps = (one + 0.2*(ps - 0.5))*hvcoord%ps0

    do k = 1, nlev
       dp1(:,:,k) = (hvcoord%hyai(k+1) - hvcoord%hyai(k))*hvcoord%ps0 + &
            &       (hvcoord%hybi(k+1) - hvcoord%hybi(k))*ps
    end do

    ! Test that for etai_ref we get the same as the usual formula.
    do j = 1, np
       do i = 1, np
          etai(i,j,:) = hvcoord%etai
       end do
    end do
    call eta_to_dp(hvcoord, ps, etai, dp2)
    nerr = nerr + assert(maxval(abs(dp2-dp1)) < 100*eps*maxval(dp1), 'eta_to_dp 1')
  end function test_eta_to_dp

  function test_deta_caas() result(nerr)
    integer, parameter :: nl = 128, nlp = nl+1
    
    real(real_kind) :: deta_ref(nlp), etam_ref(nl), deta_tol, etam(nl), &
         &             deta(nlp)
    integer :: i, k, nerr

    nerr = 0

    call random_number(deta_ref)
    deta_ref = deta_ref + 0.1
    deta_ref = deta_ref/sum(deta_ref)

    deta_tol = 10_real_kind*eps*sum(deta_ref)/size(deta_ref)
    nerr = nerr + assert(deta_tol < minval(deta_ref), 'deta_caas 1')

    ! Test: Input not touched.
    deta = deta_ref
    call deta_caas(nlp, deta_ref, deta_tol, deta)
    nerr = nerr + assert(maxval(abs(deta-deta_ref)) == zero, 'deta_caas 2')
    
    etam_ref(1) = deta_ref(1)
    do k = 2, nl
       etam_ref(k) = etam_ref(k-1) + deta_ref(k)
    end do

    ! Test: Modify one etam and only adjacent intervals change beyond eps.
    do i = 1, 2
       etam = etam_ref
       if (i == 1) then
          etam(11) = etam(11) + 1.1
       else
          etam(11) = etam(11) - 13.1
       end if
       deta(1) = etam(1)
       do k = 2, nl
          deta(k) = etam(k) - etam(k-1)
       end do
       deta(nlp) = one - etam(nl)
       nerr = nerr + assert(minval(deta) < deta_tol, 'deta_caas 3')
       call deta_caas(nlp, deta_ref, deta_tol, deta)
       nerr = nerr + assert(minval(deta) == deta_tol, 'deta_caas 4')
       nerr = nerr + assert(abs(sum(deta) - one) < 100*eps, 'deta_caas 5')
       deta = abs(deta - deta_ref)
       nerr = nerr + assert(maxval(deta(:10)) < 100*eps, 'deta_caas 6')
       nerr = nerr + assert(maxval(deta(13:)) < 100*eps, 'deta_caas 7')
    end do

    ! Test: Completely messed up levels.
    call random_number(deta)
    deta = deta - 0.5_real_kind
    if (sum(deta) < 0.1) deta = deta + (0.1 - sum(deta))/nlp
    deta = deta/sum(deta)
    call deta_caas(nlp, deta_ref, deta_tol, deta)
    nerr = nerr + assert(minval(deta) == deta_tol, 'deta_caas 8')
    nerr = nerr + assert(abs(sum(deta) - one) < 1000*eps, 'deta_caas 9')
  end function test_deta_caas

  function test_init_velocity_record() result(error)
    integer :: dtf, drf, nsub, nvel, e, error
    type (velocity_record_t) :: v

    error = 0
    dtf = 6
    drf = 2
    nsub = 3
    nvel = 4
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    nvel = 3
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    drf = 3
    nvel = 6
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    drf = 1
    nsub = 5
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    dtf = 12
    drf = 2
    nsub = 3
    nvel = -1
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    nsub = 5
    nvel = 5
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)
    dtf = 27
    drf = 3
    nsub = 51
    nvel = 99
    call init_velocity_record(1, dtf, drf, nsub, nvel, v, e)
    call test(dtf, drf, nsub, nvel, v, e)
    if (e > 0) error = 1
    call cleanup(v)

  contains
    subroutine test(dtf, drf, nsub, nvel, v, error)
      integer, intent(in) :: dtf, drf, nsub, nvel
      integer, intent(inout) :: error
      type (velocity_record_t), intent(in) :: v

      real(kind=real_kind) :: endslots(2), ys(dtf), a, x, y, y0, y1, &
           &                  xsup(2), ysup(2)
      integer :: n, e, i, k

      e = 0

      if (modulo(dtf, drf) /= 0) then
         print *, 'TESTING ERROR: dtf % drf == 0 is required:', dtf, drf
      end if

      ! Check that t_vel is monotonically increasing.
      do n = 2, v%nvel
         if (v%t_vel(n) <= v%t_vel(n-1)) e = 1
      end do

      ! Check that obs_slots does not reference end points. This should not
      ! happen b/c nvel <= navail and observations are uniformly spaced.
      do n = 1, dtf
         do i = 1, 2
            if (v%obs_slots(n,i) == 0 .or. v%obs_slots(n,i) == dtf) e = 11
         end do
      end do

      ! Check that weights sum to 1.
      ys = 0
      do n = 1, dtf
         do i = 1, 2
            if (v%obs_slots(n,i) > 0) then
               ys(v%obs_slots(n,i)) = ys(v%obs_slots(n,i)) + v%obs_wts(n,i)
            end if
         end do
      end do
      do i = 2, v%nvel-1
         if (abs(ys(i) - 1) > 1e3*eps) e = 12
      end do

      ! Test for exact interp of an affine function.
      ! Observe data forward in time.
      endslots(1) = tfn(0.d0)
      endslots(2) = tfn(real(dtf, real_kind))
      ys(:) = 0
      do n = 1, dtf
         if (modulo(n, drf) /= 0) cycle
         y = tfn(real(n, real_kind))
         do i = 1, 2
            if (v%obs_slots(n,i) == -1) cycle
            ys(v%obs_slots(n,i)) = ys(v%obs_slots(n,i)) + v%obs_wts(n,i)*y
         end do
      end do
      ! Use the data backward in time.
      do n = 1, nsub
         ! Each segment orders the data forward in time. Thus, data are always
         ! ordered forward in time but used backward.
         do i = 1, 2
            k = nsub - n + (i-1)
            xsup(i) = (k*v%t_vel(v%nvel))/nsub
            k = v%run_step(n+1-i)
            if (k == 2) then
               y0 = endslots(1)
            else
               y0 = ys(k-1)
            end if
            if (k == v%nvel) then
               y1 = endslots(2)
            else
               y1 = ys(k)
            end if
            ysup(i) = ((v%t_vel(k) - xsup(i))*y0 + (xsup(i) - v%t_vel(k-1))*y1) / &
                 &    (v%t_vel(k) - v%t_vel(k-1))
         end do
         do i = 0, 10
            a = real(i, real_kind)/10
            x = (1-a)*xsup(1) + a*xsup(2)
            y = (1-a)*ysup(1) + a*ysup(2)
            if (abs(y - tfn(x)) > 1e3*eps) e = 2
         end do
      end do

      if (error > 0 .or. e > 0) then
         print *, 'ERROR', error, e
         print '(a,i3,a,i3,a,i3,a,i3,a,i3)', 'dtf',dtf,' drf',drf,' nsub',nsub, &
              ' nvel',nvel,' v%nvel',v%nvel
         print '(a,99es11.3)', '  t_vel', (v%t_vel(n),n=1,v%nvel)
         do n = 1, dtf-1
            print '(3i3,2f5.2)', n, v%obs_slots(n,:), v%obs_wts(n,:)
         end do
         do n = 0, nsub
            print '(i3,i3)', n, v%run_step(n)
         end do
         error = 1
      end if
    end subroutine test

    function tfn(x) result(y)
      real(kind=real_kind), intent(in) :: x
      real(kind=real_kind) :: y

      y = 7.1*x - 11.5
    end function tfn

    subroutine cleanup(v)
      type (velocity_record_t), intent(inout) :: v
      deallocate(v%t_vel, v%obs_slots, v%obs_wts, v%run_step, v%vel, v%dp)
    end subroutine cleanup
  end function test_init_velocity_record

  subroutine sl_unittest(par, hvcoord)
    use kinds, only: iulog

    type (parallel_t), intent(in) :: par
    type (hvcoord_t), intent(in) :: hvcoord

    integer :: n(6)

    n(1) = test_lagrange()
    n(2) = test_reconstruct_and_limit_dp()
    n(3) = test_deta_caas()
    n(4) = test_linterp()
    n(5) = test_eta_to_dp(hvcoord)
    n(6) = test_init_velocity_record()

    if (sum(n) > 0 .and. par%masterproc) then
       write(iulog,'(a,6i2)') 'COMPOSE> sl_unittest FAIL ', n
    end if
  end subroutine sl_unittest

end module sl_advection
