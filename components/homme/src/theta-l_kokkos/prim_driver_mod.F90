#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module prim_driver_mod

  use prim_driver_base,     only: deriv1, smooth_topo_datasets
  use prim_cxx_driver_base, only: prim_init1, prim_finalize
  use kinds,                only : real_kind
  use dimensions_mod,       only : qsize, nelemd, np, qsize
  use element_mod,          only : element_t
  use prim_driver_base,     only : deriv1, smooth_topo_datasets
  use prim_cxx_driver_base, only : prim_init1, prim_finalize
  use physical_constants,   only : scale_factor, laplacian_rigid_factor
  use hybrid_mod,           only : hybrid_t
  use hybvcoord_mod,        only : hvcoord_t
  use derivative_mod,       only : derivative_t
  use time_mod,             only : timelevel_t
  
  implicit none

  public :: prim_init2
  public :: prim_run_subcycle
  public :: prim_init_elements_views
  public :: prim_init_grid_views
  public :: prim_init_geopotential_views
  public :: prim_init_state_views
  public :: prim_init_ref_states_views
  public :: prim_init_diags_views

  type, private :: PrescribedWind_t
     type (element_t), pointer :: elem(:)
     type (hybrid_t) :: hybrid
     type (hvcoord_t) :: hvcoord
     type (derivative_t) :: deriv
     integer :: nets, nete
  end type PrescribedWind_t

  type (PrescribedWind_t), private :: prescribed_wind_args

contains

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use prim_driver_base, only : deriv1, prim_init2_base => prim_init2
    use prim_state_mod,   only : prim_printstate
    use theta_f2c_mod,    only : initialize_dp3d_from_ps_c
    use control_mod,      only : prescribed_wind
    !
    ! Inputs
    !
    type (element_t),   intent(inout)         :: elem(:)
    type (hybrid_t),    intent(in)            :: hybrid
    type (TimeLevel_t), intent(inout)         :: tl       ! time level struct
    type (hvcoord_t),   intent(inout), target :: hvcoord  ! hybrid vertical coordinate struct
    integer,            intent(in)            :: nets     ! starting thread element number (private)
    integer,            intent(in)            :: nete     ! ending thread element number   (private)

    ! Call the base version of prim_init2
    call prim_init2_base(elem,hybrid,nets,nete,tl,hvcoord)

    ! Init the c data structures
    call prim_create_c_data_structures(tl,hvcoord,elem(1)%mp)

    !Init the kokkos functors (and their boundary exchanges)
    call prim_init_kokkos_functors ()

    ! Init the kokkos views
    call prim_init_elements_views (elem)

    ! Initialize dp3d from ps_v
    call initialize_dp3d_from_ps_c ()

    if (prescribed_wind == 1) then
       call init_standalone_test(elem,deriv1,hybrid,hvcoord,tl,nets,nete)
    end if
  end subroutine prim_init2

  subroutine prim_create_c_data_structures (tl, hvcoord, mp)
    use iso_c_binding, only : c_loc, c_ptr, C_NULL_CHAR
    use theta_f2c_mod, only : init_reference_element_c, init_simulation_params_c, &
                              init_time_level_c, init_hvcoord_c, init_elements_c
    use time_mod,      only : TimeLevel_t, nsplit
    use hybvcoord_mod, only : hvcoord_t
    use control_mod,   only : limiter_option, rsplit, qsplit, tstep_type, statefreq,   &
                              nu, nu_p, nu_q, nu_s, nu_div, nu_top, vert_remap_q_alg,  &
                              hypervis_order, hypervis_subcycle, hypervis_subcycle_tom,&
                              hypervis_scaling,                                        &
                              ftype, prescribed_wind, use_moisture, disable_diagnostics,   &
                              use_cpstar, transport_alg, theta_hydrostatic_mode,       &
                              dcmip16_mu, theta_advect_form, test_case,                &
                              MAX_STRING_LEN, dt_remap_factor, dt_tracer_factor,       &
                              pgrad_correction, dp3d_thresh, vtheta_thresh,            &
                              internal_diagnostics_level
    !
    ! Input(s)
    !
    type (TimeLevel_t),       intent(in) :: tl
    type (hvcoord_t), target, intent(in) :: hvcoord
    real (kind=real_kind),    intent(in) :: mp(np,np)
    !
    ! Local(s)
    !
    integer :: ie
    real (kind=real_kind), target :: dvv (np,np), elem_mp(np,np)
    type (c_ptr) :: hybrid_am_ptr, hybrid_ai_ptr, hybrid_bm_ptr, hybrid_bi_ptr
    character(len=MAX_STRING_LEN), target :: test_name

    integer :: disable_diagnostics_int, theta_hydrostatic_mode_int, use_moisture_int

    ! Initialize the C++ reference element structure (i.e., pseudo-spectral deriv matrix and ref element mass matrix)
    dvv = deriv1%dvv
    elem_mp = mp
    call init_reference_element_c(c_loc(dvv),c_loc(elem_mp))

    ! Fill the simulation params structures in C++
    test_name = TRIM(test_case) // C_NULL_CHAR

    disable_diagnostics_int = 0
    if (disable_diagnostics) disable_diagnostics_int = 1
    use_moisture_int = 0
    if (use_moisture) use_moisture_int = 1
    theta_hydrostatic_mode_int = 0
    if (theta_hydrostatic_mode) theta_hydrostatic_mode_int = 1

    call init_simulation_params_c (vert_remap_q_alg, limiter_option, rsplit, qsplit, tstep_type,  &
                                   qsize, statefreq, nu, nu_p, nu_q, nu_s, nu_div, nu_top,        &
                                   hypervis_order, hypervis_subcycle, hypervis_subcycle_tom,      &
                                   hypervis_scaling,                                              &
                                   dcmip16_mu, ftype, theta_advect_form,                          &
                                   prescribed_wind,                                               &
                                   use_moisture_int,                                              &
                                   disable_diagnostics_int,                                       &
                                   use_cpstar,                                                    &
                                   transport_alg,                                                 &
                                   theta_hydrostatic_mode_int,                                    &
                                   c_loc(test_name),                                              &
                                   dt_remap_factor, dt_tracer_factor,                             &
                                   scale_factor, laplacian_rigid_factor,                          &
                                   nsplit,                                                        &
                                   pgrad_correction,                                              &
                                   dp3d_thresh, vtheta_thresh, internal_diagnostics_level)

    ! Initialize time level structure in C++
    call init_time_level_c(tl%nm1, tl%n0, tl%np1, tl%nstep, tl%nstep0)

    ! Initialize the hybrid vertical coordinate in C++
    hybrid_am_ptr = c_loc(hvcoord%hyam)
    hybrid_ai_ptr = c_loc(hvcoord%hyai)
    hybrid_bm_ptr = c_loc(hvcoord%hybm)
    hybrid_bi_ptr = c_loc(hvcoord%hybi)
    call init_hvcoord_c (hvcoord%ps0,hybrid_am_ptr,hybrid_ai_ptr,hybrid_bm_ptr,hybrid_bi_ptr)

    ! Initialize the C++ elements structure
    call init_elements_c (nelemd)

  end subroutine prim_create_c_data_structures

  subroutine prim_init_grid_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use theta_f2c_mod, only : init_elements_2d_c
    use control_mod,   only : geometry
    use coordinate_systems_mod, only : change_coordinates, cartesian3D_t
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)
    !
    ! Local(s)
    !
    real (kind=real_kind), target, dimension(np,np,2,2)     :: elem_D, elem_Dinv, elem_metinv, elem_tensorvisc
    real (kind=real_kind), target, dimension(np,np)         :: elem_mp, elem_fcor, elem_spheremp
    real (kind=real_kind), target, dimension(np,np)         :: elem_rspheremp, elem_metdet
    real (kind=real_kind), target, dimension(np,np,3,2)     :: elem_vec_sph2cart

    type (c_ptr) :: elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr
    type (c_ptr) :: elem_spheremp_ptr, elem_rspheremp_ptr
    type (c_ptr) :: elem_metdet_ptr, elem_metinv_ptr
    type (c_ptr) :: elem_tensorvisc_ptr, elem_vec_sph2cart_ptr

    type (cartesian3D_t) :: sphere_cart
    real (kind=real_kind) :: sphere_cart_vec(3,np,np), sphere_latlon_vec(2,np,np)

    integer :: ie, i, j
    logical :: is_sphere

    elem_D_ptr            = c_loc(elem_D)
    elem_Dinv_ptr         = c_loc(elem_Dinv)
    elem_fcor_ptr         = c_loc(elem_fcor)
    elem_spheremp_ptr     = c_loc(elem_spheremp)
    elem_rspheremp_ptr    = c_loc(elem_rspheremp)
    elem_metdet_ptr       = c_loc(elem_metdet)
    elem_metinv_ptr       = c_loc(elem_metinv)
    elem_tensorvisc_ptr   = c_loc(elem_tensorvisc)
    elem_vec_sph2cart_ptr = c_loc(elem_vec_sph2cart)

    is_sphere = trim(geometry) /= 'plane'

    do ie=1,nelemd
      elem_D            = elem(ie)%D
      elem_Dinv         = elem(ie)%Dinv
      elem_fcor         = elem(ie)%fcor
      elem_spheremp     = elem(ie)%spheremp
      elem_rspheremp    = elem(ie)%rspheremp
      elem_metdet       = elem(ie)%metdet
      elem_metinv       = elem(ie)%metinv
      elem_tensorvisc   = elem(ie)%tensorVisc
      elem_vec_sph2cart = elem(ie)%vec_sphere2cart
      do j = 1,np
         do i = 1,np
            if (is_sphere) then
               sphere_cart = change_coordinates(elem(ie)%spherep(i,j))
               sphere_cart_vec(1,i,j) = sphere_cart%x
               sphere_cart_vec(2,i,j) = sphere_cart%y
               sphere_cart_vec(3,i,j) = sphere_cart%z
            else
               sphere_cart_vec(1,i,j) = elem(ie)%spherep(i,j)%lon
               sphere_cart_vec(2,i,j) = elem(ie)%spherep(i,j)%lat
               sphere_cart_vec(3,i,j) = 0
            end if
            sphere_latlon_vec(1,i,j) = elem(ie)%spherep(i,j)%lat
            sphere_latlon_vec(2,i,j) = elem(ie)%spherep(i,j)%lon
         end do
      end do
      call init_elements_2d_c (ie-1,                                      &
                               elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr,  &
                               elem_spheremp_ptr, elem_rspheremp_ptr,     &
                               elem_metdet_ptr, elem_metinv_ptr,          &
                               elem_tensorvisc_ptr, elem_vec_sph2cart_ptr,&
                               sphere_cart_vec, sphere_latlon_vec)
    enddo
  end subroutine prim_init_grid_views

  subroutine prim_init_geopotential_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use theta_f2c_mod, only : init_geopotential_c
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)
    !
    ! Local(s)
    !
    real (kind=real_kind), target, dimension(np,np)         :: elem_state_phis
    real (kind=real_kind), target, dimension(np,np,2)       :: elem_gradphis

    type (c_ptr) :: elem_state_phis_ptr, elem_gradphis_ptr

    integer :: ie

    elem_state_phis_ptr   = c_loc(elem_state_phis)
    elem_gradphis_ptr     = c_loc(elem_gradphis)

    do ie=1,nelemd
      elem_state_phis   = elem(ie)%state%phis
      elem_gradphis     = elem(ie)%derived%gradphis
      call init_geopotential_c (ie-1, elem_state_phis_ptr, elem_gradphis_ptr)
    enddo
  end subroutine prim_init_geopotential_views

  subroutine prim_init_state_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use element_state, onlY : elem_state_dp3d, elem_state_phinh_i, elem_state_ps_v, &
                              elem_state_qdp, elem_state_v,                         &
                              elem_state_vtheta_dp, elem_state_w_i
    use theta_f2c_mod, only : init_elements_states_c
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)
    !
    ! Local(s)
    !

    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr
    type (c_ptr) :: elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_state_Qdp_ptr

    elem_state_v_ptr         = c_loc(elem_state_v)
    elem_state_w_i_ptr       = c_loc(elem_state_w_i)
    elem_state_vtheta_dp_ptr = c_loc(elem_state_vtheta_dp)
    elem_state_phinh_i_ptr   = c_loc(elem_state_phinh_i)
    elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
    elem_state_Qdp_ptr       = c_loc(elem_state_Qdp)
    elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
    call init_elements_states_c (elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr,   &
                                 elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr, &
                                 elem_state_Qdp_ptr)
  end subroutine prim_init_state_views

  subroutine prim_init_ref_states_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use element_state, onlY : elem_theta_ref, elem_dp_ref, elem_phi_ref
    use theta_f2c_mod, only : init_reference_states_c
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)
    !
    ! Local(s)
    !
    type (c_ptr) :: elem_theta_ref_ptr, elem_dp_ref_ptr, elem_phi_ref_ptr

    elem_theta_ref_ptr = c_loc(elem_theta_ref)
    elem_dp_ref_ptr    = c_loc(elem_dp_ref)
    elem_phi_ref_ptr   = c_loc(elem_phi_ref)
    call init_reference_states_c (elem_theta_ref_ptr, elem_dp_ref_ptr, elem_phi_ref_ptr)
  end subroutine prim_init_ref_states_views

  subroutine prim_init_diags_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use element_state, onlY : elem_accum_iener, elem_accum_kener, elem_accum_pener, &
                              elem_accum_q1mass, elem_accum_qmass, elem_accum_qvar, &
                              elem_state_q
    use theta_f2c_mod, only : init_diagnostics_c
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)
    !
    ! Local(s)
    !
    type (c_ptr) :: elem_accum_iener_ptr, elem_accum_kener_ptr, elem_accum_pener_ptr
    type (c_ptr) :: elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr
    type (c_ptr) :: elem_state_q_ptr

    elem_state_q_ptr         = c_loc(elem_state_q)
    elem_accum_qvar_ptr      = c_loc(elem_accum_qvar)
    elem_accum_qmass_ptr     = c_loc(elem_accum_qmass)
    elem_accum_q1mass_ptr    = c_loc(elem_accum_q1mass)
    elem_accum_iener_ptr     = c_loc(elem_accum_iener)
    elem_accum_kener_ptr     = c_loc(elem_accum_kener)
    elem_accum_pener_ptr     = c_loc(elem_accum_pener)
    call init_diagnostics_c (elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr, &
                             elem_accum_q1mass_ptr, elem_accum_iener_ptr,                 &
                             elem_accum_kener_ptr, elem_accum_pener_ptr)
  end subroutine prim_init_diags_views

  subroutine prim_init_elements_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem (:)

    ! Initialize the grid-related views in C++
    call prim_init_grid_views (elem)

    ! Initialize phis/gradphis views in C++
    call prim_init_geopotential_views (elem)

    ! Initialize the 3d states views in C++
    call prim_init_state_views (elem)

    ! Initialize the reference states in C++
    call prim_init_ref_states_views (elem)

    ! Initialize the diagnostics arrays in C++
    call prim_init_diags_views (elem)
  end subroutine prim_init_elements_views

  subroutine prim_init_kokkos_functors (allocate_buffer)
    !todo-repo-unification Remove the use of c_bool here. It's used in purely
    ! F90 code.
    use iso_c_binding, only : c_int, c_bool
    use theta_f2c_mod, only : init_functors_c, init_boundary_exchanges_c
    !
    ! Optional Input
    !
    logical(kind=c_bool), intent(in), optional :: allocate_buffer  ! Whether functor memory buffer should be allocated internally
    integer(kind=c_int) :: ab
    ! Initialize the C++ functors in the C++ context
    ! If no argument allocate_buffer is present,
    ! let Homme internally allocate buffers
    ab = 1
    if (present(allocate_buffer)) then
       if (.not. allocate_buffer) ab = 0
    end if
    call init_functors_c (ab)

    ! Initialize boundary exchange structure in C++
    call init_boundary_exchanges_c ()

  end subroutine prim_init_kokkos_functors

  subroutine prim_run_subcycle(elem, hybrid, nets, nete, dt, single_column, tl, hvcoord, nsplit_iteration)
    use iso_c_binding,  only : c_int, c_ptr, c_loc
    use control_mod,    only : qsplit, rsplit, statefreq, disable_diagnostics, &
                               dt_remap_factor, dt_tracer_factor
    use dimensions_mod, only : nelemd
    use element_state,  only : elem_state_v, elem_state_w_i, elem_state_vtheta_dp,     &
                               elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v,   &
                               elem_state_Qdp, elem_state_Q, elem_derived_omega_p,     &
                               elem_derived_FM, elem_derived_FVTheta, elem_derived_FT, &
                               elem_derived_FPHI, elem_derived_FQ
    use hybrid_mod,     only : hybrid_t
    use hybvcoord_mod,  only : hvcoord_t
    use kinds,          only : real_kind
    use time_mod,       only : timelevel_t, nextOutputStep, nsplit, TimeLevel_Qdp
    use control_mod,    only : statefreq, prescribed_wind
    use parallel_mod,   only : abortmp
    use perf_mod,       only : t_startf, t_stopf
    use prim_state_mod, only : prim_printstate
    use theta_f2c_mod,  only : prim_run_subcycle_c, cxx_push_results_to_f90
    use theta_f2c_mod,  only : push_forcing_to_c, sync_diagnostics_to_host_c
    !
    ! Inputs
    !
    type (element_t) ,    intent(inout) :: elem(:)
    type (hybrid_t),      intent(in)    :: hybrid                       ! distributed parallel structure (shared)
    type (hvcoord_t),     intent(in)    :: hvcoord                      ! hybrid vertical coordinate struct
    integer,              intent(in)    :: nets                         ! starting thread element number (private)
    integer,              intent(in)    :: nete                         ! ending thread element number   (private)
    real(kind=real_kind), intent(in)    :: dt                           ! "timestep dependent" timestep
    logical,              intent(in)    :: single_column
    type (TimeLevel_t),   intent(inout) :: tl
    integer,              intent(in)    :: nsplit_iteration             !  = 1 .. nsplit
    !
    ! Locals
    !
    logical :: compute_diagnostics
    integer (kind=c_int) :: nstep_end, nstep_c, nm1_c, n0_c, np1_c
    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr
    type (c_ptr) :: elem_state_dp3d_ptr, elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_omega_p_ptr
    integer :: n0_qdp, np1_qdp
    real(kind=real_kind) :: dt_remap, dt_q, eta_ave_w
    logical :: compute_forcing_and_push_to_c, push_to_f

    if (nsplit<1) then
      call abortmp ('nsplit_is less than 1.')
    endif
    if (nsplit_iteration < 1 .or. nsplit_iteration > nsplit) then
      call abortmp ('nsplit_iteration out of range')
    endif
    if (nets/=1 .or. nete/=nelemd) then
      call abortmp ('We don''t allow to call C routines from a horizontally threaded region')
    endif

    nstep_end = tl%nstep + qsplit
    if (rsplit>0) then
       nstep_end = tl%nstep + qsplit*rsplit  ! nstep at end of this routine
    endif
    compute_diagnostics   = .false.
    if (MODULO(nstep_end,statefreq)==0 .or. (tl%nstep <= tl%nstep0+(nstep_end-tl%nstep) )) then
      compute_diagnostics = .true.
    endif
    if (disable_diagnostics) then
      compute_diagnostics = .false.
    endif

    compute_forcing_and_push_to_c = is_push_to_c_required(nsplit_iteration)

    dt_q = dt*dt_tracer_factor
    if (dt_remap_factor == 0) then
       dt_remap = dt
    else
       dt_remap = dt*dt_remap_factor
    end if

    call TimeLevel_Qdp(tl, dt_tracer_factor, n0_qdp, np1_qdp)

    ! Test forcing is only for standalone Homme (and only for some tests/configurations)
    if (compute_forcing_and_push_to_c) then
      call compute_test_forcing_f(elem,hybrid,hvcoord,tl%n0,n0_qdp,max(dt_q,dt_remap),nets,nete,tl)
      call t_startf('push_to_cxx')
      call push_forcing_to_c(elem_derived_FM,   elem_derived_FVTheta, elem_derived_FT, &
                             elem_derived_FPHI, elem_derived_FQ)
      call t_stopf('push_to_cxx')
    end if
    if (prescribed_wind == 1) then
       call init_prescribed_wind_subcycle(elem,nets,nete,tl)
    end if

    call prim_run_subcycle_c(dt,nstep_c,nm1_c,n0_c,np1_c,nextOutputStep,nsplit_iteration)

    ! Set final timelevels from C into Fortran structure
    tl%nstep = nstep_c
    tl%nm1   = nm1_c + 1
    tl%n0    = n0_c  + 1
    tl%np1   = np1_c + 1

    push_to_f = is_push_to_f_required(tl,statefreq,nextOutputStep,compute_diagnostics,nsplit_iteration)

    if (push_to_f) then
      ! Set pointers to states
      elem_state_v_ptr         = c_loc(elem_state_v)
      elem_state_w_i_ptr       = c_loc(elem_state_w_i)
      elem_state_vtheta_dp_ptr = c_loc(elem_state_vtheta_dp)
      elem_state_phinh_i_ptr   = c_loc(elem_state_phinh_i)
      elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
      elem_state_Qdp_ptr       = c_loc(elem_state_Qdp)
      elem_state_Q_ptr         = c_loc(elem_state_Q)
      elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
      elem_derived_omega_p_ptr = c_loc(elem_derived_omega_p)

      ! Copy cxx arrays back to f90 structures
      call t_startf('push_to_f90')
      call cxx_push_results_to_f90(elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr,   &
                                   elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr, &
                                   elem_state_Qdp_ptr, elem_state_Q_ptr, elem_derived_omega_p_ptr)
      call t_stopf('push_to_f90')
    endif

    ! Print some diagnostic information

    if (compute_diagnostics) then
       call t_startf('sync_diag_to_host')
       call sync_diagnostics_to_host_c()
       call t_stopf('sync_diag_to_host')
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle


!the next 2 routines have logic for push to/from F and for forcing routine
!
!STANDALONE HOMME there are 3 cases:
!
!performance:
! (no forcing, no push to c) -> (subcycle) -> (no push to f)
!
!test without forcing (can be performant if output is only at the end):
! (no forcing, no push to c) -> (subcycle) -> (push to f only for output/diagnostics)
!
!test with forcing (not performant):
! (always forcing and push to c) -> (subcycle) -> (always push to f)

  function is_push_to_c_required(nsplit_iter) result(compute_forcing_and_push_to_c)

    use control_mod, only: test_with_forcing

    integer,              intent(in) :: nsplit_iter

    logical :: compute_forcing_and_push_to_c

    compute_forcing_and_push_to_c = .false.

    ! Scream already computes all forcing using the same pointers
    ! stored in Hommexx, so the forcing is already up to date
#if defined(HOMMEXX_BENCHMARK_NOFORCING)

#elif defined(SCREAM)

#elif defined(CAM)
    if (nsplit_iter == 1) then
       compute_forcing_and_push_to_c = .true.
    endif
#else
    if(test_with_forcing)then
       compute_forcing_and_push_to_c = .true.
    endif
#endif

  end function is_push_to_c_required

  function is_push_to_f_required(tl,statefreq,nextOutputStep,compute_diagnostics,nsplit_iter) &
       result(push_to_f)

    use control_mod, only : test_with_forcing, prescribed_wind
    use time_mod,    only : timelevel_t, nsplit

    type (TimeLevel_t),   intent(in) :: tl
    integer,              intent(in) :: statefreq, nextOutputStep, nsplit_iter
    logical,              intent(in) :: compute_diagnostics

    logical                          :: push_to_f, time_for_homme_output

    push_to_f = .false.
    time_for_homme_output = &
         (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep .or. compute_diagnostics)

#ifdef HOMMEXX_BENCHMARK_NOFORCING
!standalone homme, only benchmarks
    push_to_f = .false.

#elif defined(SCREAM)
!SCREAM run, only compute_diagnostics
    push_to_f = compute_diagnostics

#elif defined(CAM)
!CAM run, push at the end of nsplit loop
    if (nsplit_iter == nsplit) then
       push_to_f = .true.
    endif

!CAM also needs some of homme output
    !if (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep .or. compute_diagnostics) then
    if ( time_for_homme_output ) then
       push_to_f = .true.
    endif

#else
!standalone homme, not benchmarks
!output
    !if (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep .or. compute_diagnostics) then
    if ( time_for_homme_output ) then
       push_to_f = .true.
    endif

!push for standalone homme with forcing, test_with_forcing=false
!for most standalone homme tests
    if (test_with_forcing) then
       push_to_f = .true.
    endif

    ! In principle this shouldn't be needed, but there are roundoff-level errors
    ! that develop if this isn't true.
    if (prescribed_wind == 1) push_to_f = .true.
#endif

  end function is_push_to_f_required

  subroutine init_standalone_test(elem,deriv,hybrid,hvcoord,tl,nets,nete)
    ! set_prescribed_wind takes hvcoord as intent(inout) because it modifies it
    ! in the first call. In the C++ dycore init, we need hvcoord already
    ! established. This routine and set_prescribed_wind_f takes care of this
    ! detail.
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use element_mod,      only : element_t
    use derivative_mod,   only : derivative_t
#if !defined(CAM) && !defined(SCREAM)
    use test_mod,         only : set_test_initial_conditions
#endif

    type (element_t),      intent(inout), target  :: elem(:)
    type (derivative_t),   intent(in)             :: deriv
    type (hybrid_t),       intent(in)             :: hybrid
    type (hvcoord_t),      intent(inout)          :: hvcoord
    type (TimeLevel_t)   , intent(in)             :: tl
    integer              , intent(in)             :: nets
    integer              , intent(in)             :: nete

#if !defined(CAM) && !defined(SCREAM)
    ! Already called in prim_driver_base::prim_init2:
    !   call set_test_initial_conditions(elem,deriv,hybrid,hvcoord,tl,nets,nete)
    ! Also already taken care of:
    !   call push_test_state_to_c_wrapper()

    ! Save arguments for the C++-F90 bridge for prescribed winds.
    prescribed_wind_args%elem => elem
    prescribed_wind_args%hybrid = hybrid
    prescribed_wind_args%hvcoord = hvcoord
    prescribed_wind_args%deriv = deriv
    prescribed_wind_args%nets = nets
    prescribed_wind_args%nete = nete
#endif
  end subroutine init_standalone_test

  subroutine compute_test_forcing_f(elem,hybrid,hvcoord,nt,ntQ,dt,nets,nete,tl)
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use element_mod,      only : element_t
#if !defined(CAM) && !defined(SCREAM)
    use test_mod,       only : compute_test_forcing
#endif
    implicit none
    type(element_t),     intent(inout) :: elem(:)                            ! element array
    type(hybrid_t),      intent(in)    :: hybrid                             ! hybrid parallel structure
    type(hvcoord_t),     intent(in)    :: hvcoord
    real(kind=real_kind),intent(in)    :: dt
    integer,             intent(in)    :: nets,nete,nt,ntQ
    type(TimeLevel_t),   intent(in)    :: tl

#if !defined(CAM) && !defined(SCREAM)
    call compute_test_forcing(elem,hybrid,hvcoord,nt,ntQ,dt,nets,nete,tl)
#endif
  end subroutine compute_test_forcing_f

  subroutine push_test_state_to_c_wrapper()
#if !defined(CAM) && !defined(SCREAM)
    use iso_c_binding, only : c_ptr, c_loc
    use perf_mod,      only : t_startf, t_stopf
    use theta_f2c_mod, only : push_test_state_to_c
    use element_state, only : elem_state_v, elem_state_w_i, elem_state_vtheta_dp,     &
                              elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v,   &
                              elem_derived_eta_dot_dpdn, elem_derived_vn0

    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr
    type (c_ptr) :: elem_state_dp3d_ptr, elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_eta_dot_dpdn_ptr, elem_derived_vn0_ptr

    call t_startf('push_to_cxx')
    elem_state_v_ptr         = c_loc(elem_state_v)
    elem_state_w_i_ptr       = c_loc(elem_state_w_i)
    elem_state_vtheta_dp_ptr = c_loc(elem_state_vtheta_dp)
    elem_state_phinh_i_ptr   = c_loc(elem_state_phinh_i)
    elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
    elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
    elem_derived_vn0_ptr     = c_loc(elem_derived_vn0)
    elem_derived_eta_dot_dpdn_ptr = c_loc(elem_derived_eta_dot_dpdn)
    call push_test_state_to_c(elem_state_ps_v_ptr, elem_state_dp3d_ptr, &
         elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr, elem_state_v_ptr, &
         elem_state_w_i_ptr, elem_derived_eta_dot_dpdn_ptr, elem_derived_vn0_ptr)
    call t_stopf('push_to_cxx')
#endif
  end subroutine push_test_state_to_c_wrapper

  subroutine init_prescribed_wind_subcycle(elem, nets, nete, tl)
    ! Set the derived values used in tracer transport on the F90 side even
    ! though most of the work is done on the C++ side. This is needed because
    ! set_prescribed_wind accumulates certain derived quantities during
    ! prim_advance_exp that get repeatedly copied from F90 to C++. Here we
    ! initialize values for accumulation.
    !   In summary: Call this before entering the prim_run_subcycle loop.
    
    use prim_driver_base, only: set_tracer_transport_derived_values
    
    type (element_t), intent(inout) :: elem(:)
    integer, intent(in) :: nets, nete
    type (timelevel_t) :: tl

    call set_tracer_transport_derived_values(elem, nets, nete, tl)
  end subroutine init_prescribed_wind_subcycle

  subroutine set_prescribed_wind_f_bridge(n0, np1, nstep, dt) bind(c)
    ! This routine is called from the C++ prim_advance_exp implementation inside
    ! the prim_run_subcycle loop.
    
    use iso_c_binding, only: c_int, c_double
    
    integer(c_int), value, intent(in) :: n0, np1, nstep
    real(c_double), value, intent(in) :: dt

    type (TimeLevel_t) :: tl

    ! Only these fields need to be valid.
    tl%n0 = n0+1
    tl%np1 = np1+1
    tl%nstep = nstep

    call set_prescribed_wind_f(prescribed_wind_args%elem, prescribed_wind_args%deriv, &
         prescribed_wind_args%hybrid, prescribed_wind_args%hvcoord, dt, tl, &
         prescribed_wind_args%nets, prescribed_wind_args%nete)
  end subroutine set_prescribed_wind_f_bridge

  subroutine set_prescribed_wind_f(elem,deriv,hybrid,hvcoord,dt,tl,nets,nete)
    ! Here we finally can compute the prescribed wind in F90 and then push the
    ! data to C++.
    
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use element_mod,      only : element_t
    use derivative_mod,   only : derivative_t
#if !defined(CAM) && !defined(SCREAM)
    use control_mod,      only : qsplit
    use test_mod,         only : set_prescribed_wind
#endif

    type (element_t),      intent(inout), target  :: elem(:)
    type (derivative_t),   intent(in)             :: deriv
    type (hvcoord_t),      intent(in)             :: hvcoord
    type (hybrid_t),       intent(in)             :: hybrid
    real (kind=real_kind), intent(in)             :: dt
    type (TimeLevel_t)   , intent(in)             :: tl
    integer              , intent(in)             :: nets
    integer              , intent(in)             :: nete

#if !defined(CAM) && !defined(SCREAM)
    type (hvcoord_t) :: hv

    real(kind=real_kind) :: eta_ave_w

    ! We need to set up an hvcoord_t that can be passed as intent(inout), even
    ! though at this point, it won't be changed in the set_prescribed_wind call.
    hv = hvcoord

    eta_ave_w = 1d0/qsplit
    call set_prescribed_wind(elem,deriv,hybrid,hv,dt,tl,nets,nete,eta_ave_w)

    call push_test_state_to_c_wrapper()
#endif
  end subroutine set_prescribed_wind_f

end module
