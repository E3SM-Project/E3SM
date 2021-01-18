#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module prim_driver_mod

  use kinds,                only : real_kind
  use dimensions_mod,       only : qsize, nelemd, np, qsize
  use element_mod,          only : element_t
  use prim_driver_base,     only : deriv1, smooth_topo_datasets
  use prim_cxx_driver_base, only : prim_init1, prim_finalize

  implicit none

  public :: prim_init2
  public :: prim_run_subcycle
  public :: prim_init_elements_views
  public :: prim_init_grid_views
  public :: prim_init_geopotential_views
  public :: prim_init_state_views
  public :: prim_init_ref_states_views
  public :: prim_init_diags_views

contains

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use time_mod,         only : timelevel_t
    use prim_driver_base, only : deriv1, prim_init2_base => prim_init2
    use prim_state_mod,   only : prim_printstate
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
  end subroutine prim_init2

  subroutine prim_create_c_data_structures (tl, hvcoord, mp)
    use iso_c_binding, only : c_loc, c_ptr, c_bool, C_NULL_CHAR
    use theta_f2c_mod, only : init_reference_element_c, init_simulation_params_c, &
                              init_time_level_c, init_hvcoord_c, init_elements_c
    use time_mod,      only : TimeLevel_t
    use hybvcoord_mod, only : hvcoord_t
    use control_mod,   only : limiter_option, rsplit, qsplit, tstep_type, statefreq,  &
                              nu, nu_p, nu_q, nu_s, nu_div, nu_top, vert_remap_q_alg, &
                              hypervis_order, hypervis_subcycle, hypervis_scaling,    &
                              ftype, prescribed_wind, moisture, disable_diagnostics,  &
                              use_cpstar, transport_alg, theta_hydrostatic_mode,      &
                              dcmip16_mu, theta_advect_form, test_case, MAX_STRING_LEN
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
    logical (kind=c_bool) :: use_semi_lagrange_transport
    real (kind=real_kind), target :: dvv (np,np), elem_mp(np,np)
    type (c_ptr) :: hybrid_am_ptr, hybrid_ai_ptr, hybrid_bm_ptr, hybrid_bi_ptr
    character(len=MAX_STRING_LEN), target :: test_name

    ! Initialize the C++ reference element structure (i.e., pseudo-spectral deriv matrix and ref element mass matrix)
    dvv = deriv1%dvv
    elem_mp = mp
    call init_reference_element_c(c_loc(dvv),c_loc(elem_mp))

    ! Fill the simulation params structures in C++
    use_semi_lagrange_transport = transport_alg > 0
    test_name = TRIM(test_case) // C_NULL_CHAR
    call init_simulation_params_c (vert_remap_q_alg, limiter_option, rsplit, qsplit, tstep_type,  &
                                   qsize, statefreq, nu, nu_p, nu_q, nu_s, nu_div, nu_top,        &
                                   hypervis_order, hypervis_subcycle, hypervis_scaling,           &
                                   dcmip16_mu, ftype, theta_advect_form,                          &
                                   LOGICAL(prescribed_wind==1,c_bool),                            &
                                   LOGICAL(moisture/="dry",c_bool),                               &
                                   LOGICAL(disable_diagnostics,c_bool),                           &
                                   LOGICAL(use_cpstar==1,c_bool),                                 &
                                   LOGICAL(use_semi_lagrange_transport,c_bool),                   &
                                   LOGICAL(theta_hydrostatic_mode,c_bool),                        &
                                   c_loc(test_name))

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

    integer :: ie

    elem_D_ptr            = c_loc(elem_D)
    elem_Dinv_ptr         = c_loc(elem_Dinv)
    elem_fcor_ptr         = c_loc(elem_fcor)
    elem_spheremp_ptr     = c_loc(elem_spheremp)
    elem_rspheremp_ptr    = c_loc(elem_rspheremp)
    elem_metdet_ptr       = c_loc(elem_metdet)
    elem_metinv_ptr       = c_loc(elem_metinv)
    elem_tensorvisc_ptr   = c_loc(elem_tensorvisc)
    elem_vec_sph2cart_ptr = c_loc(elem_vec_sph2cart)

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
      call init_elements_2d_c (ie-1,                                      &
                               elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr,  &
                               elem_spheremp_ptr, elem_rspheremp_ptr,     &
                               elem_metdet_ptr, elem_metinv_ptr,          &
                               elem_tensorvisc_ptr, elem_vec_sph2cart_ptr)
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

  subroutine prim_init_kokkos_functors ()
    use theta_f2c_mod, only : init_functors_c, init_boundary_exchanges_c

    ! Initialize the C++ functors in the C++ context
    call init_functors_c ()

    ! Initialize boundary exchange structure in C++
    call init_boundary_exchanges_c ()

  end subroutine prim_init_kokkos_functors

  subroutine prim_run_subcycle(elem, hybrid, nets, nete, dt, single_column, tl, hvcoord,nsubstep)
    use iso_c_binding,  only : c_int, c_ptr, c_loc
    use control_mod,    only : qsplit, rsplit, statefreq, disable_diagnostics
    use dimensions_mod, only : nelemd
    use element_state,  only : elem_state_v, elem_state_w_i, elem_state_vtheta_dp,     &
                               elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v,   &
                               elem_state_Qdp, elem_state_Q, elem_derived_omega_p,     &
                               elem_derived_FM, elem_derived_FVTheta, elem_derived_FT, &
                               elem_derived_FPHI, elem_derived_FQ
    use hybrid_mod,     only : hybrid_t
    use hybvcoord_mod,  only : hvcoord_t
    use kinds,          only : real_kind
    use time_mod,       only : timelevel_t, nextOutputStep, nsplit
    use control_mod,    only : statefreq
    use parallel_mod,   only : abortmp
    use perf_mod,       only : t_startf, t_stopf
    use prim_state_mod, only : prim_printstate
    use theta_f2c_mod,  only : prim_run_subcycle_c, cxx_push_results_to_f90
#ifndef SCREAM
    use theta_f2c_mod,  only : push_forcing_to_c
#endif
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
    integer,              intent(in)    :: nsubstep                     ! nsubstep = 1 .. nsplit
    !
    ! Locals
    !
    logical :: compute_diagnostics
    integer (kind=c_int) :: nstep_end, nstep_c, nm1_c, n0_c, np1_c
    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr, elem_state_phinh_i_ptr
    type (c_ptr) :: elem_state_dp3d_ptr, elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
    type (c_ptr) :: elem_derived_omega_p_ptr

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

#ifndef SCREAM
    ! Scream already computes all forcing using the same pointers
    ! stored in Hommexx, so the forcing is already up to date
    call t_startf('push_to_cxx')
    call push_forcing_to_c(elem_derived_FM,   elem_derived_FVTheta, elem_derived_FT, &
                           elem_derived_FPHI, elem_derived_FQ)
    call t_stopf('push_to_cxx')
#endif

    call prim_run_subcycle_c(dt,nstep_c,nm1_c,n0_c,np1_c,nextOutputStep)

    ! Set final timelevels from C into Fortran structure
    tl%nstep = nstep_c
    tl%nm1   = nm1_c + 1
    tl%n0    = n0_c  + 1
    tl%np1   = np1_c + 1

    if (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep .or. compute_diagnostics) then
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
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle

  subroutine setup_element_pointers (elem)
    use element_state,  only : allocate_element_arrays, elem_state_v, elem_state_w_i, elem_state_vtheta_dp, &
                               elem_state_phinh_i, elem_state_dp3d, elem_state_ps_v, elem_state_phis,       & 
                               elem_state_Qdp, elem_state_Q, elem_derived_omega_p,                          &
                               elem_accum_pener, elem_accum_kener, elem_accum_iener,                        &
                               elem_accum_qvar, elem_accum_qmass, elem_accum_q1mass
    !
    ! Inputs
    !
    type (element_t), intent(inout) :: elem(:)
    !
    ! Locals
    !
    integer :: ie

    call allocate_element_arrays(nelemd)

    do ie=1,nelemd
      elem(ie)%state%v         => elem_state_v(:,:,:,:,:,ie)
      elem(ie)%state%w_i       => elem_state_w_i(:,:,:,:,ie)
      elem(ie)%state%vtheta_dp => elem_state_vtheta_dp(:,:,:,:,ie)
      elem(ie)%state%phinh_i   => elem_state_phinh_i(:,:,:,:,ie)
      elem(ie)%state%dp3d      => elem_state_dp3d(:,:,:,:,ie)
      elem(ie)%state%ps_v      => elem_state_ps_v(:,:,:,ie)
      elem(ie)%state%Q         => elem_state_Q(:,:,:,:,ie)
      elem(ie)%state%Qdp       => elem_state_Qdp(:,:,:,:,:,ie)
      elem(ie)%state%phis      => elem_state_phis(:,:,ie)
      elem(ie)%derived%omega_p => elem_derived_omega_p(:,:,:,ie)

      elem(ie)%accum%KEner     => elem_accum_KEner    (:,:,:,ie)
      elem(ie)%accum%PEner     => elem_accum_PEner    (:,:,:,ie)
      elem(ie)%accum%IEner     => elem_accum_IEner    (:,:,:,ie)
      elem(ie)%accum%Qvar      => elem_accum_Qvar     (:,:,:,:,ie)
      elem(ie)%accum%Qmass     => elem_accum_Qmass    (:,:,:,:,ie)
      elem(ie)%accum%Q1mass    => elem_accum_Q1mass   (:,:,:,ie)
    enddo
  end subroutine setup_element_pointers

end module
