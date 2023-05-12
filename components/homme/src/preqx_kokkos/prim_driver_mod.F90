#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
module prim_driver_mod

  use kinds,                only : real_kind
  use dimensions_mod,       only : qsize, nelemd, np, qsize
  use element_mod,          only : element_t
  use prim_driver_base,     only : deriv1, smooth_topo_datasets
  use prim_cxx_driver_base, only : prim_init1, prim_finalize
  use physical_constants,   only : scale_factor, laplacian_rigid_factor

  implicit none

  public :: prim_init2
  public :: prim_init_elements_views
  public :: prim_init_kokkos_functors
  public :: prim_run_subcycle

contains

  subroutine prim_init2(elem, hybrid, nets, nete, tl, hvcoord)
    use hybrid_mod,       only : hybrid_t
    use hybvcoord_mod,    only : hvcoord_t
    use element_mod,      only : element_t
    use time_mod,         only : timelevel_t
    use prim_driver_base, only : prim_init2_base => prim_init2

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

    ! Init the kokkos functors (and their bounday exchanges)
    call prim_init_kokkos_functors ()

    ! Init the kokkos views
    call prim_init_elements_views (elem)
  end subroutine prim_init2

  subroutine prim_create_c_data_structures (tl, hvcoord, mp)
    use iso_c_binding,    only : c_loc, c_ptr, c_bool, C_NULL_CHAR
    use time_mod,         only : TimeLevel_t
    use hybvcoord_mod,    only : hvcoord_t
    use prim_driver_base, only : deriv1
    use control_mod,      only : vert_remap_q_alg, use_cpstar, transport_alg,           &
                                 tstep_type, statefreq, rsplit, qsplit, ftype,          &
                                 prescribed_wind, limiter_option, disable_diagnostics,  &
                                 nu, nu_p, nu_q, nu_s, nu_div, nu_top, moisture,        &
                                 hypervis_order, hypervis_scaling, hypervis_subcycle,   &
                                 dt_remap_factor, dt_tracer_factor
    use preqx_f2c_mod,    only : init_reference_element_c, init_simulation_params_c, &
                                 init_hvcoord_c, init_time_level_c, init_elements_c
    !
    ! Input(s)
    !
    type (TimeLevel_t),       intent(in) :: tl
    type (hvcoord_t), target, intent(in) :: hvcoord
    real(kind=real_kind),     intent(in) :: mp(np,np)
    !
    ! Local(s)
    !
    real (kind=real_kind), dimension(np,np), target :: dvv, elem_mp

    type (c_ptr) :: hybrid_am_ptr, hybrid_ai_ptr, hybrid_bm_ptr, hybrid_bi_ptr
    
    ! Initialize the C++ reference element structure (i.e., pseudo-spectral deriv matrix and ref element mass matrix)
    dvv = deriv1%dvv
    elem_mp = mp
    call init_reference_element_c(c_loc(dvv),c_loc(elem_mp))

    ! Fill the simulation params structures in C++
    call init_simulation_params_c (vert_remap_q_alg, limiter_option, rsplit, qsplit, tstep_type,  &
                                   qsize, statefreq, nu, nu_p, nu_q, nu_s, nu_div, nu_top,        &
                                   hypervis_order, hypervis_subcycle, hypervis_scaling, ftype,    &
                                   LOGICAL(prescribed_wind==1,c_bool),                            &
                                   LOGICAL(moisture/="dry",c_bool),                               &
                                   LOGICAL(disable_diagnostics,c_bool),                           &
                                   LOGICAL(use_cpstar==1,c_bool),                                 &
                                   transport_alg, dt_remap_factor, dt_tracer_factor,              &
                                   scale_factor, laplacian_rigid_factor)

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

  subroutine prim_init_elements_views (elem)
    use iso_c_binding, only : c_ptr, c_loc
    use element_mod,   only : element_t
    use element_state, only : elem_state_v, elem_state_temp, elem_state_dp3d,       &
                              elem_state_Q, elem_state_Qdp, elem_state_ps_v,        &
                              elem_accum_qvar, elem_accum_qmass, elem_accum_q1mass, &
                              elem_accum_kener, elem_accum_pener, elem_accum_iener, &
                              elem_accum_iener_wet
    use preqx_f2c_mod, only : init_elements_2d_c, init_elements_states_c, init_diagnostics_c
    !
    ! Input(s)
    !
    type (element_t), intent(in) :: elem(:)
    !
    ! Local(s)
    !
    type (c_ptr) :: elem_D_ptr, elem_Dinv_ptr
    type (c_ptr) :: elem_spheremp_ptr, elem_rspheremp_ptr, elem_metdet_ptr
    type (c_ptr) :: elem_metinv_ptr, elem_tensorvisc_ptr, elem_vec_sph2cart_ptr
    type (c_ptr) :: elem_state_phis_ptr, elem_fcor_ptr
    type (c_ptr) :: elem_state_ps_v_ptr, elem_state_dp3d_ptr
    type (c_ptr) :: elem_state_v_ptr, elem_state_temp_ptr
    type (c_ptr) :: elem_state_q_ptr, elem_state_Qdp_ptr
    type (c_ptr) :: elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr
    type (c_ptr) :: elem_accum_iener_ptr, elem_accum_iener_wet_ptr, elem_accum_kener_ptr, elem_accum_pener_ptr

    real (kind=real_kind), target, dimension(np,np,2,2)     :: elem_D, elem_Dinv, elem_metinv, elem_tensorvisc
    real (kind=real_kind), target, dimension(np,np)         :: elem_spheremp, elem_rspheremp, elem_metdet
    real (kind=real_kind), target, dimension(np,np)         :: elem_state_phis, elem_fcor
    real (kind=real_kind), target, dimension(np,np,3,2)     :: elem_vec_sph2cart
    integer :: ie

    ! Initialize the 2d element arrays in C++
    elem_D_ptr            = c_loc(elem_D)
    elem_Dinv_ptr         = c_loc(elem_Dinv)
    elem_fcor_ptr         = c_loc(elem_fcor)
    elem_spheremp_ptr     = c_loc(elem_spheremp)
    elem_rspheremp_ptr    = c_loc(elem_rspheremp)
    elem_metdet_ptr       = c_loc(elem_metdet)
    elem_metinv_ptr       = c_loc(elem_metinv)
    elem_tensorvisc_ptr   = c_loc(elem_tensorvisc)
    elem_vec_sph2cart_ptr = c_loc(elem_vec_sph2cart)
    elem_state_phis_ptr   = c_loc(elem_state_phis)

    do ie=1,nelemd
      elem_D            = elem(ie)%D
      elem_Dinv         = elem(ie)%Dinv
      elem_fcor         = elem(ie)%fcor
      elem_spheremp     = elem(ie)%spheremp
      elem_rspheremp    = elem(ie)%rspheremp
      elem_metdet       = elem(ie)%metdet
      elem_metinv       = elem(ie)%metinv
      elem_state_phis   = elem(ie)%state%phis
      elem_tensorvisc   = elem(ie)%tensorVisc
      elem_vec_sph2cart = elem(ie)%vec_sphere2cart
      call init_elements_2d_c (ie-1,                                      &
                               elem_D_ptr, elem_Dinv_ptr, elem_fcor_ptr,  &
                               elem_spheremp_ptr, elem_rspheremp_ptr,     &
                               elem_metdet_ptr, elem_metinv_ptr,          &
                               elem_state_phis_ptr,                       &
                               elem_tensorvisc_ptr, elem_vec_sph2cart_ptr)
    enddo

    ! Initialize the 3d element arrays in C++
    elem_state_v_ptr    = c_loc(elem_state_v)
    elem_state_temp_ptr = c_loc(elem_state_temp)
    elem_state_dp3d_ptr = c_loc(elem_state_dp3d)
    elem_state_Qdp_ptr  = c_loc(elem_state_Qdp)
    elem_state_ps_v_ptr = c_loc(elem_state_ps_v)
    call init_elements_states_c (elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr,   &
                                 elem_state_Qdp_ptr, elem_state_ps_v_ptr)

    ! Initialize the diagnostics arrays in C++
    elem_state_q_ptr         = c_loc(elem_state_q)
    elem_accum_qvar_ptr      = c_loc(elem_accum_qvar)
    elem_accum_qmass_ptr     = c_loc(elem_accum_qmass)
    elem_accum_q1mass_ptr    = c_loc(elem_accum_q1mass)
    elem_accum_iener_ptr     = c_loc(elem_accum_iener)
    elem_accum_iener_wet_ptr = c_loc(elem_accum_iener_wet)
    elem_accum_kener_ptr     = c_loc(elem_accum_kener)
    elem_accum_pener_ptr     = c_loc(elem_accum_pener)
    call init_diagnostics_c (elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr,           &
                             elem_accum_q1mass_ptr, elem_accum_iener_ptr, elem_accum_iener_wet_ptr, &
                             elem_accum_kener_ptr, elem_accum_pener_ptr)
  end subroutine prim_init_elements_views

  subroutine prim_init_kokkos_functors ()
    use preqx_f2c_mod, only : init_functors_c, init_boundary_exchanges_c

    ! Initialize the C++ functors in the C++ context
    call init_functors_c ()

    ! Initialize boundary exchange structure in C++
    call init_boundary_exchanges_c ()

  end subroutine prim_init_kokkos_functors

  subroutine prim_run_subcycle(elem, hybrid, nets, nete, dt, single_column, tl, hvcoord, nsplit_iteration)
    use iso_c_binding,  only : c_int, c_ptr, c_loc
    use control_mod,    only : qsplit, rsplit, statefreq
    use dimensions_mod, only : nelemd
    use element_mod,    only : element_t
    use element_state,  only : elem_state_v, elem_state_temp, elem_state_dp3d, &
                               elem_state_Qdp, elem_state_Q, elem_state_ps_v, elem_derived_omega_p
    use hybrid_mod,     only : hybrid_t
    use hybvcoord_mod,  only : hvcoord_t
    use kinds,          only : real_kind
    use time_mod,       only : timelevel_t, nextOutputStep, nsplit
    use control_mod,    only : statefreq
    use parallel_mod,   only : abortmp
    use perf_mod,       only: t_startf, t_stopf
    use prim_state_mod, only: prim_printstate
    interface
      subroutine prim_run_subcycle_c(tstep,nstep,nm1,n0,np1,next_output_step,nsplit_iteration) bind(c)
        use iso_c_binding, only: c_int, c_double
        !
        ! Inputs
        !
        integer(kind=c_int),  intent(in) :: nstep, nm1, n0, np1, next_output_step, nsplit_iteration
        real (kind=c_double), intent(in) :: tstep
      end subroutine prim_run_subcycle_c


      subroutine cxx_push_results_to_f90(elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr, &
                                         elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr,  &
                                         elem_derived_omega_p_ptr) bind(c)
        use iso_c_binding , only : c_ptr
        !
        ! Inputs
        !
        type (c_ptr),          intent(in) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
        type (c_ptr),          intent(in) :: elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
        type (c_ptr),          intent(in) :: elem_derived_omega_p_ptr
      end subroutine cxx_push_results_to_f90
    end interface
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
    integer,              intent(in)    :: nsplit_iteration             ! nsplit_iteration = 1 .. nsplit
    !
    ! Locals
    !
    logical :: compute_diagnostics
    integer (kind=c_int) :: nstep_end, nstep_c, nm1_c, n0_c, np1_c
    type (c_ptr) :: elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr
    type (c_ptr) :: elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr
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

    call prim_run_subcycle_c(dt,nstep_c,nm1_c,n0_c,np1_c,nextOutputStep,nsplit_iteration)

    ! Set final timelevels from C into Fortran structure
    tl%nstep = nstep_c
    tl%nm1   = nm1_c + 1
    tl%n0    = n0_c  + 1
    tl%np1   = np1_c + 1

    if (MODULO(tl%nstep,statefreq)==0 .or. tl%nstep >= nextOutputStep) then
      ! Set pointers to states
      elem_state_v_ptr         = c_loc(elem_state_v)
      elem_state_temp_ptr      = c_loc(elem_state_temp)
      elem_state_dp3d_ptr      = c_loc(elem_state_dp3d)
      elem_state_Qdp_ptr       = c_loc(elem_state_Qdp)
      elem_state_Q_ptr         = c_loc(elem_state_Q)
      elem_state_ps_v_ptr      = c_loc(elem_state_ps_v)
      elem_derived_omega_p_ptr = c_loc(elem_derived_omega_p)

      ! Copy cxx arrays back to f90 structures
      call t_startf('push_to_f90')
      call cxx_push_results_to_f90(elem_state_v_ptr, elem_state_temp_ptr, elem_state_dp3d_ptr, &
                                   elem_state_Qdp_ptr, elem_state_Q_ptr, elem_state_ps_v_ptr, &
                                   elem_derived_omega_p_ptr)
      call t_stopf('push_to_f90')
    endif

    ! Print some diagnostic information

    if (compute_diagnostics) then
       call prim_printstate(elem, tl, hybrid,hvcoord,nets,nete)
    end if

  end subroutine prim_run_subcycle

end module
