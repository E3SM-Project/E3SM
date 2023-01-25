module theta_f2c_mod

interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !          Initialization routines              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Copies simulation parameters to C++ structures
  subroutine init_simulation_params_c (remap_alg, limiter_option, rsplit, qsplit, time_step_type,    &
                                       qsize, state_frequency, nu, nu_p, nu_q, nu_s, nu_div, nu_top, &
                                       hypervis_order, hypervis_subcycle, hypervis_subcycle_tom,     &
                                       hypervis_scaling,                                             &
                                       dcmip16_mu, ftype, theta_adv_form, prescribed_wind, moisture, &
                                       disable_diagnostics, use_cpstar, transport_alg,               &
                                       theta_hydrostatic_mode, test_case_name, dt_remap_factor,      &
                                       dt_tracer_factor, rearth, nsplit, pgrad_correction,           &
                                       dp3d_thresh, vtheta_thresh) bind(c)

    use iso_c_binding, only: c_int, c_bool, c_double, c_ptr
    !
    ! Inputs
    !
    integer(kind=c_int),  intent(in) :: remap_alg, limiter_option, rsplit, qsplit, time_step_type, nsplit
    integer(kind=c_int),  intent(in) :: dt_remap_factor, dt_tracer_factor, transport_alg
    integer(kind=c_int),  intent(in) :: state_frequency, qsize
    real(kind=c_double),  intent(in) :: nu, nu_p, nu_q, nu_s, nu_div, nu_top, hypervis_scaling, dcmip16_mu, &
                                        rearth, dp3d_thresh, vtheta_thresh
    integer(kind=c_int),  intent(in) :: hypervis_order, hypervis_subcycle, hypervis_subcycle_tom
    integer(kind=c_int),  intent(in) :: ftype, theta_adv_form
    logical(kind=c_bool), intent(in) :: prescribed_wind, moisture, disable_diagnostics, use_cpstar
    logical(kind=c_bool), intent(in) :: theta_hydrostatic_mode, pgrad_correction
    type(c_ptr), intent(in) :: test_case_name
  end subroutine init_simulation_params_c

  ! Creates element structures in C++
  subroutine init_elements_c (nelemd) bind(c)
    use iso_c_binding, only: c_int
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: nelemd
  end subroutine init_elements_c

  ! Initialize hybrid vertical coordinate in C++ from f90 values
  subroutine init_hvcoord_c (ps0,hybrid_am_ptr,hybrid_ai_ptr,hybrid_bm_ptr,hybrid_bi_ptr) bind(c)
    use iso_c_binding, only: c_double, c_ptr
    !
    ! Inputs
    !
    real (kind=c_double),  intent(in) :: ps0
    type (c_ptr),          intent(in) :: hybrid_am_ptr, hybrid_ai_ptr
    type (c_ptr),          intent(in) :: hybrid_bm_ptr, hybrid_bi_ptr
  end subroutine init_hvcoord_c

  ! Initialize time level structure in C++
  subroutine init_time_level_c(nm1,n0,np1,nstep,nstep0) bind(c)
    use iso_c_binding, only: c_int
    !
    ! Inputs
    !
    integer(kind=c_int), intent(in) :: nm1, n0, np1, nstep, nstep0
  end subroutine init_time_level_c

  ! Copies constant geometry arrays (e.g., metric, jacobian,...) from f90 arrays into C++ views
  subroutine init_elements_2d_c (ie, D_ptr, Dinv_ptr, elem_fcor_ptr,      &
                                 elem_spheremp_ptr, elem_rspheremp_ptr,   &
                                 elem_metdet_ptr, elem_metinv_ptr,        &
                                 tensorvisc_ptr, vec_sph2cart_ptr,        &
                                 sphere_cart_vec, sphere_latlon_vec) bind(c)
    use iso_c_binding, only: c_int, c_ptr, c_double
    use dimensions_mod, only : np
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ie
    type (c_ptr) , intent(in) :: D_ptr, Dinv_ptr, elem_fcor_ptr
    type (c_ptr) , intent(in) :: elem_spheremp_ptr, elem_rspheremp_ptr
    type (c_ptr) , intent(in) :: elem_metdet_ptr, elem_metinv_ptr
    type (c_ptr) , intent(in) :: tensorvisc_ptr, vec_sph2cart_ptr
    real (kind=c_double), intent(in) :: sphere_cart_vec(3,np,np), sphere_latlon_vec(2,np,np)
  end subroutine init_elements_2d_c

  ! Copies geopotential from f90 arrays to C++ views
  subroutine init_geopotential_c (ie, phis_ptr, gradphis_ptr) bind(c)
    use iso_c_binding, only: c_int, c_ptr
    !
    ! Inputs
    !
    integer (kind=c_int), intent(in) :: ie
    type (c_ptr) , intent(in) :: phis_ptr, gradphis_ptr
  end subroutine init_geopotential_c

  ! Initializes C++ diagnostics arrays with ptrs provided from f90
  subroutine init_diagnostics_c (elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr, &
                                 elem_accum_q1mass_ptr, elem_accum_iener_ptr,                 &
                                 elem_accum_kener_ptr, elem_accum_pener_ptr) bind(c)
    use iso_c_binding, only: c_ptr
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr
    type (c_ptr), intent(in) :: elem_accum_q1mass_ptr, elem_accum_iener_ptr
    type (c_ptr), intent(in) :: elem_accum_kener_ptr, elem_accum_pener_ptr
  end subroutine init_diagnostics_c

  ! Copies states from f90 arrays into C++ views
  subroutine init_elements_states_c (elem_state_v_ptr, elem_state_w_i, elem_state_vtheta_dp_ptr, &
                                     elem_state_phinh_i, elem_state_dp3d_ptr,                    &
                                     elem_state_Qdp_ptr, elem_state_ps_v_ptr) bind(c)
    use iso_c_binding, only: c_ptr
    !
    ! Inputs
    !
    type (c_ptr) :: elem_state_v_ptr, elem_state_w_i, elem_state_vtheta_dp_ptr
    type (c_ptr) :: elem_state_phinh_i, elem_state_dp3d_ptr
    type (c_ptr) :: elem_state_Qdp_ptr, elem_state_ps_v_ptr
  end subroutine init_elements_states_c

  ! Copies reference states from f90 arrays into C++ views
  subroutine init_reference_states_c (elem_theta_ref_ptr, elem_dp_ref_ptr, elem_phi_ref_ptr) bind(c)
    use iso_c_binding, only: c_ptr
    !
    ! Inputs
    !
    type (c_ptr) :: elem_theta_ref_ptr, elem_dp_ref_ptr, elem_phi_ref_ptr
  end subroutine init_reference_states_c

  ! Initialize SEM reference element structures (mass and pseudo-spectral deriv matrices)
  subroutine init_reference_element_c (deriv_ptr, mass_ptr) bind(c)
    use iso_c_binding, only: c_ptr
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: deriv_ptr, mass_ptr
  end subroutine init_reference_element_c

  ! Create C++ functors
  subroutine init_functors_c (allocate_buffer) bind(c)
  use iso_c_binding, only: c_bool
  !
  ! Inputs
  !
  logical(kind=c_bool), intent(in) :: allocate_buffer
  end subroutine init_functors_c

  ! Initialize C++ boundary exchange structures
  subroutine init_boundary_exchanges_c () bind(c)
  end subroutine init_boundary_exchanges_c

  ! Initialize dp3d from ps_v and hybrid v coord
  subroutine initialize_dp3d_from_ps_c () bind(c)
  end subroutine initialize_dp3d_from_ps_c

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !               Run-time routines               !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Copy forcing from f90 arrays into C++ views
  subroutine push_forcing_to_c(FM, FVTheta, FT, FPHI, FQ) bind(c)
    use iso_c_binding, only: c_double
    use dimensions_mod, only: np, nlev, nlevp, nelemd, qsize_d
    !
    ! Inputs
    !
     real (kind=c_double), intent(in) :: FM(np,np,3,nlev,nelemd)
     real (kind=c_double), intent(in) :: FVTheta(np,np,nlev,nelemd), FT(np,np,nlev,nelemd)
     real (kind=c_double), intent(in) :: FPHI(np,np,nlevp,nelemd), FQ(np,np,nlev,qsize_d,nelemd)
  end subroutine push_forcing_to_c

  ! Run dycore for a full atm timesteps
  subroutine prim_run_subcycle_c(tstep,nstep,nm1,n0,np1,next_output_step,nsplit_iter) bind(c)
    use iso_c_binding, only: c_int, c_double
    !
    ! Inputs
    !
    integer(kind=c_int),  intent(inout) :: nstep, nm1, n0, np1
    integer(kind=c_int),  intent(in)    :: next_output_step, nsplit_iter
    real (kind=c_double), intent(in)    :: tstep
  end subroutine prim_run_subcycle_c

  ! Copy results from C++ views back to f90 arrays
  subroutine cxx_push_results_to_f90(elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr,   &
                                     elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr, &
                                     elem_state_Qdp_ptr, elem_state_Q_ptr, elem_derived_omega_p_ptr) bind(c)
    use iso_c_binding, only: c_ptr
    !
    ! Inputs
    !
    type (c_ptr), intent(in) :: elem_state_v_ptr, elem_state_w_i_ptr, elem_state_vtheta_dp_ptr
    type (c_ptr), intent(in) :: elem_state_phinh_i_ptr, elem_state_dp3d_ptr, elem_state_ps_v_ptr
    type (c_ptr), intent(in) :: elem_state_Qdp_ptr, elem_state_Q_ptr, elem_derived_omega_p_ptr
  end subroutine cxx_push_results_to_f90

  subroutine push_test_state_to_c( &
       ! state
       ps_v, dp3d, vtheta_dp, phinh_i, v, w_i, &
       ! derived
       eta_dot_dpdn, vn0) bind(c)
    use iso_c_binding, only: c_ptr
    
    type (c_ptr), intent(in) :: ps_v, dp3d, vtheta_dp, phinh_i, v, w_i, eta_dot_dpdn, vn0
  end subroutine push_test_state_to_c
end interface

end module theta_f2c_mod
