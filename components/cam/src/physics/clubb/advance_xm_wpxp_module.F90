!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_xm_wpxp_module

  ! Description:
  ! Contains the CLUBB advance_xm_wpxp_module scheme.

  ! References:
  ! None
  !-----------------------------------------------------------------------

  implicit none

  private ! Default scope

  public  :: advance_xm_wpxp

  private :: xm_wpxp_lhs, & 
             xm_wpxp_rhs, & 
             xm_wpxp_solve, & 
             xm_wpxp_clipping_and_stats, &
             xm_term_ta_lhs, & 
             xm_term_ta_lhs_all, & 
             wpxp_term_tp_lhs, & 
             wpxp_term_tp_lhs_all, & 
             wpxp_terms_ac_pr2_lhs, & 
             wpxp_terms_ac_pr2_lhs_all, &
             wpxp_term_pr1_lhs, & 
             wpxp_term_pr1_lhs_all, & 
             wpxp_terms_bp_pr3_rhs, &
             wpxp_terms_bp_pr3_rhs_all, &
             xm_correction_wpxp_cl, &
             damp_coefficient, &
             diagnose_upxp, &
             error_prints_xm_wpxp

  ! Parameter Constants
  integer, parameter, private :: & 
    nsub = 2, & ! Number of subdiagonals in the LHS matrix
    nsup = 2, & ! Number of superdiagonals in the LHS matrix
    xm_wpxp_thlm = 1,   & ! Named constant for thlm and wpthlp solving
    xm_wpxp_rtm = 2,    & ! Named constant for rtm and wprtp solving
    xm_wpxp_scalar = 3, & ! Named constant for sclrm and wpsclrp solving
    xm_wpxp_um = 4,     & ! Named constant for optional um and upwp solving
    xm_wpxp_vm = 5        ! Named constant for optional vm and vpwp solving

  contains

  !=============================================================================
  subroutine advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                              Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, &
                              tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, rtm_forcing, &
                              wprtp_forcing, rtm_ref, wp2thlp, thlpthvp, &
                              thlm_forcing, wpthlp_forcing, thlm_ref, &
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                              invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
                              w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
                              mixt_frac_zm, l_implemented, em, wp2sclrp, &
                              sclrpthvp, sclrm_forcing, sclrp2, exner, rcm, &
                              p_in_Pa, thvm, Cx_fnc_Richardson, &
                              pdf_implicit_coefs_terms, &
                              um_forcing, vm_forcing, ug, vg, wpthvp, &
                              fcor, um_ref, vm_ref, up2, vp2, &
                              uprcp, vprcp, rc_coef, & 
                              rtm, wprtp, thlm, wpthlp, &
                              sclrm, wpsclrp, um, upwp, vm, vpwp )

    ! Description:
    ! Advance the mean and flux terms by one timestep.

    ! References:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_eqns
    !
    ! Eqn. 16 & 17 on p. 3546 of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.

    ! See Also
    ! ``Equations for CLUBB'' Section 5:
    !   /Implicit solutions for the means and fluxes/
    !-----------------------------------------------------------------------

    use parameters_tunable, only:  & 
        C6rt,  & ! Variable(s)
        C6rtb,  & 
        C6rtc,  & 
        C6thl,  & 
        C6thlb,  & 
        C6thlc, & 
        C7,  & 
        C7b,  & 
        C7c,  & 
        c_K6,  & 
        C6rt_Lscale0, &
        C6thl_Lscale0, &
        C7_Lscale0, &
        wpxp_L_thresh

    use constants_clubb, only:  & 
        fstderr, &  ! Constant
        rt_tol, &
        thl_tol, &
        w_tol, &
        w_tol_sqd, &
        thl_tol_mfl, &
        rt_tol_mfl, &
        max_mag_correlation, &
        one, &
        one_half, &
        zero, &
        zero_threshold, &
        eps, &
        ep1

    use parameters_model, only: & 
        sclr_dim, &  ! Variable(s)
        sclr_tol, &
        ts_nudge

    use grid_class, only: & 
        gr,   & ! Variable(s)
        ddzt    ! Procedure(s)

    use grid_class, only: &
        zm2zt, & ! Procedure(s)
        zt2zm

    use model_flags, only: &
        l_clip_semi_implicit,          & ! Variable(s)
        l_use_C7_Richardson,           &
        l_explicit_turbulent_adv_wpxp, &
        l_upwind_wpxp_ta,              &
        l_predict_upwp_vpwp,           &
        l_uv_nudge

    use mono_flux_limiter, only: &
        calc_turb_adv_range ! Procedure(s)

    use pdf_closure_module, only: &
        iiPDF_new,  & ! Variable(s)
        iiPDF_ADG1, &
        iiPDF_type

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use turbulent_adv_pdf, only: &
        sgn_turbulent_velocity    ! Procedure(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constants

    use stats_type_utilities, only: &
        stat_begin_update, & ! Procedure(s)
        stat_end_update, &
        stat_update_var

    use stats_variables, only: & 
        stats_zt, &
        stats_zm, &
        irtm_matrix_condt_num, &  ! Variables
        ithlm_matrix_condt_num, &
        irtm_sdmp, &
        ithlm_sdmp, & 
        ium_sdmp, &
        ivm_sdmp, &
        ium_ndg, &
        ivm_ndg, &
        ium_ref, &
        ivm_ref, &
        ium_gf, &
        ium_cf, &
        ivm_gf, &
        ivm_cf, &
        ium_f, &
        ivm_f, &
        iupwp_pr4, &
        ivpwp_pr4, &
        iupthlp, &
        iuprtp,  &
        ivpthlp, &
        ivprtp,  &
        iupthvp, &
        iuprcp,  &
        ivpthvp, &
        ivprcp,  &
        iC7_Skw_fnc, &
        iC6rt_Skw_fnc, &
        iC6thl_Skw_fnc, &
        icoef_wp2rtp_implicit, &
        iterm_wp2rtp_explicit, &
        icoef_wp2thlp_implicit, &
        iterm_wp2thlp_explicit, &
        l_stats_samp

    use sponge_layer_damping, only: &
        rtm_sponge_damp_settings, &
        thlm_sponge_damp_settings, &
        uv_sponge_damp_settings, &
        rtm_sponge_damp_profile, &
        thlm_sponge_damp_profile, &
        uv_sponge_damp_profile, &
        sponge_damp_xm ! Procedure(s)

    use advance_helper_module, only: &
        compute_Cx_fnc_Richardson ! Procedure

    implicit none

    ! External
    intrinsic :: exp, sqrt

    ! Parameter Constants
    logical, parameter :: &
      l_iter = .true. ! True when the means and fluxes are prognosed

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      sigma_sqd_w,     & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,           & ! w wind component on momentum levels      [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      Lscale,          & ! Turbulent mixing length                  [m]
      em,              & ! Turbulent Kinetic Energy (TKE)           [m^2/s^2]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on momentum levels    [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels     [m/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      Kh_zm,           & ! Eddy diffusivity on momentum levels
      tau_C6_zm,       & ! Time-scale tau on momentum levels applied to C6 term [s]
      Skw_zm,          & ! Skewness of w on momentum levels         [-]
      wp2rtp,          & ! <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
      rtm_ref,         & ! rtm for nudging                          [kg/kg]
      wp2thlp,         & ! <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
      thlpthvp,        & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,    & ! th_l forcing (thermodynamic levels)      [K/s]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)      [K/s^2]
      thlm_ref,        & ! thlm for nudging                         [K]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on moment. levs. [K]
      ! Added for clipping by Vince Larson 29 Sep 2007
      rtp2,            & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2,           & ! th_l'^2 (momentum levels)                [K^2]
      ! End of Vince Larson's addition.
      w_1_zm,          & ! Mean w (1st PDF component)              [m/s]
      w_2_zm,          & ! Mean w (2nd PDF component)              [m/s]
      varnce_w_1_zm,   & ! Variance of w (1st PDF component)       [m^2/s^2]
      varnce_w_2_zm,   & ! Variance of w (2nd PDF component)       [m^2/s^2]
      mixt_frac_zm      ! Weight of 1st PDF component (Sk_w dependent) [-]

    logical, intent(in) ::  & 
      l_implemented      ! Flag for CLUBB being implemented in a larger model.

    ! Additional variables for passive scalars
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: & 
      wp2sclrp,      & ! <w'^2 sclr'> (thermodynamic levels)   [Units vary]
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      exner,           & ! Exner function                            [-]
      rcm,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,         & ! Air pressure                              [Pa]
      thvm,            & ! Virutal potential temperature             [K]
      Cx_fnc_Richardson  ! Cx_fnc computed from Richardson_num       [-]

    type(implicit_coefs_terms), dimension(gr%nz), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      um_forcing, & ! <u> forcing term (thermodynamic levels)      [m/s^2]
      vm_forcing, & ! <v> forcing term (thermodynamic levels)      [m/s^2]
      ug,         & ! <u> geostrophic wind (thermodynamic levels)  [m/s]
      vg,         & ! <v> geostrophic wind (thermodynamic levels)  [m/s]
      wpthvp        ! <w'thv'> (momentum levels)                   [m/s K]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      uprcp,              & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp,              & ! < v' r_c' >              [(m kg)/(s kg)]
      rc_coef               ! Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]

     real( kind = core_rknd ), intent(in) ::  &
      fcor          ! Coriolis parameter                           [s^-1]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      um_ref, & ! Reference u wind component for nudging       [m/s]
      vm_ref, & ! Reference v wind component for nudging       [m/s]
      up2,    & ! Variance of the u wind component             [m^2/s^2]
      vp2       ! Variance of the v wind component             [m^2/s^2]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      um,   & ! <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
      upwp, & ! <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      vm,   & ! <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
      vpwp    ! <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]

    ! Local variables
    real( kind = core_rknd ), dimension(nsup+nsub+1,2*gr%nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

    ! Additional variables for passive scalars
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: & 
      wpsclrp_forcing    ! <w'sclr'> forcing (momentum levels)  [m/s{un vary}]

    ! Eddy Diffusion for wpthlp and wprtp.
    real( kind = core_rknd ), dimension(gr%nz) :: Kw6  ! wpxp eddy diff. [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      a1,    & ! a_1 (momentum levels); See eqn. 24 in `Equations for CLUBB' [-]
      a1_zt    ! a_1 interpolated to thermodynamic levels                    [-]

    ! Variables for turbulent advection of predictive variances and covariances.

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2rtp_implicit_zm, & ! coef_wp2rtp_implicit interp. to m-levs. [m/s]
      term_wp2rtp_explicit_zm    ! term_wp2rtp_expl intrp m-levs [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2thlp_implicit_zm, & ! coef_wp2thlp_implicit interp. m-levs.  [m/s]
      term_wp2thlp_explicit_zm    ! term_wp2thlp_expl interp. m-levs [m^2/s^2 K]

    ! <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'> + term_wp2sclrp_explicit
    real ( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      coef_wp2sclrp_implicit, & ! Coef. that is multiplied by <w'sclr'>    [m/s]
      term_wp2sclrp_explicit    ! Term that is on the RHS    [m^2/s^2(un. vary)]

    real ( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      coef_wp2sclrp_implicit_zm, & ! coef_wp2sclrp_implicit interp. m-levs [m/s]
      term_wp2sclrp_explicit_zm    ! term_wp2sclrp_expl intrp zm [m^2/s^2(un v)]

    ! Sign of turbulent velocity (used for "upwind" turbulent advection)
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      sgn_t_vel_wprtp,  & ! Sign of the turbulent velocity for <w'rt'>       [-]
      sgn_t_vel_wpthlp    ! Sign of the turbulent velocity for <w'thl'>      [-]

    real ( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      sgn_t_vel_wpsclrp    ! Sign of the turbulent velocity for <w'sclr'>    [-]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.

    ! <w'^2 u'> = coef_wp2up_implicit * <u'w'> + term_wp2up_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2up_implicit, & ! Coefficient that is multiplied by <u'w'>  [m/s]
      term_wp2up_explicit    ! Term that is on the RHS               [m^3/s^3]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2up_implicit_zm, & ! coef_wp2up_implicit interp. to m-levs. [m/s]
      term_wp2up_explicit_zm    ! term_wp2up_expl interp. to m-levs. [m^3/s^3]

    ! <w'^2 v'> = coef_wp2vp_implicit * <v'w'> + term_wp2vp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2vp_implicit, & ! Coefficient that is multiplied by <v'w'>  [m/s]
      term_wp2vp_explicit    ! Term that is on the RHS               [m^3/s^3]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2vp_implicit_zm, & ! coef_wp2vp_implicit interp. to m-levs. [m/s]
      term_wp2vp_explicit_zm    ! term_wp2vp_expl interp. to m-levs. [m^3/s^3]

    ! Sign of turbulent velocity (used for "upwind" turbulent advection)
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      sgn_t_vel_upwp, & ! Sign of the turbulent velocity for <u'w'>       [-]
      sgn_t_vel_vpwp    ! Sign of the turbulent velocity for <v'w'>       [-]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      um_tndcy,     & ! <u> forcing term + coriolis (thermo levs)        [m/s^2]
      vm_tndcy,     & ! <v> forcing term + coriolis (thermo levs)        [m/s^2]
      upwp_forcing, & ! <u'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
      vpwp_forcing, & ! <v'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
      upthvp,       & ! <u'thv'> (momentum levels)                       [m/s K]
      vpthvp,       & ! <v'thv'> (momentum levels)                       [m/s K]
      upthlp,       & ! eastward horz turb flux of theta_l (mom levs)    [m/s K]
      vpthlp,       & ! northward horz turb flux of theta_l (mom levs)   [m/s K]
      uprtp,        & ! eastward horz turb flux of tot water (mom levs)  [m/s kg/kg]
      vprtp           ! northward horz turb flux of tot water (mom levs) [m/s kg/kg]

    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(gr%nz) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.

    ! Variables used for clipping of w'x' due to correlation
    ! of w with x, such that:
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
    ! -1 <= corr_(w,x) <= 1.
    real( kind = core_rknd ), dimension(gr%nz) :: & 
      wpxp_upper_lim, & ! Keeps correlations from becoming greater than 1.
      wpxp_lower_lim    ! Keeps correlations from becoming less than -1.

    real( kind = core_rknd ), dimension(gr%nz) :: &
      zeros_vector  ! Array of zeros, of the size of a vertical profile [-]

    real( kind = core_rknd ), allocatable, dimension(:,:) :: & 
      rhs,      & ! Right-hand sides of band diag. matrix. (LAPACK)
      rhs_save, & ! Saved Right-hand sides of band diag. matrix. (LAPACK)
      solution    ! solution vectors of band diag. matrix. (LAPACK)

    ! Constant parameters as a function of Skw.

    integer :: &
      nrhs         ! Number of RHS vectors

    real( kind = core_rknd ) :: rcond

    ! Saved values of predictive fields, prior to being advanced, for use in
    ! print statements in case of fatal error.
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      rtm_old,    & ! Saved value of r_t        [kg/kg]
      wprtp_old,  & ! Saved value of w'r_t'     [(kg/kg) m/s]
      thlm_old,   & ! Saved value of th_l       [K]
      wpthlp_old    ! Saved value of w'th_l'    [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) ::  & 
      sclrm_old, wpsclrp_old !                  [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      um_old,   & ! Saved value of <u>       [m/s]
      upwp_old, & ! Saved value of <u'w'>    [m^2/s^2]
      vm_old,   & ! Saved value of <v>       [m/s]
      vpwp_old    ! Saved value of <v'w'>    [m^2/s^2]

    ! Indices
    integer :: i, k

    !---------------------------------------------------------------------------

    ! ----- Begin Code -----
    if ( l_clip_semi_implicit &
         .or. ( ( iiPDF_type == iiPDF_new ) &
                .and. ( .not. l_explicit_turbulent_adv_wpxp ) ) ) then
       nrhs = 1
    else
       if ( l_predict_upwp_vpwp ) then
          nrhs = 4+sclr_dim
       else
          nrhs = 2+sclr_dim
       endif
    endif

    ! Allocate rhs and solution vector
    allocate( rhs(2*gr%nz,nrhs) )
    allocate( rhs_save(2*gr%nz,nrhs) )
    allocate( solution(2*gr%nz,nrhs) )

    ! This is initialized solely for the purpose of avoiding a compiler
    ! warning about uninitialized variables.
    zeros_vector = zero

    ! Save values of predictive fields to be printed in case of crash.
    rtm_old = rtm
    wprtp_old = wprtp
    thlm_old = thlm
    wpthlp_old = wpthlp
    if ( sclr_dim > 0 ) then
       sclrm_old = sclrm
       wpsclrp_old = wpsclrp
    endif ! sclr_dim > 0
    if ( l_predict_upwp_vpwp ) then
       um_old = um
       upwp_old = upwp
       vm_old = vm
       vpwp_old = vpwp
    endif ! l_predict_upwp_vpwp

    ! Compute C6 and C7 as a function of Skw
    ! The if...then is just here to save compute time
    if ( abs(C6rt-C6rtb) > abs(C6rt+C6rtb)*eps/2 ) then
      C6rt_Skw_fnc(1:gr%nz) = C6rtb + (C6rt-C6rtb) & 
        *EXP( -one_half * (Skw_zm(1:gr%nz)/C6rtc)**2 )
    else
      C6rt_Skw_fnc(1:gr%nz) = C6rtb
    endif

    if ( abs(C6thl-C6thlb) > abs(C6thl+C6thlb)*eps/2 ) then
      C6thl_Skw_fnc(1:gr%nz) = C6thlb + (C6thl-C6thlb) & 
        *EXP( -one_half * (Skw_zm(1:gr%nz)/C6thlc)**2 )
    else
      C6thl_Skw_fnc(1:gr%nz) = C6thlb
    endif

    ! Compute C7_Skw_fnc
    if ( l_use_C7_Richardson ) then
      ! New formulation based on Richardson number
      C7_Skw_fnc = Cx_fnc_Richardson
    else
      if ( abs(C7-C7b) > abs(C7+C7b)*eps/2 ) then
        C7_Skw_fnc(1:gr%nz) = C7b + (C7-C7b) & 
          *EXP( -one_half * (Skw_zm(1:gr%nz)/C7c)**2 )
      else
        C7_Skw_fnc(1:gr%nz) = C7b
      endif

      ! Damp C7 as a function of Lscale in stably stratified regions
      C7_Skw_fnc = damp_coefficient( C7, C7_Skw_fnc, &
                                     C7_Lscale0, wpxp_L_thresh, Lscale )
    end if ! l_use_C7_Richardson

    ! Damp C6 as a function of Lscale in stably stratified regions
    C6rt_Skw_fnc = damp_coefficient( C6rt, C6rt_Skw_fnc, &
                                     C6rt_Lscale0, wpxp_L_thresh, Lscale )

    C6thl_Skw_fnc = damp_coefficient( C6thl, C6thl_Skw_fnc, &
                                      C6thl_Lscale0, wpxp_L_thresh, Lscale )

    !        C6rt_Skw_fnc = C6rt
    !        C6thl_Skw_fnc = C6thl
    !        C7_Skw_fnc = C7

    if ( l_stats_samp ) then

      call stat_update_var( iC7_Skw_fnc, C7_Skw_fnc, stats_zm )
      call stat_update_var( iC6rt_Skw_fnc, C6rt_Skw_fnc, stats_zm )
      call stat_update_var( iC6thl_Skw_fnc, C6thl_Skw_fnc, stats_zm )

    end if

    if ( clubb_at_least_debug_level( 0 ) ) then
      ! Assertion check for C7_Skw_fnc
      if ( any( C7_Skw_fnc(:) > one ) .or. any( C7_Skw_fnc(:) < zero ) ) then
        write(fstderr,*) "The C7_Skw_fnc variable is outside the valid range"
        err_code = clubb_fatal_error
        return
      end if
    end if

    ! Define the Coefficent of Eddy Diffusivity for the wpthlp and wprtp.
    ! Kw6 is used for wpthlp and wprtp, which are located on momentum levels.
    ! Kw6 is located on thermodynamic levels.
    ! Kw6 = c_K6 * Kh_zt

    Kw6(1:gr%nz) = c_K6 * Kh_zt(1:gr%nz)

    ! Find the number of grid levels, both upwards and downwards, that can
    ! have an effect on the central thermodynamic level during the course of
    ! one time step due to turbulent advection.  This is used as part of the
    ! monotonic turbulent advection scheme.
    call calc_turb_adv_range( dt, w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, & ! In
                              mixt_frac_zm, &  ! In
                              low_lev_effect, high_lev_effect ) ! Out

    ! Set up the implicit coefficients and explicit terms for turbulent
    ! advection of <w'rt'>, <w'thl'>, and <w'sclr'>.
    if ( l_explicit_turbulent_adv_wpxp ) then

       ! The turbulent advection of <w'x'> is handled explicitly.

       ! The <w'rt'> turbulent advection term is entirely explicit, as
       ! term_wp2rtp_explicit is equal to <w'^2 rt'> as calculated using PDF
       ! parameters, which is general for any PDF type.  The value of
       ! <w'^2 rt'> is calculated on thermodynamic levels.  The value of
       ! coef_wp2rtp_implicit is always 0.
       coef_wp2rtp_implicit = zero
       term_wp2rtp_explicit = wp2rtp

       if ( l_upwind_wpxp_ta ) then

          ! Interpolate term_wp2rtp_explicit to momentum levels as
          ! term_wp2rtp_explicit_zm.  The value of coef_wp2rtp_implicit_zm is
          ! always 0.
          coef_wp2rtp_implicit_zm = zero
          term_wp2rtp_explicit_zm = zt2zm( term_wp2rtp_explicit )

          ! Calculate the sign of the turbulent velocity for <w'rt'>.
          sgn_t_vel_wprtp &
          = sgn_turbulent_velocity( term_wp2rtp_explicit_zm, wprtp )

       endif ! l_upwind_wpxp_ta

       ! The <w'thl'> turbulent advection term is entirely explicit, as
       ! term_wp2thlp_explicit is equal to <w'^2 thl'> as calculated using PDF
       ! parameters, which is general for any PDF type.  The value of
       ! <w'^2 thl'> is calculated on thermodynamic levels.  The value of
       ! coef_wp2thlp_implicit is always 0.
       coef_wp2thlp_implicit = zero
       term_wp2thlp_explicit = wp2thlp

       if ( l_upwind_wpxp_ta ) then

          ! Interpolate term_wp2thlp_explicit to momentum levels as
          ! term_wp2thlp_explicit_zm.  The value of coef_wp2thlp_implicit_zm is
          ! always 0.
          coef_wp2thlp_implicit_zm = zero
          term_wp2thlp_explicit_zm = zt2zm( term_wp2thlp_explicit )

          ! Calculate the sign of the turbulent velocity for <w'thl'>.
          sgn_t_vel_wpthlp &
          = sgn_turbulent_velocity( term_wp2thlp_explicit_zm, wpthlp )

       endif ! l_upwind_wpxp_ta

       do i = 1, sclr_dim, 1

          ! The <w'sclr'> turbulent advection term is entirely explicit, as
          ! term_wp2sclrp_explicit is equal to <w'^2 sclr'> as calculated
          ! using PDF parameters, which is general for any PDF type.  The value
          ! of <w'^2 sclr'> is calculated on thermodynamic levels.  The value of
          ! coef_wp2sclrp_implicit is always 0.
          coef_wp2sclrp_implicit(:,i) = zero
          term_wp2sclrp_explicit(:,i) = wp2sclrp(:,i)

          if ( l_upwind_wpxp_ta ) then

             ! Interpolate term_wp2sclrp_explicit to momentum levels as
             ! term_wp2sclrp_explicit_zm.  The value of
             ! coef_wp2sclrp_implicit_zm is always 0.
             coef_wp2sclrp_implicit_zm(:,i) = zero
             term_wp2sclrp_explicit_zm(:,i) &
             = zt2zm( term_wp2sclrp_explicit(:,i) )

             ! Calculate the sign of the turbulent velocity for <w'sclr'>.
             sgn_t_vel_wpsclrp(:,i) &
             = sgn_turbulent_velocity( term_wp2sclrp_explicit_zm(:,i), &
                                       wpsclrp(:,i) )

          endif ! l_upwind_wpxp_ta

       enddo ! i = 1, sclr_dim, 1

    else ! .not. l_explicit_turbulent_adv_xpyp

       ! The turbulent advection of <w'x'> is handled implicitly or
       ! semi-implicitly.

       if ( iiPDF_type == iiPDF_ADG1 ) then

          ! The ADG1 PDF is used.

          ! Calculate the implicit coefficients and explicit terms on
          ! thermodynamic grid levels.

          ! Calculate a_1.
          ! It is a variable that is a function of sigma_sqd_w (where
          ! sigma_sqd_w is located on momentum levels).
          a1(1:gr%nz) = one / ( one - sigma_sqd_w(1:gr%nz) )

          ! Interpolate a_1 from momentum levels to thermodynamic levels.  This
          ! will be used for the <w'x'> turbulent advection (ta) term.
          a1_zt = max( zm2zt( a1 ), zero_threshold )   ! Positive def. quantity

          ! Implicit coefficient on <w'rt'> in <w'^2 rt'> equation.
          coef_wp2rtp_implicit = a1_zt * wp3_on_wp2_zt

          ! The <w'rt'> turbulent advection term is entirely implicit, as
          ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'>.  The value of
          ! term_wp2rtp_explicit is always 0.
          term_wp2rtp_explicit = zero

          if ( l_upwind_wpxp_ta ) then

             ! Calculate coef_wp2rtp_implicit on momentum levels as
             ! coef_wp2rtp_implicit_zm.  The value of term_wp2rtp_explicit_zm is
             ! always 0.
             coef_wp2rtp_implicit_zm = a1 * wp3_on_wp2
             term_wp2rtp_explicit_zm = zero

             ! For ADG1, the sign of the turbulent velocity is the sign of
             ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity
             ! is set to wp3_on_wp2.
             sgn_t_vel_wprtp = wp3_on_wp2

          endif ! l_upwind_wpxp_ta

          ! Implicit coefficient on <w'thl'> in <w'^2 thl'> equation.
          ! For ADG1, this is the same as coef_wp2rtp_implicit.
          coef_wp2thlp_implicit = coef_wp2rtp_implicit

          ! The <w'thl'> turbulent advection term is entirely implicit, as
          ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'>.  The value of
          ! term_wp2thlp_explicit is always 0.
          term_wp2thlp_explicit = zero

          if ( l_upwind_wpxp_ta ) then

             ! Calculate coef_wp2thlp_implicit on momentum levels as
             ! coef_wp2thlp_implicit_zm.  The value of term_wp2thlp_explicit_zm
             ! is always 0.
             coef_wp2thlp_implicit_zm = coef_wp2rtp_implicit_zm
             term_wp2thlp_explicit_zm = zero

             ! For ADG1, the sign of the turbulent velocity is the sign of
             ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity
             ! is set to wp3_on_wp2.
             sgn_t_vel_wpthlp = wp3_on_wp2

          endif ! l_upwind_wpxp_ta

          do i = 1, sclr_dim, 1

             ! Implicit coefficient on <w'sclr'> in <w'^2 sclr'> equation.
             ! For ADG1, this is the same as coef_wp2rtp_implicit.
             coef_wp2sclrp_implicit(:,i) = coef_wp2rtp_implicit

             ! The <w'sclr'> turbulent advection term is entirely implicit, as
             ! <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'>.  The value of
             ! term_wp2sclrp_explicit is always 0.
             term_wp2sclrp_explicit(:,i) = zero

             if ( l_upwind_wpxp_ta ) then

                ! Calculate coef_wp2sclrp_implicit on momentum levels as
                ! coef_wp2sclrp_implicit_zm.  The value of
                ! term_wp2sclrp_explicit_zm is always 0.
                coef_wp2sclrp_implicit_zm(:,i) = coef_wp2rtp_implicit_zm
                term_wp2sclrp_explicit_zm(:,i) = zero

                ! For ADG1, the sign of the turbulent velocity is the sign of
                ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
                ! velocity is set to wp3_on_wp2.
                sgn_t_vel_wpsclrp(:,i) = wp3_on_wp2

             endif ! l_upwind_wpxp_ta

          enddo ! i = 1, sclr_dim, 1

       elseif ( iiPDF_type == iiPDF_new ) then

          ! The new PDF is used.

          ! Unpack the variables coef_wp2rtp_implicit, term_wp2rtp_explicit,
          ! coef_wp2thlp_implicit, and term_wp2thlp_explicit from
          ! pdf_implicit_coefs_terms.  The PDF parameters and the resulting
          ! implicit coefficients and explicit terms are calculated on
          ! thermodynamic levels.

          ! Implicit coefficient on <w'rt'> in <w'^2 rt'> equation.
          coef_wp2rtp_implicit = pdf_implicit_coefs_terms%coef_wp2rtp_implicit

          ! Explicit (RHS) term in <w'rt'> equation.
          term_wp2rtp_explicit = pdf_implicit_coefs_terms%term_wp2rtp_explicit

          if ( l_upwind_wpxp_ta ) then

             ! Interpolate coef_wp2rtp_implicit and term_wp2rtp_explicit
             ! to momentum levels as coef_wp2rtp_implicit_zm and
             ! term_wp2rtp_explicit_zm, respectively.
             coef_wp2rtp_implicit_zm = zt2zm( coef_wp2rtp_implicit )
             term_wp2rtp_explicit_zm = zt2zm( term_wp2rtp_explicit )

             ! Calculate the sign of the turbulent velocity for <w'rt'>.
             sgn_t_vel_wprtp &
             = sgn_turbulent_velocity( coef_wp2rtp_implicit_zm * wprtp &
                                       + term_wp2rtp_explicit_zm, wprtp )

          endif ! l_upwind_wpxp_ta

          ! Implicit coefficient on <w'thl'> in <w'^2 thl'> equation.
          coef_wp2thlp_implicit &
          = pdf_implicit_coefs_terms%coef_wp2thlp_implicit

          ! Explicit (RHS) term in <w'thl'> equation.
          term_wp2thlp_explicit &
          = pdf_implicit_coefs_terms%term_wp2thlp_explicit

          if ( l_upwind_wpxp_ta ) then

             ! Interpolate coef_wp2thlp_implicit and term_wp2thlp_explicit
             ! to momentum levels as coef_wp2thlp_implicit_zm and
             ! term_wp2thlp_explicit_zm, respectively.
             coef_wp2thlp_implicit_zm = zt2zm( coef_wp2thlp_implicit )
             term_wp2thlp_explicit_zm = zt2zm( term_wp2thlp_explicit )

             ! Calculate the sign of the turbulent velocity for <w'thl'>.
             sgn_t_vel_wpthlp &
             = sgn_turbulent_velocity( coef_wp2thlp_implicit_zm * wpthlp &
                                       + term_wp2thlp_explicit_zm, wpthlp )

          endif ! l_upwind_wpxp_ta

          do i = 1, sclr_dim, 1

             ! The code for the scalar variables will be set up later.
             coef_wp2sclrp_implicit(:,i) = zero
             term_wp2sclrp_explicit(:,i) = zero

             if ( l_upwind_wpxp_ta ) then

                ! Interpolate coef_wp2sclrp_implicit and term_wp2sclrp_explicit
                ! to momentum levels as coef_wp2sclrp_implicit_zm and
                ! term_wp2sclrp_explicit_zm, respectively.
                coef_wp2sclrp_implicit_zm(:,i) &
                = zt2zm( coef_wp2sclrp_implicit(:,i) )
                term_wp2sclrp_explicit_zm(:,i) &
                = zt2zm( term_wp2sclrp_explicit(:,i) )

                ! Calculate the sign of the turbulent velocity for <w'sclr'>.
                sgn_t_vel_wpsclrp(:,i) &
                = sgn_turbulent_velocity( coef_wp2sclrp_implicit_zm(:,i) &
                                          * wpsclrp(:,i) &
                                          + term_wp2sclrp_explicit_zm(:,i), &
                                          wpsclrp(:,i) )

             endif ! l_upwind_wpxp_ta

          enddo ! i = 1, sclr_dim, 1

       endif ! iiPDF_type

    endif ! l_explicit_turbulent_adv_xpyp

    if ( l_stats_samp ) then
       call stat_update_var( icoef_wp2rtp_implicit, coef_wp2rtp_implicit, &
                             stats_zt )
       call stat_update_var( iterm_wp2rtp_explicit, term_wp2rtp_explicit, &
                             stats_zt )
       call stat_update_var( icoef_wp2thlp_implicit, coef_wp2thlp_implicit, &
                             stats_zt )
       call stat_update_var( iterm_wp2thlp_explicit, term_wp2thlp_explicit, &
                             stats_zt )
    endif

    ! Setup and decompose matrix for each variable.

    if ( l_clip_semi_implicit &
         .or. ( ( iiPDF_type == iiPDF_new ) &
                .and. ( .not. l_explicit_turbulent_adv_wpxp ) ) ) then

      ! Compute the upper and lower limits of w'r_t' at every level,
      ! based on the correlation of w and r_t, such that:
      ! corr_(w,r_t) = w'r_t' / [ sqrt(w'^2) * sqrt(r_t'^2) ];
      ! -1 <= corr_(w,r_t) <= 1.
      if ( l_clip_semi_implicit ) then
        wpxp_upper_lim =  max_mag_correlation * sqrt( wp2 * rtp2 )
        wpxp_lower_lim = -wpxp_upper_lim
      endif

      ! Compute the implicit portion of the r_t and w'r_t' equations.
      ! Build the left-hand side matrix.
      call xm_wpxp_lhs( l_iter, dt, Kh_zm, wprtp, wm_zm, wm_zt, wp2, & ! In
                        coef_wp2rtp_implicit, coef_wp2rtp_implicit_zm, & ! In
                        sgn_t_vel_wprtp, Kw6, tau_C6_zm, C7_Skw_fnc, & ! In
                        C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! In
                        invrs_rho_ds_zm, invrs_rho_ds_zt, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, &! In
                        l_implemented, em, Lscale, thlm, exner, & ! In
                        rtm, rcm, p_in_Pa, thvm, & ! In
                        lhs ) ! Out

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_rtm, l_iter, dt, rtm, wprtp, & ! In
                        rtm_forcing, wprtp_forcing, C7_Skw_fnc, & ! In
                        rtpthvp, C6rt_Skw_fnc, tau_C6_zm, & ! In
                        coef_wp2rtp_implicit, coef_wp2rtp_implicit_zm, & ! In
                        term_wp2rtp_explicit, term_wp2rtp_explicit_zm, & ! In
                        sgn_t_vel_wprtp, rho_ds_zt, & ! In
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, & ! In
                        rhs(:,1) ) ! Out

      ! Save the value of rhs, which will be overwritten with the solution as
      ! part of the solving routine.
      rhs_save = rhs

      ! Solve r_t / w'r_t'
      if ( l_stats_samp .and. irtm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                            lhs, rhs, &                 ! Intent(inout)
                            solution, rcond )           ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution )           ! Intent(out)
      endif

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Mean total water & total water flux LU decomp. failed"
            write(fstderr,*) "rtm and wprtp LHS"
            do k = 1, gr%nz
               write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k-1)
               write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k)
            enddo ! k = 1, gr%nz
            write(fstderr,*) "rtm and wprtp RHS"
            do k = 1, gr%nz
               write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                "RHS = ", rhs_save(2*k-1,1)
               write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                "RHS = ", rhs_save(2*k,1)
            enddo ! k = 1, gr%nz
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
             rtm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             rt_tol**2, rt_tol, rcond, &            ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,1), &        ! Intent(in)
             rtm, rt_tol_mfl, wprtp )               ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "rtm monotonic flux limiter:  tridag failed"
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      ! Compute the upper and lower limits of w'th_l' at every level,
      ! based on the correlation of w and th_l, such that:
      ! corr_(w,th_l) = w'th_l' / [ sqrt(w'^2) * sqrt(th_l'^2) ];
      ! -1 <= corr_(w,th_l) <= 1.
      if ( l_clip_semi_implicit ) then
        wpxp_upper_lim =  max_mag_correlation * sqrt( wp2 * thlp2 )
        wpxp_lower_lim = -wpxp_upper_lim
      endif

      ! Compute the implicit portion of the th_l and w'th_l' equations.
      ! Build the left-hand side matrix.
      call xm_wpxp_lhs( l_iter, dt, Kh_zm, wpthlp, wm_zm, wm_zt, wp2, & ! In
                        coef_wp2thlp_implicit, coef_wp2thlp_implicit_zm, & ! In
                        sgn_t_vel_wpthlp, Kw6, tau_C6_zm, C7_Skw_fnc, & ! In
                        C6thl_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! In
                        invrs_rho_ds_zm, invrs_rho_ds_zt, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, &! In
                        l_implemented, em, Lscale, thlm, exner, & ! In
                        rtm, rcm, p_in_Pa, thvm, & ! In
                        lhs ) ! Out

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_thlm, l_iter, dt, thlm, wpthlp, & ! In
                        thlm_forcing, wpthlp_forcing, C7_Skw_fnc, & ! In
                        thlpthvp, C6thl_Skw_fnc, tau_C6_zm, & ! In
                        coef_wp2thlp_implicit, coef_wp2thlp_implicit_zm, & ! In
                        term_wp2thlp_explicit, term_wp2thlp_explicit_zm, & ! In
                        sgn_t_vel_wpthlp, rho_ds_zt, & ! In
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, & ! In
                        rhs(:,1) ) ! Out

      ! Save the value of rhs, which will be overwritten with the solution as
      ! part of the solving routine.
      rhs_save = rhs

      ! Solve for th_l / w'th_l'
      if ( l_stats_samp .and. ithlm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                            lhs, rhs, &                 ! Intent(inout)
                            solution, rcond )           ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution )           ! Intent(out)
      endif

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Liquid pot. temp & thetal flux LU decomp. failed"
            write(fstderr,*) "thlm and wpthlp LHS"
            do k = 1, gr%nz
               write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k-1)
               write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k)
            enddo ! k = 1, gr%nz
            write(fstderr,*) "thlm and wpthlp RHS"
            do k = 1, gr%nz
               write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                "RHS = ", rhs_save(2*k-1,1)
               write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                "RHS = ", rhs_save(2*k,1)
            enddo ! k = 1, gr%nz
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_thlm, dt, wp2, thlp2, wm_zt,  & ! Intent(in)
             thlm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &     ! Intent(in)
             thl_tol**2, thl_tol, rcond, &           ! Intent(in)
             low_lev_effect, high_lev_effect, &      ! Intent(in)
             l_implemented, solution(:,1),  &        ! Intent(in)
             thlm, thl_tol_mfl, wpthlp )             ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "thlm monotonic flux limiter:  tridag failed" 
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      ! Solve sclrm / wpsclrp
      ! If sclr_dim is 0, then this loop will execute 0 times.
! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
      do i = 1, 0, 1
#else
      do i = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

        ! Compute the upper and lower limits of w'sclr' at every level,
        ! based on the correlation of w and sclr, such that:
        ! corr_(w,sclr) = w'sclr' / [ sqrt(w'^2) * sqrt(sclr'^2) ];
        ! -1 <= corr_(w,sclr) <= 1.
        if ( l_clip_semi_implicit ) then
          wpxp_upper_lim(:) =  max_mag_correlation * sqrt( wp2(:) * sclrp2(:,i) )
          wpxp_lower_lim(:) = -wpxp_upper_lim(:)
        endif

        ! Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
        ! using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
        wpsclrp_forcing(:,i) = zero

        ! Compute the implicit portion of the sclr and w'sclr' equations.
        ! Build the left-hand side matrix.
        call xm_wpxp_lhs( l_iter, dt, Kh_zm, wpsclrp(:,i), &  ! In
                          wm_zm, wm_zt, wp2, & ! In
                          coef_wp2sclrp_implicit(:,i), & ! In
                          coef_wp2sclrp_implicit_zm(:,i), & ! In
                          sgn_t_vel_wpsclrp(:,i), Kw6, tau_C6_zm, & ! In
                          C7_Skw_fnc, C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! In
                          invrs_rho_ds_zm, invrs_rho_ds_zt, & ! In
                          wpxp_upper_lim, wpxp_lower_lim, &! In
                          l_implemented, em, Lscale, thlm, exner, & ! In
                          rtm, rcm, p_in_Pa, thvm, & ! In
                          lhs ) ! Out

        ! Compute the explicit portion of the sclrm and w'sclr' equations.
        ! Build the right-hand side vector.
        call xm_wpxp_rhs( xm_wpxp_scalar, l_iter, dt, sclrm(:,i), & ! In
                          wpsclrp(:,i), sclrm_forcing(:,i), &
                          wpsclrp_forcing(:,i), C7_Skw_fnc, & ! In
                          sclrpthvp(:,i), C6rt_Skw_fnc, tau_C6_zm, & ! In
                          coef_wp2sclrp_implicit(:,i), & ! In
                          coef_wp2sclrp_implicit_zm(:,i), & ! In
                          term_wp2sclrp_explicit(:,i), & ! In
                          term_wp2sclrp_explicit_zm(:,i), & ! In
                          sgn_t_vel_wpsclrp(:,i), rho_ds_zt, & ! In
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                          wpxp_upper_lim, wpxp_lower_lim, & ! In
                          rhs(:,1) ) ! Out

        ! Save the value of rhs, which will be overwritten with the solution as
        ! part of the solving routine.
        rhs_save = rhs

        ! Solve for sclrm / w'sclr'
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution )           ! Intent(out)

        if ( clubb_at_least_debug_level( 0 ) ) then
           if ( err_code == clubb_fatal_error ) then   
              write(fstderr,*) "Passive scalar # ", i, " LU decomp. failed."
              write(fstderr,*) "sclrm and wpsclrp LHS"
              do k = 1, gr%nz
                 write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                  "LHS = ", lhs(1:nsup+nsub+1,2*k-1)
                 write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                  "LHS = ", lhs(1:nsup+nsub+1,2*k)
              enddo ! k = 1, gr%nz
              write(fstderr,*) "sclrm and wpsclrp RHS"
              do k = 1, gr%nz
                 write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                  "RHS = ", rhs_save(2*k-1,1)
                 write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                  "RHS = ", rhs_save(2*k,1)
              enddo ! k = 1, gr%nz
              call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                         Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                         Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                         wp2rtp, rtpthvp, rtm_forcing, &
                                         wprtp_forcing, rtm_ref, wp2thlp, &
                                         thlpthvp, thlm_forcing, &
                                         wpthlp_forcing, thlm_ref, rho_ds_zm, &
                                         rho_ds_zt, invrs_rho_ds_zm, &
                                         invrs_rho_ds_zt, thv_ds_zm, rtp2, &
                                         thlp2, w_1_zm, w_2_zm, varnce_w_1_zm, &
                                         varnce_w_2_zm, mixt_frac_zm, &
                                         l_implemented, em, wp2sclrp, &
                                         sclrpthvp, sclrm_forcing, sclrp2, &
                                         exner, rcm, p_in_Pa, thvm, &
                                         Cx_fnc_Richardson, &
                                         pdf_implicit_coefs_terms, um_forcing, &
                                         vm_forcing, ug, vg, wpthvp, fcor, &
                                         um_ref, vm_ref, up2, vp2, uprcp, &
                                         vprcp, rc_coef, rtm, wprtp, thlm, &
                                         wpthlp, sclrm, wpsclrp, um, upwp, vm, &
                                         vpwp, rtm_old, wprtp_old, thlm_old, &
                                         wpthlp_old, sclrm_old, wpsclrp_old, &
                                         um_old, upwp_old, vm_old, vpwp_old )
              return
           endif
        endif

        call xm_wpxp_clipping_and_stats &
             ( xm_wpxp_scalar, dt, wp2, sclrp2(:,i),  & ! Intent(in)
               wm_zt, sclrm_forcing(:,i),  &            ! Intent(in)
               rho_ds_zm, rho_ds_zt, &                  ! Intent(in)
               invrs_rho_ds_zm, invrs_rho_ds_zt, &      ! Intent(in)
               sclr_tol(i)**2, sclr_tol(i), rcond, &    ! Intent(in)
               low_lev_effect, high_lev_effect, &       ! Intent(in)
               l_implemented, solution(:,1),  &         ! Intent(in)
               sclrm(:,i), sclr_tol(i), wpsclrp(:,i) )  ! Intent(inout)

        if ( clubb_at_least_debug_level( 0 ) ) then
           if ( err_code == clubb_fatal_error ) then
              write(fstderr,*) "sclrm # ", i, "monotonic flux limiter: tridag failed"
              call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                         Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                         Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                         wp2rtp, rtpthvp, rtm_forcing, &
                                         wprtp_forcing, rtm_ref, wp2thlp, &
                                         thlpthvp, thlm_forcing, &
                                         wpthlp_forcing, thlm_ref, rho_ds_zm, &
                                         rho_ds_zt, invrs_rho_ds_zm, &
                                         invrs_rho_ds_zt, thv_ds_zm, rtp2, &
                                         thlp2, w_1_zm, w_2_zm, varnce_w_1_zm, &
                                         varnce_w_2_zm, mixt_frac_zm, &
                                         l_implemented, em, wp2sclrp, &
                                         sclrpthvp, sclrm_forcing, sclrp2, &
                                         exner, rcm, p_in_Pa, thvm, &
                                         Cx_fnc_Richardson, &
                                         pdf_implicit_coefs_terms, um_forcing, &
                                         vm_forcing, ug, vg, wpthvp, fcor, &
                                         um_ref, vm_ref, up2, vp2, uprcp, &
                                         vprcp, rc_coef, rtm, wprtp, thlm, &
                                         wpthlp, sclrm, wpsclrp, um, upwp, vm, &
                                         vpwp, rtm_old, wprtp_old, thlm_old, &
                                         wpthlp_old, sclrm_old, wpsclrp_old, &
                                         um_old, upwp_old, vm_old, vpwp_old )
              return
           endif
        endif

      enddo ! passive scalars

    else

      ! Simple case, where l_clip_semi_implicit is false, and if the new PDF is
      ! used, l_explicit_turbulent_adv_wpxp is enabled.

      ! Create the lhs once
      call xm_wpxp_lhs( l_iter, dt, Kh_zm, zeros_vector, wm_zm, wm_zt, wp2, & ! In
                        coef_wp2rtp_implicit, coef_wp2rtp_implicit_zm, & ! In
                        sgn_t_vel_wprtp, Kw6, tau_C6_zm, C7_Skw_fnc, & ! In
                        C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! In
                        invrs_rho_ds_zm, invrs_rho_ds_zt, & ! In
                        zeros_vector, zeros_vector, &! In
                        l_implemented, em, Lscale, thlm, exner, & ! In
                        rtm, rcm, p_in_Pa, thvm, & ! In
                        lhs ) ! Out

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_rtm, l_iter, dt, rtm, wprtp, & ! In
                        rtm_forcing, wprtp_forcing, C7_Skw_fnc, & ! In
                        rtpthvp, C6rt_Skw_fnc, tau_C6_zm, & ! In
                        coef_wp2rtp_implicit, coef_wp2rtp_implicit_zm, & ! In
                        term_wp2rtp_explicit, term_wp2rtp_explicit_zm, & ! In
                        sgn_t_vel_wprtp, rho_ds_zt, & ! In
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, & ! In
                        rhs(:,1) ) ! Out

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_thlm, l_iter, dt, thlm, wpthlp, & ! In
                        thlm_forcing, wpthlp_forcing, C7_Skw_fnc, & ! In
                        thlpthvp, C6thl_Skw_fnc, tau_C6_zm, & ! In
                        coef_wp2thlp_implicit, coef_wp2thlp_implicit_zm, & ! In
                        term_wp2thlp_explicit, term_wp2thlp_explicit_zm, & ! In
                        sgn_t_vel_wpthlp, rho_ds_zt, & ! In
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                        wpxp_upper_lim, wpxp_lower_lim, & ! In
                        rhs(:,2) ) ! Out

! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
      do i = 1, 0, 1
#else
      do i = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

        ! Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
        ! using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
        wpsclrp_forcing(:,i) = zero

        call xm_wpxp_rhs( xm_wpxp_scalar, l_iter, dt, sclrm(:,i), & ! In
                          wpsclrp(:,i), sclrm_forcing(:,i), &
                          wpsclrp_forcing(:,i), C7_Skw_fnc, & ! In
                          sclrpthvp(:,i), C6rt_Skw_fnc, tau_C6_zm, & ! In
                          coef_wp2sclrp_implicit(:,i), & ! In
                          coef_wp2sclrp_implicit_zm(:,i), & ! In
                          term_wp2sclrp_explicit(:,i), & ! In
                          term_wp2sclrp_explicit_zm(:,i), & ! In
                          sgn_t_vel_wpsclrp(:,i), rho_ds_zt, & ! In
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                          wpxp_upper_lim, wpxp_lower_lim, & ! In
                          rhs(:,2+i) ) ! Out

      enddo

      if ( l_predict_upwp_vpwp ) then

         ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.

         ! Currently, this requires the ADG1 PDF with implicit turbulent
         ! advection.
         if ( ( iiPDF_type == iiPDF_ADG1 ) &
              .and. ( .not. l_explicit_turbulent_adv_wpxp ) ) then

            ! Calculate the implicit coefficients and explicit terms on
            ! thermodynamic grid levels.

            ! Implicit coefficient on <u'w'> in <w'^2 u'> equation.
            ! For ADG1, this is the same as coef_wp2rtp_implicit.
            coef_wp2up_implicit = coef_wp2rtp_implicit

            ! The <u'w'> turbulent advection term is entirely implicit, as
            ! <w'^2 u'> = coef_wp2up_implicit * <u'w'>.  The value of
            ! term_wp2up_explicit is always 0.
            term_wp2up_explicit = zero

            if ( l_upwind_wpxp_ta ) then

               ! Calculate coef_wp2up_implicit on momentum levels as
               ! coef_wp2up_implicit_zm.  The value of term_wp2up_explicit_zm is
               ! always 0.
               coef_wp2up_implicit_zm = coef_wp2rtp_implicit_zm
               term_wp2up_explicit_zm = zero

               ! For ADG1, the sign of the turbulent velocity is the sign of
               ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
               ! velocity is set to wp3_on_wp2.
               sgn_t_vel_upwp = wp3_on_wp2

            endif ! l_upwind_wpxp_ta

            ! Implicit coefficient on <v'w'> in <w'^2 v'> equation.
            ! For ADG1, this is the same as coef_wp2rtp_implicit.
            coef_wp2vp_implicit = coef_wp2rtp_implicit

            ! The <v'w'> turbulent advection term is entirely implicit, as
            ! <w'^2 v'> = coef_wp2vp_implicit * <v'w'>.  The value of
            ! term_wp2vp_explicit is always 0.
            term_wp2vp_explicit = zero

            if ( l_upwind_wpxp_ta ) then

               ! Calculate coef_wp2vp_implicit on momentum levels as
               ! coef_wp2vp_implicit_zm.  The value of term_wp2vp_explicit_zm is
               ! always 0.
               coef_wp2vp_implicit_zm = coef_wp2rtp_implicit_zm
               term_wp2vp_explicit_zm = zero

               ! For ADG1, the sign of the turbulent velocity is the sign of
               ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
               ! velocity is set to wp3_on_wp2.
               sgn_t_vel_vpwp = wp3_on_wp2

            endif ! l_upwind_wpxp_ta

         endif ! ( iiPDF_type == iiPDF_ADG1 )
               ! .and. ( .not. l_explicit_turbulent_adv_wpxp )

         ! Coriolis term for <u> and <v>
         if ( .not. l_implemented ) then

            ! Only compute the Coriolis term if the model is running on its own,
            ! and is not part of a larger, host model.
            um_tndcy = um_forcing - fcor * ( vg - vm )
            vm_tndcy = vm_forcing + fcor * ( ug - um )

            if ( l_stats_samp ) then

               ! um or vm term gf is completely explicit; call stat_update_var.
               call stat_update_var( ium_gf, - fcor * vg, stats_zt )
               call stat_update_var( ivm_gf, fcor * ug, stats_zt )

               ! um or vm term cf is completely explicit; call stat_update_var.
               call stat_update_var( ium_cf, fcor * vm, stats_zt )
               call stat_update_var( ivm_cf, - fcor * um, stats_zt )

               ! um or vm forcing term
               call stat_update_var( ium_f, um_forcing, stats_zt )
               call stat_update_var( ivm_f, vm_forcing, stats_zt )

            endif ! l_stats_samp

         else ! implemented in a host model

            um_tndcy = zero
            vm_tndcy = zero

         endif ! .not. l_implemented

         ! Add "extra term" and optional Coriolis term for <u'w'> and <v'w'>.
         upwp_forcing = C7_Skw_fnc * wp2 * ddzt( um )
         vpwp_forcing = C7_Skw_fnc * wp2 * ddzt( vm )

         if ( l_stats_samp ) then
            call stat_update_var( iupwp_pr4, C7_Skw_fnc * wp2 * ddzt( um ), &
                                  stats_zm )
            call stat_update_var( ivpwp_pr4, C7_Skw_fnc * wp2 * ddzt( vm ), &
                                  stats_zm )
         endif ! l_stats_samp

         call diagnose_upxp( upwp, thlm, wpthlp, um, &               ! Intent(in)
                             C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc, & ! Intent(in)
                             upthlp )                                ! Intent(out)
         call diagnose_upxp( upwp, rtm, wprtp, um, &                ! Intent(in)
                             C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc, & ! Intent(in)
                             uprtp )                                ! Intent(out)
         call diagnose_upxp( vpwp, thlm, wpthlp, vm, &               ! Intent(in)
                             C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc, & ! Intent(in)
                             vpthlp )                                ! Intent(out)
         call diagnose_upxp( vpwp, rtm, wprtp, vm, &                ! Intent(in)
                             C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc, & ! Intent(in)
                             vprtp )                                ! Intent(out)

         ! Use a crude approximation for buoyancy terms <u'thv'> and <v'thv'>.
         !upthvp = upwp * wpthvp / max( wp2, w_tol_sqd )
         !vpthvp = vpwp * wpthvp / max( wp2, w_tol_sqd )
         !upthvp = 0.3_core_rknd * ( upthlp + 200.0_core_rknd * uprtp ) &
         !         + 200._core_rknd * sign( one, upwp) * sqrt( up2 * rcm**2 )
         !vpthvp = 0.3_core_rknd * ( vpthlp + 200.0_core_rknd * vprtp ) &
         !         + 200._core_rknd * sign( one, vpwp ) * sqrt( vp2 * rcm**2 )
         upthvp = upthlp + ep1 * thv_ds_zm * uprtp + rc_coef * uprcp
         vpthvp = vpthlp + ep1 * thv_ds_zm * vprtp + rc_coef * vprcp

         if ( l_stats_samp ) then
            call stat_update_var( iupthlp, upthlp, stats_zm )
            call stat_update_var( iuprtp,  uprtp,  stats_zm )
            call stat_update_var( ivpthlp, vpthlp, stats_zm )
            call stat_update_var( ivprtp,  vprtp,  stats_zm )
            call stat_update_var( iupthvp, upthvp, stats_zm )
            call stat_update_var( iuprcp,  uprcp,  stats_zm )
            call stat_update_var( ivpthvp, vpthvp, stats_zm )
            call stat_update_var( ivprcp,  vprcp,  stats_zm )
         endif ! l_stats_samp

         call xm_wpxp_rhs( xm_wpxp_um, l_iter, dt, um, upwp, & ! In
                           um_tndcy, upwp_forcing, C7_Skw_fnc, & ! In
                           upthvp, C6rt_Skw_fnc, tau_C6_zm, & ! In
                           coef_wp2up_implicit, coef_wp2up_implicit_zm, & ! In
                           term_wp2up_explicit, term_wp2up_explicit_zm, & ! In
                           sgn_t_vel_upwp, rho_ds_zt, & ! In
                           rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                           wpxp_upper_lim, wpxp_lower_lim, & ! In
                           rhs(:,3+sclr_dim) ) ! Out

         call xm_wpxp_rhs( xm_wpxp_vm, l_iter, dt, vm, vpwp, & ! In
                           vm_tndcy, vpwp_forcing, C7_Skw_fnc, & ! In
                           vpthvp, C6rt_Skw_fnc, tau_C6_zm, & ! In
                           coef_wp2vp_implicit, coef_wp2vp_implicit_zm, & ! In
                           term_wp2vp_explicit, term_wp2vp_explicit_zm, & ! In
                           sgn_t_vel_vpwp, rho_ds_zt, & ! In
                           rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                           wpxp_upper_lim, wpxp_lower_lim, & ! In
                           rhs(:,4+sclr_dim) ) ! Out

      endif ! l_predict_upwp_vpwp

      ! Save the value of rhs, which will be overwritten with the solution as
      ! part of the solving routine.
      rhs_save = rhs

      ! Solve for all fields
      if ( l_stats_samp &
           .and. ithlm_matrix_condt_num + irtm_matrix_condt_num > 0 ) then
         call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                             lhs, rhs, &                 ! Intent(inout)
                             solution, rcond )           ! Intent(out)
      else
         call xm_wpxp_solve( nrhs, &              ! Intent(in)
                             lhs, rhs, &          ! Intent(inout)
                             solution )           ! Intent(out)
      endif

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "xm & wpxp LU decomp. failed"
            write(fstderr,*) "General xm and wpxp LHS"
            do k = 1, gr%nz
               write(fstderr,*) "zt level = ", k, "height [m] = ", gr%zt(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k-1)
               write(fstderr,*) "zm level = ", k, "height [m] = ", gr%zm(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,2*k)
            enddo ! k = 1, gr%nz
            do i = 1, nrhs
               if ( i == 1 ) then
                  write(fstderr,*) "rtm and wprtp RHS"
               elseif ( i == 2 ) then
                  write(fstderr,*) "thlm and wpthlp RHS"
               else ! i > 2
                  if ( sclr_dim > 0 ) then
                     if ( i <= 2+sclr_dim ) then
                        write(fstderr,*) "sclrm and wpsclrp RHS for sclr", i-2
                     endif ! i <= 2+sclr_dim )
                  endif ! sclr_dim > 0
                  if ( l_predict_upwp_vpwp ) then
                     if ( i == 3+sclr_dim ) then
                        write(fstderr,*) "um and upwp RHS"
                     elseif ( i == 4+sclr_dim ) then
                        write(fstderr,*) "vm and vpwp RHS"
                     endif
                  endif ! l_predict_upwp_vpwp
               endif
               do k = 1, gr%nz
                  write(fstderr,*) "zt level = ", k, &
                                   "height [m] = ", gr%zt(k), &
                                   "RHS = ", rhs_save(2*k-1,i)
                  write(fstderr,*) "zm level = ", k, &
                                   "height [m] = ", gr%zm(k), &
                                   "RHS = ", rhs_save(2*k,i)
               enddo ! k = 1, gr%nz
            enddo ! i = 1, nrhs
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
             rtm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             rt_tol**2, rt_tol, rcond, &            ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,1),  &       ! Intent(in)
             rtm, rt_tol_mfl, wprtp )               ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "rtm monotonic flux limiter:  tridag failed"
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_thlm, dt, wp2, thlp2, wm_zt, & ! Intent(in)
             thlm_forcing, rho_ds_zm, rho_ds_zt, &  ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             thl_tol**2, thl_tol, rcond, &          ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,2),  &       ! Intent(in)
             thlm, thl_tol_mfl, wpthlp )            ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "thlm monotonic flux limiter:  tridag failed"
            call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                       Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                       Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                       wp2rtp, rtpthvp, rtm_forcing, &
                                       wprtp_forcing, rtm_ref, wp2thlp, &
                                       thlpthvp, thlm_forcing, wpthlp_forcing, &
                                       thlm_ref, rho_ds_zm, rho_ds_zt, &
                                       invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                       thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                       varnce_w_1_zm, varnce_w_2_zm, &
                                       mixt_frac_zm, l_implemented, em, &
                                       wp2sclrp, sclrpthvp, sclrm_forcing, &
                                       sclrp2, exner, rcm, p_in_Pa, thvm, &
                                       Cx_fnc_Richardson, &
                                       pdf_implicit_coefs_terms, um_forcing, &
                                       vm_forcing, ug, vg, wpthvp, fcor, &
                                       um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                       rc_coef, rtm, wprtp, thlm, wpthlp, &
                                       sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                       rtm_old, wprtp_old, thlm_old, &
                                       wpthlp_old, sclrm_old, wpsclrp_old, &
                                       um_old, upwp_old, vm_old, vpwp_old )
            return
         endif
      endif

! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
      do i = 1, 0, 1
#else
      do i = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

        call xm_wpxp_clipping_and_stats &
             ( xm_wpxp_scalar, dt, wp2, sclrp2(:,i), &  ! Intent(in)
               wm_zt, sclrm_forcing(:,i), &             ! Intent(in)
               rho_ds_zm, rho_ds_zt, &                  ! Intent(in)
               invrs_rho_ds_zm, invrs_rho_ds_zt, &      ! Intent(in)
               sclr_tol(i)**2, sclr_tol(i), rcond, &    ! Intent(in)
               low_lev_effect, high_lev_effect, &       ! Intent(in)
               l_implemented, solution(:,2+i),  &       ! Intent(in)
               sclrm(:,i), sclr_tol(i), wpsclrp(:,i) )  ! Intent(inout)

        if ( clubb_at_least_debug_level( 0 ) ) then
           if ( err_code == clubb_fatal_error ) then
              write(fstderr,*) "sclrm # ", i, "monotonic flux limiter: tridag failed"
              call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                         Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                         Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                         wp2rtp, rtpthvp, rtm_forcing, &
                                         wprtp_forcing, rtm_ref, wp2thlp, &
                                         thlpthvp, thlm_forcing, &
                                         wpthlp_forcing, thlm_ref, rho_ds_zm, &
                                         rho_ds_zt, invrs_rho_ds_zm, &
                                         invrs_rho_ds_zt, thv_ds_zm, rtp2, &
                                         thlp2, w_1_zm, w_2_zm, varnce_w_1_zm, &
                                         varnce_w_2_zm, mixt_frac_zm, &
                                         l_implemented, em, wp2sclrp, &
                                         sclrpthvp, sclrm_forcing, sclrp2, &
                                         exner, rcm, p_in_Pa, thvm, &
                                         Cx_fnc_Richardson, &
                                         pdf_implicit_coefs_terms, um_forcing, &
                                         vm_forcing, ug, vg, wpthvp, fcor, &
                                         um_ref, vm_ref, up2, vp2, uprcp, &
                                         vprcp, rc_coef, rtm, wprtp, thlm, &
                                         wpthlp, sclrm, wpsclrp, um, upwp, vm, &
                                         vpwp, rtm_old, wprtp_old, thlm_old, &
                                         wpthlp_old, sclrm_old, wpsclrp_old, &
                                         um_old, upwp_old, vm_old, vpwp_old )
              return
           endif
        endif

      end do ! 1..sclr_dim

      if ( l_predict_upwp_vpwp ) then

         ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.

         call xm_wpxp_clipping_and_stats &
              ( xm_wpxp_um, dt, wp2, up2, wm_zt,       & ! Intent(in)
                um_tndcy, rho_ds_zm, rho_ds_zt,        & ! Intent(in)
                invrs_rho_ds_zm, invrs_rho_ds_zt,      & ! Intent(in)
                w_tol_sqd, w_tol, rcond,               & ! Intent(in)
                low_lev_effect, high_lev_effect,       & ! Intent(in)
                l_implemented, solution(:,3+sclr_dim), & ! Intent(in)
                um, w_tol, upwp                        ) ! Intent(inout)

         if ( clubb_at_least_debug_level( 0 ) ) then
            if ( err_code == clubb_fatal_error ) then
               write(fstderr,*) "um monotonic flux limiter:  tridag failed"
               call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                          Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                          Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                          wp2rtp, rtpthvp, rtm_forcing, &
                                          wprtp_forcing, rtm_ref, wp2thlp, &
                                          thlpthvp, thlm_forcing, &
                                          wpthlp_forcing, thlm_ref, rho_ds_zm, &
                                          rho_ds_zt, invrs_rho_ds_zm, &
                                          invrs_rho_ds_zt, thv_ds_zm, rtp2, &
                                          thlp2, w_1_zm, w_2_zm, &
                                          varnce_w_1_zm, varnce_w_2_zm, &
                                          mixt_frac_zm, l_implemented, em, &
                                          wp2sclrp, sclrpthvp, sclrm_forcing, &
                                          sclrp2, exner, rcm, p_in_Pa, thvm, &
                                          Cx_fnc_Richardson, &
                                          pdf_implicit_coefs_terms, &
                                          um_forcing, vm_forcing, ug, vg, &
                                          wpthvp, fcor, um_ref, vm_ref, up2, &
                                          vp2, uprcp, vprcp, rc_coef, rtm, &
                                          wprtp, thlm, wpthlp, sclrm, wpsclrp, &
                                          um, upwp, vm, vpwp, rtm_old, &
                                          wprtp_old, thlm_old, wpthlp_old, &
                                          sclrm_old, wpsclrp_old, um_old, &
                                          upwp_old, vm_old, vpwp_old )
               return
            endif
         endif

         call xm_wpxp_clipping_and_stats &
              ( xm_wpxp_vm, dt, wp2, vp2, wm_zt,       & ! Intent(in)
                vm_tndcy, rho_ds_zm, rho_ds_zt,        & ! Intent(in)
                invrs_rho_ds_zm, invrs_rho_ds_zt,      & ! Intent(in)
                w_tol_sqd, w_tol, rcond,               & ! Intent(in)
                low_lev_effect, high_lev_effect,       & ! Intent(in)
                l_implemented, solution(:,4+sclr_dim), & ! Intent(in)
                vm, w_tol, vpwp                        ) ! Intent(inout)

         if ( clubb_at_least_debug_level( 0 ) ) then
            if ( err_code == clubb_fatal_error ) then
               write(fstderr,*) "vm monotonic flux limiter:  tridag failed"
               call error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                          Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                          Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                          wp2rtp, rtpthvp, rtm_forcing, &
                                          wprtp_forcing, rtm_ref, wp2thlp, &
                                          thlpthvp, thlm_forcing, &
                                          wpthlp_forcing, thlm_ref, rho_ds_zm, &
                                          rho_ds_zt, invrs_rho_ds_zm, &
                                          invrs_rho_ds_zt, thv_ds_zm, rtp2, &
                                          thlp2, w_1_zm, w_2_zm, &
                                          varnce_w_1_zm, varnce_w_2_zm, &
                                          mixt_frac_zm, l_implemented, em, &
                                          wp2sclrp, sclrpthvp, sclrm_forcing, &
                                          sclrp2, exner, rcm, p_in_Pa, thvm, &
                                          Cx_fnc_Richardson, &
                                          pdf_implicit_coefs_terms, &
                                          um_forcing, vm_forcing, ug, vg, &
                                          wpthvp, fcor, um_ref, vm_ref, up2, &
                                          vp2, uprcp, vprcp, rc_coef, rtm, &
                                          wprtp, thlm, wpthlp, sclrm, wpsclrp, &
                                          um, upwp, vm, vpwp, rtm_old, &
                                          wprtp_old, thlm_old, wpthlp_old, &
                                          sclrm_old, wpsclrp_old, um_old, &
                                          upwp_old, vm_old, vpwp_old )
               return
            endif
         endif

      endif ! l_predict_upwp_vpwp

    endif ! l_clip_semi_implicit &
          ! .or. ( ( iiPDF_type == iiPDF_new ) &
          !        .and. ( .not. l_explicit_turbulent_adv_wpxp ) )

    ! De-allocate memory
    deallocate( rhs, rhs_save, solution )

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( irtm_sdmp, rtm / dt, stats_zt )
       endif

       rtm(1:gr%nz) = sponge_damp_xm( dt, gr%zt, rtm_ref(1:gr%nz), &
                                      rtm(1:gr%nz), rtm_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( irtm_sdmp, rtm / dt, stats_zt )
       endif

    endif ! rtm_sponge_damp_settings%l_sponge_damping

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( ithlm_sdmp, thlm / dt, stats_zt )
       endif

       thlm(1:gr%nz) = sponge_damp_xm( dt, gr%zt, thlm_ref(1:gr%nz), &
                                       thlm(1:gr%nz), thlm_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( ithlm_sdmp, thlm / dt, stats_zt )
       endif

    endif ! thlm_sponge_damp_settings%l_sponge_damping

    if ( l_predict_upwp_vpwp ) then

       if ( uv_sponge_damp_settings%l_sponge_damping ) then

          if ( l_stats_samp ) then
             call stat_begin_update( ium_sdmp, um / dt, stats_zt )
             call stat_begin_update( ivm_sdmp, vm / dt, stats_zt )
          endif

          um(1:gr%nz) = sponge_damp_xm( dt, gr%zt, um_ref(1:gr%nz), &
                                        um(1:gr%nz), uv_sponge_damp_profile )

          vm(1:gr%nz) = sponge_damp_xm( dt, gr%zt, vm_ref(1:gr%nz), &
                                        vm(1:gr%nz), uv_sponge_damp_profile )

          if ( l_stats_samp ) then
             call stat_end_update( ium_sdmp, um / dt, stats_zt )
             call stat_end_update( ivm_sdmp, vm / dt, stats_zt )
          endif

       endif ! uv_sponge_damp_settings%l_sponge_damping

       ! Adjust um and vm if nudging is turned on.
       if ( l_uv_nudge ) then

          ! Reflect nudging in budget
          if ( l_stats_samp ) then
             call stat_begin_update( ium_ndg, um / dt, stats_zt ) 
             call stat_begin_update( ivm_ndg, vm / dt, stats_zt )
          endif
      
          um(1:gr%nz) &
          = um(1:gr%nz) - ( ( um(1:gr%nz) - um_ref(1:gr%nz) ) * (dt/ts_nudge) )
          vm(1:gr%nz) &
          = vm(1:gr%nz) - ( ( vm(1:gr%nz) - vm_ref(1:gr%nz) ) * (dt/ts_nudge) )

          ! Reflect nudging in budget
          if ( l_stats_samp ) then
             call stat_end_update( ium_ndg, um / dt, stats_zt )
             call stat_end_update( ivm_ndg, vm / dt, stats_zt )
          endif

       endif ! l_uv_nudge

       if ( l_stats_samp ) then
          call stat_update_var( ium_ref, um_ref, stats_zt )
          call stat_update_var( ivm_ref, vm_ref, stats_zt )
       endif

    endif ! l_predict_upwp_vpwp


    return

  end subroutine advance_xm_wpxp

  !=============================================================================
  subroutine xm_wpxp_lhs( l_iter, dt, Kh_zm, wpxp, wm_zm, wm_zt, wp2, & ! In
                          coef_wp2xp_implicit, coef_wp2xp_implicit_zm, & ! In
                          sgn_turbulent_vel, Kw6, tau_C6_zm, C7_Skw_fnc, & ! In
                          C6x_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! In
                          invrs_rho_ds_zm, invrs_rho_ds_zt, & ! In
                          wpxp_upper_lim, wpxp_lower_lim, &! In
                          l_implemented, em, Lscale, thlm, exner, & ! In
                          rtm, rcm, p_in_Pa, thvm, & ! In
                          lhs ) ! Out

    ! Description:
    !   Compute LHS band diagonal matrix for xm and w'x'.
    !   This subroutine computes the implicit portion of
    !   the xm and w'x' equations.
    ! 
    ! 
    ! Notes: 
    ! 
    !   Boundary conditions:
    !       The turbulent flux (wpxp) use fixed-point boundary conditions at both the
    !       upper and lower boundaries.  Therefore, anything set in the wpxp loop
    !       at both the upper and lower boundaries would be overwritten here.
    !       However, the wpxp loop does not extend to the boundary levels.  An array
    !       with a value of 1 at the main diagonal on the left-hand side and with
    !       values of 0 at all other diagonals on the left-hand side will preserve the
    !       right-hand side value at that level.  The value of xm at level k = 1,
    !       which is below the model surface, is preserved and then overwritten to
    !       match the new value of xm at level k = 2.
    ! 
    !           xm(1)  wpxp(1) ... wpxp(nzmax)
    !         [  0.0     0.0         0.0    ]
    !         [  0.0     0.0         0.0    ]
    !         [  1.0     1.0   ...   1.0    ]
    !         [  0.0     0.0         0.0    ]
    !         [  0.0     0.0         0.0    ]
    ! 
    ! 
    !   LHS turbulent advection (ta) term:
    !        An "over-implicit" weighted time step is applied to this term.
    !        The weight of the implicit portion of this term is controlled by
    !        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
    !        equation in order to balance a weight that is not equal to 1,
    !        such that:
    !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
    !        where X is the variable that is being solved for in a predictive
    !        equation (<w'x'> in this case), y(t) is the linearized portion of
    !        the term that gets treated implicitly, and RHS is the portion of
    !        the term that is always treated explicitly.  A weight of greater
    !        than 1 can be applied to make the term more numerically stable.
    ! 
    ! 
    !   xm: Left-hand side (implicit xm portion of the code).
    ! 
    !   Thermodynamic subdiagonal (lhs index: t_km1_tdiag)
    !         [ x xm(k-1,<t+1>) ]
    !   Momentum subdiagonal (lhs index: t_km1_mdiag)
    !         [ x wpxp(k-1,<t+1>) ]
    !   Thermodynamic main diagonal (lhs index: t_k_tdiag)
    !         [ x xm(k,<t+1>) ]
    !   Momentum superdiagonal (lhs index: t_k_mdiag)
    !         [ x wpxp(k,<t+1>) ]
    !   Thermodynamic superdiagonal (lhs index: t_kp1_tdiag)
    !         [ x xm(k+1,<t+1>) ]
    ! 
    ! 
    !   w'x': Left-hand side (implicit w'x' portion of the code).
    ! 
    !   Momentum subdiagonal (lhs index: m_km1_mdiag)
    !         [ x wpxp(k-1,<t+1>) ]
    !   Thermodynamic subdiagonal (lhs index: m_k_tdiag)
    !         [ x xm(k,<t+1>) ]
    !   Momentum main diagonal (lhs index: m_k_mdiag)
    !         [ x wpxp(k,<t+1>) ]
    !   Thermodynamic superdiagonal (lhs index: m_kp1_tdiag)
    !         [ x xm(k+1,<t+1>) ]
    !   Momentum superdiagonal (lhs index: m_kp1_mdiag)
    !         [ x wpxp(k+1,<t+1>) ]
    !  
    ! 
    !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
    !   Significant changes to this routine may adversely affect computational speed
    !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
    !----------------------------------------------------------------------------------

    use parameters_tunable, only:  & 
        nu6_vert_res_dep ! Variable(s)

    use grid_class, only:  & 
        gr,  & ! Variable(s)
        zm2zt, & ! Procedure(s)
        ddzt

    use constants_clubb, only: &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_clip_semi_implicit,          & ! Variable(s)
        l_upwind_wpxp_ta,              &
        l_diffuse_rtm_and_thlm,        &
        l_stability_correct_Kh_N2_zm

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use mean_adv, only: & 
        term_ma_zt_lhs, & ! Procedure(s)
        term_ma_zt_lhs_all, &
        term_ma_zm_lhs, &
        term_ma_zm_lhs_all

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs, & ! Procedure(s)
        xpyp_term_ta_pdf_lhs_all

    use diffusion, only:  & 
        diffusion_zt_lhs, &! Procedure(s)
        diffusion_zt_lhs_all, &
        diffusion_zm_lhs, &
        diffusion_zm_lhs_all

    use clip_semi_implicit, only: & 
        clip_semi_imp_lhs ! Procedure(s)

    use stats_variables, only: & 
        ztscr01, & ! Variable(s)
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        zmscr01, & 
        zmscr02, & 
        zmscr03, & 
        zmscr04, & 
        zmscr05, & 
        zmscr06, & 
        zmscr07, & 
        zmscr08, & 
        zmscr09, & 
        zmscr10, & 
        zmscr11, & 
        zmscr12, & 
        zmscr13, & 
        zmscr14, & 
        zmscr15

    use stats_variables, only: & 
        l_stats_samp, & 
        ithlm_ma, & 
        ithlm_ta, & 
        irtm_ma, & 
        irtm_ta, & 
        iwpthlp_ma, & 
        iwpthlp_ta, & 
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_dp1, & 
        iwpthlp_sicl, & 
        iwprtp_ma, & 
        iwprtp_ta, & 
        iwprtp_tp, & 
        iwprtp_ac, & 
        iwprtp_pr1, & 
        iwprtp_pr2, & 
        iwprtp_dp1, & 
        iwprtp_sicl

    use advance_helper_module, only: &
      set_boundary_conditions_lhs, & ! Procedure(s)
      calc_stability_correction

    implicit none

    ! External
    intrinsic :: min, max

    !------------------- Input Variables -------------------
    logical, intent(in) :: l_iter

    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Timestep                                  [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      wpxp,                   & ! w'x' (momentum levs) at timestep (t) [un vary]
      Kh_zm,                  & ! Eddy diffusivity on momentum levels    [m^2/s]
      Lscale,                 & ! Turbulent mixing length                    [m]
      em,                     & ! Turbulent Kinetic Energy (TKE)       [m^2/s^2]
      thlm,                   & ! th_l (thermo. levels)                      [K]
      exner,                  & ! Exner function                             [-]
      rtm,                    & ! total water mixing ratio, r_t              [-]
      rcm,                    & ! cloud water mixing ratio, r_c          [kg/kg]
      p_in_Pa,                & ! Air pressure                              [Pa]
      thvm,                   & ! Virtual potential temperature              [K]
      wm_zm,                  & ! w wind component on momentum levels      [m/s]
      wm_zt,                  & ! w wind component on thermo. levels       [m/s]
      wp2,                    & ! w'^2 (momentum levels)               [m^2/s^2]
      coef_wp2xp_implicit,    & ! Coef. on <w'x'>; <w'^2 x'> eq.; t-levs.  [m/s]
      coef_wp2xp_implicit_zm, & ! coef_wp2xp_implicit interp. to m-levs.   [m/s]
      sgn_turbulent_vel,      & ! Sign of turbulent velocity ("upwind" ta)   [-]
      Kw6,                    & ! Coef. of eddy diffusivity for w'x'     [m^2/s]
      tau_C6_zm,              & ! Time-scale tau on momentum levels          [s]
      C7_Skw_fnc,             & ! C_7 parameter with Sk_w applied            [-]
      C6x_Skw_fnc,            & ! C_6x parameter with Sk_w applied           [-]
      rho_ds_zm,              & ! Dry, static density on momentum levs. [kg/m^3]
      rho_ds_zt,              & ! Dry, static density on thermo. levs.  [kg/m^3]
      invrs_rho_ds_zm,        & ! Inv. dry, static density at m-levs.   [m^3/kg]
      invrs_rho_ds_zt,        & ! Inv. dry, static density at t-levs.   [m^3/kg]
      wpxp_upper_lim,         & ! Keeps corrs. from becoming > 1       [un vary]
      wpxp_lower_lim            ! Keeps corrs. from becoming < -1      [un vary]

    logical, intent(in) ::  & 
      l_implemented ! Flag for CLUBB being implemented in a larger model.


    !------------------- Output Variable -------------------
    real( kind = core_rknd ), intent(out), dimension(nsup+nsub+1,2*gr%nz) ::  & 
      lhs ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)


    !------------------- Local Variables -------------------
    ! Indices
    integer :: k
    integer :: k_xm, k_wpxp

    logical :: l_upper_thresh, l_lower_thresh ! flags for clip_semi_imp_lhs

    ! These variables are used to change the amount
    ! of diffusion applied towards rtm and thlm. They are only used when
    ! l_diffuse_rtm_and_thlm = .true.
    real (kind = core_rknd), dimension(gr%nz) :: &
      zero_nu, &
      Kh_N2_zm, &
      K_zm          ! Coef. of eddy diffusivity at momentum level (k)   [m^2/s]

    real (kind = core_rknd) :: &
      constant_nu ! controls the magnitude of diffusion

    real( kind = core_rknd ), dimension(3,gr%nz) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for xm
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm,    & ! Mean advection contributions to lhs
      lhs_ta_wpxp     ! Turbulent advection contributions to lhs

    real( kind = core_rknd ), dimension(2,gr%nz) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      lhs_pr1,    & ! Pressure term 1 for w'x'
      lhs_ac_pr2    ! Accumulation of w'x' and w'x' pressure term 2
    
    real( kind = core_rknd ) :: invrs_dt    ! Reciprocal of dt, used for computational efficiency


    !------------------- Begin Code -------------------


    ! Initializations/precalculations
    invrs_dt    = 1.0_core_rknd / dt
    constant_nu = 0.1_core_rknd
    lhs         = 0.0_core_rknd


    if ( l_stability_correct_Kh_N2_zm ) then
      Kh_N2_zm = Kh_zm / calc_stability_correction( thlm, Lscale, em, exner, rtm, rcm, &
                                                    p_in_Pa, thvm )
    else
      Kh_N2_zm = Kh_zm
    end if


    
    ! Calculate diffusion terms for all momentum grid level
    call diffusion_zm_lhs_all( Kw6(:), nu6_vert_res_dep(:),      & ! Intent(in)
                               gr%invrs_dzt(:), gr%invrs_dzm(:), & ! Intent(in)
                               lhs_diff_zm(:,:)                  ) ! Intent(out) 


    ! Calculate diffusion terms for all thermodynamic grid level
    if ( l_diffuse_rtm_and_thlm ) then

        K_zm(:) = rho_ds_zm(:) * ( Kh_N2_zm(:) + constant_nu )
        zero_nu(:) = 0.0_core_rknd

        call diffusion_zt_lhs_all( K_zm(:), zero_nu(:),              & ! Intent(in)
                                   gr%invrs_dzm(:), gr%invrs_dzt(:), & ! Intent(in)
                                   lhs_diff_zt(:,:)                  ) ! Intent(out) 
        do k = 2, gr%nz 
            k_xm = 2*k - 1
            lhs(1,k_xm) = lhs(1,k_xm) + invrs_rho_ds_zt(k) * lhs_diff_zt(1,k) 
            lhs(3,k_xm) = lhs(3,k_xm) + invrs_rho_ds_zt(k) * lhs_diff_zt(2,k)
            lhs(5,k_xm) = lhs(5,k_xm) + invrs_rho_ds_zt(k) * lhs_diff_zt(3,k)
        end do

    end if


    ! Calculate mean advection terms for all momentum grid level
    call term_ma_zm_lhs_all( wm_zm(:), gr%invrs_dzm(:), & ! Intent(in)
                             lhs_ma_zm(:,:)           ) ! Intent(out)


    ! Calculate mean advection terms for all momentum grid level
    if ( .not. l_implemented ) then

        call term_ma_zt_lhs_all( wm_zt(:), gr%invrs_dzt(:), gr%invrs_dzm(:), & ! Intent(in)
                                 lhs_ma_zt(:,:)                            ) ! Intent(out)

        do k = 2, gr%nz 
            k_xm = 2*k - 1
            lhs(1,k_xm) = lhs(1,k_xm) + lhs_ma_zt(1,k)
            lhs(3,k_xm) = lhs(3,k_xm) + lhs_ma_zt(2,k)
            lhs(5,k_xm) = lhs(5,k_xm) + lhs_ma_zt(3,k)
        end do

    endif


    ! Calculate turbulent advection terms of w'x' for all grid levels
    ! An "over-implicit" weighted time step is applied to this term, see notes above
    call xpyp_term_ta_pdf_lhs_all( coef_wp2xp_implicit(:),      & ! Intent(in)
                                   rho_ds_zt(:),                & ! Intent(in)
                                   invrs_rho_ds_zm(:),          & ! Intent(in)
                                   gr%invrs_dzm(:),             & ! Intent(in)
                                   l_upwind_wpxp_ta,            & ! Intent(in)
                                   sgn_turbulent_vel(:),        & ! Intent(in)
                                   coef_wp2xp_implicit_zm(:),   & ! Intent(in)
                                   rho_ds_zm(:),                & ! Intent(in)
                                   gr%invrs_dzt(:),             & ! Intent(in)
                                   lhs_ta_wpxp(:,:)             ) ! Intent(out)



    ! Calculate turbulent advection terms of xm for all grid levels
    call xm_term_ta_lhs_all( rho_ds_zm(:),       & ! Intent(in)
                             invrs_rho_ds_zt(:), & ! Intent(in)
                             gr%invrs_dzt(:),    & ! Intent(in)
                             lhs_ta_xm(:,:)      ) ! Intent(out)


    ! Calculate turbulent production terms of w'x' for all grid level
    call wpxp_term_tp_lhs_all( wp2(:),          & ! Intent(in)
                               gr%invrs_dzm(:), & ! Intent(in)
                               lhs_tp(:,:) )      ! Intent(out)


    ! Calculate pressure term 1 for w'x' for all grid level
    ! Note:  An "over-implicit" weighted time step is applied to this term.
    call wpxp_term_pr1_lhs_all( C6x_Skw_fnc(:), & ! Intent(in)
                                tau_C6_zm(:),   & ! Intent(in)
                                lhs_pr1(:)      ) ! Intent(out)             


    ! Calculate accumulation of w'x' and w'x' pressure term 2 of w'x' for all grid level
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_pr
    call wpxp_terms_ac_pr2_lhs_all( C7_Skw_fnc(:),   & ! Intent(in)
                                    wm_zt(:),        & ! Intent(in)
                                    gr%invrs_dzm(:), & ! Intent(in)
                                    lhs_ac_pr2(:)    ) ! Intent(out)
        

    ! Lower boundary for xm, lhs(:,1)
    lhs(1,1) = 0.0_core_rknd
    lhs(2,1) = 0.0_core_rknd
    lhs(3,1) = 1.0_core_rknd
    lhs(4,1) = 0.0_core_rknd
    lhs(5,1) = 0.0_core_rknd

    ! Lower boundary for w'x', lhs(:,2)
    lhs(1,2) = 0.0_core_rknd
    lhs(2,2) = 0.0_core_rknd
    lhs(3,2) = 1.0_core_rknd
    lhs(4,2) = 0.0_core_rknd
    lhs(5,2) = 0.0_core_rknd

    ! Combine terms for lhs(:,3) to lhs(:,2*gr%nz)
    do k = 2, gr%nz

        k_xm = 2*k - 1  ! xm at odd index values
        k_wpxp = 2*k    ! w'x' at even index values

        ! ---- sum xm terms ----

        lhs(2,k_xm) = lhs(2,k_xm) + lhs_ta_xm(1,k)

        lhs(3,k_xm) = lhs(3,k_xm) + invrs_dt

        lhs(4,k_xm) = lhs(4,k_xm) + lhs_ta_xm(2,k)


        ! ---- sum w'x' terms ----

        lhs(1,k_wpxp) = lhs(1,k_wpxp) + lhs_ma_zm(1,k) + lhs_diff_zm(1,k) &
                                      + gamma_over_implicit_ts * lhs_ta_wpxp(1,k)

        lhs(2,k_wpxp) = lhs(2,k_wpxp) + lhs_tp(1,k)

        lhs(3,k_wpxp) = lhs(3,k_wpxp) + lhs_ma_zm(2,k) + lhs_diff_zm(2,k) + lhs_ac_pr2(k) &
                                      + gamma_over_implicit_ts * ( lhs_ta_wpxp(2,k) + lhs_pr1(k) )

        lhs(4,k_wpxp) = lhs(4,k_wpxp) + lhs_tp(2,k)

        lhs(5,k_wpxp) = lhs(5,k_wpxp) + lhs_ma_zm(3,k) + lhs_diff_zm(3,k) &
                                      + gamma_over_implicit_ts * lhs_ta_wpxp(3,k)

    enddo

    ! Upper boundary for w'x', , lhs(:,2*gr%nz)
    ! These were set in the loop above for memory access optimization
    lhs(1,2*gr%nz) = 0.0_core_rknd
    lhs(2,2*gr%nz) = 0.0_core_rknd
    lhs(3,2*gr%nz) = 1.0_core_rknd
    lhs(4,2*gr%nz) = 0.0_core_rknd
    lhs(5,2*gr%nz) = 0.0_core_rknd


    ! LHS time tendency
    if ( l_iter ) then
        do k = 2, gr%nz-1
            k_wpxp = 2*k 
            lhs(3,k_wpxp) = lhs(3,k_wpxp) + invrs_dt
        end do
    endif


    ! LHS portion of semi-implicit clipping term.
    if ( l_clip_semi_implicit ) then
        l_upper_thresh = .true.
        l_lower_thresh = .true.

        do k = 2, gr%nz-1
            k_wpxp = 2*k 
            lhs(3,k_wpxp) = lhs(3,k_wpxp) + clip_semi_imp_lhs( dt, wpxp(k),  & 
                                               l_upper_thresh, wpxp_upper_lim(k),  & 
                                               l_lower_thresh, wpxp_lower_lim(k) )
        end do
    endif


    ! Statistics contributions
    if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wprtp or wpthlp.

        if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
            ! Note:  An "over-implicit" weighted time step is applied to this
            !        term.  A weighting factor of greater than 1 may be used to
            !        make the term more numerically stable (see note above for
            !        LHS turbulent advection (ta) term).
            do k = 2, gr%nz-1
                zmscr01(k) = - lhs_ma_zm(3,k)
                zmscr02(k) = - lhs_ma_zm(2,k)
                zmscr03(k) = - lhs_ma_zm(1,k)
            end do
        endif

        if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
            do k = 2, gr%nz-1
                zmscr04(k) = - gamma_over_implicit_ts * lhs_ta_wpxp(3,k)
                zmscr05(k) = - gamma_over_implicit_ts * lhs_ta_wpxp(2,k)
                zmscr06(k) = - gamma_over_implicit_ts * lhs_ta_wpxp(1,k)
            end do
        endif

        if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
            do k = 2, gr%nz-1
                zmscr07(k) = - lhs_tp(2,k)
                zmscr08(k) = - lhs_tp(1,k)
            end do
        endif


        if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
            ! Note:  To find the contribution of w'x' term ac, substitute 0 for the
            !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
            do k = 2, gr%nz-1
                zmscr09(k) = - wpxp_terms_ac_pr2_lhs( zero, wm_zt(k+1), &
                                                      wm_zt(k), gr%invrs_dzm(k) )
            end do
        endif


        if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for LHS turbulent
            !        advection (ta) term).
            do k = 2, gr%nz-1
                zmscr10(k) = - gamma_over_implicit_ts * lhs_pr1(k)
            end do
        endif


        if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
            ! Note:  To find the contribution of w'x' term pr2, add 1 to the
            !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
            do k = 2, gr%nz-1
                zmscr11(k) = - wpxp_terms_ac_pr2_lhs( (one+C7_Skw_fnc(k)), wm_zt(k+1), &
                                                       wm_zt(k), gr%invrs_dzm(k) )
            end do
        endif

        if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
            do k = 2, gr%nz-1
                zmscr12(k) = - lhs_diff_zm(3,k)
                zmscr13(k) = - lhs_diff_zm(2,k)
                zmscr14(k) = - lhs_diff_zm(1,k)
            end do
        endif

        if ( l_clip_semi_implicit .and. ( iwprtp_sicl > 0 .or. iwpthlp_sicl > 0 ) ) then
            l_upper_thresh = .true.
            l_lower_thresh = .true.

            do k = 2, gr%nz-1
                zmscr15(k) = - clip_semi_imp_lhs( dt, wpxp(k),  & 
                                                  l_upper_thresh, wpxp_upper_lim(k), & 
                                                  l_lower_thresh, wpxp_lower_lim(k) )
            end do
        endif

        
        ! Statistics: implicit contributions for rtm or thlm.

        if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
            if ( .not. l_implemented ) then
                do k = 2, gr%nz
                    ztscr01(k) = - lhs_ma_zt(3,k)
                    ztscr02(k) = - lhs_ma_zt(2,k)
                    ztscr03(k) = - lhs_ma_zt(1,k)
                end do
            else
                do k = 2, gr%nz
                    ztscr01(k) = zero
                    ztscr02(k) = zero
                    ztscr03(k) = zero
                end do
            endif
        endif

        if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
            do k = 2, gr%nz
                ztscr04(k) = - lhs_ta_xm(2,k)
                ztscr05(k) = - lhs_ta_xm(1,k)
            end do
        endif

    endif

    return

  end subroutine xm_wpxp_lhs

  !=============================================================================
  subroutine xm_wpxp_rhs( solve_type, l_iter, dt, xm, wpxp, & ! In
                          xm_forcing, wpxp_forcing, C7_Skw_fnc, & ! In
                          xpthvp, C6x_Skw_fnc, tau_C6_zm, & ! In
                          coef_wp2xp_implicit, coef_wp2xp_implicit_zm, & ! In
                          term_wp2xp_explicit, term_wp2xp_explicit_zm, & ! In
                          sgn_turbulent_vel, rho_ds_zt, & ! In
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, & ! In
                          wpxp_upper_lim, wpxp_lower_lim, & ! In
                          rhs ) ! Out

    ! Description:
    ! Compute RHS vector for xm and w'x'.
    ! This subroutine computes the explicit portion of
    ! the xm and w'x' equations.
    !
    ! Notes:  
    !   For LHS turbulent advection (ta) term.
    !       An "over-implicit" weighted time step is applied to this term.
    !       The weight of the implicit portion of this term is controlled by
    !       the factor gamma_over_implicit_ts (abbreviated "gamma" in the
    !       expression below).  A factor is added to the right-hand side of
    !       the equation in order to balance a weight that is not equal to 1,
    !       such that:
    !            -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
    !       where X is the variable that is being solved for in a predictive
    !       equation (<w'x'> in this case), y(t) is the linearized portion of
    !       the term that gets treated implicitly, and RHS is the portion of
    !       the term that is always treated explicitly.  A weight of greater
    !       than 1 can be applied to make the term more numerically stable.
    !
    !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
    !   Significant changes to this routine may adversely affect computational speed
    !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
    !----------------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only:  &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_clip_semi_implicit,          & ! Variable(s)
        l_upwind_wpxp_ta

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs, & ! Procedure(s)
        xpyp_term_ta_pdf_rhs, &
        xpyp_term_ta_pdf_lhs_all, & ! Procedure(s)
        xpyp_term_ta_pdf_rhs_all

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use clip_semi_implicit, only: & 
        clip_semi_imp_rhs ! Procedure(s)

    use stats_type_utilities, only: & 
        stat_update_var_pt,   & 
        stat_begin_update_pt, &
        stat_modify_pt

    use stats_variables, only: & 
        stats_zt, & ! Variable(s)
        stats_zm, & 
        irtm_forcing, & 
        ithlm_forcing, & 
        iwprtp_bp, & 
        iwprtp_pr3, &
        iwprtp_sicl, &
        iwprtp_ta, &
        iwprtp_pr1, &
        iwprtp_forcing, &
        iwpthlp_bp, & 
        iwpthlp_pr3, & 
        iwpthlp_sicl, &
        iwpthlp_ta, &
        iwpthlp_pr1, &
        iwpthlp_forcing, &
        iupwp_bp, & 
        iupwp_pr3, &
        iupwp_ta, &
        iupwp_pr1, &
        ivpwp_bp, & 
        ivpwp_pr3, &
        ivpwp_ta, &
        ivpwp_pr1, &
        l_stats_samp

    use advance_helper_module, only: &
        set_boundary_conditions_rhs    ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      solve_type  ! Variables being solved for.

    logical, intent(in) :: l_iter

    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      xm,                     & ! xm (thermodynamic levels)               [x un]
      wpxp,                   & ! <w'x'> (momentum levels)          [{x un} m/s]
      xm_forcing,             & ! xm forcings (thermodynamic levels)  [{x un}/s]
      wpxp_forcing,           & ! <w'x'> forcing (momentum levs)  [{x un} m/s^2]
      C7_Skw_fnc,             & ! C_7 parameter with Sk_w applied            [-]
      xpthvp,                 & ! x'th_v' (momentum levels)           [{x un} K]
      C6x_Skw_fnc,            & ! C_6x parameter with Sk_w applied           [-]
      tau_C6_zm,              & ! Time-scale tau on momentum levels          [s]
      coef_wp2xp_implicit,    & ! Coef. on <w'x'>; <w'^2 x'> eq.; t-levs.  [m/s]
      coef_wp2xp_implicit_zm, & ! coef_wp2xp_implicit interp. to m-levs.   [m/s]
      term_wp2xp_explicit,    & ! Term that is on the RHS        [m^2/s^2{x un}]
      term_wp2xp_explicit_zm, & ! term_wp2xp_explicit interp. zm [m^2/s^2{x un}]
      sgn_turbulent_vel,      & ! Sign of turbulent velocity ("upwind" ta)   [-]
      rho_ds_zt,              & ! Dry, static density on thermo. levels [kg/m^3]
      rho_ds_zm,              & ! Dry, static density on momentum levs  [kg/m^3]
      invrs_rho_ds_zm,        & ! Inv. dry, static density at mom. levs [m^3/kg]
      thv_ds_zm,              & ! Dry, base-state theta_v on mom. levs.      [K]
      wpxp_upper_lim,         & ! Keeps corrs. from becoming > 1       [un vary]
      wpxp_lower_lim            ! Keeps corrs. from becoming < -1      [un vary]

    ! Output Variable
    real( kind = core_rknd ), intent(out), dimension(2*gr%nz) ::  & 
      rhs  ! Right-hand side of band diag. matrix. (LAPACK)

    ! Local Variables.

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3,gr%nz) :: &
        lhs_ta_wpxp   ! Turbulent advection terms of w'x'

    real( kind = core_rknd ), dimension(gr%nz) :: &
        lhs_pr1,    & ! Pressure term 1 for w'x'
        rhs_ta,     &
        rhs_bp_pr3    ! Buoyancy production of w'x' and w'x' pressure term 3
      
    real( kind = core_rknd ) :: &
        invrs_dt

    ! Indices
    integer :: k, k_xm, k_wpxp

    integer :: & 
      ixm_f, & 
      iwpxp_bp, & 
      iwpxp_pr3, &
      iwpxp_f, &
      iwpxp_sicl, &
      iwpxp_ta, &
      iwpxp_pr1

    logical :: l_upper_thresh, l_lower_thresh ! flags for clip_semi_imp_lhs

    ! ---- Begin Code ----

    ! Initialize output array and precalculate the reciprocal of dt
    invrs_dt = 1.0_core_rknd / dt
    rhs = 0.0_core_rknd


    ! Calculate turbulent advection terms of w'x' for all grid levels
    call xpyp_term_ta_pdf_rhs_all( term_wp2xp_explicit(:),      & ! Intent(in)
                                   rho_ds_zt(:),                & ! Intent(in)
                                   invrs_rho_ds_zm(:),          & ! Intent(in)
                                   gr%invrs_dzm(:),             & ! Intent(in)
                                   l_upwind_wpxp_ta,            & ! Intent(in)
                                   sgn_turbulent_vel(:),        & ! Intent(in)
                                   term_wp2xp_explicit_zm(:),   & ! Intent(in)
                                   rho_ds_zm(:),                & ! Intent(in)
                                   gr%invrs_dzt(:),             & ! Intent(in)
                                   rhs_ta(:)                    ) ! Intent(out)


    ! Calculate turbulent advection terms of w'x' for all grid levels
    ! An "over-implicit" weighted time step is applied to this term, see notes above
    call xpyp_term_ta_pdf_lhs_all( coef_wp2xp_implicit(:),      & ! Intent(in)
                                   rho_ds_zt(:),                & ! Intent(in)
                                   invrs_rho_ds_zm(:),          & ! Intent(in)
                                   gr%invrs_dzm(:),             & ! Intent(in)
                                   l_upwind_wpxp_ta,            & ! Intent(in)
                                   sgn_turbulent_vel(:),        & ! Intent(in)
                                   coef_wp2xp_implicit_zm(:),   & ! Intent(in)
                                   rho_ds_zm(:),                & ! Intent(in)
                                   gr%invrs_dzt(:),             & ! Intent(in)
                                   lhs_ta_wpxp(:,:)             ) ! Intent(out)


    ! Calculate pressure term 1 for w'x' for all grid level
    ! Note:  An "over-implicit" weighted time step is applied to this term.
    call wpxp_term_pr1_lhs_all( C6x_Skw_fnc(:), & ! Intent(in)
                                tau_C6_zm(:),   & ! Intent(in)
                                lhs_pr1(:)      ) ! Intent(out)


    ! Calculate buoyancy production of w'x' and w'x' pressure term 3
    call wpxp_terms_bp_pr3_rhs_all( C7_Skw_fnc(:), thv_ds_zm(:), xpthvp(:), &
                                    rhs_bp_pr3(:) )


    ! Set lower boundary for xm
    rhs(1) = xm(1)

    ! Set lower boundary for w'x'
    rhs(2) = wpxp(1)

    ! Combine terms to calculate other values, rhs(3) to rhs(gr%nz-2)
    do k = 2, gr%nz-1

        k_xm   = 2*k - 1
        k_wpxp = 2*k

        ! RHS time tendency and forcings for xm
        ! Note: xm forcings include the effects of microphysics,
        !       cloud water sedimentation, radiation, and any
        !       imposed forcings on xm.
        rhs(k_xm) = rhs(k_xm) + xm(k) * invrs_dt + xm_forcing(k)


        ! Calculate rhs values for w'x' using precalculated terms
        rhs(k_wpxp) = rhs(k_wpxp) + rhs_bp_pr3(k) + wpxp_forcing(k) + rhs_ta(k) &
                                  + ( one - gamma_over_implicit_ts ) &
                                  * ( - lhs_ta_wpxp(1,k) * wpxp(k+1) &
                                      - lhs_ta_wpxp(2,k) * wpxp(k) &
                                      - lhs_ta_wpxp(3,k) * wpxp(k-1) &
                                      - lhs_pr1(k) * wpxp(k) )
    enddo

    ! Upper boundary for xm
    rhs(2*gr%nz-1) = rhs(2*gr%nz-1) + xm(gr%nz) * invrs_dt + xm_forcing(gr%nz)

    ! Upper boundary for w'x', rhs(2*gr%nz)
    rhs(2*gr%nz) = 0.0_core_rknd


    ! RHS time tendency.
    if ( l_iter ) then
        do k = 2, gr%nz-1
            k_wpxp = 2*k
            rhs(k_wpxp) = rhs(k_wpxp) + wpxp(k) * invrs_dt
        end do
    end if


    ! RHS portion of semi-implicit clipping (sicl) term.
    if ( l_clip_semi_implicit ) then
        l_upper_thresh = .true.
        l_lower_thresh = .true.

        do k = 2, gr%nz-1
            k_wpxp = 2*k
            rhs(k_wpxp) = rhs(k_wpxp) + clip_semi_imp_rhs( dt, wpxp(k), & 
                                                           l_upper_thresh, wpxp_upper_lim(k), & 
                                                           l_lower_thresh, wpxp_lower_lim(k) )
        end do

    endif


    if ( l_stats_samp ) then

        select case ( solve_type )
            case ( xm_wpxp_rtm )  ! rtm/wprtp budget terms
              ixm_f      = irtm_forcing
              iwpxp_bp   = iwprtp_bp
              iwpxp_pr3  = iwprtp_pr3
              iwpxp_f    = iwprtp_forcing
              iwpxp_sicl = iwprtp_sicl
              iwpxp_ta   = iwprtp_ta
              iwpxp_pr1  = iwprtp_pr1
            case ( xm_wpxp_thlm ) ! thlm/wpthlp budget terms
              ixm_f      = ithlm_forcing
              iwpxp_bp   = iwpthlp_bp
              iwpxp_pr3  = iwpthlp_pr3
              iwpxp_f    = iwpthlp_forcing
              iwpxp_sicl = iwpthlp_sicl
              iwpxp_ta   = iwpthlp_ta
              iwpxp_pr1  = iwpthlp_pr1
            case ( xm_wpxp_um )  ! um/upwp budget terms
              ixm_f      = 0
              iwpxp_bp   = iupwp_bp
              iwpxp_pr3  = iupwp_pr3
              iwpxp_f    = 0
              iwpxp_sicl = 0
              iwpxp_ta   = iupwp_ta
              iwpxp_pr1  = iupwp_pr1
            case ( xm_wpxp_vm )  ! vm/vpwp budget terms
              ixm_f      = 0
              iwpxp_bp   = ivpwp_bp
              iwpxp_pr3  = ivpwp_pr3
              iwpxp_f    = 0
              iwpxp_sicl = 0
              iwpxp_ta   = ivpwp_ta
              iwpxp_pr1  = ivpwp_pr1
            case default    ! this includes the sclrm case
              ixm_f      = 0
              iwpxp_bp   = 0
              iwpxp_pr3  = 0
              iwpxp_f    = 0
              iwpxp_sicl = 0
              iwpxp_ta   = 0
              iwpxp_pr1  = 0
        end select

        do k = 2, gr%nz-1

            ! Statistics: explicit contributions for wpxp.

            ! w'x' term bp is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'x' term bp, substitute 0 for the
            !        C_7 skewness function input to function wpxp_terms_bp_pr3_rhs.
            call stat_update_var_pt( iwpxp_bp, k, & 
                wpxp_terms_bp_pr3_rhs( zero, thv_ds_zm(k), xpthvp(k) ), stats_zm )


            ! w'x' term pr3 is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of w'x' term pr3, add 1 to the
            !        C_7 skewness function input to function wpxp_terms_bp_pr2_rhs.
            call stat_update_var_pt( iwpxp_pr3, k, & 
                wpxp_terms_bp_pr3_rhs( (one+C7_Skw_fnc(k)), thv_ds_zm(k), &
                                       xpthvp(k) ), &
                                     stats_zm )

            ! w'x' forcing term is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( iwpxp_f, k, wpxp_forcing(k), stats_zm )

            ! w'x' term sicl has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on clip_semi_imp_rhs.
            if ( l_clip_semi_implicit ) then
                l_upper_thresh = .true.
                l_lower_thresh = .true.

                call stat_begin_update_pt( iwpxp_sicl, k, & 
                                           -clip_semi_imp_rhs( dt, wpxp(k), & 
                                           l_upper_thresh, wpxp_upper_lim(k), & 
                                           l_lower_thresh, wpxp_lower_lim(k) ), stats_zm )
            endif

            ! <w'x'> term ta has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on
            ! xpyp_term_ta_pdf_rhs.
            call stat_begin_update_pt( iwpxp_ta, k, &
                                       -rhs_ta(k), &
                                       stats_zm )

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS
            !        contribution from "over-implicit" weighted time step for LHS
            !        turbulent advection (ta) term).
            call stat_modify_pt( iwpxp_ta, k, &
                                 + ( one - gamma_over_implicit_ts ) &
                                   * ( - lhs_ta_wpxp(1,k) * wpxp(k+1) &
                                       - lhs_ta_wpxp(2,k) * wpxp(k) &
                                       - lhs_ta_wpxp(3,k) * wpxp(k-1) ), &
                                 stats_zm )

            ! w'x' term pr1 is normally completely implicit.  However, there is a
            ! RHS contribution from the "over-implicit" weighted time step.  A
            ! weighting factor of greater than 1 may be used to make the term more
            ! numerically stable (see note above for RHS contribution from
            ! "over-implicit" weighted time step for LHS turbulent advection (ta)
            ! term).  Therefore, w'x' term pr1 has both implicit and explicit
            ! components; call stat_begin_update_pt.  Since stat_begin_update_pt
            ! automatically subtracts the value sent in, reverse the sign on the
            ! input value.
            call stat_begin_update_pt( iwpxp_pr1, k, &
                                      - ( one - gamma_over_implicit_ts )  &
                                      * ( - lhs_pr1(k) * wpxp(k) ), stats_zm )
        end do


        ! Statistics: explicit contributions for xm
        !             (including microphysics/radiation).

        ! xm forcings term is completely explicit; call stat_update_var_pt.
        do k = 2, gr%nz
        
            call stat_update_var_pt( ixm_f, k, xm_forcing(k), stats_zt )

        end do

    endif ! l_stats_samp

    return

  end subroutine xm_wpxp_rhs

  !=============================================================================
  subroutine xm_wpxp_solve( nrhs, lhs, rhs, solution, rcond )

    ! Description:
    !   Solve for xm / w'x' using the band diagonal solver.

    ! References:
    !   None
    !------------------------------------------------------------------------

    use grid_class, only: & 
      gr ! Variable(s)

    use lapack_wrap, only:  & 
      band_solve,  & ! Procedure(s)
      band_solvex

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: &
      fstderr     ! Constant(s)
    
    use error_code, only: &
      clubb_at_least_debug_level,     & ! Procedure
      err_code,                       & ! Error indicator
      clubb_no_error                    ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nrhs ! Number of rhs vectors

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nsup+nsub+1,2*gr%nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix in LAPACK storage)

    real( kind = core_rknd ), intent(inout), dimension(2*gr%nz,nrhs) ::  & 
      rhs      ! Right-hand side of band diag. matrix. (LAPACK storage)

    real( kind = core_rknd ), intent(out), dimension(2*gr%nz,nrhs) ::  & 
      solution ! Solution to band diagonal system (LAPACK storage)

    ! Output Variables
    real( kind = core_rknd ), optional, intent(out) :: &
      rcond ! Est. of the reciprocal of the condition #

    if ( present( rcond ) ) then
      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call band_solvex( "xm_wpxp", nsup, nsub, 2*gr%nz, nrhs, & 
                        lhs, rhs, solution, rcond )
    else
      ! Perform LU decomp and solve system (LAPACK)
      call band_solve( "xm_wpxp", nsup, nsub, 2*gr%nz, nrhs, & 
                       lhs, rhs, solution )
    end if

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code /= clubb_no_error ) then
            write(fstderr,*) "Error in xm_wpxp_solve"
            return
        end if
    end if

    return
  end subroutine xm_wpxp_solve

!===============================================================================
  subroutine xm_wpxp_clipping_and_stats &
             ( solve_type, dt, wp2, xp2, wm_zt, &
               xm_forcing, rho_ds_zm, rho_ds_zt, &
               invrs_rho_ds_zm, invrs_rho_ds_zt, &
               xp2_threshold, xm_threshold, rcond, &
               low_lev_effect, high_lev_effect, &
               l_implemented, solution, &
               xm, xm_tol, wpxp )

    ! Description:
    ! Clips and computes implicit stats for an artitrary xm and wpxp
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use model_flags, only: &
        l_clip_semi_implicit, & ! Variable(s)
        l_tke_aniso

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use mono_flux_limiter, only: &
        monotonic_turbulent_flux_limit ! Procedure(s)

    use pos_definite_module, only:  & 
        pos_definite_adj ! Procedure(s)

    use clip_explicit, only: & 
        clip_covar,   & ! Procedure(s)
        clip_wprtp,   & ! Variable(s)
        clip_wpthlp,  &
        clip_upwp,    &
        clip_vpwp,    &
        clip_wpsclrp

    use model_flags, only: & 
        l_pos_def, &     ! Logical for whether to apply the positive definite scheme to rtm
        l_hole_fill, &   ! Logical for whether to apply the hole filling scheme to thlm/rtm
        l_clip_turb_adv  ! Logical for whether to clip xm when wpxp is clipped

    use constants_clubb, only: &
        fstderr, & ! Constant(s)
        one, &
        zero, &
        eps

    use fill_holes, only: &
        fill_holes_vertical ! Procedure

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

    use stats_type_utilities, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_update_var_pt, & 
        stat_end_update_pt, & 
        stat_end_update,  & 
        stat_update_var, & 
        stat_modify

    use stats_variables, only: & 
        stats_zt,  & ! Variable(s)
        stats_zm, & 
        stats_sfc, & 
        irtm_ta, & 
        irtm_ma, & 
        irtm_matrix_condt_num, & 
        irtm_pd, & 
        irtm_cl,  & 
        iwprtp_ma, & 
        iwprtp_ta, & 
        iwprtp_tp, & 
        iwprtp_ac, & 
        iwprtp_pr1, & 
        iwprtp_pr2, & 
        iwprtp_dp1, & 
        iwprtp_pd, & 
        iwprtp_sicl, & 
        ithlm_ta

    use stats_variables, only: &
        ithlm_ma, & 
        ithlm_cl, & 
        ithlm_matrix_condt_num, & 
        iwpthlp_ma, & 
        iwpthlp_ta, & 
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_dp1, & 
        iwpthlp_sicl, &
        ium_ma, &
        ium_ta, &
        iupwp_ma, &
        iupwp_ta, &
        iupwp_tp, &
        iupwp_ac, &
        iupwp_pr1, &
        iupwp_pr2, &
        iupwp_dp1, &
        ivm_ma, &
        ivm_ta, &
        ivpwp_ma, &
        ivpwp_ta, &
        ivpwp_tp, &
        ivpwp_ac, &
        ivpwp_pr1, &
        ivpwp_pr2, &
        ivpwp_dp1

    use stats_variables, only: & 
        l_stats_samp, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        zmscr01, & 
        zmscr02, & 
        zmscr03, & 
        zmscr04, & 
        zmscr05, & 
        zmscr06, & 
        zmscr07, & 
        zmscr08, & 
        zmscr09, & 
        zmscr10, & 
        zmscr11, & 
        zmscr12, & 
        zmscr13, & 
        zmscr14, & 
        zmscr15

    implicit none

    ! Constant Parameters
    logical, parameter :: &
      l_mono_flux_lim = .true., &  ! Flag for monotonic turbulent flux limiter
      l_enable_relaxed_clipping = .true., & ! Flag to relax clipping
      l_first_clip_ts = .true., &
      l_last_clip_ts  = .false.

    ! Input Variables
    integer, intent(in) ::  & 
      solve_type  ! Variables being solved for.

    real( kind = core_rknd ), intent(in) ::  & 
      dt  ! Timestep   [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      xp2,             & ! x'^2 (momentum levels)                   [{xm units}^2]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      xm_forcing,      & ! xm forcings (thermodynamic levels)       [units vary]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), intent(in) :: &
      xp2_threshold, & ! Minimum allowable value of x'^2   [units vary]
      xm_threshold,  & ! Minimum allowable value of xm     [units vary]
      xm_tol,        & ! Minimum allowable deviation of xm [units vary]
      rcond ! Reciprocal of the estimated condition number (from computing A^-1)

    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(gr%nz), intent(in) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.

    logical, intent(in) :: &
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    real( kind = core_rknd ), intent(in), dimension(2*gr%nz) :: &
      solution ! The <t+1> value of xm and wpxp   [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
      xm, &     ! The mean x field  [units vary]
      wpxp      ! The flux of x     [units vary m/s]

    ! Local Variables
    integer :: & 
      solve_type_cl ! solve_type used for clipping statistics.

    character(len=10) :: &
      solve_type_str ! solve_type as a string for debug output purposes

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      xm_n ! Old value of xm for positive definite scheme     [units vary]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      wpxp_pd, xm_pd ! Change in xm and wpxp due to the pos. def. scheme

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wpxp_chnge, &  ! Net change in w'x' due to clipping       [units vary]
      xp2_relaxed    ! Value of x'^2 * clip_factor               [units vary]

    ! Indices
    integer :: &
      k, km1, kp1, &
      k_xm, k_wpxp

    integer :: & 
      ixm_ta, & 
      ixm_ma, & 
      ixm_matrix_condt_num, & 
      ixm_pd, & 
      ixm_cl, & 
      iwpxp_ma, & 
      iwpxp_ta, & 
      iwpxp_tp, & 
      iwpxp_ac, & 
      iwpxp_pr1, & 
      iwpxp_pr2, & 
      iwpxp_dp1, & 
      iwpxp_pd, & 
      iwpxp_sicl

    ! ----- Begin code ------

    select case ( solve_type )

    case ( xm_wpxp_rtm ) ! rtm/wprtp budget terms
      ixm_ta     = irtm_ta
      ixm_ma     = irtm_ma
      ixm_pd     = irtm_pd
      ixm_cl     = irtm_cl
      iwpxp_ma   = iwprtp_ma
      iwpxp_ta   = iwprtp_ta
      iwpxp_tp   = iwprtp_tp
      iwpxp_ac   = iwprtp_ac
      iwpxp_pr1  = iwprtp_pr1
      iwpxp_pr2  = iwprtp_pr2
      iwpxp_dp1  = iwprtp_dp1
      iwpxp_pd   = iwprtp_pd
      iwpxp_sicl = iwprtp_sicl

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = irtm_matrix_condt_num

    case ( xm_wpxp_thlm ) ! thlm/wpthlp budget terms
      ixm_ta     = ithlm_ta
      ixm_ma     = ithlm_ma
      ixm_pd     = 0
      ixm_cl     = ithlm_cl
      iwpxp_ma   = iwpthlp_ma
      iwpxp_ta   = iwpthlp_ta
      iwpxp_tp   = iwpthlp_tp
      iwpxp_ac   = iwpthlp_ac
      iwpxp_pr1  = iwpthlp_pr1
      iwpxp_pr2  = iwpthlp_pr2
      iwpxp_dp1  = iwpthlp_dp1
      iwpxp_pd   = 0
      iwpxp_sicl = iwpthlp_sicl

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = ithlm_matrix_condt_num

    case ( xm_wpxp_um ) ! um/upwp budget terms
      ixm_ta     = ium_ta
      ixm_ma     = ium_ma
      ixm_pd     = 0
      ixm_cl     = 0
      iwpxp_ma   = iupwp_ma
      iwpxp_ta   = iupwp_ta
      iwpxp_tp   = iupwp_tp
      iwpxp_ac   = iupwp_ac
      iwpxp_pr1  = iupwp_pr1
      iwpxp_pr2  = iupwp_pr2
      iwpxp_dp1  = iupwp_dp1
      iwpxp_pd   = 0
      iwpxp_sicl = 0

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = 0

    case ( xm_wpxp_vm ) ! vm/vpwp budget terms
      ixm_ta     = ivm_ta
      ixm_ma     = ivm_ma
      ixm_pd     = 0
      ixm_cl     = 0
      iwpxp_ma   = ivpwp_ma
      iwpxp_ta   = ivpwp_ta
      iwpxp_tp   = ivpwp_tp
      iwpxp_ac   = ivpwp_ac
      iwpxp_pr1  = ivpwp_pr1
      iwpxp_pr2  = ivpwp_pr2
      iwpxp_dp1  = ivpwp_dp1
      iwpxp_pd   = 0
      iwpxp_sicl = 0

      ! This is a diagnostic from inverting the matrix, not a budget
      ixm_matrix_condt_num = 0

    case default  ! this includes the sclrm case
      ixm_ta     = 0
      ixm_ma     = 0
      ixm_pd     = 0
      ixm_cl     = 0
      iwpxp_ma   = 0
      iwpxp_ta   = 0
      iwpxp_tp   = 0
      iwpxp_ac   = 0
      iwpxp_pr1  = 0
      iwpxp_pr2  = 0
      iwpxp_dp1  = 0
      iwpxp_pd   = 0
      iwpxp_sicl = 0

      ixm_matrix_condt_num = 0

    end select

    ! Copy result into output arrays

    do k=1, gr%nz, 1

      k_xm   = 2 * k - 1
      k_wpxp = 2 * k

      xm_n(k) = xm(k)

      xm(k)   = solution(k_xm)
      wpxp(k) = solution(k_wpxp)

    end do ! k=1..gr%nz

    ! Lower boundary condition on xm
    xm(1) = xm(2)


    if ( l_stats_samp ) then


      if ( ixm_matrix_condt_num > 0 ) then
        ! Est. of the condition number of the mean/flux LHS matrix
        call stat_update_var_pt( ixm_matrix_condt_num, 1, one / rcond, stats_sfc )
      end if


      ! The xm loop runs between k = 2 and k = gr%nz.  The value of xm at
      ! level k = 1, which is below the model surface, is simply set equal to
      ! the value of xm at level k = 2 after the solve has been completed.
      ! Thus, the statistical code will run from levels 2 through gr%nz.

      do k = 2, gr%nz

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nz )

        ! Finalize implicit contributions for xm

        ! xm term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_ma, k, & 
            ztscr01(k) * xm(km1) & 
          + ztscr02(k) * xm(k) & 
          + ztscr03(k) * xm(kp1), stats_zt )

        ! xm term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_ta, k, & 
            ztscr04(k) * wpxp(km1) & 
          + ztscr05(k) * wpxp(k), stats_zt )

      enddo ! xm loop: 2..gr%nz


      ! The wpxp loop runs between k = 2 and k = gr%nz-1.  The value of wpxp
      ! is set to specified values at both the lowest level, k = 1, and the
      ! highest level, k = gr%nz.  Thus, the statistical code will run from
      ! levels 2 through gr%nz-1.

      do k = 2, gr%nz-1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nz )

        ! Finalize implicit contributions for wpxp

        ! w'x' term ma is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_ma, k, & 
            zmscr01(k) * wpxp(km1) & 
          + zmscr02(k) * wpxp(k) & 
          + zmscr03(k) * wpxp(kp1), stats_zm )

!       if( .not. l_upwind_wpxp_ta ) then
          ! w'x' term ta is normally completely implicit.  However, due to the
          ! RHS contribution from the "over-implicit" weighted time step,
          ! w'x' term ta has both implicit and explicit components;
          ! call stat_end_update_pt.
          call stat_end_update_pt( iwpxp_ta, k, & 
              zmscr04(k) * wpxp(km1) & 
            + zmscr05(k) * wpxp(k) & 
            + zmscr06(k) * wpxp(kp1), stats_zm )
!       endif

        ! w'x' term tp is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_tp, k, & 
            zmscr07(k) * xm(k) & 
          + zmscr08(k) * xm(kp1), stats_zm )

        ! w'x' term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_ac, k, & 
            zmscr09(k) * wpxp(k), stats_zm )

        ! w'x' term pr1 is normally completely implicit.  However, due to the
        ! RHS contribution from the "over-implicit" weighted time step,
        ! w'x' term pr1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwpxp_pr1, k, & 
            zmscr10(k) * wpxp(k), stats_zm )

        ! w'x' term pr2 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_pr2, k, & 
            zmscr11(k) * wpxp(k), stats_zm )

        ! w'x' term dp1 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_dp1, k, & 
            zmscr12(k) * wpxp(km1) & 
          + zmscr13(k) * wpxp(k) & 
          + zmscr14(k) * wpxp(kp1), stats_zm )

        ! w'x' term sicl has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_clip_semi_implicit ) then
          call stat_end_update_pt( iwpxp_sicl, k, & 
              zmscr15(k) * wpxp(k), stats_zm )
        endif

      enddo ! wpxp loop: 2..gr%nz-1


    endif ! l_stats_samp


    ! Apply a monotonic turbulent flux limiter to xm/w'x'.
    if ( l_mono_flux_lim ) then
      call monotonic_turbulent_flux_limit( solve_type, dt, xm_n, &
                                           xp2, wm_zt, xm_forcing, &
                                           rho_ds_zm, rho_ds_zt, &
                                           invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                           xp2_threshold, l_implemented, &
                                           low_lev_effect, high_lev_effect, &
                                           xm, xm_tol, wpxp )
    end if ! l_mono_flux_lim

    ! Apply a flux limiting positive definite scheme if the solution
    ! for the mean field is negative and we're determining total water
    if ( solve_type == xm_wpxp_rtm .and. l_pos_def .and. any( xm < zero ) ) then

      call pos_definite_adj( dt, "zt", xm, wpxp, & 
                             xm_n, xm_pd, wpxp_pd )

    else
      ! For stats purposes
      xm_pd   = zero
      wpxp_pd = zero

    end if ! l_pos_def and solve_type == "rtm" and rtm <n+1> less than 0

    if ( l_stats_samp ) then

      call stat_update_var( iwpxp_pd, wpxp_pd(1:gr%nz), stats_zm )

      call stat_update_var( ixm_pd, xm_pd(1:gr%nz), stats_zt )

    end if

    ! Computed value before clipping
    if ( l_stats_samp ) then
      call stat_begin_update( ixm_cl, xm / dt, & ! Intent(in)
                              stats_zt )                       ! Intent(inout)
    end if

    if ( any( xm < xm_threshold ) .and. l_hole_fill &
         .and. solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm ) then

      select case ( solve_type )
      case ( xm_wpxp_rtm )
        solve_type_str = "rtm"
      case ( xm_wpxp_thlm )
        solve_type_str = "thlm"
      case default
        solve_type_str = "scalars"
      end select

      if ( clubb_at_least_debug_level( 1 ) ) then
        do k = 1, gr%nz
          if ( xm(k) < zero ) then
            write(fstderr,*) solve_type_str//" < ", xm_threshold, &
              " in advance_xm_wpxp_module at k= ", k
          end if
        end do
      end if

      call fill_holes_vertical( 2, xm_threshold, "zt", &
                                rho_ds_zt, rho_ds_zm, &
                                xm )

    endif ! any( xm < xm_threshold ) .and. l_hole_fill
          ! .and. solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm

    if ( l_stats_samp ) then
      call stat_end_update( ixm_cl, xm / dt, & ! Intent(in) 
                            stats_zt )                       ! Intent(inout)
    end if

    ! Use solve_type to find solve_type_cl, which is used
    ! in subroutine clip_covar.
    select case ( solve_type )
    case ( xm_wpxp_rtm )
      solve_type_cl = clip_wprtp
    case ( xm_wpxp_thlm )
      solve_type_cl = clip_wpthlp
    case ( xm_wpxp_um )
      solve_type_cl = clip_upwp
    case ( xm_wpxp_vm )
      solve_type_cl = clip_vpwp
    case default
      solve_type_cl = clip_wpsclrp
    end select

    ! Clipping for w'x'
    ! Clipping w'x' at each vertical level, based on the
    ! correlation of w and x at each vertical level, such that:
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
    ! -1 <= corr_(w,x) <= 1.
    ! Since w'^2, x'^2, and w'x' are updated in different places
    ! from each other, clipping for w'x' has to be done three times
    ! (three times each for w'r_t', w'th_l', and w'sclr').  This is
    ! the second instance of w'x' clipping.

    ! Compute a slightly larger value of rt'^2 for clipping purposes.  This was
    ! added to prevent a situation in which both the variance and flux are small
    ! and the simulation gets "stuck" at the rt_tol^2 value.
    ! See ticket #389 on the CLUBB TRAC for further details.
    ! -dschanen 10 Jan 2011
    if ( l_enable_relaxed_clipping ) then
      if ( solve_type == xm_wpxp_rtm ) then
        xp2_relaxed = max( 1e-7_core_rknd , xp2 )

      else if ( solve_type == xm_wpxp_thlm ) then
        xp2_relaxed = max( 0.01_core_rknd, xp2 )

      else ! This includes the passive scalars
        xp2_relaxed = max( 1e-7_core_rknd , xp2 )

      end if

    else  ! Don't relax clipping
      xp2_relaxed = xp2

    end if

    if ( solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm ) then

       call clip_covar( solve_type_cl, l_first_clip_ts, &  ! In
                        l_last_clip_ts, dt, wp2, xp2_relaxed, &  ! In
                        wpxp, wpxp_chnge ) ! In/Out

    else ! clipping for upwp or vpwp

       if ( l_tke_aniso ) then

          call clip_covar( solve_type_cl, l_first_clip_ts, &  ! In
                           l_last_clip_ts, dt, wp2, xp2, &  ! In
                           wpxp, wpxp_chnge ) ! In/Out

       else

          call clip_covar( solve_type_cl, l_first_clip_ts, &  ! In
                           l_last_clip_ts, dt, wp2, wp2, &  ! In
                           wpxp, wpxp_chnge ) ! In/Out

       endif ! l_tke_aniso

    endif ! solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm

    ! Adjusting xm based on clipping for w'x'.
    if ( any( abs(wpxp_chnge) > eps ) .and. l_clip_turb_adv ) then
      call xm_correction_wpxp_cl( solve_type, dt, wpxp_chnge, gr%invrs_dzt, &
                                  xm )
    endif


    return
  end subroutine xm_wpxp_clipping_and_stats

  !=============================================================================
  pure function xm_term_ta_lhs( rho_ds_zm, rho_ds_zmm1, &
                                invrs_rho_ds_zt, invrs_dzt ) &
  result( lhs )

    ! Description:
    ! Turbulent advection of xm:  implicit portion of the code.
    !
    ! The d(xm)/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'x' )/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * w'x'(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(xm)/dt and
    ! d(w'x')/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! While the values of xm are found on the thermodynamic levels, the values
    ! of w'x' are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  On the momentum
    ! levels, the values of rho_ds_zm are multiplied by the values of w'x'.  The
    ! derivative of (rho_ds_zm * w'x') is taken over the intermediate (central)
    ! thermodynamic level, where it is multiplied by invrs_rho_ds_zt, yielding
    ! the desired results.
    !
    ! =====rho_ds_zm=====wpxp================================== m(k)
    !
    ! ------invrs_rho_ds_zt--------d(rho_ds*wpxp)/dz----------- t(k)
    !
    ! =====rho_ds_zmm1===wpxpm1================================ m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1,    & ! Momentum superdiagonal index.
      km1_mdiag = 2       ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rho_ds_zm,       & ! Dry, static density at momentum level (k)    [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k+1)  [kg/m^3]
      invrs_rho_ds_zt, & ! Inverse dry, static density @ thermo lev (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                  [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2) :: lhs


    ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
    lhs(k_mdiag) & 
    = + invrs_rho_ds_zt * invrs_dzt * rho_ds_zm

    ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
    lhs(km1_mdiag) & 
    = - invrs_rho_ds_zt * invrs_dzt * rho_ds_zmm1


    return
  end function xm_term_ta_lhs

  !=============================================================================
  pure subroutine xm_term_ta_lhs_all( rho_ds_zm, &
                                      invrs_rho_ds_zt, &
                                      invrs_dzt, &
                                      lhs_ta_xm )
  ! Description:
  !     This subroutine serves the same function as xm_term_ta_lhs (above), but
  !     calculated terms for all grid levels at once rather than one at a time.
  !     This was done so that this code could be vectorized and thereby sped up
  !     by the compiler. See clubb:ticket:834 for more information.
  ! 
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only:  & 
        gr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zm,       & ! Dry, static density at momentum level (k)    [kg/m^3]
      invrs_rho_ds_zt, & ! Inverse dry, static density @ thermo lev (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                  [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: &
        lhs_ta_xm

    ! Local Variable
    integer :: k

    ! Set lower boundary condition to 0
    lhs_ta_xm(1,1) = 0.0_core_rknd
    lhs_ta_xm(2,1) = 0.0_core_rknd

    ! Calculate all other grid levels
    do k = 2, gr%nz 

        ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
        lhs_ta_xm(1,k) = + invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k)

        ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
        lhs_ta_xm(2,k) = - invrs_rho_ds_zt(k) * invrs_dzt(k) * rho_ds_zm(k-1)

    end do

    return

  end subroutine xm_term_ta_lhs_all

  !=============================================================================
  pure function wpxp_term_tp_lhs( wp2, invrs_dzm ) & 
  result( lhs )

    ! Description:
    ! Turbulent production of w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains a turbulent production term:
    !
    ! - w'^2 d(xm)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w'^2 * d( xm(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of xm being used is from the
    ! next timestep, which is being advanced to in solving the d(w'x')/dt and
    ! d(xm)/dt equations.
    !
    ! This term is discretized as follows:
    !
    ! The values of xm are found on thermodynamic levels, while the values of
    ! w'^2 are found on momentum levels.  The derivative of xm is taken over the
    ! intermediate (central) momentum level, where it is multiplied by w'^2,
    ! yielding the desired result.
    !
    ! ---------------------------xmp1-------------------------- t(k+1)
    !
    ! ==========wp2=====================d(xm)/dz=============== m(k)
    !
    ! ---------------------------xm---------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag = 2         ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      wp2,      & ! w'^2(k)                       [m^2/s^2]
      invrs_dzm   ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2) :: lhs


    ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
    lhs(kp1_tdiag) & 
    = + wp2 * invrs_dzm

    ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
    lhs(k_tdiag) & 
    = - wp2 * invrs_dzm


    return
  end function wpxp_term_tp_lhs
    
  !=============================================================================
  pure subroutine wpxp_term_tp_lhs_all( wp2, &
                                        invrs_dzm, & 
                                        lhs_tp )
  ! Description:
  !     This subroutine serves the same function as wpxp_term_tp_lhs (above), but
  !     calculated terms for all grid levels at once rather than one at a time.
  !     This was done so that this code could be vectorized and thereby sped up
  !     by the compiler. See clubb:ticket:834 for more information.
  ! 
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only:  & 
        gr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wp2,      & ! w'^2(k)                       [m^2/s^2]
      invrs_dzm   ! Inverse of grid spacing (k)   [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,gr%nz), intent(out) :: &
        lhs_tp

    ! Local Variable
    integer :: k

    ! Set lower boundary to 0
    lhs_tp(1,1) = 0.0_core_rknd
    lhs_tp(2,1) = 0.0_core_rknd

    ! Calculate all other grid levels
    do k = 2, gr%nz-1

        ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
        lhs_tp(1,k) = + wp2(k) * invrs_dzm(k)

        ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
        lhs_tp(2,k) = - wp2(k) * invrs_dzm(k)

    end do

    ! Set upper boundary to 0
    lhs_tp(1,gr%nz) = 0.0_core_rknd
    lhs_tp(2,gr%nz) = 0.0_core_rknd


    return
  end subroutine wpxp_term_tp_lhs_all

  !=============================================================================
  pure function wpxp_terms_ac_pr2_lhs( C7_Skw_fnc,  & 
                                       wm_ztp1, wm_zt, invrs_dzm ) & 
  result( lhs )

    ! Description:
    ! Accumulation of w'x' and w'x' pressure term 2:  implicit portion of the
    ! code.
    !
    ! The d(w'x')/dt equation contains an accumulation term:
    !
    ! - w'x' dw/dz;
    !
    ! and pressure term 2:
    !
    ! + C_7 w'x' dw/dz.
    !
    ! Both the w'x' accumulation term and pressure term 2 are completely
    ! implicit.  The accumulation term and pressure term 2 are combined and
    ! solved together as:
    !
    ! - ( 1 - C_7 ) * w'x'(t+1) * dw/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! The terms are discretized as follows:
    !
    ! The values of w'x' are found on momentum levels, while the values of wm_zt
    ! (mean vertical velocity on thermodynamic levels) are found on
    ! thermodynamic levels.  The vertical derivative of wm_zt is taken over the
    ! intermediate (central) momentum level.  It is then multiplied by w'x'
    ! (implicitly calculated at timestep (t+1)) and the coefficients to yield
    ! the desired results.
    !
    ! -------wm_ztp1------------------------------------------- t(k+1)
    !
    ! ===============d(wm_zt)/dz============wpxp=============== m(k)
    !
    ! -------wm_zt--------------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)             [-]
      wm_ztp1,     & ! w wind component on thermodynamic level (k+1)   [m/s]
      wm_zt,       & ! w wind component on thermodynamic level (k)     [m/s]
      invrs_dzm      ! Inverse of grid spacing (k)                     [1/m]


    ! Return Variable
    real( kind = core_rknd ) :: lhs


    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs = ( one - C7_Skw_fnc ) * invrs_dzm * ( wm_ztp1 - wm_zt )


    return
  end function wpxp_terms_ac_pr2_lhs

  !====================================================================================
  pure subroutine wpxp_terms_ac_pr2_lhs_all( C7_Skw_fnc, & 
                                             wm_zt,      &
                                             invrs_dzm,  &
                                             lhs_ac_pr2  ) 
  ! Description:
  !     This subroutine serves the same function as wpxp_terms_ac_pr2_lhs (above), but
  !     calculated terms for all grid levels at once rather than one at a time.
  !     This was done so that this code could be vectorized and thereby sped up
  !     by the compiler. See clubb:ticket:834 for more information.
  ! 
  !-------------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only:  & 
        gr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)             [-]
      wm_zt,       & ! w wind component on thermodynamic level (k)     [m/s]
      invrs_dzm      ! Inverse of grid spacing (k)                     [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_ac_pr2

    ! Local Variable
    integer :: k


    ! Set lower boundary to 0
    lhs_ac_pr2(1) = 0.0_core_rknd

    do k = 2, gr%nz-1 

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs_ac_pr2(k) = ( 1.0_core_rknd - C7_Skw_fnc(k) ) * invrs_dzm(k) * ( wm_zt(k+1) - wm_zt(k) )

    end do

    lhs_ac_pr2(gr%nz) = 0.0_core_rknd


    return

  end subroutine wpxp_terms_ac_pr2_lhs_all

  !=============================================================================
  pure function wpxp_term_pr1_lhs( C6x_Skw_fnc, tau_C6_zm ) & 
  result( lhs )

    ! Description
    ! Pressure term 1 for w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains pressure term 1:
    !
    ! - ( C_6 / tau_m ) w'x'.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - ( C_6 / tau_m ) w'x'(t+1)
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! The values of w'x' are found on the momentum levels.  The values of the
    ! C_6 skewness function and time-scale tau_m are also found on the momentum
    ! levels.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C6x_Skw_fnc,  & ! C_6x parameter with Sk_w applied (k)                    [-]
      tau_C6_zm       ! Time-scale tau at momentum level (k) applied to C6 term [s]

    ! Return Variable
    real( kind = core_rknd ) :: lhs


    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs = C6x_Skw_fnc / tau_C6_zm


    return
  end function wpxp_term_pr1_lhs

  !=====================================================================================
  pure subroutine wpxp_term_pr1_lhs_all( C6x_Skw_fnc, &
                                         tau_C6_zm, &
                                         lhs_pr1 )
  ! Description:
  !     This subroutine serves the same function as wpxp_term_pr1_lhs (above), but
  !     calculated terms for all grid levels at once rather than one at a time.
  !     This was done so that this code could be vectorized and thereby sped up
  !     by the compiler. See clubb:ticket:834 for more information.
  ! 
  !-------------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only:  & 
        gr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C6x_Skw_fnc,  & ! C_6x parameter with Sk_w applied (k)                    [-]
      tau_C6_zm       ! Time-scale tau at momentum level (k) applied to C6 term [s]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
        lhs_pr1

    ! Local Variable
    integer :: k

    ! Set lower boundary to 0
    lhs_pr1(1) = 0.0_core_rknd

    ! Calculate all other grid levels
    do k = 2, gr%nz-1

        ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
        lhs_pr1(k) = C6x_Skw_fnc(k) / tau_C6_zm(k)
    
    end do

    ! Set upper boundary to 0
    lhs_pr1(gr%nz) = 0.0_core_rknd

    return

  end subroutine wpxp_term_pr1_lhs_all

  !=============================================================================
  pure function wpxp_terms_bp_pr3_rhs( C7_Skw_fnc, thv_ds_zm, xpthvp ) &
  result( rhs )

    ! Description:
    ! Buoyancy production of w'x' and w'x' pressure term 3:  explicit portion of
    ! the code.
    !
    ! The d(w'x')/dt equation contains a buoyancy production term:
    !
    ! + (g/thv_ds) x'th_v';
    !
    ! and pressure term 3:
    !
    ! - C_7 (g/thv_ds) x'th_v'.
    !
    ! Both the w'x' buoyancy production term and pressure term 3 are completely
    ! explicit.  The buoyancy production term and pressure term 3 are combined
    ! and solved together as:
    !
    ! + ( 1 - C_7 ) * (g/thv_ds) * x'th_v'.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: & ! Constants(s) 
        grav, & ! Gravitational acceleration [m/s^2]
        one

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)      [-]
      thv_ds_zm,   & ! Dry, base-state theta_v on mom. lev. (k) [K]
      xpthvp         ! x'th_v'(k)                               [K {xm units}]

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    rhs = ( grav / thv_ds_zm ) * ( one - C7_Skw_fnc ) * xpthvp


    return
  end function wpxp_terms_bp_pr3_rhs

  !=====================================================================================
  pure subroutine wpxp_terms_bp_pr3_rhs_all( C7_Skw_fnc, thv_ds_zm, xpthvp, &
                                             rhs_bp_pr3 )
  ! Description:
  !     This subroutine serves the same function as wpxp_term_pr1_lhs (above), but
  !     calculated terms for all grid levels at once rather than one at a time.
  !     This was done so that this code could be vectorized and thereby sped up
  !     by the compiler. See clubb:ticket:834 for more information.
  ! 
  !-------------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: & ! Constants(s) 
        grav     ! Gravitational acceleration [m/s^2]

    use grid_class, only:  & 
        gr

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied (k)      [-]
      thv_ds_zm,   & ! Dry, base-state theta_v on mom. lev. (k) [K]
      xpthvp         ! x'th_v'(k)                               [K {xm units}]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
        rhs_bp_pr3

    ! Local Variable
    integer :: k
    

    ! Set lower boundary to 0
    rhs_bp_pr3(1) = 0.0_core_rknd

    ! Calculate for all other grid levels
    do k = 2, gr%nz-1
        rhs_bp_pr3(k) = ( grav / thv_ds_zm(k) ) * ( 1.0_core_rknd - C7_Skw_fnc(k) ) * xpthvp(k)
    end do
    
    ! Set upper boundary to 0
    rhs_bp_pr3(gr%nz) = 0.0_core_rknd

    return

  end subroutine wpxp_terms_bp_pr3_rhs_all

  !=============================================================================
  subroutine xm_correction_wpxp_cl( solve_type, dt, wpxp_chnge, invrs_dzt, &
                                    xm )

    ! Description:
    ! Corrects the value of xm if w'x' needed to be clipped, for xm is partially
    ! based on the derivative of w'x' with respect to altitude.
    !
    ! The time-tendency equation for xm is:
    !
    ! d(xm)/dt = -w d(xm)/dz - d(w'x')/dz + d(xm)/dt|_ls;
    !
    ! where d(xm)/dt|_ls is the rate of change of xm over time due to radiation,
    ! microphysics, and/or any other large-scale forcing(s).
    !
    ! The time-tendency equation for xm is solved in conjunction with the
    ! time-tendency equation for w'x'.  Both equations are solved together in a
    ! semi-implicit manner.  However, after both equations have been solved (and
    ! thus both xm and w'x' have been advanced to the next timestep with
    ! timestep index {t+1}), the value of covariance w'x' may be clipped at any
    ! level in order to prevent the correlation of w and x from becoming greater
    ! than 1 or less than -1.
    !
    ! The correlation between w and x is:
    !
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ].
    !
    ! The correlation must always have a value between -1 and 1, such that:
    !
    ! -1 <= corr_(w,x) <= 1.
    !
    ! Therefore, there is an upper limit on w'x', such that:
    !
    ! w'x' <=  [ sqrt(w'^2) * sqrt(x'^2) ];
    !
    ! and a lower limit on w'x', such that:
    !
    ! w'x' >= -[ sqrt(w'^2) * sqrt(x'^2) ].
    !
    ! The aforementioned time-tendency equation for xm is based on the value of
    ! w'x' without being clipped (w'x'{t+1}_unclipped), such that:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz + d(xm{t})/dt|_ls;
    !
    ! where the both the mean advection term, -w d(xm{t+1})/dz, and the
    ! turbulent advection term, -d(w'x'{t+1}_unclipped)/dz, are solved
    ! completely implicitly.  The xm forcing term, +d(xm{t})/dt|_ls, is solved
    ! completely explicitly.
    !
    ! However, if w'x' needs to be clipped after being advanced one timestep,
    ! then xm needs to be altered to reflect the fact that w'x' has a different
    ! value than the value used while both were being solved together.  Ideally,
    ! the xm time-tendency equation that should be used is:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_clipped)/dz + d(xm{t})/dt|_ls.
    !
    ! However, w'x'{t+1}_clipped isn't known until after the w'x' and xm
    ! equations have been solved together.  However, a proper adjuster can be
    ! applied to xm through the use of the following relationship:
    !
    ! w'x'{t+1}_clipped = w'x'{t+1}_unclipped + w'x'{t+1}_amount_clipped;
    !
    ! at any given vertical level.
    !
    ! When the expression above is substituted into the preceeding xm
    ! time-tendency equation, the resulting equation for xm time-tendency is:
    !
    ! d(xm)/dt = -w d(xm{t+1})/dz - d(w'x'{t+1}_unclipped)/dz
    !               - d(w'x'{t+1}_amount_clipped)/dz + d(xm{t})/dt|_ls.
    !
    ! Thus, the resulting xm time-tendency equation is the same as the original
    ! xm time-tendency equation, but with added adjuster term:
    !
    ! -d(w'x'{t+1}_amount_clipped)/dz.
    !
    ! Since the adjuster term needs to be applied after xm has already been
    ! solved, it needs to be multiplied by the timestep length and added on to
    ! xm{t+1}, such that:
    !
    ! xm{t+1}_after_adjustment =
    !    xm{t+1}_before_adjustment + ( -d(w'x'{t+1}_amount_clipped)/dz ) * dt.
    !
    ! The adjuster term is discretized as follows:
    !
    ! The values of w'x' are located on the momentum levels.  Thus, the values
    ! of w'x'_amount_clipped are also located on the momentum levels.  The
    ! values of xm are located on the thermodynamic levels.  The derivatives
    ! (d/dz) of w'x'_amount_clipped are taken over the intermediate
    ! thermodynamic levels, where they are applied to xm.
    !
    ! =======wpxp_amount_clipped=============================== m(k)
    !
    ! -----------------------------d(wpxp_amount_clipped)/dz--- t(k)
    !
    ! =======wpxpm1_amount_clipped============================= m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! Note:  The results of this xm adjustment are highly dependent on the
    !        numerical stability and the smoothness of the w'^2 and x'^2 fields.
    !        An unstable "sawtooth" profile for w'^2 and/or x'^2 causes an
    !        unstable "sawtooth" profile for the upper and lower limits on w'x'.
    !        In turn, this causes an unstable "sawtooth" profile for
    !        w'x'_amount_clipped.  Taking the derivative of that such a "noisy"
    !        field and applying the results to xm causes the xm field to become
    !        more "noisy" and unstable.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr  ! Variable(s); gr%nz only.

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        l_stats_samp, & ! Variable(s)
        stats_zt, &
        ithlm_tacl, &
        irtm_tacl

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      solve_type    ! Variable that is being solved for.

    real( kind = core_rknd ), intent(in) :: &
      dt            ! Model timestep                            [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wpxp_chnge, & ! Amount of change in w'x' due to clipping  [m/s {xm units}]
      invrs_dzt     ! Inverse of grid spacing                   [1/m]

    ! Input/Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      xm            ! xm (thermodynamic levels)                 [{xm units}]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xm_tndcy_wpxp_cl ! d(xm)/dt due to clipping of w'x'       [{xm units}/s]

    integer :: k    ! Array index

    integer :: ixm_tacl  ! Statistical index


    select case ( solve_type )
    case ( xm_wpxp_rtm )
      ixm_tacl = irtm_tacl
    case ( xm_wpxp_thlm )
      ixm_tacl = ithlm_tacl
    case default
      ixm_tacl = 0
    end select

    ! Adjusting xm based on clipping for w'x'.
    ! Loop over all thermodynamic levels between the second-lowest and the
    ! highest.
    do k = 2, gr%nz, 1
      xm_tndcy_wpxp_cl(k) = - invrs_dzt(k) * ( wpxp_chnge(k) - wpxp_chnge(k-1) )
      xm(k) = xm(k) + xm_tndcy_wpxp_cl(k) * dt
    enddo

    if ( l_stats_samp ) then
      ! The adjustment to xm due to turbulent advection term clipping
      ! (xm term tacl) is completely explicit; call stat_update_var.
      call stat_update_var( ixm_tacl, xm_tndcy_wpxp_cl, stats_zt )
    endif


    return

  end subroutine xm_correction_wpxp_cl


  !=============================================================================
  pure function damp_coefficient( coefficient, Cx_Skw_fnc, max_coeff_value, &
                                  threshold, Lscale ) &
    result( damped_value )

    ! Description:
    ! Damps a given coefficient linearly based on the value of Lscale.
    ! For additional information see CLUBB ticket #431.

    use constants_clubb, only: &
        one_hundred  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only: & 
        gr ! Variable(s)

    implicit none

    ! Input variables
    real( kind = core_rknd ), intent(in) :: &
      coefficient,      &   ! The coefficient to be damped
      max_coeff_value,  &   ! Maximum value the damped coefficient should have
      threshold             ! Value of Lscale below which the damping should occur

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Lscale,           &   ! Current value of Lscale
      Cx_Skw_fnc            ! Initial skewness function before damping

    ! Local variables
    real( kind = core_rknd ), parameter :: &
      ! Added to prevent large damping at low altitudes where Lscale is small
      altitude_threshold = one_hundred  ! Altitude above which damping should occur

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: damped_value

    damped_value = Cx_Skw_fnc

    where( Lscale < threshold .and. gr%zt > altitude_threshold)
      damped_value = max_coeff_value &
                     + ( ( coefficient - max_coeff_value ) / threshold ) &
                       * Lscale
    end where

    return

  end function damp_coefficient
  !-----------------------------------------------------------------------

  !=====================================================================================
  pure subroutine diagnose_upxp( ypwp, xm, wpxp, ym, &
                                 C6x_Skw_fnc, tau_C6_zm, C7_Skw_fnc, &
                                 ypxp )
    ! Description:
    !   Diagnose turbulent horizontal flux of a conserved scalar.
    !
    ! References:
    !   Eqn. 7 of Andre et al. (1978)
    !   Eqn. 4 of Bougeault et al. (1981)
    !   github issue #841
    !
    !-------------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: & ! Constants(s)
        one     ! 1.0_core_rknd

    use grid_class, only:  &
        gr, & ! Variable
        ddzt  ! Procedure

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      ypwp,        & ! momentum flux component, either upwp or vpwp  [m^2/s^2]
      xm,          & ! grid-mean conserved thermodynamic variable, either thlm or rtm [varies]
      wpxp,        & ! vertical scalar flux, either wpthlp or wprtp [varies]
      ym,          & ! grid-mean velocity component, either um or vm [m/s]
      C6x_Skw_fnc, & ! C_6 pressure parameter with effects of Sk_w incorporated (k)  [-]
      tau_C6_zm,   & ! Time-scale tau on momentum levels applied to C6 term [s]
      C7_Skw_fnc     ! C_7 pressure parameter with effects of Sk_w incorporated (k)  [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
        ypxp        ! horizontal flux of a conserved scalar, either upthlp, uprtp, vpthlp, or vprtp

    ypxp = ( tau_C6_zm / C6x_Skw_fnc ) * &
              ( - ypwp * ddzt( xm ) - (one - C7_Skw_fnc ) * ( wpxp * ddzt( ym ) ) )

    return

  end subroutine diagnose_upxp

  !=============================================================================
  subroutine error_prints_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                   Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                   Kh_zt, Kh_zm, tau_C6_zm, Skw_zm, &
                                   wp2rtp, rtpthvp, rtm_forcing, &
                                   wprtp_forcing, rtm_ref, wp2thlp, &
                                   thlpthvp, thlm_forcing, wpthlp_forcing, &
                                   thlm_ref, rho_ds_zm, rho_ds_zt, &
                                   invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                   thv_ds_zm, rtp2, thlp2, w_1_zm, w_2_zm, &
                                   varnce_w_1_zm, varnce_w_2_zm, &
                                   mixt_frac_zm, l_implemented, em, &
                                   wp2sclrp, sclrpthvp, sclrm_forcing, &
                                   sclrp2, exner, rcm, p_in_Pa, thvm, &
                                   Cx_fnc_Richardson, &
                                   pdf_implicit_coefs_terms, um_forcing, &
                                   vm_forcing, ug, vg, wpthvp, fcor, &
                                   um_ref, vm_ref, up2, vp2, uprcp, vprcp, &
                                   rc_coef, rtm, wprtp, thlm, wpthlp, &
                                   sclrm, wpsclrp, um, upwp, vm, vpwp, &
                                   rtm_old, wprtp_old, thlm_old, &
                                   wpthlp_old, sclrm_old, wpsclrp_old, &
                                   um_old, upwp_old, vm_old, vpwp_old )

    ! Description:
    ! Prints values of model fields when fatal errors (LU decomp.) occur.
    ! All field that are passed into and out of subroutine advance_xm_wpxp are
    ! printed.  If additional fields are added to the call to subroutine
    ! advance_xm_wpxp, they should also be added here.

    use constants_clubb, only: &
        fstderr    ! Variable(s)

    use grid_class, only: &
        gr    ! Variable Type(s)

    use model_flags, only: &
        l_predict_upwp_vpwp    ! Variable(s)

    use parameters_model, only: &
        sclr_dim    ! Variable(s)

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      sigma_sqd_w,     & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,           & ! w wind component on momentum levels      [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      Lscale,          & ! Turbulent mixing length                  [m]
      em,              & ! Turbulent Kinetic Energy (TKE)           [m^2/s^2]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on momentum levels    [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels     [m/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      Kh_zm,           & ! Eddy diffusivity on momentum levels
      tau_C6_zm,       & ! Time-scale tau on momentum levels applied to C6 term [s]
      Skw_zm,          & ! Skewness of w on momentum levels         [-]
      wp2rtp,          & ! <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
      rtm_ref,         & ! rtm for nudging                          [kg/kg]
      wp2thlp,         & ! <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
      thlpthvp,        & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,    & ! th_l forcing (thermodynamic levels)      [K/s]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)      [K/s^2]
      thlm_ref,        & ! thlm for nudging                         [K]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on moment. levs. [K]
      ! Added for clipping by Vince Larson 29 Sep 2007
      rtp2,            & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2,           & ! th_l'^2 (momentum levels)                [K^2]
      ! End of Vince Larson's addition.
      w_1_zm,          & ! Mean w (1st PDF component)              [m/s]
      w_2_zm,          & ! Mean w (2nd PDF component)              [m/s]
      varnce_w_1_zm,   & ! Variance of w (1st PDF component)       [m^2/s^2]
      varnce_w_2_zm,   & ! Variance of w (2nd PDF component)       [m^2/s^2]
      mixt_frac_zm      ! Weight of 1st PDF component (Sk_w dependent) [-]

    logical, intent(in) ::  & 
      l_implemented      ! Flag for CLUBB being implemented in a larger model.

    ! Additional variables for passive scalars
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: & 
      wp2sclrp,      & ! <w'^2 sclr'> (thermodynamic levels)   [Units vary]
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      exner,           & ! Exner function                            [-]
      rcm,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,         & ! Air pressure                              [Pa]
      thvm,            & ! Virutal potential temperature             [K]
      Cx_fnc_Richardson  ! Cx_fnc computed from Richardson_num       [-]

    type(implicit_coefs_terms), dimension(gr%nz), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      um_forcing, & ! <u> forcing term (thermodynamic levels)      [m/s^2]
      vm_forcing, & ! <v> forcing term (thermodynamic levels)      [m/s^2]
      ug,         & ! <u> geostrophic wind (thermodynamic levels)  [m/s]
      vg,         & ! <v> geostrophic wind (thermodynamic levels)  [m/s]
      wpthvp        ! <w'thv'> (momentum levels)                   [m/s K]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      uprcp,   & ! < u' r_c' >                                  [(m kg)/(s kg)]
      vprcp,   & ! < v' r_c' >                                  [(m kg)/(s kg)]
      rc_coef    ! Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]

     real( kind = core_rknd ), intent(in) ::  &
      fcor          ! Coriolis parameter                           [s^-1]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      um_ref, & ! Reference u wind component for nudging       [m/s]
      vm_ref, & ! Reference v wind component for nudging       [m/s]
      up2,    & ! Variance of the u wind component             [m^2/s^2]
      vp2       ! Variance of the v wind component             [m^2/s^2]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      um,   & ! <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
      upwp, & ! <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      vm,   & ! <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
      vpwp    ! <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]

    ! Saved values of predictive fields, prior to being advanced, for use in
    ! print statements in case of fatal error.
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      rtm_old,    & ! Saved value of r_t        [kg/kg]
      wprtp_old,  & ! Saved value of w'r_t'     [(kg/kg) m/s]
      thlm_old,   & ! Saved value of th_l       [K]
      wpthlp_old    ! Saved value of w'th_l'    [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(in) ::  & 
      sclrm_old, wpsclrp_old !                  [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      um_old,   & ! Saved value of <u>       [m/s]
      upwp_old, & ! Saved value of <u'w'>    [m^2/s^2]
      vm_old,   & ! Saved value of <v>       [m/s]
      vpwp_old    ! Saved value of <v'w'>    [m^2/s^2]


    write(fstderr,*) "Error in advance_xm_wpxp", new_line('c')

    write(fstderr,*) "Intent(in)", new_line('c')

    write(fstderr,*) "gr%zt = ", gr%zt, new_line('c')
    write(fstderr,*) "gr%zm = ", gr%zm, new_line('c')
    write(fstderr,*) "dt = ", dt, new_line('c')
    write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w, new_line('c')
    write(fstderr,*) "wm_zm = ", wm_zm, new_line('c')
    write(fstderr,*) "wm_zt = ", wm_zt, new_line('c')
    write(fstderr,*) "wp2 = ", wp2, new_line('c')
    write(fstderr,*) "Lscale = ", Lscale, new_line('c')
    write(fstderr,*) "wp3_on_wp2 = ", wp3_on_wp2, new_line('c')
    write(fstderr,*) "wp3_on_wp2_zt = ", wp3_on_wp2_zt, new_line('c')
    write(fstderr,*) "Kh_zt = ", Kh_zt, new_line('c')
    write(fstderr,*) "Kh_zm = ", Kh_zm, new_line('c')
    write(fstderr,*) "tau_C6_zm = ", tau_C6_zm, new_line('c')
    write(fstderr,*) "Skw_zm = ", Skw_zm, new_line('c')
    write(fstderr,*) "wp2rtp = ", wp2rtp, new_line('c')
    write(fstderr,*) "rtpthvp = ", rtpthvp, new_line('c')
    write(fstderr,*) "rtm_forcing = ", rtm_forcing, new_line('c')
    write(fstderr,*) "wprtp_forcing = ", wprtp_forcing, new_line('c')
    write(fstderr,*) "rtm_ref = ", rtm_ref, new_line('c')
    write(fstderr,*) "wp2thlp = ", wp2thlp, new_line('c')
    write(fstderr,*) "thlpthvp = ", thlpthvp, new_line('c')
    write(fstderr,*) "thlm_forcing = ", thlm_forcing, new_line('c')
    write(fstderr,*) "wpthlp_forcing = ", wpthlp_forcing, new_line('c')
    write(fstderr,*) "thlm_ref = ", thlm_ref, new_line('c')
    write(fstderr,*) "rho_ds_zm = ", rho_ds_zm, new_line('c')
    write(fstderr,*) "rho_ds_zt = ", rho_ds_zt, new_line('c')
    write(fstderr,*) "invrs_rho_ds_zm = ", invrs_rho_ds_zm, new_line('c')
    write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt, new_line('c')
    write(fstderr,*) "thv_ds_zm = ", thv_ds_zm, new_line('c')
    write(fstderr,*) "rtp2 = ", rtp2, new_line('c')
    write(fstderr,*) "thlp2 = ", thlp2, new_line('c')
    write(fstderr,*) "w_1_zm = ", w_1_zm, new_line('c')
    write(fstderr,*) "w_2_zm = ", w_2_zm, new_line('c')
    write(fstderr,*) "varnce_w_1_zm = ", varnce_w_1_zm, new_line('c')
    write(fstderr,*) "varnce_w_2_zm = ", varnce_w_2_zm, new_line('c')
    write(fstderr,*) "mixt_frac_zm = ", mixt_frac_zm, new_line('c')
    write(fstderr,*) "l_implemented = ", l_implemented, new_line('c')
    write(fstderr,*) "em = ", em, new_line('c')
    write(fstderr,*) "exner = ", exner, new_line('c')
    write(fstderr,*) "rcm = ", rcm, new_line('c')
    write(fstderr,*) "p_in_Pa = ", p_in_Pa, new_line('c')
    write(fstderr,*) "thvm = ", thvm, new_line('c')
    write(fstderr,*) "Cx_fnc_Richardson = ", Cx_fnc_Richardson, new_line('c')
    write(fstderr,*) "pdf_implicit_coefs_terms = ", pdf_implicit_coefs_terms, &
                     new_line('c')
     
    if ( sclr_dim > 0 )  then
       write(fstderr,*) "sclrp2 = ", sclrp2, new_line('c')
       write(fstderr,*) "wp2sclrp = ", wp2sclrp, new_line('c')
       write(fstderr,*) "sclrpthvp = ", sclrpthvp, new_line('c')
       write(fstderr,*) "sclrm_forcing = ", sclrm_forcing, new_line('c')
    endif

    if ( l_predict_upwp_vpwp ) then
       write(fstderr,*) "um_forcing = ", um_forcing, new_line('c')
       write(fstderr,*) "vm_forcing = ", vm_forcing, new_line('c')
       write(fstderr,*) "ug = ", ug, new_line('c')
       write(fstderr,*) "vg = ", vg, new_line('c')
       write(fstderr,*) "wpthvp = ", wpthvp, new_line('c')
       write(fstderr,*) "fcor = ", fcor, new_line('c')
       write(fstderr,*) "um_ref = ", um_ref, new_line('c')
       write(fstderr,*) "vm_ref = ", vm_ref, new_line('c')
       write(fstderr,*) "up2 = ", up2, new_line('c')
       write(fstderr,*) "vp2 = ", vp2, new_line('c')
       write(fstderr,*) "uprcp = ", uprcp, new_line('c')
       write(fstderr,*) "vprcp = ", vprcp, new_line('c')
       write(fstderr,*) "rc_coef = ",  rc_coef, new_line('c')
    endif ! l_predict_upwp_vpwp

    write(fstderr,*) "Intent(inout)", new_line('c')
     
    write(fstderr,*) "rtm (pre-solve) = ", rtm_old, new_line('c')
    write(fstderr,*) "rtm = ", rtm, new_line('c')
    write(fstderr,*) "wprtp (pre-solve) = ", wprtp_old, new_line('c')
    write(fstderr,*) "wprtp = ", wprtp, new_line('c')
    write(fstderr,*) "thlm (pre-solve) = ", thlm_old, new_line('c')
    write(fstderr,*) "thlm = ", thlm, new_line('c')
    write(fstderr,*) "wpthlp (pre-solve) =", wpthlp_old, new_line('c')
    write(fstderr,*) "wpthlp =", wpthlp, new_line('c')

    if ( sclr_dim > 0 )  then
       write(fstderr,*) "sclrm (pre-solve) = ", sclrm_old, new_line('c')
       write(fstderr,*) "sclrm = ", sclrm, new_line('c')
       write(fstderr,*) "wpsclrp (pre-solve) = ", wpsclrp_old, new_line('c')
       write(fstderr,*) "wpsclrp = ", wpsclrp, new_line('c')
    endif

    if ( l_predict_upwp_vpwp ) then
       write(fstderr,*) "um (pre-solve) = ", um_old, new_line('c')
       write(fstderr,*) "um = ", um, new_line('c')
       write(fstderr,*) "upwp (pre-solve) = ",  upwp_old, new_line('c')
       write(fstderr,*) "upwp = ",  upwp, new_line('c')
       write(fstderr,*) "vm (pre-solve) = ", vm_old, new_line('c')
       write(fstderr,*) "vm = ", vm, new_line('c')
       write(fstderr,*) "vpwp (pre-solve) = ",  vpwp_old, new_line('c')
       write(fstderr,*) "vpwp = ",  vpwp, new_line('c')
    endif ! l_predict_upwp_vpwp


    return

  end subroutine error_prints_xm_wpxp

  !=============================================================================

end module advance_xm_wpxp_module
