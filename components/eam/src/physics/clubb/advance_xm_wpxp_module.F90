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
             wpxp_term_tp_lhs, & 
             wpxp_terms_ac_pr2_lhs, & 
             wpxp_term_pr1_lhs, & 
             wpxp_terms_bp_pr3_rhs, &
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
  subroutine advance_xm_wpxp( nz, ngrdcol, gr, dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                              Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, &
                              invrs_tau_C6_zm, tau_max_zm, Skw_zm, wp2rtp, rtpthvp, &
                              rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp, &
                              thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref, &
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                              invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
                              w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
                              mixt_frac_zm, l_implemented, em, wp2sclrp, &
                              sclrpthvp, sclrm_forcing, sclrp2, exner, rcm, &
                              p_in_Pa, thvm, Cx_fnc_Richardson, &
                              ice_supersat_frac, &
                              pdf_implicit_coefs_terms, &
                              um_forcing, vm_forcing, ug, vg, wpthvp, &
                              fcor, um_ref, vm_ref, up2, vp2, &
                              uprcp, vprcp, rc_coef, &
                              clubb_params, nu_vert_res_dep, &
                              iiPDF_type, &
                              l_predict_upwp_vpwp, &
                              l_diffuse_rtm_and_thlm, &
                              l_stability_correct_Kh_N2_zm, &
                              l_godunov_upwind_wpxp_ta, &
                              l_upwind_xm_ma, &
                              l_uv_nudge, &
                              l_tke_aniso, &
                              l_diag_Lscale_from_tau, &
                              l_use_C7_Richardson, &
                              l_brunt_vaisala_freq_moist, &
                              l_use_thvm_in_bv_freq, &
                              l_lmm_stepping, &
                              l_enable_relaxed_clipping, &
                              order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, &
                              stats_zt, stats_zm, stats_sfc, &
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

    use parameter_indices, only: &
        nparams,             & ! Variable(s)
        iC6rt,               &
        iC6rtb,              &
        iC6rtc,              &
        iC6thl,              &
        iC6thlb,             &
        iC6thlc,             &
        iC6rt_Lscale0,       &
        iC6thl_Lscale0,      &
        iC7,                 &
        iC7b,                &
        iC7c,                &
        iC7_Lscale0,         &
        ic_K6,               &
        iwpxp_L_thresh,      &
        ialtitude_threshold, &
        iC_uu_shr

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use constants_clubb, only:  & 
        fstderr, &  ! Constant
        one, &
        one_half, &
        zero, &
        eps

    use parameters_model, only: & 
        sclr_dim, &  ! Variable(s)
        ts_nudge

    use grid_class, only: & 
        grid, & ! Type
        ddzt    ! Procedure(s)

    use grid_class, only: &
        zm2zt, & ! Procedure(s)
        zt2zm

    use model_flags, only: &
        iiPDF_new,                     & ! Variable(s)
        l_explicit_turbulent_adv_wpxp

    use mono_flux_limiter, only: &
        calc_turb_adv_range ! Procedure(s)

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
        irtm_sdmp, &
        ithlm_sdmp, & 
        ium_sdmp, &
        ivm_sdmp, &
        ium_ndg, &
        ivm_ndg, &
        ium_ref, &
        ivm_ref, &
        iC7_Skw_fnc, &
        iC6rt_Skw_fnc, &
        iC6thl_Skw_fnc, &
        iC6_term, &
        l_stats_samp

    use sponge_layer_damping, only: &
        rtm_sponge_damp_settings, &
        thlm_sponge_damp_settings, &
        uv_sponge_damp_settings, &
        rtm_sponge_damp_profile, &
        thlm_sponge_damp_profile, &
        uv_sponge_damp_profile, &
        sponge_damp_xm ! Procedure(s)

    use stats_type, only: stats ! Type

    implicit none

    ! -------------------- Input Variables --------------------
    
    integer, intent(in) :: &
      nz, &
      ngrdcol
      
    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: & 
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
      invrs_tau_C6_zm, & ! Inverse time-scale on mom. levels applied to C6 term [1/s]
      tau_max_zm,      & ! Max. allowable eddy dissipation time scale on m-levs  [s]
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
      mixt_frac_zm       ! Weight of 1st PDF component (Sk_w dependent) [-]

    logical, intent(in) ::  & 
      l_implemented      ! Flag for CLUBB being implemented in a larger model.

    ! Additional variables for passive scalars
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz,sclr_dim) :: & 
      wp2sclrp,      & ! <w'^2 sclr'> (thermodynamic levels)   [Units vary]
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) ::  &
      exner,            & ! Exner function                            [-]
      rcm,              & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,          & ! Air pressure                              [Pa]
      thvm,             & ! Virutal potential temperature             [K]
      Cx_fnc_Richardson,& ! Cx_fnc computed from Richardson_num       [-]
      ice_supersat_frac

    type(implicit_coefs_terms), dimension(ngrdcol), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      um_forcing, & ! <u> forcing term (thermodynamic levels)      [m/s^2]
      vm_forcing, & ! <v> forcing term (thermodynamic levels)      [m/s^2]
      ug,         & ! <u> geostrophic wind (thermodynamic levels)  [m/s]
      vg,         & ! <v> geostrophic wind (thermodynamic levels)  [m/s]
      wpthvp        ! <w'thv'> (momentum levels)                   [m/s K]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      uprcp,              & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp,              & ! < v' r_c' >              [(m kg)/(s kg)]
      rc_coef               ! Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]

     real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      fcor          ! Coriolis parameter                           [s^-1]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      um_ref, & ! Reference u wind component for nudging       [m/s]
      vm_ref, & ! Reference v wind component for nudging       [m/s]
      up2,    & ! Variance of the u wind component             [m^2/s^2]
      vp2       ! Variance of the v wind component             [m^2/s^2]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), dimension(ngrdcol), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diffuse_rtm_and_thlm,       & ! This flag determines whether or not we want CLUBB to do
                                      ! diffusion on rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! This flag determines whether or not we want CLUBB to apply
                                      ! a stability correction
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms.
                                      ! It affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                                      ! (u'^2 + v'^2 + w'^2)
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_enable_relaxed_clipping       ! Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats

    integer, intent(in) :: &
      order_xm_wpxp, &
      order_xp2_xpyp, &
      order_wp2_wp3

    ! -------------------- Input/Output Variables --------------------

    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  & 
      um,   & ! <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
      upwp, & ! <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      vm,   & ! <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
      vpwp    ! <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
 
    ! -------------------- Local Variables --------------------

    ! Parameter Constants
    logical, parameter :: &
      l_iter = .true. ! True when the means and fluxes are prognosed

    real( kind = core_rknd ) ::  &
      C6rt,               & ! CLUBB tunable parameter C6rt
      C6rtb,              & ! CLUBB tunable parameter C6rtb
      C6rtc,              & ! CLUBB tunable parameter C6rtc
      C6thl,              & ! CLUBB tunable parameter C6thl
      C6thlb,             & ! CLUBB tunable parameter C6thlb
      C6thlc,             & ! CLUBB tunable parameter C6thlc
      C6rt_Lscale0,       & ! CLUBB tunable parameter C6rt_Lscale0
      C6thl_Lscale0,      & ! CLUBB tunable parameter C6thl_Lscale0
      C7,                 & ! CLUBB tunable parameter C7
      C7b,                & ! CLUBB tunable parameter C7b
      C7c,                & ! CLUBB tunable parameter C7c
      C7_Lscale0,         & ! CLUBB tunable parameter C7_Lscale0
      c_K6,               & ! CLUBB tunable parameter c_K6
      altitude_threshold, & ! CLUBB tunable parameter altitude_threshold
      wpxp_L_thresh         ! CLUBB tunable parameter wpxp_L_thresh

    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc, C6_term

    ! Eddy Diffusion for wpthlp and wprtp.
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: Kw6  ! wpxp eddy diff. [m^2/s]

    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(ngrdcol,nz) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.

    ! Constant parameters as a function of Skw.

    integer :: &
      nrhs         ! Number of RHS vectors

    ! Saved values of predictive fields, prior to being advanced, for use in
    ! print statements in case of fatal error.
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      rtm_old,    & ! Saved value of r_t        [kg/kg]
      wprtp_old,  & ! Saved value of w'r_t'     [(kg/kg) m/s]
      thlm_old,   & ! Saved value of th_l       [K]
      wpthlp_old    ! Saved value of w'th_l'    [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) ::  & 
      sclrm_old,   & ! Saved value of sclr      [units vary]
      wpsclrp_old    ! Saved value of wpsclrp   [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(ngrdcol,nz) ::  & 
      um_old,   & ! Saved value of <u>       [m/s]
      upwp_old, & ! Saved value of <u'w'>    [m^2/s^2]
      vm_old,   & ! Saved value of <v>       [m/s]
      vpwp_old    ! Saved value of <v'w'>    [m^2/s^2]
      
    ! LHS/RHS terms
    real( kind = core_rknd ), dimension(3,ngrdcol,nz) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for w'x'
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm       ! Mean advection contributions to lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz) :: & 
      lhs_ta_wprtp,  & ! w'r_t' turbulent advection contributions to lhs  
      lhs_ta_wpthlp, & ! w'thl' turbulent advection contributions to lhs
      lhs_ta_wpup,   & ! w'u' turbulent advection contributions to lhs
      lhs_ta_wpvp      ! w'v' turbulent advection contributions to lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz,sclr_dim) :: & 
      lhs_ta_wpsclrp    ! w'sclr' turbulent advection contributions to lhs
     
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      rhs_ta_wprtp,  & ! w'r_t' turbulent advection contributions to rhs  
      rhs_ta_wpthlp, & ! w'thl' turbulent advection contributions to rhs
      rhs_ta_wpup,   & ! w'u' turbulent advection contributions to rhs
      rhs_ta_wpvp      ! w'v' turbulent advection contributions to rhs
      
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: & 
      rhs_ta_wpsclrp    ! w'sclr' turbulent advection contributions to rhs

    real( kind = core_rknd ), dimension(2,ngrdcol,nz) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
    
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      lhs_ac_pr2,     & ! Accumulation of w'x' and w'x' pressure term 2
      lhs_pr1_wprtp,  & ! Pressure term 1 for w'r_t' for all grid levels
      lhs_pr1_wpthlp, & ! Pressure term 1 for w'thl' for all grid levels
      lhs_pr1_wpsclrp   ! Pressure term 1 for w'sclr' for all grid levels
      
    logical :: &
      l_scalar_calc   ! True if sclr_dim > 0
      
    integer :: i, j, k

    ! -------------------- Begin Code --------------------
    
    ! Check whether the passive scalars are present.
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if
      
    if ( iiPDF_type == iiPDF_new .and. ( .not. l_explicit_turbulent_adv_wpxp ) ) then
       nrhs = 1
    else
      if ( l_predict_upwp_vpwp ) then
        nrhs = 4+sclr_dim
      else
        nrhs = 2+sclr_dim
      endif
    end if

    ! Save values of predictive fields to be printed in case of crash.
    if ( l_lmm_stepping ) then
      
      rtm_old(:,:)    = rtm(:,:)
      wprtp_old(:,:)  = wprtp(:,:)
      thlm_old(:,:)   = thlm(:,:)
      wpthlp_old(:,:) = wpthlp(:,:)
          
      if ( sclr_dim > 0 ) then
        sclrm_old(:,:,:) = sclrm(:,:,:)
        wpsclrp_old(:,:,:) = wpsclrp(:,:,:)
      end if ! sclr_dim > 0
       
      if ( l_predict_upwp_vpwp ) then
        um_old(:,:) = um(:,:)
        upwp_old(:,:) = upwp(:,:)
        vm_old(:,:) = vm(:,:)
        vpwp_old(:,:) = vpwp(:,:)
      end if ! l_predict_upwp_vpwp
       
    end if ! l_lmm_stepping

    ! Unpack CLUBB tunable parameters
    C6rt = clubb_params(iC6rt)
    C6thl = clubb_params(iC6thl)
    altitude_threshold = clubb_params(ialtitude_threshold)
    wpxp_L_thresh = clubb_params(iwpxp_L_thresh)

    if ( .not. l_diag_Lscale_from_tau ) then

      ! Unpack CLUBB tunable parameters
      C6rtb = clubb_params(iC6rtb)
      C6rtc = clubb_params(iC6rtc)
      C6thlb = clubb_params(iC6thlb)
      C6thlc = clubb_params(iC6thlc)
      C6rt_Lscale0 = clubb_params(iC6rt_Lscale0)
      C6thl_Lscale0 = clubb_params(iC6thl_Lscale0)

      ! Compute C6 as a function of Skw
      ! The if...then is just here to save compute time
      if ( abs(C6rt-C6rtb) > abs(C6rt+C6rtb)*eps/2 ) then
        C6rt_Skw_fnc(:,:) = C6rtb + ( C6rt - C6rtb ) & 
                                    * exp( -one_half * (Skw_zm(:,:)/C6rtc)**2 )
      else
        C6rt_Skw_fnc(:,:) = C6rtb
      end if

      if ( abs(C6thl-C6thlb) > abs(C6thl+C6thlb)*eps/2 ) then
        C6thl_Skw_fnc(:,:) = C6thlb + ( C6thl - C6thlb ) & 
                                      * exp( -one_half * (Skw_zm(:,:)/C6thlc)**2 )
      else
        C6thl_Skw_fnc(:,:) = C6thlb
      end if

      ! Damp C6 as a function of Lscale in stably stratified regions
      call damp_coefficient( nz, ngrdcol, gr, C6rt, C6rt_Skw_fnc, &
                             C6rt_Lscale0, altitude_threshold, &
                             wpxp_L_thresh, Lscale, &
                             C6rt_Skw_fnc )

      call damp_coefficient( nz, ngrdcol, gr, C6thl, C6thl_Skw_fnc, &
                             C6thl_Lscale0, altitude_threshold, &
                             wpxp_L_thresh, Lscale, &
                             C6thl_Skw_fnc )

    else ! l_diag_Lscale_from_tau

      C6rt_Skw_fnc(:,:) = C6rt
      C6thl_Skw_fnc(:,:) = C6thl

    endif ! .not. l_diag_Lscale_from_tau

    ! Compute C7_Skw_fnc
    if ( l_use_C7_Richardson ) then

      ! New formulation based on Richardson number
      C7_Skw_fnc(:,:) = Cx_fnc_Richardson(:,:)

    else

      ! Unpack CLUBB tunable parameters
      C7 = clubb_params(iC7)
      C7b = clubb_params(iC7b)
      C7c = clubb_params(iC7c)
      C7_Lscale0 = clubb_params(iC7_Lscale0)

      ! Compute C7 as a function of Skw
      if ( abs(C7-C7b) > abs(C7+C7b)*eps/2 ) then
        C7_Skw_fnc(:,:) = C7b + ( C7 - C7b ) * exp( -one_half * (Skw_zm(:,:)/C7c)**2 )
      else
        C7_Skw_fnc(:,:) = C7b
      endif

      ! Damp C7 as a function of Lscale in stably stratified regions
      call damp_coefficient( nz, ngrdcol, gr, C7, C7_Skw_fnc, &
                             C7_Lscale0, altitude_threshold, &
                             wpxp_L_thresh, Lscale, &
                             C7_Skw_fnc )

    end if ! l_use_C7_Richardson
    
    
    if ( l_stats_samp ) then

      do i = 1, ngrdcol
        call stat_update_var( iC7_Skw_fnc, C7_Skw_fnc(i,:), & ! intent(in)
                              stats_zm(i) )                 ! intent(inout)
        call stat_update_var( iC6rt_Skw_fnc, C6rt_Skw_fnc(i,:), & ! intent(in)
                              stats_zm(i) )                     ! intent(inout
        call stat_update_var( iC6thl_Skw_fnc, C6thl_Skw_fnc(i,:), & ! intent(in)
                              stats_zm(i) )                       ! intent(inout)
      end do

    end if

    if ( clubb_at_least_debug_level( 0 ) ) then
      ! Assertion check for C7_Skw_fnc
      if ( any( C7_Skw_fnc(:,:) > one ) .or. any( C7_Skw_fnc(:,:) < zero ) ) then
        write(fstderr,*) "The C7_Skw_fnc variable is outside the valid range"
        err_code = clubb_fatal_error
        return
      end if
    end if

    ! Define the Coefficent of Eddy Diffusivity for the wpthlp and wprtp.
    ! Kw6 is used for wpthlp and wprtp, which are located on momentum levels.
    ! Kw6 is located on thermodynamic levels.
    ! Kw6 = c_K6 * Kh_zt
    c_K6 = clubb_params(ic_K6)
    Kw6(:,:) = c_K6 * Kh_zt(:,:)

    ! Find the number of grid levels, both upwards and downwards, that can
    ! have an effect on the central thermodynamic level during the course of
    ! one time step due to turbulent advection.  This is used as part of the
    ! monotonic turbulent advection scheme.
    call calc_turb_adv_range( nz, ngrdcol, gr, dt, &
                              w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, & ! intent(in)
                              mixt_frac_zm, &  ! intent(in)
                              stats_zm, & ! intent(inout)
                              low_lev_effect, high_lev_effect ) ! intent(out)

    
    ! Calculate 1st pressure terms for w'r_t', w'thl', and w'sclr'. 
    call wpxp_term_pr1_lhs( nz, ngrdcol, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc, & ! Intent(in)
                            invrs_tau_C6_zm, l_scalar_calc,                       & ! Intent(in)
                            lhs_pr1_wprtp, lhs_pr1_wpthlp, lhs_pr1_wpsclrp )        ! Intent(out)
    
    C6_term(:,:) = C6rt_Skw_fnc(:,:) * invrs_tau_C6_zm(:,:)

    if ( l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_update_var( iC6_term, C6_term(i,:), & ! intent(in)
                              stats_zm(i) )           ! intent(inout)
      end do
    end if

    call  calc_xm_wpxp_ta_terms( nz, ngrdcol, gr, wp2rtp, &  ! intent(in)
                                 wp2thlp, wp2sclrp, & ! intent(in)
                                 rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm, & ! intent(in)
                                 sigma_sqd_w, wp3_on_wp2_zt, & ! intent(in)
                                 pdf_implicit_coefs_terms, & ! intent(in)
                                 iiPDF_type, & ! intent(in)
                                 l_explicit_turbulent_adv_wpxp, l_predict_upwp_vpwp, & ! intent(in)
                                 l_scalar_calc, & ! intent(in)
                                 l_godunov_upwind_wpxp_ta, & ! intent(in)
                                 stats_zt, & ! intent(inout)
                                 lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpup, & ! intent(out)
                                 lhs_ta_wpvp, lhs_ta_wpsclrp, & ! intent(out)
                                 rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup, & ! intent(out)
                                 rhs_ta_wpvp, rhs_ta_wpsclrp ) ! intent(out)

    ! Calculate various terms that are the same between all LHS matricies
    call calc_xm_wpxp_lhs_terms( nz, ngrdcol, gr, Kh_zm, wm_zm, wm_zt, wp2,        & ! In
                                 Kw6, C7_Skw_fnc, invrs_rho_ds_zt,                 & ! In
                                 invrs_rho_ds_zm, rho_ds_zt,                       & ! In
                                 rho_ds_zm, l_implemented, em,                     & ! In
                                 Lscale, thlm, exner, rtm, rcm, p_in_Pa, thvm,     & ! In
                                 ice_supersat_frac,                                & ! In
                                 clubb_params, nu_vert_res_dep,                    & ! In
                                 l_diffuse_rtm_and_thlm,                           & ! In
                                 l_stability_correct_Kh_N2_zm,                     & ! In
                                 l_upwind_xm_ma,                                   & ! In
                                 l_brunt_vaisala_freq_moist,                       & ! In
                                 l_use_thvm_in_bv_freq,                            & ! In
                                 lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,   & ! Out
                                 lhs_tp, lhs_ta_xm, lhs_ac_pr2 ) ! Out

    ! Setup and decompose matrix for each variable.

    if ( ( iiPDF_type == iiPDF_new ) .and. ( .not. l_explicit_turbulent_adv_wpxp ) ) then

      ! LHS matrices are unique, multiple band solves required
      call solve_xm_wpxp_with_multiple_lhs( nz, ngrdcol, gr, dt, l_iter, nrhs, wm_zt, wp2,  & ! In
                                            rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,  & ! In
                                            thlm_forcing,   wpthlp_forcing, rho_ds_zm,      & ! In
                                            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,    & ! In
                                            thv_ds_zm, rtp2, thlp2, l_implemented,          & ! In
                                            sclrpthvp, sclrm_forcing, sclrp2,               & ! In
                                            low_lev_effect, high_lev_effect, C7_Skw_fnc,    & ! In
                                            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, & ! In
                                            lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpsclrp,    & ! In
                                            rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpsclrp,    & ! In
                                            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,   & ! In
                                            lhs_pr1_wpthlp, lhs_pr1_wpsclrp,                & ! In
                                            l_predict_upwp_vpwp,                            & ! In
                                            l_diffuse_rtm_and_thlm,                         & ! In
                                            l_upwind_xm_ma,                                 & ! In
                                            l_tke_aniso,                                    & ! In
                                            l_enable_relaxed_clipping,                      & ! In
                                            order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3,   & ! In
                                            stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                                            rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp )        ! Out
    else
        
      ! LHS matrices are equivalent, only one solve required
      call solve_xm_wpxp_with_single_lhs( nz, ngrdcol, gr, dt, l_iter, nrhs, wm_zt, wp2,   & ! In 
                                          invrs_tau_C6_zm, tau_max_zm,                     & ! In
                                          rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp,   & ! In
                                          thlm_forcing, wpthlp_forcing, rho_ds_zm,         & ! In
                                          rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,     & ! In
                                          thv_ds_zm, rtp2, thlp2, l_implemented,           & ! In
                                          sclrpthvp, sclrm_forcing, sclrp2, um_forcing,    & ! In
                                          vm_forcing, ug, vg, uprcp, vprcp, rc_coef, fcor, & ! In
                                          up2, vp2,                                        & ! In
                                          low_lev_effect, high_lev_effect,                 & ! In
                                          C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc,         & ! In
                                          lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,  & ! In
                                          lhs_ta_wprtp,                                    & ! In
                                          rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup,        & ! In
                                          rhs_ta_wpvp, rhs_ta_wpsclrp,                     & ! In
                                          lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp,    & ! In
                                          lhs_pr1_wpthlp, lhs_pr1_wpsclrp,                 & ! In
                                          clubb_params(iC_uu_shr),                         & ! In
                                          l_predict_upwp_vpwp,                             & ! In
                                          l_diffuse_rtm_and_thlm,                          & ! In
                                          l_upwind_xm_ma,                                  & ! In
                                          l_tke_aniso,                                     & ! In
                                          l_enable_relaxed_clipping,                       & ! In
                                          order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3,    & ! In
                                          stats_zt, stats_zm, stats_sfc,                   & ! In
                                          rtm, wprtp, thlm, wpthlp,                        & ! Out
                                          sclrm, wpsclrp, um, upwp, vm,vpwp )                ! Out
    end if ! ( ( iiPDF_type == iiPDF_new ) .and. ( .not. l_explicit_turbulent_adv_wpxp ) )

    if ( l_lmm_stepping ) then
      
      thlm(:,:) = one_half * ( thlm_old(:,:) + thlm(:,:) )
      rtm(:,:) = one_half * ( rtm_old(:,:) + rtm(:,:) )
      wpthlp(:,:) = one_half * ( wpthlp_old(:,:) + wpthlp(:,:) ) 
      wprtp(:,:) = one_half * ( wprtp_old(:,:) + wprtp(:,:) )
      
      if ( sclr_dim > 0 ) then
         sclrm(:,:,:) = one_half * ( sclrm_old(:,:,:) + sclrm(:,:,:) )
         wpsclrp(:,:,:) = one_half * ( wpsclrp_old(:,:,:) + wpsclrp(:,:,:) )
      endif ! sclr_dim > 0
      
      if ( l_predict_upwp_vpwp ) then
        um(:,:) = one_half * ( um_old(:,:) + um(:,:) )
        vm(:,:) = one_half * ( vm_old(:,:) + vm(:,:) )
        upwp(:,:) = one_half * ( upwp_old(:,:) + upwp(:,:) )
        vpwp(:,:) = one_half * ( vpwp_old(:,:) + vpwp(:,:) )    
      end if ! l_predict_upwp_vpwp 
      
    end if ! l_lmm_stepping

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        do i = 1, ngrdcol
          call error_prints_xm_wpxp( gr(i), dt, sigma_sqd_w(i,:), wm_zm(i,:), wm_zt(i,:), wp2(i,:), & ! intent(in)
                                     Lscale(i,:), wp3_on_wp2(i,:), wp3_on_wp2_zt(i,:), & ! intent(in)
                                     Kh_zt(i,:), Kh_zm(i,:), invrs_tau_C6_zm(i,:), Skw_zm(i,:), & ! intent(in)
                                     wp2rtp(i,:), rtpthvp(i,:), rtm_forcing(i,:), & ! intent(in)
                                     wprtp_forcing(i,:), rtm_ref(i,:), wp2thlp(i,:), & ! intent(in)
                                     thlpthvp(i,:), thlm_forcing(i,:), & ! intent(in)
                                     wpthlp_forcing(i,:), thlm_ref(i,:), rho_ds_zm(i,:), & ! intent(in)
                                     rho_ds_zt(i,:), invrs_rho_ds_zm(i,:), & ! intent(in)
                                     invrs_rho_ds_zt(i,:), thv_ds_zm(i,:), rtp2(i,:), & ! intent(in)
                                     thlp2(i,:), w_1_zm(i,:), w_2_zm(i,:), & ! intent(in)
                                     varnce_w_1_zm(i,:), varnce_w_2_zm(i,:), & ! intent(in)
                                     mixt_frac_zm(i,:), l_implemented, em(i,:), & ! intent(in)
                                     wp2sclrp(i,:,:), sclrpthvp(i,:,:), sclrm_forcing(i,:,:), & ! intent(in) 
                                     sclrp2(i,:,:), exner(i,:), rcm(i,:), p_in_Pa(i,:), thvm(i,:), & ! intent(in)
                                     Cx_fnc_Richardson(i,:), & ! intent(in)
                                     pdf_implicit_coefs_terms(i), & ! intent(in)
                                     um_forcing(i,:), vm_forcing(i,:), ug(i,:), vg(i,:), & ! intent(in)
                                     wpthvp(i,:), fcor(i), um_ref(i,:), vm_ref(i,:), up2(i,:), & ! intent(in)
                                     vp2(i,:), uprcp(i,:), vprcp(i,:), rc_coef(i,:), rtm(i,:), & ! intent(in)
                                     wprtp(i,:), thlm(i,:), wpthlp(i,:), sclrm(i,:,:), wpsclrp(i,:,:), & ! intent(in)
                                     um(i,:), upwp(i,:), vm(i,:), vpwp(i,:), rtm_old(i,:), & ! intent(in)
                                     wprtp_old(i,:), thlm_old(i,:), wpthlp_old(i,:), & ! intent(in)
                                     sclrm_old(i,:,:), wpsclrp_old(i,:,:), um_old(i,:), & ! intent(in)
                                     upwp_old(i,:), vm_old(i,:), vpwp_old(i,:), & ! intent(in)
                                     l_predict_upwp_vpwp, l_lmm_stepping ) ! intent(in)
        end do
      end if
    end if

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_begin_update( gr(i), irtm_sdmp, rtm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )             ! intent(inout)
        end do
      end if

      do i = 1, ngrdcol
        rtm(i,:) = sponge_damp_xm( gr(i), dt, gr(i)%zt, rtm_ref(i,:), &
                                   rtm(i,:), rtm_sponge_damp_profile )
      end do

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_end_update( gr(i), irtm_sdmp, rtm(i,:) / dt, & ! intent(in)
                                stats_zt(i) )             ! intent(inout)
        end do
      end if

    endif ! rtm_sponge_damp_settings%l_sponge_damping

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_begin_update( gr(i), ithlm_sdmp, thlm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )               ! intent(inout)
        end do
      end if

      do i = 1, ngrdcol
        thlm(i,:) = sponge_damp_xm( gr(i), dt, gr(i)%zt, thlm_ref(i,:), &
                                    thlm(i,:), thlm_sponge_damp_profile )
      end do

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_end_update( gr(i), ithlm_sdmp, thlm(i,:) / dt, & ! intent(in)
                                stats_zt(i) )               ! intent(inout)
        end do
      end if

    end if ! thlm_sponge_damp_settings%l_sponge_damping

    if ( l_predict_upwp_vpwp ) then

      if ( uv_sponge_damp_settings%l_sponge_damping ) then

        if ( l_stats_samp ) then
          do i = 1, ngrdcol
             call stat_begin_update( gr(i), ium_sdmp, um(i,:) / dt, & ! intent(in)
                                     stats_zt(i) )           ! intent(inout)
             call stat_begin_update( gr(i), ivm_sdmp, vm(i,:) / dt, & ! intent(in)
                                     stats_zt(i) )           ! intent(inout)
          end do
        end if

        do i = 1, ngrdcol
          um(i,:) = sponge_damp_xm( gr(i), dt, gr(i)%zt, um_ref(i,:), &
                                    um(i,:), uv_sponge_damp_profile )
        end do
        
        do i = 1, ngrdcol
          vm(i,:) = sponge_damp_xm( gr(i), dt, gr(i)%zt, vm_ref(i,:), &
                                    vm(i,:), uv_sponge_damp_profile )
        end do

        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_end_update( gr(i), ium_sdmp, um(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
            call stat_end_update( gr(i), ivm_sdmp, vm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
          end do
        end if

      end if ! uv_sponge_damp_settings%l_sponge_damping

      ! Adjust um and vm if nudging is turned on.
      if ( l_uv_nudge ) then

        ! Reflect nudging in budget
        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_begin_update( gr(i), ium_ndg, um(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )          ! intent(inout)
            call stat_begin_update( gr(i), ivm_ndg, vm(i,:) / dt, & ! intent(in)
                                    stats_zt(i) )          ! intent(inout)
          end do
        end if
      
        um(:,:) = um(:,:) - ( ( um(:,:) - um_ref(:,:) ) * (dt/ts_nudge) )
        vm(:,:) = vm(:,:) - ( ( vm(:,:) - vm_ref(:,:) ) * (dt/ts_nudge) )

        ! Reflect nudging in budget
        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            call stat_end_update( gr(i), ium_ndg, um(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )          ! intent(inout)
            call stat_end_update( gr(i), ivm_ndg, vm(i,:) / dt, & ! intent(in)
                                  stats_zt(i) )          ! intent(inout)
          end do
        end if

      end if ! l_uv_nudge

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_update_var( ium_ref, um_ref(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)
          call stat_update_var( ivm_ref, vm_ref(i,:), & ! intent(in)
                                stats_zt(i) )         ! intent(inout)
        end do
      end if

    end if ! l_predict_upwp_vpwp

    return

  end subroutine advance_xm_wpxp

  !======================================================================================
  subroutine xm_wpxp_lhs( nz, ngrdcol, gr, l_iter, dt, wpxp, wm_zt, C7_Skw_fnc, & ! In
                          wpxp_upper_lim, wpxp_lower_lim,                       & ! In
                          l_implemented, lhs_diff_zm, lhs_diff_zt,              & ! In
                          lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm,         & ! In
                          lhs_tp, lhs_pr1, lhs_ac_pr2,                          & ! In
                          l_diffuse_rtm_and_thlm,                               & ! In
                          lhs )                                                   ! Out
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
    !----------------------------------------------------------------------------------
    
    use grid_class, only:  & 
        grid, & ! Type
        zm2zt, & ! Procedure(s)
        ddzt

    use constants_clubb, only: &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use clip_semi_implicit, only: & 
        clip_semi_imp_lhs ! Procedure(s)

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

    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    real( kind = core_rknd ), intent(in) ::  & 
      dt    ! Timestep                                  [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: & 
      wpxp,                   & ! w'x' (momentum levs) at timestep (t) [un vary]
      wm_zt,                  & ! w wind component on thermo. levels       [m/s]
      C7_Skw_fnc,             & ! C_7 parameter with Sk_w applied            [-]
      wpxp_upper_lim,         & ! Keeps corrs. from becoming > 1       [un vary]
      wpxp_lower_lim            ! Keeps corrs. from becoming < -1      [un vary]

    logical, intent(in) ::  & 
      l_implemented, & ! Flag for CLUBB being implemented in a larger model.
      l_iter

    !------------------- Output Variable -------------------
    real( kind = core_rknd ), intent(out), dimension(nsup+nsub+1,ngrdcol,2*nz) ::  & 
      lhs ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    !------------------- Local Variables -------------------
    ! Indices
    integer :: k
    integer :: k_xm, k_wpxp

    logical :: l_upper_thresh, l_lower_thresh ! flags for clip_semi_imp_lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for xm
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm,    & ! Mean advection contributions to lhs
      lhs_ta_wpxp     ! Turbulent advection contributions to lhs

    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(in) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      lhs_ac_pr2, & ! Accumulation of w'x' and w'x' pressure term 2
      lhs_pr1       ! Pressure term 1 for w'x'

    logical, intent(in) :: &
      l_diffuse_rtm_and_thlm ! This flag determines whether or not we want CLUBB to do diffusion
                             ! on rtm and thlm
      
    real (kind = core_rknd) :: &
      invrs_dt
        
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zero_vector    ! Vector of 0s
      
    integer :: i

    !------------------- Begin Code -------------------
    
    ! Initializations/precalculations
    invrs_dt = 1.0_core_rknd / dt

    ! Lower boundary for xm, lhs(:,1)
    do i = 1, ngrdcol
      lhs(1,i,1) = 0.0_core_rknd
      lhs(2,i,1) = 0.0_core_rknd
      lhs(3,i,1) = 1.0_core_rknd
      lhs(4,i,1) = 0.0_core_rknd
      lhs(5,i,1) = 0.0_core_rknd
    end do

    ! Lower boundary for w'x', lhs(:,2)
    do i = 1, ngrdcol
      lhs(1,i,2) = 0.0_core_rknd
      lhs(2,i,2) = 0.0_core_rknd
      lhs(3,i,2) = 1.0_core_rknd
      lhs(4,i,2) = 0.0_core_rknd
      lhs(5,i,2) = 0.0_core_rknd
    end do

    ! Combine xm and w'x' terms into LHS
    do k = 2, nz

      k_xm = 2*k - 1  ! xm at odd index values
      k_wpxp = 2*k    ! w'x' at even index values
        
      do i = 1, ngrdcol

        ! ---- sum xm terms ----
        
        lhs(1,i,k_xm) = zero

        lhs(2,i,k_xm) = lhs_ta_xm(1,i,k)

        lhs(3,i,k_xm) = invrs_dt

        lhs(4,i,k_xm) = lhs_ta_xm(2,i,k)
        
        lhs(5,i,k_xm) = zero

        ! ---- sum w'x' terms ----

        lhs(1,i,k_wpxp) = lhs_ma_zm(1,i,k) + lhs_diff_zm(1,i,k) &
                          + gamma_over_implicit_ts * lhs_ta_wpxp(1,i,k)

        lhs(2,i,k_wpxp) = lhs_tp(1,i,k)

        lhs(3,i,k_wpxp) = lhs_ma_zm(2,i,k) + lhs_diff_zm(2,i,k) + lhs_ac_pr2(i,k) &
                          + gamma_over_implicit_ts * ( lhs_ta_wpxp(2,i,k) + lhs_pr1(i,k) )

        lhs(4,i,k_wpxp) = lhs_tp(2,i,k)

        lhs(5,i,k_wpxp) = lhs_ma_zm(3,i,k) + lhs_diff_zm(3,i,k) &
                        + gamma_over_implicit_ts * lhs_ta_wpxp(3,i,k)
                        
      end do
    end do

    ! Upper boundary for w'x', , lhs(:,2*gr%nz)
    ! These were set in the loop above for simplicity, so they must be set properly here
    do i = 1, ngrdcol
      lhs(1,i,2*nz) = 0.0_core_rknd
      lhs(2,i,2*nz) = 0.0_core_rknd
      lhs(3,i,2*nz) = 1.0_core_rknd
      lhs(4,i,2*nz) = 0.0_core_rknd
      lhs(5,i,2*nz) = 0.0_core_rknd
    end do
    
    ! LHS time tendency
    if ( l_iter ) then
      do k = 2, nz-1
        k_wpxp = 2*k 
        do i = 1, ngrdcol
          lhs(3,i,k_wpxp) = lhs(3,i,k_wpxp) + invrs_dt
        end do
      end do
    end if
    
    ! Calculate diffusion terms for all thermodynamic grid level
    if ( l_diffuse_rtm_and_thlm ) then
      do k = 2, nz 
        k_xm = 2*k - 1
        do i = 1, ngrdcol
          lhs(1,i,k_xm) = lhs(1,i,k_xm) + lhs_diff_zt(1,i,k) 
          lhs(3,i,k_xm) = lhs(3,i,k_xm) + lhs_diff_zt(2,i,k)
          lhs(5,i,k_xm) = lhs(5,i,k_xm) + lhs_diff_zt(3,i,k)
        end do
      end do
    end if
    
    ! Calculate mean advection terms for all momentum grid level
    if ( .not. l_implemented ) then
      do k = 2, nz 
        k_xm = 2*k - 1
        do i = 1, ngrdcol
          lhs(1,i,k_xm) = lhs(1,i,k_xm) + lhs_ma_zt(1,i,k)
          lhs(3,i,k_xm) = lhs(3,i,k_xm) + lhs_ma_zt(2,i,k)
          lhs(5,i,k_xm) = lhs(5,i,k_xm) + lhs_ma_zt(3,i,k)
        end do
      end do
    end if

    return

  end subroutine xm_wpxp_lhs

  !=============================================================================================
  subroutine calc_xm_wpxp_lhs_terms( nz, ngrdcol, gr, Kh_zm, wm_zm, wm_zt, wp2,         & ! In
                                     Kw6, C7_Skw_fnc, invrs_rho_ds_zt,                  & ! In
                                     invrs_rho_ds_zm, rho_ds_zt,                        & ! In
                                     rho_ds_zm, l_implemented, em,                      & ! In
                                     Lscale, thlm, exner, rtm, rcm, p_in_Pa, thvm,      & ! In
                                     ice_supersat_frac,                                 & ! In
                                     clubb_params, nu_vert_res_dep,                     & ! In
                                     l_diffuse_rtm_and_thlm,                            & ! In
                                     l_stability_correct_Kh_N2_zm,                      & ! In
                                     l_upwind_xm_ma,                                    & ! In
                                     l_brunt_vaisala_freq_moist,                        & ! In
                                     l_use_thvm_in_bv_freq,                             & ! In
                                     lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm,    & ! Out
                                     lhs_tp, lhs_ta_xm, lhs_ac_pr2 )                      ! Out
    ! Description:
    !   Calculate various xm and w'x' terms. These are general terms that are the same
    !   for multiple LHS matrices, so we save computations by calculating them once
    !   here, then reusing them where needed.
    !
    !-------------------------------------------------------------------------------------------
    
    use grid_class, only:  & 
        grid, & ! Type
        zm2zt, & ! Procedure(s)
        zt2zm

    use parameter_indices, only: &
        nparams, & ! Variable(s)
        ilambda0_stability_coef

    use parameters_tunable, only: &
        nu_vertical_res_dep    ! Type(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use advance_helper_module, only: &
        calc_stability_correction
      
    use mean_adv, only: & 
        term_ma_zt_lhs, &
        term_ma_zm_lhs

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs

    use diffusion, only:  & 
        diffusion_zt_lhs, &
        diffusion_zm_lhs

    use constants_clubb, only: &
        zero_threshold, &
        zero

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    !------------------- Input Variables -------------------
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: & 
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
      Kw6,                    & ! Coef. of eddy diffusivity for w'x'     [m^2/s]
      C7_Skw_fnc,             & ! C_7 parameter with Sk_w applied            [-]
      rho_ds_zm,              & ! Dry, static density on momentum levs. [kg/m^3]
      rho_ds_zt,              &
      invrs_rho_ds_zm,        &
      invrs_rho_ds_zt,        &  ! Inv. dry, static density at t-levs.   [m^3/kg]
      ice_supersat_frac

    logical, intent(in) ::  & 
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), dimension(ngrdcol), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    logical, intent(in) :: &
      l_diffuse_rtm_and_thlm,       & ! This flag determines whether or not we want CLUBB to do
                                      ! diffusion on rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! This flag determines whether or not we want CLUBB to apply
                                      ! a stability correction
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms.
                                      ! It affects rtm, thlm, sclrm, um and vm.
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq           ! Use thvm in the calculation of Brunt-Vaisala frequency
      
    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(out) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for xm
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm       ! Mean advection contributions to lhs

    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(out) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: & 
      lhs_ac_pr2    ! Accumulation of w'x' and w'x' pressure term 2
      
      
    !------------------- Local Variables -------------------
    real (kind = core_rknd), dimension(ngrdcol,nz) :: &
      Kh_N2_zm, &
      K_zm, &      ! Coef. of eddy diffusivity at momentum level (k)   [m^2/s]
      K_zt, &      ! Eddy diffusivity coefficient, thermo. levels [m2/s]
      Kw6_zm       ! Eddy diffusivity coefficient, momentum levels [m2/s]

    real (kind = core_rknd) :: &
      constant_nu ! controls the magnitude of diffusion
      
    integer :: i

    !------------------- Begin Code -------------------
    
    ! Initializations/precalculations
    constant_nu = 0.1_core_rknd
    Kw6_zm      = max( zt2zm( nz, ngrdcol, gr, Kw6 ), zero_threshold )
 
    ! Calculate turbulent advection terms of xm for all grid levels
    call xm_term_ta_lhs( nz, ngrdcol, gr,            & ! Intent(in)
                         rho_ds_zm, invrs_rho_ds_zt, & ! Intent(in)
                         lhs_ta_xm )                   ! Intent(out) 
    
                                   
    ! Calculate turbulent production terms of w'x' for all grid level
    call wpxp_term_tp_lhs( nz, ngrdcol, gr, wp2, & ! Intent(in)
                           lhs_tp )                ! Intent(out)

    ! Calculate accumulation of w'x' and w'x' pressure term 2 of w'x' for all grid level
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wpxp_pr
    do i = 1, ngrdcol
      call wpxp_terms_ac_pr2_lhs( gr(i), C7_Skw_fnc(i,:), wm_zt(i,:), gr(i)%invrs_dzm(:), & ! Intent(in)
                                  lhs_ac_pr2(i,:)                       ) ! Intent(out)
    end do        

    ! Calculate diffusion terms for all momentum grid level
    do i = 1, ngrdcol
      call diffusion_zm_lhs( gr(i), Kw6(i,:), Kw6_zm(i,:), nu_vert_res_dep(i)%nu6, & ! Intent(in)
                             gr(i)%invrs_dzt(:), gr(i)%invrs_dzm(:), & ! Intent(in)
                             invrs_rho_ds_zm(i,:), rho_ds_zt(i,:), & ! Intent(in)
                             lhs_diff_zm(:,i,:)                  ) ! Intent(out)    
    end do        
                              
    ! Calculate mean advection terms for all momentum grid level
    do i = 1, ngrdcol
      call term_ma_zm_lhs( gr(i), wm_zm(i,:), gr(i)%invrs_dzm(:), & ! Intent(in)
                           lhs_ma_zm(:,i,:)             ) ! Intent(out) 
    end do        
                               
    ! Calculate diffusion terms for all thermodynamic grid level
    if ( l_diffuse_rtm_and_thlm ) then
        
        if ( l_stability_correct_Kh_N2_zm ) then
          do i = 1, ngrdcol
            Kh_N2_zm(i,:) = Kh_zm(i,:) / &
              calc_stability_correction( gr(i), thlm(i,:), Lscale(i,:), em(i,:), exner(i,:), rtm(i,:), rcm(i,:), &
                                         p_in_Pa(i,:), thvm(i,:), ice_supersat_frac(i,:), &
                                         clubb_params(ilambda0_stability_coef), &
                                         l_brunt_vaisala_freq_moist, &
                                         l_use_thvm_in_bv_freq )
          end do
        else
          Kh_N2_zm(:,:) = Kh_zm(:,:)
        end if

        K_zm(:,:) = Kh_N2_zm(:,:) + constant_nu
        K_zt(:,:) = max( zm2zt ( nz, ngrdcol, gr, K_zm ), zero_threshold )

        do i = 1, ngrdcol
          call diffusion_zt_lhs( gr(i), K_zm(i,:), K_zt(i,:), zero,       & ! Intent(in)
                                 gr(i)%invrs_dzm(:), gr(i)%invrs_dzt(:), & ! Intent(in)
                                 invrs_rho_ds_zt(i,:), rho_ds_zm(i,:), & ! intent(in)
                                 lhs_diff_zt(:,i,:)                  ) ! Intent(out)
        end do
        
    end if        
                             
    ! Calculate mean advection terms for all thermodynamic grid level
    if ( .not. l_implemented ) then
      do i = 1, ngrdcol
        call term_ma_zt_lhs( gr(i), wm_zt(i,:), gr(i)%invrs_dzt(:), gr(i)%invrs_dzm(:), & ! Intent(in)
                             l_upwind_xm_ma,                             & ! Intent(in)
                             lhs_ma_zt(:,i,:)                              ) ! Intent(out)
      end do
    end if    
     
    return

  end subroutine calc_xm_wpxp_lhs_terms

  !=============================================================================
  subroutine xm_wpxp_rhs( nz, ngrdcol, gr, solve_type, l_iter, dt, xm, wpxp, & ! In
                          xm_forcing, wpxp_forcing, C7_Skw_fnc, & ! In
                          xpthvp, rhs_ta, thv_ds_zm, & ! In
                          lhs_pr1, lhs_ta_wpxp, & ! In
                          stats_zt, stats_zm, & ! intent(inout)
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
        grid ! Type

    use constants_clubb, only:  &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs, & ! Procedure(s)
        xpyp_term_ta_pdf_rhs

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use clip_semi_implicit, only: & 
        clip_semi_imp_rhs ! Procedure(s)

    use stats_type_utilities, only: & 
        stat_update_var,      & ! Procedure(s)
        stat_update_var_pt,   & 
        stat_begin_update_pt, &
        stat_modify_pt

    use stats_variables, only: & 
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

    use stats_type, only: stats ! Type

    implicit none

    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol 

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
  
    integer, intent(in) :: & 
      solve_type  ! Variables being solved for.

    logical, intent(in) :: l_iter

    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                  [s]
      
    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: &
      lhs_ta_wpxp   ! Turbulent advection terms of w'x'

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      xm,                     & ! xm (thermodynamic levels)               [x un]
      wpxp,                   & ! <w'x'> (momentum levels)          [{x un} m/s]
      xm_forcing,             & ! xm forcings (thermodynamic levels)  [{x un}/s]
      wpxp_forcing,           & ! <w'x'> forcing (momentum levs)  [{x un} m/s^2]
      C7_Skw_fnc,             & ! C_7 parameter with Sk_w applied            [-]
      xpthvp,                 & ! x'th_v' (momentum levels)           [{x un} K]
      thv_ds_zm,              & ! Dry, base-state theta_v on mom. levs.      [K]
      lhs_pr1,                & ! Pressure term 1 for w'x'
      rhs_ta

    !------------------- InOut Variables -------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm

    !------------------- Output Variable -------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,2*nz) ::  & 
      rhs  ! Right-hand side of band diag. matrix. (LAPACK)

    !------------------- Local Variables -------------------

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
        rhs_bp_pr3, & ! Buoyancy production of w'x' and w'x' pressure term 3
        rhs_bp,     & ! Buoyancy production of w'x' (stats only)
        rhs_pr3       ! w'x' pressure term 3 (stats only)
      
    real( kind = core_rknd ) :: &
        invrs_dt

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zero_vector    ! Vector of 0s

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
      
    integer :: i

    !------------------- Begin Code -------------------

    ! Initialize output array and precalculate the reciprocal of dt
    invrs_dt = 1.0_core_rknd / dt    
                                  
    ! Calculate buoyancy production of w'x' and w'x' pressure term 3
    call wpxp_terms_bp_pr3_rhs( nz, ngrdcol, gr, C7_Skw_fnc, thv_ds_zm, xpthvp, & ! intent(in)
                                rhs_bp_pr3 )                                      ! intent(out)
                            
    do i = 1, ngrdcol
      ! Set lower boundary for xm
      rhs(i,1) = xm(i,1)

      ! Set lower boundary for w'x'
      rhs(i,2) = wpxp(i,1)
    end do

    ! Combine terms to calculate other values, rhs(3) to rhs(gr%nz-2)
    do k = 2, nz-1

      k_xm   = 2*k - 1
      k_wpxp = 2*k

      do i = 1, ngrdcol
        ! RHS time tendency and forcings for xm
        ! Note: xm forcings include the effects of microphysics,
        !       cloud water sedimentation, radiation, and any
        !       imposed forcings on xm.
        rhs(i,k_xm) =  xm(i,k) * invrs_dt + xm_forcing(i,k)


        ! Calculate rhs values for w'x' using precalculated terms
        rhs(i,k_wpxp) =  rhs_bp_pr3(i,k) + wpxp_forcing(i,k) + rhs_ta(i,k) &
                                  + ( one - gamma_over_implicit_ts ) &
                                  * ( - lhs_ta_wpxp(1,i,k) * wpxp(i,k+1) &
                                      - lhs_ta_wpxp(2,i,k) * wpxp(i,k) &
                                      - lhs_ta_wpxp(3,i,k) * wpxp(i,k-1) &
                                      - lhs_pr1(i,k) * wpxp(i,k) )
      end do
    end do

    do i = 1, ngrdcol
      ! Upper boundary for xm
      rhs(i,2*nz-1) = xm(i,nz) * invrs_dt + xm_forcing(i,nz)

      ! Upper boundary for w'x', rhs(2*gr%nz)
      rhs(i,2*nz) = 0.0_core_rknd
    end do

    ! RHS time tendency.
    if ( l_iter ) then
      do k = 2, nz-1
        k_wpxp = 2*k
        do i = 1, ngrdcol
          rhs(i,k_wpxp) = rhs(i,k_wpxp) + wpxp(i,k) * invrs_dt
        end do
      end do
    end if
    

    if ( l_stats_samp ) then

      zero_vector = zero

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

      ! Statistics: explicit contributions for wpxp.

      ! w'x' term bp is completely explicit; call stat_update_var.
      ! Note:  To find the contribution of w'x' term bp, substitute 0 for the
      !        C_7 skewness function input to function wpxp_terms_bp_pr3_rhs.
      call wpxp_terms_bp_pr3_rhs( nz, ngrdcol, gr, zero_vector, thv_ds_zm, xpthvp, & ! intent(in)
                                  rhs_bp )                          ! intent(out)
                                    
      do i = 1, ngrdcol
        call stat_update_var( iwpxp_bp, rhs_bp(i,:), & ! intent(in)
                              stats_zm(i) )          ! intent(inout)
      end do

      ! w'x' term pr3 is completely explicit; call stat_update_var.
      ! Note:  To find the contribution of w'x' term pr3, add 1 to the
      !        C_7 skewness function input to function wpxp_terms_bp_pr2_rhs.
      call wpxp_terms_bp_pr3_rhs( nz, ngrdcol, gr, (one+C7_Skw_fnc), thv_ds_zm, xpthvp, & ! intent(in)
                                  rhs_pr3 ) ! intent(out)
                                  
      do i = 1, ngrdcol
        call stat_update_var( iwpxp_pr3, rhs_pr3(i,:), & ! intent(in)
                              stats_zm(i) )            ! intent(inout)
      end do

      do k = 2, nz-1
        do i = 1, ngrdcol

          ! w'x' forcing term is completely explicit; call stat_update_var_pt.
          call stat_update_var_pt( iwpxp_f, k, wpxp_forcing(i,k), & ! intent(in)
                                   stats_zm(i) )                     ! intent(inout)


          ! <w'x'> term ta has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on
          ! xpyp_term_ta_pdf_rhs.
          call stat_begin_update_pt( iwpxp_ta, k, -rhs_ta(i,k), & ! intent(in) 
                                     stats_zm(i) )                 ! intent(inout)

          ! Note:  An "over-implicit" weighted time step is applied to this term.
          !        A weighting factor of greater than 1 may be used to make the
          !        term more numerically stable (see note above for RHS
          !        contribution from "over-implicit" weighted time step for LHS
          !        turbulent advection (ta) term).
          call stat_modify_pt( iwpxp_ta, k, & ! intent(in)
                               + ( one - gamma_over_implicit_ts ) &
                                 * ( - lhs_ta_wpxp(1,i,k) * wpxp(i,k+1) &
                                     - lhs_ta_wpxp(2,i,k) * wpxp(i,k) &
                                     - lhs_ta_wpxp(3,i,k) * wpxp(i,k-1) ), & ! intent(in)
                               stats_zm(i) ) ! intent(inout)

          ! w'x' term pr1 is normally completely implicit.  However, there is a
          ! RHS contribution from the "over-implicit" weighted time step.  A
          ! weighting factor of greater than 1 may be used to make the term more
          ! numerically stable (see note above for RHS contribution from
          ! "over-implicit" weighted time step for LHS turbulent advection (ta)
          ! term).  Therefore, w'x' term pr1 has both implicit and explicit
          ! components; call stat_begin_update_pt.  Since stat_begin_update_pt
          ! automatically subtracts the value sent in, reverse the sign on the
          ! input value.
          call stat_begin_update_pt( iwpxp_pr1, k, & ! intent(in)
                                    - ( one - gamma_over_implicit_ts )  &
                                    * ( - lhs_pr1(i,k) * wpxp(i,k) ), & ! intent(in)
                                     stats_zm(i) ) ! intent(inout)
        end do
      end do

      
      ! Statistics: explicit contributions for xm
      !             (including microphysics/radiation).

      ! xm forcings term is completely explicit; call stat_update_var_pt.
      do k = 2, nz
        do i = 1, ngrdcol
          call stat_update_var_pt( ixm_f, k, xm_forcing(i,k), & ! intent(in)
                                   stats_zt(i) )                ! intent(inout)
        end do
      end do

    endif ! l_stats_samp

    return

  end subroutine xm_wpxp_rhs
  
  !=============================================================================================
  subroutine calc_xm_wpxp_ta_terms( nz, ngrdcol, gr, wp2rtp, &
                                    wp2thlp, wp2sclrp, &
                                    rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm, &
                                    sigma_sqd_w, wp3_on_wp2_zt, &
                                    pdf_implicit_coefs_terms, &
                                    iiPDF_type, &
                                    l_explicit_turbulent_adv_wpxp, l_predict_upwp_vpwp, &
                                    l_scalar_calc, &
                                    l_godunov_upwind_wpxp_ta, &
                                    stats_zt, & 
                                    lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpup, &
                                    lhs_ta_wpvp, lhs_ta_wpsclrp, &
                                    rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup, &
                                    rhs_ta_wpvp, rhs_ta_wpsclrp )
  !
  ! Description: This subroutine calculates the turbulent advection terms for 
  !              the left and right hand side matrices. Solutions may be entirely
  !              explicit, entirely implicit, or mixed between, depending on 
  !              various flags and the PDF type. 
  !---------------------------------------------------------------------------------------------
                                    
    use grid_class, only: &
        grid, & ! Type
        zt2zm,  & ! Procedure(s)
        zm2zt
      
    use clubb_precision, only: &
        core_rknd  ! Variable(s)
      
    use constants_clubb, only: &
        one, &
        zero, &
        zero_threshold
      
    use parameters_model, only: &
        sclr_dim  ! Number of passive scalar variables
      
    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs, &  ! Procedures
        xpyp_term_ta_pdf_lhs_godunov, &
        xpyp_term_ta_pdf_rhs, &
        sgn_turbulent_velocity
      
    use model_flags, only: &
        iiPDF_ADG1,       & ! Integer constants
        iiPDF_new,        &
        iiPDF_new_hybrid
      
    use stats_variables, only: &
        l_stats_samp,             & ! Logical constant
        icoef_wp2rtp_implicit, &
        iterm_wp2rtp_explicit, &
        icoef_wp2thlp_implicit, &
        iterm_wp2thlp_explicit
      
    use stats_type_utilities, only: & 
        stat_update_var   ! Procedure(s)

    use stats_type, only: stats ! Type

    implicit none 
    
    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    type(implicit_coefs_terms), dimension(ngrdcol), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]
                                
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wp2rtp, &
      wp2thlp, &
      rho_ds_zt, &
      invrs_rho_ds_zm, &                  
      rho_ds_zm, &           
      sigma_sqd_w, &     
      wp3_on_wp2_zt
      
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(in) :: &
      wp2sclrp
      
    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_explicit_turbulent_adv_wpxp, &
      l_scalar_calc, &
      l_predict_upwp_vpwp

    logical, intent(in) :: &
      l_godunov_upwind_wpxp_ta    ! This flag determines whether we want to use an upwind
                                  ! differencing approximation rather than a centered 
                                  ! differencing for turbulent advection terms. 
                                  ! It affects  wpxp only.

    logical, parameter :: &
      l_dummy_false = .false. ! This flag is set to false in order to replace the flag
                              ! passed into the xpyp_term_ta_pdf_rhs subroutine.
                              ! This stems from removing the l_upwind_wpxp_ta flag.
                              ! More information on this can be found on issue #926
                              ! on the clubb repository.
                              
    !------------------- Inout Variables -------------------
    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt
      
    !------------------- Output Variables -------------------
        
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(out) :: &
      lhs_ta_wprtp, &
      lhs_ta_wpthlp, &
      lhs_ta_wpup, &
      lhs_ta_wpvp
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz,sclr_dim), intent(out) :: &
      lhs_ta_wpsclrp
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rhs_ta_wprtp, &
      rhs_ta_wpthlp, &
      rhs_ta_wpup, &
      rhs_ta_wpvp
      
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(out) :: &
      rhs_ta_wpsclrp
    
    !------------------- Local Variables -------------------

    ! Variables for turbulent advection of predictive variances and covariances.

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      coef_wp2rtp_implicit_zm, & ! coef_wp2rtp_implicit interp. to m-levs. [m/s]
      term_wp2rtp_explicit_zm    ! term_wp2rtp_expl intrp m-levs [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      coef_wp2thlp_implicit_zm, & ! coef_wp2thlp_implicit interp. m-levs.  [m/s]
      term_wp2thlp_explicit_zm    ! term_wp2thlp_expl interp. m-levs [m^2/s^2 K]

    ! <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'> + term_wp2sclrp_explicit
    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      term_wp2sclrp_explicit    ! Term that is on the RHS    [m^2/s^2(un. vary)]

    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      term_wp2sclrp_explicit_zm    ! term_wp2sclrp_expl intrp zm [m^2/s^2(un v)]

    ! Sign of turbulent velocity (used for "upwind" turbulent advection)
    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sgn_t_vel_wprtp,  & ! Sign of the turbulent velocity for <w'rt'>       [-]
      sgn_t_vel_wpthlp    ! Sign of the turbulent velocity for <w'thl'>      [-]

    real ( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sgn_t_vel_wpsclrp    ! Sign of the turbulent velocity for <w'sclr'>    [-]
    
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      a1, &
      a1_zt
      
    integer :: i, j, k
    
    !------------------- Begin Code -------------------
    
    ! Set up the implicit coefficients and explicit terms for turbulent
    ! advection of <w'rt'>, <w'thl'>, and <w'sclr'>.
    if ( l_explicit_turbulent_adv_wpxp ) then

      ! The turbulent advection of <w'x'> is handled explicitly
       
      ! The turbulent advection of <w'x'> is handled explicitly, the
      ! terms are calculated only for the RHS matrices. The 
      ! term_wp2xp_explicit terms are equal to <w'x'> as calculated using PDF
      ! parameters, which are general for any PDF type. The values of
      ! <w'x'> are calculated on thermodynamic levels.
       
      ! These coefficients only need to be set if stats output is on
      if ( l_stats_samp ) then
        coef_wp2rtp_implicit(:,:)  = zero
        coef_wp2thlp_implicit(:,:) = zero
      end if
       
      ! The turbulent advection terms are handled entirely explicitly. Thus the LHS
      ! terms can be set to zero.
      lhs_ta_wprtp(:,:,:) = zero
      lhs_ta_wpthlp(:,:,:) = zero
       
      if ( l_scalar_calc ) then
        lhs_ta_wpsclrp(:,:,:,:) = zero
      end if
       
      term_wp2rtp_explicit(:,:)  = wp2rtp(:,:)
      term_wp2thlp_explicit(:,:) = wp2thlp(:,:)
      
      ! Calculate the RHS turbulent advection term for <w'r_t'>
      do i = 1, ngrdcol
        call xpyp_term_ta_pdf_rhs( gr(i), term_wp2rtp_explicit(i,:),       & ! Intent(in)
                                   rho_ds_zt(i,:),                  & ! Intent(in)
                                   invrs_rho_ds_zm(i,:),            & ! Intent(in)
                                   gr(i)%invrs_dzm(:),               & ! Intent(in)
                                   l_dummy_false,                 & ! Intent(in)
                                   sgn_t_vel_wprtp(i,:),            & ! Intent(in)
                                   term_wp2rtp_explicit_zm(i,:),    & ! Intent(in)
                                   rho_ds_zm(i,:),                  & ! Intent(in)
                                   gr(i)%invrs_dzt(:),               & ! Intent(in)
                                   rhs_ta_wprtp(i,:)                ) ! Intent(out)
      end do
       
      ! Calculate the RHS turbulent advection term for <w'thl'>
      do i = 1, ngrdcol
        call xpyp_term_ta_pdf_rhs( gr(i), term_wp2thlp_explicit(i,:),      & ! Intent(in)
                                   rho_ds_zt(i,:),                  & ! Intent(in)
                                   invrs_rho_ds_zm(i,:),            & ! Intent(in)
                                   gr(i)%invrs_dzm(:),               & ! Intent(in)
                                   l_dummy_false,                 & ! Intent(in)
                                   sgn_t_vel_wpthlp(i,:),           & ! Intent(in)
                                   term_wp2thlp_explicit_zm(i,:),   & ! Intent(in)
                                   rho_ds_zm(i,:),                  & ! Intent(in)
                                   gr(i)%invrs_dzt(:),               & ! Intent(in)
                                   rhs_ta_wpthlp(i,:)               ) ! Intent(out)  
      end do
                                     
      do j = 1, sclr_dim, 1
        
        term_wp2sclrp_explicit(:,:) = wp2sclrp(:,:,j)
        
        ! Calculate the RHS turbulent advection term for <w'thl'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_rhs( gr(i), term_wp2sclrp_explicit(i,:),     & ! Intent(in)
                                     rho_ds_zt(i,:),                  & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),            & ! Intent(in)
                                     gr(i)%invrs_dzm(:),               & ! Intent(in)
                                     l_dummy_false,                 & ! Intent(in)
                                     sgn_t_vel_wpsclrp(i,:),          & ! Intent(in)
                                     term_wp2sclrp_explicit_zm(i,:),  & ! Intent(in)
                                     rho_ds_zm(i,:),                  & ! Intent(in)
                                     gr(i)%invrs_dzt(:),               & ! Intent(in)
                                     lhs_ta_wpsclrp(:,i,:,j)            ) ! Intent(out)  
        end do

      end do ! i = 1, sclr_dim, 1

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
        a1(:,:) = one / ( one - sigma_sqd_w(:,:) )

        ! Interpolate a_1 from momentum levels to thermodynamic levels.  This
        ! will be used for the <w'x'> turbulent advection (ta) term.
        a1_zt(:,:) = max( zm2zt( nz, ngrdcol, gr, a1 ), zero_threshold )   ! Positive def. quantity
        
        coef_wp2rtp_implicit(:,:) = a1_zt(:,:) * wp3_on_wp2_zt(:,:)
        coef_wp2thlp_implicit(:,:) = coef_wp2rtp_implicit(:,:)
       
        if ( .not. l_godunov_upwind_wpxp_ta ) then
 
          ! Calculate the LHS turbulent advection term for <w'r_t'>
          do i = 1, ngrdcol
            call xpyp_term_ta_pdf_lhs( gr(i), coef_wp2rtp_implicit(i,:),     & ! Intent(in)
                                       rho_ds_zt(i,:),                & ! Intent(in)
                                       invrs_rho_ds_zm(i,:),          & ! Intent(in)
                                       gr(i)%invrs_dzm(:),             & ! Intent(in)
                                       l_dummy_false,               & ! Intent(in)
                                       sgn_t_vel_wprtp(i,:),          & ! Intent(in)
                                       coef_wp2rtp_implicit_zm(i,:),  & ! Intent(in)
                                       rho_ds_zm(i,:),                & ! Intent(in)
                                       gr(i)%invrs_dzt(:),             & ! Intent(in)
                                       lhs_ta_wprtp(:,i,:)           ) ! Intent(out)
          end do
 
        else

          ! Godunov-like method for the vertical discretization of ta term  
          coef_wp2rtp_implicit(:,:) = a1_zt(:,:) * wp3_on_wp2_zt(:,:)
          coef_wp2thlp_implicit(:,:) = coef_wp2rtp_implicit(:,:)

          do i = 1, ngrdcol
            call xpyp_term_ta_pdf_lhs_godunov( gr(i), coef_wp2rtp_implicit(i,:),     & ! Intent(in)
                                               invrs_rho_ds_zm(i,:),          & ! Intent(in)
                                               gr(i)%invrs_dzm(:),             & ! Intent(in)
                                               rho_ds_zm(i,:),                & ! Intent(in)
                                               lhs_ta_wprtp(:,i,:)           )  ! Intent(out)
          end do
      
        end if
 
        ! For ADG1, the LHS turbulent advection terms for 
        ! <w'r_t'>, <w'thl'>, <w'sclr'> are all equal
        lhs_ta_wpthlp(:,:,:) = lhs_ta_wprtp(:,:,:)
        
        if ( l_scalar_calc ) then
          do j = 1, sclr_dim
            lhs_ta_wpsclrp(:,:,:,j) = lhs_ta_wprtp(:,:,:)
          end do
        end if
        
        if ( l_stats_samp ) then
          term_wp2rtp_explicit(:,:) = zero
          term_wp2thlp_explicit(:,:) = zero
        end if

        ! The <w'r_t'>, <w'thl'>, <w'sclr'> turbulent advection terms are entirely implicit.
        ! Set the RHS turbulent advection terms to 0
        rhs_ta_wprtp(:,:) = zero
        rhs_ta_wpthlp(:,:) = zero
        rhs_ta_wpsclrp(:,:,:) = zero
        
        if ( l_predict_upwp_vpwp ) then
            
          ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.
          ! These terms are equal to the <w'r_t'> terms as well in this case
          lhs_ta_wpup(:,:,:) = lhs_ta_wprtp(:,:,:)
          lhs_ta_wpvp(:,:,:) = lhs_ta_wprtp(:,:,:)
          
          ! The <w'u'> and <w'v'> turbulent advection terms are entirely implicit.
          ! Set the RHS turbulent advection terms to 0
          rhs_ta_wpup(:,:) = zero
          rhs_ta_wpvp(:,:) = zero

        endif  

      elseif ( iiPDF_type == iiPDF_new ) then

        ! The new PDF is used.

        ! Unpack the variables coef_wp2rtp_implicit, term_wp2rtp_explicit,
        ! coef_wp2thlp_implicit, and term_wp2thlp_explicit from
        ! pdf_implicit_coefs_terms.  The PDF parameters and the resulting
        ! implicit coefficients and explicit terms are calculated on
        ! thermodynamic levels.
        do i = 1, ngrdcol
          coef_wp2rtp_implicit(i,:)  = pdf_implicit_coefs_terms(i)%coef_wp2rtp_implicit(:)
          coef_wp2thlp_implicit(i,:) = pdf_implicit_coefs_terms(i)%coef_wp2thlp_implicit(:)
          term_wp2rtp_explicit(i,:)  = pdf_implicit_coefs_terms(i)%term_wp2rtp_explicit(:)
          term_wp2thlp_explicit(i,:) = pdf_implicit_coefs_terms(i)%term_wp2thlp_explicit(:)
        end do


        ! Calculate the LHS turbulent advection term for <w'rt'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_lhs( gr(i), coef_wp2rtp_implicit(i,:),      & ! Intent(in)
                                     rho_ds_zt(i,:),                 & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),           & ! Intent(in)
                                     gr(i)%invrs_dzm(:),              & ! Intent(in)
                                     l_dummy_false,                & ! Intent(in)
                                     sgn_t_vel_wprtp(i,:),           & ! Intent(in)
                                     coef_wp2rtp_implicit_zm(i,:),   & ! Intent(in)
                                     rho_ds_zm(i,:),                 & ! Intent(in)
                                     gr(i)%invrs_dzt(:),              & ! Intent(in)
                                     lhs_ta_wprtp(:,i,:)             ) ! Intent(out) 
        end do
                                       
        ! Calculate the RHS turbulent advection term for <w'rt'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_rhs( gr(i), term_wp2rtp_explicit(i,:),      & ! Intent(in)
                                     rho_ds_zt(i,:),                 & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),           & ! Intent(in)
                                     gr(i)%invrs_dzm(:),              & ! Intent(in)
                                     l_dummy_false,                & ! Intent(in)
                                     sgn_t_vel_wprtp(i,:),           & ! Intent(in)
                                     term_wp2rtp_explicit_zm(i,:),   & ! Intent(in)
                                     rho_ds_zm(i,:),                 & ! Intent(in)
                                     gr(i)%invrs_dzt(:),              & ! Intent(in)
                                     rhs_ta_wprtp(i,:)           ) ! Intent(out)
        end do

        
        ! Calculate the LHS turbulent advection term for <w'thl'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_lhs( gr(i), coef_wp2thlp_implicit(i,:),      & ! Intent(in)
                                     rho_ds_zt(i,:),                  & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),            & ! Intent(in)
                                     gr(i)%invrs_dzm(:),               & ! Intent(in)
                                     l_dummy_false,                 & ! Intent(in)
                                     sgn_t_vel_wpthlp(i,:),           & ! Intent(in)
                                     coef_wp2thlp_implicit_zm(i,:),   & ! Intent(in)
                                     rho_ds_zm(i,:),                  & ! Intent(in)
                                     gr(i)%invrs_dzt(:),               & ! Intent(in)
                                     lhs_ta_wpthlp(:,i,:)             ) ! Intent(out) 
        end do
      
        ! Calculate the RHS turbulent advection term for <w'thl'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_rhs( gr(i), term_wp2thlp_explicit(i,:),      & ! Intent(in)
                                     rho_ds_zt(i,:),                  & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),            & ! Intent(in)
                                     gr(i)%invrs_dzm(:),               & ! Intent(in)
                                     l_dummy_false,                 & ! Intent(in)
                                     sgn_t_vel_wpthlp(i,:),           & ! Intent(in)
                                     term_wp2thlp_explicit_zm(i,:),   & ! Intent(in)
                                     rho_ds_zm(i,:),                  & ! Intent(in)
                                     gr(i)%invrs_dzt(:),               & ! Intent(in)
                                     rhs_ta_wpthlp(i,:)               ) ! Intent(out)
        end do

        ! The code for the scalar variables will be set up later.
        lhs_ta_wpsclrp(:,:,:,:) = zero
        rhs_ta_wpsclrp(:,:,:) = zero
      
      elseif ( iiPDF_type == iiPDF_new_hybrid ) then

        ! The new hybrid PDF is used.

        ! Unpack the variable coef_wp2rtp_implicit from the structure
        ! pdf_implicit_coefs_terms.  The values of coef_wp2thlp_implicit,
        ! coef_wp2up_implicit, coef_wp2vp_implict, and coef_wp2sclrp_implicit
        ! are all equal to coef_wp2rtp_implicit.  The PDF parameters and the
        ! resulting implicit coefficients are calculated on thermodynamic
        ! levels.
        do i = 1, ngrdcol
          coef_wp2rtp_implicit(i,:) = pdf_implicit_coefs_terms(i)%coef_wp2rtp_implicit(:)
          coef_wp2thlp_implicit(i,:) = coef_wp2rtp_implicit(i,:)
        end do
        

        ! Calculate the LHS turbulent advection term for <w'rt'>
        do i = 1, ngrdcol
          call xpyp_term_ta_pdf_lhs( gr(i), coef_wp2rtp_implicit(i,:),      & ! Intent(in)
                                     rho_ds_zt(i,:),                 & ! Intent(in)
                                     invrs_rho_ds_zm(i,:),           & ! Intent(in)
                                     gr(i)%invrs_dzm(:),              & ! Intent(in)
                                     l_dummy_false,                & ! Intent(in)
                                     sgn_t_vel_wprtp(i,:),           & ! Intent(in)
                                     coef_wp2rtp_implicit_zm(i,:),   & ! Intent(in)
                                     rho_ds_zm(i,:),                 & ! Intent(in)
                                     gr(i)%invrs_dzt(:),              & ! Intent(in)
                                     lhs_ta_wprtp(:,i,:)             ) ! Intent(out) 
        end do
                                       
        ! For the new hybrid PDF, the LHS turbulent advection terms for 
        ! <w'r_t'>, <w'thl'>, and <w'sclr'> are all the same.
        lhs_ta_wpthlp(:,:,:) = lhs_ta_wprtp(:,:,:)
        
        if ( l_scalar_calc ) then
          do j = 1, sclr_dim
            lhs_ta_wpsclrp(:,:,:,j) = lhs_ta_wprtp(:,:,:)
          end do
        end if
        
        if ( l_stats_samp ) then
          term_wp2rtp_explicit(:,:) = zero
          term_wp2thlp_explicit(:,:) = zero
        end if

        ! The <w'r_t'>, <w'thl'>, <w'sclr'> turbulent advection terms are
        ! entirely implicit.  Set the RHS turbulent advection terms to 0
        rhs_ta_wprtp(:,:) = zero
        rhs_ta_wpthlp(:,:) = zero
        rhs_ta_wpsclrp(:,:,:) = zero
        
        if ( l_predict_upwp_vpwp ) then
            
          ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.
          ! These terms are equal to the <w'r_t'> terms as well in this case
          lhs_ta_wpup(:,:,:) = lhs_ta_wprtp(:,:,:)
          lhs_ta_wpvp(:,:,:) = lhs_ta_wprtp(:,:,:)
          
          ! The <w'u'> and <w'v'> turbulent advection terms are entirely
          ! implicit.  Set the RHS turbulent advection terms to 0
          rhs_ta_wpup(:,:) = zero
          rhs_ta_wpvp(:,:) = zero

        endif  

      endif ! iiPDF_type

    endif ! l_explicit_turbulent_adv_xpyp
      
    if ( l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_update_var( icoef_wp2rtp_implicit, coef_wp2rtp_implicit(i,:), & ! intent(in)
                              stats_zt(i) )                                     ! intent(inout)
        call stat_update_var( iterm_wp2rtp_explicit, term_wp2rtp_explicit(i,:), & ! intent(in)
                              stats_zt(i) )                                     ! intent(inout)
        call stat_update_var( icoef_wp2thlp_implicit, coef_wp2thlp_implicit(i,:), & ! intent(in)
                              stats_zt(i) )                                       ! intent(inout)
        call stat_update_var( iterm_wp2thlp_explicit, term_wp2thlp_explicit(i,:), & ! intent(in)
                              stats_zt(i) )                                       ! intent(inout)
      end do
    endif
    
  end subroutine calc_xm_wpxp_ta_terms
  
  !==========================================================================================
  subroutine solve_xm_wpxp_with_single_lhs( nz, ngrdcol, gr, dt, l_iter, nrhs, wm_zt, wp2, &
                                            invrs_tau_C6_zm, tau_max_zm, &
                                            rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp, &
                                            thlm_forcing, wpthlp_forcing, rho_ds_zm, &
                                            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                            thv_ds_zm, rtp2, thlp2, l_implemented, &
                                            sclrpthvp, sclrm_forcing, sclrp2, um_forcing, &
                                            vm_forcing, ug, vg, uprcp, vprcp, rc_coef, fcor, &
                                            up2, vp2, &
                                            low_lev_effect, high_lev_effect, &
                                            C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc, &
                                            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, &
                                            lhs_ta_wpxp, &
                                            rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpup, &
                                            rhs_ta_wpvp, rhs_ta_wpsclrp, &
                                            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp, &
                                            lhs_pr1_wpthlp, lhs_pr1_wpsclrp, &
                                            C_uu_shr, &
                                            l_predict_upwp_vpwp, &
                                            l_diffuse_rtm_and_thlm, &
                                            l_upwind_xm_ma, &
                                            l_tke_aniso, &
                                            l_enable_relaxed_clipping, &
                                            order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, &
                                            stats_zt, stats_zm, stats_sfc, & 
                                            rtm, wprtp, thlm, wpthlp, &
                                            sclrm, wpsclrp, um, upwp, vm, vpwp )
    !            
    ! Description: This subroutine solves all xm_wpxp when all the LHS matrices are equal.
    !              The LHS matrices being equivalent allows for only a single solve, rather
    !              than a seperate solve for each field. 
    !----------------------------------------------------------------------------------------
    
    use grid_class, only: & 
        grid, & ! Type
        ddzt    ! Procedure(s)
      
    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constants
      
    use stats_type_utilities, only: & 
        stat_update_var   ! Procedure(s)
      
    use stats_variables, only: & 
        irtm_matrix_condt_num, &  ! Variables
        ithlm_matrix_condt_num, &
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
        ivpthvp, &
        l_stats_samp
        
    use parameters_model, only: & 
        sclr_dim, &  ! Variable(s)
        sclr_tol

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use constants_clubb, only:  & 
        fstderr, &  ! Constant
        rt_tol, &
        thl_tol, &
        w_tol, &
        w_tol_sqd, &
        thl_tol_mfl, &
        rt_tol_mfl, &
        zero, &
        one, &
        ep1

    use stats_type, only: stats ! Type

    implicit none
    
    ! ------------------- Input Variables -------------------

    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: & 
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      invrs_tau_C6_zm, & ! Inverse tau on momentum levels applied to C6 term [1/s]
      tau_max_zm,      & ! Max. allowable eddy dissipation time scale on m-levs [s]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
      thlpthvp,        & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,    & ! th_l forcing (thermodynamic levels)      [K/s]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)      [K/s^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on moment. levs. [K]
      rtp2,            & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2              ! th_l'^2 (momentum levels)                [K^2]

    logical, intent(in) ::  & 
      l_implemented, &      ! Flag for CLUBB being implemented in a larger model.
      l_iter

    ! Additional variables for passive scalars
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz,sclr_dim) :: & 
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      um_forcing, & ! <u> forcing term (thermodynamic levels)      [m/s^2]
      vm_forcing, & ! <v> forcing term (thermodynamic levels)      [m/s^2]
      ug,         & ! <u> geostrophic wind (thermodynamic levels)  [m/s]
      vg            ! <v> geostrophic wind (thermodynamic levels)  [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      uprcp,              & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp,              & ! < v' r_c' >              [(m kg)/(s kg)]
      rc_coef               ! Coefficient on X'r_c' in X'th_v' equation [K/(kg/kg)]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
      fcor          ! Coriolis parameter                           [s^-1]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      up2,    & ! Variance of the u wind component             [m^2/s^2]
      vp2       ! Variance of the v wind component             [m^2/s^2]

    ! LHS/RHS terms
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for w'x'
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm       ! Mean advection contributions to lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_ta_wpxp    ! w'r_t' turbulent advection contributions to lhs  
     
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      rhs_ta_wprtp,  & ! w'r_t' turbulent advection contributions to rhs  
      rhs_ta_wpthlp, & ! w'thl' turbulent advection contributions to rhs
      rhs_ta_wpup,   & ! w'u' turbulent advection contributions to rhs
      rhs_ta_wpvp      ! w'v' turbulent advection contributions to rhs
      
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(in) :: & 
      rhs_ta_wpsclrp    ! w'sclr' turbulent advection contributions to rhs

    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(in) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      lhs_ac_pr2,     & ! Accumulation of w'x' and w'x' pressure term 2
      lhs_pr1_wprtp,  & ! Pressure term 1 for w'r_t' for all grid levels
      lhs_pr1_wpthlp, & ! Pressure term 1 for w'thl' for all grid levels
      lhs_pr1_wpsclrp   ! Pressure term 1 for w'sclr' for all grid levels
      
    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(ngrdcol,nz), intent(in) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

    integer, intent(in) :: &
      nrhs         ! Number of RHS vectors

    real( kind = core_rknd ), intent(in) ::  &
      C_uu_shr    ! CLUBB tunable parameter C_uu_shr

    logical, intent(in) :: &
      l_predict_upwp_vpwp,       & ! Flag to predict <u'w'> and <v'w'> along
                                   ! with <u> and <v> alongside the advancement
                                   ! of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>,
                                   ! and <w'sclr'> in subroutine advance_xm_wpxp.
                                   ! Otherwise, <u'w'> and <v'w'> are still
                                   ! approximated by eddy diffusivity when <u>
                                   ! and <v> are advanced in subroutine
                                   ! advance_windm_edsclrm.
      l_diffuse_rtm_and_thlm,    & ! This flag determines whether or not we want
                                   ! CLUBB to do diffusion on rtm and thlm
      l_upwind_xm_ma,            & ! This flag determines whether we want to use
                                   ! an upwind differencing approximation rather
                                   ! than a centered differencing for turbulent
                                   ! or mean advection terms. It affects rtm,
                                   ! thlm, sclrm, um and vm.
      l_tke_aniso,               & ! For anisotropic turbulent kinetic energy,
                                   ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_enable_relaxed_clipping    ! Flag to relax clipping on wpxp in
                                   ! xm_wpxp_clipping_and_stats
      
    integer, intent(in) :: &
      order_xm_wpxp, &
      order_xp2_xpyp, &
      order_wp2_wp3

    ! ------------------- Input/Output Variables -------------------
    
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]
      
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  & 
      um,   & ! <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
      upwp, & ! <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      vm,   & ! <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
      vpwp    ! <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]

    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! ------------------- Local Variables -------------------
    
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,2*nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    ! Additional variables for passive scalars
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: & 
      wpsclrp_forcing    ! <w'sclr'> forcing (momentum levels)  [m/s{un vary}]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      um_tndcy,     & ! <u> forcing term + coriolis (thermo levs)        [m/s^2]
      vm_tndcy,     & ! <v> forcing term + coriolis (thermo levs)        [m/s^2]
      upwp_forcing, & ! <u'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
      vpwp_forcing, & ! <v'w'> extra RHS pressure term (mom levs)        [m^2/s^3]
      upthvp,       & ! <u'thv'> (momentum levels)                       [m/s K]
      vpthvp,       & ! <v'thv'> (momentum levels)                       [m/s K]
      upthlp,       & ! eastward horz turb flux of theta_l (mom levs)    [m/s K]
      vpthlp,       & ! northward horz turb flux of theta_l (mom levs)   [m/s K]
      uprtp,        & ! eastward horz turb flux of tot water (mom levs)  [m/s kg/kg]
      vprtp,        & ! northward horz turb flux of tot water (mom levs) [m/s kg/kg]
      tau_C6_zm       ! Time-scale tau on momentum levels applied to C6 term [s]

    real( kind = core_rknd ), dimension(ngrdcol,2*nz,nrhs) :: & 
      rhs,      & ! Right-hand sides of band diag. matrix. (LAPACK)
      rhs_save, & ! Saved Right-hand sides of band diag. matrix. (LAPACK)
      solution    ! solution vectors of band diag. matrix. (LAPACK)

    ! Constant parameters as a function of Skw.

    real( kind = core_rknd ), dimension(ngrdcol) :: rcond
      
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zeros_vector
      
    integer :: i, k, j
    
    ! ------------------- Begin Code -------------------
    
    ! This is initialized solely for the purpose of avoiding a compiler
    ! warning about uninitialized variables.
    zeros_vector(:,:) = zero
    
    ! Simple case, where the new PDF is
    ! used, l_explicit_turbulent_adv_wpxp is enabled.
        
    ! Create the lhs once
    call xm_wpxp_lhs( nz, ngrdcol, gr, l_iter, dt, zeros_vector, wm_zt, C7_Skw_fnc, & ! In
                      zeros_vector, zeros_vector,                                   & ! In
                      l_implemented, lhs_diff_zm, lhs_diff_zt,                      & ! In
                      lhs_ma_zm, lhs_ma_zt, lhs_ta_wpxp, lhs_ta_xm,                 & ! In
                      lhs_tp, lhs_pr1_wprtp, lhs_ac_pr2,                            & ! In
                      l_diffuse_rtm_and_thlm,                                       & ! In
                      lhs )                                                           ! Out

    ! Compute the explicit portion of the r_t and w'r_t' equations.
    ! Build the right-hand side vector.
    call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_rtm, l_iter, dt, rtm, wprtp,  & ! In
                      rtm_forcing, wprtp_forcing, C7_Skw_fnc,                & ! In
                      rtpthvp, rhs_ta_wprtp, thv_ds_zm,                      & ! In
                      lhs_pr1_wprtp, lhs_ta_wpxp,                            & ! In
                      stats_zt, stats_zm,                                    & ! Inout
                      rhs(:,:,1) )                                             ! Out
                        
    ! Compute the explicit portion of the th_l and w'th_l' equations.
    ! Build the right-hand side vector.
    call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_thlm, l_iter, dt, thlm, wpthlp,  & ! In
                      thlm_forcing, wpthlp_forcing, C7_Skw_fnc,                 & ! In
                      thlpthvp, rhs_ta_wpthlp, thv_ds_zm,                       & ! In
                      lhs_pr1_wpthlp, lhs_ta_wpxp,                              & ! In
                      stats_zt, stats_zm,                                       & ! Inout
                      rhs(:,:,2) )                                                ! Out

! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
    do j = 1, 0, 1
#else
    do j = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

      ! Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
      ! using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
      wpsclrp_forcing(:,:,j) = zero
      
      call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_scalar, l_iter, dt, sclrm(:,:,j), wpsclrp(:,:,j), & ! In
                        sclrm_forcing(:,:,j),                                             & ! In
                        wpsclrp_forcing(:,:,j), C7_Skw_fnc,                               & ! In
                        sclrpthvp(:,:,j), rhs_ta_wpsclrp(:,:,j), thv_ds_zm,               & ! In
                        lhs_pr1_wpsclrp, lhs_ta_wpxp,                                     & ! In
                        stats_zt, stats_zm,                                               & ! Inout
                        rhs(:,:,2+j) )                                                      ! Out
    end do

    if ( l_predict_upwp_vpwp ) then

      ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.
      ! Currently, this requires the ADG1 PDF with implicit turbulent advection.
      ! l_explicit_turbulent_adv_wpxp = false
      ! and ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_new_hybrid )

      ! Coriolis term for <u> and <v>
      if ( .not. l_implemented ) then

        ! Only compute the Coriolis term if the model is running on its own,
        ! and is not part of a larger, host model.
        do k = 1, nz
          do i = 1, ngrdcol
            um_tndcy(i,k) = um_forcing(i,k) - fcor(i) * ( vg(i,k) - vm(i,k) )
            vm_tndcy(i,k) = vm_forcing(i,k) + fcor(i) * ( ug(i,k) - um(i,k) )
          end do
        end do

        if ( l_stats_samp ) then
          do i = 1, ngrdcol
            ! um or vm term gf is completely explicit; call stat_update_var.
            call stat_update_var( ium_gf, - fcor(i) * vg(i,:), & ! intent(in)
                                  stats_zt(i) )             ! intent(inout)
            call stat_update_var( ivm_gf, fcor(i) * ug(i,:), & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)

            ! um or vm term cf is completely explicit; call stat_update_var.
            call stat_update_var( ium_cf, fcor(i) * vm(i,:), & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
            call stat_update_var( ivm_cf, - fcor(i) * um(i,:), & ! intent(in)
                                  stats_zt(i) )             ! intent(inout)

            ! um or vm forcing term
            call stat_update_var( ium_f, um_forcing(i,:), & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
            call stat_update_var( ivm_f, vm_forcing(i,:), & ! intent(in)
                                  stats_zt(i) )           ! intent(inout)
          end do
        endif ! l_stats_samp

      else ! implemented in a host model

        do k = 1, nz
          do i = 1, ngrdcol
            um_tndcy(i,k) = zero
            vm_tndcy(i,k) = zero
          end do
        end do

      end if ! .not. l_implemented

      ! Add "extra term" and optional Coriolis term for <u'w'> and <v'w'>.
      upwp_forcing(:,:) = C_uu_shr * wp2(:,:) * ddzt( nz, ngrdcol, gr(:), um(:,:) )
      vpwp_forcing(:,:) = C_uu_shr * wp2(:,:) * ddzt( nz, ngrdcol, gr(:), vm(:,:) )

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_update_var( iupwp_pr4, C_uu_shr * wp2(i,:) * ddzt( gr(i), um(i,:) ), & ! intent(in)
                                stats_zm(i) )                                ! intent(inout)
          call stat_update_var( ivpwp_pr4, C_uu_shr * wp2(i,:) * ddzt( gr(i), vm(i,:) ), & ! intent(in)
                                stats_zm(i) )                                ! intent(inout)
        end do
      end if ! l_stats_samp

      ! need tau_C6_zm for these calls
      tau_C6_zm(:,:) = min ( one / invrs_tau_C6_zm(:,:), tau_max_zm(:,:) )

      call diagnose_upxp( nz, ngrdcol, gr, upwp, thlm, wpthlp, um,  & ! Intent(in)
                          C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,     & ! Intent(in)
                          upthlp )                                    ! Intent(out)

      call diagnose_upxp( nz, ngrdcol, gr, upwp, rtm, wprtp, um,  & ! Intent(in)
                          C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,    & ! Intent(in)
                          uprtp )                                   ! Intent(out)

      call diagnose_upxp( nz, ngrdcol, gr, vpwp, thlm, wpthlp, vm,  & ! Intent(in)
                          C6thl_Skw_fnc, tau_C6_zm, C7_Skw_fnc,     & ! Intent(in)
                          vpthlp )                                    ! Intent(out)

      call diagnose_upxp( nz, ngrdcol, gr, vpwp, rtm, wprtp, vm,  & ! Intent(in)
                          C6rt_Skw_fnc, tau_C6_zm, C7_Skw_fnc,    & ! Intent(in)
                          vprtp )                                   ! Intent(out)

      ! Use a crude approximation for buoyancy terms <u'thv'> and <v'thv'>.
      !upthvp = upwp * wpthvp / max( wp2, w_tol_sqd )
      !vpthvp = vpwp * wpthvp / max( wp2, w_tol_sqd )
      !upthvp = 0.3_core_rknd * ( upthlp + 200.0_core_rknd * uprtp ) &
      !         + 200._core_rknd * sign( one, upwp) * sqrt( up2 * rcm**2 )
      !vpthvp = 0.3_core_rknd * ( vpthlp + 200.0_core_rknd * vprtp ) &
      !         + 200._core_rknd * sign( one, vpwp ) * sqrt( vp2 * rcm**2 )
      upthvp(:,:) = upthlp(:,:) + ep1 * thv_ds_zm(:,:) * uprtp(:,:) + rc_coef(:,:) * uprcp(:,:)
      vpthvp(:,:) = vpthlp(:,:) + ep1 * thv_ds_zm(:,:) * vprtp(:,:) + rc_coef(:,:) * vprcp(:,:)

      if ( l_stats_samp ) then
        do i = 1, ngrdcol
          call stat_update_var( iupthlp, upthlp(i,:), & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
          call stat_update_var( iuprtp,  uprtp(i,:),  & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
          call stat_update_var( ivpthlp, vpthlp(i,:), & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
          call stat_update_var( ivprtp,  vprtp(i,:),  & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
          call stat_update_var( iupthvp, upthvp(i,:), & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
          call stat_update_var( ivpthvp, vpthvp(i,:), & ! intent(in)
                                stats_zm(i) )         ! intent(inout)
        end do
      end if ! l_stats_samp

      call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_um, l_iter, dt, um, upwp,  & ! In
                        um_tndcy, upwp_forcing, C7_Skw_fnc,                 & ! In
                        upthvp, rhs_ta_wpup, thv_ds_zm,                     & ! In
                        lhs_pr1_wprtp, lhs_ta_wpxp,                         & ! In
                        stats_zt, stats_zm,                                 & ! Inout
                        rhs(:,:,3+sclr_dim) )                                 ! Out

      call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_vm, l_iter, dt, vm, vpwp,  & ! In
                        vm_tndcy, vpwp_forcing, C7_Skw_fnc,                 & ! In
                        vpthvp, rhs_ta_wpvp, thv_ds_zm,                     & ! In
                        lhs_pr1_wprtp, lhs_ta_wpxp,                         & ! In
                        stats_zt, stats_zm,                                 & ! Inout
                        rhs(:,:,4+sclr_dim) )                                 ! Out

    endif ! l_predict_upwp_vpwp

    ! Save the value of rhs, which will be overwritten with the solution as
    ! part of the solving routine.
    rhs_save = rhs

    ! Solve for all fields
    if ( l_stats_samp .and. ithlm_matrix_condt_num + irtm_matrix_condt_num > 0 ) then
       call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &                     ! Intent(in)
                          lhs, rhs, &                 ! Intent(inout)
                          solution, rcond )           ! Intent(out)
    else
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &              ! Intent(in)
                          lhs, rhs, &          ! Intent(inout)
                          solution )           ! Intent(out)
    end if
    

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
         
        write(fstderr,*) "xm & wpxp LU decomp. failed"
        write(fstderr,*) "General xm and wpxp LHS"
        
        do k = 1, nz
          do i = 1, ngrdcol
            write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                             "LHS = ", lhs(1:nsup+nsub+1,i,2*k-1)
            write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                             "LHS = ", lhs(1:nsup+nsub+1,i,2*k)
          end do
        end do ! k = 1, nz
          
        do j = 1, nrhs
          if ( j == 1 ) then
            write(fstderr,*) "rtm and wprtp RHS"
          elseif ( j == 2 ) then
            write(fstderr,*) "thlm and wpthlp RHS"
          else ! j > 2
            if ( sclr_dim > 0 ) then
              if ( j <= 2+sclr_dim ) then
                write(fstderr,*) "sclrm and wpsclrp RHS for sclr", j-2
              end if ! j <= 2+sclr_dim )
            end if ! sclr_dim > 0
            if ( l_predict_upwp_vpwp ) then
              if ( j == 3+sclr_dim ) then
                write(fstderr,*) "um and upwp RHS"
              elseif ( j == 4+sclr_dim ) then
                write(fstderr,*) "vm and vpwp RHS"
              end if
            end if ! l_predict_upwp_vpwp
          end if
          do k = 1, nz
            do i = 1, ngrdcol
              write(fstderr,*) "grid col = ",i,"zt level = ", k, &
                               "height [m] = ", gr(i)%zt(k), &
                               "RHS = ", rhs_save(i,2*k-1,j)
              write(fstderr,*) "grid col = ",i,"zm level = ", k, &
                               "height [m] = ", gr(i)%zm(k), &
                               "RHS = ", rhs_save(i,2*k,j)
            end do
          end do ! k = 1, nz
        end do ! j = 1, nrhs
        return
      end if
    end if
    
    call xm_wpxp_clipping_and_stats( nz, ngrdcol, &   ! Intent(in)
           gr, xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
           rtm_forcing, rho_ds_zm, rho_ds_zt, &       ! Intent(in)
           invrs_rho_ds_zm, invrs_rho_ds_zt, &        ! Intent(in)
           rt_tol**2, rt_tol, rcond, &                ! Intent(in)
           low_lev_effect, high_lev_effect, &         ! Intent(in)
           lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp, &       ! Intent(in)
           lhs_diff_zm, C7_Skw_fnc, &                 ! Intent(in)
           lhs_tp, lhs_ta_xm, lhs_pr1_wprtp, &        ! Intent(in)
           l_implemented, solution(:,:,1),  &         ! Intent(in)
           l_predict_upwp_vpwp, &                     ! Intent(in)
           l_upwind_xm_ma, &                          ! Intent(in)
           l_tke_aniso, &                             ! Intent(in)
           l_enable_relaxed_clipping, &               ! Intent(in)
           order_xm_wpxp, order_xp2_xpyp, &           ! Intent(in)
           order_wp2_wp3, &                           ! Intent(in)
           stats_zt, stats_zm, stats_sfc, &           ! intent(inout)
           rtm, rt_tol_mfl, wprtp )                   ! Intent(inout)

    if ( clubb_at_least_debug_level( 0 ) ) then
       if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "rtm monotonic flux limiter:  tridag failed"
          return
       end if
    end if

    call xm_wpxp_clipping_and_stats( nz, ngrdcol, &   ! Intent(in)
           gr, xm_wpxp_thlm, dt, wp2, thlp2, wm_zt, & ! Intent(in)
           thlm_forcing, rho_ds_zm, rho_ds_zt, &      ! Intent(in)
           invrs_rho_ds_zm, invrs_rho_ds_zt, &        ! Intent(in)
           thl_tol**2, thl_tol, rcond, &              ! Intent(in)
           low_lev_effect, high_lev_effect, &         ! Intent(in)
           lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp, &       ! Intent(in)
           lhs_diff_zm, C7_Skw_fnc, &                 ! Intent(in)
           lhs_tp, lhs_ta_xm, lhs_pr1_wprtp, &        ! Intent(in)
           l_implemented, solution(:,:,2),  &         ! Intent(in)
           l_predict_upwp_vpwp, &                     ! Intent(in)
           l_upwind_xm_ma, &                          ! Intent(in)
           l_tke_aniso, &                             ! Intent(in)
           l_enable_relaxed_clipping, &               ! Intent(in)
           order_xm_wpxp, order_xp2_xpyp, &           ! Intent(in)
           order_wp2_wp3, &                           ! Intent(in)
           stats_zt, stats_zm, stats_sfc, &           ! intent(inout)
           thlm, thl_tol_mfl, wpthlp )                ! Intent(inout)

    if ( clubb_at_least_debug_level( 0 ) ) then
       if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "thlm monotonic flux limiter:  tridag failed"
          return
       end if
    end if

! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
    do j = 1, 0, 1
#else
    do j = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15
      call xm_wpxp_clipping_and_stats( nz, ngrdcol, &               ! Intent(in)
             gr, xm_wpxp_scalar, dt, wp2, sclrp2(:,:,j), wm_zt, &   ! Intent(in)
             sclrm_forcing(:,:,j), &                                ! Intent(in)
             rho_ds_zm, rho_ds_zt, &                                ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &                    ! Intent(in)
             sclr_tol(j)**2, sclr_tol(j), rcond, &                  ! Intent(in)
             low_lev_effect, high_lev_effect, &                     ! Intent(in)
             lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp, &                   ! Intent(in)
             lhs_diff_zm, C7_Skw_fnc, &                             ! Intent(in)
             lhs_tp, lhs_ta_xm, lhs_pr1_wprtp, &                    ! Intent(in)
             l_implemented, solution(:,:,2+j),  &                   ! Intent(in)
             l_predict_upwp_vpwp, &                                 ! Intent(in)
             l_upwind_xm_ma, &                                      ! Intent(in)
             l_tke_aniso, &                                         ! Intent(in)
             l_enable_relaxed_clipping, &                           ! Intent(in)
             order_xm_wpxp, order_xp2_xpyp, &                       ! Intent(in)
             order_wp2_wp3, &                                       ! Intent(in)
             stats_zt, stats_zm, stats_sfc, &                       ! intent(inout)
             sclrm(:,:,j), sclr_tol(j), wpsclrp(:,:,j) )            ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "sclrm # ", j, "monotonic flux limiter: tridag failed"
            return
         end if
      end if

    end do ! 1..sclr_dim

    if ( l_predict_upwp_vpwp ) then

      ! Predict <u> and <u'w'>, as well as <v> and <v'w'>.
      call xm_wpxp_clipping_and_stats( nz, ngrdcol,   & ! Intent(in)
            gr, xm_wpxp_um, dt, wp2, up2, wm_zt,      & ! Intent(in)
            um_tndcy, rho_ds_zm, rho_ds_zt,           & ! Intent(in)
            invrs_rho_ds_zm, invrs_rho_ds_zt,         & ! Intent(in)
            w_tol_sqd, w_tol, rcond,                  & ! Intent(in)
            low_lev_effect, high_lev_effect,          & ! Intent(in)
            lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,        & ! Intent(in)
            lhs_diff_zm, C7_Skw_fnc,                  & ! Intent(in)
            lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,         & ! Intent(in)
            l_implemented, solution(:,:,3+sclr_dim),  & ! Intent(in)
            l_predict_upwp_vpwp,                      & ! Intent(in)
            l_upwind_xm_ma,                           & ! Intent(in)
            l_tke_aniso,                              & ! Intent(in)
            l_enable_relaxed_clipping,                & ! Intent(in)
            order_xm_wpxp, order_xp2_xpyp,            & ! Intent(in)
            order_wp2_wp3,                            & ! Intent(in)
            stats_zt, stats_zm, stats_sfc,            & ! intent(inout)
            um, w_tol, upwp                           ) ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "um monotonic flux limiter:  tridag failed"
          return
        end if
      end if

      call xm_wpxp_clipping_and_stats( nz, ngrdcol,   & ! Intent(in)
            gr, xm_wpxp_vm, dt, wp2, vp2, wm_zt,      & ! Intent(in)
            vm_tndcy, rho_ds_zm, rho_ds_zt,           & ! Intent(in)
            invrs_rho_ds_zm, invrs_rho_ds_zt,         & ! Intent(in)
            w_tol_sqd, w_tol, rcond,                  & ! Intent(in)
            low_lev_effect, high_lev_effect,          & ! Intent(in)
            lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp,        & ! Intent(in)
            lhs_diff_zm, C7_Skw_fnc,                  & ! Intent(in)
            lhs_tp, lhs_ta_xm, lhs_pr1_wprtp,         & ! Intent(in)
            l_implemented, solution(:,:,4+sclr_dim),  & ! Intent(in)
            l_predict_upwp_vpwp,                      & ! Intent(in)
            l_upwind_xm_ma,                           & ! Intent(in)
            l_tke_aniso,                              & ! Intent(in)
            l_enable_relaxed_clipping,                & ! Intent(in)
            order_xm_wpxp, order_xp2_xpyp,            & ! Intent(in)
            order_wp2_wp3,                            & ! Intent(in)
            stats_zt, stats_zm, stats_sfc,            & ! intent(inout)
            vm, w_tol, vpwp )                           ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "vm monotonic flux limiter:  tridag failed"
          return
        end if
      end if

    end if ! l_predict_upwp_vpwp
    
  end subroutine solve_xm_wpxp_with_single_lhs
  
  !==========================================================================================

  subroutine solve_xm_wpxp_with_multiple_lhs( nz, ngrdcol, gr, dt, l_iter, nrhs, wm_zt, wp2, &
                                            rtpthvp, rtm_forcing, wprtp_forcing, thlpthvp, &
                                            thlm_forcing,   wpthlp_forcing, rho_ds_zm, &
                                            rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                            thv_ds_zm, rtp2, thlp2, l_implemented, &
                                            sclrpthvp, sclrm_forcing, sclrp2, &
                                            low_lev_effect, high_lev_effect, C7_Skw_fnc, &
                                            lhs_diff_zm, lhs_diff_zt, lhs_ma_zt, lhs_ma_zm, &
                                            lhs_ta_wprtp, lhs_ta_wpthlp, lhs_ta_wpsclrp, &
                                            rhs_ta_wprtp, rhs_ta_wpthlp, rhs_ta_wpsclrp, &
                                            lhs_tp, lhs_ta_xm, lhs_ac_pr2, lhs_pr1_wprtp, &
                                            lhs_pr1_wpthlp, lhs_pr1_wpsclrp, &
                                            l_predict_upwp_vpwp, &
                                            l_diffuse_rtm_and_thlm, &
                                            l_upwind_xm_ma, &
                                            l_tke_aniso, &
                                            l_enable_relaxed_clipping, &
                                            order_xm_wpxp, order_xp2_xpyp, order_wp2_wp3, &
                                            stats_zt, stats_zm, stats_sfc, & 
                                            rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp )
    !            
    ! Description: This subroutine solves all xm_wpxp when all the LHS matrices are NOT equal.
    !              This means multiple solves are required, one for each unique LHS.
    !
    !----------------------------------------------------------------------------------------
    
    use grid_class, only: & 
        grid, & ! Type
        ddzt    ! Procedure(s)
      
    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constants
      
    use stats_type_utilities, only: & 
        stat_update_var   ! Procedure(s)
      
    use stats_variables, only: & 
        irtm_matrix_condt_num, &  ! Variables
        ithlm_matrix_condt_num, &
        l_stats_samp
        
    use parameters_model, only: & 
        sclr_dim, &  ! Variable(s)
        sclr_tol
        
    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use constants_clubb, only:  & 
        fstderr, &  ! Constant
        rt_tol, &
        thl_tol, &
        thl_tol_mfl, &
        rt_tol_mfl, &
        zero

    use stats_type, only: stats ! Type

    implicit none
    
    ! ------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr
    
    real( kind = core_rknd ), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: & 
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
      thlpthvp,        & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,    & ! th_l forcing (thermodynamic levels)      [K/s]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)      [K/s^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on moment. levs. [K]
      rtp2,            & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2              ! th_l'^2 (momentum levels)                [K^2]

    logical, intent(in) ::  & 
      l_implemented, &      ! Flag for CLUBB being implemented in a larger model.
      l_iter

    ! Additional variables for passive scalars
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz,sclr_dim) :: & 
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    ! LHS/RHS terms
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_diff_zt,  & ! Diffusion term for w'x'
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm       ! Mean advection contributions to lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_ta_wprtp,   & ! w'r_t' turbulent advection contributions to lhs
      lhs_ta_wpthlp     ! w'thl' turbulent advection contributions to lhs
      
    real( kind = core_rknd ), dimension(3,ngrdcol,nz,sclr_dim), intent(in) :: & 
      lhs_ta_wpsclrp    ! w'sclr' turbulent advection contributions to lhs  
     
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      rhs_ta_wprtp,  & ! w'r_t' turbulent advection contributions to rhs  
      rhs_ta_wpthlp    ! w'thl' turbulent advection contributions to rhs
      
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim), intent(in) :: & 
      rhs_ta_wpsclrp    ! w'sclr' turbulent advection contributions to rhs

    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(in) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      lhs_ac_pr2,     & ! Accumulation of w'x' and w'x' pressure term 2
      lhs_pr1_wprtp,  & ! Pressure term 1 for w'r_t' for all grid levels
      lhs_pr1_wpthlp, & ! Pressure term 1 for w'thl' for all grid levels
      lhs_pr1_wpsclrp   ! Pressure term 1 for w'sclr' for all grid levels
      
    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(ngrdcol,nz), intent(in) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      C7_Skw_fnc

    integer, intent(in) :: &
      nrhs         ! Number of RHS vectors

    logical, intent(in) :: &
      l_predict_upwp_vpwp,       & ! Flag to predict <u'w'> and <v'w'> along
                                   ! with <u> and <v> alongside the advancement
                                   ! of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>,
                                   ! and <w'sclr'> in subroutine advance_xm_wpxp.
                                   ! Otherwise, <u'w'> and <v'w'> are still
                                   ! approximated by eddy diffusivity when <u>
                                   ! and <v> are advanced in subroutine
                                   ! advance_windm_edsclrm.
      l_diffuse_rtm_and_thlm,    & ! This flag determines whether or not we want
                                   ! CLUBB to do diffusion on rtm and thlm
      l_upwind_xm_ma,            & ! This flag determines whether we want to use
                                   ! an upwind differencing approximation rather
                                   ! than a centered differencing for turbulent
                                   ! or mean advection terms. It affects rtm,
                                   ! thlm, sclrm, um and vm.
      l_tke_aniso,               & ! For anisotropic turbulent kinetic energy,
                                   ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_enable_relaxed_clipping    ! Flag to relax clipping on wpxp in
                                   ! xm_wpxp_clipping_and_stats
      
    integer, intent(in) :: &
      order_xm_wpxp, &
      order_xp2_xpyp, &
      order_wp2_wp3

    ! ------------------- Input/Output Variables -------------------

    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc
    
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]
      
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! ------------------- Local Variables -------------------
    
    real( kind = core_rknd ), dimension(nsup+nsub+1,ngrdcol,2*nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    real( kind = core_rknd ), dimension(ngrdcol,2*nz,nrhs) :: & 
      rhs,      & ! Right-hand sides of band diag. matrix. (LAPACK)
      rhs_save, & ! Saved Right-hand sides of band diag. matrix. (LAPACK)
      solution    ! solution vectors of band diag. matrix. (LAPACK)
      
    ! Additional variables for passive scalars
    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: & 
      wpsclrp_forcing    ! <w'sclr'> forcing (momentum levels)  [m/s{un vary}]

    ! Variables used for clipping of w'x' due to correlation
    ! of w with x, such that:
    ! corr_(w,x) = w'x' / [ sqrt(w'^2) * sqrt(x'^2) ];
    ! -1 <= corr_(w,x) <= 1.
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      wpxp_upper_lim, & ! Keeps correlations from becoming greater than 1.
      wpxp_lower_lim    ! Keeps correlations from becoming less than -1.

    ! Constant parameters as a function of Skw.

    real( kind = core_rknd ), dimension(ngrdcol) :: rcond
      
    integer :: i, j, k
      
    ! ------------------- Begin Code -------------------

    ! Compute the implicit portion of the r_t and w'r_t' equations.
    ! Build the left-hand side matrix.                 
    call xm_wpxp_lhs( nz, ngrdcol, gr, l_iter, dt, wprtp, wm_zt, C7_Skw_fnc,  & ! In
                      wpxp_upper_lim, wpxp_lower_lim,                         & ! In
                      l_implemented, lhs_diff_zm, lhs_diff_zt,                & ! In
                      lhs_ma_zm, lhs_ma_zt, lhs_ta_wprtp, lhs_ta_xm,          & ! In
                      lhs_tp, lhs_pr1_wprtp, lhs_ac_pr2,                      & ! In
                      l_diffuse_rtm_and_thlm,                                 & ! In
                      lhs )                                                     ! Out

    ! Compute the explicit portion of the r_t and w'r_t' equations.
    ! Build the right-hand side vector.
    call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_rtm, l_iter, dt, rtm, wprtp,   & ! In
                      rtm_forcing, wprtp_forcing, C7_Skw_fnc,                 & ! In
                      rtpthvp, rhs_ta_wprtp, thv_ds_zm,                       & ! In
                      lhs_pr1_wprtp, lhs_ta_wprtp,                            & ! In
                      stats_zt, stats_zm,                                     & ! Inout
                      rhs(:,:,1) )                                              ! Out

    ! Save the value of rhs, which will be overwritten with the solution as
    ! part of the solving routine.
    rhs_save = rhs

    ! Solve r_t / w'r_t'
    if ( l_stats_samp .and. irtm_matrix_condt_num > 0 ) then
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &                     ! Intent(in)
                          lhs, rhs, &                 ! Intent(inout)
                          solution, rcond )           ! Intent(out)
    else
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &              ! Intent(in)
                          lhs, rhs, &          ! Intent(inout)
                          solution )           ! Intent(out)
    end if

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        do i = 1, ngrdcol
          write(fstderr,*) "Mean total water & total water flux LU decomp. failed"
          write(fstderr,*) "rtm and wprtp LHS"
          do k = 1, nz
            write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                             "LHS = ", lhs(1:nsup+nsub+1,i,2*k-1)
            write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                             "LHS = ", lhs(1:nsup+nsub+1,i,2*k)
          end do ! k = 1, nz
          write(fstderr,*) "rtm and wprtp RHS"
          do k = 1, nz
            write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                             "RHS = ", rhs_save(i,2*k-1,1)
            write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                             "RHS = ", rhs_save(i,2*k,1)
          end do ! k = 1, nz
        end do
        return
      end if
    end if

    call xm_wpxp_clipping_and_stats( nz, ngrdcol, &   ! Intent(in)
           gr, xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
           rtm_forcing, rho_ds_zm, rho_ds_zt, &       ! Intent(in)
           invrs_rho_ds_zm, invrs_rho_ds_zt, &        ! Intent(in)
           rt_tol**2, rt_tol, rcond, &                ! Intent(in)
           low_lev_effect, high_lev_effect, &         ! Intent(in)
           lhs_ma_zt, lhs_ma_zm, lhs_ta_wprtp, &      ! Intent(in)
           lhs_diff_zm, C7_Skw_fnc, &                 ! Intent(in)
           lhs_tp, lhs_ta_xm, lhs_pr1_wprtp, &        ! Intent(in)
           l_implemented, solution(:,:,1), &          ! Intent(in)
           l_predict_upwp_vpwp, &                     ! Intent(in)
           l_upwind_xm_ma, &                          ! Intent(in)
           l_tke_aniso, &                             ! Intent(in)
           l_enable_relaxed_clipping, &               ! Intent(in)
           order_xm_wpxp, order_xp2_xpyp, &           ! Intent(in)
           order_wp2_wp3, &                           ! Intent(in)
           stats_zt, stats_zm, stats_sfc, &           ! intent(inout)
           rtm, rt_tol_mfl, wprtp )                   ! Intent(inout)

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "rtm monotonic flux limiter:  tridag failed"
        return
      end if
    end if
      
    ! Compute the implicit portion of the th_l and w'th_l' equations.
    ! Build the left-hand side matrix.
    call xm_wpxp_lhs( nz, ngrdcol, gr, l_iter, dt, wpthlp, wm_zt, C7_Skw_fnc,   & ! In
                      wpxp_upper_lim, wpxp_lower_lim,                           & ! In
                      l_implemented, lhs_diff_zm, lhs_diff_zt,                  & ! In
                      lhs_ma_zm, lhs_ma_zt, lhs_ta_wpthlp, lhs_ta_xm,           & ! In
                      lhs_tp, lhs_pr1_wpthlp, lhs_ac_pr2,                       & ! In
                      l_diffuse_rtm_and_thlm,                                   & ! In
                      lhs )                                                       ! Out

    ! Compute the explicit portion of the th_l and w'th_l' equations.
    ! Build the right-hand side vector.
    call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_thlm, l_iter, dt, thlm, wpthlp,  & ! In
                      thlm_forcing, wpthlp_forcing, C7_Skw_fnc,                 & ! In
                      thlpthvp, rhs_ta_wpthlp, thv_ds_zm,                       & ! In
                      lhs_pr1_wpthlp, lhs_ta_wpthlp,                            & ! In
                      stats_zt, stats_zm,                                       & ! Inout
                      rhs(:,:,1) )                                                ! Out

    ! Save the value of rhs, which will be overwritten with the solution as
    ! part of the solving routine.
    rhs_save = rhs

    ! Solve for th_l / w'th_l'
    if ( l_stats_samp .and. ithlm_matrix_condt_num > 0 ) then
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &                 ! Intent(in)
                          lhs, rhs, &                 ! Intent(inout)
                          solution, rcond )           ! Intent(out)
    else
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &          ! Intent(in)
                          lhs, rhs, &          ! Intent(inout)
                          solution )           ! Intent(out)
    end if

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        do i = 1, ngrdcol
          write(fstderr,*) "Liquid pot. temp & thetal flux LU decomp. failed"
          write(fstderr,*) "thlm and wpthlp LHS"
          do k = 1, nz
             write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                              "LHS = ", lhs(1:nsup+nsub+1,i,2*k-1)
             write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                              "LHS = ", lhs(1:nsup+nsub+1,i,2*k)
          end do ! k = 1, nz
          write(fstderr,*) "thlm and wpthlp RHS"
          do k = 1, nz
             write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                              "RHS = ", rhs_save(i,2*k-1,1)
             write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                              "RHS = ", rhs_save(i,2*k,1)
          end do ! k = 1, nz
        end do
        return
      end if
    end if

    call xm_wpxp_clipping_and_stats( nz, ngrdcol, &     ! Intent(in)
           gr, xm_wpxp_thlm, dt, wp2, thlp2, wm_zt,  &  ! Intent(in)
           thlm_forcing, rho_ds_zm, rho_ds_zt, &        ! Intent(in)
           invrs_rho_ds_zm, invrs_rho_ds_zt, &          ! Intent(in)
           thl_tol**2, thl_tol, rcond, &                ! Intent(in)
           low_lev_effect, high_lev_effect, &           ! Intent(in)
           lhs_ma_zt, lhs_ma_zm, lhs_ta_wpthlp, &       ! Intent(in)
           lhs_diff_zm, C7_Skw_fnc, &                   ! Intent(in)
           lhs_tp, lhs_ta_xm, lhs_pr1_wpthlp, &         ! Intent(in)
           l_implemented, solution(:,:,1),  &           ! Intent(in)
           l_predict_upwp_vpwp, &                       ! Intent(in)
           l_upwind_xm_ma, &                            ! Intent(in)
           l_tke_aniso, &                               ! Intent(in)
           l_enable_relaxed_clipping, &                 ! Intent(in)
           order_xm_wpxp, order_xp2_xpyp, &             ! Intent(in)
           order_wp2_wp3, &                             ! Intent(in)
           stats_zt, stats_zm, stats_sfc, &             ! intent(inout)
           thlm, thl_tol_mfl, wpthlp )                  ! Intent(inout)

    if ( clubb_at_least_debug_level( 0 ) ) then
      if ( err_code == clubb_fatal_error ) then
        write(fstderr,*) "thlm monotonic flux limiter:  tridag failed" 
        return
      end if
    end if

    ! Solve sclrm / wpsclrp
    ! If sclr_dim is 0, then this loop will execute 0 times.
! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
    do j = 1, 0, 1
#else
    do j = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

      ! Set <w'sclr'> forcing to 0 unless unless testing the wpsclrp code
      ! using wprtp or wpthlp (then use wprtp_forcing or wpthlp_forcing).
      wpsclrp_forcing(:,:,j) = zero
      
      ! Compute the implicit portion of the sclr and w'sclr' equations.
      ! Build the left-hand side matrix.
      call xm_wpxp_lhs( nz, ngrdcol, gr, l_iter, dt, wpsclrp(:,:,j), wm_zt, C7_Skw_fnc,   & ! In
                        wpxp_upper_lim, wpxp_lower_lim,                                   & ! In
                        l_implemented, lhs_diff_zm, lhs_diff_zt,                          & ! In
                        lhs_ma_zm, lhs_ma_zt, lhs_ta_wpsclrp(:,:,:,j), lhs_ta_xm,         & ! In
                        lhs_tp, lhs_pr1_wpsclrp, lhs_ac_pr2,                              & ! In
                        l_diffuse_rtm_and_thlm,                                           & ! In
                        lhs )                                                               ! Out

      ! Compute the explicit portion of the sclrm and w'sclr' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( nz, ngrdcol, gr, xm_wpxp_scalar, l_iter, dt, sclrm(:,:,j), wpsclrp(:,:,j),            & ! In
                        sclrm_forcing(:,:,j),                                             & ! In
                        wpsclrp_forcing(:,:,j), C7_Skw_fnc,                               & ! In
                        sclrpthvp(:,:,j), rhs_ta_wpsclrp(:,:,j), thv_ds_zm,               & ! In
                        lhs_pr1_wpsclrp, lhs_ta_wpsclrp(:,:,:,j),                         & ! In
                        stats_zt, stats_zm,                                               & ! Inout
                        rhs(:,:,1) )                                                        ! Out

      ! Save the value of rhs, which will be overwritten with the solution as
      ! part of the solving routine.
      rhs_save = rhs

      ! Solve for sclrm / w'sclr'
      call xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &              ! Intent(in)
                          lhs, rhs, &          ! Intent(inout)
                          solution )           ! Intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then   
          do i = 1, ngrdcol
            write(fstderr,*) "Passive scalar # ", j, " LU decomp. failed."
            write(fstderr,*) "sclrm and wpsclrp LHS"
            do k = 1, nz
               write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,i,2*k-1)
               write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                                "LHS = ", lhs(1:nsup+nsub+1,i,2*k)
            end do ! k = 1, nz
            write(fstderr,*) "sclrm and wpsclrp RHS"
            do k = 1, nz
               write(fstderr,*) "grid col = ",i,"zt level = ", k, "height [m] = ", gr(i)%zt(k), &
                                "RHS = ", rhs_save(i,2*k-1,1)
               write(fstderr,*) "grid col = ",i,"zm level = ", k, "height [m] = ", gr(i)%zm(k), &
                                "RHS = ", rhs_save(i,2*k,1)
            end do ! k = 1, nz
          end do
          return
        end if
      end if
      
      call xm_wpxp_clipping_and_stats( nz, ngrdcol, &               ! Intent(in)
             gr, xm_wpxp_scalar, dt, wp2, sclrp2(:,:,j), wm_zt,   & ! Intent(in)
             sclrm_forcing(:,:,j),  &                               ! Intent(in)
             rho_ds_zm, rho_ds_zt, &                                ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &                    ! Intent(in)
             sclr_tol(j)**2, sclr_tol(j), rcond, &                  ! Intent(in)
             low_lev_effect, high_lev_effect, &                     ! Intent(in)
             lhs_ma_zt, lhs_ma_zm, lhs_ta_wpsclrp(:,:,:,j), &       ! Intent(in)
             lhs_diff_zm, C7_Skw_fnc, &                             ! Intent(in)
             lhs_tp, lhs_ta_xm, lhs_pr1_wpsclrp, &                  ! Intent(in)
             l_implemented, solution(:,:,1),  &                     ! Intent(in)
             l_predict_upwp_vpwp, &                                 ! Intent(in)
             l_upwind_xm_ma, &                                      ! Intent(in)
             l_tke_aniso, &                                         ! Intent(in)
             l_enable_relaxed_clipping, &                           ! Intent(in)
             order_xm_wpxp, order_xp2_xpyp, &                       ! Intent(in)
             order_wp2_wp3, &                                       ! Intent(in)
             stats_zt, stats_zm, stats_sfc, &                       ! intent(inout)
             sclrm(:,:,j), sclr_tol(j), wpsclrp(:,:,j) )            ! Intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "sclrm # ", j, "monotonic flux limiter: tridag failed"
          return
        end if
      end if

    end do ! passive scalars
    
  end subroutine solve_xm_wpxp_with_multiple_lhs

  !=============================================================================
  subroutine xm_wpxp_solve( nz, ngrdcol, gr, nrhs, &
                            lhs, rhs, &
                            solution, rcond )

    ! Description:
    !   Solve for xm / w'x' using the band diagonal solver.

    ! References:
    !   None
    !------------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

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

    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    integer, intent(in) :: &
      nrhs ! Number of rhs vectors

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nsup+nsub+1,ngrdcol,2*nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix in LAPACK storage)

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,2*nz,nrhs) ::  & 
      rhs      ! Right-hand side of band diag. matrix. (LAPACK storage)

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,2*nz,nrhs) ::  & 
      solution ! Solution to band diagonal system (LAPACK storage)

    ! Output Variables
    real( kind = core_rknd ), optional, dimension(ngrdcol), intent(out) :: &
      rcond ! Est. of the reciprocal of the condition #
      
    ! Local Variables
    integer :: i

    if ( present( rcond ) ) then
      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      do i = 1, ngrdcol
        call band_solvex( "xm_wpxp", nsup, nsub, 2*nz, nrhs, & ! intent(in) 
                          lhs(:,i,:), rhs(i,:,:), & ! intent(inout)
                          solution(i,:,:), rcond(i) ) ! intent(out)
      end do
    else
      ! Perform LU decomp and solve system (LAPACK)
      do i = 1, ngrdcol
        call band_solve( "xm_wpxp", nsup, nsub, 2*nz, nrhs, lhs(:,i,:), & ! intent(in) 
                         rhs(i,:,:), & ! intent(inout)
                         solution(i,:,:) ) ! intent(out)
      end do
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
  subroutine xm_wpxp_clipping_and_stats( &
               nz, ngrdcol, gr, solve_type, dt, wp2, xp2, wm_zt, &
               xm_forcing, rho_ds_zm, rho_ds_zt, &
               invrs_rho_ds_zm, invrs_rho_ds_zt, &
               xp2_threshold, xm_threshold, rcond, &
               low_lev_effect, high_lev_effect, &
               lhs_ma_zt, lhs_ma_zm, lhs_ta_wpxp, &
               lhs_diff_zm, C7_Skw_fnc, &
               lhs_tp, lhs_ta_xm, lhs_pr1, &
               l_implemented, solution, &
               l_predict_upwp_vpwp, &
               l_upwind_xm_ma, &
               l_tke_aniso, &
               l_enable_relaxed_clipping, &
               order_xm_wpxp, order_xp2_xpyp, &
               order_wp2_wp3, &
               stats_zt, stats_zm, stats_sfc, & 
               xm, xm_tol, wpxp )

    ! Description:
    ! Clips and computes implicit stats for an artitrary xm and wpxp
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

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
        eps, &
        gamma_over_implicit_ts

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
        l_stats_samp, & 
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

    use stats_type, only: stats ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Constant Parameters
    logical, parameter :: &
      l_mono_flux_lim = .true.!, &  ! Flag for monotonic turbulent flux limiter
!      l_first_clip_ts = .true., &
!      l_last_clip_ts  = .false.

    logical :: &
      l_first_clip_ts, &
      l_last_clip_ts
      
    integer, intent(in) ::  & 
      solve_type  ! Variables being solved for.

    real( kind = core_rknd ), intent(in) ::  & 
      dt  ! Timestep   [s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) ::  & 
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
      xm_tol           ! Minimum allowable deviation of xm [units vary]
      
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      rcond ! Reciprocal of the estimated condition number (from computing A^-1)

    ! Variables used as part of the monotonic turbulent advection scheme.
    ! Find the lowermost and uppermost grid levels that can have an effect
    ! on the central thermodynamic level during the course of a time step,
    ! due to the effects of turbulent advection only.
    integer, dimension(ngrdcol,nz), intent(in) ::  &
      low_lev_effect, & ! Index of the lowest level that has an effect.
      high_lev_effect   ! Index of the highest level that has an effect.
    
    real( kind = core_rknd ), dimension(3,ngrdcol,nz), intent(in) :: & 
      lhs_diff_zm,  & ! Diffusion term for w'x'
      lhs_ma_zt,    & ! Mean advection contributions to lhs
      lhs_ma_zm,    & ! Mean advection contributions to lhs
      lhs_ta_wpxp     ! Turbulent advection contributions to lhs
      
    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(in) :: & 
      lhs_tp,     & ! Turbulent production terms of w'x'
      lhs_ta_xm     ! Turbulent advection terms of xm
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      lhs_pr1       ! Pressure term 1 for w'x'
      
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      C7_Skw_fnc
      
    logical, intent(in) :: &
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,2*nz) :: &
      solution ! The <t+1> value of xm and wpxp   [units vary]

    logical, intent(in) :: &
      l_predict_upwp_vpwp,       & ! Flag to predict <u'w'> and <v'w'> along
                                   ! with <u> and <v> alongside the advancement
                                   ! of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>,
                                   ! and <w'sclr'> in subroutine advance_xm_wpxp.
                                   ! Otherwise, <u'w'> and <v'w'> are still
                                   ! approximated by eddy diffusivity when <u>
                                   ! and <v> are advanced in subroutine
                                   ! advance_windm_edsclrm.
      l_upwind_xm_ma,            & ! This flag determines whether we want to use
                                   ! an upwind differencing approximation rather
                                   ! than a centered differencing for turbulent
                                   ! or mean advection terms. It affects rtm,
                                   ! thlm, sclrm, um and vm.
      l_tke_aniso,               & ! For anisotropic turbulent kinetic energy,
                                   ! i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_enable_relaxed_clipping    ! Flag to relax clipping on wpxp in
                                   ! xm_wpxp_clipping_and_stats

    integer, intent(in) :: &
      order_xm_wpxp, &
      order_xp2_xpyp, &
      order_wp2_wp3

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) :: & 
      xm, &     ! The mean x field  [units vary]
      wpxp      ! The flux of x     [units vary m/s]

    type (stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! Local Variables
    integer :: & 
      solve_type_cl ! solve_type used for clipping statistics.

    character(len=10) :: &
      solve_type_str ! solve_type as a string for debug output purposes

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      xm_n ! Old value of xm for positive definite scheme     [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      wpxp_pd, xm_pd ! Change in xm and wpxp due to the pos. def. scheme

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wpxp_chnge, &  ! Net change in w'x' due to clipping       [units vary]
      xp2_relaxed    ! Value of x'^2 * clip_factor               [units vary]
      
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
      zero_vector, &
      wpxp_ac, &
      wpxp_pr2

    ! Indices
    integer :: &
      k, i, km1, kp1, &
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

    do k=1, nz, 1

      k_xm   = 2 * k - 1
      k_wpxp = 2 * k
      
      do i = 1, ngrdcol

        xm_n(i,k) = xm(i,k)

        xm(i,k)   = solution(i,k_xm)
        wpxp(i,k) = solution(i,k_wpxp)
        
      end do

    end do ! k=1..nz

    ! Lower boundary condition on xm
    do i = 1, ngrdcol
      xm(i,1) = xm(i,2)
    end do


    if ( l_stats_samp ) then
      
      do i = 1, ngrdcol

        if ( ixm_matrix_condt_num > 0 ) then
          ! Est. of the condition number of the mean/flux LHS matrix
          call stat_update_var_pt( ixm_matrix_condt_num, 1, one / rcond(i), & ! intent(in)
                                   stats_sfc(i) )                             ! intent(inout)
        end if

        ! The xm loop runs between k = 2 and k = nz.  The value of xm at
        ! level k = 1, which is below the model surface, is simply set equal to
        ! the value of xm at level k = 2 after the solve has been completed.
        ! Thus, the statistical code will run from levels 2 through nz.

        do k = 2, nz

          km1 = max( k-1, 1 )
          kp1 = min( k+1, nz )

          ! Finalize implicit contributions for xm

          ! xm term ma is completely implicit; call stat_update_var_pt.
          if ( .not. l_implemented ) then
            call stat_update_var_pt( ixm_ma, k, & ! intent(in)
                (-lhs_ma_zt(3,i,k)) * xm(i,km1) & 
              + (-lhs_ma_zt(2,i,k)) * xm(i,k) & 
              + (-lhs_ma_zt(1,i,k)) * xm(i,kp1), & ! intent(in)
                stats_zt(i) ) ! intent(inout)
          end if

          ! xm term ta is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( ixm_ta, k, & ! intent(in)
              (-lhs_ta_xm(2,i,k)) * wpxp(i,km1) & 
            + (-lhs_ta_xm(1,i,k)) * wpxp(i,k), & ! intent(in)
              stats_zt(i) ) ! intent(inout)

        enddo ! xm loop: 2..nz

        ! The wpxp loop runs between k = 2 and k = nz-1.  The value of wpxp
        ! is set to specified values at both the lowest level, k = 1, and the
        ! highest level, k = nz.  Thus, the statistical code will run from
        ! levels 2 through nz-1.
        
        zero_vector(i,:) = 0.0_core_rknd
        
        ! Note:  To find the contribution of w'x' term ac,
        !        substitute 0 for the C_7 skewness function input
        !        to function wpxp_terms_ac_pr2_lhs.
        call wpxp_terms_ac_pr2_lhs( gr(i), zero_vector(i,:), wm_zt(i,:), gr(i)%invrs_dzm, & ! intent(in)
                                    wpxp_ac(i,:) ) ! intent(out)

        ! Note:  To find the contribution of w'x' term pr2,
        !        add 1 to the C_7 skewness function input
        !        to function wpxp_terms_ac_pr2_lhs.
        call wpxp_terms_ac_pr2_lhs( gr(i), (one+C7_Skw_fnc(i,:)), wm_zt(i,:), gr(i)%invrs_dzm, & ! intent(in)
                                    wpxp_pr2(i,:) ) ! intent(out)

        do k = 2, nz-1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, nz )

          ! Finalize implicit contributions for wpxp

          ! w'x' term ma is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( iwpxp_ma, k, & ! intent(in)
              (-lhs_ma_zm(3,i,k)) * wpxp(i,km1) & 
            + (-lhs_ma_zm(2,i,k)) * wpxp(i,k) & 
            + (-lhs_ma_zm(1,i,k)) * wpxp(i,kp1), & ! intent(in)
              stats_zm(i) ) ! intent(inout)


            call stat_end_update_pt( iwpxp_ta, k, & ! intent(in)
                (-gamma_over_implicit_ts*lhs_ta_wpxp(3,i,k)) * wpxp(i,km1) & 
              + (-gamma_over_implicit_ts*lhs_ta_wpxp(2,i,k)) * wpxp(i,k) & 
              + (-gamma_over_implicit_ts*lhs_ta_wpxp(1,i,k)) * wpxp(i,kp1), & ! intent(in)
                stats_zm(i) ) ! intent(inout)

          ! w'x' term tp is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( iwpxp_tp, k, & ! intent(in)
              (-lhs_tp(2,i,k)) * xm(i,k) & 
            + (-lhs_tp(1,i,k)) * xm(i,kp1), & ! intent(in)
              stats_zm(i) ) ! intent(inout)

          ! w'x' term ac is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( iwpxp_ac, k, & ! intent(in)
                                   -wpxp_ac(i,k) * wpxp(i,k), & ! intent(in)
                                   stats_zm(i) ) ! intent(inout)

          ! w'x' term pr1 is normally completely implicit.  However, due to the
          ! RHS contribution from the "over-implicit" weighted time step,
          ! w'x' term pr1 has both implicit and explicit components;
          ! call stat_end_update_pt.
          ! Note:  An "over-implicit" weighted time step is applied to this term.
          !        A weighting factor of greater than 1 may be used to make the
          !        term more numerically stable (see note above for LHS turbulent
          !        advection (ta) term).
          call stat_end_update_pt( iwpxp_pr1, k, & ! intent(in) 
              (-gamma_over_implicit_ts*lhs_pr1(i,k)) * wpxp(i,k), & ! intent(in)
              stats_zm(i) ) ! intent(inout)
              
          call stat_update_var_pt( iwpxp_pr2, k, & ! intent(in) 
                                  -wpxp_pr2(i,k) * wpxp(i,k), & ! intent(in)
                                   stats_zm(i) ) ! intent(inout)

          ! w'x' term dp1 is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( iwpxp_dp1, k, & ! intent(in)
              (-lhs_diff_zm(3,i,k)) * wpxp(i,km1) & 
            + (-lhs_diff_zm(2,i,k)) * wpxp(i,k) & 
            + (-lhs_diff_zm(1,i,k)) * wpxp(i,kp1), & ! intent(in)
              stats_zm(i) ) ! intent(inout)

        end do ! wpxp loop: 2..nz-1
        
      end do

    end if ! l_stats_samp


    ! Apply a monotonic turbulent flux limiter to xm/w'x'.
    if ( l_mono_flux_lim ) then
      call monotonic_turbulent_flux_limit( nz, ngrdcol, gr, solve_type, dt, xm_n, & ! intent(in)
                                           xp2, wm_zt, xm_forcing, & ! intent(in)
                                           rho_ds_zm, rho_ds_zt, & ! intent(in)
                                           invrs_rho_ds_zm, invrs_rho_ds_zt, & ! intent(in)
                                           xp2_threshold, xm_tol, l_implemented, & ! intent(in)
                                           low_lev_effect, high_lev_effect, & ! intent(in)
                                           l_upwind_xm_ma, & ! intent(in)
                                           stats_zt, stats_zm, & ! intent(inout)
                                           xm, wpxp ) ! intent(inout)
    end if ! l_mono_flux_lim

    ! Apply a flux limiting positive definite scheme if the solution
    ! for the mean field is negative and we're determining total water
    if ( solve_type == xm_wpxp_rtm .and. l_pos_def ) then
      do i = 1, ngrdcol
        if ( any( xm(i,:) < zero ) ) then
            call pos_definite_adj( gr(i), dt, "zt", xm(i,:), & ! intent(in) 
                                   wpxp(i,:), xm_n(i,:), & ! intent(inout)
                                   xm_pd(i,:), wpxp_pd(i,:) ) ! intent(out)
        end if
      end do
    else
      ! For stats purposes
      xm_pd(:,:)   = zero
      wpxp_pd(:,:) = zero

    end if ! l_pos_def and solve_type == "rtm" and rtm <n+1> less than 0

    if ( l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_update_var( iwpxp_pd, wpxp_pd(i,1:nz), & ! intent(in)
                              stats_zm(i) )                ! intent(inout)

        call stat_update_var( ixm_pd, xm_pd(i,1:nz), & ! intent(in)
                              stats_zt(i) )            ! intent(inout)
                          
        ! Computed value before clipping    
        call stat_begin_update( gr(i), ixm_cl, xm(i,:) / dt, & ! Intent(in)
                                stats_zt(i) )                  ! Intent(inout)
      end do
    end if
    
    if ( solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm .and. l_hole_fill ) then 

      do i = 1, ngrdcol
        
        if ( any( xm(i,:) < xm_threshold) ) then

          if ( clubb_at_least_debug_level( 3 ) ) then
            
            select case ( solve_type )
            case ( xm_wpxp_rtm )
              solve_type_str = "rtm"
            case ( xm_wpxp_thlm )
              solve_type_str = "thlm"
            case default
              solve_type_str = "scalars"
            end select
            
            do k = 1, nz
              if ( xm(i,k) < zero ) then
                write(fstderr,*) solve_type_str//" < ", xm_threshold, &
                  " in advance_xm_wpxp_module at k= ", k
              end if
            end do
          end if

          call fill_holes_vertical( gr(i), 2, xm_threshold, "zt", & ! intent(in)
                                    rho_ds_zt(i,:), rho_ds_zm(i,:), & ! intent(in)
                                    xm(i,:) ) ! intent(inout)

          ! Hole filling does not affect the below ground level, perform a blunt clipping
          ! here on that level to prevent small values of xm(1)
          xm(i,1) = max( xm(i,1), xm_tol )

        endif ! any( xm < xm_threshold ) .and. l_hole_fill
              ! .and. solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm
      end do
      
    end if

    if ( l_stats_samp ) then
      do i = 1, ngrdcol
        call stat_end_update( gr(i), ixm_cl, xm(i,:) / dt, & ! Intent(in) 
                              stats_zt(i) )                       ! Intent(inout)
      end do                        
    end if

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
        xp2_relaxed(:,:) = max( 1e-7_core_rknd , xp2(:,:) )

      else if ( solve_type == xm_wpxp_thlm ) then
        xp2_relaxed(:,:) = max( 0.01_core_rknd, xp2(:,:) )

      else ! This includes the passive scalars
        xp2_relaxed(:,:) = max( 1e-7_core_rknd , xp2(:,:) )

      end if

    else  ! Don't relax clipping
      xp2_relaxed(:,:) = xp2(:,:)

    end if

    if ( order_xm_wpxp < order_wp2_wp3 .and. order_xm_wpxp < order_xp2_xpyp ) then
       l_first_clip_ts = .true.
       l_last_clip_ts = .false.
    elseif ( order_xm_wpxp > order_wp2_wp3 .and. order_xm_wpxp > order_xp2_xpyp ) then
       l_first_clip_ts = .false.
       l_last_clip_ts = .true.
    else
       l_first_clip_ts = .false.
       l_last_clip_ts = .false.
    endif
    
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

    if ( solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm ) then
      do i = 1, ngrdcol
        call clip_covar( gr(i), solve_type_cl, l_first_clip_ts, &  ! In
                         l_last_clip_ts, dt, wp2(i,:), xp2_relaxed(i,:), &  ! In
                         l_predict_upwp_vpwp, & ! In
                         stats_zm(i), & ! intent(inout)
                         wpxp(i,:), wpxp_chnge(i,:) ) ! In/Out
      end do
    else ! clipping for upwp or vpwp

      if ( l_tke_aniso ) then
        do i = 1, ngrdcol
          call clip_covar( gr(i), solve_type_cl, l_first_clip_ts, &  ! In
                           l_last_clip_ts, dt, wp2(i,:), xp2(i,:), &  ! In
                           l_predict_upwp_vpwp, & ! In
                           stats_zm(i), & ! intent(inout)
                           wpxp(i,:), wpxp_chnge(i,:) ) ! In/Out
        end do
      else
        do i = 1, ngrdcol
          call clip_covar( gr(i), solve_type_cl, l_first_clip_ts, &  ! In
                           l_last_clip_ts, dt, wp2(i,:), wp2(i,:), &  ! In
                           l_predict_upwp_vpwp, & ! In
                           stats_zm(i), & ! intent(inout)
                           wpxp(i,:), wpxp_chnge(i,:) ) ! In/Out
        end do
       end if ! l_tke_aniso
    end if ! solve_type /= xm_wpxp_um .and. solve_type /= xm_wpxp_vm

    ! Adjusting xm based on clipping for w'x'.
    if ( l_clip_turb_adv ) then
      do i = 1, ngrdcol
        if ( any( abs(wpxp_chnge(i,:)) > eps ) ) then
          call xm_correction_wpxp_cl( gr(i), solve_type, dt, wpxp_chnge(i,:), gr(i)%invrs_dzt, & ! intent(in)
                                      stats_zt(i), & ! intent(inout)
                                      xm(i,:) ) ! intent(inout)
        endif
      end do
    end if

    return
  end subroutine xm_wpxp_clipping_and_stats

  !=============================================================================
  pure subroutine xm_term_ta_lhs( nz, ngrdcol, gr, &
                                  rho_ds_zm, invrs_rho_ds_zt, &
                                  lhs_ta_xm )

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
    ! =====rho_ds_zm=====wpxp================================== m(k-1)
    !
    ! The vertical indices m(k), t(k), and m(k-1) correspond with altitudes
    ! zm(k), zt(k), and zm(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: &
        zero    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      k_mdiag   = 1,    & ! Momentum superdiagonal index.
      km1_mdiag = 2       ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      rho_ds_zm,       & ! Dry, static density at momentum levels     [kg/m^3]
      invrs_rho_ds_zt    ! Inverse dry, static density at thermo levs [m^3/kg]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(out) :: &
      lhs_ta_xm    ! LHS coefficient of xm turbulent advection  [1/m]

    ! Local Variable
    integer :: i, k    ! Vertical level index


    ! Set lower boundary condition to 0
    do i = 1, ngrdcol
      lhs_ta_xm(k_mdiag,i,1)   = zero
      lhs_ta_xm(km1_mdiag,i,1) = zero
    end do

    ! Calculate term at all other grid levels.
    do k = 2, nz 
      do i = 1, ngrdcol

        ! Momentum superdiagonal [ x wpxp(k,<t+1>) ]
        lhs_ta_xm(k_mdiag,i,k) = + invrs_rho_ds_zt(i,k) * gr(i)%invrs_dzt(k) * rho_ds_zm(i,k)

        ! Momentum subdiagonal [ x wpxp(k-1,<t+1>) ]
        lhs_ta_xm(km1_mdiag,i,k) = - invrs_rho_ds_zt(i,k) * gr(i)%invrs_dzt(k) * rho_ds_zm(i,k-1)
      end do
    end do ! k = 2, gr%nz 

    return

  end subroutine xm_term_ta_lhs

  !=============================================================================
  pure subroutine wpxp_term_tp_lhs( nz, ngrdcol, gr, wp2, & 
                                    lhs_tp )

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
    ! ---------------------------xm---------------------------- t(k+1)
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

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: &
        zero    ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag = 2         ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      wp2  ! w'^2                       [m^2/s^2]

    ! Return Variable
    real( kind = core_rknd ), dimension(2,ngrdcol,nz), intent(out) :: &
        lhs_tp    ! LHS coefficient of xm for turbulent production  [1/s]

    ! Local Variable
    integer :: i, k  ! Vertical level index


    ! Set lower boundary to 0
    do i = 1, ngrdcol
      lhs_tp(1,i,1) = zero
      lhs_tp(2,i,1) = zero
    end do

    ! Calculate term at all interior grid levels.
    do k = 2, nz-1
      do i = 1, ngrdcol

       ! Thermodynamic superdiagonal [ x xm(k+1,<t+1>) ]
       lhs_tp(kp1_tdiag,i,k) = + wp2(i,k) * gr(i)%invrs_dzm(k)

       ! Thermodynamic subdiagonal [ x xm(k,<t+1>) ]
       lhs_tp(k_tdiag,i,k)   = - wp2(i,k) * gr(i)%invrs_dzm(k)
       
      end do
    end do ! k = 2, gr%nz-1

    ! Set upper boundary to 0
    do i = 1, ngrdcol
      lhs_tp(1,i,nz) = 0.0_core_rknd
      lhs_tp(2,i,nz) = 0.0_core_rknd
    end do

    return

  end subroutine wpxp_term_tp_lhs
    
  !=============================================================================
  pure subroutine wpxp_terms_ac_pr2_lhs( gr, C7_Skw_fnc, wm_zt, invrs_dzm, &
                                         lhs_ac_pr2  ) 

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
    ! -------wm_zt--------------------------------------------- t(k+1)
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

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied              [-]
      wm_zt,       & ! w wind component on thermodynamic levels     [m/s]
      invrs_dzm      ! Inverse of grid spacing                      [1/m]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      lhs_ac_pr2    ! LHS coefficient of accumulation and pressure term 2 [1/s]

    ! Local Variable
    integer :: k    ! Vertical level index


    ! Set lower boundary to 0
    lhs_ac_pr2(1) = zero

    ! Calculate term at all interior grid levels.
    do k = 2, gr%nz-1

       ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
       lhs_ac_pr2(k) &
       = ( one - C7_Skw_fnc(k) ) * invrs_dzm(k) * ( wm_zt(k+1) - wm_zt(k) )

    enddo ! k = 2, gr%nz-1 

    ! Set upper boundary to 0
    lhs_ac_pr2(gr%nz) = zero


    return

  end subroutine wpxp_terms_ac_pr2_lhs

  !=============================================================================
  pure subroutine wpxp_term_pr1_lhs( nz, ngrdcol, C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc, &
                                     invrs_tau_C6_zm, l_scalar_calc, &
                                     lhs_pr1_wprtp, lhs_pr1_wpthlp, &
                                     lhs_pr1_wpsclrp )

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
    !
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      C6rt_Skw_fnc,  & ! C_6rt parameter with Sk_w applied        [-]
      C6thl_Skw_fnc, & ! C_6thl parameter with Sk_w applied       [-]
      C7_Skw_fnc,    & ! C_7 parameter with Sk_w applied          [-]
      invrs_tau_C6_zm  ! Inverse time-scale tau at momentum levels   [1/s]

    logical, intent(in) :: &
      l_scalar_calc   ! True if sclr_dim > 0

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: & 
      lhs_pr1_wprtp,   & ! LHS coefficient for w'r_t' pressure term 1   [1/s]
      lhs_pr1_wpthlp,  & ! LHS coefficient for w'thl' pressure term 1   [1/s]
      lhs_pr1_wpsclrp    ! LHS coefficient for w'sclr' pressure term 1  [1/s]

    ! Set lower boundary to 0
    lhs_pr1_wprtp(:,1) = zero

    ! Momentum main diagonals: [ x wpxp(k,<t+1>) ]
    lhs_pr1_wprtp(:,2:nz-1) = C6rt_Skw_fnc(:,2:nz-1) * invrs_tau_C6_zm(:,2:nz-1)

    ! Set upper boundary to 0
    lhs_pr1_wprtp(:,nz) = zero
    
    ! Set lower boundary to 0
    lhs_pr1_wpthlp(:,1) = zero

    ! Momentum main diagonals: [ x wpxp(k,<t+1>) ]
    lhs_pr1_wpthlp(:,2:nz-1) = C6thl_Skw_fnc(:,2:nz-1) * invrs_tau_C6_zm(:,2:nz-1)

    ! Set upper boundary to 0
    lhs_pr1_wpthlp(:,nz) = zero
        
    if ( l_scalar_calc ) then

       ! Set lower boundary to 0
       lhs_pr1_wpsclrp(:,1) = zero

       ! Momentum main diagonals: [ x wpxp(k,<t+1>) ]
       lhs_pr1_wpsclrp(:,2:nz-1) = C7_Skw_fnc(:,2:nz-1) * invrs_tau_C6_zm(:,2:nz-1)

       ! Set upper boundary to 0
       lhs_pr1_wpsclrp(:,nz) = zero

    endif ! l_scalar_calc

    return

  end subroutine wpxp_term_pr1_lhs

  !=============================================================================
  pure subroutine wpxp_terms_bp_pr3_rhs( nz, ngrdcol, gr, C7_Skw_fnc, thv_ds_zm, xpthvp, &
                                         rhs_bp_pr3 )

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

    use grid_class, only:  & 
        grid ! Type

    use constants_clubb, only: & ! Constants(s) 
        grav, & ! Gravitational acceleration [m/s^2]
        one,  &
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: & 
      C7_Skw_fnc,  & ! C_7 parameter with Sk_w applied       [-]
      thv_ds_zm,   & ! Dry, base-state theta_v on mom. levs. [K]
      xpthvp         ! x'th_v'                               [K {xm units}]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
        rhs_bp_pr3    ! RHS portion of bouyancy prod and pressure term 3

    ! Local Variable
    integer :: i, k    ! Vertical level index
    

    ! Set lower boundary to 0
    do i = 1, ngrdcol
      rhs_bp_pr3(i,1) = zero
    end do

    ! Calculate term at all interior grid levels.
    do k = 2, nz-1
      do i = 1, ngrdcol
        rhs_bp_pr3(i,k) = ( grav / thv_ds_zm(i,k) ) * ( one - C7_Skw_fnc(i,k) ) * xpthvp(i,k)
      end do
    end do ! k = 2, nz-1

    ! Set upper boundary to 0
    do i = 1, ngrdcol
      rhs_bp_pr3(i,nz) = zero
    end do

    return

  end subroutine wpxp_terms_bp_pr3_rhs

  !=============================================================================
  subroutine xm_correction_wpxp_cl( gr, solve_type, dt, wpxp_chnge, invrs_dzt, &
                                    stats_zt, & 
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
    ! =======wpxp_amount_clipped=============================== m(k-1)
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
        grid ! Type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        l_stats_samp, & ! Variable(s)
        ithlm_tacl, &
        irtm_tacl

    use stats_type, only: stats ! Type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_zt

    type (grid), target, intent(in) :: gr

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
      call stat_update_var( ixm_tacl, xm_tndcy_wpxp_cl, & ! intent(in)
                            stats_zt )                    ! intent(inout)
    endif


    return

  end subroutine xm_correction_wpxp_cl


  !=============================================================================
  subroutine damp_coefficient( nz, ngrdcol, gr, coefficient, Cx_Skw_fnc, &
                               max_coeff_value, altitude_threshold, &
                               threshold, Lscale, &
                               damped_value )

    ! Description:
    ! Damps a given coefficient linearly based on the value of Lscale.
    ! For additional information see CLUBB ticket #431.

    use grid_class, only: & 
        grid ! Type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input variables
    real( kind = core_rknd ), intent(in) :: &
      coefficient,        & ! The coefficient to be damped
      max_coeff_value,    & ! Maximum value the damped coefficient should have
      altitude_threshold, & ! Minimum altitude where damping should occur 
      threshold             ! Value of Lscale below which the damping should occur

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Lscale,           &   ! Current value of Lscale
      Cx_Skw_fnc            ! Initial skewness function before damping

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: damped_value
    
    ! Local Variables
    integer :: i, k
    
    do k = 1, nz
      do i = 1, ngrdcol
        
        if ( Lscale(i,k) < threshold .and. gr(i)%zt(k) > altitude_threshold ) then
          damped_value(i,k) = max_coeff_value &
                              + ( ( coefficient - max_coeff_value ) / threshold ) &
                                * Lscale(i,k)
        else
          damped_value(i,k) = Cx_Skw_fnc(i,k)
        end if
        
      end do
    end do
    
    return

  end subroutine damp_coefficient
  !-----------------------------------------------------------------------

  !=====================================================================================
  pure subroutine diagnose_upxp( nz, ngrdcol, gr, ypwp, xm, wpxp, ym, &
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
        grid, & ! Type
        ddzt  ! Procedure

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      ypwp,        & ! momentum flux component, either upwp or vpwp  [m^2/s^2]
      xm,          & ! grid-mean conserved thermodynamic variable, either thlm or rtm [varies]
      wpxp,        & ! vertical scalar flux, either wpthlp or wprtp [varies]
      ym,          & ! grid-mean velocity component, either um or vm [m/s]
      C6x_Skw_fnc, & ! C_6 pressure parameter with effects of Sk_w incorporated (k)  [-]
      tau_C6_zm,   & ! Time-scale tau on momentum levels applied to C6 term [s]
      C7_Skw_fnc     ! C_7 pressure parameter with effects of Sk_w incorporated (k)  [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
        ypxp        ! horizontal flux of a conserved scalar, either upthlp, uprtp, vpthlp, or vprtp

    
      ypxp = ( tau_C6_zm / C6x_Skw_fnc ) &
             * ( - ypwp * ddzt( nz, ngrdcol, gr, xm ) &
                 - (one - C7_Skw_fnc ) * ( wpxp * ddzt( nz, ngrdcol, gr, ym ) ) )
              
    return

  end subroutine diagnose_upxp

  !=============================================================================
  subroutine error_prints_xm_wpxp( gr, dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                                   Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                                   Kh_zt, Kh_zm, invrs_tau_C6_zm, Skw_zm, &
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
                                   um_old, upwp_old, vm_old, vpwp_old, &
                                   l_predict_upwp_vpwp, l_lmm_stepping )

    ! Description:
    ! Prints values of model fields when fatal errors (LU decomp.) occur.
    ! All field that are passed into and out of subroutine advance_xm_wpxp are
    ! printed.  If additional fields are added to the call to subroutine
    ! advance_xm_wpxp, they should also be added here.

    use constants_clubb, only: &
        fstderr    ! Variable(s)

    use grid_class, only: &
        grid ! Type

    use parameters_model, only: &
        sclr_dim    ! Variable(s)

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

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
      invrs_tau_C6_zm, & ! Inverse time-scale tau on momentum levels applied to C6 term [1/s]
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

    type(implicit_coefs_terms), intent(in) :: &
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
      sclrm_old,   & ! Saved value of sclrm     [units vary]
      wpsclrp_old    ! Saved value of wpsclrp   [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      um_old,   & ! Saved value of <u>       [m/s]
      upwp_old, & ! Saved value of <u'w'>    [m^2/s^2]
      vm_old,   & ! Saved value of <v>       [m/s]
      vpwp_old    ! Saved value of <v'w'>    [m^2/s^2]

    logical, intent(in) :: &
      l_predict_upwp_vpwp, & ! Flag to predict <u'w'> and <v'w'> along with <u>
                             ! and <v> alongside the advancement of <rt>,
                             ! <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'>
                             ! in subroutine advance_xm_wpxp.  Otherwise, <u'w'>
                             ! and <v'w'> are still approximated by eddy
                             ! diffusivity when <u> and <v> are advanced in
                             ! subroutine advance_windm_edsclrm.
      l_lmm_stepping         ! Apply Linear Multistep Method (LMM) Stepping


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
    write(fstderr,*) "invrs_tau_C6_zm = ", invrs_tau_C6_zm, new_line('c')
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
    write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp2rtp_implicit = ", &
                     pdf_implicit_coefs_terms%coef_wp2rtp_implicit, &
                     new_line('c')
    write(fstderr,*) "pdf_implicit_coefs_terms%term_wp2rtp_explicit = ", &
                     pdf_implicit_coefs_terms%term_wp2rtp_explicit, &
                     new_line('c')
    write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp2thlp_implicit = ", &
                     pdf_implicit_coefs_terms%coef_wp2thlp_implicit, &
                     new_line('c')
    write(fstderr,*) "pdf_implicit_coefs_terms%term_wp2thlp_explicit = ", &
                     pdf_implicit_coefs_terms%term_wp2thlp_explicit, &
                     new_line('c')
     
    if ( sclr_dim > 0 )  then
       write(fstderr,*) "sclrp2 = ", sclrp2, new_line('c')
       write(fstderr,*) "wp2sclrp = ", wp2sclrp, new_line('c')
       write(fstderr,*) "sclrpthvp = ", sclrpthvp, new_line('c')
       write(fstderr,*) "sclrm_forcing = ", sclrm_forcing, new_line('c')
       write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp2sclrp_implicit = ", &
                        pdf_implicit_coefs_terms%coef_wp2sclrp_implicit, &
                        new_line('c')
       write(fstderr,*) "pdf_implicit_coefs_terms%term_wp2sclrp_explicit = ", &
                        pdf_implicit_coefs_terms%term_wp2sclrp_explicit, &
                        new_line('c')
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
       write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp2up_implicit = ", &
                        pdf_implicit_coefs_terms%coef_wp2up_implicit, &
                        new_line('c')
       write(fstderr,*) "pdf_implicit_coefs_terms%term_wp2up_explicit = ", &
                        pdf_implicit_coefs_terms%term_wp2up_explicit, &
                        new_line('c')
       write(fstderr,*) "pdf_implicit_coefs_terms%coef_wp2vp_implicit = ", &
                        pdf_implicit_coefs_terms%coef_wp2vp_implicit, &
                        new_line('c')
       write(fstderr,*) "pdf_implicit_coefs_terms%term_wp2vp_explicit = ", &
                        pdf_implicit_coefs_terms%term_wp2vp_explicit, &
                        new_line('c')
    endif ! l_predict_upwp_vpwp

    write(fstderr,*) "Intent(inout)", new_line('c')
     
    if ( l_lmm_stepping ) &
       write(fstderr,*) "rtm (pre-solve) = ", rtm_old, new_line('c')
    write(fstderr,*) "rtm = ", rtm, new_line('c')
    if ( l_lmm_stepping )  &
       write(fstderr,*) "wprtp (pre-solve) = ", wprtp_old, new_line('c')
    write(fstderr,*) "wprtp = ", wprtp, new_line('c')
    if ( l_lmm_stepping ) &
       write(fstderr,*) "thlm (pre-solve) = ", thlm_old, new_line('c')
    write(fstderr,*) "thlm = ", thlm, new_line('c')
    if ( l_lmm_stepping ) &
       write(fstderr,*) "wpthlp (pre-solve) =", wpthlp_old, new_line('c')
    write(fstderr,*) "wpthlp =", wpthlp, new_line('c')

    if ( sclr_dim > 0 )  then
       if ( l_lmm_stepping ) &
          write(fstderr,*) "sclrm (pre-solve) = ", sclrm_old, new_line('c')
       write(fstderr,*) "sclrm = ", sclrm, new_line('c')
       if ( l_lmm_stepping ) &
          write(fstderr,*) "wpsclrp (pre-solve) = ", wpsclrp_old, new_line('c')
       write(fstderr,*) "wpsclrp = ", wpsclrp, new_line('c')
    endif

    if ( l_predict_upwp_vpwp ) then
       if ( l_lmm_stepping ) &
          write(fstderr,*) "um (pre-solve) = ", um_old, new_line('c')
       write(fstderr,*) "um = ", um, new_line('c')
       if ( l_lmm_stepping ) &
          write(fstderr,*) "upwp (pre-solve) = ",  upwp_old, new_line('c')
       write(fstderr,*) "upwp = ",  upwp, new_line('c')
       if ( l_lmm_stepping ) &
          write(fstderr,*) "vm (pre-solve) = ", vm_old, new_line('c')
       write(fstderr,*) "vm = ", vm, new_line('c')
       if ( l_lmm_stepping ) &
          write(fstderr,*) "vpwp (pre-solve) = ",  vpwp_old, new_line('c')
       write(fstderr,*) "vpwp = ",  vpwp, new_line('c')
    endif ! l_predict_upwp_vpwp

    return

  end subroutine error_prints_xm_wpxp

  !=============================================================================

end module advance_xm_wpxp_module
