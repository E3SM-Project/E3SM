!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module model_flags

! Description:
!   Various model options that can be toggled off and on as desired.

! References:
!   None
!-------------------------------------------------------------------------------

  implicit none

  public :: setup_model_flags, clubb_config_flags_type, set_default_clubb_config_flags, &
            initialize_clubb_config_flags_type, print_clubb_config_flags

  private ! Default Scope

  ! Options for the two component normal (double Gaussian) PDF type to use for
  ! the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
  ! two-component PDF.
  integer, parameter, public :: &
    iiPDF_ADG1 = 1,       & ! ADG1 PDF
    iiPDF_ADG2 = 2,       & ! ADG2 PDF
    iiPDF_3D_Luhar = 3,   & ! 3D Luhar PDF
    iiPDF_new = 4,        & ! new PDF
    iiPDF_TSDADG = 5,     & ! new TSDADG PDF
    iiPDF_LY93 = 6,       & ! Lewellen and Yoh (1993)
    iiPDF_new_hybrid = 7    ! new hybrid PDF

  ! Options for the placement of the call to CLUBB's PDF.
  integer, parameter, public :: &
    ipdf_pre_advance_fields   = 1,   & ! Call before advancing predictive fields
    ipdf_post_advance_fields  = 2,   & ! Call after advancing predictive fields
    ipdf_pre_post_advance_fields = 3   ! Call both before and after advancing
                                       ! predictive fields

  logical, parameter, public ::  & 
    l_pos_def            = .false., & ! Flux limiting positive definite scheme on rtm
    l_hole_fill          = .true.,  & ! Hole filling pos def scheme on wp2,up2,rtp2,etc
    l_clip_turb_adv      = .false., & ! Corrects thlm/rtm when w'th_l'/w'r_t' is clipped
    l_sat_mixrat_lookup  = .false.    ! Use a lookup table for mixing length
                                      ! saturation vapor pressure calculations

  logical, parameter, public :: &
#ifdef BYTESWAP_IO
    l_byteswap_io = .true.,   & ! Don't use the native byte ordering in GrADS output
#else
    l_byteswap_io = .false.,  & ! Use the native byte ordering in GrADS output
#endif
    l_gamma_Skw   = .true.      ! Use a Skw dependent gamma parameter

  logical, parameter, public :: &
    l_use_boussinesq = .false.  ! Flag to use the Boussinesq form of the
                                ! predictive equations.  The predictive
                                ! equations are anelastic by default.

  logical, parameter, public :: &
    ! Flag to use explicit turbulent advection in the wp3 predictive equation.
    l_explicit_turbulent_adv_wp3 = .false.,  &
    ! Flag to use explicit turbulent advection in the wpxp predictive equation.
    l_explicit_turbulent_adv_wpxp = .false., &
    ! Flag to use explicit turbulent advection in the xp2 and xpyp predictive
    ! equations.
    l_explicit_turbulent_adv_xpyp = .false.

  ! Flag to advance xp3 using a simplified version of the d(xp3)/dt predictive
  ! equation or calculate it using a steady-state approximation.  When the flag
  ! is turned off, the Larson and Golaz (2005) ansatz to calculate xp3 after
  ! calculating Skx using the ansatz.
  logical, parameter, public :: &
    l_advance_xp3 = .false.

  logical, parameter, public :: &
    l_morr_xp2_mc = .false. ! Flag to include the effects of rain evaporation
                            ! on rtp2 and thlp2.  The moister (rt_1 or rt_2)
                            ! and colder (thl_1 or thl_2) will be fed into
                            ! the morrison microphys, and rain evaporation will
                            ! be allowed to increase variances

  logical, parameter, public :: &
    l_evaporate_cold_rcm = .false. ! Flag to evaporate cloud water at temperatures
                                   ! colder than -37C.  This is to be used for 
                                   ! Morrison microphysics, to prevent excess ice

  logical, parameter, public :: &
    l_cubic_interp = .false.      ! Flag to convert grid points with cubic monotonic
                                  ! spline interpolation as opposed to linear interpolation.

  logical, parameter, public :: &
    l_upwind_Kh_dp_term = .false.

  ! These are the integer constants that represent the various saturation
  ! formulas. To add a new formula, add an additional constant here,
  ! add the logic to check the strings for the new formula in clubb_core and
  ! this module, and add logic in saturation to call the proper function--
  ! the control logic will be based on these named constants.

  integer, parameter, public :: &
    saturation_bolton = 1, & ! Constant for Bolton approximations of saturation
    saturation_gfdl   = 2, & ! Constant for the GFDL approximation of saturation
    saturation_flatau = 3    ! Constant for Flatau approximations of saturation

  !-----------------------------------------------------------------------------
  ! Options that can be changed at runtime 
  ! The default values are chosen below and overwritten if desired by the user
  !----------------------------------------------------------------------------- 

  ! Flag to use high accuracy for the parabolic cylinder function
  logical, public :: &
    l_high_accuracy_parab_cyl_fnc = .false. 

!$omp threadprivate(l_high_accuracy_parab_cyl_fnc)

  logical, public :: & 
    l_quintic_poly_interp = .false. ! Use a quintic polynomial in mono_cubic_interp

!$omp threadprivate(l_quintic_poly_interp)

  ! Use to determine whether a host model has already applied the surface flux,
  ! to avoid double counting.
  logical, public :: &
    l_host_applies_sfc_fluxes = .false.

!$omp threadprivate(l_host_applies_sfc_fluxes)

  integer, public :: &
    saturation_formula = saturation_flatau ! Integer that stores the saturation formula to be used

!$omp threadprivate(saturation_formula)

  logical, parameter, public :: &
    l_silhs_rad = .false.    ! Resolve radiation over subcolumns using SILHS

#ifdef GFDL
  logical, public :: &
     I_sat_sphum       ! h1g, 2010-06-15
!$omp threadprivate( I_sat_sphum )
#endif

  ! Derived type to hold all configurable CLUBB flags
  type clubb_config_flags_type
    ! Note that all flags in this data type might also be modified in 
    ! corresponding .in files which would supercede the default values!
    
    integer :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term, &! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Peskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp in
                                      ! xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds           ! Code to linearize PBL winds

  end type clubb_config_flags_type

  contains

!===============================================================================
  subroutine setup_model_flags & 
             ( l_host_applies_sfc_fluxes_in, & 
               saturation_formula_in &
#ifdef GFDL
               ,  I_sat_sphum_in   &  ! h1g, 2010-06-15
#endif
                )

! Description:
!   Setup flags that influence the numerics, etc. of CLUBB core

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    logical, intent(in) ::  & 
      l_host_applies_sfc_fluxes_in

    character(len=*), intent(in) :: &
      saturation_formula_in

#ifdef GFDL
    logical, intent(in) ::  & 
      I_sat_sphum_in           ! h1g, 2010-06-15
#endif

    !---- Begin Code ----

    ! Logicals

    l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes_in

    l_high_accuracy_parab_cyl_fnc = .false.

    ! Integers

    ! Set up the saturation formula value
    select case ( trim( saturation_formula_in ) )
    case ( "bolton", "Bolton" )
      saturation_formula = saturation_bolton

    case ( "flatau", "Flatau" )
      saturation_formula = saturation_flatau

    case ( "gfdl", "GFDL" )
      saturation_formula = saturation_gfdl

      ! Add new saturation formulas after this.
    end select

#ifdef GFDL
    I_sat_sphum = I_sat_sphum_in  ! h1g, 2010-06-15
#endif
    return
  end subroutine setup_model_flags

!===============================================================================
  subroutine set_default_clubb_config_flags( iiPDF_type, &
                                             ipdf_call_placement, &
                                             l_use_precip_frac, &
                                             l_predict_upwp_vpwp, &
                                             l_min_wp2_from_corr_wx, &
                                             l_min_xp2_from_corr_wx, &
                                             l_C2_cloud_frac, &
                                             l_diffuse_rtm_and_thlm, &
                                             l_stability_correct_Kh_N2_zm, &
                                             l_calc_thlp2_rad, &
                                             l_upwind_xpyp_ta, &
                                             l_upwind_xm_ma, &
                                             l_uv_nudge, &
                                             l_rtm_nudge, &
                                             l_tke_aniso, &
                                             l_vert_avg_closure, &
                                             l_trapezoidal_rule_zt, &
                                             l_trapezoidal_rule_zm, &
                                             l_call_pdf_closure_twice, &
                                             l_standard_term_ta, &
                                             l_partial_upwind_wp3, &
                                             l_godunov_upwind_wpxp_ta, &
                                             l_godunov_upwind_xpyp_ta, &
                                             l_use_cloud_cover, &
                                             l_diagnose_correlations, &
                                             l_calc_w_corr, &
                                             l_const_Nc_in_cloud, &
                                             l_fix_w_chi_eta_correlations, &
                                             l_stability_correct_tau_zm, &
                                             l_damp_wp2_using_em, &
                                             l_do_expldiff_rtm_thlm, &
                                             l_Lscale_plume_centered, &
                                             l_diag_Lscale_from_tau, &
                                             l_use_C7_Richardson, &
                                             l_use_C11_Richardson, &
                                             l_use_shear_Richardson, &
                                             l_brunt_vaisala_freq_moist, &
                                             l_use_thvm_in_bv_freq, &
                                             l_rcm_supersat_adj, &
                                             l_damp_wp3_Skw_squared, &
                                             l_prescribed_avg_deltaz, &
                                             l_lmm_stepping, &
                                             l_e3sm_config, &
                                             l_vary_convect_depth, &
                                             l_use_tke_in_wp3_pr_turb_term, &
                                             l_use_tke_in_wp2_wp3_K_dfsn, &
                                             l_smooth_Heaviside_tau_wpxp, &
                                             l_enable_relaxed_clipping, &
                                             l_linearize_pbl_winds )

! Description:
!   Sets all CLUBB flags to a default setting.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Output variables
    integer, intent(out) :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical, intent(out) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Peskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp in
                                      ! xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds           ! Code to linearize PBL winds

!-----------------------------------------------------------------------
    ! Begin code
    ! WARNING: THE DEFAULT VALUES OF THE FLAGS BELOW MAY BE OVERWRITTEN
    !    BY NAMELIST VALUES FROM, E.G., configurable_clubb_flags_nl!!!

    iiPDF_type = iiPDF_ADG1
    ipdf_call_placement = ipdf_post_advance_fields
    l_use_precip_frac = .true.
    l_predict_upwp_vpwp = .true.
    l_min_wp2_from_corr_wx = .true.
    l_min_xp2_from_corr_wx = .true.
    l_C2_cloud_frac = .false.
    l_diffuse_rtm_and_thlm = .false.
    l_stability_correct_Kh_N2_zm = .false.
    l_calc_thlp2_rad = .true.
    l_upwind_xpyp_ta = .true.
    l_upwind_xm_ma = .true.
    l_uv_nudge = .false.
    l_rtm_nudge = .false.
    l_tke_aniso = .true.
    l_vert_avg_closure  = .false.
    l_trapezoidal_rule_zt = .false.
    l_trapezoidal_rule_zm = .false.
    l_call_pdf_closure_twice = .false.
    l_standard_term_ta = .false.
    l_partial_upwind_wp3 = .false.
    l_godunov_upwind_wpxp_ta = .false.
    l_godunov_upwind_xpyp_ta = .false.
    l_use_cloud_cover = .false.
    l_diagnose_correlations = .false.
    l_calc_w_corr = .false.
    l_const_Nc_in_cloud = .false.
    l_fix_w_chi_eta_correlations = .true.
    l_stability_correct_tau_zm = .false.
    l_damp_wp2_using_em = .true.
    l_do_expldiff_rtm_thlm = .false.
    l_Lscale_plume_centered = .false.
    l_diag_Lscale_from_tau = .true.
    l_use_C7_Richardson = .true.
    l_use_C11_Richardson = .false.
    l_use_shear_Richardson = .false.
    l_brunt_vaisala_freq_moist = .false.
    l_use_thvm_in_bv_freq = .false.
    l_rcm_supersat_adj = .true.
    l_damp_wp3_Skw_squared = .true.
#ifdef GFDL
    l_prescribed_avg_deltaz = .true.
#else
    l_prescribed_avg_deltaz = .false.
#endif
    l_lmm_stepping = .false.
    l_e3sm_config = .false.
    l_vary_convect_depth = .false.
    l_use_tke_in_wp3_pr_turb_term = .false.
    l_use_tke_in_wp2_wp3_K_dfsn = .false.
    l_smooth_Heaviside_tau_wpxp = .false.
    l_enable_relaxed_clipping = .false.
    l_linearize_pbl_winds = .false.

    return
  end subroutine set_default_clubb_config_flags

!===============================================================================
  subroutine initialize_clubb_config_flags_type( iiPDF_type,          &
                                                 ipdf_call_placement, &
                                                 l_use_precip_frac, &
                                                 l_predict_upwp_vpwp, &
                                                 l_min_wp2_from_corr_wx, &
                                                 l_min_xp2_from_corr_wx, &
                                                 l_C2_cloud_frac, &
                                                 l_diffuse_rtm_and_thlm, &
                                                 l_stability_correct_Kh_N2_zm, &
                                                 l_calc_thlp2_rad, &
                                                 l_upwind_xpyp_ta, &
                                                 l_upwind_xm_ma, &
                                                 l_uv_nudge, &
                                                 l_rtm_nudge, &
                                                 l_tke_aniso, &
                                                 l_vert_avg_closure, &
                                                 l_trapezoidal_rule_zt, &
                                                 l_trapezoidal_rule_zm, &
                                                 l_call_pdf_closure_twice, &
                                                 l_standard_term_ta, &
                                                 l_partial_upwind_wp3, &
                                                 l_godunov_upwind_wpxp_ta, &
                                                 l_godunov_upwind_xpyp_ta, &
                                                 l_use_cloud_cover, &
                                                 l_diagnose_correlations, &
                                                 l_calc_w_corr, &
                                                 l_const_Nc_in_cloud, &
                                                 l_fix_w_chi_eta_correlations, &
                                                 l_stability_correct_tau_zm, &
                                                 l_damp_wp2_using_em, &
                                                 l_do_expldiff_rtm_thlm, &
                                                 l_Lscale_plume_centered, &
                                                 l_diag_Lscale_from_tau, &
                                                 l_use_C7_Richardson, &
                                                 l_use_C11_Richardson, &
                                                 l_use_shear_Richardson, &
                                                 l_brunt_vaisala_freq_moist, &
                                                 l_use_thvm_in_bv_freq, &
                                                 l_rcm_supersat_adj, &
                                                 l_damp_wp3_Skw_squared, &
                                                 l_prescribed_avg_deltaz, &
                                                 l_lmm_stepping, &
                                                 l_e3sm_config, &
                                                 l_vary_convect_depth, &
                                                 l_use_tke_in_wp3_pr_turb_term, &
                                                 l_use_tke_in_wp2_wp3_K_dfsn, &
                                                 l_smooth_Heaviside_tau_wpxp, &
                                                 l_enable_relaxed_clipping, &
                                                 l_linearize_pbl_winds, &
                                                 clubb_config_flags )

! Description:
!   Initialize the clubb_config_flags_type.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Peskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp in
                                      ! xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds           ! Code to linearize PBL winds

    ! Output variables
    type(clubb_config_flags_type), intent(out) :: &
      clubb_config_flags            ! Derived type holding all configurable CLUBB flags

!-----------------------------------------------------------------------
    ! Begin code

    clubb_config_flags%iiPDF_type = iiPDF_type
    clubb_config_flags%ipdf_call_placement = ipdf_call_placement
    clubb_config_flags%l_use_precip_frac = l_use_precip_frac
    clubb_config_flags%l_predict_upwp_vpwp = l_predict_upwp_vpwp
    clubb_config_flags%l_min_wp2_from_corr_wx = l_min_wp2_from_corr_wx
    clubb_config_flags%l_min_xp2_from_corr_wx = l_min_xp2_from_corr_wx
    clubb_config_flags%l_C2_cloud_frac = l_C2_cloud_frac
    clubb_config_flags%l_diffuse_rtm_and_thlm = l_diffuse_rtm_and_thlm
    clubb_config_flags%l_stability_correct_Kh_N2_zm = l_stability_correct_Kh_N2_zm
    clubb_config_flags%l_calc_thlp2_rad = l_calc_thlp2_rad
    clubb_config_flags%l_upwind_xpyp_ta = l_upwind_xpyp_ta
    clubb_config_flags%l_upwind_xm_ma = l_upwind_xm_ma
    clubb_config_flags%l_uv_nudge = l_uv_nudge
    clubb_config_flags%l_rtm_nudge = l_rtm_nudge
    clubb_config_flags%l_tke_aniso = l_tke_aniso
    clubb_config_flags%l_vert_avg_closure  = l_vert_avg_closure
    clubb_config_flags%l_trapezoidal_rule_zt = l_trapezoidal_rule_zt
    clubb_config_flags%l_trapezoidal_rule_zm = l_trapezoidal_rule_zm
    clubb_config_flags%l_call_pdf_closure_twice = l_call_pdf_closure_twice
    clubb_config_flags%l_standard_term_ta = l_standard_term_ta
    clubb_config_flags%l_partial_upwind_wp3 = l_partial_upwind_wp3
    clubb_config_flags%l_godunov_upwind_wpxp_ta = l_godunov_upwind_wpxp_ta
    clubb_config_flags%l_godunov_upwind_xpyp_ta = l_godunov_upwind_xpyp_ta
    clubb_config_flags%l_use_cloud_cover = l_use_cloud_cover
    clubb_config_flags%l_diagnose_correlations = l_diagnose_correlations
    clubb_config_flags%l_calc_w_corr = l_calc_w_corr
    clubb_config_flags%l_const_Nc_in_cloud = l_const_Nc_in_cloud
    clubb_config_flags%l_fix_w_chi_eta_correlations = l_fix_w_chi_eta_correlations
    clubb_config_flags%l_stability_correct_tau_zm = l_stability_correct_tau_zm
    clubb_config_flags%l_damp_wp2_using_em = l_damp_wp2_using_em
    clubb_config_flags%l_do_expldiff_rtm_thlm = l_do_expldiff_rtm_thlm
    clubb_config_flags%l_Lscale_plume_centered = l_Lscale_plume_centered
    clubb_config_flags%l_diag_Lscale_from_tau = l_diag_Lscale_from_tau
    clubb_config_flags%l_use_C7_Richardson = l_use_C7_Richardson
    clubb_config_flags%l_use_C11_Richardson = l_use_C11_Richardson
    clubb_config_flags%l_use_shear_Richardson = l_use_shear_Richardson
    clubb_config_flags%l_brunt_vaisala_freq_moist = l_brunt_vaisala_freq_moist
    clubb_config_flags%l_use_thvm_in_bv_freq = l_use_thvm_in_bv_freq
    clubb_config_flags%l_rcm_supersat_adj = l_rcm_supersat_adj
    clubb_config_flags%l_damp_wp3_Skw_squared = l_damp_wp3_Skw_squared
    clubb_config_flags%l_prescribed_avg_deltaz = l_prescribed_avg_deltaz
    clubb_config_flags%l_lmm_stepping = l_lmm_stepping
    clubb_config_flags%l_e3sm_config = l_e3sm_config
    clubb_config_flags%l_vary_convect_depth = l_vary_convect_depth
    clubb_config_flags%l_use_tke_in_wp3_pr_turb_term = l_use_tke_in_wp3_pr_turb_term
    clubb_config_flags%l_use_tke_in_wp2_wp3_K_dfsn = l_use_tke_in_wp2_wp3_K_dfsn
    clubb_config_flags%l_smooth_Heaviside_tau_wpxp = l_smooth_Heaviside_tau_wpxp
    clubb_config_flags%l_enable_relaxed_clipping = l_enable_relaxed_clipping
    clubb_config_flags%l_linearize_pbl_winds = l_linearize_pbl_winds

    return
  end subroutine initialize_clubb_config_flags_type

!===============================================================================
  subroutine print_clubb_config_flags( iunit, clubb_config_flags )

! Description:
!   Prints the clubb_config_flags.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iunit ! The file to write to

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

!-----------------------------------------------------------------------
    ! Begin code

    write(iunit,*) "iiPDF_type = ", clubb_config_flags%iiPDF_type
    write(iunit,*) "ipdf_call_placement = ", clubb_config_flags%ipdf_call_placement
    write(iunit,*) "l_use_precip_frac = ", clubb_config_flags%l_use_precip_frac
    write(iunit,*) "l_predict_upwp_vpwp = ", clubb_config_flags%l_predict_upwp_vpwp
    write(iunit,*) "l_min_wp2_from_corr_wx = ", clubb_config_flags%l_min_wp2_from_corr_wx
    write(iunit,*) "l_min_xp2_from_corr_wx = ", clubb_config_flags%l_min_xp2_from_corr_wx
    write(iunit,*) "l_C2_cloud_frac = ", clubb_config_flags%l_C2_cloud_frac
    write(iunit,*) "l_diffuse_rtm_and_thlm = ", clubb_config_flags%l_diffuse_rtm_and_thlm
    write(iunit,*) "l_stability_correct_Kh_N2_zm = ", &
                   clubb_config_flags%l_stability_correct_Kh_N2_zm
    write(iunit,*) "l_calc_thlp2_rad = ", clubb_config_flags%l_calc_thlp2_rad
    write(iunit,*) "l_upwind_xpyp_ta = ", clubb_config_flags%l_upwind_xpyp_ta
    write(iunit,*) "l_upwind_xm_ma = ", clubb_config_flags%l_upwind_xm_ma
    write(iunit,*) "l_uv_nudge = ", clubb_config_flags%l_uv_nudge
    write(iunit,*) "l_rtm_nudge = ", clubb_config_flags%l_rtm_nudge
    write(iunit,*) "l_tke_aniso = ", clubb_config_flags%l_tke_aniso
    write(iunit,*) "l_vert_avg_closure = ", clubb_config_flags%l_vert_avg_closure
    write(iunit,*) "l_trapezoidal_rule_zt = ", clubb_config_flags%l_trapezoidal_rule_zt
    write(iunit,*) "l_trapezoidal_rule_zm = ", clubb_config_flags%l_trapezoidal_rule_zm
    write(iunit,*) "l_call_pdf_closure_twice = ", clubb_config_flags%l_call_pdf_closure_twice
    write(iunit,*) "l_standard_term_ta = ", clubb_config_flags%l_standard_term_ta
    write(iunit,*) "l_partial_upwind_wp3 = ", clubb_config_flags%l_partial_upwind_wp3
    write(iunit,*) "l_godunov_upwind_wpxp_ta = ", clubb_config_flags%l_godunov_upwind_wpxp_ta
    write(iunit,*) "l_godunov_upwind_xpyp_ta = ", clubb_config_flags%l_godunov_upwind_xpyp_ta
    write(iunit,*) "l_use_cloud_cover = ", clubb_config_flags%l_use_cloud_cover
    write(iunit,*) "l_diagnose_correlations = ", clubb_config_flags%l_diagnose_correlations
    write(iunit,*) "l_calc_w_corr = ", clubb_config_flags%l_calc_w_corr
    write(iunit,*) "l_const_Nc_in_cloud = ", clubb_config_flags%l_const_Nc_in_cloud
    write(iunit,*) "l_fix_w_chi_eta_correlations = ", clubb_config_flags%l_fix_w_chi_eta_correlations
    write(iunit,*) "l_stability_correct_tau_zm = ", clubb_config_flags%l_stability_correct_tau_zm
    write(iunit,*) "l_damp_wp2_using_em = ", clubb_config_flags%l_damp_wp2_using_em
    write(iunit,*) "l_do_expldiff_rtm_thlm = ", clubb_config_flags%l_do_expldiff_rtm_thlm
    write(iunit,*) "l_Lscale_plume_centered = ", clubb_config_flags%l_Lscale_plume_centered
    write(iunit,*) "l_diag_Lscale_from_tau = ", clubb_config_flags%l_diag_Lscale_from_tau
    write(iunit,*) "l_use_C7_Richardson = ", clubb_config_flags%l_use_C7_Richardson
    write(iunit,*) "l_use_C11_Richardson = ", clubb_config_flags%l_use_C11_Richardson
    write(iunit,*) "l_use_shear_Richardson = ", clubb_config_flags%l_use_shear_Richardson
    write(iunit,*) "l_brunt_vaisala_freq_moist = ", clubb_config_flags%l_brunt_vaisala_freq_moist
    write(iunit,*) "l_use_thvm_in_bv_freq = ", clubb_config_flags%l_use_thvm_in_bv_freq
    write(iunit,*) "l_rcm_supersat_adj = ", clubb_config_flags%l_rcm_supersat_adj
    write(iunit,*) "l_damp_wp3_Skw_squared = ", clubb_config_flags%l_damp_wp3_Skw_squared
    write(iunit,*) "l_prescribed_avg_deltaz = ", clubb_config_flags%l_prescribed_avg_deltaz
    write(iunit,*) "l_lmm_stepping = ", clubb_config_flags%l_lmm_stepping
    write(iunit,*) "l_e3sm_config = ", clubb_config_flags%l_e3sm_config
    write(iunit,*) "l_vary_convect_depth", clubb_config_flags%l_vary_convect_depth
    write(iunit,*) "l_use_tke_in_wp3_pr_turb_term = ", clubb_config_flags%l_use_tke_in_wp3_pr_turb_term
    write(iunit,*) "l_use_tke_in_wp2_wp3_K_dfsn = ", clubb_config_flags%l_use_tke_in_wp2_wp3_K_dfsn
    write(iunit,*) "l_smooth_Heaviside_tau_wpxp = ", clubb_config_flags%l_smooth_Heaviside_tau_wpxp
    write(iunit,*) "l_enable_relaxed_clipping = ", clubb_config_flags%l_enable_relaxed_clipping
    write(iunit,*) "l_linearize_pbl_winds = ", clubb_config_flags%l_linearize_pbl_winds

    return
  end subroutine print_clubb_config_flags

end module model_flags
