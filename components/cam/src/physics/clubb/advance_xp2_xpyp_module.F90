!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module advance_xp2_xpyp_module

  ! Description:
  ! Contains the subroutine advance_xp2_xpyp and ancillary functions.
  !-----------------------------------------------------------------------

  implicit none

  public :: advance_xp2_xpyp, &
            update_xp2_mc

  private :: xp2_xpyp_lhs, & 
             xp2_xpyp_solve, & 
             xp2_xpyp_uv_rhs, & 
             xp2_xpyp_rhs, & 
             xp2_xpyp_implicit_stats, & 
             term_tp, & 
             term_dp1_lhs, & 
             term_dp1_rhs, & 
             term_pr1, & 
             term_pr2

  private    ! Set default scope

  ! Private named constants to avoid string comparisons
  integer, parameter, private :: &
    xp2_xpyp_rtp2 = 1, &     ! Named constant for rtp2 solves
    xp2_xpyp_thlp2 = 2, &    ! Named constant for thlp2 solves
    xp2_xpyp_rtpthlp = 3, &  ! Named constant for rtpthlp solves
    xp2_xpyp_up2_vp2 = 4, &  ! Named constant for up2_vp2 solves
    xp2_xpyp_up2 = 5, &      ! Named constant for up2 solves
    xp2_xpyp_vp2 = 6, &      ! Named constant for vp2 solves
    xp2_xpyp_scalars = 7, &  ! Named constant for scalar solves
    xp2_xpyp_sclrp2 = 8, &   ! Named constant for sclrp2 solves
    xp2_xpyp_sclrprtp = 9, & ! Named constant for sclrprtp solves
    xp2_xpyp_sclrpthlp = 10  ! Named constant for sclrpthlp solves

  contains

  !=============================================================================
  subroutine advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp, thlm,        & ! In
                               wpthlp, wpthvp, um, vm, wp2, wp2_zt,    & ! In
                               wp3, upwp, vpwp, sigma_sqd_w, Skw_zm,   & ! In
                               wprtp2, wpthlp2, wprtpthlp,             & ! In
                               Kh_zt, rtp2_forcing, thlp2_forcing,     & ! In
                               rtpthlp_forcing, rho_ds_zm, rho_ds_zt,  & ! In
                               invrs_rho_ds_zm, thv_ds_zm, cloud_frac, & ! In
                               Lscale, wp3_on_wp2, wp3_on_wp2_zt,      & ! In
                               pdf_implicit_coefs_terms,               & ! In
                               l_iter, dt,                             & ! In
                               sclrm, wpsclrp,                         & ! In
                               wpsclrp2, wpsclrprtp, wpsclrpthlp,      & ! In
                               wp2_splat,                              & ! In
                               rtp2, thlp2, rtpthlp, up2, vp2,         & ! Inout
                               sclrp2, sclrprtp, sclrpthlp)              ! Inout

    ! Description:
    ! Prognose scalar variances, scalar covariances, and horizontal turbulence components.

    ! References:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xpyp_eqns
    !
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_vp2_eqns
    !  
    !   Eqn. 13, 14, 15  on p. 3545 of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !     Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.

    ! See also:
    !   ``Equations for CLUBB'', Section 4:
    !   /Steady-state solution for the variances/
    !-----------------------------------------------------------------------

    use constants_clubb, only: & 
        w_tol_sqd,  & ! Constant(s)
        rt_tol, & 
        thl_tol, &
        max_mag_correlation_flux, &
        cloud_frac_min, &
        fstderr, &
        one, &
        two_thirds, &
        one_half, &
        one_third, &
        zero, &
        zero_threshold

    use model_flags, only: & 
        l_hole_fill, &    ! logical constants
        l_single_C2_Skw, &
        l_explicit_turbulent_adv_xpyp, &
        l_upwind_xpyp_ta, &
        l_min_xp2_from_corr_wx, &
        l_C2_cloud_frac

    use parameters_tunable, only: &
        C2rt,     & ! Variable(s)
        C2thl,    &
        C2rtthl,  &
        c_K2,     &
        nu2_vert_res_dep, &
        c_K9,     &
        nu9_vert_res_dep, &
        beta,     &
        C4,       &
        C14,      &
        C5,       &
        C2,       &
        C2b,      &
        C2c

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        sclr_tol

    use grid_class, only: & 
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm

    use pdf_closure_module, only: &
        iiPDF_ADG1, & ! Variable(s)
        iiPDF_new,  &
        iiPDF_type

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use turbulent_adv_pdf, only: &
        sgn_turbulent_velocity    ! Procedure(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use clip_explicit, only: & 
        clip_covar,  & ! Procedure(s)
        clip_variance, &
        clip_variance_level, &
        clip_sclrp2, &
        clip_sclrprtp, &
        clip_sclrpthlp

    use sponge_layer_damping, only: &
        up2_vp2_sponge_damp_settings, & ! Variable(s)
        up2_vp2_sponge_damp_profile,  &
        sponge_damp_xp2                 ! Procedure(s)
      
    use stats_type_utilities, only: &
        stat_begin_update, & ! Procedure(s)
        stat_end_update,   &
        stat_modify,       &
        stat_update_var

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constants

    use stats_variables, only: &
        stats_zm,                 & ! Variable(s)
        stats_zt,                 &
        icoef_wprtp2_implicit,    &
        iterm_wprtp2_explicit,    &
        icoef_wpthlp2_implicit,   &
        iterm_wpthlp2_explicit,   &
        icoef_wprtpthlp_implicit, &
        iterm_wprtpthlp_explicit, &
        irtp2_cl,                 &
        iup2_sdmp,                &
        ivp2_sdmp,                &
        iup2_cl,                  &
        ivp2_cl,                  &
        l_stats_samp

    use array_index, only: &
        iisclr_rt, &
        iisclr_thl

    implicit none

    ! Intrinsic functions
    intrinsic :: &
      exp, sqrt, min

    ! Constant parameters
    logical, parameter :: &
      l_clip_large_rtp2 = .true. ! Clip rtp2 to be < rtm^2 * coef

    real( kind = core_rknd ), parameter :: &
      rtp2_clip_coef = one_half ! Coefficient appled the clipping threshold on rtp2 [-]

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      tau_zm,          & ! Time-scale tau on momentum levels     [s]
      wm_zm,           & ! w-wind component on momentum levels   [m/s]
      rtm,             & ! Total water mixing ratio (t-levs)     [kg/kg]
      wprtp,           & ! <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
      thlm,            & ! Liquid potential temp. (t-levs)       [K]
      wpthlp,          & ! <w'th_l'> (momentum levels)           [(m K)/s]
      wpthvp,          & ! <w'th_v'> (momentum levels)           [(m K)/s]
      um,              & ! u wind (thermodynamic levels)         [m/s]
      vm,              & ! v wind (thermodynamic levels)         [m/s]
      wp2,             & ! <w'^2> (momentum levels)              [m^2/s^2]
      wp2_zt,          & ! <w'^2> interpolated to thermo. levels [m^2/s^2]
      wp3,             & ! <w'^3> (thermodynamic levels)         [m^3/s^3]
      upwp,            & ! <u'w'> (momentum levels)              [m^2/s^2]
      vpwp,            & ! <v'w'> (momentum levels)              [m^2/s^2]
      sigma_sqd_w,     & ! sigma_sqd_w (momentum levels)         [-]
      Skw_zm,          & ! Skewness of w on momentum levels      [-]
      wprtp2,          & ! <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
      wpthlp2,         & ! <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
      wprtpthlp,       & ! <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
      Kh_zt,           & ! Eddy diffusivity on thermo. levels    [m^2/s]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
      rho_ds_zm,       & ! Dry, static density on momentum levs. [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ mom. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on mom. levs. [K]
      cloud_frac,      & ! Cloud fraction (thermodynamic levels) [-]
      Lscale,          & ! Mixing length                         [m]
      wp3_on_wp2,      & ! Smoothed version of <w'^3>/<w'^2> zm  [m/s]
      wp3_on_wp2_zt      ! Smoothed version of <w'^3>/<w'^2> zt  [m/s]

    type(implicit_coefs_terms), dimension(gr%nz), intent(in) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    logical, intent(in) :: l_iter ! Whether variances are prognostic

    real( kind = core_rknd ), intent(in) :: &
      dt             ! Model timestep                                [s]

    ! Passive scalar input
    real( kind = core_rknd ), intent(in), dimension(gr%nz, sclr_dim) ::  & 
      sclrm,       & ! Mean value; pass. scalar (t-levs.) [{sclr units}]
      wpsclrp,     & ! <w'sclr'> (momentum levels)        [m/s{sclr units}]
      wpsclrp2,    & ! <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
      wpsclrprtp,  & ! <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
      wpsclrpthlp    ! <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wp2_splat    ! Gustiness tendency for wp2 equation

    ! Input/Output variables
    ! An attribute of (inout) is also needed to import the value of the variances
    ! at the surface.  Brian.  12/18/05.
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      rtp2,    & ! <r_t'^2>                      [(kg/kg)^2]
      thlp2,   & ! <th_l'^2>                     [K^2]
      rtpthlp, & ! <r_t'th_l'>                   [(kg K)/kg]
      up2,     & ! <u'^2>                        [m^2/s^2]
      vp2        ! <v'^2>                        [m^2/s^2]

    ! Passive scalar output
    real( kind = core_rknd ), intent(inout), dimension(gr%nz, sclr_dim) ::  & 
      sclrp2, sclrprtp, sclrpthlp

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: & 
      C2sclr_1d, C2rt_1d, C2thl_1d, C2rtthl_1d, &
      C4_C14_1d     ! Parameters C4 and C14 combined for simplicity

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      a1 ! a_1 (momentum levels); See eqn. 24 in `Equations for CLUBB' [-]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      upwp_zt,    & ! <u'w'> interpolated to thermodynamic levels    [m^2/s^2]
      vpwp_zt,    & ! <v'w'> interpolated to thermodynamic levels    [m^2/s^2]
      wpsclrp_zt    ! <w'sclr'> interp. to thermo. levels  [m/s {sclrm units}]

    real( kind = core_rknd ) :: & 
      threshold     ! Minimum value for variances                   [units vary]

    real( kind = core_rknd ), dimension(3,gr%nz) ::  & 
      lhs ! Tridiagonal matrix

    real( kind = core_rknd ), dimension(gr%nz,1) :: & 
      rhs ! RHS vector of tridiagonal matrix

    real( kind = core_rknd ), dimension(gr%nz,2) :: & 
      uv_rhs,    &! RHS vectors of tridiagonal system for up2/vp2
      uv_solution ! Solution to the tridiagonal system for up2/vp2

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim*3) ::  & 
      sclr_rhs,   & ! RHS vectors of tridiagonal system for the passive scalars
      sclr_solution ! Solution to tridiagonal system for the passive scalars

    ! Variables for turbulent advection of predictive variances and covariances.

    ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtp2_implicit, & ! Coefficient that is multiplied by <rt'^2>  [m/s]
      term_wprtp2_explicit    ! Term that is on the RHS           [m/s(kg/kg)^2]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtp2_implicit_zm, & ! coef_wprtp2_implicit interp. to m-levs  [m/s]
      term_wprtp2_explicit_zm    ! term_wprtp2_expl interp m-levs [m/s(kg/kg)^2]

    ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wpthlp2_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpthlp2_implicit, & ! Coef. that is multiplied by <thl'^2>      [m/s]
      term_wpthlp2_explicit    ! Term that is on the RHS               [m/s K^2]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpthlp2_implicit_zm, & ! coef_wpthlp2_implicit interp. m-levs   [m/s]
      term_wpthlp2_explicit_zm    ! term_wpthlp2_expl interp to m-levs [m/s K^2]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtpthlp_implicit_zm, & ! coef_wprtpthlp_impl interp. m-levs   [m/s]
      term_wprtpthlp_explicit_zm    ! term_wprtpthlp_expl intrp zm [m/s(kg/kg)K]

    ! CLUBB does not produce a PDF for horizontal wind components u and v.
    ! However, turbulent advection of the variances of the horizontal wind
    ! components, <u'^2> and <v'^2>, is still handled by equations of the form:
    ! <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit; and
    ! <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit.
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpup2_wpvp2_implicit, & ! Coef. that is mult. by <u'^2>/<v'^2>  [m/s]
      term_wpup2_explicit,       & ! Term that is on the RHS (<u'^2>)  [m^3/s^3]
      term_wpvp2_explicit          ! Term that is on the RHS (<v'^2>)  [m^3/s^3]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpup2_wpvp2_implicit_zm, & ! coef_wpup2_wpvp2_impl intrp m-levs [m/s]
      term_wpup2_explicit_zm,       & ! term_wpup2_expl interp. m-levs [m^3/s^3]
      term_wpvp2_explicit_zm          ! term_wpvp2_expl interp. m-levs [m^3/s^3]

    ! Sign of turbulent velocity (used for "upwind" turbulent advection)
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      sgn_t_vel_rtp2,    & ! Sign of the turbulent velocity for <rt'^2>      [-]
      sgn_t_vel_thlp2,   & ! Sign of the turbulent velocity for <thl'^2>     [-]
      sgn_t_vel_rtpthlp, & ! Sign of the turbulent velocity for <rt'thl'>    [-]
      sgn_t_vel_up2_vp2    ! Sign of the turbulent vel. for <u'^2>/<v'^2>    [-]

    ! <w'sclr'x'> = coef_wpsclrpxp_implicit*<sclr'x'> + term_wpsclrpxp_explicit;
    ! where x is sclr, rt, or thl.
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpsclrpxp_implicit,   & ! Coef. that is multiplied by <sclr'x'> [m/s]
      term_wpsclrp2_explicit,    & ! Term that is on the RHS        [units vary]
      term_wpsclrprtp_explicit,  & ! Term that is on the RHS        [units vary]
      term_wpsclrpthlp_explicit    ! Term that is on the RHS        [units vary]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wpsclrpxp_implicit_zm,   & ! coef_wpsclrpxp_impl interp. m-levs [m/s]
      term_wpsclrp2_explicit_zm,    & ! term_wpsclrp2_expl interp zm   [un vary]
      term_wpsclrprtp_explicit_zm,  & ! term_wpsclrprtp_expl interp zm [un vary]
      term_wpsclrpthlp_explicit_zm    ! term_wpsclrpthlp_expl intrp zm [un vary]

    ! Sign of turbulent velocity (used for "upwind" turbulent advection)
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      sgn_t_vel_sclrpxp,   & ! Sign of the turbulent velocity for <sclr'x'>  [-]
      sgn_t_vel_sclrp2,    & ! Sign of the turbulent velocity for <sclr'^2>  [-]
      sgn_t_vel_sclrprtp,  & ! Sign of the turbulent velocity for <sclr'rt'> [-]
      sgn_t_vel_sclrpthlp    ! Sign of the turbulent vel. for <sclr'thl'>    [-]

    ! Eddy Diffusion for Variances and Covariances.
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      Kw2, & ! For rtp2, thlp2, rtpthlp, and passive scalars  [m^2/s]
      Kw9    ! For up2 and vp2                                [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      a1_zt,     & ! a_1 interpolated to thermodynamic levels      [-]
      wprtp_zt,  & ! w'r_t' interpolated to thermodynamic levels   [(kg/kg) m/s]
      wpthlp_zt    ! w'th_l' interpolated to thermodyamnic levels  [K m/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rtpthlp_chnge    ! Net change in r_t'th_l' due to clipping [(kg/kg) K]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      sclrprtp_chnge,  & ! Net change in sclr'r_t' due to clipping  [units vary]
      sclrpthlp_chnge    ! Net change in sclr'th_l' due to clipping [units vary]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      sclrp2_forcing,    & ! <sclr'^2> forcing (momentum levels)    [units vary]
      sclrprtp_forcing,  & ! <sclr'r_t'> forcing (momentum levels)  [units vary]
      sclrpthlp_forcing    ! <sclr'th_l'> forcing (momentum levels) [units vary]

    logical :: l_scalar_calc, l_first_clip_ts, l_last_clip_ts

    ! Minimum value of cloud fraction to multiply C2 by to calculate the
    ! adjusted value of C2.
    real( kind = core_rknd ), parameter ::  & 
      min_cloud_frac_mult = 0.10_core_rknd

    ! Loop indices
    integer :: i, k

    !---------------------------- Begin Code ----------------------------------

    if ( clubb_at_least_debug_level( 0 ) ) then
      ! Assertion check for C5
      if ( C5 > one .or. C5 < zero ) then
        write(fstderr,*) "The C5 variable is outside the valid range"
        err_code = clubb_fatal_error
        return
      end if
    end if

    if ( l_single_C2_Skw ) then

       ! Use a single value of C2 for all equations.
       C2rt_1d = C2b + ( C2 - C2b ) * exp( -one_half * ( Skw_zm / C2c )**2 )

       if ( l_C2_cloud_frac ) then
          do k = 1, gr%nz, 1
             if ( cloud_frac(k) >= cloud_frac_min ) then
                C2rt_1d(k) &
                = C2rt_1d(k) * max( min_cloud_frac_mult, cloud_frac(k) )
             endif ! cloud_frac(k) >= cloud_frac_min
          enddo ! k = 1, gr%nz, 1
       endif ! l_C2_cloud_frac

       C2thl_1d   = C2rt_1d
       C2rtthl_1d = C2rt_1d

       C2sclr_1d  = C2rt_1d

    else ! .not. l_single_C2_Skw

       ! Use 3 different values of C2 for rtp2, thlp2, rtpthlp.
       if ( l_C2_cloud_frac ) then

          do k = 1, gr%nz, 1
             if ( cloud_frac(k) >= cloud_frac_min ) then
                C2rt_1d(k) = C2rt * max( min_cloud_frac_mult, cloud_frac(k) )
                C2thl_1d(k) = C2thl * max( min_cloud_frac_mult, cloud_frac(k) )
                C2rtthl_1d(k) = C2rtthl &
                                * max( min_cloud_frac_mult, cloud_frac(k) )
             else ! cloud_frac(k) < cloud_frac_min
                C2rt_1d(k)    = C2rt
                C2thl_1d(k)   = C2thl
                C2rtthl_1d(k) = C2rtthl
             endif ! cloud_frac(k) >= cloud_frac_min
          enddo ! k = 1, gr%nz, 1

       else

          C2rt_1d    = C2rt
          C2thl_1d   = C2thl
          C2rtthl_1d = C2rtthl

       endif ! l_C2_cloud_frac

       C2sclr_1d = C2rt  ! Use rt value for now

    endif ! l_single_C2_Skw

    ! Combine C4 and C14 for simplicity
    C4_C14_1d(1:gr%nz) = ( two_thirds * C4 ) + ( one_third * C14 )

    ! Are we solving for passive scalars as well?
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if

    ! Define a_1 (located on momentum levels).
    ! It is a variable that is a function of sigma_sqd_w (where sigma_sqd_w is
    ! located on the momentum levels).  This will be used for the turbulent
    ! advection (ta) terms for the ADG1 PDF.  This will also be used for the
    ! turbulent advection of <u'^2> and <v'^2>, regardless of which PDF type or
    ! turbulent advection option is used.
    a1(1:gr%nz) = one / ( one - sigma_sqd_w(1:gr%nz) )

    ! Interpolate a_1 from the momentum levels to the thermodynamic levels.
    a1_zt = max( zm2zt( a1 ), zero_threshold ) ! Positive definite quantity

    ! Set up the implicit coefficients and explicit terms for turbulent
    ! advection of <rt'^2>, <thl'^2>, and <rt'thl'>.
    if ( l_explicit_turbulent_adv_xpyp ) then

       ! The turbulent advection of <x'y'> is handled explicitly.

       ! The <rt'^2> turbulent advection term is entirely explicit, as
       ! term_wprtp2_explicit is equal to <w'rt'^2> as calculated using PDF
       ! parameters, which is general for any PDF type.  The value of
       ! <w'rt'^2> is calculated on thermodynamic levels.  The value of
       ! coef_wprtp2_implicit is always 0.
       coef_wprtp2_implicit = zero
       term_wprtp2_explicit = wprtp2

       if ( l_upwind_xpyp_ta ) then

          ! Interpolate term_wprtp2_explicit to momentum levels as
          ! term_wprtp2_explicit_zm.  The value of coef_wprtp2_implicit_zm is
          ! always 0.
          coef_wprtp2_implicit_zm = zero
          term_wprtp2_explicit_zm = zt2zm( term_wprtp2_explicit )

          ! Calculate the sign of the turbulent velocity for <rt'^2>.
          sgn_t_vel_rtp2 &
          = sgn_turbulent_velocity( term_wprtp2_explicit_zm, rtp2 )

       endif ! l_upwind_xpyp_ta

       ! The <thl'^2> turbulent advection term is entirely explicit, as
       ! term_wpthlp2_explicit is equal to <w'thl'^2> as calculated using PDF
       ! parameters, which is general for any PDF type.  The value of
       ! <w'thl'^2> is calculated on thermodynamic levels.  The value of
       ! coef_wpthlp2_implicit is always 0.
       coef_wpthlp2_implicit = zero
       term_wpthlp2_explicit = wpthlp2

       if ( l_upwind_xpyp_ta ) then

          ! Interpolate term_wpthlp2_explicit to momentum levels as
          ! term_wpthlp2_explicit_zm.  The value of coef_wpthlp2_implicit_zm is
          ! always 0.
          coef_wpthlp2_implicit_zm = zero
          term_wpthlp2_explicit_zm = zt2zm( term_wpthlp2_explicit )

          ! Calculate the sign of the turbulent velocity for <thl'^2>.
          sgn_t_vel_thlp2 &
          = sgn_turbulent_velocity( term_wpthlp2_explicit_zm, thlp2 )

       endif ! l_upwind_xpyp_ta

       ! The <rt'thl'> turbulent advection term is entirely explicit, as
       ! term_wprtpthlp_explicit is equal to <w'rt'thl'> as calculated using PDF
       ! parameters, which is general for any PDF type.  The value of
       ! <w'rt'thl'> is calculated on thermodynamic levels.  The value of
       ! coef_wprtpthlp_implicit is always 0.
       coef_wprtpthlp_implicit = zero
       term_wprtpthlp_explicit = wprtpthlp

       if ( l_upwind_xpyp_ta ) then

          ! Interpolate term_wprtpthlp_explicit to momentum levels as
          ! term_wprtpthlp_explicit_zm.  The value of coef_wprtpthlp_implicit_zm
          ! is always 0.
          coef_wprtpthlp_implicit_zm = zero
          term_wprtpthlp_explicit_zm = zt2zm( term_wprtpthlp_explicit )

          ! Calculate the sign of the turbulent velocity for <rt'thl'>.
          sgn_t_vel_rtpthlp &
          = sgn_turbulent_velocity( term_wprtpthlp_explicit_zm, rtpthlp )

       endif ! l_upwind_xpyp_ta

    else ! .not. l_explicit_turbulent_adv_xpyp

       ! The turbulent advection of <x'y'> is handled implicitly or
       ! semi-implicitly.

       if ( iiPDF_type == iiPDF_ADG1 ) then

          ! The ADG1 PDF is used.

          ! Calculate the implicit coefficients and explicit terms on
          ! thermodynamic grid levels.

          ! Interpolate <w'rt'> and <w'thl'> from the momentum levels to the
          ! thermodynamic levels.  These will be used for the turbulent
          ! advection terms for the ADG1 PDF.
          wprtp_zt  = zm2zt( wprtp )
          wpthlp_zt = zm2zt( wpthlp )

          ! Implicit coefficient on <rt'^2> in <w'rt'^2> equation.
          coef_wprtp2_implicit = one_third * beta * a1_zt * wp3_on_wp2_zt

          ! Explicit (RHS) term in <w'rt'^2> equation.
          term_wprtp2_explicit &
          = ( one - one_third * beta ) * a1_zt**2 * wprtp_zt**2 &
            * wp3_on_wp2_zt / wp2_zt

          if ( l_upwind_xpyp_ta ) then

             ! Calculate coef_wprtp2_implicit and term_wprtp2_explicit on
             ! momentum levels as coef_wprtp2_implicit_zm and
             ! term_wprtp2_explicit_zm, respectively.
             coef_wprtp2_implicit_zm = one_third * beta * a1 * wp3_on_wp2

             term_wprtp2_explicit_zm &
             = ( one - one_third * beta ) * a1**2 * wprtp**2 * wp3_on_wp2 / wp2

             ! For ADG1, the sign of the turbulent velocity is the sign of
             ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity
             ! is set to wp3_on_wp2.
             sgn_t_vel_rtp2 = wp3_on_wp2

          endif ! l_upwind_xpyp_ta

          ! Implicit coefficient on <thl'^2> in <w'thl'^2> equation.
          ! For ADG1, this is the same as coef_wprtp2_implicit.
          coef_wpthlp2_implicit = coef_wprtp2_implicit

          ! Explicit (RHS) term in <w'thl'^2> equation.
          term_wpthlp2_explicit &
          = ( one - one_third * beta ) * a1_zt**2 * wpthlp_zt**2 &
            * wp3_on_wp2_zt / wp2_zt

          if ( l_upwind_xpyp_ta ) then

             ! Calculate coef_wpthlp2_implicit and term_wpthlp2_explicit on
             ! momentum levels as coef_wpthlp2_implicit_zm and
             ! term_wpthlp2_explicit_zm, respectively.
             coef_wpthlp2_implicit_zm = coef_wprtp2_implicit_zm

             term_wpthlp2_explicit_zm &
             = ( one - one_third * beta ) * a1**2 * wpthlp**2 * wp3_on_wp2 / wp2

             ! For ADG1, the sign of the turbulent velocity is the sign of
             ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity
             ! is set to wp3_on_wp2.
             sgn_t_vel_thlp2 = wp3_on_wp2

          endif ! l_upwind_xpyp_ta

          ! Implicit coefficient on <rt'thl'> in <w'rt'thl'> equation.
          ! For ADG1, this is the same as coef_wprtp2_implicit.
          coef_wprtpthlp_implicit = coef_wprtp2_implicit

          ! Explicit (RHS) term in <w'rt'thl'> equation.
          term_wprtpthlp_explicit &
          = ( one - one_third * beta ) * a1_zt**2 * wprtp_zt * wpthlp_zt &
            * wp3_on_wp2_zt / wp2_zt

          if ( l_upwind_xpyp_ta ) then

             ! Calculate coef_wprtpthlp_implicit and term_wprtpthlp_explicit
             ! on momentum levels as coef_wprtpthlp_implicit_zm and
             ! term_wprtpthlp_explicit_zm, respectively.
             coef_wprtpthlp_implicit_zm = coef_wprtp2_implicit_zm

             term_wprtpthlp_explicit_zm &
             = ( one - one_third * beta ) * a1**2 * wprtp * wpthlp &
               * wp3_on_wp2 / wp2

             ! For ADG1, the sign of the turbulent velocity is the sign of
             ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity
             ! is set to wp3_on_wp2.
             sgn_t_vel_rtpthlp = wp3_on_wp2

          endif ! l_upwind_xpyp_ta

       elseif ( iiPDF_type == iiPDF_new ) then

          ! The new PDF is used.

          ! Unpack the variables coef_wprtp2_implicit, coef_wpthlp2_implicit,
          ! coef_wprtpthlp_implicit, and term_wprtpthlp_explicit from
          ! pdf_implicit_coefs_terms.  The PDF parameters and the resulting
          ! implicit coefficients and explicit terms are calculated on
          ! thermodynamic levels.

          ! Implicit coefficient on <rt'^2> in <w'rt'^2> equation.
          coef_wprtp2_implicit = pdf_implicit_coefs_terms%coef_wprtp2_implicit

          ! The <rt'^2> turbulent advection term is entirely implicit, as
          ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>.  The value of
          ! term_wprtp2_explicit is always 0.
          term_wprtp2_explicit = zero

          if ( l_upwind_xpyp_ta ) then

             ! Interpolate coef_wprtp2_implicit to momentum levels as
             ! coef_wprtp2_implicit_zm.  The value of term_wprtp2_explicit_zm is
             ! always 0.
             coef_wprtp2_implicit_zm = zt2zm( coef_wprtp2_implicit )
             term_wprtp2_explicit_zm = zero

             ! Calculate the sign of the turbulent velocity for <rt'^2>.
             sgn_t_vel_rtp2 &
             = sgn_turbulent_velocity( coef_wprtp2_implicit_zm * rtp2, rtp2 )

          endif ! l_upwind_xpyp_ta

          ! Implicit coefficient on <thl'^2> in <w'thl'^2> equation.
          coef_wpthlp2_implicit = pdf_implicit_coefs_terms%coef_wpthlp2_implicit

          ! The <thl'^2> turbulent advection term is entirely implicit, as
          ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>.  The value of
          ! term_wpthlp2_explicit is always 0.
          term_wpthlp2_explicit = zero

          if ( l_upwind_xpyp_ta ) then

             ! Interpolate coef_wpthlp2_implicit to momentum levels as
             ! coef_wpthlp2_implicit_zm.  The value of term_wpthlp2_explicit_zm
             ! is always 0.
             coef_wpthlp2_implicit_zm = zt2zm( coef_wpthlp2_implicit )
             term_wpthlp2_explicit_zm = zero

             ! Calculate the sign of the turbulent velocity for <thlp'^2>.
             sgn_t_vel_thlp2 &
             = sgn_turbulent_velocity( coef_wpthlp2_implicit_zm * thlp2, thlp2 )

          endif ! l_upwind_xpyp_ta

          ! Implicit coefficient on <rt'thl'> in <w'rt'thl'> equation.
          coef_wprtpthlp_implicit &
          = pdf_implicit_coefs_terms%coef_wprtpthlp_implicit

          ! Explicit (RHS) term in <w'rt'thl'> equation.
          term_wprtpthlp_explicit &
          = pdf_implicit_coefs_terms%term_wprtpthlp_explicit

          if ( l_upwind_xpyp_ta ) then

             ! Interpolate coef_wprtpthlp_implicit and term_wprtpthlp_explicit
             ! to momentum levels as coef_wprtpthlp_implicit_zm and
             ! term_wprtpthlp_explicit_zm, respectively.
             coef_wprtpthlp_implicit_zm = zt2zm( coef_wprtpthlp_implicit )
             term_wprtpthlp_explicit_zm = zt2zm( term_wprtpthlp_explicit )

             ! Calculate the sign of the turbulent velocity for <rt'thl'>.
             sgn_t_vel_rtpthlp &
             = sgn_turbulent_velocity( coef_wprtpthlp_implicit_zm * rtpthlp &
                                       + term_wprtpthlp_explicit_zm, rtpthlp )

          endif ! l_upwind_xpyp_ta

       endif ! iiPDF_type

    endif ! l_explicit_turbulent_adv_xpyp

    if ( l_stats_samp ) then
       call stat_update_var( icoef_wprtp2_implicit, coef_wprtp2_implicit, &
                             stats_zt )
       call stat_update_var( iterm_wprtp2_explicit, term_wprtp2_explicit, &
                             stats_zt )
       call stat_update_var( icoef_wpthlp2_implicit, coef_wpthlp2_implicit, &
                             stats_zt )
       call stat_update_var( iterm_wpthlp2_explicit, term_wpthlp2_explicit, &
                             stats_zt )
       call stat_update_var( icoef_wprtpthlp_implicit, &
                             coef_wprtpthlp_implicit, stats_zt )
       call stat_update_var( iterm_wprtpthlp_explicit, &
                             term_wprtpthlp_explicit, stats_zt )
    endif ! l_stats_samp

    ! Set up the implicit coefficients and explicit terms for turbulent
    ! advection of <u'^2> and <v'^2>.

    ! CLUBB does not produce a PDF for horizontal wind components u and v.
    ! However, the code for the ADG1 PDF is still used to handle the turbulent
    ! advection for the variances of the horizontal wind components, <u'^2>
    ! and <v'^2>.  The ADG1 code is used regardless of which PDF type or
    ! turbulent advection option is used.  The implicit coefficients and
    ! explicit terms are calculated on thermodynamic grid levels.

    ! Interpolate <u'w'> and <v'w'> from the momentum levels to the
    ! thermodynamic levels.  These will be used for the turbulent advection
    ! terms in each equation.
    upwp_zt = zm2zt( upwp )
    vpwp_zt = zm2zt( vpwp )

    ! Implicit coefficient on <u'^2> or <v'^2> in <w'u'^2> or <w'v'^2> equation.
    coef_wpup2_wpvp2_implicit = one_third * beta * a1_zt * wp3_on_wp2_zt

    ! Explicit (RHS) term in <w'u'^2> equation.
    term_wpup2_explicit &
    = ( one - one_third * beta ) * a1_zt**2 * upwp_zt**2 &
      * wp3_on_wp2_zt / wp2_zt

    ! Explicit (RHS) term in <w'v'^2> equation.
    term_wpvp2_explicit &
    = ( one - one_third * beta ) * a1_zt**2 * vpwp_zt**2 &
      * wp3_on_wp2_zt / wp2_zt

    if ( l_upwind_xpyp_ta ) then

       ! Calculate coef_wpup2_wpvp2_implicit, term_wpup2_explicit, and
       ! term_wpvp2_explicit on momentum levels as coef_wpup2_wpvp2_implicit_zm,
       ! term_wpup2_explicit_zm, and term_wpvp2_explicit_zm, respectively.
       coef_wpup2_wpvp2_implicit_zm = one_third * beta * a1 * wp3_on_wp2

       term_wpup2_explicit_zm &
       = ( one - one_third * beta ) * a1**2 * upwp**2 * wp3_on_wp2 / wp2

       term_wpvp2_explicit_zm &
       = ( one - one_third * beta ) * a1**2 * vpwp**2 * wp3_on_wp2 / wp2

       ! For ADG1, the sign of the turbulent velocity is the sign of
       ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity is set
       ! to wp3_on_wp2.
       sgn_t_vel_up2_vp2 = wp3_on_wp2

    endif ! l_upwind_xpyp_ta

    ! Define the Coefficent of Eddy Diffusivity for the variances
    ! and covariances.
    do k = 1, gr%nz, 1

      ! Kw2 is used for variances and covariances rtp2, thlp2, rtpthlp, and
      ! passive scalars.  The variances and covariances are located on the
      ! momentum levels.  Kw2 is located on the thermodynamic levels.
      ! Kw2 = c_K2 * Kh_zt
      Kw2(k) = c_K2 * Kh_zt(k)

      ! Kw9 is used for variances up2 and vp2.  The variances are located on
      ! the momentum levels.  Kw9 is located on the thermodynamic levels.
      ! Kw9 = c_K9 * Kh_zt
      Kw9(k) = c_K9 * Kh_zt(k)

    enddo

    !!!!!***** r_t'^2 *****!!!!!

    ! Implicit contributions to term rtp2
    call xp2_xpyp_lhs( dt, l_iter, & ! In
                       coef_wprtp2_implicit, coef_wprtp2_implicit_zm, & ! In
                       sgn_t_vel_rtp2, tau_zm, wm_zm, Kw2, & ! In
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                       C2rt_1d, nu2_vert_res_dep, & ! In
                       lhs ) ! Out


    call xp2_xpyp_rhs( xp2_xpyp_rtp2, dt, l_iter, & ! In
                       coef_wprtp2_implicit, coef_wprtp2_implicit_zm, & ! In
                       term_wprtp2_explicit, term_wprtp2_explicit_zm, & ! In
                       sgn_t_vel_rtp2, wprtp, wprtp, & ! In
                       rtm, rtm, rtp2, rtp2_forcing, & ! In
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                       C2rt_1d, tau_zm, rt_tol**2, & ! In
                       rhs ) ! Out

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_rtp2, 1, & ! Intent(in)
                         rhs, lhs, rtp2 )    ! Intent(inout)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_rtp2, rtp2 ) ! Intent(in)
    end if

    !!!!!***** th_l'^2 *****!!!!!

    ! Implicit contributions to term thlp2
    call xp2_xpyp_lhs( dt, l_iter, & ! In
                       coef_wpthlp2_implicit, coef_wpthlp2_implicit_zm, & ! In
                       sgn_t_vel_thlp2, tau_zm, wm_zm, Kw2, & ! In
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                       C2thl_1d, nu2_vert_res_dep, & ! In
                       lhs ) ! Out

    ! Explicit contributions to thlp2
    call xp2_xpyp_rhs( xp2_xpyp_thlp2, dt, l_iter, & ! In
                       coef_wpthlp2_implicit, coef_wpthlp2_implicit_zm, & ! In
                       term_wpthlp2_explicit, term_wpthlp2_explicit_zm, & ! In
                       sgn_t_vel_thlp2, wpthlp, wpthlp, & ! In
                       thlm, thlm, thlp2, thlp2_forcing, & ! In
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                       C2thl_1d, tau_zm, thl_tol**2, & ! In
                       rhs ) ! Out

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_thlp2, 1, & ! Intent(in)
                         rhs, lhs, thlp2 )    ! Intent(inout)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_thlp2, thlp2 ) ! Intent(in)
    end if


    !!!!!***** r_t'th_l' *****!!!!!

    ! Implicit contributions to term rtpthlp
    call xp2_xpyp_lhs( dt, l_iter, & ! In
                       coef_wprtpthlp_implicit, coef_wprtpthlp_implicit_zm, &!In
                       sgn_t_vel_rtpthlp, tau_zm, wm_zm, Kw2, & ! In
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                       C2rtthl_1d, nu2_vert_res_dep, & ! In
                       lhs ) ! Out

    ! Explicit contributions to rtpthlp
    call xp2_xpyp_rhs( xp2_xpyp_rtpthlp, dt, l_iter, & ! In
                       coef_wprtpthlp_implicit, coef_wprtpthlp_implicit_zm, &!In
                       term_wprtpthlp_explicit, term_wprtpthlp_explicit_zm, &!In
                       sgn_t_vel_rtpthlp, wprtp, wpthlp, & ! In
                       rtm, thlm, rtpthlp, rtpthlp_forcing, & ! In
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                       C2rtthl_1d, tau_zm, zero_threshold, & ! In
                       rhs ) ! Out

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_rtpthlp, 1, & ! Intent(in)
                         rhs, lhs, rtpthlp )    ! Intent(inout)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_rtpthlp, rtpthlp ) ! Intent(in)
    end if


    !!!!!***** u'^2 / v'^2 *****!!!!!

    ! Implicit contributions to term up2/vp2
    call xp2_xpyp_lhs( dt, l_iter, & ! In
                       coef_wpup2_wpvp2_implicit, & ! In
                       coef_wpup2_wpvp2_implicit_zm, & ! In
                       sgn_t_vel_up2_vp2, tau_zm, wm_zm, Kw9, & ! In
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                       C4_C14_1d, nu9_vert_res_dep, & ! In
                       lhs ) ! Out

    ! Explicit contributions to up2
    call xp2_xpyp_uv_rhs( xp2_xpyp_up2, dt, l_iter, & ! In
                          coef_wpup2_wpvp2_implicit, & ! In
                          coef_wpup2_wpvp2_implicit_zm, & ! In
                          term_wpup2_explicit, term_wpup2_explicit_zm, & ! In
                          sgn_t_vel_up2_vp2, wp2, wp2_zt, wpthvp, & ! In
                          Lscale, C4_C14_1d, tau_zm,  & ! In
                          um, vm, upwp, vpwp, up2, vp2, & ! In
                          rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                          thv_ds_zm, C4, C5, C14, wp2_splat, & ! In
                          uv_rhs(:,1) ) ! Out

    ! Explicit contributions to vp2
    call xp2_xpyp_uv_rhs( xp2_xpyp_vp2, dt, l_iter, & ! In
                          coef_wpup2_wpvp2_implicit, & ! In
                          coef_wpup2_wpvp2_implicit_zm, & ! In
                          term_wpvp2_explicit, term_wpvp2_explicit_zm, & ! In
                          sgn_t_vel_up2_vp2, wp2, wp2_zt, wpthvp, & ! In
                          Lscale, C4_C14_1d, tau_zm,  & ! In
                          vm, um, vpwp, upwp, vp2, up2, & ! In
                          rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                          thv_ds_zm, C4, C5, C14, wp2_splat, & ! In
                          uv_rhs(:,2) ) ! Out

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_up2_vp2, 2, & ! Intent(in)
                         uv_rhs, lhs,         & ! Intent(inout)
                         uv_solution )          ! Intent(out)

    up2(1:gr%nz) = uv_solution(1:gr%nz,1)
    vp2(1:gr%nz) = uv_solution(1:gr%nz,2)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_up2, up2 ) ! Intent(in)
      call xp2_xpyp_implicit_stats( xp2_xpyp_vp2, vp2 ) ! Intent(in)
    end if


    ! Apply the positive definite scheme to variances
    if ( l_hole_fill ) then
      call pos_definite_variances( xp2_xpyp_rtp2, dt, rt_tol**2, & ! Intent(in)
                                   rho_ds_zm, rho_ds_zt, &         ! Intent(in)
                                   rtp2 )                          ! Intent(inout)
      call pos_definite_variances( xp2_xpyp_thlp2, dt, thl_tol**2, & ! Intent(in)
                                   rho_ds_zm, rho_ds_zt, &           ! Intent(in)
                                   thlp2 )                           ! Intent(inout)
      call pos_definite_variances( xp2_xpyp_up2, dt, w_tol_sqd, & ! Intent(in)
                                   rho_ds_zm, rho_ds_zt, &        ! Intent(in)
                                   up2 )                          ! Intent(inout)
      call pos_definite_variances( xp2_xpyp_vp2, dt, w_tol_sqd, & ! Intent(in)
                                   rho_ds_zm, rho_ds_zt, &        ! Intent(in)
                                   vp2 )                          ! Intent(inout)
    endif


    ! Clipping for r_t'^2

    !threshold = zero_threshold
    !
    !where ( wp2 >= w_tol_sqd ) &
    !   threshold = rt_tol*rt_tol

    ! The value of rtp2 is not allowed to become smaller than the threshold
    ! value of rt_tol^2.  Additionally, that threshold value may be boosted at
    ! any grid level in order to keep the overall correlation of w and rt
    ! between the values of -max_mag_correlation_flux and
    ! max_mag_correlation_flux by boosting rtp2 rather than by limiting the
    ! magnitude of wprtp.
    if ( l_min_xp2_from_corr_wx ) then

       ! The overall correlation of w and rt is:
       !
       ! corr_w_rt = wprtp / ( sqrt( wp2 ) * sqrt( rtp2 ) ).
       !
       ! Squaring both sides, the equation becomes:
       !
       ! corr_w_rt^2 = wprtp^2 / ( wp2 * rtp2 ).
       !
       ! Using max_mag_correlation_flux for the correlation and then solving for
       ! the minimum of rtp2, the equation becomes:
       !
       ! rtp2|_min = wprtp^2 / ( wp2 * max_mag_correlation_flux^2 ).
       do k = 1, gr%nz, 1

          threshold &
          = max( rt_tol**2, &
                 wprtp(k)**2 / ( wp2(k) * max_mag_correlation_flux**2 ) )

          call clip_variance_level( xp2_xpyp_rtp2, dt, threshold, k, & ! In
                                    rtp2(k) )                          ! In/out

       enddo ! k = 1, gr%nz, 1

    else

       ! Consider only the minimum tolerance threshold value for rtp2.
       threshold = rt_tol**2

       call clip_variance( xp2_xpyp_rtp2, dt, threshold, & ! Intent(in)
                           rtp2 )                          ! Intent(inout)

    endif ! l_min_xp2_from_corr_wx

    ! Special clipping on the variance of rt to prevent a large variance at
    ! higher altitudes.  This is done because we don't want the PDF to extend
    ! into the negative, and found that for latin hypercube sampling a large
    ! variance aloft leads to negative samples of total water.
    ! -dschanen 8 Dec 2010
    if ( l_clip_large_rtp2 ) then
    
      ! This overwrites stats clipping data from clip_variance
      if ( l_stats_samp ) then
        call stat_modify( irtp2_cl, -rtp2 / dt, stats_zm )
      endif
      
      do k = 1, gr%nz
        threshold = rtp2_clip_coef * rtm(k)**2
        if ( rtp2(k) > threshold ) then
          rtp2(k) = threshold
        end if
      end do ! k = 1..gr%nz
      
      if ( l_stats_samp ) then
        call stat_modify( irtp2_cl, rtp2 / dt, stats_zm )
      endif
      
    end if ! l_clip_large_rtp2
 

    ! Clipping for th_l'^2

    !threshold = zero_threshold
    !
    !where ( wp2 >= w_tol_sqd ) &
    !   threshold = thl_tol*thl_tol

    ! The value of thlp2 is not allowed to become smaller than the threshold
    ! value of thl_tol^2.  Additionally, that threshold value may be boosted at
    ! any grid level in order to keep the overall correlation of w and theta-l
    ! between the values of -max_mag_correlation_flux and
    ! max_mag_correlation_flux by boosting thlp2 rather than by limiting the
    ! magnitude of wpthlp.
    if ( l_min_xp2_from_corr_wx ) then

       ! The overall correlation of w and theta-l is:
       !
       ! corr_w_thl = wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) ).
       !
       ! Squaring both sides, the equation becomes:
       !
       ! corr_w_thl^2 = wpthlp^2 / ( wp2 * thlp2 ).
       !
       ! Using max_mag_correlation_flux for the correlation and then solving for
       ! the minimum of thlp2, the equation becomes:
       !
       ! thlp2|_min = wpthlp^2 / ( wp2 * max_mag_correlation_flux^2 ).
       do k = 1, gr%nz, 1

          threshold &
          = max( thl_tol**2, &
                 wpthlp(k)**2 / ( wp2(k) * max_mag_correlation_flux**2 ) )

          call clip_variance_level( xp2_xpyp_thlp2, dt, threshold, k, & ! In
                                    thlp2(k) )                          ! In/out

       enddo ! k = 1, gr%nz, 1

    else

       ! Consider only the minimum tolerance threshold value for thlp2.
       threshold = thl_tol**2

       call clip_variance( xp2_xpyp_thlp2, dt, threshold, & ! Intent(in)
                           thlp2 )                          ! Intent(inout)

    endif ! l_min_xp2_from_corr_wx


    ! Clipping for u'^2

    ! Clip negative values of up2
    threshold = w_tol_sqd
    call clip_variance( xp2_xpyp_up2, dt, threshold, & ! Intent(in)
                        up2 )                          ! Intent(inout)

    ! Clip excessively large values of up2
    if ( l_stats_samp ) then
      ! Store previous value in order to calculate clipping
      call stat_modify( iup2_cl, -up2 / dt, &   ! Intent(in)
                              stats_zm )             ! Intent(inout)
    endif

    where (up2 > 1000._core_rknd)
      up2 = 1000._core_rknd
    end where

    if ( l_stats_samp ) then
      ! Store final value in order to calculate clipping
      call stat_modify( iup2_cl, up2 / dt, &   ! Intent(in)
                              stats_zm )             ! Intent(inout)
    endif

    ! Clipping for v'^2

    ! Clip negative values of vp2
    threshold = w_tol_sqd
    call clip_variance( xp2_xpyp_vp2, dt, threshold, & ! Intent(in)
                        vp2 )                          ! Intent(inout)

    if ( l_stats_samp ) then
      ! Store previous value in order to calculate clipping
      call stat_modify( ivp2_cl, -vp2 / dt, &   ! Intent(in)
                              stats_zm )             ! Intent(inout)
    endif

    where (vp2 > 1000._core_rknd)
      vp2 = 1000._core_rknd
    end where

    if ( l_stats_samp ) then
      ! Store final value in order to calculate clipping
      call stat_modify( ivp2_cl, vp2 / dt, &   ! Intent(in)
                              stats_zm )             ! Intent(inout)
    endif

    ! When selected, apply sponge damping after up2 and vp2 have been advanced.
    if ( up2_vp2_sponge_damp_settings%l_sponge_damping ) then

       if ( l_stats_samp ) then
          call stat_begin_update( iup2_sdmp, up2 / dt, stats_zm )
          call stat_begin_update( ivp2_sdmp, vp2 / dt, stats_zm )
       endif

       up2 = sponge_damp_xp2( dt, gr%zm, up2, w_tol_sqd, &
                              up2_vp2_sponge_damp_profile )

       vp2 = sponge_damp_xp2( dt, gr%zm, vp2, w_tol_sqd, &
                              up2_vp2_sponge_damp_profile )

       if ( l_stats_samp ) then
          call stat_end_update( iup2_sdmp, up2 / dt, stats_zm )
          call stat_end_update( ivp2_sdmp, vp2 / dt, stats_zm )
       endif

    endif ! up2_vp2_sponge_damp_settings%l_sponge_damping


    ! Clipping for r_t'th_l'
    ! Clipping r_t'th_l' at each vertical level, based on the
    ! correlation of r_t and th_l at each vertical level, such that:
    ! corr_(r_t,th_l) = r_t'th_l' / [ sqrt(r_t'^2) * sqrt(th_l'^2) ];
    ! -1 <= corr_(r_t,th_l) <= 1.
    ! Since r_t'^2, th_l'^2, and r_t'th_l' are all computed in the
    ! same place, clipping for r_t'th_l' only has to be done once.
    l_first_clip_ts = .true.
    l_last_clip_ts = .true.
    call clip_covar( xp2_xpyp_rtpthlp, l_first_clip_ts,  & ! Intent(in)
                     l_last_clip_ts, dt, rtp2, thlp2,  &  ! Intent(in)
                     rtpthlp, rtpthlp_chnge )     ! Intent(inout)

    if ( l_scalar_calc ) then

      ! Implicit contributions to passive scalars

      !!!!!***** sclr'^2, sclr'r_t', sclr'th_l' *****!!!!!

      ! Set up the implicit coefficients for turbulent advection of <sclr'^2>,
      ! <sclr'rt'>, and <sclr'thl'>.
      if ( l_explicit_turbulent_adv_xpyp ) then

         ! The turbulent advection of <sclr'x'> is handled explicitly.

         ! The value of coef_wpsclrpxp_implicit is always 0.
         coef_wpsclrpxp_implicit = zero

         if ( l_upwind_xpyp_ta ) then

            ! The value of coef_wpsclrpxp_implicit_zm is always 0.
            coef_wpsclrpxp_implicit_zm = zero

            ! The sign of the turbulent velocity doesn't matter for the implicit
            ! component of turbulent advection.
            sgn_t_vel_sclrpxp = one 

         endif ! l_upwind_xpyp_ta

      else ! .not. l_explicit_turbulent_adv_xpyp

         ! The turbulent advection of <sclr'x'> is handled implicitly or
         ! semi-implicitly.

         if ( iiPDF_type == iiPDF_ADG1 ) then

            ! The ADG1 PDF is used.

            ! Calculate the implicit coefficients on thermodynamic grid levels.

            ! Implicit coefficient on <sclr'x'> in <w'sclr'x'> equation.
            coef_wpsclrpxp_implicit = one_third * beta * a1_zt * wp3_on_wp2_zt

            if ( l_upwind_xpyp_ta ) then

               ! Calculate coef_wpsclrpxp_implicit on momentum levels as
               ! coef_wpsclrpxp_implicit_zm.
               coef_wpsclrpxp_implicit_zm = one_third * beta * a1 * wp3_on_wp2

               ! For ADG1, the sign of the turbulent velocity is the sign of
               ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
               ! velocity is set to wp3_on_wp2.
               sgn_t_vel_sclrpxp = wp3_on_wp2

            endif ! l_upwind_xpyp_ta

         elseif ( iiPDF_type == iiPDF_new ) then

            ! The new PDF is used.

            ! The code for the scalar variables will be set up later.
            coef_wpsclrpxp_implicit = zero

            if ( l_upwind_xpyp_ta ) then
               coef_wpsclrpxp_implicit_zm = zero
               sgn_t_vel_sclrpxp = one
            endif ! l_upwind_xpyp_ta

         endif ! iiPDF_type

      endif ! l_explicit_turbulent_adv_xpyp

      call xp2_xpyp_lhs( dt, l_iter, & ! In
                         coef_wpsclrpxp_implicit, & ! In
                         coef_wpsclrpxp_implicit_zm, & ! In
                         sgn_t_vel_sclrpxp, tau_zm, wm_zm, Kw2, & ! In
                         rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                         C2sclr_1d, nu2_vert_res_dep, & ! In
                         lhs ) ! Out


      ! Explicit contributions to passive scalars

      do i = 1, sclr_dim, 1

        ! Set up the explicit terms for turbulent advection of <sclr'^2>.
        ! <sclr'rt'>, and <sclr'thl'>.
        if ( l_explicit_turbulent_adv_xpyp ) then

           ! The turbulent advection of <sclr'x'> is handled explicitly.

           ! The <sclr'x'> turbulent advection term is entirely explicit, as
           ! term_wpsclrpxp_explicit is equal to <w'sclr'x'> as calculated using
           ! PDF parameters, which is general for any PDF type.  The value of
           ! <w'sclr'x'> is calculated on thermodynamic levels.
           term_wpsclrp2_explicit = wpsclrp2(:,i)
           term_wpsclrprtp_explicit = wpsclrprtp(:,i)
           term_wpsclrpthlp_explicit = wpsclrpthlp(:,i)

           if ( l_upwind_xpyp_ta ) then

              ! Interpolate term_wpsclrpxp_explicit to momentum levels as
              ! term_wpsclrpxp_explicit_zm.
              term_wpsclrp2_explicit_zm = zt2zm( term_wpsclrp2_explicit )
              term_wpsclrprtp_explicit_zm = zt2zm( term_wpsclrprtp_explicit )
              term_wpsclrpthlp_explicit_zm = zt2zm( term_wpsclrpthlp_explicit )

              ! Calculate the sign of the turbulent velocity for <sclr'x'>.
              sgn_t_vel_sclrp2 &
              = sgn_turbulent_velocity( term_wpsclrp2_explicit_zm, sclrp2(:,i) )

              sgn_t_vel_sclrprtp &
              = sgn_turbulent_velocity( term_wpsclrprtp_explicit_zm, &
                                        sclrprtp(:,i) )

              sgn_t_vel_sclrpthlp &
              = sgn_turbulent_velocity( term_wpsclrpthlp_explicit_zm, &
                                        sclrpthlp(:,i) )

           endif ! l_upwind_xpyp_ta

        else ! .not. l_explicit_turbulent_adv_xpyp

           ! The turbulent advection of <sclr'x'> is handled implicitly or
           ! semi-implicitly.

           if ( iiPDF_type == iiPDF_ADG1 ) then

              ! The ADG1 PDF is used.

              ! Calculate the explicit terms on thermodynamic grid levels.

              ! Interpolate <w'sclr'> from the momentum levels to the
              ! thermodynamic levels.  It will be used for the turbulent
              ! advection terms for the ADG1 PDF.
              wpsclrp_zt = zm2zt( wpsclrp(:,i) )

              ! Explicit (RHS) term in <w'sclr'x'> equation.
              term_wpsclrp2_explicit &
              = ( one - one_third * beta ) * a1_zt**2 * wpsclrp_zt**2 &
                * wp3_on_wp2_zt / wp2_zt

              term_wpsclrprtp_explicit &
              = ( one - one_third * beta ) * a1_zt**2 * wpsclrp_zt * wprtp_zt &
                * wp3_on_wp2_zt / wp2_zt

              term_wpsclrpthlp_explicit &
              = ( one - one_third * beta ) * a1_zt**2 * wpsclrp_zt * wpthlp_zt &
                * wp3_on_wp2_zt / wp2_zt

              if ( l_upwind_xpyp_ta ) then

                 ! Calculate term_wpsclrpxp_explicit on momentum levels as
                 ! term_wpsclrpxp_explicit_zm.
                 term_wpsclrp2_explicit_zm &
                 = ( one - one_third * beta ) * a1**2 * wpsclrp(:,i)**2 &
                   * wp3_on_wp2 / wp2

                 term_wpsclrprtp_explicit_zm &
                 = ( one - one_third * beta ) * a1**2 * wpsclrp(:,i) * wprtp &
                   * wp3_on_wp2 / wp2

                 term_wpsclrpthlp_explicit_zm &
                 = ( one - one_third * beta ) * a1**2 * wpsclrp(:,i) * wpthlp &
                   * wp3_on_wp2 / wp2

                 ! For ADG1, the sign of the turbulent velocity is the sign of
                 ! <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
                 ! velocity is set to wp3_on_wp2.
                 sgn_t_vel_sclrp2 = wp3_on_wp2
                 sgn_t_vel_sclrprtp = wp3_on_wp2
                 sgn_t_vel_sclrpthlp = wp3_on_wp2

              endif ! l_upwind_xpyp_ta

           elseif ( iiPDF_type == iiPDF_new ) then

              ! The new PDF is used.

              ! The code for the scalar variables will be set up later.
              term_wpsclrp2_explicit = zero
              term_wpsclrprtp_explicit = zero
              term_wpsclrpthlp_explicit = zero

              if ( l_upwind_xpyp_ta ) then

                 term_wpsclrp2_explicit_zm = zero
                 term_wpsclrprtp_explicit_zm = zero
                 term_wpsclrpthlp_explicit_zm = zero

                 sgn_t_vel_sclrp2 = one
                 sgn_t_vel_sclrprtp = one
                 sgn_t_vel_sclrpthlp = one

              endif ! l_upwind_xpyp_ta

           endif ! iiPDF_type

        endif ! l_explicit_turbulent_adv_xpyp

        ! Forcing for <sclr'^2>.
        sclrp2_forcing = zero

        !!!!!***** sclr'^2 *****!!!!!

        call xp2_xpyp_rhs( xp2_xpyp_sclrp2, dt, l_iter, & ! In
                           coef_wpsclrpxp_implicit, & ! In
                           coef_wpsclrpxp_implicit_zm, & ! In
                           term_wpsclrp2_explicit, & ! In
                           term_wpsclrp2_explicit_zm, & ! In
                           sgn_t_vel_sclrp2, wpsclrp(:,i), wpsclrp(:,i), & ! In
                           sclrm(:,i), sclrm(:,i), & ! In
                           sclrp2(:,i), sclrp2_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                           C2sclr_1d, tau_zm, sclr_tol(i)**2, & ! In
                           sclr_rhs(:,i) ) ! Out

        !!!!!***** sclr'r_t' *****!!!!!
        if ( i == iisclr_rt ) then
           ! In this case we're trying to emulate rt'^2 with sclr'rt', so we
           ! handle this as we would a variance, even though generally speaking
           ! the scalar is not rt
           sclrprtp_forcing = rtp2_forcing
           threshold = rt_tol**2
        else
           sclrprtp_forcing = zero
           threshold = zero_threshold
        endif

        call xp2_xpyp_rhs( xp2_xpyp_sclrprtp, dt, l_iter, & ! In
                           coef_wpsclrpxp_implicit, & ! In
                           coef_wpsclrpxp_implicit_zm, & ! In
                           term_wpsclrprtp_explicit, & ! In
                           term_wpsclrprtp_explicit_zm, & ! In
                           sgn_t_vel_sclrprtp, wpsclrp(:,i), wprtp, & ! In
                           sclrm(:,i), rtm, sclrprtp(:,i), & ! In
                           sclrprtp_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                           C2sclr_1d, tau_zm, threshold, & ! In
                           sclr_rhs(:,i+sclr_dim) ) ! Out

        !!!!!***** sclr'th_l' *****!!!!!

        if ( i == iisclr_thl ) then
           ! In this case we're trying to emulate thl'^2 with sclr'thl', so we
           ! handle this as we did with sclr_rt, above.
           sclrpthlp_forcing = thlp2_forcing
           threshold = thl_tol**2
        else
           sclrpthlp_forcing = zero
           threshold = zero_threshold
        endif

        call xp2_xpyp_rhs( xp2_xpyp_sclrpthlp, dt, l_iter, & ! In
                           coef_wpsclrpxp_implicit, & ! In
                           coef_wpsclrpxp_implicit_zm, & ! In
                           term_wpsclrpthlp_explicit, & ! In
                           term_wpsclrpthlp_explicit_zm, & ! In
                           sgn_t_vel_sclrpthlp, wpsclrp(:,i), wpthlp, & ! In
                           sclrm(:,i), thlm, sclrpthlp(:,i), & ! In
                           sclrpthlp_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                           C2sclr_1d, tau_zm, threshold, & ! In
                           sclr_rhs(:,i+2*sclr_dim) ) ! Out

      enddo ! 1..sclr_dim


      ! Solve the tridiagonal system

      call xp2_xpyp_solve( xp2_xpyp_scalars, 3*sclr_dim, &   ! Intent(in)
                           sclr_rhs, lhs, sclr_solution )    ! Intent(inout)

      sclrp2(:,1:sclr_dim) = sclr_solution(:,1:sclr_dim)

      sclrprtp(:,1:sclr_dim) = sclr_solution(:,sclr_dim+1:2*sclr_dim)

      sclrpthlp(:,1:sclr_dim) = sclr_solution(:,2*sclr_dim+1:3*sclr_dim)

      ! Apply hole filling algorithm to the scalar variance terms
      if ( l_hole_fill ) then
        do i = 1, sclr_dim, 1
          call pos_definite_variances( xp2_xpyp_sclrp2, dt, sclr_tol(i)**2, & ! Intent(in)
                                       rho_ds_zm, rho_ds_zt, &                ! Intent(in)
                                       sclrp2(:,i) )                          ! Intent(inout)
          if ( i == iisclr_rt ) then
            ! Here again, we do this kluge here to make sclr'rt' == rt'^2
            call pos_definite_variances( xp2_xpyp_sclrprtp, dt, sclr_tol(i)**2, & ! Intent(in)
                                         rho_ds_zm, rho_ds_zt, &                  ! Intent(in)
                                         sclrprtp(:,i) )                          ! Intent(inout)
          end if
          if ( i == iisclr_thl ) then
            ! As with sclr'rt' above, but for sclr'thl'
            call pos_definite_variances( xp2_xpyp_sclrpthlp, dt, sclr_tol(i)**2, & ! Intent(in)
                                         rho_ds_zm, rho_ds_zt, &                   ! Intent(in)
                                         sclrpthlp(:,i) )                          ! Intent(inout)
          end if
        enddo
      endif


      ! Clipping for sclr'^2
      do i = 1, sclr_dim, 1

!      threshold = zero_threshold
!
!      where ( wp2 >= w_tol_sqd ) &
!         threshold = sclr_tol(i)*sclr_tol(i)

        threshold = sclr_tol(i)**2

        call clip_variance( clip_sclrp2, dt, threshold, & ! Intent(in)
                            sclrp2(:,i) )                 ! Intent(inout)

      enddo


      ! Clipping for sclr'r_t'
      ! Clipping sclr'r_t' at each vertical level, based on the
      ! correlation of sclr and r_t at each vertical level, such that:
      ! corr_(sclr,r_t) = sclr'r_t' / [ sqrt(sclr'^2) * sqrt(r_t'^2) ];
      ! -1 <= corr_(sclr,r_t) <= 1.
      ! Since sclr'^2, r_t'^2, and sclr'r_t' are all computed in the
      ! same place, clipping for sclr'r_t' only has to be done once.
      do i = 1, sclr_dim, 1

        if  ( i == iisclr_rt ) then
          ! Treat this like a variance if we're emulating rt
          threshold = sclr_tol(i) * rt_tol

          call clip_variance( clip_sclrprtp, dt, threshold, & ! Intent(in)
                              sclrprtp(:,i) )                 ! Intent(inout)
        else
          l_first_clip_ts = .true.
          l_last_clip_ts = .true.
          call clip_covar( clip_sclrprtp, l_first_clip_ts,  &            ! Intent(in) 
                           l_last_clip_ts, dt, sclrp2(:,i), rtp2(:), &  ! Intent(in)
                           sclrprtp(:,i), sclrprtp_chnge(:,i) ) ! Intent(inout)
        end if
      enddo


      ! Clipping for sclr'th_l'
      ! Clipping sclr'th_l' at each vertical level, based on the
      ! correlation of sclr and th_l at each vertical level, such that:
      ! corr_(sclr,th_l) = sclr'th_l' / [ sqrt(sclr'^2) * sqrt(th_l'^2) ];
      ! -1 <= corr_(sclr,th_l) <= 1.
      ! Since sclr'^2, th_l'^2, and sclr'th_l' are all computed in the
      ! same place, clipping for sclr'th_l' only has to be done once.
      do i = 1, sclr_dim, 1
        if ( i == iisclr_thl ) then
          ! As above, but for thl
          threshold = sclr_tol(i) * thl_tol
          call clip_variance( clip_sclrpthlp, dt, threshold, & ! Intent(in)
                              sclrpthlp(:,i) )                 ! Intent(inout)
        else
          l_first_clip_ts = .true.
          l_last_clip_ts = .true.
          call clip_covar( clip_sclrpthlp, l_first_clip_ts,  &            ! Intent(in) 
                           l_last_clip_ts, dt, sclrp2(:,i), thlp2(:), &   ! Intent(in) 
                           sclrpthlp(:,i), sclrpthlp_chnge(:,i) ) ! Intent(inout)
        end if
      enddo

    endif ! l_scalar_calc

    if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then

          write(fstderr,*) "Error in advance_xp2_xpyp"

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "tau_zm = ", tau_zm
          write(fstderr,*) "wm_zm = ", wm_zm
          write(fstderr,*) "rtm = ", rtm
          write(fstderr,*) "wprtp = ", wprtp
          write(fstderr,*) "thlm = ", thlm
          write(fstderr,*) "wpthlp = ", wpthlp
          write(fstderr,*) "wpthvp = ", wpthvp
          write(fstderr,*) "um = ", um
          write(fstderr,*) "vm = ", vm
          write(fstderr,*) "wp2 = ", wp2
          write(fstderr,*) "wp3 = ", wp3
          write(fstderr,*) "upwp = ", upwp
          write(fstderr,*) "vpwp = ", vpwp
          write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
          write(fstderr,*) "Skw_zm = ", Skw_zm
          write(fstderr,*) "Kh_zt = ", Kh_zt
          write(fstderr,*) "rtp2_forcing = ", rtp2_forcing
          write(fstderr,*) "thlp2_forcing = ", thlp2_forcing
          write(fstderr,*) "rtpthlp_forcing = ", rtpthlp_forcing
          write(fstderr,*) "rho_ds_zm = ", rho_ds_zm
          write(fstderr,*) "rho_ds_zt = ", rho_ds_zt
          write(fstderr,*) "invrs_rho_ds_zm = ", invrs_rho_ds_zm
          write(fstderr,*) "thv_ds_zm = ", thv_ds_zm
          write(fstderr,*) "wp2_zt = ", wp2_zt

          do i = 1, sclr_dim
            write(fstderr,*) "sclrm = ", i, sclrm(:,i)
            write(fstderr,*) "wpsclrp = ", i, wpsclrp(:,i)
          enddo

          write(fstderr,*) "Intent(In/Out)"

          write(fstderr,*) "rtp2 = ", rtp2
          write(fstderr,*) "thlp2 = ", thlp2
          write(fstderr,*) "rtpthlp = ", rtpthlp
          write(fstderr,*) "up2 = ", up2
          write(fstderr,*) "vp2 = ", vp2

          do i = 1, sclr_dim
            write(fstderr,*) "sclrp2 = ", i, sclrp2(:,i)
            write(fstderr,*) "sclrprtp = ", i, sclrprtp(:,i)
            write(fstderr,*) "sclrthlp = ", i, sclrpthlp(:,i)
          enddo

        endif
    end if

    return
  end subroutine advance_xp2_xpyp

  !=============================================================================
  subroutine xp2_xpyp_lhs( dt, l_iter, & ! In
                           coef_wpxpyp_implicit, coef_wpxpyp_implicit_zm, & ! In
                           sgn_turbulent_vel, tau_zm, wm_zm, Kw, & ! In
                           rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, Cn, nu, & ! In
                           lhs ) ! Out

  ! Description:
  !     Compute LHS tridiagonal matrix for a variance or covariance term
  ! 
  ! Notes:
  !     For LHS turbulent advection terms:
  !         An "over-implicit" weighted time step is applied to this term.
  !         The weight of the implicit portion of this term is controlled by
  !         the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  !         expression below).  A factor is added to the right-hand side of
  !         the equation in order to balance a weight that is not equal to 1,
  !         such that:
  !              -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  !         where X is the variable that is being solved for in a predictive
  !         equation (<x'^2> or <x'y'> in this case), y(t) is the linearized
  !         portion of the term that gets treated implicitly, and RHS is the
  !         portion of the term that is always treated explicitly.  A weight
  !         of greater than 1 can be applied to make the term more
  !         numerically stable.
  ! 
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

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs_all, & ! Procedure(s)
        xpyp_term_ta_pdf_lhs

    use mean_adv, only:  & 
        term_ma_zm_lhs_all, & ! Procedure(s)
        term_ma_zm_lhs

    use diffusion, only:  & 
        diffusion_zm_lhs_all, & ! Procedure(s)
        diffusion_zm_lhs

    use model_flags, only: &
        l_upwind_xpyp_ta    ! Variable(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use advance_helper_module, only: &
        set_boundary_conditions_lhs    ! Procedure(s)

    
    use stats_variables, only: &
        l_stats_samp, &
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
        irtp2_ma, &
        irtp2_ta, &
        irtp2_dp1, &
        irtp2_dp2, &
        ithlp2_ma, &
        ithlp2_ta, &
        ithlp2_dp1, &
        ithlp2_dp2, &
        irtpthlp_ma, &
        irtpthlp_ta, &
        irtpthlp_dp1, &
        irtpthlp_dp2, &
        iup2_ma, &
        iup2_ta, &
        iup2_dp2, &
        ivp2_ma, &
        ivp2_ta, &
        ivp2_dp2

    implicit none

    !------------------- Input Variables -------------------
    real( kind = core_rknd ), intent(in) :: & 
      dt ! Timestep length                             [s]

    logical, intent(in) :: & 
      l_iter  ! Whether the variances are prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      coef_wpxpyp_implicit,    & ! Coef. of <x'y'> in <w'x'y'> eq.; t-levs [m/s]
      coef_wpxpyp_implicit_zm, & ! coef_wpxpyp_implicit interp. to m-levs  [m/s]
      sgn_turbulent_vel,       & ! Sign of turbulent velocity ("upwind" ta)  [-]
      tau_zm,                  & ! Time-scale tau on momentum levels         [s]
      wm_zm,                   & ! w wind component on momentum levels     [m/s]
      Kw,                      & ! Coef. of eddy diffusivity (all vars.) [m^2/s]
      rho_ds_zt,               & ! Dry, static density on thermo. levs. [kg/m^3]
      rho_ds_zm,               & ! Dry, static density on m-levs.       [kg/m^3]
      invrs_rho_ds_zm,         & ! Inv. dry, static density on m-levs.  [m^3/kg]
      Cn,                      & ! Coefficient C_n                           [-]
      nu                         ! Background const. coef. of eddy diff. [m^2/s]

    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(3,gr%nz), intent(out) :: & 
      lhs         ! Implicit contributions to the term

    !---------------- Local Variables -------------------
    real( kind = core_rknd ), dimension(3,gr%nz) :: & 
      lhs_diff, & ! Diffusion contributions to lhs, dissipation term 2
      lhs_ma, & ! Mean advection contributions to lhs
      lhs_ta    ! Turbulent advection contributions to lhs
    
    real( kind = core_rknd ), dimension(gr%nz) :: &
      lhs_dp1   ! LHS dissipation term 1

    ! Array indices
    integer :: k, low_bound, high_bound

    !---------------- Begin Code -------------------

    ! LHS dissipation term 1 (dp1). An "over-implicit" weighted time 
    ! step is applied to this term (and to pressure term 1 for u'^2 and v'^2).
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xp2_dp
    do k = 2, gr%nz-1
        lhs_dp1(k) = term_dp1_lhs( Cn(k), tau_zm(k) ) * gamma_over_implicit_ts
    enddo ! k=2..gr%nz-1

    ! Calculate LHS eddy diffusion term: dissipation term 2 (dp2).
    call diffusion_zm_lhs_all( Kw(:), nu(:),                     & ! Intent(in)
                               gr%invrs_dzt(:), gr%invrs_dzm(:), & ! Intent(in)
                               lhs_diff(:,:)                     ) ! Intent(out) 

    ! Calculate LHS mean advection (ma) term.
    call term_ma_zm_lhs_all( wm_zm(:), gr%invrs_dzm(:), & ! Intent(in)
                             lhs_ma(:,:)              ) ! Intent(out)
                             

    ! Calculate LHS turbulent advection (ta) terms
    call xpyp_term_ta_pdf_lhs_all( coef_wpxpyp_implicit(:),             & ! Intent(in)
                                   rho_ds_zt(:),                        & ! Intent(in)
                                   invrs_rho_ds_zm(:),                  & ! Intent(in)
                                   gr%invrs_dzm(:),                     & ! Intent(in)
                                   l_upwind_xpyp_ta,                    & ! Intent(in)
                                   sgn_turbulent_vel(:),                & ! Intent(in)
                                   coef_wpxpyp_implicit_zm(:),          & ! Intent(in)
                                   rho_ds_zm(:),                        & ! Intent(in)
                                   gr%invrs_dzt(:),                     & ! Intent(in)
                                   lhs_ta(:,:)                        ) ! Intent(out)

    ! Combine all lhs terms into lhs, should be fully vectorized
    do k = 2, gr%nz-1

      lhs(1,k) = lhs_diff(1,k) + lhs_ma(1,k) + lhs_ta(1,k) * gamma_over_implicit_ts

      lhs(2,k) = lhs_diff(2,k) + lhs_ma(2,k) + lhs_ta(2,k) * gamma_over_implicit_ts &
                                               + lhs_dp1(k)
        
      lhs(3,k) = lhs_diff(3,k) + lhs_ma(3,k) + lhs_ta(3,k) * gamma_over_implicit_ts

    enddo ! k=2..gr%nz-1


    ! LHS time tendency.
    if ( l_iter ) then
        do k = 2, gr%nz-1
            lhs(2,k) = lhs(2,k) + ( one / dt )
        end do
    endif

    if ( l_stats_samp ) then

        ! Statistics: implicit contributions for rtp2, thlp2,
        !             rtpthlp, up2, or vp2.

        if ( irtp2_dp1 + ithlp2_dp1 + irtpthlp_dp1  > 0 ) then

            ! Note:  The statistical implicit contribution to term dp1
            !        (as well as to term pr1) for up2 and vp2 is recorded
            !        in xp2_xpyp_uv_rhs because up2 and vp2 use a special
            !        dp1/pr1 combined term.
            ! Note:  An "over-implicit" weighted time step is applied to this
            !        term.  A weighting factor of greater than 1 may be used to
            !        make the term more numerically stable (see note above for
            !        LHS turbulent advection (ta) term).
            do k = 2, gr%nz-1
                zmscr01(k) = -lhs_dp1(k)
            end do

        endif

        if ( irtp2_dp2 + ithlp2_dp2 + irtpthlp_dp2 + iup2_dp2 + ivp2_dp2 > 0 ) then

            do k = 2, gr%nz-1
                zmscr02(k) = -lhs_diff(3,k)
                zmscr03(k) = -lhs_diff(2,k)
                zmscr04(k) = -lhs_diff(1,k)
            end do

        endif

        if ( irtp2_ta + ithlp2_ta + irtpthlp_ta + iup2_ta + ivp2_ta > 0 ) then

            ! Note:  An "over-implicit" weighted time step is applied to this
            !        term.  A weighting factor of greater than 1 may be used to
            !        make the term more numerically stable (see note above for
            !        LHS turbulent advection (ta) term).
            do k = 2, gr%nz-1
                zmscr05(k) = -gamma_over_implicit_ts * lhs_ta(3,k)
                zmscr06(k) = -gamma_over_implicit_ts * lhs_ta(2,k)
                zmscr07(k) = -gamma_over_implicit_ts * lhs_ta(1,k)
            end do

        endif

        if ( irtp2_ma + ithlp2_ma + irtpthlp_ma + iup2_ma + ivp2_ma > 0 ) then

            do k = 2, gr%nz-1
                zmscr08(k) = -lhs_ma(3,k)
                zmscr09(k) = -lhs_ma(2,k)
                zmscr10(k) = -lhs_ma(1,k)
            end do

        endif

    end if


    ! Boundary Conditions
    ! These are set so that the surface_varnce value of the variances and
    ! covariances can be used at the lowest boundary and the values of those
    ! variables can be set to their respective threshold minimum values at the
    ! top boundary.  Fixed-point boundary conditions are used for both the
    ! variances and the covariances.
    low_bound = 1
    high_bound = gr%nz

    call set_boundary_conditions_lhs( 2, low_bound, high_bound, lhs )

    return

  end subroutine xp2_xpyp_lhs

  !=============================================================================
  subroutine xp2_xpyp_solve( solve_type, nrhs, rhs, lhs, xapxbp )

    ! Description:
    ! Solve a tridiagonal system
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one  ! Constant(s)

    use lapack_wrap, only:  & 
        tridag_solve,  & ! Variable(s)
        tridag_solvex !, &
!        band_solve

    use grid_class, only: & 
        gr ! Variable(s)

    use stats_type_utilities, only: & 
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only: & 
        stats_sfc, &  ! Derived type
        irtp2_matrix_condt_num, & ! Stat index Variables
        ithlp2_matrix_condt_num, & 
        irtpthlp_matrix_condt_num, & 
        iup2_vp2_matrix_condt_num, & 
        l_stats_samp  ! Logical

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    ! Input variables
    integer, intent(in) :: &
      nrhs  ! Number of right hand side vectors

    integer, intent(in) ::  & 
      solve_type ! Variable(s) description

    ! Input/Ouput variables
    real( kind = core_rknd ), dimension(gr%nz,nrhs), intent(inout) :: & 
      rhs  ! Explicit contributions to x variance/covariance term [units vary]

    real( kind = core_rknd ), dimension(3,gr%nz), intent(inout) :: & 
      lhs  ! Implicit contributions to x variance/covariance term [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz,nrhs), intent(out) ::  & 
      xapxbp ! Computed value of the variable(s) at <t+1> [units vary]

    ! Local variables
    real( kind = core_rknd ) :: rcond  ! Est. of the reciprocal of the condition # on the matrix

    integer ::  ixapxbp_matrix_condt_num ! Stat index

    character(len=10) :: &
      solve_type_str ! solve_type in string format for debug output purposes

    ! --- Begin Code ---

    select case ( solve_type )
      !------------------------------------------------------------------------
      ! Note that these are diagnostics from inverting the matrix, not a budget
      !------------------------------------------------------------------------
    case ( xp2_xpyp_rtp2 )
      ixapxbp_matrix_condt_num  = irtp2_matrix_condt_num
      solve_type_str = "rtp2"
    case ( xp2_xpyp_thlp2 )
      ixapxbp_matrix_condt_num  = ithlp2_matrix_condt_num
      solve_type_str = "thlp2"
    case ( xp2_xpyp_rtpthlp )
      ixapxbp_matrix_condt_num  = irtpthlp_matrix_condt_num
      solve_type_str = "rtpthlp"
    case ( xp2_xpyp_up2_vp2 )
      ixapxbp_matrix_condt_num  = iup2_vp2_matrix_condt_num
      solve_type_str = "up2_vp2"
    case default
      ! No condition number is setup for the passive scalars
      ixapxbp_matrix_condt_num  = 0
      solve_type_str = "scalar"
    end select

    if ( l_stats_samp .and. ixapxbp_matrix_condt_num > 0 ) then
      call tridag_solvex & 
           ( solve_type_str, gr%nz, nrhs, &                                        ! Intent(in) 
             lhs(kp1_mdiag,:), lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  & ! Intent(inout)
             xapxbp(:,1:nrhs), rcond )                                             ! Intent(out)

      ! Est. of the condition number of the variance LHS matrix
      call stat_update_var_pt( ixapxbp_matrix_condt_num, 1, one / rcond, &  ! Intent(in)
                               stats_sfc )                                  ! Intent(inout)

    else
      call tridag_solve & 
           ( solve_type_str, gr%nz, nrhs, lhs(kp1_mdiag,:),  &      ! Intent(in)
             lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  &    ! Intent(inout)
             xapxbp(:,1:nrhs) )                                     ! Intent(out)
    end if

    return
  end subroutine xp2_xpyp_solve

  !=============================================================================
  subroutine xp2_xpyp_implicit_stats( solve_type, xapxbp )

    ! Description:
    ! Finalize implicit contributions for r_t'^2, th_l'^2, r_t'th_l',
    ! u'^2, and v'^2.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: &
      gr ! Derived type variable

    use stats_type_utilities, only: & 
      stat_end_update_pt, & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: &
      stats_zm,  & ! Variable(s) 
      irtp2_dp1, &
      irtp2_dp2, &
      irtp2_ta, &
      irtp2_ma, &
      ithlp2_dp1, &
      ithlp2_dp2, &
      ithlp2_ta, &
      ithlp2_ma, &
      irtpthlp_dp1, &
      irtpthlp_dp2, &
      irtpthlp_ta, &
      irtpthlp_ma, &
      iup2_dp1, &
      iup2_dp2, &
      iup2_ta, &
      iup2_ma, &
      iup2_pr1, &
      ivp2_dp1

    use stats_variables, only: &
      ivp2_dp2, &
      ivp2_ta, &
      ivp2_ma, &
      ivp2_pr1, &
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
      zmscr11

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, min, trim

    ! Input variables
    integer, intent(in) ::  & 
      solve_type ! Variable(s) description

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      xapxbp ! Computed value of the variable at <t+1> [units vary]

    ! Local variables
    integer :: k, kp1, km1 ! Array indices

    ! Budget indices
    integer :: & 
      ixapxbp_dp1, & 
      ixapxbp_dp2, & 
      ixapxbp_ta, & 
      ixapxbp_ma, & 
      ixapxbp_pr1

    ! --- Begin Code ---

    select case ( solve_type )
    case ( xp2_xpyp_rtp2 )
      ixapxbp_dp1 = irtp2_dp1
      ixapxbp_dp2 = irtp2_dp2
      ixapxbp_ta  = irtp2_ta
      ixapxbp_ma  = irtp2_ma
      ixapxbp_pr1 = 0

    case ( xp2_xpyp_thlp2 )
      ixapxbp_dp1 = ithlp2_dp1
      ixapxbp_dp2 = ithlp2_dp2
      ixapxbp_ta  = ithlp2_ta
      ixapxbp_ma  = ithlp2_ma
      ixapxbp_pr1 = 0

    case ( xp2_xpyp_rtpthlp )
      ixapxbp_dp1 = irtpthlp_dp1
      ixapxbp_dp2 = irtpthlp_dp2
      ixapxbp_ta  = irtpthlp_ta
      ixapxbp_ma  = irtpthlp_ma
      ixapxbp_pr1 = 0

    case ( xp2_xpyp_up2 )
      ixapxbp_dp1 = iup2_dp1
      ixapxbp_dp2 = iup2_dp2
      ixapxbp_ta  = iup2_ta
      ixapxbp_ma  = iup2_ma
      ixapxbp_pr1 = iup2_pr1

    case ( xp2_xpyp_vp2 )
      ixapxbp_dp1 = ivp2_dp1
      ixapxbp_dp2 = ivp2_dp2
      ixapxbp_ta  = ivp2_ta
      ixapxbp_ma  = ivp2_ma
      ixapxbp_pr1 = ivp2_pr1

    case default ! No budgets are setup for the passive scalars
      ixapxbp_dp1 = 0
      ixapxbp_dp2 = 0
      ixapxbp_ta  = 0
      ixapxbp_ma  = 0
      ixapxbp_pr1 = 0

    end select

    do k = 2, gr%nz-1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! x'y' term dp1 has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_dp1, k, &          ! Intent(in)
                               zmscr01(k) * xapxbp(k), &  ! Intent(in)
                               stats_zm )                 ! Intent(inout)

      ! x'y' term dp2 is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_dp2, k, &            ! Intent(in)
                                 zmscr02(k) * xapxbp(km1) & ! Intent(in)
                               + zmscr03(k) * xapxbp(k) & 
                               + zmscr04(k) * xapxbp(kp1), &
                               stats_zm )                   ! Intent(inout)

      ! x'y' term ta has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_ta, k, &              ! Intent(in)
                                 zmscr05(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr06(k) * xapxbp(k) &  
                               + zmscr07(k) * xapxbp(kp1), &
                               stats_zm )                    ! Intent(inout)

      ! x'y' term ma is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_ma, k, &              ! Intent(in)
                                 zmscr08(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr09(k) * xapxbp(k) & 
                               + zmscr10(k) * xapxbp(kp1), &
                               stats_zm )                    ! Intent(inout)

      ! x'y' term pr1 has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_pr1, k, &          ! Intent(in)
                               zmscr11(k) * xapxbp(k), &  ! Intent(in)
                               stats_zm )                 ! Intent(inout)

    end do ! k=2..gr%nz-1

    return
  end subroutine xp2_xpyp_implicit_stats

  !==================================================================================
  subroutine xp2_xpyp_uv_rhs( solve_type, dt, l_iter, & ! In
                              coef_wpxp2_implicit, coef_wpxp2_implicit_zm, &! In
                              term_wpxp2_explicit, term_wpxp2_explicit_zm, &! In
                              sgn_turbulent_vel, wp2, wp2_zt, wpthvp, & ! In
                              Lscale, C4_C14_1d, tau_zm,  & ! In
                              xam, xbm, wpxap, wpxbp, xap2, xbp2, & ! In
                              rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! In
                              thv_ds_zm, C4, C5, C14, wp2_splat, & ! In
                              rhs ) ! Out

  ! Description:
  !     Explicit contributions to u'^2 or v'^2
  ! 
  ! Notes:
  !    For LHS turbulent advection (ta) term:
  !        An "over-implicit" weighted time step is applied to this term.
  !        The weight of the implicit portion of this term is controlled by
  !        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  !        expression below).  A factor is added to the right-hand side of
  !        the equation in order to balance a weight that is not equal to 1,
  !        such that:
  !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  !        where X is the variable that is being solved for in a predictive
  !        equation (x'^2 or x'y' in this case), y(t) is the linearized
  !        portion of the term that gets treated implicitly, and RHS is the
  !        portion of the term that is always treated explicitly.  A weight
  !        of greater than 1 can be applied to make the term more
  !        numerically stable.
  ! 
  ! 
  !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  !   Significant changes to this routine may adversely affect computational speed
  !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  !----------------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only:  & 
        gamma_over_implicit_ts, & ! Constant(s)
        w_tol_sqd, &
        one, &
        two_thirds, &
        one_third, &
        zero

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs_all, & ! Procedure(s)
        xpyp_term_ta_pdf_lhs,     &
        xpyp_term_ta_pdf_rhs_all, &
        xpyp_term_ta_pdf_rhs

    use model_flags, only: &
        l_upwind_xpyp_ta ! Constant(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_begin_update_pt, & ! Procedure(s)
        stat_update_var_pt, &
        stat_modify_pt

    use stats_variables, only: & 
        ivp2_ta,  & ! Variable(s)
        ivp2_tp, & 
        ivp2_dp1, & 
        ivp2_pr1, & 
        ivp2_pr2, & 
        ivp2_splat, & 
        iup2_ta, & 
        iup2_tp, & 
        iup2_dp1, & 
        iup2_pr1, & 
        iup2_pr2, & 
        iup2_splat, & 
        stats_zm, & 
        zmscr01, & 
        zmscr11, & 
        l_stats_samp

    implicit none

    !------------------- Input Variables -------------------
    integer, intent(in) :: solve_type

    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                              [s]

    logical, intent(in) :: & 
      l_iter  ! Whether x is prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      coef_wpxp2_implicit,    & ! Coef. of <x'^2> in <w'x'^2> eq.; t-levs. [m/s]
      coef_wpxp2_implicit_zm, & ! coef_wpxp2_implicit interp. to m-levs.   [m/s]
      term_wpxp2_explicit,    & ! RHS term: <w'x'^2> eq.; t-levs.      [m^3/s^3]
      term_wpxp2_explicit_zm, & ! term_wpxp2_explicit interp. m-levs.  [m^3/s^3]
      sgn_turbulent_vel,      & ! Sign of turbulent velocity ("upwind" ta)   [-]
      wp2,                    & ! w'^2 (momentum levels)              [m^2/s^2]
      wp2_zt,                 & ! w'^2 interp. to thermo. levels      [m^2/s^2]
      wpthvp,                 & ! w'th_v' (momentum levels)             [K m/s]
      Lscale,                 & ! Mixing Length                             [m]
      C4_C14_1d,              & ! Combination of model params. C_4 and C_14 [-]
      tau_zm,                 & ! Time-scale tau on momentum levels         [s]
      xam,                    & ! x_am (thermodynamic levels)             [m/s]
      xbm,                    & ! x_bm (thermodynamic levels)             [m/s]
      wpxap,                  & ! w'x_a' (momentum levels)            [m^2/s^2]
      wpxbp,                  & ! w'x_b' (momentum levels)            [m^2/s^2]
      xap2,                   & ! x_a'^2 (momentum levels)            [m^2/s^2]
      xbp2,                   & ! x_b'^2 (momentum levels)            [m^2/s^2]
      rho_ds_zt,              & ! Dry, static density on thermo. levs. [kg/m^3]
      rho_ds_zm,              & ! Dry, static density on m-levs.       [kg/m^3]
      invrs_rho_ds_zm,        & ! Inv. dry, static density on m-levs.  [m^3/kg]
      thv_ds_zm,              & ! Dry, base-state theta_v on momentum levs. [K]
      wp2_splat    ! Tendency of <w'^2> due to splatting of eddies  [m^2/s^3]

    real( kind = core_rknd ), intent(in) :: & 
      C4,  & ! Model parameter C_4                         [-]
      C5,  & ! Model parameter C_5                         [-]
      C14    ! Model parameter C_{14}                      [-]

    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      rhs    ! Explicit contributions to x variance/covariance terms

    !---------------- Local Variables -------------------

    ! Array indices
    integer :: k

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3,gr%nz) :: lhs_ta

    real( kind = core_rknd ), dimension(gr%nz) :: rhs_ta  ! Turbulent advection terms

    real( kind = core_rknd ) :: tmp

    integer :: & 
      ixapxbp_ta, & 
      ixapxbp_tp, & 
      ixapxbp_dp1, & 
      ixapxbp_pr1, & 
      ixapxbp_pr2, &
      ixapxbp_splat

    !----------------------------- Begin Code ----------------------------------

    select case ( solve_type )
    case ( xp2_xpyp_vp2 )
      ixapxbp_ta  = ivp2_ta
      ixapxbp_tp  = ivp2_tp
      ixapxbp_dp1 = ivp2_dp1
      ixapxbp_pr1 = ivp2_pr1
      ixapxbp_pr2 = ivp2_pr2
      ixapxbp_splat = ivp2_splat
    case ( xp2_xpyp_up2 )
      ixapxbp_ta  = iup2_ta
      ixapxbp_tp  = iup2_tp
      ixapxbp_dp1 = iup2_dp1
      ixapxbp_pr1 = iup2_pr1
      ixapxbp_pr2 = iup2_pr2
      ixapxbp_splat = iup2_splat
    case default ! No budgets for passive scalars
      ixapxbp_ta  = 0
      ixapxbp_tp  = 0
      ixapxbp_dp1 = 0
      ixapxbp_pr1 = 0
      ixapxbp_pr2 = 0
      ixapxbp_splat = 0
    end select

    
    ! Calculate RHS turbulent advection (ta) terms.
    call xpyp_term_ta_pdf_rhs_all( term_wpxp2_explicit(:),            & ! Intent(in)
                                   rho_ds_zt(:),                      & ! Intent(in)
                                   invrs_rho_ds_zm(:),                & ! Intent(in)
                                   gr%invrs_dzm(:),                   & ! Intent(in)
                                   l_upwind_xpyp_ta,                  & ! Intent(in)
                                   sgn_turbulent_vel(:),              & ! Intent(in)
                                   term_wpxp2_explicit_zm(:),         & ! Intent(in)
                                   rho_ds_zm(:),                      & ! Intent(in)
                                   gr%invrs_dzt(:),                   & ! Intent(in)
                                   rhs_ta(:)                        ) ! Intent(out

    ! Calculate RHS contribution from "over-implicit" weighted time step
    ! for LHS turbulent advection (ta) term. See notes above
    call xpyp_term_ta_pdf_lhs_all( coef_wpxp2_implicit(:),              & ! Intent(in)
                                   rho_ds_zt(:),                        & ! Intent(in)
                                   invrs_rho_ds_zm(:),                  & ! Intent(in)
                                   gr%invrs_dzm(:),                     & ! Intent(in)
                                   l_upwind_xpyp_ta,                    & ! Intent(in)
                                   sgn_turbulent_vel(:),                & ! Intent(in)
                                   coef_wpxp2_implicit_zm(:),           & ! Intent(in)
                                   rho_ds_zm(:),                        & ! Intent(in)
                                   gr%invrs_dzt(:),                     & ! Intent(in)
                                   lhs_ta(:,:)                        ) ! Intent(out)

    ! Vertical compression of eddies causes gustiness (increase in up2 and vp2)
    ! Add half the contribution to up2 and half to vp2
    rhs(2:gr%nz-1) = rhs_ta(2:gr%nz-1) - 0.5_core_rknd*wp2_splat(2:gr%nz-1)

    ! Finish RHS calc with vectorizable loop, functions are in source file and should
    ! be inlined with an -O2 or above compiler optimization flag
    do k = 2, gr%nz-1, 1

        rhs(k) = rhs(k) + ( one - gamma_over_implicit_ts ) &
                             * ( - lhs_ta(1,k) * xap2(k+1) &
                                 - lhs_ta(2,k) * xap2(k) &
                                 - lhs_ta(3,k) * xap2(k-1) )

        ! RHS turbulent production (tp) term.
        ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_pr 
        rhs(k) = rhs(k) + ( one - C5 ) * term_tp( xam(k+1), xam(k), xam(k+1), xam(k), & 
                                                wpxap(k), wpxap(k), gr%invrs_dzm(k) )

        ! RHS pressure term 1 (pr1) (and dissipation term 1 (dp1)).
        rhs(k) = rhs(k) + term_pr1( C4, C14, xbp2(k), wp2(k), tau_zm(k) )

        ! RHS contribution from "over-implicit" weighted time step
        ! for LHS dissipation term 1 (dp1) and pressure term 1 (pr1).
        rhs(k) = rhs(k) + ( one - gamma_over_implicit_ts ) &
                        * ( - term_dp1_lhs( C4_C14_1d(k), tau_zm(k) ) * xap2(k) )

        ! RHS pressure term 2 (pr2).
        rhs(k) = rhs(k) + term_pr2( C5, thv_ds_zm(k), wpthvp(k), wpxap(k), wpxbp(k), &
                                  xam, xbm, gr%invrs_dzm(k), k+1, k, &
                                  Lscale(k+1), Lscale(k), wp2_zt(k+1), wp2_zt(k) )
    enddo ! k=2..gr%nz-1


    ! RHS time tendency.
    if ( l_iter ) then
        do k = 2, gr%nz-1
            rhs(k) = rhs(k) + one/dt * xap2(k)
        end do
    endif


    if ( l_stats_samp ) then

        ! Statistics: explicit contributions for up2 or vp2.

        ! x'y' term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_ta_ADG1_rhs.

        do k = 2, gr%nz-1

            call stat_begin_update_pt( ixapxbp_ta, k, &    ! Intent(in)
                                       -rhs_ta(k), &
                                       stats_zm )          ! Intent(inout)

            call stat_modify_pt( ixapxbp_ta, k,  &          ! Intent(in)
                                 + ( one - gamma_over_implicit_ts )  & ! Intent(in)
                                   * ( - lhs_ta(1,k) * xap2(k+1) &
                                       - lhs_ta(2,k) * xap2(k) &
                                       - lhs_ta(3,k) * xap2(k-1) ), &
                                 stats_zm )                 ! Intent(inout)

            if ( ixapxbp_dp1 > 0 ) then

              tmp  &
              = gamma_over_implicit_ts  &
              * term_dp1_lhs( two_thirds*C4, tau_zm(k) )
              zmscr01(k) = -tmp
              call stat_begin_update_pt( ixapxbp_dp1, k, & ! Intent(in)
                   -term_pr1( C4, zero, xbp2(k), wp2(k), tau_zm(k) ), & ! Intent(in)
                                         stats_zm )        ! Intent(inout)

              lhs_ta(1,k)  &
              = term_dp1_lhs( two_thirds*C4, tau_zm(k) )
              call stat_modify_pt( ixapxbp_dp1, k, &        ! Intent(in)
                    + ( one - gamma_over_implicit_ts )  &   ! Intent(in)
                    * ( - lhs_ta(1,k) * xap2(k) ),  & ! Intent(in)
                                         stats_zm )         ! Intent(inout)

            endif

            if ( ixapxbp_pr1 > 0 ) then
              tmp  &
              = gamma_over_implicit_ts  &
              * term_dp1_lhs( one_third*C14, tau_zm(k) )
              zmscr11(k) = -tmp
              call stat_begin_update_pt( ixapxbp_pr1, k, & ! Intent(in)  
                   -term_pr1( zero, C14, xbp2(k), wp2(k), tau_zm(k) ), &! Intent(in)
                                         stats_zm )        ! Intent(inout)

              lhs_ta(1,k)  &
              = term_dp1_lhs( one_third*C14, tau_zm(k) )
              call stat_modify_pt( ixapxbp_pr1, k, &        ! Intent(in)
                    + ( one - gamma_over_implicit_ts )  &   ! Intent(in)
                    * ( - lhs_ta(1,k) * xap2(k) ),  & ! Intent(in)
                                         stats_zm )         ! Intent(inout)

            endif

            ! x'y' term pr2 is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( ixapxbp_pr2, k, & ! Intent(in)
                   term_pr2( C5, thv_ds_zm(k), wpthvp(k), wpxap(k), wpxbp(k), & ! In
                             xam, xbm, gr%invrs_dzm(k), k+1, k, &  
                             Lscale(k+1), Lscale(k), wp2_zt(k+1), wp2_zt(k) ), &
                   stats_zm )                          ! Intent(inout)

            ! x'y' term tp is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( ixapxbp_tp, k, & ! Intent(in) 
                  ( one - C5 ) &                      ! Intent(in)
                   * term_tp( xam(k+1), xam(k), xam(k+1), xam(k), &
                              wpxap(k), wpxap(k), gr%invrs_dzm(k) ), & 
                                     stats_zm )       ! Intent(inout)

            ! Vertical compression of eddies.
            call stat_update_var_pt( ixapxbp_splat, k, & ! Intent(in) 
                         -0.5_core_rknd * wp2_splat(k),  & ! Intent(in)
                                     stats_zm )       ! Intent(inout)

        end do

      endif ! l_stats_samp


    ! Boundary Conditions
    ! These are set so that the surface_varnce value of u'^2 or v'^2 can be
    ! used at the lowest boundary and the values of those variables can be
    ! set to their respective threshold minimum values at the top boundary.
    ! Fixed-point boundary conditions are used for the variances.

    rhs(1) = xap2(1)
    ! The value of u'^2 or v'^2 at the upper boundary will be set to the
    ! threshold minimum value of w_tol_sqd.
    rhs(gr%nz) = w_tol_sqd

    return
  end subroutine xp2_xpyp_uv_rhs

  !=============================================================================
  subroutine xp2_xpyp_rhs( solve_type, dt, l_iter, & ! In
                           coef_wpxpyp_implicit, coef_wpxpyp_implicit_zm, & ! In
                           term_wpxpyp_explicit, term_wpxpyp_explicit_zm, & ! In
                           sgn_turbulent_vel, wpxap, wpxbp, & ! In
                           xam, xbm, xapxbp, xpyp_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                           Cn, tau_zm, threshold, & ! In
                           rhs ) ! Out

  ! Description:
  !   Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', sclr'r_t',
  !   sclr'th_l', or sclr'^2.
  ! 
  !   Note:  
  !     For LHS turbulent advection term:
  !        An "over-implicit" weighted time step is applied to this term.
  !        The weight of the implicit portion of this term is controlled by
  !        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  !        expression below).  A factor is added to the right-hand side of
  !        the equation in order to balance a weight that is not equal to 1,
  !        such that:
  !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  !        where X is the variable that is being solved for in a predictive
  !        equation (x'^2 or x'y' in this case), y(t) is the linearized
  !        portion of the term that gets treated implicitly, and RHS is the
  !        portion of the term that is always treated explicitly.  A weight
  !        of greater than 1 can be applied to make the term more
  !        numerically stable.
  ! 
  !   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  !   Significant changes to this routine may adversely affect computational speed
  !       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  !----------------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only: &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use turbulent_adv_pdf, only: &
        xpyp_term_ta_pdf_lhs_all, & ! Procedure(s)
        xpyp_term_ta_pdf_lhs,     &
        xpyp_term_ta_pdf_rhs_all, &
        xpyp_term_ta_pdf_rhs

    use model_flags, only: &
        l_upwind_xpyp_ta    ! Constant(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use stats_type_utilities, only: & 
        stat_begin_update_pt, & ! Procedure(s)
        stat_update_var_pt, &
        stat_modify_pt

    use stats_variables, only: & 
        irtp2_ta,      & ! Variable(s)
        irtp2_tp,      & 
        irtp2_dp1,     &
        irtp2_forcing, &
        ithlp2_ta,      & 
        ithlp2_tp,      &
        ithlp2_dp1,     &
        ithlp2_forcing, & 
        irtpthlp_ta,      & 
        irtpthlp_tp1,     & 
        irtpthlp_tp2,     &
        irtpthlp_dp1,     &
        irtpthlp_forcing, & 
        stats_zm, & 
        l_stats_samp
  
    use advance_helper_module, only: &
        set_boundary_conditions_rhs    ! Procedure(s)

    implicit none

    !------------------- Input Variables -------------------
    integer, intent(in) :: solve_type

    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                              [s]

    logical, intent(in) :: & 
      l_iter   ! Whether x is prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      coef_wpxpyp_implicit,    & ! Coef. of <x'y'> in <w'x'y'> eq.; t-levs [m/s]
      coef_wpxpyp_implicit_zm, & ! coef_wpxpyp_implicit interp. to m-levs  [m/s]
      term_wpxpyp_explicit,    & ! RHS term: <w'x'y'> eq.; zt  [m/s{x un}{y un}]
      term_wpxpyp_explicit_zm, & ! term_wpxpyp_expl interp. zm [m/s{x un}{y un}]
      sgn_turbulent_vel,       & ! Sign of turbulent velocity ("upwind" ta)  [-]
      wpxap,                   & ! w'x_a' (momentum levels)     [m/s{x_a units}]
      wpxbp,                   & ! w'x_b' (momentum levels)     [m/s{x_b units}]
      xam,                     & ! x_am (thermodynamic levels)     [{x_a units}]
      xbm,                     & ! x_bm (thermodynamic levels)     [{x_b units}]
      xapxbp,                  & ! x_a'x_b' (m-levs)          [{x_a un}{x_b un}]
      xpyp_forcing,            & ! <x'y'> forcing (m-levs)      [{x un}{x un}/s]
      rho_ds_zm,               & ! Dry, static density on momentum levs [kg/m^3]
      rho_ds_zt,               & ! Dry, static density on thermo. levs. [kg/m^3]
      invrs_rho_ds_zm,         & ! Inv. dry, static density on m-levs.  [m^3/kg]
      tau_zm,                  & ! Time-scale tau on momentum levels         [s]
      Cn                         ! Coefficient C_n                           [-]

    real( kind = core_rknd ), intent(in) :: &
      threshold    ! Smallest allowable mag. value for x_a'x_b'  [{x_am units}
                   !                                              *{x_bm units}]

    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      rhs     ! Explicit contributions to x variance/covariance terms

    !---------------- Local Variables -------------------
    real( kind = core_rknd ) :: &
      turbulent_prod, & ! Turbulent production term  [{x_a units}*{x_b units}/s]
      xp2_mc_limiter    ! Largest allowable (negative) mc effect

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rhs_ta

    logical :: &
      l_clip_large_neg_mc = .false.  ! Flag to clip excessively large mc values.

    ! Fractional amount of xp2 that microphysics (mc) is allowed to deplete,
    ! where 0 <= mc_xp2_deplete_frac <= 1.  When mc_xp2_deplete_frac has a value
    ! of 0, microphysics is not allowed to deplete xp2 at all, and when
    ! mc_xp2_deplete_frac has a value of 1, microphysics is allowed to deplete
    ! xp2 down to the threshold value of x_tol^2.  This is used only when the
    ! l_clip_large_neg_mc flag is enabled.
    real( kind = core_rknd ), parameter :: &
      mc_xp2_deplete_frac = 0.75_core_rknd

    ! Array indices
    integer :: k, k_low, k_high

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3,gr%nz) :: lhs_ta

    integer :: & 
      ixapxbp_ta, & 
      ixapxbp_tp, & 
      ixapxbp_tp1, & 
      ixapxbp_tp2, &
      ixapxbp_dp1, &
      ixapxbp_f

    !------------------------------ Begin Code ---------------------------------

    select case ( solve_type )
    case ( xp2_xpyp_rtp2 )
      ixapxbp_ta  = irtp2_ta
      ixapxbp_tp  = irtp2_tp
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = irtp2_dp1
      ixapxbp_f   = irtp2_forcing
    case ( xp2_xpyp_thlp2 )
      ixapxbp_ta  = ithlp2_ta
      ixapxbp_tp  = ithlp2_tp
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = ithlp2_dp1
      ixapxbp_f   = ithlp2_forcing
    case ( xp2_xpyp_rtpthlp )
      ixapxbp_ta  = irtpthlp_ta
      ixapxbp_tp  = 0
      ixapxbp_tp1 = irtpthlp_tp1
      ixapxbp_tp2 = irtpthlp_tp2
      ixapxbp_dp1 = irtpthlp_dp1
      ixapxbp_f   = irtpthlp_forcing
    case default ! No budgets for passive scalars
      ixapxbp_ta  = 0
      ixapxbp_tp  = 0
      ixapxbp_tp1 = 0
      ixapxbp_tp2 = 0
      ixapxbp_dp1 = 0
      ixapxbp_f   = 0
    end select

    ! Calculate RHS turbulent advection (ta) term
    call xpyp_term_ta_pdf_rhs_all( term_wpxpyp_explicit(:),           & ! Intent(in)
                                   rho_ds_zt(:),                      & ! Intent(in)
                                   invrs_rho_ds_zm(:),                & ! Intent(in)
                                   gr%invrs_dzm(:),                   & ! Intent(in)
                                   l_upwind_xpyp_ta,                  & ! Intent(in)
                                   sgn_turbulent_vel(:),              & ! Intent(in)
                                   term_wpxpyp_explicit_zm(:),        & ! Intent(in)
                                   rho_ds_zm(:),                      & ! Intent(in)
                                   gr%invrs_dzt(:),                   & ! Intent(in)
                                   rhs_ta(:)                        ) ! Intent(out
    
    ! RHS contribution from "over-implicit" weighted time step
    ! for LHS turbulent advection (ta) term. See notes above
    call xpyp_term_ta_pdf_lhs_all( coef_wpxpyp_implicit(:),             & ! Intent(in)
                                   rho_ds_zt(:),                        & ! Intent(in)
                                   invrs_rho_ds_zm(:),                  & ! Intent(in)
                                   gr%invrs_dzm(:),                     & ! Intent(in)
                                   l_upwind_xpyp_ta,                    & ! Intent(in)
                                   sgn_turbulent_vel(:),                & ! Intent(in)
                                   coef_wpxpyp_implicit_zm(:),          & ! Intent(in)
                                   rho_ds_zm(:),                        & ! Intent(in)
                                   gr%invrs_dzt(:),                     & ! Intent(in)
                                   lhs_ta(:,:)                        ) ! Intent(out)

    ! Finish RHS calc with vectorizable loop, functions are in source file and should
    ! be inlined with an -O2 or above compiler optimization flag
    do k = 2, gr%nz-1

      rhs(k) = rhs_ta(k) + ( one - gamma_over_implicit_ts ) &
                           * ( - lhs_ta(1,k) * xapxbp(k+1) &
                               - lhs_ta(2,k) * xapxbp(k) &
                               - lhs_ta(3,k) * xapxbp(k-1) )

      ! RHS turbulent production (tp) term.
      rhs(k) = rhs(k) + term_tp( xam(k+1), xam(k), xbm(k+1), xbm(k), &
                                 wpxbp(k), wpxap(k), gr%invrs_dzm(k) )

      ! RHS dissipation term 1 (dp1)
      rhs(k) = rhs(k) + term_dp1_rhs( Cn(k), tau_zm(k), threshold )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS dissipation term 1 (dp1).
      rhs(k) = rhs(k)  + ( one - gamma_over_implicit_ts ) &
                       * ( - term_dp1_lhs( Cn(k), tau_zm(k) ) * xapxbp(k) )
    end do

    ! RHS <x'y'> forcing.
    ! Note: <x'y'> forcing includes the effects of microphysics on <x'y'>.
    if ( l_clip_large_neg_mc &
         .and. ( solve_type == xp2_xpyp_rtp2 .or. solve_type == xp2_xpyp_thlp2 ) ) then

        do k = 2, gr%nz-1

            turbulent_prod = term_tp( xam(k+1), xam(k), xbm(k+1), xbm(k), &
                                      wpxbp(k), wpxap(k), gr%invrs_dzm(k) )

            ! Limit the variance-depleting effects of excessively large
            ! microphysics terms on rtp2 and thlp2 in order to reduce oscillations.
            if ( turbulent_prod >= zero ) then

                ! Microphysics is allowed to deplete turbulent production and a
                ! fraction of the variance, determined by mc_xp2_deplete_frac, down
                ! to the threshold.
                xp2_mc_limiter = - mc_xp2_deplete_frac &
                               * ( xapxbp(k) - threshold ) / dt &
                               - turbulent_prod
            else
            
                ! Microphysics is allowed to deplete a fraction of the variance,
                ! determined by mc_xp2_deplete_frac, down to the threshold.
                xp2_mc_limiter = - mc_xp2_deplete_frac &
                               * ( xapxbp(k) - threshold ) / dt
            endif

            ! Note: the value of xp2_mc_limiter is always negative.
            rhs(k) = rhs(k) + max( xpyp_forcing(k), xp2_mc_limiter )

        end do

    else

        do k = 2, gr%nz-1
            rhs(k) = rhs(k) + xpyp_forcing(k)
        end do

    endif 

    
    ! RHS time tendency.
    if ( l_iter ) then
        do k = 2, gr%nz-1, 1
            rhs(k) = rhs(k) + one/dt * xapxbp(k)
        end do
    endif


    if ( l_stats_samp ) then

        do k = 2, gr%nz-1

            ! Statistics: explicit contributions for rtp2, thlp2, or rtpthlp.

            ! <x'y'> term ta has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on term_ta_explicit_rhs.
            call stat_begin_update_pt( ixapxbp_ta, k, & ! Intent(in)
                                       -rhs_ta(k), &
                                       stats_zm )       ! Intent(inout)

            call stat_modify_pt( ixapxbp_ta, k, &             ! Intent(in)
                                 + ( one - gamma_over_implicit_ts ) & ! Intent(in)
                                   * ( - lhs_ta(1,k) * xapxbp(k+1) &
                                       - lhs_ta(2,k) * xapxbp(k) &
                                       - lhs_ta(3,k) * xapxbp(k-1) ), &
                                 stats_zm )                   ! Intent(inout)

            ! x'y' term dp1 has both implicit and explicit components; call
            ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
            ! subtracts the value sent in, reverse the sign on term_dp1_rhs.
            call stat_begin_update_pt( ixapxbp_dp1, k, &           ! Intent(in)
                 -term_dp1_rhs( Cn(k), tau_zm(k), threshold ), &   ! Intent(in)
                                       stats_zm )                  ! Intent(inout)

            ! Note:  An "over-implicit" weighted time step is applied to this term.
            !        A weighting factor of greater than 1 may be used to make the
            !        term more numerically stable (see note above for RHS turbulent
            !        advection (ta) term).
            lhs_ta(1,k)  &
            = term_dp1_lhs( Cn(k), tau_zm(k) )
            call stat_modify_pt( ixapxbp_dp1, k,  &         ! Intent(in)
                  + ( one - gamma_over_implicit_ts )  &     ! Intent(in)
                  * ( - lhs_ta(1,k) * xapxbp(k) ),  & ! Intent(in)
                                       stats_zm )                 ! Intent(inout)

            ! rtp2/thlp2 case (1 turbulent production term)
            ! x'y' term tp is completely explicit; call stat_update_var_pt.
            call stat_update_var_pt( ixapxbp_tp, k, &             ! Intent(in)
                  term_tp( xam(k+1), xam(k), xbm(k+1), xbm(k), &  ! Intent(in)
                           wpxbp(k), wpxap(k), gr%invrs_dzm(k) ), & 
                                     stats_zm )                         ! Intent(inout)

            ! rtpthlp case (2 turbulent production terms)
            ! x'y' term tp1 is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of x'y' term tp1, substitute 0 for all
            !        the xam inputs and the wpxbp input to function term_tp.
            call stat_update_var_pt( ixapxbp_tp1, k, &    ! Intent(in)
                  term_tp( zero, zero, xbm(k+1), xbm(k), &  ! Intent(in)
                           zero, wpxap(k), gr%invrs_dzm(k) ), &
                                     stats_zm )                 ! Intent(inout)

            ! x'y' term tp2 is completely explicit; call stat_update_var_pt.
            ! Note:  To find the contribution of x'y' term tp2, substitute 0 for all
            !        the xbm inputs and the wpxap input to function term_tp.
            call stat_update_var_pt( ixapxbp_tp2, k, &    ! Intent(in)
                  term_tp( xam(k+1), xam(k), zero, zero, &  ! Intent(in)
                           wpxbp(k), zero, gr%invrs_dzm(k) ), &
                                     stats_zm )                 ! Intent(inout)

            ! x'y' forcing term is completely explicit; call stat_update_var_pt.
            if ( l_clip_large_neg_mc &
                 .and. ( solve_type == xp2_xpyp_rtp2 &
                         .or. solve_type == xp2_xpyp_thlp2 ) ) then
               call stat_update_var_pt( ixapxbp_f, k, &
                                        max( xpyp_forcing(k), xp2_mc_limiter ), &
                                        stats_zm )
            else
               call stat_update_var_pt( ixapxbp_f, k, xpyp_forcing(k), stats_zm )
            endif
   
        enddo ! k=2..gr%nz-1

    endif ! l_stats_samp


    ! Boundary Conditions
    ! These are set so that the surface_varnce value of rtp2, thlp2, or rtpthlp
    ! (or sclrp2, sclrprtp, or sclrpthlp) can be used at the lowest boundary and the
    ! values of those variables can be set to their respective threshold minimum
    ! values (which is 0 in the case of the covariances) at the top boundary.
    ! Fixed-point boundary conditions are used for both the variances and the
    ! covariances.

    k_low = 1
    k_high = gr%nz

    ! The value of the field at the upper boundary will be set to it's threshold
    ! minimum value, as contained in the variable 'threshold'.
    call set_boundary_conditions_rhs( xapxbp(1), k_low, threshold, k_high, rhs(:) )

    return
  end subroutine xp2_xpyp_rhs

  !=============================================================================
  pure function term_tp( xamp1, xam, xbmp1, xbm,  & 
                         wpxbp, wpxap, invrs_dzm ) & 
  result( rhs )

    ! Description:
    ! Turbulent production of x_a'x_b':  explicit portion of the code.
    !
    ! The d(x_a'x_b')/dt equation contains a turbulent production term:
    !
    ! - w'x_b' d(x_am)/dz - w'x_a' d(x_bm)/dz.
    !
    ! This term is solved for completely explicitly and is discretized as
    ! follows:
    !
    ! The values of w'x_a' and w'x_b' are found on the momentum levels, whereas
    ! the values of x_am and x_bm are found on the thermodynamic levels.  The
    ! derivatives of both x_am and x_bm are taken over the intermediate
    ! (central) momentum level.  All of the remaining mathematical operations
    ! take place at the central momentum level, yielding the desired result.
    !
    ! ---------xamp1------------xbmp1-------------------------- t(k+1)
    !
    ! ===wpxap======d(xam)/dz=========d(xbm)/dz===wpxbp======== m(k)
    !
    ! ---------xam--------------xbm---------------------------- t(k)
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

    ! Input variables
    real( kind = core_rknd ), intent(in) :: & 
      xam,       & ! x_am(k)                     [{x_am units}]
      xamp1,     & ! x_am(k+1)                   [{x_am units}]
      xbm,       & ! x_bm(k)                     [{x_bm units}]
      xbmp1,     & ! x_bm(k+1)                   [{x_bm units}]
      wpxbp,     & ! w'x_b'(k)                   [m/s {x_bm units}]
      wpxap,     & ! w'x_a'(k)                   [m/s {x_am units}]
      invrs_dzm    ! Inverse of grid spacing (k) [1/m]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs & 
    = - wpxbp * invrs_dzm * ( xamp1 - xam ) & 
      - wpxap * invrs_dzm * ( xbmp1 - xbm )

    return
  end function term_tp

  !=============================================================================
  pure function term_dp1_lhs( Cn, tau_zm )  & 
  result( lhs )

    ! Description:
    ! Dissipation term 1 for x_a'x_b':  implicit portion of the code.
    !
    ! The d(x_a'x_b')/dt equation contains dissipation term 1:
    !
    ! - ( C_n / tau_zm ) x_a'x_b'.
    !
    ! For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
    ! are the same variable), the term is damped to a certain positive
    ! threshold, such that:
    !
    ! - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
    !
    ! However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
    ! is part of pressure term 1 and is handled as part of function 'term_pr1'.
    ! Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
    ! 'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.
    !
    ! For cases where x_a'x_b' is a covariance (in other words, where x_a and
    ! x_b are different variables), threshold is set to 0, and the expression
    ! reverts to the form found in the first equation.
    !
    ! This term is broken into implicit and explicit portions.  The equations
    ! for u'^2, v'^2, and any covariances only include the implicit portion.
    ! The implicit portion of this term is:
    !
    ! - ( C_n / tau_zm ) x_a'x_b'(t+1).
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the sign is reversed and the leading "-" in front of the term
    !        is changed to a "+".
    !
    ! The timestep index (t+1) means that the value of x_a'x_b' being used is
    ! from the next timestep, which is being advanced to in solving the
    ! d(x_a'x_b')/dt equation.
    !
    ! The values of x_a'x_b' are found on momentum levels.  The values of
    ! time-scale tau_zm are also found on momentum levels.
    !
    ! Note:  For equations that use pressure term 1 (such as the equations for
    !        u'^2 and v'^2), C_n = ( 2*C_4 + C_14 ) / 3; which combines the
    !        implicit contributions for dissipation term 1 and pressure term 1
    !        into one expression.  Otherwise, C_n = C_2.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      Cn,    & ! Coefficient C_n                       [-]
      tau_zm   ! Time-scale tau at momentum levels (k) [s]

    ! Return Variable
    real( kind = core_rknd ) :: lhs

    ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
    lhs  & 
    = + Cn / tau_zm

    return
  end function term_dp1_lhs

  !=============================================================================
  pure function term_dp1_rhs( Cn, tau_zm, threshold ) &
  result( rhs )

    ! Description:
    ! Dissipation term 1 for x_a'x_b':  explicit portion of the code.
    !
    ! The d(x_a'x_b')/dt equation contains dissipation term 1:
    !
    ! - ( C_n / tau_zm ) x_a'x_b'.
    !
    ! For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
    ! are the same variable), the term is damped to a certain positive
    ! threshold, such that:
    !
    ! - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
    !
    ! However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
    ! is part of pressure term 1 and is handled as part of function 'term_pr1'.
    ! Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
    ! 'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.
    !
    ! For cases where x_a'x_b' is a covariance (in other words, where x_a and
    ! x_b are different variables), threshold is set to 0, and the expression
    ! reverts to the form found in the first equation.
    !
    ! This term is broken into implicit and explicit portions.  The equations
    ! for u'^2, v'^2, and any covariances only include the implicit portion.
    ! The explicit portion of this term is:
    !
    ! + ( C_n / tau_zm ) * threshold.
    !
    ! The values of time-scale tau_zm and the threshold are found on the
    ! momentum levels.
    !
    ! Note:  The equations that use pressure term 1 (such as the equations for
    !        u'^2 and v'^2) do not call this function.  Thus, within this
    !        function, C_n = C_2.

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Cn,       & ! Coefficient C_n                               [-]
      tau_zm,   & ! Time-scale tau at momentum levels (k)         [s]
      threshold   ! Minimum allowable magnitude value of x_a'x_b' [units vary]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs  & 
    = + ( Cn / tau_zm ) * threshold

    return
  end function term_dp1_rhs

  !=============================================================================
  pure function term_pr1( C4, C14, xbp2, wp2, tau_zm ) & 
  result( rhs )

    ! Description:
    ! Pressure term 1 for x_a'x_b':  explicit portion of the code.
    !
    ! Note:  Pressure term 1 is only used when x_a'x_b' is either u'^2 or v'^2.
    !        For the following description, pressure term 2 for u'^2 is used as
    !        the example.  Pressure term 2 for v'^2 is the same as pressure
    !        term 2 for u'^2, except that the v'^2 and u'^2 variables are
    !        switched.
    !
    ! The d(u'^2)/dt equation contains dissipation term 1:
    !
    ! - ( C_4 / tau_zm ) * ( u'^2 - (2/3)*em );
    !
    ! where em = (1/2) * ( u'^2 + v'^2 + w'^2 );
    !
    ! and with the substitution applied, dissipation term 1 becomes:
    !
    ! - ( C_4 / tau_zm ) * ( u'^2 - (1/3) * ( u'^2 + v'^2 + w'^2 ) ).
    !
    ! The d(u'^2)/dt equation also contains pressure term 1:
    !
    ! - (2/3) * epsilon;
    !
    ! where epsilon = C_14 * ( em / tau_zm ).
    !
    ! Additionally, since pressure term 1 is a damping term, em is damped only
    ! to it's minimum threshold value, em_min, where:
    !
    ! em_min = (1/2) * ( u'^2|_min + v'^2|_min + w'^2|_min )
    !      = (1/2) * ( w_tol^2 + w_tol^2 + w_tol^2 )
    !      = (3/2) * w_tol^2.
    !
    ! With the damping threshold applied, epsilon becomes:
    !
    ! epsilon = C_14 * ( ( em - em_min ) / tau_zm );
    !
    ! and with all substitutions applied, pressure term 1 becomes:
    !
    ! - (2/3) * ( C_14 / tau_zm )
    !         * [ (1/2) * ( u'^2 + v'^2 + w'^2 ) - (3/2) * w_tol^2 ].
    !
    ! Dissipation term 1 and pressure term 1 are combined and simplify to:
    !
    ! - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2
    !    + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2 + w'^2 )
    !    + ( C_14 / tau_zm ) * w_tol^2.
    !
    ! The combined term has both implicit and explicit components.
    ! The implicit component is:
    !
    ! - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2(t+1).
    !
    ! Note:  When the implicit term is brought over to the left-hand side,
    !        the sign is reversed and the leading "-" in front of the term
    !        is changed to a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(x_a'x_b')/dt equation.
    !
    ! The implicit component of the combined dp1 and pr1 term is solved in
    ! function "term_dp1_lhs" above, where "( 2*C_4 + C_14 ) / 3" is sent in
    ! as "C_n".
    !
    ! The explicit component of the combined dp1 and pr1 term is:
    !
    ! + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2(t) + w'^2(t) )
    ! + ( C_14 / tau_zm ) * w_tol^2;
    !
    ! and is discretized as follows:
    !
    ! The values for v'^2 and w'^2, as well as for tau_zm, are found on the
    ! momentum levels.  The mathematical operations all take place on the
    ! momentum levels, yielding the desired result.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        w_tol_sqd, & ! Constant(s)
        one_third

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C4,    & ! Model parameter C_4                         [-]
      C14,   & ! Model parameter C_14                        [-]
      xbp2,  & ! v'^2(k) (if solving for u'^2) or vice versa [m^2/s^2]
      wp2,   & ! w'^2(k)                                     [m^2/s^2]
      tau_zm   ! Time-scale tau at momentum levels (k)       [s]

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    rhs = + one_third * ( C4 - C14 ) * ( xbp2 + wp2 ) / tau_zm  &
          + ( C14 / tau_zm ) * w_tol_sqd

    return
  end function term_pr1

  !=============================================================================
  function term_pr2( C5, thv_ds_zm, wpthvp, upwp, vpwp, &
                          um, vm, invrs_dzm, kp1, k, & 
                          Lscalep1, Lscale, wp2_ztp1, wp2_zt ) &
  result( rhs )

    ! Description:
    ! Pressure term 2 for x_a'x_b':  explicit portion of the code.
    !
    ! Note:  Pressure term 2 is only used when x_a'x_b' is either u'^2 or v'^2.
    !        For the following description, pressure term 2 for u'^2 is used as
    !        the example.  Pressure term 2 for v'^2 is the exact same as
    !        pressure term 2 for u'^2.
    !
    ! The d(u'^2)/dt equation contains pressure term 2:
    !
    ! + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].
    !
    ! This term is solved for completely explicitly and is discretized as
    ! follows:
    !
    ! The values of w'th_v', u'w', and v'w' are found on the momentum levels,
    ! whereas the values of um and vm are found on the thermodynamic levels.
    ! Additionally, the values of thv_ds_zm are found on the momentum levels.
    ! The derivatives of both um and vm are taken over the intermediate
    ! (central) momentum level.  All the remaining mathematical operations take
    ! place at the central momentum level, yielding the desired result.
    !
    ! -----ump1------------vmp1-------------------------------------- t(k+1)
    !
    ! =upwp====d(um)/dz========d(vm)/dz==vpwp===thv_ds_zm==wpthvp==== m(k)
    !
    ! -----um--------------vm---------------------------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: & ! Constants 
        grav, & ! Gravitational acceleration [m/s^2]
        one, &
        two_thirds, &
        zero, &
        zero_threshold

    use grid_class, only: &
        gr ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: abs, max

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      C5,        & ! Model parameter C_5                            [-]
      thv_ds_zm, & ! Dry, base-state theta_v at momentum level (k)  [K]
      wpthvp,    & ! w'th_v'(k)                                     [m/K/s]
      upwp,      & ! u'w'(k)                                        [m^2/s^2]
      vpwp,      & ! v'w'(k)                                        [m^2/s^2]
      invrs_dzm, & ! Inverse of the grid spacing (k)                [1/m]
      Lscalep1,  & ! Mixing length (k+1)                            [m]
      Lscale,    & ! Mixing length (k)                              [m]
      wp2_ztp1,  & ! w'^2(k+1) (thermo. levels)                     [m^2/s^2]
      wp2_zt       ! w'^2(k)   (thermo. levels)                     [m^2/s^2]

    ! Note: Entire arrays of um and vm are now required rather than um and vm
    ! only at levels k and k+1.  The entire array is necessary when a vertical
    ! average calculation of d(um)/dz and d(vm)/dz is used. --ldgrant March 2010
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      um,  & ! mean zonal wind       [m/s]
      vm     ! mean meridional wind  [m/s]

    integer, intent(in) :: &
      kp1, & ! current level+1 in xp2_xpyp_uv_rhs loop
      k      ! current level in xp2_xpyp_uv_rhs loop

    ! Return Variable
    real( kind = core_rknd ) :: rhs

    ! Local Variable(s)  --ldgrant, March 2010
    real( kind = core_rknd ), parameter :: &
      ! Constants empirically determined for experimental version of term_pr2 
      ! ldgrant March 2010
      constant1 = one, &     ! [m/s]
      constant2 = 1000.0_core_rknd, &  ! [m]
      vert_avg_depth = 200.0_core_rknd ! Depth over which to average d(um)/dz and d(vm)/dz [m]

    real( kind = core_rknd ) :: &
      zt_high, & ! altitude above current altitude zt(k)         [m]
      um_high, & ! um at altitude zt_high                        [m/s]
      vm_high, & ! vm at altitude zt_high                        [m/s]
      zt_low,  & ! altitude below (or at) current altitude zt(k) [m]
      um_low,  & ! um at altitude zt_low                         [m/s]
      vm_low     ! vm at altitude zt_low                         [m/s]

    logical, parameter :: & 
      l_use_experimental_term_pr2 = .false., & ! If true, use experimental version
                                               ! of term_pr2 calculation
      l_use_vert_avg_winds = .true. ! If true, use vert_avg_depth average
    ! calculation for d(um)/dz and d(vm)/dz

    !------ Begin code ------------

    if( .not. l_use_experimental_term_pr2 ) then
      ! use original version of term_pr2

      ! As applied to w'2
      rhs = + two_thirds * C5 & 
              * ( ( grav / thv_ds_zm ) * wpthvp &
                  - upwp * invrs_dzm * ( um(kp1) - um(k) ) &
                  - vpwp * invrs_dzm * ( vm(kp1) - vm(k) ) &
                )

    else ! use experimental version of term_pr2 --ldgrant March 2010

      if( l_use_vert_avg_winds ) then
        ! We found that using a 200m running average of d(um)/dz and d(vm)/dz
        ! produces larger spikes in up2 and vp2 near the inversion for
        ! the stratocumulus cases.
        call find_endpts_for_vert_avg_winds &
             ( vert_avg_depth, k, um, vm, & ! intent(in)
               zt_high, um_high, vm_high, & ! intent(out)
               zt_low, um_low, vm_low )     ! intent(out)

      else ! Do not use a vertical average calculation for d(um)/dz and d(vm)/dz
        zt_high = gr%zt(kp1)
        um_high = um(kp1)
        vm_high = vm(kp1)

        zt_low  = gr%zt(k)
        um_low  = um(k)
        vm_low  = vm(k)
      end if ! l_use_vert_avg_winds

      ! *****NOTES on experimental version*****
      ! Leah Grant and Vince Larson eliminated the contribution from wpthvp
      ! because terms with d(wp2)/dz include buoyancy effects and seem to
      ! produce better results.
      !
      ! We also eliminated the contribution from the momentum flux terms
      ! because they didn't contribute to the results.
      !
      ! The constant1 line does not depend on shear.  This is important for
      ! up2 and vp2 generation in cases that have little shear such as FIRE.
      ! We also made the constant1 line proportional to d(Lscale)/dz to account
      ! for higher spikes in up2 and vp2 near a stronger inversion.  This
      ! increases up2 and vp2 near the inversion for the stratocumulus cases,
      ! but overpredicts up2 and vp2 near cloud base in cumulus cases such
      ! as BOMEX where d(Lscale)/dz is large.  Therefore, the d(Lscale)/dz
      ! contribution is commented out for now.
      !
      ! The constant2 line includes the possibility of shear generation of
      ! up2 and vp2, which is important for some cases.  The current functional
      ! form used is:
      !   constant2 * |d(wp2)/dz| * |d(vm)/dz|
      ! We use  |d(vm)/dz|  instead of  |d(um)/dz| + |d(vm)/dz|  here because
      ! this allows for different profiles of up2 and vp2, which occur for
      ! many cases.  In addition, we found that in buoyant cases, up2 is
      ! more related to d(vm)/dz and vp2 is more related to d(um)/dz.  This
      ! occurs if horizontal rolls are oriented in the direction of the shear
      ! vector.  However, in stably stratified cases, the opposite relation is
      ! true (horizontal rolls caused by shear are perpendicular to the shear
      ! vector).  This effect is not yet accounted for.
      !
      ! For better results, we reduced the value of C5 from 5.2 to 3.0 and
      ! changed the eddy diffusivity coefficient Kh so that it is
      ! proportional to 1.5*wp2 rather than to em.
      rhs = + two_thirds * C5 & 
              * ( constant1 * abs( wp2_ztp1 - wp2_zt ) * invrs_dzm &
                    ! * abs( Lscalep1 - Lscale ) * invrs_dzm &
                  + constant2 * abs( wp2_ztp1 - wp2_zt ) * invrs_dzm &
                    * abs( vm_high - vm_low ) / ( zt_high - zt_low ) &
                     + ( Lscalep1 + Lscale ) * zero &    
                             ! This line eliminates an Intel compiler
                )            ! warning that Lscalep1/Lscale are not
                             ! used. -meyern
    end if ! .not. l_use_experimental_term_pr2

    ! Added by dschanen for ticket #36
    ! We have found that when shear generation is zero this term will only be
    ! offset by hole-filling (up2_pd/vp2_pd) and reduces turbulence 
    ! unrealistically at lower altitudes to make up the difference.
    rhs = max( rhs, zero_threshold )

    return
  end function term_pr2

  !=============================================================================
  subroutine find_endpts_for_vert_avg_winds &
                  ( vert_avg_depth, k, um, vm, & ! intent(in)
                    zt_high, um_high, vm_high, & ! intent(out)
                    zt_low, um_low, vm_low )     ! intent(out)
    ! Description:
    ! This subroutine determines values of um and vm which are
    ! +/- [vert_avg_depth/2] m above and below the current altitude zt(k).
    ! This is for the purpose of using a running vertical average
    ! calculation of d(um)/dz and d(vm)/dz in term_pr2 (over a depth
    ! vert_avg_depth).  E.g. If a running average over 200m is desired,
    ! then this subroutine will determine the values of um and vm which
    ! are 100m above and below the current level.
    ! ldgrant March 2010
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two    ! Constant(s)

    use interpolation, only : &
        binary_search, lin_interpolate_two_points  ! Function(s)

    use grid_class, only: &
        gr ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      vert_avg_depth ! Depth over which to average d(um)/dz
    ! and d(vm)/dz in term_pr2 [m]

    integer, intent(in) :: &
      k ! current level in xp2_xpyp_uv_rhs loop

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      um,  & ! mean zonal wind       [m/s]
      vm     ! mean meridional wind  [m/s]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      zt_high, & ! current altitude zt(k) + depth [m]
      um_high, & ! um at altitude zt_high         [m/s]
      vm_high, & ! vm at altitude zt_high         [m/s]
      zt_low,  & ! current altitude zt(k) - depth [m]
      um_low,  & ! um at altitude zt_low          [m/s]
      vm_low     ! vm at altitude zt_low          [m/s]

    ! Local Variables
    real( kind = core_rknd ) :: depth ! vert_avg_depth/2 [m]

    integer :: k_high, k_low
    ! Number of levels above (below) the current level where altitude is
    ! [depth] greater (less) than the current altitude
    ! [unless zt(k) < [depth] from an upper/lower boundary]

    !------ Begin code ------------

    depth = vert_avg_depth / two

    ! Find the grid level that contains the altitude greater than or
    ! equal to the current altitude + depth
    k_high = binary_search( gr%nz, gr%zt, gr%zt(k)+depth )
    ! If the current altitude + depth is greater than the highest
    ! altitude, binary_search returns a value of -1
    if ( k_high == -1 ) k_high = gr%nz

    if ( k_high == gr%nz ) then
      ! Current altitude + depth is higher than or exactly at the top grid level.
      ! Since this is a ghost point, use the altitude at grid level nzmax-1
      k_high = gr%nz-1
      zt_high = gr%zt(k_high)
      um_high = um(k_high)
      vm_high = vm(k_high)
    else ! Do an interpolation to find um & vm at current altitude + depth.
      zt_high = gr%zt(k)+depth
      um_high = lin_interpolate_two_points( zt_high, gr%zt(k_high), gr%zt(k_high-1), &
                         um(k_high), um(k_high-1) )
      vm_high = lin_interpolate_two_points( zt_high, gr%zt(k_high), gr%zt(k_high-1), &
                         vm(k_high), vm(k_high-1) )
    end if ! k_high ...


    ! Find the grid level that contains the altitude less than or
    ! equal to the current altitude - depth
    k_low = binary_search( gr%nz, gr%zt, gr%zt(k)-depth )
    ! If the current altitude - depth is less than the lowest
    ! altitude, binary_search returns a value of -1
    if ( k_low == -1 ) k_low = 2

    if ( k_low == 2 ) then
      ! Current altitude - depth is less than or exactly at grid level 2.
      ! Since grid level 1 is a ghost point, use the altitude at grid level 2
      zt_low = gr%zt(k_low)
      um_low = um(k_low)
      vm_low = vm(k_low)
    else ! Do an interpolation to find um at current altitude - depth.
      zt_low = gr%zt(k)-depth
      um_low = lin_interpolate_two_points( zt_low, gr%zt(k_low), gr%zt(k_low-1), &
                        um(k_low), um(k_low-1) )
      vm_low = lin_interpolate_two_points( zt_low, gr%zt(k_low), gr%zt(k_low-1), &
                        vm(k_low), vm(k_low-1) )
    end if ! k_low ...

    return
  end subroutine find_endpts_for_vert_avg_winds

  !=============================================================================
  subroutine pos_definite_variances( solve_type, dt, tolerance, &
                                     rho_ds_zm, rho_ds_zt, &
                                     xp2_np1 )

    ! Description:
    ! Use the hole filling code to make a variance term positive definite
    !-----------------------------------------------------------------------

    use fill_holes, only: fill_holes_vertical
    use grid_class, only: gr
    use clubb_precision, only: core_rknd

    use stats_variables, only:  & 
        stats_zm, l_stats_samp, & 
        irtp2_pd, ithlp2_pd, iup2_pd, ivp2_pd ! variables
    use stats_type_utilities, only:  & 
        stat_begin_update, stat_end_update ! subroutines


    implicit none

    ! External
    intrinsic :: any, real, trim

    ! Input variables
    integer, intent(in) :: & 
      solve_type

    real( kind = core_rknd ), intent(in) :: & 
      dt        ! Model timestep              [s]

    real( kind = core_rknd ), intent(in) :: & 
      tolerance ! Threshold for xp2_np1       [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zm, & ! Dry, static density on momentum levels         [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermodynamic levels    [kg/m^3]

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      xp2_np1   ! Variance for <n+1>          [units vary]

    ! Local variables
    integer :: & 
      ixp2_pd

    select case( solve_type )
    case ( xp2_xpyp_rtp2 )
      ixp2_pd = irtp2_pd
    case ( xp2_xpyp_thlp2 )
      ixp2_pd = ithlp2_pd
    case ( xp2_xpyp_up2 )
      ixp2_pd = iup2_pd
    case ( xp2_xpyp_vp2 )
      ixp2_pd = ivp2_pd
    case default
      ixp2_pd = 0 ! This includes the passive scalars
    end select

    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_begin_update( ixp2_pd, xp2_np1 / dt, &   ! Intent(in)
                              stats_zm )                               ! Intent(inout)
    endif


    if ( any( xp2_np1 < tolerance ) ) then

      ! Call the hole-filling scheme.
      ! The first pass-through should draw from only two levels on either side
      ! of the hole.
      call fill_holes_vertical( 2, tolerance, "zm",   & ! Intent(in)
                              rho_ds_zt, rho_ds_zm, & ! Intent(in)
                              xp2_np1 )               ! Intent(inout)

    endif

    if ( l_stats_samp ) then
      ! Store previous value for effect of the positive definite scheme
      call stat_end_update( ixp2_pd, xp2_np1 / dt, & ! Intent(in)
                            stats_zm )                             ! Intent(inout)
    endif


    return
  end subroutine pos_definite_variances

  !============================================================================
  subroutine update_xp2_mc( nz, dt, cloud_frac, rcm, rvm, thlm,        &
                            wm, exner, rrm_evap, pdf_params,        &
                            rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                            rtpthlp_mc )
    !Description:
    !This subroutine is for use when l_morr_xp2_mc = .true.
    !The effects of rain evaporation on rtp2 and thlp2 are included by
    !assuming rain falls through the moist (cold) portion of the pdf.
    !This is accomplished by defining a precip_fraction and assuming a double
    !delta shaped pdf, such that the evaporation makes the moist component 
    !moister and the colder component colder. Calculations are done using
    !variables on the zt grid, and the outputs are on the zm grid --storer

    use pdf_parameter_module, only: pdf_parameter

    use grid_class, only: &
      zt2zm   ! Procedure(s)

    use constants_clubb, only: &
      cloud_frac_min, &  !Variables
      Cp, &
      Lv

    use clubb_precision, only: &
      core_rknd ! Variable(s)
      

    implicit none

    !input parameters
    integer, intent(in) :: nz ! Points in the Vertical        [-]

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac, &       !Cloud fraction                        [-]
      rcm, &              !Cloud water mixing ratio              [kg/kg]
      rvm, &              !Vapor water mixing ratio              [kg/kg]
      thlm, &             !Liquid potential temperature          [K]
      wm, &               !Mean vertical velocity                [m/s]
      exner, &            !Exner function                        [-]
      rrm_evap         !Evaporation of rain                   [kg/kg/s]
                          !It is expected that this variable is negative, as
                          !that is the convention in Morrison microphysics

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters

    !input/output variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]

    !local variables
    real( kind = core_rknd ), dimension(nz) :: &
      temp_rtp2, &          !Used only to calculate rtp2_mc  [(kg/kg)^2]
      temp_thlp2, &         !Used to calculate thlp2_mc      [K^2/s]
      temp_wp2,  &          !Used to calculate wpxp_mc       [m^2/s^2]
      rtp2_mc_zt, &   !Calculated on the zt grid             [(kg/kg)^2/s]
      thlp2_mc_zt, &  !Calculated on the zt grid             [(kg/kg)^2/s]
      wprtp_mc_zt, &  !Calculated on the zt grid             [m*(kg/kg)/s^2]
      wpthlp_mc_zt, & !Calcualted on the zt grid             [m*K/s^2]
      rtpthlp_mc_zt,& !Calculated on the zt grid             [K*(kg/kg)/s]
      precip_frac_double_delta, &!Precipitation fraction for a double delta [-]
      pf_const              ! ( 1 - pf )/( pf )                    [-]

    integer :: k

    ! ---- Begin Code ----

    ! Calculate precip_frac_double_delta
    precip_frac_double_delta(nz) = 0.0_core_rknd
    do k = nz-1, 1, -1
      if ( cloud_frac(k) > cloud_frac_min ) then
        precip_frac_double_delta(k) = cloud_frac(k)
      else
        precip_frac_double_delta(k) = precip_frac_double_delta(k+1)
      end if
    end do


    !pf_const is calculated so that when precip_frac_double_delta = 0, rtp2_mc and 
    !thlp2_mc will both be zero.  This also avoids a divide by zero error
    where ( precip_frac_double_delta > cloud_frac_min )
      pf_const = ( 1.0_core_rknd - precip_frac_double_delta ) / precip_frac_double_delta
    elsewhere
      pf_const = 0.0_core_rknd
    end where

    ! Include effects of rain evaporation on rtp2
    temp_rtp2 = pdf_params%mixt_frac &
                    * ( ( pdf_params%rt_1 - ( rcm + rvm ) )**2 + pdf_params%varnce_rt_1 ) &
                + ( 1.0_core_rknd - pdf_params%mixt_frac ) &
                    * ( ( pdf_params%rt_2 - ( rcm + rvm ) )**2 + pdf_params%varnce_rt_2 )

    rtp2_mc_zt = rrm_evap**2 * pf_const * dt &
                       + 2.0_core_rknd * abs(rrm_evap) * sqrt(temp_rtp2 * pf_const)
                       !use absolute value of evaporation, as evaporation will add
                       !to rt_1
    rtp2_mc = zt2zm( rtp2_mc_zt )

    !Include the effects of rain evaporation on thlp2
    temp_thlp2 = pdf_params%mixt_frac &
                    * ( ( pdf_params%thl_1 - thlm )**2 + pdf_params%varnce_thl_1 ) &
                 + ( 1.0_core_rknd - pdf_params%mixt_frac ) &
                    * ( ( pdf_params%thl_2 - thlm )**2 + pdf_params%varnce_thl_2 )

    thlp2_mc_zt = ( rrm_evap * Lv / ( Cp * exner) )**2 &
                       * pf_const * dt &
                       + 2.0_core_rknd * abs(rrm_evap) * Lv / ( Cp * exner ) &
                       * sqrt(temp_thlp2 * pf_const)
    
    thlp2_mc = zt2zm( thlp2_mc_zt )

    ! Include effects of rain evaporation on other moments (wprtp, wpthlp, and 
    ! rtpthlp - added 07/13 rstorer

    temp_wp2 = pdf_params%mixt_frac &
                  * ( ( pdf_params%w_1 - wm )**2 + pdf_params%varnce_w_1 ) &
               + ( 1.0_core_rknd - pdf_params%mixt_frac ) &
                  * ( ( pdf_params%w_2 - wm )**2 + pdf_params%varnce_w_2 )

    wprtp_mc_zt = abs(rrm_evap) * sqrt(pf_const) * sqrt(temp_wp2)

    wpthlp_mc_zt = -1.0_core_rknd * Lv / ( Cp * exner) * abs(rrm_evap) &
                                * sqrt(pf_const) * sqrt(temp_wp2)

    rtpthlp_mc_zt = -1.0_core_rknd * abs(rrm_evap) * sqrt( pf_const ) &
                              * ( ( Lv / (cp * exner ) ) * sqrt( temp_rtp2 ) &
                                + sqrt( temp_thlp2 ) ) &
                            - ( Lv / (cp * exner ) ) * pf_const &
                                * ( rrm_evap )**2 * dt

    wprtp_mc = zt2zm( wprtp_mc_zt )
    wpthlp_mc = zt2zm( wpthlp_mc_zt )
    rtpthlp_mc = zt2zm( rtpthlp_mc_zt )
  end subroutine update_xp2_mc

!===============================================================================

end module advance_xp2_xpyp_module
