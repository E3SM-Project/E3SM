!-----------------------------------------------------------------------
! $Id: advance_xp2_xpyp_module.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
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
             term_ta_lhs, & 
             term_ta_lhs_upwind, & 
             term_ta_rhs, & 
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
  subroutine advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp, thlm, &
                               wpthlp, wpthvp, um, vm, wp2, wp2_zt, &
                               wp3, upwp, vpwp, sigma_sqd_w, Skw_zm, &
                               Kh_zt, rtp2_forcing, thlp2_forcing, &
                               rtpthlp_forcing, rho_ds_zm, rho_ds_zt, &
                               invrs_rho_ds_zm, thv_ds_zm, &
                               Lscale, wp3_on_wp2, wp3_on_wp2_zt, &
                               l_iter, dt, &
                               sclrm, wpsclrp, & 
                               rtp2, thlp2, rtpthlp, up2, vp2,  & 
                               err_code, & 
                               sclrp2, sclrprtp, sclrpthlp )

    ! Description:
    ! Subprogram to diagnose variances by solving steady-state equations

    ! References:
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
      w_tol_sqd, & 
      fstderr, &
      one, &
      two_thirds, &
      one_half, &
      one_third, &
      zero, &
      zero_threshold

    use model_flags, only: & 
      l_hole_fill, &    ! logical constants
      l_single_C2_Skw

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
      gr,  & ! Variable(s)
      zm2zt ! Procedure(s)

    use clubb_precision, only:  & 
      core_rknd ! Variable(s)

    use clip_explicit, only: & 
      clip_covar,  & ! Procedure(s)
      clip_variance, &
      clip_sclrp2, &
      clip_sclrprtp, &
      clip_sclrpthlp
      
    use stats_type_utilities, only: & 
      stat_modify

    use error_code, only:  & 
      clubb_no_error,  & ! Variable(s)
      clubb_var_out_of_range, &
      clubb_singular_matrix

    use error_code, only:  & 
      fatal_error,    & ! Procedure(s)
      clubb_at_least_debug_level

    use stats_variables, only: & 
      stats_zm, &
      irtp2_cl, & 
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
      Kh_zt,           & ! Eddy diffusivity on thermo. levels    [m^2/s]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
      rho_ds_zm,       & ! Dry, static density on momentum levs. [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ mom. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on mom. levs. [K]
      Lscale,          & ! Mixing length                         [m]
      wp3_on_wp2,      & ! Smoothed version of <w'^3>/<w'^2> zm  [m/s]
      wp3_on_wp2_zt      ! Smoothed version of <w'^3>/<w'^2> zt  [m/s]

    logical, intent(in) :: l_iter ! Whether variances are prognostic

    real( kind = core_rknd ), intent(in) :: &
      dt             ! Model timestep                                [s]

    ! Passive scalar input
    real( kind = core_rknd ), intent(in), dimension(gr%nz, sclr_dim) ::  & 
      sclrm, wpsclrp

    ! Input/Output variables
    ! An attribute of (inout) is also needed to import the value of the variances
    ! at the surface.  Brian.  12/18/05.
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      rtp2,    & ! <r_t'^2>                      [(kg/kg)^2]
      thlp2,   & ! <th_l'^2>                     [K^2]
      rtpthlp, & ! <r_t'th_l'>                   [(kg K)/kg]
      up2,     & ! <u'^2>                        [m^2/s^2]
      vp2        ! <v'^2>                        [m^2/s^2]

    ! Output variable for singular matrices
    integer, intent(inout) :: err_code

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

    integer, dimension(5+1) :: & 
      err_code_array ! Array containing the error codes for each variable

    ! Eddy Diffusion for Variances and Covariances.
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      Kw2,      & ! For rtp2, thlp2, rtpthlp, and passive scalars  [m^2/s]
      Kw9         ! For up2 and vp2                                [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      a1_zt,      & ! a_1 interpolated to thermodynamic levels       [-]
      wprtp_zt,   & ! w'r_t' interpolated to thermodynamic levels    [(kg/kg) m/s]
      wpthlp_zt     ! w'th_l' interpolated to thermodyamnic levels   [K m/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rtpthlp_chnge  ! Net change in r_t'th_l' due to clipping [(kg/kg) K]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      sclrprtp_chnge, & ! Net change in sclr'r_t' due to clipping  [{units vary}]
      sclrpthlp_chnge   ! Net change in sclr'th_l' due to clipping [{units vary}]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      sclrp2_forcing,    & ! <sclr'^2> forcing (momentum levels)    [units vary]
      sclrprtp_forcing,  & ! <sclr'r_t'> forcing (momentum levels)  [units vary]
      sclrpthlp_forcing    ! <sclr'th_l'> forcing (momentum levels) [units vary]

    logical :: l_scalar_calc, l_first_clip_ts, l_last_clip_ts

    ! Loop indices
    integer :: i, k

    !---------------------------- Begin Code ----------------------------------

    if ( clubb_at_least_debug_level( 2 ) ) then
      ! Assertion check for C5
      if ( C5 > one .or. C5 < zero ) then
        write(fstderr,*) "The C5 variable is outside the valid range"
        err_code = clubb_var_out_of_range
        return
      end if
    end if

    if ( l_single_C2_Skw ) then
      ! Use a single value of C2 for all equations.
      C2rt_1d(1:gr%nz)  & 
      = C2b + (C2-C2b) *exp( -one_half * (Skw_zm(1:gr%nz)/C2c)**2 )

      C2thl_1d   = C2rt_1d
      C2rtthl_1d = C2rt_1d

      C2sclr_1d  = C2rt_1d
    else
      ! Use 3 different values of C2 for rtp2, thlp2, rtpthlp.
      C2rt_1d(1:gr%nz)    = C2rt
      C2thl_1d(1:gr%nz)   = C2thl
      C2rtthl_1d(1:gr%nz) = C2rtthl

      C2sclr_1d(1:gr%nz)  = C2rt  ! Use rt value for now
    end if

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
    ! located on the momentum levels).
    a1(1:gr%nz) = one / ( one - sigma_sqd_w(1:gr%nz) )


    ! Interpolate a_1, w'r_t', w'th_l', u'w', and v'w' from the momentum levels
    ! to the thermodynamic levels.  These will be used for the turbulent
    ! advection (ta) terms in each equation.
    a1_zt     = max( zm2zt( a1 ), zero_threshold )   ! Positive definite quantity
    wprtp_zt  = zm2zt( wprtp )
    wpthlp_zt = zm2zt( wpthlp )
    upwp_zt   = zm2zt( upwp )
    vpwp_zt   = zm2zt( vpwp )

    ! Initialize tridiagonal solutions to valid

    err_code_array(:) = clubb_no_error


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
    call xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, &  ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2, &          ! Intent(in)
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, &  ! Intent(in)
                       C2rt_1d, nu2_vert_res_dep, beta, &        ! Intent(in)
                       lhs )                                     ! Intent(out)


    call xp2_xpyp_rhs( xp2_xpyp_rtp2, dt, l_iter, a1, a1_zt, &  ! Intent(in)
                       wp2_zt, wprtp, wprtp_zt, wp3_on_wp2, &   ! Intent(in)
                       wp3_on_wp2_zt, wprtp, wprtp_zt, &        ! Intent(in)
                       rtm, rtm, rtp2, rtp2_forcing, &          ! Intent(in)
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! Intent(in)
                       C2rt_1d, tau_zm, rt_tol**2, beta, &      ! Intent(in)
                       rhs )                                    ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_rtp2, 1, & ! Intent(in)
                         rhs, lhs, rtp2, &   ! Intent(inout)
                         err_code_array(1) ) ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_rtp2, rtp2 ) ! Intent(in)
    end if

    !!!!!***** th_l'^2 *****!!!!!

    ! Implicit contributions to term thlp2
    call xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, & ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2, &         ! Intent(in)
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! Intent(in)
                       C2thl_1d, nu2_vert_res_dep, beta, &      ! Intent(in)
                       lhs )                                    ! Intent(out)

    ! Explicit contributions to thlp2
    call xp2_xpyp_rhs( xp2_xpyp_thlp2, dt, l_iter, a1, a1_zt, & ! Intent(in)
                       wp2_zt, wpthlp, wpthlp_zt, wp3_on_wp2, & ! Intent(in)
                       wp3_on_wp2_zt, wpthlp, wpthlp_zt, &      ! Intent(in)
                       thlm, thlm, thlp2, thlp2_forcing, &      ! Intent(in)
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! Intent(in)
                       C2thl_1d, tau_zm, thl_tol**2, beta, &    ! Intent(in)
                       rhs )                                    ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_thlp2, 1, & ! Intent(in)
                         rhs, lhs, thlp2, &   ! Intent(inout)
                         err_code_array(2) )  ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_thlp2, thlp2 ) ! Intent(in)
    end if


    !!!!!***** r_t'th_l' *****!!!!!

    ! Implicit contributions to term rtpthlp
    call xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, & ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw2, &         ! Intent(in)
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! Intent(in)
                       C2rtthl_1d, nu2_vert_res_dep, beta, &    ! Intent(in)
                       lhs )                                    ! Intent(out)

    ! Explicit contributions to rtpthlp
    call xp2_xpyp_rhs( xp2_xpyp_rtpthlp, dt, l_iter, a1, a1_zt, &  ! Intent(in)
                       wp2_zt, wprtp, wprtp_zt, wp3_on_wp2, &      ! Intent(in)
                       wp3_on_wp2_zt, wpthlp, wpthlp_zt, &         ! Intent(in)
                       rtm, thlm, rtpthlp, rtpthlp_forcing, &      ! Intent(in)
                       rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &    ! Intent(in)
                       C2rtthl_1d, tau_zm, zero_threshold, beta, & ! Intent(in)
                       rhs )                                       ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_rtpthlp, 1, & ! Intent(in)
                         rhs, lhs, rtpthlp, &   ! Intent(inout)
                         err_code_array(3) )    ! Intent(out)

    if ( l_stats_samp ) then
      call xp2_xpyp_implicit_stats( xp2_xpyp_rtpthlp, rtpthlp ) ! Intent(in)
    end if


    !!!!!***** u'^2 / v'^2 *****!!!!!

    ! Implicit contributions to term up2/vp2
    call xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, & ! Intent(in)
                       a1, a1_zt, tau_zm, wm_zm, Kw9,  &        ! Intent(in)
                       rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! Intent(in)
                       C4_C14_1d, nu9_vert_res_dep, beta, &     ! Intent(in)
                       lhs )                                    ! Intent(out)

    ! Explicit contributions to up2
    call xp2_xpyp_uv_rhs( xp2_xpyp_up2, dt, l_iter, a1, a1_zt, wp2, & ! Intent(in)
                          wp2_zt, wpthvp, Lscale, wp3_on_wp2_zt, &  ! Intent(in)
                          wp3_on_wp2, C4_C14_1d, tau_zm, &          ! Intent(in)
                          um, vm, upwp, upwp_zt, vpwp, vpwp_zt, &   ! Intent(in)
                          up2, vp2, rho_ds_zt, invrs_rho_ds_zm, &   ! Intent(in)
                          rho_ds_zm, & ! Intent(in)
                          thv_ds_zm, C4, C5, C14, beta, &           ! Intent(in)
                          uv_rhs(:,1) )                             ! Intent(out)

    ! Explicit contributions to vp2
    call xp2_xpyp_uv_rhs( xp2_xpyp_vp2, dt, l_iter, a1, a1_zt, wp2, & ! Intent(in)
                          wp2_zt, wpthvp, Lscale, wp3_on_wp2_zt, &  ! Intent(in)
                          wp3_on_wp2, C4_C14_1d, tau_zm, &          ! Intent(in)
                          vm, um, vpwp, vpwp_zt, upwp, upwp_zt, &   ! Intent(in)
                          vp2, up2, rho_ds_zt, invrs_rho_ds_zm, &   ! Intent(in)
                          rho_ds_zm, & ! Intent(in)
                          thv_ds_zm, C4, C5, C14, beta, &           ! Intent(in)
                          uv_rhs(:,2) )                             ! Intent(out)

    ! Solve the tridiagonal system
    call xp2_xpyp_solve( xp2_xpyp_up2_vp2, 2,       & ! Intent(in)
                         uv_rhs, lhs,               & ! Intent(inout)
                         uv_solution, err_code_array(4) )  ! Intent(out)

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

    threshold = rt_tol**2

    call clip_variance( xp2_xpyp_rtp2, dt, threshold, & ! Intent(in)
                        rtp2 )                          ! Intent(inout)

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

    threshold = thl_tol**2

    call clip_variance( xp2_xpyp_thlp2, dt, threshold, & ! Intent(in)
                        thlp2 )                          ! Intent(inout)


    ! Clipping for u'^2

    !threshold = zero_threshold
    threshold = w_tol_sqd

    call clip_variance( xp2_xpyp_up2, dt, threshold, & ! Intent(in)
                        up2 )                          ! Intent(inout)


    ! Clipping for v'^2

    !threshold = zero_threshold
    threshold = w_tol_sqd

    call clip_variance( xp2_xpyp_vp2, dt, threshold, & ! Intent(in)
                        vp2 )                          ! Intent(inout)


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

      call xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, & ! Intent(in)
                         a1, a1_zt, tau_zm, wm_zm, Kw2,         & ! Intent(in)
                         rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm, & ! Intent(in)
                         C2sclr_1d, nu2_vert_res_dep, beta,     & ! Intent(in)
                         lhs )                                    ! Intent(out)


      ! Explicit contributions to passive scalars

      do i = 1, sclr_dim, 1

        ! Interpolate w'sclr' from momentum levels to thermodynamic
        ! levels.  These will be used for the turbulent advection (ta)
        ! terms in each equation.
        wpsclrp_zt = zm2zt( wpsclrp(:,i) )

        ! Forcing for <sclr'^2>.
        sclrp2_forcing = zero

        !!!!!***** sclr'^2 *****!!!!!

        call xp2_xpyp_rhs( xp2_xpyp_sclrp2, dt, l_iter, a1, a1_zt, & ! In
                           wp2_zt, wpsclrp(:,i), wpsclrp_zt, wp3_on_wp2, & ! In
                           wp3_on_wp2_zt, wpsclrp(:,i), wpsclrp_zt, & ! In
                           sclrm(:,i), sclrm(:,i), sclrp2(:,i), sclrp2_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &  ! In
                           C2sclr_1d, tau_zm, sclr_tol(i)**2, beta, & ! In
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

        call xp2_xpyp_rhs( xp2_xpyp_sclrprtp, dt, l_iter, a1, a1_zt, & ! In
                           wp2_zt, wpsclrp(:,i), wpsclrp_zt, wp3_on_wp2, & ! In
                           wp3_on_wp2_zt, wprtp, wprtp_zt, & ! In
                           sclrm(:,i), rtm, sclrprtp(:,i), sclrprtp_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &  ! In
                           C2sclr_1d, tau_zm, threshold, beta, &  ! In
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

        call xp2_xpyp_rhs( xp2_xpyp_sclrpthlp, dt, l_iter, a1, a1_zt, & ! In
                           wp2_zt, wpsclrp(:,i), wpsclrp_zt, wp3_on_wp2, & ! In
                           wp3_on_wp2_zt, wpthlp, wpthlp_zt, & ! In
                           sclrm(:,i), thlm, sclrpthlp(:,i), sclrpthlp_forcing, & ! In
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & ! In
                           C2sclr_1d, tau_zm, threshold, beta, & ! In
                           sclr_rhs(:,i+2*sclr_dim) ) ! Out


      enddo ! 1..sclr_dim


      ! Solve the tridiagonal system

      call xp2_xpyp_solve( xp2_xpyp_scalars, 3*sclr_dim, &   ! Intent(in)
                           sclr_rhs, lhs, sclr_solution, &  ! Intent(inout)
                           err_code_array(6) )                   ! Intent(out)

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


    ! Check for singular matrices and bad LAPACK arguments
    if ( any( fatal_error( err_code_array ) ) ) then
      err_code = clubb_singular_matrix
    end if

    if ( fatal_error( err_code ) .and.  & 
         clubb_at_least_debug_level( 1 ) ) then

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

    return
  end subroutine advance_xp2_xpyp

  !=============================================================================
  subroutine xp2_xpyp_lhs( dt, l_iter, wp3_on_wp2_zt, wp3_on_wp2, &
                           a1, a1_zt, tau_zm, wm_zm, Kw,  &
                           rho_ds_zt, rho_ds_zm, invrs_rho_ds_zm,  &
                           Cn, nu, beta, lhs )

    ! Description:
    ! Compute LHS tridiagonal matrix for a variance or covariance term

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only:  &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_upwind_xpyp_ta ! Constant(s)

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use diffusion, only:  & 
        diffusion_zm_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zm_lhs ! Procedure(s)

    use stats_variables, only: & 
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
        l_stats_samp, &
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

    use advance_helper_module, only: set_boundary_conditions_lhs


    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      dt ! Timestep length                             [s]

    logical, intent(in) :: & 
      l_iter  ! Whether the variances are prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wp3_on_wp2,      & ! Smoothed w'^3 / w'^2 (moment. levels)       [m/s]
      wp3_on_wp2_zt,   & ! Smoothed w'^3 / w'^2 (thermo. levels)       [m/s]
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)      [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels    [-]
      tau_zm,          & ! Time-scale tau on momentum levels           [s]
      wm_zm,           & ! w wind component on momentum levels         [m/s]
      Kw,              & ! Coefficient of eddy diffusivity (all vars.) [m^2/s]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levs.  [m^3/kg]
      Cn                 ! Coefficient C_n                             [-]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      nu         ! Background constant coef. of eddy diff.        [-]
    real( kind = core_rknd ), intent(in) :: &
      beta       ! Constant model parameter beta                  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(3,gr%nz), intent(out) :: & 
      lhs ! Implicit contributions to the term

    ! Local Variables

    ! Array indices
    integer :: k, kp1, km1, low_bound, high_bound

    real( kind = core_rknd ), dimension(3) :: & 
      tmp

    ! Initialize LHS matrix to 0.
    lhs = zero

    ! Setup LHS of the tridiagonal system
    do k = 2, gr%nz-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! LHS mean advection (ma) term.
      lhs(kp1_mdiag:km1_mdiag,k) & 
      = lhs(kp1_mdiag:km1_mdiag,k) & 
      + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )

      ! LHS turbulent advection (ta) term.
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
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
      if ( .not. l_upwind_xpyp_ta ) then
        lhs(kp1_mdiag:km1_mdiag,k)  & 
        = lhs(kp1_mdiag:km1_mdiag,k)  &
        + gamma_over_implicit_ts  &
        * term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                       rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                       a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
      else
        lhs(kp1_mdiag:km1_mdiag,k)  & 
        = lhs(kp1_mdiag:km1_mdiag,k)  &
        + gamma_over_implicit_ts  &
        * term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                              wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                              gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                              invrs_rho_ds_zm(k), &
                              rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
      end if

      ! LHS dissipation term 1 (dp1)
      ! (combined with pressure term 1 (pr1) for u'^2 and v'^2).
      ! Note:  An "over-implicit" weighted time step is applied to this term
      !        (and to pressure term 1 for u'^2 and v'^2).
      lhs(k_mdiag,k)  & 
      = lhs(k_mdiag,k)  &
      + gamma_over_implicit_ts  &
      * term_dp1_lhs( Cn(k), tau_zm(k) )

      ! LHS eddy diffusion term: dissipation term 2 (dp2).
      lhs(kp1_mdiag:km1_mdiag,k)  & 
      = lhs(kp1_mdiag:km1_mdiag,k)  & 
      + diffusion_zm_lhs( Kw(k), Kw(kp1), nu, & 
                          gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                          gr%invrs_dzm(k), k )

      ! LHS time tendency.
      if ( l_iter ) then
        lhs(k_mdiag,k) = lhs(k_mdiag,k) + ( one / dt )
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
          tmp(1)  &
          = gamma_over_implicit_ts  &
          * term_dp1_lhs( Cn(k), tau_zm(k) )
          zmscr01(k) = -tmp(1)
        endif

        if ( irtp2_dp2 + ithlp2_dp2 + irtpthlp_dp2 +  & 
             iup2_dp2 + ivp2_dp2 > 0 ) then
          tmp(1:3) & 
          = diffusion_zm_lhs( Kw(k), Kw(kp1), nu, & 
                              gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                              gr%invrs_dzm(k), k )
          zmscr02(k) = -tmp(3)
          zmscr03(k) = -tmp(2)
          zmscr04(k) = -tmp(1)
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( irtp2_ta + ithlp2_ta + irtpthlp_ta + & 
             iup2_ta + ivp2_ta > 0 ) then
          if ( .not. l_upwind_xpyp_ta ) then
            tmp(1:3)  &
            = gamma_over_implicit_ts  &
            * term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                           rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                           a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
          else
              tmp(1:3)  &
              = gamma_over_implicit_ts  &
              * term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                    wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                    gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                    invrs_rho_ds_zm(k), &
                                    rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
          end if

          zmscr05(k) = -tmp(3)
          zmscr06(k) = -tmp(2)
          zmscr07(k) = -tmp(1)
        endif

        if ( irtp2_ma + ithlp2_ma + irtpthlp_ma + & 
             iup2_ma + ivp2_ma > 0 ) then
          tmp(1:3) & 
          = term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )
          zmscr08(k) = -tmp(3)
          zmscr09(k) = -tmp(2)
          zmscr10(k) = -tmp(1)
        endif

      endif ! l_stats_samp

    enddo ! k=2..gr%nz-1


    ! Boundary Conditions
    ! These are set so that the surface_varnce value of the variances and
    ! covariances can be used at the lowest boundary and the values of those
    ! variables can be set to their respective threshold minimum values at the
    ! top boundary.  Fixed-point boundary conditions are used for both the
    ! variances and the covariances.
    low_bound = 1
    high_bound = gr%nz

    call set_boundary_conditions_lhs( k_mdiag, low_bound, high_bound, lhs )

    return

  end subroutine xp2_xpyp_lhs

  !=============================================================================
  subroutine xp2_xpyp_solve( solve_type, nrhs, rhs, lhs, xapxbp, err_code )

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

    use error_code, only: &
        clubb_no_error ! Constant

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

    integer, intent(out) :: & 
      err_code ! Returns an error code in the event of a singular matrix

    ! Local variables
    real( kind = core_rknd ) :: rcond  ! Est. of the reciprocal of the condition # on the matrix

    integer ::  ixapxbp_matrix_condt_num ! Stat index

    character(len=10) :: &
      solve_type_str ! solve_type in string format for debug output purposes

    ! --- Begin Code ---
    err_code = clubb_no_error  ! Initialize to the value for no errors

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
           ( solve_type_str, gr%nz, nrhs, &                                      ! Intent(in) 
             lhs(kp1_mdiag,:), lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  & ! Intent(inout)
             xapxbp(:,1:nrhs), rcond, err_code )                                   ! Intent(out)

      ! Est. of the condition number of the variance LHS matrix
      call stat_update_var_pt( ixapxbp_matrix_condt_num, 1, one / rcond, &  ! Intent(in)
                               stats_sfc )                          ! Intent(inout)

    else
      call tridag_solve & 
           ( solve_type_str, gr%nz, nrhs, lhs(kp1_mdiag,:),  &    ! Intent(in)
             lhs(k_mdiag,:), lhs(km1_mdiag,:), rhs(:,1:nrhs),  &    ! Intent(inout)
             xapxbp(:,1:nrhs), err_code )                           ! Intent(out)
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
      call stat_end_update_pt( ixapxbp_dp1, k, &            ! Intent(in)
                               zmscr01(k) * xapxbp(k), &  ! Intent(in)
                               stats_zm )                         ! Intent(inout)

      ! x'y' term dp2 is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_dp2, k, &            ! Intent(in)
                                 zmscr02(k) * xapxbp(km1) & ! Intent(in)
                               + zmscr03(k) * xapxbp(k) & 
                               + zmscr04(k) * xapxbp(kp1), &
                               stats_zm )                         ! Intent(inout)

      ! x'y' term ta has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_ta, k, &              ! Intent(in)
                                 zmscr05(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr06(k) * xapxbp(k) &  
                               + zmscr07(k) * xapxbp(kp1), &
                               stats_zm )                          ! Intent(inout)

      ! x'y' term ma is completely implicit; call stat_update_var_pt.
      call stat_update_var_pt( ixapxbp_ma, k, &              ! Intent(in)
                                 zmscr08(k) * xapxbp(km1) &  ! Intent(in)
                               + zmscr09(k) * xapxbp(k) & 
                               + zmscr10(k) * xapxbp(kp1), &
                               stats_zm )                          ! Intent(inout)

      ! x'y' term pr1 has both implicit and explicit components;
      ! call stat_end_update_pt.
      call stat_end_update_pt( ixapxbp_pr1, k, &            ! Intent(in)
                               zmscr11(k) * xapxbp(k), &  ! Intent(in)
                               stats_zm )                         ! Intent(inout)

    end do ! k=2..gr%nz-1

    return
  end subroutine xp2_xpyp_implicit_stats

  !=============================================================================
  subroutine xp2_xpyp_uv_rhs( solve_type, dt, l_iter, a1, a1_zt, wp2, &
                              wp2_zt, wpthvp, Lscale, wp3_on_wp2_zt, &
                              wp3_on_wp2, C4_C14_1d, tau_zm,  & 
                              xam, xbm, wpxap, wpxap_zt, wpxbp, wpxbp_zt, &
                              xap2, xbp2, rho_ds_zt, invrs_rho_ds_zm, &
                              rho_ds_zm, &
                              thv_ds_zm, C4, C5, C14, beta, &
                              rhs )

    ! Description:
    ! Explicit contributions to u'^2 or v'^2
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only:  & 
        gamma_over_implicit_ts, & ! Constant(s)
        w_tol_sqd, &
        one, &
        two_thirds, &
        one_third, &
        zero

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
        iup2_ta, & 
        iup2_tp, & 
        iup2_dp1, & 
        iup2_pr1, & 
        iup2_pr2, & 
        stats_zm, & 
        zmscr01, & 
        zmscr11, & 
        l_stats_samp

    implicit none

    ! Input Variables
    integer, intent(in) :: solve_type

    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                              [s]

    logical, intent(in) :: & 
      l_iter  ! Whether x is prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)      [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels    [-]
      wp2,             & ! w'^2 (momentum levels)                      [m^2/s^2]
      wp2_zt,          & ! w'^2 interpolated to thermodynamic levels   [m^2/s^2]
      wpthvp,          & ! w'th_v' (momentum levels)                   [K m/s]
      Lscale,          & ! Mixing Length                               [m]
      wp3_on_wp2,      & ! Smoothed w'^3 / w'^2 on momentum levels     [m/s]
      wp3_on_wp2_zt,   & ! Smoothed w'^3 / w'^2 on thermo. levels      [m/s]
      C4_C14_1d,       & ! Combination of model params. C_4 and C_14   [-]
      tau_zm,          & ! Time-scale tau on momentum levels           [s]
      xam,             & ! x_am (thermodynamic levels)                 [m/s]
      xbm,             & ! x_bm (thermodynamic levels)                 [m/s]
      wpxap,           & ! w'x_a' (momentum levels)                    [m^2/s^2]
      wpxap_zt,        & ! w'x_a' interpolated to thermodynamic levels [m^2/s^2]
      wpxbp,           & ! w'x_b' (momentum levels)                    [m^2/s^2]
      wpxbp_zt,        & ! w'x_b' interpolated to thermodynamic levels [m^2/s^2]
      xap2,            & ! x_a'^2 (momentum levels)                    [m^2/s^2]
      xbp2,            & ! x_b'^2 (momentum levels)                    [m^2/s^2]
      rho_ds_zt,       & ! Dry, static density on thermodynamic levels [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels      [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levs.  [m^3/kg]
      thv_ds_zm          ! Dry, base-state theta_v on momentum levels  [K]

    real( kind = core_rknd ), intent(in) :: & 
      C4,              & ! Model parameter C_4                         [-]
      C5,              & ! Model parameter C_5                         [-]
      C14,             & ! Model parameter C_{14}                      [-]
      beta               ! Model parameter beta                        [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz,1), intent(out) :: & 
      rhs    ! Explicit contributions to x variance/covariance terms

    ! Local Variables

    ! Array indices
    integer :: k, kp1, km1

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3) :: lhs_fnc_output

    real( kind = core_rknd ) :: tmp

    integer :: & 
      ixapxbp_ta, & 
      ixapxbp_tp, & 
      ixapxbp_dp1, & 
      ixapxbp_pr1, & 
      ixapxbp_pr2

    !----------------------------- Begin Code ----------------------------------

    select case ( solve_type )
    case ( xp2_xpyp_vp2 )
      ixapxbp_ta  = ivp2_ta
      ixapxbp_tp  = ivp2_tp
      ixapxbp_dp1 = ivp2_dp1
      ixapxbp_pr1 = ivp2_pr1
      ixapxbp_pr2 = ivp2_pr2
    case ( xp2_xpyp_up2 )
      ixapxbp_ta  = iup2_ta
      ixapxbp_tp  = iup2_tp
      ixapxbp_dp1 = iup2_dp1
      ixapxbp_pr1 = iup2_pr1
      ixapxbp_pr2 = iup2_pr2
    case default ! No budgets for passive scalars
      ixapxbp_ta  = 0
      ixapxbp_tp  = 0
      ixapxbp_dp1 = 0
      ixapxbp_pr1 = 0
      ixapxbp_pr2 = 0
    end select


    ! Initialize RHS vector to 0.
    rhs = zero

    do k = 2, gr%nz-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! RHS turbulent advection (ta) term.
      rhs(k,1)  & 
      = rhs(k,1)  & 
      + term_ta_rhs( wp2_zt(kp1), wp2_zt(k),  &
                     wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                     rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                     a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k),  &
                     wpxap_zt(kp1), wpxap_zt(k), gr%invrs_dzm(k), beta )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS turbulent advection (ta) term.
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
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
      if ( .not. l_upwind_xpyp_ta ) then
        lhs_fnc_output(1:3)  &
        = term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                       rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                       a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
      else
        lhs_fnc_output(1:3)  &
        = term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                              wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                              gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                              invrs_rho_ds_zm(k), &
                              rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
      end if

      rhs(k,1)  &
      = rhs(k,1)  &
      + ( one - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * xap2(kp1)  &
          - lhs_fnc_output(2) * xap2(k)  &
          - lhs_fnc_output(3) * xap2(km1) )

      ! RHS turbulent production (tp) term.
      rhs(k,1)  &
      = rhs(k,1)  &
      + ( one - C5 )  & 
         * term_tp( xam(kp1), xam(k), xam(kp1), xam(k), & 
                    wpxap(k), wpxap(k), gr%invrs_dzm(k) )

      ! RHS pressure term 1 (pr1) (and dissipation term 1 (dp1)).
      rhs(k,1)  &
      = rhs(k,1)  &
      + term_pr1( C4, C14, xbp2(k), wp2(k), tau_zm(k) )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS dissipation term 1 (dp1) and pressure term 1 (pr1).
      !
      ! Note:  An "over-implicit" weighted time step is applied to these terms.
      lhs_fnc_output(1)  &
      = term_dp1_lhs( C4_C14_1d(k), tau_zm(k) )
      rhs(k,1)  &
      = rhs(k,1)  &
      + ( one - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * xap2(k) )

      ! RHS pressure term 2 (pr2).
      rhs(k,1)  &
      = rhs(k,1)  &
      + term_pr2( C5, thv_ds_zm(k), wpthvp(k), wpxap(k), wpxbp(k), &
                  xam, xbm, gr%invrs_dzm(k), kp1, k, &
                  Lscale(kp1), Lscale(k), wp2_zt(kp1), wp2_zt(k) )

      ! RHS time tendency.
      if ( l_iter ) then
        rhs(k,1) = rhs(k,1) + one/dt * xap2(k)
      endif

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for up2 or vp2.

        ! x'y' term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_ta_rhs.
        call stat_begin_update_pt( ixapxbp_ta, k, &                 ! Intent(in) 
        -term_ta_rhs( wp2_zt(kp1), wp2_zt(k), &   ! Intent(in)
                      wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                      rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k), &
                      a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                      wpxap_zt(kp1), wpxap_zt(k), gr%invrs_dzm(k), beta ), &
                                   stats_zm )                             ! Intent(inout)

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) term).
        if ( .not. l_upwind_xpyp_ta ) then
          lhs_fnc_output(1:3)  &
          = term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                         rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                         a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
        else ! turbulent advection is using an upwind discretization
          lhs_fnc_output(1:3)  &
          = term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                invrs_rho_ds_zm(k), &
                                rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
        end if ! ~l_upwind_xpyp_ta

        call stat_modify_pt( ixapxbp_ta, k,  &          ! Intent(in)
              + ( one - gamma_over_implicit_ts )  &     ! Intent(in)
              * ( - lhs_fnc_output(1) * xap2(kp1)  &
                  - lhs_fnc_output(2) * xap2(k)  &
                  - lhs_fnc_output(3) * xap2(km1) ),  &
                             stats_zm )                       ! Intent(inout)

        if ( ixapxbp_dp1 > 0 ) then
          ! Note:  The function term_pr1 is the explicit component of a
          !        semi-implicit solution to dp1 and pr1.
          ! Record the statistical contribution of the implicit component of
          ! term dp1 for up2 or vp2.  This will overwrite anything set
          ! statistically in xp2_xpyp_lhs for this term.
          ! Note:  To find the contribution of x'y' term dp1, substitute
          !        (2/3)*C_4 for the C_n input to function term_dp1_lhs.
          ! Note:  An "over-implicit" weighted time step is applied to this
          !        term.  A weighting factor of greater than 1 may be used to
          !        make the term more numerically stable (see note above for
          !        RHS turbulent advection (ta) term).
          tmp  &
          = gamma_over_implicit_ts  &
          * term_dp1_lhs( two_thirds*C4, tau_zm(k) )
          zmscr01(k) = -tmp
          ! Statistical contribution of the explicit component of term dp1 for
          ! up2 or vp2.
          ! x'y' term dp1 has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_pr1.
          ! Note:  To find the contribution of x'y' term dp1, substitute 0 for
          !        the C_14 input to function term_pr1.
          call stat_begin_update_pt( ixapxbp_dp1, k, &              ! Intent(in)
               -term_pr1( C4, zero, xbp2(k), wp2(k), tau_zm(k) ), &  ! Intent(in)
                                     stats_zm )                           ! Intent(inout)

          ! Note:  An "over-implicit" weighted time step is applied to this
          !        term.  A weighting factor of greater than 1 may be used to
          !        make the term more numerically stable (see note above for
          !        RHS turbulent advection (ta) term).
          lhs_fnc_output(1)  &
          = term_dp1_lhs( two_thirds*C4, tau_zm(k) )
          call stat_modify_pt( ixapxbp_dp1, k, &        ! Intent(in)
                + ( one - gamma_over_implicit_ts )  &   ! Intent(in)
                * ( - lhs_fnc_output(1) * xap2(k) ),  & ! Intent(in)
                                     stats_zm )               ! Intent(inout)

        endif

        if ( ixapxbp_pr1 > 0 ) then
          ! Note:  The function term_pr1 is the explicit component of a
          !        semi-implicit solution to dp1 and pr1.
          ! Statistical contribution of the implicit component of term pr1 for
          ! up2 or vp2.
          ! Note:  To find the contribution of x'y' term pr1, substitute
          !        (1/3)*C_14 for the C_n input to function term_dp1_lhs.
          ! Note:  An "over-implicit" weighted time step is applied to this
          !        term.  A weighting factor of greater than 1 may be used to
          !        make the term more numerically stable (see note above for
          !        RHS turbulent advection (ta) term).
          tmp  &
          = gamma_over_implicit_ts  &
          * term_dp1_lhs( one_third*C14, tau_zm(k) )
          zmscr11(k) = -tmp
          ! Statistical contribution of the explicit component of term pr1 for
          ! up2 or vp2.
          ! x'y' term pr1 has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_pr1.
          ! Note:  To find the contribution of x'y' term pr1, substitute 0 for
          !        the C_4 input to function term_pr1.
          call stat_begin_update_pt( ixapxbp_pr1, k, &                 ! Intent(in)  
               -term_pr1( zero, C14, xbp2(k), wp2(k), tau_zm(k) ), &    ! Intent(in)
                                     stats_zm )                              ! Intent(inout)

          ! Note:  An "over-implicit" weighted time step is applied to this
          !        term.  A weighting factor of greater than 1 may be used to
          !        make the term more numerically stable (see note above for
          !        RHS turbulent advection (ta) term).
          lhs_fnc_output(1)  &
          = term_dp1_lhs( one_third*C14, tau_zm(k) )
          call stat_modify_pt( ixapxbp_pr1, k, &        ! Intent(in)
                + ( one - gamma_over_implicit_ts )  &   ! Intent(in)
                * ( - lhs_fnc_output(1) * xap2(k) ),  & ! Intent(in)
                                     stats_zm )               ! Intent(inout)

        endif

        ! x'y' term pr2 is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_pr2, k, &                            ! Intent(in)
               term_pr2( C5, thv_ds_zm(k), wpthvp(k), wpxap(k), wpxbp(k), &   ! Intent(in)
                         xam, xbm, gr%invrs_dzm(k), kp1, k, &  
                         Lscale(kp1), Lscale(k), wp2_zt(kp1), wp2_zt(k) ), &
               stats_zm )                                                           ! Intent(inout)

        ! x'y' term tp is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_tp, k, &                  ! Intent(in) 
              ( one - C5 ) &                                       ! Intent(in)
               * term_tp( xam(kp1), xam(k), xam(kp1), xam(k), &
                          wpxap(k), wpxap(k), gr%invrs_dzm(k) ), & 
                                 stats_zm )                              ! Intent(inout)

      endif ! l_stats_samp

    enddo ! k=2..gr%nz-1


    ! Boundary Conditions
    ! These are set so that the surface_varnce value of u'^2 or v'^2 can be
    ! used at the lowest boundary and the values of those variables can be
    ! set to their respective threshold minimum values at the top boundary.
    ! Fixed-point boundary conditions are used for the variances.

    rhs(1,1) = xap2(1)
    ! The value of u'^2 or v'^2 at the upper boundary will be set to the
    ! threshold minimum value of w_tol_sqd.
    rhs(gr%nz,1) = w_tol_sqd

    return
  end subroutine xp2_xpyp_uv_rhs

  !=============================================================================
  subroutine xp2_xpyp_rhs( solve_type, dt, l_iter, a1, a1_zt, &
                           wp2_zt, wpxap, wpxap_zt, wp3_on_wp2, &
                           wp3_on_wp2_zt, wpxbp, wpxbp_zt, &
                           xam, xbm, xapxbp, xapxbp_forcing, &
                           rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, & 
                           Cn, tau_zm, threshold, beta, & 
                           rhs )

    ! Description:
    ! Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', sclr'r_t',
    ! sclr'th_l', or sclr'^2.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only: &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_upwind_xpyp_ta ! Constant(s)

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
  
    use advance_helper_module, only: set_boundary_conditions_rhs

    implicit none

    ! Input Variables
    integer, intent(in) :: solve_type

    real( kind = core_rknd ), intent(in) :: & 
      dt                 ! Model timestep                              [s]

    logical, intent(in) :: & 
      l_iter   ! Whether x is prognostic (T/F)

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      a1,              & ! sigma_sqd_w term a_1 (momentum levels)      [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels    [-]
      wp2_zt,          & ! w'^2 interpolated to thermodynamic levels   [m^2/s^2]
      wpxap,           & ! w'x_a' (momentum levels)                    [m/s {x_am units}]
      wpxap_zt,        & ! w'x_a' interpolated to thermodynamic levels [m/s {x_am units}]
      wp3_on_wp2,      & ! w'^3 / w'^2 on momentum levels              [m/s]
      wp3_on_wp2_zt,   & ! w'^3 / w'^2 on thermodynamic levels         [m/s]
      wpxbp,           & ! w'x_b' (momentum levels)                    [m/s {x_bm units}]
      wpxbp_zt,        & ! w'x_b' interpolated to thermodynamic levels [m/s {x_bm units}]
      xam,             & ! x_am (thermodynamic levels)                 [{x_am units}]
      xbm,             & ! x_bm (thermodynamic levels)                 [{x_bm units}]
      xapxbp,          & ! x_a'x_b' (momentum levels)                  [{x_am units}*{x_bm units}]
      xapxbp_forcing,  & ! x_a'x_b' forcing (momentum levels)          [{x_am units}*{x_bm units}/s]
      rho_ds_zm,       & ! Dry, static density on moment. levels       [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels       [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density on momentum levs.  [m^3/kg]
      tau_zm,          & ! Time-scale tau on momentum levels           [s]
      Cn                 ! Coefficient C_n                             [-]

    real( kind = core_rknd ), intent(in) :: &
      threshold,       & ! Smallest allowable mag. value for x_a'x_b'  [{x_am units}
                         !                                              *{x_bm units}]
      beta               ! Model parameter beta                        [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz,1), intent(out) :: & 
      rhs     ! Explicit contributions to x variance/covariance terms

    ! Local Variables

    ! Array indices
    integer :: k, kp1, km1, k_low, k_high

    ! For "over-implicit" weighted time step.
    ! This vector holds output from the LHS (implicit) portion of a term at a
    ! given vertical level.  This output is weighted and applied to the RHS.
    ! This is used if the implicit portion of the term is "over-implicit", which
    ! means that the LHS contribution is given extra weight (>1) in order to
    ! increase numerical stability.  A weighted factor must then be applied to
    ! the RHS in order to balance the weight.
    real( kind = core_rknd ), dimension(3) :: lhs_fnc_output

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


    ! Initialize RHS vector to 0.
    rhs = zero

    do k = 2, gr%nz-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! RHS turbulent advection (ta) term.
      rhs(k,1)  & 
      = rhs(k,1)  & 
      + term_ta_rhs( wp2_zt(kp1), wp2_zt(k),  &
                     wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                     rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                     a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k),  &
                     wpxap_zt(kp1), wpxap_zt(k), gr%invrs_dzm(k), beta )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS turbulent advection (ta) term.
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
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
      if ( .not. l_upwind_xpyp_ta ) then
        lhs_fnc_output(1:3)  &
        = term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                       rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                       a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
      else
        lhs_fnc_output(1:3)  &
        = term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                              wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                              gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                              invrs_rho_ds_zm(k), &
                              rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
      endif

      rhs(k,1)  &
      = rhs(k,1)  &
      + ( one - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * xapxbp(kp1)  &
          - lhs_fnc_output(2) * xapxbp(k)  &
          - lhs_fnc_output(3) * xapxbp(km1) )

      ! RHS turbulent production (tp) term.
      rhs(k,1)  & 
      = rhs(k,1)  & 
      + term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), & 
                 wpxbp(k), wpxap(k), gr%invrs_dzm(k) )

      ! RHS dissipation term 1 (dp1)
      rhs(k,1) &
      = rhs(k,1) + term_dp1_rhs( Cn(k), tau_zm(k), threshold )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS dissipation term 1 (dp1).
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs_fnc_output(1)  &
      = term_dp1_lhs( Cn(k), tau_zm(k) )
      rhs(k,1)  &
      = rhs(k,1)  &
      + ( one - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * xapxbp(k) )

      ! RHS time tendency.
      if ( l_iter ) then
        rhs(k,1) = rhs(k,1) + one/dt * xapxbp(k)
      endif

      ! RHS <x'y'> forcing.
      ! Note: <x'y'> forcing includes the effects of microphysics on <x'y'>.
      rhs(k,1) = rhs(k,1) + xapxbp_forcing(k)


      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for rtp2, thlp2, or rtpthlp.

        ! x'y' term ta has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_ta_rhs.
        call stat_begin_update_pt( ixapxbp_ta, k, &                ! Intent(in) 
        -term_ta_rhs( wp2_zt(kp1), wp2_zt(k), &  ! Intent(in)
                      wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                      rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k), &
                      a1_zt(kp1), a1(k), a1_zt(k), wpxbp_zt(kp1), wpxbp_zt(k), &
                      wpxap_zt(kp1), wpxap_zt(k), gr%invrs_dzm(k), beta ), &
                                   stats_zm )                            ! Intent(inout)

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) term).
        if ( .not. l_upwind_xpyp_ta ) then
          lhs_fnc_output(1:3)  &
          = term_ta_lhs( wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                         rho_ds_zt(kp1), rho_ds_zt(k), invrs_rho_ds_zm(k),  &
                         a1_zt(kp1), a1(k), a1_zt(k), gr%invrs_dzm(k), beta, k )
        else
          lhs_fnc_output(1:3)  &
          = term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                invrs_rho_ds_zm(k), &
                                rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1), beta )
        end if
        call stat_modify_pt( ixapxbp_ta, k,  &            ! Intent(in)
              + ( one - gamma_over_implicit_ts )  &       ! Intent(in)
              * ( - lhs_fnc_output(1) * xapxbp(kp1)  &
                  - lhs_fnc_output(2) * xapxbp(k)  &
                  - lhs_fnc_output(3) * xapxbp(km1) ),  &
                             stats_zm )                         ! Intent(inout)

        ! x'y' term dp1 has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on term_dp1_rhs.
        call stat_begin_update_pt( ixapxbp_dp1, k, &           ! Intent(in)
             -term_dp1_rhs( Cn(k), tau_zm(k), threshold ), &   ! Intent(in)
                                   stats_zm )                        ! Intent(inout)

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for RHS turbulent
        !        advection (ta) term).
        lhs_fnc_output(1)  &
        = term_dp1_lhs( Cn(k), tau_zm(k) )
        call stat_modify_pt( ixapxbp_dp1, k,  &         ! Intent(in)
              + ( one - gamma_over_implicit_ts )  &     ! Intent(in)
              * ( - lhs_fnc_output(1) * xapxbp(k) ),  & ! Intent(in)
                                   stats_zm )                 ! Intent(inout)

        ! rtp2/thlp2 case (1 turbulent production term)
        ! x'y' term tp is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_tp, k, &             ! Intent(in)
              term_tp( xam(kp1), xam(k), xbm(kp1), xbm(k), &  ! Intent(in)
                       wpxbp(k), wpxap(k), gr%invrs_dzm(k) ), & 
                                 stats_zm )                         ! Intent(inout)

        ! rtpthlp case (2 turbulent production terms)
        ! x'y' term tp1 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of x'y' term tp1, substitute 0 for all
        !        the xam inputs and the wpxbp input to function term_tp.
        call stat_update_var_pt( ixapxbp_tp1, k, &    ! Intent(in)
              term_tp( zero, zero, xbm(kp1), xbm(k), &  ! Intent(in)
                       zero, wpxap(k), gr%invrs_dzm(k) ), &
                                 stats_zm )                 ! Intent(inout)

        ! x'y' term tp2 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of x'y' term tp2, substitute 0 for all
        !        the xbm inputs and the wpxap input to function term_tp.
        call stat_update_var_pt( ixapxbp_tp2, k, &    ! Intent(in)
              term_tp( xam(kp1), xam(k), zero, zero, &  ! Intent(in)
                       wpxbp(k), zero, gr%invrs_dzm(k) ), &
                                 stats_zm )                 ! Intent(inout)

        ! x'y' forcing term is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixapxbp_f, k, xapxbp_forcing(k), stats_zm )

      endif ! l_stats_samp

    enddo ! k=2..gr%nz-1


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
    call set_boundary_conditions_rhs( &
            xapxbp(1), k_low, threshold, k_high, &
            rhs(:,1) )

    return
  end subroutine xp2_xpyp_rhs

  !=============================================================================
  pure function term_ta_lhs( wp3_on_wp2_ztp1, wp3_on_wp2_zt, &
                             rho_ds_ztp1, rho_ds_zt, invrs_rho_ds_zm,  &
                             a1_ztp1, a1, a1_zt, invrs_dzm, beta, level )  & 
  result( lhs )

    ! Description:
    ! Turbulent advection of x_a'x_b':  implicit portion of the code.
    !
    ! The d(x_a'x_b')/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'x_a'x_b' )/dz.
    !
    ! A substitution is made in order to close the turbulent advection term,
    ! such that:
    !
    ! w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
    !              + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
    !                * w'x_a' * w'x_b';
    !
    ! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
    ! advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d [ rho_ds * { (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
    !                    + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
    !                      * w'x_a' * w'x_b' } ]
    !     / dz;
    !
    ! which produces an implicit and an explicit portion of this term.  The
    ! implicit portion of this term is:
    !
    ! - (1/rho_ds)
    !   * d [ rho_ds * (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ]
    !     / dz.
    !
    ! Since (1/3)*beta is a constant, it can be pulled outside of the
    ! derivative.  The implicit portion of this term becomes:
    !
    ! - (1/3)*beta/rho_ds
    !   * d [ rho_ds * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'(t+1) ] / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of x_a'x_b' being used is
    ! from the next timestep, which is being advanced to in solving the
    ! d(x_a'x_b')/dt equation.
    !
    ! The implicit portion of this term is discretized as follows:
    !
    ! The values of x_a'x_b' are found on the momentum levels, as are the values
    ! of w'^2 and a_1.  The values of w'^3 are found on the thermodynamic
    ! levels.  Additionally, the values of rho_ds_zt are found on the
    ! thermodynamic levels, and the values of invrs_rho_ds_zm are found on the
    ! momentum levels.  The variables x_a'x_b', w'^2, and a_1 are each
    ! interpolated to the intermediate thermodynamic levels.  The values of the
    ! mathematical expression (called F here) within the dF/dz term are computed
    ! on the thermodynamic levels.  Then the derivative (d/dz) of the
    ! expression (F) is taken over the central momentum level, where it is
    ! multiplied by (1/3)*beta and by invrs_rho_ds_zm, yielding the desired
    ! result.  In this function, the values of F are as follows:
    !
    ! F = rho_ds_zt * a_1(t) * ( w'^3(t) / w'^2(t) ) * x_a'x_b'(t+1);
    !
    ! where the timestep index (t) stands for the index of the current timestep.
    !
    !
    ! ==a1p1========wp2p1========xapxbpp1================================ m(k+1)
    !
    ! ----a1(interp)--wp2(interp)--xapxbp(interp)--wp3p1---rho_ds_ztp1--- t(k+1)
    !
    ! ==a1==========wp2==========xapxbp=======dF/dz====invrs_rho_ds_zm=== m(k)
    !
    ! ----a1(interp)--wp2(interp)--xapxbp(interp)--wp3-----rho_ds_zt----- t(k)
    !
    ! ==a1m1========wp2m1========xapxbpm1================================ m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & ! gr%weights_zm2zt
        gr ! Variable(s)

    use constants_clubb, only: &
        one_third  ! Constant(s)

    use model_flags, only:  &
        l_standard_term_ta

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1, & ! Index for upper momentum level grid weight.
      m_below = 2    ! Index for lower momentum level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      wp3_on_wp2_ztp1, & ! w'^3 / w'^2 (k+1)                           [m/s]
      wp3_on_wp2_zt,   & ! w'^3 / w'^2 (k)                             [m/s]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)  [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum lev (k) [m^3/kg]
      a1_ztp1,         & ! a_1 interpolated to thermo. level (k+1)     [-]
      a1,              & ! a_1(k)                                      [-]
      a1_zt,           & ! a_1 interpolated to thermo. level (k)       [-]
      invrs_dzm,       & ! Inverse of grid spacing                     [1/m]
      beta               ! Model parameter                             [-]

    integer, intent(in) :: & 
      level ! Central momentum level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      tkp1, & ! Thermodynamic level directly above central momentum level.
      tk      ! Thermodynamic level directly below central momentum level.


    ! Thermodynamic level (k+1) is between momentum level (k+1)
    ! and momentum level (k).
    tkp1 = level + 1

    ! Thermodynamic level (k) is between momentum level (k)
    ! and momentum level (k-1).
    tk = level

    if ( l_standard_term_ta ) then

      ! The turbulent advection term is discretized normally, in accordance
      ! with the model equations found in the documentation and the description
      ! listed above.

      ! Momentum superdiagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag)  &
      = + one_third * beta  &
          * invrs_rho_ds_zm  &
            * invrs_dzm  &
              * rho_ds_ztp1 * a1_ztp1  &
              * wp3_on_wp2_ztp1 &
              * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag)  &
      = + one_third * beta  &
          * invrs_rho_ds_zm  &
            * invrs_dzm  &
              * (   rho_ds_ztp1 * a1_ztp1  &
                    * wp3_on_wp2_ztp1 &
                    * gr%weights_zm2zt(m_below,tkp1)  &
                  - rho_ds_zt * a1_zt  &
                    * wp3_on_wp2_zt &
                    * gr%weights_zm2zt(m_above,tk)  &
                )

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag)  &
      = - one_third * beta  &
          * invrs_rho_ds_zm  &
            * invrs_dzm  &
              * rho_ds_zt * a1_zt  &
              * wp3_on_wp2_zt &
              * gr%weights_zm2zt(m_below,tk)

    else

      ! Brian tried a new discretization for the turbulent advection term, for
      ! which the implicit portion of the term is:
      !  - (1/rho_ds)
      !    * d [ rho_ds * a_1 * (1/3)*beta * ( w'^3 / w'^2 ) * x_a'x_b' ] / dz.
      ! In order to help stabilize x_a'x_b', a_1 has been pulled outside the
      ! derivative.

      ! Momentum superdiagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag)  & 
      = + one_third * beta  &
          * invrs_rho_ds_zm * a1  &
            * invrs_dzm  &
              * rho_ds_ztp1  &
              * wp3_on_wp2_ztp1 & 
              * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag)  & 
      = + one_third * beta  &
          * invrs_rho_ds_zm * a1  &
            * invrs_dzm  &
              * (   rho_ds_ztp1  &
                    * wp3_on_wp2_ztp1 &
                    * gr%weights_zm2zt(m_below,tkp1)  &
                  - rho_ds_zt  &
                    * wp3_on_wp2_zt &
                    * gr%weights_zm2zt(m_above,tk)  &
                )

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag)  & 
      = - one_third * beta  &
          * invrs_rho_ds_zm * a1  &
            * invrs_dzm  &
              * rho_ds_zt  &
              * wp3_on_wp2_zt &
              * gr%weights_zm2zt(m_below,tk)

      ! End of Brian's a1 change.  14 Feb 2008.

    endif


    return
  end function term_ta_lhs

  !-----------------------------------------------------------------------------
  pure function term_ta_lhs_upwind( a1_zm, a1_zm_p1, a1_zm_m1, &
                                    wp3_on_wp2_p1, wp3_on_wp2, wp3_on_wp2_m1, &
                                    invrs_dzt, invrs_dzt_p1, &
                                    invrs_rho_ds_zm, & 
                                    rho_ds_zm_p1, rho_ds_zm, rho_ds_zm_m1, beta ) &
  result( lhs )

  ! Description:
  !   Turbulent advection of x_a'x_b' using an upwind differencing
  !   approximation rather than a centered difference.
  ! References:
  !   None
  !-----------------------------------------------------------------------------

    use constants_clubb, only: &
        one_third, & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      a1_zm,             & ! a_1(k) on momentum levels                     [-] 
      a1_zm_p1,          & ! a_1(k+1) on momentum levels                   [-]
      a1_zm_m1,          & ! a_1(k-1) on momentum levels                   [-]
      wp3_on_wp2_p1,     & ! Smoothed wp3 / wp2 on moment. levels (k+1)    [m/s]
      wp3_on_wp2,        & ! Smoothed wp3 / wp2 on moment. levels (k)      [m/s]
      wp3_on_wp2_m1,     & ! Smoothed wp3 / wp2 on moment. levels (k-1)    [m/s]
      invrs_dzt,         & ! Inverse of grid spacing (k)                   [1/m]
      invrs_dzt_p1,      & ! Inverse of grid spacing (k+1)                 [1/m]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum lev (k)   [m^3/kg]
      rho_ds_zm,         & ! Density of air (k)                            [kg/m^3]
      rho_ds_zm_p1,      & ! Density of air (k+1)                          [kg/m^3]
      rho_ds_zm_m1,      & ! Density of air (k-1)                          [kg/m^3]
      beta                 ! Model parameter                               [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs


    if ( wp3_on_wp2 > zero ) then

      ! Momentum main diagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag) = zero

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag)  & 
      = + one_third * beta  &
          * invrs_dzt * invrs_rho_ds_zm  &
          * rho_ds_zm * a1_zm * wp3_on_wp2

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag)  & 
      = - one_third * beta &
          * invrs_dzt * invrs_rho_ds_zm  &
          * rho_ds_zm_m1 * a1_zm_m1 * wp3_on_wp2_m1

    else ! "Wind" is blowing downward

      ! Momentum main diagonal: [ x xapxbp(k+1,<t+1>) ]
      lhs(kp1_mdiag) & 
      = + one_third * beta &
        * invrs_dzt_p1 * invrs_rho_ds_zm &
        * rho_ds_zm_p1 * a1_zm_p1 * wp3_on_wp2_p1

      ! Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
      lhs(k_mdiag) & 
      = - one_third * beta &
        * invrs_dzt_p1 * invrs_rho_ds_zm & 
        * rho_ds_zm * a1_zm * wp3_on_wp2

      ! Momentum subdiagonal: [ x xapxbp(k-1,<t+1>) ]
      lhs(km1_mdiag) = zero

    end if

    return
  end function term_ta_lhs_upwind

  !=============================================================================
  pure function term_ta_rhs( wp2_ztp1, wp2_zt, &
                             wp3_on_wp2_ztp1, wp3_on_wp2_zt,  &
                             rho_ds_ztp1, rho_ds_zt, invrs_rho_ds_zm,  &
                             a1_ztp1, a1, a1_zt, wpxbp_ztp1, wpxbp_zt,  &
                             wpxap_ztp1, wpxap_zt, invrs_dzm, beta )  &
  result( rhs )

    ! Description:
    ! Turbulent advection of x_a'x_b':  explicit portion of the code.
    !
    ! The d(x_a'x_b')/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'x_a'x_b' )/dz.
    !
    ! A substitution is made in order to close the turbulent advection term,
    ! such that:
    !
    ! w'x_a'x_b' = (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
    !              + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
    !                * w'x_a' * w'x_b';
    !
    ! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
    ! advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d [ rho_ds * { (1/3)*beta * a_1 * ( w'^3 / w'^2 ) * x_a'x_b'
    !                    + (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
    !                      * w'x_a' * w'x_b' } ]
    !     / dz;
    !
    ! which produces an implicit and an explicit portion of this term.  The
    ! explicit portion of this term is:
    !
    ! - (1/rho_ds)
    !   * d [ rho_ds * (1-(1/3)*beta) * (a_1)^2 * ( w'^3 / (w'^2)^2 )
    !         * w'x_a' * w'x_b' ] / dz.
    !
    ! Since (1-(1/3)*beta) is a constant, it can be pulled outside of the
    ! derivative.  The explicit portion of this term becomes:
    !
    ! - (1-(1/3)*beta)/rho_ds
    !   * d [ rho_ds * (a_1)^2 * ( w'^3 / (w'^2)^2 ) * w'x_a' * w'x_b' ] / dz.
    !
    ! The explicit portion of this term is discretized as follows:
    !
    ! The values of w'x_a', w'x_b', w'^2, and a_1 are found on the momentum
    ! levels.  The values of w'^3 are found on the thermodynamic levels.
    ! Additionally, the values of rho_ds_zt are found on the thermodynamic
    ! levels, and the values of invrs_rho_ds_zm are found on the momentum
    ! levels.  The variables w'x_a', w'x_b', w'^2, and a_1 are each interpolated
    ! to the intermediate thermodynamic levels.  The values of the mathematical
    ! expression (called F here) within the dF/dz term are computed on the
    ! thermodynamic levels.  Then the derivative (d/dz) of the expression (F) is
    ! taken over the central momentum level, where it is multiplied by
    ! (1-(1/3)*beta), and by invrs_rho_ds_zm, yielding the desired result.  In
    ! this function, the values of F are as follows:
    !
    ! F = rho_ds_zt * ( a_1(t) )^2 * ( w'^3(t) / ( w'^2(t) )^2 )
    !     * w'x_a'(t) * w'x_b'(t);
    !
    ! where the timestep index (t) stands for the index of the current timestep.
    !
    !
    ! =a1p1=======wp2p1=======wpxapp1=======wpxbpp1========================= m(k+1)
    !
    ! -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3p1-rho_ds_ztp1- t(k+1)
    !
    ! =a1=========wp2=========wpxap=========wpxbp===dF/dz===invrs_rho_ds_zm= m(k)
    !
    ! -a1(interp)-wp2(interp)-wpxap(interp)-wpxbp(interp)-wp3---rho_ds_zt--- t(k)
    !
    ! =a1m1=======wp2m1=======wpxapm1=======wpxbpm1========================= m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, & ! Constant(s)
        one_third

    use model_flags, only:  &
        l_standard_term_ta

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max

    ! Input variables
    real( kind = core_rknd ), intent(in) :: & 
      wp2_ztp1,        & ! w'^2 interpolated to thermo. level (k+1)   [m^2/s^2]
      wp2_zt,          & ! w'^2 interpolated to thermo. level (k)     [m^2/s^2]
      wp3_on_wp2_ztp1, & ! Smoothed w'^3 / w'^2 on thermo. level (k+1)[m^2/s^2]
      wp3_on_wp2_zt,   & ! Smoothed w'^3 / w'^2 on thermo. level (k)  [m^2/s^2]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ mome. lev (k)   [m^3/kg]
      a1_ztp1,         & ! a_1 interpolated to thermo. level (k+1)    [-]
      a1,              & ! a_1(k)                                     [-]
      a1_zt,           & ! a_1 interpolated to thermo. level (k)      [-]
      wpxbp_ztp1,      & ! w'x_b' interpolated to thermo. level (k+1) [m/s {x_bm units}]
      wpxbp_zt,        & ! w'x_b' interpolated to thermo. level (k)   [m/s {x_bm units}]
      wpxap_ztp1,      & ! w'x_a' interpolated to thermo. level (k+1) [m/s {x_am units}]
      wpxap_zt,        & ! w'x_a' interpolated to thermo. level (k)   [m/s {x_am units}]
      invrs_dzm,       & ! Inverse of grid spacing                    [1/m]
      beta               ! Model parameter                            [-]

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    if ( l_standard_term_ta ) then

      ! The turbulent advection term is discretized normally, in accordance
      ! with the model equations found in the documentation and the description
      ! listed above.

      rhs  &
      = - ( one - one_third * beta )  &
          * invrs_rho_ds_zm  &
            * invrs_dzm  &
              * (   rho_ds_ztp1 * a1_ztp1**2  &
                    * wp3_on_wp2_ztp1 / wp2_ztp1 &
                    * wpxap_ztp1 * wpxbp_ztp1  &
                  - rho_ds_zt * a1_zt**2  &
                    * wp3_on_wp2_zt / wp2_zt  &
                    * wpxap_zt * wpxbp_zt  &
                )

    else

      ! Brian tried a new discretization for the turbulent advection term, for
      ! which the explicit portion of the term is:
      !  - (1/rho_ds)
      !    * d [ rho_ds * (a_1)^2 * (1-(1/3)*beta) * ( w'^3 / (w'^2)^2 )
      !          * w'x_a' * w'x_b' ] / dz.
      ! In order to help stabilize x_a'x_b', (a_1)^2 has been pulled outside
      ! the derivative.

      rhs & 
      = - ( one - one_third * beta )  &
          * invrs_rho_ds_zm * a1**2  &
            * invrs_dzm  &
              * (   rho_ds_ztp1  &
                    * wp3_on_wp2_ztp1 / wp2_ztp1 &
                    * wpxap_ztp1 * wpxbp_ztp1  &
                  - rho_ds_zt  &
                    * wp3_on_wp2_zt / wp2_zt &
                    * wpxap_zt * wpxbp_zt  &
                )

      ! End of Brian's a1 change.  14 Feb 2008.

    endif


    return
  end function term_ta_rhs

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
        two  ! Constant(s)

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
    else if ( gr%zt(k_high) == gr%zt(k)+depth ) then
      ! Current altitude + depth falls exactly on another grid level.
      ! In this case, no interpolation is necessary.
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
    else if ( gr%zt(k_low) == gr%zt(k)-depth ) then
      ! Current altitude - depth falls exactly on another grid level.
      ! In this case, no interpolation is necessary.
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

    type(pdf_parameter), dimension(nz), intent(in) :: &
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
    else where
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
