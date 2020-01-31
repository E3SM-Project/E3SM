!-----------------------------------------------------------------------
! $Id: advance_xm_wpxp_module.F90 6146 2013-04-05 18:02:22Z raut@uwm.edu $
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
             wpxp_term_ta_lhs, & 
             wpxp_term_tp_lhs, & 
             wpxp_terms_ac_pr2_lhs, & 
             wpxp_term_pr1_lhs, & 
             wpxp_terms_bp_pr3_rhs, &
             xm_correction_wpxp_cl, &
             damp_coefficient

  ! Parameter Constants
  integer, parameter, private :: & 
    nsub = 2, & ! Number of subdiagonals in the LHS matrix
    nsup = 2, & ! Number of superdiagonals in the LHS matrix
    xm_wpxp_thlm = 1, &   ! Named constant for thlm solving
    xm_wpxp_rtm = 2, &    ! Named constant for rtm solving
    xm_wpxp_scalar = 3    ! Named constant for scalar solving

  contains

  !=============================================================================
  subroutine advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                              Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, &
                              tau_zm, Skw_zm, rtpthvp, rtm_forcing, &
                              wprtp_forcing, rtm_ref, thlpthvp, &
                              thlm_forcing, wpthlp_forcing, thlm_ref, &
                              rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                              invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
                              w1_zm, w2_zm, varnce_w1_zm, varnce_w2_zm, &
                              mixt_frac_zm, l_implemented, &
                              sclrpthvp, sclrm_forcing, sclrp2, &
                              rtm, wprtp, thlm, wpthlp, &
                              err_code, &
                              sclrm, wpsclrp )

    ! Description:
    ! Advance the mean and flux terms by one timestep.

    ! References:
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
        thl_tol_mfl, &
        rt_tol_mfl, &
        max_mag_correlation, &
        one, &
        one_half, &
        zero, &
        zero_threshold

    use parameters_model, only: & 
        sclr_dim, &  ! Variable(s)
        sclr_tol

    use grid_class, only: & 
        gr ! Variable(s)

    use grid_class, only: &
        zm2zt, & ! Procedure(s)
        zt2zm

    use model_flags, only: &
        l_clip_semi_implicit ! Variable(s)

    use mono_flux_limiter, only: &
        calc_turb_adv_range ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use error_code, only:  & 
        clubb_at_least_debug_level, & ! Procedure(s)
        reportError, &
        fatal_error

    use error_code, only:  & 
      clubb_var_out_of_range ! Constant(s)

    use stats_type, only: &
        stat_begin_update, & ! Procedure(s)
        stat_end_update, &
        stat_update_var

    use stats_variables, only: & 
        zt, &
        zm, &
        irtm_matrix_condt_num, &  ! Variables
        ithlm_matrix_condt_num, &
        irtm_sdmp, ithlm_sdmp, & 
        l_stats_samp, &
        iC7_Skw_fnc, &
        iC6rt_Skw_fnc, &
        iC6thl_Skw_fnc, &
        l_stats_samp

    use sponge_layer_damping, only: &
        rtm_sponge_damp_settings, &
        thlm_sponge_damp_settings, &
        rtm_sponge_damp_profile, &
        thlm_sponge_damp_profile, &
        sponge_damp_xm ! Procedure(s)

    implicit none

    ! External
    intrinsic :: exp, sqrt

    ! Parameter Constants
    logical, parameter :: &
      l_iter = .true. ! True when the means and fluxes are prognosed

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      sigma_sqd_w,     & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,           & ! w wind component on momentum levels      [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      Lscale,          & ! Turbulent mixing length                  [m]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on momentum levels    [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels     [m/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      tau_zm,          & ! Time-scale tau on momentum levels        [s]
      Skw_zm,          & ! Skewness of w on momentum levels         [-]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)       [(kg/kg)/s^2]
      rtm_ref,         & ! rtm for nudging                          [kg/kg]
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
      w1_zm,           & ! Mean w (1st PDF component)                   [m/s]
      w2_zm,           & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w1_zm,    & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w2_zm,    & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm       ! Weight of 1st PDF component (Sk_w dependent) [-]

    logical, intent(in) ::  & 
      l_implemented      ! Flag for CLUBB being implemented in a larger model.


    ! Additional variables for passive scalars
    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) ::  & 
      sclrpthvp, sclrm_forcing,  & !                           [Units vary]
      sclrp2                       ! For clipping Vince Larson [Units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  & 
      rtm,       & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,     & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,      & ! th_l (liquid water potential temperature) [K]
      wpthlp       ! w'th_l'                                   [K m/s]

    integer, intent(inout) :: err_code ! Error code for the model's status

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Local variables
    real( kind = core_rknd ), dimension(nsup+nsub+1,2*gr%nz) :: & 
      lhs  ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      C6rt_Skw_fnc, C6thl_Skw_fnc, C7_Skw_fnc

    ! Eddy Diffusion for wpthlp and wprtp.
    real( kind = core_rknd ), dimension(gr%nz) :: Kw6   ! wpxp eddy diff. [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      a1,     & ! a_1 (momentum levels); See eqn. 24 in `Equations for CLUBB' [-]
      a1_zt     ! a_1 interpolated to thermodynamic levels                    [-]

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

    real( kind = core_rknd ), dimension(gr%nz) :: dummy_1d ! Unreferenced array

    real( kind = core_rknd ), allocatable, dimension(:,:) :: & 
      rhs,     &! Right-hand sides of band diag. matrix. (LAPACK)
      solution  ! solution vectors of band diag. matrix. (LAPACK)

    ! Constant parameters as a function of Skw.

    integer :: &
      nrhs, &          ! Number of RHS vectors
      err_code_xm_wpxp ! Error code

    real( kind = core_rknd ) :: rcond

    ! Indices
    integer :: i

    !---------------------------------------------------------------------------

    ! ----- Begin Code -----
    if ( l_clip_semi_implicit ) then
      nrhs = 1
    else
      nrhs = 2+sclr_dim
    endif

    ! Allocate rhs and solution vector
    allocate( rhs(2*gr%nz,nrhs) )
    allocate( solution(2*gr%nz,nrhs) )

    ! This is initialized solely for the purpose of avoiding a compiler
    ! warning about uninitialized variables.
    dummy_1d = zero

    ! Compute C6 and C7 as a function of Skw
    ! The if...then is just here to save compute time
    if ( C6rt /= C6rtb ) then
      C6rt_Skw_fnc(1:gr%nz) = C6rtb + (C6rt-C6rtb) & 
        *EXP( -one_half * (Skw_zm(1:gr%nz)/C6rtc)**2 )
    else
      C6rt_Skw_fnc(1:gr%nz) = C6rtb
    endif

    if ( C6thl /= C6thlb ) then
      C6thl_Skw_fnc(1:gr%nz) = C6thlb + (C6thl-C6thlb) & 
        *EXP( -one_half * (Skw_zm(1:gr%nz)/C6thlc)**2 )
    else
      C6thl_Skw_fnc(1:gr%nz) = C6thlb
    endif

    if ( C7 /= C7b ) then
      C7_Skw_fnc(1:gr%nz) = C7b + (C7-C7b) & 
        *EXP( -one_half * (Skw_zm(1:gr%nz)/C7c)**2 )
    else
      C7_Skw_fnc(1:gr%nz) = C7b
    endif

    ! Damp C6 and C7 as a function of Lscale in stably stratified regions
    C7_Skw_fnc = damp_coefficient( C7, C7_Skw_fnc, &
                                   C7_Lscale0, wpxp_L_thresh, Lscale )
    C6rt_Skw_fnc = damp_coefficient( C6rt, C6rt_Skw_fnc, &
                                     C6rt_Lscale0, wpxp_L_thresh, Lscale )
    C6thl_Skw_fnc = damp_coefficient( C6thl, C6thl_Skw_fnc, &
                                      C6thl_Lscale0, wpxp_L_thresh, Lscale )

    !        C6rt_Skw_fnc = C6rt
    !        C6thl_Skw_fnc = C6thl
    !        C7_Skw_fnc = C7

    if ( l_stats_samp ) then

      call stat_update_var( iC7_Skw_fnc, C7_Skw_fnc, zm )
      call stat_update_var( iC6rt_Skw_fnc, C6rt_Skw_fnc, zm )
      call stat_update_var( iC6thl_Skw_fnc, C6thl_Skw_fnc, zm )

    end if

    if ( clubb_at_least_debug_level( 2 ) ) then
      ! Assertion check for C7_Skw_fnc
      if ( any( C7_Skw_fnc(:) > one ) .or. any( C7_Skw_fnc(:) < zero ) ) then
        write(fstderr,*) "The C7_Skw_fnc variable is outside the valid range"
        err_code = clubb_var_out_of_range
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
    call calc_turb_adv_range( dt, w1_zm, w2_zm, varnce_w1_zm, varnce_w2_zm, & ! In
                              mixt_frac_zm, &  ! In
                              low_lev_effect, high_lev_effect ) ! Out


    ! Define a_1 (located on momentum levels).
    ! It is a variable that is a function of sigma_sqd_w (where sigma_sqd_w is
    ! located on momentum levels).
    a1(1:gr%nz) = one / ( one - sigma_sqd_w(1:gr%nz) )

    ! Interpolate a_1 from momentum levels to thermodynamic levels.  This will
    ! be used for the w'x' turbulent advection (ta) term.
    a1_zt  = max( zm2zt( a1 ), zero_threshold )   ! Positive definite quantity

    ! Setup and decompose matrix for each variable.

    if ( l_clip_semi_implicit ) then

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
      call xm_wpxp_lhs( l_iter, dt, wprtp, a1, a1_zt, wm_zm, wm_zt,  &  ! Intent(in)
                        wp2, wp3_on_wp2, wp3_on_wp2_zt, & ! Intent(in)
                        Kw6, tau_zm, C7_Skw_fnc,  & ! Intent(in)
                        C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! Intent(in)
                        invrs_rho_ds_zm, invrs_rho_ds_zt,  & ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, l_implemented, & ! Intent(in)
                        lhs ) ! Intent(out)

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_rtm, l_iter, dt, rtm, wprtp, &      ! Intent(in)
                        rtm_forcing, wprtp_forcing, C7_Skw_fnc, &   ! Intent(in)
                        rtpthvp, C6rt_Skw_fnc, tau_zm, a1, a1_zt, & ! Intent(in)
                        wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &     ! Intent(in)
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &    ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &           ! Intent(in)
                        rhs(:,1) )                                  ! Intent(out)

      ! Solve r_t / w'r_t'
      if ( l_stats_samp .and. irtm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                            lhs, rhs, &                 ! Intent(inout)
                            solution, err_code_xm_wpxp, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution, err_code_xm_wpxp ) ! Intent(out)
      endif

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "Mean total water & total water flux LU decomp. failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
             rtm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             rt_tol**2, rt_tol, rcond, &            ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,1), &        ! Intent(in)
             rtm, rt_tol_mfl, wprtp, & ! Intent(inout)
             err_code_xm_wpxp )  ! Intent(out)

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "rtm monotonic flux limiter:  tridag failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if


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
      call xm_wpxp_lhs( l_iter, dt, wpthlp, a1, a1_zt, wm_zm, wm_zt, & ! Intent(in)
                        wp2, wp3_on_wp2, wp3_on_wp2_zt, & !  Intent(in)
                        Kw6, tau_zm, C7_Skw_fnc, & ! Intent(in)
                        C6thl_Skw_fnc, rho_ds_zm, rho_ds_zt, & ! Intent(in)
                        invrs_rho_ds_zm, invrs_rho_ds_zt, & ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, l_implemented, & ! Intent(in)
                        lhs ) ! Intent(out)

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_thlm, l_iter, dt, thlm, wpthlp, &     ! Intent(in)
                        thlm_forcing, wpthlp_forcing, C7_Skw_fnc, &   ! Intent(in)
                        thlpthvp, C6thl_Skw_fnc, tau_zm, a1, a1_zt, & ! Intent(in)
                        wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &       ! Intent(in)
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &      ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &             ! Intent(in)
                        rhs(:,1) )                                    ! Intent(out)

      ! Solve for th_l / w'th_l'
      if ( l_stats_samp .and. ithlm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                            lhs, rhs, &                 ! Intent(inout)
                            solution, err_code_xm_wpxp, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution, err_code_xm_wpxp ) ! Intent(out)
      endif

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "Liquid pot. temp & thetal flux LU decomp. failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_thlm, dt, wp2, thlp2, wm_zt,  & ! Intent(in)
             thlm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &     ! Intent(in)
             thl_tol**2, thl_tol, rcond, &           ! Intent(in)
             low_lev_effect, high_lev_effect, &      ! Intent(in)
             l_implemented, solution(:,1),  &        ! Intent(in)
             thlm, thl_tol_mfl, wpthlp, & ! Intent(inout)
             err_code_xm_wpxp ) ! Intent(out)

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "thlm monotonic flux limiter:  tridag failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

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

        ! Compute the implicit portion of the sclr and w'sclr' equations.
        ! Build the left-hand side matrix.
        call xm_wpxp_lhs( l_iter, dt, wpsclrp(:,i), a1, a1_zt, wm_zm, wm_zt, & ! Intent(in)
                          wp2, wp3_on_wp2, wp3_on_wp2_zt, & !  Intent(in)
                          Kw6, tau_zm, C7_Skw_fnc, &  ! Intent(in)
                          C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt,  &  ! Intent(in)
                          invrs_rho_ds_zm, invrs_rho_ds_zt,  &  ! Intent(in)
                          wpxp_upper_lim, wpxp_lower_lim, l_implemented, & ! Intent(in)
                          lhs ) ! Intent(out)

        ! Compute the explicit portion of the sclrm and w'sclr' equations.
        ! Build the right-hand side vector.
        call xm_wpxp_rhs( xm_wpxp_scalar, l_iter, dt, sclrm(:,i), wpsclrp(:,i), & ! Intent(in)
                          sclrm_forcing(:,i), dummy_1d, C7_Skw_fnc, &             ! Intent(in)
                          sclrpthvp(:,i), C6rt_Skw_fnc, tau_zm, a1, a1_zt, &      ! Intent(in)
                          wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &                 ! Intent(in)
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &                ! Intent(in)
                          wpxp_upper_lim, wpxp_lower_lim, &                       ! Intent(in)
                          rhs(:,1) )                                              ! Intent(out)

        ! Solve for sclrm / w'sclr'
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution, err_code_xm_wpxp ) ! Intent(out)

        if ( fatal_error( err_code_xm_wpxp ) ) then
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Passive scalar # ", i, " LU decomp. failed."
            call reportError( err_code_xm_wpxp )
          end if

          ! Overwrite the current error status with the new fatal error
          err_code = err_code_xm_wpxp

        end if

        call xm_wpxp_clipping_and_stats &
             ( xm_wpxp_scalar, dt, wp2, sclrp2(:,i),  & ! Intent(in)
               wm_zt, sclrm_forcing(:,i),  &            ! Intent(in)
               rho_ds_zm, rho_ds_zt, &                  ! Intent(in)
               invrs_rho_ds_zm, invrs_rho_ds_zt, &      ! Intent(in)
               sclr_tol(i)**2, sclr_tol(i), rcond, &    ! Intent(in)
               low_lev_effect, high_lev_effect, &       ! Intent(in)
               l_implemented, solution(:,1),  &         ! Intent(in)
               sclrm(:,i), sclr_tol(i), wpsclrp(:,i), & ! Intent(inout)
               err_code_xm_wpxp ) ! Intent(out)

       if ( fatal_error( err_code_xm_wpxp ) ) then
         if ( clubb_at_least_debug_level( 1 ) ) then
           write(fstderr,*) "sclrm # ", i, "monotonic flux limiter: tridag failed"
           call reportError( err_code_xm_wpxp )
         end if

         ! Overwrite the current error status with the new fatal error
         err_code = err_code_xm_wpxp

       end if

      enddo ! passive scalars

    else ! Simple case, where l_clip_semi_implicit is false

      ! Create the lhs once
      call xm_wpxp_lhs( l_iter, dt, dummy_1d, a1, a1_zt, wm_zm, wm_zt, & ! Intent(in)
                        wp2, wp3_on_wp2, wp3_on_wp2_zt, & ! Intent(in)
                        Kw6, tau_zm, C7_Skw_fnc, & ! Intent(in)
                        C6rt_Skw_fnc, rho_ds_zm, rho_ds_zt,  & ! Intent(in)
                        invrs_rho_ds_zm, invrs_rho_ds_zt,  & ! Intent(in)
                        dummy_1d, dummy_1d, l_implemented,  & ! Intent(in)
                        lhs ) ! Intent(out)

      ! Compute the explicit portion of the r_t and w'r_t' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_rtm, l_iter, dt, rtm, wprtp, &      ! Intent(in)
                        rtm_forcing, wprtp_forcing, C7_Skw_fnc, &   ! Intent(in)
                        rtpthvp, C6rt_Skw_fnc, tau_zm, a1, a1_zt, & ! Intent(in)
                        wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &     ! Intent(in)
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &    ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &           ! Intent(in)
                        rhs(:,1) )                                  ! Intent(out)

      ! Compute the explicit portion of the th_l and w'th_l' equations.
      ! Build the right-hand side vector.
      call xm_wpxp_rhs( xm_wpxp_thlm, l_iter, dt, thlm, wpthlp, &     ! Intent(in)
                        thlm_forcing, wpthlp_forcing, C7_Skw_fnc, &   ! Intent(in)
                        thlpthvp, C6thl_Skw_fnc, tau_zm, a1, a1_zt, & ! Intent(in)
                        wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &       ! Intent(in)
                        rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &      ! Intent(in)
                        wpxp_upper_lim, wpxp_lower_lim, &             ! Intent(in)
                        rhs(:,2) )                                    ! Intent(out)

! ---> h1g, 2010-06-15
! scalar transport, e.g, droplet and ice number concentration
! are handled in  " advance_sclrm_Nd_module.F90 "
#ifdef GFDL
      do i = 1, 0, 1
#else
      do i = 1, sclr_dim, 1
#endif
! <--- h1g, 2010-06-15

        call xm_wpxp_rhs( xm_wpxp_scalar, l_iter, dt, sclrm(:,i), wpsclrp(:,i), & ! Intent(in)
                          sclrm_forcing(:,i), dummy_1d, C7_Skw_fnc, &             ! Intent(in)
                          sclrpthvp(:,i), C6rt_Skw_fnc, tau_zm, a1, a1_zt, &      ! Intent(in)
                          wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &                 ! Intent(in)
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &                ! Intent(in)
                          wpxp_upper_lim, wpxp_lower_lim, &                       ! Intent(in)
                          rhs(:,2+i) )                                              ! Intent(out)

      enddo

      ! Solve for all fields
      if ( l_stats_samp .and. ithlm_matrix_condt_num + irtm_matrix_condt_num > 0 ) then
        call xm_wpxp_solve( nrhs, &                     ! Intent(in)
                            lhs, rhs, &                 ! Intent(inout)
                            solution, err_code_xm_wpxp, rcond ) ! Intent(out)
      else
        call xm_wpxp_solve( nrhs, &              ! Intent(in)
                            lhs, rhs, &          ! Intent(inout)
                            solution, err_code_xm_wpxp ) ! Intent(out)
      endif

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "xm_wpxp matrix LU decomp. failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_rtm, dt, wp2, rtp2, wm_zt,  &  ! Intent(in)
             rtm_forcing, rho_ds_zm, rho_ds_zt, &   ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             rt_tol**2, rt_tol, rcond, &            ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,1),  &       ! Intent(in)
             rtm, rt_tol_mfl, wprtp,  & ! Intent(inout)
             err_code_xm_wpxp ) ! Intent(out)

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "rtm monotonic flux limiter:  tridag failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

      call xm_wpxp_clipping_and_stats &
           ( xm_wpxp_thlm, dt, wp2, thlp2, wm_zt, & ! Intent(in)
             thlm_forcing, rho_ds_zm, rho_ds_zt, &  ! Intent(in)
             invrs_rho_ds_zm, invrs_rho_ds_zt, &    ! Intent(in)
             thl_tol**2, thl_tol, rcond, &          ! Intent(in)
             low_lev_effect, high_lev_effect, &     ! Intent(in)
             l_implemented, solution(:,2),  &       ! Intent(in)
             thlm, thl_tol_mfl, wpthlp, & ! Intent(inout)
             err_code_xm_wpxp ) ! Intent(out)

      if ( fatal_error( err_code_xm_wpxp ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,'(a)') "thlm monotonic flux limiter:  tridag failed"
          call reportError( err_code_xm_wpxp )
        end if

        ! Overwrite the current error status with the new fatal error
        err_code = err_code_xm_wpxp

      end if

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
               sclrm(:,i), sclr_tol(i), wpsclrp(:,i), & ! Intent(inout)
               err_code_xm_wpxp ) ! Intent(out)

        if ( fatal_error( err_code_xm_wpxp ) ) then
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "sclrm # ", i, "monotonic flux limiter: tridag failed"
            call reportError( err_code_xm_wpxp )
          end if

          ! Overwrite the current error status with the new fatal error
          err_code = err_code_xm_wpxp

        end if

      end do ! 1..sclr_dim

    end if ! l_clip_semi_implicit

    ! De-allocate memory
    deallocate( rhs, solution )

    ! Error Report
    ! Joshua Fasching Feb 2008
    if ( fatal_error( err_code ) .and. clubb_at_least_debug_level( 1 ) ) then

      write(fstderr,*) "Error in advance_xm_wpxp"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "dt = ", dt
      write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
      write(fstderr,*) "wm_zm = ", wm_zm
      write(fstderr,*) "wm_zt = ", wm_zt
      write(fstderr,*) "wp2 = ", wp2
      write(fstderr,*) "wp3_on_wp2 = ", wp3_on_wp2
      write(fstderr,*) "wp3_on_wp2_zt = ", wp3_on_wp2_zt
      write(fstderr,*) "Kh_zt = ", Kh_zt
      write(fstderr,*) "tau_zm = ", tau_zm
      write(fstderr,*) "Skw_zm = ", Skw_zm
      write(fstderr,*) "rtpthvp = ", rtpthvp
      write(fstderr,*) "rtm_forcing = ", rtm_forcing
      write(fstderr,*) "wprtp_forcing = ", wprtp_forcing
      write(fstderr,*) "rtm_ref = ", rtm_ref
      write(fstderr,*) "thlpthvp = ", thlpthvp
      write(fstderr,*) "thlm_forcing = ", thlm_forcing
      write(fstderr,*) "wpthlp_forcing = ", wpthlp_forcing
      write(fstderr,*) "thlm_ref = ", thlm_ref
      write(fstderr,*) "rho_ds_zm = ", rho_ds_zm
      write(fstderr,*) "rho_ds_zt = ", rho_ds_zt
      write(fstderr,*) "invrs_rho_ds_zm = ", invrs_rho_ds_zm
      write(fstderr,*) "invrs_rho_ds_zt = ", invrs_rho_ds_zt
      write(fstderr,*) "thv_ds_zm = ", thv_ds_zm
      write(fstderr,*) "rtp2 = ", rtp2
      write(fstderr,*) "thlp2 = ", thlp2
      write(fstderr,*) "w1_zm = ", w1_zm
      write(fstderr,*) "w2_zm = ", w2_zm
      write(fstderr,*) "varnce_w1_zm = ", varnce_w1_zm
      write(fstderr,*) "varnce_w2_zm = ", varnce_w2_zm
      write(fstderr,*) "mixt_frac_zm = ", mixt_frac_zm
      write(fstderr,*) "l_implemented = ", l_implemented
 
      if ( sclr_dim > 0 )  then
        write(fstderr,*) "sclrp2 = ", sclrp2
        write(fstderr,*) "sclrpthvp = ", sclrpthvp
        write(fstderr,*) "sclrm_forcing = ", sclrm_forcing
      end if
    
      write(fstderr,*) "Intent(inout)"
 
      write(fstderr,*) "rtm = ", rtm
      write(fstderr,*) "wprtp = ", wprtp
      write(fstderr,*) "thlm = ", thlm
      write(fstderr,*) "wpthlp =", wpthlp

      if ( sclr_dim > 0 )  then
        write(fstderr,*) "sclrm = ", sclrm
        write(fstderr,*) "wpsclrp = ", wpsclrp
      end if

    end if ! Fatal error and debug_level >= 1

    if ( rtm_sponge_damp_settings%l_sponge_damping ) then
      if( l_stats_samp ) then
        call stat_begin_update( irtm_sdmp, rtm / real( dt, kind = core_rknd ), zt )
      end if
      rtm(1:gr%nz) = sponge_damp_xm( dt, rtm_ref(1:gr%nz), rtm(1:gr%nz), &
                                       rtm_sponge_damp_profile )

      if( l_stats_samp ) then
        call stat_end_update( irtm_sdmp, rtm / real( dt, kind = core_rknd ), zt )
      end if
    endif

    if ( thlm_sponge_damp_settings%l_sponge_damping ) then
      if( l_stats_samp ) then
        call stat_begin_update( ithlm_sdmp, thlm / real( dt, kind = core_rknd ), zt )
      end if
      thlm(1:gr%nz) = sponge_damp_xm( dt, thlm_ref(1:gr%nz), thlm(1:gr%nz), &
                                        thlm_sponge_damp_profile )
      if( l_stats_samp ) then
        call stat_end_update( ithlm_sdmp, thlm / real( dt, kind = core_rknd ), zt )
      end if
    endif

    return

  end subroutine advance_xm_wpxp

  !=============================================================================
  subroutine xm_wpxp_lhs( l_iter, dt, wpxp, a1, a1_zt, wm_zm, wm_zt,  &
                          wp2, wp3_on_wp2, wp3_on_wp2_zt, &
                          Kw6, tau_zm, C7_Skw_fnc, &
                          C6x_Skw_fnc, rho_ds_zm, rho_ds_zt,  &
                          invrs_rho_ds_zm, invrs_rho_ds_zt,  &
                          wpxp_upper_lim, wpxp_lower_lim, l_implemented,  &
                          lhs )

    ! Description:
    !   Compute LHS band diagonal matrix for xm and w'x'.
    !   This subroutine computes the implicit portion of
    !   the xm and w'x' equations.

    ! References:
    !   None
    !------------------------------------------------------------------------

    use parameters_tunable, only:  & 
        nu6_vert_res_dep ! Variable(s)

    use grid_class, only:  & 
        gr,  & ! Variable(s)
        zm2zt ! Procedure(s)

    use constants_clubb, only: &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_clip_semi_implicit, & ! Variable(s)
        l_upwind_wpxp_ta

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use diffusion, only:  & 
        diffusion_zm_lhs ! Procedure(s)

    use mean_adv, only: & 
        term_ma_zt_lhs,  & ! Procedure(s)
        term_ma_zm_lhs

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

    use advance_helper_module, only: set_boundary_conditions_lhs ! Procedure(s)


    implicit none

    ! External
    intrinsic :: min, max

    ! Constant parameters
    ! Left-hand side matrix diagonal identifiers for
    ! momentum-level variable, w'x'.
    integer, parameter ::  &
      m_kp1_mdiag = 1, & ! Momentum superdiagonal index for w'x'.
      m_kp1_tdiag = 2, & ! Thermodynamic superdiagonal index for w'x'.
      m_k_mdiag   = 3, & ! Momentum main diagonal index for w'x'.
      m_k_tdiag   = 4, & ! Thermodynamic subdiagonal index for w'x'.
      m_km1_mdiag = 5    ! Momentum subdiagonal index for w'x'.

    ! Left-hand side matrix diagonal identifiers for
    ! thermodynamic-level variable, xm.
    integer, parameter ::  &
      t_kp1_tdiag = 1, & ! Thermodynamic superdiagonal index for xm.
      t_k_mdiag   = 2, & ! Momentum superdiagonal index for xm.
      t_k_tdiag   = 3, & ! Thermodynamic main diagonal index for xm.
      t_km1_mdiag = 4, & ! Momentum subdiagonal index for xm.
      t_km1_tdiag = 5    ! Thermodynamic subdiagonal index for xm.

    ! Input variables
    logical, intent(in) :: l_iter

    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: & 
      wpxp,            & ! w'x' (momentum levels) at timestep (t)    [{xm units} m/s]
      a1,              & ! a_1 (momentum levels)                     [-]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels  [-]
      wm_zm,           & ! w wind component on momentum levels       [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels  [m/s]
      wp2,             & ! w'^2 (momentum levels)                    [m^2/s^2]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on momentum levels     [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels      [m/s]
      Kw6,             & ! Coefficient of eddy diffusivity for w'x'  [m^2/s]
      tau_zm,          & ! Time-scale tau on momentum levels         [s]
      C7_Skw_fnc,      & ! C_7 parameter with Sk_w applied           [-]
      C6x_Skw_fnc,     & ! C_6x parameter with Sk_w applied          [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs.  [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      wpxp_upper_lim,  & ! Keeps correlations from becoming > 1.     [units vary]
      wpxp_lower_lim     ! Keeps correlations from becoming < -1.    [units vary]

    logical, intent(in) ::  & 
      l_implemented ! Flag for CLUBB being implemented in a larger model.

    ! Output Variable
    real( kind = core_rknd ), intent(out), dimension(nsup+nsub+1,2*gr%nz) ::  & 
      lhs ! Implicit contributions to wpxp/xm (band diag. matrix) (LAPACK)

    ! Local Variables

    ! Indices
    integer :: k, kp1, km1
    integer :: k_xm, k_wpxp
    integer :: k_wpxp_low, k_wpxp_high

    real( kind = core_rknd ), dimension(3) :: tmp

    logical :: l_upper_thresh, l_lower_thresh ! flags for clip_semi_imp_lhs


    ! Initialize the left-hand side matrix to 0.
    lhs = zero

    ! The xm loop runs between k = 2 and k = gr%nz.  The value of xm at
    ! level k = 1, which is below the model surface, is simply set equal to the
    ! value of xm at level k = 2 after the solve has been completed.

    do k = 2, gr%nz, 1

      ! Define indices

      km1 = max( k-1, 1 )

      k_xm = 2*k - 1
      ! k_wpxp is 2*k


      !!!!!***** xm *****!!!!!

      ! xm: Left-hand side (implicit xm portion of the code).
      !
      ! Thermodynamic subdiagonal (lhs index: t_km1_tdiag)
      !         [ x xm(k-1,<t+1>) ]
      ! Momentum subdiagonal (lhs index: t_km1_mdiag)
      !         [ x wpxp(k-1,<t+1>) ]
      ! Thermodynamic main diagonal (lhs index: t_k_tdiag)
      !         [ x xm(k,<t+1>) ]
      ! Momentum superdiagonal (lhs index: t_k_mdiag)
      !         [ x wpxp(k,<t+1>) ]
      ! Thermodynamic superdiagonal (lhs index: t_kp1_tdiag)
      !         [ x xm(k+1,<t+1>) ]

      ! LHS mean advection (ma) term.
      if ( .not. l_implemented ) then

        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

      else

        lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) & 
        = lhs((/t_kp1_tdiag,t_k_tdiag,t_km1_tdiag/),k_xm) + zero

      endif

      ! LHS turbulent advection (ta) term.
      lhs((/t_k_mdiag,t_km1_mdiag/),k_xm) & 
      = lhs((/t_k_mdiag,t_km1_mdiag/),k_xm) & 
      + xm_term_ta_lhs( rho_ds_zm(k), rho_ds_zm(km1), &
                        invrs_rho_ds_zt(k), gr%invrs_dzt(k) )

      ! LHS time tendency.
      lhs(t_k_tdiag,k_xm) & 
      = lhs(t_k_tdiag,k_xm) + one / real( dt, kind = core_rknd )

      if (l_stats_samp) then

        ! Statistics: implicit contributions for rtm or thlm.

        if ( irtm_ma > 0 .or. ithlm_ma > 0 ) then
          if ( .not. l_implemented ) then
            tmp(1:3) =  & 
            + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )
            ztscr01(k) = - tmp(3)
            ztscr02(k) = - tmp(2)
            ztscr03(k) = - tmp(1)
          else
            ztscr01(k) = zero
            ztscr02(k) = zero
            ztscr03(k) = zero
          endif
        endif

        if ( irtm_ta > 0 .or. ithlm_ta > 0 ) then
          tmp(1:2) = & 
          + xm_term_ta_lhs( rho_ds_zm(k), rho_ds_zm(km1), &
                            invrs_rho_ds_zt(k), gr%invrs_dzt(k) )
          ztscr04(k) = - tmp(2)
          ztscr05(k) = - tmp(1)
        endif

      endif

    enddo ! xm loop: 2..gr%nz


    ! The wpxp loop runs between k = 2 and k = gr%nz-1.  The value of wpxp
    ! is set to specified values at both the lowest level, k = 1, and the
    ! highest level, k = gr%nz.

    do k = 2, gr%nz-1, 1

      ! Define indices

      kp1 = min( k+1, gr%nz )
      km1 = max( k-1, 1 )

      ! k_xm is 2*k - 1
      k_wpxp = 2*k


      !!!!!***** w'x' *****!!!!!

      ! w'x': Left-hand side (implicit w'x' portion of the code).
      !
      ! Momentum subdiagonal (lhs index: m_km1_mdiag)
      !         [ x wpxp(k-1,<t+1>) ]
      ! Thermodynamic subdiagonal (lhs index: m_k_tdiag)
      !         [ x xm(k,<t+1>) ]
      ! Momentum main diagonal (lhs index: m_k_mdiag)
      !         [ x wpxp(k,<t+1>) ]
      ! Thermodynamic superdiagonal (lhs index: m_kp1_tdiag)
      !         [ x xm(k+1,<t+1>) ]
      ! Momentum superdiagonal (lhs index: m_kp1_mdiag)
      !         [ x wpxp(k+1,<t+1>) ]

      ! LHS mean advection (ma) term.
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )

      ! LHS turbulent advection (ta) term.
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      !        The weight of the implicit portion of this term is controlled
      !        by the factor gamma_over_implicit_ts (abbreviated "gamma" in the
      !        the equation in order to balance a weight that is not equal to 1,
      !        such that:
      !             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
      !        where X is the variable that is being solved for in a predictive
      !        equation (w'x' in this case), y(t) is the linearized portion of
      !        the term that gets treated implicitly, and RHS is the portion of
      !        the term that is always treated explicitly (in the case of the
      !        w'x' turbulent advection term, RHS = 0).  A weight of greater
      !        than 1 can be applied to make the term more numerically stable.
      if ( .not. l_upwind_wpxp_ta ) then
        lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  & 
        = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  &
        + gamma_over_implicit_ts  &
        * wpxp_term_ta_lhs( a1_zt(kp1), a1_zt(k),  & 
                            wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                            rho_ds_zt(kp1), rho_ds_zt(k),  &
                            invrs_rho_ds_zm(k),  &
                            gr%invrs_dzm(k), k )
      else
        lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  & 
        = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp)  &
        + gamma_over_implicit_ts  &
        * wpxp_term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                   wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                   gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                   invrs_rho_ds_zm(k), &
                                   rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1) )
      end if

      ! LHS turbulent production (tp) term.
      lhs((/m_kp1_tdiag,m_k_tdiag/),k_wpxp) & 
      = lhs((/m_kp1_tdiag,m_k_tdiag/),k_wpxp) & 
      + wpxp_term_tp_lhs( wp2(k), gr%invrs_dzm(k) )

      ! LHS accumulation (ac) term and pressure term 2 (pr2).
      lhs(m_k_mdiag,k_wpxp) & 
      = lhs(m_k_mdiag,k_wpxp) & 
      + wpxp_terms_ac_pr2_lhs( C7_Skw_fnc(k),  & 
                               wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k) )

      ! LHS pressure term 1 (pr1).
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs(m_k_mdiag,k_wpxp)  & 
      = lhs(m_k_mdiag,k_wpxp)  &
      + gamma_over_implicit_ts  & 
      * wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )

      ! LHS eddy diffusion term: dissipation term 1 (dp1).
      lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      = lhs((/m_kp1_mdiag,m_k_mdiag,m_km1_mdiag/),k_wpxp) & 
      + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6_vert_res_dep, & 
                          gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                          gr%invrs_dzm(k), k )

      ! LHS time tendency.
      if ( l_iter ) then
        lhs(m_k_mdiag,k_wpxp) &
        = lhs(m_k_mdiag,k_wpxp) + one / real(dt, kind = core_rknd)
      endif

      ! LHS portion of semi-implicit clipping term.
      if ( l_clip_semi_implicit ) then
        l_upper_thresh = .true.
        l_lower_thresh = .true.

        lhs(m_k_mdiag,k_wpxp) & 
        = lhs(m_k_mdiag,k_wpxp) & 
        + clip_semi_imp_lhs( dt, wpxp(k),  & 
                             l_upper_thresh, wpxp_upper_lim(k),  & 
                             l_lower_thresh, wpxp_lower_lim(k) )

      endif

      if ( l_stats_samp ) then

        ! Statistics: implicit contributions for wprtp or wpthlp.

        if ( iwprtp_ma > 0 .or. iwpthlp_ma > 0 ) then
          tmp(1:3) = & 
          + term_ma_zm_lhs( wm_zm(k), gr%invrs_dzm(k), k )
          zmscr01(k) = - tmp(3)
          zmscr02(k) = - tmp(2)
          zmscr03(k) = - tmp(1)
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( iwprtp_ta > 0 .or. iwpthlp_ta > 0 ) then
          if ( .not. l_upwind_wpxp_ta ) then
            tmp(1:3)  &
            = gamma_over_implicit_ts  &
            * wpxp_term_ta_lhs( a1_zt(kp1), a1_zt(k),  &
                                wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                                rho_ds_zt(kp1), rho_ds_zt(k),  &
                                invrs_rho_ds_zm(k),  &
                                gr%invrs_dzm(k), k )
          else
            tmp(1:3)  &
            = gamma_over_implicit_ts  &
            * wpxp_term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                       wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                       gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                       invrs_rho_ds_zm(k), &
                                       rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1) )
          end if

          zmscr04(k) = - tmp(3)
          zmscr05(k) = - tmp(2)
          zmscr06(k) = - tmp(1)
        endif

        if ( iwprtp_tp > 0 .or. iwpthlp_tp > 0 ) then
          tmp(1:2) = & 
          + wpxp_term_tp_lhs( wp2(k), gr%invrs_dzm(k) )
          zmscr07(k) = - tmp(2)
          zmscr08(k) = - tmp(1)
        endif

        ! Note:  To find the contribution of w'x' term ac, substitute 0 for the
        !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
        if ( iwprtp_ac > 0 .or. iwpthlp_ac > 0 ) then
          zmscr09(k) =  & 
          - wpxp_terms_ac_pr2_lhs( zero, & 
                                   wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k) )
        endif

        ! Note:  An "over-implicit" weighted time step is applied to this term.
        !        A weighting factor of greater than 1 may be used to make the
        !        term more numerically stable (see note above for LHS turbulent
        !        advection (ta) term).
        if ( iwprtp_pr1 > 0 .or. iwpthlp_pr1 > 0 ) then
          zmscr10(k)  &
          = - gamma_over_implicit_ts  &
            * wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )
        endif

        ! Note:  To find the contribution of w'x' term pr2, add 1 to the
        !        C_7 skewness function input to function wpxp_terms_ac_pr2_lhs.
        if ( iwprtp_pr2 > 0 .or. iwpthlp_pr2 > 0 ) then
          zmscr11(k) = & 
          - wpxp_terms_ac_pr2_lhs( (one+C7_Skw_fnc(k)), & 
                                   wm_zt(kp1), wm_zt(k), gr%invrs_dzm(k) )
        endif

        if ( iwprtp_dp1 > 0 .or. iwpthlp_dp1 > 0 ) then
          tmp(1:3) = & 
          + diffusion_zm_lhs( Kw6(k), Kw6(kp1), nu6_vert_res_dep, & 
                              gr%invrs_dzt(kp1), gr%invrs_dzt(k), &
                              gr%invrs_dzm(k), k )
          zmscr12(k) = - tmp(3)
          zmscr13(k) = - tmp(2)
          zmscr14(k) = - tmp(1)
        endif

        if ( l_clip_semi_implicit ) then
          if ( iwprtp_sicl > 0 .or. iwpthlp_sicl > 0 ) then
            l_upper_thresh = .true.
            l_lower_thresh = .true.
            zmscr15(k) = & 
            - clip_semi_imp_lhs( dt, wpxp(k),  & 
                                 l_upper_thresh, wpxp_upper_lim(k),  & 
                                 l_lower_thresh, wpxp_lower_lim(k) )
          endif
        endif

      endif

    enddo ! wpxp loop: 2..gr%nz-1


    ! Boundary conditions

    ! The turbulent flux (wpxp) use fixed-point boundary conditions at both the
    ! upper and lower boundaries.  Therefore, anything set in the wpxp loop
    ! at both the upper and lower boundaries would be overwritten here.
    ! However, the wpxp loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side and with
    ! values of 0 at all other diagonals on the left-hand side will preserve the
    ! right-hand side value at that level.  The value of xm at level k = 1,
    ! which is below the model surface, is preserved and then overwritten to
    ! match the new value of xm at level k = 2.
    !
    !   xm(1)  wpxp(1) ... wpxp(nzmax)
    ! [  0.0     0.0         0.0    ]
    ! [  0.0     0.0         0.0    ]
    ! [  1.0     1.0   ...   1.0    ]
    ! [  0.0     0.0         0.0    ]
    ! [  0.0     0.0         0.0    ]

    ! Lower boundary
    k      = 1
    k_xm   = 2*k - 1
    k_wpxp_low = 2*k

    ! Upper boundary
    k      = gr%nz
    !k_xm is 2*k - 1
    k_wpxp_high = 2*k

    call set_boundary_conditions_lhs( m_k_mdiag, k_wpxp_low, k_wpxp_high, lhs, &
                                  t_k_tdiag, k_xm)

    return

  end subroutine xm_wpxp_lhs

  !=============================================================================
  subroutine xm_wpxp_rhs( solve_type, l_iter, dt, xm, wpxp, & 
                          xm_forcing, wpxp_forcing, C7_Skw_fnc, &
                          xpthvp, C6x_Skw_fnc, tau_zm, a1, a1_zt, &
                          wp3_on_wp2, wp3_on_wp2_zt, rho_ds_zt, &
                          rho_ds_zm, invrs_rho_ds_zm, thv_ds_zm, &
                          wpxp_upper_lim, wpxp_lower_lim, &
                          rhs )

    ! Description:
    ! Compute RHS vector for xm and w'x'.
    ! This subroutine computes the explicit portion of
    ! the xm and w'x' equations.

    ! References:
    !------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use constants_clubb, only:  &
        gamma_over_implicit_ts, & ! Constant(s)
        one, &
        zero

    use model_flags, only: &
        l_clip_semi_implicit, & ! Variable(s)
        l_upwind_wpxp_ta

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use clip_semi_implicit, only: & 
        clip_semi_imp_rhs ! Procedure(s)

    use stats_type, only: & 
        stat_update_var_pt, & 
        stat_begin_update_pt

    use stats_variables, only: & 
        zt, & ! Variable(s)
        zm, & 
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
        l_stats_samp

    use advance_helper_module, only: set_boundary_conditions_rhs

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      solve_type  ! Variables being solved for.

    logical, intent(in) :: l_iter

    real(kind=time_precision), intent(in) ::  & 
      dt                 ! Timestep                                  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      xm,              & ! xm (thermodynamic levels)                 [{xm units}]
      wpxp,            & ! <w'x'> (momentum levels)                  [{xm units} m/s]
      xm_forcing,      & ! xm forcings (thermodynamic levels)        [{xm units}/s]
      wpxp_forcing,    & ! <w'x'> forcing (momentum levels)          [{xm units} m/s^2]
      C7_Skw_fnc,      & ! C_7 parameter with Sk_w applied           [-]
      xpthvp,          & ! x'th_v' (momentum levels)                 [{xm units} K]
      C6x_Skw_fnc,     & ! C_6x parameter with Sk_w applied          [-]
      tau_zm,          & ! Time-scale tau on momentum levels         [s]
      a1_zt,           & ! a_1 interpolated to thermodynamic levels  [-]
      a1,              & ! a_1                                       [-]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on moment. levels      [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels      [m/s]
      rho_ds_zt,       & ! Dry, static density on thermo.  levels    [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on moment.  levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      wpxp_upper_lim,  & ! Keeps correlations from becoming > 1.     [units vary]
      wpxp_lower_lim     ! Keeps correlations from becoming < -1.    [units vary]

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
    real( kind = core_rknd ), dimension(3) :: lhs_fnc_output

    ! Indices
    integer :: k, km1, kp1, k_xm, k_wpxp, k_xm_low, k_wpxp_low, k_wpxp_high


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
    case default    ! this includes the sclrm case
      ixm_f      = 0
      iwpxp_bp   = 0
      iwpxp_pr3  = 0
      iwpxp_f    = 0
      iwpxp_sicl = 0
      iwpxp_ta   = 0
      iwpxp_pr1  = 0
    end select


    ! Initialize the right-hand side vector to 0.
    rhs = zero

    ! The xm loop runs between k = 2 and k = gr%nz.  The value of xm at
    ! level k = 1, which is below the model surface, is simply set equal to the
    ! value of xm at level k = 2 after the solve has been completed.

    do k = 2, gr%nz, 1

      ! Define indices

      k_xm   = 2*k - 1
      ! k_wpxp is 2*k


      !!!!!***** xm *****!!!!!

      ! xm: Right-hand side (explicit xm portion of the code).

      ! RHS time tendency.
      rhs(k_xm) = rhs(k_xm) + xm(k) / real( dt, kind = core_rknd )

      ! RHS xm forcings.
      ! Note: xm forcings include the effects of microphysics,
      !       cloud water sedimentation, radiation, and any
      !       imposed forcings on xm.
      rhs(k_xm) = rhs(k_xm) + xm_forcing(k)

      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for xm
        !             (including microphysics/radiation).

        ! xm forcings term is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_f, k, xm_forcing(k), zt )

      endif ! l_stats_samp

    enddo ! xm loop: 2..gr%nz


    ! The wpxp loop runs between k = 2 and k = gr%nz-1.  The value of wpxp
    ! is set to specified values at both the lowest level, k = 1, and the
    ! highest level, k = gr%nz.

    do k = 2, gr%nz-1, 1

      ! Define indices

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nz )

      ! k_xm is 2*k - 1
      k_wpxp = 2*k


      !!!!!***** w'x' *****!!!!!

      ! w'x': Right-hand side (explicit w'x' portion of the code).

      ! RHS buoyancy production (bp) term and pressure term 3 (pr3).
      rhs(k_wpxp) & 
      = rhs(k_wpxp) & 
      + wpxp_terms_bp_pr3_rhs( C7_Skw_fnc(k), thv_ds_zm(k), xpthvp(k) )

      ! RHS time tendency.
      if ( l_iter ) then
        rhs(k_wpxp) = rhs(k_wpxp) + wpxp(k) / real( dt, kind = core_rknd )
      end if

      ! RHS <w'x'> forcing.
      ! Note: <w'x'> forcing includes the effects of microphysics on <w'x'>.
      rhs(k_wpxp) = rhs(k_wpxp) + wpxp_forcing(k)

      ! RHS portion of semi-implicit clipping (sicl) term.
      if ( l_clip_semi_implicit ) then
        l_upper_thresh = .true.
        l_lower_thresh = .true.

        rhs(k_wpxp) & 
        = rhs(k_wpxp) & 
        + clip_semi_imp_rhs( dt, wpxp(k), & 
                             l_upper_thresh, wpxp_upper_lim(k), & 
                             l_lower_thresh, wpxp_lower_lim(k) )

      endif

      if( .not. l_upwind_wpxp_ta ) then ! Only do this when not using Upwind Differencing
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
        !        equation (w'x' in this case), y(t) is the linearized portion of
        !        the term that gets treated implicitly, and RHS is the portion of
        !        the term that is always treated explicitly (in the case of the
        !        w'x' turbulent advection term, RHS = 0).  A weight of greater
        !        than 1 can be applied to make the term more numerically stable.
        lhs_fnc_output(1:3)  &
        = wpxp_term_ta_lhs( a1_zt(kp1), a1_zt(k),  &
                            wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                            rho_ds_zt(kp1), rho_ds_zt(k),  &
                            invrs_rho_ds_zm(k),  &
                            gr%invrs_dzm(k), k )
      else
        lhs_fnc_output(1:3)  &
        = wpxp_term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                   wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                   gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                   invrs_rho_ds_zm(k), &
                                   rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1) )
      endif

      rhs(k_wpxp)  &
        = rhs(k_wpxp)  &
        + ( one - gamma_over_implicit_ts )  &
        * ( - lhs_fnc_output(1) * wpxp(kp1)  &
            - lhs_fnc_output(2) * wpxp(k)  &
            - lhs_fnc_output(3) * wpxp(km1) )

      ! RHS contribution from "over-implicit" weighted time step
      ! for LHS pressure term 1 (pr1).
      !
      ! Note:  An "over-implicit" weighted time step is applied to this term.
      lhs_fnc_output(1)  &
      = wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )
      rhs(k_wpxp)  &
      = rhs(k_wpxp)  &
      + ( one - gamma_over_implicit_ts )  &
      * ( - lhs_fnc_output(1) * wpxp(k) )


      if ( l_stats_samp ) then

        ! Statistics: explicit contributions for wpxp.

        ! w'x' term bp is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'x' term bp, substitute 0 for the
        !        C_7 skewness function input to function wpxp_terms_bp_pr3_rhs.
        call stat_update_var_pt( iwpxp_bp, k, & 
            wpxp_terms_bp_pr3_rhs( zero, thv_ds_zm(k), xpthvp(k) ), zm )

        ! w'x' term pr3 is completely explicit; call stat_update_var_pt.
        ! Note:  To find the contribution of w'x' term pr3, add 1 to the
        !        C_7 skewness function input to function wpxp_terms_bp_pr2_rhs.
        call stat_update_var_pt( iwpxp_pr3, k, & 
            wpxp_terms_bp_pr3_rhs( (one+C7_Skw_fnc(k)), thv_ds_zm(k), &
                                   xpthvp(k) ), &
                                 zm )

        ! w'x' forcing term is completely explicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_f, k, wpxp_forcing(k), zm )

        ! w'x' term sicl has both implicit and explicit components; call
        ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
        ! subtracts the value sent in, reverse the sign on clip_semi_imp_rhs.
        if ( l_clip_semi_implicit ) then
          l_upper_thresh = .true.
          l_lower_thresh = .true.
          call stat_begin_update_pt( iwpxp_sicl, k, & 
             -clip_semi_imp_rhs( dt, wpxp(k), & 
                                 l_upper_thresh, wpxp_upper_lim(k), & 
                                 l_lower_thresh, wpxp_lower_lim(k) ), zm )
        endif

        if ( l_upwind_wpxp_ta ) then ! Use upwind differencing
          lhs_fnc_output(1:3)  &
          = wpxp_term_ta_lhs_upwind( a1(k), a1(kp1), a1(km1), &
                                     wp3_on_wp2(kp1), wp3_on_wp2(k), wp3_on_wp2(km1), &
                                     gr%invrs_dzt(k), gr%invrs_dzt(kp1), &
                                     invrs_rho_ds_zm(k), &
                                     rho_ds_zm(kp1), rho_ds_zm(k), rho_ds_zm(km1) )

        else
          ! w'x' term ta is normally completely implicit.  However, there is a
          ! RHS contribution from the "over-implicit" weighted time step.  A
          ! weighting factor of greater than 1 may be used to make the term more
          ! numerically stable (see note above for RHS contribution from
          ! "over-implicit" weighted time step for LHS turbulent advection (ta)
          ! term).  Therefore, w'x' term ta has both implicit and explicit
          ! components; call stat_begin_update_pt.  Since stat_begin_update_pt
          ! automatically subtracts the value sent in, reverse the sign on the
          ! input value.
          lhs_fnc_output(1:3)  &
          = wpxp_term_ta_lhs( a1_zt(kp1), a1_zt(k),  &
                              wp3_on_wp2_zt(kp1), wp3_on_wp2_zt(k), &
                              rho_ds_zt(kp1), rho_ds_zt(k),  &
                              invrs_rho_ds_zm(k),  &
                              gr%invrs_dzm(k), k )
        endif

        call stat_begin_update_pt( iwpxp_ta, k, &
              - ( one - gamma_over_implicit_ts )  &
              * ( - lhs_fnc_output(1) * wpxp(kp1)  &
                  - lhs_fnc_output(2) * wpxp(k)  &
                  - lhs_fnc_output(3) * wpxp(km1) ), zm )

        ! w'x' term pr1 is normally completely implicit.  However, there is a
        ! RHS contribution from the "over-implicit" weighted time step.  A
        ! weighting factor of greater than 1 may be used to make the term more
        ! numerically stable (see note above for RHS contribution from
        ! "over-implicit" weighted time step for LHS turbulent advection (ta)
        ! term).  Therefore, w'x' term pr1 has both implicit and explicit
        ! components; call stat_begin_update_pt.  Since stat_begin_update_pt
        ! automatically subtracts the value sent in, reverse the sign on the
        ! input value.
        lhs_fnc_output(1)  &
        = wpxp_term_pr1_lhs( C6x_Skw_fnc(k), tau_zm(k) )
        call stat_begin_update_pt( iwpxp_pr1, k, &
              - ( one - gamma_over_implicit_ts )  &
              * ( - lhs_fnc_output(1) * wpxp(k) ), zm )


      endif ! l_stats_samp

    enddo ! wpxp loop: 2..gr%nz-1


    ! Boundary conditions

    ! The turbulent flux (wpxp) use fixed-point boundary conditions at both the
    ! upper and lower boundaries.  Therefore, anything set in the wpxp loop
    ! at both the upper and lower boundaries would be overwritten here.
    ! However, the wpxp loop does not extend to the boundary levels.  An array
    ! with a value of 1 at the main diagonal on the left-hand side and with
    ! values of 0 at all other diagonals on the left-hand side will preserve the
    ! right-hand side value at that level.  The value of xm at level k = 1,
    ! which is below the model surface, is preserved and then overwritten to
    ! match the new value of xm at level k = 2.

    ! Lower boundary
    k      = 1
    k_xm_low   = 2*k - 1
    k_wpxp_low = 2*k

    ! Upper boundary
    k      = gr%nz
    !k_xm is 2*k - 1
    k_wpxp_high = 2*k


    ! The value of xm at the lower boundary will remain the same.
    ! However, the value of xm at the lower boundary gets overwritten
    ! after the matrix is solved for the next timestep, such
    ! that xm(1) = xm(2).

    ! The value of w'x' at the lower boundary will remain the same.
    ! The surface value of w'x' is set elsewhere
    ! (case-specific information).

    ! The value of w'x' at the upper boundary will be 0.
    call set_boundary_conditions_rhs( &
            wpxp(1), k_wpxp_low, zero, k_wpxp_high, &
            rhs, &
            xm(1), k_xm_low )


  end subroutine xm_wpxp_rhs

  !=============================================================================
  subroutine xm_wpxp_solve( nrhs, lhs, rhs, solution, err_code, rcond )

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

    use error_code, only: &
      clubb_no_error ! Constant

    use clubb_precision, only: &
      core_rknd ! Variable(s)

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
    integer, intent(out) :: err_code

    real( kind = core_rknd ), optional, intent(out) :: &
      rcond ! Est. of the reciprocal of the condition #

    err_code = clubb_no_error  ! Initialize to the value for no errors

    if ( present( rcond ) ) then
      ! Perform LU decomp and solve system (LAPACK with diagnostics)
      call band_solvex( "xm_wpxp", nsup, nsub, 2*gr%nz, nrhs, & 
                        lhs, rhs, solution, rcond, err_code )


    else
      ! Perform LU decomp and solve system (LAPACK)
      call band_solve( "xm_wpxp", nsup, nsub, 2*gr%nz, nrhs, & 
                       lhs, rhs, solution, err_code )
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
               xm, xm_tol, wpxp, err_code )

    ! Description:
    ! Clips and computes implicit stats for an artitrary xm and wpxp
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use model_flags, only: &
        l_clip_semi_implicit ! Variable(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use mono_flux_limiter, only: &
        monotonic_turbulent_flux_limit ! Procedure(s)

    use pos_definite_module, only:  & 
        pos_definite_adj ! Procedure(s)

    use clip_explicit, only: & 
        clip_covar, & ! Procedure(s)
        clip_wprtp, &      ! Variable(s)
        clip_wpthlp, &
        clip_wpsclrp

    use model_flags, only: & 
        l_pos_def, &     ! Logical for whether to apply the positive definite scheme to rtm
        l_hole_fill, &   ! Logical for whether to apply the hole filling scheme to thlm/rtm
        l_clip_turb_adv  ! Logical for whether to clip xm when wpxp is clipped

    use constants_clubb, only: &
        fstderr, & ! Constant(s)
        one, &
        zero

    use fill_holes, only: &
        fill_holes_driver ! Procedure

    use error_code, only: &
        clubb_at_least_debug_level, & ! Procedure(s)
        clubb_no_error   ! Constant

    use stats_type, only: & 
        stat_begin_update,  & ! Procedure(s)
        stat_update_var_pt, & 
        stat_end_update_pt, & 
        stat_end_update,  & 
        stat_update_var, & 
        stat_modify

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        zm, & 
        sfc, & 
        irtm_ta, & 
        irtm_ma, & 
        irtm_matrix_condt_num, & 
        irtm_pd, & 
        irtm_cl,  & 
        iwprtp_bt, & 
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
        iwpthlp_bt, & 
        iwpthlp_ma, & 
        iwpthlp_ta, & 
        iwpthlp_tp, & 
        iwpthlp_ac, & 
        iwpthlp_pr1, & 
        iwpthlp_pr2, & 
        iwpthlp_dp1, & 
        iwpthlp_sicl

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

    real(kind=time_precision), intent(in) ::  & 
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

    ! Output Variable
    integer, intent(out) ::  &
      err_code  ! Returns an error code in the event of a singular matrix

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
      iwpxp_bt, & 
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
    err_code = clubb_no_error  ! Initialize to the value for no errors

    select case ( solve_type )
    case ( xm_wpxp_rtm ) ! rtm/wprtp budget terms
      ixm_ta     = irtm_ta
      ixm_ma     = irtm_ma
      ixm_pd     = irtm_pd
      ixm_cl     = irtm_cl
      iwpxp_bt   = iwprtp_bt
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
      iwpxp_bt   = iwpthlp_bt
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

    case default  ! this includes the sclrm case
      ixm_ta     = 0
      ixm_ma     = 0
      ixm_pd     = 0
      ixm_cl     = 0
      iwpxp_bt   = 0
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
        call stat_update_var_pt( ixm_matrix_condt_num, 1, one / rcond, sfc )
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
          + ztscr03(k) * xm(kp1), zt )

        ! xm term ta is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( ixm_ta, k, & 
            ztscr04(k) * wpxp(km1) & 
          + ztscr05(k) * wpxp(k), zt )

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
          + zmscr03(k) * wpxp(kp1), zm )

!       if( .not. l_upwind_wpxp_ta ) then
          ! w'x' term ta is normally completely implicit.  However, due to the
          ! RHS contribution from the "over-implicit" weighted time step,
          ! w'x' term ta has both implicit and explicit components;
          ! call stat_end_update_pt.
          call stat_end_update_pt( iwpxp_ta, k, & 
              zmscr04(k) * wpxp(km1) & 
            + zmscr05(k) * wpxp(k) & 
            + zmscr06(k) * wpxp(kp1), zm )
!       endif

        ! w'x' term tp is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_tp, k, & 
            zmscr07(k) * xm(k) & 
          + zmscr08(k) * xm(kp1), zm )

        ! w'x' term ac is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_ac, k, & 
            zmscr09(k) * wpxp(k), zm )

        ! w'x' term pr1 is normally completely implicit.  However, due to the
        ! RHS contribution from the "over-implicit" weighted time step,
        ! w'x' term pr1 has both implicit and explicit components;
        ! call stat_end_update_pt.
        call stat_end_update_pt( iwpxp_pr1, k, & 
            zmscr10(k) * wpxp(k), zm )

        ! w'x' term pr2 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_pr2, k, & 
            zmscr11(k) * wpxp(k), zm )

        ! w'x' term dp1 is completely implicit; call stat_update_var_pt.
        call stat_update_var_pt( iwpxp_dp1, k, & 
            zmscr12(k) * wpxp(km1) & 
          + zmscr13(k) * wpxp(k) & 
          + zmscr14(k) * wpxp(kp1), zm )

        ! w'x' term sicl has both implicit and explicit components;
        ! call stat_end_update_pt.
        if ( l_clip_semi_implicit ) then
          call stat_end_update_pt( iwpxp_sicl, k, & 
              zmscr15(k) * wpxp(k), zm )
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
                                           xm, xm_tol, wpxp, err_code )
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

      call stat_update_var( iwpxp_pd, wpxp_pd(1:gr%nz), zm )

      call stat_update_var( ixm_pd, xm_pd(1:gr%nz), zt )

    end if

    ! Computed value before clipping
    if ( l_stats_samp ) then
      call stat_begin_update( ixm_cl, xm / real( dt, kind = core_rknd ), & ! Intent(in)
                              zt )                       ! Intent(inout)
    end if

    if ( any( xm < xm_threshold ) .and. l_hole_fill ) then

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

      call fill_holes_driver( 2, xm_threshold, "zt", &
                              rho_ds_zt, rho_ds_zm, &
                              xm )

    end if !  any( xm < xm_threshold ) .and. l_hole_fill

    if ( l_stats_samp ) then
      call stat_end_update( ixm_cl, xm / real( dt, kind = core_rknd ), & ! Intent(in) 
                            zt )                       ! Intent(inout)
    end if

    ! Use solve_type to find solve_type_cl, which is used
    ! in subroutine clip_covar.
    select case ( solve_type )
    case ( xm_wpxp_rtm )
      solve_type_cl = clip_wprtp
    case ( xm_wpxp_thlm )
      solve_type_cl = clip_wpthlp
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

    call clip_covar( solve_type_cl, l_first_clip_ts, &  ! In
                     l_last_clip_ts, dt, wp2, xp2_relaxed, &  ! In
                     wpxp, wpxp_chnge ) ! In/Out

    ! Adjusting xm based on clipping for w'x'.
    if ( any( wpxp_chnge /= zero ) .and. l_clip_turb_adv ) then
      call xm_correction_wpxp_cl( solve_type, dt, wpxp_chnge, gr%invrs_dzt, &
                                  xm )
    endif

    if ( l_stats_samp ) then

      ! wpxp time tendency
      call stat_modify( iwpxp_bt, wpxp / real( dt, kind = core_rknd ), zm )
      ! Brian Griffin; July 5, 2008.

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
  pure function wpxp_term_ta_lhs( wp3_on_wp2_ztp1, wp3_on_wp2_zt, &
                                  a1_ztp1, a1_zt, &
                                  rho_ds_ztp1, rho_ds_zt, &
                                  invrs_rho_ds_zm, &
                                  invrs_dzm, level ) & 
  result( lhs )

    ! Description:
    ! Turbulent advection of w'x':  implicit portion of the code.
    !
    ! The d(w'x')/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * w'^2x' )/dz.
    !
    ! A substitution is made in order to close the turbulent advection term,
    ! such that:
    !
    ! w'^2x' = a_1 * ( w'^3 / w'^2 ) * w'x',
    !
    ! where a_1 is a variable that is a function of sigma_sqd_w.  The turbulent
    ! advection term becomes:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d [ rho_ds * a_1 * ( w'^3 / w'^2 ) * w'x'(t+1) ] / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed
    !        to a "+".
    !
    ! The timestep index (t+1) means that the value of w'x' being used is from
    ! the next timestep, which is being advanced to in solving the d(w'x')/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of w'x', w'^2, and a_1 are found on the momentum levels, while
    ! the values of w'^3 are found on the thermodynamic levels.  Additionally,
    ! the values of rho_ds_zt are found on the thermodynamic levels, and the
    ! values of invrs_rho_ds_zm are found on the momentum levels.  Each of the
    ! variables w'x', w'^2, and a_1 are interpolated to the intermediate
    ! thermodynamic levels.  The values of the mathematical expression (called F
    ! here) within the dF/dz term are computed on the thermodynamic levels.
    ! Then, the derivative (d/dz) of the expression (F) is taken over the
    ! central momentum level, where it is multiplied by invrs_rho_ds_zm,
    ! yielding the desired result.  In this function, the values of F are as
    ! follows:
    !
    ! F = rho_ds_zt * a_1(t) * ( w'^3(t) / w'^2(t) ) * w'x'(t+1);
    !
    ! where the timestep index (t) stands for the index of the current timestep.
    !
    !
    ! =a1p1========wp2p1========wpxpp1=================================== m(k+1)
    !
    ! -----a1(interp)---wp2(interp)---wpxp(interp)---wp3p1---rho_ds_ztp1- t(k+1)
    !
    ! =a1==========wp2==========wpxp=======invrs_rho_ds_zm=======dF/dz=== m(k)
    !
    ! -----a1(interp)---wp2(interp)---wpxp(interp)---wp3-----rho_ds_zt--- t(k)
    !
    ! =a1m1========wp2m1========wpxpm1=================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only: & 
        gr ! Variable; gr%weights_zm2zt

!   use model_flags, only:  &
!       l_standard_term_ta

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1,    & ! Momentum superdiagonal index.
      k_mdiag   = 2,    & ! Momentum main diagonal index.
      km1_mdiag = 3       ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1,    & ! Index for upper momentum level grid weight.
      m_below = 2       ! Index for lower momentum level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      wp3_on_wp2_ztp1, & ! Smoothed wp3 / wp2 on thermo. levels (k+1)  [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels (k)    [m/s]
!     a1,              & ! a_1 interpolated to thermo. level (k+1)     [-]
      a1_ztp1,         & ! a_1 interpolated to thermo. level (k+1)     [-]
      a1_zt,           & ! a_1 interpolated to thermo. level (k)       [-]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)  [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum lev (k) [m^3/kg]
      invrs_dzm          ! Inverse of grid spacing (k)                 [1/m]

    integer, intent(in) :: & 
      level ! Central momentum level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      tkp1,  & ! Thermodynamic level directly above central momentum level.
      tk       ! Thermodynamic level directly below central momentum level.

    ! Thermodynamic level (k+1) is between momentum level (k+1)
    ! and momentum level (k).
    tkp1 = level + 1

    ! Thermodynamic level (k) is between momentum level (k)
    ! and momentum level (k-1).
    tk = level

    ! Note:  The w'x' turbulent advection term, which is
    !        - (1/rho_ds) * d [ rho_ds * a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz,
    !        still keeps the a_1 term inside the derivative, unlike the w'^3
    !        equation (found in advance_wp2_wp3_module.F90) and the equations for
    !        r_t'^2, th_l'^2, r_t'th_l', u'^2, v'^2, sclr'r_t', sclr'th_l', and
    !        sclr'^2 (found in advance_xp2_xpyp_module.F90).  Brian.

!   if ( l_standard_term_ta ) then

      ! Always use the standard discretization for the w'x' turbulent advection
      ! term.  Brian.

      ! The turbulent advection term is discretized normally, in accordance
      ! with the model equations found in the documentation and the description
      ! listed above.
      ! The w'x' turbulent advection term is
      ! - (1/rho_ds) * d [ rho_ds * a_1 * ( w'^3 / w'^2 ) * w'x' ] / dz

      ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
      lhs(kp1_mdiag) & 
      = + invrs_rho_ds_zm &
          * invrs_dzm & 
            * rho_ds_ztp1 * a1_ztp1 &
            * wp3_on_wp2_ztp1 & 
            * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
      lhs(k_mdiag) & 
      = + invrs_rho_ds_zm &
          * invrs_dzm & 
            * (   rho_ds_ztp1 * a1_ztp1 &
                  * wp3_on_wp2_ztp1 & 
                  * gr%weights_zm2zt(m_below,tkp1) & 
                - rho_ds_zt * a1_zt &
                  * wp3_on_wp2_zt & 
                  * gr%weights_zm2zt(m_above,tk) & 
              )

      ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
      lhs(km1_mdiag) & 
      = - invrs_rho_ds_zm &
          * invrs_dzm & 
            * rho_ds_zt * a1_zt &
            * wp3_on_wp2_zt & 
            * gr%weights_zm2zt(m_below,tk)

!   else

      ! This discretization very similar to what Brian did for the xp2_ta terms
      ! and is intended to stabilize the simulation by pulling a1 out of the
      ! derivative. It didn't seem to work very well.  -dschanen 17 Jan 2010

      ! Momentum superdiagonal: [ x wpxp(k+1,<t+1>) ]
!     lhs(kp1_mdiag) & 
!     = + invrs_rho_ds_zm * a1 &
!         * invrs_dzm & 
!           * rho_ds_ztp1 &
!           * wp3_on_wp2_ztp1 & 
!           * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
!     lhs(k_mdiag) & 
!     = + invrs_rho_ds_zm * a1 &
!         * invrs_dzm & 
!           * (   rho_ds_ztp1 &
!                 * wp3_on_wp2_ztp1 & 
!                 * gr%weights_zm2zt(m_below,tkp1) & 
!               - rho_ds_zt &
!                 * wp3_on_wp2_zt & 
!                 * gr%weights_zm2zt(m_above,tk) & 
!             )

!     ! Momentum subdiagonal: [ x wpxp(k-1,<t+1>) ]
!     lhs(km1_mdiag) & 
!     = - invrs_rho_ds_zm * a1 &
!         * invrs_dzm & 
!           * rho_ds_zt &
!           * wp3_on_wp2_zt & 
!           * gr%weights_zm2zt(m_below,tk)

!   endif ! l_standard_term_ta


    return
  end function wpxp_term_ta_lhs

  !=============================================================================
  pure function wpxp_term_ta_lhs_upwind( a1_zm, a1_zm_p1, a1_zm_m1, &
                                         wp3_on_wp2_p1, wp3_on_wp2, wp3_on_wp2_m1, &
                                         invrs_dzt, invrs_dztkp1, &
                                         invrs_rho_ds_zm, & 
                                         rho_ds_zmp1, rho_ds_zm, rho_ds_zmm1 ) &
  result( lhs )

    ! Description:
    !   Upwind Differencing for the wpxp term
    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1,    & ! Momentum superdiagonal index.
      k_mdiag   = 2,    & ! Momentum main diagonal index.
      km1_mdiag = 3       ! Momentum subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      a1_zm,             & ! a_1(k) on momentum levels                     [-] 
      a1_zm_p1,          & ! a_1(k+1) on momentum levels                   [-]
      a1_zm_m1,          & ! a_1(k-1) on momentum levels                   [-]
      wp3_on_wp2_p1,     & ! Smoothed wp3 / wp2 on moment. levels (k+1)    [m/s]
      wp3_on_wp2,        & ! Smoothed wp3 / wp2 on moment. levels (k)      [m/s]
      wp3_on_wp2_m1,     & ! Smoothed wp3 / wp2 on moment. levels (k-1)    [m/s]
      invrs_dzt,         & ! Inverse of grid spacing (k)                   [1/m]
      invrs_dztkp1,      & ! Inverse of grid spacing (k+1)                 [1/m]
      invrs_rho_ds_zm,   & ! Inv. dry, static density @ momentum lev (k)   [m^3/kg]
      rho_ds_zm,         & ! Density of air (k)                            [kg/m^3]
      rho_ds_zmp1,       & ! Density of air (k+1)                          [kg/m^3]
      rho_ds_zmm1          ! Density of air (k-1)                          [kg/m^3]

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs


    if ( wp3_on_wp2 > zero ) then

      ! "Wind" is blowing upwards (a1_zm > 0 and wp2 > 0 always)
      lhs(kp1_mdiag) = zero

      lhs(k_mdiag) &
      = + invrs_dzt * invrs_rho_ds_zm &
          * rho_ds_zm * a1_zm * wp3_on_wp2

      lhs(km1_mdiag) & 
      = - invrs_dzt * invrs_rho_ds_zm & 
          * rho_ds_zmm1 * a1_zm_m1 * wp3_on_wp2_m1

    else ! "Wind" is blowing downward

      lhs(kp1_mdiag) & 
      = + invrs_dztkp1 * invrs_rho_ds_zm &  
          * rho_ds_zmp1 * a1_zm_p1 * wp3_on_wp2_p1

      lhs(k_mdiag) & 
      = - invrs_dztkp1 * invrs_rho_ds_zm & 
          * rho_ds_zm * a1_zm * wp3_on_wp2

      lhs(km1_mdiag) = zero

    endif


    return
  end function wpxp_term_ta_lhs_upwind

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

  !=============================================================================
  pure function wpxp_term_pr1_lhs( C6x_Skw_fnc, tau_zm ) & 
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
      C6x_Skw_fnc,  & ! C_6x parameter with Sk_w applied (k)   [-]
      tau_zm          ! Time-scale tau at momentum level (k)   [s]

    ! Return Variable
    real( kind = core_rknd ) :: lhs


    ! Momentum main diagonal: [ x wpxp(k,<t+1>) ]
    lhs = C6x_Skw_fnc / tau_zm


    return
  end function wpxp_term_pr1_lhs

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
        time_precision, &
        core_rknd

    use stats_type, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        l_stats_samp, & ! Variable(s)
        zt, &
        ithlm_tacl, &
        irtm_tacl

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      solve_type    ! Variable that is being solved for.

    real(kind=time_precision), intent(in) :: &
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
      xm(k) = xm(k) + xm_tndcy_wpxp_cl(k) * real( dt, kind = core_rknd )
    enddo

    if ( l_stats_samp ) then
      ! The adjustment to xm due to turbulent advection term clipping
      ! (xm term tacl) is completely explicit; call stat_update_var.
      call stat_update_var( ixm_tacl, xm_tndcy_wpxp_cl, zt )
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
!===============================================================================

end module advance_xm_wpxp_module
