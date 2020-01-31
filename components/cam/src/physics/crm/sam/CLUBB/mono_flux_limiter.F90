!-----------------------------------------------------------------------
! $Id: mono_flux_limiter.F90 5715 2012-02-14 00:36:17Z dschanen@uwm.edu $
!===============================================================================
module mono_flux_limiter

  implicit none

  private ! Default Scope

  public :: monotonic_turbulent_flux_limit, &
            calc_turb_adv_range

  private :: mfl_xm_lhs, &
             mfl_xm_rhs, &
             mfl_xm_solve, &
             mean_vert_vel_up_down

  ! Private named constants to avoid string comparisons
  ! NOTE: These values must match the values for xm_wpxp_thlm
  ! and xm_wpxp_rtm given in advance_xm_wpxp_module!
  integer, parameter, private :: &
    mono_flux_thlm = 1, & ! Named constant for thlm mono_flux calls
    mono_flux_rtm = 2     ! Named constant for rtm mono_flux calls

  contains

  !=============================================================================
  subroutine monotonic_turbulent_flux_limit( solve_type, dt, xm_old, &
                                             xp2, wm_zt, xm_forcing, &
                                             rho_ds_zm, rho_ds_zt, &
                                             invrs_rho_ds_zm, invrs_rho_ds_zt, &
                                             xp2_threshold, l_implemented, &
                                             low_lev_effect, high_lev_effect, &
                                             xm, xm_tol, wpxp, err_code )

    ! Description:
    ! Limits the value of w'x' and corrects the value of xm when the xm turbulent
    ! advection term is not monotonic.  A monotonic turbulent advection scheme
    ! will not create new extrema for variable x, based only on turbulent
    ! advection (not considering mean advection and xm forcings).
    !
    ! Montonic turbulent advection
    ! ----------------------------
    !
    ! A monotonic turbulent advection scheme does not allow new extrema for
    ! variable x to be created (by means of turbulent advection).  In a
    ! monotonic turbulent advection scheme, when only the effects of turbulent
    ! advection are considered (neglecting forcings and mean advection), the
    ! value of variable x at a given point should not increase above the
    ! greatest value of variable x at nearby points, nor decrease below the
    ! smallest value of variable x at nearby points.  Nearby points are points
    ! that are close enough to the given point so that the value of variable x
    ! at the given point is effected by the values of variable x at the nearby
    ! points by means of transfer by turbulent winds during a time step.  Again,
    ! a monotonic scheme insures that advection only transfers around values of
    ! variable x and does not create new extrema for variable x.  A monotonic
    ! turbulent advection scheme is useful because the turbulent advection term
    ! (w'x') may go numerically unstable, resulting in large instabilities in
    ! the mean field (xm).  A monotonic turbulent advection scheme will limit
    ! the change in xm, and also in w'x'.
    !
    ! The following example illustrates the concept of monotonic turbulent
    ! advection.  Three successive vertical grid levels are shown (k-1, k, and
    ! k+1).  Three point values of theta-l are listed at every vertical grid
    ! level.  All three vertical levels have a mean theta-l (thlm) of 288.0 K.
    ! A circulation is occuring (in the direction of the arrows) in the vertical
    ! (w wind component) and in the horizontal (u and/or v wind components),
    ! such that the mean value of vertical velocity (wmm) is 0, but there is a
    ! turbulent component such that w'^2 > 0.
    !
    ! level = k+1 || --- 287.0 K --- 288.0 K --- 289.0 K --- || thlm = 288.0
    !             ||      / \--------------------->|         ||
    !             ||       |                       |         || wmm = 0; wp2 > 0
    !             ||       |<---------------------\ /        ||
    ! level = k   || --- 288.0 K --- 288.0 K --- 288.0 K --- || thlm = 288.0
    !             ||       |<---------------------/ \        ||
    !             ||       |                       |         || wmm = 0; wp2 > 0
    !             ||      \ /--------------------->|         ||
    ! level = k-1 || --- 287.5 K --- 288.0 K --- 288.5 K --- || thlm = 288.0
    !
    ! Neglecting any contributions from thlm forcings (effects of radiation,
    ! microphysics, large-scale horizontal advection, etc.), the values of
    ! theta-l as shown will be altered by only turbulent advection.  As a side
    ! note, the contribution of mean advection will be 0 since wmm = 0.  The
    ! diagram shows that the value of theta-l at the point on the right at level
    ! k will increase.  However, the values of theta-l at the other two points
    ! at level k will remain the same.  Thus, the value of thlm at level k will
    ! become greater than 288.0 K.  In the same manner, the values of thlm at
    ! the other two vertical levels (k-1 and k+1) will become smaller than
    ! 288.0 K.  However, the monotonic turbulent advection scheme insures that
    ! any theta-l point value cannot become smaller than the smallest theta-l
    ! point value (287.0 K) or larger than the largest theta-l point value
    ! (289.0 K).  Since all theta-l point values must fall between 287.0 K and
    ! 289.0 K, the level averages of theta-l (thlm) must fall between 287.0 K
    ! and 289.0 K.  Thus, any values of the turbulent flux, w'th_l', that would
    ! cause thlm to rise above 289.0 K or fall below 287.0 K, not considering
    ! the effect of other terms on thlm (such as forcings), are faulty and need
    ! to be limited appropriately.  The values of thlm also need to be corrected
    ! appropriately.
    !
    ! Formula for the limitation of w'x' and xm
    ! -----------------------------------------
    !
    ! The equation for change in the mean field, xm, over time is:
    !
    ! d(xm)/dt = -w*d(xm)/dz - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing;
    !
    ! where w*d(xm)/dz is the mean advection component,
    ! (1/rho_ds) * d( rho_ds * w'x' )/dz is the turbulent advection component,
    ! and xm_forcing is the xm forcing component.  The d(xm)/dt time tendency
    ! component is discretized as:
    !
    ! xm(k,<t+1>)/dt = xm(k,<t>)/dt - w*d(xm)/dz 
    !                  - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing.
    !
    ! The value of xm after it has been advanced to timestep (t+1) must be in an
    ! appropriate range based on the values of xm at timestep (t), the amount of
    ! xm forcings applied over the ensuing time step, and the amount of mean
    ! advection applied over the ensuing time step.  This is exactly the same
    ! thing as saying that the value of xm(k,<t+1>), with the contribution of
    ! turbulent advection included, must fall into a certain range based on the
    ! value of xm(k,<t+1>) without the contribution of the turbulent advection
    ! component over the last time step.  The following inequality is used to
    ! limit the value of xm(k,<t+1>):
    !
    ! MIN{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
    !         - x_max_dev_low(k-1,<t>),
    !      xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !         - x_max_dev_low(k,<t>), 
    !      xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
    !         - x_max_dev_low(k+1,<t>) }
    ! <= xm(k,<t+1>) <=
    ! MAX{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
    !         + x_max_dev_high(k-1,<t>), 
    !      xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !         + x_max_dev_high(k,<t>), 
    !      xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
    !         + x_max_dev_high(k+1,<t>) };
    !
    ! where x_max_dev_low is the absolute value of the deviation from the mean
    ! of the smallest point value of variable x at the given vertical level and
    ! timestep; and where x_max_dev_high is the deviation from the mean of the
    ! largest point value of variable x at the given vertical level and
    ! timestep.  For example, at vertical level (k+1) and timestep (t):
    !
    ! x_max_dev_low(k+1,<t>)  = | MIN( x(k+1,<t>) ) - xm(k+1,<t>) |;
    ! x_max_dev_high(k+1,<t>) = MAX( x(k+1,<t>) ) - xm(k+1,<t>).
    !
    ! The inequality shown above only takes into account values from the central
    ! level, one-level-below the central level, and one-level-above the central
    ! level.  This is the minimal amount of vertical levels that can have their
    ! values taken into consideration.  Any vertical level that can have it's
    ! properties advect to the given level during the course of a single time
    ! step can be taken into consideration.  However, only three levels will be
    ! considered in this example for the sake of simplicity.
    !
    ! The inequality will be written in more simple terms:
    !
    ! xm_lower_lim_allowable(k) <= xm(k,<t+1>) <= xm_upper_lim_allowable(k).
    !
    ! The inequality can now be related to the turbulent flux, w'x'(k,<t+1>),
    ! through a substitution that is made for xm(k,<t+1>), such that:
    !
    ! xm(k,<t+1>) = xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !               - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k).
    !
    ! The inequality becomes:
    !
    ! xm_lower_lim_allowable(k)
    ! <=
    !    xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !    - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k)
    ! <=
    ! xm_upper_lim_allowable(k).
    !
    ! The inequality is rearranged, and the turbulent advection term,
    ! d(w'x')/dz, is discretized:
    !
    ! xm_lower_lim_allowable(k)
    ! - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ]
    ! <=
    !    - dt * (1/rho_ds_zt(k))
    !           * invrs_dzt(k)
    !             * [   rho_ds_zm(k) * w'x'(k,<t+1>)
    !                 - rho_ds_zm(k-1) * w'x'(k-1,<t+1>) ]
    ! <=
    ! xm_upper_lim_allowable(k)
    ! - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ];
    !
    ! where invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).
    !
    ! Multiplying the inequality by -rho_ds_zt(k)/(dz*invrs_dzt(k)):
    !
    ! rho_ds_zt(k)/(dz*invrs_dzt(k))
    ! * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !     - xm_lower_lim_allowable(k) ]
    ! >=
    !    rho_ds_zm(k) * w'x'(k,<t+1>) - rho_ds_zm(k-1) * w'x'(k-1,<t+1>)
    ! >=
    ! rho_ds_zt(k)/(dz*invrs_dzt(k))
    ! * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !     - xm_upper_lim_allowable(k) ].
    !
    ! Note:  The inequality symbols have been flipped due to multiplication
    !        involving a (-) sign.
    !
    ! Adding rho_ds_zm(k-1) * w'x'(k-1,<t+1>) to the inequality:
    !
    ! rho_ds_zt(k)/(dz*invrs_dzt(k))
    ! * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !     - xm_lower_lim_allowable(k) ]
    ! + rho_ds_zm(k-1) * w'x'(k-1,<t+1>)
    ! >= rho_ds_zm(k) * w'x'(k,<t+1>) >=
    ! rho_ds_zt(k)/(dz*invrs_dzt(k))
    ! * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) 
    !     - xm_upper_lim_allowable(k) ]
    ! + rho_ds_zm(k-1) * w'x'(k-1,<t+1>).
    !
    ! The inequality is then rearranged to be based around w'x'(k,<t+1>):
    !
    ! (1/rho_ds_zm(k))
    ! * [ rho_ds_zt(k)/(dt*invrs_dzt(k)) 
    !     * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
    !         - xm_lower_lim_allowable(k) }
    !     + rho_ds_zm(k-1) * w'x'(k-1,<t+1>) ]
    ! >=   w'x'(k,<t+1>)   >=
    ! (1/rho_ds_zm(k))
    ! * [ rho_ds_zt(k)/(dt*invrs_dzt(k))
    !     * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) 
    !         - xm_upper_lim_allowable(k) }
    !     + rho_ds_zm(k-1) * w'x'(k-1,<t+1>) ].
    !
    ! The values of w'x' are found on the momentum levels, while the values of
    ! xm are found on the thermodynamic levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of rho_ds_zt
    ! are found on the thermodynamic levels.  The inequality is applied to
    ! w'x'(k,<t+1>) from vertical levels 2 through the second-highest level
    ! (gr%nz-1).  The value of w'x' at level 1 is a set surface (or lowest
    ! level) flux.  The value of w'x' at the highest level is also a set value,
    ! and therefore is not altered.
    !
    ! Approximating maximum and minimum values of x at any given vertical level
    ! -------------------------------------------------------------------------
    !
    ! The CLUBB code provides means, variances, and covariances for certain
    ! variables at all vertical levels.  However, there is no way to find the
    ! maximum or minimum point value of any variable on any vertical level.
    ! Without that information, x_max_dev_low and x_max_dev_high can't be found,
    ! and the inequality above is useless.  However, there is a way to
    ! approximate the maximum and minimum point values at any given vertical
    ! level.  The maximum and minimum point values can be approximated through
    ! the use of the variance, x'^2.
    !
    ! Just as the mean value of x, which is xm, and the turbulent flux of x,
    ! which is w'x', are known, so is the variance of x, which is x'^2.  The
    ! standard deviation of x is the square root of the variance of x.  The
    ! distribution of x along the horizontal plane (at vertical level k) is
    ! approximated to be the sum of two normal (or Gaussian) distributions.
    ! Most of the values in a normal distribution are found within 2 standard
    ! deviations from the mean.  Thus, the maximum point value of x along the
    ! horizontal plance at any vertical level can be approximated as:
    ! xm + 2*sqrt(x'^2).  Likewise, the minimum value of x along the horizontal
    ! plane at any vertical level can be approximated as:  xm - 2*sqrt(x'^2).
    !
    ! The values of x'^2 are found on the momentum levels.  The values of xm
    ! are found on the thermodynamic levels.  Thus, the values of x'^2 are
    ! interpolated to the thermodynamic levels in order to find the maximum
    ! and minimum point values of variable x.
    !
    ! The one downfall of this method is that instabilities can arise in the
    ! model where unphysically large values of x'^2 are produced.  Thus, this
    ! allows for an unphysically large deviation of xm from its values at the
    ! previous time step due to turbulent advection.  Thus, for purposes of
    ! determining the maximum and minimum point values of x, a upper limit
    ! is placed on x'^2, in order to limit the standard deviation of x.  This
    ! limit is only applied in this subroutine, and is not applied to x'^2
    ! elsewhere in the model code.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt  ! Procedure(s)

    use constants_clubb, only: &    
        zero_threshold, &
        eps, &
        fstderr

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use error_code, only:  &
        fatal_error, &  ! Procedure(s)
        clubb_no_error   ! Constant
        
    use fill_holes, only: &
        vertical_integral ! Procedure(s)

    use stats_type, only:  &
        stat_begin_update,  & ! Procedure(s)
        stat_end_update,  &
        stat_update_var

    use stats_variables, only:  &
        zm,  & ! Variable(s)
        zt,  &
        iwprtp_mfl,  &
        irtm_mfl,  &
        iwpthlp_mfl,  &
        ithlm_mfl,  &
        ithlm_old, &
        ithlm_without_ta, &
        ithlm_mfl_min, &
        ithlm_mfl_max, &
        irtm_old, &
        irtm_without_ta, &
        irtm_mfl_min, &
        irtm_mfl_max, &
        ithlm_enter_mfl, &
        ithlm_exit_mfl, &
        irtm_enter_mfl, &
        irtm_exit_mfl, &
        iwpthlp_mfl_min, &
        iwpthlp_mfl_max, &
        iwpthlp_entermfl, &
        iwpthlp_exit_mfl, &
        iwprtp_mfl_min, &
        iwprtp_mfl_max, &
        iwprtp_enter_mfl, &
        iwprtp_exit_mfl, &
        l_stats_samp

    implicit none

    ! Constant Parameters

    ! Flag for using a semi-implicit, tridiagonal method to solve for xm(t+1)
    ! when xm(t+1) needs to be changed.
    logical, parameter :: l_mfl_xm_imp_adj = .true.

    ! Input Variables
    integer, intent(in) ::  & 
      solve_type  ! Variables being solved for.

    real(kind=time_precision), intent(in) ::  &
      dt          ! Model timestep length                           [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      xm_old,          & ! xm at previous time step (thermo. levs.) [units vary]
      xp2,             & ! x'^2 (momentum levels)                   [units vary]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      xm_forcing,      & ! xm forcings (thermodynamic levels)       [units vary]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), intent(in) ::  &
      xp2_threshold, &   ! Lower limit of x'^2                      [units vary]
      xm_tol             ! Lower limit of maxdev                    [units vary]

    logical, intent(in) :: &
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    integer, dimension(gr%nz), intent(in) ::  &
      low_lev_effect, & ! Index of lowest level that has an effect (for lev. k)
      high_lev_effect   ! Index of highest level that has an effect (for lev. k)

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      xm,  &      ! xm at current time step (thermodynamic levels)  [units vary]
      wpxp        ! w'x' (momentum levels)                          [units vary]

    ! Output Variable
    integer, intent(out) ::  &
      err_code  ! Returns an error code in the event of a singular matrix

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xp2_zt,          &      ! x'^2 interpolated to thermodynamic levels  [units vary]
      xm_enter_mfl,    &      ! xm as it enters the MFL                    [units vary]
      xm_without_ta,   &      ! Value of xm without turb. adv. contrib.    [units vary]
      wpxp_net_adjust, &      ! Net amount of adjustment needed on w'x'    [units vary]      
      dxm_dt_mfl_adjust       ! Rate of change of adjustment to xm         [units vary]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      min_x_allowable_lev, & ! Smallest usuable value of x at lev k [units vary]
      max_x_allowable_lev, & ! Largest usuable value of x at lev k  [units vary]
      min_x_allowable, & ! Smallest usuable x within k +/- num_levs [units vary]
      max_x_allowable, & ! Largest usuable x within k +/- num_levs  [units vary]
      wpxp_mfl_max, & ! Upper limit on w'x'(k)                [units vary]
      wpxp_mfl_min    ! Lower limit on w'x'(k)                [units vary]

    real( kind = core_rknd ) ::  &
      max_xp2,             & ! Maximum allowable x'^2                        [units vary]
      stnd_dev_x,          & ! Standard deviation of x                       [units vary]
      max_dev,             & ! Determines approximate upper/lower limit of x [units vary]
      m_adv_term,          & ! Contribution of mean advection to d(xm)/dt    [units vary]
      xm_density_weighted, & ! Density weighted xm at domain top             [units vary]
      xm_adj_coef,         & ! Coeffecient to eliminate spikes at domain top [units vary]
      xm_vert_integral,    & ! Vertical integral of xm                       [units_vary]
      dz                     ! zm grid spacing at top of domain              [m]

    real( kind = core_rknd ), dimension(3,gr%nz) ::  &
      lhs_mfl_xm  ! Left hand side of tridiagonal matrix

    real( kind = core_rknd ), dimension(gr%nz) ::  &
      rhs_mfl_xm  ! Right hand side of tridiagonal matrix equation

    integer ::  &
      k, km1  ! Array indices

!    integer, parameter :: &
!      num_levs = 10  ! Number of levels above and below level k to look for
!                     ! maxima and minima of variable x.

    integer :: &
      low_lev, & ! Lowest level (from level k) to look for x minima and maxima
      high_lev   ! Highest level (from level k) to look for x minima and maxima

    integer ::  &
      iwpxp_mfl,  &
      ixm_mfl

    !--- Begin Code ---
    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Default Initialization required due to G95 compiler warning
    max_xp2 = 0.0_core_rknd
    dz = 0.0_core_rknd

    select case( solve_type )
    case ( mono_flux_rtm )  ! rtm/wprtp
       iwpxp_mfl = iwprtp_mfl
       ixm_mfl   = irtm_mfl
       max_xp2   = 5.0e-6_core_rknd
    case ( mono_flux_thlm ) ! thlm/wpthlp
       iwpxp_mfl = iwpthlp_mfl
       ixm_mfl   = ithlm_mfl
       max_xp2   = 5.0_core_rknd
    case default    ! passive scalars are involved
       iwpxp_mfl = 0
       ixm_mfl   = 0
       max_xp2   = 5.0_core_rknd
    end select


    if ( l_stats_samp ) then
       call stat_begin_update( iwpxp_mfl, wpxp / real( dt, kind = core_rknd ), zm )
       call stat_begin_update( ixm_mfl, xm / real( dt, kind = core_rknd ), zt )
    endif
    if ( l_stats_samp .and. solve_type == mono_flux_thlm ) then
       call stat_update_var( ithlm_enter_mfl, xm, zt )
       call stat_update_var( ithlm_old, xm_old, zt )
       call stat_update_var( iwpthlp_entermfl, xm, zm )
    elseif ( l_stats_samp .and. solve_type == mono_flux_rtm ) then
       call stat_update_var( irtm_enter_mfl, xm, zt )
       call stat_update_var( irtm_old, xm_old, zt )
       call stat_update_var( iwprtp_enter_mfl, xm, zm )
    endif

    ! Initialize arrays.
    wpxp_net_adjust = 0.0_core_rknd
    dxm_dt_mfl_adjust = 0.0_core_rknd

    ! Store the value of xm as it enters the mfl
    xm_enter_mfl = xm

    ! Interpolate x'^2 to thermodynamic levels.
    xp2_zt = max( zm2zt( xp2 ), xp2_threshold )

    ! Place an upper limit on xp2_zt.
    ! For purposes of this subroutine, an upper limit has been placed on the
    ! variance, x'^2.  This does not effect the value of x'^2 anywhere else in
    ! the model code.  The upper limit is a reasonable upper limit.  This is
    ! done to prevent unphysically large standard deviations caused by numerical
    ! instabilities in the x'^2 profile.
    xp2_zt = min( xp2_zt, max_xp2 )

    ! Find the maximum and minimum usuable values of variable x at each
    ! vertical level.  Start from level 2, which is the first level above
    ! the ground (or above the model surface).  This computation needs to be
    ! performed for all vertical levels above the ground (or model surface).
    do k = 2, gr%nz, 1

       km1 = max( k-1, 1 )
       !kp1 = min( k+1, gr%nz )

       ! Standard deviation is the square root of the variance.
       stnd_dev_x = sqrt( xp2_zt(k) )

       ! Most values are found within +/- 2 standard deviations from the mean.
       ! Use +/- 2 standard deviations from the mean as the maximum/minimum
       ! values.
       ! max_dev = 2.0_core_rknd*stnd_dev_x

       ! Set a minimum on max_dev
       max_dev = max(2.0_core_rknd * stnd_dev_x, xm_tol)

       ! Calculate the contribution of the mean advection term:
       ! m_adv_term = -wm_zt(k)*d(xm)/dz|_(k).
       ! Note:  mean advection is not applied to xm at level gr%nz.
       !if ( .not. l_implemented .and. k < gr%nz ) then
       !   tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k )
       !   m_adv_term = - tmp(1) * xm(kp1)  &
       !                - tmp(2) * xm(k)  &
       !                - tmp(3) * xm(km1)
       !else
       !   m_adv_term = 0.0_core_rknd
       !endif

       ! Shut off to avoid using new, possibly corrupt mean advection term
       m_adv_term = 0.0_core_rknd

       ! Find the value of xm without the contribution from the turbulent
       ! advection term.
       ! Note:  the contribution of xm_forcing at level gr%nz should be 0.
       xm_without_ta(k) = xm_old(k) + real( dt, kind = core_rknd )*xm_forcing(k) &
                          + real( dt, kind = core_rknd )*m_adv_term

       ! Find the minimum usuable value of variable x at each vertical level.
       ! Since variable x must be one of theta_l, r_t, or a scalar, all of
       ! which are positive definite quantities, the value must be >= 0.
       min_x_allowable_lev(k)  &
       = max( xm_without_ta(k) - max_dev, zero_threshold )

       ! Find the maximum usuable value of variable x at each vertical level.
       max_x_allowable_lev(k) = xm_without_ta(k) + max_dev

    enddo

    ! Boundary condition on xm_without_ta    
    k = 1
    xm_without_ta(k) = xm(k)
    min_x_allowable_lev(k) = min_x_allowable_lev(k+1)
    max_x_allowable_lev(k) = max_x_allowable_lev(k+1)

    ! Find the maximum and minimum usuable values of x that can effect the value
    ! of x at level k.  Then, find the upper and lower limits of w'x'.  Reset
    ! the value of w'x' if it is outside of those limits, and store the amount
    ! of adjustment that was needed to w'x'.
    ! The values of w'x' at level 1 and at level gr%nz are set values and
    ! are not altered.
    do k = 2, gr%nz-1, 1

       km1 = max( k-1, 1 )

       low_lev  = max( low_lev_effect(k), 2 )
       high_lev = min( high_lev_effect(k), gr%nz )
       !low_lev  = max( k-num_levs, 2 )
       !high_lev = min( k+num_levs, gr%nz )

       ! Find the smallest value of all relevant level minima for variable x.
       min_x_allowable(k) = minval( min_x_allowable_lev(low_lev:high_lev) )

       ! Find the largest value of all relevant level maxima for variable x.
       max_x_allowable(k) = maxval( max_x_allowable_lev(low_lev:high_lev) )

       ! Find the upper limit for w'x' for a monotonic turbulent flux.
       wpxp_mfl_max(k)  &
       = invrs_rho_ds_zm(k)  &
                  * (   ( rho_ds_zt(k) / (real( dt, kind = core_rknd )*gr%invrs_dzt(k)) )  &
                        * ( xm_without_ta(k) - min_x_allowable(k) )  &
                      + rho_ds_zm(km1) * wpxp(km1)  )

       ! Find the lower limit for w'x' for a monotonic turbulent flux.
       wpxp_mfl_min(k)  &
       = invrs_rho_ds_zm(k)  &
                  * (   ( rho_ds_zt(k) / (real( dt, kind = core_rknd )*gr%invrs_dzt(k)) )  &
                        * ( xm_without_ta(k) - max_x_allowable(k) )  &
                      + rho_ds_zm(km1) * wpxp(km1)  )

       if ( wpxp(k) > wpxp_mfl_max(k) ) then

          ! This block of print statements can be uncommented for debugging.
          !print *, "k = ", k
          !print *, "wpxp too large (mfl)"
          !print *, "xm(t) = ", xm_old(k)
          !print *, "xm(t+1) entering mfl = ", xm(k)
          !print *, "xm(t+1) without ta = ", xm_without_ta(k)
          !print *, "max x allowable = ", max_x_allowable(k)
          !print *, "min x allowable = ", min_x_allowable(k)
          !print *, "1/rho_ds_zm(k) = ", invrs_rho_ds_zm(k)
          !print *, "rho_ds_zt(k) = ", rho_ds_zt(k)
          !print *, "rho_ds_zt(k)*(delta_zt/dt) = ",  &
          !             real( rho_ds_zt(k) / (dt*gr%invrs_dzt(k)) )
          !print *, "xm without ta - min x allow = ",  &
          !             xm_without_ta(k) - min_x_allowable(k)
          !print *, "rho_ds_zm(km1) = ", rho_ds_zm(km1)
          !print *, "wpxp(km1) = ", wpxp(km1)
          !print *, "rho_ds_zm(km1) * wpxp(km1) = ", rho_ds_zm(km1) * wpxp(km1)
          !print *, "wpxp upper lim = ", wpxp_mfl_max(k)
          !print *, "wpxp before adjustment = ", wpxp(k)

          ! Determine the net amount of adjustment needed for w'x'.
          wpxp_net_adjust(k) = wpxp_mfl_max(k) - wpxp(k)

          ! Reset the value of w'x' to the upper limit allowed by the
          ! monotonic flux limiter.
          wpxp(k) = wpxp_mfl_max(k)

       elseif ( wpxp(k) < wpxp_mfl_min(k) ) then

          ! This block of print statements can be uncommented for debugging.
          !print *, "k = ", k
          !print *, "wpxp too small (mfl)"
          !print *, "xm(t) = ", xm_old(k)
          !print *, "xm(t+1) entering mfl = ", xm(k)
          !print *, "xm(t+1) without ta = ", xm_without_ta(k)
          !print *, "max x allowable = ", max_x_allowable(k)
          !print *, "min x allowable = ", min_x_allowable(k)
          !print *, "1/rho_ds_zm(k) = ", invrs_rho_ds_zm(k)
          !print *, "rho_ds_zt(k) = ", rho_ds_zt(k)
          !print *, "rho_ds_zt(k)*(delta_zt/dt) = ",  &
          !             real( rho_ds_zt(k) / (dt*gr%invrs_dzt(k)) )
          !print *, "xm without ta - max x allow = ",  &
          !             xm_without_ta(k) - max_x_allowable(k)
          !print *, "rho_ds_zm(km1) = ", rho_ds_zm(km1)
          !print *, "wpxp(km1) = ", wpxp(km1)
          !print *, "rho_ds_zm(km1) * wpxp(km1) = ", rho_ds_zm(km1) * wpxp(km1)
          !print *, "wpxp lower lim = ", wpxp_mfl_min(k)
          !print *, "wpxp before adjustment = ", wpxp(k)

          ! Determine the net amount of adjustment needed for w'x'.
          wpxp_net_adjust(k) = wpxp_mfl_min(k) - wpxp(k)

          ! Reset the value of w'x' to the lower limit allowed by the
          ! monotonic flux limiter.
          wpxp(k) = wpxp_mfl_min(k)

       ! This block of code can be uncommented for debugging.
       !else
       !
       !   ! wpxp(k) is okay.
       !   if ( wpxp_net_adjust(km1) /= 0.0_core_rknd ) then
       !      print *, "k = ", k
       !      print *, "wpxp is in an acceptable range (mfl)"
       !      print *, "xm(t) = ", xm_old(k)
       !      print *, "xm(t+1) entering mfl = ", xm(k)
       !      print *, "xm(t+1) without ta = ", xm_without_ta(k)
       !      print *, "max x allowable = ", max_x_allowable(k)
       !      print *, "min x allowable = ", min_x_allowable(k)
       !      print *, "1/rho_ds_zm(k) = ", invrs_rho_ds_zm(k)
       !      print *, "rho_ds_zt(k) = ", rho_ds_zt(k)
       !      print *, "rho_ds_zt(k)*(delta_zt/dt) = ",  &
       !                   real( rho_ds_zt(k) / (dt*gr%invrs_dzt(k)) )
       !      print *, "xm without ta - min x allow = ",  &
       !                   xm_without_ta(k) - min_x_allowable(k)
       !      print *, "xm without ta - max x allow = ",  &
       !                   xm_without_ta(k) - max_x_allowable(k)
       !      print *, "rho_ds_zm(km1) = ", rho_ds_zm(km1)
       !      print *, "wpxp(km1) = ", wpxp(km1)
       !      print *, "rho_ds_zm(km1) * wpxp(km1) = ",  &
       !                   rho_ds_zm(km1) * wpxp(km1)
       !      print *, "wpxp upper lim = ", wpxp_mfl_max(k)
       !      print *, "wpxp lower lim = ", wpxp_mfl_min(k)
       !      print *, "wpxp (stays the same) = ", wpxp(k)
       !   endif
       !
       endif

    enddo

    ! Boundary conditions
    min_x_allowable(1) = 0._core_rknd
    max_x_allowable(1) = 0._core_rknd

    min_x_allowable(gr%nz) = 0._core_rknd
    max_x_allowable(gr%nz) = 0._core_rknd

    wpxp_mfl_min(1) = 0._core_rknd
    wpxp_mfl_max(1) = 0._core_rknd

    wpxp_mfl_min(gr%nz) = 0._core_rknd
    wpxp_mfl_max(gr%nz) = 0._core_rknd

    if ( l_stats_samp .and. solve_type == mono_flux_thlm ) then
       call stat_update_var( ithlm_without_ta, xm_without_ta, zt )
       call stat_update_var( ithlm_mfl_min, min_x_allowable, zt )
       call stat_update_var( ithlm_mfl_max, max_x_allowable, zt )
       call stat_update_var( iwpthlp_mfl_min, wpxp_mfl_min, zm )
       call stat_update_var( iwpthlp_mfl_max, wpxp_mfl_max, zm )
    elseif ( l_stats_samp .and. solve_type == mono_flux_rtm ) then
       call stat_update_var( irtm_without_ta, xm_without_ta, zt )
       call stat_update_var( irtm_mfl_min, min_x_allowable, zt )
       call stat_update_var( irtm_mfl_max, max_x_allowable, zt )
       call stat_update_var( iwprtp_mfl_min, wpxp_mfl_min, zm )
       call stat_update_var( iwprtp_mfl_max, wpxp_mfl_max, zm )
    endif


    if ( any( wpxp_net_adjust(:) /= 0.0_core_rknd ) ) then

       ! Reset the value of xm to compensate for the change to w'x'.

       if ( l_mfl_xm_imp_adj ) then

          ! A tridiagonal matrix is used to semi-implicitly re-solve for the
          ! values of xm at timestep index (t+1).

          ! Set up the left-hand side of the tridiagonal matrix equation.
          call mfl_xm_lhs( dt, wm_zt, l_implemented, &
                           lhs_mfl_xm )

          ! Set up the right-hand side of tridiagonal matrix equation.
          call mfl_xm_rhs( dt, xm_old, wpxp, xm_forcing, &
                           rho_ds_zm, invrs_rho_ds_zt, &
                           rhs_mfl_xm )

          ! Solve the tridiagonal matrix equation.
          call mfl_xm_solve( solve_type, lhs_mfl_xm, rhs_mfl_xm,  &
                             xm, err_code )

          ! Check for errors
          if ( fatal_error( err_code ) ) return

       else  ! l_mfl_xm_imp_adj = .false.

          ! An explicit adjustment is made to the values of xm at timestep
          ! index (t+1), which is based upon the array of the amounts of w'x'
          ! adjustments.

          do k = 2, gr%nz, 1

             km1 = max( k-1, 1 )

             ! The rate of change of the adjustment to xm due to the monotonic
             ! flux limiter.
             dxm_dt_mfl_adjust(k)  &
             = - invrs_rho_ds_zt(k)  &
                 * gr%invrs_dzt(k)  &
                   * (   rho_ds_zm(k) * wpxp_net_adjust(k)  &
                       - rho_ds_zm(km1) * wpxp_net_adjust(km1) )

             ! The net change to xm due to the monotonic flux limiter is the
             ! rate of change multiplied by the time step length.  Add the
             ! product to xm to find the new xm resulting from the monotonic
             ! flux limiter.
             xm(k) = xm(k) + dxm_dt_mfl_adjust(k) * real( dt, kind = core_rknd )

          enddo

          ! Boundary condition on xm
          xm(1) = xm(2)

       endif  ! l_mfl_xm_imp_adj

       ! This code can be uncommented for debugging.
       !do k = 1, gr%nz, 1
       !   print *, "k = ", k, "xm(t) = ", xm_old(k), "new xm(t+1) = ", xm(k)
       !enddo

       !Ensure there are no spikes at the top of the domain
       if (abs( xm(gr%nz) - xm_enter_mfl(gr%nz) ) > 10._core_rknd * xm_tol) then
          dz = gr%zm(gr%nz) - gr%zm(gr%nz - 1)

          xm_density_weighted = rho_ds_zt(gr%nz) &
                                * (xm(gr%nz) - xm_enter_mfl(gr%nz)) &
                                * dz

          xm_vert_integral &
          = vertical_integral  &
              ( ((gr%nz - 1) - 2 + 1), rho_ds_zt(2:gr%nz - 1), &
                xm(2:gr%nz - 1), gr%invrs_dzt(2:gr%nz - 1) )

          !Check to ensure the vertical integral is not zero to avoid a divide
          !by zero error
          if (xm_vert_integral < eps) then
             write(fstderr,*) "Vertical integral of xm is zero;", & 
                              "mfl will remove spike at top of domain,", &
                              "but it will not conserve xm."

             !Remove the spike at the top of the domain
             xm(gr%nz) = xm_enter_mfl(gr%nz)      
          else
             xm_adj_coef = xm_density_weighted / xm_vert_integral

             !xm_adj_coef can not be smaller than -1
             if (xm_adj_coef < -0.99_core_rknd) then
                write(fstderr,*) "xm_adj_coef in mfl less than -0.99, " &
                                 // "mx_adj_coef set to -0.99"
                xm_adj_coef = -0.99_core_rknd
             endif

             !Apply the adjustment
             xm = xm * (1._core_rknd + xm_adj_coef)

             !Remove the spike at the top of the domain
             xm(gr%nz) = xm_enter_mfl(gr%nz)

             !This code can be uncommented to ensure conservation
             !if (abs(sum(rho_ds_zt(2:gr%nz) * xm(2:gr%nz) / gr%invrs_dzt(2:gr%nz)) - & 
             !    sum(rho_ds_zt(2:gr%nz) * xm_enter_mfl(2:gr%nz) / gr%invrs_dzt(2:gr%nz)))&
             !    > (1000 * xm_tol)) then
             !   write(fstderr,*) "NON-CONSERVATION in MFL", trim( solve_type ), &
             !      abs(sum(rho_ds_zt(2:gr%nz) * xm(2:gr%nz) / gr%invrs_dzt(2:gr%nz)) - &
             !       sum(rho_ds_zt(2:gr%nz) * xm_enter_mfl(2:gr%nz) / &
             !              gr%invrs_dzt(2:gr%nz)))
             !
             !   write(fstderr,*) "XM_ENTER_MFL=", xm_enter_mfl 
             !   write(fstderr,*) "XM_AFTER_SPIKE_REMOVAL", xm 
             !   write(fstderr,*) "XM_TOL", xm_tol
             !   write(fstderr,*) "XM_ADJ_COEF", xm_adj_coef   
             !endif

          endif ! xm_vert_integral < eps
       endif ! spike at domain top

    endif ! any( wpxp_net_adjust(:) /= 0.0_core_rknd )


    if ( l_stats_samp ) then

       call stat_end_update( iwpxp_mfl, wpxp / real( dt, kind = core_rknd ), zm )

       call stat_end_update( ixm_mfl, xm / real( dt, kind = core_rknd ), zt )

       if ( solve_type == mono_flux_thlm ) then
          call stat_update_var( ithlm_exit_mfl, xm, zt )
          call stat_update_var( iwpthlp_exit_mfl, xm, zm )
       elseif ( solve_type == mono_flux_rtm ) then
          call stat_update_var( irtm_exit_mfl, xm, zt )
          call stat_update_var( iwprtp_exit_mfl, xm, zm )
       endif

    endif


    return
  end subroutine monotonic_turbulent_flux_limit

  !=============================================================================
  subroutine mfl_xm_lhs( dt, wm_zt, l_implemented, &
                         lhs )

    ! Description:
    ! This subroutine is part of the process of re-solving for xm at timestep
    ! index (t+1).  This is done because the original solving process produced
    ! values outside of what is deemed acceptable by the monotonic flux limiter.
    ! Unlike the original formulation for advancing xm one timestep, which
    ! combines w'x' and xm in a band-diagonal solver, this formulation uses a
    ! tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    ! is known.
    !
    ! Subroutine mfl_xm_lhs sets up the left-hand side of the matrix equation.

    use grid_class, only: & 
        gr  ! Variable(s)

    use mean_adv, only: & 
        term_ma_zt_lhs ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real(kind=time_precision), intent(in) ::  &
      dt     ! Model timestep length                      [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm_zt  ! w wind component on thermodynamic levels   [m/s]

    logical, intent(in) :: &
      l_implemented   ! Flag for CLUBB being implemented in a larger model.

    ! Output Variables
    real( kind = core_rknd ), dimension(3,gr%nz), intent(out) ::  & 
      lhs    ! Left hand side of tridiagonal matrix

    ! Local Variables
    integer :: k, km1  ! Array index


    !-----------------------------------------------------------------------

    ! Initialize the left-hand side matrix to 0.
    lhs = 0.0_core_rknd


    ! The xm loop runs between k = 2 and k = gr%nz.  The value of xm at
    ! level k = 1, which is below the model surface, is simply set equal to the
    ! value of xm at level k = 2 after the solve has been completed.

    ! Setup LHS of the tridiagonal system
    do k = 2, gr%nz, 1

       km1 = max( k-1,1 )

       ! LHS xm mean advection (ma) term.
       if ( .not. l_implemented ) then

          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) &
          + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

       else

          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) + 0.0_core_rknd

       endif

       ! LHS xm time tendency.
       lhs(k_tdiag,k) &
       = lhs(k_tdiag,k) + 1.0_core_rknd / real( dt, kind = core_rknd )

    enddo ! xm loop: 2..gr%nz

    ! Boundary conditions.

    ! Lower boundary
    k = 1
    lhs(:,k)       = 0.0_core_rknd
    lhs(k_tdiag,k) = 1.0_core_rknd

    return
  end subroutine mfl_xm_lhs

  !=============================================================================
  subroutine mfl_xm_rhs( dt, xm_old, wpxp, xm_forcing, &
                         rho_ds_zm, invrs_rho_ds_zt, &
                         rhs )

    ! Description:
    ! This subroutine is part of the process of re-solving for xm at timestep
    ! index (t+1).  This is done because the original solving process produced
    ! values outside of what is deemed acceptable by the monotonic flux limiter.
    ! Unlike the original formulation for advancing xm one timestep, which
    ! combines w'x' and xm in a band-diagonal solver, this formulation uses a
    ! tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    ! is known.
    !
    ! Subroutine mfl_xm_rhs sets up the right-hand side of the matrix equation.

    use grid_class, only: & 
        gr  ! Variable(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) ::  &
      dt                 ! Model timestep length                    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      xm_old,          & ! xm; timestep (t) (thermodynamic levels)  [units vary]
      wpxp,            & ! w'x'; timestep (t+1); limited (m-levs.)  [units vary]
      xm_forcing,      & ! xm forcings (thermodynamic levels)       [units vary]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      rhs         ! Right hand side of tridiagonal matrix equation

    ! Local Variables
    integer :: k, km1  ! Array indices

    !-----------------------------------------------------------------------

    ! Initialize the right-hand side vector to 0.
    rhs = 0.0_core_rknd


    ! The xm loop runs between k = 2 and k = gr%nz.  The value of xm at
    ! level k = 1, which is below the model surface, is simply set equal to the
    ! value of xm at level k = 2 after the solve has been completed.

    do k = 2, gr%nz, 1

       ! Define indices
       km1 = max( k-1, 1 )

       ! RHS xm time tendency.
       rhs(k) = rhs(k) + xm_old(k) / real( dt, kind = core_rknd )

       ! RHS xm turbulent advection (ta) term.
       ! Note:  Normally, the turbulent advection (ta) term is treated
       !        implicitly when advancing xm one timestep, as both xm and w'x'
       !        are advanced together from timestep index (t) to timestep
       !        index (t+1).  However, in this case, both xm and w'x' have
       !        already been advanced one timestep.  However, w'x'(t+1) has been
       !        limited after the fact, and therefore it's values at timestep
       !        index (t+1) are known.  Thus, in re-solving for xm(t+1), the
       !        derivative of w'x'(t+1) can be placed on the right-hand side of
       !        the d(xm)/dt equation.
       rhs(k) &
       = rhs(k) &
       - invrs_rho_ds_zt(k)  &
         * gr%invrs_dzt(k)  &
           * ( rho_ds_zm(k) * wpxp(k) - rho_ds_zm(km1) * wpxp(km1) )

       ! RHS xm forcings.
       ! Note: xm forcings include the effects of microphysics,
       !       cloud water sedimentation, radiation, and any
       !       imposed forcings on xm.
       rhs(k) = rhs(k) + xm_forcing(k)

    enddo ! xm loop: 2..gr%nz

    ! Boundary conditions

    ! Lower Boundary
    k = 1
    ! The value of xm at the lower boundary will remain the same.  However, the
    ! value of xm at the lower boundary gets overwritten after the matrix is
    ! solved for the next timestep, such that xm(1) = xm(2).
    rhs(k) = xm_old(k)

    return
  end subroutine mfl_xm_rhs

  !=============================================================================
  subroutine mfl_xm_solve( solve_type, lhs, rhs,  &
                           xm, err_code )

    ! Description:
    ! This subroutine is part of the process of re-solving for xm at timestep
    ! index (t+1).  This is done because the original solving process produced
    ! values outside of what is deemed acceptable by the monotonic flux limiter.
    ! Unlike the original formulation for advancing xm one timestep, which
    ! combines w'x' and xm in a band-diagonal solver, this formulation uses a
    ! tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    ! is known.
    !
    ! Subroutine mfl_xm_solve solves the tridiagonal matrix equation for xm at
    ! timestep index (t+1).

    use grid_class, only: &
        gr  ! Variable(s)

    use lapack_wrap, only:  & 
        tridag_solve  ! Procedure(s)

    use error_code, only:  &
        fatal_error, &  ! Procedure(s)
        clubb_no_error   ! Constant

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    integer, intent(in) ::  & 
      solve_type  ! Variables being solved for.

    real( kind = core_rknd ), dimension(3,gr%nz), intent(inout) ::  & 
      lhs  ! Left hand side of tridiagonal matrix

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      rhs  ! Right hand side of tridiagonal matrix equation

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      xm   ! Value of variable being solved for at timestep (t+1)   [units vary]

    integer, intent(out) ::  &
      err_code  ! Returns an error code in the event of a singular matrix

    ! Local variable
    character(len=10) :: &
      solve_type_str ! solve_type as a string for debug output purposes

    !-----------------------------------------------------------------------

    err_code = clubb_no_error  ! Initialize to the value for no errors

    select case( solve_type )
    case ( mono_flux_rtm )
      solve_type_str = "rtm"
    case ( mono_flux_thlm )
      solve_type_str = "thlm"
    case default
      solve_type_str = "scalars"
    end select

    ! Solve for xm at timestep index (t+1) using the tridiagonal solver.
    call tridag_solve & 
         ( solve_type_str, gr%nz, 1, lhs(kp1_tdiag,:),  &  ! Intent(in)
           lhs(k_tdiag,:), lhs(km1_tdiag,:), rhs,  &         ! Intent(inout)
           xm, err_code )                                    ! Intent(out)

    ! Check for errors
    if ( fatal_error( err_code ) ) return

    ! Boundary condition on xm
    xm(1) = xm(2)

    return
  end subroutine mfl_xm_solve

  !=============================================================================
  subroutine calc_turb_adv_range( dt, w1_zm, w2_zm, varnce_w1_zm, varnce_w2_zm, &
                                  mixt_frac_zm, &
                                  low_lev_effect, high_lev_effect )

    ! Description:
    ! Calculates the lowermost and uppermost thermodynamic grid levels that can
    ! effect the base (or central) thermodynamic level through the effects of
    ! turbulent advection over the course of one time step.  This is used as
    ! part of the monotonic turbulent advection scheme.
    !
    ! One method is to use the vertical velocity at each level to determine the
    ! amount of time that it takes to travel across that particular grid level.
    ! The method is to keep on advancing one grid level until either (a) the 
    ! total sum of time taken reaches or exceeds the model time step length,
    ! (b) the top or bottom of the model is reached, or (c) a level is reached
    ! where the vertical velocity component (with turbulence included) is
    ! oriented completely opposite of the direction of travel towards the base
    ! (or central) thermodynamic level.  An example of situation (c) would be,
    ! while starting from a higher altitude and searching downward for all
    ! upward vertical velocity components, encountering a strong downdraft
    ! where the vertical velocity at every single point is oriented downward.
    ! Such a situation would occur when the mean vertical velocity (wm_zm)
    ! exceeds any turbulent component (w') that would be oriented upwards.
    !
    ! Another method is to simply set the thickness (in meters) of the layer
    ! that turbulent advection is allowed to act over, for purposes of the 
    ! monotonic turbulent advection scheme.  The lowermost and uppermost
    ! grid level that can effect the base (or central) thermodynamic level
    ! is computed based on the thickness and altitude of each level.
    
    ! References:
    !-----------------------------------------------------------------------
    
    use grid_class, only:  &
        gr  ! Variable(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    implicit none
   
    ! Constant parameters 
    logical, parameter ::  &
      l_constant_thickness = .false.  ! Toggle constant or variable thickness.

    real( kind = core_rknd ), parameter ::  &
      const_thick = 150.0_core_rknd  ! Constant thickness value               [m]

    ! Input Variables
    real(kind=time_precision), intent(in) ::  &
      dt ! Model timestep length                       [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      w1_zm,        & ! Mean w (1st PDF component)                   [m/s]
      w2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm    ! Weight of 1st PDF component (Sk_w dependent) [-]

    ! Output Variables
    integer, dimension(gr%nz), intent(out) ::  &
      low_lev_effect, & ! Index of lowest level that has an effect (for lev. k)
      high_lev_effect   ! Index of highest level that has an effect (for lev. k)

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      vert_vel_up,  & ! Average upwards vertical velocity component   [m/s]
      vert_vel_down   ! Average downwards vertical velocity component [m/s]

    real(kind=time_precision) ::  &
      dt_one_grid_lev, & ! Amount of time to travel one grid box           [s]
      dt_all_grid_levs   ! Running count of amount of time taken to travel [s]

    integer :: k, j

    ! ---- Begin Code ----

    if ( l_constant_thickness ) then ! thickness is a constant value.

       ! The value of w'x' may only be altered between levels 3 and gr%nz-2.
       do k = 3, gr%nz-2, 1

          ! Compute the number of levels that effect the central thermodynamic
          ! level through upwards motion (traveling from lower levels to reach
          ! the central thermodynamic level).

          ! Start with the index of the thermodynamic level immediately below
          ! the central thermodynamic level.
          j = k - 1

          do ! loop downwards until answer is found.

             if ( gr%zt(k) - gr%zt(j) >= const_thick ) then

                ! Stop, the current grid level is the lowest level that can
                ! be considered.
                low_lev_effect(k) = j

                exit

             else

                ! Thermodynamic level 1 cannot be considered because it is
                ! located below the surface or below the bottom of the model.
                ! The lowest level that can be considered is thermodynamic
                ! level 2.
                if ( j == 2 ) then

                   ! The current level (level 2) is the lowest level that can
                   ! be considered.
                   low_lev_effect(k) = j

                   exit

                else

                   ! Increment to the next vertical level down.
                   j = j - 1

                endif

             endif

          enddo ! downwards loop


          ! Compute the number of levels that effect the central thermodynamic
          ! level through downwards motion (traveling from higher levels to
          ! reach the central thermodynamic level).

          ! Start with the index of the thermodynamic level immediately above
          ! the central thermodynamic level.
          j = k + 1

          do ! loop upwards until answer is found.

             if ( gr%zt(j) - gr%zt(k) >= const_thick ) then

                ! Stop, the current grid level is the highest level that can
                ! be considered.
                high_lev_effect(k) = j

                exit

             else

                ! The highest level that can be considered is thermodynamic
                ! level gr%nz.
                if ( j == gr%nz ) then

                   ! The current level (level gr%nz) is the highest level
                   ! that can be considered.
                   high_lev_effect(k) = j

                   exit

                else

                   ! Increment to the next vertical level up.
                   j = j + 1

                endif

             endif

          enddo ! upwards loop

       enddo ! k = 3, gr%nz-2


    else ! thickness based on vertical velocity and time step length.

       ! Find the average upwards vertical velocity and the average downwards
       ! vertical velocity.
       ! Note:  A level that has all vertical wind moving downwards will have a
       !        vert_vel_up value that is 0, and vice versa.
       call mean_vert_vel_up_down( w1_zm, w2_zm, varnce_w1_zm, varnce_w2_zm, & !  In
                                   mixt_frac_zm, 0.0_core_rknd,  & ! In
                                   vert_vel_down, vert_vel_up )

       ! The value of w'x' may only be altered between levels 3 and gr%nz-2.
       do k = 3, gr%nz-2, 1

          ! Compute the number of levels that effect the central thermodynamic
          ! level through upwards motion (traveling from lower levels to reach
          ! the central thermodynamic level).

          ! Start with the index of the thermodynamic level immediately below
          ! the central thermodynamic level.
          j = k - 1

          ! Initialize the overall delta t counter to 0.
          dt_all_grid_levs = 0.0_time_precision

          do ! loop downwards until answer is found.

             ! Continue if there is some component of upwards vertical velocity.
             if ( vert_vel_up(j) > 0.0_core_rknd ) then

                ! Compute the amount of time it takes to travel one grid level
                ! upwards:  delta_t = delta_z / vert_vel_up.
                dt_one_grid_lev = real( (1.0_core_rknd/gr%invrs_dzm(j)) / vert_vel_up(j), &
                                        kind=time_precision )

                ! Total time elapsed for crossing all grid levels that have been
                ! passed, thus far.
                dt_all_grid_levs = dt_all_grid_levs + dt_one_grid_lev

                ! Stop if has taken more than one model time step (overall) to
                ! travel the entire extent of the current vertical grid level.
                if ( dt_all_grid_levs >= dt ) then

                   ! The current level is the lowest level that can be
                   ! considered.
                   low_lev_effect(k) = j

                   exit

                ! Continue if the total elapsed time has not reached or exceeded
                ! one model time step.
                else

                   ! Thermodynamic level 1 cannot be considered because it is
                   ! located below the surface or below the bottom of the model.
                   ! The lowest level that can be considered is thermodynamic
                   ! level 2.
                   if ( j == 2 ) then

                      ! The current level (level 2) is the lowest level that can
                      ! be considered.
                      low_lev_effect(k) = j

                      exit

                   else

                      ! Increment to the next vertical level down.
                      j = j - 1

                   endif

                endif

             ! Stop if there isn't a component of upwards vertical velocity.
             else

                ! The current level cannot be considered.  The lowest level that
                ! can be considered is one-level-above the current level.
                low_lev_effect(k) = j + 1

                exit

             endif

          enddo ! downwards loop


          ! Compute the number of levels that effect the central thermodynamic
          ! level through downwards motion (traveling from higher levels to
          ! reach the central thermodynamic level).

          ! Start with the index of the thermodynamic level immediately above
          ! the central thermodynamic level.
          j = k + 1

          ! Initialize the overall delta t counter to 0.
          dt_all_grid_levs = 0.0_time_precision

          do ! loop upwards until answer is found.

             ! Continue if there is some component of downwards vertical velocity.
             if ( vert_vel_down(j-1) < 0.0_core_rknd ) then

                ! Compute the amount of time it takes to travel one grid level
                ! downwards:  delta_t = - delta_z / vert_vel_down.
                ! Note:  There is a (-) sign in front of delta_z because the
                !        distance traveled is downwards.  Since vert_vel_down
                !        has a negative value, dt_one_grid_lev will be a
                !        positive value.
                dt_one_grid_lev = real( -(1.0_core_rknd/gr%invrs_dzm(j-1)) / vert_vel_down(j-1), &
                                        kind=time_precision )

                ! Total time elapsed for crossing all grid levels that have been
                ! passed, thus far.
                dt_all_grid_levs = real( dt_all_grid_levs + dt_one_grid_lev, kind=time_precision )

                ! Stop if has taken more than one model time step (overall) to
                ! travel the entire extent of the current vertical grid level.
                if ( dt_all_grid_levs >= dt ) then

                   ! The current level is the highest level that can be
                   ! considered.
                   high_lev_effect(k) = j

                   exit

                ! Continue if the total elapsed time has not reached or exceeded
                ! one model time step.
                else

                   ! The highest level that can be considered is thermodynamic
                   ! level gr%nz.
                   if ( j == gr%nz ) then

                      ! The current level (level gr%nz) is the highest level
                      ! that can be considered.
                      high_lev_effect(k) = j

                      exit

                   else

                      ! Increment to the next vertical level up.
                      j = j + 1

                   endif

                endif

             ! Stop if there isn't a component of downwards vertical velocity.
             else

                ! The current level cannot be considered.  The highest level
                ! that can be considered is one-level-below the current level.
                high_lev_effect(k) = j - 1

                exit

             endif

          enddo  ! upwards loop

       enddo ! k = 3, gr%nz-2

    endif ! l_constant_thickness


    ! Information for levels 1, 2, gr%nz-1, and gr%nz is not needed.
    ! However, set the values at these levels for purposes of not having odd
    ! values in the arrays.
    low_lev_effect(1)  = 1
    high_lev_effect(1) = 1
    low_lev_effect(2)  = 2
    high_lev_effect(2) = 2
    low_lev_effect(gr%nz-1)  = gr%nz-1
    high_lev_effect(gr%nz-1) = gr%nz
    low_lev_effect(gr%nz)    = gr%nz
    high_lev_effect(gr%nz)   = gr%nz


    return
  end subroutine calc_turb_adv_range

  !=============================================================================
  subroutine mean_vert_vel_up_down( w1_zm, w2_zm, varnce_w1_zm, varnce_w2_zm, &
                                    mixt_frac_zm, w_ref, &
                                    mean_w_down, mean_w_up )

    ! Description
    ! The values of vertical velocity, along a horizontal plane at any given
    ! vertical level, are not allowed by CLUBB to be uniform.  In other words,
    ! there must be some variance in vertical velocity.  This subroutine
    ! calculates the mean of all values of vertical velocity, at any given
    ! vertical level, that are greater than a certain reference velocity.  This
    ! subroutine also calculates the mean of all values of vertical velocity, at
    ! any given vertical level, that are less than a certain reference velocity.
    ! The reference velocity is usually 0 m/s, in which case this subroutine
    ! calculates the average positive (upward) velocity and the average negative
    ! (downward) velocity.  However, the reference velocity may be other values,
    ! such as wm_zm, which is the overall mean vertical velocity.  If the
    ! reference velocity is wm_zm, this subroutine calculates the average of all
    ! values of w that are on the positive ("upward") side of the mean and the
    ! average of all values of w that are on the negative ("downward") side of
    ! the mean.  These mean positive and negative vertical velocities are useful
    ! in determining how long, on average, it takes a parcel of air, being
    ! driven by subgrid updrafts or downdrafts, to traverse the length of the
    ! vertical grid level.
    !
    ! Method
    ! ------
    !
    ! The CLUBB model uses a joint PDF of vertical velocity, liquid water
    ! potential temperature, and total water mixing ratio to determine subgrid
    ! variability.
    !
    ! The values of vertical velocity, w, along an undefined horizontal plane
    ! at any vertical level, are considered to approximately follow a
    ! distribution that is a mixture of two normal (or Gaussian) distributions.
    ! The values of w that are a part of the 1st normal distribution are
    ! referred to as w1, and the values of w that are part of the 2nd normal
    ! distribution are referred to as w2.  Note that these distributions
    ! overlap, and there are many values of w that are found in both w1 and w2.
    !
    ! The probability density function (PDF) for w, P(w), is:
    !
    ! P(w) = mixt_frac*P(w1) + (1-mixt_frac)*P(w2);
    !
    ! where "mixt_frac" is the weight of the 1st normal distribution, and P(w1) and
    ! P(w2) are the equations for the 1st and 2nd normal distributions,
    ! respectively:
    !
    ! P(w1) = 1 / ( sigma_w1 * sqrt(2*PI) ) 
    !         * EXP[ -(w1-mu_w1)^2 / (2*sigma_w1^2) ]; and
    !
    ! P(w2) = 1 / ( sigma_w2 * sqrt(2*PI) ) 
    !         * EXP[ -(w2-mu_w2)^2 / (2*sigma_w2^2) ].
    !
    ! The mean of the 1st normal distribution is mu_w1, and the standard
    ! deviation of the 1st normal distribution is sigma_w1.  The mean of the
    ! 2nd normal distribution is mu_w2, and the standard deviation of the 2nd
    ! normal distribution is sigma_w2.
    !
    ! The average value of w, distributed according to the probability
    ! distribution, between limits alpha and beta, is:
    !
    ! <w|_(alpha:beta)> = INT(alpha:beta) w P(w) dw.
    !
    ! The average value of w over a certain domain is used to determine the
    ! average positive and negative (as compared to the reference velocity)
    ! values of w at any vertical level.
    !
    ! Average Negative Vertical Velocity
    ! ----------------------------------
    !
    ! The average of all values of w in the distribution that are below the
    ! reference velocity, w|_ref, is the mean value of w over the domain
    ! -inf <= w <= w|_ref, such that:
    !
    ! <w|_(-inf:w|_ref)> = INT(-inf:w|_ref) w P(w) dw.
    !                    = mixt_frac * INT(-inf:w|_ref) w1 P(w1) dw1
    !                      + (1-mixt_frac) * INT(-inf:w|_ref) w2 P(w2) dw2.
    !
    ! For each normal distribution in the mixture of normal distribution, i
    ! (where "i" can be 1 or 2):
    !
    ! INT(-inf:w|_ref) wi P(wi) dwi =
    !   - ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
    !   + mu_wi * (1/2)*[ 1 + erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];
    !
    ! where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
    ! the standard deviations of w for the ith normal distribution, and erf( )
    ! is the error function.
    !
    ! The mean of all values of w <= w|_ref is:
    !
    ! <w|_(-inf:w|_ref)> =
    ! mixt_frac * { - ( sigma_w1 / sqrt(2*PI) ) 
    !                 * EXP[ -(w|_ref-mu_w1)^2 / (2*sigma_w1^2) ]
    !               + mu_w1 * (1/2)
    !                 *[1 + erf( (w|_ref-mu_w1) / (sqrt(2)*sigma_w1) )] }
    ! + (1-mixt_frac) * { - ( sigma_w2 / sqrt(2*PI) ) 
    !                       * EXP[ -(w|_ref-mu_w2)^2 / (2*sigma_w2^2) ]
    !                     + mu_w2 * (1/2)
    !                       *[1 + erf( (w|_ref-mu_w2) / (sqrt(2)*sigma_w2) )] }.
    !
    ! Average Positive Vertical Velocity
    ! ----------------------------------
    !
    ! The average of all values of w in the distribution that are above the
    ! reference velocity, w|_ref, is the mean value of w over the domain
    ! w|_ref <= w <= inf, such that:
    !
    ! <w|_(w|_ref:inf)> = INT(w|_ref:inf) w P(w) dw.
    !                   = mixt_frac * INT(w|_ref:inf) w1 P(w1) dw1
    !                     + (1-mixt_frac) * INT(w|_ref:inf) w2 P(w2) dw2.
    !
    ! For each normal distribution in the mixture of normal distribution, i
    ! (where "i" can be 1 or 2):
    !
    ! INT(w|_ref:inf) wi P(wi) dwi =
    !     ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
    !   + mu_wi * (1/2)*[ 1 - erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];
    !
    ! where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
    ! the standard deviations of w for the ith normal distribution, and erf( )
    ! is the error function.
    !
    ! The mean of all values of w >= w|_ref is:
    !
    ! <w|_(w|_ref:inf)> =
    ! mixt_frac * {   ( sigma_w1 / sqrt(2*PI) ) 
    !                * EXP[ -(w|_ref-mu_w1)^2 / (2*sigma_w1^2) ]
    !               + mu_w1 * (1/2)
    !                 *[1 - erf( (w|_ref-mu_w1) / (sqrt(2)*sigma_w1) )] }
    ! + (1-mixt_frac) * {   ( sigma_w2 / sqrt(2*PI) ) 
    !                      * EXP[ -(w|_ref-mu_w2)^2 / (2*sigma_w2^2) ]
    !                     + mu_w2 * (1/2)
    !                       *[1 - erf( (w|_ref-mu_w2) / (sqrt(2)*sigma_w2) )] }.
    !
    ! Special Limitations:
    ! --------------------
    !
    ! A normal distribution has a domain from -inf to inf.  However, the mixture
    ! of normal distributions is an approximation of the distribution of values
    ! of w along a horizontal plane at any given vertical level.  Vertical
    ! velocity, w, has absolute minimum and maximum values (that cannot be
    ! predicted by the PDF).  The absolute maximum and minimum for each normal
    ! distribution is most likely found within 2 or 3 standard deviations of the
    ! mean for the relevant normal distribution.  In other words, for each
    ! normal distribution in the mixture of normal distributions, all the values
    ! of w are found within 2 or 3 standard deviations on both sides of the
    ! mean.  Therefore, if one (or both) of the normal distributions has a mean
    ! that is more than 3 standard deviations away from the reference velocity,
    ! then that entire w distribution is found on ONE side of the reference
    ! velocity.
    !
    ! Therefore:
    !
    ! a) where mu_wi + 3*sigma_wi <= w|_ref:
    !
    !       The entire ith normal distribution of w is on the negative side of
    !       w|_ref; and
    !
    !       INT(-inf:w|_ref) wi P(wi) dwi = mu_wi; and
    !       INT(inf:w|_ref) wi P(wi) dwi = 0.
    !
    ! b) where mu_wi - 3*sigma_wi >= w|_ref:
    !
    !       The entire ith normal distribution of w is on the positive side of
    !       w|_ref; and
    !
    !       INT(-inf:w|_ref) wi P(wi) dwi = 0; and
    !       INT(inf:w|_ref) wi P(wi) dwi = mu_wi.
    !
    ! Note:  A value of 3 standard deviations above and below the mean of the
    !        ith normal distribution was chosen for the approximate maximum and
    !        minimum values of the ith normal distribution because 99.7% of
    !        values in a normal distribution are found within 3 standard
    !        deviations from the mean (compared to 95.4% for 2 standard
    !        deviations).  The value of 3 standard deviations provides for a
    !        reasonable estimate of the absolute maximum and minimum of w, while
    !        covering a great majority of the normal distribution.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        gr,  & ! Variable(s)
        zt2zm  ! Procedure(s)

    use constants_clubb, only: &
        sqrt_2pi, &
        sqrt_2

    use anl_erf, only:  & 
        erf ! Procedure(s)
            ! The error function

    use stats_type, only:  &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only:  &
        zm,  & ! Variable(s)
        imean_w_up, &
        imean_w_down, &
        l_stats_samp

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      w1_zm,        & ! Mean w (1st PDF component)                   [m/s]
      w2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm    ! Weight of 1st PDF component (Sk_w dependent) [-]

    real( kind = core_rknd ), intent(in) ::  &
      w_ref          ! Reference velocity, w|_ref (normally = 0)   [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mean_w_down, & ! Overall mean w (<= w|_ref)                  [m/s]
      mean_w_up      ! Overall mean w (>= w|_ref)                  [m/s]

    ! Local Variables

    real( kind = core_rknd ) :: &
      sigma_w1, & ! Standard deviation of w for 1st normal distribution    [m/s]
      sigma_w2, & ! Standard deviation of w for 2nd normal distribution    [m/s]
      mean_w_down_1st, & ! Mean w (<= w|_ref) from 1st normal distribution [m/s]
      mean_w_down_2nd, & ! Mean w (<= w|_ref) from 2nd normal distribution [m/s]
      mean_w_up_1st, &   ! Mean w (>= w|_ref) from 1st normal distribution [m/s]
      mean_w_up_2nd, &   ! Mean w (>= w|_ref) from 2nd normal distribution [m/s]
      exp_cache, & ! Cache of exponential calculations to reduce runtime
      erf_cache    ! Cache of error function calculations to reduce runtime

    integer :: k  ! Vertical loop index

    ! ---- Begin Code ----

    ! Loop over momentum levels from 2 to gr%nz-1.  Levels 1 and gr%nz
    ! are not needed.
    do k = 2, gr%nz-1, 1

       ! Standard deviation of w for the 1st normal distribution.
       sigma_w1 = sqrt( varnce_w1_zm(k) )

       ! Standard deviation of w for the 2nd normal distribution.
       sigma_w2 = sqrt( varnce_w2_zm(k) )


       ! Contributions from the 1st normal distribution.
       if ( w1_zm(k) + 3._core_rknd*sigma_w1 <= w_ref ) then

          ! The entire 1st normal is on the negative side of w|_ref.
          mean_w_down_1st = w1_zm(k)
          mean_w_up_1st   = 0.0_core_rknd

       elseif ( w1_zm(k) - 3._core_rknd*sigma_w1 >= w_ref ) then

          ! The entire 1st normal is on the positive side of w|_ref.
          mean_w_down_1st = 0.0_core_rknd
          mean_w_up_1st   = w1_zm(k)

       else

          ! The exponential calculation is pulled out as it is reused in both
          ! equations. This should save one calculation of the
          ! exp( -(w_ref-w1_zm(k))**2 ... etc. part of the formula.
          ! ~~EIHoppe//20090618
          exp_cache = exp( -(w_ref-w1_zm(k))**2 / (2.0_core_rknd*sigma_w1**2) ) 

          ! Added cache of the error function calculations.
          ! This should save one calculation of the erf(...) part
          ! of the formula.
          ! ~~EIHoppe//20090623
          erf_cache = erf( (w_ref-w1_zm(k)) / (sqrt_2*sigma_w1) )

          ! The 1st normal has values on both sides of w_ref.
          mean_w_down_1st =  &
             - (sigma_w1/sqrt_2pi)  &
!             * exp( -(w_ref-w1_zm(k))**2 / (2.0_core_rknd*sigma_w1**2) )  &
             * exp_cache  &
!             + w1(k) * 0.5_core_rknd*( 1.0_core_rknd + erf( (w_ref-w1(k)) / (sqrt_2*sigma_w1) ) )
             + w1_zm(k) * 0.5_core_rknd*( 1.0_core_rknd + erf_cache)

          mean_w_up_1st =  &
             + (sigma_w1/sqrt_2pi)  &
!             * exp( -(w_ref-w1(k))**2 / (2.0_core_rknd*sigma_w1**2) )  &
             * exp_cache  &
!             + w1(k) * 0.5_core_rknd*( 1.0_core_rknd - erf( (w_ref-w1(k)) / (sqrt_2*sigma_w1) ) )
             + w1_zm(k) * 0.5_core_rknd*( 1.0_core_rknd - erf_cache)

          ! /EIHoppe changes

       endif


       ! Contributions from the 2nd normal distribution.
       if ( w2_zm(k) + 3._core_rknd*sigma_w2 <= w_ref ) then

          ! The entire 2nd normal is on the negative side of w|_ref.
          mean_w_down_2nd = w2_zm(k)
          mean_w_up_2nd   = 0.0_core_rknd

       elseif ( w2_zm(k) - 3._core_rknd*sigma_w2 >= w_ref ) then

          ! The entire 2nd normal is on the positive side of w|_ref.
          mean_w_down_2nd = 0.0_core_rknd
          mean_w_up_2nd   = w2_zm(k)

       else

          ! The exponential calculation is pulled out as it is reused in both
          ! equations. This should save one calculation of the
          ! exp( -(w_ref-w1(k))**2 ... etc. part of the formula.
          ! ~~EIHoppe//20090618
          exp_cache = exp( -(w_ref-w2_zm(k))**2 / (2.0_core_rknd*sigma_w2**2) ) 

          ! Added cache of the error function calculations.
          ! This should save one calculation of the erf(...) part
          ! of the formula.
          ! ~~EIHoppe//20090623
          erf_cache = erf( (w_ref-w2_zm(k)) / (sqrt_2*sigma_w2) )

          ! The 2nd normal has values on both sides of w_ref.
          mean_w_down_2nd =  &
             - (sigma_w2/sqrt_2pi)  &
!            * exp( -(w_ref-w2_zm(k))**2 / (2.0_core_rknd*sigma_w2**2) )  &
             * exp_cache  &
!            + w2_zm(k) * 0.5_core_rknd*( 1.0_core_rknd + erf( (w_ref-w2(k)) / (sqrt_2*sigma_w2) ) )
             + w2_zm(k) * 0.5_core_rknd*( 1.0_core_rknd + erf_cache)

          mean_w_up_2nd =  &
             + (sigma_w2/sqrt_2pi)  &
!            * exp( -(w_ref-w2(k))**2 / (2.0_core_rknd*sigma_w2**2) )  &
             * exp_cache  &
!            + w2(k) * 0.5_core_rknd*( 1.0_core_rknd - erf( (w_ref-w2(k)) / (sqrt_2*sigma_w2) ) )
             + w2_zm(k) * 0.5_core_rknd*( 1.0_core_rknd - erf_cache)

          ! /EIHoppe changes

       endif

       ! Overall mean of downwards w.
       mean_w_down(k) = mixt_frac_zm(k) * mean_w_down_1st  &
                        + ( 1.0_core_rknd - mixt_frac_zm(k) ) * mean_w_down_2nd

       ! Overall mean of upwards w.
       mean_w_up(k) = mixt_frac_zm(k) * mean_w_up_1st  &
                      + ( 1.0_core_rknd - mixt_frac_zm(k) ) * mean_w_up_2nd

       if ( l_stats_samp ) then

          call stat_update_var_pt( imean_w_up, k, mean_w_up(k), zm )

          call stat_update_var_pt( imean_w_down, k, mean_w_down(k), zm )

       endif ! l_stats_samp

    enddo ! k = 2, gr%nz, 1


    return
  end subroutine mean_vert_vel_up_down

!===============================================================================

end module mono_flux_limiter
