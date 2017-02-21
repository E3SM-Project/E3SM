!-------------------------------------------------------------------------
! $Id: advance_helper_module.F90 8133 2016-06-12 19:12:36Z raut@uwm.edu $
!===============================================================================
module advance_helper_module

! Description:
!   This module contains helper methods for the advance_* modules.
!------------------------------------------------------------------------

  implicit none

  public :: &
    set_boundary_conditions_lhs, &
    set_boundary_conditions_rhs, &
    calc_stability_correction,   &
    calc_brunt_vaisala_freq_sqd, &
    compute_Cx_fnc_Richardson

  private ! Set Default Scope

  contains

  !---------------------------------------------------------------------------
  subroutine set_boundary_conditions_lhs( diag_index, low_bound, high_bound, lhs, &
                                      diag_index2, low_bound2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a left-hand side LAPACK matrix.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables
    integer, intent(in) :: &
      diag_index, low_bound, high_bound ! boundary indexes for the first variable

    ! Input / Output Variables
    real( kind = core_rknd ), dimension(:,:), intent(inout) :: &
      lhs ! left hand side of the LAPACK matrix equation

    ! Optional Input Variables
    integer, intent(in), optional :: &
      diag_index2, low_bound2, high_bound2 ! boundary indexes for the second variable

    ! --------------------- BEGIN CODE ----------------------

    if ( ( present( low_bound2 ) .or. present( high_bound2 ) ) .and. &
         ( .not. present( diag_index2 ) ) ) then

      stop "Boundary index provided without diag_index."

    end if

    ! Set the lower boundaries for the first variable
    lhs(:,low_bound) = 0.0_core_rknd
    lhs(diag_index,low_bound) = 1.0_core_rknd

    ! Set the upper boundaries for the first variable
    lhs(:,high_bound) = 0.0_core_rknd
    lhs(diag_index,high_bound) = 1.0_core_rknd

    ! Set the lower boundaries for the second variable, if it is provided
    if ( present( low_bound2 ) ) then

      lhs(:,low_bound2) = 0.0_core_rknd
      lhs(diag_index2,low_bound2) = 1.0_core_rknd

    end if

    ! Set the upper boundaries for the second variable, if it is provided
    if ( present( high_bound2 ) ) then

      lhs(:,high_bound2) = 0.0_core_rknd
      lhs(diag_index2,high_bound2) = 1.0_core_rknd

    end if

    return
  end subroutine set_boundary_conditions_lhs

  !--------------------------------------------------------------------------
  subroutine set_boundary_conditions_rhs( &
               low_value, low_bound, high_value, high_bound, &
               rhs, &
               low_value2, low_bound2, high_value2, high_bound2 )

  ! Description:
  !   Sets the boundary conditions for a right-hand side LAPACK vector.
  !
  ! References:
  !   none
  !---------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Exernal 
    intrinsic :: present

    ! Input Variables

    ! The values for the first variable
    real( kind = core_rknd ), intent(in) :: low_value, high_value

    ! The bounds for the first variable
    integer, intent(in) :: low_bound, high_bound

    ! Input / Output Variables

    ! The right-hand side vector
    real( kind = core_rknd ), dimension(:), intent(inout) :: rhs

    ! Optional Input Variables

    ! The values for the second variable
    real( kind = core_rknd ), intent(in), optional :: low_value2, high_value2

    ! The bounds for the second variable
    integer, intent(in), optional :: low_bound2, high_bound2


    ! -------------------- BEGIN CODE ------------------------

    ! Stop execution if a boundary was provided without a value
    if ( (present( low_bound2 ) .and. (.not. present( low_value2 ))) .or. &
         (present( high_bound2 ) .and. (.not. present( high_value2 ))) ) then

      stop "Boundary condition provided without value."

    end if

    ! Set the lower and upper bounds for the first variable
    rhs(low_bound) = low_value
    rhs(high_bound) = high_value

    ! If a lower bound was given for the second variable, set it
    if ( present( low_bound2 ) ) then
      rhs(low_bound2) = low_value2
    end if

    ! If an upper bound was given for the second variable, set it
    if ( present( high_bound2 ) ) then
      rhs(high_bound2) = high_value2
    end if

    return
  end subroutine set_boundary_conditions_rhs

  !===============================================================================
  function calc_stability_correction( thlm, Lscale, em, exner, rtm, rcm, p_in_Pa, &
                                      cloud_frac, thvm ) &
    result ( stability_correction )
  !
  ! Description:
  !   Stability Factor
  !
  ! References:
  !
  !--------------------------------------------------------------------

    use parameters_tunable, only: &
      lambda0_stability_coef ! Variable(s)

    use constants_clubb, only: &
      zero    ! Constant(s)

    use grid_class, only:  &
      gr, & ! Variable(s)
      zt2zm    ! Procedure(s)

    use clubb_precision, only:  &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      Lscale,          & ! Turbulent mixing length                   [m]
      em,              & ! Turbulent Kinetic Energy (TKE)            [m^2/s^2]
      thlm,            & ! th_l (thermo. levels)                     [K]
      exner,           & ! Exner function                            [-]
      rtm,             & ! total water mixing ratio, r_t             [kg/kg]
      rcm,             & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,         & ! Air pressure                              [Pa]
      cloud_frac,      & ! Cloud fraction                            [-]
      thvm               ! Virtual potential temperature             [K]

    ! Result
    real( kind = core_rknd ), dimension(gr%nz) :: &
      stability_correction

    real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd, & !  []
      lambda0_stability

    !------------ Begin Code --------------
    brunt_vaisala_freq_sqd = calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, &
                                                          cloud_frac, thvm )
    lambda0_stability = merge( lambda0_stability_coef, zero, brunt_vaisala_freq_sqd > zero )

    stability_correction = 1.0_core_rknd &
      + min( lambda0_stability * brunt_vaisala_freq_sqd * zt2zm( Lscale )**2 / em, 3.0_core_rknd )

    return
  end function calc_stability_correction

  !===============================================================================
  function calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, cloud_frac, thvm ) &
    result( brunt_vaisala_freq_sqd )

  ! Description:
  !   Calculate the Brunt-Vaisala frequency squared, N^2.

  ! References:
  !   ?
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Konstant

    use constants_clubb, only: &
      grav, & ! Constant(s)
      cloud_frac_min, &
      Lv, Cp, Rd, ep, &
      one

    use parameters_model, only: & 
      T0 ! Variable! 

    use grid_class, only: &
      gr,     & ! Variable
      ddzt,   &  ! Procedure(s)
      zt2zm

    use T_in_K_module, only: &
      thlm2T_in_K ! Procedure

    use saturation, only: &
      sat_mixrat_liq ! Procedure

    use model_flags, only: &
      l_brunt_vaisala_freq_moist, & ! Variable(s)
      l_use_thvm_in_bv_freq

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      thlm,    &  ! th_l (thermo. levels)              [K]
      exner,   &  ! Exner function                     [-]
      rtm,     &  ! total water mixing ratio, r_t      [kg/kg]
      rcm,     &  ! cloud water mixing ratio, r_c      [kg/kg]
      p_in_Pa, &  ! Air pressure                       [Pa]
      cloud_frac, & ! Cloud fraction                   [-]
      thvm        ! Virtual potential temperature      [K]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd ! Brunt-Vaisala frequency squared, N^2 [1/s^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      T_in_K, T_in_K_zm, rsat, rsat_zm, thm, thm_zm, ddzt_thlm, &
      ddzt_thm, ddzt_rsat, ddzt_rtm, thvm_zm, ddzt_thvm

    integer :: k

  !---------------------------------------------------------------------
    !----- Begin Code -----
    ddzt_thlm = ddzt( thlm )
    thvm_zm = zt2zm( thvm )
    ddzt_thvm = ddzt( thvm )

    if ( l_brunt_vaisala_freq_moist ) then
      ! These parameters are needed to compute the moist Brunt-Vaisala
      ! frequency.
      T_in_K = thlm2T_in_K( thlm, exner, rcm )
      T_in_K_zm = zt2zm( T_in_K )
      rsat = sat_mixrat_liq( p_in_Pa, T_in_K )
      rsat_zm = zt2zm( rsat )
      ddzt_rsat = ddzt( rsat )
      thm = thlm + Lv/(Cp*exner) * rcm
      thm_zm = zt2zm( thm )
      ddzt_thm = ddzt( thm )
      ddzt_rtm = ddzt( rtm )
    end if

    do k=1, gr%nz

      if ( .not. l_brunt_vaisala_freq_moist ) then

        ! Dry Brunt-Vaisala frequency
        if ( l_use_thvm_in_bv_freq ) then
          brunt_vaisala_freq_sqd(k) = ( grav / thvm_zm(k) ) * ddzt_thvm(k)
        else
          brunt_vaisala_freq_sqd(k) = ( grav / T0 ) * ddzt_thlm(k)
        end if

      else ! l_brunt_vaisala_freq_moist

        ! In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and Klemp (1982)
        brunt_vaisala_freq_sqd(k) = &
          grav * ( ((one + Lv*rsat_zm(k) / (Rd*T_in_K_zm(k))) / &
            (one + ep*(Lv**2)*rsat_zm(k)/(Cp*Rd*T_in_K_zm(k)**2))) * &
            ( (one/thm_zm(k) * ddzt_thm(k)) + (Lv/(Cp*T_in_K_zm(k)))*ddzt_rsat(k)) - &
            ddzt_rtm(k) )

      end if ! .not. l_brunt_vaisala_freq_moist

    end do ! k=1, gr%nz

    return
  end function calc_brunt_vaisala_freq_sqd

!===============================================================================
  function compute_Cx_fnc_Richardson( thlm, um, vm, em, Lscale, exner, rtm, &
                                      rcm, p_in_Pa, cloud_frac, thvm, rho_ds_zm ) &
    result( Cx_fnc_Richardson )

  ! Description:
  !   Compute Cx as a function of the Richardson number

  ! References:
  !   cam:ticket:59
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd  ! Konstant

    use grid_class, only: &
      gr,   & ! Variable
      ddzt, & ! Procedure(s)
      zt2zm

    use constants_clubb, only: &
      one_fourth, &     ! Constant(s)
      one_third,  &
      one,        &
      five

    use interpolation, only: &
      linear_interp_factor ! Procedure

    use stats_variables, only: &
      iRichardson_num, &    ! Variable(s)
      ibrunt_vaisala_freq_sqd, &
      ishear_sqd, &
      stats_zm,       &
      l_stats_samp

    use stats_type_utilities, only: &
      stat_update_var      ! Procedure

    implicit none

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      Richardson_num_divisor_threshold = 1.0e-6_core_rknd, &
      Richardson_num_min = one_fourth, &
      Richardson_num_max = 400._core_rknd,       &
      Cx_min            = one_third,   &
      Cx_max            = 0.95_core_rknd,         &
      Cx_fnc_Richardson_below_ground_value = one

    logical, parameter :: &
      l_Cx_fnc_Richardson_vert_avg = .false.,& ! Vertically average Cx_fnc_Richardson over a
                                        !  distance of Lscale
      l_Richardson_vert_avg = .true. , & ! Vertically average Richardson_num over a
                                         !  distance of Lscale
      l_use_shear_turb_freq_sqd = .false.! Use turb_freq_sqd and shear_sqd in denominator of
                                         !  Richardson_num

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      thlm,    & ! th_l (liquid water potential temperature)      [K]
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      em,      & ! Turbulent Kinetic Energy (TKE)                 [m^2/s^2]
      Lscale,  & ! Turbulent mixing length                        [m]
      exner,   & ! Exner function                                 [-]
      rtm,     & ! total water mixing ratio, r_t                  [kg/kg]
      rcm,     & ! cloud water mixing ratio, r_c                  [kg/kg]
      p_in_Pa, & ! Air pressure                                   [Pa]
      cloud_frac, & ! Cloud fraction                              [-]
      thvm,    & ! Virtual potential temperature                  [K]
      rho_ds_zm  ! Dry static density on momentum levels          [kg/m^3]


    ! Output Variable
    real( kind = core_rknd), dimension(gr%nz) :: &
      Cx_fnc_Richardson

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd, &
      Richardson_num, &
      dum_dz, dvm_dz, &
      shear_sqd, &
      turb_freq_sqd, &
      Lscale_zm

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    brunt_vaisala_freq_sqd = calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, &
                                                          cloud_frac, thvm )

    ! Statistics sampling
    if ( l_stats_samp ) then

      ! NOTE: This is a kludgy place to sample brunt_vaisala_freq_sqd, because
      ! it is used in multiple places, and depending on CLUBB parameters, it
      ! could be computed in another place and not here. In the future, we
      ! should compute brunt_vaisala_freq_sqd once, and pass it around
      ! everywhere. This will save on computational expense as well.
      call stat_update_var( ibrunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd, stats_zm )

    end if ! l_stats_samp

    Lscale_zm = zt2zm( Lscale )

    if ( l_use_shear_turb_freq_sqd ) then
      ! Calculate shear_sqd
      dum_dz = ddzt( um )
      dvm_dz = ddzt( vm )
      shear_sqd = dum_dz**2 + dvm_dz**2

      turb_freq_sqd = em / Lscale_zm**2
      Richardson_num = brunt_vaisala_freq_sqd / max( shear_sqd, turb_freq_sqd, &
                                                     Richardson_num_divisor_threshold )

      if ( l_stats_samp ) &
        call stat_update_var( ishear_sqd, shear_sqd, stats_zm )
    else
      Richardson_num = brunt_vaisala_freq_sqd / Richardson_num_divisor_threshold
    end if

    if ( l_Richardson_vert_avg ) then
      ! Clip below-min values of Richardson_num
      Richardson_num = max( Richardson_num, Richardson_num_min )

      Richardson_num = Lscale_width_vert_avg( Richardson_num, Lscale_zm, rho_ds_zm, &
                                              Richardson_num_max )
    end if

    ! Cx_fnc_Richardson is interpolated based on the value of Richardson_num
    where ( Richardson_num <= Richardson_num_min )
      Cx_fnc_Richardson = Cx_min
    else where ( Richardson_num >= Richardson_num_max )
      Cx_fnc_Richardson = Cx_max
    else where
      ! Linear interpolation
      Cx_fnc_Richardson = &
        linear_interp_factor( (Richardson_num-Richardson_num_min) / &
                              (Richardson_num_max-Richardson_num_min), Cx_max, Cx_min )
    end where

    if ( l_Cx_fnc_Richardson_vert_avg ) then
      Cx_fnc_Richardson = Lscale_width_vert_avg( Cx_fnc_Richardson, Lscale_zm, rho_ds_zm, &
                                                 Cx_fnc_Richardson_below_ground_value )
    end if

    ! On some compilers, roundoff error can result in Cx_fnc_Richardson being
    ! slightly outside the range [0,1]. Thus, it is clipped here.
    Cx_fnc_Richardson = max( 0.0_core_rknd, min( 1.0_core_rknd, Cx_fnc_Richardson ) )

    ! Stats sampling
    if ( l_stats_samp ) then
      call stat_update_var( iRichardson_num, Richardson_num, stats_zm )
    end if

  end function compute_Cx_fnc_Richardson
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  function Lscale_width_vert_avg( var_profile, Lscale_zm, rho_ds_zm, var_below_ground_value )

  ! Description:
  !   Averages a profile with a running mean of width Lscale_zm

  ! References:
  !   cam:ticket:59

    use clubb_precision, only: &
      core_rknd ! Precision

    use grid_class, only: &
      gr ! Variable

    use fill_holes, only: &
      vertical_avg ! Procedure

    use constants_clubb, only: &
      one_half

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      var_profile, &      ! Profile on momentum levels
      Lscale_zm, &        ! Lscale on momentum levels
      rho_ds_zm           ! Dry static energy on momentum levels!

    real( kind = core_rknd ), intent(in) :: &
      var_below_ground_value ! Value to use below ground

    ! Result Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Lscale_width_vert_avg ! Vertically averaged profile (on momentum levels)

    ! Local Variables
    integer :: k, k_avg_lower, k_avg_upper, k_inner_loop

    real( kind = core_rknd ), dimension(gr%nz) :: &
      one_half_avg_width

    real( kind = core_rknd ), dimension(:), allocatable :: &
      rho_ds_zm_virtual, &
      var_profile_virtual, &
      invrs_dzm_virtual

    integer :: n_virtual_levels, n_below_ground_levels

  !----------------------------------------------------------------------
    !----- Begin Code -----
    one_half_avg_width = max( Lscale_zm, 500.0_core_rknd )

    outer_vert_loop: do k=1, gr%nz

      !-----------------------------------------------------------------------
      ! Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
      !-----------------------------------------------------------------------

      k_avg_upper = k

      inner_vert_loop_upward: do k_inner_loop=k+1, gr%nz
        if ( gr%zm(k_inner_loop) - gr%zm(k) <= one_half_avg_width(k) ) then
          ! Include this height level in the average.
          k_avg_upper = k_inner_loop
        else
          ! Do not include this level in the average. No point in searching further.
          exit
        end if
      end do inner_vert_loop_upward

      k_avg_lower = k
      inner_vert_loop_downward: do k_inner_loop=k-1, 1, -1
        if ( gr%zm(k) - gr%zm(k_inner_loop) <= one_half_avg_width(k) ) then
          ! Include this height level in the average.
          k_avg_lower = k_inner_loop
        else
          ! Do not include this level in the average. No point in searching further.
          exit
        end if
      end do inner_vert_loop_downward

      ! Compute the number of levels below ground to include.
      if ( k_avg_lower > 1 ) then
        ! k=1, the lowest "real" level, is not included in the average, so no
        ! below-ground levels should be included.
        n_below_ground_levels = 0
      else
        ! The number of below-ground levels included is equal to the distance
        ! below the lowest level spanned by one_half_avg_width(k)
        ! divided by the distance between vertical levels below ground; the
        ! latter is assumed to be the same as the distance between the first and
        ! second vertical levels.
        n_below_ground_levels = int( ( one_half_avg_width(k)-(gr%zm(k)-gr%zm(1)) ) / &
                                        ( gr%zm(2)-gr%zm(1) ) )
      end if

      ! Prepare the virtual levels!
      n_virtual_levels = k_avg_upper-k_avg_lower+n_below_ground_levels+1
      allocate( rho_ds_zm_virtual(n_virtual_levels), var_profile_virtual(n_virtual_levels), &
                invrs_dzm_virtual(n_virtual_levels) )

      ! All vertical levels have rho_ds_zm and invrs_dzm_virtual equal to the
      ! values at k=1. The value of var_profile at k=1 is given as an argument
      ! to this function.
      if ( n_below_ground_levels > 0 ) then
        rho_ds_zm_virtual(1:n_below_ground_levels) = rho_ds_zm(1)
        var_profile_virtual(1:n_below_ground_levels) = var_below_ground_value
        invrs_dzm_virtual(1:n_below_ground_levels) = gr%invrs_dzm(1)
      end if

      ! Set up the above-ground virtual levels.
      rho_ds_zm_virtual(n_below_ground_levels+1:n_virtual_levels) = &
        rho_ds_zm(k_avg_lower:k_avg_upper)
      var_profile_virtual(n_below_ground_levels+1:n_virtual_levels) = &
        var_profile(k_avg_lower:k_avg_upper)
      invrs_dzm_virtual(n_below_ground_levels+1:n_virtual_levels) = &
        gr%invrs_dzm(k_avg_lower:k_avg_upper)

      ! Finally, compute the average.
      Lscale_width_vert_avg(k) = vertical_avg( n_virtual_levels, rho_ds_zm_virtual, &
                                               var_profile_virtual, invrs_dzm_virtual )

      deallocate( rho_ds_zm_virtual, var_profile_virtual, invrs_dzm_virtual )

    end do outer_vert_loop

    return
  end function Lscale_width_vert_avg

end module advance_helper_module
