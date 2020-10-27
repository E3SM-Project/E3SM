!-------------------------------------------------------------------------
! $Id: advance_helper_module.F90 8769 2018-08-11 22:54:36Z vlarson@uwm.edu $
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
    compute_Cx_fnc_Richardson, &
    term_wp2_splat, term_wp3_splat

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
  function calc_stability_correction( thlm, Lscale, em, exner, rtm, rcm, &
                                      p_in_Pa, thvm ) &
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
      thvm               ! Virtual potential temperature             [K]

    ! Result
    real( kind = core_rknd ), dimension(gr%nz) :: &
      stability_correction

    real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd, & !  []
      lambda0_stability

    !------------ Begin Code --------------
    call calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, thvm, &
                                      brunt_vaisala_freq_sqd )

    lambda0_stability = merge( lambda0_stability_coef, zero, brunt_vaisala_freq_sqd > zero )

    stability_correction = 1.0_core_rknd &
      + min( lambda0_stability * brunt_vaisala_freq_sqd * zt2zm( Lscale )**2 / em, 3.0_core_rknd )

    return
  end function calc_stability_correction

  !===============================================================================
  subroutine calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, thvm, &
                                           brunt_vaisala_freq_sqd )

  ! Description:
  !   Calculate the Brunt-Vaisala frequency squared, N^2.

  ! References:
  !   ?
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Konstant

    use constants_clubb, only: &
      grav, & ! Constant(s)
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
      thvm        ! Virtual potential temperature      [K]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
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

    if ( .not. l_brunt_vaisala_freq_moist ) then

        ! Dry Brunt-Vaisala frequency
        if ( l_use_thvm_in_bv_freq ) then

          brunt_vaisala_freq_sqd(:) = ( grav / thvm_zm(:) ) * ddzt_thvm(:)

        else

          brunt_vaisala_freq_sqd(:) = ( grav / T0 ) * ddzt_thlm(:)
    
        end if

    else ! l_brunt_vaisala_freq_moist

        T_in_K = thlm2T_in_K( thlm, exner, rcm )
        T_in_K_zm = zt2zm( T_in_K )
        rsat = sat_mixrat_liq( p_in_Pa, T_in_K )
        rsat_zm = zt2zm( rsat )
        ddzt_rsat = ddzt( rsat )
        thm = thlm + Lv/(Cp*exner) * rcm
        thm_zm = zt2zm( thm )
        ddzt_thm = ddzt( thm )
        ddzt_rtm = ddzt( rtm )

        do k=1, gr%nz

            ! In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and Klemp (1982)
            brunt_vaisala_freq_sqd(k) = &
              grav * ( ((one + Lv*rsat_zm(k) / (Rd*T_in_K_zm(k))) / &
                (one + ep*(Lv**2)*rsat_zm(k)/(Cp*Rd*T_in_K_zm(k)**2))) * &
                ( (one/thm_zm(k) * ddzt_thm(k)) + (Lv/(Cp*T_in_K_zm(k)))*ddzt_rsat(k)) - &
                ddzt_rtm(k) )

        end do ! k=1, gr%nz

    end if ! .not. l_brunt_vaisala_freq_moist

    return

  end subroutine calc_brunt_vaisala_freq_sqd

!===============================================================================
  subroutine compute_Cx_fnc_Richardson( thlm, um, vm, em, Lscale, exner, rtm, &
                                        rcm, p_in_Pa, thvm, rho_ds_zm, &
                                        Cx_fnc_Richardson )

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
      one

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
      thvm,    & ! Virtual potential temperature                  [K]
      rho_ds_zm  ! Dry static density on momentum levels          [kg/m^3]


    ! Output Variable
    real( kind = core_rknd), dimension(gr%nz), intent(out) :: &
      Cx_fnc_Richardson

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      brunt_vaisala_freq_sqd, &
      Richardson_num, &
      dum_dz, dvm_dz, &
      shear_sqd, &
      turb_freq_sqd, &
      Lscale_zm

    real ( kind = core_rknd ), dimension(gr%nz) :: &
        invrs_min_max_diff, &
        invrs_num_div_thresh

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    call calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, thvm, &
                                      brunt_vaisala_freq_sqd )

    invrs_min_max_diff = 1.0_core_rknd / ( Richardson_num_max - Richardson_num_min )
    invrs_num_div_thresh = 1.0_core_rknd / Richardson_num_divisor_threshold

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
      Richardson_num = brunt_vaisala_freq_sqd * invrs_num_div_thresh
    end if

    if ( l_Richardson_vert_avg ) then
      ! Clip below-min values of Richardson_num
      Richardson_num = max( Richardson_num, Richardson_num_min )

      Richardson_num = Lscale_width_vert_avg( Richardson_num, Lscale_zm, rho_ds_zm, &
                                              Richardson_num_max )
    end if

    ! Cx_fnc_Richardson is interpolated based on the value of Richardson_num
    Cx_fnc_Richardson = linear_interp_factor( &
                        (Richardson_num-Richardson_num_min) * invrs_min_max_diff , Cx_max, Cx_min )

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

  end subroutine compute_Cx_fnc_Richardson
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
    integer :: &
        k, i,        & ! Loop variable
        k_avg_lower, &
        k_avg_upper

    real( kind = core_rknd ), dimension(gr%nz) :: &
      one_half_avg_width, &
      numer_terms, &
      denom_terms

    integer :: n_below_ground_levels

    real( kind = core_rknd ) :: & 
      numer_integral, & ! Integral in the numerator (see description)
      denom_integral    ! Integral in the denominator (see description)

  !----------------------------------------------------------------------

    !----- Begin Code -----

    one_half_avg_width = max( Lscale_zm, 500.0_core_rknd )

    ! Pre calculate numerator and denominator terms
    do k=1, gr%nz
        numer_terms(k) = rho_ds_zm(k) * gr%dzm(k) * var_profile(k)
        denom_terms(k) = rho_ds_zm(k) * gr%dzm(k)
    end do

    k_avg_upper = 2
    k_avg_lower = 1

    ! For every grid level
    do k=1, gr%nz

        !-----------------------------------------------------------------------
        ! Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
        ! 
        !     k_avg_upper and k_avg_lower are saved each loop iteration, this 
        !     improves computational efficiency since their values are likely
        !     within one or two grid levels of where they were last found to 
        !     be. This is because one_half_avg_width does not change drastically
        !     from one grid level to the next. Thus, less searching is required
        !     by allowing the search to start at a value that is close to the
        !     desired value and allowing each value to increment or decrement 
        !     as needed. 
        !-----------------------------------------------------------------------


        ! Determine if k_avg_upper needs to increment or decrement
        if ( gr%zm(k_avg_upper) - gr%zm(k) > one_half_avg_width(k) ) then

            ! k_avg_upper is too large, decrement it
            do while ( gr%zm(k_avg_upper) - gr%zm(k) > one_half_avg_width(k) )
                k_avg_upper = k_avg_upper - 1
            end do

        elseif ( k_avg_upper < gr%nz ) then

            ! k_avg_upper is too small, increment it
            do while ( gr%zm(k_avg_upper+1) - gr%zm(k) <= one_half_avg_width(k) )

                k_avg_upper = k_avg_upper + 1

                if ( k_avg_upper == gr%nz ) exit

            end do

        end if


        ! Determine if k_avg_lower needs to increment or decrement
        if ( gr%zm(k) - gr%zm(k_avg_lower) > one_half_avg_width(k) ) then

            ! k_avg_lower is too small, increment it
            do while ( gr%zm(k) - gr%zm(k_avg_lower) > one_half_avg_width(k) )

                k_avg_lower = k_avg_lower + 1

            end do

        elseif ( k_avg_lower > 1 ) then

            ! k_avg_lower is too large, decrement it
            do while ( gr%zm(k) - gr%zm(k_avg_lower-1) <= one_half_avg_width(k) )

                k_avg_lower = k_avg_lower - 1

                if ( k_avg_lower == 1 ) exit

            end do 

        end if


        ! Compute the number of levels below ground to include.
        if ( k_avg_lower > 1 ) then

            ! k=1, the lowest "real" level, is not included in the average, so no
            ! below-ground levels should be included.
            n_below_ground_levels = 0

            numer_integral = 0.0_core_rknd
            denom_integral = 0.0_core_rknd

        else

            ! The number of below-ground levels included is equal to the distance
            ! below the lowest level spanned by one_half_avg_width(k)
            ! divided by the distance between vertical levels below ground; the
            ! latter is assumed to be the same as the distance between the first and
            ! second vertical levels.
            n_below_ground_levels = int( ( one_half_avg_width(k)-(gr%zm(k)-gr%zm(1)) ) / &
                                        ( gr%zm(2)-gr%zm(1) ) )

            numer_integral = n_below_ground_levels * denom_terms(1) * var_below_ground_value
            denom_integral = n_below_ground_levels * denom_terms(1)

        end if

            
        ! Add numerator and denominator terms for all above-ground levels
        do i = k_avg_lower, k_avg_upper

            numer_integral = numer_integral + numer_terms(i)
            denom_integral = denom_integral + denom_terms(i)

        end do

        Lscale_width_vert_avg(k) = numer_integral / denom_integral

    end do

    return

  end function Lscale_width_vert_avg

 !============================================================================
  subroutine term_wp2_splat( C_wp2_splat, nz, dt, wp2, wp2_zt, tau_zm, &
                             wp2_splat )


  ! Description:
  !   This subroutine computes the (negative) tendency of wp2 due
  !   to "splatting" of eddies, e.g., near the ground or a Sc inversion.
  !   term_splat is intended to be added to the right-hand side of 
  !   the wp2 equation, and -0.5*term_splat is intended to be added to each 
  !   of the up2 and vp2 equations.  The functional form of term splat is
  !
  !   term_splat \propto - w'2 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2

    ! Included Modules
    use grid_class, only: &
      ddzt   ! Procedure(s)

    use clubb_precision, only: & 
      core_rknd 

    use constants_clubb, only: &
      five  ! Constant(s)

    implicit none 

    ! Input Variables
    integer, intent(in) :: & 
      nz          ! Number of vertical levels                    [-]

    real( kind = core_rknd ), intent(in) :: & 
      C_wp2_splat, &          ! Tuning parameter                    [-]
      dt                      ! CLUBB computational time step       [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: & 
      wp2,     &  ! Variance of vertical velocity on the momentum grid [m^2/s^2]
      wp2_zt,  &  ! Variance of vertical velocity on the thermodynamic grid [m^2/s^2]
      tau_zm     ! Turbulent time scale on the momentum grid               [s]

    ! Output Variable
    real( kind = core_rknd ), dimension(nz), intent(out) :: & 
      wp2_splat     ! Tendency of <w'^2> due to splatting of eddies (on zm grid) [m^2/s^3]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: & 
      d_sqrt_wp2_dz    ! d/dz( sqrt( w'2 ) )                 [1/s]

    ! ---- Begin Code ----

    d_sqrt_wp2_dz = ddzt( sqrt( wp2_zt ) )
    ! The splatting term is clipped so that the incremental change doesn't exceed 5 times the
    !   value of wp2 itself.  This prevents spikes in wp2 from being propagated to up2 and vp2.
    !   However, it does introduce undesired dependence on the time step.
    !   Someday we may wish to treat this term using a semi-implicit discretization.
    wp2_splat = - wp2 * min( five/dt, C_wp2_splat * tau_zm * d_sqrt_wp2_dz**2 )
    !wp2_splat = - C_wp2_splat * wp2 * 900.0_core_rknd * d_sqrt_wp2_dz**2

  end subroutine term_wp2_splat

  !============================================================================
  subroutine term_wp3_splat( C_wp2_splat, nz, dt, wp2, wp3, tau_zt, &
                             wp3_splat )

  ! Description:
  !   This subroutine computes the damping of wp3 due
  !   to "splatting" of eddies, e.g., as they approach the ground or a Sc inversion.
  !   term_wp3_splat is intended to be added to the right-hand side of 
  !   the wp3 equation.  The functional form of wp3_splat is
  !
  !   wp3_splat \propto - w'3 * (turbulent time scale) * ( d/dz( sqrt(w'2) ) )^2
  !
  !   If the coefficient on wp3_splat is at least 1.5 times greater than the 
  !   coefficient on wp2_splat, then skewness will be damped, promoting
  !   more stratiform layers.

    ! Included Modules
    use grid_class, only: &
      ddzm   ! Procedure(s)

    use clubb_precision, only: & 
      core_rknd 

    use constants_clubb, only: &
      three, &  ! Constant(s)
      five

    implicit none 

    ! Input Variables
    integer, intent(in) :: & 
      nz          ! Number of vertical levels                    [-]

    real( kind = core_rknd ), intent(in) :: & 
      C_wp2_splat, &          ! Tuning parameter                    [-]
      dt                      ! CLUBB computational time step       [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: & 
      wp2,     &  ! Variance of vertical velocity on the momentum grid [m^2/s^2]
      wp3,     &  ! Third moment of vertical velocity on the momentum grid [m^3/s^3]
      tau_zt      ! Turbulent time scale on the thermal grid               [s]

    ! Output Variable
    real( kind = core_rknd ), dimension(nz), intent(out) :: & 
      wp3_splat     ! Tendency of <w'^3> due to splatting of eddies (on zt grid) [m^3/s^4]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: & 
      d_sqrt_wp2_dz    ! d/dz( sqrt( w'2 ) )                 [1/s]

    ! ---- Begin Code ----

    d_sqrt_wp2_dz = ddzm( sqrt( wp2 ) )
    ! The splatting term is clipped so that the incremental change doesn't exceed 5 times the
    !   value of wp2 itself.  This prevents spikes in wp2 from being propagated to up2 and vp2.
    !   However, it does introduce undesired dependence on the time step.
    !   Someday we may wish to treat this term using a semi-implicit discretization.
    wp3_splat = - wp3 * min( five/dt, three * C_wp2_splat * tau_zt * d_sqrt_wp2_dz**2 )
    !wp3_splat = - three * C_wp2_splat * wp3 * 900._core_rknd * d_sqrt_wp2_dz**2

  end subroutine term_wp3_splat

end module advance_helper_module
