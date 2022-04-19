!===============================================================================
module calc_pressure

  implicit none

  public :: update_pressure, & ! Procedure(s)
            init_pressure,   &
            calculate_thvm

  private ! default scope

  contains

  !=============================================================================
  subroutine update_pressure( gr, thlm, rtm, rcm, rho_ds_zt, thv_ds_zt, &
                              p_in_Pa, exner, &
                              p_in_Pa_zm, exner_zm )

    ! Description:
    ! Updates pressure according to the hydrostatic approximation.  Combining
    ! the moist equation of state and the hydrostatic approximation, the change
    ! of pressure with respect to height can be calculated based on theta_v,
    ! such that:
    !
    ! dp/dz = - p * grav / ( Rd * theta_v * exner );
    !
    ! where exner = ( p / p0 )^(Rd/Cp);
    !
    ! and where p0 is a reference pressure of 100000 Pa.
    !
    ! The integral equation is set up to integrate over p on the left-hand side
    ! and integrate over z on the right-hand side.  The equation is:
    !
    ! INT(p1:p2) p^(Rd/Cp-1) dp
    ! = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.
    !
    ! The value of mean theta_v (thvm) is calculated at each thermodynamic grid
    ! level, and linear interpolation is used in the integral equation for all
    ! altitudes in-between successive thermodynamic levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! The integrals are solved, and the results for pressure can be rewritten
    ! in terms of exner, such that:
    !
    ! exner2 - exner1
    !   | - ( grav / Cp )
    !   |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    ! = | where thvm2 /= thvm1;
    !   |
    !   | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    ! 
    ! The value of pressure (exner) can be calculated using the above equation
    ! at all vertical levels once the value of pressure is known at one level.
    !
    ! The model fields are usually least variant over time at the top of the
    ! model, so the value of pressure is calculated at the uppermost
    ! thermodynamic level, and then the values of pressure (exner) can be
    ! calculated at all other vertical levels (both thermodynamic levels and
    ! momentum levels).  The following equation is used to calculate pressure
    ! at the uppermost thermodynamic level:
    !
    ! p = ( ( rho_dry * Rd / p0^(Rd/Cp) )
    !       * thm * ( 1 + (Rv/Rd) * rvm ) )^(1/(1-Rd/Cp));
    !
    ! where rho_dry is the density of dry air, where rvm = rtm - rcm, and where:
    !
    ! thm = thlm + ( Lv / ( Cp * exner|_old ) ) * rcm;
    !
    ! and exner|_old is exner from the previous timestep.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one,   & ! 1
        Rd,    & ! Gas constant for dry air                    [J/(kg K)]
        Lv,    & ! Latent heat of vaporization                 [J/kg]
        Cp,    & ! Specific heat of dry air                    [J/(kg K)]
        kappa, & ! Rd/Cp                                       [-]
        ep2,   & ! Rv/Rd (Rv is gas constant for water vapor)  [-]
        p0,    & ! Reference pressure of 100000 Pa             [Pa]
        grav     ! Acceleration of gravity (9.81 m/s^2)        [m/s^2]

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      thlm,      & ! Mean liquid water potential temperature         [K]
      rtm,       & ! Mean total water mixing ratio                   [kg/kg]
      rcm,       & ! Mean cloud water mixing ratio                   [kg/kg]
      rho_ds_zt, & ! Dry, state base-state density (thermo. levels)  [kg/m^3]
      thv_ds_zt    ! Reference theta_v on thermodynamic levels       [K]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      p_in_Pa,    & ! Pressure (thermodynamic levels)        [Pa]
      exner         ! Exner function (thermodynamic levels)  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      p_in_Pa_zm, & ! Pressure on momentum levels            [Pa]
      exner_zm      ! Exner function on momentum levels      [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      thvm,    & ! Mean theta_v (thermodynamic levels)                [K]
      thvm_zm    ! Mean theta_v interpolated to momentum grid levels  [K]

    real( kind = core_rknd ) :: &
      thm_nz, & ! Theta at the uppermost thermodynamic grid level    [K]
      rvm_nz    ! Water vapor mixing ratio; uppermost thermo. level  [kg/kg]

    real( kind = core_rknd ), parameter :: &
      g_ov_Cp = grav / Cp, &     ! g / Cp  [K/m]
      invrs_kappa = one / kappa  ! 1 / kappa

    ! Flag to calculate pressure and exner on momentum levels.  Otherwise,
    ! linear interpolation of exner will be used.
    logical, parameter :: &
      l_calc_p_exner_m_levs = .true.

    integer :: k  ! Vertical level index


    ! Calculate thvm on thermodynamic grid levels.
    thvm = calculate_thvm( thlm, rtm, rcm, exner, thv_ds_zt )

    ! Interpolate thvm to momentum grid levels.
    thvm_zm = zt2zm( gr, thvm )

    ! Calculate mean theta (thm) at the uppermost thermodynamic vertical grid
    ! level.
    thm_nz = thlm(gr%nz) + ( Lv / ( Cp * exner(gr%nz) ) ) * rcm(gr%nz)

    ! Calculate mean water vapor mixing ratio (rvm) at the uppermost
    ! thermodynamic vertical grid level.
    rvm_nz = rtm(gr%nz) - rcm(gr%nz)

    ! Update pressure at the uppermost thermodynamic grid level.
    p_in_Pa(gr%nz) &
    = ( ( rho_ds_zt(gr%nz) * Rd / p0**kappa ) &
        * thm_nz * ( one + ep2 * rvm_nz ) )**(one/(one-kappa))

    ! Calculate exner at the uppermost thermodynamic grid level.
    exner(gr%nz) = ( p_in_Pa(gr%nz) / p0 )**kappa

    ! Calculate exner at the uppermost momentum grid level, which is located
    ! above the uppermost thermodynamic grid level.
    ! exner2
    ! = exner1
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   - | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    if ( l_calc_p_exner_m_levs ) then

       if ( abs( thvm(gr%nz) - thvm_zm(gr%nz) ) &
            > epsilon( thvm ) * thvm(gr%nz) ) then

          exner_zm(gr%nz) &
          = exner(gr%nz) &
            - g_ov_Cp * ( gr%zm(gr%nz) - gr%zt(gr%nz) ) &
              / ( thvm_zm(gr%nz) - thvm(gr%nz) ) &
              * log( thvm_zm(gr%nz) / thvm(gr%nz) )

       else ! thvm(k+1) = thvm_zm(k)

          exner_zm(gr%nz) &
          = exner(gr%nz) &
            - g_ov_Cp * ( gr%zm(gr%nz) - gr%zt(gr%nz) ) / thvm_zm(gr%nz)

       endif

    else ! .not. l_calc_p_exner_m_levs

       ! Interpolate exner to momentum levels
       exner_zm(gr%nz) = zt2zm( gr, exner, gr%nz )

    endif ! l_calc_p_exner_m_levs

    ! Calculate pressure on the uppermost momentum level.
    p_in_Pa_zm(gr%nz) = p0 * exner_zm(gr%nz)**invrs_kappa

    ! Calculate exner at all other thermodynamic and momentum grid levels,
    ! which are all located below the uppermost thermodynamic grid level.
    ! exner1
    ! = exner2
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   + | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    do k = gr%nz-1, 2, -1

       ! Calculate exner on thermodynamic levels.
       if ( abs( thvm(k+1) - thvm(k) ) > epsilon( thvm ) * thvm(k+1) ) then

          exner(k) &
          = exner(k+1) &
            + g_ov_Cp * ( gr%zt(k+1) - gr%zt(k) ) / ( thvm(k+1) - thvm(k) ) &
              * log( thvm(k+1) / thvm(k) )

       else ! thvm(k+1) = thvm(k)

          exner(k) = exner(k+1) + g_ov_Cp * ( gr%zt(k+1) - gr%zt(k) ) / thvm(k)

       endif

       if ( l_calc_p_exner_m_levs ) then

          ! Calculate exner on momentum levels.
          if ( abs( thvm(k+1) - thvm_zm(k) ) &
               > epsilon( thvm ) * thvm(k+1) ) then

             exner_zm(k) &
             = exner(k+1) &
               + g_ov_Cp * ( gr%zt(k+1) - gr%zm(k) ) &
                 / ( thvm(k+1) - thvm_zm(k) ) &
                 * log( thvm(k+1) / thvm_zm(k) )

          else ! thvm(k+1) = thvm_zm(k)

             exner_zm(k) &
             = exner(k+1) + g_ov_Cp * ( gr%zt(k+1) - gr%zm(k) ) / thvm_zm(k)

          endif

       else ! .not. l_calc_p_exner_m_levs

          ! Interpolate exner to momentum levels
          exner_zm(k) = zt2zm( gr, exner, k )

       endif ! l_calc_p_exner_m_levs

    enddo ! k = gr%nz-1, 2, -1

#ifdef MKL
    ! MKL VML functions available. vdpowx(n,a,b,y) computes a(1:n)^b = y(1:n)
    ! This temporarily store exner(_zm)**invrs_kappa in p_in_Pa(_zm), before
    ! multiplying p_in_Pa(_zm) by p0 to complete the calculation.

    ! Calculate pressure on thermodynamic levels
    call vdpowx( gr%nz-2, exner(2:gr%nz-1), invrs_kappa, p_in_Pa(2:gr%nz-1) )
    p_in_Pa(2:gr%nz-1) =  p_in_Pa(2:gr%nz-1) * p0

    ! Calculate pressure on momentum levels.
    call vdpowx( gr%nz-2, exner_zm(2:gr%nz-1), invrs_kappa, p_in_Pa_zm(2:gr%nz-1) )
    p_in_Pa_zm(2:gr%nz-1) =  p_in_Pa_zm(2:gr%nz-1) * p0
#else
    ! Calculate pressure on thermodynamic levels.
    p_in_Pa(2:gr%nz-1) = p0 * exner(2:gr%nz-1)**invrs_kappa

    ! Calculate pressure on momentum levels.
    p_in_Pa_zm(2:gr%nz-1) = p0 * exner_zm(2:gr%nz-1)**invrs_kappa
#endif

    ! Calculate exner the model lower boundary or surface.
    ! exner1
    ! = exner2
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   + | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    if ( abs( thvm(2) - thvm_zm(1) ) > epsilon( thvm ) * thvm(2) ) then

       exner_zm(1) &
       = exner(2) &
         + g_ov_Cp * ( gr%zt(2) - gr%zm(1) ) / ( thvm(2) - thvm_zm(1) ) &
           * log( thvm(2) / thvm_zm(1) )

     else ! thvm(k+1) = thvm_zm(k)

       exner_zm(1) &
       = exner(2) + g_ov_Cp * ( gr%zt(2) - gr%zm(1) ) / thvm_zm(1)

    endif

    ! Calculate pressure at the model lower boundary or surface.
    p_in_Pa_zm(1) = p0 * exner_zm(1)**invrs_kappa

    ! For the lowest thermodynamic level, which is below the model lower
    ! boundary, set pressure and exner to the pressure and exner found at the
    ! model lower boundary.
    p_in_Pa(1) = p_in_Pa_zm(1)
    exner(1) = exner_zm(1)


    return

  end subroutine update_pressure

  !=============================================================================
  subroutine init_pressure( gr, thvm, p_sfc, &
                            p_in_Pa, exner, p_in_Pa_zm, exner_zm )

    ! Description:
    ! Calculates the initial pressure according to the hydrostatic
    ! approximation.  Combining the moist equation of state and the hydrostatic
    ! approximation, the change of pressure with respect to height can be
    ! calculated based on theta_v, such that:
    !
    ! dp/dz = - p * grav / ( Rd * theta_v * exner );
    !
    ! where exner = ( p / p0 )^(Rd/Cp);
    !
    ! and where p0 is a reference pressure of 100000 Pa.
    !
    ! The integral equation is set up to integrate over p on the left-hand side
    ! and integrate over z on the right-hand side.  The equation is:
    !
    ! INT(p1:p2) p^(Rd/Cp-1) dp
    ! = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.
    !
    ! The value of mean theta_v (thvm) is calculated at each thermodynamic grid
    ! level, and linear interpolation is used in the integral equation for all
    ! altitude in-between successive thermodynamic levels, such that:
    !
    ! thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.
    !
    ! The integrals are solved, and the results for pressure can be rewritten
    ! in terms of exner, such that:
    !
    ! exner2 - exner1
    !   | - ( grav / Cp )
    !   |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    ! = | where thvm2 /= thvm1;
    !   |
    !   | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    ! 
    ! The value of pressure (exner) can be calculated using the above equation
    ! at all vertical levels once the value of pressure is known at one level.
    ! Since the surface pressure is known at the initial time, that allows
    ! pressure to be calculated for the rest of the vertical profile.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid, & ! Type
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one,   & ! 1
        Cp,    & ! Specific heat of dry air                    [J/(kg K)]
        kappa, & ! Rd/Cp                                       [-]
        p0,    & ! Reference pressure of 100000 Pa             [Pa]
        grav     ! Acceleration of gravity (9.81 m/s^2)        [m/s^2]

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      thvm    ! Mean theta_v (thermodynamic levels)                [K]

    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Surface pressure                                   [Pa]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      p_in_Pa,    & ! Pressure (thermodynamic levels)        [Pa]
      exner,      & ! Exner function (thermodynamic levels)  [-]
      p_in_Pa_zm, & ! Pressure on momentum levels            [Pa]
      exner_zm      ! Exner function on momentum levels      [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      thvm_zm    ! Mean theta_v interpolated to momentum grid levels  [K]

    real( kind = core_rknd ), parameter :: &
      g_ov_Cp = grav / Cp  ! g / Cp  [K/m]

    integer :: k  ! Vertical level index


    ! The pressure (and exner) at the lowest momentum level is the surface
    ! pressure (and exner based on the surface pressure).
    p_in_Pa_zm(1) = p_sfc
    exner_zm(1) = ( p_sfc / p0 )**kappa

    ! Set the pressure (and exner) at the lowest thermodynamic level, which is
    ! below the model lower boundary, to their values at the model lower
    ! boundary or surface.
    p_in_Pa(1) = p_in_Pa_zm(1)
    exner(1) = exner_zm(1)

    ! Interpolate theta_v to momentum levels.
    thvm_zm = zt2zm( gr, thvm )

    ! Calculate exner at all other thermodynamic and momentum grid levels.
    ! exner2
    ! = exner1
    !     | ( grav / Cp )
    !     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    !   - | where thvm2 /= thvm1;
    !     |
    !     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).

    ! Calculate exner at thermodynamic level 2 (first thermodynamic grid level
    ! that is above the lower boundary).
    if ( abs( thvm(2) - thvm_zm(1) ) > epsilon( thvm ) * thvm(2) ) then

       exner(2) &
       = exner_zm(1) &
         - g_ov_Cp * ( gr%zt(2) - gr%zm(1) ) / ( thvm(2) - thvm_zm(1) ) &
           * log( thvm(2) / thvm_zm(1) )

    else ! thvm(2) = thvm_zm(1)

       exner(2) = exner_zm(1) - g_ov_Cp * ( gr%zt(2) - gr%zm(1) ) / thvm(2)

    endif

    ! Calculate pressure on thermodynamic levels.
    p_in_Pa(2) = p0 * exner(2)**(one/kappa)

    ! Loop over all other thermodynamic vertical grid levels.
    do k = 3, gr%nz, 1

       ! Calculate exner on thermodynamic levels.
       if ( abs( thvm(k) - thvm(k-1) ) > epsilon( thvm ) * thvm(k) ) then

          exner(k) &
          = exner(k-1) &
            - g_ov_Cp * ( gr%zt(k) - gr%zt(k-1) ) / ( thvm(k) - thvm(k-1) ) &
              * log( thvm(k) / thvm(k-1) )

       else ! thvm(k+1) = thvm(k)

          exner(k) = exner(k-1) - g_ov_Cp * ( gr%zt(k) - gr%zt(k-1) ) / thvm(k)

       endif

       ! Calculate pressure on thermodynamic levels.
       p_in_Pa(k) = p0 * exner(k)**(one/kappa)

    enddo ! k = 2, gr%nz, 1

    ! Loop over all momentum grid levels.
    do k = 2, gr%nz, 1

       ! Calculate exner on momentum levels.
       if ( abs( thvm(k) - thvm_zm(k) ) > epsilon( thvm ) * thvm(k) ) then

          exner_zm(k) &
          = exner(k) &
            - g_ov_Cp * ( gr%zm(k) - gr%zt(k) ) / ( thvm_zm(k) - thvm(k) ) &
              * log( thvm_zm(k) / thvm(k) )

       else ! thvm(k) = thvm_zm(k)

          exner_zm(k) &
          = exner(k) - g_ov_Cp * ( gr%zm(k) - gr%zt(k) ) / thvm_zm(k)

       endif

       ! Calculate pressure on momentum levels.
       p_in_Pa_zm(k) = p0 * exner_zm(k)**(one/kappa)

    enddo ! k = 2, gr%nz, 1


    return

  end subroutine init_pressure

  !=============================================================================
  elemental function calculate_thvm( thlm, rtm, rcm, exner, thv_ds_zt ) &
  result( thvm )

    ! Description:
    ! Calculates mean theta_v based on a linearized approximation to the theta_v
    ! equation.  This equation also includes liquid water loading.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        Lv,  & ! Latent Heat of Vaporizaion    [J/kg]
        Cp,  & ! Specific Heat of Dry Air      [J/(kg K)]
        ep1, & ! Rv/Rd - 1                     [-]
        ep2    ! Rv/Rd                         [-]

    use clubb_precision, only: &
        core_rknd

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thlm,      & ! Mean theta_l (thermodynamic levels)          [K]
      rtm,       & ! Mean total water (thermodynamic levels)      [kg/kg]
      rcm,       & ! Mean cloud water (thermodynamic levels)      [kg/kg]
      exner,     & ! Exner function (thermodynamic levels)        [-]
      thv_ds_zt    ! Reference theta_v on thermodynamic levels    [K]

    ! Return Variable
    real( kind = core_rknd ) :: &
      thvm    ! Mean theta_v (thermodynamic levels)    [K]


    ! Calculate mean theta_v
    thvm = thlm + ep1 * thv_ds_zt * rtm &
           + ( Lv / ( Cp * exner ) - ep2 * thv_ds_zt ) * rcm


    return

  end function calculate_thvm
  !=============================================================================

end module calc_pressure
