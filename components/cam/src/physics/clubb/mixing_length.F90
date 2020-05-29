!-----------------------------------------------------------------------
! $Id: mixing_length.F90 8664 2018-05-10 20:21:35Z huebler@uwm.edu $
!===============================================================================
module mixing_length

  implicit none

  private ! Default Scope

  public :: compute_mixing_length

  contains

  !=============================================================================
  subroutine compute_mixing_length( thvm, thlm, &
                             rtm, em, Lscale_max, p_in_Pa, &
                             exner, thv_ds, mu, l_implemented, &
                             Lscale, Lscale_up, Lscale_down )
    ! Description:
    !   Larson's 5th moist, nonlocal length scale
    ! 
    ! References:
    !   Section 3b ( /Eddy length formulation/ ) of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    ! 
    ! Notes: 
    ! 
    !   The equation for the rate of change of theta_l and r_t of the parcel with
    !   respect to height, due to entrainment, is:
    ! 
    !           d(thl_par)/dz = - mu * ( thl_parcel - thl_environment );
    ! 
    !           d(rt_par)/dz = - mu * ( rt_parcel - rt_environment );
    ! 
    !   where mu is the entrainment rate,
    !   such that:
    ! 
    !           mu = (1/m)*(dm/dz);
    ! 
    !   where m is the mass of the parcel.  The value of mu is set to be a
    !   constant.
    ! 
    !   The differential equations are solved for given the boundary condition 
    !   and given the fact that the value of thl_environment and rt_environment 
    !   are treated as changing linearly for a parcel of air from one grid level
    !   to the next.
    ! 
    !   For the special case where entrainment rate, mu, is set to 0,
    !   thl_parcel and rt_parcel remain constant
    ! 
    ! 
    !   The equation for Lscale_up is:
    ! 
    !       INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
    ! 
    !   and for Lscale_down
    !   
    !       INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
    ! 
    !   where thv_par is theta_v of the parcel, thvm is the mean
    !   environmental value of theta_v, z_i is the altitude that the parcel
    !   started from, and em is the mean value of TKE at
    !   altitude z_i (which gives the parcel its initial boost)
    ! 
    !   The increment of CAPE (convective air potential energy) for any two 
    !   successive vertical levels is:
    !   
    !       Upwards:
    !           CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
    !
    !       Downwards: 
    !           CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz
    ! 
    !   Thus, the derivative of CAPE with respect to height is:
    ! 
    !           dCAPE/dz = g * ( thv_par - thvm ) / thvm.
    ! 
    !   A purely trapezoidal rule is used between levels, and is considered
    !   to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be 
    !   of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
    !   where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )
    ! 
    !   The integral is evaluated to find the CAPE increment between two
    !   successive vertical levels.  The result either adds to or depletes
    !   from the total amount of energy that keeps the parcel ascending/descending.
    ! 
    ! 
    ! IMPORTANT NOTE:
    !   This subroutine has been optimized by adding precalculations, rearranging
    !   equations to avoid divides, and modifying the algorithm entirely.
    !       -Gunther Huebler
    ! 
    !   The algorithm previously used looped over every grid level, following a
    !   a parcel up from its initial grid level to its max. The very nature of this 
    !   algorithm is an N^2 
    !--------------------------------------------------------------------------------

    ! mu = (1/M) dM/dz > 0.  mu=0 for no entrainment.
    ! Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
    ! When mu was fixed, we used the value mu = 6.e-4

    use constants_clubb, only:  &  ! Variable(s)
        Cp,             & ! Dry air specific heat at constant pressure [J/kg/K]
        Rd,             & ! Dry air gas constant                       [J/kg/K]
        ep,             & ! Rd / Rv                                    [-]
        ep1,            & ! (1-ep)/ep                                  [-]
        ep2,            & ! 1/ep                                       [-]
        Lv,             & ! Latent heat of vaporiztion                 [J/kg/K]
        grav,           & ! Gravitational acceleration                 [m/s^2]
        fstderr,        &
        zero_threshold, &
        eps,            &
        em_min

    use parameters_tunable, only:  &  ! Variable(s)
        lmin    ! Minimum value for Lscale                         [m]

    use grid_class, only:  & 
        gr,  & ! Variable(s)
        zm2zt ! Procedure(s)

    use numerical_check, only:  & 
        length_check ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  & 
      zlmin = 0.1_core_rknd, & ! Minimum value for Lscale [m]
      Lscale_sfclyr_depth = 500._core_rknd ! [m]

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      thvm,    & ! Virtual potential temp. on themodynamic level  [K]
      thlm,    & ! Liquid potential temp. on themodynamic level   [K]
      rtm,     & ! Total water mixing ratio on themodynamic level [kg/kg]
      em,      & ! em = 3/2 * w'^2; on momentum level             [m^2/s^2]
      exner,   & ! Exner function on thermodynamic level          [-]
      p_in_Pa, & ! Pressure on thermodynamic level                [Pa]
      thv_ds     ! Dry, base-state theta_v on thermodynamic level [K]
    ! Note:  thv_ds used as a reference theta_l here

    real( kind = core_rknd ), intent(in) :: &
      Lscale_max ! Maximum allowable value for Lscale             [m]

    real( kind = core_rknd ), intent(in) :: &
      mu  ! mu Fractional extrainment rate per unit altitude      [1/m]

    logical, intent(in) :: &
      l_implemented ! Flag for CLUBB being implemented in a larger model

    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  & 
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    ! Local Variables

    integer :: i, j

    real( kind = core_rknd ) :: tke, CAPE_incr

    real( kind = core_rknd ) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    ! Temporary arrays to store calculations to speed runtime
    real( kind = core_rknd ), dimension(gr%nz) :: &
        exp_mu_dzm, &
        invrs_dzm_on_mu, &
        grav_on_thvm, &
        thl_par_j_precalc, &
        rt_par_j_precalc, &
        tke_i, &
        thl_par_1, &
        rt_par_1, &
        tl_par_1, &
        rsatl_par_1, &
        s_par_1, &
        rc_par_1, &
        thv_par_1, &
        dCAPE_dz_1, &
        CAPE_incr_1, &
        Lv_coef, &
        entrain_coef

    ! Minimum value for Lscale that will taper off with height
    real( kind = core_rknd ) :: lminh

    ! Parcel quantities at grid level j
    real( kind = core_rknd ) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

    ! Used in latent heating calculation
    real( kind = core_rknd ) :: tl_par_j, rsatl_par_j, s_par_j

    ! Variables to make L nonlocal
    real( kind = core_rknd ) :: Lscale_up_max_alt, Lscale_down_min_alt

    ! Variables used to precalculate values
    real( kind = core_rknd ) :: &
        Lv2_coef, &
        tl_par_j_sqd, &
        invrs_dCAPE_diff, &
        invrs_Lscale_sfclyr_depth

    ! ---- Begin Code ----

    !---------- Mixing length computation ----------------------------------

    if( abs(mu) < eps ) then
        write(fstderr,*) "Entrainment rate mu cannot be 0"
        stop "Fatal error in subroutine compute_mixing_length"
    end if


    ! Initialize arrays and precalculate values for computational efficiency
    do i = 1, gr%nz

        ! Initialize up and down arrays
        Lscale_up(i) = zlmin
        Lscale_down(i) = zlmin

        ! Precalculate values to avoid unnecessary calculations later
        exp_mu_dzm(i) = exp( -mu * gr%dzm(i) )
        invrs_dzm_on_mu(i) = ( gr%invrs_dzm(i) ) / mu
        grav_on_thvm(i) = grav / thvm(i)
        Lv_coef(i) = Lv / ( exner(i) * cp ) - ep2 * thv_ds(i)
        entrain_coef(i) = ( 1.0_core_rknd - exp_mu_dzm(i) ) * invrs_dzm_on_mu(i)

    end do

    ! Avoid uninitialized memory (these values are not used in Lscale)
    Lscale_up(1)   = 0.0_core_rknd
    Lscale_down(1) = 0.0_core_rknd

    ! Precalculations of single values to avoid unnecessary calculations later
    Lv2_coef = ep * Lv**2 / ( Rd * cp )
    invrs_Lscale_sfclyr_depth = 1.0_core_rknd / Lscale_sfclyr_depth


    ! Calculate initial turbulent kinetic energy for each grid level
    tke_i = zm2zt( em )
    

    ! ---------------- Upwards Length Scale Calculation ----------------

    ! Precalculate values for upward Lscale, these are useful only if a parcel can rise
    ! more than one level. They are used in the equations that calculate thl and rt 
    ! recursively for a parcel as it ascends
    do j = 2, gr%nz-1

        thl_par_j_precalc(j) = thlm(j) - thlm(j-1) * exp_mu_dzm(j-1)  &
                               - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)

        rt_par_j_precalc(j) = rtm(j) - rtm(j-1) * exp_mu_dzm(j-1)  &
                              - ( rtm(j) - rtm(j-1) ) * entrain_coef(j-1)
    end do


    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain 
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations. 

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    do j = 3, gr%nz

        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j-1) ) * entrain_coef(j-1)

        tl_par_1(j) = thl_par_1(j) * exner(j)

        rt_par_1(j) = rtm(j) - ( rtm(j) - rtm(j-1) ) * entrain_coef(j-1)
       
    end do


    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental
    rsatl_par_1(3:) = compute_rsat_parcel( p_in_Pa(3:), tl_par_1(3:) )


    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    do j = 3, gr%nz

        tl_par_j_sqd = tl_par_1(j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )  
        ! 
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )  
        s_par_1(j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )

        rc_par_1(j) = max( s_par_1(j), zero_threshold )

        ! theta_v of entraining parcel at grid level j 
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        
        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(j) = 0.5_core_rknd * dCAPE_dz_1(j) * gr%dzm(j-1)

    end do

    
    ! Calculate Lscale_up for each grid level. If the TKE from a parcel has not been completely
    ! exhausted by the initial change then continue the exhaustion calculations here for a single
    ! grid level at a time until the TKE is exhausted. 

    Lscale_up_max_alt = 0._core_rknd    ! Set initial max value for Lscale_up to 0
    do i = 2, gr%nz-2

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i) + CAPE_incr_1(i+1) > 0.0_core_rknd ) then

            ! Calculate new TKE for parcel
            tke = tke_i(i) + CAPE_incr_1(i+1)

            ! Set j to 2 levels above current Lscale_up level, this is because we've already
            ! determined that the parcel can rise at least 1 full level
            j = i + 2

            ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j = thl_par_1(i+1)
            rt_par_j  = rt_par_1(i+1)
            dCAPE_dz_j_minus_1 = dCAPE_dz_1(i+1)


            ! Continue change in TKE calculations until it is exhausted or the max grid
            ! level has been reached. j is the next grid level above the level that can
            ! be reached for a parcel starting at level i. If TKE is exhausted in this loop
            ! that means the parcel starting at i cannot reach level j, but has reached j-1
            do while ( j < gr%nz )
           
                ! thl, rt of parcel are conserved except for entrainment
                ! 
                ! The values of thl_env and rt_env are treated as changing linearly for a parcel 
                ! of air ascending from level j-1 to level j

                ! theta_l of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j-1)


                ! r_t of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
                rt_par_j = rt_par_j_precalc(j) + rt_par_j * exp_mu_dzm(j-1)
                

                ! Include effects of latent heating on Lscale_up 6/12/00
                ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
                ! Probably should use properties of bump 1 in Gaussian, not mean!!!

                tl_par_j = thl_par_j*exner(j)

                rsatl_par_j = compute_rsat_parcel( p_in_Pa(j), tl_par_j )

                tl_par_j_sqd = tl_par_j**2

                ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
                !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
                ! and SD's beta (eqn. 8),
                !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
                ! 
                ! Simplified by multiplying top and bottom by tl^2 to avoid a 
                ! divide and precalculating ep * Lv**2 / ( Rd * cp )  
                s_par_j = ( rt_par_j - rsatl_par_j ) * tl_par_j_sqd &
                          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

                rc_par_j = max( s_par_j, zero_threshold )

                ! theta_v of entraining parcel at grid level j 
                thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j  &
                            + Lv_coef(j) * rc_par_j

                ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )

                ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
                ! Trapezoidal estimate between grid levels j and j-1
                CAPE_incr = 0.5_core_rknd * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * gr%dzm(j-1)

                ! Exit loop early if tke has been exhaused between level j and j+1
                if ( tke + CAPE_incr <= 0.0_core_rknd ) then
                    exit
                end if

                ! Save previous dCAPE value for next cycle
                dCAPE_dz_j_minus_1 = dCAPE_dz_j

                ! Caclulate new TKE and increment j
                tke = tke + CAPE_incr
                j = j + 1

            enddo


            ! Add full grid level thickness for each grid level that was passed without the TKE
            ! being exhausted, difference between starting level (i) and last level passed (j-1)
            Lscale_up(i) = Lscale_up(i) + gr%zt(j-1) - gr%zt(i)


            if ( j < gr%nz ) then

                ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
                ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.
        
                if ( abs( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_minus_1 ) * eps ) then

                    ! Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0
                    ! Find the remaining distance z - z_0 that it takes to 
                    ! exhaust the remaining TKE

                    Lscale_up(i) = Lscale_up(i) + ( - tke / dCAPE_dz_j )

                else

                    ! Case used for most scenarios where dCAPE/dz|_(z_1) /= dCAPE/dz|_(z_0)
                    ! Find the remaining distance z - z_0 that it takes to exhaust the
                    ! remaining TKE (tke_i), using the quadratic formula (only the
                    ! negative (-) root works in this scenario).
                    invrs_dCAPE_diff = 1.0_core_rknd / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )

                    Lscale_up(i) = Lscale_up(i) &
                                   - dCAPE_dz_j_minus_1 * invrs_dCAPE_diff * gr%dzm(j-1)  &
                                   - sqrt( dCAPE_dz_j_minus_1**2 &
                                            - 2.0_core_rknd * tke * gr%invrs_dzm(j-1) & 
                                              * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) ) &
                                     * invrs_dCAPE_diff  * gr%dzm(j-1)
                endif

            end if

        else    ! TKE for parcel at level (i) was exhaused before one full grid level
            
            ! Find the remaining distance z - z_0 that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula. Simplified 
            ! since dCAPE_dz_j_minus_1 = 0.0
            Lscale_up(i) = Lscale_up(i) - sqrt( - 2.0_core_rknd * tke_i(i) &
                                                  * gr%dzm(i) * dCAPE_dz_1(i+1) ) &
                                          / dCAPE_dz_1(i+1)  
        endif


        ! If a parcel at a previous grid level can rise past the parcel at this grid level
        ! then this one should also be able to rise up to that height. This feature insures 
        ! that the profile of Lscale_up will be smooth, thus reducing numerical instability.
        if ( gr%zt(i) + Lscale_up(i) < Lscale_up_max_alt ) then

            ! A lower starting parcel can ascend higher than this one, set height to the max
            ! that any lower starting parcel can ascend to
            Lscale_up(i) = Lscale_up_max_alt - gr%zt(i)
        else 
            
            ! This parcel can ascend higher than any below it, save final height
            Lscale_up_max_alt = Lscale_up(i) + gr%zt(i)
        end if


    end do


    ! ---------------- Downwards Length Scale Calculation ----------------

    ! Precalculate values for downward Lscale, these are useful only if a parcel can descend
    ! more than one level. They are used in the equations that calculate thl and rt 
    ! recursively for a parcel as it descends
    do j = 2, gr%nz-1

        thl_par_j_precalc(j) = thlm(j) - thlm(j+1) * exp_mu_dzm(j)  &
                               - ( thlm(j) - thlm(j+1) ) * entrain_coef(j)

        rt_par_j_precalc(j) = rtm(j) - rtm(j+1) * exp_mu_dzm(j)  &
                              - ( rtm(j) - rtm(j+1) ) * entrain_coef(j)
    end do


    ! Calculate the initial change in TKE for each level. This is done for computational
    ! efficiency, it helps because there will be at least one calculation for each grid level,
    ! meaning the first one can be done for every grid level and therefore the calculations can
    ! be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain 
    ! how many more iterations should be done for each individual grid level, and calculating
    ! one change in TKE for each level until all are exhausted will result in many unnessary
    ! and expensive calculations. 

    ! Calculate initial thl, tl, and rt for parcels at each grid level
    do j = 2, gr%nz-1

        thl_par_1(j) = thlm(j) - ( thlm(j) - thlm(j+1) )  * entrain_coef(j)

        tl_par_1(j) = thl_par_1(j) * exner(j)

        rt_par_1(j) = rtm(j) - ( rtm(j) - rtm(j+1) ) * entrain_coef(j)
       
    end do


    ! Caclculate initial rsatl for parcels at each grid level, this function is elemental
    rsatl_par_1(2:) = compute_rsat_parcel( p_in_Pa(2:), tl_par_1(2:) )


    ! Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
    do j = 2, gr%nz-1

        tl_par_j_sqd = tl_par_1(j)**2

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
        !                           s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
        ! and SD's beta (eqn. 8),
        !                           beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )  
        ! 
        ! Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
        ! ep * Lv**2 / ( Rd * cp )  
        s_par_1(j) = ( rt_par_1(j) - rsatl_par_1(j) ) * tl_par_j_sqd &
                     / ( tl_par_j_sqd + Lv2_coef * rsatl_par_1(j) )

        rc_par_1(j) = max( s_par_1(j), zero_threshold )

        ! theta_v of entraining parcel at grid level j 
        thv_par_1(j) = thl_par_1(j) + ep1 * thv_ds(j) * rt_par_1(j) + Lv_coef(j) * rc_par_1(j)

        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        dCAPE_dz_1(j) = grav_on_thvm(j) * ( thv_par_1(j) - thvm(j) )

        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
        ! Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
        CAPE_incr_1(j) = 0.5_core_rknd * dCAPE_dz_1(j) * gr%dzm(j)

    end do


    ! Calculate Lscale_down for each grid level. If the TKE from a parcel has not been completely
    ! exhausted by the initial change then continue the exhaustion calculations here for a single
    ! grid level at a time until the TKE is exhausted. 

    Lscale_down_min_alt = gr%zt(gr%nz)  ! Set initial min value for Lscale_down to max zt
    do i = gr%nz, 3, -1

        ! If the initial turbulent kinetic energy (tke) has not been exhausted for this grid level
        if ( tke_i(i) - CAPE_incr_1(i-1) > 0.0_core_rknd ) then

            ! Calculate new TKE for parcel
            tke = tke_i(i) - CAPE_incr_1(i-1)

            ! Set j to 2 levels below current Lscale_down level, this is because we've already
            ! determined that the parcel can descend at least 1 full level
            j = i - 2

            ! Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j = thl_par_1(i-1)
            rt_par_j = rt_par_1(i-1)
            dCAPE_dz_j_plus_1 = dCAPE_dz_1(i-1)


            ! Continue change in TKE calculations until it is exhausted or the min grid
            ! level has been reached. j is the next grid level below the level that can
            ! be reached for a parcel starting at level i. If TKE is exhausted in this loop
            ! that means the parcel starting at i cannot sink to level j, but can sink to j+1
            do while ( j >= 2 )

                ! thl, rt of parcel are conserved except for entrainment
                ! 
                ! The values of thl_env and rt_env are treated as changing linearly for a parcel 
                ! of air descending from level j to level j-1

                ! theta_l of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(thl_par)/dz = - mu * ( thl_par - thl_env )
                thl_par_j = thl_par_j_precalc(j) + thl_par_j * exp_mu_dzm(j)


                ! r_t of the parcel starting at grid level i, and currenly
                ! at grid level j
                ! 
                ! d(rt_par)/dz = - mu * ( rt_par - rt_env )
                rt_par_j = rt_par_j_precalc(j) + rt_par_j * exp_mu_dzm(j)


                ! Include effects of latent heating on Lscale_up 6/12/00
                ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
                ! Probably should use properties of bump 1 in Gaussian, not mean!!!

                tl_par_j = thl_par_j*exner(j)

                rsatl_par_j = compute_rsat_parcel( p_in_Pa(j), tl_par_j )

                tl_par_j_sqd = tl_par_j**2

                ! s from Lewellen and Yoh 1993 (LY) eqn. 1 
                !                         s = ( rt - rsatl ) / ( 1 + beta * rsatl )  
                ! and SD's beta (eqn. 8),
                !                         beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
                ! 
                ! Simplified by multiplying top and bottom by tl^2 to avoid a 
                ! divide and precalculating ep * Lv**2 / ( Rd * cp )  
                s_par_j = (rt_par_j - rsatl_par_j) * tl_par_j_sqd &
                          / ( tl_par_j_sqd + Lv2_coef * rsatl_par_j )

                rc_par_j = max( s_par_j, zero_threshold )

                ! theta_v of entraining parcel at grid level j 
                thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j + Lv_coef(j) * rc_par_j

                ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
                dCAPE_dz_j = grav_on_thvm(j) * ( thv_par_j - thvm(j) )

                ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
                ! Trapezoidal estimate between grid levels j+1 and j
                CAPE_incr = 0.5_core_rknd * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * gr%dzm(j)

                ! Exit loop early if tke has been exhaused between level j+1 and j
                if ( tke - CAPE_incr <= 0.0_core_rknd ) then
                    exit
                endif

                ! Save previous dCAPE value for next cycle
                dCAPE_dz_j_plus_1 = dCAPE_dz_j

                ! Caclulate new TKE and increment j
                tke = tke - CAPE_incr
                j = j - 1

            enddo

            ! Add full grid level thickness for each grid level that was passed without the TKE
            ! being exhausted, difference between starting level (i) and last level passed (j+1)
            Lscale_down(i) = Lscale_down(i) + gr%zt(i) - gr%zt(j+1)

            
            if ( j >= 2 ) then
        
                ! Loop terminated early, meaning TKE was completely exhaused at grid level j.
                ! Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.

                if ( abs( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) * 2 <= &
                     abs( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) * eps ) then

                    ! Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0
                    ! Find the remaining distance z_0 - z that it takes to 
                    ! exhaust the remaining TKE

                    Lscale_down(i) = Lscale_down(i) + ( tke / dCAPE_dz_j )

                else

                    ! Case used for most scenarios where dCAPE/dz|_(z_(-1)) /= dCAPE/dz|_(z_0)
                    ! Find the remaining distance z_0 - z that it takes to exhaust the
                    ! remaining TKE (tke_i), using the quadratic formula (only the
                    ! negative (-) root works in this scenario) -- however, the
                    ! negative (-) root is divided by another negative (-) factor,
                    ! which results in an overall plus (+) sign in front of the
                    ! square root term in the equation below).
                    invrs_dCAPE_diff = 1.0_core_rknd / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )

                    Lscale_down(i) = Lscale_down(i) &
                                     - dCAPE_dz_j_plus_1 * invrs_dCAPE_diff * gr%dzm(j)  &
                                     + sqrt( dCAPE_dz_j_plus_1**2 &
                                             + 2.0_core_rknd * tke * gr%invrs_dzm(j)  &
                                               * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                                       * invrs_dCAPE_diff * gr%dzm(j)
                endif

            end if

        else    ! TKE for parcel at level (i) was exhaused before one full grid level

            ! Find the remaining distance z_0 - z that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula. Simplified 
            ! since dCAPE_dz_j_plus_1 = 0.0
            Lscale_down(i) = Lscale_down(i) + sqrt( 2.0_core_rknd * tke_i(i) &
                                                    * gr%dzm(i-1) * dCAPE_dz_1(i-1) ) &
                                              / dCAPE_dz_1(i-1) 
        endif

        ! If a parcel at a previous grid level can descend past the parcel at this grid level
        ! then this one should also be able to descend down to that height. This feature insures
        ! that the profile of Lscale_down will be smooth, thus reducing numerical instability.
        if ( gr%zt(i) - Lscale_down(i) > Lscale_down_min_alt ) then
            Lscale_down(i) = gr%zt(i) - Lscale_down_min_alt
        else 
            Lscale_down_min_alt = gr%zt(i) - Lscale_down(i)
        end if

    end do


    ! ---------------- Final Lscale Calculation ----------------

    do i = 2, gr%nz, 1

        ! Make lminh a linear function starting at value lmin at the bottom
        ! and going to zero at 500 meters in altitude.
        if( l_implemented ) then

            ! Within a host model, increase mixing length in 500 m layer above *ground*
            lminh = max( zero_threshold, Lscale_sfclyr_depth - ( gr%zt(i) - gr%zm(1) ) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        else

            ! In standalone mode, increase mixing length in 500 m layer above *mean sea level*
            lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i) ) &
                    * lmin * invrs_Lscale_sfclyr_depth
        end if

        Lscale_up(i)    = max( lminh, Lscale_up(i) )
        Lscale_down(i)  = max( lminh, Lscale_down(i) )

        ! When L is large, turbulence is weakly damped
        ! When L is small, turbulence is strongly damped
        ! Use a geometric mean to determine final Lscale so that L tends to become small 
        ! if either Lscale_up or Lscale_down becomes small. 
        Lscale(i) = sqrt( Lscale_up(i) * Lscale_down(i) )

    enddo

    ! Set the value of Lscale at the upper and lower boundaries.
    Lscale(1) = Lscale(2)
    Lscale(gr%nz) = Lscale(gr%nz-1)

    ! Vince Larson limited Lscale to allow host
    ! model to take over deep convection.  13 Feb 2008.
    Lscale = min( Lscale, Lscale_max )

 
    ! Ensure that no Lscale values are NaN
    if ( clubb_at_least_debug_level( 1 ) ) then

        call length_check( Lscale, Lscale_up, Lscale_down )

        if ( err_code == clubb_fatal_error ) then

          write(fstderr,*) "Errors in compute_mixing_length subroutine"

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "thvm = ", thvm
          write(fstderr,*) "thlm = ", thlm
          write(fstderr,*) "rtm = ", rtm
          write(fstderr,*) "em = ", em
          write(fstderr,*) "exner = ", exner
          write(fstderr,*) "p_in_Pa = ", p_in_Pa
          write(fstderr,*) "thv_ds = ", thv_ds

          write(fstderr,*) "Intent(out)"

          write(fstderr,*) "Lscale = ", Lscale
          write(fstderr,*) "Lscale_up = ", Lscale_up
          write(fstderr,*) "Lscale_down = ", Lscale_down
     
        endif ! Fatal error

    end if

    return

  end subroutine compute_mixing_length

!===============================================================================

  elemental function compute_rsat_parcel( p_in_Pa, tl_par ) result( rsatl_par )

  ! Description:
  !   Computes rsat for a parcel

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Precision

    use model_flags, only: &
      l_sat_mixrat_lookup ! Variable(s)

    use saturation, only:  & 
      sat_mixrat_liq, & ! Procedure(s)
      sat_mixrat_liq_lookup, &
      sat_mixrat_ice

    use constants_clubb, only: &
      zero, &
      one, &
      T_freeze_K, &
      fstderr

    use error_code, only: &
      clubb_at_least_debug_level

    implicit none

    ! Parameters
    logical, parameter :: &
      l_include_ice = .false. ! Include ice in calculation of rsat_par

    real( kind = core_rknd ), parameter :: &
      T_all_ice = 233.15_core_rknd ! Temperature at which only ice is included in calculation [K]

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      p_in_Pa, &    ! Pressure at this grid level           [Pa]
      tl_par        ! tl at this grid level                 [K]

    ! Result variable
    real( kind = core_rknd ) :: &
      rsatl_par

    real( kind = core_rknd ) :: &
      sat_ice_ratio, &    ! Ratio of interpolation between sat_mixrat_liq and sat_mixrat_ice.
                          ! sat_ice_ratio=1 ---> rsatl_par = sat_mixrat_ice
      sat_mixrat_liq_res, &
      sat_mixrat_ice_res

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    if ( l_include_ice ) then
      if ( tl_par >= T_freeze_K ) then
        sat_ice_ratio = zero
      else if ( tl_par <= T_all_ice ) then
        sat_ice_ratio = one
      else
        ! Linear interpolation
        sat_ice_ratio = ( T_freeze_K - tl_par ) / ( T_freeze_K - T_all_ice )
      end if
    else
      sat_ice_ratio = zero
    end if ! l_include_ice

    if ( sat_ice_ratio < one ) then
      ! Include liquid.
      if ( l_sat_mixrat_lookup ) then
        sat_mixrat_liq_res = sat_mixrat_liq_lookup( p_in_Pa, tl_par )
      else
        sat_mixrat_liq_res = sat_mixrat_liq( p_in_Pa, tl_par )
      end if
    else
      sat_mixrat_liq_res = zero
    end if

    if ( sat_ice_ratio > zero ) then
      ! Include ice.
      sat_mixrat_ice_res = sat_mixrat_ice( p_in_Pa, tl_par )
    else
      sat_mixrat_ice_res = zero
    end if

    rsatl_par = sat_mixrat_liq_res + ( sat_mixrat_ice_res - sat_mixrat_liq_res ) * sat_ice_ratio

    return
  end function compute_rsat_parcel

!===============================================================================

end module mixing_length
