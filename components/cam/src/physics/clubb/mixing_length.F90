!-----------------------------------------------------------------------
! $Id: mixing_length.F90 7226 2014-08-19 15:52:41Z betlej@uwm.edu $
!===============================================================================
module mixing_length

  implicit none

  private ! Default Scope

  public :: compute_length

  contains

  !=============================================================================
  subroutine compute_length( thvm, thlm, rtm, em, Lscale_max, &
                             p_in_Pa, exner, thv_ds, mu, l_implemented, &
                             err_code, &
                             Lscale, Lscale_up, Lscale_down )
    ! Description:
    ! Larson's 5th moist, nonlocal length scale

    ! References:
    ! Section 3b ( /Eddy length formulation/ ) of
    ! ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    ! Method and Model Description'' Golaz, et al. (2002)
    ! JAS, Vol. 59, pp. 3540--3551.

    !-----------------------------------------------------------------------

    ! mu = (1/M) dM/dz > 0.  mu=0 for no entrainment.
    ! Siebesma recommends mu=2e-3, although most schemes use mu=1e-4
    ! When mu was fixed, we used the value mu = 6.e-4

    use constants_clubb, only:  &  ! Variable(s)
        Cp,            & ! Dry air specific heat at constant pressure [J/kg/K]
        Rd,            & ! Dry air gas constant                       [J/kg/K]
        ep,            & ! Rd / Rv                                    [-]
        ep1,           & ! (1-ep)/ep                                  [-]
        ep2,           & ! 1/ep                                       [-]
        Lv,            & ! Latent heat of vaporiztion                 [J/kg/K]
        grav,          & ! Gravitational acceleration                 [m/s^2]
        fstderr,       &
        zero_threshold

    use parameters_tunable, only:  &  ! Variable(s)
        lmin    ! Minimum value for Lscale                         [m]

    use grid_class, only:  & 
        gr,  & ! Variable(s)
        zm2zt ! Procedure(s)

    use numerical_check, only:  & 
        length_check ! Procedure(s)

    use saturation, only:  & 
        sat_mixrat_liq, & ! Procedure(s)
        sat_mixrat_liq_lookup

    use error_code, only:  & 
        clubb_at_least_debug_level, & ! Procedure(s)
        fatal_error

    use error_code, only:  & 
        clubb_no_error ! Constant

    use model_flags, only: &
        l_sat_mixrat_lookup ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

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

    ! Output Variables
    integer, intent(inout) :: & 
      err_code

    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  & 
      Lscale,    & ! Mixing length      [m]
      Lscale_up, & ! Mixing length up   [m]
      Lscale_down  ! Mixing length down [m]

    ! Local Variables

    integer :: i, j, &
      err_code_Lscale

    real( kind = core_rknd ) :: tke_i, CAPE_incr

    real( kind = core_rknd ) :: dCAPE_dz_j, dCAPE_dz_j_minus_1, dCAPE_dz_j_plus_1

    ! Temporary arrays to store calculations to speed runtime
    real( kind = core_rknd ), dimension(gr%nz) :: exp_mu_dzm, invrs_dzm_on_mu

    ! Minimum value for Lscale that will taper off with height
    real( kind = core_rknd ) :: lminh

    ! Parcel quantities at grid level j
    real( kind = core_rknd ) :: thl_par_j, rt_par_j, rc_par_j, thv_par_j

    ! Used in latent heating calculation
    real( kind = core_rknd ) :: tl_par_j, rsatl_par_j, beta_par_j, & 
            s_par_j

    ! Parcel quantities at grid level j-1
    real( kind = core_rknd ) :: thl_par_j_minus_1, rt_par_j_minus_1

    ! Parcel quantities at grid level j+1
    real( kind = core_rknd ) :: thl_par_j_plus_1, rt_par_j_plus_1

    ! Variables to make L nonlocal
    real( kind = core_rknd ) :: Lscale_up_max_alt, Lscale_down_min_alt

    ! ---- Begin Code ----

    err_code_Lscale = clubb_no_error

    !---------- Mixing length computation ----------------------------------

    ! Avoid uninitialized memory (these values are not used in Lscale)
    ! -dschanen 12 March 2008
    Lscale_up(1)   = 0.0_core_rknd
    Lscale_down(1) = 0.0_core_rknd

    ! Initialize exp_mu_dzm--sets each exp_mu_dzm value to its corresponding
    !   exp(-mu/gr%invrs_dzm) value. In theory, this saves 11 computations of
    !   exp(-mu/gr%invrs_dzm) used below.
    ! ~~EIHoppe//20090615
    exp_mu_dzm(:)  = exp( -mu/gr%invrs_dzm(:) )

    ! Initialize invrs_dzm_on_mu -- sets each invrs_dzm_on_mu value to its
    ! corresponding (gr%invrs_dzm/mu) value. This will save computations of
    ! this value below.
    ! ~EIHoppe//20100728
    invrs_dzm_on_mu(:) = (gr%invrs_dzm(:))/mu

    !!!!! Compute Lscale_up for every vertical level.

    ! Upwards loop

    Lscale_up_max_alt = 0._core_rknd
    do i = 2, gr%nz, 1

      tke_i = zm2zt( em, i )   ! TKE interpolated to thermodynamic level

      Lscale_up(i) = zlmin
      j = i + 1

      thl_par_j_minus_1 = thlm(i)
      rt_par_j_minus_1  = rtm(i)
      dCAPE_dz_j_minus_1 = 0.0_core_rknd

      do while ((tke_i > 0._core_rknd) .and. (j < gr%nz))

        ! thl, rt of parcel are conserved except for entrainment

        ! theta_l of the parcel at grid level j.
        !
        ! The equation for the rate of change of theta_l of the parcel with
        ! respect to height, due to entrainment, is:
        !
        ! d(thl_par)/dz = - mu * ( thl_par - thl_env );
        !
        ! where thl_par is theta_l of the parcel, thl_env is theta_l of the
        ! ambient (or environmental) air, and mu is the entrainment rate,
        ! such that:
        !
        ! mu = (1/m)*(dm/dz);
        !
        ! where m is the mass of the parcel.  The value of mu is set to be a
        ! constant.
        !
        ! The differential equation is solved for thl_par_j (thl_par at
        ! height gr%zt(j)) given the boundary condition thl_par_j_minus_1
        ! (thl_par at height gr%zt(j-1)), and given the fact that the value
        ! of thl_env is treated as changing linearly for a parcel of air
        ! ascending from level j-1 (where thl_env has the value thlm(j-1)) to
        ! level j (where thl_env has the value thlm(j)).
        !
        ! For the special case where entrainment rate, mu, is set to 0,
        ! thl_par remains constant as the parcel ascends.

        if ( mu /= 0.0_core_rknd ) then

          ! The ascending parcel is entraining at rate mu.

          ! Calculation changed to use pre-calculated exp(-mu/gr%invrs_dzm)
          ! values. ~~EIHoppe//20090615

          ! Calculation changed to use pre-calculated mu/gr%invrs_dzm values.
          ! ~EIHoppe//20100728

          thl_par_j = thlm(j) - thlm(j-1)*exp_mu_dzm(j-1)  &
                      - ( 1.0_core_rknd - exp_mu_dzm(j-1))  &
                        * ( (thlm(j) - thlm(j-1))  &
                        * invrs_dzm_on_mu(j-1) ) &
!                               / (mu/gr%invrs_dzm(j-1)) )  &
                      + thl_par_j_minus_1 * exp_mu_dzm(j-1)

        else

          ! The ascending parcel is not entraining.

          thl_par_j = thl_par_j_minus_1

        endif

        ! r_t of the parcel at grid level j.
        !
        ! The equation for the rate of change of r_t of the parcel with
        ! respect to height, due to entrainment, is:
        !
        ! d(rt_par)/dz = - mu * ( rt_par - rt_env );
        !
        ! where rt_par is r_t of the parcel, rt_env is r_t of the ambient (or
        ! environmental) air, and mu is the entrainment rate, such that:
        !
        ! mu = (1/m)*(dm/dz);
        !
        ! where m is the mass of the parcel.  The value of mu is set to be a
        ! constant.
        !
        ! The differential equation is solved for rt_par_j (rt_par at height
        ! gr%zt(j)) given the boundary condition rt_par_j_minus_1 (rt_par at
        ! height gr%zt(j-1)), and given the fact that the value of rt_env is
        ! treated as changing linearly for a parcel of air ascending from
        ! level j-1 (where rt_env has the value rtm(j-1)) to level j (where
        ! rt_env has the value rtm(j)).
        !
        ! For the special case where entrainment rate, mu, is set to 0,
        ! rt_par remains constant as the parcel ascends.

        if ( mu /= 0.0_core_rknd ) then

          ! The ascending parcel is entraining at rate mu.

          ! Calculation changed to use pre-calculated exp(-mu/gr%invrs_dzm)
          ! values. ~~EIHoppe//20090615

          ! Calculation changed to use pre-calculated mu/gr%invrs_dzm values.
          ! ~EIHoppe//20100728

          rt_par_j = rtm(j) - rtm(j-1)*exp_mu_dzm(j-1)  &
                     - ( 1.0_core_rknd - exp_mu_dzm(j-1))  &
                       * ( (rtm(j) - rtm(j-1)) &
                        * invrs_dzm_on_mu(j-1) ) &
!                          / (mu/gr%invrs_dzm(j-1)) )  &
                     + rt_par_j_minus_1 * exp_mu_dzm(j-1)

        else

          ! The ascending parcel is not entraining.

          rt_par_j = rt_par_j_minus_1

        endif

        ! Include effects of latent heating on Lscale_up 6/12/00
        ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
        ! Probably should use properties of bump 1 in Gaussian, not mean!!!

        ! Calculate r_c of the parcel at grid level j based on the values of
        ! theta_l of the parcel and r_t of the parcel at grid level j.
        tl_par_j = thl_par_j*exner(j)
        if ( l_sat_mixrat_lookup ) then
          rsatl_par_j = sat_mixrat_liq_lookup( p_in_Pa(j), tl_par_j )
        else
          rsatl_par_j = sat_mixrat_liq( p_in_Pa(j), tl_par_j )
        end if
        ! SD's beta (eqn. 8)
        beta_par_j = ep*(Lv/(Rd*tl_par_j))*(Lv/(cp*tl_par_j))
        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        s_par_j = (rt_par_j-rsatl_par_j)/(1._core_rknd+beta_par_j*rsatl_par_j)
        rc_par_j = max( s_par_j, zero_threshold )

        ! theta_v of entraining parcel at grid level j.
        thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j  &
                    + ( Lv / (exner(j)*cp) - ep2 * thv_ds(j) ) * rc_par_j

        ! Lscale_up and CAPE increment.
        !
        ! The equation for Lscale_up is:
        !
        ! INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
        !
        ! where thv_par is theta_v of the parcel, thvm is the mean
        ! environmental value of theta_v, z_i is the altitude that the parcel
        ! started its ascent from, and em is the mean value of TKE at
        ! altitude z_i (which gives the parcel its initial upward boost).
        !
        ! The increment of CAPE for any two successive vertical levels (z_0
        ! and z_1, such that z_0 < z_1, and where z_0 is gr%zt(j-1) and z_1
        ! is gr%zt(j)) is:
        !
        ! CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz.
        !
        ! Thus, the derivative of CAPE with respect to height is:
        !
        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        !
        ! A purely trapezoidal rule is used between levels z_0 and z_1, such
        ! that dCAPE/dz is evaluated at levels z_0 and z_1, and is considered
        ! to vary linearly at all altitudes z_0 <= z <= z_1.  Thus, dCAPE/dz
        ! is considered to be of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
        ! where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 ).
        !
        ! The integral is evaluated to find the CAPE increment between two
        ! successive vertical levels.  The result either adds to or depletes
        ! from the total amount of energy that keeps the parcel ascending.

        dCAPE_dz_j = ( grav/thvm(j) ) * ( thv_par_j - thvm(j) )

        CAPE_incr = 0.5_core_rknd * ( dCAPE_dz_j + dCAPE_dz_j_minus_1 )  &
                          / gr%invrs_dzm(j-1)

        if ( tke_i + CAPE_incr > 0.0_core_rknd ) then

          ! The total amount of CAPE increment has not exhausted the initial
          ! TKE (plus any additions by CAPE increments due to upward
          ! buoyancy) that boosted and carried the parcel upward.  The
          ! thickness of the full grid level is added to Lscale_up.

          Lscale_up(i) = Lscale_up(i) + gr%zt(j) - gr%zt(j-1)

        else

          ! The total amount of CAPE increment has exhausted the initial TKE
          ! (plus any additions by CAPE increments due to upward buoyancy)
          ! that boosted and carried the parcel upward.  Add the thickness
          ! z - z_0 (where z_0 < z <= z_1) to Lscale_up.  The calculation of
          ! Lscale_up is complete.

          if ( dCAPE_dz_j == dCAPE_dz_j_minus_1 ) then

            ! Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0,
            ! thus making factor A (above) equal to 0.  Find the remaining
            ! distance z - z_0 that it takes to exhaust the remaining TKE
            ! (tke_i).

            Lscale_up(i)  &
            = Lscale_up(i)  &
            + ( - tke_i / dCAPE_dz_j )

          else

            ! Case used for most scenarios where dCAPE/dz|_(z_1)
            ! /= dCAPE/dz|_(z_0), thus making factor A /= 0.  Find the
            ! remaining distance z - z_0 that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula (only the
            ! negative (-) root works in this scenario).

            Lscale_up(i)  &
            = Lscale_up(i)  &
            + ( - dCAPE_dz_j_minus_1 /  &
                 ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) )  &
                 / gr%invrs_dzm(j-1)  &
            - sqrt( dCAPE_dz_j_minus_1**2  &
                    - 2.0_core_rknd * tke_i * gr%invrs_dzm(j-1)  &
                      * ( dCAPE_dz_j - dCAPE_dz_j_minus_1 ) )  &
                 / ( dCAPE_dz_j - dCAPE_dz_j_minus_1 )  &
                 / gr%invrs_dzm(j-1)

          endif

        endif

        ! Reset values for use during the next vertical level up.

        thl_par_j_minus_1 = thl_par_j
        rt_par_j_minus_1 = rt_par_j
        dCAPE_dz_j_minus_1 = dCAPE_dz_j

        tke_i = tke_i + CAPE_incr
        j = j + 1

      enddo

      ! Make Lscale_up nonlocal
      !
      ! This code makes the value of Lscale_up nonlocal.  Thus, if a parcel
      ! starting from a lower altitude can ascend to altitude
      ! Lscale_up_max_alt, then a parcel starting from a higher altitude should
      ! also be able to ascend to at least altitude Lscale_up_max_alt, even if
      ! the local result of Lscale_up for the parcel that started at a higher
      ! altitude is not sufficient for the parcel to reach altitude
      ! Lscale_up_max_alt.
      !
      ! For example, if it was found that a parcel starting at an altitude of
      ! 100 m. ascended to an altitude of 2100 m. (an Lscale_up value of
      ! 2000 m.), then a parcel starting at an altitude of 200 m. should also
      ! be able to ascend to an altitude of at least 2100 m.  If Lscale_up
      ! was found to be only 1800 m. for the parcel starting at 200 m.
      ! (resulting in the parcel only being able to ascend to an altitude of
      ! 2000 m.), then this code will overwrite the 1800 m. value with a
      ! Lscale_up value of 1900 m. (so that the parcel reaches an altitude of
      ! 2100 m.).
      !
      ! This feature insures that the profile of Lscale_up will be very smooth,
      ! thus reducing numerical instability in the model.

      Lscale_up_max_alt = max( Lscale_up_max_alt, Lscale_up(i)+gr%zt(i) )

      if ( ( gr%zt(i) + Lscale_up(i) ) < Lscale_up_max_alt ) then
        Lscale_up(i) = Lscale_up_max_alt - gr%zt(i)
      endif

    enddo


    !!!!! Compute Lscale_down for every vertical level.

    ! Do it again for downwards particle motion.
    ! For now, do not include latent heat

    ! Chris Golaz modification to include effects on latent heating
    ! on Lscale_down

    Lscale_down_min_alt = gr%zt(gr%nz)
    do i = gr%nz, 2, -1

      tke_i = zm2zt( em, i )   ! TKE interpolated to thermodynamic level

      Lscale_down(i) = zlmin
      j = i - 1

      thl_par_j_plus_1 = thlm(i)
      rt_par_j_plus_1 = rtm(i)
      dCAPE_dz_j_plus_1 = 0.0_core_rknd

      do while ( (tke_i > 0._core_rknd) .and. (j >= 2) )

        ! thl, rt of parcel are conserved except for entrainment

        ! theta_l of the parcel at grid level j.
        !
        ! The equation for the rate of change of theta_l of the parcel with
        ! respect to height, due to entrainment, is:
        !
        ! d(thl_par)/dz = - mu * ( thl_par - thl_env );
        !
        ! where thl_par is theta_l of the parcel, thl_env is theta_l of the
        ! ambient (or environmental) air, and mu is the entrainment rate,
        ! such that:
        !
        ! mu = (1/m)*(dm/dz);
        !
        ! where m is the mass of the parcel.  The value of mu is set to be a
        ! constant.
        !
        ! NOTE:  For an entraining, descending parcel, parcel mass will
        !        increase as height decreases.  Thus dm/dz < 0, and therefore
        !        mu < 0.  However, in the equation for thl_par_j, mu is always
        !        multiplied by the delta_z factor ( gr%zt(j) - gr%zt(j+1) ),
        !        which always has the propery delta_z < 0 for a descending
        !        parcel.  Thus, mu*delta_z > 0, just as for an entraining,
        !        ascending parcel.  Therefore, the same general form of the
        !        entrainment equation (only with differing grid level indices)
        !        can be used for both the ascending and descending parcels.
        !
        ! The differential equation is solved for thl_par_j (thl_par at
        ! height gr%zt(j)) given the boundary condition thl_par_j_plus_1
        ! (thl_par at height gr%zt(j+1)), and given the fact that the value
        ! of thl_env is treated as changing linearly for a parcel of air
        ! descending from level j+1 (where thl_env has the value thlm(j+1)) to
        ! level j (where thl_env has the value thlm(j)).
        !
        ! For the special case where entrainment rate, mu, is set to 0,
        ! thl_par remains constant as the parcel descends.

        if ( mu /= 0.0_core_rknd ) then

          ! The descending parcel is entraining at rate mu.

          ! Calculation changed to use pre-calculated exp(-mu/gr%invrs_dzm)
          ! values. ~~EIHoppe//20090615

          ! Calculation changed to use pre-calculated mu/gr%invrs_dzm values.
          ! ~EIHoppe//20100728

          thl_par_j = thlm(j) - thlm(j+1)*exp_mu_dzm(j)  &
                      - ( 1.0_core_rknd - exp_mu_dzm(j))  &
                        * ( (thlm(j) - thlm(j+1)) &
                        * invrs_dzm_on_mu(j) ) &
!                          / (mu/gr%invrs_dzm(j)) )  &
                      + thl_par_j_plus_1 * exp_mu_dzm(j)

        else

          ! The descending parcel is not entraining.

          thl_par_j = thl_par_j_plus_1

        endif

        ! r_t of the parcel at grid level j.
        !
        ! The equation for the rate of change of r_t of the parcel with
        ! respect to height, due to entrainment, is:
        !
        ! d(rt_par)/dz = - mu * ( rt_par - rt_env );
        !
        ! where rt_par is r_t of the parcel, rt_env is r_t of the ambient (or
        ! environmental) air, and mu is the entrainment rate, such that:
        !
        ! mu = (1/m)*(dm/dz);
        !
        ! where m is the mass of the parcel.  The value of mu is set to be a
        ! constant.
        !
        ! NOTE:  For an entraining, descending parcel, parcel mass will
        !        increase as height decreases.  Thus dm/dz < 0, and therefore
        !        mu < 0.  However, in the equation for rt_par_j, mu is always
        !        multiplied by the delta_z factor ( gr%zt(j) - gr%zt(j+1) ),
        !        which always has the propery delta_z < 0 for a descending
        !        parcel.  Thus, mu*delta_z > 0, just as for an entraining,
        !        ascending parcel.  Therefore, the same general form of the
        !        entrainment equation (only with differing grid level indices)
        !        can be used for both the ascending and descending parcels.
        !
        ! The differential equation is solved for rt_par_j (rt_par at height
        ! gr%zt(j)) given the boundary condition rt_par_j_plus_1 (rt_par at
        ! height gr%zt(j+1)), and given the fact that the value of rt_env is
        ! treated as changing linearly for a parcel of air descending from
        ! level j+1 (where rt_env has the value rtm(j+1)) to level j (where
        ! rt_env has the value rtm(j)).
        !
        ! For the special case where entrainment rate, mu, is set to 0,
        ! rt_par remains constant as the parcel descends.

        if ( mu /= 0.0_core_rknd ) then

          ! The descending parcel is entraining at rate mu.

          ! Calculation changed to use pre-calculated exp(-mu/gr%invrs_dzm)
          ! values. ~~EIHoppe//20090615

          ! Calculation changed to use pre-calculated mu/gr%invrs_dzm values.
          ! ~EIHoppe//20100728

          rt_par_j = rtm(j) - rtm(j+1)*exp_mu_dzm(j)  &
                     - ( 1.0_core_rknd - exp_mu_dzm(j) )  &
                       * ( (rtm(j) - rtm(j+1)) &
!                         / (mu/gr%invrs_dzm(j)) )  &
                        * invrs_dzm_on_mu(j) ) &
                     + rt_par_j_plus_1 * exp_mu_dzm(j)

        else

          ! The descending parcel is not entraining.

          rt_par_j = rt_par_j_plus_1

        endif

        ! Include effects of latent heating on Lscale_down
        ! Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
        ! Probably should use properties of bump 1 in Gaussian, not mean!!!

        ! Calculate r_c of the parcel at grid level j based on the values of
        ! theta_l of the parcel and r_t of the parcel at grid level j.
        tl_par_j = thl_par_j*exner(j)
        if ( l_sat_mixrat_lookup ) then
          rsatl_par_j = sat_mixrat_liq_lookup( p_in_Pa(j), tl_par_j )
        else
          rsatl_par_j = sat_mixrat_liq( p_in_Pa(j), tl_par_j )
        end if
        ! SD's beta (eqn. 8)
        beta_par_j = ep*(Lv/(Rd*tl_par_j))*(Lv/(cp*tl_par_j))
        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        s_par_j = (rt_par_j-rsatl_par_j)/(1._core_rknd+beta_par_j*rsatl_par_j)
        rc_par_j = max( s_par_j, zero_threshold )

        ! theta_v of the entraining parcel at grid level j.
        thv_par_j = thl_par_j + ep1 * thv_ds(j) * rt_par_j &
                    + ( Lv / (exner(j)*cp) - ep2 * thv_ds(j) ) * rc_par_j

        ! Lscale_down and CAPE increment.
        !
        ! The equation for Lscale_down (where Lscale_down is the absolute
        ! value of downward distance) is:
        !
        ! INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
        !
        ! where thv_par is theta_v of the parcel, thvm is the mean
        ! environmental value of theta_v, z_i is the altitude that the parcel
        ! started its descent from, and em is the mean value of TKE at
        ! altitude z_i (which gives the parcel its initial downward boost).
        !
        ! The increment of CAPE for any two successive vertical levels (z_0
        ! and z_(-1), such that z_(-1) < z_0, and where z_0 is gr%zt(j+1) and
        ! z_(-1) is gr%zt(j)) is:
        !
        ! CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz.
        !
        ! Thus, the derivative of CAPE with respect to height is:
        !
        ! dCAPE/dz = g * ( thv_par - thvm ) / thvm.
        !
        ! A purely trapezoidal rule is used between levels z_(-1) and z_0,
        ! such that dCAPE/dz is evaluated at levels z_(-1) and z_0, and is
        ! considered to vary linearly at all altitudes z_(-1) <= z <= z_0.
        ! Thus, dCAPE/dz is considered to be of the form:
        ! A * (z-zo) + dCAPE/dz|_(z_0), where
        ! A = ( dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) ) / ( z_(-1) - z_0 ).
        !
        ! The integral is evaluated to find the CAPE increment between two
        ! successive vertical levels.  The result either adds to or depletes
        ! from the total amount of energy that keeps the parcel descending.

        dCAPE_dz_j = ( grav/thvm(j) ) * ( thv_par_j - thvm(j) )

        CAPE_incr = 0.5_core_rknd * ( dCAPE_dz_j + dCAPE_dz_j_plus_1 ) / gr%invrs_dzm(j)

        if ( tke_i - CAPE_incr > 0.0_core_rknd ) then

          ! The total amount of CAPE increment has not exhausted the initial
          ! TKE (plus any additions by CAPE increments due to downward
          ! buoyancy) that boosted and carried the parcel downward.  The
          ! thickness of the full grid level is added to Lscale_down.

          Lscale_down(i) = Lscale_down(i) + gr%zt(j+1) - gr%zt(j)

        else

          ! The total amount of CAPE increment has exhausted the initial TKE
          ! (plus any additions by CAPE increments due to downward buoyancy)
          ! that boosted and carried the parcel downward.  Add the thickness
          ! z_0 - z (where z_(-1) <= z < z_0) to Lscale_down.  The
          ! calculation of Lscale_down is complete.

          if ( dCAPE_dz_j == dCAPE_dz_j_plus_1 ) then

            ! Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0,
            ! thus making factor A (above) equal to 0.  Find the remaining
            ! distance z_0 - z that it takes to exhaust the remaining TKE
            ! (tke_i).

            Lscale_down(i)  &
            = Lscale_down(i)  &
            + ( tke_i / dCAPE_dz_j )

          else

            ! Case used for most scenarios where dCAPE/dz|_(z_(-1))
            ! /= dCAPE/dz|_(z_0), thus making factor A /= 0.  Find the
            ! remaining distance z_0 - z that it takes to exhaust the
            ! remaining TKE (tke_i), using the quadratic formula (only the
            ! negative (-) root works in this scenario -- however, the
            ! negative (-) root is divided by another negative (-) factor,
            ! which results in an overall plus (+) sign in front of the
            ! square root term in the equation below).

            Lscale_down(i)  &
            = Lscale_down(i)  &
            + ( - dCAPE_dz_j_plus_1 /  &
                 ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                 / gr%invrs_dzm(j)  &
            + sqrt( dCAPE_dz_j_plus_1**2  &
                    + 2.0_core_rknd * tke_i * gr%invrs_dzm(j)  &
                      * ( dCAPE_dz_j - dCAPE_dz_j_plus_1 ) )  &
                 / ( dCAPE_dz_j - dCAPE_dz_j_plus_1 )  &
                 / gr%invrs_dzm(j)

          endif

        endif

        ! Reset values for use during the next vertical level down.

        thl_par_j_plus_1 = thl_par_j
        rt_par_j_plus_1  = rt_par_j
        dCAPE_dz_j_plus_1 = dCAPE_dz_j

        tke_i = tke_i - CAPE_incr
        j = j - 1

      enddo

      ! Make Lscale_down nonlocal
      !
      ! This code makes the value of Lscale_down nonlocal.  Thus, if a parcel
      ! starting from a higher altitude can descend to altitude
      ! Lscale_down_min_alt, then a parcel starting from a lower altitude
      ! should also be able to descend to at least altitude
      ! Lscale_down_min_alt, even if the local result of Lscale_down for the
      ! parcel that started at a lower altitude is not sufficient for the
      ! parcel to reach altitude Lscale_down_min_alt.
      !
      ! For example, if it was found that a parcel starting at an altitude of
      ! 1100 m. descended to an altitude of 100 m. (an Lscale_down value of
      ! 1000 m.), then a parcel starting at an altitude of 1000 m. should also
      ! be able to descend to an altitude of at least 100 m.  If Lscale_down
      ! was found to be only 800 m. for the parcel starting at 1000 m.
      ! (resulting in the parcel only being able to descend to an altitude of
      ! 200 m.), then this code will overwrite the 800 m. value with a
      ! Lscale_down value of 900 m. (so that the parcel reaches an altitude of
      ! 100 m.).
      !
      ! This feature insures that the profile of Lscale_down will be very
      ! smooth, thus reducing numerical instability in the model.

      Lscale_down_min_alt = min( Lscale_down_min_alt, gr%zt(i)-Lscale_down(i) )

      if ( (gr%zt(i)-Lscale_down(i)) > Lscale_down_min_alt ) then
        Lscale_down(i) = gr%zt(i) - Lscale_down_min_alt
      endif

    enddo


    !!!!! Compute Lscale for every vertical level.

    do i = 2, gr%nz, 1

      ! The equation for Lscale is:
      !
      ! Lscale = sqrt( Lscale_up * Lscale_down ).

      ! Make lminh a linear function starting at value lmin at the bottom
      ! and going to zero at 500 meters in altitude.
      ! -dschanen 27 April 2007
      if( l_implemented ) then
        ! Within a host model, increase mixing length in 500 m layer above *ground*
        lminh = max( zero_threshold, Lscale_sfclyr_depth - (gr%zt(i) - gr%zm(1)) ) &
                        * ( lmin / Lscale_sfclyr_depth )
      else
        ! In standalone mode, increase mixing length in 500 m layer above *mean sea level*
        lminh = max( zero_threshold, Lscale_sfclyr_depth - gr%zt(i) ) &
                 * ( lmin / Lscale_sfclyr_depth )
      end if

      Lscale_up(i)    = max( lminh, Lscale_up(i) )
      Lscale_down(i)  = max( lminh, Lscale_down(i) )

      Lscale(i) = sqrt( Lscale_up(i)*Lscale_down(i) )

    enddo

    ! Set the value of Lscale at the upper and lower boundaries.
    Lscale(1) = Lscale(2)
    Lscale(gr%nz) = Lscale(gr%nz-1)

    ! Vince Larson limited Lscale to allow host
    ! model to take over deep convection.  13 Feb 2008.

    !Lscale = min( Lscale, 1e5 )
    Lscale = min( Lscale, Lscale_max )

    if( clubb_at_least_debug_level( 2 ) ) then

      ! Ensure that the output from this subroutine is valid.
      call length_check( Lscale, Lscale_up, Lscale_down, err_code_Lscale )
      ! Joshua Fasching January 2008

      ! Error Reporting
      ! Joshua Fasching February 2008

      if ( fatal_error( err_code_Lscale ) ) then

        write(fstderr,*) "Errors in length subroutine"

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

        ! Overwrite the last error code with this new fatal error
        err_code = err_code_Lscale

      endif ! Fatal error

    endif ! clubb_debug_level

    return

  end subroutine compute_length

!===============================================================================

end module mixing_length
