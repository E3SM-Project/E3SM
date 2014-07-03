!-----------------------------------------------------------------------
! $Id: hydrostatic_module.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
!===============================================================================
module hydrostatic_module

  implicit none

  private ! Default Scope

  public :: hydrostatic, &
            inverse_hydrostatic

  private :: calc_exner_const_thvm, &
             calc_exner_linear_thvm, &
             calc_z_linear_thvm

  contains

!===============================================================================
  subroutine hydrostatic( thvm, p_sfc, &
                          p_in_Pa, p_in_Pa_zm, &
                          exner, exner_zm, &
                          rho, rho_zm )

    ! Description:
    ! This subroutine integrates the hydrostatic equation.
    !
    ! The hydrostatic equation is of the form:
    !
    ! dp/dz = - rho * grav.
    !
    ! This equation can be re-written in terms of d(exner)/dz, such that:
    !
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) } / { R_d * rho } ] * d(exner)/dz
    ! = - grav / C_p;
    !
    ! which can also be expressed as:
    !
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) } / { R_d * rho_d * ( 1 + r_v + r_c ) } ]
    !    * d(exner)/dz
    ! = - grav / C_p.
    !
    ! Furthermore, the moist equation of state can be written as:
    !
    ! theta =
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) }
    !   / { R_d * rho_d * ( 1 + (R_v/R_d)*r_v ) } ].
    !
    ! The relationship between theta and theta_v (including water vapor and
    ! cloud water) is:
    !
    ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
    !
    ! which, when substituted into the above equation, changes the equation of
    ! state to:
    !
    ! theta_v =
    ! [ { p0^(R_d/C_p) * p^(C_v/C_p) }
    !   / { R_d * rho_d * ( 1 + r_v + r_c ) } ].
    !
    ! This equation is substituted into the d(exner)/dz form of the hydrostatic
    ! equation, resulting in:
    !
    ! theta_v * d(exner)/dz = - grav / C_p;
    !
    ! which can be re-written as:
    !
    ! d(exner)/dz = - grav / ( C_p * theta_v ).
    !
    ! This subroutine integrates the above equation to solve for exner, such
    ! that:
    !
    ! INT(exner_1:exner_2) d(exner) =
    ! - ( grav / C_p ) * INT(z_1:z_2) ( 1 / theta_v ) dz.
    !
    !
    ! The resulting value of exner is used to calculate pressure.  Then, the
    ! values of pressure, exner, and theta_v can be used to calculate density.

    ! References:
    !
    !------------------------------------------------------------------------

    use constants_clubb, only: & 
        kappa,  & ! Variable(s)
        p0, & 
        Rd, &
        zero_threshold

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt,  & ! Procedure(s)
        zt2zm

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      p_sfc    ! Pressure at the surface                     [Pa]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      thvm    ! Virtual potential temperature               [K]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  & 
      p_in_Pa,    & ! Pressure (thermodynamic levels)         [Pa]
      p_in_Pa_zm, & ! Pressure on momentum levels             [Pa]
      exner,      & ! Exner function (thermodynamic levels)   [-]
      exner_zm,   & ! Exner function on momentum levels       [-]
      rho,        & ! Density (thermodynamic levels)          [kg/m^3]
      rho_zm        ! Density on momentum levels              [kg/m^3]

    !  Local Variables
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      thvm_zm       ! Theta_v interpolated to momentum levels  [K]

    real( kind = core_rknd ) :: &
      dthvm_dz      ! Constant d(thvm)/dz between successive levels   [K/m]

    integer :: k

    ! Interpolate thvm from thermodynamic to momentum levels.  Linear
    ! interpolation is used, except for the uppermost momentum level, where a
    ! linear extension is used.  Since thvm is considered to either be constant
    ! or vary linearly over the depth of a grid level, this interpolation is
    ! consistent with the rest of this code.
    thvm_zm = zt2zm( thvm )

    ! Exner is defined on thermodynamic grid levels except for the value at
    ! index 1.  Since thermodynamic level 1 is below the surface, it is
    ! disregarded, and the value of exner(1) corresponds to surface value, which
    ! is actually at momentum level 1.
    exner(1) = ( p_sfc/p0 )**kappa
    exner_zm(1) = ( p_sfc/p0 )**kappa

    ! Consider the value of exner at thermodynamic level (2) to be based on
    ! a constant thvm between thermodynamic level (2) and momentum level (1),
    ! which is the surface or model lower boundary.  Since thlm(1) is set equal
    ! to thlm(2), the values of thvm are considered to be basically constant
    ! near the ground.
    exner(2) &
    = calc_exner_const_thvm( thvm(2), gr%zt(2), gr%zm(1), exner(1) )

    ! Given the value of exner at thermodynamic level k-1, and considering
    ! thvm to vary linearly between its values at thermodynamic levels k
    ! and k-1, the value of exner can be found at thermodynamic level k,
    ! as well as at intermediate momentum level k-1.
    do k = 3, gr%nz

      dthvm_dz = gr%invrs_dzm(k-1) * ( thvm(k) - thvm(k-1) )

      if ( dthvm_dz /= 0.0_core_rknd ) then

        exner(k) &
        = calc_exner_linear_thvm( thvm(k-1), dthvm_dz, &
                                  gr%zt(k-1), gr%zt(k), exner(k-1) )

        exner_zm(k-1) &
        = calc_exner_linear_thvm( thvm(k-1), dthvm_dz, &
                                  gr%zt(k-1), gr%zm(k-1), exner(k-1) )

      else ! dthvm_dz = 0

        exner(k) &
        = calc_exner_const_thvm &
             ( thvm(k), gr%zt(k), gr%zt(k-1), exner(k-1) )

        exner_zm(k-1) &
        = calc_exner_const_thvm &
             ( thvm(k), gr%zm(k-1), gr%zt(k-1), exner(k-1) )

      endif

    enddo ! k = 3, gr%nz

    ! Find the value of exner_zm at momentum level gr%nz by using a linear
    ! extension of thvm from the two thermodynamic level immediately below
    ! momentum level gr%nz.
    dthvm_dz = ( thvm_zm(gr%nz) - thvm(gr%nz) ) &
               / ( gr%zm(gr%nz) - gr%zt(gr%nz) )

    if ( dthvm_dz /= 0.0_core_rknd ) then

      exner_zm(gr%nz) &
      = calc_exner_linear_thvm &
           ( thvm(gr%nz), dthvm_dz, &
             gr%zt(gr%nz), gr%zm(gr%nz), exner(gr%nz) )

    else ! dthvm_dz = 0

      exner_zm(gr%nz) &
      = calc_exner_const_thvm &
           ( thvm(gr%nz), gr%zm(gr%nz), gr%zt(gr%nz), exner(gr%nz) )

    endif

    ! Calculate pressure based on the values of exner.

    do k = 1, gr%nz
      p_in_Pa(k) = p0 * exner(k)**( 1._core_rknd/kappa )
      p_in_Pa_zm(k) = p0 * exner_zm(k)**( 1._core_rknd/kappa )
    enddo

    ! Calculate density based on pressure, exner, and thvm.

    do k = 1, gr%nz
      rho(k) = p_in_Pa(k) / ( Rd * thvm(k) * exner(k) )
      rho_zm(k) = p_in_Pa_zm(k) / ( Rd * thvm_zm(k) * exner_zm(k) )
    enddo


    return
  end subroutine hydrostatic

!===============================================================================
  subroutine inverse_hydrostatic( p_sfc, zm_init, nlevels, thvm, exner, &
                                  z )

    ! Description:
    ! Subprogram to integrate the inverse of hydrostatic equation

    ! References:
    !
    !------------------------------------------------------------------------

    use constants_clubb, only: &
        p0,     & ! Constant(s)
        kappa,  &
        fstderr

    use interpolation, only: &
        binary_search ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      p_sfc,    & ! Pressure at the surface      [Pa]
      zm_init    ! Altitude at the surface      [m]

    integer, intent(in) ::  &
      nlevels  ! Number of levels in the sounding [-]

    real( kind = core_rknd ), intent(in), dimension(nlevels) ::  & 
      thvm,  & ! Virtual potential temperature   [K]
      exner    ! Exner function                  [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(nlevels) ::  & 
      z        ! Height                    [m]

    !  Local Variables
    integer :: k

    real( kind = core_rknd ), dimension(nlevels) ::  &
      ref_z_snd  ! Altitude minus altitude of the lowest sounding level  [m]

    real( kind = core_rknd ), dimension(nlevels) ::  &
      exner_reverse_array  ! Array of exner snd. values in reverse order [-]

    real( kind = core_rknd ) ::  &
      exner_sfc,    & ! Value of exner at the surface                    [-]
      ref_z_sfc,    & ! Alt. diff between surface and lowest snd. level  [m]
      z_snd_bottom, & ! Altitude of the bottom of the input sounding     [m]
      dthvm_dexner    ! Constant rate of change of thvm with respect to
    ! exner between sounding levels k-1 and k          [K]

    integer ::  &
      rev_low_idx, &
      low_idx, &
      high_idx


    ! Variable ref_z_sfc is initialized to 0.0 to avoid a compiler warning.
    ref_z_sfc = 0.0_core_rknd

    ! The variable ref_z_snd is the altitude of each sounding level compared to
    ! the altitude of the lowest sounding level.  Thus, the value of ref_z_snd
    ! at sounding level 1 is 0.  The lowest sounding level may or may not be
    ! right at the surface, and therefore an adjustment may be required to find
    ! the actual altitude above ground.
    ref_z_snd(1) = 0.0_core_rknd

    do k = 2, nlevels

      ! The value of thvm is given at two successive sounding levels.  For
      ! purposes of achieving a quality estimate of altitude at each pressure
      ! sounding level, the value of thvm is considered to vary linearly
      ! with respect to exner between two successive sounding levels.  Thus,
      ! there is a constant d(thvm)/d(exner) between the two successive
      ! sounding levels.  If thvm is constant, then d(thvm)/d(exner) is 0.
      dthvm_dexner = ( thvm(k) - thvm(k-1) ) / ( exner(k) - exner(k-1) )

      ! Calculate the value of the reference height at sounding level k, based
      ! the value of thvm at sounding level k-1, the constant value of
      ! d(thvm)/d(exner), the value of exner at sounding levels k-1 and k, and
      ! the reference altitude at sounding level k-1.
      ref_z_snd(k) &
      = calc_z_linear_thvm( thvm(k-1), dthvm_dexner, &
                            exner(k-1), exner(k), ref_z_snd(k-1) )

    enddo

    ! Find the actual (above ground) altitude of the sounding levels from the
    ! reference altitudes.

    ! The pressure at the surface (or model lower boundary), p_sfc, is found at
    ! the altitude of the surface (or model lower boundary), zm_init.

    ! Find the value of exner at the surface from the pressure at the surface.
    exner_sfc = ( p_sfc / p0 )**kappa

    ! Find the value of exner_sfc compared to the values of exner in the exner
    ! sounding profile.

    if ( exner_sfc < exner(nlevels) ) then

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is less than all the
      ! values of exner in the sounding (and thus the surface is located above
      ! all the levels of the sounding), then there is insufficient information
      ! to run the model.  Stop the run.

      write(fstderr,*) "The entire sounding is below the model surface."
      stop

    elseif ( exner_sfc > exner(1) ) then

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is greater than all the
      ! values of exner in the sounding (and thus the surface is located below
      ! all the levels of the sounding), use a linear extension of thvm to find
      ! thvm at the surface.  Thus, d(thvm)/d(exner) is the same as its value
      ! between sounding levels 1 and 2.  If the surface is so far below the
      ! sounding that gr%zt(2) is below the first sounding level, the code in
      ! subroutine read_sounding (found in sounding.F90) will stop the run.

      ! Calculate the appropriate d(thvm)/d(exner).
      dthvm_dexner = ( thvm(2) - thvm(1) ) / ( exner(2) - exner(1) )

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_z_linear_thvm( thvm(1), dthvm_dexner, &
                            exner(1), exner_sfc, ref_z_snd(1) )

    else  ! exner(nlevels) < exner_sfc < exner(1)

      ! Since the values of exner decrease monotonically with height (and thus
      ! with sounding level), the value of exner_sfc is between two values of
      ! exner (at some levels k-1 and k) in the sounding, and the value of
      ! d(thvm)/d(exner) is the same as between those two levels in the above
      ! calculation.

      ! The value of exner_sfc is between two levels of the exner sounding.
      ! Find the index of the lower level.

      ! In order to use the binary search, the array must be sorted from least
      ! value to greatest value.  Since exner decreases with altitude (and
      ! vertical level), the array that is sent to function binary_search must
      ! be the exact reverse of exner.
      ! Thus, exner(1) becomes exner_reverse_array(nlevels), exner(nlevels)
      ! becomes exner_reverse_array(1), etc.
      do k = 1, nlevels, 1
        exner_reverse_array(k) = exner(nlevels-k+1)
      enddo
      ! The output from the binary search yields the first value in the
      ! exner_reverse_array that is greater than or equal to exner_sfc.  Thus,
      ! in regards to the regular exner array, this is the reverse index of
      ! the lower sounding level for exner_sfc.  For example, if exner_sfc
      ! is found between exner(1) and exner(2), the binary search for exner_sfc
      ! in regards to exner_reverse_index will return a value of nlevels.
      ! Once the actual lower level index is calculated, the result will be 1.
      rev_low_idx = binary_search( nlevels, exner_reverse_array, exner_sfc )

      ! Find the lower level index for the regular exner profile from the
      ! lower level index for the reverse exner profile.
      low_idx = nlevels - rev_low_idx + 1

      ! Find the index of the upper level.
      high_idx = low_idx + 1

      ! Calculate the appropriate d(thvm)/d(exner).
      dthvm_dexner = ( thvm(high_idx) - thvm(low_idx) )  &
                       /  ( exner(high_idx) - exner(low_idx) )

      ! Calculate the difference between the altitude of the surface (or model
      ! lower boundary) and the altitude of the lowest level of the sounding.
      ref_z_sfc  &
      = calc_z_linear_thvm( thvm(low_idx), dthvm_dexner, &
                            exner(low_idx), exner_sfc, ref_z_snd(low_idx) )

    endif  ! exner_sfc

    ! Find the altitude of the bottom of the sounding.
    z_snd_bottom = zm_init - ref_z_sfc

    ! Calculate the sounding altitude profile based
    ! on z_snd_bottom and ref_z_snd.
    do k = 1, nlevels, 1
      z(k) = z_snd_bottom + ref_z_snd(k)
    enddo


    return
  end subroutine inverse_hydrostatic

!===============================================================================
  pure function calc_exner_const_thvm( thvm, z_2, z_1, exner_1 ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a constant thvm over the depth of the
    ! level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * thvm).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) (1/thvm) dz.
    !
    ! Since thvm is considered to be a constant over the depth of the layer
    ! between z_1 and z_2, the equation can be written as:
    !
    ! INT(exner_1:exner_2) d(exner) = - grav / ( Cp * thvm ) INT(z_1:z_2) dz.
    !
    ! Solving the integral:
    !
    ! exner_2 = exner_1 - [ grav / ( Cp * thvm ) ] * ( z_2 - z_1 ).

    ! References:
    !-------------------------------------------------------------------

    use constants_clubb, only: &
        grav,  & ! Gravitational acceleration                  [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure [J/(kg*K)]

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thvm,    & ! Constant value of thvm over the layer    [K]
      z_2,     & ! Altitude at the top of the layer         [m]
      z_1,     & ! Altitude at the bottom of the layer      [m]
      exner_1    ! Exner at the bottom of the layer         [-]

    ! Return Variable
    real( kind = core_rknd ) :: exner_2  ! Exner at the top of the layer        [-]

    ! Calculate exner at top of the layer.
    exner_2 = exner_1 - ( grav / ( Cp * thvm ) ) * ( z_2 - z_1 )

    return
  end function calc_exner_const_thvm

!===============================================================================
  pure function calc_exner_linear_thvm( thvm_km1, dthvm_dz, &
                                        z_km1, z_2, exner_km1  ) &
  result( exner_2 )

    ! Description:
    ! This function solves for exner at a level, given exner at another level,
    ! the altitudes of both levels, and a value of thvm that is considered to
    ! vary linearly over the depth of the level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * thvm).
    !
    ! This equation is integrated to solve for exner, such that:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) (1/thvm) dz.
    !
    ! The value of thvm is considered to vary linearly (with respect to height)
    ! over the depth of the level (resulting in a constant d(thvm)/dz over the
    ! depth of the level).  The entire level between z_1 and z_2 must be
    ! encompassed between two levels with two known values of thvm.  The value
    ! of thvm at the upper level (z_up) is called thvm_up, and the value of thvm
    ! at the lower level (z_low) is called thvm_low.  Again, the values of thvm
    ! at all interior altitudes, z_low <= z_1 < z <= z_2 <= z_up, behave
    ! linearly between thvm_low and thvm_up, such that:
    !
    ! thvm(z)
    ! = [ ( thvm_up - thvm_low ) / ( z_up - z_low ) ] * ( z - z_low)
    !   + thvm_low
    ! = [ d(thvm)/dz ] * ( z - z_low ) + thvm_low
    ! = C_a*z + C_b;
    !
    ! where:
    !
    ! C_a
    ! = ( thvm_up - thvm_low ) / ( z_up - z_low )
    ! = d(thvm)/dz;
    !
    ! and:
    !
    ! C_b
    ! = thvm_low - [ ( thvm_up - thvm_low ) / ( z_up - z_low ) ] * z_low
    ! = thvm_low - [ d(thvm)/dz ] * z_low.
    !
    ! The integral becomes:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) [ 1 / ( C_a*z + C_b ) ] dz.
    !
    ! Performing a u-substitution ( u = C_a*z + C_b ), the equation becomes:
    !
    ! INT(exner_1:exner_2) d(exner)
    ! = - ( grav / Cp ) * ( 1 / C_a ) INT(z=z_1:z=z_2) (1/u) du.
    !
    ! Solving the integral, and then re-substituting for u:
    !
    ! exner_2 = exner_1
    !           - ( grav / Cp ) * ( 1 / C_a )
    !             * ln [ ( C_a*z_2 + C_b ) / ( C_a*z_1 + C_b ) ].
    !
    ! Re-substituting for C_a and C_b:
    !
    ! exner_2
    ! = exner_1
    !   - ( grav / Cp ) * ( 1 / {d(thvm)/dz} )
    !     * ln [   ( {d(thvm)/dz}*z_2 + thvm_low - {d(thvm)/dz}*z_low )
    !            / ( {d(thvm)/dz}*z_1 + thvm_low - {d(thvm)/dz}*z_low ) ].
    !
    ! This equation is used to calculate exner_2 using exner_1, which is at the
    ! same level as z_1.  Furthermore, thvm_low and z_low are taken from the
    ! same level as z_1 and exner_1.  Thus, z_1 = z_low.  Therefore:
    !
    ! exner_2
    ! = exner_low
    !   - ( grav / Cp ) * ( 1 / {d(thvm)/dz} )
    !     * ln [ ( thvm_low + {d(thvm)/dz}*(z_2-z_low) ) / thvm_low ].
    !
    ! Considering either a thermodynamic or sounding level k-1 as the low level
    ! in the integration, and that thvm varies linearly between level k-1 and
    ! level k:
    !
    ! exner_2
    ! = exner(k-1)
    !   - ( grav / Cp ) * ( 1 / {d(thvm)/dz} )
    !     * ln [ ( thvm(k-1) + {d(thvm)/dz}*(z_2-z(k-1)) ) / thvm(k-1) ];
    !
    ! where:
    !
    ! d(thvm)/dz = ( thvm(k) - thvm(k-1) ) / ( z(k) - z(k-1) );
    !
    ! and where z(k-1) < z_2 <= z(k); and {d(thvm)/dz} /= 0.  If the value of
    ! {d(thvm)/dz} is 0, then thvm is considered to be a constant over the depth
    ! of the level.  The appropriate equation is found in pure function
    ! calc_exner_const_thvm.

    ! References:
    !-------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        grav,  & ! Gravitational acceleration                   [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure  [J/(kg*K)]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thvm_km1,  & ! Value of thvm at level k-1                     [K]
      dthvm_dz,  & ! Constant d(thvm)/dz between levels k-1 and k   [K/m]
      z_km1,     & ! Altitude at level k-1                          [m]
      z_2,       & ! Altitude at the top of the layer               [m]
      exner_km1    ! Exner at level k-1                             [-]

    ! Return Variable
    real( kind = core_rknd ) :: exner_2 ! Exner at the top of the layer                 [-]

    ! Calculate exner at the top of the layer.
    exner_2  &
    = exner_km1 &
      - ( grav / Cp ) * ( 1.0_core_rknd / dthvm_dz )  &
        * log(  ( thvm_km1 + dthvm_dz * ( z_2 - z_km1 ) )  /  thvm_km1  )

    return
  end function calc_exner_linear_thvm

!===============================================================================
  pure function calc_z_linear_thvm( thvm_km1, dthvm_dexner, &
                                    exner_km1, exner_2, z_km1 ) &
  result( z_2 )

    ! Description:
    ! This function solves for z (altitude) at a level, given altitude at
    ! another level, the values of exner at both levels, and a value of thvm
    ! that is considered to vary linearly over the depth of the level.
    !
    ! The derivative of exner is given by the following equation:
    !
    ! d(exner)/dz = - grav / (Cp * thvm).
    !
    ! This equation is integrated to solve for z, such that:
    !
    ! INT(exner_1:exner_2) thvm d(exner) = - ( grav / Cp ) INT(z_1:z_2) dz.
    !
    ! The value of thvm is considered to vary linearly (with respect to exner)
    ! over the depth of the level (resulting in a constant d(thvm)/d(exner) over
    ! the depth of the level).  The entire level between exner_1 and exner_2
    ! must be encompassed between two levels with two known values of thvm.  The
    ! value of thvm at the upper level (exner_up) is called thvm_up, and the
    ! value of thvm at the lower level (exner_low) is called thvm_low.  Again,
    ! the values of thvm at all interior exner levels,
    ! exner_low >= exner_1 > exner >= exner_2 >= exner_up, behave linearly
    ! between thvm_low and thvm_up, such that:
    !
    ! thvm(exner)
    ! = [ ( thvm_up - thvm_low ) / ( exner_up - exner_low ) ]
    !     * ( exner - exner_low )
    !   + thvm_low
    ! = [ d(thvm)/d(exner) ] * ( exner - exner_low ) + thvm_low
    ! = C_a*z + C_b;
    !
    ! where:
    !
    ! C_a
    ! = ( thvm_up - thvm_low ) / ( exner_up - exner_low )
    ! = d(thvm)/d(exner);
    !
    ! and:
    !
    ! C_b
    ! = thvm_low
    !   - [ ( thvm_up - thvm_low ) / ( exner_up - exner_low ) ] * exner_low
    ! = thvm_low - [ d(thvm)/d(exner) ] * exner_low.
    !
    ! The integral becomes:
    !
    ! INT(exner_1:exner_2) ( C_a*exner + C_b ) d(exner)
    ! = - ( grav / Cp ) INT(z_1:z_2) dz.
    !
    ! Solving the integral:
    !
    ! z_2
    ! = z_1
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner_1}^2 )
    !          + ( thvm_low - {d(thvm)/d(exner)} * exner_low )
    !            * ( exner_2 - exner_1 )  ].
    !
    ! This equation is used to calculate z_2 using z_1, which is at the same
    ! level as exner_1.  Furthermore, thvm_low and exner_low are taken from the
    ! same level as exner_1 and z_1.  Thus, exner_1 = exner_low.  Therefore:
    !
    ! z_2
    ! = z_low
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner_low}^2 )
    !          + ( thvm_low - {d(thvm)/d(exner)} * exner_low )
    !            * ( exner_2 - exner_low )  ].
    !
    ! Considering a sounding level k-1 as the low level in the integration, and
    ! that thvm varies linearly (with respect to exner) between level k-1 and
    ! level k:
    !
    ! z_2
    ! = z(k-1)
    !   - ( Cp / grav )
    !     * [    (1/2) * {d(thvm)/d(exner)} * ( {exner_2}^2 - {exner(k-1)}^2 )
    !          + ( thvm(k-1) - {d(thvm)/d(exner)} * exner(k-1) )
    !            * ( exner_2 - exner(k-1) )  ];
    !
    ! where:
    !
    ! d(thvm)/d(exner)
    ! = ( thvm(k) - thvm(k-1) ) / ( exner(k) - exner(k-1) );
    !
    ! and where exner(k-1) > exner_2 >= exner(k).  If the value of
    ! d(thvm)/d(exner) is 0, then thvm is considered to be a constant over the
    ! depth of the level, and the equation will reduce to:
    !
    ! z_2 = z(k-1) - ( Cp / grav ) * thvm(k-1) * ( exner_2 - exner(k-1) ).
    !
    !
    ! IMPORTANT NOTE:
    !
    ! CLUBB is an altitude-based model.  All linear interpolations (and
    ! extensions) are based on considering a variable to change linearly with
    ! respect to altitude, rather than with respect to exner.  An exception is
    ! made here to calculate the altitude of a sounding level based on a
    ! sounding given in terms of a pressure coordinate rather than a height
    ! coordinate.  After the altitude of the sounding level has been calculated,
    ! the values of the sounding variables are interpolated onto the model grid
    ! linearly with respect to altitude.  Therefore, considering a variable to
    ! change linearly with respect to exner is not consistent with the rest of
    ! the model code, but provides for a better estimation of the altitude of
    ! the sounding levels (than simply considering thvm to be constant over the
    ! depth of the sounding level).

    ! References:
    !-------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        grav,  & ! Gravitational acceleration                   [m/s^2]
        Cp       ! Specific heat of dry air at const. pressure  [J/(kg*K)]

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      thvm_km1,     & ! Value of thvm at sounding level k-1                 [K]
      dthvm_dexner, & ! Constant d(thvm)/d(exner) between levels k-1 and k  [K]
      exner_km1,    & ! Value of exner at sounding level k-1                [-]
      exner_2,      & ! Value of exner at the top of the layer              [-]
      z_km1           ! Altitude at sounding level k-1                      [m]

    ! Return Variable
    real( kind = core_rknd ) :: z_2       ! Altitude at the top of the layer                    [m]

    ! Calculate z_2 at the top of the layer.
    z_2  &
    = z_km1  &
      - ( Cp / grav )  &
        * (   0.5_core_rknd * dthvm_dexner * ( exner_2**2 - exner_km1**2 )  &
            + ( thvm_km1 - dthvm_dexner * exner_km1 )  &
              * ( exner_2 - exner_km1 )  &
          )

    return
  end function calc_z_linear_thvm

!===============================================================================

end module hydrostatic_module
