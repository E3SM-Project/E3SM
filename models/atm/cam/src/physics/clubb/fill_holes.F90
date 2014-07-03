!-----------------------------------------------------------------------
! $Id: fill_holes.F90 5623 2012-01-17 17:55:26Z connork@uwm.edu $
!===============================================================================
module fill_holes

  implicit none

  public :: fill_holes_driver, &
            vertical_avg, &
            vertical_integral

  private :: fill_holes_multiplicative

  private ! Set Default Scope

  contains

  !=============================================================================
  subroutine fill_holes_driver( num_pts, threshold, field_grid, &
                                rho_ds, rho_ds_zm, &
                                field )

    ! Description:
    ! This subroutine clips values of 'field' that are below 'threshold' as much
    ! as possible (i.e. "fills holes"), but conserves the total integrated mass
    ! of 'field'.  This prevents clipping from acting as a spurious source.
    !
    ! Mass is conserved by reducing the clipped field everywhere by a constant
    ! multiplicative coefficient.
    !
    ! This subroutine does not guarantee that the clipped field will exceed
    ! threshold everywhere; blunt clipping is needed for that.

    ! References:
    !   ``Numerical Methods for Wave Equations in Geophysical Fluid
    !     Dynamics'', Durran (1999), p. 292.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
       gr ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      num_pts  ! The number of points on either side of the hole;
               ! Mass is drawn from these points to fill the hole.  []

    real( kind = core_rknd ), intent(in) :: & 
      threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                 ! fall                           [Units vary; same as field]

    character(len=2), intent(in) :: & 
      field_grid ! The grid of the field, either zt or zm

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      rho_ds,    & ! Dry, static density on thermodynamic levels    [kg/m^3]
      rho_ds_zm    ! Dry, static density on momentum levels         [kg/m^3]

    ! Input/Output variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      field  ! The field (e.g. wp2) that contains holes [Units same as threshold]

    ! Local Variables
    integer :: & 
      k,             & ! Loop index for absolute grid level              []
      begin_idx,     & ! Lower grid level of local hole-filling range    []
      end_idx,       & ! Upper grid level of local hole-filling range    []
      upper_hf_level   ! Upper grid level of global hole-filling range   []

    !-----------------------------------------------------------------------

    ! Check whether any holes exist in the entire profile.
    ! The lowest level (k=1) should not be included, as the hole-filling scheme
    ! should not alter the set value of 'field' at the surface (for momentum
    ! level variables), or consider the value of 'field' at a level below the
    ! surface (for thermodynamic level variables).  For momentum level variables
    ! only, the hole-filling scheme should not alter the set value of 'field' at
    ! the upper boundary level (k=gr%nz).

    if ( field_grid == "zt" ) then
      ! 'field' is on the zt (thermodynamic level) grid
      upper_hf_level = gr%nz
    elseif ( field_grid == "zm" )  then
      ! 'field' is on the zm (momentum level) grid
      upper_hf_level = gr%nz-1
    endif

    if ( any( field( 2:upper_hf_level ) < threshold ) ) then

      ! Make one pass up the profile, filling holes as much as we can using
      ! nearby mass.
      ! The lowest level (k=1) should not be included in the loop, as the
      ! hole-filling scheme should not alter the set value of 'field' at the
      ! surface (for momentum level variables), or consider the value of
      ! 'field' at a level below the surface (for thermodynamic level
      ! variables).  For momentum level variables only, the hole-filling scheme
      ! should not alter the set value of 'field' at the upper boundary
      ! level (k=gr%nz).
      do k = 2+num_pts, upper_hf_level-num_pts, 1

        begin_idx = k - num_pts
        end_idx   = k + num_pts

        if ( any( field( begin_idx:end_idx ) < threshold ) ) then

          ! 'field' is on the zt (thermodynamic level) grid
          if ( field_grid == "zt" ) then
            call fill_holes_multiplicative &
                    ( begin_idx, end_idx, threshold, &
                      rho_ds(begin_idx:end_idx), gr%invrs_dzt(begin_idx:end_idx), &
                      field(begin_idx:end_idx) )
                      
          ! 'field' is on the zm (momentum level) grid
          elseif ( field_grid == "zm" )  then
            call fill_holes_multiplicative &
                    ( begin_idx, end_idx, threshold, &
                      rho_ds_zm(begin_idx:end_idx), gr%invrs_dzm(begin_idx:end_idx), &
                      field(begin_idx:end_idx) )
          endif

        endif

      enddo

      ! Fill holes globally, to maximize the chance that all holes are filled.
      ! The lowest level (k=1) should not be included, as the hole-filling
      ! scheme should not alter the set value of 'field' at the surface (for
      ! momentum level variables), or consider the value of 'field' at a level
      ! below the surface (for thermodynamic level variables).  For momentum
      ! level variables only, the hole-filling scheme should not alter the set
      ! value of 'field' at the upper boundary level (k=gr%nz).
      if ( any( field( 2:upper_hf_level ) < threshold ) ) then

        ! 'field' is on the zt (thermodynamic level) grid
        if ( field_grid == "zt" ) then
          call fill_holes_multiplicative &
                 ( 2, upper_hf_level, threshold, &
                   rho_ds(2:upper_hf_level), gr%invrs_dzt(2:upper_hf_level), &
                   field(2:upper_hf_level) )
                   
        ! 'field' is on the zm (momentum level) grid
        elseif ( field_grid == "zm" )  then
            call fill_holes_multiplicative &
                 ( 2, upper_hf_level, threshold, &
                   rho_ds_zm(2:upper_hf_level), gr%invrs_dzm(2:upper_hf_level), &
                   field(2:upper_hf_level) )
        endif

      endif

    endif  ! End overall check for existence of holes

    return

  end subroutine fill_holes_driver

  !=============================================================================
  subroutine fill_holes_multiplicative &
                 ( begin_idx, end_idx, threshold, &
                   rho, invrs_dz, &
                   field )

    ! Description:
    ! This subroutine clips values of 'field' that are below 'threshold' as much
    ! as possible (i.e. "fills holes"), but conserves the total integrated mass
    ! of 'field'.  This prevents clipping from acting as a spurious source.
    !
    ! Mass is conserved by reducing the clipped field everywhere by a constant
    ! multiplicative coefficient.
    !
    ! This subroutine does not guarantee that the clipped field will exceed
    ! threshold everywhere; blunt clipping is needed for that.

    ! References:
    ! ``Numerical Methods for Wave Equations in Geophysical Fluid
    ! Dynamics", Durran (1999), p. 292.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      begin_idx, & ! The beginning index (e.g. k=2) of the range of hole-filling 
      end_idx      ! The end index (e.g. k=gr%nz) of the range of hole-filling

    real( kind = core_rknd ), intent(in) :: & 
      threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not fall
                 !                              [Units vary; same as field]

    real( kind = core_rknd ), dimension(end_idx-begin_idx+1), intent(in) ::  & 
      rho,     &  ! Dry, static density on either thermodynamic or momentum levels   [kg/m^3]
      invrs_dz    ! Reciprocal of thermodynamic or momentum level thickness depending on whether
                  ! we're on zt or zm grid.

    ! Input/Output variable
    real( kind = core_rknd ), dimension(end_idx-begin_idx+1), intent(inout) ::  & 
      field  ! The field (e.g. wp2) that contains holes
             !                                  [Units same as threshold]

    ! Local Variables
    real( kind = core_rknd ), dimension(end_idx-begin_idx+1)  ::  & 
      field_clipped  ! The raw field (e.g. wp2) that contains no holes
                     !                          [Units same as threshold]

    real( kind = core_rknd ) ::  & 
      field_avg,  &        ! Vertical average of field [Units of field]
      field_clipped_avg, & ! Vertical average of clipped field [Units of field]
      mass_fraction        ! Coefficient that multiplies clipped field
                           ! in order to conserve mass.                      []

    !-----------------------------------------------------------------------

    ! Compute the field's vertical average, which we must conserve.
    field_avg = vertical_avg( (end_idx-begin_idx+1), rho, &
                                  field, invrs_dz )

    ! Clip small or negative values from field.
    if ( field_avg >= threshold ) then
      ! We know we can fill in holes completely
      field_clipped = max( threshold, field )
    else
      ! We can only fill in holes partly;
      ! to do so, we remove all mass above threshold.
      field_clipped = min( threshold, field )
    endif

    ! Compute the clipped field's vertical integral.
    ! clipped_total_mass >= original_total_mass
    field_clipped_avg = vertical_avg( (end_idx-begin_idx+1), rho, &
                                      field_clipped, invrs_dz )

    ! If the difference between the field_clipped_avg and the threshold is so
    ! small that it falls within numerical round-off, return to the parent
    ! subroutine without altering the field in order to avoid divide-by-zero
    ! error.
    !if ( abs(field_clipped_avg - threshold)  &
    !      < threshold*epsilon(threshold) ) then
    if ( abs(field_clipped_avg - threshold) == 0.0_core_rknd ) then
      return
    endif

    ! Compute coefficient that makes the clipped field have the same mass as the
    ! original field.  We should always have mass_fraction > 0.
    mass_fraction = ( field_avg - threshold ) / & 
                          ( field_clipped_avg - threshold )

    ! Output normalized, filled field
    field = mass_fraction * ( field_clipped - threshold )  & 
                 + threshold


    return

  end subroutine fill_holes_multiplicative

  !=============================================================================
  function vertical_avg( total_idx, rho_ds, &
                             field, invrs_dz )

    ! Description:
    ! Computes the density-weighted vertical average of a field.
    !
    ! The average value of a function, f, over a set domain, [a,b], is
    ! calculated by the equation:
    !
    ! f_avg = ( INT(a:b) f*g ) / ( INT(a:b) g );
    !
    ! as long as f is continous and g is nonnegative and integrable.  Therefore,
    ! the density-weighted (by dry, static, base-static density) vertical
    ! average value of any model field, x, is calculated by the equation:
    !
    ! x_avg|_z = ( INT(z_bot:z_top) x rho_ds dz )
    !            / ( INT(z_bot:z_top) rho_ds dz );
    !
    ! where z_bot is the bottom of the vertical domain, and z_top is the top of
    ! the vertical domain.
    !
    ! This calculation is done slightly differently depending on whether x is a
    ! thermodynamic-level or a momentum-level variable.
    !
    ! Thermodynamic-level computation:
    
    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given thermodynamic level, x and rho_ds are
    ! both thermodynamic-level variables, and delta_z(k) = zm(k) - zm(k-1).  The
    ! indices k_bot and k_top are the indices of the respective lower and upper
    ! thermodynamic levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) thermodynamic level is below ground (or below the
    ! official lower boundary at the first momentum level), so it should not
    ! count in a vertical average, whether that vertical average is used for
    ! the hole-filling scheme or for statistical purposes. Begin no lower
    ! than level k=2, which is the first thermodynamic level above ground (or
    ! above the model lower boundary).
    !
    ! For cases where hole-filling over the entire (global) vertical domain
    ! is desired, or where statistics over the entire (global) vertical
    ! domain are desired, the lower (thermodynamic-level) index of k = 2 and
    ! the upper (thermodynamic-level) index of k = gr%nz, means that the
    ! overall vertical domain will be gr%zm(gr%nz) - gr%zm(1).
    !
    !
    ! Momentum-level computation:
    !
    ! For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the
    ! numerator integral, is calculated as:
    !
    ! SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k);
    !
    ! where k is the index of the given momentum level, x and rho_ds are both
    ! momentum-level variables, and delta_z(k) = zt(k+1) - zt(k).  The indices
    ! k_bot and k_top are the indices of the respective lower and upper momentum
    ! levels involved in the integration.
    !
    ! Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral,
    ! is calculated as:
    !
    ! SUM(k_bot:k_top) rho_ds(k) delta_z(k).
    !
    ! The first (k=1) momentum level is right at ground level (or right at
    ! the official lower boundary).  The momentum level variables that call
    ! the hole-filling scheme have set values at the surface (or lower
    ! boundary), and those set values should not be changed.  Therefore, the
    ! vertical average (for purposes of hole-filling) should not include the
    ! surface level (or lower boundary level).  For hole-filling purposes,
    ! begin no lower than level k=2, which is the second momentum level above
    ! ground (or above the model lower boundary).  Likewise, the value at the
    ! model upper boundary (k=gr%nz) is also set for momentum level
    ! variables.  That value should also not be changed.
    !
    ! However, this function is also used to keep track (for statistical
    ! purposes) of the vertical average of certain variables.  In that case,
    ! the vertical average needs to be taken over the entire vertical domain
    ! (level 1 to level gr%nz).
    !
    !
    ! In both the thermodynamic-level computation and the momentum-level
    ! computation, the numerator integral is divided by the denominator integral
    ! in order to find the average value (over the vertical domain) of x.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      total_idx ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds, & ! Dry, static density on either thermodynamic or momentum levels    [kg/m^3]
      field,  & ! The field (e.g. wp2) to be vertically averaged                    [Units vary]
      invrs_dz  ! Reciprocal of thermodynamic or momentum level thickness           [1/m]
                ! depending on whether we're on zt or zm grid.
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = 1.

    ! Output variable
    real( kind = core_rknd ) :: & 
      vertical_avg  ! Vertical average of field    [Units of field]

    ! Local variables
    real( kind = core_rknd ) :: & 
      numer_integral, & ! Integral in the numerator (see description)
      denom_integral    ! Integral in the denominator (see description)
      
    real( kind = core_rknd ), dimension(total_idx) :: &
      denom_field       ! When computing the vertical integral in the denominator
                        ! there is no field variable, so create a "dummy" variable
                        ! with value of 1 to pass as an argument

    !-----------------------------------------------------------------------
    
    ! Fill array with 1's (see variable description)
    denom_field = 1.0_core_rknd

    ! Initializing vertical_avg to avoid a compiler warning.
    vertical_avg = 0.0_core_rknd
    
     
    ! Compute the numerator integral.
    ! Multiply the variable 'field' at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    ! Note:  The level thickness at level k is the distance between either
    !        momentum level k and momentum level k-1, or
    !        thermodynamic level k+1 and thermodynamic level k, depending
    !        on which field grid is being analyzed. Thus, 1.0/invrs_dz(k)
    !        is the level thickness for level k.
    ! Note:  The values of 'field' and rho_ds are passed into this function
    !        so that field(1) and rho_ds(1) are actually 'field' and rho_ds
    !        at the level k = 1.
       
    numer_integral = vertical_integral( total_idx, rho_ds(1:total_idx), &
                                            field(1:total_idx), invrs_dz(1:total_idx) )
    
    ! Compute the denominator integral.
    ! Multiply rho_ds at level k by the level thickness
    ! at level k.  Then, sum over all vertical levels.
    denom_integral = vertical_integral( total_idx, rho_ds(1:total_idx), &
                                            denom_field(1:total_idx), invrs_dz(1:total_idx) )

    ! Find the vertical average of 'field'.
    vertical_avg = numer_integral / denom_integral

    return
  end function vertical_avg

  !=============================================================================
  pure function vertical_integral( total_idx, rho_ds, &
                                       field, invrs_dz )

    ! Description:
    ! Computes the vertical integral. rho_ds, field, and invrs_dz must all be
    ! of size total_idx and should all start at the same index.
    ! 
    
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: & 
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      invrs_dz   ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = begin_idx.

    ! Local variables
    real( kind = core_rknd ) :: &
      vertical_integral ! Integral in the numerator (see description)

    !-----------------------------------------------------------------------

    !  Assertion checks: that begin_idx <= gr%nz - 1
    !                    that end_idx   >= 2
    !                    that begin_idx <= end_idx


    ! Initializing vertical_integral to avoid a compiler warning.
    vertical_integral = 0.0_core_rknd

    ! Compute the integral.
    ! Multiply the field at level k by rho_ds at level k and by
    ! the level thickness at level k.  Then, sum over all vertical levels.
    ! Note:  The values of the field and rho_ds are passed into this function
    !        so that field(1) and rho_ds(1) are actually the field and rho_ds
    !        at level k_start.
    vertical_integral = sum( field * rho_ds / invrs_dz )

    return
  end function vertical_integral

!===============================================================================

end module fill_holes
