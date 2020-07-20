!-----------------------------------------------------------------------
! $Id: fill_holes.F90 7315 2014-09-30 20:49:54Z schemena@uwm.edu $
!===============================================================================
module fill_holes

  implicit none

  public :: fill_holes_driver, &
            fill_holes_vertical, &
            hole_filling_hm_one_lev, &
            fill_holes_hydromet, &
            fill_holes_wv, &
            vertical_avg, &
            vertical_integral, &
            setup_stats_indices

  private :: fill_holes_multiplicative

  private ! Set Default Scope

  contains

  !=============================================================================
  subroutine fill_holes_vertical( num_draw_pts, threshold, field_grid, &
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
      num_draw_pts  ! The number of points on either side of the hole;
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
      do k = 2+num_draw_pts, upper_hf_level-num_draw_pts, 1

        begin_idx = k - num_draw_pts
        end_idx   = k + num_draw_pts

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

  end subroutine fill_holes_vertical

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

  subroutine hole_filling_hm_one_lev( num_hm_fill, hm_one_lev, & ! Intent(in)
                                   hm_one_lev_filled ) ! Intent(out)

  ! Description:
  ! Fills holes between same-phase (i.e. either liquid or frozen) hydrometeors for
  ! one height level.
  !
  ! Warning: Do not input hydrometeors of different phases, e.g. liquid and frozen.
  ! Otherwise heat will not be conserved.
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: num_hm_fill ! number of hydrometeors involved

    real(kind = core_rknd), dimension(num_hm_fill), intent(in) :: hm_one_lev

    ! Output Variables
    real(kind = core_rknd), dimension(num_hm_fill), intent(out) :: hm_one_lev_filled

    ! Local Variables
    integer :: num_neg_hm ! number of holes

    real(kind = core_rknd) :: &
      total_hole, & ! Size of the hole ( missing mass, less than 0 )
      total_mass    ! Total mass to fill the hole
      ! total mass of water substance = total_mass + total_hole

    integer :: i ! loop iterator

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Initialization
    hm_one_lev_filled = 0._core_rknd
    total_hole = 0._core_rknd
    total_mass = 0._core_rknd
    num_neg_hm = 0

    ! Determine the total size of the hole and the number of neg. hydrometeors
    ! and the total mass of hole filling material
    do i=1, num_hm_fill
!       print *, "hm_one_lev(",i,") = ", hm_one_lev(i)
       if ( hm_one_lev(i) < zero ) then
          total_hole = total_hole + hm_one_lev(i) ! less than zero
          num_neg_hm = num_neg_hm + 1
       else
          total_mass = total_mass + hm_one_lev(i)
       endif

    enddo

!    print *, "total_hole = ", total_hole
!    print *, "total_mass = ", total_mass
!    print *, "num_neg_hm = ", num_neg_hm

    ! There is no water substance at all to fill the hole
    if ( total_mass == zero ) then

       if ( clubb_at_least_debug_level(2) ) then
          print *, "Warning: One level hole filling was not successful! total_mass = 0"
       endif

       hm_one_lev_filled = hm_one_lev

       return
    endif

    ! Fill the holes and adjust the remaining quantities:
    ! hm_filled(i) = 0, if hm(i) < 0
    ! or
    ! hm_filled(i) = (1 + total_hole/total_mass)*hm(i), if hm(i) > 0
    do i=1, num_hm_fill

       ! if there is not enough material, fill the holes partially with all the material available
       if ( abs(total_hole) > total_mass ) then

          if ( clubb_at_least_debug_level(2) ) then
             print *, "Warning: One level hole was not able to fill holes completely!" // &
                      " The holes were filled partially. |total_hole| > total_mass"
          endif

          hm_one_lev_filled(i) = min(hm_one_lev(i), zero) * ( one + total_mass / total_hole )

       else ! fill holes completely
          hm_one_lev_filled(i) = max(hm_one_lev(i), zero) * ( one + total_hole / total_mass )

       endif

    enddo

    ! Assertion checks (water substance conservation, non-negativity)
    if ( clubb_at_least_debug_level( 2 ) ) then

       if ( sum( hm_one_lev ) /= sum(hm_one_lev_filled) ) then
          print *, "Warning: Hole filling was not conservative!"
       endif

       if ( any( hm_one_lev_filled < zero ) ) then
          print *, "Warning: Hole filling failed! A hole could not be filled."
       endif

    endif

    return

  end subroutine hole_filling_hm_one_lev
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine fill_holes_hydromet( nz, hydromet_dim, hydromet, & ! Intent(in)
                                  hydromet_filled ) ! Intent(out)

  ! Description:
  ! Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  ! The hole filling conserves water substance between all same-phase (frozen or liquid)
  ! hydrometeors at each height level.
  !
  ! Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  !
  ! Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use array_index, only: &
        l_frozen_hm, & ! Variable(s)
        l_mix_rat_hm

    use constants_clubb, only: &
        zero

    implicit none

    ! Input Variables
    integer, intent(in) :: hydromet_dim, nz

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_filled

    ! Local Variables
    integer :: i,j ! Loop iterators

    integer :: num_frozen_hm ! Number of frozen hydrometeor mixing ratios

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      hydromet_frozen,       & ! Frozen hydrometeor mixing ratios
      hydromet_frozen_filled   ! Frozen hydrometeor mixing ratios after hole filling

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Determine the number of frozen hydrometeor mixing ratios
    num_frozen_hm = 0
    do i=1,hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          num_frozen_hm = num_frozen_hm + 1
       endif
    enddo

    ! Allocation
    allocate( hydromet_frozen(nz,num_frozen_hm) )
    allocate( hydromet_frozen_filled(nz,num_frozen_hm) )

    ! Determine frozen hydrometeor mixing ratios
    j = 1
    do i = 1,hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          hydromet_frozen(:,j) = hydromet(:,i)
          j = j+1
       endif
    enddo

    ! Fill holes for the frozen hydrometeors
    do i=1,nz
       if ( any( hydromet_frozen(i,:) < zero ) ) then
          call hole_filling_hm_one_lev( num_frozen_hm, hydromet_frozen(i,:), & ! Intent(in)
                                     hydromet_frozen_filled(i,:) ) ! Intent(out)
       else
          hydromet_frozen_filled(i,:) = hydromet_frozen(i,:)
       endif
    enddo

    ! Setup the filled hydromet array
    j = 1
    do i=1, hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          hydromet_filled(:,i) = hydromet_frozen_filled(:,j)
          j = j+1
       else
          hydromet_filled(:,i) = hydromet(:,i)
       endif
    enddo

    !!! Here we could do the same hole filling for all the liquid phase hydrometeors

    return
  end subroutine fill_holes_hydromet
  !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
  subroutine fill_holes_wv( nz, dt, exner, hydromet_name, & ! Intent(in)
                            rvm_mc, thlm_mc, hydromet )! Intent(inout)

  ! Description:
  ! Fills holes using the cloud water mixing ratio from the current height level.
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero_threshold, &
        Lv, &
        Ls, &
        Cp

    implicit none

    ! Input Variables
    integer, intent(in) :: nz

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Timestep         [s]

    character(len=10), intent(in) :: hydromet_name

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner        ! Exner function                            [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      hydromet, &  ! Hydrometeor array                         [units vary]
      rvm_mc, &
      thlm_mc

    ! Local Variables
    integer :: k ! Loop iterator

    real( kind = core_rknd ) :: rvm_clip_tndcy
  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do k = 2, nz, 1

       if ( hydromet(k) < zero_threshold ) then

          ! Set rvm_clip_tndcy to the time tendency applied to vapor and removed
          ! from the hydrometeor.
          rvm_clip_tndcy = hydromet(k) / dt

          ! Adjust the tendency rvm_mc accordingly
          rvm_mc(k) = rvm_mc(k) + rvm_clip_tndcy

          ! Adjust the tendency of thlm_mc according to whether the
          ! effect is an evaporation or sublimation tendency.
          select case ( trim( hydromet_name ) )
          case( "rrm" )
             thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * ( Lv / ( Cp*exner(k) ) )
          case( "rim", "rsm", "rgm" )
             thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * ( Ls / ( Cp*exner(k) ) )
          case default
             stop "Fatal error in microphys_driver"
          end select

          ! Set the mixing ratio to 0
          hydromet(k) = zero_threshold

       endif ! hydromet(k,i) < 0

    enddo ! k = 2..gr%nz

    return
  end subroutine fill_holes_wv
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine fill_holes_driver( nz, dt, hydromet_dim,        & ! Intent(in)
                                l_fill_holes_hm,             & ! Intent(in)
                                rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
                                thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

  ! Description:
  ! Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  ! The hole filling conserves water substance between all same-phase (frozen or liquid)
  ! hydrometeors at each height level.
  !
  ! Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  !
  ! Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use grid_class, only: &
        gr  ! Variable(s)

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use constants_clubb, only: &
        pi,              &
        four_thirds,     &
        one,             &
        zero,            &
        zero_threshold,  &
        Lv,              &
        Ls,              &
        Cp,              &
        rho_lw,          &
        rho_ice,         &
        fstderr

    use array_index, only: &
        hydromet_list, & ! Names of the hydrometeor species
        hydromet_tol

    use array_index, only: &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm

    use index_mapping, only: &
        Nx2rx_hm_idx, & ! Procedure(s)
        mvr_hm_max

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use stats_type_utilities, only: &
        stat_begin_update, & ! Subroutines
        stat_end_update

    use stats_variables, only: &
        stats_zt, &  ! Variables
        l_stats_samp

    implicit none

    intrinsic :: trim

    ! Input Variables
    integer, intent(in) :: hydromet_dim, nz

    logical, intent(in) :: l_fill_holes_hm

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Timestep         [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zm, & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner  ! Exner function                                       [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz, hydromet_dim), intent(inout) :: &
      hydromet

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water            [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp. [K/s]

    ! Local Variables
    integer :: i, k ! Loop iterators

    real( kind = core_rknd ), dimension(nz, hydromet_dim) :: &
      hydromet_filled   ! Frozen hydrometeor mixing ratios after hole filling

    character( len = 10 ) :: hydromet_name

    real( kind = core_rknd ) :: &
      Nxm_min_coef, & ! Coefficient for min. mean value of a concentration [1/kg]
      max_velocity    ! Maximum sedimentation velocity                     [m/s]

    integer :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
               ixrm_bt, ixrm_mc

    logical :: l_hole_fill = .true.

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Start stats output for the _hf variables (changes in the hydromet array
    ! due to fill_holes_hydromet and fill_holes_vertical)
    if ( l_stats_samp ) then

       do i = 1, hydromet_dim

          ! Set up the stats indices for hydrometeor at index i
          call setup_stats_indices( i,                           & ! Intent(in)
                                    ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                    ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                    max_velocity )                 ! Intent(inout)

          call stat_begin_update( ixrm_hf, hydromet(:,i) &
                                           / dt, stats_zt )

       enddo ! i = 1, hydromet_dim

    endif ! l_stats_samp

    ! If we're dealing with negative hydrometeors, we first try to fill the
    ! holes proportionally from other same-phase hydrometeors at each height
    ! level.
    if ( any( hydromet < zero_threshold ) .and. l_fill_holes_hm ) then

       call fill_holes_hydromet( nz, hydromet_dim, hydromet, & ! Intent(in)
                                 hydromet_filled ) ! Intent(out)

       hydromet = hydromet_filled

    endif ! any( hydromet < zero ) .and. l_fill_holes_hm

    hydromet_filled = zero

    do i = 1, hydromet_dim

      ! Set up the stats indices for hydrometeor at index i
      call setup_stats_indices( i,                           & ! Intent(in)
                                ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                max_velocity )                 ! Intent(inout)

      ! Print warning message if any hydrometeor species has a value < 0.
      if ( clubb_at_least_debug_level( 1 ) ) then
         if ( any( hydromet(:,i) < zero_threshold ) ) then

            hydromet_name = hydromet_list(i)

            do k = 1, nz
               if ( hydromet(k,i) < zero_threshold ) then
                  write(fstderr,*) trim( hydromet_name ) //" < ", &
                                   zero_threshold, &
                                   " in fill_holes_driver at k= ", k
               endif ! hydromet(k,i) < 0
            enddo ! k = 1, nz
         endif ! hydromet(:,i) < 0       
      endif ! clubb_at_least_debug_level( 1 )


      ! Store the previous value of the hydrometeor for the effect of the
      ! hole-filling scheme.
!      if ( l_stats_samp ) then
!         call stat_begin_update( ixrm_hf, hydromet(:,i) &
!                                          / dt, stats_zt )
!      endif

      ! If we're dealing with a mixing ratio and hole filling is enabled,
      ! then we apply the hole filling algorithm
      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            ! Apply the hole filling algorithm
            call fill_holes_vertical( 2, zero_threshold, "zt", &
                                      rho_ds_zt, rho_ds_zm, &
                                      hydromet(:,i) )

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the
      ! hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_end_update( ixrm_hf, hydromet(:,i) &
                                        / dt, stats_zt )
      endif

      ! Store the previous value of the hydrometeor for the effect of the water
      ! vapor hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_begin_update( ixrm_wvhf, hydromet(:,i) &
                                            / dt, stats_zt )
      endif

      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            ! If the hole filling algorithm failed, then we attempt to fill
            ! the missing mass with water vapor mixing ratio.
            ! We noticed this is needed for ASEX A209, particularly if Latin
            ! hypercube sampling is enabled.  -dschanen 11 Nov 2010
            call fill_holes_wv( nz, dt, exner, hydromet_name, & ! Intent(in)
                                rvm_mc, thlm_mc, hydromet(:,i) )   ! Intent(out)

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the water vapor
      ! hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_end_update( ixrm_wvhf, hydromet(:,i) &
                                          / dt, stats_zt )
      endif

      ! Clipping for hydrometeor mixing ratios.
      if ( l_mix_rat_hm(i) ) then

         ! Store the previous value of the hydrometeor for the effect of
         ! clipping.
         if ( l_stats_samp ) then
            call stat_begin_update( ixrm_cl, &
                                    hydromet(:,i) &
                                    / dt, &
                                    stats_zt )
         endif

         if ( any( hydromet(:,i) < zero_threshold ) ) then

            ! Clip any remaining negative values of precipitating hydrometeor
            ! mixing ratios to 0.
            where ( hydromet(:,i) < zero_threshold )
               hydromet(:,i) = zero_threshold
            end where

         endif ! hydromet(:,i) < 0

         ! Eliminate very small values of mean precipitating hydrometeor mixing
         ! ratios by setting them to 0.
         do k = 2, gr%nz, 1

            if ( hydromet(k,i) <= hydromet_tol(i) ) then

               rvm_mc(k) &
               = rvm_mc(k) &
                 + ( hydromet(k,i) / dt )

               if ( .not. l_frozen_hm(i) ) then

                  ! Rain water mixing ratio
   
                  thlm_mc(k) &
                  = thlm_mc(k) &
                    - ( Lv / ( Cp * exner(k) ) ) &
                      * ( hydromet(k,i) / dt )

               else ! Frozen hydrometeor mixing ratio

                  thlm_mc(k) &
                  = thlm_mc(k) &
                    - ( Ls / ( Cp * exner(k) ) ) &
                      * ( hydromet(k,i) / dt )

               endif ! l_frozen_hm(i)

               hydromet(k,i) = zero

            endif ! hydromet(k,i) <= hydromet_tol(i)

         enddo ! k = 2, gr%nz, 1


         ! Enter the new value of the hydrometeor for the effect of clipping.
         if ( l_stats_samp ) then
            call stat_end_update( ixrm_cl, hydromet(:,i) &
                                           / dt, stats_zt )
         endif

      endif ! l_mix_rat_hm(i)

    enddo ! i = 1, hydromet_dim, 1

    ! Clipping for hydrometeor concentrations.
    do i = 1, hydromet_dim

      if ( .not. l_mix_rat_hm(i) ) then

         ! Set up the stats indices for hydrometeor at index i
         call setup_stats_indices( i,                           & ! Intent(in)
                                   ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                   ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                   max_velocity )                 ! Intent(inout)

         ! Store the previous value of the hydrometeor for the effect of
         ! clipping.
         if ( l_stats_samp ) then
            call stat_begin_update( ixrm_cl, &
                                    hydromet(:,i) &
                                    / dt, &
                                    stats_zt )
         endif

         if ( .not. l_frozen_hm(i) ) then

            ! Clipping for mean rain drop concentration, <Nr>.
            ! When mean rain water mixing ratio, <rr>, is found at a grid level,
            ! mean rain drop concentration must be at least a minimum value so
            ! that average rain drop mean volume radius stays within an upper
            ! bound.  Otherwise, mean rain drop concentration is 0.

            ! The minimum mean rain drop concentration is given by:
            !
            ! <Nr> = <rr> / ( (4/3) * pi * rho_lw * mvr_rain_max^3 ).

            Nxm_min_coef &
            = one / ( four_thirds * pi * rho_lw * mvr_hm_max(i)**3 )

         else ! l_frozen_hm(i)

            ! Clipping for mean frozen hydrometeor concentration, <Nx>.
            ! When mean frozen hydrometeor mixing ratio, <rx>, is found at a
            ! grid level, mean frozen hydrometeor concentration must be at least
            ! a minimum value so that average frozen hydrometeor mean volume
            ! radius stays within an upper bound.  Otherwise, mean frozen
            ! hydrometeor concentration is 0.

            ! The minimum mean frozen hydrometeor concentration is given by:
            !
            ! <Nx> = <rx> / ( (4/3) * pi * rho_ice * mvr_x_max^3 ).

            Nxm_min_coef &
            = one / ( four_thirds * pi * rho_ice * mvr_hm_max(i)**3 )

         endif ! .not. l_frozen_hm(i)

         ! Loop over vertical levels and increase hydrometeor concentrations
         ! when necessary.
         do k = 2, gr%nz, 1

            if ( hydromet(k,Nx2rx_hm_idx(i)) > zero ) then

               ! Hydrometeor mixing ratio, <rx>, is found at the grid level.
               hydromet(k,i) &
               = max( hydromet(k,i), &
                      Nxm_min_coef * hydromet(k,Nx2rx_hm_idx(i)) )

            else ! <rx> = 0

               hydromet(k,i) = zero

            endif ! hydromet(k,Nx2rx_hm_idx(i)) > 0

         enddo ! k = 2, gr%nz, 1

         ! Enter the new value of the hydrometeor for the effect of clipping.
         if ( l_stats_samp ) then
            call stat_end_update( ixrm_cl, hydromet(:,i) &
                                           / dt, stats_zt )
         endif

      endif ! .not. l_mix_rat_hm(i)

    enddo ! i = 1, hydromet_dim, 1


    return

  end subroutine fill_holes_driver

  !-----------------------------------------------------------------------
  subroutine setup_stats_indices( ihm,                         & ! Intent(in)
                                  ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                  ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                  max_velocity )                 ! Intent(inout)

  ! Description:
  !
  ! Determines the stats output indices depending on the hydrometeor.

  ! Attention: hydromet_list needs to be set up before this routine is called.
  !
  ! Bogus example
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------


    use array_index, only: &
        hydromet_list  ! Names of the hydrometeor species

    use stats_variables, only: &
        irrm_bt,   & ! Variable(s)
        irrm_mc,   &
        irrm_hf,   &
        irrm_wvhf, &
        irrm_cl,   &
        irim_bt,   &
        irim_mc,   &
        irim_hf,   &
        irim_wvhf, &
        irim_cl,   &
        irgm_bt,   &
        irgm_mc,   &
        irgm_hf,   &
        irgm_wvhf, &
        irgm_cl,   &
        irsm_bt,   &
        irsm_mc,   &
        irsm_hf,   &
        irsm_wvhf, &
        irsm_cl

    use stats_variables, only: &
        iNrm_bt, & ! Variable(s)
        iNrm_mc, &
        iNrm_cl, &
        iNim_bt, &
        iNim_cl, &
        iNim_mc, &
        iNsm_bt, &
        iNsm_cl, &
        iNsm_mc, &
        iNgm_bt, &
        iNgm_cl, &
        iNgm_mc, &
        iNcm_bt, &
        iNcm_cl, &
        iNcm_mc

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero

    implicit none

    ! Input Variables
    integer, intent(in) :: ihm

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      max_velocity ! Maximum sedimentation velocity [m/s]

    integer, intent(inout) :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
                              ixrm_bt, ixrm_mc

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Initializing max_velocity in order to avoid a compiler warning.
    ! Regardless of the case, it will be reset in the 'select case'
    ! statement immediately below.
    max_velocity = zero

    select case ( trim( hydromet_list(ihm) ) )
      case ( "rrm" )
        ixrm_bt   = irrm_bt
        ixrm_hf   = irrm_hf
        ixrm_wvhf = irrm_wvhf
        ixrm_cl   = irrm_cl
        ixrm_mc   = irrm_mc

        max_velocity = -9.1_core_rknd ! m/s

      case ( "rim" )
        ixrm_bt   = irim_bt
        ixrm_hf   = irim_hf
        ixrm_wvhf = irim_wvhf
        ixrm_cl   = irim_cl
        ixrm_mc   = irim_mc

        max_velocity = -1.2_core_rknd ! m/s

      case ( "rsm" )
        ixrm_bt   = irsm_bt
        ixrm_hf   = irsm_hf
        ixrm_wvhf = irsm_wvhf
        ixrm_cl   = irsm_cl
        ixrm_mc   = irsm_mc

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "rgm" )
        ixrm_bt   = irgm_bt
        ixrm_hf   = irgm_hf
        ixrm_wvhf = irgm_wvhf
        ixrm_cl   = irgm_cl
        ixrm_mc   = irgm_mc

        max_velocity = -20._core_rknd ! m/s

      case ( "Nrm" )
        ixrm_bt   = iNrm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNrm_cl
        ixrm_mc   = iNrm_mc

        max_velocity = -9.1_core_rknd ! m/s

      case ( "Nim" )
        ixrm_bt   = iNim_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNim_cl
        ixrm_mc   = iNim_mc

        max_velocity = -1.2_core_rknd ! m/s

      case ( "Nsm" )
        ixrm_bt   = iNsm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNsm_cl
        ixrm_mc   = iNsm_mc

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "Ngm" )
        ixrm_bt   = iNgm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNgm_cl
        ixrm_mc   = iNgm_mc

        max_velocity = -20._core_rknd ! m/s

      case ( "Ncm" )
        ixrm_bt   = iNcm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNcm_cl
        ixrm_mc   = iNcm_mc

        ! Use the rain water limit, since Morrison has no explicit limit on
        ! cloud water.  Presumably these numbers are never large.
        ! -dschanen 28 Sept 2009
        max_velocity = -9.1_core_rknd ! m/s

      case default
        ixrm_bt   = 0
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = 0
        ixrm_mc   = 0

        max_velocity = -9.1_core_rknd ! m/s

    end select


    return

  end subroutine setup_stats_indices
  !-----------------------------------------------------------------------

end module fill_holes
