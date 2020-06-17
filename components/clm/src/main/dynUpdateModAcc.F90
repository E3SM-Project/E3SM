module dynUpdateModAcc

  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod, only : bounds_type
  use dynColumnStateUpdaterMod, only : column_state_updater_type
  use dynPatchStateUpdaterMod , only : patch_state_updater_type
  use VegetationType       , only : veg_pp
  use ColumnType           , only : col_pp

  public

contains
  !-----------------------------------------------------------------------
  subroutine update_column_state_acc(column_state_updater, bounds, &
       vals_input, vals_input_valid, has_prognostic_state, &
       fractional_area_old, fractional_area_new, &
       var, non_conserved_mass, adjustment)
    !
    ! !DESCRIPTION:
    ! Do the work of updating a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(column_state_updater_type), intent(in) :: column_state_updater
    type(bounds_type), intent(in) :: bounds

    ! value used as input for each column, if that column is shrinking (can differ from
    ! var if we're doing special handling of that column)
    real(r8), intent(in) :: vals_input( bounds%begc: )

    ! whether each item in vals_input_valid is valid. An entry can be invalid if there is
    ! was no way to derive an input for that column.
    logical, intent(in) :: vals_input_valid( bounds%begc: )

    ! whether each column simulates the given variable (which, among other things,
    ! determines whether it can accept mass of this variable)
    logical, intent(in) :: has_prognostic_state( bounds%begc: )

    ! Fraction of each column over which the state variable applies, for both the old and
    ! new subgrid weights
    real(r8), intent(in) :: fractional_area_old( bounds%begc: )
    real(r8), intent(in) :: fractional_area_new( bounds%begc: )

    ! column-level variable of interest, updated in-place
    real(r8), intent(inout) :: var( bounds%begc: )

    ! mass lost (per unit of grid cell area) from each grid cell, if doing special
    ! handling that leads mass to not be conserved; this can happen due to growing columns
    ! where has_prognostic_state is false, or due to shrinking columns where
    ! has_prognostic_state is false but vals_input /= 0. Positive denotes mass lost from
    ! the grid cell, negative denotes mass gained by the grid cell.
    real(r8), intent(inout) :: non_conserved_mass( bounds%begg: )

    ! Apparent state adjustment in each column
    real(r8), optional, intent(inout) :: adjustment( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    integer  :: c, g

    ! whether vals_input /= var in the columns where it should be equal
    logical :: bad_vals_input(bounds%begc:bounds%endc)

    ! fractional area lost from this column
    real(r8) :: area_lost

    ! area-weighted amount lost from a given column
    real(r8) :: area_weighted_loss

    ! area-weighted amount lost from decreasing weights, in each grid cell
    ! ((mass-per-unit-area) * (fractional area lost))
    real(r8) :: total_loss_grc(bounds%begg:bounds%endg)

    ! total area lost by columns decreasing in area (fractional area of grid cell)
    ! Note that this should exactly equal the total area gained by columns increasing in area
    real(r8) :: total_area_lost_grc(bounds%begg:bounds%endg)

    ! amount of state gain needed per unit area of area gained
    real(r8) :: gain_per_unit_area_grc(bounds%begg:bounds%endg)

    ! mass gained by a given column ((mass-gained-per-unit-area) * (fractional area gained))
    real(r8) :: mass_gained

    ! value in a column before update
    real(r8) :: val_old

    !character(len=:), allocatable :: message

    !character(len=*), parameter :: subname = 'update_column_state'
    !-----------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ! Error-checking on inputs
    ! ------------------------------------------------------------------------
    !!SHR_ASSERT_ALL((ubound(var) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!SHR_ASSERT_ALL((ubound(vals_input) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!SHR_ASSERT_ALL((ubound(has_prognostic_state) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!SHR_ASSERT_ALL((ubound(fractional_area_old) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!SHR_ASSERT_ALL((ubound(fractional_area_new) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!SHR_ASSERT_ALL((ubound(non_conserved_mass) == (/bounds%endg/)), errMsg(sourcefile, __LINE__))
    !!if (present(adjustment)) then
    !!   SHR_ASSERT_ALL((ubound(adjustment) == (/bounds%endc/)), errMsg(sourcefile, __LINE__))
    !!end if

    ! For the sake of conservation - including the calculation of non_conserved_mass - we
    ! assume that vals_input == var wherever has_prognostic_state is .true. We ensure that
    ! is the case here.
    where(has_prognostic_state .and. vals_input_valid)
       bad_vals_input = (vals_input /= var)
    elsewhere
       ! where has_prognostic_state is false, vals_input can be anything
       bad_vals_input = .false.
    end where

    ! ------------------------------------------------------------------------
    ! Begin main work
    ! ------------------------------------------------------------------------

    ! Determine the total mass loss for each grid cell, along with the gross area loss
    ! (which should match the gross area gain)
    total_loss_grc(bounds%begg:bounds%endg) = 0._r8
    total_area_lost_grc(bounds%begg:bounds%endg) = 0._r8
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       if (column_state_updater%area_gained_col(c) < 0._r8) then
          if (.not. vals_input_valid(c)) then
             print *, ' ERROR: shrinking column without valid input value'
          end if
          area_lost = -1._r8 * column_state_updater%area_gained_col(c)
          total_area_lost_grc(g) = total_area_lost_grc(g) + area_lost
          area_weighted_loss = area_lost * vals_input(c) * fractional_area_old(c)
          total_loss_grc(g) = total_loss_grc(g) + area_weighted_loss
          if (.not. has_prognostic_state(c)) then
             ! If a column doesn't model this state variable, then its vals_input value is
             ! really some fictitious quantity. So we track how much of this fictitious
             ! quantity we added to the system.
             non_conserved_mass(g) = non_conserved_mass(g) - area_weighted_loss
          end if
       end if
    end do

    ! Determine the mass loss per unit area for each grid cell. We essentially lump all of
    ! the loss together in a "loss" pool in each grid cell, so that we can then
    ! distribute that loss amongst the growing columns.
    do g = bounds%begg, bounds%endg
       if (total_area_lost_grc(g) > 0._r8) then
          gain_per_unit_area_grc(g) = total_loss_grc(g) / total_area_lost_grc(g)
       else
          gain_per_unit_area_grc(g) = 0._r8
       end if
    end do

    ! Distribute gain to growing columns
    do c = bounds%begc, bounds%endc
       g = col_pp%gridcell(c)
       if (column_state_updater%area_gained_col(c) > 0._r8) then
          mass_gained = column_state_updater%area_gained_col(c) * gain_per_unit_area_grc(g)
          if (has_prognostic_state(c)) then
             val_old = var(c)

             ! Need to make sure fractional_area_new /= 0 to avoid divide-by-zero. Note
             ! that fractional_area_new == 0 can only happen if both
             ! fractional_area_old(c) == 0 and the fractional_areas of the shrinking
             ! columns were all 0 - in which case the value of var is irrelevant for
             ! conservation purposes.
             if (fractional_area_new(c) /= 0._r8) then
                var(c) = (column_state_updater%cwtgcell_old(c) * var(c) * fractional_area_old(c) + mass_gained) / &
                     (column_state_updater%cwtgcell_new(c) * fractional_area_new(c))
             end if

             if (present(adjustment)) then
                adjustment(c) = var(c) * fractional_area_new(c) - &
                     val_old * fractional_area_old(c)
             end if
          else
             non_conserved_mass(g) = non_conserved_mass(g) + mass_gained
          end if
       end if
    end do

  end subroutine update_column_state_acc

  !-----------------------------------------------------------------------
  subroutine update_column_state_no_special_handling_acc(column_state_updater, bounds, clump_index, &
       var, fractional_area_old, fractional_area_new, adjustment)
    !
    ! !DESCRIPTION:
    ! Adjust the values of a column-level state variable due to changes in subgrid
    ! weights.
    !
    ! This method does no special handling of any columns.
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(column_state_updater_type), intent(in) :: column_state_updater
    type(bounds_type), intent(in) :: bounds

    ! Index of clump on which we're currently operating. Note that this implies that this
    ! routine must be called from within a clump loop.
    integer, intent(in) :: clump_index

    real(r8), intent(inout) :: var( bounds%begc: ) ! column-level variable

    ! Fraction of each column over which the state variable applies. See module-level
    ! documentation for details. You must provide both old & new fractional areas, or
    ! neither: it is invalid to provide just one. Fractional areas should be valid for all
    ! columns, and fractional_area_new should have been computed based on a call to
    ! update_column_state_no_special_handling: code that works with these fractional areas
    ! is not able to do any special handling.
    real(r8), optional, intent(in) :: fractional_area_old( bounds%begc: )
    real(r8), optional, intent(in) :: fractional_area_new( bounds%begc: )

    ! Apparent state adjustment in each column
    real(r8), optional, intent(out) :: adjustment( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: vals_input(bounds%begc:bounds%endc)
    logical  :: vals_input_valid(bounds%begc:bounds%endc)
    logical  :: has_prognostic_state(bounds%begc:bounds%endc)
    real(r8) :: non_conserved_mass(bounds%begg:bounds%endg)
    real(r8), parameter :: conservation_tolerance = 1.e-12_r8
    integer  :: g

    !character(len=*), parameter :: subname = 'update_column_state_no_special_handling'
    !-----------------------------------------------------------------------


    ! Even if there's no work to be done, need to zero out adjustment, since it's
    ! intent(out), and caller may expect it to return in a reasonable state.
    if (present(adjustment)) then
       adjustment(bounds%begc:bounds%endc) = 0._r8
    end if

    if (column_state_updater%any_changes(clump_index)) then

       vals_input(bounds%begc:bounds%endc) = var(bounds%begc:bounds%endc)
       vals_input_valid(bounds%begc:bounds%endc) = .true.
       has_prognostic_state(bounds%begc:bounds%endc) = .true.
       non_conserved_mass(bounds%begg:bounds%endg) = 0._r8

       ! explicit bounds not needed on any of these arguments - and specifying explicit
       ! bounds defeats some later bounds checking (for fractional_area_old and
       ! fractional_area_new)
       call update_column_state_with_optional_fractions_acc( column_state_updater,&
            bounds = bounds, &
            vals_input = vals_input, &
            vals_input_valid = vals_input_valid, &
            has_prognostic_state = has_prognostic_state, &
            var = var, &
            non_conserved_mass = non_conserved_mass, &
            fractional_area_old = fractional_area_old, &
            fractional_area_new = fractional_area_new, &
            adjustment = adjustment)

       ! Since there is no special handling in this routine, the non_conserved_mass variable
       ! should not have any accumulation. We allow for roundoff-level accumulation in case
       ! non-conserved mass is determined in a way that is prone to roundoff-level errors.
       !err_msg = subname//': ERROR: failure to conserve mass when using no special handling'
       do g = bounds%begg, bounds%endg
         non_conserved_mass(g) = abs(non_conserved_mass(g))
       end do
       if(sum(non_conserved_mass(bounds%begg:bounds%endg)) < conservation_tolerance) Then
         print *,  'ERROR: failure to conserve mass when using no special handling'
       end if
    end if

  end subroutine update_column_state_no_special_handling_acc


  !-----------------------------------------------------------------------
  subroutine update_column_state_with_optional_fractions_acc(column_state_updater, bounds, &
       vals_input, vals_input_valid, has_prognostic_state, &
       var, non_conserved_mass, &
       fractional_area_old, fractional_area_new, &
       adjustment)
    !
    ! !DESCRIPTION:
    ! Intermediate routine between the public routines and the real work routine
    ! (update_column_state). This routine determines the fractional areas to use in the
    ! call to update_column_state, and then does the call to update_column_state.
    !$acc routine seq
    ! !USES:
    !
    ! !ARGUMENTS:
    type(column_state_updater_type), intent(in) :: column_state_updater
    type(bounds_type), intent(in) :: bounds

    ! value used as input for each column, if that column is shrinking (can differ from
    ! var if we're doing special handling of that column)
    real(r8), intent(in) :: vals_input( bounds%begc: )

    ! whether each item in vals_input_valid is valid. An entry can be invalid if there is
    ! was no way to derive an input for that column.
    logical, intent(in) :: vals_input_valid( bounds%begc: )

    ! whether each column simulates the given variable (which, among other things,
    ! determines whether it can accept mass of this variable)
    logical, intent(in) :: has_prognostic_state( bounds%begc: )

    ! column-level variable of interest, updated in-place
    real(r8), intent(inout) :: var( bounds%begc: )

    ! mass lost (per unit of grid cell area) from each grid cell, if doing special
    ! handling that leads mass to not be conserved; this can happen due to growing columns
    ! where has_prognostic_state is false, or due to shrinking columns where
    ! has_prognostic_state is false but vals_input /= 0. Positive denotes mass lost from
    ! the grid cell, negative denotes mass gained by the grid cell.
    real(r8), intent(inout) :: non_conserved_mass( bounds%begg: )

    ! Fraction of each column over which the state variable applies. See module-level
    ! documentation for details. You must provide both old & new fractional areas, or
    ! neither: it is invalid to provide just one. Fractional areas should be valid for all
    ! columns, and fractional_area_new should have been computed based on a call to
    ! update_column_state_no_special_handling: code that works with these fractional areas
    ! is not able to do any special handling.
    real(r8), optional, intent(in) :: fractional_area_old( bounds%begc: )
    real(r8), optional, intent(in) :: fractional_area_new( bounds%begc: )

    ! Apparent state adjustment in each column
    real(r8), optional, intent(inout) :: adjustment( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: my_fractional_area_old(bounds%begc:bounds%endc)
    real(r8) :: my_fractional_area_new(bounds%begc:bounds%endc)

    !-----------------------------------------------------------------------

    if (present(fractional_area_old) .and. .not. present(fractional_area_new)) then
       print *,' ERROR: If fractional_area_old is provided, then fractional_area_new must be provided, too'
    end if

    if (present(fractional_area_new) .and. .not. present(fractional_area_old)) then
       print *,' ERROR: If fractional_area_new is provided, then fractional_area_old must be provided, too'
    end if

    if (present(fractional_area_old)) then
       my_fractional_area_old(bounds%begc:bounds%endc) = fractional_area_old(bounds%begc:bounds%endc)
    else
       my_fractional_area_old(bounds%begc:bounds%endc) = 1._r8
    end if

    if (present(fractional_area_new)) then
       my_fractional_area_new(bounds%begc:bounds%endc) = fractional_area_new(bounds%begc:bounds%endc)
    else
       my_fractional_area_new(bounds%begc:bounds%endc) = 1._r8
    end if

    call update_column_state_acc(column_state_updater,&
         bounds = bounds, &
         vals_input = vals_input(bounds%begc:bounds%endc), &
         vals_input_valid     = vals_input_valid(bounds%begc:bounds%endc), &
         has_prognostic_state = has_prognostic_state(bounds%begc:bounds%endc), &
         fractional_area_old  = my_fractional_area_old(bounds%begc:bounds%endc), &
         fractional_area_new  = my_fractional_area_new(bounds%begc:bounds%endc), &
         var = var(bounds%begc:bounds%endc), &
         non_conserved_mass = non_conserved_mass(bounds%begc:bounds%endc), &
         adjustment = adjustment)

  end subroutine update_column_state_with_optional_fractions_acc

  !-----------------------------------------------------------------------
  function old_weight_was_zeroAcc(patch_state_updater, bounds)
    !
    ! !DESCRIPTION:
    ! Returns a patch-level logical array that is true wherever the patch weight was zero
    ! prior to weight updates
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(in) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    logical :: old_weight_was_zero(bounds%begp:bounds%endp)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: p
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       old_weight_was_zero(p) = (patch_state_updater%pwtgcell_old(p) == 0._r8)
    end do

  end function old_weight_was_zeroAcc

  !-----------------------------------------------------------------------
  function patch_grewAcc(patch_state_updater, bounds)
    !
    ! !DESCRIPTION:
    ! Returns a patch-level logical array that is true wherever the patch grew in this
    ! time step
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(in) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    logical :: patch_grew(bounds%begp:bounds%endp)  ! function result
    !
    ! !LOCAL VARIABLES:
    integer :: p
    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       patch_grew(p) = (patch_state_updater%dwt(p) > 0._r8)
    end do

  end function patch_grewAcc

  !-----------------------------------------------------------------------
  subroutine update_patch_stateAcc(patch_state_updater, bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       var, flux_out_col_area, flux_out_grc_area, &
       seed, seed_addition)
    !
    ! !DESCRIPTION:
    ! Update a patch-level state variable and compute associated fluxes based on changing
    ! patch areas
    !
    ! For growing patches, this subroutine adjusts var. For shrinking patches, this
    ! subroutine accumulates flux in flux_out_col_area.
    !
    ! Changes are only made within the given filter. Note that this filter should include
    ! inactive as well as active patches, so that it includes patches that just became
    ! inactive.
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(in) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer, intent(in) :: filterp_with_inactive(:) ! patch filter that includes inactive points
    real(r8), intent(inout) :: var( bounds%begp: ) ! patch-level state variable

    ! Accumulated flux from shrinking areas, expressed as mass per unit area COLUMN, using
    ! the OLD column weight. (The use of the old column weight is appropriate if the
    ! fluxes from shrinking patches are applied to column-level state variables BEFORE
    ! doing column-level state adjustments via the column_state_updater). For shrinking
    ! areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_col_area( bounds%begp: )

    ! Accumulated flux from shrinking areas, expressed as mass per unit area GRIDCELL. For
    ! shrinking areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_grc_area( bounds%begp: )

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed( bounds%begp: )

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition( bounds%begp: )
    !
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p, c

    !-----------------------------------------------------------------------

    if (present(seed_addition)) then
       if (.not. present(seed)) then
          print *, ' ERROR: seed_addition can only be provided if seed is provided'
       end if
    end if

    do fp = 1, num_filterp_with_inactive
       p = filterp_with_inactive(fp)
       c = veg_pp%column(p)

       if (patch_state_updater%dwt(p) > 0._r8) then
          var(p) = var(p) * patch_state_updater%growing_old_fraction(p)
          if (present(seed)) then
             var(p) = var(p) + seed(p) * patch_state_updater%growing_new_fraction(p)
             if (present(seed_addition)) then
                seed_addition(p) = seed_addition(p) + seed(p) * patch_state_updater%dwt(p)
             end if
          end if

       else if (patch_state_updater%dwt(p) < 0._r8) then
          if (present(flux_out_grc_area)) then
             flux_out_grc_area(p) = flux_out_grc_area(p) + var(p) * patch_state_updater%dwt(p)
          end if
          if (present(flux_out_col_area)) then
             ! No need to check for divide by 0 here: If dwt < 0 then we must have had
             ! cwtgcell_old > 0.
             flux_out_col_area(p) = flux_out_col_area(p) + &
                  var(p) * (patch_state_updater%dwt(p) / patch_state_updater%cwtgcell_old(c))
          end if
       end if
    end do

  end subroutine update_patch_stateAcc

  !-----------------------------------------------------------------------
  subroutine update_patch_state_partition_flux_by_typeAcc(patch_state_updater, bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       flux1_fraction_by_pft_type, &
       var, flux1_out, flux2_out, &
       seed, seed_addition)
    !
    ! !DESCRIPTION:
    ! Update a patch-level state variable and compute associated fluxes based on
    ! changing patch areas, with flux out partitioned into two fluxes. Partitioning is
    ! based on pft type.
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(in) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_filterp_with_inactive ! number of points in filterp_with_inactive
    integer, intent(in) :: filterp_with_inactive(:) ! patch filter that includes inactive points
    real(r8), intent(in) :: flux1_fraction_by_pft_type( 0: ) ! fraction of flux that goes into flux1_out, indexed by pft type
    real(r8), intent(inout) :: var( bounds%begp: ) ! patch-level state variable

    ! Accumulated fluxes from shrinking areas. For shrinking areas, these are given as
    ! NEGATIVE quantities. Even though these are patch-level arrays, they are expressed
    ! as mass per unit area GRIDCELL (so these are equivalent to the flux_out_grc_area
    ! argument in the main update_patch_state routine).
    real(r8), intent(inout) :: flux1_out( bounds%begp: )
    real(r8), intent(inout) :: flux2_out( bounds%begp: )

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed( bounds%begp: )

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer :: fp, p
    real(r8) :: total_flux_out(bounds%begp:bounds%endp)
    real(r8) :: my_flux1_fraction
    !-----------------------------------------------------------------------


    total_flux_out(bounds%begp:bounds%endp) = 0._r8
    call update_patch_stateAcc(patch_state_updater, bounds, &
       num_filterp_with_inactive, filterp_with_inactive, &
       var, flux_out_grc_area = total_flux_out, &
       seed = seed, seed_addition = seed_addition)

    do fp = 1, num_filterp_with_inactive
       p = filterp_with_inactive(fp)
       my_flux1_fraction = flux1_fraction_by_pft_type(veg_pp%itype(p))
       flux1_out(p) = flux1_out(p) + total_flux_out(p) * my_flux1_fraction
       flux2_out(p) = flux2_out(p) + total_flux_out(p) * (1._r8 - my_flux1_fraction)
    end do

  end subroutine update_patch_state_partition_flux_by_typeAcc
  !-----------------------------------------------------------------------
  subroutine patch_set_new_weightsAcc(patch_state_updater, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights after dyn subgrid updates
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(inout) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p

    !-----------------------------------------------------------------------

    do p = bounds%begp, bounds%endp
       patch_state_updater%pwtgcell_new(p) = veg_pp%wtgcell(p)
       patch_state_updater%dwt(p)   = patch_state_updater%pwtgcell_new(p) &
                                    - patch_state_updater%pwtgcell_old(p)
       if (patch_state_updater%dwt(p) > 0._r8) then
          patch_state_updater%growing_old_fraction(p) = &
              patch_state_updater%pwtgcell_old(p)/ patch_state_updater%pwtgcell_new(p)
          patch_state_updater%growing_new_fraction(p) = &
              patch_state_updater%dwt(p)/ patch_state_updater%pwtgcell_new(p)
       else
          ! These values are unused in this case, but set them to something reasonable for
          ! safety. (We could set them to NaN, but that requires a more expensive
          ! subroutine call, using the shr_infnan_mod infrastructure.)
          patch_state_updater%growing_old_fraction(p) = 1._r8
          patch_state_updater%growing_new_fraction(p) = 0._r8
       end if
    end do

  end subroutine patch_set_new_weightsAcc

  !-----------------------------------------------------------------------
  subroutine column_set_new_weightsAcc(column_state_updater, bounds, clump_index)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights after dyn subgrid updates
    !
    ! !USES:
    !$acc routine seq
    ! !ARGUMENTS:
    type(column_state_updater_type), intent(inout) :: column_state_updater
    type(bounds_type), intent(in) :: bounds

    ! Index of clump on which we're currently operating. Note that this implies that this
    ! routine must be called from within a clump loop.
    integer,   intent(in) :: clump_index

    !
    ! !LOCAL VARIABLES:
    integer :: c
    !-----------------------------------------------------------------------

    column_state_updater%any_changes(clump_index) = .false.

    do c = bounds%begc, bounds%endc
       column_state_updater%cwtgcell_new(c)     = col_pp%wtgcell(c)
       column_state_updater%area_gained_col(c)  = &
          column_state_updater%cwtgcell_new(c) - column_state_updater%cwtgcell_old(c)
       if (column_state_updater%area_gained_col(c) /= 0._r8) then
          column_state_updater%any_changes(clump_index) = .true.
       end if
    end do

  end subroutine column_set_new_weightsAcc

end module dynUpdateModAcc
