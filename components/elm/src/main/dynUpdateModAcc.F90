module dynUpdateModAcc

   use shr_kind_mod , only : r8 => shr_kind_r8
   use decompMod, only : bounds_type
   use dynColumnStateUpdaterMod, only : column_state_updater_type
   use dynPatchStateUpdaterMod , only : patch_state_updater_type 
   use GridcellType    , only : grc_pp
   use LandunitType    , only : lun_pp
   use VegetationType       , only : veg_pp
   use ColumnType           , only : col_pp
   use elm_varcon      , only : ispval
   use dynColumnTemplateMod, only : TEMPLATE_NONE_FOUND
   
   implicit none 
   public
   public :: update_column_state_acc
   public :: update_column_state_no_special_handling_acc 
   public :: old_weight_was_zeroAcc 
   public :: patch_grewAcc
   public :: update_patch_stateAcc 
   public :: update_patch_state_partition_flux_by_typeAcc
   public :: patch_set_new_weightsAcc
   public :: column_set_new_weightsAcc
   public :: set_old_patch_weightsAcc 
   public :: set_old_column_weightsAcc

contains

  !-----------------------------------------------------------------------
  subroutine update_column_state_acc(column_state_updater, bounds, &
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

          area_lost = -1._r8 * column_state_updater%area_gained_col(c)
          total_area_lost_grc(g) = total_area_lost_grc(g) + area_lost
          area_weighted_loss = area_lost * var(c) 
          total_loss_grc(g) = total_loss_grc(g) + area_weighted_loss
         
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
         val_old = var(c)
         !
         ! Need to make sure fractional_area_new /= 0 to avoid divide-by-zero. Note
         ! that fractional_area_new == 0 can only happen if both
         ! fractional_area_old(c) == 0 and the fractional_areas of the shrinking
         ! columns were all 0 - in which case the value of var is irrelevant for
         ! conservation purposes.
         var(c) = (column_state_updater%cwtgcell_old(c) *    var(c) + mass_gained) / &
              (column_state_updater%cwtgcell_new(c))

         adjustment(c) = var(c) - val_old 
       end if
    end do

  end subroutine update_column_state_acc

  !-----------------------------------------------------------------------
  subroutine update_column_state_no_special_handling_acc(column_state_updater, bounds, clump_index, &
       var, adjustment)
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

    ! Apparent state adjustment in each column
    real(r8),  intent(out) :: adjustment( bounds%begc: )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: non_conserved_mass!(bounds%begg:bounds%endg)
    real(r8), parameter :: conservation_tolerance = 1.e-12_r8
    integer  :: c,g

    ! fractional area lost from this column
    real(r8) :: area_lost

    ! area-weighted amount lost from a given column
    real(r8) :: area_weighted_loss

    ! area-weighted amount lost from decreasing weights, in each grid cell
    ! ((mass-per-unit-area) * (fractional area lost))
    real(r8) :: total_loss_grc!(bounds%begg:bounds%endg)

    ! total area lost by columns decreasing in area (fractional area of grid cell)
    ! Note that this should exactly equal the total area gained by columns increasing in area
    real(r8) :: total_area_lost_grc!(bounds%begg:bounds%endg)

    ! amount of state gain needed per unit area of area gained
    real(r8) :: gain_per_unit_area_grc!(bounds%begg:bounds%endg)

    ! mass gained by a given column ((mass-gained-per-unit-area) * (fractional area gained))
    real(r8) :: mass_gained

    ! value in a column before update
    real(r8) :: val_old

    !character(len=*), parameter :: subname = 'update_column_state_no_special_handling'
    !-----------------------------------------------------------------------

   ! Even if there's no work to be done, need to zero out adjustment, since it's
   ! intent(out), and caller may expect it to return in a reasonable state.
     
   ! adjustment(bounds%begc:bounds%endc) = 0._r8


      ! non_conserved_mass(bounds%begg:bounds%endg) = 0._r8

      ! call update_column_state_acc(column_state_updater, bounds, &
      !          var, non_conserved_mass, adjustment)
      ! ------------------------------------------------------------------------
      ! Begin main work
      ! ------------------------------------------------------------------------

      ! Determine the total mass loss for each grid cell, along with the gross area loss
      ! (which should match the gross area gain)
      total_loss_grc = 0._r8    !(bounds%begg:bounds%endg) = 0._r8
      total_area_lost_grc = 0._r8 !(bounds%begg:bounds%endg) = 0._r8
      do c = bounds%begc, bounds%endc
         g = col_pp%gridcell(c)
         if (column_state_updater%area_gained_col(c) < 0._r8) then
  
            area_lost = -1._r8 * column_state_updater%area_gained_col(c)
            total_area_lost_grc = total_area_lost_grc + area_lost
            area_weighted_loss = area_lost * var(c) 
            total_loss_grc = total_loss_grc + area_weighted_loss
           
         end if
      end do
      ! Determine the mass loss per unit area for each grid cell. We essentially lump all of
      ! the loss together in a "loss" pool in each grid cell, so that we can then
      ! distribute that loss amongst the growing columns.
      if (total_area_lost_grc > 0._r8) then
         gain_per_unit_area_grc = total_loss_grc / total_area_lost_grc
      else
         gain_per_unit_area_grc = 0._r8
      end if

      ! Distribute gain to growing columns
      do c = bounds%begc, bounds%endc
         g = col_pp%gridcell(c)
         adjustment(c) = 0._r8
         if (column_state_updater%area_gained_col(c) > 0._r8) then
            mass_gained = column_state_updater%area_gained_col(c) * gain_per_unit_area_grc
            val_old = var(c)
            !
            ! Need to make sure fractional_area_new /= 0 to avoid divide-by-zero. Note
            ! that fractional_area_new == 0 can only happen if both
            ! fractional_area_old(c) == 0 and the fractional_areas of the shrinking
            ! columns were all 0 - in which case the value of var is irrelevant for
            ! conservation purposes.
            var(c) = (column_state_updater%cwtgcell_old(c) * var(c) + mass_gained) / &
               (column_state_updater%cwtgcell_new(c))

            adjustment(c) = var(c) - val_old 
         end if
       end do
      ! Since there is no special handling in this routine, the non_conserved_mass variable
      ! should not have any accumulation. We allow for roundoff-level accumulation in case
      ! non-conserved mass is determined in a way that is prone to roundoff-level errors.
      !err_msg = subname//': ERROR: failure to conserve mass when using no special handling'
      non_conserved_mass = abs(non_conserved_mass)
      if(non_conserved_mass > conservation_tolerance) print *, "Error Mass not conserved "
      ! if(sum(non_conserved_mass(bounds%begg:bounds%endg)) < conservation_tolerance) Then
      !   print *,  'ERROR: failure to conserve mass when using no special handling'
      ! end if

  end subroutine update_column_state_no_special_handling_acc

  !-----------------------------------------------------------------------
  function old_weight_was_zeroAcc(patch_state_updater, bounds) result(old_weight_was_zero)
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
  function patch_grewAcc(patch_state_updater, bounds) result(patch_grew)
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
  subroutine update_patch_stateAcc(patch_state_updater,p,c, &
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
    integer , value  , intent(in) :: p,c
    real(r8), intent(inout) :: var   ! patch-level state variable

    ! Accumulated flux from shrinking areas, expressed as mass per unit area COLUMN, using
    ! the OLD column weight. (The use of the old column weight is appropriate if the
    ! fluxes from shrinking patches are applied to column-level state variables BEFORE
    ! doing column-level state adjustments via the column_state_updater). For shrinking
    ! areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_col_area

    ! Accumulated flux from shrinking areas, expressed as mass per unit area GRIDCELL. For
    ! shrinking areas, this is given as a NEGATIVE quantity. Often you will provide one of
    ! flux_out_col_area or flux_out_grc_area, but it is okay to provide both, or to
    ! provide neither if you don't need to track the flux out from this state variable.
    real(r8), intent(inout), optional :: flux_out_grc_area

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition
    !
    !
    ! !LOCAL VARIABLES:

    !-----------------------------------------------------------------------

   !  if (present(seed_addition)) then
   !     if (.not. present(seed)) then
   !        print *, ' ERROR: seed_addition can only be provided if seed is provided'
   !     end if
   !  end if

      if (patch_state_updater%dwt(p) > 0._r8) then
          var  = var  * patch_state_updater%growing_old_fraction(p)
      if (present(seed)) then
         var  = var  + seed  * patch_state_updater%growing_new_fraction(p)
         if (present(seed_addition)) then
            seed_addition  = seed_addition  + seed  * patch_state_updater%dwt(p)
         end if
       end if
       !
    else if (patch_state_updater%dwt(p) < 0._r8) then
      if (present(flux_out_grc_area)) then
         flux_out_grc_area  = flux_out_grc_area  + var  * patch_state_updater%dwt(p)
      end if
      if (present(flux_out_col_area)) then
         ! No need to check for divide by 0 here: If dwt < 0 then we must have had
         ! cwtgcell_old > 0.
         flux_out_col_area  = flux_out_col_area  + &
              var  * (patch_state_updater%dwt(p) / patch_state_updater%cwtgcell_old(c))
      end if
    end if

  end subroutine update_patch_stateAcc

  !-----------------------------------------------------------------------
  subroutine update_patch_state_partition_flux_by_typeAcc(patch_state_updater, &
       p,c,flux1_out_dest, flux1_fraction_by_pft_type, &
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
    integer   , value, intent(in) :: p,c
    character(len=1), intent(in) :: flux1_out_dest  ! flux1_out to column or grid (default)
    real(r8), intent(in) :: flux1_fraction_by_pft_type( 0: ) ! fraction of flux that goes into flux1_out, indexed by pft type
    real(r8), intent(inout) :: var!( bounds%begp: ) ! patch-level state variable

    ! Accumulated fluxes from shrinking areas. For shrinking areas, these are given as
    ! NEGATIVE quantities. Even though these are patch-level arrays, they are expressed
    ! as mass per unit area GRIDCELL (so these are equivalent to the flux_out_grc_area
    ! argument in the main update_patch_state routine).
    real(r8), intent(inout) :: flux1_out!( bounds%begp: )
    real(r8), intent(inout) :: flux2_out!( bounds%begp: )

    ! If provided, this gives some 'seed' amount added to the state in the area into
    ! which each growing patch grows. The value is ignored for patches that are either
    ! constant or shrinking in area.
    real(r8), intent(in), optional :: seed!( bounds%begp: )

    ! If provided, this accumulates the amount of seed added to each patch. This gives
    ! seed(p) * dwt(p). This can only be provided if seed is provided. Even though this is
    ! a patch-level array, it is expressed as mass per unit area GRIDCELL.
    real(r8), intent(inout), optional :: seed_addition!( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    real(r8) :: total_flux_out!(bounds%begp:bounds%endp)
    real(r8) :: my_flux1_fraction
    !-----------------------------------------------------------------------


    total_flux_out  = 0._r8
    if (flux1_out_dest=='c') then
       call update_patch_stateAcc(patch_state_updater, p, c, &
          var, flux_out_col_area = total_flux_out, &
          seed = seed, seed_addition = seed_addition)
    else
       call update_patch_stateAcc(patch_state_updater,p,c, &
          var, flux_out_grc_area = total_flux_out, &
          seed = seed, seed_addition = seed_addition)
    end if

    my_flux1_fraction = flux1_fraction_by_pft_type(veg_pp%itype(p))
    flux1_out = flux1_out + total_flux_out * my_flux1_fraction
    flux2_out = flux2_out + total_flux_out * (1._r8 - my_flux1_fraction)

  end subroutine update_patch_state_partition_flux_by_typeAcc
  !-----------------------------------------------------------------------
  subroutine patch_set_new_weightsAcc(patch_state_updater, bounds)
    !
    ! !DESCRIPTION:
    ! Set subgrid weights after dyn subgrid updates
    !
    ! !USES:
    ! !ARGUMENTS:
    type(patch_state_updater_type), intent(inout) :: patch_state_updater
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p

    !-----------------------------------------------------------------------
    !$acc parallel loop independent gang vector default(present) 
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

  !-----------------------------------------------------------------------
  subroutine set_old_patch_weightsAcc(this, bounds)
   !
   ! !DESCRIPTION:
   ! Set subgrid weights before dyn subgrid updates
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   !$acc routine seq
   type(patch_state_updater_type), intent(inout) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: p
   integer :: c

   ! character(len=*), parameter :: subname = 'set_old_weights'
   !-----------------------------------------------------------------------

   do p = bounds%begp, bounds%endp
      c = veg_pp%column(p)
      this%pwtgcell_old(p) = veg_pp%wtgcell(p)
      this%cwtgcell_old(c) = col_pp%wtgcell(c)
   end do

 end subroutine set_old_patch_weightsAcc

 !-----------------------------------------------------------------------
 subroutine set_old_column_weightsAcc(this, bounds)
   !
   ! !DESCRIPTION:
   ! Set subgrid weights before dyn subgrid updates
   !
   ! !USES:
   ! !ARGUMENTS:
   !$acc routine seq 
   use landunit_varcon, only : istsoil

   type(column_state_updater_type) , intent(inout) :: this
   type(bounds_type)                , intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   logical :: found  ! whether a suitable template column has been found
   integer :: g,l,c, c_l  ! indices of grid cell, landunit, column
   ! character(len=*), parameter :: subname = 'set_old_weights'
   !-----------------------------------------------------------------------
   do c = bounds%begc, bounds%endc
     this%cwtgcell_old(c) = col_pp%wtgcell(c)
   end do

   do c = bounds%begc, bounds%endc
      found = .false.
      g = col_pp%gridcell(c)
      l = grc_pp%landunit_indices(istsoil, g)

      ! If this landunit exists on this grid cell...
      if (l /= ispval) then
         ! Loop through columns on this landunit; stop if as soon as we find an active
         ! column: that will serve as the template
         c_l = lun_pp%coli(l)
         do while (.not. found .and. c_l <= lun_pp%colf(l))
            if (col_pp%active(c_l)) then
               found = .true.
            else
               c_l = c_l + 1
            end if
         end do
      end if

      if (found) then
         this%natveg_template_col(c) = c_l
      else
         this%natveg_template_col(c) = TEMPLATE_NONE_FOUND
      end if
      ! call template_col_from_natveg_array(bounds, col_pp%(bounds%begc:bounds%endc), &
      !        this%natveg_template_col(bounds%begc:bounds%endc))
   end do 

 end subroutine set_old_column_weightsAcc


end module dynUpdateModAcc
