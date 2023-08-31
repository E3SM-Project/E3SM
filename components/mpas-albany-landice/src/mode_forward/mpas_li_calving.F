! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_calving
!
!> \brief MPAS land ice calving scheme
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This module contains several options for calving ice.
!
!-----------------------------------------------------------------------

module li_calving

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log

   use li_setup
   use li_mask
   use li_constants

   implicit none
   private

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_calve_ice, li_restore_calving_front, li_apply_front_ablation_velocity, &
             li_calculate_damage, li_finalize_damage_after_advection, li_flood_fill

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************
   contains
!***********************************************************************


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_calve_ice
!
!> \brief MPAS land ice calving scheme
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This routine contains several options for calving ice:
!> (0) Do nothing
!> (1) Calve all floating ice
!> (2) Calve ice based on a topographic threshold
!> (3) Calve ice based on an ice thickness threshold
!-----------------------------------------------------------------------

   subroutine li_calve_ice(domain, err)

      use li_advection

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool

      ! calving-relevant config options
      character (len=StrKIND), pointer :: config_calving, &
                                          config_front_mass_bal_grounded
      logical, pointer :: config_print_calving_info, config_data_calving
      logical, pointer :: config_apply_calving_mask
      real(kind=RKIND), pointer :: config_calving_timescale

      integer, pointer :: nCells
      integer, pointer :: config_number_of_blocks

      real (kind=RKIND), pointer :: deltat  !< time step (s)

      integer, dimension(:), pointer :: &
           indexToCellID       ! list of global cell IDs

      real (kind=RKIND) ::  &
           calvingFraction ! fraction of ice that calves in each column; depends on calving_timescale

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography       ! bed topography (negative below sea level)

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThickness    ! thickness of ice that calves (computed in this subroutine)
                               ! typically the entire ice thickness, but will be a fraction of the thickness
                               ! if calving_timescale > dt

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThicknessFromThreshold   ! thickness of ice that calves from threshold processes (computed in this subroutine)

      real (kind=RKIND), dimension(:), pointer :: calvingVelocity

      type (field1dReal), pointer :: originalThicknessField

      real (kind=RKIND), dimension(:), pointer :: originalThickness

      integer :: iCell

      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)
      call mpas_pool_get_config(liConfigs, 'config_apply_calving_mask', config_apply_calving_mask)
      call mpas_pool_get_config(liConfigs, 'config_front_mass_bal_grounded', config_front_mass_bal_grounded)
      call mpas_pool_get_config(liConfigs, 'config_calving_timescale', config_calving_timescale)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_data_calving', config_data_calving)
      call mpas_pool_get_config(liConfigs, 'config_number_of_blocks', config_number_of_blocks)

      ! Zero calvingThickness here instead of or in addition to in individual subroutines.
      ! This is necessary when using damage threshold calving with other calving
      ! routines. Some individual routines still set calvingThickness to zero, but
      ! this is redundant. li_apply_front_ablation_velocity will zero
      ! calvingThickness, so any routines that use that should not be applied
      ! after other calving routines.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         calvingThickness(:) = 0.0_RKIND
         calvingThicknessFromThreshold(:) = 0.0_RKIND
         calvingVelocity(:) = 0.0_RKIND
         block => block % next
      end do

      ! Get deltat from first block (same on all blocks)
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)

      ! Update mask and geometry before calling any calving routines.  May not be necessary, but best be safe.
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)
      call li_update_geometry(geometryPool)

      ! based on the calving timescale, set the fraction of ice that calves
      if (config_calving_timescale > 0.0_RKIND) then
         calvingFraction = min(deltat/config_calving_timescale, 1.0_RKIND)
         !WHL - debug
         if (config_print_calving_info) then
            call mpas_log_write('Calving a fraction of the ice thickness at each timestep')
            call mpas_log_write('deltat (yr) = $r, calvingFraction = $r', realArgs=(/deltat/scyr, calvingFraction/))
         endif
      else
         calvingFraction = 1.0_RKIND   ! calve the entire thickness in eligible columns
      endif

      if (config_print_calving_info) then
         call mpas_log_write('Do ice calving, option = ' // trim(config_calving))
         call mpas_log_write('Calving timscale (yr) = $r', realArgs=(/config_calving_timescale / scyr/))
      endif

      ! In data calving mode we need to calculate calving flux but not have it be applied.
      ! However, the eigencalving method requires multiple applications of the calvingThickness
      ! to the thickness.  So the simplest method to apply data calving is to store the old
      ! thickness and then set it back when we are done.
      if (config_data_calving) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'workCell2',  originalThicknessField)
         call mpas_allocate_scratch_field(originalThicknessField, single_block_in = .false.)
         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_array(geometryPool, 'thickness', thickness)
            call mpas_pool_get_array(scratchPool, 'workCell2', originalThickness)

            ! Store old thickness here
            originalThickness(:) = thickness(:)

            block => block % next
         end do
      endif
      
      ! compute calvingThickness based on the calving_config option
      if (trim(config_calving) == 'none') then

         ! Do nothing, but continue through rest of routine

      elseif (trim(config_calving) == 'thickness_threshold') then

         call thickness_calving(domain, calvingFraction, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'floating') then

         call floating_calving(domain, calvingFraction, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'topographic_threshold') then

         call topographic_calving(domain, calvingFraction, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'eigencalving') then

         call eigencalving(domain, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'damagecalving') then
         ! Does not currently support multiple blocks per processor.
         if (config_number_of_blocks /= 0) then
            call mpas_log_write("damagecalving requires that config_number_of_blocks=0", MPAS_LOG_ERR)
            err_tmp = 1
         endif
         err = ior(err,err_tmp)
         call damage_calving(domain, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'von_Mises_stress') then

         call von_Mises_calving(domain, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'specified_calving_velocity') then

         call specified_calving_velocity(domain, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'ismip6_retreat') then

         if (trim(config_front_mass_bal_grounded) .ne. 'none') then
            call mpas_log_write("config_front_mass_bal_grounded must be set to 'none' '// &
                                'when config_calving = 'ismip6_retreat'.", MPAS_LOG_ERR)
            err = 1
         else
            call ismip6_retreat(domain, err_tmp)
            err = ior(err, err_tmp)
         end if

      else

         call mpas_log_write("Invalid option for config_calving specified: " // trim(config_calving), MPAS_LOG_ERR)
         err = 1

      endif

      ! Consider mask calving as a possible additional step
      ! Mask calving can occur by itself or in conjunction with a physical calving law
      if (config_apply_calving_mask) then
         call mask_calving(domain, err_tmp)
         err = ior(err, err_tmp)
      endif

      ! now also remove any icebergs
      call remove_icebergs(domain)

      ! Final operations after calving has been applied.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! In data calving mode we just calculate what should be calved but don't actually calve it.
         ! So set thickness back to original value.
         if (config_data_calving) then
            call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
            call mpas_pool_get_array(scratchPool, 'workCell2', originalThickness)
            thickness(:) = originalThickness(:)
         endif

         ! Optionally, print a list of cells with calving
         if (config_print_calving_info) then
            call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
            call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

            call mpas_log_write(' ')
            call mpas_log_write('Global cell ID, bedTopography, calvingThickness:')
            do iCell = 1, nCells
               if (calvingThickness(iCell) > 0.0_RKIND) then
                  call mpas_log_write("$i $r $r", intArgs=(/indexToCellID(iCell)/), &
                     realArgs=(/bedTopography(iCell), calvingThickness(iCell)/))
               endif
            enddo
         endif   ! config_print_calving_info

         ! Update mask and geometry
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)
         call li_update_geometry(geometryPool)

         block => block % next
      end do

      if (config_data_calving) then
         call mpas_deallocate_scratch_field(originalThicknessField, single_block_in=.false.)
      endif

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_calve_ice.", MPAS_LOG_ERR)
      endif

   !--------------------------------------------------------------------
    end subroutine li_calve_ice


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine li_restore_calving_front
!
!> \brief MPAS land ice restore the calving front
!> \author William Lipscomb
!> \date   September 2015
!> \details
!>  This routine restores the calving front to its initial position:
!> (1) It removes any floating ice that has advanced beyond the initial front.
!> (2) It adds back a thin layer of ice wherever the ice has retreated from
!>     the initial front.
!-----------------------------------------------------------------------

   subroutine li_restore_calving_front(domain, err)

      use li_thermal, only: li_init_linear_temperature_in_column
      use li_advection

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: velocityPool

      integer, pointer :: nCells, nVertLevels

      logical, pointer :: &
           config_print_calving_info, &
           config_restore_calving_front_prevent_retreat

      real (kind=RKIND), pointer ::   &
           config_sea_level,          &
           config_dynamic_thickness

      integer, dimension(:), pointer :: &
           cellMask, &            ! bit mask describing whether ice is floating, dynamically active, etc.
           indexToCellID

      real(kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma    ! vertical sigma coordinate at layer midpoints

      ! geometry and calving fields
      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography,    & ! elevation of the bed
           calvingVelocity,  &
           calvingThickness, & ! thickness of ice that calves
                               ! > 0 for cells below sea level that were initially ice-free and now have ice
           restoreThickness    ! thickness of ice that is added to restore the calving front to its initial position
                               ! > 0 for cells below sea level that were initially ice-covered and now have very thin or no ice

      real (kind=RKIND) ::  &
           restoreThicknessMin  ! small thickness to which ice is restored should it fall below this thickness

      ! thermal fields
      ! These are needed to initialize the temperature profile in restored columns.
      real (kind=RKIND), dimension(:,:), pointer :: &
           temperature,           &   ! interior ice temperature
           waterFrac                  ! interior water fraction

      real (kind=RKIND), dimension(:), pointer :: &
           surfaceAirTemperature, &   ! surface air temperature
           surfaceTemperature,    &   ! surface ice temperature
           basalTemperature           ! basal ice temperature

      real (kind=RKIND), dimension(:), pointer :: surfaceSpeed

      integer :: iCell, err_tmp

      !WHL - debug
      logical, parameter :: circular_shelf = .false.
      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow
      integer :: k


      ! No longer a need to remove icebergs here, because li_calve_ice is now always called first
      ! and remove_icebergs is always called from calving routine.

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)  ! required for cellMask computation
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         ! get required fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

         ! get required fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'restoreThickness', restoreThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)

         ! get required fields from the thermal pool
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)

         call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)

         ! get config variables
         call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
         call mpas_pool_get_config(liConfigs, 'config_restore_calving_front_prevent_retreat', config_restore_calving_front_prevent_retreat)

         if (config_print_calving_info) then
            call mpas_log_write('Restore calving front')
            call mpas_log_write('max thickness (m) = $r', realArgs=(/maxval(thickness)/))

            !WHL - debug - for circular shelf test case
!            if (circular_shelf) then
!               call mpas_log_write('Initial ice thickness'
!               do iRow = nRows, 1, -1
!                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                     write(stdoutUnit,'(a3)',advance='no') '    '
!                  endif
!!!                 do i = 1, nCellsPerRow
!                  do i = nCellsPerRow/2 - 2, nCellsPerRow
!                     iCell = (iRow-1)*nCellsPerRow + i
!!!                    write(stdoutUnit,'(i5)',advance='no') iCell
!                     write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
!                  enddo
!                  write(stdoutUnit,*) ' '
!               enddo
!            endif   ! circular_shelf

         endif

         ! set restoreThicknessMin
         ! It should be less than config_dynamic_thickness so that the restored ice remains dynamically inactive,
         !  even with a certain amount of natural variability.
         ! It should also be large enough to permit stable thermal calculations.
         ! For now, setting it to 1/10 of config_dynamic_thickness

         restoreThicknessMin = 0.1_RKIND * config_dynamic_thickness

         ! calculate masks - so we know where the calving front was located initially
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! initialize
         restoreThickness = 0.0_RKIND

         if (config_restore_calving_front_prevent_retreat) then
            do iCell = 1, nCells
               if (bedTopography(iCell) < config_sea_level) then
                  ! The bed is below sea level; test for calving-front advance and retreat.
                  if (li_mask_is_initial_ice(cellMask(iCell)) .and. thickness(iCell) < restoreThicknessMin) then
                     ! Ice was present in this cell initially, but now is either very thin or absent;
                     !  reset the thickness to restoreThicknessMin.
                     ! Note: Mass is not conserved.
                     !       Save the difference (restoreThicknessMin - thickness) so as to keep track of energy non-conservation.

                     if (config_print_calving_info) then
                        call mpas_log_write('Restore ice: indexToCellID=$i, thickness=$r', intArgs=(/indexToCellID(iCell)/), realArgs=(/thickness(iCell)/))
                     endif

                     restoreThickness(iCell) = restoreThicknessMin - thickness(iCell)
                     thickness(iCell) = restoreThicknessMin

                     ! Initialize a linear temperature profile in the column
                     ! Note: Energy is not conserved.

                     call li_init_linear_temperature_in_column(&
                          nVertLevels,                   &
                          layerCenterSigma,              &
                          thickness(iCell),              &
                          surfaceAirTemperature(iCell),  &
                          temperature(:,iCell),          &
                          waterFrac(:,iCell),            &
                          surfaceTemperature(iCell),     &
                          basalTemperature(iCell))
                  endif ! retreated
               endif ! bed below sea level
            enddo ! nCells
         endif ! if preventing retreat

         ! Now check for marine regions that have advanced
         do iCell = 1, nCells
            if (bedTopography(iCell) < config_sea_level) then
               ! The bed is below sea level; test for calving-front advance and retreat.
               if (.not.li_mask_is_initial_ice(cellMask(iCell)) .and. thickness(iCell) > 0.0_RKIND) then

                  ! This cell was initially ice-free but now has ice.
                  ! Remove the ice and add it to calvingThickness.

                  if (config_print_calving_info) then
                     call mpas_log_write('Remove ice:  indexToCellID=$i, thickness=$r', intArgs=(/indexToCellID(iCell)/), realArgs=(/thickness(iCell)/))
                  endif

                  calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
                  thickness(iCell) = 0.0_RKIND
                  ! Assign flow speed to calving speed - this allows calculation of licalvf in ISMIP6 postprocessing
                  calvingVelocity(iCell) = surfaceSpeed(iCell)
               endif   ! li_mask_is_initial_ice
            endif    ! bedTopography < config_sea_level
         enddo   ! iCell

         block => block % next
      enddo

      ! Update mask and geometry
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         call li_update_geometry(geometryPool)

         block => block % next
      end do

      if (config_print_calving_info) then
         call mpas_log_write('Restored the initial calving front')

         !WHL - debug - for circular shelf test case
!         if (circular_shelf) then
!            write(stdoutUnit,*) 'Final ice thickness'
!            do iRow = nRows, 1, -1
!               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stdoutUnit,'(a3)',advance='no') '    '
!               endif
!!!               do i = 1, nCellsPerRow
!               do i = nCellsPerRow/2 - 2, nCellsPerRow
!                  iCell = (iRow-1)*nCellsPerRow + i
!!!                  write(stdoutUnit,'(i5)',advance='no') iCell
!                  write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
!               enddo
!               write(stdoutUnit,*) ' '
!            enddo
!         endif  ! circular shelf

      endif

      ! === error check
      if (err > 0) then
          call mpas_log_write("An error has occurred in li_restore_calving_front.", MPAS_LOG_ERR)
      endif


    end subroutine li_restore_calving_front

!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine thickness_calving
!
!> \brief Calve ice based on thickness threshold
!> \author William Lipscomb, moved to separate subroutine by Matt Hoffman Feb. 2018
!> \date   September 2015
!> \details  calve ice thinner than a thickness threshold
!> Note: The thickness-threshold option is different from the others.
!>       For the other options, we look at each cell and determine whether it meets the calving-law criteria
!>        (e.g., ice is floating, or the topography lies below a given level).
!>        If a cell meets the criteria and lies in the calving domain (e.g., at the margin), it is calved.
!>       For the thickness-threshold option, ice thinner than config_calving_thickness is calved,
!>        but only if it lies beyond a protected ring of thin ice at the floating margin.
!>       The reason for this more complicated approach is that we do not want to remove all floating ice
!>         thinner than the calving thickness, because then we would remove thin ice that has just
!>         been advected from active cells at the margin, and the calving front would be unable to advance.
!>       By protecting a ring of inactive ice (thickness < config_calving_thickness) at the margin,
!>        we allow ice in these cells to thicken and become active, thus advancing the calving front.
!>       The calving front retreats when active floating ice thins to become inactive, removing protection
!>        from previously protected cells.
!>
!>      Specifically, the rules are as follows:
!>      - Mark cells as active-for-calving if either (1) grounded, with thickness > config_dynamic_thickness
!>        or (2) floating, with thickness > config_calving_thickness.
!>      - Mark cells as ocean if not land and not active.
!>      - Mark cells as lying on the inactive margin if not active, but with an active neighbor.
!>      - Calve ice in ocean cells that are not on the inactive margin.

!-----------------------------------------------------------------------
   subroutine thickness_calving(domain, calvingFraction, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: calvingFraction !< fraction of possible ice to calve on this timestep

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block

      real(kind=RKIND), pointer ::  &
                                   config_calving_thickness, &
                                   config_sea_level, &
                                   config_dynamic_thickness,  &
                                   config_ice_density,  &
                                   config_ocean_density

      logical, pointer :: config_print_calving_info

      type (field1dInteger), pointer :: activeForCalvingMaskField
      integer, dimension(:), pointer :: activeForCalvingMask   ! = 1 for grounded cells thicker than config_dynamic_thickness;
                                                               ! = 1 for floating cells thicker than config_calving_thickness;
                                                               ! = 0 elsewhere

      type (field1dInteger), pointer :: inactiveMarginMaskField
      integer, dimension(:), pointer :: inactiveMarginMask
      ! = 1 for inactive cells (thin or no ice) that have 1 or more active neighbors

      type (field1dInteger), pointer :: oceanMaskField
      integer, dimension(:), pointer :: oceanMask   ! = 1 for cells that are not land and do not have active ice
                                                    ! includes floating cells with inactive icea

      type (mpas_pool_type), pointer :: meshPool, geometryPool, velocityPool, scratchPool

      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell

      integer, dimension(:,:), pointer :: &
           cellsOnCell         ! list of cells that neighbor each cell

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography       ! bed topography (negative below sea level)

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThickness    ! thickness of ice that calves (computed in this subroutine)

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThicknessFromThreshold    ! thickness of ice that calves (computed in this subroutine)

      integer, pointer :: nCells
      integer :: iCell, iCellOnCell, iCellNeighbor
      real (kind=RKIND) :: flotationThickness  ! thickness at which marine-based ice starts to float


      err = 0

      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)


      if (config_print_calving_info) then
         call mpas_log_write('Calving thickness (m) = $r', realArgs=(/config_calving_thickness/))
         call mpas_log_write('Dynamic thickness (m) = $r', realArgs=(/config_dynamic_thickness/))
      endif

      ! Make sure config_calving_thickness > config_dynamic_thickness.
      ! Otherwise the algorithm will not work.
      if (config_calving_thickness < config_dynamic_thickness) then
         call mpas_log_write('config_calving_thickness (m) = $r', realArgs=(/config_calving_thickness/), messageType=MPAS_LOG_ERR)
         call mpas_log_write('config_dynamic_thickness (m) = $r', realArgs=(/config_dynamic_thickness/), messageType=MPAS_LOG_ERR)
         call mpas_log_write('Must have config_calving_thickness > config_dynamic_thickness', messageType=MPAS_LOG_ERR)
         err = 1
      endif

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)  ! required for cellMask computation
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! get required fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

         ! get required fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)

         ! get scratch fields for calving
         ! 'true' flag means to allocate the field for a single block
         call mpas_pool_get_field(scratchPool, 'iceCellMask', activeForCalvingMaskField)
         call mpas_allocate_scratch_field(activeForCalvingMaskField, .true.)
         activeForCalvingMask => activeForCalvingMaskField % array

         call mpas_pool_get_field(scratchPool, 'iceCellMask2',  inactiveMarginMaskField)
         call mpas_allocate_scratch_field(inactiveMarginMaskField, .true.)
         inactiveMarginMask => inactiveMarginMaskField % array

         call mpas_pool_get_field(scratchPool, 'iceCellMask3', oceanMaskField)
         call mpas_allocate_scratch_field(oceanMaskField, .true.)
         oceanMask => oceanMaskField % array

         ! Initialize
         calvingThickness = 0.0_RKIND

         ! Identify cells that are active-for-calving:
         ! (1) Grounded ice with thickness > config_dynamic_thickness
         ! (2) Floating ice with thickness > config_calving_thickness

         activeForCalvingMask(:) = 0

         do iCell = 1, nCells

            if (bedTopography(iCell) >= config_sea_level) then   ! land cell
               if (thickness(iCell) > config_dynamic_thickness) then    ! active for calving
                  activeForCalvingMask(iCell) = 1
               endif
            else    ! marine cell, topography below sea level
               flotationThickness = (config_ocean_density/config_ice_density) * (config_sea_level - bedTopography(iCell))
               if (thickness(iCell) < flotationThickness) then   ! floating
                  if (thickness(iCell) > config_calving_thickness) then
                     activeForCalvingMask(iCell) = 1
                  endif
               else   ! grounded marine ice
                  if (thickness(iCell) > config_dynamic_thickness) then    ! active for calving
                     activeForCalvingMask(iCell) = 1
                  endif
               endif   ! floating or grounded
            endif   ! land or marine

         enddo   ! iCell

         ! Identify cells that are inactive but border active-for-calving cells

         inactiveMarginMask(:) = 0

         do iCell = 1, nCells
            if (activeForCalvingMask(iCell) == 0) then   ! inactive

               ! check whether any neighbor cells are active
               do iCellOnCell = 1, nEdgesOnCell(iCell)
                  iCellNeighbor = cellsOnCell(iCellOnCell,iCell)
                  if (activeForCalvingMask(iCellNeighbor) == 1) then  ! neighbor cell is active
                     inactiveMarginMask(iCell) = 1
                     exit
                  endif
               enddo   ! iCellOnCell

            endif   ! inactive
         enddo    ! iCell

         ! Identify ocean cells (not land and not active ice, but including inactive floating ice)

         where (bedTopography < config_sea_level .and. activeForCalvingMask == 0)
            oceanMask = 1
         elsewhere
            oceanMask = 0
         endwhere

!         if (config_print_calving_info) then
!
!            write(stdoutUnit,*) 'Active-for-calving mask'
!            do iRow = nRows, 1, -1
!               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stdoutUnit,'(a3)',advance='no') '    '
!               endif
!!            do i = 1, nCellsPerRow
!               do i = nCellsPerRow/2 - 2, nCellsPerRow
!                  iCell = (iRow-1)*nCellsPerRow + i
!!               write(stdoutUnit,'(i5)',advance='no') iCell
!                  write(stdoutUnit,'(i8)',advance='no') activeForCalvingMask(iCell)
!               enddo
!               write(stdoutUnit,*) ' '
!            enddo
!
!            write(stdoutUnit,*) 'Inactive margin mask'
!            do iRow = nRows, 1, -1
!               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                  write(stdoutUnit,'(a3)',advance='no') '    '
!               endif
!!            do i = 1, nCellsPerRow
!               do i = nCellsPerRow/2 - 2, nCellsPerRow
!                  iCell = (iRow-1)*nCellsPerRow + i
!!               write(stdoutUnit,'(i5)',advance='no') iCell
!                  write(stdoutUnit,'(i8)',advance='no') inactiveMarginMask(iCell)
!               enddo
!               write(stdoutUnit,*) ' '
!            enddo
!
!         endif  ! config_print_calving_info

         ! Calve ice in ocean cells that are not on the protected inactive margin

         where (oceanMask == 1 .and. inactiveMarginMask == 0 .and. thickness > 0.0_RKIND)
            calvingThickness = thickness * calvingFraction
         endwhere
         calvingThicknessFromThreshold = calvingThickness

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

      ! clean up
      call mpas_deallocate_scratch_field(activeForCalvingMaskField, .true.)
      call mpas_deallocate_scratch_field(inactiveMarginMaskField, .true.)
      call mpas_deallocate_scratch_field(oceanMaskField, .true.)

   end subroutine thickness_calving


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine floating_calving
!
!> \brief Calve any ice that is floating
!> \author William Lipscomb, moved to separate subroutine by Matt Hoffman Feb. 2018
!> \date   September 2015
!> \details
!-----------------------------------------------------------------------
   subroutine floating_calving(domain, calvingFraction, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: calvingFraction !< fraction of possible ice to calve on this timestep

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: thickness
      integer, dimension(:), pointer :: cellMask

      err = 0

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         calvingThickness = 0.0_RKIND

         ! Note: The floating_ice mask includes all floating ice, both inactive and active
         where (li_mask_is_floating_ice(cellMask))
            calvingThickness = thickness * calvingFraction
         endwhere
         calvingThicknessFromThreshold = calvingThickness

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

   end subroutine floating_calving

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine remove_small_islands
!
!> \brief Remove very small islands that lead to velocity solver problems
!> \author Matthew Hoffman
!> \date   Summer 2018
!> \details This routine finds and eliminates very small islands that lead to
!> unrealistic velocities in the Albany velocity solver.  Specifically, this
!> finds one- and two-cell masses of ice that are surrounded by open ocean
!> and eliminates them by sending them to the calving flux.
!-----------------------------------------------------------------------

   subroutine remove_small_islands(meshPool, geometryPool)
      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input: Geometry pool

      logical, pointer :: config_remove_small_islands
      real(kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, pointer :: nCellsSolve
      integer :: iCell, jCell, n, nIceNeighbors, nIceNeighbors2, neighborWithIce
      integer :: nOpenOceanNeighbors, nOpenOceanNeighbors2

      call mpas_pool_get_config(liConfigs, 'config_remove_small_islands', config_remove_small_islands)
      if (.not. config_remove_small_islands) then
         return  ! skip this entire routine if disabled
      endif

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

      do iCell = 1, nCellsSolve
         if (li_mask_is_ice(cellMask(iCell))) then ! might as well do for both grounded or floating
                                                   ! (1 or 2 cell floating masses are icebergs)
            nIceNeighbors = 0
            nOpenOceanNeighbors = 0
            do n = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(n, iCell)
               if (li_mask_is_ice(cellMask(jCell))) then
                  nIceNeighbors = nIceNeighbors + 1
                  neighborWithIce = jCell
               endif
               if (.not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level) then
                  nOpenOceanNeighbors = nOpenOceanNeighbors + 1
               endif
            enddo
            if ((nIceNeighbors == 0) .and. (nOpenOceanNeighbors == nEdgesOnCell(iCell))) then
               ! If this is a single cell of ice surrounded by open ocean, kill this location
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
               calvingThicknessFromThreshold(iCell) = calvingThicknessFromThreshold(iCell) + thickness(iCell)
               thickness(iCell) = 0.0_RKIND
            elseif (nIceNeighbors == 1) then
               ! check if this neighbor has any additional neighbors with ice
               nIceNeighbors2 = 0
               nOpenOceanNeighbors2 = 0
               do n = 1, nEdgesOnCell(neighborWithIce)
                  jCell = cellsOnCell(n, neighborWithIce)
                  if (li_mask_is_ice(cellMask(jCell))) then
                     nIceNeighbors2 = nIceNeighbors2 + 1
                  endif
                  if (.not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level) then
                     nOpenOceanNeighbors2 = nOpenOceanNeighbors2 + 1
                  endif
               enddo
               if ((nIceNeighbors2 == 1) .and. (nOpenOceanNeighbors2 == nEdgesOnCell(iCell)-1)) then
                  ! <- only neighbor with ice must have been iCell
                  ! kill both cells
                  calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
                  calvingThicknessFromThreshold(iCell) = calvingThicknessFromThreshold(iCell) + thickness(iCell)
                  thickness(iCell) = 0.0_RKIND
                  calvingThickness(neighborWithIce) = calvingThickness(neighborWithIce) + thickness(neighborWithIce)
                  calvingThicknessFromThreshold(neighborWithIce) = calvingThicknessFromThreshold(neighborWithIce) + thickness(neighborWithIce)
                  thickness(neighborWithIce) = 0.0_RKIND
               endif

            endif ! check on nIceNeighbors

         endif ! check if iCell has ice
      end do ! loop over cells

   end subroutine remove_small_islands


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine topographic_calving
!
!> \brief Calve any floating ice existing where ocean bathymetry is deeper than a threshold
!> \author William Lipscomb, moved to separate subroutine by Matt Hoffman Feb. 2018
!> \date   September 2015
!> \details calve ice where the bed topography lies below a threshold depth
!-----------------------------------------------------------------------
   subroutine topographic_calving(domain, calvingFraction, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      real (kind=RKIND), intent(in) :: calvingFraction !< fraction of possible ice to calve on this timestep

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool
      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold
      real(kind=RKIND), pointer :: config_calving_topography
      real(kind=RKIND), pointer :: config_sea_level
      logical, pointer :: config_print_calving_info
      real (kind=RKIND), dimension(:), pointer :: bedTopography, thickness
      integer, dimension(:), pointer :: cellMask

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_topography', config_calving_topography)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      if (config_print_calving_info) then
         call mpas_log_write('Calving topographic threshold (m) = $r', realArgs=(/config_calving_topography/))
      endif

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         calvingThickness = 0.0_RKIND

         where ( (li_mask_is_floating_ice(cellMask)) .and. (bedTopography < config_calving_topography + config_sea_level) )
            calvingThickness = thickness * calvingFraction
         endwhere
         calvingThicknessFromThreshold = calvingThickness

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

   end subroutine topographic_calving


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine eigencalving
!
!> \brief Calve ice from the calving front based on a strain rate threshold
!> \author Matthew Hoffman
!> \date   Feb. 2018
!> \details Calve using eigencalving scheme in which calving rate is taken
!>  proportional to the product of principle strain rates, if both extensional,
!>  and zero otherwise. Described in detail in:
!>  Levermann, A., T. Albrecht, R. Winkelmann, M. A. Martin, M. Haseloff, and
!>  I. Joughin (2012), Kinematic first-order calving law implies potential for
!>  abrupt ice-shelf retreat, Cryosph., 6(2), 273–286, doi:10.5194/tc-6-273-2012.
!-----------------------------------------------------------------------
   subroutine eigencalving(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool
      real(kind=RKIND), pointer :: config_calving_eigencalving_parameter_scalar_value
      character (len=StrKIND), pointer :: config_calving_eigencalving_parameter_source
      character (len=StrKIND), pointer :: config_damage_calving_method
      real(kind=RKIND), pointer :: config_damage_calving_threshold
      logical, pointer :: config_print_calving_info
      real(kind=RKIND), pointer :: config_calving_thickness
      real (kind=RKIND), dimension(:), pointer :: eigencalvingParameter
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), pointer :: calvingCFLdt
      real (kind=RKIND), pointer :: dtCalvingCFLratio
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: calvingFrontMask
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, pointer :: nCells
      integer :: iCell, jCell, iNeighbor
      logical :: dynamicNeighbor
      real(kind=RKIND) :: calvingSubtotal
      integer :: err_tmp
      logical :: applyToGrounded, applyToFloating, applyToGroundingLine
      real (kind=RKIND), pointer :: totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume

      err = 0

      ! Logical arrays needed for apply_ablation_velocity
      applyToGrounded = .true.
      applyToFloating = .true.
      applyToGroundingLine = .false.

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_scalar_value', &
               config_calving_eigencalving_parameter_scalar_value)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_source', &
              config_calving_eigencalving_parameter_source)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_threshold', config_damage_calving_threshold)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_method', config_damage_calving_method)

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! get fields
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'calvingFrontMask', calvingFrontMask)
         call mpas_pool_get_array(geometryPool, 'eigencalvingParameter', eigencalvingParameter)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         call mpas_pool_get_array(velocityPool, 'eMax', eMax)
         call mpas_pool_get_array(velocityPool, 'eMin', eMin)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
         call mpas_pool_get_array(geometryPool, 'dtCalvingCFLratio', dtCalvingCFLratio)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedCalvedVolume', totalRatebasedCalvedVolume)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedUncalvedVolume', totalRatebasedUncalvedVolume)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

         ! get parameter value
         if (trim(config_calving_eigencalving_parameter_source) == 'scalar') then
            eigencalvingParameter = config_calving_eigencalving_parameter_scalar_value
         elseif (trim(config_calving_eigencalving_parameter_source) == 'data') then
            ! do nothing - use value from input file
         else
            err = 1
            call mpas_log_write("Invalid value specified for option config_calving_eigencalving_parameter_source" // &
                  config_calving_eigencalving_parameter_source, MPAS_LOG_ERR)
         endif

         if (config_print_calving_info) then
            call mpas_log_write("eigencalvingParameter (m s) value range: Min=$r, Max=$r", &
                    realArgs=(/minval(eigencalvingParameter), maxval(eigencalvingParameter)/))
         endif

         calvingVelocity(:) = 0.0_RKIND
         ! First calculate the front retreat rate (Levermann eq. 1)
         calvingVelocity(:) = eigencalvingParameter(:) * max(0.0_RKIND, eMax(:)) * max(0.0_RKIND, eMin(:)) ! m/s

         call mpas_log_write("calling li_apply_front_ablation_velocity from eigencalving")
         ! Convert calvingVelocity to calvingThickness
         call li_apply_front_ablation_velocity(meshPool, geometryPool, velocityPool, &
                                              calvingThickness, calvingVelocity, applyToGrounded, &
                                              applyToFloating, applyToGroundingLine, domain, calvingCFLdt, dtCalvingCFLratio, &
                                              totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume, err_tmp)
         err = ior(err, err_tmp)

         if ( trim(config_damage_calving_method) == 'none' ) then
            ! do nothing
         elseif ( trim(config_damage_calving_method) == 'threshold' ) then
            ! remove ice exceeding damage threshold
            call mpas_log_write('config_damage_calving_method == threshold; &
                                 removing ice with damage > $r', realArgs=(/config_damage_calving_threshold/))

            call apply_calving_damage_threshold(meshPool, geometryPool, scratchPool, domain, err_tmp)
            err = ior(err, err_tmp)
         elseif ( trim(config_damage_calving_method) == 'calving_rate' ) then
             call mpas_log_write('config_damage_calving_method == calving_rate &
                                  is not supported with config_calving == eigencalving', MPAS_LOG_ERR)
             err = 1
             return
         else
             call mpas_log_write('Invalid setting for config_damage_calving_method', MPAS_LOG_ERR)
             err = 1
             return
         endif
         ! Update halos on calvingThickness or faceMeltingThickness before
         ! applying it.
         ! Testing seemed to indicate this is not necessary, but I don't
         ! understand
         ! why not, so leaving it.
         ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'calvingThickness')
         call mpas_timer_stop("halo updates")
         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! Now also remove thin floating, dynamic ice (based on chosen thickness threshold) after mask is updated.

         !where ((li_mask_is_floating_ice(cellMask) .and. li_mask_is_dynamic_ice(cellMask) .and. &
         ! thickness<config_calving_thickness))
         ! The above commented version removed thin ice anywhere on the shelf.  In theory this seemed preferable,
         ! but in practice it has the tendency to create 'holes' in relatively stagnant areas of ice shelves along
         ! the grounding line.  Once these holes developed they would grow and eventually collapse the shelf from behind.

         ! This criteria below only remove too-thin ice at the new calving front,
         ! meaning just one 'row' of cells per timestep.  This could be expanded to continue
         ! removing ice backward until all connected too-thin ice has been removed.
         ! Tests of the current implementation show reasonable behavior.
         do iCell = 1, nCells
            if (calvingFrontMask(iCell) == 1 .and. thickness(iCell) < config_calving_thickness) then
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
               thickness(iCell) = 0.0_RKIND
            endif
         enddo
         ! TODO: global reduce & reporting on amount of calving generated in this step

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! remove abandoned floating ice (i.e. icebergs) and add it to the calving flux
         ! Defined as: floating ice (dynamic or non-dynamic) that is not adjacent to dynamic ice (floating or grounded)
         ! This won't necessarily find all abandoned ice, but in practice it does a pretty good job at general cleanup
         calvingSubtotal = 0.0_RKIND
         do iCell = 1, nCells
           if (li_mask_is_floating_ice(cellMask(iCell))) then
              ! check neighbors for dynamic ice (floating or grounded)
              dynamicNeighbor = .false.
              do iNeighbor = 1, nEdgesOnCell(iCell)
                 jCell = cellsOnCell(iNeighbor, iCell)
                 if (li_mask_is_dynamic_ice(cellMask(jCell))) dynamicNeighbor = .true.
              enddo
              if (.not. dynamicNeighbor) then  ! calve this ice
                 calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
                 thickness(iCell) = 0.0_RKIND
                 calvingSubtotal = calvingSubtotal + calvingThickness(iCell) * areaCell(iCell)
              endif
           endif
         enddo
         ! TODO: global reduce & reporting on amount of calving generated in this step

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

   end subroutine eigencalving


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine specified_calving_velocity
!
!> \brief use a specified calving velocity given by input data
!> \author Matthew Hoffman
!> \date   Feb. 2018
!> \details we can specify the calving velocity by i) a constant value
!> given by config_calving_velocity_const in namelist and ii) source data
!> specified by calvingVelocityData
!-----------------------------------------------------------------------
   subroutine specified_calving_velocity(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool
      logical, pointer :: config_print_calving_info
      real(kind=RKIND), pointer :: config_calving_thickness
      real(kind=RKIND), pointer :: config_calving_velocity_const
      character (len=StrKIND), pointer :: config_calving_specified_source
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: calvingVelocityData
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), pointer :: calvingCFLdt
      real (kind=RKIND), pointer :: dtCalvingCFLratio
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: calvingFrontMask
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, pointer :: nCells
      integer :: iCell, jCell, iNeighbor
      logical :: dynamicNeighbor
      real(kind=RKIND) :: calvingSubtotal
      real (kind=RKIND), pointer :: totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume
      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_calving_velocity_const', config_calving_velocity_const)
      call mpas_pool_get_config(liConfigs, 'config_calving_specified_source', config_calving_specified_source)

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! get fields
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'calvingFrontMask', calvingFrontMask)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         call mpas_pool_get_array(geometryPool, 'calvingVelocityData', calvingVelocityData)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
         call mpas_pool_get_array(geometryPool, 'dtCalvingCFLratio', dtCalvingCFLratio)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedCalvedVolume', totalRatebasedCalvedVolume)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedUncalvedVolume', totalRatebasedUncalvedVolume)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)


         ! get parameter value
         if (trim(config_calving_specified_source) == 'const') then
             calvingVelocity = config_calving_velocity_const
         elseif (trim(config_calving_specified_source) == 'data') then
             calvingVelocity = calvingVelocityData
         else
            err = ior(err, 1)
            call mpas_log_write("Invalid value specified for option config_calving_specified_source" // &
                  config_calving_specified_source, MPAS_LOG_ERR)
         endif

         if (config_print_calving_info) then
            call mpas_log_write("calvingVelocity(m s) value range on this processor: Min=$r, Max=$r", &
                    realArgs=(/minval(calvingVelocity), maxval(calvingVelocity)/))
         endif

         ! Convert calvingVelocity to calvingThickness
         call li_apply_front_ablation_velocity(meshPool, geometryPool, velocityPool, calvingThickness, calvingVelocity, &
                                              applyToGrounded=.true., applyToFloating=.true., applyToGroundingLine=.false., &
                                              domain=domain, maxDt=calvingCFLdt, CFLratio=dtCalvingCFLratio, &
                                              totalAblatedVolume=totalRatebasedCalvedVolume, &
                                              totalUnablatedVolume=totalRatebasedUncalvedVolume, &
                                              err=err_tmp)
         err = ior(err, err_tmp)

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! Now also remove thin floating, dynamic ice (based on chosen thickness threshold) after mask is updated.

         !where ((li_mask_is_floating_ice(cellMask) .and. li_mask_is_dynamic_ice(cellMask) .and. &
         ! thickness<config_calving_thickness))
         ! The above commented version removed thin ice anywhere on the shelf.  In theory this seemed preferable,
         ! but in practice it has the tendency to create 'holes' in relatively stagnant areas of ice shelves along
         ! the grounding line.  Once these holes developed they would grow and eventually collapse the shelf from behind.

         ! This criteria below only remove too-thin ice at the new calving front,
         ! meaning just one 'row' of cells per timestep.  This could be expanded to continue
         ! removing ice backward until all connected too-thin ice has been removed.
         ! Tests of the current implementation show reasonable behavior.
         do iCell = 1, nCells
            if (calvingFrontMask(iCell) == 1 .and. thickness(iCell) < config_calving_thickness) then
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
               thickness(iCell) = 0.0_RKIND
            endif
         enddo
         ! TODO: global reduce & reporting on amount of calving generated in this step

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! remove abandoned floating ice (i.e. icebergs) and add it to the calving flux
         ! Defined as: floating ice (dynamic or non-dynamic) that is not adjacent to dynamic ice (floating or grounded)
         ! This won't necessarily find all abandoned ice, but in practice it does a pretty good job at general cleanup
         calvingSubtotal = 0.0_RKIND
         do iCell = 1, nCells
           if (li_mask_is_floating_ice(cellMask(iCell))) then
              ! check neighbors for dynamic ice (floating or grounded)
              dynamicNeighbor = .false.
              do iNeighbor = 1, nEdgesOnCell(iCell)
                 jCell = cellsOnCell(iNeighbor, iCell)
                 if (li_mask_is_dynamic_ice(cellMask(jCell))) dynamicNeighbor = .true.
              enddo
              if (.not. dynamicNeighbor) then  ! calve this ice
                 calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
                 thickness(iCell) = 0.0_RKIND
                 calvingSubtotal = calvingSubtotal + calvingThickness(iCell) * areaCell(iCell)
              endif
           endif
         enddo
         ! TODO: global reduce & reporting on amount of calving generated in this step

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

   end subroutine specified_calving_velocity


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine von_Mises_calving
!
!> \brief Apply calving law based on von Mises stress to grounded marine margins
!> \author Trevor Hillebrand
!> \date   April 2020
!> \details Uses the calving law described by Morlighem et al. (2016):
!> "Modeling of Store Gletscher’s calving dynamics, West Greenland, in response
!> to ocean thermal forcing" to calculate calvingVelocity at grounded tidewater
!glacier margins
!-----------------------------------------------------------------------
   subroutine von_Mises_calving(domain, err)

      use li_diagnostic_vars

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, jCell, iNeighbor, nGroundedNeighbors, err_tmp
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool, &
                                        velocityPool, scratchPool, thermalPool
      real (kind=RKIND), pointer :: config_grounded_von_Mises_threshold_stress, &
                                    config_floating_von_Mises_threshold_stress, &
                                    config_flowLawExponent, config_calving_speed_limit, &
                                    config_damage_calving_threshold
      character (len=StrKIND), pointer :: config_grounded_von_Mises_threshold_stress_source, &
                                          config_floating_von_Mises_threshold_stress_source
      character (len=StrKIND), pointer :: config_flowParamA_calculation
      character (len=StrKIND), pointer :: config_damage_calving_method
      logical, pointer :: config_use_Albany_flowA_eqn_for_vM
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin, &
                                        calvingVelocity, thickness, &
                                        xvelmean, yvelmean, calvingThickness, &
                                        floatingVonMisesThresholdStress, &
                                        groundedVonMisesThresholdStress, &
                                        surfaceSpeed
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA, &
                                        temperature, layerThickness
      real (kind=RKIND), pointer :: config_default_flowParamA
      real (kind=RKIND), pointer :: calvingCFLdt
      real (kind=RKIND), pointer :: dtCalvingCFLratio
      integer, pointer :: nCells
      integer, dimension(:), pointer :: nCellsArray  ! array giving index to owned cells, plus each halo level
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: indexToCellID
      real (kind=RKIND), dimension(:), pointer :: vonMisesStress, damage
      real (kind=RKIND), pointer :: totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume
      logical :: applyToGrounded, applyToFloating, applyToGroundingLine
      logical :: badCell, badNeighbor
      integer :: kCell, jNeighbor, nGoodNeighbors
      real (kind=RKIND) :: vMsum
      real(kind=RKIND), dimension(:), allocatable :: vonMisesStressOld

      err = 0

      ! Logical arrays needed for li_apply_front_ablation_velocity
      applyToGrounded = .true.
      applyToFloating = .true.
      applyToGroundingLine = .false.

      call mpas_pool_get_config(liConfigs, 'config_grounded_von_Mises_threshold_stress_source', config_grounded_von_Mises_threshold_stress_source)
      call mpas_pool_get_config(liConfigs, 'config_floating_von_Mises_threshold_stress_source', config_floating_von_Mises_threshold_stress_source)
      call mpas_pool_get_config(liConfigs, 'config_grounded_von_Mises_threshold_stress', config_grounded_von_Mises_threshold_stress)
      call mpas_pool_get_config(liConfigs, 'config_floating_von_Mises_threshold_stress', config_floating_von_Mises_threshold_stress)
      call mpas_pool_get_config(liConfigs, 'config_calving_speed_limit', config_calving_speed_limit)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_threshold', config_damage_calving_threshold)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_method', config_damage_calving_method)
      call mpas_pool_get_config(liConfigs, 'config_flowParamA_calculation', config_flowParamA_calculation)

      !call mpas_pool_get_config(liConfigs, 'config_default_flowParamA',
      !config_default_flowParamA) ! REMOVE THIS ONCE YOU CAN GET A FROM
      !ALBANY!!!!!
      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', config_flowLawExponent)
      call mpas_pool_get_config(liConfigs, 'config_use_Albany_flowA_eqn_for_vM', config_use_Albany_flowA_eqn_for_vM)

      block => domain % blocklist

      do while (associated(block))
         ! get pools
          call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
          call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
          call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
          call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         ! get fields
          call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
          call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
          call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
          call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
          call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness)
          call mpas_pool_get_array(velocityPool, 'eMax', eMax)
          call mpas_pool_get_array(velocityPool, 'eMin', eMin)
          call mpas_pool_get_array(velocityPool, 'vonMisesStress', vonMisesStress)
          call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
          call mpas_pool_get_array(velocityPool, 'xvelmean', xvelmean)
          call mpas_pool_get_array(velocityPool, 'yvelmean', yvelmean)
          call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)
          call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
          call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
          call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
          call mpas_pool_get_array(geometryPool, 'dtCalvingCFLratio', dtCalvingCFLratio)
          call mpas_pool_get_array(geometryPool, 'thickness', thickness)
          call mpas_pool_get_array(geometryPool, 'groundedVonMisesThresholdStress', groundedVonMisesThresholdStress)
          call mpas_pool_get_array(geometryPool, 'floatingVonMisesThresholdStress', floatingVonMisesThresholdStress)
          call mpas_pool_get_array(geometryPool, 'totalRatebasedCalvedVolume', totalRatebasedCalvedVolume)
          call mpas_pool_get_array(geometryPool, 'totalRatebasedUncalvedVolume', totalRatebasedUncalvedVolume)
          call mpas_pool_get_array(thermalPool, 'temperature', temperature)

          ! get parameter value and check that values are valid
          if (trim(config_grounded_von_Mises_threshold_stress_source) == 'scalar') then
                groundedVonMisesThresholdStress(:) = config_grounded_von_Mises_threshold_stress
          elseif (trim(config_grounded_von_Mises_threshold_stress_source) == 'data') then
                ! do nothing - use value from input file
          else
                err = 1
                call mpas_log_write("Invalid value specified for option config_grounded_von_Mises_threshold_stress_source" // &
                      config_grounded_von_Mises_threshold_stress_source, MPAS_LOG_ERR)
          endif

          if ( minval(groundedVonMisesThresholdStress(1:nCells)) <= 0.0_RKIND ) then
                err = 1
                call mpas_log_write("groundedVonMisesThresholdStress must be >0.0", MPAS_LOG_ERR)
          endif

          if (trim(config_floating_von_Mises_threshold_stress_source) == 'scalar') then
                floatingVonMisesThresholdStress(:) = config_floating_von_Mises_threshold_stress
          elseif (trim(config_floating_von_Mises_threshold_stress_source) == 'data') then
                ! do nothing - use value from input file
          else
                err = 1
                call mpas_log_write("Invalid value specified for option config_floating_von_Mises_threshold_stress_source" // &
                      config_floating_von_Mises_threshold_stress_source, MPAS_LOG_ERR)
          endif

          if ( minval(floatingVonMisesThresholdStress(:)) < 0.0_RKIND ) then
                err = 1
                call mpas_log_write("floatingVonMisesThresholdStress must be >=0.0", MPAS_LOG_ERR)
          endif

          ! If von Mises threshold stresses contain invalid values, do not
          ! continue.
          if ( err == 1 ) then
             return
          endif

          vonMisesStress(:) = 0.0_RKIND
          
          ! get flowParamA from MPAS or use Albany-like equation
          if ( config_use_Albany_flowA_eqn_for_vM ) then
          !calculate Albany-type flowParamA
             call mpas_log_write("config_use_Albany_flowA_eqn_for_vM not yet supported", MPAS_LOG_ERR)
             err = 1
          else
             if (trim(config_flowParamA_calculation) /= "PB1982") then
                call mpas_log_write("If using von Mises calving, config_flowParamA_calculation must be set to 'PB1982'", MPAS_LOG_ERR)
                err = 1
                return
             endif
             call li_calculate_flowParamA(meshPool, temperature, thickness,flowParamA,err) ! Get MPAS flowParamA
          endif

          !Using a depth-averaged ice viscosity parameter B_depthAvg
          !=sum(layerThickness(:,iCell) *
          !flowParamA(:,iCell)**(-1.0_RKIND/config_flowLawExponent), dim=1) /
          !thickness(iCell)
          ! Calculate effective von  Mises stress.
          calvingVelocity(:) = 0.0_RKIND

          do iCell = 1,nCells
             if ( thickness(iCell) > 0.0_RKIND ) then
                vonMisesStress(iCell) = sqrt(3.0_RKIND) * sum(layerThickness(:,iCell) * &
                              flowParamA(:,iCell)**(-1.0_RKIND/config_flowLawExponent),dim=1) / &
                              thickness(iCell) * ( 0.5_RKIND * ( (max(0.0_RKIND, eMax(iCell)))**2.0_RKIND + &
                              (max(0.0_RKIND, eMin(iCell)))**2.0_RKIND) )**(1.0_RKIND / (2.0_RKIND * config_flowLawExponent))
             endif
          enddo

          ! halo update on vonMisesStress before extrap
          call mpas_timer_start("halo updates")
          call mpas_dmpar_field_halo_exch(domain, 'vonMisesStress')
          call mpas_timer_stop("halo updates")

          ! the row of cells that are adjacent to invalid velocity will have garbage strain rates and vM stress
          ! so extrapolate into those cells from their valid neighbors
          do iCell = 1, nCellsArray(2)  ! this is owned cells plus first halo
             if (surfaceSpeed(iCell) == 0.0_RKIND) then
                vonMisesStress(iCell) = 0.0_RKIND  ! remove the vM calculation where there is no velocity - garbage can occur here
             else
                ! only clean up cells that have a valid velocity - cells outside the solution can be left alone
                badCell = .false.
                do iNeighbor = 1, nEdgesOnCell(iCell)
                   jCell = cellsOnCell(iNeighbor, iCell)
                   if ((jCell < nCells + 1) .and. surfaceSpeed(jCell) == 0.0_RKIND) then
                      ! if one neighbor has a 0 velo, the stress here is bad
                      ! Ignore nonexistent neighbors (nCells+1)
                      badCell = .true.
                      exit ! no need to keep looping over neighbors
                   endif
                enddo
                if (badCell) then
                   vMsum = 0.0_RKIND
                   nGoodNeighbors = 0
                   do iNeighbor = 1, nEdgesOnCell(iCell)
                      jCell = cellsOnCell(iNeighbor, iCell)
                      if (jCell < nCells + 1) then  ! don't check nonexistent neighbors
                         badNeighbor = .false.
                         do jNeighbor = 1, nEdgesOnCell(jCell)
                            kCell = cellsOnCell(jNeighbor, jCell)
                            if ((kCell < nCells + 1) .and. (surfaceSpeed(kCell) == 0.0_RKIND)) then
                               ! if one neighbor of this neighbor has a 0 velo, this neighbor is not to be used
                               ! Ignore nonexistent neighbors
                               badNeighbor = .true.
                               exit
                            endif
                         enddo
                         if (.not. badNeighbor) then
                            ! if this neighbor has a usable stress, accumulate its value for averaging
                            vMsum = vMsum + vonMisesStress(jCell)
                            nGoodNeighbors = nGoodNeighbors + 1
                         endif
                      endif
                   enddo
                   if (nGoodNeighbors > 0) then
                      vonMisesStress(iCell) = vMsum / real(nGoodNeighbors, kind=RKIND)
                   else
                      call mpas_log_write('Extrapolation failed for von mises stress at cell $i', MPAS_LOG_WARN, intArgs=(/indexToCellID(iCell)/))
                   endif
                endif ! if badCell
             endif ! if velo>0
          enddo ! loop over nCells

          ! no halo update needed yet because we included the first halo above

          ! Now extrapolate one additional row forward for good measure
          allocate(vonMisesStressOld(nCells+1)) ! temp var for extrap
          vonMisesStressOld(:) = vonMisesStress(:)
          do iCell = 1, nCellsArray(1)  ! just owned cells (equivalent to nCellsSolve)
             ! Find cells needing extrap
             if (vonMisesStressOld(iCell) == 0.0) then
                vMsum = 0.0_RKIND
                nGoodNeighbors = 0
                do iNeighbor = 1, nEdgesOnCell(iCell)
                   jCell = cellsOnCell(iNeighbor, iCell)
                   ! Use this neighbor if it has a nonzero vM value
                   if ((jCell < nCells + 1) .and. vonMisesStressOld(jCell) > 0.0_RKIND) then
                      vMsum = vMsum + vonMisesStressOld(jCell)
                      nGoodNeighbors = nGoodNeighbors + 1
                   endif
                enddo
                if (nGoodNeighbors > 0) then
                   vonMisesStress(iCell) = vMsum / real(nGoodNeighbors, kind=RKIND)
                endif
             endif ! if vMOld=0
          enddo ! loop over nCells
          deallocate(vonMisesStressOld)

          ! another halo update on vonMisesStress after extrap
          call mpas_timer_start("halo updates")
          call mpas_dmpar_field_halo_exch(domain, 'vonMisesStress')
          call mpas_timer_stop("halo updates")

          do iCell = 1,nCells
             ! Calculate calving velocity for grounded cells at marine margin
             if ( .not. li_mask_is_floating_ice(cellMask(iCell)) ) then
                calvingVelocity(iCell) = min(sqrt(xvelmean(iCell)**2.0_RKIND + yvelmean(iCell)**2.0_RKIND) * &
                                         vonMisesStress(iCell) / groundedVonMisesThresholdStress(iCell), config_calving_speed_limit)
             ! If config_floating_von_Mises_threshold_stress is not 0.0, calculate
             ! calvingVelocity. If config_floating_von_Mises_threshold_stress is
             ! 0.0, remove floating ice in loop below.
             elseif ( li_mask_is_floating_ice(cellMask(iCell)) .and. config_floating_von_Mises_threshold_stress .ne. 0.0_RKIND) then
                calvingVelocity(iCell) = min(sqrt(xvelmean(iCell)**2 + yvelmean(iCell)**2) * &
                                         vonMisesStress(iCell) / floatingVonMisesThresholdStress(iCell), config_calving_speed_limit)
             endif
          enddo

          call mpas_timer_start("halo updates")
          call mpas_dmpar_field_halo_exch(domain, 'calvingVelocity')
          call mpas_timer_stop("halo updates")

          call mpas_log_write("calling li_apply_front_ablation_velocity from von Mises stress calving routine")
         ! Convert calvingVelocity to calvingThickness
          call li_apply_front_ablation_velocity(meshPool, geometryPool,velocityPool, &
                                              calvingThickness, calvingVelocity, applyToGrounded, &
                                              applyToFloating, applyToGroundingLine, domain, &
                                              calvingCFLdt, dtCalvingCFLratio, &
                                              totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume, err_tmp)
          err = ior(err, err_tmp)
         ! Update halos on calvingThickness or faceMeltingThickness before
         ! applying it.
         ! Testing seemed to indicate this is not necessary, but I don't
         ! understand
         ! why not, so leaving it.
         ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'calvingThickness')
         call mpas_timer_stop("halo updates")

          ! If floating VM threshold is zero, remove floating ice. Remove
          ! non-dynamic cells adjacent to floating ice, but leave
          ! non-dynamic cells adjacent to grounded ice.
          ! TODO: Make an exception for the case of floating ice surrounded by grounded ice
          ! (i.e., floating on a subglacial lake) using a flood fill routine
          if ( config_floating_von_Mises_threshold_stress .eq. 0.0_RKIND ) then
             do iCell = 1,nCells
                if ( li_mask_is_floating_ice(cellMask(iCell)) .and. &
                     li_mask_is_dynamic_ice(cellMask(iCell)) ) then
                     calvingThickness(iCell) = thickness(iCell)
                elseif ( li_mask_is_floating_ice(cellMask(icell)) .and. &
                    (.not. li_mask_is_dynamic_ice(cellMask(iCell))) ) then
                    nGroundedNeighbors = 0
                    do iNeighbor = 1, nEdgesOnCell(iCell)
                       jCell = cellsOnCell(iNeighbor, iCell)
                       if ( li_mask_is_grounded_ice(cellMask(jCell)) ) then
                           nGroundedNeighbors = nGroundedNeighbors + 1
                       endif
                    enddo
                    if ( nGroundedNeighbors == 0 ) then
                       calvingThickness(iCell) = thickness(iCell)
                    endif
                endif
             enddo
          endif

          if ( trim(config_damage_calving_method) == 'none' ) then
             ! do nothing
          elseif ( trim(config_damage_calving_method) == 'threshold' ) then
             ! remove ice exceeding damage threshold
             call mpas_log_write('config_damage_calving_method == threshold; &
                                  removing ice with damage > $r', realArgs=(/config_damage_calving_threshold/))
            
             call apply_calving_damage_threshold(meshPool, geometryPool, scratchPool, domain, err_tmp)
             err = ior(err, err_tmp)
          elseif ( trim(config_damage_calving_method) == 'calving_rate' ) then
              call mpas_log_write('config_damage_calving_method == calving_rate &
                                   is not supported with config_calving == von_Mises_stress', MPAS_LOG_ERR)
              err = 1
              return
          else
              call mpas_log_write('Invalid setting for config_damage_calving_method', MPAS_LOG_ERR)
              err = 1
              return
          endif

          ! Update halos on calvingThickness or faceMeltingThickness before
          ! applying it.
          ! Testing seemed to indicate this is not necessary, but I don't
          ! understand
          ! why not, so leaving it.
          ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
          call mpas_timer_start("halo updates")
          call mpas_dmpar_field_halo_exch(domain, 'calvingThickness')
          call mpas_timer_stop("halo updates")

          ! === apply calving ===
          thickness(:) = thickness(:) - calvingThickness(:)

          ! update mask
          call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
          err = ior(err, err_tmp)

          call remove_small_islands(meshPool, geometryPool)
          
          block => block % next

      enddo ! associated(block)

   end subroutine von_Mises_calving

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine ismip6_retreat
!
!> /brief Use ISMIP6 retreat parameterization for Greenland (Slater et al., 2019, 2020).
!> /author Trevor Hillebrand
!> /date November 2021
!> /details This routine applies the ISMIP6 glacier retreat parateterization
!> based on ocean thermal forcing and subglacial discharge. This routine uses multiple
!> time-levels fields provided by ISMIP6. Once we use E3SM fields, this will need to be restructured.
!> The routine searches for an input stream named 'ismip6-gis', from which it reads ismip6Runoff and 
!> ismip6_2dThermalForcing fields from the previous forcing time whenever the mostRecentAccessTime
!> attribute of the ismip6-gis stream is updated. It then calculates calvingVelocity based on
!> the change in ismip6Runoff and ismip6_2dThermalForcing between the two most recent time-levels, and
!> calls apply_front_ablation_velocity to calculate a calvingThickness. config_front_mass_bal_grounded
!> must be set to 'none' when using this routine.
!> Slater, D. A., Straneo, F., Felikson, D., Little, C. M., Goelzer, H.,
!> Fettweis, X., & Holte, J. (2019).
!> Estimating Greenland tidewater glacier retreat driven by submarine melting.
!> The Cryosphere, 13(9), 2489-2509.
!> https://doi.org/10.5194/tc-13-2489-2019
!> Slater, D. A., Felikson, D., Straneo, F., Goelzer, H., Little, C. M.,
!> Morlighem, M., et al. (2020).
!> Twenty-first century ocean forcing of the Greenland ice sheet for modelling
!> of sea level contribution.
!> The Cryosphere, 14(3), 985-1008. https://doi.org/10.5194/tc-14-985-2020
   subroutine ismip6_retreat(domain, err) 

      use li_diagnostic_vars
      use mpas_timekeeping
      use mpas_stream_manager
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: & 
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: iCell, jCell, iEdge, iNeighbor, err_tmp
      type (block_type), pointer :: block
      integer, dimension(:,:), pointer :: cellsOnCell
      type (mpas_pool_type), pointer :: geometryPool, meshPool, &
                                        velocityPool
      real (kind=RKIND), pointer :: config_ismip6_retreat_k, seaLevel
      logical, pointer :: config_do_restart
      real (kind=RKIND), dimension(:), pointer ::  &
                                        calvingVelocity, thickness, &
                                        xvelmean, yvelmean, calvingThickness, &
                                        bedTopography
      real (kind=RKIND), pointer :: calvingCFLdt
      real (kind=RKIND), pointer :: dtCalvingCFLratio
      integer, pointer :: nCells, timestepNumber
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      real (kind=RKIND) :: waterDepth
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: TFocean, ismip6Runoff
      real (kind=RKIND), dimension(:), pointer :: TFoceanPrevious, ismip6RunoffPrevious
      real (kind=RKIND), dimension(:), pointer :: TFoceanCurrent, ismip6RunoffCurrent
      real (kind=RKIND), dimension(:), allocatable :: submergedArea
      real (kind=RKIND) :: deltatForcing ! time between forcing updates
      type (MPAS_Time_Type), save :: forcingTime, forcingTimeOld
      character (len=StrKIND), pointer :: forcingTimeStamp
      character(len=StrKIND) :: forcingTimeOldStamp
      type (MPAS_stream_list_type), pointer :: stream_cursor
      logical :: streamFound ! used to throw an error if required stream is not found
      real (kind=RKIND), pointer :: totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume

      call mpas_pool_get_config(liConfigs, 'config_ismip6_retreat_k', config_ismip6_retreat_k)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', seaLevel)
      call mpas_pool_get_config(liConfigs, 'config_do_restart', config_do_restart)

      if (config_ismip6_retreat_k .ge. 0.0) then
         call mpas_log_write('Error: config_ismip6_retreat_k should be negative, but is >= 0.0.')
         err = 1
      endif

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
      ! get fields
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'forcingTimeStamp', forcingTimeStamp)
      call mpas_pool_get_array(meshPool, 'timestepNumber', timestepNumber)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
      call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
      call mpas_pool_get_array(geometryPool, 'dtCalvingCFLratio', dtCalvingCFLratio)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(velocityPool, 'xvelmean', xvelmean)
      call mpas_pool_get_array(velocityPool, 'yvelmean', yvelmean)
      call mpas_pool_get_array(geometryPool, 'ismip6_2dThermalForcing', TFocean)
      call mpas_pool_get_array(geometryPool, 'ismip6Runoff', ismip6Runoff)
      call mpas_pool_get_array(geometryPool, 'ismip6_2dThermalForcingPrevious', TFoceanPrevious)
      call mpas_pool_get_array(geometryPool, 'ismip6RunoffPrevious', ismip6RunoffPrevious)
      call mpas_pool_get_array(geometryPool, 'ismip6_2dThermalForcingCurrent', TFoceanCurrent)
      call mpas_pool_get_array(geometryPool, 'ismip6RunoffCurrent', ismip6RunoffCurrent)
      call mpas_pool_get_array(geometryPool, 'totalRatebasedCalvedVolume', totalRatebasedCalvedVolume)
      call mpas_pool_get_array(geometryPool, 'totalRatebasedUncalvedVolume', totalRatebasedUncalvedVolume)
      
      ! submergedArea used in runoff unit conversion 
      allocate(submergedArea(nCells+1))

      streamFound = .false. ! Changed to true if ismip6-gis stream is found, otherwise throws error
      stream_cursor => domain % streamManager % streams % head
      do while (associated(stream_cursor))
         ! Input stream with forcings must be called 'ismip6-gis', or this will throw an error.
         if ( trim(stream_cursor % name) == 'ismip6-gis' .and. (stream_cursor % valid) ) then
            streamFound = .true.
            if ( (timestepNumber == 0) .or. (config_do_restart .and. (timestepNumber == 1)) ) then
               ! On the first timestep of a cold start OR restart initialize all needed fields by forcing a read of
               ! the previous forcing data on the first timestep.  This could be potentially be an unnecessary read on the
               ! initial time of a restart run, but that is a small cost for simplifying logic and making the initial time of both
               ! cold starts and restarts behave the same.  The potentially extra read does not affect the algorithm.

               ! Use forcing fields and time from most recent time in file, which was read on init.
               call mpas_get_time(stream_cursor%mostRecentAccessTime, dateTimeString=forcingTimeStamp, ierr=err_tmp)
               err = ior(err, err_tmp)

               TFoceanCurrent = TFocean
               ismip6RunoffCurrent = ismip6Runoff

               ! Force a read of this stream, and use the second most recent time in the file
               call mpas_stream_mgr_read(domain % streamManager, streamID = stream_cursor % name, rightNow = .true., &
                                when = forcingTimeStamp, whence = MPAS_STREAM_LATEST_STRICTLY_BEFORE, saveActualWhen= .true., ierr=err_tmp)
               err = ior(err, err_tmp)

               call mpas_get_time(stream_cursor%mostRecentAccessTime, dateTimeString=forcingTimeOldStamp, ierr=err_tmp)
               err = ior(err, err_tmp)

               TFoceanPrevious = TFocean
               ismip6RunoffPrevious = ismip6Runoff

               call mpas_log_write("  * Forced a read of input stream 'ismip6-gis'" // &
                                         " from time: " // trim(forcingTimeOldStamp))

               call mpas_set_time(forcingTime, dateTimeString=forcingTimeStamp, ierr=err_tmp)
               err = ior(err, err_tmp)
               call mpas_set_time(forcingTimeOld, dateTimeString=forcingTimeOldStamp, ierr=err_tmp)
               err = ior(err, err_tmp)

               ! Do not overwrite forcing fields with values from previous time
               ismip6Runoff = ismip6RunoffCurrent
               TFocean = TFoceanCurrent
            elseif (forcingTime < stream_cursor % mostRecentAccessTime) then
               TFoceanPrevious = TFoceanCurrent
               ismip6RunoffPrevious = ismip6RunoffCurrent
               TFoceanCurrent = TFocean
               ismip6RunoffCurrent = ismip6Runoff

               forcingTimeOld = forcingTime
               forcingTime = stream_cursor % mostRecentAccessTime
               call mpas_get_time(forcingTime, dateTimeString=forcingTimeStamp, ierr=err_tmp)
               err = ior(err, err_tmp)
               call mpas_get_time(forcingTimeOld, dateTimeString=forcingTimeOldStamp, ierr=err_tmp)
               err = ior(err, err_tmp)
               call mpas_log_write("  * ismip6 retreat forcings have been updated:" // &
                                " forcingTime: " // trim(forcingTimeStamp) //  &
                                "; forcingTimeOld: " // trim(forcingTimeOldStamp))
            else
               call mpas_log_write(" * No new forcing times for ismip6 retreat this timestep.")
            endif
            exit ! We have processed the stream we wanted.
         endif
         stream_cursor => stream_cursor % next
      enddo ! end loop over stream_cursors

      ! check that the ismip6-gis stream was found
      if (.not. streamFound) then
         call mpas_log_write('Input stream ismip6-gis is required for config_calving = ismip6_retreat, ' // &
                             'but was not found.', MPAS_LOG_ERR)
         err = 1
      endif
      ! Get submerged area of each cell to convert from kg m^{-2} s^{-1} to m^3 s^{-1}
      ! TODO: ensure that this unit conversion is still appropriate when we
      ! eventually use this subroutine with fields passed from the E3SM coupler
      ! rather than external ISMIP6 forcings.
      submergedArea(:) = 0.0_RKIND
      do iCell = 1, nCells
         if ( li_mask_is_dynamic_margin(cellMask(iCell)) .and. bedTopography(iCell) < seaLevel ) then
            waterDepth = max(0.0_RKIND, seaLevel - bedTopography(iCell))
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               ! sum up length of edges adjacent to open ocean or non-dynamic cells
               if ( (.not. li_mask_is_dynamic_ice(cellMask(jCell))) &
                       .and. (bedTopography(jCell) < seaLevel) ) then
                  submergedArea(iCell) = submergedArea(iCell) + dvEdge(iEdge) * waterDepth
               endif
            enddo
         endif
      enddo

      ! Get the time between forcings to calculate retreat rate.
      call mpas_get_timeInterval(forcingTime - forcingTimeOld, dt=deltatForcing)
      call mpas_log_write("ismip6 retreat deltatForcing = $r", realArgs = (/deltatForcing/))

      ! This if-statement is probably unnecessary, but this ensures there is no divide-by-zero error
      ! if somehow deltatForcing = 0.0.
      if (forcingTime .ne. forcingTimeOld) then
         do iCell=1, nCells
            ! For now, hard-coding factor of 4.0 to roughly convert mean annual to mean summer runoff.
            ! Factor of 1000.0 is freshwater density.
            ! TODO: ensure that this scaling factor is still appropriate when we
            ! eventually use this subroutine with fields passed from the E3SM coupler
            ! rather than external ISMIP6 forcings.
            calvingVelocity(iCell) = max(-1.0_RKIND * config_ismip6_retreat_k * ( ( (4.0_RKIND * submergedArea(iCell) &
                                         * ismip6Runoff(iCell) / 1000.0_RKIND)**0.4_RKIND * TFocean(iCell) ) -  &
                                         ( (4.0_RKIND * submergedArea(iCell) * ismip6RunoffPrevious(iCell) / 1000.0_RKIND)**0.4_RKIND &
                                         * TFoceanPrevious(iCell) ) ) / deltatForcing + &
                                         sqrt(xvelmean(iCell)**2.0_RKIND + yvelmean(iCell)**2.0_RKIND), 0.0_RKIND)
         enddo
      else ! throw an error if the forcing times are the same
          call mpas_log_write("Forcings used in ismip6_retreat have the same timestamp. Check xtime in the input file(s).")
          err = 1
      endif

      ! Testing indicates this halo update is necessary before calling li_apply_front_ablation_velocity.
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'calvingVelocity')
      call mpas_timer_stop("halo updates")

      call mpas_log_write("calling li_apply_front_ablation_velocity from ismip6 retreat routine")
      ! Convert calvingVelocity to calvingThickness
      call li_apply_front_ablation_velocity(meshPool, geometryPool,velocityPool, &
                                              calvingThickness, calvingVelocity, applyToGrounded=.true., &
                                              applyToFloating=.true., applyToGroundingLine=.false., &
                                              domain=domain, maxDt=calvingCFLdt, CFLratio=dtCalvingCFLratio, &
                                              totalAblatedVolume=totalRatebasedCalvedVolume, &
                                              totalUnablatedVolume=totalRatebasedUncalvedVolume, &
                                              err=err_tmp)
      err = ior(err, err_tmp)

      ! Update halos on calvingThickness before applying it.
      ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'calvingThickness')
      call mpas_timer_stop("halo updates")

      ! === apply calving ===
      thickness(:) = thickness(:) - calvingThickness(:)

      ! update mask
      call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
      err = ior(err, err_tmp)
      call remove_small_islands(meshPool, geometryPool)

      deallocate(submergedArea)

   end subroutine ismip6_retreat
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine li_apply_front_ablation_velocity
!
!> \brief Convert a calving or melting velocity to an ice thickness removal
!> \author Matthew Hoffman and Trevor Hillebrand
!> \date   March 2020
!> \details This routine takes a calving or front melting velocity and converts it to a
!> thickness removal.  The basic idea is to multiply the calving velocity
!> by the edge length and height to get a volume flux.  However, the details
!> are complicated by the unstructured Voronoi mesh and the use of partially-filled
!> non-dynamic floating grid cells.
!> The algorithm calculates the required calving volume flux based on edges, and then
!> applies the volume removal to grid cells, converting that to a thickness removal using areaCell.
!> This process occurs in phases.
!>
!> The first phase handles non-dynamic floating cells. This is done by calculating the required calving flux
!> at the edges between those non-dynamic floating cells and the open ocean.
!> The calving velocity and ice thickness are copied from the upstream dynamic
!> floating cells, because these quantities do not have valid values on the non-dynamic cells.
!> After calculating the calving flux determined at each of these edges, the non-dyamic floating cells
!> are looped over, and volume is marked for removal based on the calving flux on the edges of each cell.
!>
!> The final phase handles dynamic floating cells that are eligible for calving.
!> There are two ways this can happen - either they are adjacent to the open ocean (no buffer of
!> non-dynamic cells present), or they are adjacent to non-dynamic floating cells that had their
!> entire volume removed by calving, but there still is calving demand remaining.
!> The first case is handled easily by calculating calving flux on edges between dynamic floating cells
!> and open ocean using the calving velocity and ice thickness of the dynamic cell itself.
!> In the second case, "leftover" required calving flux is distributed from a non-dynamic cell that
!> has been completely "drained" to any edges between that non-dynamic cell and dynamic cells upstream.
!> This is done by weighting the flux to be distrbuted by the length of each such edge in the direction
!> perpendicular to the calving flux.
!>
!> The location at which ablation is controlled by applyToGrounded, applyToFloating, and applyToGroundingLine,
!> which are set by the calling routine.
!> If applyToGrounded = .true., ablation is applied to grounded marine margin cells, which do not have dynamic ice shelf neighbors
!> If applyToGroundingLine = .true., ablation is applied to grounding line, regardless of whether there is an ice shelf or not.
!> If there in dynamic floating ice, ablation is applied to the floating cell adjacent to grounding line;
!> otherwise it is applied to the last grounded cell.
!> If applyToFloating = .true., ablation is applied to dynamic floating margin cells
!> The output of this routine is calvingThickness or faceMeltingThickness, which then needs to be applied to thickness
!> by the calling routine.

!> IMPORTANT: The calling routine may need to perform a halo update on calving or melting velocity
!> before the call. Testing shows this is necessary for von Mises calving, but not for eigencalving, for instance.
!> The calling routine should always perform a halo update on calving or melting thickness after the call and before applying ablation.
!-----------------------------------------------------------------------

   subroutine li_apply_front_ablation_velocity(meshPool, geometryPool, velocityPool, ablationThickness, ablationVelocity, &
                   applyToGrounded, applyToFloating, applyToGroundingLine, domain, maxDt, CFLratio, &
                   totalAblatedVolume, totalUnablatedVolume, err)

      use ieee_arithmetic, only : ieee_is_nan

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), pointer, intent(in) :: velocityPool !< Input: velocity pool
      logical, intent(in) :: applyToFloating, applyToGrounded, &
                             applyToGroundingLine
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain          !< Input: domain object
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input: geometry pool
      real (kind=RKIND), dimension(:), intent(inout) :: ablationVelocity
      real (kind=RKIND), dimension(:), intent(out) :: ablationThickness
      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      real (kind=RKIND), optional :: maxDt
          !< Output: an approximation of the max allowable dt based on a CFL-like condition calculated from ablation velocity
      real (kind=RKIND), optional :: CFLratio
          !< Output: the ratio of the actual timestep being applied to the maximum allowable timestep from the CFL-like condition
      real (kind=RKIND), optional :: totalAblatedVolume
          !< Output: the total ablated volume
      real (kind=RKIND), optional :: totalUnablatedVolume
          !< Output: the total unablated volume
      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nEdges, nEdgesSolve, nCells, nCellsSolve, maxEdges
      integer :: iEdge, iCell, jCell, kCell, iNeighbor, jNeighbor
      integer :: nEmptyNeighbors, nGroundedNeighbors, counter, nTwoCellsBack, nOneCellBack
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: lowerSurface
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: cellMask, edgeMask
      integer, dimension(:), pointer :: frontAblationMask
      integer, dimension(:), pointer :: openOceanMask
      real (kind=RKIND), dimension(:), pointer :: calvingThickness, faceMeltingThickness
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: config_calving_error_threshold
      logical, pointer :: config_distribute_unablatedVolumeDynCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity, faceMeltSpeed
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructY
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      real (kind=RKIND), dimension(:), allocatable :: thicknessForAblation
      real (kind=RKIND), dimension(:), allocatable :: uvelForAblation, vvelForAblation
      integer, dimension(:), allocatable :: oneCellBackList, twoCellsBackList
      real (kind=RKIND), dimension(:), pointer :: requiredAblationVolumeNonDynEdge, requiredAblationVolumeNonDynCell
      real (kind=RKIND), dimension(:), pointer :: requiredAblationVolumeDynEdge, requiredAblationVolumeDynCell
      real (kind=RKIND), dimension(:), pointer :: ablatedVolumeNonDynCell, ablatedVolumeDynCell
      real (kind=RKIND), dimension(:), pointer :: unablatedVolumeNonDynCell, unablatedVolumeDynCell
      real (kind=RKIND), dimension(:), allocatable :: cellVolume
      real(kind=RKIND) :: ablationSubtotal1, ablationSubtotal2, ablationSubtotal3
      real(kind=RKIND) :: thkSum, ablationVelSum, uvelSum, vvelSum
      integer :: thkCount
      real(kind=RKIND) :: removeVolumeHere
      real(kind=RKIND) :: volumeAvailableToPass !< a temp. var. for accounting purposes that indicates how much volume
                          !< to be calved is 'leftover' after trying to calve non-dynamic cells and can be transferred to
                          !< neighboring dynamic cells
      real(kind=RKIND) :: ablateLengthCell, ablateLengthEdge
      real(kind=RKIND), dimension(6) :: localInfo, globalInfo
      real(kind=RKIND) :: edgeLengthScaling
      real(kind=RKIND), parameter :: ablationSmallThk = 1.0e-8_RKIND ! in meters, a small thickness threshold
      real(kind=RKIND) :: minOfMaxAllowableDt
      integer :: err_tmp, nDynNeighbors

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_calving_error_threshold', config_calving_error_threshold)
      call mpas_pool_get_config(liConfigs, 'config_distribute_unablatedVolumeDynCell', config_distribute_unablatedVolumeDynCell)

      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
      call mpas_pool_get_array(geometryPool, 'faceMeltSpeed', faceMeltSpeed)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)
      call mpas_pool_get_array(geometryPool, 'unablatedVolumeNonDynCell', unablatedVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'unablatedVolumeDynCell', unablatedVolumeDynCell)
      call mpas_pool_get_array(geometryPool, 'ablatedVolumeNonDynCell', ablatedVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'ablatedVolumeDynCell', ablatedVolumeDynCell)
      call mpas_pool_get_array(geometryPool, 'frontAblationMask', frontAblationMask)
      call mpas_pool_get_array(geometryPool, 'openOceanMask', openOceanMask)
      call mpas_pool_get_array(geometryPool, 'requiredAblationVolumeNonDynEdge', requiredAblationVolumeNonDynEdge)
      call mpas_pool_get_array(geometryPool, 'requiredAblationVolumeDynEdge', requiredAblationVolumeDynEdge)
      call mpas_pool_get_array(geometryPool, 'requiredAblationVolumeNonDynCell', requiredAblationVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'requiredAblationVolumeDynCell', requiredAblationVolumeDynCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)


      ! Step 0: Define the frontAblationMask that will be used to determine where ablation should be applied.
      ! This mask can be for 1. a floating calving front, 2. a grounded marine margin (no dynamic floating ice shelf attached),
      ! 3. along the grounding line (whether or not there is floating ice adjacent), or any combination of the three.
      ! The mask includes both the dynamic AND non-dynamic cells along the prescribed interface type(s).
      frontAblationMask(:) = 0

      ! First identify the dynamic cells for each mask type.

      if ( applyToGrounded .or. applyToGroundingLine ) then
         do iCell = 1, nCells
         ! If applyToGrounded: Define marine marginal cells as those that are (1) at the ice
         ! margin, (2) have at least one neighboring cell without ice, (3) contain
         ! grounded ice, and (4) have bed topography below sea level.
         ! OR is adjacent to an inactive floating margin cell
            if (li_mask_is_grounded_ice(cellMask(iCell)) &
               ! GL here means cell is grounded but has a neighbor that is floating or ocean.
               ! Note that as of Oct 2022, this is no longer the general
               ! definition of the grounding line throughout the code.
               .and. bedTopography(iCell) < config_sea_level &
               .and. li_mask_is_dynamic_ice(cellMask(iCell)) ) then

               ! Check if neighboring cells contain ice and have bed topo below sea level
               nEmptyNeighbors = 0
               do iEdge = 1, nEdgesOnCell(iCell)
                  iNeighbor = cellsOnCell(iEdge, iCell)
                  if ( (((thickness(iNeighbor) == 0.0_RKIND) .and. bedTopography(iNeighbor) < config_sea_level)) &
                     !< these previous conditions are open ocean
                     .or. &
                     !> these following conditions are inactive floating margin cells
                     ((li_mask_is_floating_ice(cellMask(iNeighbor)) &
                     .and. li_mask_is_margin(cellMask(iNeighbor)) &
                     .and. (.not. li_mask_is_dynamic_ice(cellMask(iNeighbor)))))) then
                     nEmptyNeighbors = nEmptyNeighbors + 1
                  endif
               enddo
               if (nEmptyNeighbors > 0) then
                  frontAblationMask(iCell) = 1
               endif ! nEmptyNeighbors
            endif ! cell is GL etc.
         enddo ! nCells loop
      endif ! applyToGrounded

      if ( applyToGroundingLine ) then
         ! If applyToGroundingLine: Define marine marginal cells as those that
         ! are at the margin (defined by li_mask_is_grounding_line) if no
         ! floating ice or are the first floating cell if there is floating ice.
         ! The case of no floating ice was handled above already, so only looking for
         ! GL locations with floating ice here.
         do iCell = 1,nCells
            if ( li_mask_is_grounding_line(cellMask(iCell)) ) then
               do iEdge = 1, nEdgesOnCell(iCell)
                  iNeighbor = cellsOnCell(iEdge, iCell)
                  if (li_mask_is_floating_ice(cellMask(iNeighbor)) .and. li_mask_is_dynamic_ice(cellMask(iNeighbor))) then
                     frontAblationMask(iNeighbor) = 1
                  endif
               enddo
            endif
         enddo
      endif

      if ( applyToFloating ) then
         call find_open_ocean(domain, openOceanMask, err_tmp)
         err = ior(err, err_tmp)
         do iCell = 1,nCells
            ! If applyToFloating: marine marginal cells are dynamic floating margin cells
            if ( li_mask_is_floating_ice(cellMask(iCell)) .and. &
               li_mask_is_dynamic_margin(cellMask(iCell)) ) then
               do iEdge = 1, nEdgesOnCell(iCell)
                  iNeighbor = cellsOnCell(iEdge, iCell)
                  ! only mark this cell if it neighbors the open ocean (don't apply calving to holes in ice shelf)
                  if (openOceanMask(iNeighbor) == 1) then
                     frontAblationMask(iCell) = 1
                     exit ! no need to keep searching for neighbors
                  endif
               enddo
            endif
         enddo
      endif

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'frontAblationMask')
      call mpas_timer_stop("halo updates")

      ! Now extend the mask forward to include non-dynamic floating neighbors (appropriate for any of the mask types)
      ! Use value 2 to indicate the non-dynamic cells
      do iCell = 1,nCells
         if (frontAblationMask(iCell) == 1) then
            do iEdge = 1, nEdgesOnCell(iCell)
               iNeighbor = cellsOnCell(iEdge, iCell)
               if ( li_mask_is_floating_ice(cellMask(iNeighbor)) .and. &
                  .not. li_mask_is_dynamic_ice(cellMask(iNeighbor)) ) then
                  frontAblationMask(iNeighbor) = 2 ! Use 2 for these non-dynamic floating neighbors
               endif
            enddo
         endif
      enddo
      ! Make a second pass to include potentially 'stranded' nondynamic floating cells
      do iCell = 1, nCells
         if (frontAblationMask(iCell) == 2) then
            ! Check if this non-dynamic masked cell has any non-dynamic, floating, but unmasked neighbors we need to include
            do iEdge = 1, nEdgesOnCell(iCell)
               iNeighbor = cellsOnCell(iEdge, iCell)
               if ((.not. li_mask_is_dynamic_ice(cellMask(iNeighbor))) .and. li_mask_is_floating_ice(cellMask(iNeighbor)) &
                  .and. frontAblationMask(iNeighbor) == 0) then
                  frontAblationMask(iNeighbor) = 3 ! Use 3 for these extended non-dynamic neighbors
               endif
            enddo
         endif
      enddo

      ! Calculate an approximate ablation CFL limiting dt
      if (present(maxDt)) then
         minOfMaxAllowableDt = 1.0e16_RKIND ! Initialize to large number
         do iEdge = 1, nEdgesSolve
            if ((frontAblationMask(cellsOnEdge(1,iEdge)) > 0) .or. (frontAblationMask(cellsOnEdge(2,iEdge)) > 0)) then
               minOfMaxAllowableDt = min(minOfMaxAllowableDt, &
                    0.5_RKIND * dcEdge(iEdge) / &
                    (max(ablationVelocity(cellsOnEdge(1,iEdge)), ablationVelocity(cellsOnEdge(2,iEdge))) + 1.0E-18_RKIND) )
            endif
         enddo
      endif

      ! Init fields for accounting
      ablationThickness(:) = 0.0_RKIND
      ablatedVolumeNonDynCell(:) = 0.0_RKIND
      ablatedVolumeDynCell(:) = 0.0_RKIND
      unablatedVolumeNonDynCell(:) = 0.0_RKIND
      unablatedVolumeDynCell(:) = 0.0_RKIND
      allocate(cellVolume(nCells+1))
      cellVolume(:) = areaCell(:) * thickness(:)
      requiredAblationVolumeNonDynEdge(:) = 0.0_RKIND
      requiredAblationVolumeNonDynCell(:) = 0.0_RKIND
      requiredAblationVolumeDynEdge(:) = 0.0_RKIND
      requiredAblationVolumeDynCell(:) = 0.0_RKIND
      allocate(thicknessForAblation(nCells+1))
      thicknessForAblation = thickness
      allocate(uvelForAblation(nCells+1))
      uvelForAblation(:) = uReconstructX(1,:)
      allocate(vvelForAblation(nCells+1))
      vvelForAblation(:) = uReconstructY(1,:)
      allocate(oneCellBackList(maxEdges+1))
      oneCellBackList(:) = 0
      allocate(twoCellsBackList((maxEdges+1)**2))
      twoCellsBackList(:) = 0

      ! 1. Calculate ablation rate for all non-dynamic cells by working with their ocean-going edges
      ablationSubtotal1 = 0.0_RKIND
      do iCell = 1, nCells
         if ((frontAblationMask(iCell)>0) .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell))) ) then
            ! a1. Translate the calving front height based on dynamic cells to the non-dynamic locations
            ! find mean (or min?) of thickness in dynamic neighbors
            ablationVelSum = 0.0_RKIND
            thkSum = 0.0_RKIND
            uvelSum = 0.0_RKIND
            vvelSum = 0.0_RKIND
            thkCount = 0
            nGroundedNeighbors = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if ( li_mask_is_dynamic_ice(cellMask(jCell)) .and. li_mask_is_floating_ice(cellMask(jCell)) ) then
                  thkSum = thkSum + thicknessForAblation(jCell)
                  uvelSum = uvelSum + uvelForAblation(jCell)
                  vvelSum = vvelSum + vvelforAblation(jCell)
                  ablationVelSum = ablationVelSum + ablationVelocity(jCell)
                  thkCount = thkCount + 1
               else if ( li_mask_is_grounding_line(cellMask(jCell)) .or. li_mask_is_grounded_ice(cellMask(jCell)) ) then
                  nGroundedNeighbors = nGroundedNeighbors + 1
               endif
            enddo

            ! if floating non-dynamic cell has *any* grounded neighbors, use ice thickness and
            ! velocity from two cells back instead of the one-cell-back (neighboring) locations above
            ! Note the two-cell back averaging is independent of whether those locations
            ! are floating are grounded.
            ! The method implemented here gives good results for Humboldt Glacier (realistic grounded margin),
            ! but as of 3/16/21 it has not been tested at a realistic case that contains both grounded and floating
            ! margins (e.g. an ice shelf in an embayment with grounded ice along the lateral margins).  If unexpected
            ! behavior occurs in such a case in the future, an alternative implementation here might work better:
            ! Use 2-back values to populate 1-back locations that are grounded, and then use 1-back values to calculate
            ! the value at each primary cell location.  This would be a mix of true-1-back values for dynamic floating cells
            ! and the corresponding 2-back values at grounded cells.
            if ( nGroundedNeighbors > 0 ) then
               thkSum = 0.0_RKIND
               uvelSum = 0.0_RKIND
               vvelSum = 0.0_RKIND
               ablationVelSum = 0.0_RKIND
               thkCount = 0
               nOneCellBack = 0
               oneCellBackList(:) = 0
               twoCellsBackList(:) = 0
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  jCell = cellsOnCell(iNeighbor, iCell)
                  if ( li_mask_is_ice(cellMask(jCell)) ) then
                        nOneCellBack = nOneCellBack + 1
                        oneCellBackList(nOneCellBack) = jCell
                  endif
               enddo
               if ( nOneCellBack > 0 ) then
               ! loop over oneCellBackList and add neighbors
                  counter = 0
                  do iNeighbor = 1, nOneCellBack
                     jCell = oneCellBackList(iNeighbor)
                     ! create twoCellsBackList so we don't add one cell multiple
                     ! times
                     nTwoCellsBack = 0
                     do jNeighbor = 1, nEdgesOnCell(jCell)
                        kCell = cellsOnCell(jNeighbor, jCell)
                        if ( li_mask_is_dynamic_ice(cellMask(kCell)) .and. &
                           (.not. any(oneCellBackList==kCell) ) .and. &
                           (.not. any(twoCellsBackList==kCell) ) ) then
                             counter = counter+1
                             twoCellsBackList(counter) = kCell
                        endif
                     enddo
                     nTwoCellsBack = nTwoCellsBack + counter
                  enddo
                  ! loop through twoCellsBackList and average thickness and
                  ! velocities
                  do kCell = 1, nTwoCellsBack
                     thkSum = thkSum + thicknessForAblation(twoCellsBackList(kCell))
                     uvelSum = uvelSum + uvelForAblation(twoCellsBackList(kCell))
                     vvelSum = vvelSum + vvelforAblation(twoCellsBackList(kCell))
                     ablationVelSum = ablationVelSum + ablationVelocity(twoCellsBackList(kCell))
                     thkCount = thkCount + 1
                  enddo
               endif
            endif
            if (thkCount == 0) then
               !call mpas_log_write("Found a stranded non-dynamic floating
               !cell: cell $i with thickness=$r m.", MPAS_LOG_WARN, &
               !   intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
            else
               thicknessForAblation(iCell) = thkSum / real(thkCount, kind=RKIND)
               ablationVelocity(iCell) = ablationVelSum / real(thkCount, kind=RKIND)
               uvelForAblation(iCell) = uvelSum / real(thkCount, kind=RKIND)
               vvelForAblation(iCell) = vvelSum / real(thkCount, kind=RKIND)
            endif
         endif
      enddo


          ! a2.Set thicknessForAblation and ablationVelocity of stranded non-dynamic cells to
          ! average of non-dynamic neighbors. Stranded non-dynamic cells are
          ! identified as those with frontAblationMask == 3
      do iCell = 1, nCells
         if ( frontAblationMask(iCell) == 3 ) then
            ablationVelSum = 0.0_RKIND
            thkSum = 0.0_RKIND
            uvelSum = 0.0_RKIND
            vvelSum = 0.0_RKIND
            thkCount = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if ( li_mask_is_floating_ice(cellMask(jCell)) &
                    .and. (.not. li_mask_is_dynamic_ice(cellMask(jCell))) ) then
                  thkSum = thkSum + thicknessForAblation(jCell)
                  uvelSum = uvelSum + uvelForAblation(jCell)
                  vvelSum = vvelSum + vvelForAblation(jCell)
                  ablationVelSum = ablationVelSum + ablationVelocity(jCell)
                  thkCount = thkCount + 1
               endif
            enddo
            if (thkCount == 0) then
                !call mpas_log_write("Found a stranded non-dynamic floating cell: cell $i with thickness=$r m.", MPAS_LOG_WARN, &
               !   intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
            else
                thicknessForAblation(iCell) = thkSum / real(thkCount, kind=RKIND)
                ablationVelocity(iCell) = ablationVelSum / real(thkCount, kind=RKIND)
                uvelForAblation(iCell) = uvelSum / real(thkCount, kind=RKIND)
                vvelForAblation(iCell) = vvelSum / real(thkCount, kind=RKIND)
            endif
         endif
      enddo

            ! b. Translate the ablationVelocity and thickness from non-dynamic cells to their ocean-going edges
            !    to calculate the calving volume on those edges and this cell
      do iCell = 1, nCells
         if (frontAblationMask(iCell)>0 .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell)))) then

            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_margin(edgeMask(iEdge)) .and. & !< edge is a margin
                   .not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level & ! ensure margin is w/ocn
                   ) then
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uvelForAblation(iCell), vvelForAblation(iCell))
                  requiredAblationVolumeNonDynEdge(iEdge) = ablationVelocity(iCell) * &
                          edgeLengthScaling * dvEdge(iEdge) * thicknessForAblation(iCell) * deltat
                  requiredAblationVolumeNonDynCell(iCell) = requiredAblationVolumeNonDynCell(iCell) + &
                          requiredAblationVolumeNonDynEdge(iEdge) ! Keep running total
               endif
            enddo

            ! c. Apply ablationThickness here
            removeVolumeHere = min(cellVolume(iCell), requiredAblationVolumeNonDynCell(iCell))  ! Don't use more than available
            ablationThickness(iCell) = removeVolumeHere / areaCell(iCell)
            ablatedVolumeNonDynCell(iCell) = removeVolumeHere
            unablatedVolumeNonDynCell(iCell) = requiredAblationVolumeNonDynCell(iCell) - removeVolumeHere
            if (bedTopography(iCell) >= config_sea_level) then
               unablatedVolumeNonDynCell(iCell) = 0.0_RKIND ! Ignore leftover potential calving when we hit sea level
            endif
            cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere
            if (iCell <= nCellsSolve) ablationSubtotal1 = ablationSubtotal1 + removeVolumeHere
         endif
      enddo
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'unablatedVolumeNonDynCell')
      call mpas_timer_stop("halo updates")
      !call mpas_log_write("Done calculating calving for nondynamic floating cells. Removed $r m^3", realArgs=(/calvingSubtotal1/))



      ! 2. Calculate calving for dynamic cells
      ablationSubtotal2 = 0.0_RKIND

      ! 2a1. for grounded dynamic margins, don't use the thickness there because it may be unrealistically thin due to ablation
      ! from previous time steps.  Instead copy forward the thickness from the row cells behind.
      do iCell = 1, nCells
         if (frontAblationMask(iCell)>0 .and. li_mask_is_margin(cellMask(iCell)) .and. &
                 li_mask_is_grounded_ice(cellMask(iCell))) then
            thkSum = 0.0_RKIND
            thkCount = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(cellMask(jCell)) .and. li_mask_is_grounded_ice(cellMask(jCell)) .and. &
                       .not. li_mask_is_margin(cellMask(jCell))) then
                  thkSum = thkSum + thicknessForAblation(jCell)
                  thkCount = thkCount + 1
               endif
            enddo
            if (thkCount == 0) then
               !call mpas_log_write("Found a stranded non-dynamic floating cell: cell $i with thickness=$r m.", MPAS_LOG_WARN, &
               !   intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
            else
               thicknessForAblation(iCell) = thkSum / real(thkCount, kind=RKIND)
            endif
         endif
      enddo

      ! a. Calculate calving on dynamic margin edges
      do iCell = 1, nCells
         if ( frontAblationMask(iCell)>0 .and. li_mask_is_dynamic_ice(cellMask(iCell)) & ! a dynamic cell in the mask...
            .and. li_mask_is_margin(cellMask(iCell)) ) then                            ! that is at the edge of the ice
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if ((.not. li_mask_is_ice(cellMask(jCell))) .and. (bedTopography(jCell) < config_sea_level)) then
                  iEdge = edgesOnCell(iNeighbor, iCell)
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uvelForAblation(iCell), vvelForAblation(iCell))
                  requiredAblationVolumeDynEdge(iEdge) = ablationVelocity(iCell) * &
                          edgeLengthScaling * dvEdge(iEdge) * thicknessForAblation(iCell) * deltat
               endif
            enddo
         endif
      enddo
      ! b. copy ablation remaining in non-dynamic cells to dynamic edges
      ! Assume height and velocity are uniform, but edge length is not.
      do iCell = 1, nCells
         if ( (.not. li_mask_is_dynamic_ice(cellMask(iCell))) .and. &
                 unablatedVolumeNonDynCell(iCell) > 0.0_RKIND) then
            ! This is a non-dynamic location that still has calving to offer - a location where calving needs to be propagated

            !call mpas_log_write("Passing calving from non-dynamic cell $i.  $r available to pass", &
            !        realArgs=(/unablatedVolumeNonDynCell(iCell)/), intArgs=(/iCell/))

            volumeAvailableToPass = unablatedVolumeNonDynCell(iCell)

            ! Find total length of interface with dynamic cells at this cell
            ablateLengthCell = 0.0_RKIND
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(edgeMask(iEdge))) then
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uvelForAblation(iCell), vvelForAblation(iCell))
                  ablateLengthEdge = edgeLengthScaling * dvEdge(iEdge)
                  ablateLengthCell = ablateLengthCell + ablateLengthEdge
               endif
            enddo

            ! Now that we know ablateLengthCell, pass along the required calving volume relative to the interface length
            if ( ablateLengthCell > 0.0_RKIND ) then
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(iNeighbor, iCell)
                  jCell = cellsOnCell(iNeighbor, iCell)
                  if (li_mask_is_dynamic_ice(edgeMask(iEdge))) then
                     edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uvelForAblation(iCell), vvelForAblation(iCell))
                     ablateLengthEdge = edgeLengthScaling * dvEdge(iEdge)
                     if (requiredAblationVolumeDynEdge(iEdge) > 0.0_RKIND) then
                        call mpas_log_write("Unexpectedly found a dynamic edge that already has calving assigned to it." // &
                              "  There is a flaw in li_apply_front_ablation_velocity that needs to be fixed!", MPAS_LOG_ERR)
                        err_tmp = 1
                        err = ior(err, err_tmp)
                     endif
                     requiredAblationVolumeDynEdge(iEdge) = ablateLengthEdge / ablateLengthCell * volumeAvailableToPass
                     unablatedVolumeNonDynCell(iCell) = unablatedVolumeNonDynCell(iCell) - requiredAblationVolumeDynEdge(iEdge)
                     !call mpas_log_write("   Passed calving $r from non-dynamic cell $i to dynamic cell $i", &
                     !   realArgs=(/requiredAblationVolumeDynEdge(iEdge)/), intArgs=(/iCell, jCell/)) ! Note: change to global IDs
                     if (bedTopography(jCell) >= config_sea_level) then
                        requiredAblationVolumeDynEdge(iEdge) = 0.0_RKIND ! Don't pass this further if it would go above sea level
                     endif
                  endif
               enddo
            endif
         endif
      enddo
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'requiredAblationVolumeDynEdge')
      call mpas_timer_stop("halo updates")

      ! c. Now apply ablation to each cell
      do iCell = 1, nCells
         if ((li_mask_is_dynamic_ice(cellMask(iCell))) .and. (bedTopography(iCell) < config_sea_level)) then
            ! Can loop over all dyn cells - only ones with calving on their edges will have calving applied.
            ! Note, we ignore any leftover calving/melting if it would be applied to cells where the bed is above sea level.
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               ! No need to check what type of edge this is - only required edges are nonzero
               requiredAblationVolumeDynCell(iCell) = requiredAblationVolumeDynCell(iCell) + &
                          requiredAblationVolumeDynEdge(iEdge) ! Keep running total
            enddo
            ! c. Apply calving here
            removeVolumeHere = min(cellVolume(iCell), requiredAblationVolumeDynCell(iCell))  ! Don't use more than available
            ablationThickness(iCell) = removeVolumeHere / areaCell(iCell)
            ablatedVolumeDynCell(iCell) = removeVolumeHere
            unablatedVolumeDynCell(iCell) = requiredAblationVolumeDynCell(iCell) - removeVolumeHere
            cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere
            if (iCell <= nCellsSolve) ablationSubtotal2 = ablationSubtotal2 + removeVolumeHere
         endif
      enddo

      ! Clean up to account for roundoff level errors that can occur
      do iCell = 1, nCells
         if (abs(ablationThickness(iCell) - thickness(iCell)) < ablationSmallThk) then
            ablationThickness(iCell) = thickness(iCell)
         endif
      enddo

      if ( config_distribute_unablatedVolumeDynCell ) then
         ! 3. Distribute unablatedVolumeDynCell among neighboring cells. This is
         ! necesssary when a cell is fully depleted but still has not ablated as
         ! much ice as required by the ablation parameterization. This still may
         ! result in large error for very high ablation rates, which should be
         ! prevented by limiting the timestep based on the ablation velocity.

         ! a. Calculate volume to be distributed evenly among neighboring dynamic cells
         do iCell = 1, nCells
            if ( (unablatedVolumeDynCell(iCell) > 0.0_RKIND) ) then
               nDynNeighbors = 0
               ! Count neighbors between which to evenly distribute
               ! unablatedVolumeDynCell
               do iEdge = 1, nEdgesOnCell(iCell)
                  iNeighbor = cellsOnCell(iEdge, iCell)
                  if ((li_mask_is_dynamic_ice(cellMask(iNeighbor))) .and. (bedTopography(iNeighbor) < config_sea_level) &
                       .and. ( cellVolume(iNeighbor) > (ablationSmallThk * areaCell(iNeighbor)) ) ) then
                     nDynNeighbors = nDynNeighbors + 1
                  endif
               enddo
               if ( nDynNeighbors > 0 ) then
                  ! Now distribute unablatedVolumeDynCell between those neighbors
                  do iEdge = 1, nEdgesOnCell(iCell)
                     iNeighbor = cellsOnCell(iEdge, iCell)
                     if ((li_mask_is_dynamic_ice(cellMask(iNeighbor))) .and. (bedTopography(iNeighbor) < config_sea_level) &
                       .and. ( cellVolume(iNeighbor) > (ablationSmallThk * areaCell(iNeighbor)) ) ) then
                        unablatedVolumeDynCell(iNeighbor) = unablatedVolumeDynCell(iNeighbor) + &
                                                                unablatedVolumeDynCell(iCell) / nDynNeighbors
                     endif
                  enddo
                  ! This cell now has no more unablatedVolume
                  unablatedVolumeDynCell(iCell) = 0.0_RKIND
               endif
            endif
         enddo

         ! b. Apply ablation
         do iCell = 1, nCells
             if ((li_mask_is_dynamic_ice(cellMask(iCell))) .and. (bedTopography(iCell) < config_sea_level) &
                  .and. ( cellVolume(iCell) > (ablationSmallThk * areaCell(iCell)) ) &
                  .and. (unablatedVolumeDynCell(iCell) > 0.0_RKIND) ) then
                  removeVolumeHere = min(cellVolume(iCell), unablatedVolumeDynCell(iCell))  ! Don't use more than available
                  ablationThickness(iCell) = ablationThickness(iCell) + removeVolumeHere / areaCell(iCell) ! add to ablationThickness calculated in 2c
                  ablatedVolumeDynCell(iCell) = ablatedVolumeDynCell(iCell) + removeVolumeHere
                  unablatedVolumeDynCell(iCell) = unablatedVolumeDynCell(iCell) - removeVolumeHere
                  cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere
                  ! For debugging it may be helpful to add a separate metric for
                  ! this step.
                  if (iCell <= nCellsSolve) ablationSubtotal2 = ablationSubtotal2 + removeVolumeHere
              endif
         enddo
         !call mpas_log_write("Done calculating ablation for dynamic floating cells. Removed $r m^3", realArgs=(/calvingSubtotal2/))

         ! c. Clean up to account for roundoff level errors that can occur
         do iCell = 1, nCells
            if (abs(ablationThickness(iCell) - thickness(iCell)) < ablationSmallThk) then
               ablationThickness(iCell) = thickness(iCell)
            endif
         enddo
      endif ! config_distribute_unablatedVolumeDynCell

      ! 4. Clean up after applying ablation velocity from parameterization
      ! a. Zap any stranded ice.  This only needs to be considered for floating ice.
      ablationSubtotal3 = 0.0_RKIND
      do iCell = 1, nCells
         if (li_mask_is_floating_ice(cellMask(iCell)) .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell)))) then
            thkCount = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(cellMask(jCell))) then
                  thkCount = thkCount + 1
               endif
            enddo
            if (thkCount == 0) then
               ablationThickness(iCell) = thickness(iCell)
               if (iCell <= nCellsSolve) ablationSubtotal3 = ablationSubtotal3 + ablationThickness(iCell)
            endif
         endif
         if (ieee_is_nan(calvingThickness(iCell))) then
            call mpas_log_write("NaN detected in calvingThickness at cell $i", MPAS_LOG_ERR, intArgs=(/iCell/))
            err_tmp = 1
            err = ior(err, err_tmp)
         endif
      enddo

      ! b. Clean up to account for roundoff level errors that can occur
      do iCell = 1, nCells
         if (abs(ablationThickness(iCell) - thickness(iCell)) < ablationSmallThk) then
            ablationThickness(iCell) = thickness(iCell)
         endif
      enddo

      ! Clean up to unablated volume arrays before checking for errors
      do iCell = 1, nCells
         ! Eliminate leftover potental calving from cells above sea level to avoid erroneous warnings/errors
         if (bedTopography(iCell) >= config_sea_level) then
            unablatedVolumeNonDynCell(iCell) = 0.0_RKIND
            unablatedVolumeDynCell(iCell) = 0.0_RKIND
         endif
         ! Clean up round off.  Because this is volume, threshold can be kind of large
         if (unablatedVolumeNonDynCell(iCell) < (ablationSmallThk * areaCell(iCell))) then
            unablatedVolumeNonDynCell(iCell) = 0.0_RKIND
         endif
         if (unablatedVolumeDynCell(iCell) < (ablationSmallThk * areaCell(iCell))) then
            unablatedVolumeDynCell(iCell) = 0.0_RKIND
         endif
      enddo

      ! End of routine accounting/reporting
      localInfo(1) = ablationSubtotal1
      localInfo(2) = ablationSubtotal2
      localInfo(3) = ablationSubtotal3
      localInfo(4) = sum(ablationThickness(1:nCellsSolve) * areaCell(1:nCellsSolve))
      localInfo(5) = sum(unablatedVolumeNonDynCell(1:nCellsSolve))
      localInfo(6) = sum(unablatedVolumeDynCell(1:nCellsSolve))
      ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
      call mpas_dmpar_sum_real_array(domain % dminfo, 6, localInfo, globalInfo)
      if (present(maxDt)) then
         ! Do this global reduce at the same time as the previous for efficiency
         call mpas_dmpar_min_real(domain % dminfo, minOfMaxAllowableDt, maxDt)
         if (present(CFLratio)) then
            CFLratio = deltat / maxDt
            if (CFLratio > 1.0_RKIND) then
               call mpas_log_write("Ratio of dt to calving CFL dt exceeds 1: $r", MPAS_LOG_WARN, realArgs = (/CFLratio/))
            endif
         endif
      endif
      call mpas_log_write("== Ablation complete. Total calved                    = $r", realArgs = (/globalInfo(4)/))
      call mpas_log_write("==                   Ablated from non-dynamic cells   = $r", realArgs = (/globalInfo(1)/))
      call mpas_log_write("==                   Ablated from dynamic cells       = $r", realArgs = (/globalInfo(2)/))
      call mpas_log_write("==                   Stranded floating cells deleted = $r", realArgs = (/globalInfo(3)/))
      call mpas_log_write("== Unablated volumes: Non-dynamic cells=$r; Dynamic cells=$r", realArgs=(/globalInfo(5),globalInfo(6)/))
      if (((globalInfo(5) + globalInfo(6)) / (globalInfo(4) + 1.0e-30_RKIND) > 0.001_RKIND) .and. &
          (globalInfo(4) > 1000.0_RKIND**2)) then ! Include some small amount of total calving for comparison
         call mpas_log_write("Failed to ablate $r m^3. ($r% of total ablated)", MPAS_LOG_WARN, &
                 realArgs=(/globalInfo(5) + globalInfo(6), &
                 100.0_RKIND * (globalInfo(5) + globalInfo(6)) / (globalInfo(4)+1.0e-30_RKIND)/))
      endif
      if (((globalInfo(5) + globalInfo(6)) / (globalInfo(4) + 1.0e-30_RKIND) > config_calving_error_threshold) .and. &
          (globalInfo(4) > 1000.0_RKIND**2)) then ! Include some small amount of total calving for comparison
         call mpas_log_write("Failed to ablate an amount greater than $r% of the ice ablated.  " // &
                 "Try reducing time step or li_apply_front_ablation_velocity may need improvements.", &
                 MPAS_LOG_ERR, realArgs=(/config_calving_error_threshold * 100.0_RKIND/))
         err_tmp = 1
         err = ior(err, err_tmp)
      endif
      if (globalInfo(3) / (globalInfo(4) + 1.0e-30_RKIND) > 0.01_RKIND) then
         ! If stranded ice deletion results in more than a small amount of the total calving flux, this routine will not
         ! be accurate.  If this error gets triggered, either try using a smaller timstep or this routine needs improving.
         call mpas_log_write("Deleting stranded floating cells accounts for more than 1% of ablation loss." // &
                 "  Try using a smaller timestep or li_apply_front_ablation_velocity may need improvements for this simulation.")
         err_tmp = 1
         err = ior(err, err_tmp)
      endif

      if (present(totalAblatedVolume)) then
         totalAblatedVolume = globalInfo(4)
      endif
      if (present(totalUnablatedVolume)) then
         totalUnablatedVolume = globalInfo(5) + globalInfo(6)
      endif

      deallocate(cellVolume)
      deallocate(thicknessForAblation)
      deallocate(uvelForAblation)
      deallocate(vvelForAblation)
      deallocate(oneCellBackList)
      deallocate(twoCellsBackList)
      call mpas_log_write("Finished with li_apply_front_ablation_velocity")

    end subroutine li_apply_front_ablation_velocity


    ! Helper function for subroutine li_apply_front_ablation_velocity
    ! Calculates the amount to scale an edge length based on the orientation of the edge with the surface velocity
    function scale_edge_length(angleEdgeHere, u, v)
       real(kind=RKIND), intent(in) :: angleEdgeHere
       real(kind=RKIND), intent(in) :: u
       real(kind=RKIND), intent(in) :: v
       real(kind=RKIND) :: scale_edge_length

       real(kind=RKIND) :: mag

       mag = sqrt(u**2 + v**2)
       if (mag == 0.0_RKIND) mag = 1.0_RKIND
       scale_edge_length =  abs(u/mag * cos(angleEdgeHere) + v/mag * sin(angleEdgeHere))  ! dot product of unit vectors
       !scale_edge_length = abs(cos(angleEdgeHere - atan2(v, u)))
    end function scale_edge_length


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine find_open_ocean
!
!> \brief Generate a mask of open ocean (i.e. ignore 'holes' in ice shelf)
!-----------------------------------------------------------------------
   subroutine find_open_ocean(domain, openOceanMask, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, dimension(:), intent(out) :: openOceanMask
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool, scratchPool
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell
      type (field1dInteger), pointer :: seedMaskField, growMaskField
      integer, dimension(:), pointer :: seedMask, growMask, cellMask
      real (kind=RKIND), dimension(:), pointer :: bedTopography, thickness
      real (kind=RKIND), pointer :: config_sea_level
      integer, pointer :: nCells
      integer :: iCell, iNeighbor, jCell

      err = 0

      block => domain % blocklist

      call mpas_pool_get_subpool(domain % blocklist % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
      call mpas_allocate_scratch_field(seedMaskField, single_block_in = .true.)
      seedMask => seedMaskField % array
      call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
      call mpas_allocate_scratch_field(growMaskField, single_block_in = .true.)
      growMask => growMaskField % array

      seedMask(:) = 0
      growMask(:) = 0
      do iCell = 1, nCells
         if ((bedTopography(iCell) < config_sea_level) .and. ( .not. li_mask_is_dynamic_ice(cellMask(iCell)) ) ) then
            growMask(iCell) = 1  ! includes floating non-dynamic cells
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (jCell == nCells + 1) then
                  seedMask(iCell) = 1  ! this is the seed locations - open ocean along domain boundary
                  exit
               endif
            enddo
         endif
      enddo

      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'seedMask')
      call mpas_timer_stop("halo updates")

      call li_flood_fill(seedMask, growMask, domain)
      openOceanMask = seedMask

      call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.true.)
      call mpas_deallocate_scratch_field(growMaskField, single_block_in=.true.)
   end subroutine find_open_ocean


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine damagecalving
!
!> \brief Calve ice from the calving front based on the Bassis-Ma (2015) damage theory.
!> We can use two options: "calving_rate" or "threshold".
!> If we use "calving_rate", the calving front retreats by removing ice at a speed of k * damage,
!> where k is a calving parameter. If we use "threshold", the ice at the calving front will be completely
!> gone if the damage is above the threshold.
!> \author Tong Zhang
!> \date   May. 2019
!> \details Connect calving with the damage model
!> Bassis, Jeremy N., and Y. Ma. "Evolution of basal crevasses links ice shelf stability to ocean forcing."
!> Earth and Planetary Science Letters 409 (2015): 203-211.

!-----------------------------------------------------------------------
   subroutine damage_calving(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool
      logical, pointer :: config_print_calving_info
      real (kind=RKIND), pointer :: config_sea_level
      real(kind=RKIND), pointer :: config_calving_thickness
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), pointer :: calvingCFLdt
      real (kind=RKIND), pointer :: dtCalvingCFLratio
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin
      real (kind=RKIND), dimension(:), pointer :: damage
      integer, dimension(:), pointer :: calvingFrontMask
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, pointer :: nCells
      integer :: iCell, jCell, iNeighbor
      real(kind=RKIND) :: cellCalvingFrontLength, cellCalvingFrontHeight
      integer :: err_tmp
      logical :: dynamicNeighbor
      real(kind=RKIND), pointer :: config_damagecalvingParameter
      real(kind=RKIND) :: calvingSubtotal
      character (len=StrKIND), pointer :: config_damage_calving_method
      real(kind=RKIND), pointer :: config_damage_calving_threshold
      real (kind=RKIND), pointer :: totalRatebasedCalvedVolume, totalRatebasedUncalvedVolume

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_damagecalvingParameter', config_damagecalvingParameter)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_method', config_damage_calving_method)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_threshold', config_damage_calving_threshold)

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get fields
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'damage', damage)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         call mpas_pool_get_array(geometryPool, 'calvingFrontMask', calvingFrontMask)
         call mpas_pool_get_array(geometryPool, 'calvingCFLdt', calvingCFLdt)
         call mpas_pool_get_array(geometryPool, 'dtCalvingCFLratio', dtCalvingCFLratio)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedCalvedVolume', totalRatebasedCalvedVolume)
         call mpas_pool_get_array(geometryPool, 'totalRatebasedUncalvedVolume', totalRatebasedUncalvedVolume)


         call calculate_calving_front_mask(meshPool, geometryPool, calvingFrontMask)

         if (trim(config_damage_calving_method) == 'calving_rate') then
             if (config_damage_calving_threshold >= 1.0_RKIND .or. config_damage_calving_threshold < 0.0_RKIND) then
                call mpas_log_write("Invalid value of config_damage_calving_threshold specified for 'calving_rate' option.  " // &
                        "Valid values are in the range [0,1).", MPAS_LOG_ERR)
                err = ior(err, err_tmp)
             endif
             ! First calculate the front retreat rate (m/s)
             calvingVelocity(:) = config_damagecalvingParameter * &
                max(0.0_RKIND, (damage(:) - config_damage_calving_threshold) / (1.0_RKIND - config_damage_calving_threshold)) &
                * real(li_mask_is_floating_ice_int(cellMask(:)), kind=RKIND) ! calculate only for floating ice
             call li_apply_front_ablation_velocity(meshPool, geometryPool, velocityPool, calvingThickness, calvingVelocity, &
                                              applyToGrounded=.false., applyToFloating=.true., applyToGroundingLine=.false., &
                                              domain=domain, maxDt=calvingCFLdt, CFLratio=dtCalvingCFLratio, &
                                              totalAblatedVolume=totalRatebasedCalvedVolume, &
                                              totalUnablatedVolume=totalRatebasedUncalvedVolume, &
                                              err=err_tmp)
             err = ior(err, err_tmp)
         elseif (trim(config_damage_calving_method) == 'threshold') then
             call apply_calving_damage_threshold(meshPool, geometryPool, scratchPool, domain, err_tmp)
             err = ior(err, err_tmp)
         else
             call mpas_log_write("Unknown value for config_damage_calving_method was specified!", MPAS_LOG_ERR)
         endif

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! Now also remove thin floating, dynamic ice (based on chosen thickness threshold) after mask is updated.
         ! This criteria below only remove too-thin ice at the new calving front,
         ! meaning just one 'row' of cells per timestep.  This could be expanded to continue
         ! removing ice backward until all connected too-thin ice has been removed.
         ! Tests of the current implementation show reasonable behavior.
         do iCell = 1, nCells
            if (calvingFrontMask(iCell) == 1 .and. thickness(iCell) < config_calving_thickness) then
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
               thickness(iCell) = 0.0_RKIND
            endif
         enddo

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! remove abandoned floating ice (i.e. icebergs) and add it to the calving flux
         ! Defined as: floating ice (dynamic or non-dynamic) that is not adjacent to dynamic ice (floating or grounded)
         ! This won't necessarily find all abandoned ice, but in practice it does a pretty good job at general cleanup
         calvingSubtotal = 0.0_RKIND
         do iCell = 1, nCells
           if (li_mask_is_floating_ice(cellMask(iCell))) then
              ! check neighbors for dynamic ice (floating or grounded)
              dynamicNeighbor = .false.
              do iNeighbor = 1, nEdgesOnCell(iCell)
                 jCell = cellsOnCell(iNeighbor, iCell)
                 if (li_mask_is_dynamic_ice(cellMask(jCell))) dynamicNeighbor = .true.
              enddo
              if (.not. dynamicNeighbor) then  ! calve this ice
                 calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
                 thickness(iCell) = 0.0_RKIND
                 calvingSubtotal = calvingSubtotal + calvingThickness(iCell) * areaCell(iCell)
              endif
           endif
         enddo
         ! TODO: global reduce & reporting on amount of calving generated in this step

         call remove_small_islands(meshPool, geometryPool)

         block => block % next

      enddo

   end subroutine damage_calving

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine li_calculate_damage
!
!> \author Tong Zhang
!> \date   May. 2019
!> \details calulate the damage tracer for floating ice shelves (damage at grounded ice is constrained as 0)
!> Bassis, Jeremy N., and Y. Ma. "Evolution of basal crevasses links ice shelf stability to ocean forcing."
!> Earth and Planetary Science Letters 409 (2015): 203-211.

!-----------------------------------------------------------------------

   subroutine li_calculate_damage(domain, rkTend, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag
      real (kind=RKIND), dimension(:,:), intent(out) :: rkTend
      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: thermalPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool

      real(kind=RKIND), pointer :: config_damage_preserve_threshold
      character (len=StrKIND), pointer :: config_damage_gl_setting
      real(kind=RKIND), pointer :: config_default_flowParamA
      real(kind=RKIND), pointer :: config_flowLawExponent
      logical, pointer :: config_print_calving_info

      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: eMax
      real (kind=RKIND), dimension(:), pointer :: tauMax, tauMin
      real (kind=RKIND), dimension(:), pointer :: floatingBasalMassBal
      real (kind=RKIND), dimension(:), pointer :: damage
      real (kind=RKIND), dimension(:), pointer :: ddamagedt
      real (kind=RKIND), dimension(:), pointer :: s0, nstar
      real (kind=RKIND), dimension(:), pointer :: damageSource
      real (kind=RKIND), dimension(:), pointer :: damageNye
      real (kind=RKIND), dimension(:), pointer :: damageMax
      real (kind=RKIND), pointer ::  &
           config_ice_density,            & ! ice density
           config_ocean_density             ! ocean density

      real (kind=RKIND), dimension(:), pointer :: principalStrainRateRatio

      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltat !< time step (s)
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, pointer :: nCells, nVertLevels
      integer :: iCell, jCell, iNeighbor, n_damage_downstream, k
      real(kind=RKIND) :: damage_downstream
      real(kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY
      real(kind=RKIND), dimension(6) :: localMinInfo, localMaxInfo, globalMinInfo, globalMaxInfo

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', config_flowLawExponent)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_damage_preserve_threshold', config_damage_preserve_threshold)
      call mpas_pool_get_config(liConfigs, 'config_damage_gl_setting', config_damage_gl_setting)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get fields
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'damage', damage)
         call mpas_pool_get_array(geometryPool, 'ddamagedt', ddamagedt)
         call mpas_pool_get_array(geometryPool, 's0', s0)
         call mpas_pool_get_array(geometryPool, 'nstar', nstar)
         call mpas_pool_get_array(geometryPool, 'damageSource', damageSource)
         call mpas_pool_get_array(geometryPool, 'damageNye', damageNye)
         call mpas_pool_get_array(geometryPool, 'damageMax', damageMax)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)

         call mpas_pool_get_array(velocityPool, 'eMax', eMax)
         call mpas_pool_get_array(velocityPool, 'tauMax', tauMax)
         call mpas_pool_get_array(velocityPool, 'tauMin', tauMin)
         call mpas_pool_get_array(velocityPool, 'principalStrainRateRatio', principalStrainRateRatio)

         where ( thickness == 0.0_RKIND .or. tauMax == 0.0_RKIND )
             principalStrainRateRatio = 0.0_RKIND
             s0 = 0.0_RKIND
             nstar = 0.0_RKIND
         elsewhere
             ! Compute the hydrostatic to tensile stress ratio (equation 24 from Bassis and Ma, 2015, EPSL, 409, C)
             ! In Bassis & Ma this is beta, defined as ratio of principal strain rates.  Here we calculate as ratio of
             ! deviatoric stresses, which is equivalent.
             principalStrainRateRatio = tauMin / tauMax

             s0 = &
             config_ice_density*(config_ocean_density-config_ice_density)*gravity*thickness/(2.0_RKIND*tauMax*config_ocean_density)

             ! Compute the effective flow law exponent (equation 11 from Bassis and Ma, 2015, EPSL, 409, C)
             nstar(:) = 4.0_RKIND * config_flowLawExponent * &
                     (1.0_RKIND + principalStrainRateRatio(:) + principalStrainRateRatio(:)**2) / &
                     (4.0_RKIND * (1.0_RKIND + principalStrainRateRatio(:) + principalStrainRateRatio(:)**2) + &
                      3.0_RKIND * (config_flowLawExponent - 1.0_RKIND )* principalStrainRateRatio(:)**2)
         endwhere


         do iCell = 1, nCells
         if (thickness(iCell) == 0.0_RKIND) then
             damageSource(iCell) = 0.0_RKIND
         else
             damageSource(iCell) = nstar(iCell) * (1.0_RKIND - s0(iCell)) * eMax(iCell) - floatingBasalMassBal(iCell) &
                     / config_ice_density / thickness(iCell)
             ! RHS of equation 26 from Bassis and Ma (2015, EPSL, 409, C). Note that the last set of terms after the minus sign,
             ! mdot / rho_i / H, differs from the mdot/H in the paper because our mdot has units of kg/m^2/s and to convert that
             ! back to the right units (1/sec) you need to multiply by a factor of m^2/kg (1/rho_i * 1/H gives that factor).
         endif
         enddo

         !! Commented out lines immediately below here are for use with manufacured experiments
         !ddamagedt(:) = damageSource(:)*damage(:) - 0.44/2.0 * (uReconstructX(1,:)/1000.0*(1+epsmax(:)*seconds)* &
         ! exp(-uReconstructX(1,:)*seconds/1000.0) + damageSource(:) * (1+exp(-uReconstructX(1,:)*seconds/1000.0)))
         !ddamagedt(:) = damageSource(:)*damage(:) - 0.44/2.0 * (uReconstructX(1,:)/1000.0*exp(-uReconstructX(1,:)*seconds/1000.0)+&
         !    damageSource(:) * (1+exp(-uReconstructX(1,:)*seconds/1000.0)))

         ddamagedt(:) = damageSource(:) * damage(:)
         ! Update damage tendency for RK time integration
         do k = 1, nVertLevels
            rkTend(k,:) = rkTend(k,:) + ddamagedt(:)
         enddo

         damage(:) = damage(:) + ddamagedt(:) * deltat

         do iCell = 1, nCells
             if (thickness(iCell) == 0.0_RKIND) then
                 damageNye(iCell) = 0.0_RKIND
             else
                 damageNye(iCell) = (2.0_RKIND + principalStrainRateRatio(iCell)) * tauMax(iCell) / &
                         (gravity * (config_ocean_density - config_ice_density) * thickness(iCell))
             endif
         enddo

         do iCell = 1, nCells
            if (damage(iCell) > config_damage_preserve_threshold) then
                damageMax(iCell) = damage(iCell)
            end if
         end do
         ! save the damageMax value for restoring later if that option is chosen (i.e., do NOT allowing for healing to occur)

         do iCell = 1, nCells
            if (damage(iCell) < 0.0_RKIND) then
                damage(iCell) = 0.0_RKIND
            end if
            if (damage(iCell) < damageNye(iCell)) then
                damage(iCell) = damageNye(iCell)
            end if
            if (damage(iCell) > 1.0_RKIND) then
                damage(iCell) = 1.0_RKIND
            end if
         end do
         ! damage is always larger than the Nye value for initialization of damage evolution.

         do iCell = 1, nCells
            if ((li_mask_is_grounded_ice(cellMask(iCell))) .or. (thickness(iCell) .eq. 0.0_RKIND)) then
                damage(iCell) = 0.0_RKIND
            end if
         end do
         ! the damage of grounded ice is kept as 0, as the strain rate calculation is only valid for ice shelf

        ! Options for initializing damage where ice goes afloat:
        if (trim(config_damage_gl_setting) == 'extrapolate') then
           ! set the damage for a cell at the GL to the average damage value for its neighbor cells downstream on the ice shelf
           do iCell = 1, nCells
              if (li_mask_is_grounding_line(cellMask(iCell))) then

                  damage_downstream = 0.0_RKIND
                  n_damage_downstream = 0
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                      jCell = cellsOnCell(iNeighbor, iCell)
                      if (li_mask_is_floating_ice(cellMask(jCell))) then

                          damage_downstream = damage_downstream + damage(jCell)
                          n_damage_downstream =  n_damage_downstream + 1
                      end if
                  end do

                  if (n_damage_downstream == 0) then
                      damage(iCell) = 0.0_RKIND
                  else
                      damage(iCell) = damage_downstream / n_damage_downstream
                  end if

              end if
           end do
        elseif (trim(config_damage_gl_setting) == "nye") then
            ! set the damage at the first floating cells to the Nye value
            do iCell = 1, nCells
                if (li_mask_is_grounding_line(cellMask(iCell))) then
                    do iNeighbor = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(iNeighbor, iCell)
                        if (li_mask_is_floating_ice(cellMask(jCell))) then
                            damage(jCell) = damageNye(jCell)
                        end if
                    end do
                end if
            end do
        else
            call mpas_log_write("Unknown value specified for config_damage_gl_setting!", MPAS_LOG_ERR)
            err = ior(err, 1)
        endif

        if (config_print_calving_info) then
           ! End of routine accounting/reporting
           localMinInfo(1) = minval(damageSource)
           localMaxInfo(1) = maxval(damageSource)
           localMinInfo(2) = minval(s0)
           localMaxInfo(2) = maxval(s0)
           localMinInfo(3) = minval(nstar)
           localMaxInfo(3) = maxval(nstar)
           localMinInfo(4) = minval(principalStrainRateRatio)
           localMaxInfo(4) = maxval(principalStrainRateRatio)
           localMinInfo(5) = minval(tauMax)
           localMaxInfo(5) = maxval(tauMax)
           localMinInfo(6) = minval(damage)
           localMaxInfo(6) = maxval(damage)
           ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
           call mpas_dmpar_min_real_array(domain % dminfo, 6, localMinInfo, globalMinInfo)
           call mpas_dmpar_max_real_array(domain % dminfo, 6, localMaxInfo, globalMaxInfo)

           call mpas_log_write("damageSource value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(1), globalMaxInfo(1)/))
           call mpas_log_write("s0 value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(2), globalMaxInfo(2)/))
           call mpas_log_write("nstar value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(3), globalMaxInfo(3)/))
           call mpas_log_write("principalStrainRateRatio value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(4), globalMaxInfo(4)/))
           call mpas_log_write("tauMax value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(5), globalMaxInfo(5)/))
           call mpas_log_write("damage value range: Min=$r, Max=$r", &
                        realArgs=(/globalMinInfo(6), globalMaxInfo(6)/))
        end if
        ! temporary message output in the log file


         block => block % next

      enddo

   end subroutine li_calculate_damage

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine li_finalize_damage_after_advection
!
!> \author Tong Zhang
!> \date   Oct. 2019
!> \details This routine finalizes the damage calculation after advection by several ways: 1) set the damage at grounding line
!>          using an extrapolation method or an Nye method.
!>          2) preserve the damage to the original value before advection if it exceeds some large number
!>          (config_damage_preserve_threshold). For example, if config_damage_preserve_threshold = 0.5,
!>          then the localtions with damage > 0.5 on the ice shelf will not heal.
!>          Thus, if we turn on this function, we do not allow the damage to heal if damage > config_damage_preserve_threshold.
!>           3) coupling rheology to damage if the option config_damage_rheology_coupling is set to True in namelist
!-----------------------------------------------------------------------

   subroutine li_finalize_damage_after_advection(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool
      real(kind=RKIND), pointer :: config_damage_stiffness_min
      logical, pointer :: config_damage_rheology_coupling
      logical, pointer :: config_preserve_damage
      logical, pointer :: config_print_calving_info
      character (len=StrKIND), pointer :: config_damage_gl_setting
      real (kind=RKIND), dimension(:), pointer :: damage
      real (kind=RKIND), dimension(:), pointer :: damageMax
      real (kind=RKIND), dimension(:), pointer :: damageNye
      real (kind=RKIND), dimension(:), pointer :: stiffnessFactor


      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, pointer :: nCells
      integer :: iCell, jCell, iNeighbor, n_damage_downstream
      real(kind=RKIND) :: damage_downstream
      integer :: err_tmp

      err = 0


      call mpas_pool_get_config(liConfigs, 'config_damage_stiffness_min', config_damage_stiffness_min)
      call mpas_pool_get_config(liConfigs, 'config_damage_rheology_coupling', config_damage_rheology_coupling)
      call mpas_pool_get_config(liConfigs, 'config_damage_gl_setting', config_damage_gl_setting)
      call mpas_pool_get_config(liConfigs, 'config_preserve_damage', config_preserve_damage)

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         ! get fields
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'damage', damage)
         call mpas_pool_get_array(geometryPool, 'damageMax', damageMax)
         call mpas_pool_get_array(geometryPool, 'damageNye', damageNye)
         call mpas_pool_get_array(velocityPool, 'stiffnessFactor', stiffnessFactor)

         ! make sure masks are up to date.  May not be necessary, but safer to do anyway.
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         if (config_preserve_damage) then
             do iCell = 1, nCells
                if (damage(iCell) < damageMax(iCell)) then
                    damage(iCell) = damageMax(iCell)
                end if
             end do
         endif
         ! put the damageMax value back to preserve the damage value (no heal)

         where (damage < 0.0_RKIND)
             damage = 0.0_RKIND
         end where

         where (damage > 1.0_RKIND)
             damage = 1.0_RKIND
         end where


         do iCell = 1, nCells
            if (li_mask_is_grounded_ice(cellMask(iCell)) .or. .not. li_mask_is_ice(cellMask(iCell))) then
                damage(iCell) = 0.0_RKIND
            end if
         end do

         if (trim(config_damage_gl_setting) == 'extrapolate') then
             do iCell = 1, nCells
                if (li_mask_is_grounding_line(cellMask(iCell))) then
                    damage_downstream = 0.0_RKIND
                    n_damage_downstream = 0
                    do iNeighbor = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(iNeighbor, iCell)
                        if (li_mask_is_floating_ice(cellMask(jCell))) then

                            damage_downstream = damage_downstream + damage(jCell)
                            n_damage_downstream =  n_damage_downstream + 1
                        end if
                    end do

                    damage(iCell) = damage_downstream

                    if (n_damage_downstream == 0) then
                        damage(iCell) = 0.0_RKIND
                    else
                        damage(iCell) = damage_downstream / n_damage_downstream
                    end if

                end if
            end do
            !set the damage at GL to the mean damage value for its neighboring floating cells

        elseif (trim(config_damage_gl_setting) == 'nye') then
             do iCell = 1, nCells
                if (li_mask_is_grounding_line(cellMask(iCell))) then
                    do iNeighbor = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(iNeighbor, iCell)
                        if (li_mask_is_floating_ice(cellMask(jCell))) then
                            damage(jCell) = damageNye(jCell)
                        end if
                    end do
                end if
            end do
            ! always set the damage at the first floating cells to the Nye value
        endif

        ! Limit range of values again after doing extrapolation
        where (damage < 0.0_RKIND)
            damage = 0.0_RKIND
        end where

        where (damage > 1.0_RKIND)
            damage = 1.0_RKIND
        end where

        if (config_damage_rheology_coupling) then
            do iCell = 1, nCells
                if (li_mask_is_floating_ice(cellMask(iCell))) then
                    stiffnessFactor(iCell) = 1.0_RKIND - damage(iCell)
                    if (stiffnessFactor(iCell) < config_damage_stiffness_min) then
                        stiffnessFactor(iCell) = config_damage_stiffness_min
                    end if
                end if
            end do
        end if


        block => block % next

      enddo

   end subroutine li_finalize_damage_after_advection


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine apply_calving_damage_threshold
!
!> \brief Calve any ice that is damaged beyond a specified threshold
!> \author Tong Zhang, Matt Hoffman, Trevor Hillebrand
!> \date   Nov. 2020, March 2021, April 2021
!> \details This routine specified floating ice to be calved wherever the damage
!> value exceeds a specified threshold, assuming the ice is connected to the calving
!> front.
!-----------------------------------------------------------------------
   subroutine apply_calving_damage_threshold(meshPool, geometryPool, scratchPool, domain, err)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), pointer, intent(in) :: scratchPool !< Input: scratch pool

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain          !< Input/Output: domain object
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input: geometry pool

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (field1dInteger), pointer :: seedMaskField
      type (field1dInteger), pointer :: growMaskField
      integer, dimension(:), pointer :: seedMask, growMask !masks to pass to flood-fill routine
      real (kind=RKIND), dimension(:), pointer :: calvingThickness !< Output: the applied calving rate as a thickness
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold !< Output: the applied calving thickness from threshold processes
      real (kind=RKIND), dimension(:), pointer :: thickness, damage
      real(kind=RKIND), pointer :: config_damage_calving_threshold
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: cellMask
      integer, pointer :: nCells
      integer :: iCell, iNeighbor, jCell

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_damage_calving_threshold', config_damage_calving_threshold)

      ! get fields
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
      call mpas_pool_get_array(geometryPool, 'damage', damage)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      ! Allocate scratch fields for flood-fill
      call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
      call mpas_allocate_scratch_field(seedMaskField, single_block_in = .true.)
      seedMask => seedMaskField % array
      seedMask(:) = 0

      call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
      call mpas_allocate_scratch_field(growMaskField, single_block_in = .true.)
      growMask => growMaskField % array
      growMask(:) = 0

      ! define seed and grow masks for flood fill.
      where ( li_mask_is_dynamic_margin(cellMask) .and. (damage .ge. config_damage_calving_threshold) )
              seedMask = 1
      end where

      where ( seedMask == 0 .and. (damage .ge. config_damage_calving_threshold) )
              growMask = 1
      end where

      call li_flood_fill(seedMask, growMask, domain)

      ! Remove ice from flood-filled mask, and any neighboring non-dynamic ice
      do iCell = 1, nCells
         if ( seedMask(iCell) == 1 ) then
            ! First remove ice from "thin" neighbors
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_floating_ice(cellMask(jCell)) .and. .not. li_mask_is_dynamic_ice(cellMask(jCell))) then
                  ! this is a thin neighbor - remove the whole cell volume
                  calvingThickness(jCell) = thickness(jCell)
                  calvingThicknessFromThreshold(jCell) = calvingThicknessFromThreshold(jCell) + thickness(jCell)
               endif
            enddo

            calvingThickness(iCell) = thickness(iCell)
            calvingThicknessFromThreshold(iCell) = calvingThicknessFromThreshold(iCell) + thickness(iCell)
         endif ! if cell is calving margin

      enddo ! cell loop
      ! deallocate scratch fields used for flood fill
      call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.true.)
      call mpas_deallocate_scratch_field(growMaskField, single_block_in=.true.)

   end subroutine apply_calving_damage_threshold


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine mask_calving
!
!> \brief Calve ice based on a mask that is input
!> \author Matthew Hoffman
!> \date   July 2019
!> \details  This routine uses the field calvingMask to eliminate floating
!> ice.  Locations marked with 0 will not be calved.  Locations with all
!> all other values will be calved.  Note that time-varying values
!> can be input.
!-----------------------------------------------------------------------
   subroutine mask_calving(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: &
         domain          !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold
      integer, dimension(:), pointer :: calvingMask
      integer, dimension(:), pointer :: cellMask
      integer, pointer :: nCells, nCellsSolve
      integer :: iCell
      integer :: localMaskCellCount, globalMaskCellCount
      integer :: err_tmp

      err = 0

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         ! get fields
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
         call mpas_pool_get_array(geometryPool, 'calvingMask', calvingMask)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         localMaskCellCount = 0

         ! === apply calving ===
         do iCell = 1, nCells
            if (li_mask_is_floating_ice(cellMask(iCell)) .and. (calvingMask(iCell) >= 1)) then
               !call mpas_log_write("Found masked cell at $i", intArgs=(/iCell/))
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)
               calvingThicknessFromThreshold(iCell) = calvingThicknessFromThreshold(iCell) + thickness(iCell)
               thickness(iCell) = 0.0_RKIND
               if (iCell <= nCellsSolve) localMaskCellCount = localMaskCellCount + 1
            endif
         enddo

         call mpas_dmpar_sum_int(domain % dminfo, localMaskCellCount, globalMaskCellCount)
         call mpas_log_write("Mask calving applied.  Removed $i floating cells with mask.", intArgs=(/globalMaskCellCount/))

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         call remove_small_islands(meshPool, geometryPool)

         block => block % next
      enddo

   end subroutine mask_calving

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine calculate_calving_front_mask
!
!> \brief  Calculate mask indicating position of effective calving front
!> \author Matthew Hoffman
!> \date   Feb. 2018
!> \details Mmake mask for effective calving front.
!> This is last dynamic floating cell, but also make sure it has a neighbor that is open ocean or thin floating ice.
!-----------------------------------------------------------------------
   subroutine calculate_calving_front_mask(meshPool, geometryPool, calvingFrontMask)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), intent(in) :: geometryPool !< Input: geometry pool

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, dimension(:), intent(out) :: calvingFrontMask !< Output: calving front mask

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer, pointer :: nCells
      integer :: iCell, iNeighbor, jCell, jNeighbor, kCell
      logical :: oceanNeighborCell, oceanNeighborEdge
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell ! list of edges that neighbor each cell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), pointer :: config_sea_level

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      ! get fields
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

      !calvingFrontMaskEdge = 0 !initialize
      calvingFrontMask = 0 !initialize

      do iCell = 1, nCells
         if ( (li_mask_is_floating_ice(cellMask(iCell))) .and. (li_mask_is_dynamic_margin(cellMask(iCell))) ) then
            oceanNeighborCell = .false.
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               oceanNeighborEdge = .false.
               if (li_mask_is_floating_ice(cellMask(jCell)) .and. .not. li_mask_is_dynamic_ice(cellMask(jCell))) then
                  ! make sure this neighbor is adjacent to open ocean (and not thin floating ice up against the coast)
                  do jNeighbor = 1, nEdgesOnCell(jCell)
                     kCell = cellsOnCell(jNeighbor, jCell)
                     if (.not. li_mask_is_ice(cellMask(kCell)) .and. bedTopography(kCell) < config_sea_level) then
                        oceanNeighborEdge = .true. ! iCell neighbors thin ice that in turn neighbors open ocean
                     endif
                  enddo
               endif
               if (.not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level) then
                  oceanNeighborEdge = .true. ! this is an open ocean neighbor
               endif
               !if (oceanNeighborEdge) then
               !   calvingFrontMaskEdge(edgesOnCell(iNeighbor, iCell)) = 1
               !endif
               oceanNeighborCell = (oceanNeighborCell .or. oceanNeighborEdge)
            enddo
            if (oceanNeighborCell) then
               calvingFrontMask(iCell) = 1
            endif
         endif
      enddo

   end subroutine calculate_calving_front_mask


   subroutine remove_icebergs(domain)
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain !< Input/Output: domain object


      ! Local variables
      logical, pointer :: config_remove_icebergs
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool

      type (field1dInteger), pointer :: seedMaskField
      type (field1dInteger), pointer :: growMaskField
      integer, dimension(:), pointer :: seedMask, growMask !masks to pass to flood-fill routine
      !integer, dimension(:), allocatable :: seedMaskOld !mask of where icebergs are removed, for debugging

      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: calvingThicknessFromThreshold
      real (kind=RKIND), dimension(:), pointer :: thickness
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell

      integer, pointer :: nCells, nCellsSolve
      integer :: iCell, jCell, n
      integer :: newMaskCountLocal, newMaskCountLocalAccum, newMaskCountGlobal
      integer :: err_tmp, err
      integer :: globalLoopCount, localLoopCount
      integer :: localIcebergCellCount, globalIcebergCellCount

      err = 0


      call mpas_pool_get_config(liConfigs, 'config_remove_icebergs', config_remove_icebergs)
      if (.not. config_remove_icebergs) then
              return  ! skip this entire routine if disabled
      endif

      call mpas_timer_start("iceberg detection")
      call mpas_log_write("Iceberg-detection flood-fill begin.")
      
      ! Allocate needed scratch fields
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
      call mpas_allocate_scratch_field(seedMaskField, single_block_in = .false.)
      seedMask => seedMaskField % array

      call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
      call mpas_allocate_scratch_field(growMaskField, single_block_in = .false.)
      growMask => growMaskField % array

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      
      seedMask(:) = 0
      growMask(:) = 0

      call mpas_log_write("Iceberg-detection flood-fill: allocated.")
      ! First mark grounded ice to initialize flood fill mask
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(geometryPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

         ! make sure masks are up to date.  May not be necessary, but safer to
         ! do anyway.
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         call mpas_log_write("Iceberg-detection flood-fill: updated masks.")
         newMaskCountLocal = 0
         do iCell = 1, nCellsSolve
            if (li_mask_is_grounded_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell))) then
               seedMask(iCell) = 1
               newMaskCountLocal = newMaskCountLocal + 1
            endif
         enddo
         call mpas_log_write("Initialized $i cells to local mask", intArgs=(/newMaskCountLocal/))

         block => block % next
      end do

      call mpas_log_write("Iceberg-detection flood-fill initialization complete.")
      
      where ( (seedMask == 0) .and. li_mask_is_floating_ice(cellMask(:)) .and. li_mask_is_dynamic_ice(cellMask(:)) )
             growMask = 1
      endwhere
      call li_flood_fill(seedMask, growMask, domain)

      ! Add floating non-dynamic fringe, but exclude dynamic ice isolated by
      ! non-dynamic ice, which can cause velocity solver to fail to converge.
      ! This is a bit of an expensive solution to this problem
      ! that may only occur once in a while, so we might need to revisit this.
      
      where ( li_mask_is_floating_ice(cellMask(:)) .and. .not. li_mask_is_dynamic_ice(cellMask(:)) )
             growMask = 1
      elsewhere
             growMask = 0
      endwhere
      call li_flood_fill(seedMask, growMask, domain)
      
      ! Now remove any ice that was not flood-filled - these are icebergs
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingThicknessFromThreshold', calvingThicknessFromThreshold)
         call mpas_pool_get_dimension(geometryPool, 'nCells', nCells)
         call mpas_pool_get_dimension(geometryPool, 'nCellsSolve', nCellsSolve)
         !allocate(seedMaskOld(nCells+1)) ! debug: make this a mask of where icebergs were removed
         !seedMaskOld(:) = 0

         localIcebergCellCount = 0
         do iCell = 1, nCellsSolve
            if (seedMask(iCell) == 0 .and. li_mask_is_floating_ice(cellMask(iCell))) then
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)  ! remove any remaining ice here
               calvingThicknessFromThreshold(iCell) = calvingThicknessFromThreshold(iCell) + thickness(iCell)  ! remove any remaining ice here
               thickness(iCell) = 0.0_RKIND
               localIcebergCellCount = localIcebergCellCount + 1
               !seedMaskOld(iCell) = 1 ! debug: make this a mask of where icebergs were removed
            endif
         enddo

         block => block % next
      end do

      ! Calculate total iceberg cells removed (not needed but useful for log - note requires extra reduce!)
      call mpas_dmpar_sum_int(domain % dminfo, localIcebergCellCount, globalIcebergCellCount)

      ! update halo - probably not needed?
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_timer_stop("halo updates")

      ! clean up
      !deallocate(seedMaskOld) !un-comment for debugging
      call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.false.)
      call mpas_deallocate_scratch_field(growMaskField, single_block_in=.false.)
      call mpas_log_write("Iceberg-detection flood-fill complete. Removed $i iceberg cells.", intArgs=(/globalIcebergCellCount/))
      call mpas_timer_stop("iceberg detection")
   end subroutine remove_icebergs

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine li_flood_fill
!
!> \brief  Flood fill routine to differentiate ice sheet from icebergs, critically damaged regions, etc. 
!> The calling routine defines the seedMask to specify the initial masked region, and the growMask 
!> to specify the potential region to be filled by the flood-fill routine. The flood-fill routine
!> then adds the contiguous areas of the growMask to the seedMask, which is passed back to the 
!> calling routine for application. Any integer fields in Registry can be used for seedMask and
!> growMask, but there are existing variables with those names that can be used.
!> \author Trevor Hillebrand and Matthew Hoffman
!> \date   March 2021
!> \details 
!-----------------------------------------------------------------------
   subroutine li_flood_fill(seedMask, growMask, domain)
      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain !< Input/Output: domain object
      integer, dimension(:), intent(inout) :: seedMask
      integer, dimension(:), intent(in) :: growMask
      ! Local variables
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool

      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:), allocatable :: seedMaskOld

      integer, pointer :: nCells, nCellsSolve
      integer :: iCell, jCell, n
      integer :: newMaskCountLocal, newMaskCountLocalAccum, newMaskCountGlobal
      integer :: err_tmp, err
      integer :: globalLoopCount, localLoopCount

      err = 0

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      allocate(seedMaskOld(nCells+1))
      seedMaskOld(:) = 0

      !call mpas_log_write("Flood-fill: allocated.")

      ! First mark grounded ice to initialize flood fill mask
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(geometryPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

         ! make sure masks are up to date.  May not be necessary, but safer to
         ! do anyway.
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         !call mpas_log_write("Flood-fill: updated masks.")
         newMaskCountLocal = sum(seedMask)
         !call mpas_log_write("Initialized $i cells to local mask", intArgs=(/newMaskCountLocal/))

         block => block % next
      end do

      !call mpas_log_write("Flood-fill initialization complete.")

      ! Outer loop over processors (should also have a loop over blocks)
      ! Inner loop over cells on that processor

      ! Initialize global mask count
      call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocal, newMaskCountGlobal)
      globalLoopCount = 0

      call mpas_log_write("Flood fill initialized $i cells to global seedMask", intArgs=(/newMaskCountGlobal/))

      do while (newMaskCountGlobal > 0)
         globalLoopCount = globalLoopCount + 1
         !call mpas_log_write("  Starting global processor loop $i", intArgs=(/globalLoopCount/))
         ! First Update halos
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'growMask')
         call mpas_dmpar_field_halo_exch(domain, 'seedMask')
         call mpas_timer_stop("halo updates")

         ! initialize counter of cells locally updated during this outer loop
         newMaskCountLocalAccum = 0

         ! Now update (advance) mask locally

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

            ! initialize local loop
            localLoopCount = 0
            newMaskCountLocal = 1  ! need to make sure we enter the loop
            do while (newMaskCountLocal > 0)
               localLoopCount = localLoopCount + 1
               !call mpas_log_write("    Starting local cell loop $i", intArgs=(/localLoopCount/))

               ! initialize
               newMaskCountLocal = 0
               seedMaskOld(:) = seedMask(:)

               do iCell = 1, nCellsSolve ! this gives owned cells only
                  if ( growMask(iCell)>0 ) then
                     ! If it has a marked neighbor, then add it to the mask
                     do n = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(n, iCell)
                        if ( (seedMaskOld(jCell) == 1) .and. (seedMask(iCell) .ne. 1) ) then
                           seedMask(iCell) = 1
                           newMaskCountLocal = newMaskCountLocal + 1
                           exit ! skip the rest of this do-loop - no need to check additional neighbors
                        endif
                     enddo
                  endif ! if not already marked
               enddo ! loop over cells

               ! Accumulate cells added locally until we do the next global
               ! reduce
               newMaskCountLocalAccum = newMaskCountLocalAccum + newMaskCountLocal
               !call mpas_log_write("    Added $i new cells to local mask", intArgs=(/newMaskCountLocal/))
            enddo ! local mask loop

            block => block % next
         end do

         ! update count of cells added to mask globally
         call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocalAccum, newMaskCountGlobal)
         call mpas_log_write("  Added $i new cells to global mask", intArgs=(/newMaskCountGlobal/))

         if (globalLoopCount>200) then
            call mpas_log_write("Too many global loops!", MPAS_LOG_ERR)
         endif
      end do ! global loop
      deallocate(seedMaskOld)

      call mpas_log_write("Flood fill complete.")

   end subroutine li_flood_fill

end module li_calving


