! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   public :: li_calve_ice, li_restore_calving_front, calculate_damage, restore_damage

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
      character (len=StrKIND), pointer :: config_calving
      logical, pointer :: config_print_calving_info, config_data_calving
      real(kind=RKIND), pointer :: config_calving_timescale

      integer, pointer :: nCells

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

      type (field1dReal), pointer :: originalThicknessField

      real (kind=RKIND), dimension(:), pointer :: originalThickness

      integer :: iCell

      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)
      call mpas_pool_get_config(liConfigs, 'config_calving_timescale', config_calving_timescale)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_data_calving', config_data_calving)

      if (trim(config_calving) == 'none') then
         return ! do nothing
      endif

      ! Get deltat from first block (same on all blocks)
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)

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
      if (trim(config_calving) == 'thickness_threshold') then

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

         call damage_calving(domain, err_tmp)
         err = ior(err, err_tmp)

      elseif (trim(config_calving) == 'mask') then

         call mask_calving(domain, err_tmp)
         err = ior(err, err_tmp)

      else

         call mpas_log_write("Invalid option for config_calving specified: " // trim(config_calving), MPAS_LOG_ERR)
         err = 1

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

      integer, pointer :: nCellsSolve, nVertLevels

      logical, pointer :: &
           config_print_calving_info

      real (kind=RKIND), pointer ::   &
           config_sea_level,          &
           config_dynamic_thickness

      integer, dimension(:), pointer :: &
           cellMask            ! bit mask describing whether ice is floating, dynamically active, etc.

      real(kind=RKIND), dimension(:), pointer :: &
           layerCenterSigma    ! vertical sigma coordinate at layer midpoints

      ! geometry and calving fields
      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography,    & ! elevation of the bed
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

      integer :: iCell, err_tmp

      !WHL - debug
      logical, parameter :: circular_shelf = .false.
      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow
      integer :: k


      ! first remove any icebergs - do it first so restore-calving can put back thin ice in those places
      call remove_icebergs(domain)

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
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         ! get required fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

         ! get required fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'restoreThickness', restoreThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         ! get required fields from the thermal pool
         call mpas_pool_get_array(thermalPool, 'temperature', temperature)
         call mpas_pool_get_array(thermalPool, 'waterFrac', waterFrac)
         call mpas_pool_get_array(thermalPool, 'surfaceAirTemperature', surfaceAirTemperature)
         call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
         call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)

         ! get config variables
         call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)

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
         calvingThickness = 0.0_RKIND
         restoreThickness = 0.0_RKIND

         ! loop over locally owned cells
         do iCell = 1, nCellsSolve

            if (bedTopography(iCell) < config_sea_level) then

               ! The bed is below sea level; test for calving-front advance and retreat.

               if (li_mask_is_initial_ice(cellMask(iCell)) .and. thickness(iCell) < restoreThicknessMin) then

                  ! Ice was present in this cell initially, but now is either very thin or absent;
                  !  reset the thickness to restoreThicknessMin.
                  ! Note: Mass is not conserved.
                  !       Save the difference (restoreThicknessMin - thickness) so as to keep track of energy non-conservation.

                  if (config_print_calving_info) then
                     call mpas_log_write('Restore ice: iCell=$i, thickness=$r', intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
                  endif

                  restoreThickness(iCell) = restoreThicknessMin - thickness(iCell)
                  thickness(iCell) = restoreThicknessMin

                  ! Initial a linear temperature profile in the column
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

               elseif (.not.li_mask_is_initial_ice(cellMask(iCell)) .and. thickness(iCell) > 0.0_RKIND) then

                  ! This cell was initially ice-free but now has ice.
                  ! Remove the ice and add it to calvingThickness.

                  if (config_print_calving_info) then
                     call mpas_log_write('Remove ice:  iCell=$i, thickness=$r', intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
                  endif

                  calvingThickness(iCell) = thickness(iCell)
                  thickness(iCell) = 0.0_RKIND

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
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         calvingThickness = 0.0_RKIND

         ! Note: The floating_ice mask includes all floating ice, both inactive and active
         where (li_mask_is_floating_ice(cellMask))
            calvingThickness = thickness * calvingFraction
         endwhere

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
               thickness(iCell) = 0.0_RKIND
            elseif (nIceNeighbors == 1) then
               ! check if this neighbor has any additional neighbors with ice
               nIceNeighbors2 = 0
               nOpenOceanNeighbors = 0
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
                  thickness(iCell) = 0.0_RKIND
                  calvingThickness(neighborWithIce) = calvingThickness(neighborWithIce) + thickness(neighborWithIce)
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
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         calvingThickness = 0.0_RKIND

         where ( (li_mask_is_floating_ice(cellMask)) .and. (bedTopography < config_calving_topography + config_sea_level) )
            calvingThickness = thickness * calvingFraction
         endwhere

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
      logical, pointer :: config_print_calving_info
      real(kind=RKIND), pointer :: config_calving_thickness
      real (kind=RKIND), dimension(:), pointer :: eigencalvingParameter
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
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

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_scalar_value', &
               config_calving_eigencalving_parameter_scalar_value)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_source', &
              config_calving_eigencalving_parameter_source)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)

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

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         calvingVelocity(:) = 0.0_RKIND
         ! First calculate the front retreat rate (Levermann eq. 1)
         !calvingVelocity(:) = eigencalvingParameter(:) * max(0.0_RKIND, eMax(:)) * max(0.0_RKIND, eMin(:)) & ! m/s
         !      * real(li_mask_is_floating_ice_int(cellMask(:)), kind=RKIND)
                  ! calculate only for floating ice - map of "potential" calving rate
         calvingVelocity(:) = config_calving_eigencalving_parameter_scalar_value
        

         ! Convert calvingVelocity to calvingThickness
         call apply_calving_velocity(meshPool, geometryPool, velocityPool, domain, err)


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
!    routine damagecalving
!
!> \brief Calve ice from the calving front based on the Bassis-Ma (2015) damage theory 
!> \author Tong Zhang
!> \date   May. 2019
!> \details Calving using damage model
!> Bassis, Jeremy N., and Y. Ma. "Evolution of basal crevasses links ice shelf stability to ocean forcing." 
!> Earth and Planetary Science Letters 409 (2015): 203-211.
! TODO need to properly connect damage tracer to calving flux

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
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin
      real (kind=RKIND), dimension(:), pointer :: damage
      integer, dimension(:), pointer :: cellMask
      type (field1dInteger), pointer :: calvingFrontMaskField
      integer, dimension(:), pointer :: calvingFrontMask
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

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_damagecalvingParameter', config_damagecalvingParameter)

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

         call mpas_pool_get_field(scratchPool, 'iceCellMask2',  calvingFrontMaskField)
         call mpas_allocate_scratch_field(calvingFrontMaskField, .true.)
         calvingFrontMask => calvingFrontMaskField % array

         !call calculate_calving_front_mask(meshPool, geometryPool, calvingFrontMask)

         calvingVelocity(:) = 0.0_RKIND
         ! First calculate the front retreat rate (Levermann eq. 1)
         calvingVelocity(:) = config_damagecalvingParameter * damage(:) & ! m/s
               * real(li_mask_is_floating_ice_int(cellMask(:)), kind=RKIND)
                  ! calculate only for floating ice - map of "potential" calving rate

         !call remove_calving_ice(meshPool, geometryPool, velocityPool, scratchPool, err)
         call apply_calving_velocity(meshPool, geometryPool, velocityPool, domain, err)

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

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
!    routine calculate_damage
!
!> \author Tong Zhang
!> \date   May. 2019
!> \details calulate the damage tracer 
!> Bassis, Jeremy N., and Y. Ma. "Evolution of basal crevasses links ice shelf stability to ocean forcing." 
!> Earth and Planetary Science Letters 409 (2015): 203-211.

!-----------------------------------------------------------------------

   subroutine calculate_damage(domain, err)

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
      type (mpas_pool_type), pointer :: thermalPool 
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: scratchPool
      real(kind=RKIND), pointer :: config_calving_eigencalving_parameter_scalar_value
      real(kind=RKIND), pointer :: config_damage_restore_threshold, config_damage_calving_threshold
      logical, pointer :: config_use_constant_A
      real(kind=RKIND), pointer :: config_default_flowParamA 
      character (len=StrKIND), pointer :: config_calving_eigencalving_parameter_source
      logical, pointer :: config_print_calving_info
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: lowerSurface
      real (kind=RKIND), dimension(:), pointer :: eigencalvingParameter
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), dimension(:), pointer :: requiredCalvingVolumeRate
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin, exx, exy, eyx, eyy
      real (kind=RKIND), dimension(:,:), pointer :: flowParamA
      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), dimension(:), pointer :: floatingBasalMassBal
      real (kind=RKIND), dimension(:), pointer :: damage
      real (kind=RKIND), dimension(:), pointer :: ddamagedt
      real (kind=RKIND), dimension(:), pointer :: szero
      real (kind=RKIND), dimension(:), pointer :: source
      real (kind=RKIND), dimension(:), pointer :: damage_nye
      real (kind=RKIND), dimension(:), pointer :: damage_max
      real (kind=RKIND), dimension(:), pointer :: stiffnessFactor
      real (kind=RKIND), pointer ::  &
           config_ice_density,            & ! ice density
           config_ocean_density             ! ocean density

      real (kind=RKIND), save :: rhoi    ! ice density (kg m^{-3}), copied from config_ice_density
      real (kind=RKIND), save :: rhoo    ! ocean density (kg m^{-3}), copied from config_ocean_density

      type (field1dReal), pointer :: meanFlowParamAField
      real, dimension(:), pointer :: meanFlowParamA
      type (field1dReal), pointer :: viscField
      real, dimension(:), pointer :: visc
      type (field1dReal), pointer :: sigmayxField
      real, dimension(:), pointer :: sigmayx
      type (field1dReal), pointer :: sigmaxyField
      real, dimension(:), pointer :: sigmaxy
      type (field1dReal), pointer :: sigmaxxField
      real, dimension(:), pointer :: sigmaxx
      type (field1dReal), pointer :: sigmayyField
      real, dimension(:), pointer :: sigmayy
      type (field1dReal), pointer :: sigmaminField
      real, dimension(:), pointer :: sigmamin
      type (field1dReal), pointer :: sigmamaxField
      real, dimension(:), pointer :: sigmamax
      type (field1dReal), pointer :: alphaField
      real, dimension(:), pointer :: alpha
      type (field1dReal), pointer :: epsmaxField
      real, dimension(:), pointer :: epsmax
      type (field1dReal), pointer :: nstarField
      real, dimension(:), pointer :: nstar
      !type (field1dReal), pointer :: szeroField
      !real, dimension(:), pointer :: szero
      !type (field1dReal), pointer :: sourceField
      !real, dimension(:), pointer :: source
      !type (field1dReal), pointer :: ddamagedtField
      !real, dimension(:), pointer :: ddamagedt
      type (field1dReal), pointer :: forceField
      real, dimension(:), pointer :: force
      type (field1dReal), pointer :: meanTField
      real, dimension(:), pointer :: meanT

      integer, dimension(:), pointer :: cellMask
      type (field1dInteger), pointer :: calvingFrontMaskField
      integer, dimension(:), pointer :: calvingFrontMask
      real (kind=RKIND), pointer :: deltat, daysSinceStart !< time step (s)
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, dvEdge
      integer, pointer :: nCells, nVertLevels
      integer :: iCell, jCell, iNeighbor, n_damage_downstream
      real(kind=RKIND) :: damage_downstream, damage_tmp, cellCalvingFrontLength, cellCalvingFrontHeight, seconds, R, A01, A02, Q1, Q2
      real(kind=RKIND), dimension(:,:), pointer :: uReconstructX, uReconstructY
      character (len=StrKIND), pointer :: xtime
      integer :: err_tmp, loc

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_scalar_value', &
               config_calving_eigencalving_parameter_scalar_value)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_source', config_calving_eigencalving_parameter_source)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_damage_restore_threshold', config_damage_restore_threshold)
      call mpas_pool_get_config(liConfigs, 'config_damage_calving_threshold', config_damage_calving_threshold)
      call mpas_pool_get_config(liConfigs, 'config_default_flowParamA', config_default_flowParamA)
      call mpas_pool_get_config(liConfigs, 'config_use_constant_A', config_use_constant_A)

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

         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'daysSinceStart', daysSinceStart)
         call mpas_pool_get_array(meshPool, 'xtime', xtime)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
         call mpas_pool_get_array(geometryPool, 'eigencalvingParameter', eigencalvingParameter)
         call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeRate', requiredCalvingVolumeRate)
         call mpas_pool_get_array(geometryPool, 'damage', damage)
         call mpas_pool_get_array(geometryPool, 'ddamagedt', ddamagedt)
         call mpas_pool_get_array(geometryPool, 'szero', szero)
         call mpas_pool_get_array(geometryPool, 'source', source)
         call mpas_pool_get_array(geometryPool, 'damage_nye', damage_nye)
         call mpas_pool_get_array(geometryPool, 'damage_max', damage_max)

         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
         call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
         call mpas_pool_get_array(velocityPool, 'eMax', eMax)
         call mpas_pool_get_array(velocityPool, 'eMin', eMin)
         call mpas_pool_get_array(velocityPool, 'exx', exx)
         call mpas_pool_get_array(velocityPool, 'exy', exy)
         call mpas_pool_get_array(velocityPool, 'eyx', eyx)
         call mpas_pool_get_array(velocityPool, 'eyy', eyy)
         call mpas_pool_get_array(velocityPool, 'flowParamA', flowParamA)
         call mpas_pool_get_array(velocityPool, 'stiffnessFactor', stiffnessFactor)

         call mpas_pool_get_array(thermalPool, 'temperature', temperature)

         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)

         call mpas_pool_get_field(scratchPool, 'iceCellMask2',  calvingFrontMaskField)
         call mpas_allocate_scratch_field(calvingFrontMaskField, .true.)
         calvingFrontMask => calvingFrontMaskField % array

         call mpas_pool_get_field(scratchPool, 'meanFlowParamAField',  meanFlowParamAField)
         call mpas_allocate_scratch_field(meanFlowParamAField, .true.)
         meanFlowParamA => meanFlowParamAField % array

         call mpas_pool_get_field(scratchPool, 'viscField',  viscField)
         call mpas_allocate_scratch_field(viscField, .true.)
         visc => viscField % array

         call mpas_pool_get_field(scratchPool, 'sigmaxyField', sigmaxyField )
         call mpas_allocate_scratch_field(sigmaxyField, .true.)
         sigmaxy => sigmaxyField % array

         call mpas_pool_get_field(scratchPool, 'sigmayxField', sigmayxField )
         call mpas_allocate_scratch_field(sigmayxField, .true.)
         sigmayx => sigmayxField % array

         call mpas_pool_get_field(scratchPool, 'sigmaxxField', sigmaxxField )
         call mpas_allocate_scratch_field(sigmaxxField, .true.)
         sigmaxx => sigmaxxField % array

         call mpas_pool_get_field(scratchPool, 'sigmayyField', sigmayyField )
         call mpas_allocate_scratch_field(sigmayyField, .true.)
         sigmayy => sigmayyField % array

         call mpas_pool_get_field(scratchPool, 'sigmaminField', sigmaminField )
         call mpas_allocate_scratch_field(sigmaminField, .true.)
         sigmamin => sigmaminField % array

         call mpas_pool_get_field(scratchPool, 'sigmamaxField', sigmamaxField )
         call mpas_allocate_scratch_field(sigmamaxField, .true.)
         sigmamax => sigmamaxField % array

         call mpas_pool_get_field(scratchPool, 'alphaField', alphaField )
         call mpas_allocate_scratch_field(alphaField, .true.)
         alpha => alphaField % array

         call mpas_pool_get_field(scratchPool, 'epsmaxField', epsmaxField )
         call mpas_allocate_scratch_field(epsmaxField, .true.)
         epsmax => epsmaxField % array

         call mpas_pool_get_field(scratchPool, 'nstarField', nstarField )
         call mpas_allocate_scratch_field(nstarField, .true.)
         nstar => nstarField % array

         !call mpas_pool_get_field(scratchPool, 'szeroField', szeroField )
         !call mpas_allocate_scratch_field(szeroField, .true.)
         !szero => szeroField % array

         !call mpas_pool_get_field(scratchPool, 'sourceField', sourceField )
         !call mpas_allocate_scratch_field(sourceField, .true.)
         !source => sourceField % array

         !call mpas_pool_get_field(scratchPool, 'ddamagedtField', ddamagedtField )
         !call mpas_allocate_scratch_field(ddamagedtField, .true.)
         !ddamagedt => ddamagedtField % array

         call mpas_pool_get_field(scratchPool, 'forceField', forceField )
         call mpas_allocate_scratch_field(forceField, .true.)
         force => forceField % array

         call mpas_pool_get_field(scratchPool, 'meanTField', meanTField )
         call mpas_allocate_scratch_field(meanTField, .true.)
         meanT => meanTField % array

         call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
         call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
         rhoi = config_ice_density
         rhoo = config_ocean_density

         meanT(:) = sum(temperature(:,:), dim=1)/nVertLevels

         call li_compute_gradient_2d(meshPool, uReconstructX(1,:), exx, exy, err_tmp)
         err = ior(err, err_tmp)

         call li_compute_gradient_2d(meshPool, uReconstructY(1,:), eyx, eyy, err_tmp)
         err = ior(err, err_tmp)

         call li_calculate_flowParamA(meshPool, temperature, thickness, flowParamA, err_tmp)
         err = ior(err, err_tmp)


         meanFlowParamA(:) = sum(flowParamA(:,:), dim=1)/nVertLevels
         call mpas_log_write("meanFlowParamA: Min=$r, Max=$r", &
                    realArgs=(/minval(meanFlowParamA), maxval(meanFlowParamA)/))
         if (config_use_constant_A) then
             visc(:) = 0.5*stiffnessFactor(:)*config_default_flowParamA**(-1.0/3.0)*(sqrt(exx(:)**2+eyy(:)**2 + exx(:)*eyy(:) + 0.25*(exy(:)+eyx(:))**2+1.0e-30))**((1.0-3.0)/3.0)
         else
             visc(:) = 0.5*stiffnessFactor(:)*meanFlowParamA(:)**(-1.0/3.0)*(sqrt(exx(:)**2+eyy(:)**2 + exx(:)*eyy(:) + 0.25*(exy(:)+eyx(:))**2))**((1.0-3.0)/3.0)
         end if
         ! recover effective viscosity from vel. and temp. fields

         sigmaxy(:) = 2*visc(:)*(exy(:)+eyx(:))/2.0
         sigmayx(:) = 2*visc(:)*(eyx(:)+exy(:))/2.0

         sigmaxx(:) = 2*visc(:)*exx(:)
         sigmayy(:) = 2*visc(:)*eyy(:)
         !deviatoric stress

         sigmamin(:) = (sigmaxx(:)+sigmayy(:))/2.0-sqrt(((sigmaxx(:)-sigmayy(:))/2.0)**2+((sigmaxy(:)+sigmayx(:))/2)**2 + 1.0e-30)
         sigmamax(:) = (sigmaxx(:)+sigmayy(:))/2.0+sqrt(((sigmaxx(:)-sigmayy(:))/2.0)**2+((sigmaxy(:)+sigmayx(:))/2)**2 + 1.0e-30)
         !deviatoric principal stress

         alpha(:) = sigmamin(:) / (sigmamax(:)+1.0e-30)
         epsmax(:) = sigmamax(:) / (2.0d0*visc(:)+1.0e-30)

         nstar(:) = 12.0d0*(1.0d0+alpha(:)+alpha(:)**2)/(4.0d0*(1.0d0+alpha(:)+alpha(:)**2)+6.0d0*alpha(:)**2+1.0e-30) 
          ! Compute the effective flow law exponent

         szero(:) = rhoi*(rhoo-rhoi)*gravity*thickness(:)/(2.0d0*sigmamax(:)*rhoo+1.0e-30)
          ! Compute the hydrostatic to tensile stress ratio

         !do iCell = 1, nCells
         !    if (li_mask_is_floating_ice(cellMask(iCell))) then
         !        floatingBasalMassBal(iCell) = 0.2/31536000.0*tanh((lowerSurface(iCell)-bedTopography(iCell))/75.0) * max(-100-lowerSurface(iCell),0.0)
         !    end if
         !end do
         ! Temporary change of BMB for MISMIP+ Ice1r

         source(:) = nstar(:) * (1.0d0-szero(:))*epsmax(:)-floatingBasalMassBal(:)/rhoi/(thickness(:)+1.0e-30)

         seconds = daysSinceStart*24.0*3600.0

         ! below is the force for the manufacured experiments
         !force(:) = -0.44/2.0 * (uReconstructX(1,:)/1000.0*(1+epsmax(:)*seconds)*exp(-uReconstructX(1,:)*seconds/1000.0) +    &
         !    source(:) * (1+exp(-uReconstructX(1,:)*seconds/1000.0)))
         !force(:) = -0.44/2.0 * (uReconstructX(1,:)/1000.0*exp(-uReconstructX(1,:)*seconds/1000.0) +    &
         !    source(:) * (1+exp(-uReconstructX(1,:)*seconds/1000.0)))

         ! set force to zero if not manufactured
         force(:) = 0.0

         ddamagedt(:) = source(:) * damage(:) + force(:)

         damage(:) = damage(:) + ddamagedt(:) * deltat

         !damage(:) =  0.5
         ! temporarily set damage to a constant value, 0.8

         damage_nye(:) = (2.0+alpha(:)) * sigmamax(:) / (gravity*(rhoo-rhoi)*thickness(:)+1e-30)
         !damage_nye(:) = 0.5
         ! temporarily set damage_nye to a constant value, 0.8

         do iCell = 1, nCells
            if (damage(iCell) > config_damage_restore_threshold) then
                damage_max(iCell) = damage(iCell)
            end if
         end do

         do iCell = 1, nCells
            if (damage(iCell) < 0.0_RKIND) then
                damage(iCell) = 0.0_RKIND
            end if
            if (damage(iCell) < damage_nye(iCell)) then
                damage(iCell) = damage_nye(iCell)
            end if
            if (damage(iCell) > 1.0_RKIND) then
                damage(iCell) = 1.0_RKIND
            end if
         end do


         !do iCell = 1, nCells
         !   if ((li_mask_is_grounded_ice(cellMask(iCell))) .or. (.not.(li_mask_is_dynamic_ice(cellMask(iCell))))) then
         !       damage(iCell) = 0.0_RKIND
         !   end if
         !end do

         do iCell = 1, nCells
            if (li_mask_is_grounded_ice(cellMask(iCell))) then
                damage(iCell) = 0.0_RKIND
            end if
         end do

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
        ! set the damage for a cell at the GL to the average damage value for its neighbor cells downstream on the ice shelf 
                        

        if (.true.) then

             call mpas_log_write("source value range: Min=$r, Max=$r", &
                        realArgs=(/minval(source), maxval(source)/))
             call mpas_log_write("szero  value range: Min=$r, Max=$r", &
                        realArgs=(/minval(szero), maxval(szero)/))
             call mpas_log_write("nstar value range: Min=$r, Max=$r", &
                        realArgs=(/minval(nstar), maxval(nstar)/))
             call mpas_log_write("deltat value range: Min=$r, Max=$r", &
                        realArgs=(/deltat, deltat/))
             call mpas_log_write("sigmamax value range: Min=$r, Max=$r", &
                        realArgs=(/minval(sigmamax), maxval(sigmamax)/))
             call mpas_log_write("damage value range: Min=$r, Max=$r", &
                        realArgs=(/minval(damage), maxval(damage)/))
             call mpas_log_write("A value range: Min=$r, Max=$r", &
                        realArgs=(/minval(meanFlowParamA), maxval(meanFlowParamA)/))
             call mpas_log_write("alpha value range: Min=$r, Max=$r", &
                        realArgs=(/minval(alpha), maxval(alpha)/))
        end if
        ! temporary message output in the log file


         call mpas_deallocate_scratch_field(calvingFrontMaskField, .true.)
         call mpas_deallocate_scratch_field(meanFlowParamAField, .true.)
         call mpas_deallocate_scratch_field(viscField, .true.)
         call mpas_deallocate_scratch_field(sigmaxxField, .true.)
         call mpas_deallocate_scratch_field(sigmayyField, .true.)
         call mpas_deallocate_scratch_field(sigmaminField, .true.)
         call mpas_deallocate_scratch_field(sigmamaxField, .true.)
         call mpas_deallocate_scratch_field(alphaField, .true.)
         call mpas_deallocate_scratch_field(epsmaxField, .true.)
         call mpas_deallocate_scratch_field(nstarField, .true.)
         !call mpas_deallocate_scratch_field(szeroField, .true.)
         !call mpas_deallocate_scratch_field(sourceField, .true.)
         !call mpas_deallocate_scratch_field(ddamagedtField, .true.)
         call mpas_deallocate_scratch_field(forceField, .true.)
         call mpas_deallocate_scratch_field(meanTField, .true.)
         ! call mpas_deallocate_scratch_field(damageField, .true.)

         block => block % next
             
      enddo

   end subroutine calculate_damage 

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine restore_damage
!
!> \author Tong Zhang
!> \date   Oct. 2019
!> \details restore the damage to the original value before advection if it exceeds some large number
!-----------------------------------------------------------------------

   subroutine restore_damage(domain, err)

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
      logical, pointer :: config_print_calving_info
      real (kind=RKIND), dimension(:), pointer :: damage
      real (kind=RKIND), dimension(:), pointer :: damage_max
      real (kind=RKIND), dimension(:), pointer :: stiffnessFactor
      real (kind=RKIND), pointer ::  &
           config_ice_density,            & ! ice density
           config_ocean_density             ! ocean density

      real (kind=RKIND), save :: rhoi    ! ice density (kg m^{-3}), copied from config_ice_density
      real (kind=RKIND), save :: rhoo    ! ocean density (kg m^{-3}), copied from config_ocean_density

      integer, dimension(:), pointer :: cellMask
      type (field1dInteger), pointer :: calvingFrontMaskField
      integer, dimension(:), pointer :: calvingFrontMask
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, dvEdge
      integer, pointer :: nCells, nVertLevels
      integer :: iCell, jCell, iNeighbor, n_damage_downstream
      real(kind=RKIND) :: damage_downstream, damage_tmp, cellCalvingFrontLength, cellCalvingFrontHeight, seconds
      integer :: err_tmp, loc

      err = 0


      call mpas_pool_get_config(liConfigs, 'config_damage_stiffness_min', config_damage_stiffness_min)
      call mpas_pool_get_config(liConfigs, 'config_damage_rheology_coupling', config_damage_rheology_coupling)

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
         call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'damage', damage)
         call mpas_pool_get_array(geometryPool, 'damage_max', damage_max)
         call mpas_pool_get_array(velocityPool, 'stiffnessFactor', stiffnessFactor)


         do iCell = 1, nCells
            if (damage(iCell) < damage_max(iCell)) then
                damage(iCell) = damage_max(iCell)
            end if
         end do
         ! put the damage_max value back to restore the damage value

         where (damage < 0.0_RKIND)
             damage = 0.0_RKIND
         end where

         where (damage > 1.0_RKIND)
             damage = 1.0_RKIND
         end where


         !do iCell = 1, nCells
         !   if ((li_mask_is_grounded_ice(cellMask(iCell))) .or. (.not.(li_mask_is_dynamic_ice(cellMask(iCell))))) then
         !       damage(iCell) = 0.0_RKIND
         !   end if
         !end do

         do iCell = 1, nCells
            if (li_mask_is_grounded_ice(cellMask(iCell))) then
                damage(iCell) = 0.0_RKIND
            end if
         end do

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

        if (config_damage_rheology_coupling) then
            do iCell = 1, nCells
                if (li_mask_is_floating_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell))) then
                    stiffnessFactor(iCell) = 1.0_RKIND - damage(iCell)
                    if (stiffnessFactor(iCell) < config_damage_stiffness_min) then
                        stiffnessFactor(iCell) = config_damage_stiffness_min
                    end if
                end if
            end do
        end if


         block => block % next
             
      enddo

   end subroutine restore_damage 

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine apply_calving_velocity
!
!> \brief Convert a calving velocity to an ice thickness removal
!> \author Matthew Hoffman
!> \date   March 2020
!> \details This routine takes a calving velocity and converts it to a
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
!> The output of this routine is calvingThickness, which then needs to be applied to thickness
!> by the calling routine.
!-----------------------------------------------------------------------

   subroutine apply_calving_velocity(meshPool, geometryPool, velocityPool, domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), pointer, intent(in) :: velocityPool !< Input: velocity pool

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain          !< Input: domain object
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input: geometry pool

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nEdges, nCells, nCellsSolve
      integer :: iEdge, iCell, jCell, iNeighbor
      integer :: nEmptyNeighbors
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      integer, dimension(:), pointer :: cellMask, edgeMask
      integer, dimension(:), pointer :: groundedMarineMarginMask
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:), pointer :: angleEdge
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:), pointer :: requiredCalvingVolumeNonDynEdge, requiredCalvingVolumeDynEdge
      real (kind=RKIND), dimension(:), pointer :: requiredCalvingVolumeNonDynCell, requiredCalvingVolumeDynCell
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructX
      real (kind=RKIND), dimension(:,:), pointer :: uReconstructY
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      real (kind=RKIND), dimension(:), allocatable :: thicknessForCalving
      real (kind=RKIND), dimension(:), pointer :: calvedVolumeNonDynCell, calvedVolumeDynCell
      real (kind=RKIND), dimension(:), pointer :: uncalvedVolumeNonDynCell, uncalvedVolumeDynCell
      real (kind=RKIND), dimension(:), allocatable :: cellVolume
      real(kind=RKIND) :: calvingSubtotal1, calvingSubtotal2, calvingSubtotal3
      real(kind=RKIND) :: thkSum
      integer :: thkCount
      real(kind=RKIND) :: removeVolumeHere
      real(kind=RKIND) :: volumeAvailableToPass !< a temp. var. for accounting purposes that indicates how much volume
                          !< to be calved is 'leftover' after trying to calve non-dynamic cells and can be transferred to
                          !< neighboring dynamic cells
      real(kind=RKIND) :: calvLengthCell, calvLengthEdge
      real(kind=RKIND), dimension(6) :: localInfo, globalInfo
      real(kind=RKIND) :: edgeLengthScaling
      real(kind=RKIND), parameter :: calvingSmallThk = 1.0e-8 ! in meters, a small thickness threshold
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'uncalvedVolumeNonDynCell', uncalvedVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'uncalvedVolumeDynCell', uncalvedVolumeDynCell)
      call mpas_pool_get_array(geometryPool, 'calvedVolumeNonDynCell', calvedVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'calvedVolumeDynCell', calvedVolumeDynCell)
      call mpas_pool_get_array(geometryPool, 'groundedMarineMarginMask', groundedMarineMarginMask)
      call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeNonDynEdge', requiredCalvingVolumeNonDynEdge)
      call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeDynEdge', requiredCalvingVolumeDynEdge)
      call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeNonDynCell', requiredCalvingVolumeNonDynCell)
      call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeDynCell', requiredCalvingVolumeDynCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)

      groundedMarineMarginMask(:) = 0
      do iCell = 1, nCells
         ! Define grounded marine marginal cells as those that are (1) at the ice
         ! margin, (2) have at least one neighboring cell without ice, (3) contain
         ! grounded ice, and (4) have bed topography below sea level.
         ! OR is adjacent to an inactive floating margin cell

         ! Check if neighboring cells contain ice and have bed topo below sea level
         nEmptyNeighbors = 0
         do iEdge = 1, nEdgesOnCell(iCell)
           iNeighbor = cellsOnCell(iEdge, iCell)
           if ( (((thickness(iNeighbor) == 0.0_RKIND) &
              .and. bedTopography(iNeighbor) < config_sea_level)) &
              .or. ((li_mask_is_floating_ice(cellMask(iNeighbor)) &
              .and. li_mask_is_margin(cellMask(iNeighbor)) &
              .and. (.not. li_mask_is_dynamic_ice(cellMask(iNeighbor)))))) then
               nEmptyNeighbors = nEmptyNeighbors + 1
           endif
         enddo

         if ( nEmptyNeighbors > 0 &
           .and. li_mask_is_grounded_ice(cellMask(iCell)) &
           .and. bedTopography(iCell) < config_sea_level &
           .and. li_mask_is_dynamic_ice(cellMask(iCell)) ) then

            groundedMarineMarginMask(iCell) = 1
         else
            groundedMarineMarginMask(iCell) = 0
         endif
      enddo


      ! Init fields for accounting
      calvingThickness(:) = 0.0_RKIND
      calvedVolumeNonDynCell(:) = 0.0_RKIND
      calvedVolumeDynCell(:) = 0.0_RKIND
      uncalvedVolumeNonDynCell(:) = 0.0_RKIND
      uncalvedVolumeDynCell(:) = 0.0_RKIND
      allocate(cellVolume(nCells+1))
      cellVolume(:) = areaCell(:) * thickness(:)
      requiredCalvingVolumeNonDynEdge(:) = 0.0_RKIND
      requiredCalvingVolumeNonDynCell(:) = 0.0_RKIND
      requiredCalvingVolumeDynEdge(:) = 0.0_RKIND
      requiredCalvingVolumeDynCell(:) = 0.0_RKIND
      allocate(thicknessForCalving(nCells+1))
      thicknessForCalving = 0.0_RKIND



      ! 1. Calculate calving rate for all non-dynamic cells by working with their ocean-going edges
      calvingSubtotal1 = 0.0_RKIND
      do iCell = 1, nCellsSolve
         if (li_mask_is_floating_ice(cellMask(iCell)) .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell)))) then
            ! a. Translate the calving front height based on dynamic cells to the non-dynamic locations
            ! find mean (or min?) of thickness in dynamic neighbors
            thkSum = 0.0_RKIND
            thkCount = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(cellMask(jCell))) then
                  thkSum = thkSum + thickness(jCell)
                  thkCount = thkCount + 1
               endif
            enddo
            if (thkCount == 0) then
               !call mpas_log_write("Found a stranded non-dynamic floating cell: cell $i with thickness=$r m.", MPAS_LOG_WARN, &
               !   intArgs=(/iCell/), realArgs=(/thickness(iCell)/))
            else
               thicknessForCalving(iCell) = thkSum / real(thkCount, kind=RKIND)
            endif
            ! b. Translate the calvingVelocity and thickness from non-dynamic cells to their ocean-going edges
            !    to calculate the calving volume on those edges and this cell
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_floating_ice(edgeMask(iEdge)) .and. li_mask_is_margin(edgeMask(iEdge)) .and. &
                   .not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level & ! ensure margin is w/ocn
                   ) then
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uReconstructX(1,iCell), uReconstructY(1,iCell))
                  requiredCalvingVolumeNonDynEdge(iEdge) = calvingVelocity(iCell) * &
                          edgeLengthScaling * dvEdge(iEdge) * thicknessForCalving(iCell) * deltat
                  requiredCalvingVolumeNonDynCell(iCell) = requiredCalvingVolumeNonDynCell(iCell) + &
                          requiredCalvingVolumeNonDynEdge(iEdge) ! Keep running total
               endif
            enddo
            ! c. Calculate calvingThickness here
            removeVolumeHere = min(cellVolume(iCell), requiredCalvingVolumeNonDynCell(iCell))  ! Don't use more than available
            calvingThickness(iCell) = removeVolumeHere / areaCell(iCell)
            calvedVolumeNonDynCell(iCell) = removeVolumeHere
            uncalvedVolumeNonDynCell(iCell) = requiredCalvingVolumeNonDynCell(iCell) - removeVolumeHere
            cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere
            calvingSubtotal1 = calvingSubtotal1 + removeVolumeHere
         endif
      enddo
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'uncalvedVolumeNonDynCell')
      call mpas_timer_stop("halo updates")
      !call mpas_log_write("Done calculating calving for non-dynamic floating cells. Removed $r m^3", realArgs=(/calvingSubtotal1/))


      ! 2. Calculate calving for dynamic cells
      calvingSubtotal2 = 0.0_RKIND

      ! a. Calculate calving on dynamic margin edges
      do iCell = 1, nCellsSolve
         if ( ( (li_mask_is_floating_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell)) ) & ! floating dynamic
                .or.  (groundedMarineMarginMask(iCell) == 1) ) &
              .and. li_mask_is_margin(cellMask(iCell)) ) then
            ! only consider dynamic margins adjacent to open ocean (need to exclude dynamic margins adjacent to land)
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if ((.not. li_mask_is_ice(cellMask(jCell))) .and. (bedTopography(jCell) <= config_sea_level)) then
                  iEdge = edgesOnCell(iNeighbor, iCell)
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uReconstructX(1,iCell), uReconstructY(1,iCell))
                  requiredCalvingVolumeDynEdge(iEdge) = calvingVelocity(iCell) * &
                          edgeLengthScaling * dvEdge(iEdge) * thickness(iCell) * deltat
               endif
            enddo
         endif
      enddo

      ! b. copy calving remaining in non-dynamic cells to dynamic edges
      ! Assume height and velocity are uniform, but edge length is not.
      do iCell = 1, nCells
         if ( (.not. li_mask_is_dynamic_ice(cellMask(iCell))) .and. &
                 uncalvedVolumeNonDynCell(iCell) > 0.0_RKIND) then
            ! This is a non-dynamic location that still has calving to offer - a location where calving needs to be propagated

            !call mpas_log_write("Passing calving from non-dynamic cell $i.  $r available to pass", &
            !        realArgs=(/uncalvedVolumeNonDynCell(iCell)/), intArgs=(/iCell/))

            volumeAvailableToPass = uncalvedVolumeNonDynCell(iCell)

            ! Find total length of interface with dynamic cells at this cell (i.e. calculate calvLengthCell)
            calvLengthCell = 0.0_RKIND
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(edgeMask(iEdge))) then
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uReconstructX(1,iCell), uReconstructY(1,iCell))
                  calvLengthEdge = edgeLengthScaling * dvEdge(iEdge)
                  calvLengthCell = calvLengthCell + calvLengthEdge
               endif
            enddo

            ! Now that we know calvLengthCell, pass along the required calving volume relative to the interface length
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(edgeMask(iEdge))) then
                  edgeLengthScaling = scale_edge_length(angleEdge(iEdge), uReconstructX(1,iCell), uReconstructY(1,iCell))
                  calvLengthEdge = edgeLengthScaling * dvEdge(iEdge)
                  if (requiredCalvingVolumeDynEdge(iEdge) > 0.0_RKIND) then
                     call mpas_log_write("Unexpectedly found a dynamic edge that already has calving assigned to it." // &
                           "  There is a flaw in apply_calving_velocity that needs to be fixed!", MPAS_LOG_ERR)
                     err_tmp = 1
                     err = ior(err, err_tmp)
                  endif
                  requiredCalvingVolumeDynEdge(iEdge) = calvLengthEdge / calvLengthCell * volumeAvailableToPass
                  uncalvedVolumeNonDynCell(iCell) = uncalvedVolumeNonDynCell(iCell) - requiredCalvingVolumeDynEdge(iEdge)
                  !call mpas_log_write("   Passed calving $r from non-dynamic cell $i to dynamic cell $i", &
                  !        realArgs=(/requiredCalvingVolumeDynEdge(iEdge)/), intArgs=(/iCell, jCell/))
               endif
            enddo
         endif
      enddo
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'requiredCalvingVolumeDynEdge')
      call mpas_timer_stop("halo updates")

      ! c. Now calculate calvingThickness for each cell
      do iCell = 1, nCellsSolve
         if ((li_mask_is_dynamic_ice(cellMask(iCell))) .and. (bedTopography(iCell) <= config_sea_level)) then
            ! Can loop over all dyn cells - only ones with calving on their edges will have calving applied.
            do iNeighbor = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(iNeighbor, iCell)
               ! No need to check what type of edge this is - only required edges are nonzero
               requiredCalvingVolumeDynCell(iCell) = requiredCalvingVolumeDynCell(iCell) + &
                          requiredCalvingVolumeDynEdge(iEdge) ! Keep running total
            enddo
            ! c. Apply calving here
            removeVolumeHere = min(cellVolume(iCell), requiredCalvingVolumeDynCell(iCell))  ! Don't use more than available
            calvingThickness(iCell) = removeVolumeHere / areaCell(iCell)
            calvedVolumeDynCell(iCell) = removeVolumeHere
            uncalvedVolumeDynCell(iCell) = requiredCalvingVolumeDynCell(iCell) - removeVolumeHere
            cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere
            calvingSubtotal2 = calvingSubtotal2 + removeVolumeHere
         endif
      enddo
      !call mpas_log_write("Done calculating calving for dynamic floating cells. Removed $r m^3", realArgs=(/calvingSubtotal2/))


      ! Clean up to account for roundoff level errors that can occur
      do iCell = 1, nCellsSolve
         if (abs(calvingThickness(iCell) - thickness(iCell)) < calvingSmallThk) then
            calvingThickness(iCell) = thickness(iCell)
         endif
      enddo

      ! Clean up - zap any stranded ice
      calvingSubtotal3 = 0.0_RKIND
      do iCell = 1, nCellsSolve
         if (li_mask_is_floating_ice(cellMask(iCell)) .and. (.not. li_mask_is_dynamic_ice(cellMask(iCell)))) then
            thkCount = 0
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_dynamic_ice(cellMask(jCell))) then
                  thkCount = thkCount + 1
               endif
            enddo
            if (thkCount == 0) then
               calvingThickness(iCell) = thickness(iCell)
               calvingSubtotal3 = calvingSubtotal3 + calvingThickness(iCell)
            endif
         endif
         if (isnan(areaCell(iCell))) then
            call mpas_log_write("NaN detected in calvingThickness at cell $i", MPAS_LOG_ERR, intArgs=(/iCell/))
            err_tmp = 1
            err = ior(err, err_tmp)
         endif
      enddo


      ! End of routine accounting/reporting
      localInfo(1) = calvingSubtotal1
      localInfo(2) = calvingSubtotal2
      localInfo(3) = calvingSubtotal3
      localInfo(4) = sum(calvingThickness(1:nCellsSolve) * areaCell(1:nCellsSolve))
      localInfo(5) = sum(uncalvedVolumeNonDynCell(1:nCellsSolve))
      localInfo(6) = sum(uncalvedVolumeDynCell(1:nCellsSolve))
      ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
      call mpas_dmpar_sum_real_array(domain % dminfo, 6, localInfo, globalInfo)
      call mpas_log_write("== Calving complete. Total calved                    = $r", realArgs = (/globalInfo(4)/))
      call mpas_log_write("==                   Calved from non-dynamic cells   = $r", realArgs = (/globalInfo(1)/))
      call mpas_log_write("==                   Calved from dynamic cells       = $r", realArgs = (/globalInfo(2)/))
      call mpas_log_write("==                   Stranded floating cells deleted = $r", realArgs = (/globalInfo(3)/))
      call mpas_log_write("== Uncalved volumes: Non-dynamic cells=$r; Dynamic cells=$r", realArgs=(/globalInfo(5),globalInfo(6)/))
      if (((globalInfo(5) + globalInfo(6)) / (globalInfo(4) + 1.0e-30_RKIND) > 0.001_RKIND) .and. &
          (globalInfo(4) > 1000.0_RKIND**2)) then ! Include some small amount of total calving for comparison
         call mpas_log_write("Calving scheme failed to calve $r m^3. ($r% of total calved)", MPAS_LOG_WARN, &
                 realArgs=(/globalInfo(5) + globalInfo(6), &
                 100.0_RKIND * (globalInfo(5) + globalInfo(6)) / (globalInfo(4)+1.0e-30_RKIND)/))
      endif
      if (((globalInfo(5) + globalInfo(6)) / (globalInfo(4) + 1.0e-30_RKIND) > 0.1_RKIND) .and. &
          (globalInfo(4) > 1000.0_RKIND**2)) then ! Include some small amount of total calving for comparison
         call mpas_log_write("Calving scheme failed to calve an amount greater than 10% of the ice calved.  " // &
                 "Try reducing time step or apply_calving_velocity may need improvements.", &
                 MPAS_LOG_ERR, realArgs=(/globalInfo(5) + globalInfo(6), &
                 100.0_RKIND * (globalInfo(5) + globalInfo(6)) / (globalInfo(4)+1.0e-30_RKIND)/))
         err_tmp = 1
         err = ior(err, err_tmp)
      endif
      if (globalInfo(3) / (globalInfo(4) + 1.0e-30_RKIND) > 0.01_RKIND) then
         ! If stranded ice deletion results in more than a small amount of the total calving flux, this routine will not
         ! be accurate.  If this error gets triggered, either try using a smaller timstep or this routine needs improving.
         call mpas_log_write("Deleting stranded floating cells accounts for more than 1% of calving loss." // &
                 "  Try using a smaller timestep or apply_calving_velocity may need improvements for this simulation.")
         err_tmp = 1
         err = ior(err, err_tmp)
      endif


      ! Update halos on calvingThickness before applying it.
      ! Testing seemed to indicate this is not necessary, but I don't understand why not, so leaving it.
      ! NOTE: THIS WILL NOT WORK ON MULTIPLE BLOCKS PER PROCESSOR
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'calvingThickness')
      call mpas_timer_stop("halo updates")

      deallocate(cellVolume)
      deallocate(thicknessForCalving)

    end subroutine apply_calving_velocity


    ! Helper function for subroutine apply_calving_velocity
    ! Calculates the amount to scale an edge length based on the orientation of the edge with the surface velocity
    function scale_edge_length(angleEdgeHere, u, v)
       real(kind=RKIND), intent(in) :: angleEdgeHere
       real(kind=RKIND), intent(in) :: u
       real(kind=RKIND), intent(in) :: v
       real(kind=RKIND) :: scale_edge_length

       real(kind=RKIND) :: mag

       mag = sqrt(u**2 + v**2)
       if (mag == 0.0_RKIND) mag = 1.0_RKIND
       !scale_edge_length =  abs(u/mag * cos(-1.0*angleEdgeHere) + v/mag * sin(-1.0*angleEdgeHere))  ! dot product of unit vectors
       scale_edge_length = abs(cos(angleEdgeHere - atan2(v, u)))
    end function scale_edge_length



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
      integer, dimension(:), pointer :: calvingMask
      integer, dimension(:), pointer :: cellMask
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
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'calvingMask', calvingMask)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)

         calvingThickness = 0.0_RKIND

         ! === apply calving ===
         where (li_mask_is_floating_ice(cellMask) .and. (calvingMask >= 1))
            calvingThickness = thickness
            thickness = 0.0_RKIND
         end where

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

      type (field1dInteger), pointer :: contiguousDynamicIceField
      type (field1dInteger), pointer :: contiguousDynamicIceOldField

      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: thickness
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell

      integer, dimension(:), pointer :: contiguousDynamicIce, contiguousDynamicIceOld
      integer, pointer :: nCells, nCellsSolve
      integer, dimension(:), pointer :: nCellsArray
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
      call mpas_pool_get_field(geometryPool, 'contiguousDynamicIce',  contiguousDynamicIceField)
      call mpas_allocate_scratch_field(contiguousDynamicIceField, single_block_in = .false.)

      call mpas_pool_get_field(geometryPool, 'contiguousDynamicIceOld',  contiguousDynamicIceOldField)
      call mpas_allocate_scratch_field(contiguousDynamicIceOldField, single_block_in = .false.)

      call mpas_log_write("Iceberg-detection flood-fill: allocated.")
      ! First mark grounded ice to initialize flood fill mask
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'contiguousDynamicIce', contiguousDynamicIce)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(geometryPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)

         ! make sure masks are up to date.  May not be necessary, but safer to do anyway.
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         call mpas_log_write("Iceberg-detection flood-fill: updated masks.")
         contiguousDynamicIce(:) = 0  ! initialize
         newMaskCountLocal = 0
         do iCell = 1, nCellsSolve
            if (li_mask_is_grounded_ice(cellMask(iCell)) .and. li_mask_is_dynamic_ice(cellMask(iCell))) then
               contiguousDynamicIce(iCell) = 1
               newMaskCountLocal = newMaskCountLocal + 1
            endif
         enddo
         call mpas_log_write("Initialized $i cells to local mask", intArgs=(/newMaskCountLocal/))

         block => block % next
      end do

      call mpas_log_write("Iceberg-detection flood-fill initialization complete.")

      ! Outer loop over processors (should also have a loop over blocks)
      ! Inner loop over cells on that processor

      ! Initialize global mask count
      call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocal, newMaskCountGlobal)
      globalLoopCount = 0

      call mpas_log_write("Initialized $i cells to global mask", intArgs=(/newMaskCountGlobal/))

      do while (newMaskCountGlobal > 0)
         globalLoopCount = globalLoopCount + 1
         call mpas_log_write("  Starting global processor loop $i", intArgs=(/globalLoopCount/))
         ! First Update halos
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'contiguousDynamicIce')
         call mpas_timer_stop("halo updates")

         ! initialize counter of cells locally updated during this outer loop
         newMaskCountLocalAccum = 0

         ! Now update (advance) mask locally

         block => domain % blocklist
         do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
            call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
            call mpas_pool_get_array(geometryPool, 'contiguousDynamicIce', contiguousDynamicIce)
            call mpas_pool_get_array(geometryPool, 'contiguousDynamicIceOld', contiguousDynamicIceOld)
            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
            call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)


            ! initialize local loop
            localLoopCount = 0
            newMaskCountLocal = 1  ! need to make sure we enter the loop
            do while (newMaskCountLocal > 0)
               localLoopCount = localLoopCount + 1
               call mpas_log_write("    Starting local cell loop $i", intArgs=(/localLoopCount/))

               ! initialize
               newMaskCountLocal = 0
               contiguousDynamicIceOld(:) = contiguousDynamicIce(:)

               do iCell = 1, nCellsArray(1) ! this gives owned cells only
                  if (contiguousDynamicIceOld(iCell) == 0 .and. &   ! this cell not yet marked
                     li_mask_is_dynamic_ice(cellMask(iCell))) then  ! and is dynamic
                     ! If it has a marked neighbor, then add it to the mask
                     do n = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(n, iCell)
                        if (contiguousDynamicIceOld(jCell) == 1) then
                           contiguousDynamicIce(iCell) = 1
                           newMaskCountLocal = newMaskCountLocal + 1
                           exit ! skip the rest of this do-loop - no need to check additional neighbors
                        endif
                     enddo
                  endif ! if not already marked
               enddo ! loop over cells

               ! Accumulate cells added locally until we do the next global reduce
               newMaskCountLocalAccum = newMaskCountLocalAccum + newMaskCountLocal
               call mpas_log_write("    Added $i new cells to local mask", intArgs=(/newMaskCountLocal/))
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


      ! Now remove any ice that was not flood-filled - these are icebergs
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'contiguousDynamicIce', contiguousDynamicIce)
         call mpas_pool_get_array(geometryPool, 'contiguousDynamicIceOld', contiguousDynamicIceOld)
         call mpas_pool_get_dimension(geometryPool, 'nCells', nCells)
         call mpas_pool_get_dimension(geometryPool, 'nCellsSolve', nCellsSolve)

         contiguousDynamicIceOld(:) = 0
         do iCell = 1, nCellsSolve
            if (contiguousDynamicIce(iCell) == 0 .and. li_mask_is_dynamic_ice(cellMask(iCell))) then
               calvingThickness(iCell) = calvingThickness(iCell) + thickness(iCell)  ! remove any remaining ice here
               thickness(iCell) = 0.0_RKIND
               contiguousDynamicIceOld(iCell) = 1 ! debug: make this a mask of where icebergs were removed
            endif
         enddo
         localIcebergCellCount = sum(contiguousDynamicIceOld)

!         where(contiguousDynamicIce == 0 .and. li_mask_is_dynamic_ice(cellMask))
!            calvingThickness = calvingThickness + thickness  ! remove any remaining ice here
!            thickness = 0.0_RKIND
!         end where

         block => block % next
      end do

      ! Calculate total iceberg cells removed (not needed but useful for log - note requires extra reduce!)
      call mpas_dmpar_sum_int(domain % dminfo, localIcebergCellCount, globalIcebergCellCount)

      ! update halo - probably not needed?
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'thickness')
      call mpas_timer_stop("halo updates")

      ! clean up
      call mpas_deallocate_scratch_field(contiguousDynamicIceField, single_block_in=.false.)
      call mpas_deallocate_scratch_field(contiguousDynamicIceOldField, single_block_in=.false.)

      call mpas_log_write("Iceberg-detection flood-fill complete. Removed $i iceberg cells.", intArgs=(/globalIcebergCellCount/))
      call mpas_timer_stop("iceberg detection")
   end subroutine remove_icebergs


end module li_calving


