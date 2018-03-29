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

   public :: li_calve_ice, li_restore_calving_front

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
      ! the calvingThickness field gets applied to the thickness state field
      ! after this if-construct
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

      else

         call mpas_log_write("Invalid option for config_calving specified: " // trim(config_calving), MPAS_LOG_ERR)
         err = 1

      endif


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
                  call mpas_log_write("$i $r $r", intArgs=(/indexToCellID(iCell)/), realArgs=(/bedTopography(iCell), calvingThickness(iCell)/))
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
           waterfrac                  ! interior water fraction

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
         call mpas_pool_get_array(thermalPool, 'waterfrac', waterfrac)
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
                       waterfrac(:,iCell),            &
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
      real (kind=RKIND), dimension(:), pointer :: calvingThickness    ! thickness of ice that calves (computed in this subroutine)
      real (kind=RKIND), dimension(:), pointer :: thickness
      integer, dimension(:), pointer :: cellMask

      err = 0

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
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

         block => block % next
      enddo

   end subroutine floating_calving


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
      type (mpas_pool_type), pointer :: geometryPool
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
      real (kind=RKIND), pointer :: config_sea_level
      real(kind=RKIND), pointer :: config_calving_thickness
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), dimension(:), pointer :: eigencalvingParameter
      real (kind=RKIND), dimension(:), pointer :: calvingThickness
      real (kind=RKIND), dimension(:), pointer :: requiredCalvingVolumeRate
      real (kind=RKIND), dimension(:), pointer :: calvingVelocity
      real (kind=RKIND), dimension(:), pointer :: eMax, eMin
      integer, dimension(:), pointer :: cellMask
      type (field1dInteger), pointer :: calvingFrontMaskField
      integer, dimension(:), pointer :: calvingFrontMask
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      integer, pointer :: nCells
      integer :: iCell, jCell, iNeighbor
      real(kind=RKIND) :: cellCalvingFrontLength, cellCalvingFrontHeight
      logical :: dynamicNeighbor
      integer :: err_tmp

      err = 0

      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_scalar_value', &
               config_calving_eigencalving_parameter_scalar_value)
      call mpas_pool_get_config(liConfigs, 'config_calving_eigencalving_parameter_source', config_calving_eigencalving_parameter_source)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)

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
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'eigencalvingParameter', eigencalvingParameter)
         call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeRate', requiredCalvingVolumeRate)
         call mpas_pool_get_array(geometryPool, 'calvingVelocity', calvingVelocity)
         call mpas_pool_get_array(velocityPool, 'eMax', eMax)
         call mpas_pool_get_array(velocityPool, 'eMin', eMin)

         call mpas_pool_get_field(scratchPool, 'iceCellMask2',  calvingFrontMaskField)
         call mpas_allocate_scratch_field(calvingFrontMaskField, .true.)
         calvingFrontMask => calvingFrontMaskField % array

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

         ! make mask for effective calving front.
         ! This is last dynamic cell, but also make sure it has a neighbor that is open ocean or thin floating ice.
         call calculate_calving_front_mask(meshPool, geometryPool, calvingFrontMask)

         calvingVelocity(:) = 0.0_RKIND
         requiredCalvingVolumeRate(:) = 0.0_RKIND
         ! First calculate the front retreat rate (Levermann eq. 1)
         calvingVelocity(:) = eigencalvingParameter(:) * max(0.0_RKIND, eMax(:)) * max(0.0_RKIND, eMin(:)) & ! m/s
               * real(li_mask_is_floating_ice_int(cellMask(:)), kind=RKIND) ! calculate only for floating ice - map of "potential" calving rate
         do iCell = 1, nCells
            if (calvingFrontMask(iCell) == 1) then

               ! convert to a volume flux per cell masking some assumptions
               ! get front length from edge lengths abutting ocean or thin ice
               ! get front height from max thickness of neighbors (since thickness near the edge could be screwy)
               cellCalvingFrontLength = 0.0_RKIND
               cellCalvingFrontHeight = 0.0_RKIND
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  jCell = cellsOnCell(iNeighbor, iCell)
                  if ( (li_mask_is_floating_ice(cellMask(jCell)) .and. .not. li_mask_is_dynamic_ice(cellMask(jCell))) .or. &  ! thin ice
                       (.not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level) ) then  ! open ocean
                       cellCalvingFrontLength = cellCalvingFrontLength + dvEdge(edgesOnCell(iNeighbor, iCell))
                  endif
                  cellCalvingFrontHeight = max(cellCalvingFrontHeight, thickness(jCell))
               enddo

               requiredCalvingVolumeRate(iCell) = calvingVelocity(iCell) * cellCalvingFrontLength * cellCalvingFrontHeight ! m^3/s
            endif
         enddo

         call distribute_calving_flux(meshPool, geometryPool, scratchPool, calvingFrontMask, err_tmp)
         err = ior(err, err_tmp)

         ! === apply calving ===
         thickness(:) = thickness(:) - calvingThickness(:)

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)
         call calculate_calving_front_mask(meshPool, geometryPool, calvingFrontMask)


         ! Now also remove thin floating, dynamic ice (based on chosen thickness threshold) after mask is updated.

         !where ((li_mask_is_floating_ice(cellMask) .and. li_mask_is_dynamic_ice(cellMask) .and. thickness<config_calving_thickness))
         ! The above commented version removed thin ice anywhere on the shelf.  In theory this seemed preferable,
         ! but in practice it has the tendency to create 'holes' in relatively stagnant areas of ice shelves along
         ! the grounding line.  Once these holes developed they would grow and eventually collapse the shelf from behind.

         ! This criteria below only remove too-thin ice at the new calving front,
         ! meaning just one 'row' of cells per timestep.  This could be expanded to continue
         ! removing ice backward until all connected too-thin ice has been removed.
         ! Tests of the current implementation show reasonable behavior.
         where (calvingFrontMask == 1 .and. thickness<config_calving_thickness)
             calvingThickness = calvingThickness + thickness
             thickness = 0.0_RKIND
         end where

         ! update mask
         call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
         err = ior(err, err_tmp)

         ! remove abandoned floating ice (i.e. icebergs) and add it to the calving flux
         ! Defined as: floating ice (dynamic or non-dynamic) that is not adjacent to dynamic ice (floating or grounded)
         ! This won't necessarily find all abandoned ice, but in practice it does a pretty good job at general cleanup
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
               endif
            endif
         enddo

         call mpas_deallocate_scratch_field(calvingFrontMaskField, .true.)

         block => block % next
      enddo

   end subroutine eigencalving


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!    routine distribute_calving_flux
!
!> \brief Apply a specified calving flux along the calving front
!> \author Matthew Hoffman
!> \date   Feb. 2018
!> \details Applies a specified calving flux along the calving front.
!> The calving volume needs to be distributed in three ways:
!> 1. We need to first remove any "thin" ice in front of this cell
!> 2. Then we remove ice from this cell
!> 3. If there is still additional ice to be removed, we have to recursively
!>    remove ice "inland" of this cell until the total required mass is removed
!-----------------------------------------------------------------------
   subroutine distribute_calving_flux(meshPool, geometryPool, scratchPool, calvingFrontMask, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh pool
      type (mpas_pool_type), pointer, intent(in) :: scratchPool !< Input: scratch pool
      integer, dimension(:), intent(in) :: calvingFrontMask

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer, intent(inout) :: geometryPool !< Input: geometry pool

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
!      logical, pointer :: config_print_calving_info
      real (kind=RKIND), dimension(:), pointer :: requiredCalvingVolumeRate  !< Input: the specified calving flux at calving front margin cells
      real (kind=RKIND), dimension(:), pointer :: calvingThickness !< Output: the applied calving rate as a thickness
      real (kind=RKIND), dimension(:), pointer :: thickness
      real (kind=RKIND), dimension(:), pointer :: areaCell
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), pointer :: deltat  !< time step (s)
      integer, pointer :: nCells
      integer :: iCell, iNeighbor, jCell
      real(kind=RKIND) :: volumeLeft
      real(kind=RKIND) :: removeVolumeHere
      type (field1dReal), pointer :: cellVolumeField
      real(kind=RKIND), dimension(:), pointer :: cellVolume
      real(kind=RKIND), dimension(:), pointer :: uncalvedVolume
      integer :: inwardNeighbors
      integer :: uncalvedCount
      real(kind=RKIND) :: uncalvedTotal

      err = 0

      ! get fields
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
      call mpas_pool_get_array(geometryPool, 'requiredCalvingVolumeRate', requiredCalvingVolumeRate)
      call mpas_pool_get_array(geometryPool, 'uncalvedVolume', uncalvedVolume)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      call mpas_pool_get_field(scratchPool, 'workCell',  cellVolumeField)
      call mpas_allocate_scratch_field(cellVolumeField, .true.)
      cellVolume => cellVolumeField % array

      calvingThickness(:) = 0.0_RKIND
      uncalvedVolume(:) = 0.0_RKIND

      cellVolume(:) = areaCell(:) * thickness(:)

      ! The calving volume needs to be distributed in three ways:
      ! 1. We need to first remove any "thin" ice in front of this cell
      ! 2. Then we remove ice from this cell
      ! 3. If there is still additional ice to be removed, we have to recursively
      !    remove ice "inland" of this cell until the total required mass is removed
      do iCell = 1, nCells
         if (calvingFrontMask(iCell) == 1) then
            volumeLeft = requiredCalvingVolumeRate(iCell) * deltat  ! units m^3

            ! First remove ice from "thin" neighbors
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_floating_ice(cellMask(jCell)) .and. .not. li_mask_is_dynamic_ice(cellMask(jCell))) then
                  ! this is a thin neighbor - remove as much ice from here as we can  TODO: could distribute this evenly amongst neighbors
                  removeVolumeHere = min(volumeLeft, cellVolume(jCell)) ! how much we want to remove here
                  calvingThickness(jCell) = calvingThickness(jCell) + removeVolumeHere / areaCell(jCell) ! apply to the field that will be used, in thickness units
                  cellVolume(jCell) = cellVolume(jCell) - removeVolumeHere ! update accounting on cell volume
                  volumeLeft = volumeLeft - removeVolumeHere ! update accounting on how much left to distribute from current iCell
               endif
            enddo

            if (volumeLeft > 0.0_RKIND) then
               ! Now remove ice from iCell
               removeVolumeHere = min(volumeLeft, cellVolume(iCell))
               calvingThickness(iCell) = calvingThickness(iCell) + removeVolumeHere / areaCell(iCell) ! apply to the field that will be used in thickness units
               cellVolume(iCell) = cellVolume(iCell) - removeVolumeHere ! update accounting on cell volume
               volumeLeft = volumeLeft - removeVolumeHere ! update accounting on how much left to distribute from current iCell
            endif

            if (volumeLeft > 0.0_RKIND) then
               ! Now remove ice from neighbors inward on shelf
               ! first count up how many there are
               inwardNeighbors = 0
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  jCell = cellsOnCell(iNeighbor, iCell)
                  if (li_mask_is_floating_ice(cellMask(jCell)) .and. li_mask_is_dynamic_ice(cellMask(jCell)) &
                          .and. (.not. li_mask_is_dynamic_margin(jCell))) then
                     inwardNeighbors = inwardNeighbors + 1
                  endif
               enddo
               ! Now distribute the flux evenly amongst the neighbors
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  jCell = cellsOnCell(iNeighbor, iCell)
                  if (li_mask_is_floating_ice(cellMask(jCell)) .and. li_mask_is_dynamic_ice(cellMask(jCell)) &
                          .and. (.not. li_mask_is_dynamic_margin(jCell))) then
                     ! this is thick neighbor that is not itself a margin - remove as much ice from here as we can  
                     removeVolumeHere = min(volumeLeft / real(inwardNeighbors, kind=RKIND), cellVolume(jCell)) ! how much we want to remove here
                     calvingThickness(jCell) = calvingThickness(jCell) + removeVolumeHere / areaCell(jCell) ! apply to the field that will be used in thickness units
                     cellVolume(jCell) = cellVolume(jCell) - removeVolumeHere ! update accounting on cell volume
                     volumeLeft = volumeLeft - removeVolumeHere ! update accounting on how much left to distribute from current iCell
                  endif
               enddo
               !TODO: need to recursively distribute across neighbors until fully depleted :(
            endif

            ! If we didn't calve enough ice, record that to allow assessment of how bad that is.
            if (volumeLeft > 0.0_RKIND) then
               uncalvedVolume(iCell) = 0.0_RKIND
            endif

         endif ! if cell is calving margin
      enddo ! cell loop

      if (maxval(uncalvedVolume) > 0.0_RKIND) then

         uncalvedTotal = sum(uncalvedVolume)
         uncalvedCount = count(uncalvedVolume > 0.0_RKIND)

         call mpas_log_write("distribute_calving_flux failed to distribute all required ice - ice was left after depleting all neighbors." &
                 // "  Search needs to be expanded to neighbors' neighbors.")
         call mpas_log_write("   On this processor: $i cells contain uncalved ice, for a total uncalved volume of $r m^3 ($r%).", &
                 MPAS_LOG_WARN, intArgs=(/uncalvedCount/), &
                 realArgs=(/uncalvedTotal, 100.0_RKIND * uncalvedTotal/(requiredCalvingVolumeRate(iCell) * deltat)/))
      endif

      call mpas_deallocate_scratch_field(cellVolumeField, .true.)

   end subroutine distribute_calving_flux


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
      integer, dimension(:) :: calvingFrontMask !< Output: calving front mask

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer, pointer :: nCells
      integer :: iCell, iNeighbor, jCell, jNeighbor, kCell
      logical :: oceanNeighbor
      integer, dimension(:), pointer :: nEdgesOnCell ! number of cells that border each cell
      integer, dimension(:,:), pointer :: cellsOnCell ! list of cells that neighbor each cell
      integer, dimension(:), pointer :: cellMask
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), pointer :: config_sea_level

      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

      ! get fields
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

      calvingFrontMask = 0 !initialize

      do iCell = 1, nCells
         if ( (li_mask_is_floating_ice(cellMask(iCell))) .and. (li_mask_is_dynamic_margin(cellMask(iCell))) ) then
            oceanNeighbor = .false.
            do iNeighbor = 1, nEdgesOnCell(iCell)
               jCell = cellsOnCell(iNeighbor, iCell)
               if (li_mask_is_floating_ice(cellMask(jCell)) .and. .not. li_mask_is_dynamic_ice(cellMask(jCell))) then
                  ! make sure this neighbor is adjacent to open ocean (and not thin floating ice up against the coast)
                  do jNeighbor = 1, nEdgesOnCell(jCell)
                     kCell = cellsOnCell(jNeighbor, jCell)
                     if (.not. li_mask_is_ice(cellMask(kCell)) .and. bedTopography(kCell) < config_sea_level) then
                        oceanNeighbor = .true. ! iCell neighbors thin ice that in turn neighbors open ocean
                     endif
                  enddo
               endif
               if (.not. li_mask_is_ice(cellMask(jCell)) .and. bedTopography(jCell) < config_sea_level) then
                  oceanNeighbor = .true. ! this is an open ocean neighbor
               endif
            enddo
            if (oceanNeighbor) then
               calvingFrontMask(iCell) = 1
            endif
         endif
      enddo

   end subroutine calculate_calving_front_mask


end module li_calving


