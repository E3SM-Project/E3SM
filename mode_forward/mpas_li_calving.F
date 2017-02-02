! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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
   use mpas_constants
   use mpas_dmpar
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

      type (dm_info), pointer :: dminfo
      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: velocityPool

      ! calving-relevant config options
      character (len=StrKIND), pointer :: config_calving
      logical, pointer :: config_print_calving_info, config_data_calving

      real(kind=RKIND), pointer :: config_calving_topography, &
                                   config_calving_thickness, &
                                   config_calving_timescale, &
                                   config_sea_level, &
                                   config_dynamic_thickness,  &
                                   config_ice_density,  &
                                   config_ocean_density

      ! fields for calving masks
      ! The calvingLaw mask is used for the floating and topographic_threshold options
      ! The activeForCalving, inactiveMargin and ocean masks are used for the more complex thickness_threshold option

      type (field1dInteger), pointer :: calvingLawMaskField
      integer, dimension(:), pointer :: calvingLawMask     ! = 1 where calving-law criterion is satisfied, else = 0

      type (field1dInteger), pointer :: activeForCalvingMaskField
      integer, dimension(:), pointer :: activeForCalvingMask   ! = 1 for grounded cells thicker than config_dynamic_thickness;
                                                               ! = 1 for floating cells thicker than config_calving_thickness;
                                                               ! = 0 elsewhere

      type (field1dInteger), pointer :: inactiveMarginMaskField
      integer, dimension(:), pointer :: inactiveMarginMask
      ! = 1 for inactive cells (thin or no ice) that have 1 or more active neighbors

      type (field1dInteger), pointer :: oceanMaskField
      integer, dimension(:), pointer :: oceanMask   ! = 1 for cells that are not land and do not have active ice
                                                    ! includes floating cells with inactive ice

      integer, pointer :: nCells

      real (kind=RKIND), pointer :: deltat  !< time step (s)

      integer, dimension(:), pointer :: &
           nEdgesOnCell,     & ! number of cells that border each cell
           cellMask            ! bit mask describing whether ice is floating, dynamically active, etc.

      integer, dimension(:), pointer :: &
           indexToCellID       ! list of global cell IDs

      integer, dimension(:,:), pointer :: &
           cellsOnCell         ! list of cells that neighbor each cell

      real (kind=RKIND) ::  &
           calvingFraction,  & ! fraction of ice that calves in each column; depends on calving_timescale
           flotationThickness  ! thickness at which marine-based ice starts to float

      real (kind=RKIND), dimension(:), pointer :: &
           thickness,        & ! ice thickness
           bedTopography       ! bed topography (negative below sea level)

      real (kind=RKIND), dimension(:), pointer :: &
           calvingThickness    ! thickness of ice that calves (computed in this subroutine)
                               ! typically the entire ice thickness, but will be a fraction of the thickness
                               ! if calving_timescale > dt

      integer :: iCell, iCellOnCell, iCellNeighbor

      integer :: err_tmp

      !WHL - debug
      logical, parameter :: circular_shelf = .false.
      integer, parameter :: ncellsPerRow = 40
      integer, parameter :: nRows = 46
      integer :: i, iRow

      err = 0

      ! get config options
      dminfo => domain % dminfo
      call mpas_pool_get_config(liConfigs, 'config_calving', config_calving)
      call mpas_pool_get_config(liConfigs, 'config_data_calving', config_data_calving)
      call mpas_pool_get_config(liConfigs, 'config_calving_thickness', config_calving_thickness)
      call mpas_pool_get_config(liConfigs, 'config_calving_topography', config_calving_topography)
      call mpas_pool_get_config(liConfigs, 'config_calving_timescale', config_calving_timescale)
      call mpas_pool_get_config(liConfigs, 'config_print_calving_info', config_print_calving_info)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)

      ! Get deltat from first block (same on all blocks)
      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)

      ! based on the calving timescale, set the fraction of ice that calves
      if (config_calving_timescale > 0.0_RKIND) then
         calvingFraction = min(deltat/config_calving_timescale, 1.0_RKIND)
         !WHL - debug
         if (config_print_calving_info) then
            write(stdoutUnit,*) 'Calving a fraction of the ice thickness at each timestep'
            write(stdoutUnit,*) 'deltat (yr), calvingFraction =', deltat/scyr, calvingFraction
         endif
      else
         calvingFraction = 1.0_RKIND   ! calve the entire thickness in eligible columns
      endif

      ! block loop
      block => domain % blocklist
      do while (associated(block))

         ! get pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)  ! required for cellMask computation
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

         ! get dimensions
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         ! get required fields from the mesh pool
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)  ! diagnostic only

         ! get required fields from the geometry pool
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         if (config_print_calving_info) then
            write(stdoutUnit,*) 'Do ice calving, option = ', trim(config_calving)
            write(stdoutUnit,*) 'Calving timscale (yr) = ', config_calving_timescale / scyr

            !WHL - debug - for circular shelf test case
            if (circular_shelf) then
               write(stdoutUnit,*) ' '
               write(stdoutUnit,*) 'Ice thickness before calving'
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                     write(stdoutUnit,'(a3)',advance='no') '    '
                  endif
!!                 do i = 1, nCellsPerRow
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!!                    write(stdoutUnit,'(i5)',advance='no') iCell
                     write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
                  enddo
                  write(stdoutUnit,*) ' '
               enddo
            endif   ! circular_shelf

         endif  ! config_print_calving_info

         ! initialize
         calvingThickness = 0.0_RKIND

         ! compute calvingThickness based on the calving_config option
         ! the calvingThickness field gets applied to the thickness state field
         ! after this if-construct

         if (trim(config_calving) == 'none') then

            ! do nothing

         elseif (trim(config_calving) == 'thickness_threshold') then

            ! calve ice thinner than a thickness threshold

            ! Note: The thickness-threshold option is different from the others.
            !       For the other options, we look at each cell and determine whether it meets the calving-law criteria
            !        (e.g., ice is floating, or the topography lies below a given level).
            !        If a cell meets the criteria and lies in the calving domain (e.g., at the margin), it is calved.
            !       For the thickness-threshold option, ice thinner than config_calving_thickness is calved,
            !        but only if it lies beyond a protected ring of thin ice at the floating margin.
            !       The reason for this more complicated approach is that we do not want to remove all floating ice
            !         thinner than the calving thickness, because then we would remove thin ice that has just
            !         been advected from active cells at the margin, and the calving front would be unable to advance.
            !       By protecting a ring of inactive ice (thickness < config_calving_thickness) at the margin,
            !        we allow ice in these cells to thicken and become active, thus advancing the calving front.
            !       The calving front retreats when active floating ice thins to become inactive, removing protection
            !        from previously protected cells.
            !
            !      Specifically, the rules are as follows:
            !      - Mark cells as active-for-calving if either (1) grounded, with thickness > config_dynamic_thickness
            !        or (2) floating, with thickness > config_calving_thickness.
            !      - Mark cells as ocean if not land and not active.
            !      - Mark cells as lying on the inactive margin if not active, but with an active neighbor.
            !      - Calve ice in ocean cells that are not on the inactive margin.

            if (config_print_calving_info) then
               write(stdoutUnit,*) 'Calving thickness (m) =', config_calving_thickness
               write(stdoutUnit,*) 'Dynamic thickness (m) =', config_dynamic_thickness
            endif

            ! Make sure config_calving_thickness > config_dynamic_thickness.
            ! Otherwise the algorithm will not work.

            if (config_calving_thickness < config_dynamic_thickness) then
               write(stderrUnit,*) 'ERROR: Must have config_calving_thickness > config_dynamic_thickness'
               write(stderrUnit,*) 'config_calving_thickness (m) =', config_calving_thickness
               write(stderrUnit,*) 'config_dynamic_thickness (m) =', config_dynamic_thickness
!!               call mpas_dmpar_global_abort('MPAS-landice: Aborting with calving error')
            endif

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

            if (config_print_calving_info) then

!               write(stdoutUnit,*) 'Active-for-calving mask'
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                     write(stdoutUnit,'(a3)',advance='no') '    '
                  endif
!!               do i = 1, nCellsPerRow
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!!                  write(stdoutUnit,'(i5)',advance='no') iCell
!                     write(stdoutUnit,'(i8)',advance='no') activeForCalvingMask(iCell)
                  enddo
!                  write(stdoutUnit,*) ' '
               enddo

!               write(stdoutUnit,*) 'Inactive margin mask'
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
!                     write(stdoutUnit,'(a3)',advance='no') '    '
                  endif
!!               do i = 1, nCellsPerRow
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!!                  write(stdoutUnit,'(i5)',advance='no') iCell
!                     write(stdoutUnit,'(i8)',advance='no') inactiveMarginMask(iCell)
                  enddo
!                  write(stdoutUnit,*) ' '
               enddo

            endif  ! config_print_calving_info

            ! Calve ice in ocean cells that are not on the protected inactive margin

            where (oceanMask == 1 .and. inactiveMarginMask == 0 .and. thickness > 0.0_RKIND)
               calvingThickness = thickness * calvingFraction
            endwhere

            ! clean up
            call mpas_deallocate_scratch_field(activeForCalvingMaskField, .true.)
            call mpas_deallocate_scratch_field(inactiveMarginMaskField, .true.)
            call mpas_deallocate_scratch_field(oceanMaskField, .true.)

         else   ! other calving options (floating, topographic_threshold)

            ! get scratch fields for calving

            call mpas_pool_get_field(scratchPool, 'iceCellMask',  calvingLawMaskField)
            call mpas_allocate_scratch_field(calvingLawMaskField, .true.)
            calvingLawMask => calvingLawMaskField % array

            if (trim(config_calving) == 'floating') then

               ! calculate masks - so we know where the ice is floating
               call li_calculate_mask(meshPool, velocityPool, geometryPool, err_tmp)
               err = ior(err, err_tmp)

               ! calve floating ice
               ! Note: The floating_ice mask includes all floating ice, both inactive and active

               where (li_mask_is_floating_ice(cellMask))
                  calvingLawMask = 1
               elsewhere
                  calvingLawMask = 0
               endwhere

            elseif (trim(config_calving) == 'topographic_threshold') then

               ! calve ice where the bed topography lies below a threshold depth

               if (config_print_calving_info) then
                  write(stdoutUnit,*) 'Calving topographic threshold (m) =', config_calving_topography
               endif

               where (bedTopography < config_calving_topography + config_sea_level)
                  calvingLawMask = 1
               elsewhere
                  calvingLawMask = 0
               endwhere

               !Note: Could add other calving options here (e.g., damage-based calving)

            endif   ! floating or topographic_threshold

            ! Calve ice where specified by the calving law.
            ! If the calving timescale > 0, then will calve only a fraction of the thickness in each calving cell.

            where (calvingLawMask == 1)

               calvingThickness = thickness * calvingFraction

            endwhere

            call mpas_deallocate_scratch_field(calvingLawMaskField, .true.)

         endif   ! config_calving

         ! Now that calvingThickness has been calculated, apply it to the thickness field,
         ! unless we are in a data calving mode.  In data calving mode we just calculate
         ! what should be calved but don't actually calve it.
         if (.not. config_data_calving) then
            thickness = thickness - calvingThickness
         endif


         ! Optionally, print a list of cells with calving

         if (config_print_calving_info) then

            write(stdoutUnit,*) ' '
            write(stdoutUnit,*) 'Global cell ID, bedTopography, calvingThickness:'
            do iCell = 1, nCells
               if (calvingThickness(iCell) > 0.0_RKIND) then
                  write(stdoutUnit,*) indexToCellID(iCell), bedTopography(iCell), calvingThickness(iCell)
               endif
            enddo

            !WHL - debug - for circular shelf test case
            if (circular_shelf) then
               write(stdoutUnit,*) ' '
               write(stdoutUnit,*) 'Ice thickness after calving'
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                     write(stdoutUnit,'(a3)',advance='no') '    '
                  endif
!!                 do i = 1, nCellsPerRow
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!!                    write(stdoutUnit,'(i5)',advance='no') iCell
                     write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
                  enddo
                  write(stdoutUnit,*) ' '
               enddo
            endif   ! circular_shelf

         endif   ! config_print_calving_info

         block => block % next
      enddo


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_calve_ice."
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

      type (dm_info), pointer :: dminfo
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
            write(stdoutUnit,*) 'Restore calving front'
            write(stdoutUnit,*) 'max thickness (m) =', maxval(thickness)

            !WHL - debug - for circular shelf test case
            if (circular_shelf) then
               write(stdoutUnit,*) 'Initial ice thickness'
               do iRow = nRows, 1, -1
                  if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                     write(stdoutUnit,'(a3)',advance='no') '    '
                  endif
!!                 do i = 1, nCellsPerRow
                  do i = nCellsPerRow/2 - 2, nCellsPerRow
                     iCell = (iRow-1)*nCellsPerRow + i
!!                    write(stdoutUnit,'(i5)',advance='no') iCell
                     write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
                  enddo
                  write(stdoutUnit,*) ' '
               enddo
            endif   ! circular_shelf

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
                     write(stdoutUnit,*) 'Restore ice: iCell, thickness =', iCell, thickness(iCell)
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
                     write(stdoutUnit,*) 'Remove ice:  iCell, thickness =', iCell, thickness(iCell)
                  endif

                  calvingThickness(iCell) = thickness(iCell)
                  thickness(iCell) = 0.0_RKIND

               endif   ! li_mask_is_initial_ice

            endif    ! bedTopography < config_sea_level

         enddo   ! iCell

         block => block % next
      enddo

      if (config_print_calving_info) then
         write(stdoutUnit,*) 'Restored the initial calving front'

         !WHL - debug - for circular shelf test case
         if (circular_shelf) then
            write(stdoutUnit,*) 'Final ice thickness'
            do iRow = nRows, 1, -1
               if (mod(iRow,2) == 0) then  ! indent for even-numbered rows
                  write(stdoutUnit,'(a3)',advance='no') '    '
               endif
!!               do i = 1, nCellsPerRow
               do i = nCellsPerRow/2 - 2, nCellsPerRow
                  iCell = (iRow-1)*nCellsPerRow + i
!!                  write(stdoutUnit,'(i5)',advance='no') iCell
                  write(stdoutUnit,'(f8.2)',advance='no') thickness(iCell)
               enddo
               write(stdoutUnit,*) ' '
            enddo
         endif  ! circular shelf

      endif

      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_restore_calving_front."
      endif


    end subroutine li_restore_calving_front

!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************

! No private subroutines so far

end module li_calving


