! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_statistics
!
!> \MPAS land ice global and local statistics
!> \author William Lipscomb
!> \date   22 January 2015
!> \details
!>  This module contains routines for computing glocal and local
!>   statistics and other diagnostic info.
!>  It is based on a similar module in CISM.
!
!-----------------------------------------------------------------------

module li_statistics

   use mpas_derived_types
   use mpas_constants
   use mpas_dmpar
   use mpas_timer
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

   public :: li_compute_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************
   contains

!***********************************************************************
!
!  routine li_compute_statistics
!
!> \brief   Computes global and local statistics
!> \author  William Lipscomb
!> \date    22 January 2015
!> \details
!>  This routine computes global statistics for the full domain, along with
!>  local statistics and diagnostic info for a user-specified grid cell.
!>
!-----------------------------------------------------------------------

   subroutine li_compute_statistics(domain, itimestep)

     use li_diagnostic_vars, only: li_calculate_layerThickness

     implicit none

     ! Input/output arguments
     type (domain_type), intent(inout) :: domain !< Input/Output: domain object
     integer, intent(in) :: itimestep   !< Input: current time step counter

     ! Local variables

     type (block_type),     pointer :: block
     type (dm_info),        pointer :: dminfo

     ! pools
     type (mpas_pool_type), pointer :: meshPool
     type (mpas_pool_type), pointer :: geometryPool
     type (mpas_pool_type), pointer :: velocityPool
     type (mpas_pool_type), pointer :: thermalPool
     type (mpas_pool_type), pointer :: scratchPool

     ! mesh dimensions
     integer, pointer :: nVertLevels, nVertInterfaces, nCellsSolve, nEdgesSolve

     ! mesh arrays
     integer, dimension(:),   pointer :: indexToCellID, indexToEdgeID
     integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
     real (kind=RKIND), dimension(:), pointer ::  areaCell
     real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
     real (kind=RKIND), dimension(:), pointer :: bedTopography, sfcMassBal

     ! state variables
     character (len=StrKIND),             pointer :: xtime
     integer,           dimension(:),     pointer :: cellMask
     integer,           dimension(:),     pointer :: edgeMask
     real (kind=RKIND), dimension(:),     pointer :: upperSurface
     real (kind=RKIND), dimension(:),     pointer :: thickness
     real (kind=RKIND), dimension(:),     pointer :: surfaceTemperature, basalTemperature
     real (kind=RKIND), dimension(:,:),   pointer :: layerThickness, normalVelocity
     real (kind=RKIND), dimension(:,:),   pointer :: temperature

     ! scratch variables
     type (field1dInteger), pointer :: iceCellMaskField
     integer, dimension(:), pointer :: iceCellMask

     type (field1dInteger), pointer :: iceEdgeMaskField
     integer, dimension(:), pointer :: iceEdgeMask

     type (field2dReal), pointer :: workLevelCellField
     real (kind=RKIND), dimension(:,:), pointer :: workLevelCell

     ! config variables
     real (kind=RKIND), pointer :: rhoi    ! ice density (kg/m^3)
     integer, pointer :: statsCellID       ! global ID of cell for which we write stats/diagnostics

     ! other pointers
!     integer, pointer :: indexTemperature

     ! work variables and arrays
     real (kind=RKIND) :: localSum, localMin, localMax, localVertSumMin, localVertSumMax
     real (kind=RKIND) :: localAbsoluteMax
     integer :: localMinlocElement, localMinlocLevel
     integer :: localMaxlocElement, localMaxlocLevel
     integer :: localVertSumMinlocElement, localVertSumMaxlocElement
     real (kind=RKIND), dimension(:), allocatable :: workLevel

     ! sums and max/mins on local processor
     real (kind=RKIND) :: iceAreaSum, iceVolumeSum, iceEnergySum
     real (kind=RKIND) :: thicknessMax, thicknessMin
     real (kind=RKIND) :: temperatureMax, temperatureMin
     real (kind=RKIND) :: velocityMax, basalVelocityMax

     ! cells, edges and levels where max/min values are located (cell/edge indices are global)
     integer :: thicknessMinlocCell, thicknessMaxlocCell
     integer :: temperatureMinlocCell, temperatureMinlocLevel
     integer :: temperatureMaxlocCell, temperatureMaxlocLevel
     integer :: velocityMaxlocEdge, velocityMaxlocLevel
     integer :: basalVelocityMaxlocEdge

     ! global sums and max/mins
     real (kind=RKIND) :: globalIceAreaSum, globalIceVolumeSum, globalIceEnergySum
     real (kind=RKIND) :: globalThicknessMax, globalThicknessMin, globalThicknessMean
     real (kind=RKIND) :: globalTemperatureMax, globalTemperatureMin, globalTemperatureMean
     real (kind=RKIND) :: globalVelocityMax, globalBasalVelocityMax

     ! diagnostic info for user-specified grid cell
     integer :: diagnosticCellLocal, diagnosticBlockIDLocal, diagnosticProcIDLocal
     real (kind=RKIND) :: diagnosticUpperSurfaceLocal, diagnosticThicknessLocal, diagnosticBedTopographyLocal
     real (kind=RKIND) :: diagnosticSfcMassBalLocal, diagnosticSurfaceTemperatureLocal, diagnosticBasalTemperatureLocal

     integer :: diagnosticCell, diagnosticBlockID, diagnosticProcID
     real (kind=RKIND) :: diagnosticUpperSurface, diagnosticThickness, diagnosticBedTopography
     real (kind=RKIND) :: diagnosticSfcMassBal, diagnosticSurfaceTemperature, diagnosticBasalTemperature
     real (kind=RKIND), dimension(:), allocatable :: diagnosticSpeed, diagnosticTemperature

     integer :: iCell, iCell1, iCell2, iEdge, kLevel
     integer :: proc

     block => domain % blocklist
     dminfo => domain % dminfo

     ! initialize statistics for this processor
     ! If we have > 1 block/proc, these are summed or reduced over the block

     iceAreaSum   =  0.0_RKIND
     iceVolumeSum =  0.0_RKIND
     iceEnergySum =  0.0_RKIND
     thicknessMax =     -huge(0.0_RKIND)
     thicknessMin =      huge(0.0_RKIND)
     temperatureMax =   -huge(0.0_RKIND)
     temperatureMin =    huge(0.0_RKIND)
     velocityMax =      -huge(0.0_RKIND)
     basalVelocityMax = -huge(0.0_RKIND)

     thicknessMinlocCell = 0
     thicknessMaxlocCell = 0
     temperatureMinlocCell = 0
     temperatureMinlocLevel = 0
     velocityMaxlocEdge = 0
     velocityMaxlocLevel = 0
     basalVelocityMaxlocEdge = 0

     ! initialize info for diagnostic grid cell
     ! These values will be overwritten on the processor owning this grid cell
     diagnosticCell = 0
     diagnosticBlockID = 0
     diagnosticProcID = 0
     diagnosticUpperSurface = 0.0_RKIND
     diagnosticThickness = 0.0_RKIND
     diagnosticBedTopography = 0.0_RKIND
     diagnosticSfcMassBal = 0.0_RKIND
     diagnosticSurfaceTemperature = 0.0_RKIND
     diagnosticBasalTemperature = 0.0_RKIND

     do while (associated(block))

        ! pools

        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
        call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
        call mpas_pool_get_subpool(block % structs, 'thermal', thermalPool)
        call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

        ! mesh dimensions
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        ! mesh arrays
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
        call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

        ! LI mesh arrays
        call mpas_pool_get_array(meshPool, 'xtime', xtime)
        call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)

        ! Geometry arrays
        call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
        call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
        call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
        call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
        call mpas_pool_get_array(geometryPool, 'thickness', thickness, timeLevel = 1)
        call mpas_pool_get_array(geometryPool, 'layerThickness', layerThickness, timeLevel = 1)
        call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)

        ! scratch arrays
        call mpas_pool_get_field(scratchPool, 'iceCellMask', iceCellMaskField)
        call mpas_pool_get_field(scratchPool, 'iceEdgeMask', iceEdgeMaskField)
        call mpas_pool_get_field(scratchPool, 'workLevelCell', workLevelCellField)

        ! velocity arrays
        call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)

        ! thermal arrays
        call mpas_pool_get_array(thermalPool, 'temperature', temperature)
        call mpas_pool_get_array(thermalPool, 'surfaceTemperature', surfaceTemperature)
        call mpas_pool_get_array(thermalPool, 'basalTemperature', basalTemperature)

        ! config settings
        call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
        call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', statsCellID)

        call mpas_allocate_scratch_field(iceCellMaskField, .true.)
        iceCellMask => iceCellMaskField % array
        call mpas_allocate_scratch_field(iceEdgeMaskField, .true.)
        iceEdgeMask => iceEdgeMaskField % array
        call mpas_allocate_scratch_field(workLevelCellField, .true.)
        workLevelCell => workLevelCellField % array

        ! compute ice cell mask (= 1 for cells where ice is present, else = 0)
        ! Note: Global sums are taken only over cells with mask = 1

        where (li_mask_is_ice(cellMask))
           iceCellMask = 1
        elsewhere
           iceCellMask = 0
        endwhere

        ! compute ice edge mask (= 1 for edges of cells where ice is present, else = 0)

        where (li_mask_is_ice(edgeMask))
           iceEdgeMask = 1
        elsewhere
           iceEdgeMask = 0
        endwhere

        ! given thickness, compute layerThickness
        call li_calculate_layerThickness(meshPool, thickness, layerThickness)

        ! Compute statistics on local block

        ! max and min ice thickness
        ! optional maxloc/minloc arguments give the local cell IDs where max/mins are located

        call li_compute_field_local_stats(dminfo,                           &
                                          nVertLevels,     nCellsSolve,     &
                                          layerThickness,                   &
                                          iceCellMask(1:nCellsSolve),       &
                                          localSum,                         &
                                          localMin,        localMax,        &
                                          localVertSumMin, localVertSumMax, &
                                          localVertSumMinlocElement = localVertSumMinlocElement, &
                                          localVertSumMaxlocElement = localVertSumMaxlocElement)

        if (localVertSumMin < thicknessMin) then
           thicknessMin = localVertSumMin
           thicknessMinlocCell = indexToCellID(localVertSumMinlocElement)
        endif

        if (localVertSumMax > thicknessMax) then
           thicknessMax = localVertSumMax
           thicknessMaxlocCell = indexToCellID(localVertSumMaxlocElement)
        endif

        ! work array whose vertical sum = 1 everywhere
        workLevelCell(1,:) = 1.0_RKIND
        workLevelCell(2:nVertLevels,:) = 0.0_RKIND

        ! total ice area
        call li_compute_field_area_weighted_local_stats                       &
                                          (dminfo,                            &
                                           nVertLevels,       nCellsSolve,    &
                                           areaCell(1:nCellsSolve),           &
                                           workLevelCell(:,1:nCellsSolve),    &
                                           iceCellMask(1:nCellsSolve),        &
                                           localSum,                          &
                                           localMin,          localMax,       &
                                           localVertSumMin,   localVertSumMax)

        iceAreaSum = iceAreaSum + localSum

        ! work array with a value of 1 in each layer
        workLevelCell(:,:) = 1.0_RKIND

        ! total ice volume
        call li_compute_field_volume_weighted_local_stats                                 &
                                          (dminfo,                                        &
                                           nVertLevels,                     nCellsSolve,  &
                                           areaCell(1:nCellsSolve),                       &
                                           layerThickness(:,1:nCellsSolve),               &
                                           workLevelCell(:,1:nCellsSolve),                &
                                           iceCellMask(1:nCellsSolve),                    &
                                           localSum,                                      &
                                           localMin,                        localMax,     &
                                           localVertSumMin,                 localVertSumMax)

        iceVolumeSum = iceVolumeSum + localSum

        ! max and min temperature
        ! optional maxloc/minloc arguments give the local cell IDs where max/mins are located

        call li_compute_field_local_stats(dminfo,                                    &
                                          nVertLevels,        nCellsSolve,           &
                                          temperature(:,1:nCellsSolve),              &
                                          iceCellMask(1:nCellsSolve),                &
                                          localSum,                                  &
                                          localMin,           localMax,              &
                                          localVertSumMin,    localVertSumMax,       &
                                          localMinlocElement, localMinlocLevel,      &
                                          localMaxlocElement, localMaxlocLevel)

        if (localMin < temperatureMin) then
           temperatureMin = localMin
           temperatureMinlocCell  = indexToCellID(localMinlocElement)
           temperatureMinlocLevel = localMinlocLevel
        endif

        if (localMax > temperatureMax) then
           temperatureMax = localMax
           temperatureMaxlocCell  = indexToCellID(localMaxlocElement)
           temperatureMaxlocLevel = localMaxlocLevel
        endif

        ! total ice energy (relative to 0 deg C)
        !TODO - Compute ice energy differently if using the enthalpy scheme

        call li_compute_field_volume_weighted_local_stats                                 &
                                          (dminfo,                                        &
                                           nVertLevels,                     nCellsSolve,  &
                                           areaCell(1:nCellsSolve),                       &
                                           layerThickness(:,1:nCellsSolve),               &
                                           temperature(:,1:nCellsSolve),                  &
                                           iceCellMask(1:nCellsSolve),                    &
                                           localSum,                                      &
                                           localMin,                        localMax,     &
                                           localVertSumMin,                 localVertSumMax)

        iceEnergySum = iceEnergySum + localSum*rhoi*cp_ice

        ! normal velocity at cell edges; find the maximum magnitude

        call li_compute_field_local_stats(dminfo,                               &
                                          nVertInterfaces,        nEdgesSolve,  &
                                          normalVelocity(:,1:nEdgesSolve),      &
                                          iceEdgeMask(1:nEdgesSolve),           &
                                          localSum,                             &
                                          localMin,           localMax,         &
                                          localVertSumMin,    localVertSumMax,  &
                                          localMinlocElement, localMinlocLevel, &
                                          localMaxlocElement, localMaxlocLevel)

        localAbsoluteMax = max(localMax, -localMin)   ! take maximum magnitude, independent of sign

        if (localAbsoluteMax > velocityMax) then
           velocityMax = localAbsoluteMax
           if (localMax > abs(localMin)) then
              velocityMaxlocEdge = indexToEdgeID(localMaxlocElement)
              velocityMaxlocLevel = localMaxlocLevel
           else
              velocityMaxlocEdge = indexToEdgeId(localMinlocElement)
              velocityMaxlocLevel = localMinlocLevel
           endif
        endif

        ! basal velocity at cell edges; find the maximum magnitude
        ! Note: If velocity is located at layer midpoints, this is actually the
        !       velocity in the lowest layer

        call li_compute_field_local_stats(dminfo,                                    &
                                          1,                  nEdgesSolve,           &
                                          normalVelocity(nVertInterfaces,1:nEdgesSolve), &
                                          iceEdgeMask(1:nEdgesSolve),                &
                                          localSum,                                  &
                                          localMin,           localMax,              &
                                          localVertSumMin,    localVertSumMax,       &
                                          localMinlocElement, localMinlocLevel,      &
                                          localMaxlocElement, localMaxlocLevel)

        localAbsoluteMax = max(localMax, -localMin)   ! take maximum magnitude, independent of sign

        if (localAbsoluteMax > basalVelocityMax) then
           basalVelocityMax = localAbsoluteMax
           if (localMax > abs(localMin)) then
              basalVelocityMaxlocEdge = indexToEdgeID(localMaxlocElement)
           else
              basalVelocityMaxlocEdge = indexToEdgeID(localMinlocElement)
           endif
        endif

        ! allocate and initialize some diagnostic arrays if not done already
        if (.not. allocated(diagnosticSpeed)) then
           allocate(diagnosticSpeed(nVertInterfaces))
           diagnosticSpeed(:) = 0.0_RKIND
        endif

        if (.not. allocated(diagnosticTemperature)) then
           allocate(diagnosticTemperature(nVertLevels))
           diagnosticTemperature(:) = 0.0_RKIND
        endif

        ! Determine whether the user-specified diagnostic cell is on this block
        ! If so, then set some diagnostics to be broadcast later to the head processor

        do iCell = 1, nCellsSolve
           if (indexToCellId(iCell) == statsCellID) then  ! this is the diagnostic cell
              diagnosticCell = iCell
              diagnosticBlockID = block % localBlockID
              diagnosticProcID = dminfo % my_proc_id
              diagnosticUpperSurface = upperSurface(iCell)
!!              diagnosticThickness = sum(layerThickness(:,iCell))
              diagnosticThickness = thickness(iCell)
              diagnosticBedTopography = bedTopography(iCell)
              diagnosticSfcMassBal = sfcMassBal(iCell) * scyr / 1000.0_RKIND   ! convert from kg/m^2/s to m/yr
              diagnosticSurfaceTemperature = surfaceTemperature(iCell)
              diagnosticBasalTemperature = basalTemperature(iCell)
              iEdge = edgesOnCell(1,iCell)  ! arbitrarily choose edge #1 for velocity diagnostics
                                            ! alternatively, could write out the cell-center velocity
              diagnosticSpeed(:) = normalVelocity(:,iEdge) * scyr  ! convert from m/s to m/yr
              diagnosticTemperature(:) = temperature(:,iCell)
           endif
        enddo

        ! clean up
        call mpas_deallocate_scratch_field(iceCellMaskField, .true.)
        call mpas_deallocate_scratch_field(iceEdgeMaskField, .true.)
        call mpas_deallocate_scratch_field(workLevelCellField, .true.)

        block => block % next
     enddo  ! block loop

     ! Compute global statistics
     ! TODO: Reduce the number of global reductions in this subroutine?
     !       This could be done by packing quantities into arrays.

     ! global sums
     call mpas_dmpar_sum_real(dminfo, iceAreaSum, globalIceAreaSum)
     call mpas_dmpar_sum_real(dminfo, iceVolumeSum, globalIceVolumeSum)
     call mpas_dmpar_sum_real(dminfo, iceEnergySum, globalIceEnergySum)

     ! global means
     !TODO - Replace temperature mean with enthalpy mean if using enthalpy scheme

     if (globalIceAreaSum > 0.0_RKIND) then
        globalThicknessMean = globalIceVolumeSum / globalIceAreaSum
     else
        globalThicknessMean = 0.0_RKIND
     endif

     if (globalIceVolumeSum > 0.0_RKIND) then
        globalTemperatureMean = globalIceEnergySum / (globalIceVolumeSum * rhoi * cp_ice)
     else
        globalTemperatureMean = 0.0_RKIND
     endif

     ! global max/mins of state variables
     ! First determine the global max/min and the proc on which it resides
     ! Then broadcast the global max/min and its cell/edge/level location to all processors

     call mpas_dmpar_minloc_real(dminfo, thicknessMin, globalThicknessMin, proc)
     call mpas_dmpar_bcast_real (dminfo, globalThicknessMin, proc)
     call mpas_dmpar_bcast_int  (dminfo, thicknessMinlocCell, proc)

     call mpas_dmpar_maxloc_real(dminfo, thicknessMax, globalThicknessMax, proc)
     call mpas_dmpar_bcast_real (dminfo, globalThicknessMax, proc)
     call mpas_dmpar_bcast_int  (dminfo, thicknessMaxlocCell, proc)

     call mpas_dmpar_minloc_real(dminfo, temperatureMin, globalTemperatureMin, proc)
     call mpas_dmpar_bcast_real (dminfo, globalTemperatureMin, proc)
     call mpas_dmpar_bcast_int  (dminfo, temperatureMinlocCell, proc)
     call mpas_dmpar_bcast_int  (dminfo, temperatureMinlocLevel, proc)

     call mpas_dmpar_maxloc_real(dminfo, temperatureMax, globalTemperatureMax, proc)
     call mpas_dmpar_bcast_real (dminfo, globalTemperatureMax, proc)
     call mpas_dmpar_bcast_int  (dminfo, temperatureMaxlocCell, proc)
     call mpas_dmpar_bcast_int  (dminfo, temperatureMaxlocLevel, proc)

     call mpas_dmpar_maxloc_real(dminfo, velocityMax, globalVelocityMax, proc)
     call mpas_dmpar_bcast_real (dminfo, globalVelocityMax, proc)
     call mpas_dmpar_bcast_int  (dminfo, velocityMaxlocEdge, proc)
     call mpas_dmpar_bcast_int  (dminfo, velocityMaxlocLevel, proc)

     call mpas_dmpar_maxloc_real(dminfo, basalVelocityMax, globalBasalVelocityMax, proc)
     call mpas_dmpar_bcast_real (dminfo, globalBasalVelocityMax, proc)
     call mpas_dmpar_bcast_int  (dminfo, basalVelocityMaxlocEdge, proc)

     ! global reductions for user-specified diagnostic cell
     ! Note: These reductions are done with global sums rather than broadcasts.
     !       Global sums will work provided that the quantity of interest
     !        has nonzero values on only a single processor.
     diagnosticCellLocal = diagnosticCell
     diagnosticBlockIDLocal = diagnosticBlockID
     diagnosticProcIDLocal = diagnosticProcID

     diagnosticUpperSurfaceLocal = diagnosticUpperSurface
     diagnosticThicknessLocal = diagnosticThickness
     diagnosticBedTopographyLocal = diagnosticBedTopography
     diagnosticSfcMassBalLocal = diagnosticSfcMassBal
     diagnosticSurfaceTemperatureLocal = diagnosticSurfaceTemperature
     diagnosticBasalTemperatureLocal = diagnosticBasalTemperature

     call mpas_dmpar_sum_int  (dminfo, diagnosticCellLocal, diagnosticCell)
     call mpas_dmpar_sum_int  (dminfo, diagnosticBlockIDLocal, diagnosticBlockID)
     call mpas_dmpar_sum_int  (dminfo, diagnosticProcIDLocal, diagnosticProcID)

     call mpas_dmpar_sum_real (dminfo, diagnosticUpperSurfaceLocal, diagnosticUpperSurface)
     call mpas_dmpar_sum_real (dminfo, diagnosticThicknessLocal, diagnosticThickness)
     call mpas_dmpar_sum_real (dminfo, diagnosticBedTopographyLocal, diagnosticBedTopography)
     call mpas_dmpar_sum_real (dminfo, diagnosticSfcMassBalLocal, diagnosticSfcMassBal)
     call mpas_dmpar_sum_real (dminfo, diagnosticSurfaceTemperatureLocal, diagnosticSurfaceTemperature)
     call mpas_dmpar_sum_real (dminfo, diagnosticBasalTemperatureLocal, diagnosticBasalTemperature)

     allocate (workLevel(nVertInterfaces))
     call mpas_dmpar_sum_real_array(dminfo, nVertInterfaces, diagnosticSpeed, workLevel)
     diagnosticSpeed(:) = workLevel(:)
     deallocate(workLevel)

     allocate (workLevel(nVertLevels))
     call mpas_dmpar_sum_real_array(dminfo, nVertLevels, diagnosticTemperature, workLevel)
     diagnosticTemperature(:) = workLevel(:)
     deallocate(workLevel)

     ! Write global and local stats to the log file

     if (dminfo % my_proc_id == IO_NODE) then
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a60)') '------------------------------------------------------------'
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a25,a20)')          'Global statistics: time =', trim(xtime)
        write(stdoutUnit,'(a25,i8)')           '               timestep =', itimestep
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a32,e24.16)')       'Total ice area (km^2)           ',   &
                                       globalIceAreaSum*1.0e-6_RKIND       ! convert from m^2 to km^2
        write(stdoutUnit,'(a32,e24.16)')       'Total ice volume (km^3)         ',   &
                                       globalIceVolumeSum*1.0e-9_RKIND     ! convert from m^3 to km^3
        write(stdoutUnit,'(a32,e24.16)')       'Total ice energy (J)            ',   &
                                       globalIceEnergySum
        write(stdoutUnit,'(a32,f24.16,i8)')    'Max thickness (m), cell         ',   &
                                       globalThicknessMax, thicknessMaxlocCell
        write(stdoutUnit,'(a32,f24.16,i8)')    'Min thickness (m), cell         ',   &
                                       globalThicknessMin, thicknessMinlocCell
        write(stdoutUnit,'(a32,f24.16)')       'Mean thickness (m)              ',   &
                                       globalthicknessMean
        write(stdoutUnit,'(a32,f24.16,i8,i4)') 'Max temperature (C), cell, level',   &
                                       globalTemperatureMax - kelvin_to_celsius, temperatureMaxlocCell, temperatureMaxlocLevel
        write(stdoutUnit,'(a32,f24.16,i8,i4)') 'Min temperature (C), cell, level',   &
                                       globalTemperatureMin - kelvin_to_celsius , temperatureMinlocCell, temperatureMinlocLevel
        write(stdoutUnit,'(a32,f24.16)')       'Mean temperature (C)            ',   &
                                       globalTemperatureMean - kelvin_to_celsius
        write(stdoutUnit,'(a32,f24.16,i8,i4)') 'Max velocity (m/yr), edge, level',   &
                                       globalVelocityMax * scyr, velocityMaxlocEdge, velocityMaxlocLevel
        write(stdoutUnit,'(a32,f24.16,i8,i4)') 'Max basal velo (m/yr), edge     ',   &
                                       globalBasalVelocityMax * scyr, basalVelocityMaxlocEdge
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a30,i6)')  'Column diagnostics: cell ID = ', statsCellID
        write(stdoutUnit,'(a30,3i6)') 'Local cell ID, block, proc =  ', diagnosticCell, diagnosticBlockID, diagnosticProcID
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a25,f24.16)') 'Upper surface (m)        ', diagnosticUpperSurface
        write(stdoutUnit,'(a25,f24.16)') 'Thickness (m)            ', diagnosticThickness
        write(stdoutUnit,'(a25,f24.16)') 'Bed topography (m)       ', diagnosticBedTopography
        write(stdoutUnit,'(a25,f24.16)') 'Sfc mass balance (m/yr)  ', diagnosticSfcMassBal
        write(stdoutUnit,*) ' '
        write(stdoutUnit,'(a55)') 'Sigma        Ice speed (m/yr)       Ice temperature (C)'
        write(stdoutUnit,'(f6.4, a25, f24.16)') 0.0_RKIND, '------', diagnosticSurfaceTemperature - kelvin_to_celsius
        do kLevel = 1, nVertLevels
           write(stdoutUnit,'(f6.4, f25.16, f24.16)') &
                layerCenterSigma(kLevel), diagnosticSpeed(kLevel), diagnosticTemperature(kLevel) - kelvin_to_celsius
        end do
        write(stdoutUnit,'(f6.4, a25, f24.16)') 1.0_RKIND, '------', diagnosticBasalTemperature - kelvin_to_celsius
        write(stdoutUnit,*) ' '
     endif   ! my_proc_id = IO_NODE

     ! clean up
     deallocate(diagnosticTemperature)
     deallocate(diagnosticSpeed)

   end subroutine li_compute_statistics

!***********************************************************************
!
!  routine li_compute_field_local_stats
!
!> \brief   Computes statistics for a field on a single block
!> \author  William Lipscomb
!> \date    22 January 2015
!> \details
!>  This routine computes statistics (sum, max/min, vertical sum max/min)
!>  for a real array on a single block.
!
!-----------------------------------------------------------------------

   subroutine li_compute_field_local_stats(dminfo,                               &
                                           nVertLevels,        nElements,        &
                                           field,              mask,             &
                                           localSum,                             &
                                           localMin,           localMax,         &
                                           localVertSumMin,    localVertSumMax,  &
                                           localMinlocElement, localMinlocLevel, &
                                           localMaxlocElement, localMaxlocLevel, &
                                           localVertSumMinlocElement,            &
                                           localVertSumMaxlocElement)

     ! Compute field statistics without area or volume weighting

     implicit none

     ! Input/output arguments
     type (dm_info), intent(in) :: dminfo
     integer, intent(in) :: nVertLevels, nElements

     real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: &
          field    ! input field for which statistics are computed

     integer, dimension(nElements), intent(in) ::  &
          mask     ! = 0 or 1; compute stats only over region where mask = 1

     real (kind=RKIND), intent(out) :: localSum, localMin, localMax
     real (kind=RKIND), intent(out) :: localVertSumMin, localVertSumMax

     integer, intent(out), optional :: localMinlocElement, localMinlocLevel
     integer, intent(out), optional :: localMaxlocElement, localMaxlocLevel
     integer, intent(out), optional :: localVertSumMinlocElement
     integer, intent(out), optional :: localVertSumMaxlocElement

     ! Local variables
     integer :: i, k

     localSum = 0.0_RKIND
     do i = 1, nElements
        localSum = localSum + real(mask(i),RKIND) * sum(field(:,i))
     end do

     if (present(localMinlocElement) .and. present(localMinlocLevel)) then
        localMin = 1.0e34_RKIND
        localMinlocElement = 0
        localMinlocLevel = 0
        do i = 1, nElements
           do k = 1, nVertLevels
              if (field(k,i) < localMin) then
                 localMin = field(k,i)
                 localMinlocElement = i
                 localMinlocLevel = k
              endif
           enddo
        enddo
     else
        localMin = minval(field)
     endif

     if (present(localMaxlocElement) .and. present(localMaxlocLevel)) then
        localMax = -1.0e34_RKIND
        localMaxlocElement = 0
        localMaxlocLevel = 0
        do i = 1, nElements
           do k = 1, nVertLevels
              if (field(k,i) > localMax) then
                 localMax = field(k,i)
                 localMaxlocElement = i
                 localMaxlocLevel = k
              endif
           enddo
        enddo
     else
        localMax = maxval(field)
     endif

     if (present(localVertSumMinlocElement)) then
        localVertSumMin = 1.0e34_RKIND
        localVertSumMinlocElement = 0
        do i = 1, nElements
           if (sum(field(:,i)) < localVertSumMin) then
              localVertSumMin = sum(field(:,i))
              localVertSumMinlocElement = i
           endif
        enddo
     else
        localVertSumMin = minval(sum(field,1))
     endif

     if (present(localVertSumMaxlocElement)) then
        localVertSumMax = -1.0e34_RKIND
        localVertSumMaxlocElement = 0
        do i = 1, nElements
           if (sum(field(:,i)) > localVertSumMax) then
              localVertSumMax = sum(field(:,i))
              localVertSumMaxlocElement = i
           endif
        enddo
     else
        localVertSumMax = maxval(sum(field,1))
     endif

   end subroutine li_compute_field_local_stats

!***********************************************************************
!
!  routine li_compute_field_area_weighted_local_stats
!
!> \brief   Computes area-weighted statistics for a field on a single block
!> \author  William Lipscomb
!> \date    22 January 2015
!> \details
!>  This routine computes statistics (sum, max/min, vertical sum max/min)
!>  for a real array on a single block. The sum is weighted by the input
!>  field 'areas' (typically the grid cell area).
!
!-----------------------------------------------------------------------

   subroutine li_compute_field_area_weighted_local_stats(dminfo,                      &
                                                         nVertLevels,     nElements,  &
                                                         areas,                       &
                                                         field,           mask,       &
                                                         localSum,                    &
                                                         localMin,        localMax,   &
                                                         localVertSumMin, localVertSumMax)

     ! Compute field statistics weighted by area

     implicit none

     ! Input/output arguments
     type (dm_info), intent(in) :: dminfo
     integer, intent(in) :: nVertLevels, nElements

     real (kind=RKIND), dimension(nElements), intent(in) :: &
          areas    ! grid cell areas

     real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: &
          field    ! input field for which statistics are computed

     integer, dimension(nElements), intent(in) ::  &
          mask     ! = 0 or 1; compute stats only over region where mask = 1

     real (kind=RKIND), intent(out) :: localSum, localMin, localMax
     real (kind=RKIND), intent(out) :: localVertSumMin, localVertSumMax

     ! Local variables
     integer :: i

     localSum = 0.0_RKIND
     do i = 1, nElements
        localSum = localSum + real(mask(i),RKIND) * areas(i) * sum(field(:,i))
     end do

     localMin = minval(field)
     localMax = maxval(field)

     localVertSumMin = minval(sum(field,1))
     localVertSumMax = maxval(sum(field,1))

    end subroutine li_compute_field_area_weighted_local_stats

!***********************************************************************
!
!  routine li_compute_field_volume_weighted_local_stats
!
!> \brief   Computes volume-weighted statistics for a field on a single block
!> \author  William Lipscomb
!> \date    22 January 2015
!> \details
!>  This routine computes statistics (sum, max/min, vertical sum max/min)
!>  for a real array on a single block. The sum is weighted by the product
!>  of the input fields 'areas' (typically the grid cell area) and 'layerThickness'.
!
!-----------------------------------------------------------------------

    subroutine li_compute_field_volume_weighted_local_stats(dminfo,                          &
                                                            nVertLevels,     nElements,      &
                                                            areas,           layerThickness, &
                                                            field,           mask,           &
                                                            localSum,                        &
                                                            localMin,        localMax,       &
                                                            localVertSumMin, localVertSumMax)

      implicit none

      ! Input/output arguments
      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements

      real (kind=RKIND), dimension(nElements), intent(in) :: &
           areas    ! element areas

      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: &
           layerThickness    ! ice thickness in each layer

      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: &
           field    ! input field for which statistics are computed

      integer, dimension(nElements), intent(in) ::  &
           mask     ! = 0 or 1; compute stats only over region where mask = 1

      real (kind=RKIND), intent(out) :: localSum, localMin, localMax
      real (kind=RKIND), intent(out) :: localVertSumMin, localVertSumMax

      ! Local variables
      integer :: i

      localSum = 0.0_RKIND
      do i = 1, nElements
         localSum = localSum + real(mask(i),RKIND) * areas(i) * sum(layerThickness(:,i)*field(:,i))
      end do

      localMin = minval(field)
      localMax = maxval(field)

      localVertSumMin = minval(sum(layerThickness*field,1))
      localVertSumMax = maxval(sum(layerThickness*field,1))

   end subroutine li_compute_field_volume_weighted_local_stats

! The remaining code is from an older module by Matt Hoffman.
! Keeping it here for reference.
!=====================================================================================

!      ! 6. Write out your global stat to the file
!      if (dminfo % my_proc_id == IO_NODE) then
!         fileID = land_ice_get_free_unit()
!
!         if (config_write_initial_stats .and. (timeIndex == 0)) then
!             open(fileID, file='GlobalIntegrals.txt',STATUS='unknown')
!         elseif ( .not.(config_write_initial_stats) .and. (timeIndex/config_stats_interval == 1) ) then
!             open(fileID, file='GlobalIntegrals.txt',STATUS='unknown')
!         else
!             open(fileID, file='GlobalIntegrals.txt',POSITION='append')
!         endif
!!         write(fileID,'(1i0, 100es24.16)') timeIndex, timeIndex*dt, globalFluidThickness, &
!!             globalPotentialVorticity, globalPotentialEnstrophy, &
!!                        globalEnergy, globalCoriolisEnergyTendency, globalKineticEnergyTendency+globalPotentialEnergyTendency, &
!!                        globalKineticEnergy, globalPotentialEnergy
!
!      endif

!   integer function land_ice_get_free_unit()
!      implicit none

!      integer :: index
!      logical :: isOpened

!      land_ice_get_free_unit = 0
!      do index = 1,99
!         if((index /= 5) .and. (index /= 6)) then
!            inquire(unit = index, opened = isOpened)
!            if( .not. isOpened) then
!               land_ice_get_free_unit = index
!               return
!            end if
!         end if
!      end do
!   end function land_ice_get_free_unit

 end module li_statistics
