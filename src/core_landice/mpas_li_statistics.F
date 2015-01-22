
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

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_timer
   use li_setup
   use li_mask

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

   subroutine li_compute_statistics(domain, timeLevel, timeIndex)

     ! dminfo is the domain info needed for global communication
     ! state contains the state variables needed to compute global diagnostics
     ! grid contains the meta data about the grid
     ! timeIndex is the current time step counter

     implicit none

     ! Input/output arguments
     type (domain_type), intent(inout) :: domain !< Input/Output: domain information
     integer, intent(in) :: timeLevel
     integer, intent(in) :: timeIndex

     ! Local variables

     type (block_type),     pointer :: block
     type (dm_info),        pointer :: dminfo

     ! pools
     type (mpas_pool_type), pointer :: statePool
     type (mpas_pool_type), pointer :: meshPool
     type (mpas_pool_type), pointer :: diagnosticsPool

     ! mesh dimensions 
     integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve

     ! mesh arrays
     integer, dimension(:),   pointer :: indexToCellID, indexToEdgeID
     integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
     real (kind=RKIND), dimension(:), pointer ::  areaCell
     real (kind=RKIND), dimension(:), pointer :: layerCenterSigma
     real (kind=RKIND), dimension(:), pointer :: bedTopography, sfcMassBal  

     ! state variables
     character (len=StrKIND),             pointer :: xtime
     integer,           dimension(:),     pointer :: cellMask
     real (kind=RKIND), dimension(:),     pointer :: upperSurface
     real (kind=RKIND), dimension(:),     pointer :: surfaceTemperature, basalTemperature
     real (kind=RKIND), dimension(:,:),   pointer :: layerThickness, normalVelocity
     real (kind=RKIND), dimension(:,:,:), pointer :: tracers

     ! config variables
     real (kind=RKIND), pointer :: scyr    ! seconds per year
     real (kind=RKIND), pointer :: rhoi    ! ice density (kg/m^3)
     real (kind=RKIND), pointer :: shci    ! specific heat capacity of ice (J/deg/kg)
     integer, pointer :: statsCellID       ! global ID of cell for which we write stats/diagnostics

     ! masks
     integer, dimension(:), allocatable :: iceCellMask, iceEdgeMask

     ! work variables and arrays
     real (kind=RKIND) :: localSum, localMin, localMax, localVertSumMin, localVertSumMax 
     integer :: localMinlocElement, localMinlocLevel
     integer :: localMaxlocElement, localMaxlocLevel
     integer :: localVertSumMinlocElement, localVertSumMaxlocElement
     real (kind=RKIND), dimension(:), allocatable :: workArray1D
     real (kind=RKIND), dimension(:,:), allocatable :: workArray

     ! sums and max/mins on local block
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
     integer :: diagnosticCell, diagnosticBlockID, diagnosticProcID
     real (kind=RKIND) :: diagnosticUpperSurface, diagnosticThickness, diagnosticBedTopography
     real (kind=RKIND) :: diagnosticSfcMassBal, diagnosticSurfaceTemperature, diagnosticBasalTemperature
     real (kind=RKIND), dimension(:), allocatable :: diagnosticSpeed, diagnosticTemperature

     integer :: iCell, iCell1, iCell2, iEdge, iTracer, kLevel
     integer :: proc

     block => domain % blocklist
     dminfo => domain % dminfo

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
        call mpas_pool_get_subpool(block % structs, 'state', statePool)
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

        ! mesh dimensions
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)

        ! mesh arrays
        call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
        call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
        call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
        call mpas_pool_get_array(meshPool, 'layerCenterSigma', layerCenterSigma)
        !TODO - Make bedTopography and sfcMassBal state variables?
        call mpas_pool_get_array(meshPool, 'bedTopography', bedTopography)
        call mpas_pool_get_array(meshPool, 'sfcMassBal', sfcMassBal)
         
        ! state variables
        call mpas_pool_get_array(statePool, 'cellMask', cellMask, timeLevel)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
        call mpas_pool_get_array(statePool, 'upperSurface', upperSurface,   timeLevel)
        call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
        call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)
        call mpas_pool_get_array(statePool, 'surfaceTemperature', surfaceTemperature, timeLevel)
        call mpas_pool_get_array(statePool, 'basalTemperature', basalTemperature, timeLevel)
        !Note: xtime only has one time level, but stating it explicitly here to avoid confusion
        call mpas_pool_get_array(statePool, 'xtime', xtime, timeLevel=1)

        ! config settings
        call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
        call mpas_pool_get_config(liConfigs, 'config_ice_specific_heat', shci)
        call mpas_pool_get_config(liConfigs, 'config_seconds_per_year', scyr)
        call mpas_pool_get_config(liConfigs, 'config_stats_cell_ID', statsCellID)

        ! compute ice cell mask (= 1 for cells where ice is present, else = 0)
        ! Note: Global sums are taken only over cells with mask = 1

        allocate(iceCellMask(nCellsSolve))
        do iCell = 1, nCellsSolve
           if (li_mask_is_ice(cellMask(iCell))) then
              iceCellMask(iCell) = 1
           else
              iceCellMask(iCell) = 0
           endif
        enddo

        ! compute ice edge mask (= 1 for edges of cells where ice is present, else = 0)
        allocate(iceEdgeMask(nEdgesSolve))
        do iEdge = 1, nEdgesSolve
           iCell1 = cellsOnEdge(1,iEdge)                                                                                                          
           iCell2 = cellsOnEdge(2,iEdge)  
           if (li_mask_is_ice(cellMask(iCell1)) .or. li_mask_is_ice(cellMask(iCell2))) then
              iceEdgeMask(iEdge) = 1
           else
              iceEdgeMask(iEdge) = 0
           endif
        enddo

        ! Compute statistics on local block

        ! max and min ice thickness
        ! optional maxloc/minloc arguments give the local cell IDs where max/mins are located

        call li_compute_field_local_stats(dminfo,                           &
                                          nVertLevels,     nCellsSolve,     &
                                          layerThickness,  iceCellMask,     &
                                          localSum,                         &
                                          localMin,        localMax,        &
                                          localVertSumMin, localVertSumMax, &
                                          localVertSumMinlocElement = localVertSumMinlocElement, &
                                          localVertSumMaxlocElement = localVertSumMaxlocElement)

        thicknessMin = localVertSumMin
        thicknessMax = localVertSumMax

        ! global cell index for max/min values
        thicknessMinlocCell = indexToCellID(localVertSumMinlocElement)
        thicknessMaxlocCell = indexToCellID(localVertSumMaxlocElement)

        ! work array whose vertical sum = 1 everywhere

        if (.not. allocated(workArray)) allocate(workArray(nVertLevels,nCellsSolve))
        workArray(1,:) = 1.0_RKIND
        workArray(2:nVertLevels,:) = 0.0_RKIND

        ! total ice area
        call li_compute_field_area_weighted_local_stats                       &
                                          (dminfo,                            &
                                           nVertLevels,       nCellsSolve,    &
                                           areaCell(1:nCellsSolve),           &
                                           workArray,         iceCellMask,    &
                                           localSum,                          &
                                           localMin,          localMax,       &
                                           localVertSumMin,   localVertSumMax)

        iceAreaSum = localSum

        ! work array with a value of 1 in each layer
        workArray(:,:) = 1.0_RKIND

        ! total ice volume
        call li_compute_field_volume_weighted_local_stats                                 &
                                          (dminfo,                                        &
                                           nVertLevels,                     nCellsSolve,  &
                                           areaCell(1:nCellsSolve),                       &
                                           layerThickness(:,1:nCellsSolve),               & 
                                           workArray,                       iceCellMask,  &
                                           localSum,                                      &
                                           localMin,                        localMax,     &
                                           localVertSumMin,                 localVertSumMax)

        iceVolumeSum = localSum

        ! max and min temperature
        ! optional maxloc/minloc arguments give the local cell IDs where max/mins are located

        iTracer = 1    ! assume temperature is the first tracer array

        call li_compute_field_local_stats(dminfo,                                &
                                          nVertLevels,        nCellsSolve,       &
                                          tracers(iTracer,:,1:nCellsSolve),      &
                                          iceCellMask,                           &
                                          localSum,                              &
                                          localMin,           localMax,          &
                                          localVertSumMin,    localVertSumMax,   &
                                          localMinlocElement, localMinlocLevel,  &
                                          localMaxlocElement, localMaxlocLevel)

        temperatureMax = localMax
        temperatureMin = localMin

        ! global cell and level indices for max/min values
        temperatureMinlocCell  = indexToCellID(localMinlocElement)
        temperatureMinlocLevel = localMinlocLevel
        temperatureMaxlocCell  = indexToCellID(localMaxlocElement)
        temperatureMaxlocLevel = localMaxlocLevel

        ! total ice energy (relative to 0 deg C)
        !TODO - Compute ice energy differently if using the enthalpy scheme

        call li_compute_field_volume_weighted_local_stats                                 &
                                          (dminfo,                                        &
                                           nVertLevels,                     nCellsSolve,  &
                                           areaCell(1:nCellsSolve),                       &
                                           layerThickness(:,1:nCellsSolve),               & 
                                           tracers(iTracer,:,1:nCellsSolve),              &
                                           iceCellMask,                                   &
                                           localSum,                                      &
                                           localMin,                        localMax,     &
                                           localVertSumMin,                 localVertSumMax)

        iceEnergySum = localSum * rhoi * shci

        ! normal velocity at cell edges; find the maximum magnitude
        !TODO - If velocity is defined at layer interfaces, then nVertLevels -> nVertLevels + 1

        call li_compute_field_local_stats(dminfo,                               &
                                          nVertLevels,        nEdgesSolve,      &
                                          normalVelocity(:,1:nEdgesSolve),      &
                                          iceEdgeMask,                          &
                                          localSum,                             &
                                          localMin,           localMax,         &
                                          localVertSumMin,    localVertSumMax,  &
                                          localMinlocElement, localMinlocLevel, &
                                          localMaxlocElement, localMaxlocLevel)

        velocityMax = max(localMax, -localMin)

        ! global edge and level indices for max value
        if (localMax > abs(localMin)) then
           velocityMaxlocEdge = indexToEdgeID(localMaxlocElement)
           velocityMaxlocLevel = localMaxlocLevel
        else
           velocityMaxlocEdge = indexToEdgeId(localMinlocElement)
           velocityMaxlocLevel = localMinlocLevel
        endif

        ! basal velocity at cell edges; find the maximum magnitude
        ! Note: If velocity is located at layer midpoints, this is actually the 
        !       velocity in the lowest layer

        call li_compute_field_local_stats(dminfo,                              &
                                          1,               nEdgesSolve,        &
                                          normalVelocity(nVertLevels,1:nEdgesSolve), &
                                          iceEdgeMask,                         &
                                          localSum,                            &
                                          localMin,           localMax,        &
                                          localVertSumMin,    localVertSumMax, &
                                          localMinlocElement, localMinlocLevel,&
                                          localMaxlocElement, localMaxlocLevel)
         
        basalVelocityMax = max(localMax, -localMin)

        ! global edge and level indices for max value
        if (localMax > abs(localMin)) then
           basalVelocityMaxlocEdge = indexToEdgeID(localMaxlocElement)
        else
           basalVelocityMaxlocEdge = indexToEdgeID(localMinlocElement)
        endif

        ! allocate and initialize some diagnostic arrays if not done already 
        !TODO - If velocity is defined at layer interfaces, then nVertLevels -> nVertLevels + 1
        if (.not. allocated(diagnosticTemperature)) then
           allocate(diagnosticTemperature(nVertLevels))
           diagnosticTemperature(:) = 0.0_RKIND
        endif

        if (.not. allocated(diagnosticSpeed)) then
           allocate(diagnosticSpeed(nVertLevels))
           diagnosticSpeed(:) = 0.0_RKIND
        endif

        ! Determine whether the user-specified diagnostic cell is on this block
        ! If so, then set some diagnostics to be broadcast later to the head processor

        do iCell = 1, nCellsSolve
           if (indexToCellId(iCell) == statsCellID) then  ! this is the diagnostic cell
              diagnosticCell = iCell
              diagnosticBlockID = block % localBlockID
              diagnosticProcID = dminfo % my_proc_id
              diagnosticUpperSurface = upperSurface(iCell)
              diagnosticThickness = sum(layerThickness(:,iCell))
              diagnosticBedTopography = bedTopography(iCell)
              diagnosticSfcMassBal = sfcMassBal(iCell) * scyr / 1000.0_RKIND   ! convert from kg/m^2/s to m/yr
              diagnosticSurfaceTemperature = surfaceTemperature(iCell)
              diagnosticBasalTemperature = basalTemperature(iCell)
              iEdge = edgesOnCell(1,iCell)  ! arbitrarily choose edge #1 for velocity diagnostics
              diagnosticSpeed(:) = normalVelocity(:,iEdge) * scyr  ! convert from m/s to m/yr
              diagnosticTemperature(:) = tracers(1,:,iCell)  ! assume temperature is tracer #1
           endif
        enddo
        
        block => block % next
     end do  ! block loop

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
        globalTemperatureMean = globalIceEnergySum / (globalIceVolumeSum * rhoi * shci)
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
     
     call mpas_dmpar_sum_int  (dminfo, diagnosticCell, diagnosticCell)
     call mpas_dmpar_sum_int  (dminfo, diagnosticBlockID, diagnosticBlockID)
     call mpas_dmpar_sum_int  (dminfo, diagnosticProcID, diagnosticProcID)

     call mpas_dmpar_sum_real (dminfo, diagnosticUpperSurface, diagnosticUpperSurface)
     call mpas_dmpar_sum_real (dminfo, diagnosticThickness, diagnosticThickness)
     call mpas_dmpar_sum_real (dminfo, diagnosticBedTopography, diagnosticBedTopography)
     call mpas_dmpar_sum_real (dminfo, diagnosticSfcMassBal, diagnosticSfcMassBal)
     call mpas_dmpar_sum_real (dminfo, diagnosticSurfaceTemperature, diagnosticSurfaceTemperature)
     call mpas_dmpar_sum_real (dminfo, diagnosticBasalTemperature, diagnosticBasalTemperature)

     !TODO - Change to nVertLevels + 1 if velocity lives on layer interfaces
     allocate (workArray1d(nVertLevels))
     call mpas_dmpar_sum_real_array(dminfo, nVertLevels, diagnosticSpeed, workArray1d)
     diagnosticSpeed(:) = workArray1D(:)

     call mpas_dmpar_sum_real_array(dminfo, nVertLevels, diagnosticTemperature, workArray1D)
     diagnosticTemperature(:) = workArray1D(:)
     
     ! Write global and local stats to the log file
     !TODO - Change stdout (6) to another fileunit?
     
     if (dminfo % my_proc_id == IO_NODE) then
        write(6,*) ' '
        write(6,'(a60)') '------------------------------------------------------------'
        write(6,*) ' '
        write(6,'(a25,a20)')          'Global statistics: time =', trim(xtime)
        write(6,'(a25,i8)')           '               timestep =', timeIndex
        write(6,*) ' '
        write(6,'(a32,e24.16)')       'Total ice area (km^2)           ',   &
                                       globalIceAreaSum*1.0d-6       ! convert from m^2 to km^2
        write(6,'(a32,e24.16)')       'Total ice volume (km^3)         ',   &
                                       globalIceVolumeSum*1.0d-9     ! convert from m^3 to km^3
        write(6,'(a32,e24.16)')       'Total ice energy (J)            ',   &
                                       globalIceEnergySum
        write(6,'(a32,f24.16,i8)')    'Max thickness (m), cell         ',   &
                                       globalThicknessMax, thicknessMaxlocCell
        write(6,'(a32,f24.16,i8)')    'Min thickness (m), cell         ',   &
                                       globalThicknessMin, thicknessMinlocCell
        write(6,'(a32,f24.16)')       'Mean thickness (m)              ',   &
                                       globalthicknessMean
        write(6,'(a32,f24.16,i8,i4)') 'Max temperature (C), cell, level',   &
                                       globalTemperatureMax, temperatureMaxlocCell, temperatureMaxlocLevel
        write(6,'(a32,f24.16,i8,i4)') 'Min temperature (C), cell, level',   &
                                       globalTemperatureMin, temperatureMinlocCell, temperatureMinlocLevel
        write(6,'(a32,f24.16)')       'Mean temperature (C)            ',   &
                                       globalTemperatureMean
        write(6,'(a32,f24.16,i8,i4)') 'Max velocity (m/yr), edge, level',   &
                                       globalVelocityMax * scyr, velocityMaxlocEdge, velocityMaxlocLevel
        write(6,'(a32,f24.16,i8,i4)') 'Max basal velo (m/yr), edge     ',   &
                                       globalBasalVelocityMax * scyr, basalVelocityMaxlocEdge
        write(6,*) ' '
        write(6,'(a30,i6)')  'Column diagnostics: cell ID = ', statsCellID
        write(6,'(a30,3i6)') 'Local cell ID, block, proc =  ', diagnosticCell, diagnosticBlockID, diagnosticProcID
        write(6,*) ' '
        write(6,'(a25,f24.16)') 'Upper surface (m)        ', diagnosticUpperSurface
        write(6,'(a25,f24.16)') 'Thickness (m)            ', diagnosticThickness
        write(6,'(a25,f24.16)') 'Bed topography (m)       ', diagnosticBedTopography
        write(6,'(a25,f24.16)') 'Sfc mass balance (m/yr)  ', diagnosticSfcMassBal
        write(6,*) ' '
        write(6,'(a55)') 'Sigma        Ice speed (m/yr)       Ice temperature (C)'
        write(6,'(f6.4, a25, f24.16)') 0.0_RKIND, '------', diagnosticSurfaceTemperature
        do kLevel = 1, nVertLevels
           write(6,'(f6.4, f25.16, f24.16)') &
                layerCenterSigma(kLevel), diagnosticSpeed(kLevel), diagnosticTemperature(kLevel)
        end do
        write(6,'(f6.4, a25, f24.16)') 1.0_RKIND, '------', diagnosticBasalTemperature
        write(6,*) ' '
     endif   ! my_proc_id = IO_NODE

     ! clean up
     deallocate(workArray)
     deallocate(workArray1d)
     deallocate(diagnosticSpeed)
     deallocate(diagnosticTemperature)
     
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
        localMin = 1.0e34
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
        localMax = -1.0e34
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
        localVertSumMin = 1.0e34
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
        localVertSumMax = -1.0e34
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
!!         write(fileID,'(1i0, 100es24.16)') timeIndex, timeIndex*dt, globalFluidThickness, globalPotentialVorticity, globalPotentialEnstrophy, &
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
