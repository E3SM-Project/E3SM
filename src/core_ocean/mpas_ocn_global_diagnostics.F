! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_global_diagnostics
!
!> \brief MPAS ocean statistics for the full domain
!> \author Mark Petersen and Xylar Asay-Davis
!> \date   15 April 2013
!> \details
!>  This module contains routines to compute global statistics such as
!>  minimum, maximum, and mean of variables, taken over the full domain.
!
!-----------------------------------------------------------------------
module ocn_global_diagnostics

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_timer

   implicit none
   save
   public

   type (timer_node), pointer :: diagBlockTimer, diagMPITimer

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_compute_global_diagnostics
!
!> \brief MPAS ocean statistics for the full domain
!> \author Mark Petersen and Xylar Asay-Davis
!> \date   15 April 2013
!> \details
!>  This routines to computes and writes global statistics such as
!>  minimum, maximum, and mean of variables, taken over the full domain.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_global_diagnostics(domain, timeLevel, timeIndex, dt)!{{{

      ! Note: this routine assumes that there is only one block per processor. No looping
      ! is preformed over blocks.
      ! dminfo is the domain info needed for global communication
      ! state contains the state variables needed to compute global diagnostics
      ! mesh conains the meta data about the mesh
      ! timeIndex is the current time step counter
      ! dt is the duration of each time step
      !
      ! Sums of variables at vertices are not weighted by thickness (since layerThickness is not known at
      !    vertices as it is at cell centers and at edges).
      !
      ! RMS here is volume-weighted root mean square, i.e.
      ! rms = sqrt( sum( T_i^2*v_i) / sum(v_i) )
      ! where T is the field of interest
      ! and v is the volume of the cell.

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information
      integer, intent(in) :: timeIndex
      real (kind=RKIND), intent(in) :: dt

      type (block_type), pointer :: block
      type (dm_info), pointer :: dminfo
      type (state_type), pointer :: state
      type (mesh_type), pointer :: mesh
      type (diagnostics_type), pointer :: diagnostics

      integer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nCellsGlobal, nEdgesGlobal, nVerticesGlobal, iTracer
      integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: timeLevel,k,i, num_tracers, fileID
      integer :: timeYYYY, timeMM, timeDD, timeH, timeM, timeS
      character*1 timeChar
      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal, localCFL, localSum, areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal, time_days
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, tangentialVelocity, layerThicknessEdge, relativeVorticity, kineticEnergyCell, &
         normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, pressure, montgomeryPotential, vertTransportVelocityTop, vertVelocityTop, &
         lowFreqDivergence, highFreqThickness, density
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      
      real (kind=RKIND), dimension(kMaxVariables) :: sums, sumSquares, mins, maxes, averages, rms, verticalSumMins, verticalSumMaxes, reductions
      real (kind=RKIND), dimension(kMaxVariables) :: sums_tmp, sumSquares_tmp, mins_tmp, maxes_tmp, averages_tmp, verticalSumMins_tmp, verticalSumMaxes_tmp

      real (kind=RKIND), dimension(:,:), allocatable :: enstrophy, normalizedAbsoluteVorticity, workArray

      block => domain % blocklist
      dminfo => domain % dminfo

      sums = 0.0
      sumSquares = 0.0
      mins = 1.0e34
      maxes = -1.0e34
      averages = 0.0
      verticalSumMins = 1.0e34
      verticalSumMaxes = -1.0e34
      reductions = 0.0

      call mpas_timer_start("diagnostic block loop", .false., diagBlockTimer)
      do while (associated(block))
         state => block % state % time_levs(timeLevel) % state
         mesh => block % mesh
         diagnostics => block % diagnostics
         
         num_tracers = state % num_tracers

         nVertLevels = mesh % nVertLevels
         nCellsSolve = mesh % nCellsSolve
         nEdgesSolve = mesh % nEdgesSolve
         nVerticesSolve = mesh % nVerticesSolve

         areaCell          => mesh % areaCell % array
         dcEdge            => mesh % dcEdge % array
         dvEdge            => mesh % dvEdge % array
         areaTriangle      => mesh % areaTriangle % array
         maxLevelCell      => mesh % maxLevelCell % array
         maxLevelEdgeTop   => mesh % maxLevelEdgeTop % array
         maxLevelVertexBot => mesh % maxLevelVertexBot % array

         allocate(areaEdge(1:nEdgesSolve))
         areaEdge = dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve)

         layerThickness    => state % layerThickness % array
         normalVelocity    => state % normalVelocity % array
         tracers           => state % tracers % array
         lowFreqDivergence => state % lowFreqDivergence % array
         highFreqThickness => state % highFreqThickness % array

         density                          => diagnostics % density % array
         montgomeryPotential              => diagnostics % montgomeryPotential % array
         pressure                         => diagnostics % pressure % array
         relativeVorticity                => diagnostics % relativeVorticity % array
         normalizedRelativeVorticityEdge  => diagnostics % normalizedRelativeVorticityEdge % array
         normalizedPlanetaryVorticityEdge => diagnostics % normalizedPlanetaryVorticityEdge % array
         vertTransportVelocityTop         => diagnostics % vertTransportVelocityTop % array
         vertVelocityTop                  => diagnostics % vertVelocityTop % array
         tangentialVelocity               => diagnostics % tangentialVelocity % array
         layerThicknessEdge               => diagnostics % layerThicknessEdge % array
         kineticEnergyCell                => diagnostics % kineticEnergyCell % array

         allocate(workArray(nVertLevels,nCellsSolve))

         variableIndex = 0
         ! layerThickness
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! normalVelocity
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), layerThicknessEdge(:,1:nEdgesSolve), &
            normalVelocity(:,1:nEdgesSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! tangentialVelocity
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), layerThicknessEdge(:,1:nEdgesSolve), &
            tangentialVelocity(:,1:nEdgesSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! layerThicknessEdge
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), layerThicknessEdge(:,1:nEdgesSolve), &
            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! relativeVorticity
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nVerticesSolve, maxLevelVertexBot(1:nVerticesSolve), areaTriangle(1:nVerticesSolve), relativeVorticity(:,1:nVerticesSolve), &
            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! enstrophy
         allocate(enstrophy(nVertLevels,nVerticesSolve))
         enstrophy(:,:)=relativeVorticity(:,1:nVerticesSolve)**2
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nVerticesSolve, maxLevelVertexBot(1:nVerticesSolve), areaTriangle(1:nVerticesSolve), &
            enstrophy(:,:), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), &
            verticalSumMins_tmp(variableIndex), verticalSumMaxes_tmp(variableIndex))
         deallocate(enstrophy)
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! kineticEnergyCell
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            kineticEnergyCell(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! normalizedAbsoluteVorticity
         allocate(normalizedAbsoluteVorticity(nVertLevels,nEdgesSolve))
         normalizedAbsoluteVorticity(:,:) = normalizedRelativeVorticityEdge(:,1:nEdgesSolve) + normalizedPlanetaryVorticityEdge(:,1:nEdgesSolve)
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), layerThicknessEdge(:,1:nEdgesSolve), &
            normalizedAbsoluteVorticity(:,1:nEdgesSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         deallocate(normalizedAbsoluteVorticity)
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! pressure
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            pressure(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! montgomeryPotential
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            montgomeryPotential(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! vertVelocityTop vertical velocity
         variableIndex = variableIndex + 1
         workArray = vertVelocityTop(1:nVertLevels,1:nCellsSolve)
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            workArray, sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! vertTransportVelocityTop vertical velocity
         variableIndex = variableIndex + 1
         workArray = vertTransportVelocityTop(1:nVertLevels,1:nCellsSolve)
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            workArray, sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! lowFreqDivergence
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            lowFreqDivergence(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! highFreqThickness
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
            highFreqThickness(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! Tracers
         do iTracer=1,num_tracers
            variableIndex = variableIndex + 1
            workArray = Tracers(iTracer,:,1:nCellsSolve)
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
               workArray, sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
               verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         enddo
         deallocate(workArray)

         nVariables = variableIndex
         nSums = nVariables
         nMins = nVariables
         nMaxes = nVariables

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + sum(areaCell(1:nCellsSolve))

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + sum(dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve))

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + sum(areaTriangle(1:nVerticesSolve))

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + nCellsSolve

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + nEdgesSolve

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + nVerticesSolve

         localCFL = 0.0
         do elementIndex = 1,nEdgesSolve
            localCFL = max(localCFL, maxval(dt*normalVelocity(:,elementIndex)/dcEdge(elementIndex)))
         end do
         nMaxes = nMaxes + 1
         maxes(nMaxes) = localCFL

         do i = 1, nVariables
            mins(nMins+i) = min(mins(nMins+i),verticalSumMins_tmp(i))
            maxes(nMaxes+i) = max(maxes(nMaxes+i),verticalSumMaxes_tmp(i))
         end do

         nMins = nMins + nVariables
         nMaxes = nMaxes + nVariables

         block => block % next
      end do
      call mpas_timer_stop("diagnostic block loop", diagBlockTimer)
      call mpas_timer_start("diagnostics mpi", .false., diagMPITimer)

      ! global reduction of the 5 arrays (packed into 3 to minimize global communication)
      call mpas_dmpar_sum_real_array(dminfo, nSums, sums(1:nSums), reductions(1:nSums))
      sums(1:nVariables) = reductions(1:nVariables)
      areaCellGlobal = reductions(nVariables+1)
      areaEdgeGlobal = reductions(nVariables+2)
      areaTriangleGlobal = reductions(nVariables+3)
      nCellsGlobal = int(reductions(nVariables+4))
      nEdgesGlobal = int(reductions(nVariables+5))
      nVerticesGlobal = int(reductions(nVariables+6))
      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
      sumSquares(1:nVariables) = reductions(1:nVariables)

      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
      mins(1:nVariables) = reductions(1:nVariables)
      verticalSumMins(1:nVariables) = reductions(nMins-nVariables+1:nMins)

      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))
      maxes(1:nVariables) = reductions(1:nVariables)
      CFLNumberGlobal = reductions(nVariables+1)
      verticalSumMaxes(1:nVariables) = reductions(nMaxes-nVariables+1:nMaxes)

      volumeCellGlobal = sums(1)
      volumeEdgeGlobal = sums(4)


      ! compute the averages (slightly different depending on how the sum was computed)
      variableIndex = 0

      ! time, in days, using a 360 day calendar
      read (diagnostics % xtime % scalar, '(i4,5(a1,i2))'), timeYYYY, timeChar, timeMM, timeChar, timeDD, timeChar, timeH, timeChar, timeM, timeChar, timeS
      ! subtract 31.0 because calendar starts on 00-01-01
      time_days = timeYYYY*360.0 + timeMM*30.0 + timeDD + (timeH + (timeM + timeS/60.0)/60.0)/24.0 - 31.0

      ! layerThickness
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal*nVertLevels))

      ! normalVelocity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeEdgeGlobal)

      ! tangentialVelocity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeEdgeGlobal)

      ! layerThicknessEdge
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaEdgeGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaEdgeGlobal*nVertLevels))

      ! relativeVorticity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaTriangleGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaTriangleGlobal*nVertLevels))

      ! enstrophy
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaTriangleGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaTriangleGlobal*nVertLevels))

      ! kineticEnergyCell
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! normalizedAbsoluteVorticity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeEdgeGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeEdgeGlobal)

      ! pressure
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! montgomeryPotential
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! vertVelocityTop vertical velocity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! vertTransportVelocityTop vertical velocity
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! lowFreqDivergence
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! highFreqThickness
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)

      ! Tracers
      do iTracer=1,num_tracers
        variableIndex = variableIndex + 1
        averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
        rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)
      enddo

      call mpas_timer_stop("diagnostics mpi", diagMPITimer)

      ! write out the data to files
      if (dminfo % my_proc_id == IO_NODE) then
         fileID = getFreeUnit()
         open(fileID,file='stats_min.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') time_days, mins(1:nVariables)
         close (fileID)
         open(fileID,file='stats_max.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') time_days, maxes(1:nVariables)
         close (fileID)
         open(fileID,file='stats_sum.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') time_days, sums(1:nVariables)
         close (fileID)
         open(fileID,file='stats_rms.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') time_days, rms(1:nVariables)
         close (fileID)
         open(fileID,file='stats_avg.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') time_days, averages(1:nVariables)
         close (fileID)
         open(fileID,file='stats_time.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(i10,10x,a,100es24.14)') timeIndex, &
               trim(diagnostics % xtime % scalar), dt, &
               CFLNumberGlobal
         close (fileID)
         open(fileID,file='stats_colmin.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') verticalSumMins(1:nVariables)
         close (fileID)
         open(fileID,file='stats_colmax.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') verticalSumMaxes(1:nVariables)
         close (fileID)
      end if

      diagnostics % areaCellGlobal % scalar = areaCellGlobal
      diagnostics % areaEdgeGlobal % scalar = areaEdgeGlobal
      diagnostics % areaTriangleGlobal % scalar = areaTriangleGlobal

      diagnostics % volumeCellGlobal % scalar = volumeCellGlobal
      diagnostics % volumeEdgeGlobal % scalar = volumeEdgeGlobal
      diagnostics % CFLNumberGlobal % scalar = CFLNumberGlobal
      deallocate(areaEdge)

   end subroutine ocn_compute_global_diagnostics!}}}

   integer function getFreeUnit()!{{{
      implicit none

      integer :: index
      logical :: isOpened

      getFreeUnit = 0
      do index = 1,99
         if((index /= 5) .and. (index /= 6)) then
            inquire(unit = index, opened = isOpened)
            if( .not. isOpened) then
               getFreeUnit = index
               return
            end if
         end if
      end do
   end function getFreeUnit!}}}

   subroutine ocn_compute_field_local_stats(dminfo, nVertLevels, nElements, field, localSum, localMin, localMax, localVertSumMin, &!{{{
      localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      localSum = sum(field)
      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(field,1))
      localVertSumMax = maxval(sum(field,1))

   end subroutine ocn_compute_field_local_stats!}}}

   subroutine ocn_compute_field_area_weighted_local_stats(dminfo, nVertLevels, nElements, areas, field, localSum, localMin, &!{{{
      localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      integer :: elementIndex

      localSum = 0.0
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(field(:,elementIndex))
      end do

      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(field,1))
      localVertSumMax = maxval(sum(field,1))

   end subroutine ocn_compute_field_area_weighted_local_stats!}}}

   subroutine ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nElements, maxLevel, areas, field, &!{{{
      localSum, localRMS, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      integer, dimension(nElements), intent(in) :: maxLevel
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localRMS, localMin, localMax, localVertSumMin, &
      localVertSumMax

      integer :: elementIndex
      real (kind=RKIND) :: colSum, colRMS, colSumAbs

      localSum = 0.0
      localRMS = 0.0
      localMin =  1.0e34
      localMax = -1.0e34
      localVertSumMin =  1.0e34
      localVertSumMax = -1.0e34

      do elementIndex = 1, nElements
        colSum = sum(field(1:maxLevel(elementIndex),elementIndex))
        localSum = localSum + areas(elementIndex) * colSum
        colRMS = sum(field(1:maxLevel(elementIndex),elementIndex)**2)
        localRMS = localRMS + areas(elementIndex) * colRMS
        localMin = min(localMin,minval(field(1:maxLevel(elementIndex),elementIndex)))
        localMax = max(localMax,maxval(field(1:maxLevel(elementIndex),elementIndex)))
        localVertSumMin = min(localVertSumMin,colSum)
        localVertSumMax = max(localVertSumMax,colSum)
      end do

   end subroutine ocn_compute_field_area_weighted_local_stats_max_level!}}}

   subroutine ocn_compute_field_thickness_weighted_local_stats(dminfo, nVertLevels, nElements, h, field, &!{{{
      localSum, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      integer :: elementIndex

      localSum = sum(h*field)
      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(h*field,1))
      localVertSumMax = maxval(sum(h*field,1))

   end subroutine ocn_compute_field_thickness_weighted_local_stats!}}}

   subroutine ocn_compute_field_volume_weighted_local_stats(dminfo, nVertLevels, nElements, areas, layerThickness, field, &!{{{
      localSum, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: layerThickness
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localMin, localMax, localVertSumMin, &
      localVertSumMax

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      integer :: elementIndex

      localSum = 0.0
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(layerThickness(:,elementIndex)*field(:,elementIndex))
      end do

      localMin = minval(field)
      localMax = maxval(field)
      localVertSumMin = minval(sum(layerThickness*field,1))
      localVertSumMax = maxval(sum(layerThickness*field,1))

   end subroutine ocn_compute_field_volume_weighted_local_stats!}}}

   subroutine ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nElements, maxLevel, areas, layerThickness, field, &!{{{
      localSum, localRMS, localMin, localMax, localVertSumMin, localVertSumMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      integer, dimension(nElements), intent(in) :: maxLevel
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: layerThickness
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: localSum, localRMS, localMin, localMax, localVertSumMin, &
         localVertSumMax

      integer :: elementIndex
      real (kind=RKIND) :: thicknessWeightedColSum, thicknessWeightedColRMS, thicknessWeightedColSumAbs
      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      localSum = 0.0
      localRMS = 0.0
      localMin =  1.0e34
      localMax = -1.0e34
      localVertSumMin =  1.0e34
      localVertSumMax = -1.0e34

      do elementIndex = 1, nElements
        thicknessWeightedColSum = sum(layerThickness(1:maxLevel(elementIndex),elementIndex)*field(1:maxLevel(elementIndex),elementIndex))
        localSum = localSum + areas(elementIndex) * thicknessWeightedColSum
        thicknessWeightedColRMS = sum(layerThickness(1:maxLevel(elementIndex),elementIndex)*field(1:maxLevel(elementIndex),elementIndex)**2)
        localRMS = localRMS + areas(elementIndex) * thicknessWeightedColRMS
        localMin = min(localMin,minval(field(1:maxLevel(elementIndex),elementIndex)))
        localMax = max(localMax,maxval(field(1:maxLevel(elementIndex),elementIndex)))
        localVertSumMin = min(localVertSumMin,thicknessWeightedColSum)
        localVertSumMax = max(localVertSumMax,thicknessWeightedColSum)
      end do

   end subroutine ocn_compute_field_volume_weighted_local_stats_max_level!}}}

   subroutine ocn_compute_global_sum(dminfo, nVertLevels, nElements, field, globalSum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum

      real (kind=RKIND) :: localSum

      localSum = sum(field)
      call mpas_dmpar_sum_real(dminfo, localSum, globalSum)

   end subroutine ocn_compute_global_sum!}}}

   subroutine ocn_compute_area_weighted_global_sum(dminfo, nVertLevels, nElements, areas, field, globalSum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum
      
      integer :: elementIndex
      real (kind=RKIND) :: localSum

      localSum = 0.
      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * sum(field(:,elementIndex))
      end do
   
      call mpas_dmpar_sum_real(dminfo, localSum, globalSum)
       
   end subroutine ocn_compute_area_weighted_global_sum!}}}

   subroutine ocn_compute_volume_weighted_global_sum(dminfo, nVertLevels, nElements, areas, h, field, globalSum)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalSum

      real (kind=RKIND), dimension(nVertLevels, nElements) :: hTimesField

      hTimesField = h*field

      call ocn_compute_area_weighted_global_sum(dminfo, nVertLevels, nElements, areas, hTimesField, globalSum)

   end subroutine ocn_compute_volume_weighted_global_sum!}}}

   subroutine ocn_compute_global_min(dminfo, nVertLevels, nElements, field, globalMin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(field)
      call mpas_dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine ocn_compute_global_min!}}}

   subroutine ocn_compute_global_max(dminfo, nVertLevels, nElements, field, globalMax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(field)
      call mpas_dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine ocn_compute_global_max!}}}

   subroutine ocn_compute_global_vert_sum_horiz_min(dminfo, nVertLevels, nElements, field, globalMin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(sum(field,1))
      call mpas_dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine ocn_compute_global_vert_sum_horiz_min!}}}

   subroutine ocn_compute_global_vert_sum_horiz_max(dminfo, nVertLevels, nElements, field, globalMax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(sum(field,1))
      call mpas_dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine ocn_compute_global_vert_sum_horiz_max!}}}

   subroutine ocn_compute_global_vert_thickness_weighted_sum_horiz_min(dminfo, nVertLevels, nElements, h, field, globalMin)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h, field
      real (kind=RKIND), intent(out) :: globalMin

      real (kind=RKIND) :: localMin

      localMin = minval(sum(h*field,1))
      call mpas_dmpar_min_real(dminfo, localMin, globalMin)

   end subroutine ocn_compute_global_vert_thickness_weighted_sum_horiz_min!}}}

   subroutine ocn_compute_global_vert_thickness_weighted_sum_horiz_max(dminfo, nVertLevels, nElements, h, field, globalMax)!{{{

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nVertLevels, nElements
      real (kind=RKIND), dimension(nVertLevels, nElements), intent(in) :: h, field
      real (kind=RKIND), intent(out) :: globalMax

      real (kind=RKIND) :: localMax

      localMax = maxval(sum(h*field,1))
      call mpas_dmpar_max_real(dminfo, localMax, globalMax)

   end subroutine ocn_compute_global_vert_thickness_weighted_sum_horiz_max!}}}

   subroutine ocn_global_diagnostics_init(dminfo,err)!{{{

      ! Create stats_readme.txt file listing variable names

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(out) :: err
      integer :: fileID, i

      err = 0

      if (dminfo % my_proc_id == IO_NODE) then
         fileID = getFreeUnit()
         open(fileID,file='stats_readme.txt',STATUS='UNKNOWN', POSITION='rewind')

            write (fileID,'(a)') 'readme file for MPAS-Ocean global statistics'
            write (fileID,'(/,a)') 'stats_time.txt. contains: timeIndex, timestamp, dt, CFLNumberGlobal'
            write (fileID,'(/,a)') 'All other stats_*.txt. contain the following columns.  Rows correspond to timestamps in rows of stats_time.txt'
            write (fileID,'(a)')   "See user's guide for units associated with these variables."

            i=1
            write (fileID,'(i5,a)') i,'. time, in days, using a 360 day calendar'; i=i+1
            write (fileID,'(i5,a)') i,'. layerThickness'; i=i+1
            write (fileID,'(i5,a)') i,'. normalVelocity'; i=i+1
            write (fileID,'(i5,a)') i,'. tangentialVelocity'; i=i+1
            write (fileID,'(i5,a)') i,'. layerThicknessEdge'; i=i+1
            write (fileID,'(i5,a)') i,'. relativeVorticity'; i=i+1
            write (fileID,'(i5,a)') i,'. enstrophy = relativeVorticity**2'; i=i+1
            write (fileID,'(i5,a)') i,'. kineticEnergyCell'; i=i+1
            write (fileID,'(i5,a)') i,'. normalizedAbsoluteVorticity = (relative vorticity + planetary vorticity)/layer thickness'; i=i+1
            write (fileID,'(i5,a)') i,'. pressure'; i=i+1
            write (fileID,'(i5,a)') i,'. montgomeryPotential'; i=i+1
            write (fileID,'(i5,a)') i,'. vertVelocityTop vertical velocity'; i=i+1
            write (fileID,'(i5,a)') i,'. vertTransportVelocityTop vertical transport'; i=i+1
            write (fileID,'(i5,a)') i,'. lowFreqDivergence'; i=i+1
            write (fileID,'(i5,a)') i,'. highFreqThickness'; i=i+1
            write (fileID,'(i5,a)') i,'. Tracers: usually T, S, then others in remaining columns'

            write (fileID,'(/,a)') 'A chain of simple unix commands may be used to access a specific part of the data. For example,'
            write (fileID,'(a)') 'to view the last three values of column seven in the global average, use:'
            write (fileID,'(a)') "cat stats_avg.txt | awk '{print $7}' | tail -n3"
  
         close (fileID)
      endif

   end subroutine ocn_global_diagnostics_init!}}}


end module ocn_global_diagnostics
