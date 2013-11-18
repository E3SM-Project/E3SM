! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  oac_am_global_stats
!
!> \brief Driver for MPAS ocean analysis core
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  This is the driver for the 
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module oac_am_global_stats

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar

   use ocn_diagnostics_routines

   implicit none
   private
   save

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

   public :: oac_init_am_global_stats, &
             oac_compute_am_global_stats, &
             oac_write_am_global_stats, &
             oac_restart_am_global_stats, &
             oac_finalize_am_global_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amGlobalStatsTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine oac_init_am_global_stats
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_init_am_global_stats(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp

      err = 0

print *, 'oac_init_am_global_stats'

   end subroutine oac_init_am_global_stats!}}}

!***********************************************************************
!
!  routine oac_compute_am_global_stats
!
!> \brief   Computeialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_compute_am_global_stats(amGlobalStats, state_ptr, dminfo, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(in), pointer :: state_ptr
      type (dm_info), intent(in), pointer :: dminfo
      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      type (amGlobalStats_type), intent(inout) :: amGlobalStats
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp
      type (block_type), pointer :: block

      type (state_type), pointer :: state
      type (mesh_type), pointer :: mesh
      type (scratch_type), pointer :: scratch
      ! mrp remove diagnostics links later
      type (diagnostics_type), pointer :: diagnostics

      integer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, nCellsGlobal, nEdgesGlobal, nVerticesGlobal, iTracer
      integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: k,i, num_tracers, fileID
      integer :: timeYYYY, timeMM, timeDD, timeH, timeM, timeS
      character*1 timeChar
      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal, localCFL, localSum, areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal, time_days
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, tangentialVelocity, layerThicknessEdge, kineticEnergyCell, &
         normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, pressure, montgomeryPotential, vertTransportVelocityTop, vertVelocityTop, &
         lowFreqDivergence, highFreqThickness, density
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      
      real (kind=RKIND), dimension(kMaxVariables) :: sums, sumSquares, mins, maxes, averages, rms, verticalSumMins, verticalSumMaxes, reductions
      real (kind=RKIND), dimension(kMaxVariables) :: sums_tmp, sumSquares_tmp, mins_tmp, maxes_tmp, averages_tmp, verticalSumMins_tmp, verticalSumMaxes_tmp

      real (kind=RKIND), dimension(:,:), allocatable :: enstrophy, normalizedAbsoluteVorticity, workArray


      err = 0

      call mpas_timer_start("compute_am_global_stats", .false., amGlobalStatsTimer)

      sums = 0.0
      sumSquares = 0.0
      mins = 1.0e34
      maxes = -1.0e34
      averages = 0.0
      verticalSumMins = 1.0e34
      verticalSumMaxes = -1.0e34
      reductions = 0.0

      block => state_ptr % block
      do while (associated(block))
         state => block % state % time_levs(timeLevel) % state
         mesh => block % mesh
         scratch => block % scratch
         ! mrp remove diagnostics links later
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

         call mpas_allocate_scratch_field(scratch % workArray1Cell3D, .true.)
         call mpas_allocate_scratch_field(scratch % workArray2Cell3D, .true.)
         call mpas_allocate_scratch_field(scratch % workArray1Vertex3D, .true.)
         call mpas_allocate_scratch_field(scratch % workArray2Vertex3D, .true.)
         call mpas_allocate_scratch_field(scratch % workArray1Edge3D, .true.)
         call mpas_allocate_scratch_field(scratch % workArray2Edge3D, .true.)

         density                          => diagnostics % density % array
         montgomeryPotential              => diagnostics % montgomeryPotential % array
         pressure                         => diagnostics % pressure % array
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
         call ocn_relativeVorticity_circulation(scratch % workArray1Vertex3D % array, scratch % workArray2Vertex3D % array, &
            mesh, state % normalVelocity % array, err_tmp)
            err = ior(err, err_tmp)
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nVerticesSolve, maxLevelVertexBot(1:nVerticesSolve), areaTriangle(1:nVerticesSolve), &
            scratch % workArray1Vertex3D % array(:,1:nVerticesSolve), &
            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! enstrophy 
         ! Note: enstrophy uses the value in workArray1Vertex3D for relative vorticity, computed directly above.
         allocate(enstrophy(nVertLevels,nVerticesSolve))
         enstrophy(:,:)=scratch % workArray1Vertex3D % array(:,1:nVerticesSolve)**2
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
            ! note: dt not connected.  dt will not be available in post processing mode.
            !localCFL = max(localCFL, maxval(dt*normalVelocity(:,elementIndex)/dcEdge(elementIndex)))
         end do
         nMaxes = nMaxes + 1
         maxes(nMaxes) = localCFL

         do i = 1, nVariables
            mins(nMins+i) = min(mins(nMins+i),verticalSumMins_tmp(i))
            maxes(nMaxes+i) = max(maxes(nMaxes+i),verticalSumMaxes_tmp(i))
         end do

         nMins = nMins + nVariables
         nMaxes = nMaxes + nVariables

         call mpas_deallocate_scratch_field(scratch % workArray1Cell3D, .true.)
         call mpas_deallocate_scratch_field(scratch % workArray2Cell3D, .true.)
         call mpas_deallocate_scratch_field(scratch % workArray1Vertex3D, .true.)
         call mpas_deallocate_scratch_field(scratch % workArray2Vertex3D, .true.)
         call mpas_deallocate_scratch_field(scratch % workArray1Edge3D, .true.)
         call mpas_deallocate_scratch_field(scratch % workArray2Edge3D, .true.)

         block => block % next
      end do

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

      ! write out the data to files
!       if (dminfo % my_proc_id == IO_NODE) then
!          fileID = getFreeUnit()
!          open(fileID,file='stats_min.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') time_days, mins(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_max.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') time_days, maxes(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_sum.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') time_days, sums(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_rms.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') time_days, rms(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_avg.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') time_days, averages(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_time.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(i10,10x,a,100es24.14)') timeIndex, &
!                trim(diagnostics % xtime % scalar), dt, &
!                CFLNumberGlobal
!          close (fileID)
!          open(fileID,file='stats_colmin.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') verticalSumMins(1:nVariables)
!          close (fileID)
!          open(fileID,file='stats_colmax.txt',STATUS='UNKNOWN', POSITION='append')
!             write (fileID,'(100es24.14)') verticalSumMaxes(1:nVariables)
!          close (fileID)
!       end if

      amGlobalStats % minGlobalStats % array(1:nVariables) =  mins(1:nVariables)
      amGlobalStats % maxGlobalStats % array(1:nVariables) =  maxes(1:nVariables)
      amGlobalStats % sumGlobalStats % array(1:nVariables) =  sums(1:nVariables)
      amGlobalStats % rmsGlobalStats % array(1:nVariables) =  rms(1:nVariables)
      amGlobalStats % avgGlobalStats % array(1:nVariables) =  averages(1:nVariables)
      amGlobalStats % vertSumMinGlobalStats % array(1:nVariables) =  verticalSumMins(1:nVariables)
      amGlobalStats % vertSumMaxGlobalStats % array(1:nVariables) =  verticalSumMaxes(1:nVariables)


      diagnostics % areaCellGlobal % scalar = areaCellGlobal
      diagnostics % areaEdgeGlobal % scalar = areaEdgeGlobal
      diagnostics % areaTriangleGlobal % scalar = areaTriangleGlobal

      diagnostics % volumeCellGlobal % scalar = volumeCellGlobal
      diagnostics % volumeEdgeGlobal % scalar = volumeEdgeGlobal
      diagnostics % CFLNumberGlobal % scalar = CFLNumberGlobal
      deallocate(areaEdge)

      call mpas_timer_stop("am_global_stats", amGlobalStatsTimer)

   end subroutine oac_compute_am_global_stats!}}}

!***********************************************************************
!
!  routine oac_restart_am_global_stats
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_restart_am_global_stats(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp

      err = 0

   end subroutine oac_restart_am_global_stats!}}}

!***********************************************************************
!
!  routine oac_write_am_global_stats
!
!> \brief   Driver for MPAS-Ocean analysis output
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine writes all output for this MPAS-Ocean analysis member.
!>  At this time this is just a stub, and all analysis output is written
!>  to the output file specified by config_output_name.
!
!-----------------------------------------------------------------------

   subroutine oac_write_am_global_stats(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp

      err = 0

   end subroutine oac_write_am_global_stats!}}}

!***********************************************************************
!
!  routine oac_finalize_am_global_stats
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine oac_finalize_am_global_stats(err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err_tmp

      err = 0

print *, 'oac_finalize_am_global_stats'

   end subroutine oac_finalize_am_global_stats!}}}

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


end module oac_am_global_stats

! vim: foldmethod=marker
