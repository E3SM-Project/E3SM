! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_global_stats
!
!> \brief MPAS ocean analysis core member: global statistics
!> \author Mark Petersen
!> \date   November 2013
!> \details
!>  MPAS ocean analysis core member: global statistics
!
!-----------------------------------------------------------------------

module ocn_global_stats

   use mpas_grid_types
   use mpas_timer
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
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

   public :: ocn_setup_packages_global_stats, &
             ocn_init_global_stats, &
             ocn_compute_global_stats, &
             ocn_restart_global_stats, &
             ocn_finalize_global_stats

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
!  routine ocn_setup_packages_global_stats
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_global_stats(configPool, packagePool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: configPool
      type (mpas_pool_type), intent(in) :: packagePool

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

      logical, pointer :: amGlobalStatsActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amGlobalStatsActive', amGlobalStatsActive)

      ! turn on package for this analysis member
      amGlobalStatsActive = .true.

   end subroutine ocn_setup_packages_global_stats!}}}


!***********************************************************************
!
!  routine ocn_init_global_stats
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_global_stats(domain, err)!{{{

      use mpas_packages

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

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_init_global_stats!}}}

!***********************************************************************
!
!  routine ocn_compute_global_stats
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_global_stats(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      type (mpas_pool_type), pointer :: amGlobalStatsPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      integer :: err_tmp
      integer :: nCellsGlobal, nEdgesGlobal, nVerticesGlobal, iTracer
      integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: k, i, fileID
      integer :: timeYYYY, timeMM, timeDD, timeH, timeM, timeS
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_tracers
      character*1 timeChar
      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      character (len=StrKIND), pointer :: xtime

      real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal, localCFL, localSum, areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal, time_days
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, tangentialVelocity, layerThicknessEdge, relativeVorticity, kineticEnergyCell, &
         normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, pressure, montgomeryPotential, vertTransportVelocityTop, vertVelocityTop, &
         lowFreqDivergence, highFreqThickness, density
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      
      real (kind=RKIND), dimension(:), pointer :: minGlobalStats,maxGlobalStats,sumGlobalStats, averages, rms, verticalSumMins, verticalSumMaxes
      real (kind=RKIND), dimension(kMaxVariables) :: sumSquares, reductions, sums, mins, maxes
      real (kind=RKIND), dimension(kMaxVariables) :: sums_tmp, sumSquares_tmp, mins_tmp, maxes_tmp, averages_tmp, verticalSumMins_tmp, verticalSumMaxes_tmp

      real (kind=RKIND), dimension(:,:), allocatable :: enstrophy, normalizedAbsoluteVorticity, workArray

      logical, pointer :: thicknessFilterActive

      err = 0

      dminfo = domain % dminfo

      call mpas_timer_start("compute_global_stats", .false., amGlobalStatsTimer)

      call mpas_pool_get_package(ocnPackages, 'thicknessFilterActive', thicknessFilterActive)

      ! write out data to Analysis Member output
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amGlobalStats', amGlobalStatsPool)
      call mpas_pool_get_array(amGlobalStatsPool, 'minGlobalStats', minGlobalStats)
      call mpas_pool_get_array(amGlobalStatsPool, 'maxGlobalStats', maxGlobalStats)
      call mpas_pool_get_array(amGlobalStatsPool, 'sumGlobalStats', sumGlobalStats)
      call mpas_pool_get_array(amGlobalStatsPool, 'rmsGlobalStats', rms)
      call mpas_pool_get_array(amGlobalStatsPool, 'avgGlobalStats', averages)
      call mpas_pool_get_array(amGlobalStatsPool, 'vertSumMinGlobalStats', verticalSumMins)
      call mpas_pool_get_array(amGlobalStatsPool, 'vertSumMaxGlobalStats', verticalSumMaxes)

      sums = 0.0
      sumSquares = 0.0
      mins = 1.0e34
      maxes = -1.0e34
      averages = 0.0
      verticalSumMins = 1.0e34
      verticalSumMaxes = -1.0e34
      reductions = 0.0

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'amGlobalStats', amGlobalStatsPool)
         
         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)        
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)          
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)          
         call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)    
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)    
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop) 
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness)  
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity)  
         call mpas_pool_get_array(statePool, 'tracers', tracers)         
         if(thicknessFilterActive) then
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness)
         end if

         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', montgomeryPotential)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
         call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdge)
         call mpas_pool_get_array(diagnosticsPool, 'normalizedPlanetaryVorticityEdge', normalizedPlanetaryVorticityEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)

         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

         allocate(areaEdge(1:nEdgesSolve))
         areaEdge = dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve)


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
         if (thicknessFilterActive) then
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
               lowFreqDivergence(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
               verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         end if
   
         ! highFreqThickness
         variableIndex = variableIndex + 1
         if (thicknessFilterActive) then
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), layerThickness(:,1:nCellsSolve), &
               highFreqThickness(:,1:nCellsSolve), sums_tmp(variableIndex), sumSquares_tmp(variableIndex), mins_tmp(variableIndex), maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
               verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         end if

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

         deallocate(areaEdge)

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
      read (xtime, '(i4,5(a1,i2))'), timeYYYY, timeChar, timeMM, timeChar, timeDD, timeChar, timeH, timeChar, timeM, timeChar, timeS
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

      if (thicknessFilterActive) then
         ! lowFreqDivergence
         variableIndex = variableIndex + 1
         averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)
   
         ! highFreqThickness
         variableIndex = variableIndex + 1
         averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)
      else
         ! lowFreqDivergence
         variableIndex = variableIndex + 1
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
   
         ! highFreqThickness
         variableIndex = variableIndex + 1
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
      end if

      ! Tracers
      do iTracer=1,num_tracers
        variableIndex = variableIndex + 1
        averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
        rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)
      enddo

      minGlobalStats(1:nVariables) =  mins(1:nVariables)
      maxGlobalStats(1:nVariables) =  maxes(1:nVariables)
      sumGlobalStats(1:nVariables) =  sums(1:nVariables)

      call mpas_timer_stop("global_stats", amGlobalStatsTimer)

   end subroutine ocn_compute_global_stats!}}}

!***********************************************************************
!
!  routine ocn_restart_global_stats
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_global_stats(domain, err)!{{{

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

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_restart_global_stats!}}}

!***********************************************************************
!
!  routine ocn_finalize_global_stats
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    November 2013
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_global_stats(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

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

      err = 0

   end subroutine ocn_finalize_global_stats!}}}

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

end module ocn_global_stats

! vim: foldmethod=marker
