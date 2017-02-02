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

   use mpas_derived_types
   use mpas_pool_routines
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

   public :: ocn_init_global_stats, &
             ocn_compute_global_stats, &
             ocn_restart_global_stats, &
             ocn_finalize_global_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

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

      type (dm_info), pointer :: dminfo
      logical, pointer :: config_AM_globalStats_text_file
      character (len=StrKIND), pointer :: config_AM_globalStats_directory
      integer :: fileID, i

      err = 0

      call mpas_pool_get_config(domain % configs, 'config_AM_globalStats_directory', config_AM_globalStats_directory)
      call mpas_pool_get_config(domain % configs, 'config_AM_globalStats_text_file', config_AM_globalStats_text_file)
      if (config_AM_globalStats_text_file) then
      dminfo => domain % dminfo
      if (dminfo % my_proc_id == IO_NODE) then
         fileID = getFreeUnit()
         open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_readme.txt',STATUS='UNKNOWN', POSITION='rewind')

            write (fileID,'(a)') 'readme file for MPAS-Ocean global statistics'
            write (fileID,'(/,a)') 'stats_time.txt. contains: timestamp, time (in days), CFLNumberGlobal'
            write (fileID,'(/,a)') 'All other stats_*.txt. contain the following columns.  Rows correspond to timestamps ' &
                                // 'in rows of stats_time.txt'
            write (fileID,'(a)')   "See user's guide for units associated with these variables."

            i=1
            write (fileID,'(i5,a)') i,'. time, in days'; i=i+1
            write (fileID,'(i5,a)') i,'. layerThickness'; i=i+1
            write (fileID,'(i5,a)') i,'. normalVelocity'; i=i+1
            write (fileID,'(i5,a)') i,'. tangentialVelocity'; i=i+1
            write (fileID,'(i5,a)') i,'. layerThicknessEdge'; i=i+1
            write (fileID,'(i5,a)') i,'. relativeVorticity'; i=i+1
            write (fileID,'(i5,a)') i,'. enstrophy = relativeVorticity**2'; i=i+1
            write (fileID,'(i5,a)') i,'. kineticEnergyCell'; i=i+1
            write (fileID,'(i5,a)') i,'. normalizedAbsoluteVorticity = (relative vorticity + planetary vorticity)/layer ' &
                                   // 'thickness'; i=i+1
            write (fileID,'(i5,a)') i,'. pressure'; i=i+1
            write (fileID,'(i5,a)') i,'. montgomeryPotential'; i=i+1
            write (fileID,'(i5,a)') i,'. vertVelocityTop vertical velocity'; i=i+1
            write (fileID,'(i5,a)') i,'. vertAleTransportTop vertical transport'; i=i+1
            write (fileID,'(i5,a)') i,'. lowFreqDivergence'; i=i+1
            write (fileID,'(i5,a)') i,'. highFreqThickness'; i=i+1
            write (fileID,'(i5,a)') i,'. Tracers: usually T, S, then others in remaining columns'

            write (fileID,'(/,a)') 'A chain of simple unix commands may be used to access a specific part of the data. ' &
                                // 'For example,'
            write (fileID,'(a)') 'to view the last three values of column seven in the global average, use:'
            write (fileID,'(a)') "cat stats_avg.txt | awk '{print $7}' | tail -n3"

         close (fileID)
      endif

      endif

   end subroutine ocn_init_global_stats!}}}

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

      type (mpas_pool_type), pointer :: globalStatsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool

      integer :: err_tmp
      integer :: nCellsGlobal, nEdgesGlobal, nVerticesGlobal, iTracer
      integer :: iEdge, variableIndex, nVariables, nSums, nMaxes, nMins
      integer :: k, i, fileID
      integer :: timeYYYY, timeMM, timeDD, timeH, timeM, timeS
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_activeTracers

      integer, parameter :: kMaxVariables = 1024 ! this must be a little more than double the number of variables to be reduced
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot, &
         landiceMask

      character (len=StrKIND), pointer :: xtime, simulationStartTime
      type (MPAS_Time_type) :: xtime_timeType, simulationStartTime_timeType
      type (MPAS_TimeInterval_type) :: timeStep

      real (kind=RKIND) :: localCFL, localSum, dt, currentVolume
      real (kind=RKIND), pointer :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal, areaCellGlobal, &
                                    areaEdgeGlobal, areaTriangleGlobal, totalVolumeChange, netFreshwaterInput, &
                                    absoluteFreshWaterConservation, relativeFreshWaterConservation, &
                                    landIceFloatingAreaSum
      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge, areaTriangle, areaEdge, landIceFloatingArea
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, normalVelocity, tangentialVelocity, layerThicknessEdge, &
            relativeVorticity, kineticEnergyCell, normalizedRelativeVorticityEdge, &
            normalizedPlanetaryVorticityEdge, pressure, montgomeryPotential, &
            vertAleTransportTop, vertVelocityTop, lowFreqDivergence, highFreqThickness, &
            density, layerThicknessPreviousTimestep, activeTracersSurfaceFlux, &
            activeTracersPistonVelocity, activeTracersSurfaceRestoringValue

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      real (kind=RKIND), dimension(:), pointer :: evaporationFlux, snowFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux, riverRunoffFlux, iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux, landIceFreshwaterFlux
      real (kind=RKIND), dimension(:), pointer :: accumulatedLandIceMass, accumulatedLandIceHeat, &
          accumulatedLandIceFrazilMass

      real (kind=RKIND), dimension(:,:), pointer :: frazilLayerThicknessTendency

      real (kind=RKIND), pointer :: daysSinceStartOfSim
      real (kind=RKIND), dimension(:), pointer :: minGlobalStats,maxGlobalStats,sumGlobalStats, averages, rms, verticalSumMins, &
          verticalSumMaxes
      real (kind=RKIND), dimension(kMaxVariables) :: sumSquares, reductions, sums, mins, maxes
      real (kind=RKIND), dimension(kMaxVariables) :: sums_tmp, sumSquares_tmp, mins_tmp, maxes_tmp, averages_tmp, &
             verticalSumMins_tmp, verticalSumMaxes_tmp

      real (kind=RKIND), dimension(:,:), allocatable :: enstrophy, normalizedAbsoluteVorticity, workArray

      logical, pointer :: thicknessFilterActive
      logical, pointer :: config_AM_globalStats_text_file
      character (len=StrKIND), pointer :: config_AM_globalStats_directory

      real (kind=RKIND), dimension(:), allocatable :: restoringSaltFlux

      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_package(ocnPackages, 'thicknessFilterActive', thicknessFilterActive)

      ! write out data to Analysis Member output
      call mpas_pool_get_subpool(domain % blocklist % structs, 'globalStatsAM', globalStatsAMPool)
      call mpas_pool_get_array(globalStatsAMPool, 'minGlobalStats', minGlobalStats)
      call mpas_pool_get_array(globalStatsAMPool, 'maxGlobalStats', maxGlobalStats)
      call mpas_pool_get_array(globalStatsAMPool, 'sumGlobalStats', sumGlobalStats)
      call mpas_pool_get_array(globalStatsAMPool, 'rmsGlobalStats', rms)
      call mpas_pool_get_array(globalStatsAMPool, 'avgGlobalStats', averages)
      call mpas_pool_get_array(globalStatsAMPool, 'vertSumMinGlobalStats', verticalSumMins)
      call mpas_pool_get_array(globalStatsAMPool, 'vertSumMaxGlobalStats', verticalSumMaxes)
      call mpas_pool_get_array(globalStatsAMPool, 'areaCellGlobal', areaCellGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'areaEdgeGlobal', areaEdgeGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'areaTriangleGlobal', areaTriangleGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'volumeCellGlobal', volumeCellGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'volumeEdgeGlobal', volumeEdgeGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'CFLNumberGlobal', CFLNumberGlobal)
      call mpas_pool_get_array(globalStatsAMPool, 'landIceFloatingAreaSum', landIceFloatingAreaSum)

      call mpas_pool_get_array(globalStatsAMPool, 'totalVolumeChange', totalVolumeChange)
      call mpas_pool_get_array(globalStatsAMPool, 'netFreshwaterInput', netFreshwaterInput)
      call mpas_pool_get_array(globalStatsAMPool, 'absoluteFreshWaterConservation', absoluteFreshWaterConservation)
      call mpas_pool_get_array(globalStatsAMPool, 'relativeFreshWaterConservation', relativeFreshWaterConservation)

      sums = 0.0_RKIND
      sumSquares = 0.0_RKIND
      mins = 1.0e34_RKIND
      maxes = -1.0e34_RKIND
      averages = 0.0_RKIND
      verticalSumMins = 1.0e34_RKIND
      verticalSumMaxes = -1.0e34_RKIND
      reductions = 0.0_RKIND

      timeStep = mpas_get_clock_timestep(domain % clock, ierr=err_tmp)
      call mpas_get_timeInterval(timeStep, dt=dt)

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

         call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
         call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)

         call mpas_pool_get_dimension(tracersPool, 'num_activeTracers', num_activeTracers)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessPreviousTimestep, 2)
         call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         if(thicknessFilterActive) then
            call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence, 1)
            call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, 1)
         end if

         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', montgomeryPotential)
         call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
         call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdge)
         call mpas_pool_get_array(diagnosticsPool, 'normalizedPlanetaryVorticityEdge', normalizedPlanetaryVorticityEdge)
         call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
         call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
         call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
         call mpas_pool_get_array(diagnosticsPool, 'simulationStartTime',simulationStartTime)
         call mpas_pool_get_array(diagnosticsPool, 'daysSinceStartOfSim',daysSinceStartOfSim)

         call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
         call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
         call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
         call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
         call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
         call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
         call mpas_pool_get_array(forcingPool, 'frazilLayerThicknessTendency', frazilLayerThicknessTendency)
         call mpas_pool_get_array(forcingPool, 'landIceFreshwaterFlux', landIceFreshwaterFlux)
         call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

         call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)
         call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', &
                                  activeTracersPistonVelocity)
         call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
                                  activeTracersSurfaceRestoringValue)


         call mpas_pool_get_array(statePool, 'accumulatedLandIceMass', accumulatedLandIceMass, 1)
         call mpas_pool_get_array(statePool, 'accumulatedLandIceHeat', accumulatedLandIceHeat, 1)
         call mpas_pool_get_array(statePool, 'accumulatedLandIceFrazilMass', accumulatedLandIceFrazilMass, 1)

         allocate(areaEdge(1:nEdgesSolve))
         areaEdge = dcEdge(1:nEdgesSolve)*dvEdge(1:nEdgesSolve)

         allocate(landIceFloatingArea(1:nCellsSolve))
         if ( associated(landIceMask) ) then
            landIceFloatingArea = landIceMask(1:nCellsSolve)*areaCell(1:nCellsSolve)
         else
            landIceFloatingArea = 0.
         end if

         allocate(workArray(nVertLevels,nCellsSolve))

         variableIndex = 0
         ! layerThickness
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                            maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                            layerThickness(:,1:nCellsSolve), sums_tmp(variableIndex), &
                            sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                            maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
                            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! normalVelocity
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, &
                              maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), &
                              layerThicknessEdge(:,1:nEdgesSolve), &
                              normalVelocity(:,1:nEdgesSolve), sums_tmp(variableIndex), &
                              sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                              maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
                              verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! tangentialVelocity
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, &
                              maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), &
                              layerThicknessEdge(:,1:nEdgesSolve), &
                              tangentialVelocity(:,1:nEdgesSolve), &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
                              verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! layerThicknessEdge
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, &
                            maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), &
                            layerThicknessEdge(:,1:nEdgesSolve), sums_tmp(variableIndex), &
                            sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                            maxes_tmp(variableIndex), verticalSumMins_tmp(variableIndex), &
                            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! relativeVorticity
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nVerticesSolve, &
                            maxLevelVertexBot(1:nVerticesSolve), &
                            areaTriangle(1:nVerticesSolve), &
                            relativeVorticity(:,1:nVerticesSolve), &
                            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                            mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                            verticalSumMins_tmp(variableIndex), &
                            verticalSumMaxes_tmp(variableIndex))
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
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nVerticesSolve, &
                            maxLevelVertexBot(1:nVerticesSolve), &
                            areaTriangle(1:nVerticesSolve), enstrophy(:,:), &
                            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                            mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                            verticalSumMins_tmp(variableIndex), &
                            verticalSumMaxes_tmp(variableIndex))
         deallocate(enstrophy)
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! kineticEnergyCell
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                              maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                              layerThickness(:,1:nCellsSolve), &
                              kineticEnergyCell(:,1:nCellsSolve), &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
                              verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! normalizedAbsoluteVorticity
         allocate(normalizedAbsoluteVorticity(nVertLevels,nEdgesSolve))
         normalizedAbsoluteVorticity(:,:) = normalizedRelativeVorticityEdge(:,1:nEdgesSolve) &
  + normalizedPlanetaryVorticityEdge(:,1:nEdgesSolve)
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nEdgesSolve, &
                              maxLevelEdgeTop(1:nEdgesSolve), areaEdge(1:nEdgesSolve), &
                              layerThicknessEdge(:,1:nEdgesSolve), &
                              normalizedAbsoluteVorticity(:,1:nEdgesSolve), &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
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
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                              maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                              layerThickness(:,1:nCellsSolve), pressure(:,1:nCellsSolve), &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
                              verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! montgomeryPotential
         variableIndex = variableIndex + 1
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                              maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                              layerThickness(:,1:nCellsSolve), &
                              montgomeryPotential(:,1:nCellsSolve), &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
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
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                              maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                              layerThickness(:,1:nCellsSolve), workArray, &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
                              verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! vertAleTransportTop vertical velocity
         variableIndex = variableIndex + 1
         workArray = vertAleTransportTop(1:nVertLevels,1:nCellsSolve)
         call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                              maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                              layerThickness(:,1:nCellsSolve), workArray, &
                              sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                              mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                              verticalSumMins_tmp(variableIndex), &
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
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                                 maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                                 layerThickness(:,1:nCellsSolve), &
                                 lowFreqDivergence(:,1:nCellsSolve), &
                                 sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                                 mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                                 verticalSumMins_tmp(variableIndex), &
                                 verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
            verticalSumMins_tmp(variableIndex) = 0.0_RKIND
            verticalSumMaxes_tmp(variableIndex) = 0.0_RKIND
            sums(variableIndex) = 0.0_RKIND
            sumSquares(variableIndex) = 0.0_RKIND
            mins(variableIndex) = 0.0_RKIND
            maxes(variableIndex) = 0.0_RKIND
            verticalSumMins(variableIndex) = 0.0_RKIND
            verticalSumMaxes(variableIndex) = 0.0_RKIND
         end if

         ! highFreqThickness
         variableIndex = variableIndex + 1
         if (thicknessFilterActive) then
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                                 maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                                 layerThickness(:,1:nCellsSolve), &
                                 highFreqThickness(:,1:nCellsSolve), &
                                 sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                                 mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                                 verticalSumMins_tmp(variableIndex), &
                                 verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
            verticalSumMins_tmp(variableIndex) = 0.0_RKIND
            verticalSumMaxes_tmp(variableIndex) = 0.0_RKIND
            sums(variableIndex) = 0.0_RKIND
            sumSquares(variableIndex) = 0.0_RKIND
            mins(variableIndex) = 0.0_RKIND
            maxes(variableIndex) = 0.0_RKIND
            verticalSumMins(variableIndex) = 0.0_RKIND
            verticalSumMaxes(variableIndex) = 0.0_RKIND
         end if

         ! active Tracers
         do iTracer=1,num_activeTracers
            variableIndex = variableIndex + 1
            workArray = activeTracers(iTracer,:,1:nCellsSolve)
            call ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                                 maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                                 layerThickness(:,1:nCellsSolve), workArray, &
                                 sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                                 mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                                 verticalSumMins_tmp(variableIndex), &
                                 verticalSumMaxes_tmp(variableIndex))
            sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
            sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
            mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
            maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
            verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
            verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))
         enddo

         ! layerThickness from previous timestep
         variableIndex = variableIndex + 1
         call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                            maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                            layerThicknessPreviousTimestep(:,1:nCellsSolve), &
                            sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                            mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                            verticalSumMins_tmp(variableIndex), &
                            verticalSumMaxes_tmp(variableIndex))
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! frazilLayerThicknessTendency
         variableIndex = variableIndex + 1
         if ( associated(frazilLayerThicknessTendency) ) then
            call ocn_compute_field_area_weighted_local_stats_max_level(dminfo, nVertLevels, nCellsSolve, &
                               maxLevelCell(1:nCellsSolve), areaCell(1:nCellsSolve), &
                               frazilLayerThicknessTendency(:,1:nCellsSolve), &
                               sums_tmp(variableIndex), sumSquares_tmp(variableIndex), &
                               mins_tmp(variableIndex), maxes_tmp(variableIndex), &
                               verticalSumMins_tmp(variableIndex), &
                               verticalSumMaxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
            verticalSumMins_tmp(variableIndex) = 0.0_RKIND
            verticalSumMaxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = min(verticalSumMins(variableIndex), verticalSumMins_tmp(variableIndex))
         verticalSumMaxes(variableIndex) = max(verticalSumMaxes(variableIndex), verticalSumMaxes_tmp(variableIndex))

         ! evaporationFlux
         variableIndex = variableIndex + 1
         if ( associated(evaporationFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         evaporationFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


         ! rainFlux
         variableIndex = variableIndex + 1
         if ( associated(rainFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         rainFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


         ! snowFlux
         variableIndex = variableIndex + 1
         if ( associated(snowFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         snowFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


         ! seaIceFreshWaterFlux
         variableIndex = variableIndex + 1
         if ( associated(seaIceFreshWaterFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         seaIceFreshWaterFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


         ! riverRunoffFlux
         variableIndex = variableIndex + 1
         if ( associated(riverRunoffFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         riverRunoffFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if
         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


         ! iceRunoffFlux
         variableIndex = variableIndex + 1
         if ( associated(iceRunoffFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         areaCell(1:nCellsSolve), &
                         iceRunoffFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! temperature flux
         variableIndex = variableIndex + 1
         if ( associated(activeTracersSurfaceFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                   areaCell(1:nCellsSolve), &
                   activeTracersSurfaceFlux(1,1:nCellsSolve), sums_tmp(variableIndex), &
                   sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                   maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! salinity flux
         variableIndex = variableIndex + 1
         if ( associated(activeTracersSurfaceFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                   areaCell(1:nCellsSolve), &
                   activeTracersSurfaceFlux(2,1:nCellsSolve), sums_tmp(variableIndex), &
                   sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                   maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! salinity restoring flux
         variableIndex = variableIndex + 1
         if ( associated(activeTracersPistonVelocity) .and. associated(activeTracersSurfaceRestoringValue) ) then
            allocate (restoringSaltFlux(nCellsSolve))
            restoringSaltFlux = activeTracersPistonVelocity(2,1:nCellsSolve)*  &
                                activeTracersSurfaceRestoringValue(2,1:nCellsSolve)
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                   areaCell(1:nCellsSolve), &
                   restoringSaltFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                   sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                   maxes_tmp(variableIndex))
            deallocate (restoringSaltFlux)
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! landIceFreshwaterFlux
         variableIndex = variableIndex + 1
         if ( associated(landIceFreshwaterFlux) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         landIceFloatingArea, &
                         landIceFreshwaterFlux(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! accumulatedLandIceMass
         variableIndex = variableIndex + 1
         if ( associated(accumulatedLandIceMass) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         landIceFloatingArea, &
                         accumulatedLandIceMass(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! accumulatedLandIceHeat
         variableIndex = variableIndex + 1
         if ( associated(accumulatedLandIceHeat) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         landIceFloatingArea, &
                         accumulatedLandIceHeat(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)

         ! accumulatedLandIceFrazilMass
         variableIndex = variableIndex + 1
         if ( associated(accumulatedLandIceFrazilMass) ) then
            call ocn_compute_field_area_weighted_local_stats_surface(dminfo, nCellsSolve, &
                         landIceFloatingArea, &
                         accumulatedLandIceFrazilMass(1:nCellsSolve), sums_tmp(variableIndex), &
                         sumSquares_tmp(variableIndex), mins_tmp(variableIndex), &
                         maxes_tmp(variableIndex))
         else
            sums_tmp(variableIndex) = 0.0_RKIND
            sumSquares_tmp(variableIndex) = 0.0_RKIND
            mins_tmp(variableIndex) = 0.0_RKIND
            maxes_tmp(variableIndex) = 0.0_RKIND
         end if

         sums(variableIndex) = sums(variableIndex) + sums_tmp(variableIndex)
         sumSquares(variableIndex) = sumSquares(variableIndex) + sumSquares_tmp(variableIndex)
         mins(variableIndex) = min(mins(variableIndex), mins_tmp(variableIndex))
         maxes(variableIndex) = max(maxes(variableIndex), maxes_tmp(variableIndex))
         verticalSumMins(variableIndex) = mins(variableIndex)
         verticalSumMaxes(variableIndex) = maxes(variableIndex)


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

         nSums = nSums + 1
         sums(nSums) = sums(nSums) + sum(landIceFloatingArea)

         localCFL = 0.0_RKIND
         do iEdge = 1,nEdgesSolve
            localCFL = max(localCFL, maxval(dt*abs(normalVelocity(1:maxLevelEdgeTop(iEdge),iEdge))/dcEdge(iEdge)))
         end do
         nMaxes = nMaxes + 1
         maxes(nMaxes) = localCFL
         do i = 1, nVariables
            mins(nMins+i) = min(mins(nMins+i),verticalSumMins(i))
            maxes(nMaxes+i) = max(maxes(nMaxes+i),verticalSumMaxes(i))
         end do

         nMins = nMins + nVariables
         nMaxes = nMaxes + nVariables

         deallocate(areaEdge)

         deallocate(landIceFloatingArea)

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
      landIceFloatingAreaSum = reductions(nVariables+7)
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

      ! layerThickness
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal*nVertLevels))

      ! set current ocean volume to be used for volume change diagnostic
      currentVolume = sums(variableIndex)

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

      ! vertAleTransportTop vertical velocity
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

      ! active tracers
      do iTracer=1,num_activeTracers
        variableIndex = variableIndex + 1
        averages(variableIndex) = sums(variableIndex)/volumeCellGlobal
        rms(variableIndex) = sqrt(sumSquares(variableIndex)/volumeCellGlobal)
      enddo

      ! layerThickness from previous timestep
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal*nVertLevels))

      ! calculate change in ocean volume
      totalVolumeChange = currentVolume - sums(variableIndex)

      ! frazilLayerThicknessTendency
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal*nVertLevels)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal*nVertLevels))

      ! start accumulating fresh water inputs
      netFreshwaterInput = dt * sums(variableIndex)

      ! evaporationFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! rainFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! snowFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! seaIceFreshWaterFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! riverRunoffFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! iceRunoffFlux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! temperature flux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! salinity flux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! salinity restoring flux
      variableIndex = variableIndex + 1
      averages(variableIndex) = sums(variableIndex)/(areaCellGlobal)
      rms(variableIndex) = sqrt(sumSquares(variableIndex)/(areaCellGlobal))

      ! landIceFreshwaterFlux
      variableIndex = variableIndex + 1
      if (associated(landIceFreshwaterFlux)) then
         averages(variableIndex) = sums(variableIndex)/(landIceFloatingAreaSum)
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/(landIceFloatingAreaSum))
      else
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
      end if

      ! continue accumulating fresh water inputs
      netFreshwaterInput = netFreshwaterInput + sums(variableIndex) * dt/rho_sw

      ! accumulatedLandIceMass
      variableIndex = variableIndex + 1
      if (associated(accumulatedLandIceMass)) then
         averages(variableIndex) = sums(variableIndex)/(landIceFloatingAreaSum)
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/(landIceFloatingAreaSum))
      else
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
      end if

      ! accumulatedLandIceHeat
      variableIndex = variableIndex + 1
      if (associated(accumulatedLandIceHeat)) then
         averages(variableIndex) = sums(variableIndex)/(landIceFloatingAreaSum)
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/(landIceFloatingAreaSum))
      else
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
      end if

      ! accumulatedLandIceFrazilMass
      variableIndex = variableIndex + 1
      if (associated(accumulatedLandIceFrazilMass)) then
         averages(variableIndex) = sums(variableIndex)/(landIceFloatingAreaSum)
         rms(variableIndex) = sqrt(sumSquares(variableIndex)/(landIceFloatingAreaSum))
      else
         averages(variableIndex) = 0.0_RKIND
         rms(variableIndex) = 0.0_RKIND
      end if

      ! calculate fresh water conservation check quantities
      absoluteFreshWaterConservation = totalVolumeChange - netFreshwaterInput
      if (totalVolumeChange == 0.0_RKIND) then
         relativeFreshWaterConservation = 0.0_RKIND
      else
         relativeFreshWaterConservation = (totalVolumeChange - netFreshwaterInput)/totalVolumeChange
      endif

      minGlobalStats(1:nVariables) =  mins(1:nVariables)
      maxGlobalStats(1:nVariables) =  maxes(1:nVariables)
      sumGlobalStats(1:nVariables) =  sums(1:nVariables)

      call mpas_pool_get_config(domain % configs, 'config_AM_globalStats_text_file', config_AM_globalStats_text_file)
      call mpas_pool_get_config(domain % configs, 'config_AM_globalStats_directory', config_AM_globalStats_directory)
      if (config_AM_globalStats_text_file) then

         call mpas_pool_get_subpool(domain % blocklist % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)

         ! write out the data to files
         if (dminfo % my_proc_id == IO_NODE) then
            fileID = getFreeUnit()
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_min.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') daysSinceStartOfSim, mins(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_max.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') daysSinceStartOfSim, maxes(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_sum.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') daysSinceStartOfSim, sums(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_rms.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') daysSinceStartOfSim, rms(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_avg.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') daysSinceStartOfSim, averages(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_time.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(a,2es24.14)') trim(xtime), daysSinceStartOfSim, CFLNumberGlobal
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_colmin.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') verticalSumMins(1:nVariables)
            close (fileID)
            open(fileID,file=trim(config_AM_globalStats_directory)//'/stats_colmax.txt',STATUS='UNKNOWN', POSITION='append')
            write (fileID,'(100es24.14)') verticalSumMaxes(1:nVariables)
            close (fileID)
         end if

      endif

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

      localSum = 0.0_RKIND
      localRMS = 0.0_RKIND
      localMin =  1.0e34_RKIND
      localMax = -1.0e34_RKIND
      localVertSumMin =  1.0e34_RKIND
      localVertSumMax = -1.0e34_RKIND

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

   subroutine ocn_compute_field_volume_weighted_local_stats_max_level(dminfo, nVertLevels, nElements, maxLevel, areas, & !{{{
                              layerThickness, field, localSum, localRMS, localMin, &
                              localMax, localVertSumMin, localVertSumMax)

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

      localSum = 0.0_RKIND
      localRMS = 0.0_RKIND
      localMin =  1.0e34_RKIND
      localMax = -1.0e34_RKIND
      localVertSumMin =  1.0e34_RKIND
      localVertSumMax = -1.0e34_RKIND

      do elementIndex = 1, nElements
        thicknessWeightedColSum = sum(layerThickness(1:maxLevel(elementIndex),elementIndex) &
                                * field(1:maxLevel(elementIndex),elementIndex))
        localSum = localSum + areas(elementIndex) * thicknessWeightedColSum
        thicknessWeightedColRMS = sum(layerThickness(1:maxLevel(elementIndex),elementIndex) &
                                * field(1:maxLevel(elementIndex),elementIndex)**2)
        localRMS = localRMS + areas(elementIndex) * thicknessWeightedColRMS
        localMin = min(localMin,minval(field(1:maxLevel(elementIndex),elementIndex)))
        localMax = max(localMax,maxval(field(1:maxLevel(elementIndex),elementIndex)))
        localVertSumMin = min(localVertSumMin,thicknessWeightedColSum)
        localVertSumMax = max(localVertSumMax,thicknessWeightedColSum)
      end do

   end subroutine ocn_compute_field_volume_weighted_local_stats_max_level!}}}

   subroutine ocn_compute_field_area_weighted_local_stats_surface(dminfo, nElements, areas, field, &!{{{
      localSum, localRMS, localMin, localMax)

      implicit none

      type (dm_info), intent(in) :: dminfo
      integer, intent(in) :: nElements
      real (kind=RKIND), dimension(nElements), intent(in) :: areas, field
      real (kind=RKIND), intent(out) :: localSum, localRMS, localMin, localMax

      integer :: elementIndex

      localSum = 0.0_RKIND
      localRMS = 0.0_RKIND
      localMin =  1.0e34_RKIND
      localMax = -1.0e34_RKIND

      do elementIndex = 1, nElements
        localSum = localSum + areas(elementIndex) * field(elementIndex)
        localRMS = localRMS + areas(elementIndex) * field(elementIndex)**2
        localMin = min(localMin,field(elementIndex))
        localMax = max(localMax,field(elementIndex))
      end do

   end subroutine ocn_compute_field_area_weighted_local_stats_surface!}}}

end module ocn_global_stats

! vim: foldmethod=marker
