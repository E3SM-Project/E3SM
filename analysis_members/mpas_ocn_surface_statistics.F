! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_statistics
!
!> \brief MPAS ocean analysis member: surface statistics
!> \author Todd Ringler
!> \date   April 24, 2015
!> \details
!>  MPAS ocean analysis member: surface statistics
!
!-----------------------------------------------------------------------

module ocn_surface_statistics

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

   public :: ocn_setup_packages_surface_statistics, &
             ocn_init_surface_statistics, &
             ocn_compute_surface_statistics, &
             ocn_restart_surface_statistics, &
             ocn_finalize_surface_statistics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type (timer_node), pointer :: amSurfaceStatisticsTimer

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_setup_packages_surface_statistics
!
!> \brief   Set up packages for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine is intended to configure the packages for this MPAS
!>   ocean analysis member
!
!-----------------------------------------------------------------------

   subroutine ocn_setup_packages_surface_statistics(configPool, packagePool, err)!{{{

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

      logical, pointer :: amSurfaceStatisticsPkgActive

      err = 0

      call mpas_pool_get_package(packagePool, 'amSurfaceStatisticsPkgActive', amSurfaceStatisticsPkgActive)

      ! turn on package for this analysis member
      amSurfaceStatisticsPkgActive = .true.

      write(6,*) ' active surface statistics'

   end subroutine ocn_setup_packages_surface_statistics!}}}


!***********************************************************************
!
!  routine ocn_init_surface_statistics
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all initializations required for the 
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_surface_statistics(domain, err)!{{{

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

      write(6,*) ' surf stat inited'

   end subroutine ocn_init_surface_statistics!}}}

!***********************************************************************
!
!  routine ocn_compute_surface_statistics
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_surface_statistics(domain, timeLevel, err)!{{{

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

   !  integer :: nCellsGlobal, nEdgesGlobal, nVerticesGlobal
   !  integer :: elementIndex, variableIndex, nVariables, nSums, nMaxes, nMins
   !  integer :: timeYYYY, timeMM, timeDD, timeH, timeM, timeS
   !  character*1 timeChar
   !  character (len=StrKIND), pointer :: xtime
   !  real (kind=RKIND) :: volumeCellGlobal, volumeEdgeGlobal, CFLNumberGlobal, localCFL, localSum, areaCellGlobal, areaEdgeGlobal, areaTriangleGlobal, time_days
   !  real (kind=RKIND), dimension(kMaxVariables) :: sumSquares, reductions, sums, mins, maxes
   !  real (kind=RKIND), dimension(kMaxVariables) :: sums_tmp, sumSquares_tmp, mins_tmp, maxes_tmp, averages_tmp, verticalSumMins_tmp, verticalSumMaxes_tmp
   !  real (kind=RKIND), dimension(:,:), allocatable :: enstrophy, normalizedAbsoluteVorticity, workArray

!---------------------------------------

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: amSurfaceStatisticsPool
      type (mpas_pool_type), pointer :: amSurfaceStatisticsLocationPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      real (kind=RKIND), dimension(:), pointer ::   minSurfaceStatistics, &
                                                    maxSurfaceStatistics, &
                                                    avgSurfaceStatistics

      real (kind=RKIND), dimension(:,:), pointer :: locationMinSurfaceStatistics, &
                                                    locationMaxSurfaceStatistics


      ! pointers to data in pools to be analyzed
      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:,:), pointer :: kineticEnergyCell
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      ! pointers to data in mesh pool
      integer, pointer :: nCellsSolve
      real (kind=RKIND), dimension(:), pointer ::  areaCell

      ! scratch space
      type(field2DReal), pointer :: workArrayField
      real (kind=RKIND), dimension(:,:), pointer :: workArray
      type(field1DReal), pointer :: workMinField, workMaxField, workAvgField
      real (kind=RKIND), dimension(:), pointer :: workMin, workMax, workAvg

      ! local variables
      integer, parameter :: maxVariables = 128
      integer :: iCell, iVariable, iTracer, err_tmp

      ! package flag
      logical, pointer :: amSurfaceStatisticsPkgActive

        write(94,*) ' entered compute surface statistics'
        write(6,*) ' entered compute surface statistics'

      ! assume no error
      err = 0

      ! check to see if package is on
      call mpas_pool_get_package(ocnPackages, 'amSurfaceStatisticsPkgActive', amSurfaceStatisticsPkgActive)
      if ( .not. amSurfaceStatisticsPkgActive ) return

      ! set highest level pointer
      dminfo = domain % dminfo

      ! start timer
      call mpas_timer_start("compute_surface_statistics", .false., amSurfaceStatisticsTimer)

      ! get pointers to analysis member arrays
      !
      ! TDR: should this be inside the "do while block" loop?
      !
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatistics', amSurfaceStatisticsPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatisticsLocation', amSurfaceStatisticsLocationPool)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'minSurfaceStatistics', minSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'maxSurfaceStatistics', maxSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'avgSurfaceStatistics', avgSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsLocationPool, 'locationMinSurfaceStatistics', locationMinSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsLocationPool, 'locationMaxSurfaceStatistics', locationMaxSurfaceStatistics)

  ! HERE
  !   sums = 0.0
  !   sumSquares = 0.0
  !   mins = 1.0e34
  !   maxes = -1.0e34
  !   averages = 0.0
  !   verticalSumMins = 1.0e34
  !   verticalSumMaxes = -1.0e34
  !   reductions = 0.0
  ! HERE

      block => domain % blocklist
      do while (associated(block))
         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'surfaceStatisticsScratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

    ! HERE -- how to get num_tracers
    !    call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         ! get pointers to data that will be analyzed
         call mpas_pool_get_array(statePool, 'ssh', ssh)  
       ! call mpas_pool_get_array(statePool, 'tracers', tracers)         
       ! call mpas_pool_get_array(diagnosticsPool, 'density', density)
       ! call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
       ! call mpas_pool_get_array(diagnosticsPool, 'layerThickness', layerThickness)
       ! call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)

         ! get pointers to scratch variables
         call mpas_pool_get_field(scratchPool, 'workArray', workArrayField)
         call mpas_pool_get_field(scratchPool, 'workMin', workMinField)
         call mpas_pool_get_field(scratchPool, 'workMax', workMaxField)
         call mpas_pool_get_field(scratchPool, 'workAvg', workAvgField)
         call mpas_allocate_scratch_field(workArrayField, .true.)
         call mpas_allocate_scratch_field(workMinField, .true.)
         call mpas_allocate_scratch_field(workMaxField, .true.)
         call mpas_allocate_scratch_field(workAvgField, .true.)
         workArray => workArrayField % array
         workMin => workMinField % array
         workMax => workMaxField % array
         workAvg => workAvgField % array

         workArray(1,:) = areaCell(:)
         workArray(2,:) = ssh(:)

    !    call compute_surface_statistics(workArray, work2)
         workMin(1) = minval(workArray(1,:))
         workMin(2) = minval(workArray(2,:))

         workMin = 12.5

   !     call mpas_deallocate_scratch_field(workArrayField, .true.)
   !     call mpas_deallocate_scratch_field(workMinField, .true.)
   !     call mpas_deallocate_scratch_field(workMaxField, .true.)
   !     call mpas_deallocate_scratch_field(workAvgField, .true.)

         block => block % next
      end do

      ! communication
   !  call mpas_dmpar_sum_real_array(dminfo, nSums, sums(1:nSums), reductions(1:nSums))
   !  call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
   !  call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
   !  call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

      ! HERE save data that will be written to file
      minSurfaceStatistics(1) = workMin(1)
      minSurfaceStatistics(2) = workMin(2)

      call mpas_deallocate_scratch_field(workArrayField, .true.)
      call mpas_deallocate_scratch_field(workMinField, .true.)
      call mpas_deallocate_scratch_field(workMaxField, .true.)
      call mpas_deallocate_scratch_field(workAvgField, .true.)

   !  minGlobalStats(1:nVariables) =  mins(1:nVariables)
   !  maxGlobalStats(1:nVariables) =  maxes(1:nVariables)
   !  sumGlobalStats(1:nVariables) =  sums(1:nVariables)

      call mpas_timer_stop("surface_statistics", amSurfaceStatisticsTimer)

   end subroutine ocn_compute_surface_statistics!}}}

!***********************************************************************
!
!  routine ocn_restart_surface_statistics
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_surface_statistics(domain, err)!{{{

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

   end subroutine ocn_restart_surface_statistics!}}}

!***********************************************************************
!
!  routine ocn_finalize_surface_statistics
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details 
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_surface_statistics(domain, err)!{{{

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

   end subroutine ocn_finalize_surface_statistics!}}}

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

end module ocn_surface_statistics

! vim: foldmethod=marker
