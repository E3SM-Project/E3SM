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
      integer, pointer :: nCells, nCellsSolve, nSurfaceStatisticsFields
      real (kind=RKIND), dimension(:), pointer ::  areaCell

      ! scratch space
      type(field2DReal), pointer :: workArrayField
      real (kind=RKIND), dimension(:,:), pointer :: workArray
      type(field1DReal), pointer :: workMaskField, workMinField, workMaxField, workSumField
      real (kind=RKIND), dimension(:), pointer :: workMask, workMin, workMax, workSum

      ! local variables
      integer, parameter :: kTarget = 1
      integer :: nDefinedDataFields
      integer :: iCell, iVariable, iTracer, err_tmp

      ! package flag
      logical, pointer :: amSurfaceStatisticsPkgActive

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
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatistics', amSurfaceStatisticsPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'amSurfaceStatisticsLocation', amSurfaceStatisticsLocationPool)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'minSurfaceStatistics', minSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'maxSurfaceStatistics', maxSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsPool, 'avgSurfaceStatistics', avgSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsLocationPool, 'locationMinSurfaceStatistics', locationMinSurfaceStatistics)
      call mpas_pool_get_array(amSurfaceStatisticsLocationPool, 'locationMaxSurfaceStatistics', locationMaxSurfaceStatistics)

      ! get pointers to scratch variables
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceStatisticsScratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'workArray', workArrayField)
      call mpas_pool_get_field(scratchPool, 'workMask', workMaskField)
      call mpas_pool_get_field(scratchPool, 'workMin', workMinField)
      call mpas_pool_get_field(scratchPool, 'workMax', workMaxField)
      call mpas_pool_get_field(scratchPool, 'workSum', workSumField)
      call mpas_allocate_scratch_field(workArrayField, .true.)
      call mpas_allocate_scratch_field(workMaskField, .true.)
      call mpas_allocate_scratch_field(workMinField, .true.)
      call mpas_allocate_scratch_field(workMaxField, .true.)
      call mpas_allocate_scratch_field(workSumField, .true.)
      workArray => workArrayField % array
      workMask => workMaskField % array
      workMin => workMinField % array
      workMax => workMaxField % array
      workSum => workSumField % array

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
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nSurfaceStatisticsFields', nSurfaceStatisticsFields)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

         ! test to make sure the arrays are big enough
         nDefinedDataFields = size(avgSurfaceStatistics)
         if (nDefinedDataFields > nSurfaceStatisticsFields) then
              ! call error/finalize here
         endif

    ! HERE -- how to get num_tracers
    !    call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         ! get pointers to data that will be analyzed
         call mpas_pool_get_array(statePool, 'ssh', ssh)  
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
       ! call mpas_pool_get_array(statePool, 'tracers', tracers)         
       ! call mpas_pool_get_array(diagnosticsPool, 'density', density)
       ! call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
       ! call mpas_pool_get_array(diagnosticsPool, 'layerThickness', layerThickness)

         ! copy data into work array
         workArray(1,:) = areaCell(:)
         workArray(2,:) = ssh(:)
         workArray(3,:) = kineticEnergyCell(kTarget,:)

         call compute_mask(nCells, nCellsSolve, workMask)
         call compute_surface_statistics(nDefinedDataFields, nCells, workArray, workMask, workMin, workMax, workSum)

         ! temporary
    !    workMin(1) = minval(workArray(1,:))
    !    workMin(2) = minval(workArray(2,:))
    !    workMin(3) = minval(workArray(3,:))
    !    workMax(1) = maxval(workArray(1,:))
    !    workMax(2) = maxval(workArray(2,:))
    !    workMax(3) = maxval(workArray(3,:))

         write(6,*) 'workMax ', workMax(:)

         block => block % next
      end do

      ! communication
      call mpas_dmpar_sum_real_array(dminfo, nDefinedDataFields, workSum(1:nDefinedDataFields), avgSurfaceStatistics(1:nDefinedDataFields))
      call mpas_dmpar_min_real_array(dminfo, nDefinedDataFields, workMin(1:nDefinedDataFields), minSurfaceStatistics(1:nDefinedDataFields))
      call mpas_dmpar_min_real_array(dminfo, nDefinedDataFields, workMax(1:nDefinedDataFields), maxSurfaceStatistics(1:nDefinedDataFields))

      ! normalize averages
      avgSurfaceStatistics(:) = avgSurfaceStatistics(:) / avgSurfaceStatistics(1)

      write(6,*) 'avgSurfaceStatistics ', avgSurfaceStatistics(:)
      write(6,*) 'minSurfaceStatistics ', minSurfaceStatistics(:)
      write(6,*) 'maxSurfaceStatistics ', maxSurfaceStatistics(:)

      ! deallocate scratch fields
      call mpas_deallocate_scratch_field(workArrayField, .true.)
      call mpas_deallocate_scratch_field(workMaskField, .true.)
      call mpas_deallocate_scratch_field(workMinField, .true.)
      call mpas_deallocate_scratch_field(workMaxField, .true.)
      call mpas_deallocate_scratch_field(workSumField, .true.)

      ! stop timer
      call mpas_timer_stop("surface_statistics", amSurfaceStatisticsTimer)

   contains

   subroutine compute_mask(nCells, nCellsSolve, workMask)
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to 
   ! allow for min/max/avg to represent specific regions of the ocean domain
   integer, intent(in) :: nCells, nCellsSolve
   real(kind=RKIND), dimension(:), intent(out) :: workMask
   integer :: iCell

   workMask(:) = 1.0_RKIND
   do iCell=nCellsSolve+1,nCells
      workMask(iCell) = 0.0_RKIND
   enddo

   end subroutine compute_mask


   subroutine compute_surface_statistics(nDefinedDataFields, nCells, workArray, workMask, workMin, workMax, workSum)
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCells
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   real(kind=RKIND), dimension(:), intent(out) :: workMask, workMin, workMax, workSum
   integer :: iCell, iData

   workSum = 0.0
   do iCell=1,nCells
    workSum(1) = workSum(1) + workArray(1,iCell)
    workSum(2:nDefinedDataFields) = workSum(2:nDefinedDataFields) + workArray(1,iCell)*workArray(2:nDefinedDataFields,iCell)*workMask(iCell)
   enddo

   do iData=1,nDefinedDataFields
      workMin(iData) = minval(workArray(iData,:)*workMask(:))
      workMax(iData) = maxval(workArray(iData,:)*workMask(:))
   enddo

   end subroutine compute_surface_statistics

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
