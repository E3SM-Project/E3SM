! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_surface_area_weighted_averages
!
!> \brief MPAS ocean analysis member: areal-{min,max,avg} of 2D fields
!> \author Todd Ringler
!> \date   April 24, 2015
!> \details
!>  MPAS ocean analysis member: surface-area-weighted averages
!
!-----------------------------------------------------------------------

module ocn_surface_area_weighted_averages

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

   public :: ocn_init_surface_area_weighted_averages, &
             ocn_compute_surface_area_weighted_averages, &
             ocn_restart_surface_area_weighted_averages, &
             ocn_finalize_surface_area_weighted_averages

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_surface_area_weighted_averages
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_surface_area_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_init_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_compute_surface_area_weighted_averages
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_surface_area_weighted_averages(domain, timeLevel, err)!{{{

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

      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: surfaceAreaWeightedAveragesAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool

      real (kind=RKIND), dimension(:,:), pointer :: minValueWithinOceanRegion
      real (kind=RKIND), dimension(:,:), pointer :: maxValueWithinOceanRegion
      real (kind=RKIND), dimension(:,:), pointer :: avgValueWithinOceanRegion

      ! pointers to data in pools to be analyzed
      ! bulkForcing Pkg
      real (kind=RKIND), dimension(:), pointer :: latentHeatFlux
      real (kind=RKIND), dimension(:), pointer :: sensibleHeatFlux
      real (kind=RKIND), dimension(:), pointer :: longWaveHeatFluxUp
      real (kind=RKIND), dimension(:), pointer :: longWaveHeatFluxDown
      real (kind=RKIND), dimension(:), pointer :: seaIceHeatFlux
      real (kind=RKIND), dimension(:), pointer :: shortWaveHeatFlux
      real (kind=RKIND), dimension(:), pointer :: evaporationFlux
      real (kind=RKIND), dimension(:), pointer :: seaIceFreshWaterFlux
      real (kind=RKIND), dimension(:), pointer :: riverRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: iceRunoffFlux
      real (kind=RKIND), dimension(:), pointer :: rainFlux
      real (kind=RKIND), dimension(:), pointer :: snowFlux

      ! frazilIce Pkg
      real (kind=RKIND), dimension(:), pointer :: seaIceEnergy

      real (kind=RKIND), dimension(:), pointer :: surfaceThicknessFlux
      real (kind=RKIND), dimension(:,:), pointer :: activeTracersSurfaceFlux
      real (kind=RKIND), dimension(:), pointer ::  penetrativeTemperatureFlux

      real (kind=RKIND), dimension(:), pointer :: seaIceSalinityFlux
      real (kind=RKIND), dimension(:), pointer :: surfaceStressMagnitude
      real (kind=RKIND), dimension(:), pointer :: windStressZonal
      real (kind=RKIND), dimension(:), pointer :: windStressMeridional
      real (kind=RKIND), dimension(:), pointer :: atmosphericPressure
      real (kind=RKIND), dimension(:), pointer :: ssh
      real (kind=RKIND), dimension(:), pointer :: boundaryLayerDepth
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      ! pointers to data in mesh pool
      integer, pointer :: nCells, nCellsSolve, nSfcAreaWeightedAvgFields, nOceanRegions
      integer, pointer :: indexTemperature, indexSalinity
      real (kind=RKIND), dimension(:), pointer ::  areaCell, lonCell, latCell

      ! scratch space
      type(field2DReal), pointer :: workArrayField
      real (kind=RKIND), dimension(:,:), pointer :: workArray
      type(field1DReal), pointer :: workMaskField, workMinField, workMaxField, workSumField
      real (kind=RKIND), dimension(:), pointer :: workMask, workMin, workMax, workSum

      ! local variables
      integer :: iDataField, nDefinedDataFields
      integer :: iCell, iRegion, iTracer, err_tmp

      ! package flag
      logical, pointer :: activeTracersBulkRestoringPKG
      logical, pointer :: frazilIcePkgActive

      ! buffers data for message passaging
      integer :: kBuffer, kBufferLength
      real (kind=RKIND), dimension(:), allocatable :: workBufferSum, workBufferSumReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMin, workBufferMinReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMax, workBufferMaxReduced

      ! assume no error
      err = 0

      ! get status of other packages
      call mpas_pool_get_package(ocnPackages, 'activeTracersBulkRestoringPKGActive', activeTracersBulkRestoringPKG)
      call mpas_pool_get_package(ocnPackages, 'frazilIceActive', frazilIcePkgActive)

      ! set highest level pointer
      dminfo = domain % dminfo

      ! find the number of regions
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegions', nOceanRegions)

      ! find the number of data fields
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)

      ! allocate buffer for message passing
      kBuffer=0
      kBufferLength=nOceanRegions*nSfcAreaWeightedAvgFields
      allocate(workBufferSum(kBufferLength))
      allocate(workBufferMin(kBufferLength))
      allocate(workBufferMax(kBufferLength))
      allocate(workBufferSumReduced(kBufferLength))
      allocate(workBufferMinReduced(kBufferLength))
      allocate(workBufferMaxReduced(kBufferLength))
      workBufferSum=0.0_RKIND
      workBufferMin=0.0_RKIND
      workBufferMax=0.0_RKIND
      workBufferSumReduced=0.0_RKIND
      workBufferMinReduced=0.0_RKIND
      workBufferMaxReduced=0.0_RKIND

      ! loop over all ocean regions
      do iRegion=1,nOceanRegions

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAM', surfaceAreaWeightedAveragesAMPool)
      call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'minValueWithinOceanRegion', minValueWithinOceanRegion)
      call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'maxValueWithinOceanRegion', maxValueWithinOceanRegion)
      call mpas_pool_get_array(surfaceAreaWeightedAveragesAMPool, 'avgValueWithinOceanRegion', avgValueWithinOceanRegion)

      ! get pointers to scratch variables
      call mpas_pool_get_subpool(domain % blocklist % structs, 'surfaceAreaWeightedAveragesAMScratch', scratchPool)
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

      ! loop over blocks
      ! NOTE: code is not valid for multiple blocks !
      block => domain % blocklist
      do while (associated(block))
         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
         call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
         call mpas_pool_get_dimension(block % dimensions, 'nOceanRegions', nOceanRegions)
         call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
         call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)

         ! test to make sure the arrays are big enough
         nDefinedDataFields = size(avgValueWithinOceanRegion,dim=1)
         if (nDefinedDataFields > nSfcAreaWeightedAvgFields) then
             write (stderrUnit,*) 'Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields'
             write (stderrUnit,*) '     :  increase size of ocn_surface_area_weighted_averages scratch space'
             call mpas_dmpar_global_abort('MPAS-ocean: Abort: nDefinedDataFields > nSfcAreaWeightedAvgFields')
         endif

         ! get pointers to data that will be analyzed
         ! listed in the order in which the fields appear in {avg,min,max}SurfaceStatistics
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'latentHeatFlux', latentHeatFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxUp', longWaveHeatFluxUp)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'longWaveHeatFluxDown', longWaveHeatFluxDown)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceHeatFlux', seaIceHeatFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'shortWaveHeatFlux', shortWaveHeatFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'evaporationFlux', evaporationFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceFreshWaterFlux', seaIceFreshWaterFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'riverRunoffFlux', riverRunoffFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'iceRunoffFlux', iceRunoffFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
         if (frazilIcePkgActive)   call mpas_pool_get_array(forcingPool, 'seaIceEnergy', seaIceEnergy)
         call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
         call mpas_pool_get_array(tracersSurfaceFluxPool, 'activeTracersSurfaceFlux', activeTracersSurfaceFlux)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'seaIceSalinityFlux', seaIceSalinityFlux)
         call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
         if (activeTracersBulkRestoringPKG) call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)
         call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
         call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth',boundaryLayerDepth)

         ! compute mask
         call compute_mask(nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)

         ! copy data into work array
         workArray( :,:) = 0.0_RKIND
         workArray( 1,:) = workMask(:)
         workArray( 2,:) = areaCell(:)
         if (activeTracersBulkRestoringPKG) workArray( 3,:) = latentHeatFlux(:)
         if (activeTracersBulkRestoringPKG) workArray( 4,:) = sensibleHeatFlux(:)
         if (activeTracersBulkRestoringPKG) workArray( 5,:) = longWaveHeatFluxUp(:)
         if (activeTracersBulkRestoringPKG) workArray( 6,:) = longWaveHeatFluxDown(:)
         if (activeTracersBulkRestoringPKG) workArray( 7,:) = seaIceHeatFlux(:)
         if (activeTracersBulkRestoringPKG) workArray( 8,:) = shortWaveHeatFlux(:)
         if (activeTracersBulkRestoringPKG) workArray( 9,:) = evaporationFlux(:)
         if (activeTracersBulkRestoringPKG) workArray(10,:) = seaIceFreshWaterFlux(:)
         if (activeTracersBulkRestoringPKG) workArray(11,:) = riverRunoffFlux(:)
         if (activeTracersBulkRestoringPKG) workArray(12,:) = iceRunoffFlux(:)
         if (activeTracersBulkRestoringPKG) workArray(13,:) = rainFlux(:)
         if (activeTracersBulkRestoringPKG) workArray(14,:) = snowFlux(:)
         if (frazilIcePkgActive) workArray(15,:) = seaIceEnergy(:)
         workArray(16,:) = surfaceThicknessFlux(:)
         workArray(17,:) = activeTracersSurfaceFlux(indexTemperature,:)
         workArray(18,:) = activeTracersSurfaceFlux(indexSalinity,:)
         if (activeTracersBulkRestoringPKG) workArray(19,:) = seaIceSalinityFlux(:)
         workArray(20,:) = surfaceStressMagnitude(:)
         if (activeTracersBulkRestoringPKG) workArray(21,:) = windStressZonal(:)
         if (activeTracersBulkRestoringPKG) workArray(22,:) = windStressMeridional(:)
         workArray(23,:) = atmosphericPressure(:)
         workArray(24,:) = ssh(:)
         workArray(25,:) = activeTracers(indexTemperature,1,:)
         workArray(26,:) = activeTracers(indexSalinity,1,:)
         workArray(27,:) = boundaryLayerDepth(:)

         ! build net heat, salinity and fresh water budget
         ! net heat into ocean = latentHeatFlux + sensibleHeatFlux + longWaveHeatFluxUp + longWaveHeatFluxDown
         !                     + shortWaveHeatFlux + seaIceHeatFlux + (?seaIceEnergy?)
         ! net salinity into ocean = seaIceSalinityFlux
         ! net freshwater into ocean = evaporationFlux + seaIceFreshWaterFlux + riverRunoffFlux + iceRunoffFlux
         !                           + rainFlux + snowFlux + (?seaIceEnergy?)
         if (activeTracersBulkRestoringPKG) then
            workArray(28,:) =   latentHeatFlux(:) &
                              + sensibleHeatFlux(:) &
                              + longWaveHeatFluxUp(:) &
                              + longWaveHeatFluxDown(:) &
                              + shortWaveHeatFlux(:)  &
                              + seaIceHeatFlux(:)
                     !        + seaIceEnergy
            workArray(29,:) =   seaIceSalinityFlux(:)
            workArray(30,:) =   evaporationFlux(:) &
                              + seaIceFreshWaterFlux(:) &
                              + riverRunoffFlux(:) &
                              + iceRunoffFlux(:) &
                              + rainFlux(:) &
                              + snowFlux(:)
                     !        + seaIceEnergy(:)
         end if

         call compute_statistics(nDefinedDataFields, nCellsSolve, workArray, workMask, workMin, workMax, workSum)

         ! store data in buffer in order to allow only three dmpar calls
         do iDataField=1,nDefinedDataFields
            kBuffer = kBuffer+1
            workBufferSum(kBuffer) = workSum(iDataField)
            workBufferMin(kBuffer) = workMin(iDataField)
            workBufferMax(kBuffer) = workMax(iDataField)
         enddo

         block => block % next
      end do

      end do ! iRegion

      ! communication
      call mpas_dmpar_sum_real_array(dminfo, kBufferLength, workBufferSum, workBufferSumReduced )
      call mpas_dmpar_min_real_array(dminfo, kBufferLength, workBufferMin, workBufferMinReduced )
      call mpas_dmpar_max_real_array(dminfo, kBufferLength, workBufferMax, workBufferMaxReduced )

      ! unpack the buffer into intent(out) of this analysis member
      kBuffer=0
      do iRegion=1,nOceanRegions
        do iDataField=1,nDefinedDataFields
           kBuffer = kBuffer+1
           avgValueWithinOceanRegion(iDataField,iRegion)=workBufferSumReduced(kBuffer)
           minValueWithinOceanRegion(iDataField,iRegion)=workBufferMinReduced(kBuffer)
           maxValueWithinOceanRegion(iDataField,iRegion)=workBufferMaxReduced(kBuffer)
        enddo
      enddo

      ! normalize averages
      do iRegion=1,nOceanRegions
         ! normalize all field by total area
         do iDataField=3,nDefinedDataFields
            avgValueWithinOceanRegion(iDataField,iRegion) = avgValueWithinOceanRegion(iDataField,iRegion) &
                                                          / max(avgValueWithinOceanRegion(2,iRegion),1.0e-8_RKIND)
         enddo
         ! normalize total area by number of cells in region
         avgValueWithinOceanRegion(2,iRegion) = avgValueWithinOceanRegion(2,iRegion) &
                                              / max(avgValueWithinOceanRegion(1,iRegion),1.0e-8_RKIND)
      enddo

      ! deallocate scratch fields
      call mpas_deallocate_scratch_field(workArrayField, .true.)
      call mpas_deallocate_scratch_field(workMaskField, .true.)
      call mpas_deallocate_scratch_field(workMinField, .true.)
      call mpas_deallocate_scratch_field(workMaxField, .true.)
      call mpas_deallocate_scratch_field(workSumField, .true.)

      ! deallocate buffers
      deallocate(workBufferSumReduced)
      deallocate(workBufferMinReduced)
      deallocate(workBufferMaxReduced)

   contains

   subroutine compute_mask(nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)!{{{
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to
   ! allow for min/max/avg to represent specific regions of the ocean domain
   !
   ! NOTE: computes_mask is temporary. workMask should be intent(in) to this entire module !
   !
   integer, intent(in) :: nCells, nCellsSolve, iRegion
   real(kind=RKIND), dimension(:), intent(in) :: lonCell, latCell
   real(kind=RKIND), dimension(:), intent(out) :: workMask
   integer :: iCell
   real(kind=RKIND) :: dtr

   dtr = 4.0_RKIND*atan(1.0_RKIND) / 180.0_RKIND
   workMask(:) = 0.0_RKIND
   do iCell=1,nCellsSolve
      workMask(iCell) = 1.0_RKIND
   enddo

   if (iRegion.eq.1) then
      ! Arctic
      do iCell=1,nCellsSolve
        if(latCell(iCell).lt. 60.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.2) then
      ! Equatorial
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt. 15.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt.-15.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.3) then
      ! Southern Ocean
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.-50.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.4) then
      ! Nino 3
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.210.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.270.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.5) then
      ! Nino 4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.160.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.210.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   elseif (iRegion.eq.6) then
      ! Nino 3.4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.190.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.240.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
   else
      ! global (do nothing!)
   endif

   end subroutine compute_mask!}}}


   subroutine compute_statistics(nDefinedDataFields, nCellsSolve, workArray, workMask, workMin, workMax, workSum)!{{{
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCellsSolve
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   real(kind=RKIND), dimension(:), intent(in) :: workMask
   real(kind=RKIND), dimension(:), intent(out) :: workMin, workMax, workSum
   integer :: iCell, iData

   workSum = 0.0_RKIND
   do iCell=1,nCellsSolve
    workSum(1) = workSum(1) + workMask(iCell)
    workSum(2) = workSum(2) + workArray(2,iCell)*workMask(iCell)
    do iData=3,nDefinedDataFields
      workSum(iData) = workSum(iData) + workArray(2,iCell)*workArray(iData,iCell)*workMask(iCell)
    enddo
   enddo

   do iData=1,nDefinedDataFields
      workMin(iData) = minval(workArray(iData,1:nCellsSolve),workMask(1:nCellsSolve)>0.5_RKIND)
      workMax(iData) = maxval(workArray(iData,1:nCellsSolve),workMask(1:nCellsSolve)>0.5_RKIND)
   enddo

   end subroutine compute_statistics!}}}

   end subroutine ocn_compute_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_restart_surface_area_weighted_averages
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_surface_area_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_restart_surface_area_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_finalize_surface_area_weighted_averages
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_surface_area_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_finalize_surface_area_weighted_averages!}}}

end module ocn_surface_area_weighted_averages

! vim: foldmethod=marker
