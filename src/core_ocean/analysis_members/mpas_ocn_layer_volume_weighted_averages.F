! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_layer_volume_weighted_averages
!
!> \brief MPAS ocean analysis member: horizonal layer volume weighted averages at each vertical level
!> \author Todd Ringler
!> \date   April 24, 2015
!> \details
!>  MPAS ocean analysis member: layer_volume_weighted_averages
!
!-----------------------------------------------------------------------

module ocn_layer_volume_weighted_averages

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

   public :: ocn_init_layer_volume_weighted_averages, &
             ocn_compute_layer_volume_weighted_averages, &
             ocn_restart_layer_volume_weighted_averages, &
             ocn_finalize_layer_volume_weighted_averages

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_layer_volume_weighted_averages
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_layer_volume_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_init_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_compute_layer_volume_weighted_averages
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_layer_volume_weighted_averages(domain, timeLevel, err)!{{{

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
      type (mpas_pool_type), pointer :: layerVolumeWeightedAverageAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: forcingPool

      real (kind=RKIND), dimension(:,:,:), pointer :: minValueWithinOceanLayerRegion
      real (kind=RKIND), dimension(:,:,:), pointer :: maxValueWithinOceanLayerRegion
      real (kind=RKIND), dimension(:,:,:), pointer :: avgValueWithinOceanLayerRegion
      real (kind=RKIND), dimension(:,:),   pointer :: minValueWithinOceanVolumeRegion
      real (kind=RKIND), dimension(:,:),   pointer :: maxValueWithinOceanVolumeRegion
      real (kind=RKIND), dimension(:,:),   pointer :: avgValueWithinOceanVolumeRegion

      ! pointers to data in pools to be analyzed
      real (kind=RKIND), dimension(:,:),   pointer :: layerThickness
      real (kind=RKIND), dimension(:,:),   pointer :: density
      real (kind=RKIND), dimension(:,:),   pointer :: potentialDensity
      real (kind=RKIND), dimension(:,:),   pointer :: BruntVaisalaFreqTop
      real (kind=RKIND), dimension(:,:),   pointer :: velocityZonal
      real (kind=RKIND), dimension(:,:),   pointer :: velocityMeridional
      real (kind=RKIND), dimension(:,:),   pointer :: vertVelocityTop
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:),   pointer :: kineticEnergyCell
      real (kind=RKIND), dimension(:,:),   pointer :: relativeVorticityCell
      real (kind=RKIND), dimension(:,:),   pointer :: divergence

      ! pointers to data in mesh pool
      integer, pointer :: nVertLevels, nCells, nCellsSolve, nLayerVolWeightedAvgFields, nOceanRegionsTmp
      integer, pointer :: indexTemperature, indexSalinity
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer ::  areaCell, lonCell, latCell

      ! scratch space
      type(field2DReal), pointer :: workArrayField
      real (kind=RKIND), dimension(:,:), pointer :: workArray
      type(field1DReal), pointer :: workMaskField, workMinField, workMaxField, workSumField
      real (kind=RKIND), dimension(:), pointer :: workMask, workMin, workMax, workSum

      ! local variables
      integer :: iDataField, nDefinedDataFields
      integer :: iCell, iLevel, iRegion, iTracer, err_tmp

      ! package flag
      logical, pointer :: layerVolumeWeightedAverageAMPKGActive

      ! buffers data for message passaging
      integer :: kBuffer, kBufferLength
      real (kind=RKIND), dimension(:), allocatable :: workBufferSum, workBufferSumReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMin, workBufferMinReduced
      real (kind=RKIND), dimension(:), allocatable :: workBufferMax, workBufferMaxReduced

      ! assume no error
      err = 0

      ! set highest level pointer
      dminfo = domain % dminfo

      ! find the number of regions, number of data fields and number of vertical levels
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegionsTmp', nOceanRegionsTmp)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nLayerVolWeightedAvgFields', nLayerVolWeightedAvgFields)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      ! allocate buffer for message passing
      kBuffer=0
      kBufferLength=nOceanRegionsTmp*nLayerVolWeightedAvgFields*nVertLevels
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

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'layerVolumeWeightedAverageAM', layerVolumeWeightedAverageAMPool)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'minValueWithinOceanLayerRegion', minValueWithinOceanLayerRegion)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'maxValueWithinOceanLayerRegion', maxValueWithinOceanLayerRegion)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'avgValueWithinOceanLayerRegion', avgValueWithinOceanLayerRegion)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'minValueWithinOceanVolumeRegion', minValueWithinOceanVolumeRegion)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'maxValueWithinOceanVolumeRegion', maxValueWithinOceanVolumeRegion)
      call mpas_pool_get_array(layerVolumeWeightedAverageAMPool, 'avgValueWithinOceanVolumeRegion', avgValueWithinOceanVolumeRegion)

      ! loop over blocks
      ! NOTE: code is not valid for multiple blocks !
      block => domain % blocklist
      do while (associated(block))

         ! get pointers to scratch variables
         call mpas_pool_get_subpool(block % structs, 'layerVolumeWeightedAverageAMScratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'workArrayLayerVolume', workArrayField)
         call mpas_pool_get_field(scratchPool, 'workMaskLayerVolume', workMaskField)
         call mpas_pool_get_field(scratchPool, 'workMinLayerVolume', workMinField)
         call mpas_pool_get_field(scratchPool, 'workMaxLayerVolume', workMaxField)
         call mpas_pool_get_field(scratchPool, 'workSumLayerVolume', workSumField)
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

         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nLayerVolWeightedAvgFields', nLayerVolWeightedAvgFields)
         call mpas_pool_get_dimension(block % dimensions, 'nOceanRegionsTmp', nOceanRegionsTmp)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         ! test to make sure the arrays are big enough
         nDefinedDataFields = size(avgValueWithinOceanLayerRegion,dim=1)
         if (nDefinedDataFields > nLayerVolWeightedAvgFields) then
             write (stderrUnit,*) 'Abort: nDefinedDataFields > nLayerVolWeightedAvgFields'
             write (stderrUnit,*) '     :  increase size of ocn_layer_volume_weighted_averages scratch space'
             call mpas_dmpar_abort(dminfo)
         endif

         ! get pointers to data that will be analyzed
         ! listed in the order in which the fields appear in {avg,min,max}ValueWithinOceanLayerRegion
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(diagnosticsPool, 'density', density)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
         call mpas_pool_get_array(diagnosticsPool, 'velocityZonal', velocityZonal)
         call mpas_pool_get_array(diagnosticsPool, 'velocityMeridional', velocityMeridional)
         call mpas_pool_get_array(diagnosticsPool, 'vertVelocityTop', vertVelocityTop)
         call mpas_pool_get_array(statePool, 'tracers', tracers, 1)
         call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
         call mpas_pool_get_array(diagnosticsPool, 'relativeVorticityCell', relativeVorticityCell)
         call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence)

         ! initialize buffers
         workBufferSum(:) = 0.0_RKIND
         workBufferMin(:) = +1.0e20_RKIND
         workBufferMax(:) = -1.0e20_RKIND

         ! loop over all ocean regions
         do iRegion=1,nOceanRegionsTmp

         ! loop over the vertical
         do iLevel=1,nVertLevels

            ! compute mask
            call compute_mask(iLevel, maxLevelCell, nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)

            ! copy data into work array
            workArray( :,:) = 0.0
            workArray( 1,:) = workMask(:)
            workArray( 2,:) = areaCell(:)
            workArray( 3,:) = layerThickness(iLevel,:)
            workArray( 4,:) = density(iLevel,:)
            workArray( 5,:) = potentialDensity(iLevel,:)
            workArray( 6,:) = BruntVaisalaFreqTop(iLevel,:)
            workArray( 7,:) = velocityZonal(iLevel,:)
            workArray( 8,:) = velocityMeridional(iLevel,:)
            workArray( 9,:) = vertVelocityTop(iLevel,:)
            workArray(10,:) = tracers(indexTemperature,iLevel,:)
            workArray(11,:) = tracers(indexSalinity,iLevel,:)
            workArray(12,:) = kineticEnergyCell(iLevel,:)
            workArray(13,:) = relativeVorticityCell(iLevel,:)
            workArray(14,:) = divergence(iLevel,:)
            workArray(15,:) = relativeVorticityCell(iLevel,:)*relativeVorticityCell(iLevel,:)

            call compute_statistics(nDefinedDataFields, nCellsSolve, workArray, workMask, workMin, workMax, workSum)

            ! store data in buffer in order to allow only three dmpar calls
            do iDataField=1,nDefinedDataFields
              kBuffer = kBuffer+1
              workBufferSum(kBuffer) = workBufferSum(kBuffer) + workSum(iDataField)
              workBufferMin(kBuffer) = min(workBufferMin(kBuffer), workMin(iDataField))
              workBufferMax(kBuffer) = max(workBufferMax(kBuffer), workMax(iDataField))
            enddo

         enddo ! iLevel

         end do ! iRegion
         kBuffer = 0

         ! deallocate scratch fields
         call mpas_deallocate_scratch_field(workArrayField, .true.)
         call mpas_deallocate_scratch_field(workMaskField, .true.)
         call mpas_deallocate_scratch_field(workMinField, .true.)
         call mpas_deallocate_scratch_field(workMaxField, .true.)
         call mpas_deallocate_scratch_field(workSumField, .true.)

         block => block % next
      end do

      ! communication
      call mpas_dmpar_sum_real_array(dminfo, kBufferLength, workBufferSum, workBufferSumReduced )
      call mpas_dmpar_min_real_array(dminfo, kBufferLength, workBufferMin, workBufferMinReduced )
      call mpas_dmpar_max_real_array(dminfo, kBufferLength, workBufferMax, workBufferMaxReduced )

      ! unpack the buffer into intent(out) of this analysis member
      kBuffer=0
      do iRegion=1,nOceanRegionsTmp
        do iLevel=1,nVertLevels
           do iDataField=1,nDefinedDataFields
              kBuffer = kBuffer+1
              avgValueWithinOceanLayerRegion(iDataField,iLevel,iRegion)=workBufferSumReduced(kBuffer)
              minValueWithinOceanLayerRegion(iDataField,iLevel,iRegion)=workBufferMinReduced(kBuffer)
              maxValueWithinOceanLayerRegion(iDataField,iLevel,iRegion)=workBufferMaxReduced(kBuffer)
           enddo
        enddo
      enddo

      ! compute vertical sum before layer-by-layer normalization
      minValueWithinOceanVolumeRegion = 0.0_RKIND
      maxValueWithinOceanVolumeRegion = 0.0_RKIND
      avgValueWithinOceanVolumeRegion = 0.0_RKIND
      do iRegion=1,nOceanRegionsTmp
        do iDataField=1,nDefinedDataFields
          do iLevel=1,nVertLevels
             avgValueWithinOceanVolumeRegion(iDataField, iRegion) = avgValueWithinOceanVolumeRegion(iDataField, iRegion) + avgValueWithinOceanLayerRegion(iDataField,iLevel,iRegion)
          enddo
        enddo
        do iDataField=4,nDefinedDataFields
          avgValueWithinOceanVolumeRegion(iDataField, iRegion) = avgValueWithinOceanVolumeRegion(iDataField, iRegion) / max(avgValueWithinOceanVolumeRegion(3,iRegion),1.0e-8_RKIND)
        enddo
        ! normalize total region volume by total volume cell area
        avgValueWithinOceanVolumeRegion(3,iRegion) = avgValueWithinOceanVolumeRegion(3,iRegion) / max(avgValueWithinOceanVolumeRegion(2,iRegion),1.0e-8_RKIND)
        ! normalize total volume cell area by total number of cells
        avgValueWithinOceanVolumeRegion(2,iRegion) = avgValueWithinOceanVolumeRegion(2,iRegion) / max(avgValueWithinOceanVolumeRegion(1,iRegion),1.0e-8_RKIND)
      enddo

      ! find min/max with region volume
      do iRegion=1,nOceanRegionsTmp
        do iDataField=1,nDefinedDataFields
           minValueWithinOceanVolumeRegion(iDataField, iRegion) = minval(minValueWithinOceanLayerRegion(iDataField,:,iRegion))
           maxValueWithinOceanVolumeRegion(iDataField, iRegion) = maxval(minValueWithinOceanLayerRegion(iDataField,:,iRegion))
        enddo
      enddo

      ! normalize averages layer-by-layer
      do iRegion=1,nOceanRegionsTmp
      do iLevel=1,nVertLevels
         ! normalize all field by total volume in each layer
         do iDataField=4,nDefinedDataFields
            avgValueWithinOceanLayerRegion(iDataField,iLevel,iRegion) = avgValueWithinOceanLayerRegion(iDataField,iLevel,iRegion) / max(avgValueWithinOceanLayerRegion(3,iLevel,iRegion),1.0e-8_RKIND)
         enddo
         ! normalize total layer volume by layer area
         avgValueWithinOceanLayerRegion(3,iLevel,iRegion) = avgValueWithinOceanLayerRegion(3,iLevel,iRegion) / max(avgValueWithinOceanLayerRegion(2,iLevel,iRegion),1.0e-8_RKIND)
         ! normalize total layer area by number of cells in region
         avgValueWithinOceanLayerRegion(2,iLevel,iRegion) = avgValueWithinOceanLayerRegion(2,iLevel,iRegion) / max(avgValueWithinOceanLayerRegion(1,iLevel,iRegion),1.0e-8_RKIND)
      enddo
      enddo

      ! deallocate buffers
      deallocate(workBufferSumReduced)
      deallocate(workBufferMinReduced)
      deallocate(workBufferMaxReduced)

   contains

   subroutine compute_mask(iLevel, maxLevelCell, nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to
   ! allow for min/max/avg to represent specific regions of the ocean domain
   !
   ! NOTE: computes_mask is temporary. workMask should be intent(in) to this entire module !
   !
   integer, intent(in) :: iLevel, nCells, nCellsSolve, iRegion
   integer, intent(in), dimension(:) :: maxLevelCell
   real(kind=RKIND), dimension(:), intent(in) :: lonCell, latCell
   real(kind=RKIND), dimension(:), intent(out) :: workMask
   integer :: iCell
   real(kind=RKIND) :: dtr

   dtr = 4.0_RKIND*atan(1.0_RKIND) / 180.0_RKIND 
   workMask(:) = 0.0_RKIND
   do iCell=1,nCellsSolve
      if(iLevel.le.maxLevelCell(iCell)) workMask(iCell) = 1.0_RKIND
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

   end subroutine compute_mask


   subroutine compute_statistics(nDefinedDataFields, nCellsSolve, workArray, workMask, workMin, workMax, workSum)
   ! this subroutines does the actual summing, min, max, masking ect
   ! this hides the messy code from the high-level subroutine

   integer, intent(in) :: nDefinedDataFields, nCellsSolve
   real(kind=RKIND), dimension(:,:), intent(in) :: workArray
   real(kind=RKIND), dimension(:), intent(in)   :: workMask
   real(kind=RKIND), dimension(:), intent(out)  :: workMin, workMax, workSum
   integer :: iCell, iDataField
   real(kind=RKIND) :: cellMask, cellArea, cellVolume

   workSum = 0.0
   do iCell=1,nCellsSolve
    cellMask   = workMask(iCell)                                           ! mask
    cellArea   = cellMask * workArray(2,iCell)                             ! area
    cellVolume = cellArea * workArray(3,iCell)                             ! volume
    workSum(1) = workSum(1) + cellMask                                     ! 0/1 mask sum
    workSum(2) = workSum(2) + cellArea                                     ! area sum
    workSum(3) = workSum(3) + cellVolume                                   ! volume sum
    do iDataField=4,nDefinedDataFields
      workSum(iDataField) = workSum(iDataField) + cellVolume*workArray(iDataField,iCell)  ! volume-weighted sum
    enddo
   enddo

   do iDataField=1,nDefinedDataFields
      workMin(iDataField) = minval(workArray(iDataField,1:nCellsSolve),workMask(1:nCellsSolve)>0.5_RKIND)
      workMax(iDataField) = maxval(workArray(iDataField,1:nCellsSolve),workMask(1:nCellsSolve)>0.5_RKIND)
   enddo

   end subroutine compute_statistics

   end subroutine ocn_compute_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_restart_layer_volume_weighted_averages
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_layer_volume_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_restart_layer_volume_weighted_averages!}}}

!***********************************************************************
!
!  routine ocn_finalize_layer_volume_weighted_averages
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    April 24, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_layer_volume_weighted_averages(domain, err)!{{{

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

   end subroutine ocn_finalize_layer_volume_weighted_averages!}}}

end module ocn_layer_volume_weighted_averages

! vim: foldmethod=marker
