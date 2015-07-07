! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_water_mass_census
!
!> \brief MPAS ocean analysis member: water mass census
!> \author Todd Ringler
!> \date   May 12, 2015
!> \details
!>  MPAS ocean analysis member: water_mass_census
!>  This analysis member sorts the ocean water volume based on it
!>  temperature and salinity.
!
!-----------------------------------------------------------------------

module ocn_water_mass_census

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

   public :: ocn_init_water_mass_census, &
             ocn_compute_water_mass_census, &
             ocn_restart_water_mass_census, &
             ocn_finalize_water_mass_census

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_water_mass_census
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    May 10, 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_water_mass_census(domain, err)!{{{

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

   end subroutine ocn_init_water_mass_census!}}}

!***********************************************************************
!
!  routine ocn_compute_water_mass_census
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    May 10, 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_water_mass_census(domain, timeLevel, err)!{{{

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
      type (mpas_pool_type), pointer :: waterMassCensusAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: diagnosticsPool

      real (kind=RKIND), dimension(:,:,:), pointer :: waterMassFractionalDistribution
      real (kind=RKIND), dimension(:,:,:), pointer :: potentialDensityOfTSDiagram
      real (kind=RKIND), dimension(:,:,:), pointer :: zPositionOfTSDiagram
      real (kind=RKIND), dimension(:,:),   pointer :: waterMassCensusTemperatureValues
      real (kind=RKIND), dimension(:,:),   pointer :: waterMassCensusSalinityValues

      ! pointers to data in pools required for T/S water mass census
      real (kind=RKIND), dimension(:,:),   pointer :: layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers
      real (kind=RKIND), dimension(:,:),   pointer :: potentialDensity
      real (kind=RKIND), dimension(:,:),   pointer :: zMid

      ! pointers to data in mesh pool
      ! (note: nOceanRegionsTmpCensus, lonCell, latCell to be removed when region mask is intent(in))
      integer, pointer :: nCells, nCellsSolve, nOceanRegionsTmpCensus
      integer, pointer :: index_temperature, index_salinity
      integer, pointer :: nTemperatureBins, nSalinityBins
      integer, pointer :: nTemperatureBinsP1, nSalinityBinsP1
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer ::  areaCell, lonCell, latCell

      ! local variables
      integer :: iCell, iRegion, iLevel, iTracer, iTemperatureBin, iSalinityBin, err_tmp
      real (kind=RKIND), pointer :: minTemperature, maxTemperature
      real (kind=RKIND), pointer :: minSalinity, maxSalinity
      real (kind=RKIND) :: deltaTemperature, deltaSalinity, temperature, salinity, density, zPosition, volume

      ! package flag
      logical, pointer :: waterMassCensusAMPKGActive

      ! buffers data for message passaging
      integer :: kBuffer, kBufferLength
      real (kind=RKIND), dimension(:), allocatable :: workBufferSum, workBufferSumReduced
      ! (note: regionMask will (soon) be intent(in) and workMask can then be elimated)
      real (kind=RKIND), dimension(:), allocatable :: workMask
      real (kind=RKIND), dimension(:,:), allocatable :: regionMask

      ! assume no error
      err = 0

      ! set highest level pointer
      dminfo = domain % dminfo

      ! find the number of regions, number of data fields and number of vertical levels
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nOceanRegionsTmpCensus', nOceanRegionsTmpCensus)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTemperatureBins', nTemperatureBins)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nSalinityBins', nSalinityBins)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nTemperatureBinsP1', nTemperatureBinsP1)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nSalinityBinsP1', nSalinityBinsP1)

      ! allocate buffer for message passing
      kBuffer=0
      kBufferLength=3*nOceanRegionsTmpCensus*nTemperatureBins*nSalinityBins
      allocate(workBufferSum(kBufferLength))
      allocate(workBufferSumReduced(kBufferLength))
      workBufferSum=0.0_RKIND
      workBufferSumReduced=0.0_RKIND

      ! all code below will go away when regionMask is intent(in)
      ! allocate region mask and fill array
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nCells', nCells)
      allocate(workMask(nCells))
      allocate(regionMask(nOceanRegionsTmpCensus,nCells))
      block => domain % blocklist
      do while (associated(block))
         ! get pointers to pools
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         do iRegion=1,nOceanRegionsTmpCensus
           call compute_mask(maxLevelCell, nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)
           regionMask(iRegion,:) = workMask(:)
         enddo
         block => block % next
      enddo
      ! all code above will go away when regionMask is intent(in)

      ! get pointers to analysis member arrays
      call mpas_pool_get_subpool(domain % blocklist % structs, 'waterMassCensusAM', waterMassCensusAMPool)
      call mpas_pool_get_array(waterMassCensusAMPool, 'waterMassCensusTemperatureValues', waterMassCensusTemperatureValues)
      call mpas_pool_get_array(waterMassCensusAMPool, 'waterMassCensusSalinityValues', waterMassCensusSalinityValues)
      call mpas_pool_get_array(waterMassCensusAMPool, 'waterMassFractionalDistribution', waterMassFractionalDistribution)
      call mpas_pool_get_array(waterMassCensusAMPool, 'potentialDensityOfTSDiagram', potentialDensityOfTSDiagram)
      call mpas_pool_get_array(waterMassCensusAMPool, 'zPositionOfTSDiagram', zPositionOfTSDiagram)

      ! get run-time configure variables
      call mpas_pool_get_config(domain % configs, 'config_AM_waterMassCensus_minTemperature', minTemperature)
      call mpas_pool_get_config(domain % configs, 'config_AM_waterMassCensus_maxTemperature', maxTemperature)
      call mpas_pool_get_config(domain % configs, 'config_AM_waterMassCensus_minSalinity', minSalinity)
      call mpas_pool_get_config(domain % configs, 'config_AM_waterMassCensus_maxSalinity', maxSalinity)

      do iRegion=1,nOceanRegionsTmpCensus
        ! compute temperature and salinity domains
        ! (note: the ability to have different t/s domains for different regions is not yet built out)
        deltaTemperature = (maxTemperature-minTemperature)/nTemperatureBins
        do iTemperatureBin=1,nTemperatureBinsP1
          waterMassCensusTemperatureValues(iTemperatureBin,iRegion) = minTemperature + deltaTemperature*(iTemperatureBin-1)
        enddo
        deltaSalinity = (maxSalinity-minSalinity)/nSalinityBins
        do iSalinityBin=1,nSalinityBinsP1
          waterMassCensusSalinityValues(iSalinityBin,iRegion) = minSalinity + deltaSalinity*(iSalinityBin-1)
        enddo
      enddo ! iRegion

      ! initialize intent(out) of this analysis member
      waterMassFractionalDistribution(:,:,:)=0.0_RKIND
      potentialDensityOfTSDiagram(:,:,:)=0.0_RKIND
      zPositionOfTSDiagram(:,:,:)=0.0_RKIND

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))
         ! get pointers to pools
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         ! get pointers to mesh
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nOceanRegionsTmpCensus', nOceanRegionsTmpCensus)
         call mpas_pool_get_dimension(statePool, 'index_temperature', index_temperature)
         call mpas_pool_get_dimension(statePool, 'index_salinity', index_salinity)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         ! get pointers to data needed for analysis
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
         call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

         ! loop over and bin all data
         do iCell=1,nCellsSolve
         do iLevel=1,maxLevelCell(iCell)

           ! make copies of data for convienence
           temperature = tracers(index_temperature,iLevel,iCell)
           salinity = tracers(index_salinity,iLevel,iCell)
           density = potentialDensity(iLevel,iCell)
           zPosition = zMid(iLevel,iCell)
           volume = layerThickness(iLevel,iCell) * areaCell(iCell)

           ! find temperature bin, cycle if bin is out of range
           iTemperatureBin = int((temperature-minTemperature)/deltaTemperature) + 1
           if (iTemperatureBin < 1) cycle
           if (iTemperatureBin > nTemperatureBins) cycle

           ! find salinity bin, cycle if bin is out of range
           iSalinityBin = int((salinity-minSalinity)/deltaSalinity) + 1
           if (iSalinityBin < 1) cycle
           if (iSalinityBin > nSalinityBins) cycle

           do iRegion=1,nOceanRegionsTmpCensus
             ! add volume into water mass census array for each region
             waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion) = &
                       waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion)  &
                       + volume * regionMask(iRegion,iCell)
             potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = &
                       potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion)  &
                       + density * volume * regionMask(iRegion,iCell)
             zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = &
                       zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion)  &
                       + zPosition * volume * regionMask(iRegion,iCell)
           enddo

         enddo   ! iLevel
         enddo   ! iCell

         block => block % next
      end do   ! block loop

      ! store data in buffer in order to allow only one dmpar calls
      kBuffer=0
      do iTemperatureBin=1,nTemperatureBins
      do iSalinityBin=1,nSalinityBins
      do iRegion=1,nOceanRegionsTmpCensus
        kBuffer = kBuffer+1
        workBufferSum(kBuffer) = waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion)
        kBuffer = kBuffer+1
        workBufferSum(kBuffer) = potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion)
        kBuffer = kBuffer+1
        workBufferSum(kBuffer) = zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion)
      enddo
      enddo
      enddo

      ! communication
      call mpas_dmpar_sum_real_array(dminfo, kBufferLength, workBufferSum, workBufferSumReduced )

      ! unpack the buffer into intent(out) of this analysis member
      kBuffer=0
      do iTemperatureBin=1,nTemperatureBins
      do iSalinityBin=1,nSalinityBins
      do iRegion=1,nOceanRegionsTmpCensus
        kBuffer = kBuffer+1
        waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion) = workBufferSumReduced(kBuffer)
        kBuffer = kBuffer+1
        potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = workBufferSumReduced(kBuffer)
        kBuffer = kBuffer+1
        zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = workBufferSumReduced(kBuffer)
      enddo
      enddo
      enddo

      ! normalize potentialDensityOfTSDiagram by volume in each T,S bin
      do iTemperatureBin=1,nTemperatureBins
      do iSalinityBin=1,nSalinityBins
      do iRegion=1,nOceanRegionsTmpCensus
        potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = &
          potentialDensityOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) / &
          max(waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion), 1.0e-8_RKIND)
        zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) = &
          zPositionOfTSDiagram(iTemperatureBin,iSalinityBin,iRegion) / &
          max(waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion), 1.0e-8_RKIND)
      enddo
      enddo
      enddo

      ! use workBufferSum as workspace to find total volume for each region
      workBufferSum = 0.0_RKIND
      do iTemperatureBin=1,nTemperatureBins
      do iSalinityBin=1,nSalinityBins
      do iRegion=1,nOceanRegionsTmpCensus
        workBufferSum(iRegion) = workBufferSum(iRegion) + waterMassFractionalDistribution(iTemperatureBin,iSalinityBin,iRegion)
      enddo
      enddo
      enddo

      ! use this sum to convert waterMassFractionalDistribution from total volume to fractional volume
      do iRegion=1,nOceanRegionsTmpCensus
        waterMassFractionalDistribution(:,:,iRegion) = waterMassFractionalDistribution(:,:,iRegion) / max(workBufferSum(iRegion),1.0e-8_RKIND)
      enddo

      ! deallocate buffers
      deallocate(workBufferSum)
      deallocate(workBufferSumReduced)
      deallocate(regionMask)
      deallocate(workMask)

   contains

   subroutine compute_mask(maxLevelCell, nCells, nCellsSolve, iRegion, lonCell, latCell, workMask)
   ! this subroutines produces a 0/1 mask that is multiplied with workArray to
   ! allow for min/max/avg to represent specific regions of the ocean domain
   !
   ! NOTE: computes_mask is temporary. workMask should be intent(in) to this entire module !
   !
   integer, intent(in) :: nCells, nCellsSolve, iRegion
   integer, intent(in), dimension(:) :: maxLevelCell
   real(kind=RKIND), dimension(:), intent(in) :: lonCell, latCell
   real(kind=RKIND), dimension(:), intent(out) :: workMask
   integer :: iCell
   real(kind=RKIND) :: dtr

   dtr = 4.0*atan(1.0) / 180.0_RKIND
   workMask(:) = 0.0_RKIND
   do iCell=1,nCellsSolve
      workMask(iCell) = 1.0_RKIND
   enddo

   if (iRegion.eq.1) then
      ! Arctic
      do iCell=1,nCellsSolve
        if(latCell(iCell).lt. 60.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Arctic ', sum(workMask)
   elseif (iRegion.eq.2) then
      ! Equatorial
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt. 15.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt.-15.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Equatorial ', sum(workMask)
   elseif (iRegion.eq.3) then
      ! Southern Ocean
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.-50.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Southern Ocean ', sum(workMask)
   elseif (iRegion.eq.4) then
      ! Nino 3
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.210.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.270.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Nino 3 ', sum(workMask)
   elseif (iRegion.eq.5) then
      ! Nino 4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.160.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.210.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Nino 4 ', sum(workMask)
   elseif (iRegion.eq.6) then
      ! Nino 3.4
      do iCell=1,nCellsSolve
        if(latCell(iCell).gt.  5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(latCell(iCell).lt. -5.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).lt.190.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
        if(lonCell(iCell).gt.240.0_RKIND*dtr) workMask(iCell) = 0.0_RKIND
      enddo
      write(6,*) ' Nino 3.4 ', sum(workMask)
   else
      ! global (do nothing!)
      write(6,*) ' Global ', sum(workMask)
   endif

   end subroutine compute_mask

   end subroutine ocn_compute_water_mass_census!}}}

!***********************************************************************
!
!  routine ocn_restart_water_mass_census
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    May 10, 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_water_mass_census(domain, err)!{{{

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

   end subroutine ocn_restart_water_mass_census!}}}

!***********************************************************************
!
!  routine ocn_finalize_water_mass_census
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Todd Ringler
!> \date    May 10, 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_water_mass_census(domain, err)!{{{

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

   end subroutine ocn_finalize_water_mass_census!}}}

end module ocn_water_mass_census

! vim: foldmethod=marker
