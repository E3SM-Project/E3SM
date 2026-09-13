! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_li_regional_stats
!
!> \brief MPAS land ice analysis mode member: mpas_li_regional_stats
!> \author Stephen Price
!> \date   4-29-2016
!> \details
!>
!>
!-----------------------------------------------------------------------
module li_regional_stats

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_global_sum_mod, only: mpas_global_sum_nfld

   use li_mask
   use li_constants

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

   public :: li_init_regional_stats, &
             li_compute_regional_stats, &
             li_restart_regional_stats, &
             li_finalize_regional_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_init_regional_stats
!
!> \brief   Initialize MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_init_regional_stats(domain, memberName, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character (len=*), intent(in) :: memberName

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

   end subroutine li_init_regional_stats

!***********************************************************************
!
!  routine li_compute_regional_stats
!
!> \brief   Compute MPAS-Land Ice analysis member
!> \author  S. Price
!> \date   4-29-2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_compute_regional_stats(domain, memberName, timeLevel, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel
      character (len=*), intent(in) :: memberName

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
      type (mpas_pool_type), pointer :: regionalStatsAMPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: regionalStatsAM
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: regionsPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: hydroPool

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), pointer ::  config_ice_density
      real (kind=RKIND), pointer ::  deltat
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  dvEdge
      real (kind=RKIND), dimension(:), pointer ::  dcEdge
      real (kind=RKIND), dimension(:), pointer ::  thickness
      real (kind=RKIND), dimension(:), pointer ::  bedTopography
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBalApplied
      real (kind=RKIND), dimension(:), pointer ::  groundedSfcMassBalApplied
      real (kind=RKIND), dimension(:), pointer ::  basalMassBalApplied
      real (kind=RKIND), dimension(:), pointer ::  groundedBasalMassBalApplied
      real (kind=RKIND), dimension(:), pointer ::  floatingBasalMassBalApplied
      real (kind=RKIND), dimension(:), pointer ::  calvingThickness
      real (kind=RKIND), dimension(:), pointer ::  faceMeltingThickness
      real (kind=RKIND), dimension(:), pointer ::  surfaceSpeed
      real (kind=RKIND), dimension(:), pointer ::  basalSpeed
      real (kind=RKIND), dimension(:), pointer ::  fluxAcrossGroundingLine
      real (kind=RKIND), dimension(:), pointer ::  groundedToFloatingThickness
      real (kind=RKIND), dimension(:,:), pointer ::  normalVelocity

      real (kind=RKIND), dimension(:), pointer ::  waterThickness
      real (kind=RKIND), dimension(:), pointer ::  basalMeltInput
      real (kind=RKIND), dimension(:), pointer ::  externalWaterInput
      real (kind=RKIND), dimension(:), pointer ::  channelMelt
      real (kind=RKIND), dimension(:), pointer ::  waterFlux
      real (kind=RKIND), dimension(:), pointer ::  channelDischarge
      real (kind=RKIND), dimension(:), pointer ::  waterPressure

      ! config options needed
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoi ! config_ice_density
      real (kind=RKIND), pointer :: rhow ! config_ocean_density
      real (kind=RKIND), pointer :: bedBumpMax ! config_SGH_bed_roughness_max
      logical, pointer :: config_SGH

      integer, dimension(:,:), pointer :: regionCellMasks
      integer, dimension(:,:), pointer :: regionEdgeMasks
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: edgeMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: hydroTerrestrialMarginMask
      integer, pointer :: nRegions, nRegionGroups !, maxRegionsInGroup  !! maxRegionsInGroup not needed / used yet
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels
      integer :: k, iCell, iEdge
      integer :: upwindCell
      integer :: iRegion, iGroup, iVar

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), dimension(:), pointer ::  regionalIceArea, regionalIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalVolumeAboveFloatation
      real (kind=RKIND), dimension(:), pointer ::  regionalGroundedIceArea, regionalGroundedIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalFloatingIceArea, regionalFloatingIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalIceThicknessMin, regionalIceThicknessMax, regionalIceThicknessMean
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSfcMassBal, regionalSumBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundedSfcMassBal, regionalSumFloatingSfcMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundedBasalMassBal, regionalSumFloatingBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumCalvingFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalSumFaceMeltingFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundingLineFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundingLineMigrationFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgNetAccumulation
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgGroundedBasalMelt
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgSubshelfMelt
      real (kind=RKIND), dimension(:), pointer ::  regionalSurfaceSpeedMax
      real (kind=RKIND), dimension(:), pointer ::  regionalBasalSpeedMax
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSubglacialWaterVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalSumBasalMeltInput
      real (kind=RKIND), dimension(:), pointer ::  regionalSumExternalWaterInput
      real (kind=RKIND), dimension(:), pointer ::  regionalSumChannelMelt
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSubglacialLakeVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSubglacialLakeArea
      real (kind=RKIND), dimension(:), pointer ::  regionalSumDistWaterFluxMarineMargin
      real (kind=RKIND), dimension(:), pointer ::  regionalSumDistWaterFluxTerrestrialMargin
      real (kind=RKIND), dimension(:), pointer ::  regionalSumChnlWaterFluxMarineMargin
      real (kind=RKIND), dimension(:), pointer ::  regionalSumChnlWaterFluxTerrestrialMargin
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgFlotationFraction

      ! storage for min/max reductions over blocks (kept as block-local scalars per region,
      ! since min/max are exact/associative operations and do not need a reproducible sum)
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionThickMin, blockRegionThickMax
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionMaxSurfaceSpeed
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionMaxBasalSpeed

      ! Raw (region-independent) per-cell and per-edge contributions to the sum-type
      ! reductions, computed once per cell/edge (outside the region loop). For each region,
      ! these are multiplied by that region's cell/edge mask and reduced with the reproducible
      ! mpas_global_sum_nfld routine, instead of accumulating a block-local scalar per region
      ! and reducing it with a plain MPI_Allreduce (mpas_dmpar_sum_real_array), which is not
      ! decomposition-independent due to floating-point non-associativity.
      ! NOTE: MALI supports only a single block per MPI task, so these arrays are sized
      ! and filled for the single block on this task.
      integer, parameter :: nCellSumVars = 22, nEdgeSumVars = 6
      real (kind=RKIND), dimension(:,:), allocatable :: cellRawContribution, edgeRawContribution
      real (kind=RKIND), dimension(:,:), allocatable :: cellSumContribution, edgeSumContribution
      real (kind=RKIND), dimension(nCellSumVars) :: cellGlobalSums
      real (kind=RKIND), dimension(nEdgeSumVars) :: edgeGlobalSums
      integer, dimension(:), allocatable :: edgeUpwindCell  ! upwind cell for each GL edge, used to assign GL flux to a region

      ! column indices into cellRawContribution/cellSumContribution/cellGlobalSums
      integer, parameter :: cIceArea = 1, cIceVolume = 2, cGroundedIceArea = 3, cGroundedIceVolume = 4, &
                             cFloatingIceArea = 5, cFloatingIceVolume = 6, cSfcMassBal = 7, cGroundedSfcMassBal = 8, &
                             cFloatingSfcMassBal = 9, cBasalMassBal = 10, cGroundedBasalMassBal = 11, &
                             cFloatingBasalMassBal = 12, cCalvingFlux = 13, cFaceMeltingFlux = 14, cVAF = 15, &
                             cGLMigrationFlux = 16, cSubglacialWaterVolume = 17, cBasalMeltInput = 18, &
                             cExternalWaterInput = 19, cLakeVolume = 20, cLakeArea = 21, cFlotationFraction = 22

      ! column indices into edgeRawContribution/edgeSumContribution/edgeGlobalSums
      integer, parameter :: eGLflux = 1, eChannelMelt = 2, eGLMeltFlux = 3, eTerrestrialMeltFlux = 4, &
                             eChannelGLMeltFlux = 5, eChannelTerrestrialMeltFlux = 6
      ! local variable needed calculating numerator of regionalAvgFlotationFraction
      real (kind=RKIND), dimension(:), allocatable:: regionalSumFlotationFraction

      ! local variables
      real (kind=RKIND) :: fluxSign

      ! variables for processing stats
      integer, parameter :: kMaxVariables = 43 ! Increase if number of stats increase
      integer :: nVars
      real (kind=RKIND), dimension(kMaxVariables) :: reductions, mins, maxes

      err = 0

      dminfo = domain % dminfo

      ! Get needed configs
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhow)
      call mpas_pool_get_config(liConfigs, 'config_SGH_bed_roughness_max', bedBumpMax)
      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

         !! get dimensions assocated with region masks
         call mpas_pool_get_dimension(block % dimensions, 'nRegions', nRegions)
         call mpas_pool_get_dimension(block % dimensions, 'nRegionGroups', nRegionGroups)

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'regions', regionsPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)
!        call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
         call mpas_pool_get_array(geometryPool, 'sfcMassBalApplied', sfcMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'groundedSfcMassBalApplied', groundedSfcMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'basalMassBalApplied', basalMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'groundedBasalMassBalApplied', groundedBasalMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBalApplied', floatingBasalMassBalApplied)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'faceMeltingThickness', faceMeltingThickness)
         call mpas_pool_get_array(geometryPool, 'groundedToFloatingThickness', groundedToFloatingThickness)
         call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)
         call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
         call mpas_pool_get_array(velocityPool, 'fluxAcrossGroundingLine', fluxAcrossGroundingLine)
         call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
         if (config_SGH) then
            call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
            call mpas_pool_get_array(hydroPool, 'waterThickness', waterThickness)
            call mpas_pool_get_array(hydroPool, 'basalMeltInput', basalMeltInput)
            call mpas_pool_get_array(hydroPool, 'externalWaterInput', externalWaterInput)
            call mpas_pool_get_array(hydroPool, 'channelMelt', channelMelt)
            call mpas_pool_get_array(hydroPool, 'hydroMarineMarginMask', hydroMarineMarginMask)
            call mpas_pool_get_array(hydroPool, 'hydroTerrestrialMarginMask', hydroTerrestrialMarginMask)
            call mpas_pool_get_array(hydroPool, 'waterFlux', waterFlux)
            call mpas_pool_get_array(hydroPool, 'channelDischarge', channelDischarge)
            call mpas_pool_get_array(hydroPool, 'waterPressure', waterPressure)
         endif 

         ! get region cell masks from regionMasks.nc input file
         call mpas_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)
         if (config_SGH) then
             ! regionEdgeMasks is required if calculating regional stats with SGH, 
             ! in the fututre will probably be required for standalone regional stats
             ! but this is used a stopgap to forestall backwards compatibility problems
            call mpas_pool_get_array(regionsPool, 'regionEdgeMasks', regionEdgeMasks)
         endif
         ! allocate & initialize sums over blocks to 0
         allocate(blockRegionThickMin(nRegions)); allocate(blockRegionThickMax(nRegions))
         allocate(blockRegionMaxBasalSpeed(nRegions)); allocate(blockRegionMaxSurfaceSpeed(nRegions))
         if (config_SGH) then
            allocate(regionalSumFlotationFraction(nRegions))
         endif

         blockRegionThickMin = 10000.0_RKIND; blockRegionThickMax = 0.0_RKIND
         blockRegionMaxBasalSpeed = 0.0_RKIND; blockRegionMaxSurfaceSpeed = 0.0_RKIND

         if (config_SGH) then
            regionalSumFlotationFraction = 0.0_RKIND
         endif

         ! Compute raw (region-independent) per-cell/per-edge contributions once; these will be
         ! multiplied by each region's mask and reduced (reproducibly) after the cell/edge loops.
         if (.not. allocated(cellRawContribution)) allocate(cellRawContribution(nCellsSolve, nCellSumVars))
         if (.not. allocated(edgeRawContribution)) allocate(edgeRawContribution(nEdgesSolve, nEdgeSumVars))
         if (.not. allocated(cellSumContribution)) allocate(cellSumContribution(nCellsSolve, nCellSumVars))
         if (.not. allocated(edgeSumContribution)) allocate(edgeSumContribution(nEdgesSolve, nEdgeSumVars))
         if (.not. allocated(edgeUpwindCell)) allocate(edgeUpwindCell(nEdgesSolve))
         cellRawContribution(:,:) = 0.0_RKIND
         edgeRawContribution(:,:) = 0.0_RKIND
         edgeUpwindCell(:) = 0

         do iCell = 1,nCellsSolve      ! loop over cells

            ! regional sums of ice area over cells with ice (m^2)
            cellRawContribution(iCell, cIceArea) = areaCell(iCell) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND)

            ! regional sums of ice volume over cells with ice (m^3)
            cellRawContribution(iCell, cIceVolume) = areaCell(iCell) * thickness(iCell) &
               * real(li_mask_is_ice_int(cellMask(iCell)),RKIND)

            ! regional sums of ice volume above floatation over cells with ice (m^3)
            cellRawContribution(iCell, cVAF) = real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
               * areaCell(iCell) * ( thickness(iCell) + (rhow / rhoi) &
               * min(0.0_RKIND, (bedTopography(iCell) - config_sea_level)) )

            ! regional sums of grounded ice area over cells with ice (m^2)
            cellRawContribution(iCell, cGroundedIceArea) = areaCell(iCell) &
                * real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND)

            ! regional sums of grounded ice volume over cells with ice (m^3)
            cellRawContribution(iCell, cGroundedIceVolume) = areaCell(iCell) &
              * thickness(iCell) * real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND)

            ! regional sums of floating ice area over cells with ice (m^2)
            cellRawContribution(iCell, cFloatingIceArea) = areaCell(iCell) &
                * real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND)

            ! regional sums of floating ice volume over cells with ice (m^3)
            cellRawContribution(iCell, cFloatingIceVolume) = areaCell(iCell) &
              * thickness(iCell) * real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND)

            ! regional sum of sfc mass balance (kg yr^{-1})
            cellRawContribution(iCell, cSfcMassBal) = areaCell(iCell) * sfcMassBalApplied(iCell) * scyr

            ! regional sum of sfc mass balance on grounded ice (kg yr^{-1})
            cellRawContribution(iCell, cGroundedSfcMassBal) = areaCell(iCell) * groundedSfcMassBalApplied(iCell) * scyr

            ! regional sum of sfc mass balance on floatingice (kg yr^{-1})
            cellRawContribution(iCell, cFloatingSfcMassBal) = &
                areaCell(iCell) * (sfcMassBalApplied(iCell) - groundedSfcMassBalApplied(iCell)) * scyr

            ! regional sum of basal mass balance (kg yr^{-1})
            cellRawContribution(iCell, cBasalMassBal) = areaCell(iCell) * basalMassBalApplied(iCell) * scyr

            ! regional sum of floating basal mass balance (kg yr^{-1})
            cellRawContribution(iCell, cFloatingBasalMassBal) = areaCell(iCell) * floatingBasalMassBalApplied(iCell) * scyr

            ! regional sum of grounded basal mass balance (kg yr^{-1})
            cellRawContribution(iCell, cGroundedBasalMassBal) = areaCell(iCell) * groundedBasalMassBalApplied(iCell) * scyr

            ! regional sum of mass lass due do calving (kg yr^{-1})
            cellRawContribution(iCell, cCalvingFlux) = &
              calvingThickness(iCell) * areaCell(iCell) * config_ice_density / ( deltat / scyr )

            ! regional sum of mass lass due do face-melting (kg yr^{-1})
            cellRawContribution(iCell, cFaceMeltingFlux) = &
              faceMeltingThickness(iCell) * areaCell(iCell) * config_ice_density / ( deltat / scyr )

            ! GL migration flux
            cellRawContribution(iCell, cGLMigrationFlux) = &
               groundedToFloatingThickness(iCell) * areaCell(iCell) * rhoi / (deltat / scyr)

            !! Subglacial Hydrology Calculations
            if (config_SGH) then

               ! Subglacial Water Volume
               cellRawContribution(iCell, cSubglacialWaterVolume) = waterThickness(iCell) * areaCell(iCell)

               if (li_mask_is_grounded_ice(cellMask(iCell))) then
                  ! Basal melt input
                  cellRawContribution(iCell, cBasalMeltInput) = basalMeltInput(iCell) * areaCell(iCell)
               endif

               ! External water input
               cellRawContribution(iCell, cExternalWaterInput) = &
                                                            real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) * &
                                                            externalWaterInput(iCell) * areaCell(iCell)

               ! Lake Volume
               if (waterThickness(iCell) > bedBumpMax) then
                  cellRawContribution(iCell, cLakeVolume) = (waterThickness(iCell) - bedBumpMax) * areaCell(iCell)
               endif

               ! Lake Area
               if (waterThickness(iCell) > bedBumpMax) then
                   cellRawContribution(iCell, cLakeArea) = areaCell(iCell)
               endif

               ! Area-weighted flotation fraction for grounded ice
               if (li_mask_is_grounded_ice(cellMask(iCell))) then
                   cellRawContribution(iCell, cFlotationFraction) = &
                     ( waterPressure(iCell) / rhoi / gravity / thickness(iCell) ) * areaCell(iCell)
               endif
            endif

           do iRegion = 1,nRegions     ! loop over regions, for min/max only (exact, region-restricted)

            ! max, min thickness values (m)
            if( ( thickness(iCell) * real(regionCellMasks(iRegion,iCell),RKIND) * &
               real(li_mask_is_ice_int(cellMask(iCell)),RKIND) ) > blockRegionThickMax(iRegion) ) then
               blockRegionThickMax(iRegion) = thickness(iCell)
            endif
            if( ( thickness(iCell) * real(regionCellMasks(iRegion,iCell),RKIND) * &
               real(li_mask_is_ice_int(cellMask(iCell)),RKIND) ) < blockRegionThickMin(iRegion) .and. &
               ( thickness(iCell) * real(regionCellMasks(iRegion,iCell),RKIND) &
               * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) ) > 0.0_RKIND )then
               blockRegionThickMin(iRegion) = thickness(iCell)
            endif

            ! max basal and sfc speeds (m/yr)
            if( ( surfaceSpeed(iCell) * real(regionCellMasks(iRegion,iCell),RKIND) * &
               real(li_mask_is_ice_int(cellMask(iCell)),RKIND) ) > blockRegionMaxSurfaceSpeed(iRegion) ) then
               blockRegionMaxSurfaceSpeed(iRegion) = surfaceSpeed(iCell)
            endif
            if( ( basalSpeed(iCell) * real(regionCellMasks(iRegion,iCell),RKIND) * &
               real(li_mask_is_ice_int(cellMask(iCell)),RKIND) ) > blockRegionMaxBasalSpeed(iRegion) ) then
               blockRegionMaxBasalSpeed(iRegion) = basalSpeed(iCell)
            endif

           end do ! end loop over regions
         end do ! end loop over cells


         do iEdge = 1,nEdgesSolve      ! loop over edges

            if (li_mask_is_grounding_line(edgeMask(iEdge))) then
               ! Determine upwind cell - we are choosing to assign edge quantities to the region of the upwind cell
               if (normalVelocity(1,iEdge) > 0.0_RKIND) then  ! using surface velo to determine upwind cell
                  ! upwind cell is cell 1 on this edge
                  upwindCell = cellsOnEdge(1,iEdge)
               else
                  ! upwind cell is cell 2 on this edge
                  upwindCell = cellsOnEdge(2,iEdge)
               endif

               ! record the raw GL flux and its upwind cell; the region assignment (based on
               ! regionCellMasks(iRegion, upwindCell)) is applied per-region below
               edgeUpwindCell(iEdge) = upwindCell
               edgeRawContribution(iEdge, eGLflux) = fluxAcrossGroundingLine(iEdge) * dvEdge(iEdge) * scyr * config_ice_density
                  ! convert from m^3/s to kg/yr
            end if ! if GL

            ! assign the SGH stats based on edge values using the regionEdgeMask 
            if (config_SGH) then

               ! Channel Melt
               edgeRawContribution(iEdge, eChannelMelt) = abs(channelMelt(iEdge) * dcEdge(iEdge))

               ! Meltwater Flux across the grounding line
               edgeRawContribution(iEdge, eGLMeltFlux) = &
                 abs(hydroMarineMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Flux across terrestrial margins
               edgeRawContribution(iEdge, eTerrestrialMeltFlux) = &
                 abs(hydroTerrestrialMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Discharge in channels across grounding line
               edgeRawContribution(iEdge, eChannelGLMeltFlux) = &
                 abs(hydroMarineMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)

               ! Meltwater discharge in channels across terrestrial margin
               edgeRawContribution(iEdge, eChannelTerrestrialMeltFlux) = &
                 abs( hydroTerrestrialMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)
            endif ! if SGH is on
         end do ! end loop over edges

         block => block % next

      end do    ! end loop over blocks

      ! --- Perform Reductions ---
      ! For each type of reduction (sum, min, max), set up an array so we only need
      ! a single reduction for each type, rather than a reduction for each variable
      ! since these global communications can be expensive on many processors.
      ! For now, the mapping of variable names to indices is handled manually,
      ! but if the number of variables grows, we could consider using a pool (a dictionary)
      ! or some other strategy for automating the relationships.
      ! Once the reduction is complete, stick the reduced value into the globalStats pool member.
      ! Even though some (most?) variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.

      ! ---------------------------
      ! compute sums (and means) over all procs
      ! ---------------------------

      do iRegion = 1,nRegions     ! loop over regions

        reductions = 0.0_RKIND

        ! Mask the region-independent raw contributions by this region's cell/edge mask,
        ! then perform a single reproducible reduction across all cell (and all edge) variables.
        ! (GL flux is masked by the region of its GL edge's upwind cell, not by regionEdgeMasks.)
        do iVar = 1, nCellSumVars
           cellSumContribution(1:nCellsSolve,iVar) = cellRawContribution(1:nCellsSolve,iVar) * &
              real(regionCellMasks(iRegion,1:nCellsSolve),RKIND)
        end do
        edgeSumContribution(:,:) = 0.0_RKIND
        do iEdge = 1, nEdgesSolve
           if (edgeUpwindCell(iEdge) > 0) then
              edgeSumContribution(iEdge,eGLflux) = edgeRawContribution(iEdge,eGLflux) * &
                 real(regionCellMasks(iRegion,edgeUpwindCell(iEdge)),RKIND)
           end if
        end do
        if (config_SGH) then
           do iVar = 1, nEdgeSumVars
              if (iVar /= eGLflux) then
                 edgeSumContribution(1:nEdgesSolve,iVar) = edgeRawContribution(1:nEdgesSolve,iVar) * &
                    real(regionEdgeMasks(iRegion,1:nEdgesSolve),RKIND)
              end if
           end do
        end if

        cellGlobalSums(:) = mpas_global_sum_nfld(cellSumContribution, dminfo % comm)
        edgeGlobalSums(:) = mpas_global_sum_nfld(edgeSumContribution, dminfo % comm)

        reductions(1) = cellGlobalSums(cIceArea)
        reductions(2) = cellGlobalSums(cIceVolume)
        reductions(3) = cellGlobalSums(cGroundedIceArea)
        reductions(4) = cellGlobalSums(cGroundedIceVolume)
        reductions(5) = cellGlobalSums(cFloatingIceArea)
        reductions(6) = cellGlobalSums(cFloatingIceVolume)
        reductions(7) = cellGlobalSums(cSfcMassBal)
        reductions(8) = cellGlobalSums(cGroundedSfcMassBal)
        reductions(9) = cellGlobalSums(cFloatingSfcMassBal)
        reductions(10) = cellGlobalSums(cBasalMassBal)
        reductions(11) = cellGlobalSums(cGroundedBasalMassBal)
        reductions(12) = cellGlobalSums(cFloatingBasalMassBal)
        reductions(13) = cellGlobalSums(cCalvingFlux)
        reductions(14) = cellGlobalSums(cFaceMeltingFlux)
        reductions(15) = cellGlobalSums(cVAF)
        reductions(16) = edgeGlobalSums(eGLflux)
        reductions(17) = cellGlobalSums(cGLMigrationFlux)
        if (config_SGH) then
           reductions(18) = cellGlobalSums(cSubglacialWaterVolume)
           reductions(19) = cellGlobalSums(cBasalMeltInput)
           reductions(20) = cellGlobalSums(cExternalWaterInput)
           reductions(21) = edgeGlobalSums(eChannelMelt)
           reductions(22) = cellGlobalSums(cLakeVolume)
           reductions(23) = cellGlobalSums(cLakeArea)
           reductions(24) = edgeGlobalSums(eGLMeltFlux)
           reductions(25) = edgeGlobalSums(eTerrestrialMeltFlux)
           reductions(26) = edgeGlobalSums(eChannelGLMeltFlux)
           reductions(27) = edgeGlobalSums(eChannelTerrestrialMeltFlux)
           reductions(28) = cellGlobalSums(cFlotationFraction)
        endif

        block => domain % blocklist
        do while (associated(block))

           call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

           ! get values from regional stats pool
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceArea', regionalIceArea)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceVolume', regionalIceVolume)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalVolumeAboveFloatation', regionalVolumeAboveFloatation)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalGroundedIceArea', regionalGroundedIceArea)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalGroundedIceVolume', regionalGroundedIceVolume)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalFloatingIceArea', regionalFloatingIceArea)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalFloatingIceVolume', regionalFloatingIceVolume)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMean', regionalIceThicknessMean)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumSfcMassBal', regionalSumSfcMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgNetAccumulation', regionalAvgNetAccumulation)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundedSfcMassBal', regionalSumGroundedSfcMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumFloatingSfcMassBal', regionalSumFloatingSfcMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumBasalMassBal', regionalSumBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundedBasalMassBal', regionalSumGroundedBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgGroundedBasalMelt', regionalAvgGroundedBasalMelt)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumFloatingBasalMassBal', regionalSumFloatingBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgSubshelfMelt', regionalAvgSubshelfMelt)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumCalvingFlux', regionalSumCalvingFlux)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumFaceMeltingFlux', regionalSumFaceMeltingFlux)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundingLineFlux', regionalSumGroundingLineFlux)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundingLineMigrationFlux', &
                   regionalSumGroundingLineMigrationFlux)
           if (config_SGH) then
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumSubglacialWaterVolume', regionalSumSubglacialWaterVolume)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumBasalMeltInput', regionalSumBasalMeltInput)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumExternalWaterInput', regionalSumExternalWaterInput)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumChannelMelt', regionalSumChannelMelt)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumSubglacialLakeVolume', regionalSumSubglacialLakeVolume)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumSubglacialLakeArea', regionalSumSubglacialLakeArea)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumDistWaterFluxMarineMargin',regionalSumDistWaterFluxMarineMargin)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumDistWaterFluxTerrestrialMargin', regionalSumDistWaterFluxTerrestrialMargin)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumChnlWaterFluxMarineMargin',regionalSumChnlWaterFluxMarineMargin)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumChnlWaterFluxTerrestrialMargin', regionalSumChnlWaterFluxTerrestrialMargin)
              call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgFlotationFraction', regionalAvgFlotationFraction)
           endif

           regionalIceArea(iRegion) = reductions(1)
           regionalIceVolume(iRegion) = reductions(2)
           regionalGroundedIceArea(iRegion) = reductions(3)
           regionalGroundedIceVolume(iRegion) = reductions(4)
           regionalFloatingIceArea(iRegion) = reductions(5)
           regionalFloatingIceVolume(iRegion) = reductions(6)
           regionalSumSfcMassBal(iRegion) = reductions(7)
           regionalSumGroundedSfcMassBal(iRegion) = reductions(8)
           regionalSumFloatingSfcMassBal(iRegion) = reductions(9)
           regionalSumBasalMassBal(iRegion) = reductions(10)
           regionalSumGroundedBasalMassBal(iRegion) = reductions(11)
           regionalSumFloatingBasalMassBal(iRegion) = reductions(12)
           regionalSumCalvingFlux(iRegion) = reductions(13)
           regionalSumFaceMeltingFlux(iRegion) = reductions(14)
           regionalVolumeAboveFloatation(iRegion) = reductions(15)
           regionalSumGroundingLineFlux(iRegion) = reductions(16)
           regionalSumGroundingLineMigrationFlux(iRegion) = reductions(17)
           if (config_SGH) then
              regionalSumSubglacialWaterVolume(iRegion) = reductions(18)
              regionalSumBasalMeltInput(iRegion) = reductions(19)
              regionalSumExternalWaterInput(iRegion) = reductions(20)
              regionalSumChannelMelt(iRegion) = reductions(21)
              regionalSumSubglacialLakeVolume(iRegion) = reductions(22)
              regionalSumSubglacialLakeArea(iRegion) = reductions(23)
              regionalSumDistWaterFluxMarineMargin(iRegion) = reductions(24)
              regionalSumDistWaterFluxTerrestrialMargin(iRegion) = reductions(25)
              regionalSumChnlWaterFluxMarineMargin(iRegion) = reductions(26)
              regionalSumChnlWaterFluxTerrestrialMargin(iRegion) = reductions(27)
              regionalSumFlotationFraction(iRegion) = reductions(28)
           endif

           if (regionalIceArea(iRegion) > 0.0_RKIND) then
              regionalIceThicknessMean(iRegion) = regionalIceVolume(iRegion) / regionalIceArea(iRegion)
              regionalAvgNetAccumulation(iRegion) = regionalSumSfcMassBal(iRegion) / regionalIceArea(iRegion) / rhoi
           else
              regionalIceThicknessMean(iRegion) = 0.0_RKIND
              regionalAvgNetAccumulation(iRegion) = 0.0_RKIND
           endif
           if (regionalGroundedIceArea(iRegion) > 0.0_RKIND) then
              regionalAvgGroundedBasalMelt(iRegion) = -1.0_RKIND * regionalSumGroundedBasalMassBal(iRegion) / &
                 regionalGroundedIceArea(iRegion) / rhoi
           else
              regionalAvgGroundedBasalMelt(iRegion) = 0.0_RKIND
           endif
           if (regionalFloatingIceArea(iRegion) > 0.0_RKIND) then
              regionalAvgSubshelfMelt(iRegion) = -1.0_RKIND * regionalSumFloatingBasalMassBal(iRegion) / &
                 regionalFloatingIceArea(iRegion) / rhoi
           else
              regionalAvgSubshelfMelt(iRegion) = 0.0_RKIND
           endif

           if (config_SGH) then
              if (regionalGroundedIceArea(iRegion) > 0.0_RKIND) then
                  regionalAvgFlotationFraction(iRegion) = regionalSumFlotationFraction(iRegion) / regionalGroundedIceArea(iRegion)
              else
                  regionalAvgFlotationFraction(iRegion) = 0.0_RKIND
              endif
           endif

        block => block % next
        end do

        ! ---------------------------
        ! compute mins
        ! ---------------------------
        mins = 0.0_RKIND
        reductions = 0.0_RKIND
        mins(1) = blockRegionThickMin(iRegion)
        nVars = 1
        call mpas_dmpar_min_real_array(dminfo, nVars, mins(1:nVars), reductions(1:nVars))

        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMin', regionalIceThicknessMin)
           regionalIceThicknessMin(iRegion) = reductions(1)
           block => block % next
        end do

        ! ---------------------------
        ! compute maxes
        ! ---------------------------
        maxes = 0.0_RKIND
        reductions = 0.0_RKIND
        maxes(1) = blockRegionThickMax(iRegion)
        maxes(2) = blockRegionMaxSurfaceSpeed(iRegion)* scyr  ! convert units to match Registry
        maxes(3) = blockRegionMaxBasalSpeed(iRegion)* scyr    ! convert units to match Registry
        nVars = 3
        call mpas_dmpar_max_real_array(dminfo, nVars, maxes(1:nVars), reductions(1:nVars))

        block => domain % blocklist
        do while (associated(block))
           call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMax', regionalIceThicknessMax)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSurfaceSpeedMax', regionalSurfaceSpeedMax)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalBasalSpeedMax', regionalBasalSpeedMax)
           regionalIceThicknessMax(iRegion) = reductions(1)
           regionalSurfaceSpeedMax(iRegion) = reductions(2)
           regionalBasalSpeedMax(iRegion) = reductions(3)
           block => block % next
        end do

      end do    ! loop over regions

      ! deallocate storage for sums over blocks
      deallocate(blockRegionThickMin); deallocate(blockRegionThickMax)
      deallocate(blockRegionMaxBasalSpeed); deallocate(blockRegionMaxSurfaceSpeed)
      deallocate(cellRawContribution); deallocate(edgeRawContribution)
      deallocate(cellSumContribution); deallocate(edgeSumContribution)
      deallocate(edgeUpwindCell)
      if (config_SGH) then
         ! local variable needed for calculating numerator of regionalAvgFlotationFraction
         deallocate(regionalSumFlotationFraction)
      endif 

   end subroutine li_compute_regional_stats

!***********************************************************************
!
!  routine li_restart_regional_stats
!
!> \brief   Save restart for MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_restart_regional_stats(domain, memberName, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character (len=*), intent(in) :: memberName

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

   end subroutine li_restart_regional_stats

!***********************************************************************
!
!  routine li_finalize_regional_stats
!
!> \brief   Finalize MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_finalize_regional_stats(domain, memberName, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character (len=*), intent(in) :: memberName

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

   end subroutine li_finalize_regional_stats

end module li_regional_stats

! vim: foldmethod=marker
