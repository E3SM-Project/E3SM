! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_li_global_stats
!
!> \brief MPAS land ice analysis mode member: mpas_li_global_stats
!> \author Stephen Price
!> \date   8-30-2015
!> \details
!>
!>
!-----------------------------------------------------------------------
module li_global_stats

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

   public :: li_init_global_stats, &
             li_compute_global_stats, &
             li_restart_global_stats, &
             li_finalize_global_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_init_global_stats
!
!> \brief   Initialize MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_init_global_stats(domain, memberName, err)

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

   end subroutine li_init_global_stats

!***********************************************************************
!
!  routine li_compute_global_stats
!
!> \brief   Compute MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_compute_global_stats(domain, memberName, timeLevel, err)

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

      type (mpas_pool_type), pointer :: globalStatsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: globalStatsAM
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: velocityPool
      type (mpas_pool_type), pointer :: hydroPool

      ! arrays, vars needed from other pools for calculations here
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

      real (kind=RKIND), dimension(:), pointer ::  waterThickness
      real (kind=RKIND), dimension(:), pointer ::  basalMeltInput
      real (kind=RKIND), dimension(:), pointer ::  externalWaterInput
      real (kind=RKIND), dimension(:), pointer ::  channelMelt
      real (kind=RKIND), dimension(:), pointer ::  waterFlux
      real (kind=RKIND), dimension(:), pointer ::  channelDischarge
      real (kind=RKIND), dimension(:), pointer ::  waterPressure

      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: hydroMarineMarginMask
      integer, dimension(:), pointer :: hydroTerrestrialMarginMask
      integer, pointer :: nCellsSolve
      integer, pointer :: nEdgesSolve

      ! config options needed
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoi ! config_ice_density
      real (kind=RKIND), pointer :: rhow ! config_ocean_density
      real (kind=RKIND), pointer :: bedBumpMax ! config_SGH_bed_roughness_max
      logical, pointer :: config_SGH

      ! Local counters
      integer :: k, iCell, iEdge, nCellsGrounded

      ! scalars to be calculated here from global reductions
      real (kind=RKIND), pointer ::  totalIceArea, totalIceVolume
      real (kind=RKIND), pointer ::  volumeAboveFloatation
      real (kind=RKIND), pointer ::  groundedIceArea, groundedIceVolume
      real (kind=RKIND), pointer ::  floatingIceArea, floatingIceVolume
      real (kind=RKIND), pointer ::  iceThicknessMax, iceThicknessMin, iceThicknessMean
      real (kind=RKIND), pointer ::  totalSfcMassBal, totalBasalMassBal
      real (kind=RKIND), pointer ::  totalGroundedSfcMassBal, totalFloatingSfcMassBal
      real (kind=RKIND), pointer ::  totalFaceMeltingFlux
      real (kind=RKIND), pointer ::  totalGroundedBasalMassBal, totalFloatingBasalMassBal
      real (kind=RKIND), pointer ::  avgNetAccumulation
      real (kind=RKIND), pointer ::  avgGroundedBasalMelt
      real (kind=RKIND), pointer ::  avgSubshelfMelt
      real (kind=RKIND), pointer ::  totalCalvingFlux
      real (kind=RKIND), pointer ::  groundingLineFlux
      real (kind=RKIND), pointer ::  groundingLineMigrationFlux
      real (kind=RKIND), pointer ::  surfaceSpeedMax
      real (kind=RKIND), pointer ::  basalSpeedMax
      real (kind=RKIND), pointer ::  totalSubglacialWaterVolume
      real (kind=RKIND), pointer ::  totalBasalMeltInput
      real (kind=RKIND), pointer ::  totalExternalWaterInput
      real (kind=RKIND), pointer ::  totalChannelMelt
      real (kind=RKIND), pointer ::  totalSubglacialLakeVolume
      real (kind=RKIND), pointer ::  totalSubglacialLakeArea
      real (kind=RKIND), pointer ::  totalDistWaterFluxMarineMargin
      real (kind=RKIND), pointer ::  totalDistWaterFluxTerrestrialMargin
      real (kind=RKIND), pointer ::  totalChnlWaterFluxMarineMargin
      real (kind=RKIND), pointer ::  totalChnlWaterFluxTerrestrialMargin
      real (kind=RKIND), pointer ::  avgFlotationFraction

      ! scalar reductions over all blocks on this processor
      ! (used only for the exact min/max reductions, which remain unchanged)
      real (kind=RKIND) ::  blockThickMin, blockThickMax
      real (kind=RKIND) ::  blockMaxSurfaceSpeed
      real (kind=RKIND) ::  blockMaxBasalSpeed

      ! Per-cell and per-edge contributions to the sum-type reductions above.
      ! Rather than accumulating a block-local scalar and reducing it with a plain
      ! MPI_Allreduce (mpas_dmpar_sum_real_array), which is not decomposition-independent
      ! due to floating-point non-associativity, we build the raw per-cell/per-edge
      ! contribution arrays here and reduce them once (after the block loop) with the
      ! reproducible mpas_global_sum_nfld routine.
      ! NOTE: MALI supports only a single block per MPI task, so these arrays are sized
      ! and filled for the single block on this task.
      integer, parameter :: nCellSumVars = 22, nEdgeSumVars = 6
      real (kind=RKIND), dimension(:,:), allocatable :: cellSumContribution, edgeSumContribution
      real (kind=RKIND), dimension(nCellSumVars) :: cellGlobalSums
      real (kind=RKIND), dimension(nEdgeSumVars) :: edgeGlobalSums

      ! column indices into cellSumContribution / cellGlobalSums
      integer, parameter :: cIceArea = 1, cIceVolume = 2, cGroundedIceArea = 3, cGroundedIceVolume = 4, &
                             cFloatingIceArea = 5, cFloatingIceVolume = 6, cSfcMassBal = 7, cGroundedSfcMassBal = 8, &
                             cFloatingSfcMassBal = 9, cBasalMassBal = 10, cGroundedBasalMassBal = 11, &
                             cFloatingBasalMassBal = 12, cCalvingFlux = 13, cFaceMeltingFlux = 14, cVAF = 15, &
                             cGLMigrationFlux = 16, cSubglacialWaterVolume = 17, cBasalMeltInput = 18, &
                             cExternalWaterInput = 19, cLakeVolume = 20, cLakeArea = 21, cFlotationFraction = 22

      ! column indices into edgeSumContribution / edgeGlobalSums
      integer, parameter :: eGLflux = 1, eChannelMelt = 2, eGLMeltFlux = 3, eTerrestrialMeltFlux = 4, &
                             eChannelGLMeltFlux = 5, eChannelTerrestrialMeltFlux = 6

      ! Local variables for calculations
      real (kind=RKIND) ::  totalFlotationFraction

      ! variables for processing stats
      integer, parameter :: kMaxVariables = 35 ! Increase if number of stats increase
      integer :: nVars
      real (kind=RKIND), dimension(kMaxVariables) :: reductions, mins, maxes

      err = 0

      dminfo = domain % dminfo

      ! initialize max to 0, min to large number
      blockThickMin = 100000.0_RKIND
      blockThickMax = 0.0_RKIND
      blockMaxSurfaceSpeed = 0.0_RKIND
      blockMaxBasalSpeed = 0.0_RKIND

      ! Get needed configs
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhow)
      call mpas_pool_get_config(liConfigs, 'config_SGH_bed_roughness_max', bedBumpMax)
      call mpas_pool_get_config(liConfigs, 'config_SGH', config_SGH)

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)
         call mpas_pool_get_subpool(block % structs, 'hydro', hydroPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
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

         if (.not. allocated(cellSumContribution)) allocate(cellSumContribution(nCellsSolve, nCellSumVars))
         if (.not. allocated(edgeSumContribution)) allocate(edgeSumContribution(nEdgesSolve, nEdgeSumVars))
         cellSumContribution(:,:) = 0.0_RKIND
         edgeSumContribution(:,:) = 0.0_RKIND

         if (config_SGH) then
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

         ! loop over cells
         do iCell = 1,nCellsSolve

            ! sums of ice area and volume over cells (m^2 and m^3)
            cellSumContribution(iCell, cIceArea) = real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell)
            cellSumContribution(iCell, cIceVolume) = real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            cellSumContribution(iCell, cVAF) = real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) * areaCell(iCell) * &
                ( thickness(iCell) + (rhow / rhoi) * min(0.0_RKIND, (bedTopography(iCell) - config_sea_level)) )

            cellSumContribution(iCell, cGroundedIceArea) = real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
                 * areaCell(iCell)

            cellSumContribution(iCell, cGroundedIceVolume) = real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            cellSumContribution(iCell, cFloatingIceArea) = real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell)

            cellSumContribution(iCell, cFloatingIceVolume) = real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            ! max, min thickness values (m)
            if( thickness(iCell) > blockThickMax)then
                blockThickMax = thickness(iCell)
            endif
            if( thickness(iCell) < blockThickMin .and. thickness(iCell) > 0.0_RKIND)then
                blockThickMin = thickness(iCell)
            endif

            ! SMB (kg yr^{-1})
            cellSumContribution(iCell, cSfcMassBal) = areaCell(iCell) * sfcMassBalApplied(iCell) * scyr
            cellSumContribution(iCell, cGroundedSfcMassBal) = areaCell(iCell) * groundedSfcMassBalApplied(iCell) * scyr
            cellSumContribution(iCell, cFloatingSfcMassBal) = &
               (sfcMassBalApplied(iCell) - groundedSfcMassBalApplied(iCell)) * areaCell(iCell) * scyr

            ! BMB (kg yr-1)
            cellSumContribution(iCell, cBasalMassBal) = areaCell(iCell) * basalMassBalApplied(iCell) * scyr
            cellSumContribution(iCell, cGroundedBasalMassBal) = areaCell(iCell) * groundedBasalMassBalApplied(iCell) * scyr
            cellSumContribution(iCell, cFloatingBasalMassBal) = areaCell(iCell) * floatingBasalMassBalApplied(iCell) * scyr

            ! mass loss due do calving (kg yr^{-1})
            cellSumContribution(iCell, cCalvingFlux) = calvingThickness(iCell) * &
               areaCell(iCell) * rhoi / ( deltat / scyr )

            ! mass loss due to face-melting (kg yr^{-1})
            cellSumContribution(iCell, cFaceMeltingFlux) = faceMeltingThickness(iCell) * &
               areaCell(iCell) * rhoi / ( deltat / scyr )

            ! max surface speed
            if (surfaceSpeed(iCell) > blockMaxSurfaceSpeed) then
               blockMaxSurfaceSpeed = surfaceSpeed(iCell)
            endif

            ! max basal speed
            if (basalSpeed(iCell) > blockMaxBasalSpeed) then
               blockMaxBasalSpeed = basalSpeed(iCell)
            endif

            ! GL migration flux
            cellSumContribution(iCell, cGLMigrationFlux) = groundedToFloatingThickness(iCell) * areaCell(iCell) &
               * rhoi / (deltat / scyr)  ! convert from m to kg/yr

            !! Subglacial Hydrology Calculations
            if (config_SGH) then

               ! Subglacial Water Volume
               cellSumContribution(iCell, cSubglacialWaterVolume) = waterThickness(iCell) * areaCell(iCell)

               ! Basal melt input
               cellSumContribution(iCell, cBasalMeltInput) = real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) * &
                                        basalMeltInput(iCell) * areaCell(iCell)

               ! External water input
               cellSumContribution(iCell, cExternalWaterInput) = &
                   real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) * externalWaterInput(iCell) * areaCell(iCell)

               ! Lake Volume
               if (waterThickness(iCell) > bedBumpMax) then
                  cellSumContribution(iCell, cLakeVolume) = (waterThickness(iCell) - bedBumpMax) * areaCell(iCell)
               endif

               ! Lake Area
               if (waterThickness(iCell) > bedBumpMax) then
                   cellSumContribution(iCell, cLakeArea) = areaCell(iCell)
               endif

               ! Area-weighted flotation fraction for grounded ice
               if (li_mask_is_grounded_ice(cellMask(iCell))) then
                   cellSumContribution(iCell, cFlotationFraction) = &
                      ( waterPressure(iCell) / rhoi / gravity / thickness(iCell) ) * areaCell(iCell)
               endif
            endif


         end do ! end loop over cells

         ! Loop over edges
         do iEdge = 1, nEdgesSolve

            ! Flux across GL, units = kg/yr
            edgeSumContribution(iEdge, eGLflux) = fluxAcrossGroundingLine(iEdge) * dvEdge(iEdge) &
                      * scyr * rhoi ! convert from m^2/s to kg/yr

            if (config_SGH) then

               ! Channel Melt
               edgeSumContribution(iEdge, eChannelMelt) = abs(channelMelt(iEdge) * dcEdge(iEdge))

               ! Meltwater Flux across the grounding line
               edgeSumContribution(iEdge, eGLMeltFlux) = abs(hydroMarineMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Flux across terrestrial margins
               edgeSumContribution(iEdge, eTerrestrialMeltFlux) = abs(hydroTerrestrialMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Discharge in channels across grounding line
               edgeSumContribution(iEdge, eChannelGLMeltFlux) = abs(hydroMarineMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)

               ! Meltwater discharge in channels across terrestrial margin
               edgeSumContribution(iEdge, eChannelTerrestrialMeltFlux) = abs( hydroTerrestrialMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)

            endif ! is SGH is turned on 
         end do ! end loop over edges
         
         block => block % next
      end do    ! end loop over blocks

      ! Reduce the per-cell/per-edge contribution arrays once, using the reproducible
      ! (decomposition-independent) global sum, rather than reducing block-local partial sums
      ! with a plain MPI_Allreduce (which is not reproducible due to floating-point
      ! non-associativity).
      cellGlobalSums(:) = mpas_global_sum_nfld(cellSumContribution, dminfo % comm)
      edgeGlobalSums(:) = mpas_global_sum_nfld(edgeSumContribution, dminfo % comm)
      deallocate(cellSumContribution)
      deallocate(edgeSumContribution)

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
      ! Note: cellGlobalSums and edgeGlobalSums have already been reduced above (reproducibly),
      ! so here we just unpack them into the reductions array using the same variable-to-index
      ! mapping that was previously used for the plain-sum reduction.
      reductions = 0.0_RKIND

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
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)

         call mpas_pool_get_array(globalStatsAMPool, 'totalIceArea', totalIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'totalIceVolume', totalIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'volumeAboveFloatation', volumeAboveFloatation)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceArea', floatingIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceVolume', floatingIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceArea', groundedIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceVolume', groundedIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMean', iceThicknessMean)
         call mpas_pool_get_array(globalStatsAMPool, 'totalSfcMassBal', totalSfcMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalGroundedSfcMassBal', totalGroundedSfcMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalFloatingSfcMassBal', totalFloatingSfcMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'avgNetAccumulation', avgNetAccumulation)
         call mpas_pool_get_array(globalStatsAMPool, 'totalBasalMassBal', totalBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalGroundedBasalMassBal', totalGroundedBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'avgGroundedBasalMelt', avgGroundedBasalMelt)
         call mpas_pool_get_array(globalStatsAMPool, 'totalFloatingBasalMassBal', totalFloatingBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'avgSubshelfMelt', avgSubshelfMelt)
         call mpas_pool_get_array(globalStatsAMPool, 'totalCalvingFlux', totalCalvingFlux)
         call mpas_pool_get_array(globalStatsAMPool, 'totalFaceMeltingFlux', totalFaceMeltingFlux)
         call mpas_pool_get_array(globalStatsAMPool, 'groundingLineFlux', groundingLineFlux)
         call mpas_pool_get_array(globalStatsAMPool, 'groundingLineMigrationFlux', groundingLineMigrationFlux)
         if (config_SGH) then
            call mpas_pool_get_array(globalStatsAMPool, 'totalSubglacialWaterVolume', totalSubglacialWaterVolume)
            call mpas_pool_get_array(globalStatsAMPool, 'totalBasalMeltInput', totalBasalMeltInput)
            call mpas_pool_get_array(globalStatsAMPool, 'totalExternalWaterInput', totalExternalWaterInput)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChannelMelt', totalChannelMelt)
            call mpas_pool_get_array(globalStatsAMPool, 'totalSubglacialLakeVolume', totalSubglacialLakeVolume)
            call mpas_pool_get_array(globalStatsAMPool, 'totalSubglacialLakeArea', totalSubglacialLakeArea)
            call mpas_pool_get_array(globalStatsAMPool, 'totalDistWaterFluxMarineMargin',totalDistWaterFluxMarineMargin)
            call mpas_pool_get_array(globalStatsAMPool, 'totalDistWaterFluxTerrestrialMargin', totalDistWaterFluxTerrestrialMargin)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChnlWaterFluxMarineMargin',totalChnlWaterFluxMarineMargin)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChnlWaterFluxTerrestrialMargin', totalChnlWaterFluxTerrestrialMargin)
            call mpas_pool_get_array(globalStatsAMPool, 'avgFlotationFraction', avgFlotationFraction)
         endif

         totalIceArea = reductions(1)
         totalIceVolume = reductions(2)
         groundedIceArea = reductions(3)
         groundedIceVolume = reductions(4)
         floatingIceArea = reductions(5)
         floatingIceVolume = reductions(6)
         totalSfcMassBal = reductions(7)
         totalGroundedSfcMassBal = reductions(8)
         totalFloatingSfcMassBal = reductions(9)
         totalBasalMassBal = reductions(10)
         totalGroundedBasalMassBal = reductions(11)
         totalFloatingBasalMassBal = reductions(12)
         totalCalvingFlux = reductions(13)
         totalFaceMeltingFlux = reductions(14)
         volumeAboveFloatation = reductions(15)
         groundingLineFlux = reductions(16)
         groundingLineMigrationFlux = reductions(17)
         if (config_SGH) then
            totalSubglacialWaterVolume = reductions(18)
            totalBasalMeltInput = reductions(19)
            totalExternalWaterInput = reductions(20)
            totalChannelMelt = reductions(21)
            totalSubglacialLakeVolume = reductions(22)
            totalSubglacialLakeArea = reductions(23)
            totalDistWaterFluxMarineMargin = reductions(24)
            totalDistWaterFluxTerrestrialMargin = reductions(25)
            totalChnlWaterFluxMarineMargin = reductions(26)
            totalChnlWaterFluxTerrestrialMargin = reductions(27)
            totalFlotationFraction = reductions(28)
         endif

         if (totalIceArea > 0.0_RKIND) then
            iceThicknessMean = totalIceVolume / totalIceArea
            avgNetAccumulation = totalSfcMassBal / totalIceArea / rhoi
         else
            iceThicknessMean = 0.0_RKIND
            avgNetAccumulation = 0.0_RKIND
         endif

         if (groundedIceArea > 0.0_RKIND) then
            avgGroundedBasalMelt = -1.0_RKIND * totalGroundedBasalMassBal / groundedIceArea / rhoi
         else
            avgGroundedBasalMelt = 0.0_RKIND
         endif

         if (floatingIceArea > 0.0_RKIND) then
            avgSubshelfMelt = -1.0_RKIND * totalFloatingBasalMassBal / floatingIceArea / rhoi
         else
            avgSubshelfMelt = 0.0_RKIND
         endif
         if (config_SGH) then
            if (groundedIceArea > 0.0_RKIND) then
                avgFlotationFraction = totalFlotationFraction / groundedIceArea
            else
                avgFlotationFraction = 0.0_RKIND
            endif
         endif



         block => block % next
      end do

      ! ---------------------------
      ! compute mins
      ! ---------------------------
      mins = 0.0_RKIND
      reductions = 0.0_RKIND
      mins(1) = blockThickMin
      nVars = 1
      call mpas_dmpar_min_real_array(dminfo, nVars, mins(1:nVars), reductions(1:nVars))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)

         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMin', iceThicknessMin)

         iceThicknessMin = reductions(1)

         block => block % next
      end do

      ! ---------------------------
      ! compute maxes
      ! ---------------------------
      maxes = 0.0_RKIND
      reductions = 0.0_RKIND
      maxes(1) = blockThickMax
      maxes(2) = blockMaxSurfaceSpeed * scyr  ! convert units to match Registry
      maxes(3) = blockMaxBasalSpeed * scyr  ! convert units to match Registry
      nVars = 3
      call mpas_dmpar_max_real_array(dminfo, nVars, maxes(1:nVars), reductions(1:nVars))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)

         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMax', iceThicknessMax)
         call mpas_pool_get_array(globalStatsAMPool, 'surfaceSpeedMax', surfaceSpeedMax)
         call mpas_pool_get_array(globalStatsAMPool, 'basalSpeedMax', basalSpeedMax)

         iceThicknessMax = reductions(1)
         surfaceSpeedMax = reductions(2)
         basalSpeedMax = reductions(3)

         block => block % next
      end do



   end subroutine li_compute_global_stats

!***********************************************************************
!
!  routine li_restart_global_stats
!
!> \brief   Save restart for MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_restart_global_stats(domain, memberName, err)

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

   end subroutine li_restart_global_stats

!***********************************************************************
!
!  routine li_finalize_global_stats
!
!> \brief   Finalize MPAS-Land Ice analysis member
!> \author  S. Price
!> \date    9/9/2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Land Ice analysis member.
!
!-----------------------------------------------------------------------

   subroutine li_finalize_global_stats(domain, memberName, err)

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

   end subroutine li_finalize_global_stats

end module li_global_stats

! vim: foldmethod=marker
