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
      real (kind=RKIND), pointer ::  totalLakeVolume
      real (kind=RKIND), pointer ::  totalLakeArea
      real (kind=RKIND), pointer ::  totalGLMeltFlux
      real (kind=RKIND), pointer ::  totalTerrestrialMeltFlux
      real (kind=RKIND), pointer ::  totalChannelGLMeltFlux
      real (kind=RKIND), pointer ::  totalChannelTerrestrialMeltFlux
      real (kind=RKIND), pointer ::  totalFlotationFraction
      real (kind=RKIND), pointer ::  avgFlotationFraction

      ! scalar reductions over all blocks on this processor
      real (kind=RKIND) ::  blockSumIceArea, blockSumIceVolume
      real (kind=RKIND) ::  blockSumVAF
      real (kind=RKIND) ::  blockSumGroundedIceArea, blockSumGroundedIceVolume
      real (kind=RKIND) ::  blockSumFloatingIceArea, blockSumFloatingIceVolume
      real (kind=RKIND) ::  blockThickMin, blockThickMax
      real (kind=RKIND) ::  blockSumSfcMassBal, blockSumBasalMassBal
      real (kind=RKIND) ::  blockSumGroundedSfcMassBal, blockSumFloatingSfcMassBal
      real (kind=RKIND) ::  blockSumGroundedBasalMassBal, blockSumFloatingBasalMassBal
      real (kind=RKIND) ::  blockSumCalvingFlux
      real (kind=RKIND) ::  blockSumFaceMeltingFlux
      real (kind=RKIND) ::  blockMaxSurfaceSpeed
      real (kind=RKIND) ::  blockMaxBasalSpeed
      real (kind=RKIND) ::  blockGLflux
      real (kind=RKIND) ::  blockGLMigrationFlux
      real (kind=RKIND) ::  blockSumSubglacialWaterVolume
      real (kind=RKIND) ::  blockSumBasalMeltInput
      real (kind=RKIND) ::  blockSumExternalWaterInput
      real (kind=RKIND) ::  blockSumChannelMelt
      real (kind=RKIND) ::  blockSumLakeVolume
      real (kind=RKIND) ::  blockSumLakeArea
      real (kind=RKIND) ::  blockSumGLMeltFlux
      real (kind=RKIND) ::  blockSumTerrestrialMeltFlux
      real (kind=RKIND) ::  blockSumChannelGLMeltFlux
      real (kind=RKIND) ::  blockSumChannelTerrestrialMeltFlux
      real (kind=RKIND) ::  blockSumFlotationFraction

      ! Local variables for calculations

      ! variables for processing stats
      integer, parameter :: kMaxVariables = 35 ! Increase if number of stats increase 
      integer :: nVars
      real (kind=RKIND), dimension(kMaxVariables) :: reductions, sums, mins, maxes

      err = 0

      dminfo = domain % dminfo

      ! initialize sums over blocks to 0
      blockSumIceArea = 0.0_RKIND
      blockSumIceVolume = 0.0_RKIND
      blockSumVAF = 0.0_RKIND
      blockSumGroundedIceArea = 0.0_RKIND
      blockSumGroundedIceVolume = 0.0_RKIND
      blockSumFloatingIceArea = 0.0_RKIND
      blockSumFloatingIceVolume = 0.0_RKIND
      blockSumSfcMassBal = 0.0_RKIND
      blockSumGroundedSfcMassBal = 0.0_RKIND
      blockSumFloatingSfcMassBal = 0.0_RKIND
      blockSumBasalMassBal = 0.0_RKIND
      blockSumGroundedBasalMassBal = 0.0_RKIND
      blockSumFloatingBasalMassBal = 0.0_RKIND
      blockSumCalvingFlux = 0.0_RKIND
      blockSumFaceMeltingFlux = 0.0_RKIND
      blockGLflux = 0.0_RKIND
      blockGLMigrationFlux = 0.0_RKIND
      blockSumSubglacialWaterVolume = 0.0_RKIND
      blockSumBasalMeltInput = 0.0_RKIND
      blockSumExternalWaterInput = 0.0_RKIND
      blockSumChannelMelt = 0.0_RKIND
      blockSumLakeVolume = 0.0_RKIND
      blockSumLakeArea = 0.0_RKIND
      blockSumGLMeltFlux = 0.0_RKIND
      blockSumTerrestrialMeltFlux = 0.0_RKIND
      blockSumChannelGLMeltFlux = 0.0_RKIND
      blockSumChannelTerrestrialMeltFlux = 0.0_RKIND
      blockSumFlotationFraction = 0.0_RKIND

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

            ! sums of ice area and volume over cells (m^2 and m^3)i
            blockSumIceArea = blockSumIceArea + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell)
            blockSumIceVolume = blockSumIceVolume + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            blockSumVAF = blockSumVAF + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) * areaCell(iCell) * &
                ( thickness(iCell) + (rhow / rhoi) * min(0.0_RKIND, (bedTopography(iCell) - config_sea_level)) )

            blockSumGroundedIceArea = blockSumGroundedIceArea + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
                 * areaCell(iCell)
           
            blockSumGroundedIceVolume = blockSumGroundedIceVolume + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            blockSumFloatingIceArea = blockSumFloatingIceArea + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell)
           
            blockSumFloatingIceVolume = blockSumFloatingIceVolume + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

            ! max, min thickness values (m)
            if( thickness(iCell) > blockThickMax)then
                blockThickMax = thickness(iCell)
            endif
            if( thickness(iCell) < blockThickMin .and. thickness(iCell) > 0.0_RKIND)then
                blockThickMin = thickness(iCell)
            endif

            ! SMB (kg yr^{-1})
            blockSumSfcMassBal = blockSumSfcMassBal + areaCell(iCell) * sfcMassBalApplied(iCell) * scyr
            blockSumGroundedSfcMassBal = blockSumGroundedSfcMassBal + areaCell(iCell) * groundedSfcMassBalApplied(iCell) * scyr
            blockSumFloatingSfcMassBal = blockSumFloatingSfcMassBal + &
               (sfcMassBalApplied(iCell) - groundedSfcMassBalApplied(iCell)) * areaCell(iCell) * scyr
      
            ! BMB (kg yr-1)
            blockSumBasalMassBal = blockSumBasalMassBal + areaCell(iCell) * basalMassBalApplied(iCell) * scyr
            blockSumGroundedBasalMassBal = blockSumGroundedBasalMassBal + areaCell(iCell) * groundedBasalMassBalApplied(iCell) * scyr
            blockSumFloatingBasalMassBal = blockSumFloatingBasalMassBal + areaCell(iCell) * floatingBasalMassBalApplied(iCell) * scyr

            ! mass loss due do calving (kg yr^{-1})
            blockSumCalvingFlux = blockSumCalvingFlux + calvingThickness(iCell) * &
               areaCell(iCell) * rhoi / ( deltat / scyr )

            ! mass loss due to face-melting (kg yr^{-1})
            blockSumFaceMeltingFlux = blockSumFaceMeltingFlux + faceMeltingThickness(iCell) * &
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
            blockGLMigrationFlux = blockGLMigrationFlux + groundedToFloatingThickness(iCell) * areaCell(iCell) &
               * rhoi / (deltat / scyr)  ! convert from m to kg/yr

            !! Subglacial Hydrology Calculations
            if (config_SGH) then

               ! Subglacial Water Volume
               blockSumSubglacialWaterVolume = blockSumSubglacialWaterVolume + waterThickness(iCell) * areaCell(iCell)

               ! Basal melt input
               blockSumBasalMeltInput = blockSumBasalMeltInput + basalMeltInput(iCell) * areaCell(iCell)

               ! External water input
               blockSumExternalWaterInput = blockSumExternalWaterInput + externalWaterInput(iCell) * areaCell(iCell)

               ! Lake Volume
               if (waterThickness(iCell) > bedBumpMax) then
                  blockSumLakeVolume = blockSumLakeVolume + (waterThickness(iCell) - bedBumpMax) * areaCell(iCell)
               endif
                 
               ! Lake Area
               if (waterThickness(iCell) > bedBumpMax) then
                   blockSumLakeArea = blockSumLakeArea + areaCell(iCell)
               endif

               ! Area-weighted flotation fraction for grounded ice
               if (li_mask_is_grounded_ice(cellMask(iCell))) then
                   blockSumFlotationFraction = blockSumFlotationFraction + ( waterPressure(iCell) / rhow / gravity / thickness(iCell) ) * areaCell(iCell)
               endif
            endif


         end do ! end loop over cells

         if (config_SGH) then
            ! Loop over edges
            do iEdge = 1, nEdgesSolve

               ! Flux across GL, units = kg/yr
               blockGLflux = blockGLflux + fluxAcrossGroundingLine(iEdge) * dvEdge(iEdge) &
                         * scyr * rhoi ! convert from m^2/s to kg/yr

               ! Channel Melt
               blockSumChannelMelt = blockSumChannelMelt + abs(channelMelt(iEdge) * dcEdge(iEdge))

               ! Meltwater Flux across the grounding line
               blockSumGLMeltFlux = blockSumGLMeltFlux + abs(hydroMarineMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Flux across terrestrial margins
               blockSumTerrestrialMeltFlux = blockSumTerrestrialMeltFlux + abs(hydroTerrestrialMarginMask(iEdge) * waterFlux(iEdge) * dvEdge(iEdge) * rho_water)

               ! Meltwater Discharge in channels across grounding line
               blockSumChannelGLMeltFlux = blockSumChannelGLMeltFlux + abs(hydroMarineMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)

               ! Meltwater discharge in channels across terrestrial margin
               blockSumChannelTerrestrialMeltFlux = blockSumChannelTerrestrialMeltFlux + abs( hydroTerrestrialMarginMask(iEdge) * channelDischarge(iEdge) * rho_water)

            end do ! end loop over edges
         endif
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
      sums = 0.0_RKIND
      reductions = 0.0_RKIND

      ! Fill the array with the values accumulated over all blocks on this processor.
      ! Note that the indices here need to be manually matched when unpacking the reduced
      ! values after the global reduction.
      ! Instructions: To add a new field, add a new sums(??) line with a larger index.
      !     Then update nVars to the new value.  Make sure kMaxVariables is large enough above.
      !     Then unpack the new index in the block loop into the appropriate variable.

      sums(1) = blockSumIceArea
      sums(2) = blockSumIceVolume
      sums(3) = blockSumGroundedIceArea
      sums(4) = blockSumGroundedIceVolume
      sums(5) = blockSumFloatingIceArea
      sums(6) = blockSumFloatingIceVolume
      sums(7) = blockSumSfcMassBal
      sums(8) = blockSumGroundedSfcMassBal
      sums(9) = blockSumFloatingSfcMassBal
      sums(10) = blockSumBasalMassBal
      sums(11) = blockSumGroundedBasalMassBal
      sums(12) = blockSumFloatingBasalMassBal
      sums(13) = blockSumCalvingFlux
      sums(14) = blockSumFaceMeltingFlux
      sums(15) = blockSumVAF
      if (config_SGH) then
         sums(16) = blockGLflux
         sums(17) = blockGLMigrationflux
         sums(18) = blockSumSubglacialWaterVolume
         sums(19) = blockSumBasalMeltInput
         sums(20) = blockSumExternalWaterInput
         sums(21) = blockSumChannelMelt
         sums(22) = blockSumLakeVolume
         sums(23) = blockSumLakeArea
         sums(24) = blockSumGLMeltFlux
         sums(25) = blockSumTerrestrialMeltFlux
         sums(26) = blockSumChannelGLMeltFlux
         sums(27) = blockSumChannelTerrestrialMeltFlux
         sums(28) = blockSumFlotationFraction
         nVars = 28
      else
        nVars = 15
      endif

      call mpas_dmpar_sum_real_array(dminfo, nVars, sums(1:nVars), reductions(1:nVars))

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
         if (config_SGH) then
            call mpas_pool_get_array(globalStatsAMPool, 'groundingLineFlux', groundingLineFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'groundingLineMigrationFlux', groundingLineMigrationFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'totalSubglacialWaterVolume', totalSubglacialWaterVolume)
            call mpas_pool_get_array(globalStatsAMPool, 'totalBasalMeltInput', totalBasalMeltInput)
            call mpas_pool_get_array(globalStatsAMPool, 'totalExternalWaterInput', totalExternalWaterInput)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChannelMelt', totalChannelMelt)
            call mpas_pool_get_array(globalStatsAMPool, 'totalLakeVolume', totalLakeVolume)
            call mpas_pool_get_array(globalStatsAMPool, 'totalLakeArea', totalLakeArea)
            call mpas_pool_get_array(globalStatsAMPool, 'totalGLMeltFlux',totalGLMeltFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'totalTerrestrialMeltFlux', totalTerrestrialMeltFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChannelGLMeltFlux',totalChannelGLMeltFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'totalChannelTerrestrialMeltFlux', totalChannelTerrestrialMeltFlux)
            call mpas_pool_get_array(globalStatsAMPool, 'totalFlotationFraction', totalFlotationFraction)
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
         if (config_SGH) then
            groundingLineFlux = reductions(16)
            groundingLineMigrationFlux = reductions(17)
            totalSubglacialWaterVolume = reductions(18)
            totalBasalMeltInput = reductions(19)
            totalExternalWaterInput = reductions(20)
            totalChannelMelt = reductions(21)
            totalLakevolume = reductions(22)
            totalLakeArea = reductions(23)
            totalGLMeltFlux = reductions(24)
            totalTerrestrialMeltFlux = reductions(25)
            totalChannelGLMeltFlux = reductions(26)
            totalChannelTerrestrialMeltFlux = reductions(27)
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
