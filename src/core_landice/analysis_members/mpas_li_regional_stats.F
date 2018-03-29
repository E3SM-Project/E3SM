! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   subroutine li_init_regional_stats(domain, memberName, err)!{{{

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

   end subroutine li_init_regional_stats!}}}

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

   subroutine li_compute_regional_stats(domain, memberName, timeLevel, err)!{{{

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

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), pointer ::  config_ice_density
      real (kind=RKIND), pointer ::  deltat
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  dvEdge  
      real (kind=RKIND), dimension(:), pointer ::  thickness
      real (kind=RKIND), dimension(:), pointer ::  bedTopography
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBal
      real (kind=RKIND), dimension(:), pointer ::  basalMassBal
      real (kind=RKIND), dimension(:), pointer ::  groundedBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  floatingBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  calvingThickness
      real (kind=RKIND), dimension(:), pointer ::  surfaceSpeed  
      real (kind=RKIND), dimension(:), pointer ::  basalSpeed
      real (kind=RKIND), dimension(:,:), pointer ::  layerNormalVelocity 
      real (kind=RKIND), dimension(:,:), pointer ::  layerThicknessEdge

      ! config options needed
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoi ! config_ice_density
      real (kind=RKIND), pointer :: rhow ! config_ocean_density

      integer, dimension(:,:), pointer :: regionCellMasks
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: edgeMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nRegions, nRegionGroups !, maxRegionsInGroup  !! maxRegionsInGroup not needed / used yet
      integer, pointer :: nCellsSolve, nEdgesSolve, nVertLevels
      integer :: k, iCell, iCell1, iCell2, iEdge 
      integer :: iRegion, iGroup

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), dimension(:), pointer ::  regionalIceArea, regionalIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalVolumeAboveFloatation
      real (kind=RKIND), dimension(:), pointer ::  regionalGroundedIceArea, regionalGroundedIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalFloatingIceArea, regionalFloatingIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalIceThicknessMin, regionalIceThicknessMax, regionalIceThicknessMean
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSfcMassBal, regionalSumBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundedBasalMassBal, regionalSumFloatingBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumCalvingFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalSumGroundingLineFlux
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgNetAccumulation
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgGroundedBasalMelt
      real (kind=RKIND), dimension(:), pointer ::  regionalAvgSubshelfMelt
      real (kind=RKIND), dimension(:), pointer ::  regionalSurfaceSpeedMax
      real (kind=RKIND), dimension(:), pointer ::  regionalBasalSpeedMax 

      ! storage for sums over blocks
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionIceArea, blockSumRegionIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionVAF
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionGroundedIceArea, blockSumRegionGroundedIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionFloatingIceArea, blockSumRegionFloatingIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionThickMin, blockRegionThickMax
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionSfcMassBal, blockSumRegionBasalMassBal
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionGroundedBasalMassBal, blockSumRegionFloatingBasalMassBal
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionCalvingFlux
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionGLflux
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionMaxSurfaceSpeed
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionMaxBasalSpeed

      ! local variables
      real (kind=RKIND) :: fluxSign

      ! variables for processing stats
      integer, parameter :: kMaxVariables = 16 ! Increase if number of stats increase
      integer :: nVars
      real (kind=RKIND), dimension(kMaxVariables) :: reductions, sums, mins, maxes

      err = 0

      dminfo = domain % dminfo

      ! Get needed configs
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhow)

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
!         call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

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
         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
         call mpas_pool_get_array(geometryPool, 'groundedBasalMassBal', groundedBasalMassBal)
         call mpas_pool_get_array(geometryPool, 'floatingBasalMassBal', floatingBasalMassBal)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)
         call mpas_pool_get_array(geometryPool, 'layerThicknessEdge', layerThicknessEdge)
         call mpas_pool_get_array(velocityPool, 'surfaceSpeed', surfaceSpeed)
         call mpas_pool_get_array(velocityPool, 'basalSpeed', basalSpeed)
         call mpas_pool_get_array(velocityPool, 'layerNormalVelocity', layerNormalVelocity)

         ! get region cell masks from regionMasks.nc input file
         call mpas_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)

         ! allocate & initialize sums over blocks to 0
         allocate(blockSumRegionIceArea(nRegions)); allocate(blockSumRegionIceVolume(nRegions))
         allocate(blockSumRegionVAF(nRegions))
         allocate(blockSumRegionGroundedIceArea(nRegions)); allocate(blockSumRegionGroundedIceVolume(nRegions))
         allocate(blockSumRegionFloatingIceArea(nRegions)); allocate(blockSumRegionFloatingIceVolume(nRegions))
         allocate(blockRegionThickMin(nRegions)); allocate(blockRegionThickMax(nRegions))
         allocate(blockRegionMaxBasalSpeed(nRegions)); allocate(blockRegionMaxSurfaceSpeed(nRegions))
         allocate(blockSumRegionSfcMassBal(nRegions)); allocate(blockSumRegionBasalMassBal(nRegions))
         allocate(blockSumRegionGroundedBasalMassBal(nRegions)); allocate(blockSumRegionFloatingBasalMassBal(nRegions))
         allocate(blockSumRegionCalvingFlux(nRegions)); allocate(blockSumRegionGLflux(nRegions))

         blockSumRegionIceArea = 0.0_RKIND; blockSumRegionIceVolume = 0.0_RKIND
         blockSumRegionVAF = 0.0_RKIND
         blockSumRegionGroundedIceArea = 0.0_RKIND; blockSumRegionGroundedIceVolume = 0.0_RKIND
         blockSumRegionFloatingIceArea = 0.0_RKIND; blockSumRegionFloatingIceVolume = 0.0_RKIND
         blockRegionThickMin = 10000.0_RKIND; blockRegionThickMax = 0.0_RKIND
         blockRegionMaxBasalSpeed = 0.0_RKIND; blockRegionMaxSurfaceSpeed = 0.0_RKIND
         blockSumRegionSfcMassBal = 0.0_RKIND; blockSumRegionBasalMassBal = 0.0_RKIND
         blockSumRegionGroundedBasalMassBal = 0.0_RKIND; blockSumRegionFloatingBasalMassBal = 0.0_RKIND
         blockSumRegionCalvingFlux = 0.0_RKIND; blockSumRegionGLflux = 0.0_RKIND

         do iCell = 1,nCellsSolve      ! loop over cells
 !         do iGroup = 1,nRegionGroups  ! loop over groups
           do iRegion = 1,nRegions     ! loop over regions

            ! regional sums of ice area over cells with ice (m^2)
            blockSumRegionIceArea(iRegion) = blockSumRegionIceArea(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
               * areaCell(iCell) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of ice volume over cells with ice (m^3)
            blockSumRegionIceVolume(iRegion) = blockSumRegionIceVolume(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
               * areaCell(iCell) * thickness(iCell) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of ice volume above floatation over cells with ice (m^3)
            blockSumRegionVAF(iRegion) = blockSumRegionVAF(iRegion) + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
               * real(regionCellMasks(iRegion,iCell),RKIND) * areaCell(iCell) * ( thickness(iCell) + (rhow / rhoi) &
               * min(0.0_RKIND, (bedTopography(iCell) - config_sea_level)) )

            ! regional sums of grounded ice area over cells with ice (m^2)
            blockSumRegionGroundedIceArea(iRegion) = blockSumRegionGroundedIceArea(iRegion) + &
               ( real(regionCellMasks(iRegion,iCell),RKIND) * areaCell(iCell) &
                * real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of grounded ice volume over cells with ice (m^3)
            blockSumRegionGroundedIceVolume(iRegion) = blockSumRegionGroundedIceVolume(iRegion) + &
            ( real(regionCellMasks(iRegion,iCell),RKIND) * areaCell(iCell) &
              * thickness(iCell) * real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of floating ice area over cells with ice (m^2)
            blockSumRegionFloatingIceArea(iRegion) = blockSumRegionFloatingIceArea(iRegion) + &
              ( real(regionCellMasks(iRegion,iCell),RKIND) * areaCell(iCell) &
                * real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of floating ice volume over cells with ice (m^3)
            blockSumRegionFloatingIceVolume(iRegion) = blockSumRegionFloatingIceVolume(iRegion) + &
            ( real(regionCellMasks(iRegion,iCell),RKIND) * areaCell(iCell) &
              * thickness(iCell) * real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) )

            ! regional sum of sfc mass balance (kg yr^{-1})  
            blockSumRegionSfcMassBal(iRegion) = blockSumRegionSfcMassBal(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
                * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) * areaCell(iCell) * sfcMassBal(iCell) * scyr )

            ! regional sum of basal mass balance (kg yr^{-1})
            blockSumRegionBasalMassBal(iRegion) = blockSumRegionBasalMassBal(iRegion) + &
              ( real(regionCellMasks(iRegion,iCell),RKIND) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * basalMassBal(iCell) * scyr )

            ! regional sum of grounded basal mass balance (kg yr^{-1})
            blockSumRegionFloatingBasalMassBal(iRegion) = blockSumRegionFloatingBasalMassBal(iRegion) + &
              ( real(regionCellMasks(iRegion,iCell),RKIND) * real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * floatingBasalMassBal(iCell) * scyr )

            ! regional sum of floating basal mass balance (kg yr^{-1})
            blockSumRegionGroundedBasalMassBal(iRegion) = blockSumRegionGroundedBasalMassBal(iRegion) + &
              ( real(regionCellMasks(iRegion,iCell),RKIND) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * groundedBasalMassBal(iCell) * scyr )

            ! regional sum of mass lass due do calving (kg yr^{-1})
            blockSumRegionCalvingFlux(iRegion) = blockSumRegionCalvingFlux(iRegion) + &
            ( real(regionCellMasks(iRegion,iCell),RKIND) &
              * calvingThickness(iCell) * areaCell(iCell) * config_ice_density / ( deltat / scyr ) )

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
 !         end do ! end loop over groups
         end do ! end loop over cells


        do iEdge = 1,nEdgesSolve      ! loop over edges
 !         do iGroup = 1,nRegionGroups  ! loop over groups
           do iRegion = 1,nRegions     ! loop over regions

            if (li_mask_is_grounding_line(edgeMask(iEdge))) then
               ! Determine sign of this edge relative to GL
               ! (+=grounded to floating, -=floating to grounded)
               iCell1 = cellsOnEdge(1,iEdge)
               iCell2 = cellsOnEdge(2,iEdge)
               if (li_mask_is_grounded_ice(cellMask(iCell1))) then
                  fluxSign = 1.0_RKIND
               else
                  fluxSign = -1.0_RKIND
               endif

               ! first need to determine if the edge flux is being calc. on is part of the region of interest
               ! Use cellsOnEdge indices from above to index relevant regionCellMasks values and sum
               if(regionCellMasks(iRegion,iCell1) + regionCellMasks(iRegion,iCell2) > 0) then

                 ! Loop over levels
                 do k = 1, nVertLevels
                   ! Flux across GL, units = kg/yr
                   blockSumRegionGLflux(iRegion) = blockSumRegionGLflux(iRegion) + fluxSign * layerNormalVelocity(k, iEdge) * &
                      dvEdge(iEdge) * layerThicknessEdge(k, iEdge) &
                      * scyr * config_ice_density  ! convert from m^3/s to kg/yr
                 end do ! end loop over levels

               end if ! if edge is on cell in region of interest

            end if ! if GL

           end do ! end loop over regions
 !         end do ! end loop over groups
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

        sums = 0.0_RKIND
        reductions = 0.0_RKIND
  
        sums(1) = blockSumRegionIceArea(iRegion)
        sums(2) = blockSumRegionIceVolume(iRegion)
        sums(3) = blockSumRegionGroundedIceArea(iRegion)
        sums(4) = blockSumRegionGroundedIceVolume(iRegion)
        sums(5) = blockSumRegionFloatingIceArea(iRegion)
        sums(6) = blockSumRegionFloatingIceVolume(iRegion)
        sums(7) = blockSumRegionSfcMassBal(iRegion)
        sums(8) = blockSumRegionBasalMassBal(iRegion)
        sums(9) = blockSumRegionGroundedBasalMassBal(iRegion)
        sums(10) = blockSumRegionFloatingBasalMassBal(iRegion)
        sums(11) = blockSumRegionCalvingFlux(iRegion)
        sums(12) = blockSumRegionVAF(iRegion)
        sums(13) = blockSumRegionGLflux(iRegion)
        nVars = 13

        call mpas_dmpar_sum_real_array(dminfo, nVars, sums(1:nVars), reductions(1:nVars))

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
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumBasalMassBal', regionalSumBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundedBasalMassBal', regionalSumGroundedBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgGroundedBasalMelt', regionalAvgGroundedBasalMelt)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumFloatingBasalMassBal', regionalSumFloatingBasalMassBal)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalAvgSubshelfMelt', regionalAvgSubshelfMelt)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumCalvingFlux', regionalSumCalvingFlux)
           call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumGroundingLineFlux', regionalSumGroundingLineFlux)
  
           regionalIceArea(iRegion) = reductions(1)
           regionalIceVolume(iRegion) = reductions(2)
           regionalGroundedIceArea(iRegion) = reductions(3)
           regionalGroundedIceVolume(iRegion) = reductions(4)
           regionalFloatingIceArea(iRegion) = reductions(5)
           regionalFloatingIceVolume(iRegion) = reductions(6)
           regionalSumSfcMassBal(iRegion) = reductions(7)
           regionalSumBasalMassBal(iRegion) = reductions(8)
           regionalSumGroundedBasalMassBal(iRegion) = reductions(9)
           regionalSumFloatingBasalMassBal(iRegion) = reductions(10)
           regionalSumCalvingFlux(iRegion) = reductions(11)
           regionalVolumeAboveFloatation(iRegion) = reductions(12)
           regionalSumGroundingLineFlux(iRegion) = reductions(13)

           regionalIceThicknessMean(iRegion) = regionalIceVolume(iRegion) / regionalIceArea(iRegion)
           regionalAvgNetAccumulation(iRegion) = regionalSumSfcMassBal(iRegion) / regionalIceArea(iRegion) / rhoi
           regionalAvgGroundedBasalMelt(iRegion) = -1.0_RKIND * regionalSumGroundedBasalMassBal(iRegion) / &
           regionalGroundedIceArea(iRegion) / rhoi
           regionalAvgSubshelfMelt(iRegion) = -1.0_RKIND * regionalSumFloatingBasalMassBal(iRegion) / &
           regionalFloatingIceArea(iRegion) / rhoi

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
      deallocate(blockSumRegionIceArea); deallocate(blockSumRegionIceVolume)
      deallocate(blockSumRegionVAF)
      deallocate(blockSumRegionGroundedIceArea); deallocate(blockSumRegionGroundedIceVolume)
      deallocate(blockSumRegionFloatingIceArea); deallocate(blockSumRegionFloatingIceVolume)
      deallocate(blockRegionThickMin); deallocate(blockRegionThickMax)
      deallocate(blockRegionMaxBasalSpeed); deallocate(blockRegionMaxSurfaceSpeed)
      deallocate(blockSumRegionSfcMassBal); deallocate(blockSumRegionBasalMassBal)
      deallocate(blockSumRegionGroundedBasalMassBal); deallocate(blockSumRegionFloatingBasalMassBal)
      deallocate(blockSumRegionCalvingFlux); deallocate(blockSumRegionGLflux)

   end subroutine li_compute_regional_stats!}}}

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

   subroutine li_restart_regional_stats(domain, memberName, err)!{{{

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

   end subroutine li_restart_regional_stats!}}}

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

   subroutine li_finalize_regional_stats(domain, memberName, err)!{{{

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

   end subroutine li_finalize_regional_stats!}}}

end module li_regional_stats

! vim: foldmethod=marker
