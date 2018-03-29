! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
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

   subroutine li_init_global_stats(domain, memberName, err)!{{{

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

   end subroutine li_init_global_stats!}}}

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

   subroutine li_compute_global_stats(domain, memberName, timeLevel, err)!{{{

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
      integer, dimension(:), pointer :: cellMask
      integer, dimension(:), pointer :: edgeMask
      integer, dimension(:,:), pointer :: cellsOnEdge
      integer, pointer :: nCellsSolve
      integer, pointer :: nEdgesSolve
      integer, pointer :: nVertLevels

      ! config options needed
      real (kind=RKIND), pointer :: config_sea_level
      real (kind=RKIND), pointer :: rhoi ! config_ice_density
      real (kind=RKIND), pointer :: rhow ! config_ocean_density

      ! Local counters
      integer :: k, iCell, iEdge
      integer :: iCell1, iCell2

      ! scalars to be calculated here from global reductions
      real (kind=RKIND), pointer ::  totalIceArea, totalIceVolume
      real (kind=RKIND), pointer ::  volumeAboveFloatation
      real (kind=RKIND), pointer ::  groundedIceArea, groundedIceVolume
      real (kind=RKIND), pointer ::  floatingIceArea, floatingIceVolume
      real (kind=RKIND), pointer ::  iceThicknessMax, iceThicknessMin, iceThicknessMean
      real (kind=RKIND), pointer ::  totalSfcMassBal, totalBasalMassBal
      real (kind=RKIND), pointer ::  totalGroundedBasalMassBal, totalFloatingBasalMassBal
      real (kind=RKIND), pointer ::  avgNetAccumulation
      real (kind=RKIND), pointer ::  avgGroundedBasalMelt
      real (kind=RKIND), pointer ::  avgSubshelfMelt
      real (kind=RKIND), pointer ::  totalCalvingFlux
      real (kind=RKIND), pointer ::  groundingLineFlux
      real (kind=RKIND), pointer ::  surfaceSpeedMax
      real (kind=RKIND), pointer ::  basalSpeedMax

      ! scalar reductions over all blocks on this processor
      real (kind=RKIND) ::  blockSumIceArea, blockSumIceVolume
      real (kind=RKIND) ::  blockSumVAF
      real (kind=RKIND) ::  blockSumGroundedIceArea, blockSumGroundedIceVolume
      real (kind=RKIND) ::  blockSumFloatingIceArea, blockSumFloatingIceVolume
      real (kind=RKIND) ::  blockThickMin, blockThickMax
      real (kind=RKIND) ::  blockSumSfcMassBal, blockSumBasalMassBal
      real (kind=RKIND) ::  blockSumGroundedBasalMassBal, blockSumFloatingBasalMassBal
      real (kind=RKIND) ::  blockSumCalvingFlux
      real (kind=RKIND) ::  blockMaxSurfaceSpeed
      real (kind=RKIND) ::  blockMaxBasalSpeed
      real (kind=RKIND) ::  blockGLflux

      ! Local variables for calculations
      real (kind=RKIND) :: fluxSign

      ! variables for processing stats
      integer, parameter :: kMaxVariables = 16 ! Increase if number of stats increase
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
      blockSumBasalMassBal = 0.0_RKIND
      blockSumGroundedBasalMassBal = 0.0_RKIND
      blockSumFloatingBasalMassBal = 0.0_RKIND
      blockSumCalvingFlux = 0.0_RKIND
      blockGLflux = 0.0_RKIND

      ! initialize max to 0, min to large number
      blockThickMin = 100000.0_RKIND
      blockThickMax = 0.0_RKIND
      blockMaxSurfaceSpeed = 0.0_RKIND
      blockMaxBasalSpeed = 0.0_RKIND

      ! Get needed configs
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', rhoi)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', rhow)

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'velocity', velocityPool)

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

         ! loop over cells
         do iCell = 1,nCellsSolve

            ! sums of ice area and volume over cells (m^2 and m^3)
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

            ! sfc and basal mass balance (kg yr^{-1})
            !SFP: These calculations need to be tested still
            blockSumSfcMassBal = blockSumSfcMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * sfcMassBal(iCell) * scyr
            blockSumBasalMassBal = blockSumBasalMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * basalMassBal(iCell) * scyr
            blockSumGroundedBasalMassBal = blockSumGroundedBasalMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * groundedBasalMassBal(iCell) * scyr
            blockSumFloatingBasalMassBal = blockSumFloatingBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * floatingBasalMassBal(iCell) * scyr

            ! mass lass due do calving (kg yr^{-1})
            blockSumCalvingFlux = blockSumCalvingFlux + calvingThickness(iCell) * &
               areaCell(iCell) * config_ice_density / ( deltat / scyr )

            ! max surface speed
            if (surfaceSpeed(iCell) > blockMaxSurfaceSpeed) then
               blockMaxSurfaceSpeed = surfaceSpeed(iCell)
            endif

            ! max basal speed
            if (basalSpeed(iCell) > blockMaxBasalSpeed) then
               blockMaxBasalSpeed = basalSpeed(iCell)
            endif

         end do ! end loop over cells

         ! Loop over edges
         do iEdge = 1, nEdgesSolve
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

               ! Loop over levels
               do k = 1, nVertLevels
                  ! Flux across GL, units = kg/yr
                  blockGLflux = blockGLflux + fluxSign * layerNormalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge) &
                      * scyr * config_ice_density  ! convert from m^3/s to kg/yr
               end do ! end loop over levels
            end if ! if GL

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
      sums(8) = blockSumBasalMassBal
      sums(9) = blockSumGroundedBasalMassBal
      sums(10) = blockSumFloatingBasalMassBal
      sums(11) = blockSumCalvingFlux
      sums(12) = blockSumVAF
      sums(13) = blockGLflux
      nVars = 13

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
         call mpas_pool_get_array(globalStatsAMPool, 'avgNetAccumulation', avgNetAccumulation)
         call mpas_pool_get_array(globalStatsAMPool, 'totalBasalMassBal', totalBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalGroundedBasalMassBal', totalGroundedBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'avgGroundedBasalMelt', avgGroundedBasalMelt)
         call mpas_pool_get_array(globalStatsAMPool, 'totalFloatingBasalMassBal', totalFloatingBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'avgSubshelfMelt', avgSubshelfMelt)
         call mpas_pool_get_array(globalStatsAMPool, 'totalCalvingFlux', totalCalvingFlux)
         call mpas_pool_get_array(globalStatsAMPool, 'groundingLineFlux', groundingLineFlux)

         totalIceArea = reductions(1)
         totalIceVolume = reductions(2)
         groundedIceArea = reductions(3)
         groundedIceVolume = reductions(4)
         floatingIceArea = reductions(5)
         floatingIceVolume = reductions(6)
         totalSfcMassBal = reductions(7)
         totalBasalMassBal = reductions(8)
         totalGroundedBasalMassBal = reductions(9)
         totalFloatingBasalMassBal = reductions(10)
         totalCalvingFlux = reductions(11)
         volumeAboveFloatation = reductions(12)
         groundingLineFlux = reductions(13)

         iceThicknessMean = totalIceVolume / totalIceArea
         avgNetAccumulation = totalSfcMassBal / totalIceArea / rhoi
         avgGroundedBasalMelt = -1.0_RKIND * totalGroundedBasalMassBal / groundedIceArea / rhoi
         avgSubshelfMelt = -1.0_RKIND * totalFloatingBasalMassBal / floatingIceArea / rhoi

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



   end subroutine li_compute_global_stats!}}}

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

   subroutine li_restart_global_stats(domain, memberName, err)!{{{

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

   end subroutine li_restart_global_stats!}}}

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

   subroutine li_finalize_global_stats(domain, memberName, err)!{{{

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

   end subroutine li_finalize_global_stats!}}}

end module li_global_stats

! vim: foldmethod=marker
