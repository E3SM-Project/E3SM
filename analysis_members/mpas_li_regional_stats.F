! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), pointer ::  config_ice_density
      real (kind=RKIND), pointer ::  deltat
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  thickness
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBal
      real (kind=RKIND), dimension(:), pointer ::  basalMassBal
      real (kind=RKIND), dimension(:), pointer ::  calvingThickness

      integer, dimension(:,:), pointer :: regionCellMasks
      integer, dimension(:), pointer :: cellMask
      integer, pointer :: nCellsSolve, nRegions, nRegionGroups !, maxRegionsInGroup  !! maxRegionsInGroup not needed / used yet
      integer :: k, iCell, iRegion, iGroup

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), dimension(:), pointer ::  regionalIceArea, regionalIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalGroundedIceArea, regionalGroundedIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalFloatingIceArea, regionalFloatingIceVolume
      real (kind=RKIND), dimension(:), pointer ::  regionalIceThicknessMin, regionalIceThicknessMax, regionalIceThicknessMean
      real (kind=RKIND), dimension(:), pointer ::  regionalSumSfcMassBal, regionalSumBasalMassBal
      real (kind=RKIND), dimension(:), pointer ::  regionalSumCalvingFlux


      ! storage for sums over blocks
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionIceArea, blockSumRegionIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionGroundedIceArea, blockSumRegionGroundedIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionFloatingIceArea, blockSumRegionFloatingIceVolume
      real (kind=RKIND), dimension(:), allocatable ::  blockRegionThickMin, blockRegionThickMax
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionSfcMassBal, blockSumRegionBasalMassBal
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionCalvingFlux

!      real (kind=RKIND) ::  blockSumCalvingFlux


      ! local parameters
      real (kind=RKIND), parameter :: scyr = 31536000.0_RKIND  ! seconds per 365-day year

      err = 0

      dminfo = domain % dminfo

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

         !! get dimensions assocated with region masks
         call mpas_pool_get_dimension(block % dimensions, 'nRegions', nRegions)
         call mpas_pool_get_dimension(block % dimensions, 'nRegionGroups', nRegionGroups)

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'regions', regionsPool)
         call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)

         ! get region cell masks from regionMasks.nc input file
         call mpas_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)

         ! get values from regional stats pool
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceArea', regionalIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceVolume', regionalIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalGroundedIceArea', regionalGroundedIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalGroundedIceVolume', regionalGroundedIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalFloatingIceArea', regionalFloatingIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalFloatingIceVolume', regionalFloatingIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMin', regionalIceThicknessMin)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMax', regionalIceThicknessMax)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceThicknessMean', regionalIceThicknessMean)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumSfcMassBal', regionalSumSfcMassBal)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumBasalMassBal', regionalSumBasalMassBal)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalSumCalvingFlux', regionalSumCalvingFlux)

         ! allocate & initialize sums over blocks to 0
         allocate(blockSumRegionIceArea(nRegions)); allocate(blockSumRegionIceVolume(nRegions))
         allocate(blockSumRegionGroundedIceArea(nRegions)); allocate(blockSumRegionGroundedIceVolume(nRegions))
         allocate(blockSumRegionFloatingIceArea(nRegions)); allocate(blockSumRegionFloatingIceVolume(nRegions))
         allocate(blockRegionThickMin(nRegions)); allocate(blockRegionThickMax(nRegions))
         allocate(blockSumRegionSfcMassBal(nRegions)); allocate(blockSumRegionBasalMassBal(nRegions))
         allocate(blockSumRegionCalvingFlux(nRegions))

         blockSumRegionIceArea = 0.0_RKIND; blockSumRegionIceVolume = 0.0_RKIND
         blockSumRegionGroundedIceArea = 0.0_RKIND; blockSumRegionGroundedIceVolume = 0.0_RKIND
         blockSumRegionFloatingIceArea = 0.0_RKIND; blockSumRegionFloatingIceVolume = 0.0_RKIND
         blockRegionThickMin = 10000.0_RKIND; blockRegionThickMax = 0.0_RKIND
         blockSumRegionSfcMassBal = 0.0_RKIND; blockSumRegionBasalMassBal = 0.0_RKIND
         blockSumRegionCalvingFlux = 0.0_RKIND

         do iCell = 1,nCellsSolve      ! loop over cells
 !         do iGroup = 1,nRegionGroups  ! loop over groups
           do iRegion = 1,nRegions     ! loop over regions


            ! regional sums of ice area over cells with ice (m^2)
            blockSumRegionIceArea(iRegion) = blockSumRegionIceArea(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
               * areaCell(iCell) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) )

            ! regional sums of ice volume over cells with ice (m^3)
            blockSumRegionIceVolume(iRegion) = blockSumRegionIceVolume(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
               * areaCell(iCell) * thickness(iCell) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) )

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

            ! regional sum of sfc mass balance (kg yr^{-1})  !! NOTE: still needs to be tested with realistic SMB field !!
            blockSumRegionSfcMassBal(iRegion) = blockSumRegionSfcMassBal(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
                * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) * areaCell(iCell) * sfcMassBal(iCell) * scyr )

            ! regional sum of basal mass balance (kg yr^{-1})
               !! NOTE: still needs to be tested with realistic basal mass balance field !!
            blockSumRegionBasalMassBal(iRegion) = blockSumRegionBasalMassBal(iRegion) + &
              ( real(regionCellMasks(iRegion,iCell),RKIND) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * basalMassBal(iCell) * scyr )

            ! regional sum of mass lass due do calving (kg yr^{-1})
            blockSumRegionCalvingFlux(iRegion) = blockSumRegionCalvingFlux(iRegion) + &
            ( real(regionCellMasks(iRegion,iCell),RKIND) * real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
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


           end do ! end loop over regions
 !         end do ! end loop over groups
         end do ! end loop over cells

         block => block % next

      end do    ! end loop over blocks

      ! compute sums over all procs
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceArea, regionalIceArea )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceVolume, regionalIceVolume )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionGroundedIceArea, regionalGroundedIceArea )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionGroundedIceVolume, regionalGroundedIceVolume )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionFloatingIceArea, regionalFloatingIceArea )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionFloatingIceVolume, regionalFloatingIceVolume )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionSfcMassBal, regionalSumSfcMassBal )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionBasalMassBal, regionalSumBasalMassBal )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionCalvingFlux, regionalSumCalvingFlux )

      ! find min, max, mean thickness over all procs
      call mpas_dmpar_min_real_array(dminfo, nRegions, blockRegionThickMin, regionalIceThicknessMin)
      call mpas_dmpar_max_real_array(dminfo, nRegions, blockRegionThickMax, regionalIceThicknessMax)
      regionalIceThicknessMean = regionalIceVolume / regionalIceArea


      ! deallocate storage for sums over blocks
      deallocate(blockSumRegionIceArea); deallocate(blockSumRegionIceVolume)
      deallocate(blockSumRegionGroundedIceArea); deallocate(blockSumRegionGroundedIceVolume)
      deallocate(blockSumRegionFloatingIceArea); deallocate(blockSumRegionFloatingIceVolume)
      deallocate(blockRegionThickMin); deallocate(blockRegionThickMax)
      deallocate(blockSumRegionSfcMassBal); deallocate(blockSumRegionBasalMassBal)
      deallocate(blockSumRegionCalvingFlux)

      ! Even though some variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

         ! assignment of final regionalStatsAM variables could occur here.

         block => block % next
      end do

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
