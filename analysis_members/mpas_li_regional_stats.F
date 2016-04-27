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
!> \date   8-30-2015
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
!> \date    9/9/2015    
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

      type (mpas_pool_type), pointer :: regionalStatsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: regionalStatsAM
      type (mpas_pool_type), pointer :: geometryPool
      type (mpas_pool_type), pointer :: regionsPool           !!SFP: added

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  thickness 
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBal  
      real (kind=RKIND), dimension(:), pointer ::  basalMassBal

      integer, dimension(:,:), pointer :: regionCellMasks 
      integer, dimension(:), pointer :: cellMask 
      integer, pointer :: nCellsSolve, nRegions, nRegionGroups !, maxRegionsInGroup  !! SFP: is last one needed here?
      integer :: k, iCell, iRegion, iGroup

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), pointer ::  RegionalIceArea, RegionalIceVolume
      real (kind=RKIND), pointer ::  RegionalGroundedIceArea, RegionalGroundedIceVolume  
      real (kind=RKIND), pointer ::  RegionalFloatingIceArea, RegionalFloatingIceVolume    
!      real (kind=RKIND), pointer ::  WAIStotalBasalMassBal!, WAIStotalSfcMassBal 
!      real (kind=RKIND), pointer ::  EAIStotalBasalMassBal!, EAIStotalSfcMassBal 
!      real (kind=RKIND), pointer ::  FRISBasalMassBal, RossISBasalMassBal, AmeryISBasalMassBal, PineIslandThwaitesISBasalMassBal

      ! storage for sums over blocks
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionIceArea, blockSumRegionIceVolume
!      real (kind=RKIND) ::  blockSumRegionGroundedIceArea, blockSumRegionGroundedIceVolume
!      real (kind=RKIND) ::  blockSumRegionFloatingIceArea, blockSumRegionFloatingIceVolume
!      real (kind=RKIND) ::  blockSumWAISBasalMassBal!, blockSumWAISSfcMassBal
!      real (kind=RKIND) ::  blockSumEAISBasalMassBal!, blockSumEASISfcMassBal
!      real (kind=RKIND) ::  blockSumFRISBasalMassBal, blockSumRossISBasalMassBal 
!      real (kind=RKIND) ::  blockSumAmeryISBasalMassBal, blockSumPineIslandThwaitesISBasalMassBal

      ! local parameters
      real (kind=RKIND), parameter :: scyr = 31536000.0_RKIND  ! seconds per 365-day year

      err = 0

      dminfo = domain % dminfo

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

      !! SFP: get dimensions assocated with region masks
      call mpas_pool_get_dimension(block % dimensions, 'nRegions', nRegions)
      call mpas_pool_get_dimension(block % dimensions, 'nRegionGroups', nRegionGroups)

      ! allocate & initialize sums over blocks to 0
      allocate(blockSumRegionIceArea(nRegions)); allocate(blockSumRegionIceVolume(nRegions))

      blockSumRegionIceArea = 0.0_RKIND; blockSumRegionIceVolume = 0.0_RKIND
!      blockSumRegionGroundedIceArea = 0.0_RKIND; blockSumRegionGroundedIceVolume = 0.0_RKIND
!      blockSumRegionFloatingIceArea = 0.0_RKIND; blockSumRegionFloatingIceVolume = 0.0_RKIND

!      blockSumWAISBasalMassBal = 0.0_RKIND!; blockSumWAISSfcMassBal = 0.0_RKIND
!      blockSumEAISBasalMassBal = 0.0_RKIND!; blockSumEAISSfcMassBal = 0.0_RKIND
!      blockSumFRISBasalMassBal = 0.0_RKIND; blockSumRossISBasalMassBal = 0.0_RKIND
!      blockSumAmeryISBasalMassBal = 0.0_RKIND; blockSumPineIslandThwaitesISBasalMassBal = 0.0_RKIND

!      ! loop over blocks
!      block => domain % blocklist
!      do while (associated(block))

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
!         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)

         !!SFP: get region cell masks from regionMasks.nc file
         call mpas_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)      

         ! get values from regional stats pool
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalIceArea', RegionalIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalIceVolume', RegionalIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalFloatingIceArea', RegionalFloatingIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalFloatingIceVolume', RegionalFloatingIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalGroundedIceArea', RegionalGroundedIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'RegionalGroundedIceVolume', RegionalGroundedIceVolume)
!         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalBasalMassBal', WAIStotalBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalSfcMassBal', WAIStotalSfcMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalBasalMassBal', EAIStotalBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalSfcMassBal', EAIStotalSfcMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'FRISBasalMassBal', FRISBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'RossISBasalMassBal', RossISBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'AmeryISBasalMassBal', AmeryISBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'PineIslandThwaitesISBasalMassBal', PineIslandThwaitesISBasalMassBal)

         ! loop over cells
         do iCell = 1,nCellsSolve
!         do iGroup = 1,nRegionGroups
           do iRegion = 1,nRegions

            ! sums of ice area and volume over cells (m^2 and m^3)

            ! Regional areas and volumes
            blockSumRegionIceArea(iRegion) = blockSumRegionIceArea(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
                * areaCell(iCell) ) 

            blockSumRegionIceVolume(iRegion) = blockSumRegionIceVolume(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) &
                * areaCell(iCell) * thickness(iCell) )

!            blockSumWAISGroundedIceArea = blockSumWAISGroundedIceArea + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
!                 * areaCell(iCell) !* real(mask_is_wais(iCell),RKIND) 
!            blockSumWAISGroundedIceVolume = blockSumWAISGroundedIceVolume + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * thickness(iCell) !* real(mask_is_wais(iCell),RKIND) 
!            blockSumWAISFloatingIceArea = blockSumWAISFloatingIceArea + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) !* real(mask_is_wais(iCell),RKIND) 
!            blockSumWAISFloatingIceVolume = blockSumWAISFloatingIceVolume + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * thickness(iCell) !* real(mask_is_wais(iCell),RKIND) 

!            ! sfc and basal mass balance for WAIS, EAIS, and specific shelves (kg yr^{-1})
!            !SFP: These calculations need to be tested still; placeholder masks should return int not logical
!            blockSumWAISBasalMassBal = blockSumWAISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_wais(iCell),RKIND)
!            blockSumEAISBasalMassBal = blockSumEAISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_eais(iCell),RKIND)
!            blockSumFRISBasalMassBal = blockSumFRISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_fris(iCell),RKIND)
!            blockSumRossISBasalMassBal = blockSumRossISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_ross(iCell),RKIND)
!            blockSumAmeryISBasalMassBal = blockSumAmeryISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_amery(iCell),RKIND)
!            blockSumPineIslandThwaitesISBasalMassBal = blockSumPineIslandThwaitesISBasalMassBal + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr !* real(mask_is_pigthwaites(iCell),RKIND)

!            blockSumSfcMassBal = blockSumSfcMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * sfcMassBal(iCell) * scyr


           end do ! end loop over regions
!         end do ! end loop over groups 
         end do ! end loop over cells

         ! deallocate storage for sums over blocks  
         deallocate(blockSumRegionIceArea); deallocate(blockSumRegionIceVolume)

         block => block % next

      end do    ! end loop over blocks

      ! compute sums over all procs
!      call mpas_dmpar_sum_real(dminfo, blockSumRegionIceArea, RegionalIceArea)

!      call mpas_dmpar_sum_real(dminfo, blockSumWAISIceVolume, WAISIceVolume)
!      call mpas_dmpar_sum_real(dminfo, blockSumWAISGroundedIceArea, WAISgroundedIceArea)
!      call mpas_dmpar_sum_real(dminfo, blockSumWAISGroundedIceVolume, WAISgroundedIceVolume)
!      call mpas_dmpar_sum_real(dminfo, blockSumWAISFloatingIceArea, WAISfloatingIceArea)
!      call mpas_dmpar_sum_real(dminfo, blockSumWAISFloatingIceVolume, WAISfloatingIceVolume)

!      call mpas_dmpar_sum_real(dminfo, blockSumSfcMassBal, totalSfcMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumWAISBasalMassBal, WAIStotalBasalMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumEAISBasalMassBal, EAIStotalBasalMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumFRISBasalMassBal, FRISBasalMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumRossISBasalMassBal, RossISBasalMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumAmeryISBasalMassBal, AmeryISBasalMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumPineIslandThwaitesISBasalMassBal, PineIslandThwaitesISBasalMassBal) 


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
