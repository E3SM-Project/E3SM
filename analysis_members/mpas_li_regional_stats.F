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

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  thickness 
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBal  
      real (kind=RKIND), dimension(:), pointer ::  basalMassBal

      integer, dimension(:), pointer :: cellMask 
      integer, pointer :: nCellsSolve
      integer :: k, iCell

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), pointer ::  WAIStotalIceArea, WAIStotalIceVolume
      real (kind=RKIND), pointer ::  EAIStotalIceArea, EAIStotalIceVolume
      real (kind=RKIND), pointer ::  WAISgroundedIceArea, WAISgroundedIceVolume
      real (kind=RKIND), pointer ::  EAISgroundedIceArea, EAISgroundedIceVolume
      real (kind=RKIND), pointer ::  WAISfloatingIceArea, WAISfloatingIceVolume
      real (kind=RKIND), pointer ::  EAISfloatingIceArea, EAISfloatingIceVolume
      real (kind=RKIND), pointer ::  WAIStotalBasalMassBal!, WAIStotalSfcMassBal 
      real (kind=RKIND), pointer ::  EAIStotalBasalMassBal!, EAIStotalSfcMassBal 

      ! scalar sums over blocks
      real (kind=RKIND) ::  blockSumWAISIceArea, blockSumWAISIceVolume
      real (kind=RKIND) ::  blockSumEAISIceArea, blockSumEAISIceVolume
      real (kind=RKIND) ::  blockSumWAISGroundedIceArea, blockSumWAISGroundedIceVolume
      real (kind=RKIND) ::  blockSumEAISGroundedIceArea, blockSumEAISGroundedIceVolume
      real (kind=RKIND) ::  blockSumWAISFloatingIceArea, blockSumWAISFloatingIceVolume
      real (kind=RKIND) ::  blockSumEAISFloatingIceArea, blockSumEAISFloatingIceVolume
      real (kind=RKIND) ::  blockSumWAISBasalMassBal!, blockSumWAISSfcMassBal
      real (kind=RKIND) ::  blockSumEAISBasalMassBal!, blockSumEASISfcMassBal

      ! local parameters
      real (kind=RKIND), parameter :: scyr = 31536000.0_RKIND  ! seconds per 365-day year

      err = 0

      dminfo = domain % dminfo

      ! initialize sums over blocks to 0
      blockSumWAISIceArea = 0.0_RKIND; blockSumWAISIceVolume = 0.0_RKIND
      blockSumEAISIceArea = 0.0_RKIND; blockSumEAISIceVolume = 0.0_RKIND
      blockSumWAISGroundedIceArea = 0.0_RKIND; blockSumWAISGroundedIceVolume = 0.0_RKIND
      blockSumEAISGroundedIceArea = 0.0_RKIND; blockSumWAISGroundedIceVolume = 0.0_RKIND
      blockSumWAISFloatingIceArea = 0.0_RKIND; blockSumWAISFloatingIceVolume = 0.0_RKIND
      blockSumEAISFloatingIceArea = 0.0_RKIND; blockSumEAISFloatingIceVolume = 0.0_RKIND
      blockSumWAISBasalMassBal = 0.0_RKIND!; blockSumWAISSfcMassBal = 0.0_RKIND
      blockSumEAISBasalMassBal = 0.0_RKIND!; blockSumEAISSfcMassBal = 0.0_RKIND

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

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

         ! get values from regional stats pool
         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalIceArea', WAIStotalIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalIceVolume', WAIStotalIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'WAISfloatingIceArea', WAISfloatingIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'WAISfloatingIceVolume', WAISfloatingIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'WAISgroundedIceArea', WAISgroundedIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'WAISgroundedIceVolume', WAISgroundedIceVolume)

         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalIceArea', EAIStotalIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalIceVolume', EAIStotalIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAISfloatingIceArea', EAISfloatingIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAISfloatingIceVolume', EAISfloatingIceVolume)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAISgroundedIceArea', EAISgroundedIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAISgroundedIceVolume', EAISgroundedIceVolume)


         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalBasalMassBal', WAIStotalBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'WAIStotalSfcMassBal', WAIStotalSfcMassBal)
         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalBasalMassBal', EAIStotalBasalMassBal)
!         call mpas_pool_get_array(regionalStatsAMPool, 'EAIStotalSfcMassBal', EAIStotalSfcMassBal)

         ! loop over cells
         do iCell = 1,nCellsSolve

!            ! sums of ice area and volume over cells (m^2 and m^3)
!            blockSumIceArea = blockSumIceArea + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) 
!            blockSumIceVolume = blockSumIceVolume + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * thickness(iCell) 
!
!            blockSumGroundedIceArea = blockSumGroundedIceArea + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
!                 * areaCell(iCell) 
!            blockSumGroundedIceVolume = blockSumGroundedIceVolume + real(li_mask_is_grounded_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * thickness(iCell) 
!
!            blockSumFloatingIceArea = blockSumFloatingIceArea + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) 
!            blockSumFloatingIceVolume = blockSumFloatingIceVolume + real(li_mask_is_floating_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * thickness(iCell) 

            ! sfc and basal mass balance (kg yr^{-1})
            !SFP: These calculations need to be tested still
!            blockSumSfcMassBal = blockSumSfcMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * sfcMassBal(iCell) * scyr
!            blockSumBasalMassBal = blockSumBasalMassBal + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
!                * areaCell(iCell) * basalMassBal(iCell) * scyr

         end do ! end loop over cells

         block => block % next

      end do    ! end loop over blocks

      ! compute sums over all procs
!      call mpas_dmpar_sum_real(dminfo, blockSumIceArea, totalIceArea)
!      call mpas_dmpar_sum_real(dminfo, blockSumIceVolume, totalIceVolume)
!      call mpas_dmpar_sum_real(dminfo, blockSumGroundedIceArea, groundedIceArea)
!      call mpas_dmpar_sum_real(dminfo, blockSumGroundedIceVolume, groundedIceVolume)
!      call mpas_dmpar_sum_real(dminfo, blockSumFloatingIceArea, floatingIceArea)
!      call mpas_dmpar_sum_real(dminfo, blockSumFloatingIceVolume, floatingIceVolume)
!      call mpas_dmpar_sum_real(dminfo, blockSumSfcMassBal, totalSfcMassBal) 
!      call mpas_dmpar_sum_real(dminfo, blockSumBasalMassBal, totalBasalMassBal) 

      ! find min, max, mean thickness over all procs
!      call mpas_dmpar_min_real(dminfo, blockThickMin, iceThicknessMin)
!      call mpas_dmpar_max_real(dminfo, blockThickMax, iceThicknessMax)

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
