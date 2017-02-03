! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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

      ! arrays, vars needed from other pools for calculations here
      real (kind=RKIND), pointer ::  config_ice_density
      real (kind=RKIND), pointer ::  deltat
      real (kind=RKIND), dimension(:), pointer ::  areaCell
      real (kind=RKIND), dimension(:), pointer ::  thickness
      real (kind=RKIND), dimension(:), pointer ::  sfcMassBal
      real (kind=RKIND), dimension(:), pointer ::  basalMassBal
      real (kind=RKIND), dimension(:), pointer ::  calvingThickness

      integer, dimension(:), pointer :: cellMask
      integer, pointer :: nCellsSolve
      integer :: k, iCell

      ! scalars to be calculated here from global sums
      real (kind=RKIND), pointer ::  totalIceArea, totalIceVolume
      real (kind=RKIND), pointer ::  groundedIceArea, groundedIceVolume
      real (kind=RKIND), pointer ::  floatingIceArea, floatingIceVolume
      real (kind=RKIND), pointer ::  iceThicknessMax, iceThicknessMin, iceThicknessMean
      real (kind=RKIND), pointer ::  totalSfcMassBal, totalBasalMassBal
      real (kind=RKIND), pointer ::  totalCalvingFlux

      ! scalar sums over blocks
      real (kind=RKIND) ::  blockSumIceArea, blockSumIceVolume
      real (kind=RKIND) ::  blockSumGroundedIceArea, blockSumGroundedIceVolume
      real (kind=RKIND) ::  blockSumFloatingIceArea, blockSumFloatingIceVolume
      real (kind=RKIND) ::  blockThickMin, blockThickMax
      real (kind=RKIND) ::  blockSumSfcMassBal, blockSumBasalMassBal
      real (kind=RKIND) ::  blockSumCalvingFlux

      ! local parameters
      real (kind=RKIND), parameter :: scyr = 31536000.0_RKIND  ! seconds per 365-day year

      err = 0

      dminfo = domain % dminfo

      ! initialize sums over blocks to 0
      blockSumIceArea = 0.0_RKIND
      blockSumIceVolume = 0.0_RKIND
      blockSumGroundedIceArea = 0.0_RKIND
      blockSumGroundedIceVolume = 0.0_RKIND
      blockSumFloatingIceArea = 0.0_RKIND
      blockSumFloatingIceVolume = 0.0_RKIND
      blockSumSfcMassBal = 0.0_RKIND
      blockSumBasalMassBal = 0.0_RKIND
      blockSumCalvingFlux = 0.0_RKIND

      ! initialize max, min, mean values to 0
      blockThickMin = 100000.0_RKIND
      blockThickMax = 0.0_RKIND

      ! loop over blocks
      block => domain % blocklist
      do while (associated(block))

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'deltat', deltat)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)
         call mpas_pool_get_array(geometryPool, 'calvingThickness', calvingThickness)

         ! get values from global stats pool
         call mpas_pool_get_array(globalStatsAMPool, 'totalIceArea', totalIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'totalIceVolume', totalIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceArea', floatingIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceVolume', floatingIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceArea', groundedIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceVolume', groundedIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMax', iceThicknessMax)
         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMin', iceThicknessMin)
         call mpas_pool_get_array(globalStatsAMPool, 'iceThicknessMean', iceThicknessMean)
         call mpas_pool_get_array(globalStatsAMPool, 'totalSfcMassBal', totalSfcMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalBasalMassBal', totalBasalMassBal)
         call mpas_pool_get_array(globalStatsAMPool, 'totalCalvingFlux', totalCalvingFlux)

         ! loop over cells
         do iCell = 1,nCellsSolve

            ! sums of ice area and volume over cells (m^2 and m^3)
            blockSumIceArea = blockSumIceArea + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell)
            blockSumIceVolume = blockSumIceVolume + real(li_mask_is_ice_int(cellMask(iCell)),RKIND) &
                * areaCell(iCell) * thickness(iCell)

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

            ! mass lass due do calving (kg yr^{-1})
            blockSumCalvingFlux = blockSumCalvingFlux + calvingThickness(iCell) * &
               areaCell(iCell) * config_ice_density / ( deltat / scyr )

         end do ! end loop over cells

         block => block % next

      end do    ! end loop over blocks

      ! compute sums over all procs
      call mpas_dmpar_sum_real(dminfo, blockSumIceArea, totalIceArea)
      call mpas_dmpar_sum_real(dminfo, blockSumIceVolume, totalIceVolume)
      call mpas_dmpar_sum_real(dminfo, blockSumGroundedIceArea, groundedIceArea)
      call mpas_dmpar_sum_real(dminfo, blockSumGroundedIceVolume, groundedIceVolume)
      call mpas_dmpar_sum_real(dminfo, blockSumFloatingIceArea, floatingIceArea)
      call mpas_dmpar_sum_real(dminfo, blockSumFloatingIceVolume, floatingIceVolume)
      call mpas_dmpar_sum_real(dminfo, blockSumSfcMassBal, totalSfcMassBal)
      call mpas_dmpar_sum_real(dminfo, blockSumBasalMassBal, totalBasalMassBal)
      call mpas_dmpar_sum_real(dminfo, blockSumCalvingFlux, totalCalvingFlux)

      ! find min, max, mean thickness over all procs
      call mpas_dmpar_min_real(dminfo, blockThickMin, iceThicknessMin)
      call mpas_dmpar_max_real(dminfo, blockThickMax, iceThicknessMax)
      iceThicknessMean = totalIceVolume / totalIceArea

      ! Even though some variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)

         ! assignment of final globalStatsAM variables could occur here.

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
