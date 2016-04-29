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
      type (mpas_pool_type), pointer :: regionsPool                !!SFP: added

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
      integer, pointer :: nCellsSolve, nRegions, nRegionGroups !, maxRegionsInGroup  !! SFP: is last one needed here?
      integer :: k, iCell, iRegion, iGroup

      ! scalars to be calculated here from regional sums
      real (kind=RKIND), dimension(:), pointer ::  regionalIceArea, regionalIceVolume

      ! storage for sums over blocks
      real (kind=RKIND), dimension(:), allocatable ::  blockSumRegionIceArea, blockSumRegionIceVolume

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
      allocate(blockSumRegionIceArea(nRegions)) 
      allocate(blockSumRegionIceVolume(nRegions))

      blockSumRegionIceArea = 0.0_RKIND; blockSumRegionIceVolume = 0.0_RKIND

!      ! loop over blocks
!      block => domain % blocklist
!      do while (associated(block))

         ! get structs from pools
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'regions', regionsPool)
         call mpas_pool_get_subpool(block % structs, 'regionalStatsAM', regionalStatsAMPool)

         ! get values and arrays from standard pools
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         print *, 'nCellsSolve = ', nCellsSolve

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
!         call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)
         call mpas_pool_get_array(geometryPool, 'basalMassBal', basalMassBal)

         !!SFP: get region cell masks from regionMasks.nc file
         call mpas_pool_get_array(regionsPool, 'regionCellMasks', regionCellMasks)      
         print *, 'sum( regionCellMasks(1,:) ) = ', sum( regionCellMasks(1,:) )
         print *, 'sum( regionCellMasks(2,:) ) = ', sum( regionCellMasks(2,:) )
         print *, 'size of region cell mask = ', size( regionCellMasks, 1), size( regionCellMasks, 2)

         ! get values from regional stats pool
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceArea', regionalIceArea)
         call mpas_pool_get_array(regionalStatsAMPool, 'regionalIceVolume', regionalIceVolume)

         ! loop over cells
         do iCell = 1,nCellsSolve
!         do iGroup = 1,nRegionGroups
           do iRegion = 1,nRegions

            ! sums of ice area and volume over cells (m^2 and m^3)

            ! Regional areas and volumes
            blockSumRegionIceArea(iRegion) = blockSumRegionIceArea(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) & 
               * areaCell(iCell) ) 

!!! SFP: Note that this calc. needs to only occur where there is ice ... needs to include cellMask also? Check globalStatsAM for this too !!!
            blockSumRegionIceVolume(iRegion) = blockSumRegionIceVolume(iRegion) + ( real(regionCellMasks(iRegion,iCell),RKIND) & 
               * areaCell(iCell) * thickness(iCell) )

            print *, 'nRegions = ', nRegions
            print *, 'iRegion = ', iRegion
            print *, 'iCell=', iCell
            print *, 'thickness(iCell) = ', thickness(iCell)
            print *, 'area(iCell) = ', areaCell(iCell)
            print *, 'regionCellMasks(iRegion,iCell) = ', regionCellMasks(iRegion,iCell)
            print *, 'blockSumRegionIceArea = ', blockSumRegionIceArea(iRegion)
            print *, 'blockSumRegionIceVolume = ', blockSumRegionIceVolume(iRegion)
            print *, ' '

           end do ! end loop over regions
!         end do ! end loop over groups 
         end do ! end loop over cells


         block => block % next

      end do    ! end loop over blocks

      !!SFP: debug
!            print *, 'regionCellMasks(1,:) = ', regionCellMasks(1,:)
!            print *, 'regionCellMasks(2,:) = ', regionCellMasks(2,:)


      ! compute sums over all procs
!      do iRegion = 1,nRegions

! SFP: debug
!       regionalIceArea(iRegion) = blockSumRegionIceArea(iRegion)
!       regionalIceVolume(iRegion) = blockSumRegionIceVolume(iRegion)

! SFP: sample calls from ocean analysis member - 2nd arg. is array size? !!!
      !call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceArea(iRegion), regionalIceArea(iRegion) )
      !call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceVolume(iRegion), regionalIceVolume(iRegion) )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceArea, regionalIceArea )
      call mpas_dmpar_sum_real_array(dminfo, nRegions, blockSumRegionIceVolume, regionalIceVolume )

!      end do ! end loop over regions


      ! deallocate storage for sums over blocks  
      deallocate(blockSumRegionIceArea); deallocate(blockSumRegionIceVolume)


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
