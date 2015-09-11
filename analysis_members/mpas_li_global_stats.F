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

      ! Here are some example variables which may be needed for your analysis member
!      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_tracers
      integer, pointer :: nCellsSolve
!      integer :: iTracer, k, iCell
      integer :: k, iCell
!      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

!      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge
      real (kind=RKIND), dimension(:), pointer ::  areaCell

      real (kind=RKIND), dimension(:), pointer ::  thickness 
      integer, dimension(:), pointer :: cellMask 

      ! simple 1 or 0 masks to be used here for calc. global sums over floating or grounded ice
      integer, dimension(:), pointer :: iceMask 
      integer, dimension(:), pointer :: groundedMask 
      integer, dimension(:), pointer :: floatingMask 

      ! scalars to be calculated here from global sums
      real (kind=RKIND), pointer ::  totalIceArea 
      real (kind=RKIND), pointer ::  totalIceVolume
      real (kind=RKIND), pointer ::  groundedIceArea 
      real (kind=RKIND), pointer ::  groundedIceVolume
      real (kind=RKIND), pointer ::  floatingIceArea 
      real (kind=RKIND), pointer ::  floatingIceVolume

      err = 0

      dminfo = domain % dminfo

      ! initialize scalar global sums and work masks to zero
      totalIceArea = 0.0_RKIND
      totalIceVolume = 0.0_RKIND
      groundedIceArea = 0.0_RKIND
      groundedIceVolume = 0.0_RKIND
      floatingIceArea = 0.0_RKIND
      floatingIceVolume = 0.0_RKIND
      iceMask = 0
      groundedMask = 0
      floatingMask = 0     

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'globalStatsAM', globalStatsAMPool)

         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)

         ! Here are some example variables which may be needed for your analysis member
!         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

!         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
!         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
!         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
!         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
!         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
!         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
!         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
!         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)

         call mpas_pool_get_array(globalStatsAMPool, 'totalIceArea', totalIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'totalIceVolume', totalIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceArea', floatingIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'floatingIceVolume', floatingIceVolume)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceArea', groundedIceArea)
         call mpas_pool_get_array(globalStatsAMPool, 'groundedIceVolume', groundedIceVolume)

         ! populate work masks (1 and 0 based for multiplication of area and thickness fields)
         where( cellMask == 32 ); iceMask = 1; endwhere
         where( cellMask == 4 ); floatingMask = 1; endwhere
         groundedMask = iceMask - floatingMask 

         ! Computations which are functions of nCells, nEdges, or nVertices
         ! must be placed within this block loop
         ! Here are some example loops
         do iCell = 1,nCellsSolve

!            do k = 1, maxLevelCell(iCell)
!               do iTracer = 1, num_tracers
               ! computations on tracers(iTracer,k, iCell)
!               end do
!            end do

            ! calculate total ice area and volume
            totalIceArea = totalIceArea + real( iceMask(iCell), RKIND) * areaCell(iCell) 
            totalIceVolume = totalIceVolume + real( iceMask(iCell), RKIND) * areaCell(iCell) * thickness(iCell) 

            ! debugging
!            print *, 'totalIceArea=', totalIceArea
!            print *, 'totalIceVolume=', totalIceVolume

            ! calculate grounded ice area and volume
            groundedIceArea = groundedIceArea + real( groundedMask(iCell), RKIND) * areaCell(iCell) 
            groundedIceVolume = groundedIceVolume + real( groundedMask(iCell), RKIND) * areaCell(iCell) * thickness(iCell) 

            ! calculate floating ice area and volume
            floatingIceArea = floatingIceArea + real( floatingMask(iCell), RKIND) * areaCell(iCell) 
            floatingIceVolume = floatingIceVolume + real( floatingMask(iCell), RKIND) * areaCell(iCell) * thickness(iCell) 

         end do

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      ! Here are some examples.  See mpas_oac_global_stats.F for further details.
!      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
!      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
!      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

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
