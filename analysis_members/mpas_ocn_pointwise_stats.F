! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_pointwise_stats
!
!> \brief MPAS ocean analysis mode member: pointwise_stats
!> \author Mark Petersen
!> \date   Jan 2016
!> \details
!>  MPAS ocean analysis mode member: pointwise_stats
!>
!-----------------------------------------------------------------------

module ocn_pointwise_stats

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_init_pointwise_stats, &
             ocn_compute_pointwise_stats, &
             ocn_restart_pointwise_stats, &
             ocn_finalize_pointwise_stats

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: nLocalCellPoints

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_pointwise_stats
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Jan 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_pointwise_stats(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

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

      type (mpas_pool_type), pointer :: pointwiseStatsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: pointwiseStatsAM

      integer, pointer :: nCells, nCellsSolve, nPoints, nVertLevels
      integer :: iCell, iPoint, i
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:), pointer :: indexToCellID,pointCellGlobalID, pointCellLocalID, pointMaxLevelCell, &
           indexToPointCellLocalID
      integer, dimension(:), allocatable :: pointMaxLevelCellLocal

      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nPoints', nPoints)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      allocate(pointMaxLevelCellLocal(nPoints))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_subpool(block % structs, 'pointwiseStatsAM', pointwiseStatsAMPool)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'indexToCellID',indexToCellID)
         call mpas_pool_get_array(pointwiseStatsAMPool,'pointCellGlobalID',pointCellGlobalID)
         call mpas_pool_get_array(pointwiseStatsAMPool,'pointMaxLevelCell',pointMaxLevelCell)
         call mpas_pool_get_array(pointwiseStatsAMPool,'pointCellLocalID',pointCellLocalID)
         call mpas_pool_get_array(pointwiseStatsAMPool,'indexToPointCellLocalID',indexToPointCellLocalID)

         pointCellLocalID = nCells + 1
         pointMaxLevelCellLocal = 0
         indexToPointCellLocalID = 0
         
         ! Initialize index arrays to record pointwise data.
         i = 0
         do iCell = 1,nCellsSolve
            do iPoint = 1,nPoints
               if (indexToCellID(iCell).eq.pointCellGlobalID(iPoint)) then
                  i = i + 1
                  indexToPointCellLocalID(i) = iPoint
                  pointCellLocalID(iPoint) = iCell
                  pointMaxLevelCellLocal(iPoint) = maxLevelCell(iCell)
               endif
            end do
         end do
         nLocalCellPoints = i

         block => block % next
      end do

      ! Collect data over processors by summing array.  Only owned cells are nonzero.
      call mpas_dmpar_sum_int_array(dminfo, nPoints, pointMaxLevelCellLocal, pointMaxLevelCell)

      deallocate(pointMaxLevelCellLocal)

   end subroutine ocn_init_pointwise_stats!}}}

!***********************************************************************
!
!  routine ocn_compute_pointwise_stats
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Jan 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_pointwise_stats(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

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

      type (mpas_pool_type), pointer :: pointwiseStatsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: pointwiseStatsAM

      integer, pointer :: nPoints, nVertLevels
      integer :: iPoint, i
      integer, dimension(:), pointer :: pointCellGlobalID, pointCellLocalID, indexToPointCellLocalID

      ! Doug replace these when you iterate over stream variables:
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessPointValue, layerThickness
      real (kind=RKIND), dimension(:,:), allocatable :: sumLayerThicknessPointValue

      err = 0

      dminfo = domain % dminfo

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nPoints', nPoints)
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)
      allocate(sumLayerThicknessPointValue(nVertLevels,nPoints))

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

         call mpas_pool_get_subpool(block % structs, 'pointwiseStatsAM', pointwiseStatsAMPool)
         call mpas_pool_get_dimension(block % dimensions, 'nPoints', nPoints)
         call mpas_pool_get_array(pointwiseStatsAMPool,'pointCellGlobalID',pointCellGlobalID)
         call mpas_pool_get_array(pointwiseStatsAMPool,'pointCellLocalID',pointCellLocalID)
         call mpas_pool_get_array(pointwiseStatsAMPool,'indexToPointCellLocalID',indexToPointCellLocalID)

         ! Doug: replace these when you iterate over stream variables:
         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(pointwiseStatsAMPool,'layerThicknessPointValue',layerThicknessPointValue)

         ! Doug: loop over names in stream rather than layerThickness
         layerThicknessPointValue = 0.0_RKIND
         do i = 1,nLocalCellPoints
            layerThicknessPointValue(:,indexToPointCellLocalID(i)) = &
                 layerThickness(:,pointCellLocalID(indexToPointCellLocalID(i)))
         end do
         !end do

         block => block % next
      end do

      ! Collect data over processors by summing array.  Only owned cells are nonzero.
      call mpas_dmpar_sum_real_array(dminfo, nVertLevels*nPoints, layerThicknessPointValue, sumLayerThicknessPointValue)

      layerThicknessPointValue = sumLayerThicknessPointValue

      deallocate(sumLayerThicknessPointValue)

    end subroutine ocn_compute_pointwise_stats!}}}

!***********************************************************************
!
!  routine ocn_restart_pointwise_stats
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Jan 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_pointwise_stats(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

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

   end subroutine ocn_restart_pointwise_stats!}}}

!***********************************************************************
!
!  routine ocn_finalize_pointwise_stats
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    Jan 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_pointwise_stats(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

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

   end subroutine ocn_finalize_pointwise_stats!}}}

end module ocn_pointwise_stats

! vim: foldmethod=marker
