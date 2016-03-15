! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_debug_diagnostics
!
!> \brief MPAS ocean analysis mode member: debug_diagnostics
!> \author Mark Petersen
!> \date   March 2016
!> \details
!>  MPAS ocean analysis mode member: debug_diagnostics
!>  Compute diagnostics used for debugging.
!>
!-----------------------------------------------------------------------

module ocn_debug_diagnostics

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

   public :: ocn_init_debug_diagnostics, &
             ocn_compute_debug_diagnostics, &
             ocn_restart_debug_diagnostics, &
             ocn_finalize_debug_diagnostics

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_debug_diagnostics
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_debug_diagnostics(domain, err)!{{{

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

   end subroutine ocn_init_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_compute_debug_diagnostics
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_debug_diagnostics(domain, timeLevel, err)!{{{

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

      type (mpas_pool_type), pointer :: debugDiagnosticsAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      !type (mpas_pool_type), pointer :: debugDiagnosticsAM

      ! Here are some example variables which may be needed for your analysis member
      integer :: iEdge, c1, c2, k
      integer, pointer :: nEdges
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) :: dzVert1, dzVert2, dzEdgeK, dzEdgeKp1, rx1, localMaxRx1
      real (kind=RKIND), pointer :: globalRx1Max
      real (kind=RKIND), dimension(:), pointer :: rx1MaxCell
      real (kind=RKIND), dimension(:,:), pointer :: zMid

      err = 0

      localMaxRx1 = 0.0_RKIND

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         !call mpas_pool_get_subpool(block % structs, 'debugDiagnosticsAM', debugDiagnosticsAMPool)

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

         !-----------------------------------------------------------------
         !
         ! Compute Haney number, rx1
         !
         !-----------------------------------------------------------------

         call mpas_pool_get_array(diagnosticsPool, 'rx1MaxCell', rx1MaxCell)
         call mpas_pool_get_array(diagnosticsPool, 'globalRx1Max', globalRx1Max)

        ! These could be included for edge or cell fields with depth:
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1Edge', rx1Edge)
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1Cell', rx1Cell)
        ! call mpas_pool_get_array(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge)
        ! rx1Edge(:,:) = 0.0_RKIND
        ! rx1Cell(:,:) = 0.0_RKIND
        ! rx1MaxEdge(:) = 0.0_RKIND

        rx1MaxCell(:) = 0.0_RKIND
        do iEdge = 1,nEdges
          c1 = cellsOnEdge(1,iEdge)
          c2 = cellsOnEdge(2,iEdge)
          do k = 1,maxLevelEdgeTop(iEdge)-1
            dzVert1 = zMid(k,c1)-zMid(k+1,c1)
            dzVert2 = zMid(k,c2)-zMid(k+1,c2)
            dzEdgeK = zMid(k,c2)-zMid(k,c1)
            dzEdgeKp1 = zMid(k+1,c2)-zMid(k+1,c1)

            rx1 = abs(dzEdgeK+dzEdgeKp1)/(dzVert1+dzVert2)

            rx1MaxCell(c1) = max(rx1MaxCell(c1),rx1)
            rx1MaxCell(c2) = max(rx1MaxCell(c2),rx1)

            ! These could be included for edge or cell fields with depth:
            ! rx1Edge(k,iEdge) = rx1
            ! rx1Cell(k,c1) = max(rx1Cell(k,c1),rx1)
            ! rx1Cell(k,c2) = max(rx1Cell(k,c2),rx1)
            ! rx1MaxEdge(iEdge) = max(rx1MaxEdge(iEdge),rx1)
          end do
        end do

        localMaxRx1 = max(localMaxRx1,maxval(rx1MaxCell))

        block => block % next
      end do
      call mpas_dmpar_max_real(dminfo, localMaxRx1, globalRx1Max)
      write (stdoutUnit,'(a, es10.2)') ' global max of rx1:', globalRx1Max

   end subroutine ocn_compute_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_restart_debug_diagnostics
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_debug_diagnostics(domain, err)!{{{

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

   end subroutine ocn_restart_debug_diagnostics!}}}

!***********************************************************************
!
!  routine ocn_finalize_debug_diagnostics
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Mark Petersen
!> \date    March 2016
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_debug_diagnostics(domain, err)!{{{

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

   end subroutine ocn_finalize_debug_diagnostics!}}}

end module ocn_debug_diagnostics

! vim: foldmethod=marker
