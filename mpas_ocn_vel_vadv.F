! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_vadv
!
!> \brief MPAS ocean vertical advection 
!> \author Doug Jacobsen
!> \date   15 September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencies for vertical advection.
!>
!
!-----------------------------------------------------------------------

module ocn_vel_vadv

   use mpas_grid_types
   use mpas_configure

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

   public :: ocn_vel_vadv_tend, &
             ocn_vel_vadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: velVadvOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_vadv_tend
!
!> \brief   Computes tendency term for vertical advection
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details 
!>  This routine computes the vertical advection tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_tend(grid, u, layerThicknessEdge, vertTransportVelocityTop, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         u     !< Input: Horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge,&!< Input: thickness at edge
         vertTransportVelocityTop  !< Input: Vertical velocity on top layer

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

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

      integer :: iEdge, nEdgesSolve, cell1, cell2, k
      integer :: nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND) :: vertTransportVelocityTopEdge
      real (kind=RKIND), dimension(:), allocatable :: w_dudzTopEdge

      if(.not.velVadvOn) return

      err = 0

      nVertLevels = grid % nVertLevels
      nEdgesSolve = grid % nEdgesSolve
      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgeMask => grid % edgeMask % array

      allocate(w_dudzTopEdge(nVertLevels+1))
      w_dudzTopEdge = 0.0
      do iEdge=1,nEdgesSolve
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        do k=2,maxLevelEdgeTop(iEdge)
          ! Average w from cell center to edge
          vertTransportVelocityTopEdge = 0.5*(vertTransportVelocityTop(k,cell1)+vertTransportVelocityTop(k,cell2))

          ! compute dudz at vertical interface with first order derivative.
          w_dudzTopEdge(k) = vertTransportVelocityTopEdge * (u(k-1,iEdge)-u(k,iEdge)) &
                       / (0.5*(layerThicknessEdge(k-1,iEdge) + layerThicknessEdge(k,iEdge)))
        end do
        w_dudzTopEdge(maxLevelEdgeTop(iEdge)+1) = 0.0
        ! Average w*du/dz from vertical interface to vertical middle of cell
        do k=1,maxLevelEdgeTop(iEdge)

          tend(k,iEdge) = tend(k,iEdge) - edgeMask(k, iEdge) * 0.5 * (w_dudzTopEdge(k) + w_dudzTopEdge(k+1))
        enddo
      enddo
      deallocate(w_dudzTopEdge)

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_vadv_init
!
!> \brief   Initializes ocean momentum vertical advection
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  vertical velocity advection in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_vadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! Output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0
      velVadvOn = .false.

      if (config_vert_coord_movement.ne.'isopycnal') then
          velVadvOn = .true.
      end if

      if(config_disable_vel_vadv) velVadvOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_vadv_init!}}}

!***********************************************************************

end module ocn_vel_vadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
