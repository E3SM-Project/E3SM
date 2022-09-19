! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_explicit_bottom_drag
!
!> \brief MPAS ocean explicit bottom drag
!> \author Mark Petersen
!> \date   August 2017
!> \details
!>  This module contains the routine for computing
!>  tendencies from explicit bottom drag.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_explicit_bottom_drag

   use mpas_timer

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_forcing

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

   public :: ocn_vel_forcing_explicit_bottom_drag_tend, &
             ocn_vel_forcing_explicit_bottom_drag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: &
      explicitBottomDragOff  ! on/off switch for explicit bottom drag

   real (kind=RKIND) :: &
      dragCoeff          ! drag coefficient

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_bottom_drag_tend
!
!> \brief   Computes tendency term from explicit bottom drag
!> \author  Mark Petersen
!> \date    15 August 2017
!> \details
!>  This routine computes the explicit bottom drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_bottom_drag_tend(normVelocity, &
                                 KECell, layerThickEdgeMean, tend, err) !{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normVelocity,      &!< [in] normal velocity
         KECell,            &!< [in] kinetic energy at cell
         layerThickEdgeMean  !< [in] mean layer thickness at edge

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend            !< [inout] accumulated velocity tendency

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         iEdge,      &! loop index for edge loop
         k,          &! vertical index of lowest active layer at edge
         cell1, cell2 ! neighbor cell addresses across edge

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and return if not turned on
      !*** Otherwise start timer

      err = 0
      if (explicitBottomDragOff) return
      call mpas_timer_start('vel explicit bottom drag')

      ! Explicit bottom drag term:
      ! du/dt = ... - c |u| u / h
      ! appied to bottom layer only.
      ! This term comes from the bottom boundary condition in the vertical
      ! momentum mixing, and is explicit if both |u| and u are chosen to be at
      ! time level n.

#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(cellsOnEdge, maxLevelEdgeTop, KECell, &
      !$acc            tend, normVelocity, layerThickEdgeMean) &
      !$acc    private(k, cell1, cell2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k, cell1, cell2)
#endif
      do iEdge = 1, nEdgesOwned
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        k =  maxLevelEdgeTop(iEdge)

        if (k > 0) then
           tend(k,iEdge) = tend(k,iEdge) - dragCoeff* &
                           sqrt(KECell(k,cell1) + KECell(k,cell2))* &
                           normVelocity(k,iEdge)/layerThickEdgeMean(k,iEdge)
        end if

      enddo
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop('vel explicit bottom drag')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_bottom_drag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_explicit_bottom_drag_init
!
!> \brief   Initializes ocean explicit bottom drag forcing
!> \author  Mark Petersen
!> \date    August 2017
!> \details
!>  This routine initializes quantities related to explicit bottom drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_explicit_bottom_drag_init(err)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< [out] error flag

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize return error code and set module defaults
      err = 0

      explicitBottomDragOff  = .true.
      dragCoeff = 0.0_RKIND

      !*** Reset values based on input model configuration

      if (config_use_explicit_bottom_drag) then
         explicitBottomDragOff = .false.
         dragCoeff = config_explicit_bottom_drag_coeff
      endif

      if (config_disable_vel_explicit_bottom_drag) then
         explicitBottomDragOff = .true.
         dragCoeff = 0.0_RKIND
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_bottom_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_explicit_bottom_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
