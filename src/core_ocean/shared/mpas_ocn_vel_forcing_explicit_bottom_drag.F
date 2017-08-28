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

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer

   use ocn_constants
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

   logical :: explicitBottomDragOn
   real (kind=RKIND) :: explicitBottomDragCoef

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

   subroutine ocn_vel_forcing_explicit_bottom_drag_tend(meshPool, normalVelocity, & !{{{
                                     kineticEnergyCell, layerThicknessEdge, tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,       &!< Input: velocity
         kineticEnergyCell,    &!< Input: kinetic energy at cell
         layerThicknessEdge     !< Input: thickness at edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

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

      integer :: iEdge, k, cell1, cell2, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      err = 0

      if ( .not. explicitBottomDragOn ) return

      call mpas_timer_start('vel explicit bottom drag')

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      nEdges = nEdgesArray( 1 )

      !$omp do schedule(runtime) private(k, cell1, cell2)
      do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        k =  maxLevelEdgeTop(iEdge)

        ! Explicit bottom drag term:
        ! du/dt = ... - c |u| u / h
        ! appied to bottom layer only.
        ! This term comes from the bottom boundary condition in the vertical
        ! momentum mixing, and is explicit if both |u| and u are chosen to be at
        ! time level n. 

        tend(k,iEdge) = tend(k,iEdge) - explicitBottomDragCoef * & 
           sqrt(kineticEnergyCell(k,cell1) + kineticEnergyCell(k,cell2)) * normalVelocity(k,iEdge) / layerThicknessEdge(k,iEdge)

      enddo
      !$omp end do

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

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_vel_explicit_bottom_drag
      logical, pointer :: config_use_explicit_bottom_drag
      real (kind=RKIND), pointer :: config_explicit_bottom_drag_coeff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_explicit_bottom_drag', config_use_explicit_bottom_drag)
      call mpas_pool_get_config(ocnConfigs, 'config_explicit_bottom_drag_coeff', config_explicit_bottom_drag_coeff)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_explicit_bottom_drag', config_disable_vel_explicit_bottom_drag)

      explicitBottomDragCoef = 0.0_RKIND

      if (config_use_explicit_bottom_drag) then
          explicitBottomDragOn = .true.
          explicitBottomDragCoef = config_explicit_bottom_drag_coeff
      endif

      if (config_disable_vel_explicit_bottom_drag) explicitBottomDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_explicit_bottom_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_explicit_bottom_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
