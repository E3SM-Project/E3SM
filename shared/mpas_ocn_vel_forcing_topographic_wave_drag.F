! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_forcing_topographic_wave_drag
!
!> \brief MPAS ocean topographic wave drag
!> \author Nairita Pal
!> \date   Oct 2020
!> \details
!>  This module contains the routine for computing
!>  tendencies from topographic wave drag.
!
!-----------------------------------------------------------------------

module ocn_vel_forcing_topographic_wave_drag

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer

   use ocn_constants
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

   public :: ocn_vel_forcing_topographic_wave_drag_tend, &
             ocn_vel_forcing_topographic_wave_drag_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: topographicWaveDragOn
   real (kind=RKIND) :: topographicWaveDragCoeff

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vel_forcing_topographic_wave_drag_tend
!
!> \brief   Computes tendency term from topographic wave drag
!> \author  Nairita Pal
!> \date    15 Oct 2020
!> \details
!>  This routine computes the topographic wave drag tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_topographic_wave_drag_tend(normalVelocity, & !{{{
                     kineticEnergyCell, layerThicknessEdge, topographic_wave_drag_rinv, &
                     tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity,       &!< Input: velocity
         kineticEnergyCell,    &!< Input: kinetic energy at cell
         layerThicknessEdge   !< Input: thickness at edge

        real (kind=RKIND), dimension(:), intent(in) :: topographic_wave_drag_rinv

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

      integer :: iEdge, k, cell1, cell2

      err = 0
      if ( .not. topographicWaveDragOn ) return

      call mpas_timer_start('vel topographic drag')

      !$omp do schedule(runtime) private(k, cell1, cell2)
      do iEdge = 1, nEdgesOwned
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        k =  maxLevelEdgeTop(iEdge)
        ! topographic wave drag term:
        ! du/dt = ... (1.0/topographic_wave_drag_rinv) * u
        ! topographic_wave_drag_rinv is the e-folding time which is read from a file
        tend(k,iEdge) = tend(k,iEdge) - topographicWaveDragCoeff &
           /(0.5_RKIND*(topographic_wave_drag_rinv(cell1)+topographic_wave_drag_rinv(cell2)))*normalVelocity(k,iEdge)

      enddo
      !$omp end do

      call mpas_timer_stop('vel topographic drag')

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_tend!}}}

!***********************************************************************
!
!  routine ocn_vel_forcing_topographic_wave_drag_init
!
!> \brief   Initializes ocean topographic wave drag forcing
!> \author  Nairita Pal
!> \date    Oct 2020
!> \details
!>  This routine initializes quantities related to topographic wave  drag
!>  in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_forcing_topographic_wave_drag_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_vel_topographic_wave_drag
      logical, pointer :: config_use_topographic_wave_drag
      real (kind=RKIND), pointer :: config_topographic_wave_drag_coeff

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_topographic_wave_drag', config_use_topographic_wave_drag)
      call mpas_pool_get_config(ocnConfigs, 'config_topographic_wave_drag_coeff', config_topographic_wave_drag_coeff)
      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_topographic_wave_drag', config_disable_vel_topographic_wave_drag)

      topographicWaveDragCoeff = 0.0_RKIND

      if (config_use_topographic_wave_drag) then
          topographicWaveDragOn = .true.
          topographicWaveDragCoeff = config_topographic_wave_drag_coeff
      endif

      if (config_disable_vel_topographic_wave_drag) topographicWaveDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_topographic_wave_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
