! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use ocn_diagnostics_variables
   use ocn_forcing
   use ocn_config

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

   logical :: topographicWaveDragOn, tensorScheme 
   real (kind=RKIND) :: topographicWaveDragCoeff
!   real (kind=RKIND), dimension(:), allocatable :: topographicWaveDrag

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
                     tend, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity   !< Input: velocity

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

      integer :: iEdge, k

      err = 0
      if ( .not. topographicWaveDragOn ) return

      call mpas_timer_start('vel topographic wave drag')

      ! JSL: C = pi std(H)^2 N_B / L
      ! ZAE: C = gamma H (grad(H))^2 N_B bar(N) / (8 * pi^2 * omega)
      !          N(z) = N_0 exp(z/1300)
      ! Drag: rho_0 C u
      
      if ( .not. tensorScheme) then
          !$omp do schedule(runtime) private(k)
          do iEdge = 1, nEdgesOwned
            k =  maxLevelEdgeTop(iEdge)
            ! topographic wave drag term:
            ! du/dt = ... 1/rinv * u
            ! rinv is the e-folding time use in HyCOM. Here
            ! topographic_wave_drag = 1/rinv
            if (k>0) then
              temp_twd(iEdge) = topographicWaveDrag(iEdge) * normalVelocity(k,iEdge)
              tend(k,iEdge) = tend(k,iEdge) - topographicWaveDrag(iEdge) * normalVelocity(k,iEdge)
            endif
          enddo
          !$omp end do
      else
          !$omp do schedule(runtime) private(k)
          do iEdge = 1, nEdgesOwned
            k =  maxLevelEdgeTop(iEdge)
            if (k>0) then
              temp_twd(iEdge) = topographicWaveDrag(iEdge) * (normalVelocity(k,iEdge) * normalCoeff(iEdge) + tangentialVelocity(k,iEdge) * tangentialCoeff(iEdge))
              tend(k,iEdge) = tend(k,iEdge) - topographicWaveDrag(iEdge) &
                 * (normalVelocity(k,iEdge) * normalCoeff(iEdge) + tangentialVelocity(k,iEdge) * tangentialCoeff(iEdge))
            endif
          enddo
          !$omp end do
      endif

      call mpas_timer_stop('vel topographic wave drag')

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

      real (kind=RKIND) :: topoDragCutoffDepth, topoDragCutoffWidth
      real (kind=RKIND) :: omegaM2 = 1.405189e-4_RKIND
      real (kind=RKIND), dimension(:), allocatable :: topoDragTurnOff

      err = 0


      topographicWaveDragCoeff = 0.0_RKIND

      if (config_use_topographic_wave_drag) then
          topographicWaveDragOn = .true.
          topographicWaveDragCoeff = config_topographic_wave_drag_coeff
          topoDragCutoffDepth = config_topographic_wave_drag_cutoff_depth
          topoDragCutoffWidth = config_topographic_wave_drag_cutoff_width
          topoDragTurnOff = (TANH((-bottomDepthEdge-topoDragCutoffDepth)/topoDragCutoffWidth)+1)/2
          if (config_topographic_wave_drag_scheme.EQ."JSL") then
              tensorScheme = .false.
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Jayne and St. Laurent")
              call mpas_log_write("")
              topographicWaveDrag = topographicWaveDragCoeff*topoDragTurnOff*pii*(bathy_stddev**2)*topo_buoyancy_N1B/10000.0_RKIND/bottomDepthEdge
              ! C = pi std(H)^2 Nb / L, L=10,000

          else if (config_topographic_wave_drag_scheme.EQ."ZAE") then
              tensorScheme = .false.
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Zaron and Egbert")
              call mpas_log_write("")
              topographicWaveDrag = topographicWaveDragCoeff*(bed_slope_edges**2)*topo_buoyancy_N1B*topo_buoyancy_N1V / (8.0_RKIND * 1.405189e-4_RKIND * pii**2)
!              topographicWaveDrag = topographicWaveDragCoeff*topoDragTurnOff*50.0_RKIND*bottomDepthEdge*(bed_slope_edges**2)*((5.24e-3_RKIND)*1300.0_RKIND*(1 - exp(-bottomDepthEdge/1300.0_RKIND))/bottomDepthEdge)*((5.24e-3_RKIND)*exp(-bottomDepthEdge/1300.0_RKIND))/ (8 * 1.405189e-4_RKIND * pii**2)
              ! C = G*H*(grad H)^2 * Nb * Nv / (8 pi^2 omega), take M2 omega = 1.405189e-4_RKIND
!              Nv = (5.24e-3_RKIND)*1300*(1 - exp(-bottomDepth/1300))/bottomDepth
!              Nb = (5.24e-3_RKIND)*exp(-bottomDepth/1300)
          else if (config_topographic_wave_drag_scheme.EQ."LGF") then
              tensorScheme = .true.
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Local Generation Formula")
              call mpas_log_write("")
              topographicWaveDrag = topographicWaveDragCoeff*topoDragTurnOff*(sqrt((topo_buoyancy_N1B**2 - omegaM2**2)*(topo_buoyancy_N1V**2 - omegaM2**2)))/(4*pii*omegaM2)
              normalCoeff = (lonGradEdge**2)*(cos(angleEdge)**2) + (latGradEdge**2)*(sin(angleEdge)**2) - 2*latGradEdge*lonGradEdge*sin(angleEdge)*(cos(angleEdge))
              tangentialCoeff = sin(angleEdge)*cos(angleEdge)*(lonGradEdge**2 - latGradEdge**2) + latGradEdge*lonGradEdge*(cos(angleEdge)**2 - sin(angleEdge)**2)

          else 
              call mpas_log_write("")

              call mpas_log_write("Invalid parameter for config_topopgraphic_wave_drag_scheme. It must be one of: 'JSL', 'ZAE', or 'LGF'" , MPAS_LOG_CRIT)
              call mpas_log_write("")
          end if
      endif

      if (config_disable_vel_topographic_wave_drag) topographicWaveDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_topographic_wave_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
