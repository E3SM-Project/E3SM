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
!> \author Nairita Pal, Kristin Barton
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
!> \author  Nairita Pal, Kristin Barton
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
            if (k>0) then
              temp_twd(iEdge) = topographicWaveDrag(iEdge) * normalVelocity(k,iEdge)
              tend(k,iEdge) = tend(k,iEdge) - temp_twd(iEdge)
            endif
          enddo
          !$omp end do
      else
          !$omp do schedule(runtime) private(k)
          do iEdge = 1, nEdgesOwned
            k =  maxLevelEdgeTop(iEdge)
            if (k>0) then
              temp_twd(iEdge) = topographicWaveDrag(iEdge) * (normalVelocity(k,iEdge) * normalCoeffTWD(iEdge) + tangentialVelocity(k,iEdge) * tangentialCoeffTWD(iEdge))
              tend(k,iEdge) = tend(k,iEdge) - temp_twd(iEdge)
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
!> \author  Nairita Pal, Kristin Barton
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
      real (kind=RKIND) :: kappa
      real (kind=RKIND) :: N0, L, Nb, Nbar, gam
      real (kind=RKIND) :: topoDragTurnOff

      integer :: iEdge, k, kmax

      err = 0


      topographicWaveDragCoeff = 0.0_RKIND

      if (config_use_topographic_wave_drag) then
          topographicWaveDragOn = .true.
          topographicWaveDragCoeff = config_topographic_wave_drag_coeff
          topoDragCutoffDepth = config_topographic_wave_drag_cutoff_depth
          topoDragCutoffWidth = config_topographic_wave_drag_cutoff_width
          if (config_topographic_wave_drag_scheme.EQ."JSL") then
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Jayne and St. Laurent")
              call mpas_log_write("")
              tensorScheme = .false.
              kappa = pi/10000.0_RKIND

              do iEdge = 1, nEdgesAll
                  kmax =  maxLevelEdgeTop(iEdge)
                  topographicWaveDrag(iEdge) = 0.0_RKIND
                  if (kmax > 0) then
                      k =  maxLevelEdgeTop(iEdge)
                      topographicWaveDrag(iEdge) = topographicWaveDragCoeff*kappa &
                                                 * bathy_stddev(iEdge)**2 &
                                                 * topo_buoyancy_N1B(iEdge) &
                                                 / layerThickEdgeMean(k,iEdge)
                  endif
              end do

          else if (config_topographic_wave_drag_scheme.EQ."ZAE") then
              tensorScheme = .false.
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Zaron and Egbert")
              call mpas_log_write("")

              ! Cze = Gamma H grad(H)^2 Nb Nbar / (8 pi^2 omega)
              ! Gamma = 50
              ! N(z) = N0 exp(z/1300)
              ! N0 = 5.24 x 10^-3
              ! z = 0 at surface
              ! Nb is N(z) evaluated at z = -H
              ! Nbar = 1300 N0 (1 - exp(-H/1300))/H

              N0 = 5.24e-3_RKIND
              L = 1300.0_RKIND
              gam = 50.0_RKIND
              do iEdge = 1, nEdgesAll
                  kmax =  maxLevelEdgeTop(iEdge)
                  topographicWaveDrag(iEdge) = 0.0_RKIND
                  if (kmax > 0) then
                      k = 1 ! Note: Only works for single layer case

                      Nb = N0*exp(-layerThickEdgeMean(k,iEdge)/L)
                      Nbar = N0*L*(1.0_RKIND - exp(-layerThickEdgeMean(k,iEdge)/L)) / layerThickEdgeMean(k,iEdge)

                      topographicWaveDrag(iEdge) = topographicWaveDragCoeff * gam  &
                                                 * bed_slope_edges(iEdge)**2 * Nbar * Nb &
                                                 / (8.0_RKIND * omegaM2 * pi**2)
                  endif
              enddo
          else if (config_topographic_wave_drag_scheme.EQ."LGF") then
              tensorScheme = .true.
              call mpas_log_write("")
              call mpas_log_write(" Topographic wave drag scheme is: Local Generation Formula")
              call mpas_log_write("")
              do iEdge = 1, nEdgesAll
                  kmax =  maxLevelEdgeTop(iEdge)
                  topographicWaveDrag(iEdge) = 0.0_RKIND
                  if (kmax > 0) then
                      topographicWaveDrag(iEdge) = topographicWaveDragCoeff  &
                                                 * (sqrt((topo_buoyancy_N1B(iEdge)**2 - omegaM2**2)   &
                                                       * (topo_buoyancy_N1V(iEdge)**2 - omegaM2**2))) &
                                                 / (4.0_RKIND*pi*omegaM2)
                      normalCoeffTWD(iEdge) = (lonGradEdge(iEdge)**2) &
                                         * (cos(angleEdge(iEdge))**2) + (latGradEdge(iEdge)**2)*(sin(angleEdge(iEdge))**2) &
                                         - 2.0_RKIND*latGradEdge(iEdge)*lonGradEdge(iEdge)*sin(angleEdge(iEdge))*(cos(angleEdge(iEdge)))
                      tangentialCoeffTWD(iEdge) = sin(angleEdge(iEdge)) &
                                             * cos(angleEdge(iEdge)) &
                                             * (lonGradEdge(iEdge)**2 - latGradEdge(iEdge)**2) &
                                             + latGradEdge(iEdge)*lonGradEdge(iEdge)*(cos(angleEdge(iEdge))**2 - sin(angleEdge(iEdge))**2)
                  endif
              enddo
          else 
              call mpas_log_write("")

              call mpas_log_write("Invalid parameter for config_topopgraphic_wave_drag_scheme. It must be one of: 'JSL', 'ZAE', or 'LGF'" , MPAS_LOG_CRIT)
              call mpas_log_write("")
          end if

          ! Apply ramping factor to wave drag
          do iEdge = 1, nEdgesAll
              kmax =  maxLevelEdgeTop(iEdge)
              if (kmax > 0) then
                  k = 1
                  topoDragTurnOff = (TANH((layerThickEdgeMean(k,iEdge)-topoDragCutoffDepth)/topoDragCutoffWidth)+1.0_RKIND)/2.0_RKIND
                  topographicWaveDrag(iEdge) = topoDragTurnOff * topographicWaveDrag(iEdge)
              endif
          enddo
          
      endif

      if (config_disable_vel_topographic_wave_drag) topographicWaveDragOn = .false.

   !--------------------------------------------------------------------

   end subroutine ocn_vel_forcing_topographic_wave_drag_init!}}}

!***********************************************************************

end module ocn_vel_forcing_topographic_wave_drag

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
