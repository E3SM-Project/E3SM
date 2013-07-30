! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tendency
!
!> \brief MPAS ocean tendency driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routines for computing
!>  tendency terms for the ocean primitive equations.
!
!-----------------------------------------------------------------------

module ocn_tendency

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_timer

   use ocn_tracer_advection

   use ocn_thick_hadv
   use ocn_thick_vadv

   use ocn_vel_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vmix

   use ocn_tracer_hmix
   use ocn_high_freq_thickness_hmix_del2
   use ocn_restoring

   implicit none
   private
   save

   type (timer_node), pointer :: thickHadvTimer, thickVadvTimer
   type (timer_node), pointer :: velCorTimer, velVadvTimer, velPgradTimer, velHmixTimer, velForceTimer
   type (timer_node), pointer :: tracerHadvTimer, tracerVadvTimer, tracerHmixTimer, tracerRestoringTimer

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

   public :: ocn_tend_thick, &
             ocn_tend_vel, &
             ocn_tend_tracer, &
             ocn_tend_frequency_filtered_thickness, &
             ocn_tendency_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: apply_Dhf_to_hhf, apply_hff_restore_to_hhf

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tend_thick
!
!> \brief   Computes thickness tendency
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine computes the thickness tendency for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_thick(tend, s, grid)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessEdge, vertTransportVelocityTop, tend_layerThickness, uTransport

      integer :: err

      call mpas_timer_start("ocn_tend_thick")

      uTransport  => s % uTransport % array
      vertTransportVelocityTop        => s % vertTransportVelocityTop % array
      layerThicknessEdge      => s % layerThicknessEdge % array

      tend_layerThickness      => tend % layerThickness % array
                  
      !
      ! height tendency: start accumulating tendency terms
      !
      tend_layerThickness = 0.0

      if(config_disable_thick_all_tend) return

      !
      ! height tendency: horizontal advection term -\nabla\cdot ( hu)
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3. 
      ! for explanation of divergence operator.
      !
      ! QC Comment (3/15/12): need to make sure that uTranport is the right
      ! transport velocity here.
      call mpas_timer_start("hadv", .false., thickHadvTimer)
      call ocn_thick_hadv_tend(grid, uTransport, layerThicknessEdge, tend_layerThickness, err)
      call mpas_timer_stop("hadv", thickHadvTimer)

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      call mpas_timer_start("vadv", .false., thickVadvTimer)
      call ocn_thick_vadv_tend(grid, vertTransportVelocityTop, tend_layerThickness, err)
      call mpas_timer_stop("vadv", thickVadvTimer)

      call mpas_timer_stop("ocn_tend_thick")
   
   end subroutine ocn_tend_thick!}}}

!***********************************************************************
!
!  routine ocn_tend_vel
!
!> \brief   Computes velocity tendency
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine computes the velocity tendency for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_vel(tend, s, d, grid)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (diagnostics_type), intent(in) :: d !< Input: Diagnostic information
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThicknessEdge, h, normalVelocity, density, zMid, pressure, &
        tend_normalVelocity, circulation, relativeVorticity, viscosity, kineticEnergyCell, &
        normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, &
        montgomeryPotential, vertTransportVelocityTop, divergence, vertViscTopOfEdge

      real (kind=RKIND), dimension(:,:), pointer :: normalVelocityForcing

      integer :: err

      call mpas_timer_start("ocn_tend_vel")

      normalVelocity => s % normalVelocity % array
      density         => s % density % array
      vertTransportVelocityTop        => s % vertTransportVelocityTop % array
      zMid        => s % zMid % array
      layerThicknessEdge      => s % layerThicknessEdge % array
      viscosity   => s % viscosity % array
      relativeVorticity   => s % relativeVorticity % array
      divergence  => s % divergence % array
      kineticEnergyCell          => s % kineticEnergyCell % array
      normalizedRelativeVorticityEdge     => s % normalizedRelativeVorticityEdge % array
      normalizedPlanetaryVorticityEdge     => s % normalizedPlanetaryVorticityEdge % array
      montgomeryPotential     => s % montgomeryPotential % array
      pressure    => s % pressure % array
      vertViscTopOfEdge => d % vertViscTopOfEdge % array

      tend_normalVelocity      => tend % normalVelocity % array
                  
      normalVelocityForcing => grid % normalVelocityForcing % array

      !
      ! velocity tendency: start accumulating tendency terms
      !
      tend_normalVelocity(:,:) = 0.0

      if(config_disable_vel_all_tend) return

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !

      call mpas_timer_start("coriolis", .false., velCorTimer)
      call ocn_vel_coriolis_tend(grid, normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, layerThicknessEdge, normalVelocity, kineticEnergyCell, tend_normalVelocity, err)
      call mpas_timer_stop("coriolis", velCorTimer)

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      call mpas_timer_start("vadv", .false., velVadvTimer)
      call ocn_vel_vadv_tend(grid, normalVelocity, layerThicknessEdge, vertTransportVelocityTop, tend_normalVelocity, err)
      call mpas_timer_stop("vadv", velVadvTimer)

      !
      ! velocity tendency: pressure gradient
      !
      call mpas_timer_start("pressure grad", .false., velPgradTimer)
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then
          call ocn_vel_pressure_grad_tend(grid, montgomeryPotential,  zMid, density, tend_normalVelocity, err)
      else
          call ocn_vel_pressure_grad_tend(grid, pressure, zMid, density, tend_normalVelocity, err)
      end if
      call mpas_timer_stop("pressure grad", velPgradTimer)

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla relativeVorticity )
      !   strictly only valid for config_mom_del2 == constant
      !
      call mpas_timer_start("hmix", .false., velHmixTimer)
      call ocn_vel_hmix_tend(grid, divergence, relativeVorticity, viscosity, tend_normalVelocity, err)
      call mpas_timer_stop("hmix", velHmixTimer)

      !
      ! velocity tendency: forcing and bottom drag
      !

      call mpas_timer_start("forcings", .false., velForceTimer)
      call ocn_vel_forcing_tend(grid, normalVelocity, normalVelocityForcing, layerThicknessEdge, tend_normalVelocity, err)
      call mpas_timer_stop("forcings", velForceTimer)

      !
      ! velocity tendency: vertical mixing d/dz( nu_v du/dz))
      !
      call mpas_timer_stop("ocn_tend_vel")

   end subroutine ocn_tend_vel!}}}

!***********************************************************************
!
!  routine ocn_tend_tracer
!
!> \brief   Computes tracer tendency
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details 
!>  This routine computes tracer tendencies for the ocean
!
!-----------------------------------------------------------------------
   subroutine ocn_tend_tracer(tend, s, d, grid, dt)!{{{
      implicit none

      type (tend_type), intent(inout) :: tend !< Input/Output: Tendency structure
      type (state_type), intent(in) :: s !< Input: State information
      type (diagnostics_type), intent(in) :: d !< Input: Diagnostic information
      type (mesh_type), intent(in) :: grid !< Input: Grid information
      real (kind=RKIND), intent(in) :: dt !< Input: Time step

      real (kind=RKIND), dimension(:,:), pointer :: &
        uTransport, layerThickness,vertTransportVelocityTop, layerThicknessEdge, vertDiffTopOfCell, tend_layerThickness, normalThicknessFlux
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr

      integer :: err, iEdge, k

      call mpas_timer_start("ocn_tend_tracer")

      uTransport  => s % uTransport % array
      layerThickness           => s % layerThickness % array
      vertTransportVelocityTop        => s % vertTransportVelocityTop % array
      tracers     => s % tracers % array
      layerThicknessEdge      => s % layerThicknessEdge % array
      vertDiffTopOfCell => d % vertDiffTopOfCell % array

      tend_tr     => tend % tracers % array
      tend_layerThickness      => tend % layerThickness % array

      !
      ! initialize tracer tendency (RHS of tracer equation) to zero.
      !
      tend_tr(:,:,:) = 0.0

      if(config_disable_tr_all_tend) return

      allocate(normalThicknessFlux(grid % nVertLevels, grid % nEdges+1))
      !
      ! QC Comment (3/15/12): need to make sure that uTransport is the right
      ! transport velocity for the tracer.
      do iEdge = 1, grid % nEdges
         do k = 1, grid % nVertLevels
            normalThicknessFlux(k, iEdge) = uTransport(k, iEdge) * layerThicknessEdge(k, iEdge)
         end do
      end do

      !
      ! tracer tendency: horizontal advection term -div( layerThickness \phi u)
      !

      ! Monotonoic Advection, or standard advection
      call mpas_timer_start("adv", .false., tracerHadvTimer)
      call ocn_tracer_advection_tend(tracers, normalThicknessFlux, vertTransportVelocityTop, layerThickness, layerThickness, dt, grid, tend_layerThickness, tend_tr)
      call mpas_timer_stop("adv", tracerHadvTimer)

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      call mpas_timer_start("hmix", .false., tracerHmixTimer)
      call ocn_tracer_hmix_tend(grid, layerThicknessEdge, tracers, tend_tr, err)
      call mpas_timer_stop("hmix", tracerHmixTimer)

      !
      ! add restoring to T and S in top model layer
      !
      call mpas_timer_start("restoring", .false., tracerRestoringTimer)

      call ocn_restoring_tend(grid, layerThickness, s%index_temperature, s%index_salinity, tracers, tend_tr, err)

      call mpas_timer_stop("restoring", tracerRestoringTimer)

      call mpas_timer_stop("ocn_tend_tracer")

      deallocate(normalThicknessFlux)

   end subroutine ocn_tend_tracer!}}}

!***********************************************************************
!
!  routine ocn_tend_frequency_filtered_thickness
!
!> \brief   Compute tendencies needed for frequency filtered thickness
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine compute high frequency thickness tendency and the 
!>  low freqency divergence.  It is only called when
!>  frequency_filtered_thickness is on (z-tilde)
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_frequency_filtered_thickness(tend, s, grid)!{{{

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid

      integer :: err, nCells, nVertLevels, iCell, i, k, iEdge
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: flux, invAreaCell, div_hu_btr, Dlf_restore_time_sec, hhf_restore_time_sec, &
         totalThickness
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThicknessEdge, &
         layerThickness, &
         lowFreqDivergence, highFreqThickness, &
         tend_lowFreqDivergence, tend_highFreqThickness
      real (kind=RKIND), dimension(:), allocatable:: div_hu

      call mpas_timer_start("ocn_tend_frequency_filtered_thickness")
      err = 0

      nCells      = grid % nCells
      nVertLevels = grid % nVertLevels

      nEdgesOnCell      => grid % nEdgesOnCell % array
      areaCell          => grid % areaCell % array
      edgesOnCell       => grid % edgesOnCell % array
      edgeSignOnCell    => grid % edgeSignOnCell % array
      maxLevelCell      => grid % maxLevelCell % array
      maxLevelEdgeBot   => grid % maxLevelEdgeBot % array
      dvEdge            => grid % dvEdge % array

      normalVelocity => s % normalVelocity % array
      layerThickness => s % layerThickness % array
      layerThicknessEdge => s % layerThicknessEdge % array
      lowFreqDivergence => s % lowFreqDivergence % array
      highFreqThickness => s % highFreqThickness % array
      tend_lowFreqDivergence => tend % lowFreqDivergence % array
      tend_highFreqThickness => tend % highFreqThickness % array

      allocate(div_hu(nVertLevels))

      !
      ! Low Frequency Divergence and high frequency thickness Tendency
      !
      tend_lowFreqDivergence = 0.0
      tend_highFreqThickness = 0.0

      ! Convert restore time from days to seconds
      Dlf_restore_time_sec = config_Dlf_restore_time*86400.0
      hhf_restore_time_sec = config_hhf_restore_time*86400.0
      do iCell=1,nCells
        div_hu(:) = 0.0
        div_hu_btr = 0.0
        invAreaCell = 1.0 / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          do k = 1, maxLevelEdgeBot(iEdge)
            flux = layerThicknessEdge(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
            div_hu(k) = div_hu(k) - flux
            div_hu_btr = div_hu_btr - flux
          end do
        end do

        totalThickness = sum(layerThickness(1:maxLevelCell(iCell),iCell))
        do k = 1, maxLevelCell(iCell)

           tend_lowFreqDivergence(k,iCell) = &
              -2.0*pii/Dlf_restore_time_sec &
              *(lowFreqDivergence(k,iCell)  - div_hu(k) &
                + div_hu_btr*layerThickness(k,iCell)/totalThickness)

           tend_highFreqThickness(k,iCell) = &
                apply_Dhf_to_hhf*( - div_hu(k) + div_hu_btr*layerThickness(k,iCell)/totalThickness &
                                   + lowFreqDivergence(k,iCell) ) &
              + apply_hff_restore_to_hhf*( -2.0*pii/hhf_restore_time_sec * highFreqThickness(k,iCell) )

        end do

      end do

      deallocate(div_hu)

      !
      !  high frequency thickness tendency: del2 horizontal hhf diffusion, div(\kappa_{hf} \nabla h^{hf})
      !
      call mpas_timer_start("hmix", .false., tracerHmixTimer)
      call ocn_high_freq_thickness_hmix_del2_tend(grid, highFreqThickness, tend_highFreqThickness, err)
      call mpas_timer_stop("hmix", tracerHmixTimer)

      call mpas_timer_stop("ocn_tend_frequency_filtered_thickness")

   end subroutine ocn_tend_frequency_filtered_thickness!}}}

!***********************************************************************
!
!  routine ocn_tendency_init
!
!> \brief   Initializes flags used within tendency routines.
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    4 November 2011
!> \details 
!>  This routine initializes flags related to quantities computed within
!>  other tendency routines.
!
!-----------------------------------------------------------------------
    subroutine ocn_tendency_init(err)!{{{
        integer, intent(out) :: err !< Output: Error flag

        err = 0

	if (config_apply_Dhf_to_hhf) apply_Dhf_to_hhf = 1
	if (config_apply_hff_restore_to_hhf) apply_hff_restore_to_hhf = 1

! mrp always turn to 1;
   apply_Dhf_to_hhf = 1
   apply_hff_restore_to_hhf = 1

!mrp end
    end subroutine ocn_tendency_init!}}}

!***********************************************************************

end module ocn_tendency

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
