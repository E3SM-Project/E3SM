! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
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
   use mpas_constants
   use mpas_timer

   use ocn_constants

   use ocn_tracer_advection
   use ocn_tracer_short_wave_absorption

   use ocn_thick_hadv
   use ocn_thick_vadv
   use ocn_thick_surface_flux

   use ocn_vel_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vmix

   use ocn_tracer_hmix
   use ocn_high_freq_thickness_hmix_del2
   use ocn_tracer_surface_flux

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
             ocn_tend_freq_filtered_thickness, &
             ocn_tendency_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: apply_Dhf_to_hhf, use_highFreqThick_restore

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

   subroutine ocn_tend_thick(tendPool, forcingPool, diagnosticsPool, meshPool)!{{{
      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real (kind=RKIND), dimension(:), pointer :: surfaceMassFlux
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessEdge, vertTransportVelocityTop, tend_layerThickness, uTransport, transmissionCoefficients

      integer :: err

      logical, pointer :: config_disable_thick_all_tend

      call mpas_timer_start("ocn_tend_thick")

      call mpas_pool_get_config(ocnConfigs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)

      call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)

      call mpas_pool_get_array(tendPool, 'layerThickness', tend_layerThickness)

      call mpas_pool_get_array(forcingPool, 'surfaceMassFlux', surfaceMassFlux)
      call mpas_pool_get_array(forcingPool, 'transmissionCoefficients', transmissionCoefficients)
                  
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
      call ocn_thick_hadv_tend(meshPool, uTransport, layerThicknessEdge, tend_layerThickness, err)
      call mpas_timer_stop("hadv", thickHadvTimer)

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      call mpas_timer_start("vadv", .false., thickVadvTimer)
      call ocn_thick_vadv_tend(meshPool, vertTransportVelocityTop, tend_layerThickness, err)
      call mpas_timer_stop("vadv", thickVadvTimer)

      !
      ! surface flux tendency
      !
      call mpas_timer_start("surface flux", .false.)
      call ocn_thick_surface_flux_tend(meshPool, transmissionCoefficients, layerThickness, surfaceMassFlux, tend_layerThickness, err)
      call mpas_timer_stop("surface flux")

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

   subroutine ocn_tend_vel(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, scratchPool, timeLevelIn)!{{{
      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input: Scratch structure
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state fields

      real (kind=RKIND), dimension(:), pointer :: surfaceWindStress

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThicknessEdge, normalVelocity, tangentialVelocity, density, zMid, pressure, &
        tend_normalVelocity, circulation, relativeVorticity, viscosity, kineticEnergyCell, &
        normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, &
        montgomeryPotential, vertTransportVelocityTop, divergence, vertViscTopOfEdge

      integer :: timeLevel

      integer :: err

      logical, pointer :: config_disable_vel_all_tend
      character (len=StrKIND), pointer :: config_pressure_gradient_type

      call mpas_timer_start("ocn_tend_vel")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(ocnConfigs, 'config_pressure_gradient_type', config_pressure_gradient_type)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'relativeVorticity', relativeVorticity)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedRelativeVorticityEdge', normalizedRelativeVorticityEdge)
      call mpas_pool_get_array(diagnosticsPool, 'normalizedPlanetaryVorticityEdge', normalizedPlanetaryVorticityEdge)
      call mpas_pool_get_array(diagnosticsPool, 'divergence', divergence)
      call mpas_pool_get_array(diagnosticsPool, 'viscosity', viscosity)
      call mpas_pool_get_array(diagnosticsPool, 'montgomeryPotential', montgomeryPotential)
      call mpas_pool_get_array(diagnosticsPool, 'pressure', pressure)
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)

      call mpas_pool_get_array(tendPool, 'normalVelocity', tend_normalVelocity)
                  
      call mpas_pool_get_array(forcingPool, 'surfaceWindStress', surfaceWindStress)

      !
      ! velocity tendency: start accumulating tendency terms
      !
      tend_normalVelocity(:,:) = 0.0

      if(config_disable_vel_all_tend) return

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !

      call mpas_timer_start("coriolis", .false., velCorTimer)
      call ocn_vel_coriolis_tend(meshPool, normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, layerThicknessEdge, &
         normalVelocity, kineticEnergyCell, tend_normalVelocity, err)
      call mpas_timer_stop("coriolis", velCorTimer)

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      call mpas_timer_start("vadv", .false., velVadvTimer)
      call ocn_vel_vadv_tend(meshPool, normalVelocity, layerThicknessEdge, vertTransportVelocityTop, tend_normalVelocity, err)
      call mpas_timer_stop("vadv", velVadvTimer)

      !
      ! velocity tendency: pressure gradient
      !
      call mpas_timer_start("pressure grad", .false., velPgradTimer)
      if (config_pressure_gradient_type.eq.'MontgomeryPotential') then
          call ocn_vel_pressure_grad_tend(meshPool, montgomeryPotential,  zMid, density, tend_normalVelocity, err)
      else
          call ocn_vel_pressure_grad_tend(meshPool, pressure, zMid, density, tend_normalVelocity, err)
      end if
      call mpas_timer_stop("pressure grad", velPgradTimer)

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla relativeVorticity )
      !   strictly only valid for config_mom_del2 == constant
      !
      call mpas_timer_start("hmix", .false., velHmixTimer)
      call ocn_vel_hmix_tend(meshPool, divergence, relativeVorticity, normalVelocity, tangentialVelocity, viscosity, &
         tend_normalVelocity, scratchPool, err)
      call mpas_timer_stop("hmix", velHmixTimer)

      !
      ! velocity tendency: forcing and bottom drag
      !

      call mpas_timer_start("forcings", .false., velForceTimer)
      call ocn_vel_forcing_tend(meshPool, normalVelocity, surfaceWindStress, layerThicknessEdge, tend_normalVelocity, err)
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
   subroutine ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, dt, timeLevelIn)!{{{
      implicit none

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostic information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      real (kind=RKIND), intent(in) :: dt !< Input: Time step
      integer, intent(in), optional :: timeLevelIn

      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux
      real (kind=RKIND), dimension(:,:), pointer :: &
        uTransport, layerThickness,vertTransportVelocityTop, layerThicknessEdge, vertDiffTopOfCell, tend_layerThickness, normalThicknessFlux, surfaceTracerFlux, transmissionCoefficients
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracers, tend_tr

      integer :: err, iEdge, k
      integer, pointer :: nVertLevels, nEdges, indexTemperature
      integer :: timeLevel

      logical, pointer :: config_disable_tr_all_tend

      call mpas_timer_start("ocn_tend_tracer")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'uTransport', uTransport)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertTransportVelocityTop', vertTransportVelocityTop)

      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceTracerFlux', surfaceTracerFlux)
      call mpas_pool_get_array(forcingPool, 'transmissionCoefficients', transmissionCoefficients)

      call mpas_pool_get_array(tendPool, 'tracers', tend_tr)
      call mpas_pool_get_array(tendPool, 'layerThickness', tend_layerThickness)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)

      !
      ! initialize tracer tendency (RHS of tracer equation) to zero.
      !
      tend_tr(:,:,:) = 0.0

      if(config_disable_tr_all_tend) return

      allocate(normalThicknessFlux(nVertLevels, nEdges+1))
      !
      ! QC Comment (3/15/12): need to make sure that uTransport is the right
      ! transport velocity for the tracer.
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            normalThicknessFlux(k, iEdge) = uTransport(k, iEdge) * layerThicknessEdge(k, iEdge)
         end do
      end do

      !
      ! tracer tendency: horizontal advection term -div( layerThickness \phi u)
      !

      ! Monotonoic Advection, or standard advection
      call mpas_timer_start("adv", .false., tracerHadvTimer)
      call ocn_tracer_advection_tend(tracers, normalThicknessFlux, vertTransportVelocityTop, layerThickness, layerThickness, dt, meshPool, tend_layerThickness, tend_tr)
      call mpas_timer_stop("adv", tracerHadvTimer)

      !
      ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
      !
      call mpas_timer_start("hmix", .false., tracerHmixTimer)
      call ocn_tracer_hmix_tend(meshPool, layerThicknessEdge, tracers, tend_tr, err)
      call mpas_timer_stop("hmix", tracerHmixTimer)

      !
      ! Perform forcing from surface fluxes
      !
      call mpas_timer_start("surface_flux", .false.)
      call ocn_tracer_surface_flux_tend(meshPool, transmissionCoefficients, layerThickness, surfaceTracerFlux, tend_tr, err)
      call mpas_timer_stop("surface_flux")

      !
      ! Performing shortwave absorption
      !
      call mpas_timer_start("short wave", .false.)
      call ocn_tracer_short_wave_absorption_tend(meshPool, indexTemperature, layerThickness, penetrativeTemperatureFlux, tend_tr, err)
      call mpas_timer_stop("short wave")

      call mpas_timer_stop("ocn_tend_tracer")

      deallocate(normalThicknessFlux)

   end subroutine ocn_tend_tracer!}}}

!***********************************************************************
!
!  routine ocn_tend_freq_filtered_thickness
!
!> \brief   Compute tendencies needed for frequency filtered thickness
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine compute high frequency thickness tendency and the 
!>  low freqency divergence.  It is only called when
!>  config_freq_filtered_thickness is true (z-tilde)
!
!-----------------------------------------------------------------------
   subroutine ocn_tend_freq_filtered_thickness(tendPool, statePool, diagnosticsPool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state fields

      integer :: timeLevel
      integer :: err, iCell, i, k, iEdge
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeBot, nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: flux, invAreaCell, div_hu_btr, thickness_filter_timescale_sec, highFreqThick_restore_time_sec, &
         totalThickness
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, layerThicknessEdge, &
         layerThickness, &
         lowFreqDivergence, highFreqThickness, &
         tend_lowFreqDivergence, tend_highFreqThickness
      real (kind=RKIND), dimension(:), allocatable:: div_hu

      real (kind=RKIND), pointer :: config_thickness_filter_timescale, config_highFreqThick_restore_time

      call mpas_timer_start("ocn_tend_freq_filtered_thickness")
      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_thickness_filter_timescale', config_thickness_filter_timescale)
      call mpas_pool_get_config(ocnConfigs, 'config_highFreqThick_restore_time', config_highFreqThick_restore_time)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence, timeLevel)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)

      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', tend_lowFreqDivergence)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', tend_highFreqThickness)

      allocate(div_hu(nVertLevels))

      !
      ! Low Frequency Divergence and high frequency thickness Tendency
      !
      tend_lowFreqDivergence = 0.0
      tend_highFreqThickness = 0.0

      ! Convert restore time from days to seconds
      thickness_filter_timescale_sec = config_thickness_filter_timescale*86400.0
      highFreqThick_restore_time_sec = config_highFreqThick_restore_time*86400.0
      do iCell = 1, nCells
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
              -2.0 * pii / thickness_filter_timescale_sec &
              *(lowFreqDivergence(k,iCell)  - div_hu(k) &
                + div_hu_btr * layerThickness(k,iCell) / totalThickness)

           tend_highFreqThickness(k,iCell) = &
              - div_hu(k) + div_hu_btr * layerThickness(k,iCell) / totalThickness + lowFreqDivergence(k,iCell) &
              + use_highFreqThick_restore*( -2.0 * pii / highFreqThick_restore_time_sec * highFreqThickness(k,iCell) )

        end do

      end do

      deallocate(div_hu)

      !
      !  high frequency thickness tendency: del2 horizontal hhf diffusion, div(\kappa_{hf} \nabla h^{hf})
      !
      call mpas_timer_start("hmix", .false., tracerHmixTimer)
      call ocn_high_freq_thickness_hmix_del2_tend(meshPool, highFreqThickness, tend_highFreqThickness, err)
      call mpas_timer_stop("hmix", tracerHmixTimer)

      call mpas_timer_stop("ocn_tend_freq_filtered_thickness")

   end subroutine ocn_tend_freq_filtered_thickness!}}}

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

        logical, pointer :: config_use_highFreqThick_restore

        err = 0

        call mpas_pool_get_config(ocnConfigs, 'config_use_highFreqThick_restore', config_use_highFreqThick_restore)

        if (config_use_highFreqThick_restore) then
           use_highFreqThick_restore = 1
        else
           use_highFreqThick_restore = 0
        endif

    end subroutine ocn_tendency_init!}}}

!***********************************************************************

end module ocn_tendency

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
