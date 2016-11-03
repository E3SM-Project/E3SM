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

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer
   use mpas_threading
   use ocn_diagnostics
   use ocn_constants

   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes
   use ocn_frazil_forcing

   use ocn_tracer_hmix
   use ocn_high_freq_thickness_hmix_del2
   use ocn_tracer_advection
   use ocn_tracer_short_wave_absorption
   use ocn_tracer_nonlocalflux
   use ocn_tracer_surface_restoring
   use ocn_tracer_interior_restoring
   use ocn_tracer_exponential_decay
   use ocn_tracer_ideal_age
   use ocn_tracer_TTD
   use ocn_tracer_surface_flux_to_tend
   use ocn_tracer_ecosys
   use ocn_tracer_DMS
   use ocn_tracer_MacroMolecules

   use ocn_thick_hadv
   use ocn_thick_vadv
   use ocn_thick_surface_flux

   use ocn_vel_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vmix

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
      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information
      type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information

      real (kind=RKIND), dimension(:), pointer :: surfaceThicknessFlux
      real (kind=RKIND), dimension(:), pointer :: surfaceThicknessFluxRunoff
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, layerThicknessEdge, &
         vertAleTransportTop, tend_layerThickness, normalTransportVelocity, fractionAbsorbed, fractionAbsorbedRunoff

      integer, pointer :: nCells
      integer :: err, iCell

      logical, pointer :: config_disable_thick_all_tend

      call mpas_pool_get_config(ocnConfigs, 'config_disable_thick_all_tend', config_disable_thick_all_tend)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)

      call mpas_pool_get_array(tendPool, 'layerThickness', tend_layerThickness)

      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxRunoff', surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', fractionAbsorbedRunoff)

      !
      ! height tendency: start accumulating tendency terms
      !
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         tend_layerThickness(:, iCell) = 0.0_RKIND
         surfaceThicknessFlux(iCell) = 0.0_RKIND
         surfaceThicknessFluxRunoff(iCell) = 0.0_RKIND
      end do
      !$omp end do

      if(config_disable_thick_all_tend) return

      call mpas_timer_start("ocn_tend_thick")

      ! Build suface mass flux array from bulk
      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, surfaceThicknessFlux, surfaceThicknessFluxRunoff, err)

      ! Build suface thickness flux array from land ice
      call ocn_surface_land_ice_fluxes_thick(meshPool, forcingPool, surfaceThicknessFlux, err)

      !
      ! height tendency: horizontal advection term -\nabla\cdot ( hu)
      !
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      ! for explanation of divergence operator.
      !
      ! QC Comment (3/15/12): need to make sure that uTranport is the right
      ! transport velocity here.
      call ocn_thick_hadv_tend(meshPool, normalTransportVelocity, layerThicknessEdge, tend_layerThickness, err)

      !
      ! height tendency: vertical advection term -d/dz(hw)
      !
      call ocn_thick_vadv_tend(meshPool, vertAleTransportTop, tend_layerThickness, err)

      !
      ! surface flux tendency
      !

      call ocn_thick_surface_flux_tend(meshPool, fractionAbsorbed, fractionAbsorbedRunoff, layerThickness,  &
         surfaceThicknessFlux, surfaceThicknessFluxRunoff, tend_layerThickness, err)

      !
      ! surface flux tendency
      !
      call ocn_frazil_forcing_layer_thickness(meshPool, forcingPool, tend_layerThickness, err)

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

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer :: surfaceStress, surfaceStressMagnitude, surfaceFluxAttenuationCoefficient
      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThicknessEdge, normalVelocity, tangentialVelocity, density, potentialDensity, zMid, pressure, &
        tend_normalVelocity, circulation, relativeVorticity, viscosity, kineticEnergyCell, &
        normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, &
        montgomeryPotential, vertAleTransportTop, divergence, vertViscTopOfEdge, &
        inSituThermalExpansionCoeff, inSituSalineContractionCoeff
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      integer :: timeLevel

      integer :: err, iEdge, iCell
      integer, pointer :: indexTemperature, indexSalinity, nEdges, nCells

      logical, pointer :: config_disable_vel_all_tend
      character (len=StrKIND), pointer :: config_pressure_gradient_type

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_config(ocnConfigs, 'config_disable_vel_all_tend', config_disable_vel_all_tend)
      call mpas_pool_get_config(ocnConfigs, 'config_pressure_gradient_type', config_pressure_gradient_type)

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(diagnosticsPool, 'kineticEnergyCell', kineticEnergyCell)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
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
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'tangentialVelocity', tangentialVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficient', surfaceFluxAttenuationCoefficient)

      call mpas_pool_get_array(tendPool, 'normalVelocity', tend_normalVelocity)

      call mpas_pool_get_array(forcingPool, 'surfaceStress', surfaceStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', surfaceStressMagnitude)

      !
      ! velocity tendency: start accumulating tendency terms
      !
      !$omp do schedule(runtime)
      do iEdge = 1, nEdges
         tend_normalVelocity(:, iEdge) = 0.0_RKIND
         surfaceStress(iEdge) = 0.0_RKIND
      end do
      !$omp end do

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         surfaceStressMagnitude(iCell) = 0.0_RKIND
      end do
      !$omp end do

      if(config_disable_vel_all_tend) return

      call mpas_timer_start("ocn_tend_vel")

      ! Build bulk forcing suface stress
      call ocn_surface_bulk_forcing_vel(meshPool, forcingPool, surfaceStress, surfaceStressMagnitude, err)

      ! Add top drag to suface stress
      call ocn_surface_land_ice_fluxes_vel(meshPool, diagnosticsPool, surfaceStress, surfaceStressMagnitude, err)

      !
      ! velocity tendency: nonlinear Coriolis term and grad of kinetic energy
      !
      call ocn_vel_coriolis_tend(meshPool, normalizedRelativeVorticityEdge, normalizedPlanetaryVorticityEdge, layerThicknessEdge, &
         normalVelocity, kineticEnergyCell, tend_normalVelocity, err)

      !
      ! velocity tendency: vertical advection term -w du/dz
      !
      call ocn_vel_vadv_tend(meshPool, normalVelocity, layerThicknessEdge, vertAleTransportTop, tend_normalVelocity, err)

      !
      ! velocity tendency: pressure gradient
      !
      if (config_pressure_gradient_type.eq.'Jacobian_from_TS') then
         ! only pass EOS derivatives if needed.
         call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff',inSituThermalExpansionCoeff)
         call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
         call ocn_vel_pressure_grad_tend(meshPool, pressure, montgomeryPotential, zMid, density, potentialDensity, &
              indexTemperature, indexSalinity, activeTracers, tend_normalVelocity, err, &
              inSituThermalExpansionCoeff,inSituSalineContractionCoeff)
      else
         call ocn_vel_pressure_grad_tend(meshPool, pressure, montgomeryPotential, zMid, density, potentialDensity, &
              indexTemperature, indexSalinity, activeTracers, tend_normalVelocity, err, &
              inSituThermalExpansionCoeff,inSituSalineContractionCoeff)
      endif

      !
      ! velocity tendency: del2 dissipation, \nu_2 \nabla^2 u
      !   computed as \nu( \nabla divergence + k \times \nabla relativeVorticity )
      !   strictly only valid for config_mom_del2 == constant
      !
      call ocn_vel_hmix_tend(meshPool, scratchPool, divergence, relativeVorticity, normalVelocity, tangentialVelocity, viscosity, &
         tend_normalVelocity, err)

      !
      ! velocity tendency: forcing and bottom drag
      !
      call ocn_vel_forcing_tend(meshPool, normalVelocity, surfaceFluxAttenuationCoefficient, surfaceStress, layerThicknessEdge, &
                                tend_normalVelocity, err)


      !
      ! velocity tendency: vertical mixing d/dz( nu_v du/dz))
      !
      call mpas_timer_stop("ocn_tend_vel")
      call mpas_threading_barrier()

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
   subroutine ocn_tend_tracer(tendPool, statePool, forcingPool, diagnosticsPool, meshPool, swForcingPool, scratchPool, dt, & !{{{
                              timeLevelIn)
      implicit none

      !
      ! intent in/out
      !
      type (mpas_pool_type), intent(inout) :: tendPool        !< Input/Output: Tendency structure
      type (mpas_pool_type), intent(in)    :: statePool       !< Input: State information
      type (mpas_pool_type), intent(inout) :: forcingPool     !< Input: Forcing information
      type (mpas_pool_type), intent(inout)    :: diagnosticsPool !< Input: Diagnostic information
      type (mpas_pool_type), intent(inout)    :: meshPool     !< Input: Mesh information
      type (mpas_pool_type), intent(in)    :: swForcingPool   !< Input: sw data input info
      type (mpas_pool_type), intent(in)    :: scratchPool     !< Input: Scratch information
      real (kind=RKIND), intent(in) :: dt                     !< Input: Time step
      integer, intent(in), optional :: timeLevelIn            !< Input/Optional: Time Level Indes

      !
      ! additional pools
      !
      type (mpas_pool_type), pointer :: tracersPool, tracersTendPool            ! tracers and their tendency
      type (mpas_pool_type), pointer :: tracersSurfaceFluxPool                  ! surface fluxes
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool       ! surface restoring
      type (mpas_pool_type), pointer :: tracersInteriorRestoringFieldsPool      ! interior restoring
      type (mpas_pool_type), pointer :: tracersExponentialDecayFieldsPool       ! exponential decay
      type (mpas_pool_type), pointer :: tracersIdealAgeFieldsPool               ! ideal age
      type (mpas_pool_type), pointer :: tracersTTDFieldsPool                    ! transit time distribution

      ! scalar pointers
      integer :: nTracerGroup
      integer, pointer :: nVertLevels, nEdges, nCells, nCellsSolve, indexTemperature, indexSalinity
      logical, pointer :: config_disable_tr_all_tend, config_use_cvmix_kpp
      logical, pointer :: config_use_tracerGroup, config_use_tracerGroup_surface_bulk_forcing, &
                          config_use_tracerGroup_surface_restoring, config_use_tracerGroup_interior_restoring, &
                          config_use_tracerGroup_exponential_decay, config_use_tracerGroup_idealAge_forcing, &
                          config_use_tracerGroup_ttd_forcing, config_use_surface_salinity_monthly_restoring

      ! iterator for tracer categories
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: modifiedGroupName
      character (len=StrKIND) :: modifiedConfigName
      !
      ! one dimensional pointers
      !
      real (kind=RKIND), dimension(:), pointer :: penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL
      real (kind=RKIND), dimension(:), pointer :: tracerGroupExponentialDecayRate
      integer, dimension(:), pointer :: maxLevelCell

      !
      ! two dimensional pointers
      !
      real (kind=RKIND), dimension(:,:), pointer :: tracerGroupPistonVelocity, tracerGroupSurfaceRestoringValue, &
                                                    tracerGroupIdealAgeMask, tracerGroupTTDMask

      real (kind=RKIND), dimension(:,:), pointer :: &
        normalTransportVelocity, layerThickness,vertAleTransportTop, layerThicknessEdge, vertDiffTopOfCell, &
        tend_layerThickness, normalThicknessFlux, tracerGroupSurfaceFlux, fractionAbsorbed, zMid, relativeSlopeTopOfEdge, &
        relativeSlopeTapering, relativeSlopeTaperingCell, fractionAbsorbedRunoff, tracerGroupSurfaceFluxRunoff, &
        tracerGroupSurfaceFluxRemoved, nonLocalSurfaceTracerFlux

      !
      ! three dimensional pointers
      !
      real (kind=RKIND), dimension(:,:,:), pointer :: &
        tracerGroup, tracerGroupTend, vertNonLocalFlux

      real (kind=RKIND), dimension(:,:,:), pointer :: &
        activeTracers, &  !  need T, S for ecosys
        ecosysTracers     !  need ecosys for DMS and MacroMolecules

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroupInteriorRestoringRate, tracerGroupInteriorRestoringValue

      !
      ! Field pointers
      !
      type (field2DReal), pointer :: normalThicknessFluxField

      !
      ! local integers/reals/logicals
      !
      integer :: err, iCell, iEdge, k, timeLevel, nTracersEcosys

      !
      ! set time level of optional argument is present
      !
      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! get tracers pools
      !
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)

      !
      ! get dimensions
      !
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      !
      ! get configure options
      !
      call mpas_pool_get_config(ocnConfigs, 'config_disable_tr_all_tend', config_disable_tr_all_tend)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)

      !
      ! get arrays
      !
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTopOfEdge', relativeSlopeTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTapering', relativeSlopeTapering)
      call mpas_pool_get_array(diagnosticsPool, 'relativeSlopeTaperingCell', relativeSlopeTaperingCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(diagnosticsPool, 'penetrativeTemperatureFluxOBL', penetrativeTemperatureFluxOBL)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', fractionAbsorbedRunoff)
      call mpas_pool_get_array(tendPool, 'layerThickness', tend_layerThickness)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      if(config_disable_tr_all_tend) return

      call mpas_timer_start("ocn_tend_tracer")

      !allocate(normalThicknessFlux(nVertLevels, nEdges+1))
      call mpas_pool_get_field(scratchPool, 'normalThicknessFlux', normalThicknessFluxField)
      call mpas_allocate_scratch_field(normalThicknessFluxField, .true.)
      call mpas_threading_barrier()

      normalThicknessFlux => normalThicknessFluxField % array

      !
      ! transport velocity for the tracer.
      !
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            normalThicknessFlux(k, iEdge) = normalTransportVelocity(k, iEdge) * layerThicknessEdge(k, iEdge)
         end do
      end do
      !$omp end do

      !
      ! begin iterate over tracer categories
      !
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then

            ! load configure setting for this category
            !
            modifiedConfigName = 'config_use_' // trim(groupItr % memberName)
            call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup)

            if ( config_use_tracerGroup ) then
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_surface_bulk_forcing'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_surface_bulk_forcing)
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_surface_restoring'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_surface_restoring)
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_interior_restoring'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_interior_restoring)
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_exponential_decay'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_exponential_decay)
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_idealAge_forcing'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_idealAge_forcing)
               modifiedConfigName = 'config_use_' // trim(groupItr % memberName) // '_ttd_forcing'
               call mpas_pool_get_config(ocnConfigs, modifiedConfigName, config_use_tracerGroup_ttd_forcing)


               ! Get tracer group, and other groups (tendencies, etc.)
               call mpas_pool_get_array(tracersPool, trim(groupItr % memberName), tracerGroup, timeLevel)
               nTracerGroup = size(tracerGroup, dim=1)

               ! Get Tendency array
               modifiedGroupName = trim(groupItr % memberName) // "Tend"
               call mpas_pool_get_array(tracersTendPool, trim(modifiedGroupName), tracerGroupTend)

               ! Get surface flux array
               modifiedGroupName = trim(groupItr % memberName) // "SurfaceFlux"
               call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFlux)

               ! Get Array of total surface temp/salt flux (includes thickness
               ! tendencies
               call mpas_pool_get_array(tracersSurfaceFluxPool, 'nonLocalSurfaceTracerFlux', nonLocalSurfaceTracerFlux)

               ! Get surface flux due to river runoff array
!maltrud only active tracers have runoff flux for now
               if (trim(groupItr % memberName) == 'activeTracers') then
               modifiedGroupName = trim(groupItr % memberName) // "SurfaceFluxRunoff"
               call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFluxRunoff)
               endif

               ! Get surface flux removed array to keep track of how much flux is ignored
               modifiedGroupName = trim(groupItr % memberName) // "SurfaceFluxRemoved"
               call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFluxRemoved)

               !
               ! initialize tracer surface flux and tendency to zero.
               !
               !$omp do schedule(runtime)
               do iCell = 1, nCells
                 tracerGroupTend(:,:, iCell) = 0.0_RKIND
                 tracerGroupSurfaceFlux(:, iCell) = 0.0_RKIND
               end do
               !$omp end do

               if (trim(groupItr % memberName) == 'activeTracers') then
                  !$omp do schedule(runtime)
                  do iCell = 1, nCells
                    tracerGroupSurfaceFluxRunoff(:, iCell) = 0.0_RKIND
                    tracerGroupSurfaceFluxRemoved(:, iCell) = 0.0_RKIND
                  end do
                  !$omp end do
               endif

               !
               ! fill components of surface tracer flux
               !
               if (config_use_tracerGroup_surface_bulk_forcing) then
                  call ocn_surface_bulk_forcing_tracers(meshPool, groupItr % memberName, forcingPool, tracerGroup, &
                                                        tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff, &
                                                        tracerGroupSurfaceFluxRemoved, dt, layerThickness, err)
               end if

               !
               ! compute ecosystem source-sink tendencies and net surface fluxes
               ! NOTE: must be called before ocn_tracer_surface_flux_tend
               !
               if ( trim(groupItr % memberName) == 'ecosysTracers' ) then
                  call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
                  call ocn_tracer_ecosys_compute(activeTracers, tracerGroup, forcingPool, nTracerGroup, &
                     nCellsSolve, maxLevelCell, nVertLevels, layerThickness, zMid, indexTemperature,    &
                     indexSalinity, tracerGroupTend, err)

                  call ocn_tracer_ecosys_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracerGroup, nCellsSolve, zMid, indexTemperature, indexSalinity, tracerGroupSurfaceFlux, err)
               endif

               !
               ! compute DMS source-sink tendencies and net surface fluxes
               ! NOTE: must be called before ocn_tracer_surface_flux_tend
               !
               if ( trim(groupItr % memberName) == 'DMSTracers' ) then
                  call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
                  nTracersEcosys = size(ecosysTracers, dim=1)
                  call ocn_tracer_DMS_compute(activeTracers, tracerGroup, nTracerGroup, ecosysTracers,   &
                     nTracersEcosys, forcingPool, nCellsSolve, maxLevelCell,  &
                     nVertLevels, layerThickness, indexTemperature, indexSalinity, tracerGroupTend, err)

                  call ocn_tracer_DMS_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracerGroup, nCellsSolve, zMid, indexTemperature, indexSalinity, tracerGroupSurfaceFlux, err)
               endif

               !
               ! compute MacroMolecules source-sink tendencies and net surface fluxes
               ! NOTE: must be called before ocn_tracer_surface_flux_tend
               !
               if ( trim(groupItr % memberName) == 'MacroMoleculesTracers' ) then
                  call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
                  nTracersEcosys = size(ecosysTracers, dim=1)
                  call ocn_tracer_MacroMolecules_compute(tracerGroup, nTracerGroup, ecosysTracers, nTracersEcosys, forcingPool, &
                     nCellsSolve, maxLevelCell, nVertLevels, layerThickness,  &
                     tracerGroupTend, err)

                  call ocn_tracer_MacroMolecules_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracerGroup, nCellsSolve, zMid, indexTemperature, indexSalinity, tracerGroupSurfaceFlux, err)
               endif

               !
               ! ocean surface restoring
               !
                 if (config_use_tracerGroup_surface_restoring) then
                     call mpas_timer_start("surface_restoring_" // trim(groupItr % memberName))
                     call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "PistonVelocity"
                     call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupPistonVelocity)
                     modifiedGroupName = trim(groupItr % memberName) // "SurfaceRestoringValue"
                     call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupSurfaceRestoringValue)
!maltrud monthly surface salinity restoring is a special case
                     call mpas_pool_get_config(ocnConfigs, 'config_use_surface_salinity_monthly_restoring',  &
                                                            config_use_surface_salinity_monthly_restoring)
!  within ocn_tracer_surface_restoring_compute makes sure this is only applied to salinity (using indexSalinity)
!                    if (trim(groupItr % memberName) /= 'activeTracers')   &
!                       config_use_surface_salinity_monthly_restoring = .false.
                     call ocn_tracer_surface_restoring_compute(groupItr % memberName, nTracerGroup, nCells, tracerGroup,  &
                        tracerGroupPistonVelocity, tracerGroupSurfaceRestoringValue, tracerGroupSurfaceFlux, indexSalinity, &
                        config_use_surface_salinity_monthly_restoring, err)
                     call mpas_timer_stop("surface_restoring_" // trim(groupItr % memberName))
                 endif

               ! tracer fluxes at the land-ice / ocean interface
               ! this is a flux at the top ocean surface -- so these fluxes are added into tracerGroupSurfaceFlux
               call ocn_surface_land_ice_fluxes_tracers(meshPool, groupItr % memberName, forcingPool, tracerGroupSurfaceFlux, err)

               !
               ! other additions to tracerGroupSurfaceFlux should be added here
               !

               !
               ! now begin to accumulate the RHS tracer tendencies.
               !

               !
               ! interior restoring forcing tendency
               !
                if (config_use_tracerGroup_interior_restoring) then
                     call mpas_timer_start("interior_restoring_" // trim(groupItr % memberName), .false.)
                     call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "InteriorRestoringRate"
                     call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupInteriorRestoringRate)
                     modifiedGroupName = trim(groupItr % memberName) // "InteriorRestoringValue"
                     call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupInteriorRestoringValue)
                     call ocn_tracer_interior_restoring_compute(nTracerGroup, nCells, maxLevelCell, layerThickness, &
                        tracerGroup, tracerGroupInteriorRestoringRate, tracerGroupInteriorRestoringValue, tracerGroupTend, err)
                     call mpas_timer_stop("interior_restoring_" // trim(groupItr % memberName))
                endif

               !
               ! exponential decay tendency
               !
                if (config_use_tracerGroup_exponential_decay) then
                     write (stderrUnit,'(a)') 'exponential decay not fully tested'
                     call mpas_timer_start("exponential decay " // trim(groupItr % memberName))
                     call mpas_pool_get_subpool(forcingPool, 'tracersExponentialDecayFields', tracersExponentialDecayFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "ExponentialDecayRate"
                     call mpas_pool_get_array(tracersExponentialDecayFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupExponentialDecayRate)
                     call ocn_tracer_exponential_decay_compute(nTracerGroup, nCells, maxLevelCell, layerThickness, &
                        tracerGroup, tracerGroupExponentialDecayRate, tracerGroupTend, err)

                     call mpas_timer_stop("exponential decay " // trim(groupItr % memberName))
                endif

               !
               ! ideal age forcing tendency
               !   note: ocn_tracer_ideal_age_compute resets tracers in top layer to zero
               !
                if (config_use_tracerGroup_idealAge_forcing) then
                     call mpas_timer_start("ideal age " // trim(groupItr % memberName))
                     write (stderrUnit,'(a)') 'ideal age not fully tested'
                     call mpas_pool_get_subpool(forcingPool, 'tracersIdealAgeFields', tracersIdealAgeFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "IdealAgeMask"
                     call mpas_pool_get_array(tracersIdealAgeFieldsPool, trim(modifiedGroupName), tracerGroupIdealAgeMask)
                     call ocn_tracer_ideal_age_compute(nTracerGroup, nCells, maxLevelCell, layerThickness, &
                          tracerGroupIdealAgeMask, tracerGroup, tracerGroupTend, err)

                     call mpas_timer_stop("ideal age " // trim(groupItr % memberName))
                endif

               !
               ! transit-time distribution (TTD) forcing tendency
               !   note: no tendency is actually computed in ocn_tracer_TTD_compute
               !   note: rather, tracerGroup is reset to tracerGroupTTDMask in top-most layer
               !
                if (config_use_tracerGroup_ttd_forcing) then
                     call mpas_timer_start("TTD " // trim(groupItr % memberName))
                     write (stderrUnit,'(a)') 'TTD not fully tested'
                     call mpas_pool_get_subpool(forcingPool, 'tracersTTDFields', tracersTTDFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "TTDMask"
                     call mpas_pool_get_array(tracersTTDFieldsPool, trim(modifiedGroupName), tracerGroupTTDMask)
                     call ocn_tracer_TTD_compute(nTracerGroup, nCells, maxLevelCell, layerThickness, &
                          tracerGroupTTDMask, tracerGroup, err)

                     call mpas_timer_stop("TTD " // trim(groupItr % memberName))
                endif

               !
               ! tracer tendency: horizontal advection term -div( layerThickness \phi u)
               !

               ! Monotonoic Advection, or standard advection
               call ocn_tracer_advection_tend(tracerGroup, normalThicknessFlux, vertAleTransportTop, layerThickness, &
                                              layerThickness, dt, meshPool, scratchPool, tend_layerThickness, tracerGroupTend)

               !
               ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
               !
               call ocn_tracer_hmix_tend(meshPool, scratchPool, layerThicknessEdge, zMid, tracerGroup, &
                                         relativeSlopeTopOfEdge, relativeSlopeTapering, relativeSlopeTaperingCell, &
                                         tracerGroupTend, err)

               !
               ! convert the surface tracer flux into a tracer tendency by distributing the flux across some number
               ! of surface layers
               !
               call ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, fractionAbsorbedRunoff, layerThickness, &
                                                 tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff,  &
                                                 tracerGroupTend, err)

               !
               ! Performing shortwave absorption
               !
               if ( trim(groupItr % memberName) == 'activeTracers' ) then
                  call ocn_tracer_short_wave_absorption_tend(meshPool, swForcingPool, forcingPool, indexTemperature, &
                            layerThickness, penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tracerGroupTend, err)
               endif

               !
               ! Compute tracer tendency due to non-local flux computed in KPP
               !
               if (config_use_cvmix_kpp) then
                 call mpas_timer_start("non-local flux from KPP")
                 call ocn_compute_KPP_input_fields(statePool, forcingPool, meshPool, diagnosticsPool, scratchPool, timeLevel)
                 if( trim(groupItr % memberName) == 'activeTracers' ) then
                     call ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, nonLocalSurfaceTracerFlux, tracerGroupTend, err)
                 else
                     call ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, tracerGroupSurfaceFlux, tracerGroupTend, err)
                 endif
                 call mpas_timer_stop("non-local flux from KPP")
               end if

               !
               ! Compute tracer tendency due to production/destruction of frazil ice
               !
               call ocn_frazil_forcing_tracers(meshPool, tracersPool, groupItr%memberName, forcingPool, tracerGroupTend, err)
            end if
         end if
      end do

      call mpas_threading_barrier()
      call mpas_deallocate_scratch_field(normalThicknessFluxField, .true.)

      call mpas_timer_stop("ocn_tend_tracer")

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

      !
      ! Low Frequency Divergence and high frequency thickness Tendency
      !

      ! Convert restore time from days to seconds
      thickness_filter_timescale_sec = config_thickness_filter_timescale*86400.0_RKIND
      highFreqThick_restore_time_sec = config_highFreqThick_restore_time*86400.0_RKIND

      allocate(div_hu(nVertLevels))

      !$omp do schedule(runtime) private(div_hu_btr, invAreaCell, i, iEdge, k, totalThickness, flux)
      do iCell = 1, nCells
        tend_lowFreqDivergence(:, iCell) = 0.0_RKIND
        tend_highFreqThickness(:, iCell) = 0.0_RKIND
        div_hu(:) = 0.0_RKIND
        div_hu_btr = 0.0_RKIND
        invAreaCell = 1.0_RKIND / areaCell(iCell)

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
      !$omp end do

      deallocate(div_hu)

      !
      !  high frequency thickness tendency: del2 horizontal hhf diffusion, div(\kappa_{hf} \nabla h^{hf})
      !
      call ocn_high_freq_thickness_hmix_del2_tend(meshPool, highFreqThickness, tend_highFreqThickness, err)

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
