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
   use ocn_diagnostics_variables
   use ocn_constants
   use ocn_config
   use ocn_mesh

   use ocn_surface_bulk_forcing
   use ocn_surface_land_ice_fluxes
   use ocn_frazil_forcing
   use ocn_tidal_forcing

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

   use ocn_vel_hadv_coriolis
   use ocn_vel_pressure_grad
   use ocn_vel_vadv
   use ocn_vel_hmix
   use ocn_vel_forcing
   use ocn_vmix
   use ocn_wetting_drying
   use ocn_vel_tidal_potential

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

   subroutine ocn_tend_thick(tendPool, forcingPool, meshPool)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool            !< [in] Mesh information (for pass-thru)

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         tendPool,          &!< [inout] Pool with accumulated tendencies
         forcingPool         !< [inout] Forcing information

      ! pointers for the actual output variables within the pools above
      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         surfaceThicknessFlux,     &! surface thickness flux
         surfaceThicknessFluxRunoff ! surface thickness flux from runoff

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         tendThick,        &! accumulated layer thickness tendency
         fractionAbsorbed, &! fraction of sfc flux absorbed
         fractionAbsorbedRunoff ! fraction of runoff flux absorbed

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::   &
         err,      &! internal error flag
         k, iCell   ! loop iterators for vertical, cell loops

      ! End preamble
      !-------------
      ! Begin code

      ! retrieve variables from pools
      call mpas_pool_get_array(tendPool,    'layerThickness', &
                                             tendThick)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFlux', &
                                             surfaceThicknessFlux)
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxRunoff', &
                                             surfaceThicknessFluxRunoff)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', &
                                             fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', &
                                             fractionAbsorbedRunoff)

      !
      ! layer thickness tendency: 
      !    initialize to zero and start accumulating tendency terms
      !

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCellsAll
         surfaceThicknessFlux(iCell) = 0.0_RKIND
         surfaceThicknessFluxRunoff(iCell) = 0.0_RKIND
         do k=1,nVertLevels
            tendThick(k, iCell) = 0.0_RKIND
         end do
      end do
      !$omp end do
      !$omp end parallel

      ! If turned off, return with zero fluxes, tendencies
      ! Otherwise, start time and call routines to accumulate
      if (config_disable_thick_all_tend) return
      call mpas_timer_start("ocn_tend_thick")

      ! Compute surface mass flux array from bulk forcing
      call ocn_surface_bulk_forcing_thick(meshPool, forcingPool, &
                                          surfaceThicknessFlux, &
                                          surfaceThicknessFluxRunoff, err)

      ! Compute surface thickness flux from land ice
      call ocn_surface_land_ice_fluxes_thick(meshPool, forcingPool, &
                                             surfaceThicknessFlux, err)

      !
      ! Compute horizontal advection term -\nabla\cdot ( hu)
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      ! for explanation of divergence operator.
      !
      call ocn_thick_hadv_tend(meshPool, normalTransportVelocity, &
                               layerThickEdge, tendThick, err)

      ! Compute vertical advection term -d/dz(hw)
      call ocn_thick_vadv_tend(meshPool, vertAleTransportTop, &
                               tendThick, err)

      ! Compute surface flux tendency
      call ocn_thick_surface_flux_tend(meshPool, fractionAbsorbed, &
                                       fractionAbsorbedRunoff, &
                                       surfaceThicknessFlux, &
                                       surfaceThicknessFluxRunoff, &
                                       tendThick, err)

      ! Compute contribution from frazil ice formation
      call ocn_frazil_forcing_layer_thickness(meshPool, forcingPool, &
                                              tendThick, err)

      ! Compute thickness change due to tidal forcing
      call ocn_tidal_forcing_layer_thickness(meshPool, forcingPool, &
                                             tendThick, err)

      call mpas_timer_stop("ocn_tend_thick")

   !--------------------------------------------------------------------

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

   subroutine ocn_tend_vel(tendPool, statePool, forcingPool, &
                           timeLevelIn, dt)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt             !< [in] time step

      integer, intent(in), optional :: &
         timeLevelIn    !< [in] Time level for state fields

      ! Input structures with lots of fields
      type (mpas_pool_type), intent(in) :: statePool ! replace

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      ! currently updated as pool variables - change to args?

      type (mpas_pool_type), intent(inout) :: &
         tendPool,          &!< [out] Tendency structure w/ vel tend
         forcingPool         !< [out] Forcing structure w/ sfc stresses

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         tendVel             !< [out] normal velocity tendency at edges

      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         sfcStress,         &!< [out] surface stress at edges
         sfcStressMag        !< [out] surface stress magnitude (on cell)

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::           &
         err,              &! local error flag
         iEdge, iCell, k,  &! loop indices for edge, cell and vertical
         indxTemp,         &! tracer index for temperature
         indxSalt,         &! tracer index for salinity
         timeLevel          ! time level to use for state variables

      ! various pointers for array retrievals

      integer, pointer :: indexTemperature, indexSalinity
      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         ssh               ! sea surface height

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         normalVelocity,     &! normal velocity
         layerThickness, &
         circulation

      real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
         activeTracers

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Initialize error code and determine time level

      err = 0
      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !*** Retrieve pool variables
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(statePool,   'normalVelocity', &
                                             normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool,   'layerThickness', &
                                             layerThickness, timeLevel)
      call mpas_pool_get_array(statePool,   'ssh', ssh, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                             activeTracers, timeLevel)

      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSalinity)
      indxTemp = indexTemperature
      indxSalt = indexSalinity

      call mpas_pool_get_array(tendPool, 'normalVelocity', tendVel)

      call mpas_pool_get_array(forcingPool, 'surfaceStress', &
                                             sfcStress)
      call mpas_pool_get_array(forcingPool, 'surfaceStressMagnitude', &
                                             sfcStressMag)

      !*** Transfer data to device
      !$acc enter data &
      !$acc    copyin(tendVel, sfcStress, sfcStressMag, &
      !$acc    sfcFlxAttCoeff, ssh, layerThickEdge, normalVelocity, &
      !$acc    tangentialVelocity, density, potentialDensity, zMid, &
      !$acc    pressure, layerThickness, circulation, activeTracers, &
      !$acc    relativeVorticity, kineticEnergyCell, normRelVortEdge, &
      !$acc    normPlanetVortEdge, montgomeryPotential, &
      !$acc    vertAleTransportTop, divergence, vertViscTopOfEdge, &
      !$acc    thermExpCoeff, salineContractCoeff, wettingVelocity)

      !*** Set output variables to zero before accumulating results
      
#ifdef MPAS_OPENACC
      !$acc parallel loop present(sfcStress, tendVel) private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1, nEdgesAll
         sfcStress(iEdge) = 0.0_RKIND
         do k=1,nVertLevels
            tendVel(k,iEdge) = 0.0_RKIND
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc parallel loop present(sfcStressMag)
#else
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         sfcStressMag(iCell) = 0.0_RKIND
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !*** Return early if all vel tendencies disabled
      !*** Otherwise, start timer

      if (config_disable_vel_all_tend) return
      call mpas_timer_start("ocn_tend_vel")

      ! Compute bulk forcing surface stress
      call ocn_surface_bulk_forcing_vel( &
                                    forcingPool, sfcStress, &
                                    sfcStressMag, err)

      ! Add top drag to surface stress
      call ocn_surface_land_ice_fluxes_vel( &
                                       sfcStress, sfcStressMag, err)

      ! Add nonlinear Coriolis term and horizontal advection of
      ! momentum, formulated as grad of kinetic energy
      call ocn_vel_hadv_coriolis_tend(normRelVortEdge, &
                                      normPlanetVortEdge, &
                                      layerThickEdge, normalVelocity, &
                                      kineticEnergyCell, tendVel, err)

      ! Add vertical advection term -w du/dz
      call ocn_vel_vadv_tend(normalVelocity, layerThickEdge, &
                             vertAleTransportTop, tendVel, err)

      ! Add pressure gradient
      call ocn_vel_pressure_grad_tend(ssh, pressure, &
                               montgomeryPotential, zMid, &
                               density, potentialDensity, &
                               indxTemp, indxSalt, activeTracers, &
                               thermExpCoeff, salineContractCoeff, &
                               tendVel, err)

      ! Add tidal potential (if needed) 
      call ocn_compute_tidal_potential_forcing(err)
      if (config_time_integrator == 'RK4') then
        ! for split explicit, tidal forcing is added in barotropic subcycles
        call ocn_vel_tidal_potential_tend(ssh, tendVel, err) 
      endif

      ! Add horizontal mixing
      call ocn_vel_hmix_tend(divergence, relativeVorticity, tendVel,err)

      ! Add forcing and bottom drag
      call ocn_vel_forcing_tend(normalVelocity, sfcFlxAttCoeff, &
                                sfcStress, kineticEnergyCell, &
                                layerThickEdge, tendVel, err)

      ! vertical mixing treated implicitly in a later routine
      ! adjust total velocity tendency based on wetting/drying

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) &
      !$acc    present(tendVel, wettingVelocity)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1, nEdgesAll
      do k=1,nVertLevels
         tendVel(k,iEdge) = tendVel(k,iEdge)* &
                            (1.0_RKIND - wettingVelocity(k,iEdge))
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !*** Transfer any needed data back to host and delete inputs
      !$acc exit data &
      !$acc    copyout(tendVel, sfcStress, sfcStressMag) &
      !$acc    delete( &
      !$acc    sfcFlxAttCoeff, ssh, layerThickEdge, normalVelocity, &
      !$acc    tangentialVelocity, density, potentialDensity, zMid, &
      !$acc    pressure, layerThickness, circulation, activeTracers, &
      !$acc    relativeVorticity, kineticEnergyCell, normRelVortEdge, &
      !$acc    normPlanetVortEdge, montgomeryPotential, &
      !$acc    vertAleTransportTop, divergence, vertViscTopOfEdge, &
      !$acc    thermExpCoeff, salineContractCoeff, wettingVelocity)

      ! stop the timer and exit

      call mpas_timer_stop("ocn_tend_vel")

   !--------------------------------------------------------------------

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

   subroutine ocn_tend_tracer(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, &
                              dt, activeTracersOnlyIn, timeLevelIn )!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool,          &!< [in] ocean state variables
         meshPool             !< [in] mesh information

      real (kind=RKIND), intent(in) :: &
         dt                   !< [in] time step (seconds)

      logical, intent(in), optional :: &
         activeTracersOnlyIn  !< [in] only compute for active tracers

      integer, intent(in), optional :: &
         timeLevelIn          !< [in] time index to use for state vars

      ! these variables input from shared diagnostics module:
      !    normalTransportVelocity - transport velocity across edge
      !    layerThickEdge          - layer thickness on edge

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         forcingPool,   &!< [inout] forcing data
         swForcingPool, &!< [inout] shortwave forcing data
         tendPool        !< [inout] Tendency terms for all variables

      integer :: &
         iCell, iEdge, k, n, &! loop counters for cell,edge,vert,tracer
         err,                &! internal error flag 
         timeLevel,          &! time level to use for state variables
         indxTemp, indxSalt, &! tracer index for temperature, salinity
         nTracersGroup,      &! number of tracers in each tracer group
         nTracersEcosys       ! number of ecosystem tracers

      logical :: &
         activeTracersOnly,  &! only compute tend for active tracers
         isActiveTracer,     &! is this group the active tracer group
         computeBudgets       ! compute active tracer budgets

      ! iterator and modified strings for tracer categories
      type (mpas_pool_iterator_type) :: groupItr
      character (len=StrKIND) :: &
         groupName,            &! tracer group name
         configBase,           &! base name for config options
         modifiedGroupName,    &! variations on group name for options
         modifiedConfigName     ! string to construct config var names

      ! sub-pools with various data sets
      type (mpas_pool_type), pointer :: &
         tracersPool,            &! tracer variables
         tracersTendPool,        &! tracer tendencies
         tracersSurfaceFluxPool, &! surface fluxes
         tracersSurfaceRestoringFieldsPool, &! surface restoring
         tracersInteriorRestoringFieldsPool,&! interior restoring
         tracersExponentialDecayFieldsPool, &! exponential decay
         tracersIdealAgeFieldsPool,         &! ideal age
         tracersTTDFieldsPool     ! transit time distribution

      ! pointers for retrieving data from pools
      integer, pointer :: &
         indexTemperature, indexSalinity ! tracer index for temp,salt

      logical, pointer :: & ! option configuration flags
         config_use_tracerGroup, &
         config_use_tracerGroup_surface_bulk_forcing, &
         config_use_tracerGroup_surface_restoring,    &
         config_use_tracerGroup_interior_restoring,   &
         config_use_tracerGroup_exponential_decay,    &
         config_use_tracerGroup_idealAge_forcing,     &
         config_use_tracerGroup_ttd_forcing

      real (kind=RKIND), pointer :: &
         salinity_restoring_constant_piston_velocity

      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         penetrativeTemperatureFlux, &! heat flux penetrating below sfc
         tracerGroupExponentialDecayRate ! exp decay rate for forcing

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         tracerGroupPistonVelocity,        &! forcing piston velocity
         tracerGroupSurfaceRestoringValue, &! restore sfc to this value
         tracerGroupIdealAgeMask,          &! mask to reset ideal age
         tracerGroupTTDMask                 ! mask for TTD values

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         layerThickness,               &! layer thickness
         tracerGroupSurfaceFlux,       &! tracer flux at surface
         fractionAbsorbed,             &! frac sfc flux aborbed
         fractionAbsorbedRunoff,       &! frac runoff flux aborbed
         tracerGroupSurfaceFluxRunoff, &! runoff flux
         tracerGroupSurfaceFluxRemoved,&! total sfc flux absorbed
         nonLocalSurfaceTracerFlux      ! non-local fluxes (eg KPP)

      real (kind=RKIND), dimension(:,:,:), pointer, contiguous :: &
         tracerGroup,     &! tracers in current tracer group
         tracerGroupTend, &! tendencies for current tracer group
         activeTracers,   &! retain active tracers
         ecosysTracers,   &! retain ecosys tracers for DMS, MacroMols
         tracerGroupInteriorRestoringRate, &! int restoring timescale 
         tracerGroupInteriorRestoringValue  ! int restoring value

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:), allocatable :: &
         normalThicknessFlux ! Flux of thickness through edge (m^2/s)

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !
      ! set time level of optional argument is present
      !
      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      if (present(activeTracersOnlyIn)) then
         activeTracersOnly = activeTracersOnlyIn
      else
         activeTracersOnly = .false.
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
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      !
      ! get arrays
      !
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', fractionAbsorbedRunoff)

      if(config_disable_tr_all_tend) return

      call mpas_timer_start("ocn_tend_tracer")

      allocate(normalThicknessFlux(nVertLevels, nEdgesAll+1))

      !
      ! transport velocity for the tracer.
      !
      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdgesAll
         do k = 1, nVertLevels
            normalThicknessFlux(k, iEdge) = normalTransportVelocity(k, iEdge) * layerThickEdge(k, iEdge)
         end do
      end do
      !$omp end do
      !$omp end parallel

      !
      ! begin iterate over tracer categories
      !
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
           ! Only compute tendencies for active tracers if activeTracersOnly flag is true.
           if ( .not.activeTracersOnly .or. trim(groupItr % memberName)=='activeTracers') then
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
               nTracersGroup = size(tracerGroup, dim=1)

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
!maltrud only active tracers have runoff flux for now, but we still need to associate for ALL tracers
               modifiedGroupName = trim(groupItr % memberName) // "SurfaceFluxRunoff"
               call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFluxRunoff)

               ! Get surface flux removed array to keep track of how much flux is ignored
               modifiedGroupName = trim(groupItr % memberName) // "SurfaceFluxRemoved"
               call mpas_pool_get_array(tracersSurfaceFluxPool, trim(modifiedGroupName), tracerGroupSurfaceFluxRemoved)

               !
               ! initialize tracer surface flux and tendency to zero.
               !
               !$omp parallel
               !$omp do schedule(runtime)
               do iCell = 1, nCellsAll
                 tracerGroupTend(:,:, iCell) = 0.0_RKIND
                 tracerGroupSurfaceFlux(:, iCell) = 0.0_RKIND
               end do
               !$omp end do
               !$omp end parallel

               !
               ! fill components of surface tracer flux
               !
               if (config_use_tracerGroup_surface_bulk_forcing) then
                  !$omp parallel
                  !$omp do schedule(runtime)
                  do iCell = 1, nCellsAll
                    tracerGroupSurfaceFluxRunoff(:, iCell) = 0.0_RKIND
                    tracerGroupSurfaceFluxRemoved(:, iCell) = 0.0_RKIND
                  end do
                  !$omp end do
                  !$omp end parallel

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
                     nCellsOwned, latCell, minLevelCell, maxLevelCell, nVertLevels, layerThickness, zMid, indexTemperature, &
                     indexSalinity, tracerGroupTend, err)

                  call ocn_tracer_ecosys_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracerGroup, nCellsOwned, zMid, minLevelCell, indexTemperature, indexSalinity, tracerGroupSurfaceFlux, err)
               endif

               !
               ! compute DMS source-sink tendencies and net surface fluxes
               ! NOTE: must be called before ocn_tracer_surface_flux_tend
               !
               if ( trim(groupItr % memberName) == 'DMSTracers' ) then
                  call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
                  nTracersEcosys = size(ecosysTracers, dim=1)
                  call ocn_tracer_DMS_compute(activeTracers, tracerGroup, nTracerGroup, ecosysTracers,   &
                     nTracersEcosys, forcingPool, nCellsOwned, minLevelCell, maxLevelCell,  &
                     nVertLevels, layerThickness, indexTemperature, indexSalinity, tracerGroupTend, err)

                  call ocn_tracer_DMS_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracerGroup, nCellsOwned, zMid, minLevelCell, indexTemperature, indexSalinity, tracerGroupSurfaceFlux,  &
                     tracerGroupSurfaceFluxRemoved, err)
               endif

               !
               ! compute MacroMolecules source-sink tendencies and net surface fluxes
               ! NOTE: must be called before ocn_tracer_surface_flux_tend
               !
               if ( trim(groupItr % memberName) == 'MacroMoleculesTracers' ) then
                  call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
                  nTracersEcosys = size(ecosysTracers, dim=1)
                  call ocn_tracer_MacroMolecules_compute(tracerGroup, nTracerGroup, ecosysTracers, nTracersEcosys, forcingPool, &
                     nCellsOwned, minLevelCell, maxLevelCell, nVertLevels, layerThickness,  &
                     tracerGroupTend, err)

                  call ocn_tracer_MacroMolecules_surface_flux_compute(activeTracers, tracerGroup, forcingPool,  &
                     nTracersGroup, nCellsOwned, zMid, indexTemperature, indexSalinity, tracerGroupSurfaceFlux, err)
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
                     ! Note: monthly surface salinity restoring is a special case for tracer restoring
                     call ocn_tracer_surface_restoring_compute(groupItr % memberName, nTracersGroup, nCellsAll, tracerGroup,  &
                        tracerGroupPistonVelocity, tracerGroupSurfaceRestoringValue, tracerGroupSurfaceFlux, indexSalinity, &
                        config_use_surface_salinity_monthly_restoring, &
                        config_salinity_restoring_constant_piston_velocity, err)
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
                     call ocn_tracer_interior_restoring_compute(nTracerGroup, nCellsAll, minLevelCell, maxLevelCell, layerThickness, &
                        tracerGroup, tracerGroupInteriorRestoringRate, tracerGroupInteriorRestoringValue, tracerGroupTend, err)
                     call mpas_timer_stop("interior_restoring_" // trim(groupItr % memberName))
                endif

               !
               ! exponential decay tendency
               !
                if (config_use_tracerGroup_exponential_decay) then
                     call mpas_log_write( &
                        "WARNING: exponential decay not fully tested", &
                        MPAS_LOG_WARN)
                     call mpas_timer_start("exponential decay " // trim(groupItr % memberName))
                     call mpas_pool_get_subpool(forcingPool, 'tracersExponentialDecayFields', tracersExponentialDecayFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "ExponentialDecayRate"
                     call mpas_pool_get_array(tracersExponentialDecayFieldsPool, trim(modifiedGroupName), &
                                              tracerGroupExponentialDecayRate)
                     call ocn_tracer_exponential_decay_compute(nTracerGroup, nCellsAll, minLevelCell, maxLevelCell, layerThickness, &
                        tracerGroup, tracerGroupExponentialDecayRate, tracerGroupTend, err)

                     call mpas_timer_stop("exponential decay " // trim(groupItr % memberName))
                endif

               !
               ! ideal age forcing tendency
               !   note: ocn_tracer_ideal_age_compute resets tracers in top layer to zero
               !
                if (config_use_tracerGroup_idealAge_forcing) then
                     call mpas_timer_start("ideal age " // trim(groupItr % memberName))
                     call mpas_log_write( &
                        "WARNING: ideal age not fully tested", &
                        MPAS_LOG_WARN)
                     call mpas_pool_get_subpool(forcingPool, 'tracersIdealAgeFields', tracersIdealAgeFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "IdealAgeMask"
                     call mpas_pool_get_array(tracersIdealAgeFieldsPool, trim(modifiedGroupName), tracerGroupIdealAgeMask)
                     call ocn_tracer_ideal_age_compute(nTracerGroup, nCellsAll, minLevelCell, maxLevelCell, layerThickness, &
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
                     call mpas_log_write( &
                        "WARNING: TTD not fully tested", &
                        MPAS_LOG_WARN)
                     call mpas_pool_get_subpool(forcingPool, 'tracersTTDFields', tracersTTDFieldsPool)
                     modifiedGroupName = trim(groupItr % memberName) // "TTDMask"
                     call mpas_pool_get_array(tracersTTDFieldsPool, trim(modifiedGroupName), tracerGroupTTDMask)
                     call ocn_tracer_TTD_compute(nTracersGroup, nCellsAll, maxLevelCell, layerThickness, &
                          tracerGroupTTDMask, tracerGroup, err)

                     call mpas_timer_stop("TTD " // trim(groupItr % memberName))
                endif

               !
               ! tracer tendency: horizontal advection term -div( layerThickness \phi u)
               !

               ! Monotonoic Advection, or standard advection
               ! Tendency for tracer budget is stored within the tracer adv
               ! routine
               call ocn_tracer_advection_tend(tracerGroup, normalThicknessFlux, vertAleTransportTop, layerThickness, &
                                              dt, meshPool, tracerGroupTend, &
                                              trim(groupItr % memberName))

               !
               ! tracer tendency: del2 horizontal tracer diffusion, div(h \kappa_2 \nabla \phi)
               !

               isActiveTracer = .false.
               if(trim(groupItr % memberName)=='activeTracers') then
                 isActiveTracer = .true.
                 if (config_compute_active_tracer_budgets) then
                    !$omp parallel
                    !$omp do schedule(runtime)
                    do iCell = 1, nCellsAll
                      activeTracerHorMixTendency(:,:,iCell) = tracerGroupTend(:,:,iCell)
                    enddo
                    !$omp end do
                    !$omp end parallel
                 endif
               endif
               call ocn_tracer_hmix_tend(meshPool, layerThickEdge, layerThickness, zMid, tracerGroup, &
                                         RediKappa, slopeTriadUp, slopeTriadDown, dt, isActiveTracer, &
                                         RediKappaSfcTaper, RediKappaScaling, rediLimiterCount, tracerGroupTend, err)

               if(trim(groupItr % memberName)=='activeTracers') then
                 isActiveTracer = .true.
                 if (config_compute_active_tracer_budgets) then
                    !$omp parallel
                    !$omp do schedule(runtime)
                    do iCell = 1, nCellsAll
                      activeTracerHorMixTendency(:,:,iCell) = tracerGroupTend(:,:,iCell) - &
                                  activeTracerHorMixTendency(:,:,iCell)
                    enddo
                    !$omp end do
                    !$omp end parallel
                 endif
               endif
 
               !
               ! convert the surface tracer flux into a tracer tendency by distributing the flux across some number
               ! of surface layers
               !
               if (config_compute_active_tracer_budgets) then
                  if ( trim(groupItr % memberName) == 'activeTracers' ) then
                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        activeTracerSurfaceFluxTendency(:,:,iCell) = tracerGroupTend(:,:,iCell)
                     end do
                     !$omp end do
                     !$omp end parallel
                  endif
               endif

               call ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, fractionAbsorbedRunoff, layerThickness, &
                                                 tracerGroupSurfaceFlux, tracerGroupSurfaceFluxRunoff,  &
                                                 tracerGroupTend, err)

               !
               ! Performing shortwave absorption
               !
               if ( trim(groupItr % memberName) == 'activeTracers' ) then

                  if (config_compute_active_tracer_budgets) then
                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        activeTracerSurfaceFluxTendency(:,:,iCell) = tracerGroupTend(:,:,iCell) &
                           - activeTracerSurfaceFluxTendency(:,:,iCell)
                        temperatureShortWaveTendency(:,iCell) = tracerGroupTend(indexTemperature,:,iCell)
                     end do
                     !$omp end do
                     !$omp end parallel
                  endif

                  call ocn_tracer_short_wave_absorption_tend(meshPool, swForcingPool, forcingPool, indexTemperature, &
                           layerThickness, penetrativeTemperatureFlux, penetrativeTemperatureFluxOBL, tracerGroupTend, &
                           bottomLayerShortwaveTemperatureFlux, err)

                  if (config_compute_active_tracer_budgets) then
                     !$omp parallel
                     !$omp do schedule(runtime)
                     do iCell = 1, nCellsAll
                        temperatureShortWaveTendency(:,iCell) = tracerGroupTend(indexTemperature,:,iCell) - &
                                        temperatureShortWaveTendency(:,iCell)
                     end do
                     !$omp end do
                     !$omp end parallel
                  endif
               endif

               !
               ! Compute tracer tendency due to non-local flux computed in KPP
               !
               if (config_use_cvmix_kpp) then
                  call mpas_timer_start("non-local flux from KPP")
                  call ocn_compute_KPP_input_fields(statePool, forcingPool, meshPool, timeLevel)
                  if (.not. config_cvmix_kpp_nonlocal_with_implicit_mix) then
                     if( trim(groupItr % memberName) == 'activeTracers' ) then
                        if (config_compute_active_tracer_budgets) then
                           !$omp parallel
                           !$omp do schedule(runtime)
                           do iCell = 1, nCellsAll
                              activeTracerNonLocalTendency(:,:,iCell) = tracerGroupTend(:,:,iCell)
                           end do
                           !$omp end do
                           !$omp end parallel
                        endif
                        call ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, nonLocalSurfaceTracerFlux, tracerGroupTend, err)
                        if (config_compute_active_tracer_budgets) then
                           !$omp parallel
                           !$omp do schedule(runtime)
                           do iCell = 1, nCellsAll
                           do k = 1, nVertLevels
                           do n = 1, nTracersGroup
                              activeTracerNonLocalTendency(n,k,iCell) = &
                                           tracerGroupTend(n,k,iCell) - &
                              activeTracerNonLocalTendency(n,k,iCell)
                           end do
                           end do
                           end do
                           !$omp end do
                           !$omp end parallel
                        endif
                     else
                        call ocn_tracer_nonlocalflux_tend(meshPool, vertNonLocalFlux, tracerGroupSurfaceFlux, tracerGroupTend, err)
                     endif
                  end if
                  call mpas_timer_stop("non-local flux from KPP")
               end if

               !
               ! Compute tracer tendency due to production/destruction of frazil ice
               !
               call ocn_frazil_forcing_tracers(meshPool, tracersPool, groupItr%memberName, forcingPool, tracerGroupTend, err)
            end if
           end if ! active only
         end if
      end do

      deallocate(normalThicknessFlux)

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
   subroutine ocn_tend_freq_filtered_thickness(tendPool, statePool, meshPool, timeLevelIn)!{{{

      type (mpas_pool_type), intent(inout) :: tendPool !< Input/Output: Tendency information
      type (mpas_pool_type), intent(in) :: statePool !< Input: State information
      type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
      integer, intent(in), optional :: timeLevelIn !< Input: Time level for state fields

      integer :: timeLevel
      integer :: err, iCell, i, k, iEdge
      integer, pointer :: nCells

      real (kind=RKIND) :: flux, invAreaCell, div_hu_btr, thickness_filter_timescale_sec, highFreqThick_restore_time_sec, &
         totalThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, &
         layerThickness, &
         lowFreqDivergence, highFreqThickness, &
         tend_lowFreqDivergence, tend_highFreqThickness
      real (kind=RKIND), dimension(:), allocatable:: div_hu

      call mpas_timer_start("ocn_tend_freq_filtered_thickness")

      err = 0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', lowFreqDivergence, timeLevel)
      call mpas_pool_get_array(statePool, 'highFreqThickness', highFreqThickness, timeLevel)

      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', tend_lowFreqDivergence)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', tend_highFreqThickness)

      !
      ! Low Frequency Divergence and high frequency thickness Tendency
      !

      ! Convert restore time from days to seconds
      thickness_filter_timescale_sec = config_thickness_filter_timescale*86400.0_RKIND
      highFreqThick_restore_time_sec = config_highFreqThick_restore_time*86400.0_RKIND

      allocate(div_hu(nVertLevels))

      !$omp parallel
      !$omp do schedule(runtime) private(div_hu_btr, invAreaCell, i, iEdge, k, flux, totalThickness) &
      !$omp firstprivate(div_hu)
      do iCell = 1, nCellsAll
        tend_lowFreqDivergence(:, iCell) = 0.0_RKIND
        tend_highFreqThickness(:, iCell) = 0.0_RKIND
        div_hu(:) = 0.0_RKIND
        div_hu_btr = 0.0_RKIND
        invAreaCell = 1.0_RKIND / areaCell(iCell)

        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)

          ! Note: changing loop limit to minLevelEdgeBot introduces non-bit-for-bit changes
          do k = 1, maxLevelEdgeTop(iEdge)
            flux = layerThickEdge(k, iEdge) * normalVelocity(k, iEdge) * dvEdge(iEdge) * edgeSignOnCell(i, iCell) * invAreaCell
            div_hu(k) = div_hu(k) - flux
            div_hu_btr = div_hu_btr - flux
          end do
        end do

        totalThickness = sum(layerThickness(minLevelCell(iCell):maxLevelCell(iCell),iCell))
        do k = minLevelCell(iCell), maxLevelCell(iCell)

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
      !$omp end parallel

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
