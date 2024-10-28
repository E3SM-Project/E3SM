! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use mpas_timekeeping
   use mpas_threading
   use mpas_dmpar
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
   use ocn_tracer_CFC

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
   use ocn_vel_self_attraction_loading
   use ocn_manufactured_solution

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

   subroutine ocn_tend_thick(tendPool, forcingPool)!{{{

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         tendPool,          &!< [inout] Pool with accumulated tendencies
         forcingPool         !< [inout] Forcing information

      ! pointers for the actual output variables within the pools above
      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         surfaceThicknessFlux,     &! surface thickness flux
         surfaceThicknessFluxRunoff, &! surface thickness flux from runoff
         surfaceThicknessFluxSubglacialRunoff ! surface thickness flux from runoff

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         tendThick,        &! accumulated layer thickness tendency
         fractionAbsorbed, &! fraction of sfc flux absorbed
         fractionAbsorbedRunoff, &! fraction of sfc flux absorbed
         fractionAbsorbedSubglacialRunoff ! fraction of subglacial runoff flux absorbed

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
      call mpas_pool_get_array(forcingPool, 'surfaceThicknessFluxSubglacialRunoff', &
                                             surfaceThicknessFluxSubglacialRunoff)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', &
                                             fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', &
                                             fractionAbsorbedRunoff)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedSubglacialRunoff', &
                                             fractionAbsorbedSubglacialRunoff)



      !
      ! layer thickness tendency: 
      !    initialize to zero and start accumulating tendency terms
      !
#ifdef MPAS_OPENACC
      !$acc enter data create(tendThick, surfaceThicknessFlux, surfaceThicknessFluxRunoff)

      !$acc parallel loop &
      !$acc     present(tendThick, surfaceThicknessFlux, surfaceThicknessFluxRunoff) &
      !$acc     private(k)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCellsAll
         surfaceThicknessFlux(iCell) = 0.0_RKIND
         surfaceThicknessFluxRunoff(iCell) = 0.0_RKIND
         do k=1,nVertLevels
            tendThick(k, iCell) = 0.0_RKIND
         end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      if (trim(config_subglacial_runoff_mode) == 'data') then
#ifdef MPAS_OPENACC
      !$acc enter data create(surfaceThicknessFluxSubglacialRunoff)

      !$acc parallel loop &
      !$acc     present(surfaceThicknessFluxSubglacialRunoff)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
         do iCell = 1, nCellsAll
            surfaceThicknessFluxSubglacialRunoff(iCell) = 0.0_RKIND
         end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif
      end if

      ! If turned off, return with zero fluxes, tendencies
      ! Otherwise, start time and call routines to accumulate
      if (config_disable_thick_all_tend) return
      call mpas_timer_start("ocn_tend_thick")

      ! Compute surface mass flux array from bulk forcing
      call ocn_surface_bulk_forcing_thick(forcingPool, &
                                          surfaceThicknessFlux, &
                                          surfaceThicknessFluxRunoff, err)

      if (trim(config_subglacial_runoff_mode) == 'data') then
            call ocn_surface_bulk_forcing_thick_subglacial_runoff(forcingPool, &
                                          surfaceThicknessFluxSubglacialRunoff, err)
      end if

      ! Compute surface thickness flux from land ice
      call ocn_surface_land_ice_fluxes_thick(forcingPool, &
                                             surfaceThicknessFlux, err)

      !
      ! Compute horizontal advection term -\nabla\cdot ( hu)
      ! See Ringler et al. (2010) jcp paper, eqn 19, 21, and fig. 3.
      ! for explanation of divergence operator.
      !
      call ocn_thick_hadv_tend(normalTransportVelocity, &
                               layerThickEdgeFlux, tendThick, err)

      ! Compute vertical advection term -d/dz(hw)
      call ocn_thick_vadv_tend(vertAleTransportTop, &
                               tendThick, err)

      ! Compute surface flux tendency
      call ocn_thick_surface_flux_tend(fractionAbsorbed, &
                                       fractionAbsorbedRunoff, &
                                       surfaceThicknessFlux, &
                                       surfaceThicknessFluxRunoff, &
                                       tendThick, err)
      if (trim(config_subglacial_runoff_mode) == 'data') then
            call ocn_thick_surface_flux_tend_subglacial_runoff(fractionAbsorbedSubglacialRunoff, &
                                       surfaceThicknessFluxSubglacialRunoff, &
                                       tendThick, err)
      end if

      ! Compute contribution from frazil ice formation
      call ocn_frazil_forcing_layer_thickness(forcingPool, &
                                              tendThick, err)

      ! Compute thickness change due to tidal forcing
      call ocn_tidal_forcing_layer_thickness(forcingPool, &
                                             tendThick, err)

      ! Compute thickness tendency to manufactured solution
      call ocn_manufactured_solution_tend_thick(tendThick, err)

#ifdef MPAS_OPENACC
      !$acc exit data copyout(tendThick, surfaceThicknessFlux, surfaceThicknessFluxRunoff, &
      !$acc      surfaceThicknessFluxSubglacialRunoff)
#endif

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
    ! NOTE: Domain being passed in for halo update on SAL
   subroutine ocn_tend_vel(domain, tendPool, statePool, forcingPool, &
                           timeLevelIn, dminfo, dt)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

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
      type (dm_info), intent(in) :: dminfo

      real (kind=RKIND), dimension(:), pointer, contiguous :: &
         ssh               ! sea surface height

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         normalVelocity,     &! normal velocity
         layerThickness

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
      !$acc update device(zMid, divergence, relativeVorticity, normPlanetVortEdge, &
      !$acc               density, normRelVortEdge, montgomeryPotential, pressure, &
      !$acc               thermExpCoeff, salineContractCoeff, tangentialVelocity, &
      !$acc               layerThickEdgeFlux, kineticEnergyCell, sfcFlxAttCoeff, potentialDensity, &
      !$acc               vertAleTransportTop, vertViscTopOfEdge, wettingVelocityFactor)
      !$acc enter data &
      !$acc    copyin(tendVel, sfcStress, sfcStressMag, &
      !$acc    ssh, normalVelocity, &
      !$acc    layerThickness, activeTracers)

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
                                      layerThickEdgeFlux, normalVelocity, &
                                      kineticEnergyCell, tendVel, err)

      ! Add vertical advection term -w du/dz
      call ocn_vel_vadv_tend(normalVelocity, layerThickEdgeFlux, &
                             vertAleTransportTop, tendVel, err)

      ! Add pressure gradient
      call ocn_vel_pressure_grad_tend(ssh, pressure, surfacePressure, &
                               montgomeryPotential, zMid, &
                               density, potentialDensity, &
                               indxTemp, indxSalt, activeTracers, &
                               thermExpCoeff, salineContractCoeff, &
                               tendVel, err)

      ! Compute pgf additions due to self-attraction and loading
      if(mpas_is_alarm_ringing(domain % clock, 'salComputeAlarm', ierr=err)) then 
#ifdef MPAS_DEBUG
          call mpas_log_write('       Computing SAL')
#endif
          call ocn_compute_self_attraction_loading(domain, forcingPool, dminfo, &
                                                   ssh, surfacePressure, err)
          call mpas_reset_clock_alarm(domain % clock, 'salComputeAlarm', ierr=err)
      endif

      ! Add tidal potential (if needed) 
      call ocn_compute_tidal_potential_forcing(err)
      if ( (config_time_integrator == 'RK4') &
           .or. (config_time_integrator =='LTS') &
           .or. (config_time_integrator == 'FB_LTS') ) then
        ! for split explicit, tidal forcing is added in barotropic subcycles
        call ocn_vel_tidal_potential_tend(ssh,surfacePressure, tendVel, err)
      endif

      ! Add horizontal mixing
      call ocn_vel_hmix_tend(divergence, relativeVorticity, tendVel,err)

      ! Add forcing and bottom drag
      call ocn_vel_forcing_tend(normalVelocity, sfcFlxAttCoeff, &
                                sfcStress, kineticEnergyCell, &
                                layerThickEdgeDrag, &
                                layerThickEdgeMean, tendVel, err)
 
      ! Compute tendency term for manufactured solution
      call ocn_manufactured_solution_tend_vel(tendVel, err)

      ! vertical mixing treated implicitly in a later routine
      ! adjust total velocity tendency based on wetting/drying

#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) &
      !$acc    present(tendVel, wettingVelocityFactor)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1, nEdgesAll
      do k=1,nVertLevels
         tendVel(k,iEdge) = tendVel(k,iEdge) * (1.0_RKIND - wettingVelocityFactor(k,iEdge))
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
      !$acc    ssh, normalVelocity, &
      !$acc    layerThickness, activeTracers)

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
      !    layerThickEdgeFlux      - flux-related layer thickness on edge

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
         indxCFC11,          &! tracer index for CFC11
         indxCFC12,          &! tracer index for CFC12
         nTracersGroup,      &! number of tracers in each tracer group
         nTracersEcosys       ! number of ecosystem tracers

      logical :: &
         activeTracersOnly,  &! only compute tend for active tracers
         isActiveTracer,     &! is this group the active tracer group
         computeBudgets       ! compute active tracer budgets

      ! iterator and modified strings for tracer categories
      type (mpas_pool_iterator_type) :: groupItr
      character (len=:), allocatable :: &
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
      integer, pointer :: &
         indexCFC11, indexCFC12 ! tracer indicies for CFCs

      logical, pointer :: & ! option configuration flags
         config_use_tracerGroup, &
         config_use_tracerGroup_surface_bulk_forcing, &
         config_use_tracerGroup_surface_restoring,    &
         config_use_tracerGroup_interior_restoring,   &
         config_use_tracerGroup_exponential_decay,    &
         config_use_tracerGroup_idealAge_forcing,     &
         config_use_tracerGroup_ttd_forcing

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
         fractionAbsorbed,             &! frac sfc flux absorbed
         fractionAbsorbedRunoff,       &! frac runoff flux absorbed
         fractionAbsorbedSubglacialRunoff,       &! frac subglacial runoff flux absorbed
         tracerGroupSurfaceFluxRunoff, &! runoff flux
         tracerGroupSurfaceFluxSubglacialRunoff, &! subglacial runoff flux
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

      ! return if all tracer tendencies disabled
      ! otherwise start timer and proceed
   
      if (config_disable_tr_all_tend) return
      call mpas_timer_start("ocn_tend_tracer")

      ! set time level to default 1 or override with input arg
      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      ! set flag to compute only active tracers (default false)
      if (present(activeTracersOnlyIn)) then
         activeTracersOnly = activeTracersOnlyIn
      else
         activeTracersOnly = .false.
      end if

      ! retrieve data from pools 
      call mpas_pool_get_subpool(statePool,   'tracers', &
                                               tracersPool)
      call mpas_pool_get_subpool(tendPool,    'tracersTend', &
                                               tracersTendPool)
      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', &
                                               tracersSurfaceFluxPool)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSalinity)
      indxTemp = indexTemperature ! convert to int scalar due to
      indxSalt = indexSalinity    ! problems with int scalar pointers

      call mpas_pool_get_array(statePool,   'layerThickness', &
                                             layerThickness, timeLevel)
      call mpas_pool_get_array(forcingPool, 'penetrativeTemperatureFlux', &
                                             penetrativeTemperatureFlux)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', &
                                             fractionAbsorbed)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', &
                                             fractionAbsorbedRunoff)
      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedSubglacialRunoff', &
                                             fractionAbsorbedSubglacialRunoff)

      ! allocate and transfer data not specific to tracer groups
      allocate(normalThicknessFlux(nVertLevels, nEdgesAll+1))

      !$acc enter data create(normalThicknessFlux) &
      !$acc            copyin(layerThickness)
      !TEMPORARY - once diagnotics module ported, these should already
      !            be on the device
      !$acc update device (normalTransportVelocity, layerThickEdgeFlux, &
      !$acc                vertAleTransportTop)

      ! compute transport velocity to use for all tracers
#ifdef MPAS_OPENACC
      !$acc parallel loop collapse(2) &
      !$acc    present(normalThicknessFlux, normalTransportVelocity, &
      !$acc            layerThickEdgeFlux)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1,nEdgesAll
      do k = 1,nVertLevels
         normalThicknessFlux(k,iEdge) = &
                 normalTransportVelocity(k,iEdge)* &
                 layerThickEdgeFlux(k, iEdge)
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !
      ! Start computation of tracer tendencies.
      ! Loop over tracer groups and only compute if:
      !  - the pool member is a field
      !  - the group is the active tracer group
      !  - other groups if the active tracer only flag not on
      !  - the specific use flag for the given group is enabled 
      ! We do not indent these various loops/conditions to avoid
      ! excessive indentation.
      !

      ! loop over groups
      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )

      ! make sure group member is a field
      if ( groupItr % memberType == MPAS_POOL_FIELD ) then

      ! Set active tracer flag and compute if this is the active tracer
      ! group. If not active tracer, check active tracer only flag.

      groupName = trim(groupItr%memberName)
      if (groupName=='activeTracers') then
         isActiveTracer = .true.
      else
         isActiveTracer = .false.
      endif
      if (isActiveTracer .or. (.not. activeTracersOnly)) then

      ! Finally, check config_use flag for this tracer group
      configBase = 'config_use_' // groupName
      call mpas_pool_get_config(ocnConfigs, configBase, &
                                config_use_tracerGroup)

      if ( config_use_tracerGroup ) then

         ! This group meets all requirements, now check for 
         ! requested options

         computeBudgets =  (isActiveTracer .and.  &
                            config_compute_active_tracer_budgets)

         modifiedConfigName = configBase // '_surface_bulk_forcing'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_surface_bulk_forcing)
         modifiedConfigName = configBase // '_surface_restoring'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_surface_restoring)
         modifiedConfigName = configBase // '_interior_restoring'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_interior_restoring)
         modifiedConfigName = configBase // '_exponential_decay'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_exponential_decay)
         modifiedConfigName = configBase // '_idealAge_forcing'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_idealAge_forcing)
         modifiedConfigName = configBase // '_ttd_forcing'
         call mpas_pool_get_config(ocnConfigs, modifiedConfigName, &
                            config_use_tracerGroup_ttd_forcing)

         ! Get tracers from pool, determine number of tracers in group
         call mpas_pool_get_array(tracersPool, groupName, &
                                  tracerGroup, timeLevel)
         nTracersGroup = size(tracerGroup, dim=1)

         ! Get Tendency array
         modifiedGroupName = groupName // "Tend"
         call mpas_pool_get_array(tracersTendPool, &
                                  modifiedGroupName, &
                                  tracerGroupTend)

         ! Get surface flux array
         modifiedGroupName = groupName // "SurfaceFlux"
         call mpas_pool_get_array(tracersSurfaceFluxPool, &
                                  modifiedGroupName, &
                                  tracerGroupSurfaceFlux)

         ! Get Array of total surface temp/salt flux 
         ! (includes thickness tendencies)
         call mpas_pool_get_array(tracersSurfaceFluxPool, &
                                  'nonLocalSurfaceTracerFlux', &
                                   nonLocalSurfaceTracerFlux)

         ! Get surface flux due to river runoff array
         !maltrud only active tracers have runoff flux for now,
         ! but we still need to associate for ALL tracers
         modifiedGroupName = groupName // "SurfaceFluxRunoff"
         call mpas_pool_get_array(tracersSurfaceFluxPool, &
                                  modifiedGroupName, &
                                  tracerGroupSurfaceFluxRunoff)

         if (trim(config_subglacial_runoff_mode) == 'data') then
             ! Get surface flux due to subglacial runoff array
             ! only active tracers have subglacial runoff flux for now,
             ! but we still need to associate for ALL tracers
             modifiedGroupName = groupName // "SurfaceFluxSubglacialRunoff"
             call mpas_pool_get_array(tracersSurfaceFluxPool, &
                                      modifiedGroupName, &
                                      tracerGroupSurfaceFluxSubglacialRunoff)
         end if

         ! Get surface flux removed array to keep track of how much
         ! flux is ignored
         modifiedGroupName = groupName // "SurfaceFluxRemoved"
         call mpas_pool_get_array(tracersSurfaceFluxPool, &
                                  modifiedGroupName, &
                                  tracerGroupSurfaceFluxRemoved)

         !
         ! initialize tracer surface fluxes and tendencies to zero.
         !
         !$omp parallel
         !$omp do schedule(runtime) private(k,n)
         do iCell = 1, nCellsAll
         do k=1,nVertLevels
         do n=1,nTracersGroup
            tracerGroupTend(n,k,iCell) = 0.0_RKIND
         end do
         end do
         end do
         !$omp end do

         !$omp do schedule(runtime) private(n)
         do iCell = 1, nCellsAll
         do n=1,nTracersGroup
            tracerGroupSurfaceFlux       (n,iCell) = 0.0_RKIND
            tracerGroupSurfaceFluxRunoff (n,iCell) = 0.0_RKIND
            tracerGroupSurfaceFluxRemoved(n,iCell) = 0.0_RKIND
         end do
         end do
         !$omp end do
         !$omp end parallel

         if (trim(config_subglacial_runoff_mode) == 'data') then
         !$omp parallel
         !$omp do schedule(runtime) private(n)
         do iCell = 1, nCellsAll
             do n=1,nTracersGroup
                tracerGroupSurfaceFluxSubglacialRunoff (n,iCell) = 0.0_RKIND
             end do
         end do
         !$omp end do
         !$omp end parallel
         end if

         !
         ! compute surface tracer flux from bulk forcing
         !
         if (config_use_tracerGroup_surface_bulk_forcing) then

            call ocn_surface_bulk_forcing_tracers(meshPool, groupName, &
                                  forcingPool, tracerGroup, &
                                  tracerGroupSurfaceFlux, &
                                  tracerGroupSurfaceFluxRunoff, &
                                  tracerGroupSurfaceFluxRemoved, &
                                  dt, layerThickness, err)
            if (trim(config_subglacial_runoff_mode) == 'data') then
                call ocn_surface_bulk_forcing_tracers_subglacial_runoff(meshPool, groupName, &
                                      forcingPool, &
                                      tracerGroupSurfaceFluxSubglacialRunoff, err)
            end if

         end if ! sfc bulk forcing

         !
         ! compute ecosystem source-sink tendencies and net surface fluxes
         ! NOTE: must be called before ocn_tracer_surface_flux_tend
         !
         if ( groupName == 'ecosysTracers' ) then

            call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                     activeTracers, timeLevel)
            call ocn_tracer_ecosys_compute(activeTracers, tracerGroup, &
                                           forcingPool, nTracersGroup, &
                                           nCellsOwned, latCell, minLevelCell, &
                                           maxLevelCell, nVertLevels, &
                                           layerThickness, zMid, &
                                           indxTemp, indxSalt, &
                                           tracerGroupTend, err)

            call ocn_tracer_ecosys_surface_flux_compute(activeTracers, &
                                           tracerGroup, forcingPool,  &
                                           nTracersGroup, nCellsOwned, &
                                           zMid, minLevelCell, indxTemp, &
                                           indxSalt, tracerGroupSurfaceFlux, &
                                           err)
                                           
         endif ! ecosys tracers

         !
         ! compute DMS source-sink tendencies and net surface fluxes
         ! NOTE: must be called before ocn_tracer_surface_flux_tend
         !
         if ( groupName == 'DMSTracers' ) then
            call mpas_pool_get_array(tracersPool, 'ecosysTracers', &
                                     ecosysTracers, timeLevel)
            nTracersEcosys = size(ecosysTracers, dim=1)
            call ocn_tracer_DMS_compute(activeTracers, tracerGroup, &
                                        nTracersGroup, ecosysTracers, &
                                        nTracersEcosys, forcingPool, &
                                        nCellsOwned, minLevelCell, &
                                        maxLevelCell, nVertLevels, &
                                        layerThickness, indxTemp, &
                                        indxSalt, tracerGroupTend, &
                                        err)

            call ocn_tracer_DMS_surface_flux_compute(activeTracers, &
                                        tracerGroup, forcingPool,  &
                                        nTracersGroup, nCellsOwned, &
                                        zMid, minLevelCell, indxTemp, &
                                        indxSalt, tracerGroupSurfaceFlux,  &
                                        tracerGroupSurfaceFluxRemoved, err)
         endif ! DMS

         !
         ! compute MacroMolecules source-sink tendencies and net surface fluxes
         ! NOTE: must be called before ocn_tracer_surface_flux_tend
         !
         if ( groupName == 'MacroMoleculesTracers' ) then
            call mpas_pool_get_array(tracersPool, 'ecosysTracers', &
                                     ecosysTracers, timeLevel)
            nTracersEcosys = size(ecosysTracers, dim=1)
            call ocn_tracer_MacroMolecules_compute(tracerGroup, &
                                  nTracersGroup, ecosysTracers, &
                                  nTracersEcosys, forcingPool, &
                                  nCellsOwned, minLevelCell, &
                                  maxLevelCell, nVertLevels, &
                                  layerThickness,  &
                                  tracerGroupTend, err)

            call ocn_tracer_MacroMolecules_surface_flux_compute( &
                                  activeTracers, tracerGroup, &
                                  forcingPool, nTracersGroup, &
                                  nCellsOwned, zMid, &
                                  indxTemp, indxSalt, &
                                  tracerGroupSurfaceFlux, err)

         endif ! MacroMolecules

         !
         ! compute CFC surface fluxes, there are no interior sources/sinks
         !
         if ( groupName == 'CFCTracers' ) then

            call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                     activeTracers, timeLevel)

            indxCFC11 = 0
            indxCFC12 = 0
            if (config_use_CFC11) then
              call mpas_pool_get_dimension(tracersPool, 'index_CFC11', indexCFC11)
              indxCFC11 = indexCFC11 ! convert to int scalar due to problems with int scalar pointers
            endif
            if (config_use_CFC12) then
              call mpas_pool_get_dimension(tracersPool, 'index_CFC12', indexCFC12)
              indxCFC12 = indexCFC12 ! convert to int scalar due to problems with int scalar pointers
            endif
                                                 
            call ocn_tracer_CFC_surface_flux_compute(activeTracers, &
                                        tracerGroup, forcingPool,  &
                                        nTracersGroup, nCellsOwned, &
                                        indxTemp, indxSalt, tracerGroupSurfaceFlux, indxCFC11, indxCFC12, err)
         endif ! CFC

         !
         ! ocean surface restoring
         !
         if (config_use_tracerGroup_surface_restoring) then
            call mpas_timer_start("surface_restoring_"//groupName)
            call mpas_pool_get_subpool(forcingPool, &
                                   'tracersSurfaceRestoringFields', &
                                    tracersSurfaceRestoringFieldsPool)

            modifiedGroupName = groupName // "PistonVelocity"
            call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool,&
                                     modifiedGroupName, &
                                     tracerGroupPistonVelocity)

            modifiedGroupName = groupName//"SurfaceRestoringValue"
            call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool,&
                                     modifiedGroupName, &
                                     tracerGroupSurfaceRestoringValue)

            ! Note: monthly surface salinity restoring is a special
            ! case for tracer restoring
            call ocn_tracer_surface_restoring_compute(groupName, &
                   nTracersGroup, nCellsAll, tracerGroup,  &
                   tracerGroupPistonVelocity, &
                   tracerGroupSurfaceRestoringValue, &
                   tracerGroupSurfaceFlux, indxSalt, &
                   config_use_surface_salinity_monthly_restoring, &
                   config_salinity_restoring_constant_piston_velocity, &
                   err)
            call mpas_timer_stop("surface_restoring_"//groupName)

         endif ! sfc restoring

         ! tracer fluxes at the land-ice / ocean interface
         ! this is a flux at the top ocean surface -- so these
         ! fluxes are added into tracerGroupSurfaceFlux
         call ocn_surface_land_ice_fluxes_tracers(meshPool, groupName, &
                             forcingPool, tracerGroupSurfaceFlux, err)

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
            call mpas_timer_start("interior_restoring_"//groupName,&
                                  .false.)
            call mpas_pool_get_subpool(forcingPool, &
                              'tracersInteriorRestoringFields', &
                               tracersInteriorRestoringFieldsPool)

            modifiedGroupName = groupName//"InteriorRestoringRate"
            call mpas_pool_get_array(tracersInteriorRestoringFieldsPool,&
                                     modifiedGroupName, &
                                     tracerGroupInteriorRestoringRate)

            modifiedGroupName = groupName//"InteriorRestoringValue"
            call mpas_pool_get_array(tracersInteriorRestoringFieldsPool,&
                                     modifiedGroupName, &
                                     tracerGroupInteriorRestoringValue)

            call ocn_tracer_interior_restoring_compute(nTracersGroup, &
                        nCellsAll, minLevelCell, maxLevelCell, &
                        layerThickness, tracerGroup, &
                        tracerGroupInteriorRestoringRate, &
                        tracerGroupInteriorRestoringValue, &
                        tracerGroupTend, err)

            call mpas_timer_stop("interior_restoring_"//groupName)

         endif ! interior restoring

         !
         ! exponential decay tendency
         !
         if (config_use_tracerGroup_exponential_decay) then
            call mpas_log_write( &
                        "WARNING: exponential decay not fully tested", &
                        MPAS_LOG_WARN)
            call mpas_timer_start("exponential decay "//groupName)
            call mpas_pool_get_subpool(forcingPool, &
                               'tracersExponentialDecayFields', &
                                tracersExponentialDecayFieldsPool)

            modifiedGroupName = groupName // "ExponentialDecayRate"
            call mpas_pool_get_array(tracersExponentialDecayFieldsPool,&
                                     modifiedGroupName, &
                                     tracerGroupExponentialDecayRate)

            call ocn_tracer_exponential_decay_compute(nTracersGroup, &
                        nCellsAll, minLevelCell, maxLevelCell, &
                        layerThickness, tracerGroup, &
                        tracerGroupExponentialDecayRate, &
                        tracerGroupTend, err)

            call mpas_timer_stop("exponential decay "//groupName)

         endif ! exponential decay

         !
         ! ideal age forcing tendency
         !   note: ocn_tracer_ideal_age_compute resets tracers in top
         !   layer to zero
         !
         if (config_use_tracerGroup_idealAge_forcing) then
            call mpas_timer_start("ideal age " // groupName)
            if (trim(groupName) /= 'idealAgeTracers') call mpas_log_write( &
                                                      "WARNING: ideal age not fully tested", &
                                                      MPAS_LOG_WARN)

            call mpas_pool_get_subpool(forcingPool, &
                                       'tracersIdealAgeFields', &
                                        tracersIdealAgeFieldsPool)

            modifiedGroupName = groupName // "IdealAgeMask"
            call mpas_pool_get_array(tracersIdealAgeFieldsPool, &
                                     modifiedGroupName, &
                                     tracerGroupIdealAgeMask)
            call ocn_tracer_ideal_age_compute(nTracersGroup, nCellsAll,&
                                      minLevelCell, maxLevelCell, &
                                      layerThickness, &
                                      tracerGroupIdealAgeMask, &
                                      tracerGroup, tracerGroupTend, err)

            call mpas_timer_stop("ideal age " // groupName)

         endif ! ideal age

         !
         ! transit-time distribution (TTD) forcing tendency
         !   note: no tendency is actually computed in ocn_tracer_TTD_compute
         !   note: rather, tracerGroup is reset to tracerGroupTTDMask in
         !   top-most layer
         !
         if (config_use_tracerGroup_ttd_forcing) then
            call mpas_timer_start("TTD " // groupName)
            call mpas_log_write("WARNING: TTD not fully tested", &
                                MPAS_LOG_WARN)

            call mpas_pool_get_subpool(forcingPool, &
                                       'tracersTTDFields', &
                                        tracersTTDFieldsPool)

            modifiedGroupName = groupName // "TTDMask"
            call mpas_pool_get_array(tracersTTDFieldsPool, &
                                     modifiedGroupName, &
                                     tracerGroupTTDMask)
            call ocn_tracer_TTD_compute(nTracersGroup, nCellsAll, &
                                        maxLevelCell, layerThickness, &
                                        tracerGroupTTDMask, &
                                        tracerGroup, err)

            call mpas_timer_stop("TTD " // groupName)
         endif ! TTD

         !$acc enter data &
         !$acc    copyin(tracerGroup, tracerGroupTend)

         !
         ! Accumulate tendencies for all tracers
         ! First is the horizontal advection term
         !     -div( layerThickness \phi u)
         ! Tracer budget is computed and stored within the tracer adv
         ! routine
         call ocn_tracer_advection_tend(tracerGroup, &
                            normalThicknessFlux, vertAleTransportTop, &
                            layerThickness, dt, tracerGroupTend, &
                            computeBudgets)

         !$acc update host (activeTracerHorizontalAdvectionEdgeFlux, &
         !$acc              activeTracerHorizontalAdvectionTendency, &
         !$acc              activeTracerVerticalAdvectionTopFlux, &
         !$acc              activeTracerVerticalAdvectionTendency)

         !$acc exit data &
         !$acc    copyout(tracerGroupTend) &
         !$acc    delete(tracerGroup)

         !
         ! Add horizontal tracer diffusion,
         !    eg del2 form: div(h \kappa_2 \nabla \phi)
         !

         if (computeBudgets) then
            ! compute budget contribution by taking the difference
            ! of the total tendency before and after the call
            ! here we store the before value for later subtraction
            !$omp parallel
            !$omp do schedule(runtime) private(k,n)
            do iCell = 1, nCellsAll
            do k=1,nVertLevels
            do n=1,nTracersGroup
               activeTracerHorMixTendency(n,k,iCell) = &
                          tracerGroupTend(n,k,iCell)
            enddo
            enddo
            enddo
            !$omp end do
            !$omp end parallel
         endif ! compute budgets

         call ocn_tracer_hmix_tend(meshPool, layerThickEdgeMean, &
                                   layerThickness, zMid, tracerGroup, &
                                   RediKappa, slopeTriadUp, &
                                   slopeTriadDown, limiterUp, &
                                   limiterDown, dt, isActiveTracer, &
                                   rediLimiterCount, &
                                   tracerGroupTend, err)

         if (computeBudgets) then
            ! Finish budget calculation by subtracting before, after
            !$omp parallel
            !$omp do schedule(runtime) private(k,n)
            do iCell = 1, nCellsAll
            do k=1,nVertLevels
            do n=1,nTracersGroup
               activeTracerHorMixTendency(n,k,iCell) = &
                          tracerGroupTend(n,k,iCell) - &! tendency after
               activeTracerHorMixTendency(n,k,iCell)    ! tendency before
            enddo
            enddo
            enddo
            !$omp end do
            !$omp end parallel
         endif ! compute budgets

         !
         ! convert the surface tracer flux into a tracer tendency by
         ! distributing the flux across some number of surface layers
         !
         if (computeBudgets) then
            ! as in hmix above, the budget contribution is difference
            ! between total tend before/after call - store before value
            !$omp parallel
            !$omp do schedule(runtime) private(k,n)
            do iCell = 1, nCellsAll
            do k = 1, nVertLevels
            do n = 1, nTracersGroup
               activeTracerSurfaceFluxTendency(n,k,iCell) = &
                               tracerGroupTend(n,k,iCell)
            end do
            end do
            end do
            !$omp end do
            !$omp end parallel
         endif ! compute budgets

         call ocn_tracer_surface_flux_tend(meshPool, fractionAbsorbed, &
                               fractionAbsorbedRunoff, &
                               layerThickness, tracerGroupSurfaceFlux, &
                               tracerGroupSurfaceFluxRunoff,  &
                               tracerGroupTend, err)
         if (trim(config_subglacial_runoff_mode) == 'data') then
                call ocn_tracer_surface_flux_tend_subglacial_runoff(meshPool, fractionAbsorbedSubglacialRunoff, &
                               layerThickness, tracerGroupSurfaceFluxSubglacialRunoff,  &
                               tracerGroupTend, err)
         end if

         !
         ! Compute shortwave absorption
         !

         if (computeBudgets) then
            ! finish previous budget difference and store
            ! the before value for short wave
            !$omp parallel
            !$omp do schedule(runtime) private(k,n)
            do iCell = 1, nCellsAll
            do k = 1, nVertLevels
               do n = 1, nTracersGroup
                  activeTracerSurfaceFluxTendency(n,k,iCell) = &
                                  tracerGroupTend(n,k,iCell) - &
                  activeTracerSurfaceFluxTendency(n,k,iCell)
               end do

               temperatureShortWaveTendency(k,iCell) = &
                     tracerGroupTend(indxTemp,k,iCell)
            end do
            end do
            !$omp end do
            !$omp end parallel
         endif ! compute budgets

         ! only need sw flux for temperature
         if (isActiveTracer) &
            call ocn_tracer_short_wave_absorption_tend(meshPool, &
                       swForcingPool, forcingPool, indxTemp, &
                       layerThickness, penetrativeTemperatureFlux, &
                       penetrativeTemperatureFluxOBL, tracerGroupTend, &
                       bottomLayerShortwaveTemperatureFlux, err)

         if (computeBudgets) then
            !$omp parallel
            !$omp do schedule(runtime) private(k)
            do iCell = 1, nCellsAll
            do k=1,nVertLevels
               temperatureShortWaveTendency(k,iCell) = &
                   tracerGroupTend(indxTemp,k,iCell) - &
               temperatureShortWaveTendency(k,iCell)
            end do
            end do
            !$omp end do
            !$omp end parallel
         endif ! compute budgets

         !
         ! Compute tracer tendency due to non-local flux computed in KPP
         !
         if (config_use_cvmix_kpp .or. config_use_gotm) then
            call ocn_compute_KPP_input_fields(statePool, forcingPool,&
                                                 meshPool, timeLevel)

            if (.not. config_cvmix_kpp_nonlocal_with_implicit_mix) then
               call mpas_timer_start("non-local flux from KPP")
               if (computeBudgets) then
                  !$omp parallel
                  !$omp do schedule(runtime) private(k,n)
                  do iCell = 1, nCellsAll
                  do k = 1, nVertLevels
                  do n = 1, nTracersGroup
                     activeTracerNonLocalTendency(n,k,iCell) = &
                                  tracerGroupTend(n,k,iCell)
                  end do
                  end do
                  end do
                  !$omp end do
                  !$omp end parallel
               endif ! compute budgets

               if (isActiveTracer) then
                  call ocn_tracer_nonlocalflux_tend(meshPool, &
                                          vertNonLocalFlux, &
                                          nonLocalSurfaceTracerFlux, &
                                          tracerGroupTend, err)
               else ! not active
                  call ocn_tracer_nonlocalflux_tend(meshPool, &
                                          vertNonLocalFlux, &
                                          tracerGroupSurfaceFlux, &
                                          tracerGroupTend, err)
               endif ! active tracer

               if (computeBudgets) then
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
               endif ! compute budgets
               call mpas_timer_stop("non-local flux from KPP")
            end if ! not non-local with implicit mix
         end if ! KPP

         ! Compute tracer tendency due to production/destruction of
         ! frazil ice
         call ocn_frazil_forcing_tracers(meshPool, tracersPool, &
                     groupName, forcingPool, tracerGroupTend, err)


      end if ! config use flag on
      end if ! active only
      end if ! member is field
      end do ! loop over tracer groups

      !$acc exit data delete (normalThicknessFlux, layerThickness)

      deallocate(normalThicknessFlux)

      call mpas_timer_stop("ocn_tend_tracer")

   !--------------------------------------------------------------------

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

   subroutine ocn_tend_freq_filtered_thickness(tendPool, statePool, &
                                               timeLevelIn)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool           !< [in] ocean state information

      ! actual input fields used from statePool
      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         normalVelocity,         &! normal velocity
         layerThickness,         &! layer thickness
         lowFreqDivergence,      &! low freq divergence of hu
         highFreqThickness        ! high freq thickness

      integer, intent(in), optional :: &
         timeLevelIn         !< [in] Time level for state fields

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         tendPool            !< [out] Tendency terms

      ! actual output fields from tendPool
      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         tend_lowFreqDivergence, &! tendency for low freq divergence
         tend_highFreqThickness   ! tendency for high freq thickness

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer ::     &
         i,k,iCell,  &! loop iterators for edges, vertical, cells
         iEdge,      &! index for edge of cell
         err,        &! local error flag
         timeLevel    ! time level for ocean state

      real (kind=RKIND) :: &
         flux,             &! local flux across edge
         invAreaCell,      &! inverse cell area
         div_hu_btr,       &! barotropic avg of div_hu
         totalThickness     ! total depth of column

      real (kind=RKIND) :: &
         thickness_filter_timescale_sec, &
         highFreqThick_restore_time_sec

      real (kind=RKIND), dimension(:), allocatable :: &
         div_hu             ! local temp for divergence of hu

      ! End preamble
      !-------------
      ! Begin code

      call mpas_timer_start("ocn_tend_freq_filtered_thickness")
      err = 0

      ! Set time level for state variables (default is 1)
      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      ! Retrieve variables from pools
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergence, timeLevel)
      call mpas_pool_get_array(statePool, 'highFreqThickness', &
                                           highFreqThickness, timeLevel)

      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', &
                                     tend_lowFreqDivergence)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                     tend_highFreqThickness)

      !
      ! Low Frequency Divergence and high frequency thickness Tendency
      !

      ! Convert restore time from days to seconds
      thickness_filter_timescale_sec = config_thickness_filter_timescale*86400.0_RKIND
      highFreqThick_restore_time_sec = config_highFreqThick_restore_time*86400.0_RKIND

      allocate(div_hu(nVertLevels))

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp    private(i, k, iCell, iEdge, div_hu_btr, &
      !$omp            invAreaCell, flux, totalThickness) &
      !$omp    firstprivate(div_hu)
      do iCell = 1, nCellsAll
         do k=1,nVertLevels
            tend_lowFreqDivergence(k, iCell) = 0.0_RKIND
            tend_highFreqThickness(k, iCell) = 0.0_RKIND
            div_hu(k) = 0.0_RKIND
         end do ! vert (k) loop
         div_hu_btr = 0.0_RKIND
         invAreaCell = 1.0_RKIND / areaCell(iCell)

         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)

            do k = 1, maxLevelEdgeTop(iEdge)
               flux = layerThickEdgeFlux(k,iEdge)*normalVelocity(k,iEdge)* &
                      dvEdge(iEdge)*edgeSignOnCell(i,iCell)*invAreaCell
               div_hu(k) = div_hu(k) - flux
               div_hu_btr = div_hu_btr - flux
            end do
         end do ! edges on cell

         totalThickness = sum(layerThickness(1:maxLevelCell(iCell),iCell))
         do k = 1, maxLevelCell(iCell)

            tend_lowFreqDivergence(k,iCell) = &
               -2.0*pi/thickness_filter_timescale_sec &
              *(lowFreqDivergence(k,iCell)  - div_hu(k) &
                + div_hu_btr * layerThickness(k,iCell)/totalThickness)

            tend_highFreqThickness(k,iCell) = - div_hu(k) + &
                  div_hu_btr*layerThickness(k,iCell)/totalThickness + &
                  lowFreqDivergence(k,iCell) + &
                  use_highFreqThick_restore* &
                  (-2.0*pi/highFreqThick_restore_time_sec* &
                   highFreqThickness(k,iCell) )

         end do ! vert (k) loop
      end do ! cell loop
      !$omp end do
      !$omp end parallel

      deallocate(div_hu)

      !
      !  add horizontal hhf diffusion, div(\kappa_{hf} \nabla h^{hf})
      !
      call ocn_high_freq_thickness_hmix_del2_tend(highFreqThickness, &
                                        tend_highFreqThickness, err)

      call mpas_timer_stop("ocn_tend_freq_filtered_thickness")

   !--------------------------------------------------------------------

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
