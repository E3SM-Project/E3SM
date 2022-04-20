! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tendency_benthos
!
!> \brief MPAS ocean tendency driver
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011
!> \details
!>  This module contains the routines for computing
!>  tendency terms for the ocean primitive equations.
!
!-----------------------------------------------------------------------

module ocn_tendency_benthos

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
   use ocn_tracer_CFC
   use ocn_tracer_benthos

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

   public :: ocn_tend_benthos

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
!
!  routine ocn_tend_benthos
!
!> \brief   Computes tracer tendency
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine computes tracer tendencies for the ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_tend_benthos(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, &
                              diagnosticsPool, dt, activeTracersOnlyIn, timeLevelIn )!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool,          &!< [in] ocean state variables
         meshPool             !< [in] mesh information

      real (kind=RKIND), dimension(:), pointer :: bottomDepth
      real (kind=RKIND), dimension(:), pointer :: benthosInterfaceLayerThickness

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), intent(in) :: &
         dt                   !< [in] time step (seconds)

      logical, intent(in), optional :: &
         activeTracersOnlyIn  !< [in] only compute for active tracers

      integer, intent(in), optional :: &
         timeLevelIn          !< [in] time index to use for state vars

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         forcingPool,   &!< [inout] forcing data
         swForcingPool, &!< [inout] shortwave forcing data
         tendPool,      &!< [inout] Tendency terms for all variables
         diagnosticsPool !< [inout] tracer diagnostics

      integer :: &
         iCell, iEdge, k, n, &! loop counters for cell,edge,vert,tracer
         iCellBenthic,       &! loop counter for benthic cells
         err,                &! internal error flag
         timeLevel,          &! time level to use for state variables
         indxTemp, indxSalt, &! tracer index for temperature, salinity
         indxCFC11,          &! tracer index for CFC11
         indxCFC12,          &! tracer index for CFC12
         nTracersGroup,      &! number of tracers in each tracer group
         nTracersEcosys,     &! number of ecosystem tracers
         nBenthicTracers      ! number of benthic tracers

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

      ! benthos biogeochemistry
      type (mpas_pool_type), pointer :: &
         seafloorTracersPool,  &! benthos tracers
         seafloorTendPool,     &! benthos tracer tendencies
         seafloorStoragePool,  &
         benthosAuxiliaryPool   ! holds index to benthic mesh

      ! pointers for retrieving data from pools
      integer, pointer :: &
         indexTemperature, indexSalinity, &! tracer index for temp,salt
         nBenthicVertLevels ! number of benthic vertical grid levels

      integer, pointer :: &
         indexCFC11, indexCFC12 ! tracer indicies for CFCs

      logical, pointer :: & ! option configuration flags
         config_use_tracerGroup, &
         config_use_tracerGroup_surface_bulk_forcing, &
         config_use_tracerGroup_surface_restoring,    &
         config_use_tracerGroup_interior_restoring,   &
         config_use_tracerGroup_exponential_decay,    &
         config_use_tracerGroup_idealAge_forcing,     &
         config_use_tracerGroup_ttd_forcing,          &
         config_use_benthosTracers,                   &
         config_use_ecosysTracers

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
         tracerGroupInteriorRestoringValue,&! int restoring value
         benthosTracers,  &
         benthosTracersTend

      real (kind=RKIND), dimension(:,:), pointer, contiguous :: &
         benthosPhPrev3d, &! retain ph in sediment
         benthosTracersStorage! retain losses of benthos tracers

      ! Scratch Arrays
      real (kind=RKIND), dimension(:,:), allocatable :: &
         normalThicknessFlux ! Flux of thickness through edge (m^2/s)

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! return if all tracer tendencies disabled
      ! otherwise start timer and proceed

      !if (config_disable_tr_all_tend) return
      call mpas_timer_start("ocn_tend_benthos")

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
      call mpas_pool_get_subpool(statePool,   'seafloorTracers', &
                                               seafloorTracersPool)
      call mpas_pool_get_subpool(tendPool,    'seafloorTend', &
                                               seafloorTendPool)
      call mpas_pool_get_subpool(statePool,   'seafloorStorage', &
                                               seafloorStoragePool)
      call mpas_pool_get_dimension(meshPool,  'nBenthicVertLevels', &
                                               nBenthicVertLevels)
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

      call mpas_pool_get_config(ocnConfigs, 'config_use_benthosTracers', config_use_benthosTracers)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ecosysTracers', config_use_ecosysTracers)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracers, timeLevel)
      call mpas_pool_get_array(seafloorTracersPool, 'benthosPhPrev3d', benthosPhPrev3d, timeLevel)
      call mpas_pool_get_array(seafloorTendPool, 'benthosTracersTend', benthosTracersTend)
      call mpas_pool_get_array(seafloorStoragePool, 'benthosTracersStorage', benthosTracersStorage)
      call mpas_pool_get_array(seafloorTendPool, 'benthosInterfaceLayerThickness', benthosInterfaceLayerThickness)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, timeLevel)

      nTracersEcosys = size(ecosysTracers, dim=1)
      nBenthicTracers = size(benthosTracers, dim=1)

      call mpas_timer_start("ocn_tend_tracer")
     !
     ! initialize benthos tracer tendencies to zero.
     !

     !$omp parallel
     !$omp do schedule(runtime) private(k,n)
     do iCellBenthic = 1, nBenthicCells
        do k=1,nBenthicVertLevels+1
        do n=1,nBenthicTracers
           benthosTracersTend(n,k,iCellBenthic) = 0.0_RKIND
        end do
        end do
      end do

      !$omp end do
      !$omp end parallel

      call benthos_tracer_compute(benthosTracers, benthosPhPrev3d, benthosTracersStorage, &
               activeTracers, nBenthicTracers, benthosInterfaceLayerThickness, &
               nTracersEcosys, statePool, tendPool, diagnosticsPool, forcingPool, &
               nCellsOwned, nBenthicVertLevels, bottomDepth, density, maxLevelCell, &
               nVertLevels, layerThickness, indexTemperature, indexSalinity, benthosTracersTend,&
               dt, err, config_use_ecosysTracers)

      call mpas_timer_stop("ocn_tend_tracer")

   !--------------------------------------------------------------------

   end subroutine ocn_tend_benthos!}}}

!***********************************************************************

end module ocn_tendency_benthos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
