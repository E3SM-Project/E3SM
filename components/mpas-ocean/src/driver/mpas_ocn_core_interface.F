! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
module ocn_core_interface

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_constants
   use mpas_log
   use ocn_core
   use mpas_attlist

   use ocn_forward_mode
   use ocn_analysis_mode
#ifndef EXCLUDE_INIT_MODE
   use ocn_init_mode
#endif

   private

   public :: ocn_setup_core, ocn_setup_domain

   contains
   !***********************************************************************
   !
   !  routine ocn_setup_core
   !
   !> \brief   Ocean core setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This routine is intended to setup the necessary variables within a core_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_setup_core(core)!{{{
      type (core_type), pointer :: core

      core % core_init => ocn_core_init
      core % core_run => ocn_core_run
      core % core_finalize => ocn_core_finalize
      core % define_packages => ocn_define_packages
      core % setup_packages => ocn_setup_packages
      core % setup_decompositions => ocn_setup_decompositions
      core % setup_clock => ocn_setup_clock
      core % setup_log => ocn_setup_log
      core % get_mesh_stream => ocn_get_mesh_stream
      core % setup_immutable_streams => ocn_setup_immutable_streams
      core % setup_derived_dimensions => ocn_setup_derived_dimensions
      core % setup_decomposed_dimensions => ocn_setup_decomposed_dimensions
      core % setup_block => ocn_setup_block
      core % setup_namelist => ocn_setup_namelists

      core % Conventions = 'MPAS'
      core % source = 'MPAS'
#include "../inc/core_variables.inc"

   end subroutine ocn_setup_core!}}}


   !***********************************************************************
   !
   !  routine ocn_setup_domain
   !
   !> \brief   Ocean domain setup routine
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This routine is intended to setup the necessary variables within a domain_type
   !>  for the ocean core.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_setup_domain(domain)!{{{
      type (domain_type), pointer :: domain

#include "../inc/domain_variables.inc"

   end subroutine ocn_setup_domain!}}}


   !***********************************************************************
   !
   !  function ocn_setup_packages
   !
   !> \brief   Pacakge setup routine
   !> \author  Doug Jacobsen
   !> \date    03/12/2015
   !> \details
   !>  This function is intended to correctly configure the packages for this MPAS
   !>   core. It can use any Fortran logic to properly configure packages, and it
   !>   can also make use of any namelist options. All variables in the model are
   !>   *not* allocated until after this routine is called.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_packages(configPool, packagePool, iocontext) result(ierr)!{{{

      use ocn_analysis_driver

      type (mpas_pool_type), intent(inout) :: configPool
      type (mpas_pool_type), intent(inout) :: packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext

      integer :: ierr

      integer :: err_tmp

      logical, pointer :: forwardModeActive
      logical, pointer :: analysisModeActive
      logical, pointer :: initModeActive
      logical, pointer :: thicknessFilterActive
      logical, pointer :: splitTimeIntegratorActive
      logical, pointer :: semiImplicitTimePKGActive
      logical, pointer :: splitAB2TimeIntegratorActive
      logical, pointer :: windStressBulkPKGActive
      logical, pointer :: variableBottomDragPKGActive
      logical, pointer :: tracerBudgetActive
      logical, pointer :: landIcePressurePKGActive
      logical, pointer :: dataLandIceFluxesPKGActive
      logical, pointer :: landIceFluxesPKGActive
      logical, pointer :: landIceCouplingPKGActive
      logical, pointer :: dataSubglacialRunoffFluxPKGActive
      logical, pointer :: thicknessBulkPKGActive
      logical, pointer :: frazilIceActive
      logical, pointer :: tidalForcingActive
      logical, pointer :: tidalPotentialForcingPKGActive
      logical, pointer :: vegetationDragPKGActive
      logical, pointer :: topographicWaveDragPKGActive
      logical, pointer :: inSituEOSActive
      logical, pointer :: variableShortwaveActive
      logical, pointer :: gmActive
      logical, pointer :: submesoActive
      logical, pointer :: timeVaryingAtmosphericForcingPKGActive
      logical, pointer :: timeVaryingLandIceForcingPKGActive
      logical, pointer :: gotmPKGActive
      logical, pointer :: verticalRemapPKGActive
      logical, pointer :: activeWavePKGActive
      logical, pointer :: subgridWetDryPKGActive
      logical, pointer :: scaledDISMFPKGActive

      type (mpas_pool_iterator_type) :: pkgItr
      logical, pointer :: packageActive

      logical, pointer :: tracerGroupPKGActive
      logical, pointer :: tracerGroupBulkRestoringPKGActive
      logical, pointer :: tracerGroupSurfaceRestoringPKGActive
      logical, pointer :: tracerGroupInteriorRestoringPKGActive
      logical, pointer :: tracerGroupExponentialDecayPKGActive
      logical, pointer :: tracerGroupIdealAgePKGActive
      logical, pointer :: tracerGroupTTDPKGActive

      logical, pointer :: config_use_tracerGroup
      logical, pointer :: config_use_tracerGroup_surface_bulk_forcing
      logical, pointer :: config_use_tracerGroup_surface_restoring
      logical, pointer :: config_use_tracerGroup_interior_restoring
      logical, pointer :: config_use_tracerGroup_exponential_decay
      logical, pointer :: config_use_tracerGroup_idealAge_forcing
      logical, pointer :: config_use_tracerGroup_ttd_forcing

      logical, pointer :: config_use_freq_filtered_thickness
      logical, pointer :: config_use_frazil_ice_formation
      logical, pointer :: config_use_tidal_forcing
      logical, pointer :: config_use_self_attraction_loading
      logical, pointer :: config_use_tidal_potential_forcing
      logical, pointer :: config_use_topographic_wave_drag
      logical, pointer :: config_use_GM
      logical, pointer :: config_submesoscale_enable
      logical, pointer :: config_use_Redi
      logical, pointer :: config_use_vegetation_drag
      logical, pointer :: config_use_time_varying_atmospheric_forcing
      logical, pointer :: config_use_time_varying_land_ice_forcing
      logical, pointer :: config_use_gotm
      logical, pointer :: config_use_active_wave
      logical, pointer :: config_use_subgrid_wetting_drying
      logical, pointer :: config_scale_dismf_by_removed_ice_runoff

      character (len=StrKIND), pointer :: config_time_integrator
      character (len=StrKIND), pointer :: config_ocean_run_mode
      character (len=StrKIND), pointer :: config_pressure_gradient_type
      character (len=StrKIND), pointer :: config_sw_absorption_type
      character (len=StrKIND), pointer :: config_vert_advection_method

      logical, pointer :: config_use_variable_drag
      logical, pointer :: config_use_bulk_wind_stress
      logical, pointer :: config_use_bulk_thickness_flux
      logical, pointer :: config_compute_active_tracer_budgets
      character (len=StrKIND), pointer :: config_land_ice_flux_mode
      character (len=StrKIND), pointer :: config_subglacial_runoff_mode

       type (mpas_pool_iterator_type) :: groupItr
       character (len=StrKIND) :: tracerGroupName, configName, packageName
       integer :: startIndex, strLen

      ierr = 0

      !
      ! determine the mode being used
      !
      call mpas_pool_get_package(packagePool, 'forwardModeActive', forwardModeActive)
      call mpas_pool_get_package(packagePool, 'analysisModeActive', analysisModeActive)
      call mpas_pool_get_package(packagePool, 'initModeActive', initModeActive)
      call mpas_pool_get_config(configPool, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         forwardModeActive = .true.
      endif
      if ( trim(config_ocean_run_mode) == 'analysis') then
         analysisModeActive = .true.
      endif
      if ( trim(config_ocean_run_mode) == 'init') then
         initModeActive = .true.
      endif

      !
      ! test for integration scheme
      !
      call mpas_pool_get_package(packagePool, 'splitTimeIntegratorActive', splitTimeIntegratorActive)
      call mpas_pool_get_config(configPool, 'config_time_integrator', config_time_integrator)
      if ( forwardModeActive ) then
         if (    config_time_integrator == trim('split_explicit') &
            .or. config_time_integrator == trim('unsplit_explicit') &
            .or. config_time_integrator == trim('split_implicit') &
            .or. config_time_integrator == trim('split_explicit_ab2') ) then
            splitTimeIntegratorActive = .true.
         end if
      endif
      call mpas_pool_get_package(packagePool, 'semiImplicitTimePKGActive', semiImplicitTimePKGActive)
      if ( forwardModeActive ) then
         if (config_time_integrator == trim('split_implicit') ) then
            semiImplicitTimePKGActive = .true.
         end if
      endif

      call mpas_pool_get_package(packagePool, 'splitAB2TimeIntegratorActive', splitAB2TimeIntegratorActive)
      if ( forwardModeActive ) then
         if (config_time_integrator == trim('split_explicit_ab2') ) then
            splitAB2TimeIntegratorActive = .true.
         end if
      endif

      !
      ! test for time filtering scheme
      !
      call mpas_pool_get_package(packagePool, 'thicknessFilterActive', thicknessFilterActive)
      call mpas_pool_get_config(configPool, 'config_use_freq_filtered_thickness', config_use_freq_filtered_thickness)
      if ( forwardModeActive ) then
         if (config_use_freq_filtered_thickness) then
            thicknessFilterActive = .true.
         end if
      endif

      !
      ! test for bulk forcing of layer thickness, thicknessBulkPKG
      !
      call mpas_pool_get_package(packagePool, 'thicknessBulkPKGActive', thicknessBulkPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_bulk_thickness_flux', config_use_bulk_thickness_flux)
      if ( config_use_bulk_thickness_flux ) then
         thicknessBulkPKGActive = .true.
      end if

      !
      ! test for bulk forcing of momentum by wind stress, windStressBulkPKG
      !
      call mpas_pool_get_package(packagePool, 'windStressBulkPKGActive', windStressBulkPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
      if ( config_use_bulk_wind_stress ) then
         windStressBulkPKGActive = .true.
      end if

      !
      ! test for variable bottom drag of momentum, variableBottomDragPKG
      call mpas_pool_get_package(packagePool, 'variableBottomDragPKGActive', variableBottomDragPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_variable_drag', config_use_variable_drag)
      if ( config_use_variable_drag ) then
        variableBottomDragPKGActive = .true.
      end if

      !
      ! test for tracer budget
      !
      call mpas_pool_get_package(packagePool, 'tracerBudgetActive', tracerBudgetActive)
      call mpas_pool_get_config(configPool, 'config_compute_active_tracer_budgets', config_compute_active_tracer_budgets)
      if ( config_compute_active_tracer_budgets ) then
         tracerBudgetActive = .true.
      end if

      !
      ! Test if chlorophyll, solar zenith angle, and clear sky radiation should be used
      !

      call mpas_pool_get_package(packagePool,'variableShortwaveActive',variableShortwaveActive)
      call mpas_pool_get_config(configPool,'config_sw_absorption_type',config_sw_absorption_type)
      if (trim (config_sw_absorption_type) == 'ohlmann00') then
         variableShortwaveActive = .true.
      end if

      !
      ! test for land ice pressure, landIcePressurePKG
      ! test for land ice fluxes, landIceFluxesPKG
      ! test for land ice coupling, landIceCouplingPKG
      !
      call mpas_pool_get_package(packagePool, 'landIcePressurePKGActive', landIcePressurePKGActive)
      call mpas_pool_get_package(packagePool, 'dataLandIceFluxesPKGActive', &
                                 dataLandIceFluxesPKGActive)
      call mpas_pool_get_package(packagePool, 'landIceFluxesPKGActive', landIceFluxesPKGActive)
      call mpas_pool_get_package(packagePool, 'landIceCouplingPKGActive', landIceCouplingPKGActive)
      call mpas_pool_get_config(configPool, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      if ( trim(config_land_ice_flux_mode) == 'pressure_only' ) then
         landIcePressurePKGActive = .true.
      else if ( trim(config_land_ice_flux_mode) == 'data' ) then
         landIcePressurePKGActive = .true.
         dataLandIceFluxesPKGActive = .true.
      else if ( trim(config_land_ice_flux_mode) == 'standalone' ) then
         landIcePressurePKGActive = .true.
         landIceFluxesPKGActive = .true.
      else if ( trim(config_land_ice_flux_mode) == 'coupled' ) then
         landIcePressurePKGActive = .true.
         landIceFluxesPKGActive = .true.
         landIceCouplingPKGActive = .true.
      end if

      !
      ! test for use of subglacial runoff flux, dataSubglacialRunoffFluxPKGActive
      !
      call mpas_pool_get_package(packagePool, 'dataSubglacialRunoffFluxPKGActive', dataSubglacialRunoffFluxPKGActive)
      call mpas_pool_get_config(configPool, 'config_subglacial_runoff_mode', config_subglacial_runoff_mode)
      if ( trim(config_subglacial_runoff_mode) == 'data' ) then
         dataSubglacialRunoffFluxPKGActive = .true.
      end if

      !
      ! test for use of frazil ice formation, frazilIceActive
      !
      call mpas_pool_get_package(packagePool, 'frazilIceActive', frazilIceActive)
      call mpas_pool_get_config(configPool, 'config_use_frazil_ice_formation', config_use_frazil_ice_formation)
      if (config_use_frazil_ice_formation) then
         frazilIceActive = .true.
      end if

      !
      ! test for use of tidal forcing, tidalForcingActive
      !
      call mpas_pool_get_package(packagePool, 'tidalForcingActive', tidalForcingActive)
      call mpas_pool_get_config(configPool, 'config_use_tidal_forcing', config_use_tidal_forcing)
      if (config_use_tidal_forcing) then
         tidalForcingActive = .true.
      end if

      !
      ! test for use of tidal potential forcing, tidalPotentialForcingPKGActive
      !
      call mpas_pool_get_package(packagePool, 'tidalPotentialForcingPKGActive', tidalPotentialForcingPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_tidal_potential_forcing', config_use_tidal_potential_forcing)
      if (config_use_tidal_potential_forcing) then
         tidalPotentialForcingPKGActive = .true.
      end if
      ! Set self-attraction and loading, dependent on tidal potential forcing
      call mpas_pool_get_config(configPool, 'config_use_self_attraction_loading', config_use_self_attraction_loading)

      !
      ! test for use of vegetation drag, vegetationDragPKGActive
      !
      call mpas_pool_get_package(packagePool, 'vegetationDragPKGActive', vegetationDragPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_vegetation_drag', config_use_vegetation_drag)
      if (config_use_vegetation_drag) then
         vegetationDragPKGActive = .true.
      end if

      !
      ! test for use of topographic wave drag, topographicWaveDragPKG
      !
      call mpas_pool_get_package(packagePool, 'topographicWaveDragPKGActive', topographicWaveDragPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_topographic_wave_drag', config_use_topographic_wave_drag)
      if (config_use_topographic_wave_drag) then
         topographicWaveDragPKGActive = .true.
      end if

      !
      ! test for use of gm
      !
      call mpas_pool_get_package(packagePool, 'gmActive', gmActive)
      call mpas_pool_get_config(configPool, 'config_use_GM', config_use_GM)
      call mpas_pool_get_config(configPool, 'config_use_Redi', config_use_Redi)
      if (config_use_GM.or.config_use_Redi) then
         gmActive = .true.
      end if

      !
      ! test for use of gm
      !
      call mpas_pool_get_package(packagePool, 'submesoActive', submesoActive)
      call mpas_pool_get_config(configPool, 'config_submesoscale_enable', config_submesoscale_enable)
      if (config_submesoscale_enable) then
         submesoActive = .true.
      end if

      ! test for time-varying forcing
      !
      call mpas_pool_get_package(packagePool, 'timeVaryingAtmosphericForcingPKGActive', timeVaryingAtmosphericForcingPKGActive)
      call mpas_pool_get_config(configPool, &
                               'config_use_time_varying_atmospheric_forcing', &
                                config_use_time_varying_atmospheric_forcing)
      if (config_use_time_varying_atmospheric_forcing) then
         timeVaryingAtmosphericForcingPKGActive = .true.
      endif

      call mpas_pool_get_package(packagePool, 'timeVaryingLandIceForcingPKGActive', timeVaryingLandIceForcingPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_time_varying_land_ice_forcing', config_use_time_varying_land_ice_forcing)
      if (config_use_time_varying_land_ice_forcing) then
         timeVaryingLandIceForcingPKGActive = .true.
      endif

      !
      ! test for use of GOTM
      !
      call mpas_pool_get_package(packagePool, 'gotmPKGActive', gotmPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_gotm', config_use_gotm)
      if (config_use_gotm) then
         gotmPKGActive = .true.
      end if

      !
      ! test for use of vertical remapping
      !
      call mpas_pool_get_package(packagePool, 'verticalRemapPKGActive', verticalRemapPKGActive)
      call mpas_pool_get_config(configPool, 'config_vert_advection_method', config_vert_advection_method)
      if (config_vert_advection_method == 'remap') then
         verticalRemapPKGActive = .true.
      end if

      ! test for active waves
      !
      call mpas_pool_get_package(packagePool, 'activeWavePKGActive', activeWavePKGActive)
      call mpas_pool_get_config(configPool, &
                               'config_use_active_wave', &
                                config_use_active_wave)
      if (config_use_active_wave) then
         activeWavePKGActive = .true.
      endif

      !
      ! test for use of subgrid wetting and drying
      !
      call mpas_pool_get_package(packagePool, 'subgridWetDryPKGActive', subgridWetDryPKGActive)
      call mpas_pool_get_config(configPool, 'config_use_subgrid_wetting_drying', config_use_subgrid_wetting_drying)
      if (config_use_subgrid_wetting_drying) then
         subgridWetDryPKGActive = .true.
      end if

      !
      ! test for scaling data ice-shelf melt fluxes by the running mean of removed ice runoff
      !
      call mpas_pool_get_package(packagePool, 'scaledDISMFPKGActive', scaledDISMFPKGActive)
      call mpas_pool_get_config(configPool, &
                                'config_scale_dismf_by_removed_ice_runoff', &
                                config_scale_dismf_by_removed_ice_runoff)
      if (config_scale_dismf_by_removed_ice_runoff) then
         scaledDISMFPKGActive = .true.
      end if

      !
      ! call into analysis member driver to set analysis member packages
      !
      call ocn_analysis_setup_packages(configPool, packagePool, iocontext, err_tmp)
      ierr = ior(ierr, err_tmp)


      !
      ! if in init mode, validate configuration
      !
#ifndef EXCLUDE_INIT_MODE
      if ( initModeActive ) then
         call ocn_init_mode_validate_configuration(configPool, packagePool, iocontext, ierr)
      endif
#endif

      !
      ! iterate over tracer groups
      ! each tracer group is toggled on/off using packages
      ! test each package
      !
      call mpas_pool_begin_iteration(packagePool)
      do while ( mpas_pool_get_next_member(packagePool, groupItr) )
         startIndex = index(groupItr % memberName, 'TracersPKG')
         if ( startIndex .ne. 0 ) then
            strLen = len_trim(groupItr % memberName)
            tracerGroupName = groupItr % memberName(1:strLen-9)

            configName = 'config_use_' // trim(tracerGroupName)
            call mpas_pool_get_config(configPool, configName, config_use_tracerGroup)
            if ( config_use_tracerGroup ) then
               packageName = trim(tracerGroupName) // 'PKGActive'
               call mpas_pool_get_package(packagePool, packageName, tracerGroupPKGActive)
               tracerGroupPKGActive = .true.

               configName = 'config_use_' // trim(tracerGroupName) // '_surface_bulk_forcing'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_surface_bulk_forcing)

               if ( config_use_tracerGroup_surface_bulk_forcing ) then
                  packageName = trim(tracerGroupName) // 'BulkRestoringPKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupBulkRestoringPKGActive)
                  tracerGroupBulkRestoringPKGActive = .true.
               end if

               configName = 'config_use_' // trim(tracerGroupName) // '_surface_restoring'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_surface_restoring)

               if ( config_use_tracerGroup_surface_restoring ) then
                  packageName = trim(tracerGroupName) // 'SurfaceRestoringPKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupSurfaceRestoringPKGActive)
                  tracerGroupSurfaceRestoringPKGActive = .true.
               end if

               configName = 'config_use_' // trim(tracerGroupName) // '_interior_restoring'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_interior_restoring)
               if ( config_use_tracerGroup_interior_restoring ) then
                  packageName = trim(tracerGroupName) // 'InteriorRestoringPKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupInteriorRestoringPKGActive)
                  tracerGroupInteriorRestoringPKGActive = .true.
               end if

               configName = 'config_use_' // trim(tracerGroupName) // '_exponential_decay'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_exponential_decay)
               if ( config_use_tracerGroup_exponential_decay ) then
                  packageName = trim(tracerGroupName) // 'ExponentialDecayPKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupExponentialDecayPKGActive)
                  tracerGroupExponentialDecayPKGActive = .true.
               end if

               configName = 'config_use_' // trim(tracerGroupName) // '_idealAge_forcing'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_idealAge_forcing)
               if ( config_use_tracerGroup_idealAge_forcing ) then
                  packageName = trim(tracerGroupName) // 'IdealAgePKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupIdealAgePKGActive)
                  tracerGroupIdealAgePKGActive = .true.
               end if

               configName = 'config_use_' // trim(tracerGroupName) // '_ttd_forcing'
               call mpas_pool_get_config(configPool, configName, config_use_tracerGroup_ttd_forcing)
               if ( config_use_tracerGroup_ttd_forcing ) then
                  packageName = trim(tracerGroupName) // 'TTDPKGActive'
                  call mpas_pool_get_package(packagePool, packageName, tracerGroupTTDPKGActive)
                  tracerGroupTTDPKGActive = .true.
               end if
            end if
         end if
      end do

      !
      ! test for conflicts, i.e. package settings that are inconsistent in combination
      !



      call mpas_log_write( '')
      call mpas_log_write( '  **** Summary of ocean packages ****')
      call mpas_pool_begin_iteration(packagePool)
      do while ( mpas_pool_get_next_member(packagePool, pkgItr) )

         if ( pkgItr % memberType == MPAS_POOL_PACKAGE ) then
            call mpas_pool_get_package(packagePool, pkgItr % memberName, packageActive)
            if ( packageActive ) then
               call mpas_log_write( '      ' // trim(pkgItr % memberName) // ' = ON')
            else
               call mpas_log_write( '      ' // trim(pkgItr % memberName) // ' = OFF')
            end if
         end if
      end do
      call mpas_log_write( '  ***********************************')
      call mpas_log_write( '')

   end function ocn_setup_packages!}}}


   !***********************************************************************
   !
   !  routine ocn_setup_decompositions
   !
   !> \brief   Decomposition setup routine
   !> \author  Doug Jacobsen
   !> \date    04/08/2015
   !> \details
   !>  This routine is intended to create the decomposition list within a
   !>  domain type, and register any decompositons the core wants within it.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_decompositions(decompList) result(ierr)!{{{

      use mpas_derived_types
      use mpas_decomp

      implicit none

      type (mpas_decomp_list), pointer :: decompList

      integer :: ierr
      procedure (mpas_decomp_function), pointer :: decompFunc

      ierr = 0

      call mpas_decomp_create_decomp_list(decompList)

      decompFunc => mpas_uniform_decomp

      call mpas_decomp_register_method(decompList, 'uniform', decompFunc, iErr)

      if ( iErr == MPAS_DECOMP_NOERR ) then
         iErr = 0
      end if

   end function ocn_setup_decompositions!}}}


   !***********************************************************************
   !
   !  function ocn_setup_clock
   !
   !> \brief   Pacakge setup routine
   !> \author  Michael Duda
   !> \date    6 August 2014
   !> \details
   !>  The purpose of this function is to allow the core to set up a simulation
   !>  clock that will be used by the I/O subsystem for timing reads and writes
   !>  of I/O streams.
   !>  This function is called from the superstructure after the framework
   !>  has been initialized but before any fields have been allocated and
   !>  initial fields have been read from input files. However, all namelist
   !>  options are available.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_clock(core_clock, configs) result(ierr)!{{{

      use mpas_derived_types

      implicit none

      type (MPAS_Clock_type), intent(inout) :: core_clock
      type (mpas_pool_type), intent(inout) :: configs
      integer :: ierr

      character(len=StrKIND), pointer :: config_ocean_run_mode

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' ) then
         ierr = ocn_forward_mode_setup_clock(core_clock, configs)
      else if ( trim(config_ocean_run_mode) == 'analysis' ) then
         ierr = ocn_analysis_mode_setup_clock(core_clock, configs)
#ifndef EXCLUDE_INIT_MODE
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         ierr = ocn_init_mode_setup_clock(core_clock, configs)
#endif
      end if


   end function ocn_setup_clock!}}}


   !***********************************************************************
   !
   !  function ocn_setup_log
   !
   !> \brief   Log setup routine
   !> \author  Matt Hoffman
   !> \date    14 February 2017
   !> \details
   !>  The purpose of this routine is to set up the logging manager
   !>  and allow the core to specify details of the configuration.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_log(logInfo, domain, unitNumbers) result(iErr)!{{{

      use mpas_derived_types
      use mpas_log

      implicit none

      type (mpas_log_type), intent(inout), pointer :: logInfo    !< logging information object to set up
      type (domain_type), intent(in), pointer :: domain          !< domain object to provide info for setting up log manager
      integer, dimension(2), intent(in), optional :: unitNumbers !< Fortran unit numbers to use for output and error logs
      integer :: iErr

      ! Local variables
      integer :: local_err

      iErr = 0

      ! Initialize log manager
      call mpas_log_init(logInfo, domain, unitNumbers=unitNumbers, err=local_err)
      iErr = ior(iErr, local_err)

      ! Set core specific options here
      ! (At present, there are not any.  There could eventually be choices about the file naming conventions
      !  or other settings controlling behavior.)

      ! After core has had a chance to modify log defaults, open the output log
      call mpas_log_open(err=local_err)
      iErr = ior(iErr, local_err)

   end function ocn_setup_log!}}}


   !***********************************************************************
   !
   !  function ocn_get_mesh_stream
   !
   !> \brief   Returns the name of the stream containing mesh information
   !> \author  Michael Duda
   !> \date    8 August 2014
   !> \details
   !>  This function returns the name of the I/O stream containing dimensions,
   !>  attributes, and mesh fields needed by the framework bootstrapping
   !>  routine. At the time this routine is called, only namelist options
   !>  are available.
   !
   !-----------------------------------------------------------------------
   function ocn_get_mesh_stream(configs, stream) result(ierr)!{{{

      use mpas_derived_types
      use mpas_pool_routines

      implicit none

      type (mpas_pool_type), intent(inout) :: configs
      character(len=StrKIND), intent(out) :: stream
      integer :: ierr

      character(len=StrKIND), pointer :: config_ocean_run_mode

      ierr = 0

      call mpas_pool_get_config(configs, 'config_ocean_run_mode', config_ocean_run_mode)

      if ( trim(config_ocean_run_mode) == 'forward' .or. trim(config_ocean_run_mode) == 'analysis' ) then
         write(stream,'(a)') 'mesh'
      else if ( trim(config_ocean_run_mode) == 'init' ) then
         write(stream,'(a)') 'input_init'
      end if

   end function ocn_get_mesh_stream!}}}


   !***********************************************************************
   !
   !  function ocn_setup_block
   !
   !> \brief   Ocean block setup function
   !> \author  Doug Jacobsen
   !> \date    03/18/2015
   !> \details
   !>  This function is a wrapper function to properly setup a block to be a
   !>  ocean core block.
   !
   !-----------------------------------------------------------------------
   function ocn_setup_block(block) result(iErr)!{{{
      use mpas_derived_types
      type (block_type), pointer :: block
      integer :: iErr

      iErr = 0
      call ocn_generate_structs(block, block % structs, block % dimensions, block % packages)
   end function ocn_setup_block!}}}

#include "../inc/setup_immutable_streams.inc"

#include "../inc/block_dimension_routines.inc"

#include "../inc/define_packages.inc"

#include "../inc/structs_and_variables.inc"

#include "../inc/namelist_call.inc"

#include "../inc/namelist_defines.inc"

end module ocn_core_interface

