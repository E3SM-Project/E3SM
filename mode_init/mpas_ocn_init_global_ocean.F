! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_global_ocean
!
!> \brief MPAS ocean initialize case -- Global Ocean
!> \author Doug Jacobsen
!> \date   03/04/2014
!> \details
!>  This module contains the routines for initializing the
!>  the global ocean test case
!
!-----------------------------------------------------------------------

module ocn_init_global_ocean

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_io
   use mpas_io_streams
   use mpas_stream_manager
   use mpas_timekeeping
   use mpas_dmpar

   use ocn_constants
   use ocn_equation_of_state
   use ocn_init_cell_markers
   use ocn_init_vertical_grids
   use ocn_init_interpolation
   use ocn_init_ssh_and_landIcePressure

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

   public :: ocn_init_setup_global_ocean, &
             ocn_init_validate_global_ocean

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   ! 2D input variables.
   ! Note Tracer and ecosys variables may have different original grids.
   integer :: nLatTracer, nLonTracer, nDepthTracer
   integer :: nLonEcosys, nLatEcosys, nDepthEcosys
   integer :: nDepthOutput, nTimes

   ! 3D input variables
   integer :: nLatWind, nLonWind
   integer :: nLatTopo, nLonTopo
   integer :: nLonSW, nLatSW
   integer :: nLatLandIceThk, nLonLandIceThk

   type (field1DReal) :: depthOutput
   type (field1DReal) :: tracerLat, tracerLon, tracerDepth
   type (field1DReal) :: windLat, windLon
   type (field1DReal) :: topoLat, topoLon
   type (field1DReal) :: swDataLat, swDataLon
   type (field1DReal) :: landIceThkLat, landIceThkLon

   type (field2DReal) :: topoIC, zonalWindIC, meridionalWindIC, chlorophyllIC, zenithAngleIC, clearSkyIC
   type (field2DReal) :: landIceThkIC, landIceDraftIC
   type (field2DReal) :: oceanFracIC, landIceFracIC, groundedFracIC

   type (field3DReal) :: tracerIC, ecosysForcingIC

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean
!
!> \brief   Setup for global ocean test case
!> \author  Mark Petersen, Doug Jacobsen
!> \date    08/23/2016
!> \details
!>  This routine sets up the initial conditions for the global ocean test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_global_ocean(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: meshPool, forcingPool, statePool, tracersPool, scratchPool
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration, config_sw_absorption_type, &
                                          config_global_ocean_ecosys_file, config_global_ocean_ecosys_forcing_file

      character (len=StrKIND) :: interpTracerName

      type (block_type), pointer :: block_ptr

      logical, pointer :: config_global_ocean_cull_inland_seas
      logical, pointer :: config_global_ocean_depress_by_land_ice

      logical, pointer :: on_a_sphere

      type (mpas_pool_type), pointer :: ecosysAuxiliary  ! additional ecosys fields

      real (kind=RKIND), dimension(:), pointer :: PH_PREV, PH_PREV_ALT_CO2, pocToSed
      real (kind=RKIND), dimension(:, :), pointer :: PH_PREV_3D, PH_PREV_ALT_CO2_3D
      real (kind=RKIND), dimension(:, :), pointer :: FESEDFLUX

      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND), dimension(:, :, :), pointer :: ecosysTracers, activeTracers, debugTracers, &
           DMSTracers, MacroMoleculesTracers
      integer, pointer :: nVertLevels, nCellsSolve, tracerIndex
      integer :: iCell, k, iTracer
      integer, dimension(3) :: indexField

      type (field2DReal), pointer :: interpActiveTracerField, interpEcosysTracerField, &
           interpActiveTracerSmoothField, interpEcosysTracerSmoothField

      type (field3DReal), pointer :: ecosysTracersField

      character (len=StrKIND) :: fieldName, poolName

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      call mpas_pool_get_config(domain % configs, 'config_sw_absorption_type', config_sw_absorption_type)
      call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_file', config_global_ocean_ecosys_file)
      call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_forcing_file', &
                                                   config_global_ocean_ecosys_forcing_file)

      if (trim(config_init_configuration) /= "global_ocean") return

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      ! the following are to test if each is associated before calling init routines
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
      call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', ecosysTracers, 1)
      call mpas_pool_get_array(tracersPool, 'DMSTracers', DMSTracers, 1)
      call mpas_pool_get_array(tracersPool, 'MacroMoleculesTracers', MacroMoleculesTracers, 1)

      if ( .not. on_a_sphere ) call mpas_dmpar_global_abort('MPAS-ocean: ERROR: The global ocean configuration can ' &
                                                         // 'only be applied to a spherical mesh. Exiting...')

      call mpas_pool_get_config(domain % configs, 'config_global_ocean_cull_inland_seas', &
            config_global_ocean_cull_inland_seas)
      call mpas_pool_get_config(domain % configs, 'config_global_ocean_depress_by_land_ice', &
            config_global_ocean_depress_by_land_ice)

      !***********************************************************************
      !
      !  Topography
      !
      !***********************************************************************

      write(stderrUnit,*) 'Reading depth levels.'
      call ocn_init_setup_global_ocean_read_depth_levels(domain, iErr)

      write(stderrUnit,*) 'Reading topography data.'
      call ocn_init_setup_global_ocean_read_topo(domain, iErr)
      write(stderrUnit,*) 'Interpolating topography data.'
      call ocn_init_setup_global_ocean_create_model_topo(domain, iErr)
      write(stderrUnit,*) 'Cleaning up topography IC fields'
      call ocn_init_global_ocean_destroy_topo_fields()

      !***********************************************************************
      !
      !  Land ice depression
      !
      !***********************************************************************

      if (config_global_ocean_depress_by_land_ice) then
         write(stderrUnit,*) 'Reading land ice topography data.'
         call ocn_init_setup_global_ocean_read_land_ice_topography(domain, iErr)
         write(stderrUnit,*) 'Interpolating land ice topography data.'
         call ocn_init_setup_global_ocean_interpolate_land_ice_topography(domain, iErr)
      end if

      write(stderrUnit,*) 'Initializing vertical coordinate with ssh = 0.'
      ! compute the vertical grid (layerThickness, restingThickness, maxLevelCell, zMid)
      ! based on bottomDepth and refBottomDepth and apply PBCS if requested
      call ocn_init_ssh_and_landIcePressure_vertical_grid(domain, iErr)

      if(iErr .ne. 0) then
        write(stderrUnit,*) 'ERROR: ocn_init_ssh_and_landIcePressure_vertical_grid failed.'
        call mpas_dmpar_finalize(domain % dminfo)
      end if

      !***********************************************************************
      !
      !  Shortwave data
      !
      !***********************************************************************

      if(trim(config_sw_absorption_type) == 'ohlmann00') then

         write(stderrUnit,*) 'Reading penetrating shortwave lat/lon data'
         call ocn_init_setup_global_ocean_read_swData_lat_lon(domain,iErr)
         write(stderrUnit,*) 'Interpolating penetrating shortwave data'
         call ocn_init_setup_global_ocean_interpolate_swData(domain,iErr)
         write(stderrUnit,*) 'Cleaning penetrating shortwave data'
         call ocn_init_global_ocean_destroy_swData_fields()

      endif

      !***********************************************************************
      !
      !  Active tracers (temperature and salinity)
      !
      !***********************************************************************

      write(stderrUnit,*) 'Reading tracer Lat/Lon coordinates'
      call ocn_init_setup_global_ocean_read_tracer_lat_lon(domain, iErr)

      allocate(tracerIC % attLists(1))
      allocate(tracerIC % array(nLonTracer, nLatTracer, nDepthTracer))

      call mpas_pool_get_field(scratchPool, 'interpActiveTracer', interpActiveTracerField)
      call mpas_allocate_scratch_field(interpActiveTracerField, .false.)
      call mpas_pool_get_field(scratchPool, 'interpActiveTracerSmooth', interpActiveTracerSmoothField)
      call mpas_allocate_scratch_field(interpActiveTracerSmoothField, .false.)

      interpTracerName = 'interpActiveTracer'

      fieldName = 'temperature'
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', tracerIndex)
      call ocn_init_setup_global_ocean_read_temperature(domain, iErr)
      call ocn_init_setup_global_ocean_interpolate_tracers(domain, activeTracers, tracerIndex, interpTracerName, iErr)

      fieldName = 'salinity'
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', tracerIndex)
      call ocn_init_setup_global_ocean_read_salinity(domain, iErr)
      call ocn_init_setup_global_ocean_interpolate_tracers(domain, activeTracers, tracerIndex,  interpTracerName, iErr)

      deallocate(tracerIC % array)
      deallocate(tracerIC % attLists)
      call mpas_deallocate_scratch_field(interpActiveTracerField, .false.)
      call mpas_deallocate_scratch_field(interpActiveTracerSmoothField, .false.)
      call ocn_init_global_ocean_destroy_tracer_fields()

      !***********************************************************************
      !
      !  Debug tracers
      !
      !***********************************************************************

      call mpas_pool_get_dimension(tracersPool, 'index_tracer1', tracerIndex)
      if(associated(debugTracers)) then
         debugTracers(tracerIndex,:,:) = 1.0_RKIND
      end if

      !***********************************************************************
      !
      !  Ecosystem tracers
      !
      !***********************************************************************

      if ( associated(ecosysTracers) ) then

         write(stderrUnit,*) 'Reading ecosys lat/lon data'
         call ocn_init_setup_global_ocean_read_ecosys_lat_lon(domain,iErr)

         write(stderrUnit,*) 'Reading ecosys IC.'
         allocate(tracerIC % attLists(1))
         allocate(tracerIC % array(nLonEcosys, nLatEcosys, nDepthEcosys))

         call mpas_pool_get_field(scratchPool, 'interpEcosysTracer', interpEcosysTracerField)
         call mpas_allocate_scratch_field(interpEcosysTracerField, .false.)
         call mpas_pool_get_field(scratchPool, 'interpEcosysTracerSmooth', interpEcosysTracerSmoothField)
         call mpas_allocate_scratch_field(interpEcosysTracerSmoothField, .false.)

         interpTracerName = 'interpEcosysTracer'

         call mpas_pool_get_field(tracersPool, 'ecosysTracers', ecosysTracersField, 1)
         do iTracer = 1, size( ecosysTracersField % constituentNames)
            call ocn_init_setup_global_ocean_read_ecosys(domain, ecosysTracersField % constituentNames(iTracer), &
                 config_global_ocean_ecosys_file, iErr)
            call ocn_init_setup_global_ocean_interpolate_tracers(domain, ecosysTracers, iTracer, interpTracerName, iErr)
         end do

         deallocate(tracerIC % array)
         call mpas_deallocate_scratch_field(interpEcosysTracerField, .false.)
         call mpas_deallocate_scratch_field(interpEcosysTracerSmoothField, .false.)

         !***********************************************************************
         !
         !  Ecosystem forcing tracers
         !
         !***********************************************************************

         write(stderrUnit,*) 'Reading ecosys forcing.'
         allocate(ecosysForcingIC % attLists(1))
         allocate(ecosysForcingIC % array(nLonEcosys, nLatEcosys, 1))

         fieldName = 'dust_FLUX_IN'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'IRON_FLUX_IN'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'depositionFluxNO3'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'depositionFluxNH4'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'pocToSed'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxNO3'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxPO4'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxSiO3'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxFe'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxDOC'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxDON'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxDOP'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxDIC'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'riverFluxALK'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'atmosphericCO2'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'atmosphericCO2_ALT_CO2'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'windSpeedSquared10m'; poolName = 'ecosysAuxiliary'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)

         fieldName = 'iceFraction'; poolName = 'forcing'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'atmosphericPressure'; poolName = 'forcing'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)
         fieldName = 'shortWaveHeatFlux'; poolName = 'forcing'
         call ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName,  &
              config_global_ocean_ecosys_forcing_file, iErr)
         call ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)

         write(stderrUnit,*) 'Cleaning up ecosys IC fields'
         deallocate(ecosysForcingIC % array)
         call ocn_init_global_ocean_destroy_ecosys_fields()

         !***********************************************************************
         !
         !  Ecosystem pH arrays
         !
         !***********************************************************************
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
            call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

            call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV', PH_PREV)
            call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2', PH_PREV_ALT_CO2)
            call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_3D', PH_PREV_3D)
            call mpas_pool_get_array(ecosysAuxiliary, 'PH_PREV_ALT_CO2_3D', PH_PREV_ALT_CO2_3D)
            call mpas_pool_get_array(ecosysAuxiliary, 'FESEDFLUX', FESEDFLUX)
            call mpas_pool_get_array(ecosysAuxiliary, 'pocToSed', pocToSed)

            do iCell = 1, nCellsSolve
               PH_PREV(iCell) = 8.0_RKIND
               PH_PREV_ALT_CO2(iCell) = 8.0_RKIND
               do k = 1, nVertLevels
                  PH_PREV_3D(k, iCell) = 8.0_RKIND
                  PH_PREV_ALT_CO2_3D(k, iCell) = 8.0_RKIND

                  if(maxLevelCell(iCell) == k) then
                     FESEDFLUX(k, iCell) = pocToSed(iCell)*6.8e-4_RKIND
                  else
                     FESEDFLUX(k, iCell) = 0.0_RKIND
                  end if
               end do
            end do

            block_ptr => block_ptr % next
         end do

      end if

      if ( associated(DMSTracers) ) then
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

            call mpas_pool_get_dimension(tracersPool, 'index_DMS', tracerIndex)
            indexField(1) = tracerIndex
            call mpas_pool_get_dimension(tracersPool, 'index_DMSP', tracerIndex)
            indexField(2) = tracerIndex
            do iCell = 1, nCellsSolve
               do k = 1, nVertLevels
                  DMSTracers(indexField(1), k, iCell) = 0.0_RKIND
                  DMSTracers(indexField(2), k, iCell) = 0.0_RKIND
               end do
            end do
            block_ptr => block_ptr % next
         end do
      end if

      if ( associated(MacroMoleculesTracers) ) then
         block_ptr => domain % blocklist
         do while(associated(block_ptr))
            call mpas_pool_get_subpool(domain % blocklist % structs, 'state', statePool)
            call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

            call mpas_pool_get_dimension(tracersPool, 'index_PROT', tracerIndex)
            indexField(1) = tracerIndex
            call mpas_pool_get_dimension(tracersPool, 'index_POLY', tracerIndex)
            indexField(2) = tracerIndex
            call mpas_pool_get_dimension(tracersPool, 'index_LIP', tracerIndex)
            indexField(3) = tracerIndex
            do iCell = 1, nCellsSolve
               do k = 1, nVertLevels
                  MacroMoleculesTracers(indexField(1), k, iCell) = 0.0_RKIND
                  MacroMoleculesTracers(indexField(2), k, iCell) = 0.0_RKIND
                  MacroMoleculesTracers(indexField(3), k, iCell) = 0.0_RKIND
               end do
            end do
            block_ptr => block_ptr % next
         end do
      end if

      write(stderrUnit,*) 'Reading windstress IC.'
      call ocn_init_setup_global_ocean_read_windstress(domain, iErr)
      write(stderrUnit,*) 'Interpolating windstress.'
      call ocn_init_setup_global_ocean_interpolate_windstress(domain, iErr)
      write(stderrUnit,*) 'Destroying windstress fields'
      call ocn_init_global_ocean_destroy_windstress_fields()

      if (config_global_ocean_depress_by_land_ice) then
         write(stderrUnit,*) 'Recalculating ocean layer topography due to land ice depression'
         ! compute or update the land-ice pressure (or possibly SSH), also computing density along the way
         ! If this is the initial guess, the vertical grid and activeTracers may also be recomputed based on SSH
         call ocn_init_ssh_and_landIcePressure_balance(domain, iErr)

         if(iErr .ne. 0) then
            write(stderrUnit,*) 'ERROR: ocn_init_ssh_and_landIcePressure_balance failed.'
            call mpas_dmpar_finalize(domain % dminfo)
         end if

         write(stderrUnit,*) 'Cleaning up land ice topography IC fields'
         call ocn_init_global_ocean_destroy_land_ice_topography_fields()
      end if

      write(stderrUnit,*) 'Copying restoring fields'
      ! this occurs after ocn_init_ssh_and_landIcePressure_balance because activeTracers may have been remapped
      ! to a new vertical coordinate
      call ocn_init_setup_global_ocean_interpolate_restoring(domain, iErr)

      write(stderrUnit,*) 'Compute Haney number'
      call ocn_compute_Haney_number(domain, iErr)
      if(iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: ocn_compute_Haney_number failed.'
         call mpas_dmpar_finalize(domain % dminfo)
      end if

      if (config_global_ocean_cull_inland_seas) then
         write(stderrUnit,*) 'Removing inland seas.'
         call ocn_init_setup_global_ocean_cull_inland_seas(domain, iErr)
      end if

      block_ptr => domain % blocklist
      do while (associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

         call ocn_mark_maxlevelcell(meshPool, iErr)
         block_ptr => block_ptr % next
      end do

      !--------------------------------------------------------------------

    end subroutine ocn_init_setup_global_ocean!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_topo
!
!> \brief   Read the topography IC file
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the topography IC file, including latitude and longitude
!>   information for topography data.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: topographyStream

       character (len=StrKIND), pointer :: config_global_ocean_topography_file, &
                                           config_global_ocean_topography_lat_varname, &
                                           config_global_ocean_topography_nlat_dimname, &
                                           config_global_ocean_topography_lon_varname, &
                                           config_global_ocean_topography_nlon_dimname, &
                                           config_global_ocean_topography_varname, &
                                           config_global_ocean_topography_ocean_frac_varname

       logical, pointer :: config_global_ocean_topography_latlon_degrees, config_global_ocean_topography_has_ocean_frac

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_file', &
                                 config_global_ocean_topography_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_lat_varname', &
                                 config_global_ocean_topography_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_nlat_dimname', &
                                 config_global_ocean_topography_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_lon_varname', &
                                 config_global_ocean_topography_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_nlon_dimname', &
                                 config_global_ocean_topography_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_varname', &
                                 config_global_ocean_topography_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_latlon_degrees', &
                                 config_global_ocean_topography_latlon_degrees)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_ocean_frac_varname', &
                                 config_global_ocean_topography_ocean_frac_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_has_ocean_frac', &
                                 config_global_ocean_topography_has_ocean_frac)

       ! Define stream for depth levels
       call MPAS_createStream(topographyStream, domain % iocontext, config_global_ocean_topography_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup topoLat, topoLon, and topoIC fields for stream to be read in
       topoLat % fieldName = trim(config_global_ocean_topography_lat_varname)
       topoLat % dimSizes(1) = nLatTopo
       topoLat % dimNames(1) = trim(config_global_ocean_topography_nlat_dimname)
       topoLat % isVarArray = .false.
       topoLat % isPersistent = .true.
       topoLat % isActive = .true.
       topoLat % hasTimeDimension = .false.
       topoLat % block => domain % blocklist
       allocate(topoLat % attLists(1))
       allocate(topoLat % array(nLatTopo))

       topoLon % fieldName = trim(config_global_ocean_topography_lon_varname)
       topoLon % dimSizes(1) = nLonTopo
       topoLon % dimNames(1) = trim(config_global_ocean_topography_nlon_dimname)
       topoLon % isVarArray = .false.
       topoLon % isPersistent = .true.
       topoLon % isActive = .true.
       topoLon % hasTimeDimension = .false.
       topoLon % block => domain % blocklist
       allocate(topoLon % attLists(1))
       allocate(topoLon % array(nLonTopo))

       topoIC % fieldName = trim(config_global_ocean_topography_varname)
       topoIC % dimSizes(1) = nLonTopo
       topoIC % dimSizes(2) = nLatTopo
       topoIC % dimNames(1) = trim(config_global_ocean_topography_nlon_dimname)
       topoIC % dimNames(2) = trim(config_global_ocean_topography_nlat_dimname)
       topoIC % isVarArray = .false.
       topoIC % isPersistent = .true.
       topoIC % isActive = .true.
       topoIC % hasTimeDimension = .false.
       topoIC % block => domain % blocklist
       allocate(topoIC % attLists(1))
       allocate(topoIC % array(nLonTopo, nLatTopo))

       ! Add topoLat, topoLon, and topoIC fields to stream
       call MPAS_streamAddField(topographyStream, topoLat, iErr)
       call MPAS_streamAddField(topographyStream, topoLon, iErr)
       call MPAS_streamAddField(topographyStream, topoIC, iErr)

       if(config_global_ocean_topography_has_ocean_frac) then
          oceanFracIC % fieldName = trim(config_global_ocean_topography_ocean_frac_varname)
          oceanFracIC % dimSizes(1) = nLonTopo
          oceanFracIC % dimSizes(2) = nLatTopo
          oceanFracIC % dimNames(1) = trim(config_global_ocean_topography_nlon_dimname)
          oceanFracIC % dimNames(2) = trim(config_global_ocean_topography_nlat_dimname)
          oceanFracIC % isVarArray = .false.
          oceanFracIC % isPersistent = .true.
          oceanFracIC % isActive = .true.
          oceanFracIC % hasTimeDimension = .false.
          oceanFracIC % block => domain % blocklist
          allocate(oceanFracIC % attLists(1))
          allocate(oceanFracIC % array(nLonTopo, nLatTopo))

          call MPAS_streamAddField(topographyStream, oceanFracIC, iErr)
       end if

       ! Read stream
       call MPAS_readStream(topographyStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(topographyStream)

       if (config_global_ocean_topography_latlon_degrees) then
          topoLat % array(:) = topoLat % array(:) * pii / 180.0_RKIND
          topoLon % array(:) = topoLon % array(:) * pii / 180.0_RKIND
       end if

    end subroutine ocn_init_setup_global_ocean_read_topo!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_land_ice_topography
!
!> \brief   Read the ice sheet thickness IC file
!> \author  Jeremy Fyke, Xylar Asay-Davis, Mark Petersen (modified from Doug Jacobsen code)
!> \date    06/15/2015
!> \details
!>  This routine reads the ice sheet topography IC file, including latitude and longitude
!>   information for ice sheet topography data.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_land_ice_topography(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: landIceThicknessStream

       character (len=StrKIND), pointer :: config_global_ocean_land_ice_topo_file, config_global_ocean_land_ice_topo_lat_varname, &
            config_global_ocean_land_ice_topo_draft_varname, &
            config_global_ocean_land_ice_topo_ice_frac_varname, &
            config_global_ocean_land_ice_topo_grounded_frac_varname, &
            config_global_ocean_land_ice_topo_nlat_dimname, config_global_ocean_land_ice_topo_lon_varname, &
            config_global_ocean_land_ice_topo_nlon_dimname, config_global_ocean_land_ice_topo_thickness_varname

       logical, pointer :: config_global_ocean_land_ice_topo_latlon_degrees

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_file', &
            config_global_ocean_land_ice_topo_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_lat_varname', &
            config_global_ocean_land_ice_topo_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_nlat_dimname', &
            config_global_ocean_land_ice_topo_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_lon_varname', &
            config_global_ocean_land_ice_topo_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_nlon_dimname', &
            config_global_ocean_land_ice_topo_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_thickness_varname', &
            config_global_ocean_land_ice_topo_thickness_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_draft_varname', &
            config_global_ocean_land_ice_topo_draft_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_ice_frac_varname', &
            config_global_ocean_land_ice_topo_ice_frac_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_grounded_frac_varname', &
            config_global_ocean_land_ice_topo_grounded_frac_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_land_ice_topo_latlon_degrees', &
            config_global_ocean_land_ice_topo_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(landIceThicknessStream, domain % iocontext, config_global_ocean_land_ice_topo_file, &
            MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup landIceThkLat, landIceThkLon, and landIceThkIC fields for stream to be read in
       landIceThkLat % fieldName = trim(config_global_ocean_land_ice_topo_lat_varname)
       landIceThkLat % dimSizes(1) = nLatLandIceThk
       landIceThkLat % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlat_dimname)
       landIceThkLat % isVarArray = .false.
       landIceThkLat % isPersistent = .true.
       landIceThkLat % isActive = .true.
       landIceThkLat % hasTimeDimension = .false.
       landIceThkLat % block => domain % blocklist
       allocate(landIceThkLat % attLists(1))
       allocate(landIceThkLat % array(nLatLandIceThk))

       landIceThkLon % fieldName = trim(config_global_ocean_land_ice_topo_lon_varname)
       landIceThkLon % dimSizes(1) = nLonLandIceThk
       landIceThkLon % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlon_dimname)
       landIceThkLon % isVarArray = .false.
       landIceThkLon % isPersistent = .true.
       landIceThkLon % isActive = .true.
       landIceThkLon % hasTimeDimension = .false.
       landIceThkLon % block => domain % blocklist
       allocate(landIceThkLon % attLists(1))
       allocate(landIceThkLon % array(nLonLandIceThk))

       landIceThkIC % fieldName = trim(config_global_ocean_land_ice_topo_thickness_varname)
       landIceThkIC % dimSizes(1) = nLonLandIceThk
       landIceThkIC % dimSizes(2) = nLatLandIceThk
       landIceThkIC % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlon_dimname)
       landIceThkIC % dimNames(2) = trim(config_global_ocean_land_ice_topo_nlat_dimname)
       landIceThkIC % isVarArray = .false.
       landIceThkIC % isPersistent = .true.
       landIceThkIC % isActive = .true.
       landIceThkIC % hasTimeDimension = .false.
       landIceThkIC % block => domain % blocklist
       allocate(landIceThkIC % attLists(1))
       allocate(landIceThkIC % array(nLonLandIceThk, nLatLandIceThk))

       landIceDraftIC % fieldName = trim(config_global_ocean_land_ice_topo_draft_varname)
       landIceDraftIC % dimSizes(1) = nLonLandIceThk
       landIceDraftIC % dimSizes(2) = nLatLandIceThk
       landIceDraftIC % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlon_dimname)
       landIceDraftIC % dimNames(2) = trim(config_global_ocean_land_ice_topo_nlat_dimname)
       landIceDraftIC % isVarArray = .false.
       landIceDraftIC % isPersistent = .true.
       landIceDraftIC % isActive = .true.
       landIceDraftIC % hasTimeDimension = .false.
       landIceDraftIC % block => domain % blocklist
       allocate(landIceDraftIC % attLists(1))
       allocate(landIceDraftIC % array(nLonLandIceThk, nLatLandIceThk))

       landIceFracIC % fieldName = trim(config_global_ocean_land_ice_topo_ice_frac_varname)
       landIceFracIC % dimSizes(1) = nLonLandIceThk
       landIceFracIC % dimSizes(2) = nLatLandIceThk
       landIceFracIC % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlon_dimname)
       landIceFracIC % dimNames(2) = trim(config_global_ocean_land_ice_topo_nlat_dimname)
       landIceFracIC % isVarArray = .false.
       landIceFracIC % isPersistent = .true.
       landIceFracIC % isActive = .true.
       landIceFracIC % hasTimeDimension = .false.
       landIceFracIC % block => domain % blocklist
       allocate(landIceFracIC % attLists(1))
       allocate(landIceFracIC % array(nLonLandIceThk, nLatLandIceThk))

       groundedFracIC % fieldName = trim(config_global_ocean_land_ice_topo_grounded_frac_varname)
       groundedFracIC % dimSizes(1) = nLonLandIceThk
       groundedFracIC % dimSizes(2) = nLatLandIceThk
       groundedFracIC % dimNames(1) = trim(config_global_ocean_land_ice_topo_nlon_dimname)
       groundedFracIC % dimNames(2) = trim(config_global_ocean_land_ice_topo_nlat_dimname)
       groundedFracIC % isVarArray = .false.
       groundedFracIC % isPersistent = .true.
       groundedFracIC % isActive = .true.
       groundedFracIC % hasTimeDimension = .false.
       groundedFracIC % block => domain % blocklist
       allocate(groundedFracIC % attLists(1))
       allocate(groundedFracIC % array(nLonLandIceThk, nLatLandIceThk))

       ! Add landIceThkLat, landIceThkLon, and landIceThkIC fields to stream
       call MPAS_streamAddField(landIceThicknessStream, landIceThkLat, iErr)
       call MPAS_streamAddField(landIceThicknessStream, landIceThkLon, iErr)
       call MPAS_streamAddField(landIceThicknessStream, landIceThkIC, iErr)
       call MPAS_streamAddField(landIceThicknessStream, landIceDraftIC, iErr)
       call MPAS_streamAddField(landIceThicknessStream, landIceFracIC, iErr)
       call MPAS_streamAddField(landIceThicknessStream, groundedFracIC, iErr)

       ! Read stream
       call MPAS_readStream(landIceThicknessStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(landIceThicknessStream)

       if (config_global_ocean_land_ice_topo_latlon_degrees) then
          landIceThkLat % array(:) = landIceThkLat % array(:) * pii / 180.0_RKIND
          landIceThkLon % array(:) = landIceThkLon % array(:) * pii / 180.0_RKIND
       end if

    end subroutine ocn_init_setup_global_ocean_read_land_ice_topography!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_create_model_topo
!
!> \brief   Interpolate the topography IC to MPAS mesh
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine interpolates topography data to the MPAS mesh.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_create_model_topo(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, scratchPool, verticalMeshPool, criticalPassagesPool

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, bottomDepth, bottomDepthObserved, &
            refBottomDepth, refLayerThickness, refZMid, oceanFracObserved

       integer, pointer :: nCells, nCellsSolve, nVertLevels

       type (field1DInteger), pointer :: maxLevelCellField, smoothedLevelsField
       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       integer :: iCell, coc, j, k, maxLevel
       integer :: minimum_levels

       character (len=StrKIND), pointer :: config_global_ocean_topography_method
       logical, pointer :: config_global_ocean_smooth_topography, config_global_ocean_topography_has_ocean_frac, &
                           config_global_ocean_deepen_critical_passages
       real (kind=RKIND), pointer :: config_global_ocean_minimum_depth

       logical :: isOcean

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_method', config_global_ocean_topography_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_minimum_depth', config_global_ocean_minimum_depth)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_smooth_topography', config_global_ocean_smooth_topography)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_deepen_critical_passages', &
                                 config_global_ocean_deepen_critical_passages)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_has_ocean_frac', &
                                 config_global_ocean_topography_has_ocean_frac)

       call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(meshPool, 'bottomDepthObserved', bottomDepthObserved)
          call mpas_pool_get_array(meshPool, 'oceanFracObserved', oceanFracObserved)
          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)


          do k = 1, nVertLevels
             if (refBottomDepth(k).gt.config_global_ocean_minimum_depth) then
                minimum_levels = k
                write (stdoutUnit,'(a,f8.2,2a,i5,a,f8.2,a)') 'config_global_ocean_minimum_depth=', &
                   config_global_ocean_minimum_depth,' m.  ', 'Setting minimum layer index to ', &
                   minimum_levels, ' with a bottom depth of ', refBottomDepth(k), ' m.'
                exit
             end if
          end do

          ! Record depth of the bottom of the ocean, before any alterations for modeling purposes.
          if (config_global_ocean_topography_method .eq. "nearest_neighbor") then

             call ocn_init_interpolation_nearest_horiz(topoLon % array, topoLat % array, &
                                                       topoIC % array, nLonTopo, nLatTopo, &
                                                       lonCell, latCell, bottomDepthObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             if (config_global_ocean_topography_has_ocean_frac) then
                call ocn_init_interpolation_nearest_horiz(topoLon % array, topoLat % array, &
                                                          oceanFracIC % array, nLonTopo, nLatTopo, &
                                                          lonCell, latCell, oceanFracObserved, nCells, &
                                                          inXPeriod = 2.0_RKIND * pii)
             end if

          elseif (config_global_ocean_topography_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(topoLon % array, topoLat % array, &
                                                        topoIC % array, nLonTopo, nLatTopo, &
                                                        lonCell, latCell, bottomDepthObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)

             if (config_global_ocean_topography_has_ocean_frac) then
                call ocn_init_interpolation_bilinear_horiz(topoLon % array, topoLat % array, &
                                                           oceanFracIC % array, nLonTopo, nLatTopo, &
                                                           lonCell, latCell, oceanFracObserved, nCells, &
                                                           inXPeriod = 2.0_RKIND * pii)
             end if

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_topography_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          do iCell = 1, nCells
             ! Record depth of the bottom of the ocean, before any alterations for modeling purposes.
             ! Flip the sign to positive down.
             bottomDepthObserved(iCell) = -bottomDepthObserved(iCell)
             !isOcean = bottomDepthObserved(iCell) > 0.0_RKIND
             isOcean = .true.
             if (config_global_ocean_topography_has_ocean_frac) then
                ! if there is an ocean-fraction field, mark cells that are < 50% ocean as land
                isOcean = isOcean .and. (oceanFracObserved(iCell) >= 0.5_RKIND)
             end if
             if (isOcean) then
                ! Enforce minimum depth
                bottomDepth(iCell) = max(bottomDepthObserved(iCell), refBottomDepth(minimum_levels))

                maxLevelCell(iCell) = -1
                do k = 1, nVertLevels
                   if (refBottomDepth(k) >= bottomDepth(iCell)) then
                      maxLevelCell(iCell) = k
                      exit
                   end if
                end do

                if (maxLevelCell(iCell) == -1) then
                   maxLevelCell(iCell) = nVertLevels
                   bottomDepth(iCell) = refBottomDepth( nVertLevels )
                end if

             else
                bottomDepth(iCell) = 0.0_RKIND
                maxLevelCell(iCell) = -1
             end if
          end do

          if (config_global_ocean_deepen_critical_passages) then
             call mpas_pool_get_subpool(block_ptr % structs, 'criticalPassages', criticalPassagesPool)
             call ocn_init_setup_global_ocean_deepen_critical_passages(meshPool, criticalPassagesPool, iErr)
          end if

          ! Smooth depth levels. Enforce different in maxLevelCell to only be a maximum
          ! of 1 vertical level between two neighboring cells.
          if (config_global_ocean_smooth_topography) then
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_field(scratchPool, 'smoothedLevels', smoothedLevelsField)

             call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
             call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

             call mpas_allocate_scratch_field(smoothedLevelsField, .true.)

             maxLevelCell(nCells+1) = -1
             smoothedLevelsField % array = maxLevelCell

             do iCell = 1, nCellsSolve
                maxLevel = 0
                do j = 1, nEdgesOnCell(iCell)
                   coc = cellsOnCell(j, iCell)
                   maxLevel = max(maxLevel, maxLevelCell(coc))
                end do

                if (maxLevel < maxLevelCell(iCell) ) then
                   smoothedLevelsField % array(iCell) = maxLevel + 1
                   bottomDepth(iCell) = refBottomDepth(maxLevel + 1)
                end if
             end do

             maxLevelCell(:) = smoothedLevelsField % array(:)

             call mpas_deallocate_scratch_field(smoothedLevelsField, .true.)
          end if

          block_ptr => block_ptr % next
       end do

       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_field(meshPool, 'maxLevelCell', maxLevelCellField)
       call mpas_dmpar_exch_halo_field(maxLevelCellField)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

          call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)
          call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)

          ! Compute refLayerThickness and refZMid
          call ocn_compute_layerThickness_zMid_from_bottomDepth(refLayerThickness,refZMid, &
               refBottomDepth,refBottomDepth(nVertLevels), &
               nVertLevels,nVertLevels,iErr)

          block_ptr => block_ptr % next
       end do

     end subroutine ocn_init_setup_global_ocean_create_model_topo!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_deepen_critical_passages
!
!> \brief   Deepen critical passages in model grid so that relevant seas are included.
!> \author  Xylar Asay-Davis
!> \date    5 April 2016
!> \details
!>   Deepen cirtical passages in model grid so that relevant seas are included.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_deepen_critical_passages(meshPool, criticalPassagesPool, iErr)!{{{
       type (mpas_pool_type), intent(inout) :: meshPool, criticalPassagesPool
       integer, intent(out) :: iErr

       integer, pointer :: nCells, nTransects, nVertLevels

       real(kind=RKIND), dimension(:), pointer :: bottomDepth, criticalPassageDepths, &
                                                  refBottomDepth
       integer, dimension(:,:), pointer :: criticalPassageMasks
       integer, dimension(:), pointer :: maxLevelCell

       integer, dimension(:), pointer :: criticalPassageLevel

       integer :: iCell, iTransect, k

       iErr = 0

       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

       call mpas_pool_get_dimension(criticalPassagesPool, 'nTransects', nTransects)

       call mpas_pool_get_array(criticalPassagesPool, 'transectCellMasks', criticalPassageMasks)
       call mpas_pool_get_array(criticalPassagesPool, 'depthTransects', criticalPassageDepths)
       call mpas_pool_get_array(criticalPassagesPool, 'criticalPassageLevel', criticalPassageLevel)

       do iTransect = 1,nTransects
          do k=1,nVertLevels
             criticalPassageLevel(iTransect) = nVertLevels
             if(refBottomDepth(k) > criticalPassageDepths(iTransect)) then
                criticalPassageLevel(iTransect) = k
                exit
             end if
          end do
       end do

       do iCell = 1,nCells
          do iTransect = 1,nTransects
             if(criticalPassageMasks(iTransect, iCell) == 0) cycle
             k = criticalPassageLevel(iTransect)
             if(bottomDepth(iCell) < refBottomDepth(k)) then
                bottomDepth(iCell) = refBottomDepth(k)
                maxLevelCell(iCell) = k
             end if
          end do
       end do

    end subroutine ocn_init_setup_global_ocean_deepen_critical_passages!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_land_ice_topography
!
!> \brief   Interpolate the topography IC to MPAS mesh
!> \author  Jeremy Fyke, Xylar Asay-Davis, Mark Petersen
!> \date    06/25/2014
!> \details
!>  This routine interpolates ice sheet thickness data to the MPAS mesh.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_land_ice_topography(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, forcingPool, landIceInitPool, diagnosticsPool, &
                                         statePool

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real (kind=RKIND), dimension(:), pointer :: landIceThkObserved, landIceDraftObserved, &
                                                   landIceFracObserved, landIceGroundedFracObserved

       real (kind=RKIND), dimension(:), pointer :: landIcePressure, landIceFraction, ssh, &
                                                   bottomDepth

       integer, pointer :: nCells
       integer, dimension(:), pointer :: maxLevelCell, landIceMask, modifySSHMask
       real (kind=RKIND), pointer :: config_land_ice_flux_rho_ice
       character (len=StrKIND), pointer :: config_global_ocean_topography_method, config_iterative_init_variable

       integer :: iCell

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_topography_method', &
                                 config_global_ocean_topography_method)
       call mpas_pool_get_config(domain % configs, 'config_land_ice_flux_rho_ice', config_land_ice_flux_rho_ice)
       call mpas_pool_get_config(domain % configs, 'config_iterative_init_variable', config_iterative_init_variable)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))

          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'landIceInit', landIceInitPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
          call mpas_pool_get_array(landIceInitPool, 'landIceDraftObserved', landIceDraftObserved)
          call mpas_pool_get_array(landIceInitPool, 'landIceThkObserved', landIceThkObserved)
          call mpas_pool_get_array(landIceInitPool, 'landIceFracObserved', landIceFracObserved)
          call mpas_pool_get_array(landIceInitPool, 'landIceGroundedFracObserved', landIceGroundedFracObserved)
          call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
          call mpas_pool_get_array(diagnosticsPool, 'modifySSHMask', modifySSHMask)
          call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
          call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
          call mpas_pool_get_array(forcingPool, 'landIceFraction', landIceFraction)

          if (config_global_ocean_topography_method .eq. "nearest_neighbor") then

             call ocn_init_interpolation_nearest_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                       landIceThkIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                       lonCell, latCell, landIceThkObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_nearest_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                       landIceDraftIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                       lonCell, latCell, landIceDraftObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_nearest_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                       landIceFracIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                       lonCell, latCell, landIceFracObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_nearest_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                       groundedFracIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                       lonCell, latCell, landIceGroundedFracObserved, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

          elseif (config_global_ocean_topography_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                        landIceThkIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                        lonCell, latCell, landIceThkObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_bilinear_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                        landIceDraftIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                        lonCell, latCell, landIceDraftObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_bilinear_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                        landIceFracIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                        lonCell, latCell, landIceFracObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_bilinear_horiz(landIceThkLon % array, landIceThkLat % array, &
                                                        groundedFracIC % array, nLonLandIceThk, nLatLandIceThk, &
                                                        lonCell, latCell, landIceGroundedFracObserved, nCells, &
                                                        inXPeriod = 2.0_RKIND * pii)

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_topography_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          ssh(:) = 0.0_RKIND
          landIceFraction(:) = 0.0_RKIND
          landIceMask(:) = 0
          modifySSHMask(:) = 0
          landIcePressure(:) = 0.0_RKIND
          do iCell = 1, nCells

             if(landIceFracObserved(iCell) > 0.5_RKIND) then
                landIceFraction(iCell) = landIceFracObserved(iCell)
                landIceMask(iCell) = 1
             end if

             ! nothing to do here if the cell is land
             if (maxLevelCell(iCell) <= 0) cycle

             if(config_iterative_init_variable == 'ssh') then
                ! we compute the land-ice pressure first and find out the SSH
                landIcePressure(iCell) = max(0.0_RKIND, config_land_ice_flux_rho_ice &
                   * gravity * landIceThkObserved(iCell))
                if(landIcePressure(iCell) > 0.0_RKIND) then
                   modifySSHMask(iCell) = 1
                end if
             else if(config_iterative_init_variable == 'landIcePressure' &
                .or. config_iterative_init_variable == 'landIcePressure_from_top_density') then
                ! we compute the SSH first and find out the land-ice pressure
                ssh(iCell) = min(0.0_RKIND,landIceDraftObserved(iCell))
                if(ssh(iCell) < 0.0_RKIND) then
                   modifySSHMask(iCell) = 1
                end if
             else
                write(stderrUnit,*) 'ERROR: Invalid choice of config_iterative_init_variable.'
                iErr = 1
                call mpas_dmpar_finalize(domain % dminfo)
             end if
          end do

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_interpolate_land_ice_topography!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_cull_inland_seas
!
!> \brief   Read the topography IC file
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine removes all inland seas. These are defined as isolated ocean cells.
!>   It uses a parallel version of an advancing front algorithm which might not be
!>   optimal for this purpose.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_cull_inland_seas(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: scratchPool, meshPool

       type (field1DInteger), pointer :: cullStackField, touchedCellField, oceanCellField

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, bottomDepth
       integer, dimension(:), pointer :: stack, oceanMask, touchMask
       integer, pointer :: stackSize

       integer :: iCell
       integer :: localStackSize, globalStackSize
       integer :: j, coc
       integer :: touched

       integer, pointer :: nCells, nCellsSolve, nVertLevels
       integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
       integer, dimension(:, :), pointer :: cellsOnCell

       iErr = 0

       call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

       call mpas_pool_get_field(scratchPool, 'cullStack', cullStackField)
       call mpas_pool_get_field(scratchPool, 'touchedCell', touchedCellField)
       call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)

       call mpas_allocate_scratch_field(cullStackField, .false.)
       call mpas_allocate_scratch_field(touchedCellField, .false.)
       call mpas_allocate_scratch_field(oceanCellField, .false.)

       ! Seed all deepest points for advancing front algorithm
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(scratchPool, 'cullStack', stack)
          call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
          call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
          call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

          stack(:) = 0
          oceanMask(:) = 0
          touchMask(:) = 0
          stackSize = 0

          ! Add all cells that have maxLevelCell == nVertLevels to stack
          do iCell = 1, nCellsSolve
             if (maxLevelCell(iCell) == nVertLevels) then
                stackSize = stackSize + 1
                stack(stackSize) = iCell
                touchMask(iCell) = 1
                oceanMask(iCell) = 1
             end if
          end do

          block_ptr => block_ptr % next
       end do

       ! Advancing front algorithm continues until all stacks on all processes are empty.
       globalStackSize = 1
       do while(globalStackSize /= 0)
          ! Advance front on each block with a non-zero stack until stack is empty.
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
             call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
             call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
             call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

             call mpas_pool_get_array(scratchPool, 'cullStack', stack)
             call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
             call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
             call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

             touched = 0
             do while(stackSize > 0)
                iCell = stack(stackSize)
                stackSize = stackSize - 1
                do j = 1, nEdgesOnCell(iCell)
                   coc = cellsOnCell(j, iCell)
                   if (touchMask(coc) == 0 .and. bottomDepth(coc) > 0.0_RKIND) then
                      oceanMask(coc) = 1
                      stackSize = stackSize + 1
                      stack(stackSize) = coc
                   end if
                   touchMask(coc) = 1
                   touched = touched + 1
                end do
             end do

             block_ptr => block_ptr % next
          end do

          ! Perform a halo exchange on oceanMask
          call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
          call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)
          call mpas_dmpar_exch_halo_field(oceanCellField)

          ! Check to see if any cells have been masked as ocean in the halo that have not been touched.
          ! If there are any, add them to the stack. Also, compute globalStackSize
          localStackSize = 0
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
             call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
             call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

             call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
             call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

             call mpas_pool_get_array(scratchPool, 'cullStack', stack)
             call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)
             call mpas_pool_get_array(scratchPool, 'touchedCell', touchMask)
             call mpas_pool_get_array(scratchPool, 'cullStackSize', stackSize)

             do iCell = nCellsSolve, nCells
                if (oceanMask(iCell) == 1 .and. touchMask(iCell) == 0) then
                   stackSize = stackSize + 1
                   stack(stackSize) = iCell
                   touchMask(iCell) = 1
                end if
             end do

             localStackSize = localStackSize + stackSize
             block_ptr => block_ptr % next
          end do

          call mpas_dmpar_sum_int(domain % dminfo, localStackSize, globalStackSize)
       end do

       ! Mark all cells that aren't ocean cells for removal
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(scratchPool, 'oceanCell', oceanMask)

          do iCell = 1, nCellsSolve
             if (oceanMask(iCell) == 0) then
                maxLevelCell(iCell) = -1
             end if
          end do
          block_ptr => block_ptr % next
       end do

       call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

       call mpas_pool_get_field(scratchPool, 'cullStack', cullStackField)
       call mpas_pool_get_field(scratchPool, 'touchedCell', touchedCellField)
       call mpas_pool_get_field(scratchPool, 'oceanCell', oceanCellField)

       call mpas_deallocate_scratch_field(cullStackField, .false.)
       call mpas_deallocate_scratch_field(touchedCellField, .false.)
       call mpas_deallocate_scratch_field(oceanCellField, .false.)

    end subroutine ocn_init_setup_global_ocean_cull_inland_seas!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_depth_levels
!
!> \brief   Read depth levels for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the depth levels from the temperature IC file and sets
!>  refBottomDepth accordingly
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_depth_levels(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (MPAS_Stream_type) :: depthStream

       type (mpas_pool_type), pointer :: meshPool

       character (len=StrKIND), pointer :: config_global_ocean_depth_file, config_global_ocean_depth_varname, &
                                           config_global_ocean_depth_dimname

       integer :: k

       real (kind=RKIND), pointer :: config_global_ocean_depth_conversion_factor

       real (kind=RKIND), dimension(:), pointer :: refBottomDepth

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_file', config_global_ocean_depth_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_varname', config_global_ocean_depth_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_conversion_factor', &
                                 config_global_ocean_depth_conversion_factor)

       ! Define stream for depth levels
       call MPAS_createStream(depthStream, domain % iocontext, config_global_ocean_depth_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup depth field for stream to be read in
       depthOutput % fieldName = trim(config_global_ocean_depth_varname)
       depthOutput % dimSizes(1) = nDepthOutput
       depthOutput % dimNames(1) = trim(config_global_ocean_depth_dimname)
       depthOutput % isVarArray = .false.
       depthOutput % isPersistent = .true.
       depthOutput % isActive = .true.
       depthOutput % hasTimeDimension = .false.
       depthOutput % block => domain % blocklist
       allocate(depthOutput % attLists(1))
       allocate(depthOutput % array(nDepthOutput))

       ! Add depth field to stream
       call MPAS_streamAddField(depthStream, depthOutput, iErr)

       ! Read stream
       call MPAS_readStream(depthStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(depthStream)
       depthOutput % array(:) = depthOutput % array(:) * config_global_ocean_depth_conversion_factor

       ! Set refBottomDepth depending on depth levels. And convert appropriately
       block_ptr => domain % blocklist
       do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         ! depthOutput is the mid-depth of each layer.  Convert to bottom depth.
         refBottomDepth(1) = 2.0_RKIND * depthOutput % array(1)
         do k=2,nDepthOutput
            refBottomDepth(k) = refBottomDepth(k-1) + 2*(depthOutput % array(k) - refBottomDepth(k-1))
         enddo

         block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_read_depth_levels!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_tracer_lat_lon
!
!> \brief   Read Lat/Lon for tracers in global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the latitude and longitude coordinats for tracers from the temperature IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_tracer_lat_lon(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_temperature_file, &
            config_global_ocean_tracer_lat_varname, config_global_ocean_tracer_nlat_dimname, &
            config_global_ocean_tracer_lon_varname, config_global_ocean_tracer_nlon_dimname, &
            config_global_ocean_tracer_depth_varname, config_global_ocean_tracer_ndepth_dimname

       real (kind=RKIND), pointer :: config_global_ocean_tracer_depth_conversion_factor

       logical, pointer :: config_global_ocean_tracer_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_lat_varname', &
                                 config_global_ocean_tracer_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', &
                                 config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_lon_varname', &
                                 config_global_ocean_tracer_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', &
                                 config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_depth_varname', &
                                 config_global_ocean_tracer_depth_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_ndepth_dimname', &
                                 config_global_ocean_tracer_ndepth_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_latlon_degrees', &
                                 config_global_ocean_tracer_latlon_degrees)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_depth_conversion_factor', &
                                 config_global_ocean_tracer_depth_conversion_factor)

       ! Define stream for depth levels
       call MPAS_createStream(tracerStream, domain % iocontext, config_global_ocean_temperature_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup tracerLat and tracerLon fields for stream to be read in
       tracerLat % fieldName = trim(config_global_ocean_tracer_lat_varname)
       tracerLat % dimSizes(1) = nLatTracer
       tracerLat % dimNames(1) = trim(config_global_ocean_tracer_nlat_dimname)
       tracerLat % isVarArray = .false.
       tracerLat % isPersistent = .true.
       tracerLat % isActive = .true.
       tracerLat % hasTimeDimension = .false.
       tracerLat % block => domain % blocklist
       allocate(tracerLat % attLists(1))
       allocate(tracerLat % array(nLatTracer))

       tracerLon % fieldName = trim(config_global_ocean_tracer_lon_varname)
       tracerLon % dimSizes(1) = nLonTracer
       tracerLon % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       tracerLon % isVarArray = .false.
       tracerLon % isPersistent = .true.
       tracerLon % isActive = .true.
       tracerLon % hasTimeDimension = .false.
       tracerLon % block => domain % blocklist
       allocate(tracerLon % attLists(1))
       allocate(tracerLon % array(nLonTracer))

       tracerDepth % fieldName = trim(config_global_ocean_tracer_depth_varname)
       tracerDepth % dimSizes(1) = nDepthTracer
       tracerDepth % dimNames(1) = trim(config_global_ocean_tracer_ndepth_dimname)
       tracerDepth % isVarArray = .false.
       tracerDepth % isPersistent = .true.
       tracerDepth % isActive = .true.
       tracerDepth % hasTimeDimension = .false.
       tracerDepth % block => domain % blocklist
       allocate(tracerDepth % attLists(1))
       allocate(tracerDepth % array(nDepthTracer))

       ! Add tracerLat and tracerLon fields to stream
       call MPAS_streamAddField(tracerStream, tracerLat, iErr)
       call MPAS_streamAddField(tracerStream, tracerLon, iErr)
       call MPAS_streamAddField(tracerStream, tracerDepth, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

       ! note IC tracer depth is in z coordinates, so negative
       tracerDepth % array(:) = - tracerDepth % array(:) * config_global_ocean_tracer_depth_conversion_factor

       if (config_global_ocean_tracer_latlon_degrees) then
          do iLat = 1, nLatTracer
             tracerLat % array(iLat) = tracerLat % array(iLat) * pii / 180.0_RKIND
          end do

          do iLon = 1, nLonTracer
             tracerLon % array(iLon) = tracerLon % array(iLon) * pii / 180.0_RKIND
          end do
       end if

    end subroutine ocn_init_setup_global_ocean_read_tracer_lat_lon!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_swData_lat_lon
!
!> \brief   Read Lat/Lon for swData in global ocean test case
!> \author  Luke Van Roekel
!> \date    11/16/2015
!> \details
!>  This routine reads the latitude and longitude coordinats for swData (chlorophyll,
!           clearSkyRadiation, zenithangle) from the swData IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_swData_lat_lon(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: SWStream

       character (len=StrKIND), pointer :: config_global_ocean_swData_file, config_global_ocean_swData_lat_varname, &
                                           config_global_ocean_swData_nlat_dimname, config_global_ocean_swData_lon_varname, &
                                           config_global_ocean_swData_nlon_dimname

       logical, pointer :: config_global_ocean_swData_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_lat_varname', &
                                 config_global_ocean_swData_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlat_dimname', &
                                 config_global_ocean_swData_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_lon_varname', &
                                 config_global_ocean_swData_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlon_dimname', &
                                 config_global_ocean_swData_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_latlon_degrees', &
                                 config_global_ocean_swData_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(SWStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup tracerLat and tracerLon fields for stream to be read in
       swDataLat % fieldName = trim(config_global_ocean_swData_lat_varname)
       swDataLat % dimSizes(1) = nLatSW
       swDataLat % dimNames(1) = trim(config_global_ocean_swData_nlat_dimname)
       swDataLat % isVarArray = .false.
       swDataLat % isPersistent = .true.
       swDataLat % isActive = .true.
       swDataLat % hasTimeDimension = .false.
       swDataLat % block => domain % blocklist
       allocate(swDataLat % attLists(1))
       allocate(swDataLat % array(nLatSW))

       swDataLon % fieldName = trim(config_global_ocean_swData_lon_varname)
       swDataLon % dimSizes(1) = nLonSW
       swDataLon % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
       swDataLon % isVarArray = .false.
       swDataLon % isPersistent = .true.
       swDataLon % isActive = .true.
       swDataLon % hasTimeDimension = .false.
       swDataLon % block => domain % blocklist
       allocate(swDataLon % attLists(1))
       allocate(swDataLon % array(nLonSW))

       ! Add tracerLat and tracerLon fields to stream
       call MPAS_streamAddField(SWStream, swDataLat, iErr)
       call MPAS_streamAddField(SWStream, swDataLon, iErr)

       ! Read stream
       call MPAS_readStream(SWStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(SWStream)

       if (config_global_ocean_swData_latlon_degrees) then
          do iLat = 1, nLatSW
             swDataLat % array(iLat) = swDataLat % array(iLat) * pii / 180.0_RKIND
          end do

          do iLon = 1, nLonSW
             swDataLon % array(iLon) = swDataLon % array(iLon) * pii / 180.0_RKIND
          end do
       end if

    end subroutine ocn_init_setup_global_ocean_read_swData_lat_lon!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_ecosys_lat_lon
!
!> \brief   Read Lat/Lon for ecosys in global ocean test case
!> \author  Mark Petersen
!> \date    Aug 19 2016
!> \details
!>  This routine reads the latitude and longitude coordinats for ecosys
!>  from the ecosys IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_ecosys_lat_lon(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_ecosys_file, &
            config_global_ocean_ecosys_lat_varname, config_global_ocean_ecosys_nlat_dimname, &
            config_global_ocean_ecosys_lon_varname, config_global_ocean_ecosys_nlon_dimname, &
            config_global_ocean_ecosys_depth_varname, config_global_ocean_ecosys_ndepth_dimname

       real (kind=RKIND), pointer :: config_global_ocean_ecosys_depth_conversion_factor

       logical, pointer :: config_global_ocean_ecosys_latlon_degrees

       integer :: iLat, iLon

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_file', config_global_ocean_ecosys_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_lat_varname', &
                                 config_global_ocean_ecosys_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlat_dimname', &
                                 config_global_ocean_ecosys_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_lon_varname', &
                                 config_global_ocean_ecosys_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlon_dimname', &
                                 config_global_ocean_ecosys_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_depth_varname', &
                                 config_global_ocean_ecosys_depth_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_ndepth_dimname', &
                                 config_global_ocean_ecosys_ndepth_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_latlon_degrees', &
                                 config_global_ocean_ecosys_latlon_degrees)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_depth_conversion_factor', &
                                 config_global_ocean_ecosys_depth_conversion_factor)

       ! Define stream for depth levels
       call MPAS_createStream(tracerStream, domain % iocontext, config_global_ocean_ecosys_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup tracerLat and tracerLon fields for stream to be read in
       tracerLat % fieldName = trim(config_global_ocean_ecosys_lat_varname)
       tracerLat % dimSizes(1) = nLatEcosys
       tracerLat % dimNames(1) = trim(config_global_ocean_ecosys_nlat_dimname)
       tracerLat % isVarArray = .false.
       tracerLat % isPersistent = .true.
       tracerLat % isActive = .true.
       tracerLat % hasTimeDimension = .false.
       tracerLat % block => domain % blocklist
       allocate(tracerLat % attLists(1))
       allocate(tracerLat % array(nLatEcosys))

       tracerLon % fieldName = trim(config_global_ocean_ecosys_lon_varname)
       tracerLon % dimSizes(1) = nLonEcosys
       tracerLon % dimNames(1) = trim(config_global_ocean_ecosys_nlon_dimname)
       tracerLon % isVarArray = .false.
       tracerLon % isPersistent = .true.
       tracerLon % isActive = .true.
       tracerLon % hasTimeDimension = .false.
       tracerLon % block => domain % blocklist
       allocate(tracerLon % attLists(1))
       allocate(tracerLon % array(nLonEcosys))

       tracerDepth % fieldName = trim(config_global_ocean_ecosys_depth_varname)
       tracerDepth % dimSizes(1) = nDepthEcosys
       tracerDepth % dimNames(1) = trim(config_global_ocean_ecosys_ndepth_dimname)
       tracerDepth % isVarArray = .false.
       tracerDepth % isPersistent = .true.
       tracerDepth % isActive = .true.
       tracerDepth % hasTimeDimension = .false.
       tracerDepth % block => domain % blocklist
       allocate(tracerDepth % attLists(1))
       allocate(tracerDepth % array(nDepthEcosys))

       ! Add ecosys Lat, Lon and Depth fields to stream
       call MPAS_streamAddField(tracerStream, tracerLat, iErr)
       call MPAS_streamAddField(tracerStream, tracerLon, iErr)
       call MPAS_streamAddField(tracerStream, tracerDepth, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

       ! note IC tracer depth is in z coordinates, so negative
       tracerDepth % array(:) = - tracerDepth % array(:) * config_global_ocean_ecosys_depth_conversion_factor

       if (config_global_ocean_ecosys_latlon_degrees) then
          do iLat = 1, nLatEcosys
             tracerLat % array(iLat) = tracerLat % array(iLat) * pii / 180.0_RKIND
          end do

          do iLon = 1, nLonEcosys
             tracerLon % array(iLon) = tracerLon % array(iLon) * pii / 180.0_RKIND
          end do
       end if

    end subroutine ocn_init_setup_global_ocean_read_ecosys_lat_lon!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_temperature
!
!> \brief   Read temperature ICs for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the temperature field from the temperature IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_temperature(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_temperature_file, config_global_ocean_temperature_varname, &
                                           config_global_ocean_tracer_nlon_dimname, config_global_ocean_tracer_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_file', config_global_ocean_temperature_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_temperature_varname', &
                                 config_global_ocean_temperature_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', &
                                 config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', &
                                 config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)

       ! Define stream for temperature IC
       call MPAS_createStream(tracerStream, domain % iocontext, config_global_ocean_temperature_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup temperature field for stream to be read in
       tracerIC % fieldName = trim(config_global_ocean_temperature_varname)
       tracerIC % dimSizes(1) = nLonTracer
       tracerIC % dimSizes(2) = nLatTracer
       tracerIC % dimSizes(3) = nDepthTracer
       tracerIC % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       tracerIC % dimNames(2) = trim(config_global_ocean_tracer_nlat_dimname)
       tracerIC % dimNames(3) = trim(config_global_ocean_depth_dimname)
       tracerIC % isVarArray = .false.
       tracerIC % isPersistent = .true.
       tracerIC % isActive = .true.
       tracerIC % hasTimeDimension = .false.
       tracerIC % block => domain % blocklist

       ! Add temperature field to stream
       call MPAS_streamAddField(tracerStream, tracerIC, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

    end subroutine ocn_init_setup_global_ocean_read_temperature!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_salinity
!
!> \brief   Read salinity ICs for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the salinity field from the salinity IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_salinity(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_salinity_file, config_global_ocean_salinity_varname, &
                                           config_global_ocean_tracer_nlon_dimname, config_global_ocean_tracer_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_salinity_file', config_global_ocean_salinity_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_salinity_varname', config_global_ocean_salinity_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlon_dimname', &
                                 config_global_ocean_tracer_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_nlat_dimname', &
                                 config_global_ocean_tracer_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', config_global_ocean_depth_dimname)

       ! Define stream for salinity IC
       call MPAS_createStream(tracerStream, domain % iocontext, config_global_ocean_salinity_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup salinity field for stream to be read in
       tracerIC % fieldName = trim(config_global_ocean_salinity_varname)
       tracerIC % dimSizes(1) = nLonTracer
       tracerIC % dimSizes(2) = nLatTracer
       tracerIC % dimSizes(3) = nDepthTracer
       tracerIC % dimNames(1) = trim(config_global_ocean_tracer_nlon_dimname)
       tracerIC % dimNames(2) = trim(config_global_ocean_tracer_nlat_dimname)
       tracerIC % dimNames(3) = trim(config_global_ocean_depth_dimname)
       tracerIC % isVarArray = .false.
       tracerIC % isPersistent = .true.
       tracerIC % isActive = .true.
       tracerIC % hasTimeDimension = .false.
       tracerIC % block => domain % blocklist
       allocate(tracerIC % attLists(1))
       allocate(tracerIC % array(nLonTracer, nLatTracer, nDepthTracer))

       ! Add salinity field to stream
       call MPAS_streamAddField(tracerStream, tracerIC, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

    end subroutine ocn_init_setup_global_ocean_read_salinity!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_ecosys
!
!> \brief   Read ecosys ICs for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the ecosys fields from the ecosys IC file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_ecosys(domain, fieldName, fileName, iErr)!{{{

       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       character (len=StrKIND), intent(in) :: fieldName, fileName

       type (MPAS_Stream_type) :: tracerStream

       character (len=StrKIND), pointer :: config_global_ocean_ecosys_nlon_dimname,  &
                                           config_global_ocean_ecosys_nlat_dimname, &
                                           config_global_ocean_depth_dimname

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlon_dimname', &
            config_global_ocean_ecosys_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlat_dimname', &
            config_global_ocean_ecosys_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_depth_dimname', &
            config_global_ocean_depth_dimname)

       ! Define stream for ecosys IC
       call MPAS_createStream(tracerStream, domain % iocontext, fileName, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup ecosys field for stream to be read in
       tracerIC % fieldName = trim(fieldName)
       tracerIC % dimSizes(1) = nLonEcosys
       tracerIC % dimSizes(2) = nLatEcosys
       tracerIC % dimSizes(3) = nDepthEcosys
       tracerIC % dimNames(1) = trim(config_global_ocean_ecosys_nlon_dimname)
       tracerIC % dimNames(2) = trim(config_global_ocean_ecosys_nlat_dimname)
       tracerIC % dimNames(3) = trim(config_global_ocean_depth_dimname)
       tracerIC % isVarArray = .false.
       tracerIC % isPersistent = .true.
       tracerIC % isActive = .true.
       tracerIC % hasTimeDimension = .false.
       tracerIC % block => domain % blocklist

       ! Add ecosys field to stream
       call MPAS_streamAddField(tracerStream, tracerIC, iErr)

       ! Read stream
       call MPAS_readStream(tracerStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(tracerStream)

    end subroutine ocn_init_setup_global_ocean_read_ecosys!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_ecosys_forcing
!
!> \brief   Read ecosys forcing for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine reads the ecosys forcing fields from the ecosys forcing file.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_ecosys_forcing(domain, fieldName, fileName, iErr)!{{{

       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       character (len=StrKIND), intent(in) :: fieldName, fileName

       type (MPAS_Stream_type) :: ecosysForcingStream

       character (len=StrKIND), pointer :: config_global_ocean_ecosys_nlon_dimname,  &
                                           config_global_ocean_ecosys_nlat_dimname, &
                                           config_global_ocean_ecosys_forcing_time_dimname

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlon_dimname', &
            config_global_ocean_ecosys_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_nlat_dimname', &
            config_global_ocean_ecosys_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_forcing_time_dimname', &
            config_global_ocean_ecosys_forcing_time_dimname)

       ! Define stream for ecosys forcing
       call MPAS_createStream(ecosysForcingStream, domain % iocontext, fileName, MPAS_IO_NETCDF, MPAS_IO_READ, ierr=iErr)

       ! Setup ecosys field for stream to be read in
       ecosysForcingIC % fieldName = trim(fieldName)
       ecosysForcingIC % dimSizes(1) = nLonEcosys
       ecosysForcingIC % dimSizes(2) = nLatEcosys
       nTimes = 1
       ecosysForcingIC % dimSizes(3) = nTimes
       ecosysForcingIC % dimNames(1) = trim(config_global_ocean_ecosys_nlon_dimname)
       ecosysForcingIC % dimNames(2) = trim(config_global_ocean_ecosys_nlat_dimname)
       ecosysForcingIC % dimNames(3) = trim(config_global_ocean_ecosys_forcing_time_dimname)
       ecosysForcingIC % isVarArray = .false.
       ecosysForcingIC % isPersistent = .true.
       ecosysForcingIC % isActive = .true.
       ecosysForcingIC % hasTimeDimension = .false.
       ecosysForcingIC % block => domain % blocklist

       ! Add ecosys field to stream
       call MPAS_streamAddField(ecosysForcingStream, ecosysForcingIC, iErr)

       ! Read stream
       call MPAS_readStream(ecosysForcingStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(ecosysForcingStream)

    end subroutine ocn_init_setup_global_ocean_read_ecosys_forcing!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_tracers
!
!> \brief   Interpolate tracer quantities to MPAS grid
!> \author  Mark Petersen, Doug Jacobsen, Xylar Asay-Davis
!> \date    08/23/2016
!> \details
!>  This routine interpolates the temperature/salinity data read in from the
!>  initial condition file to the MPAS grid.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_tracers(domain, tracerArray, tracerIndex, interpTracerName, iErr)!{{{

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), dimension(:, :, :), intent(inout) :: tracerArray
      integer, intent(in) :: tracerIndex
      integer, intent(out) :: iErr
      character (len=StrKIND), intent(in) :: interpTracerName

      type (block_type), pointer :: block_ptr
      type (mpas_pool_type), pointer :: meshPool, scratchPool, diagnosticsPool

      real (kind=RKIND) :: counter
      integer :: iSmooth, j, coc, iCell, k, nDepth
      integer, pointer :: nCells, nVertLevels

      integer, dimension(:), pointer :: maxLevelCell, nEdgesOnCell
      integer, dimension(:, :), pointer :: cellsOnCell

      real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
      real (kind=RKIND), dimension(:, :), pointer :: smoothedTracer

      character (len=StrKIND), pointer :: config_global_ocean_tracer_method
      integer, pointer ::  config_global_ocean_smooth_TS_iterations

      real (kind=RKIND), dimension(:), pointer :: outTracerColumn
      real (kind=RKIND), dimension(:,:), pointer :: zMid
      integer :: inKMax, outKMax

      type (field2DReal), pointer :: interpTracerField
      real (kind=RKIND), dimension(:,:), pointer :: interpTracer

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_global_ocean_tracer_method', &
           config_global_ocean_tracer_method)
      call mpas_pool_get_config(domain % configs, 'config_global_ocean_smooth_TS_iterations', &
           config_global_ocean_smooth_TS_iterations)

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

         call mpas_pool_get_array(scratchPool, trim(interpTracerName), interpTracer)

         if (config_global_ocean_tracer_method .eq. "nearest_neighbor") then
            call ocn_init_interpolation_nearest_horiz(tracerLon % array, tracerLat % array, &
                 tracerIC % array, nLonTracer, nLatTracer, &
                 lonCell, latCell, interpTracer, nCells, &
                 inXPeriod = 2.0_RKIND * pii)

         elseif (config_global_ocean_tracer_method .eq. "bilinear_interpolation") then

            call ocn_init_interpolation_bilinear_horiz(tracerLon % array, tracerLat % array, &
                 tracerIC % array, nLonTracer, nLatTracer, &
                 lonCell, latCell, interpTracer, nCells, &
                 inXPeriod = 2.0_RKIND * pii)

         else
            write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_tracer_method.'
            iErr = 1
            call mpas_dmpar_finalize(domain % dminfo)
         endif

         block_ptr => block_ptr % next
      end do

      ! Smooth the tracer
      if (config_global_ocean_smooth_TS_iterations .gt. 0) then
         call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)

         do iSmooth = 1,config_global_ocean_smooth_TS_iterations

            block_ptr => domain % blocklist
            do while(associated(block_ptr))
               call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
               call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

               call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

               call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
               call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

               call mpas_pool_get_array(scratchPool, trim(interpTracerName)//'Smooth', smoothedTracer)
               call mpas_pool_get_array(scratchPool, trim(interpTracerName), interpTracer)

               nDepth = size(smoothedTracer, 1)

               ! initialize smoothed tracer with original values
               smoothedTracer = interpTracer

               do iCell = 1, nCells

                  do k = 1, nDepth
                     ! Initialize counter to 1 because of central cell in stencil.
                     counter = 1

                     do j = 1, nEdgesOnCell(iCell)
                        coc = cellsOnCell(j, iCell)
                        if (coc<nCells+1) then
                           smoothedTracer(k, iCell) = smoothedTracer(k, iCell) + interpTracer (k, coc)
                           counter = counter + 1
                        end if

                     end do ! edgesOnCell

                     smoothedTracer(k, iCell) = smoothedTracer(k, iCell) / counter

                  end do ! k level

               end do ! iCell

               interpTracer = smoothedTracer

               block_ptr => block_ptr % next
            end do

            call mpas_pool_get_field(scratchPool, trim(interpTracerName), interpTracerField)
            call mpas_dmpar_exch_halo_field(interpTracerField)

         end do ! iSmooth

      endif

      ! reinterpolate tracers from original depths to zMid for PBCs or other modified vertical coordinates

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
         call mpas_pool_get_array(scratchPool, trim(interpTracerName), interpTracer)

         allocate(outTracerColumn(nVertLevels))

         inKMax = size(interpTracer, 1)
         do iCell = 1, nCells
            outKMax = maxLevelCell(iCell)
            if(outKMax < 1) cycle ! nothing to interpolate

            outTracerColumn(:) = 9.969209968386869e+36_RKIND
            call ocn_init_interpolation_linear_vert(tracerDepth % array(1:inKMax), &
                 interpTracer(1:inKMax,iCell), &
                 inKMax, &
                 zMid(1:outKMax,iCell), &
                 outTracerColumn(1:outKMax), &
                 outKMax, &
                 extrapolate=.false.)
            tracerArray(tracerIndex,1:outKMax,iCell) = outTracerColumn(1:outKMax)

         end do

         deallocate(outTracerColumn)

         block_ptr => block_ptr % next
      end do

    end subroutine ocn_init_setup_global_ocean_interpolate_tracers!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_ecosys_forcing
!
!> \brief   Interpolate ecosys forcing quantities to MPAS grid
!> \author  Doug Jacobsen
!> \date    03/05/2014
!> \details
!>  This routine interpolates the ecosys forcing data read in from the
!>  forcing file to the MPAS grid.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_ecosys_forcing(domain, fieldName, poolName, iErr)!{{{

       type (domain_type), intent(inout) :: domain
       character (len=StrKIND), intent(in) :: fieldName, poolName
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr
       type (mpas_pool_type), pointer :: meshPool, forcingPool, ecosysAuxiliary

       integer :: timeCounter
       integer, pointer :: nCells, nCellsSolve

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
       real (kind=RKIND), dimension(:), pointer :: ecosysForcingField

       character (len=StrKIND), pointer :: config_global_ocean_ecosys_method
       integer, pointer ::  config_global_ocean_smooth_ecosys_iterations

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_ecosys_method', &
            config_global_ocean_ecosys_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_smooth_ecosys_iterations', &
            config_global_ocean_smooth_ecosys_iterations)

       timeCounter = 1

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
          call mpas_pool_get_subpool(forcingPool, 'ecosysAuxiliary', ecosysAuxiliary)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

          if (trim(poolName) == 'ecosysAuxiliary') then
             call mpas_pool_get_array(ecosysAuxiliary, trim(fieldName), ecosysForcingField, 1)
          else if (trim(poolName) == 'forcing') then
             call mpas_pool_get_array(forcingPool, trim(fieldName), ecosysForcingField, 1)
          end if

          if (config_global_ocean_ecosys_method .eq. "nearest_neighbor") then
             call ocn_init_interpolation_nearest_horiz(tracerLon % array, tracerLat % array, &
                  ecosysForcingIC % array(:,:,timeCounter), nLonEcosys, nLatEcosys, &
                  lonCell, latCell, ecosysForcingField, nCells, &
                  inXPeriod = 2.0_RKIND * pii)

          elseif (config_global_ocean_ecosys_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(tracerLon % array, tracerLat % array, &
                  ecosysForcingIC % array(:,:,timeCounter), nLonEcosys, nLatEcosys, &
                  lonCell, latCell, ecosysForcingField, nCells, &
                  inXPeriod = 2.0_RKIND * pii)
          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_ecosys_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_interpolate_ecosys_forcing!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_restoring
!
!> \brief   Copy tracer quantities for restoring
!> \author  Doug Jacobsen, Xylar Asay-Davis
!> \date    03/05/2014
!> \details
!>  This routine copies temperature/salinity into surface and interior restoring fields.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_restoring(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr
       type (mpas_pool_type), pointer :: meshPool, statePool, tracersPool, forcingPool
       type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool, tracersInteriorRestoringFieldsPool

       real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
       real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
       real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate

       real (kind=RKIND), pointer :: config_global_ocean_piston_velocity
       real (kind=RKIND), pointer :: config_global_ocean_interior_restore_rate

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_piston_velocity', config_global_ocean_piston_velocity)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_interior_restore_rate', &
          config_global_ocean_interior_restore_rate)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
          call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

          call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
          call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
          call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
                                   activeTracersSurfaceRestoringValue, 1)
          call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', &
                                   activeTracersPistonVelocity, 1)
          call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue', &
                                   activeTracersInteriorRestoringValue, 1)
          call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringRate', &
                                   activeTracersInteriorRestoringRate, 1)

          ! set interior restoring values and rate
          if ( associated(activeTracersInteriorRestoringValue) .and. associated(activeTracers) ) then
             activeTracersInteriorRestoringValue(:, :, :) = activeTracers(:, :, :)
          end if

          if ( associated(activeTracersInteriorRestoringRate) ) then
             activeTracersInteriorRestoringRate(:, :, :) = config_global_ocean_interior_restore_rate
          end if

           ! set surface restoring values and rate
          if ( associated(activeTracersSurfaceRestoringValue) .and. associated(activeTracers) ) then
             activeTracersSurfaceRestoringValue(:, :) = activeTracers(:, 1, :)
          end if

          if ( associated(activeTracersPistonVelocity) ) then
             activeTracersPistonVelocity(:, :) = config_global_ocean_piston_velocity
          end if

          block_ptr => block_ptr % next
       end do


    end subroutine ocn_init_setup_global_ocean_interpolate_restoring!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_swData
!
!> \brief   Interpolate penetrating shortwave radiation quantities to MPAS grid
!> \author  Luke Van Roekel, Xylar Asay-Davis
!> \date    11/11/2015
!> \details
!>  This routine interpolates the penetrating swData data read in from the
!>  initial condition file to the MPAS grid.
!
!-----------------------------------------------------------------------

subroutine ocn_init_setup_global_ocean_interpolate_swData(domain, iErr)!{{{
   type (domain_type), intent(inout) :: domain
   integer, intent(out) :: iErr

   type (block_type), pointer :: block_ptr
   type (MPAS_Stream_type) :: zenithStream, chlorophyllStream, clearSkyStream
   type (mpas_pool_type), pointer :: meshPool, statePool, shortwavePool, diagnosticsPool

   type(MPAS_TimeInterval_type) :: timeStep ! time step interval
   type(MPAS_Time_Type) :: currentTime
   character(len=STRKIND) :: currentTimeStamp

   integer :: monIndex
   integer, pointer :: nCells, nCellsSolve, maxLevelCell

   integer, dimension(12), parameter :: daysInMonth     = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   real (kind=RKIND), dimension(:), pointer :: chlorophyllData, zenithAngle, clearSkyRadiation
   real (kind=RKIND), dimension(:), pointer :: latCell, lonCell

   character (len=StrKIND), pointer :: config_global_ocean_swData_method, xtime
   character (len=StrKIND), pointer :: config_global_ocean_swData_file, config_global_ocean_zenithAngle_varname, &
                                       config_global_ocean_swData_nlon_dimname, config_global_ocean_swData_nlat_dimname, &
                                       config_global_ocean_chlorophyll_varname, config_global_ocean_clearSky_varname

   iErr = 0

   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_file', config_global_ocean_swData_file)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_zenithAngle_varname', config_global_ocean_zenithAngle_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_chlorophyll_varname', config_global_ocean_chlorophyll_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_clearSky_varname', config_global_ocean_clearSky_varname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlon_dimname', config_global_ocean_swData_nlon_dimname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_nlat_dimname', config_global_ocean_swData_nlat_dimname)
   call mpas_pool_get_config(domain % configs, 'config_global_ocean_swData_method', config_global_ocean_swData_method)
   ! Define stream for zenithAngle IC
   call MPAS_createStream(zenithStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, &
                          ierr=iErr)
   call MPAS_createStream(chlorophyllStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, &
                          ierr=iErr)
   call MPAS_createStream(clearSkyStream, domain % iocontext, config_global_ocean_swData_file, MPAS_IO_NETCDF, MPAS_IO_READ, &
                          ierr=iErr)


   ! Setup zenithAngle field for stream to be read in

   zenithAngleIC % fieldName = trim(config_global_ocean_zenithAngle_varname)
   zenithAngleIC % dimSizes(1) = nLonSW
   zenithAngleIC % dimSizes(2) = nLatSW
   zenithAngleIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   zenithAngleIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)

   zenithAngleIC % isVarArray = .false.
   zenithAngleIC % isPersistent = .true.
   zenithAngleIC % isActive = .true.
   zenithAngleIC % hasTimeDimension = .false.
   zenithAngleIC % block => domain % blocklist
   allocate(zenithAngleIC % attLists(1))
   allocate(zenithAngleIC % array(nLonSW, nLatSW))

   ! Setup zenithAngle field for stream to be read in
   chlorophyllIC % fieldName = trim(config_global_ocean_chlorophyll_varname)
   chlorophyllIC % dimSizes(1) = nLonSW
   chlorophyllIC % dimSizes(2) = nLatSW
   chlorophyllIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   chlorophyllIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)
   chlorophyllIC % isVarArray = .false.
   chlorophyllIC % isPersistent = .true.
   chlorophyllIC % isActive = .true.
   chlorophyllIC % hasTimeDimension = .false.
   chlorophyllIC % block => domain % blocklist
   allocate(chlorophyllIC % attLists(1))
   allocate(chlorophyllIC % array(nLonSW, nLatSW))

   ! Setup zenithAngle field for stream to be read in
   clearSKYIC % fieldName = trim(config_global_ocean_clearSky_varname)
   clearSKYIC % dimSizes(1) = nLonSW
   clearSKYIC % dimSizes(2) = nLatSW
   clearSKYIC % dimNames(1) = trim(config_global_ocean_swData_nlon_dimname)
   clearSKYIC % dimNames(2) = trim(config_global_ocean_swData_nlat_dimname)
   clearSKYIC % isVarArray = .false.
   clearSKYIC % isPersistent = .true.
   clearSKYIC % isActive = .true.
   clearSKYIC % hasTimeDimension = .false.
   clearSKYIC % block => domain % blocklist
   allocate(clearSKYIC % attLists(1))
   allocate(clearSKYIC % array(nLonSW, nLatSW))
   ! Add chlorophyll field to stream

   call MPAS_streamAddField(zenithStream, zenithAngleIC, iErr)
   call MPAS_streamAddField(chlorophyllStream, chlorophyllIC, iErr)
   call MPAS_streamAddField(clearSkyStream, clearSKYIC, iErr)

   do monIndex=1,12
   block_ptr => domain % blocklist
   do while(associated(block_ptr))
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'shortwave', shortwavePool)
      call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics', diagnosticsPool)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

      call mpas_pool_get_array(diagnosticsPool, 'xtime', xtime)
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

      call mpas_pool_get_array(shortwavePool, 'chlorophyllData', chlorophyllData)
      call mpas_pool_get_array(shortWavePool, 'zenithAngle', zenithAngle)
      call mpas_pool_get_array(shortWavePool, 'clearSkyRadiation', clearSkyRadiation)

   ! Read stream

      call MPAS_readStream(zenithStream, monIndex, iErr)
      call MPAS_readStream(chlorophyllStream, monIndex, iErr)
      call MPAS_readStream(clearSkyStream, monIndex, iErr)


      if (config_global_ocean_swData_method .eq. "nearest_neighbor") then
         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   chlorophyllIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, chlorophyllData, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   zenithAngleIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, zenithAngle, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_nearest_horiz(swDataLon % array, swDataLat % array, &
                                                   clearSKYIC % array, nLonSW, nLatSW, &
                                                   lonCell, latCell, clearSkyRadiation, nCells, &
                                                   inXPeriod = 2.0_RKIND * pii)

      elseif (config_global_ocean_swData_method .eq. "bilinear_interpolation") then
         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    chlorophyllIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, chlorophyllData, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    zenithAngleIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, zenithAngle, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)

         call ocn_init_interpolation_bilinear_horiz(swDataLon % array, swDataLat % array, &
                                                    clearSKYIC % array, nLonSW, nLatSW, &
                                                    lonCell, latCell, clearSkyRadiation, nCells, &
                                                    inXPeriod = 2.0_RKIND * pii)
     else
         write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_swData_method.'
         iErr = 1
         call mpas_dmpar_finalize(domain % dminfo)
      endif

      block_ptr => block_ptr % next
   end do !loop on blocks

  ! increment clock with month string

   currentTime = mpas_get_clock_time(domain % clock, MPAS_NOW, iErr)
   call mpas_get_time(currentTime, dateTimeString=currentTimeStamp)

   xtime=currentTimeStamp
   call mpas_stream_mgr_write(domain % streamManager, streamID='shortwave_forcing_data_init', &
                              forceWriteNow=.true., ierr=ierr)
   call mpas_set_timeInterval(timeStep, dt=real(daysInMonth(monIndex),RKIND)*86400.0_RKIND)
   call mpas_advance_clock(domain % clock, timeStep)

   enddo  !ends loop over months

   ! Close stream
   call MPAS_closeStream(zenithStream)
   call MPAS_closeStream(chlorophyllStream)
   call MPAS_closeStream(clearSkyStream)

   ! reset mpas clock for other streams and final write

   currentTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, iErr)
   call mpas_set_clock_time(domain%clock, currentTime  , MPAS_NOW,iErr)
   currentTime = mpas_get_clock_time(domain % clock, MPAS_NOW, iErr)
   call mpas_get_time(currentTime, dateTimeString=currentTimeStamp)

   xtime=currentTimeStamp
   call mpas_stream_mgr_reset_alarms(domain%streamManager, streamID='shortwave_forcing_data_init', &
                                     direction=MPAS_STREAM_OUTPUT, ierr=ierr)

end subroutine ocn_init_setup_global_ocean_interpolate_swData!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_read_windstress
!
!> \brief   Read the windstress IC file
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details
!>  This routine reads the windstress IC file, including latitude and longitude
!>   information for windstress data.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_read_windstress(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (MPAS_Stream_type) :: windstressStream

       character (len=StrKIND), pointer :: config_global_ocean_windstress_file, config_global_ocean_windstress_lat_varname, &
                                           config_global_ocean_windstress_nlat_dimname, &
                                           config_global_ocean_windstress_lon_varname, &
                                           config_global_ocean_windstress_nlon_dimname, &
                                           config_global_ocean_windstress_zonal_varname, &
                                           config_global_ocean_windstress_meridional_varname

       logical, pointer :: config_global_ocean_windstress_latlon_degrees

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_file', config_global_ocean_windstress_file)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_lat_varname', &
                                 config_global_ocean_windstress_lat_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_nlat_dimname', &
                                 config_global_ocean_windstress_nlat_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_lon_varname', &
                                 config_global_ocean_windstress_lon_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_nlon_dimname', &
                                 config_global_ocean_windstress_nlon_dimname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_zonal_varname', &
                                 config_global_ocean_windstress_zonal_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_meridional_varname', &
                                 config_global_ocean_windstress_meridional_varname)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_latlon_degrees', &
                                 config_global_ocean_windstress_latlon_degrees)

       ! Define stream for depth levels
       call MPAS_createStream(windstressStream, domain % iocontext, config_global_ocean_windstress_file, MPAS_IO_NETCDF, &
                              MPAS_IO_READ, ierr=iErr)

       ! Setup windLat, windLon, and windIC fields for stream to be read in
       windLat % fieldName = trim(config_global_ocean_windstress_lat_varname)
       windLat % dimSizes(1) = nLatWind
       windLat % dimNames(1) = trim(config_global_ocean_windstress_nlat_dimname)
       windLat % isVarArray = .false.
       windLat % isPersistent = .true.
       windLat % isActive = .true.
       windLat % hasTimeDimension = .false.
       windLat % block => domain % blocklist
       allocate(windLat % attLists(1))
       allocate(windLat % array(nLatWind))

       windLon % fieldName = trim(config_global_ocean_windstress_lon_varname)
       windLon % dimSizes(1) = nLonWind
       windLon % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       windLon % isVarArray = .false.
       windLon % isPersistent = .true.
       windLon % isActive = .true.
       windLon % hasTimeDimension = .false.
       windLon % block => domain % blocklist
       allocate(windLon % attLists(1))
       allocate(windLon % array(nLonWind))

       zonalWindIC % fieldName = trim(config_global_ocean_windstress_zonal_varname)
       zonalWindIC % dimSizes(1) = nLonWind
       zonalWindIC % dimSizes(2) = nLatWind
       zonalWindIC % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       zonalWindIC % dimNames(2) = trim(config_global_ocean_windstress_nlat_dimname)
       zonalWindIC % isVarArray = .false.
       zonalWindIC % isPersistent = .true.
       zonalWindIC % isActive = .true.
       zonalWindIC % hasTimeDimension = .false.
       zonalWindIC % block => domain % blocklist
       allocate(zonalWindIC % attLists(1))
       allocate(zonalWindIC % array(nLonWind, nLatWind))

       meridionalWindIC % fieldName = trim(config_global_ocean_windstress_meridional_varname)
       meridionalWindIC % dimSizes(1) = nLonWind
       meridionalWindIC % dimSizes(2) = nLatWind
       meridionalWindIC % dimNames(1) = trim(config_global_ocean_windstress_nlon_dimname)
       meridionalWindIC % dimNames(2) = trim(config_global_ocean_windstress_nlat_dimname)
       meridionalWindIC % isVarArray = .false.
       meridionalWindIC % isPersistent = .true.
       meridionalWindIC % isActive = .true.
       meridionalWindIC % hasTimeDimension = .false.
       meridionalWindIC % block => domain % blocklist
       allocate(meridionalWindIC % attLists(1))
       allocate(meridionalWindIC % array(nLonWind, nLatWind))

       ! Add windLat, windLon, and windIC fields to stream
       call MPAS_streamAddField(windstressStream, windLat, iErr)
       call MPAS_streamAddField(windstressStream, windLon, iErr)
       call MPAS_streamAddField(windstressStream, zonalWindIC, iErr)
       call MPAS_streamAddField(windstressStream, meridionalWindIC, iErr)

       ! Read stream
       call MPAS_readStream(windstressStream, 1, iErr)

       ! Close stream
       call MPAS_closeStream(windstressStream)

       if (config_global_ocean_windstress_latlon_degrees) then
          windLat % array(:) = windLat % array(:) * pii / 180.0_RKIND
          windLon % array(:) = windLon % array(:) * pii / 180.0_RKIND
       end if

    end subroutine ocn_init_setup_global_ocean_read_windstress!}}}

!***********************************************************************
!
!  routine ocn_init_setup_global_ocean_interpolate_windstress
!
!> \brief   Interpolate the windstress IC to MPAS mesh
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details
!>  This routine interpolates windstress data to the MPAS mesh.
!
!-----------------------------------------------------------------------

    subroutine ocn_init_setup_global_ocean_interpolate_windstress(domain, iErr)!{{{
       type (domain_type), intent(inout) :: domain
       integer, intent(out) :: iErr

       type (block_type), pointer :: block_ptr

       type (mpas_pool_type), pointer :: meshPool, forcingPool

       real (kind=RKIND), dimension(:), pointer :: latCell, lonCell, windStressZonal, windStressMeridional

       integer, pointer :: nCells

       logical, pointer :: config_use_bulk_wind_stress
       character (len=StrKIND), pointer :: config_global_ocean_windstress_method
       real (kind=RKIND), pointer :: config_global_ocean_windstress_conversion_factor

       iErr = 0

       call mpas_pool_get_config(domain % configs, 'config_use_bulk_wind_stress', config_use_bulk_wind_stress)
       if (.not.config_use_bulk_wind_stress) then
          write(stderrUnit,'(A)') ' WARNING: wind stress not initialized because config_use_bulk_wind_stress = .false.'
          return
       endif

       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_method', config_global_ocean_windstress_method)
       call mpas_pool_get_config(domain % configs, 'config_global_ocean_windstress_conversion_factor', &
                                 config_global_ocean_windstress_conversion_factor)

       block_ptr => domain % blocklist
       do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'latCell', latCell)
          call mpas_pool_get_array(meshPool, 'lonCell', lonCell)

          call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
          call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)

          if (config_global_ocean_windstress_method .eq. "nearest_neighbor") then
             call ocn_init_interpolation_nearest_horiz(windLon % array, windLat % array, &
                                                       zonalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressZonal, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

             call ocn_init_interpolation_nearest_horiz(windLon % array, windLat % array, &
                                                       meridionalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressMeridional, nCells, &
                                                       inXPeriod = 2.0_RKIND * pii)

          elseif (config_global_ocean_windstress_method .eq. "bilinear_interpolation") then
             call ocn_init_interpolation_bilinear_horiz(windLon % array, windLat % array, &
                                                       zonalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressZonal, nCells, &
                                                       inXPeriod = 2.0_RKIND*pii)

             call ocn_init_interpolation_bilinear_horiz(windLon % array, windLat % array, &
                                                       meridionalWindIC % array, nLonWind, nLatWind, &
                                                       lonCell, latCell, windStressMeridional, nCells, &
                                                       inXPeriod = 2.0_RKIND*pii)

          else
             write(stderrUnit,*) 'ERROR: Invalid choice of config_global_ocean_windstress_method.'
             iErr = 1
             call mpas_dmpar_finalize(domain % dminfo)
          endif

          windStressZonal(:) = windStressZonal(:) * config_global_ocean_windstress_conversion_factor
          windStressMeridional(:) = windStressMeridional(:) * config_global_ocean_windstress_conversion_factor

          block_ptr => block_ptr % next
       end do

    end subroutine ocn_init_setup_global_ocean_interpolate_windstress!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_tracer_fields
!
!> \brief   Tracer field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine destroys the fields that were created to hold tracer
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_tracer_fields()!{{{
        deallocate(tracerLat % array)
        deallocate(tracerLon % array)
    end subroutine ocn_init_global_ocean_destroy_tracer_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_topo_fields
!
!> \brief   Topography field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details
!>  This routine destroys the fields that were created to hold topography
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_topo_fields()!{{{
        deallocate(topoIC % array)
        if(associated(oceanFracIC % array)) deallocate(oceanFracIC % array)
        deallocate(topoLat % array)
        deallocate(topoLon % array)
    end subroutine ocn_init_global_ocean_destroy_topo_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_land_ice_topography_fields
!
!> \brief   Topography field cleanup routine
!> \author  Jeremy Fyke, Xylar Asay-Davis, Mark Petersen
!> \date    06/23/2015
!> \details
!>  This routine destroys the fields created to hold land ice topography
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_land_ice_topography_fields()!{{{
        deallocate(landIceThkIC % array)
        deallocate(landIceDraftIC % array)
        deallocate(landIceThkLat % array)
        deallocate(landIceThkLon % array)
    end subroutine ocn_init_global_ocean_destroy_land_ice_topography_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_windstress_fields
!
!> \brief   Windstress field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/07/2014
!> \details
!>  This routine destroys the fields that were created to hold windstress
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_windstress_fields()!{{{
        deallocate(zonalWindIC % array)
        deallocate(meridionalWindIC % array)
        deallocate(windLat % array)
        deallocate(windLon % array)
    end subroutine ocn_init_global_ocean_destroy_windstress_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_swData_fields
!
!> \brief   penetrating shortwave data fields cleanup routine
!> \author  Luke Van Roekel
!> \date    11/11/2015
!> \details
!>  This routine destroys the fields that were created to hold penetrating sw radiation data
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_swData_fields()!{{{
        deallocate(chlorophyllIC % array)
        deallocate(zenithAngleIC % array)
        deallocate(clearSKYIC % array)
    end subroutine ocn_init_global_ocean_destroy_swData_fields!}}}

!***********************************************************************
!
!  routine ocn_init_global_ocean_destroy_ecosys_fields
!
!> \brief   Ecosys field cleanup routine
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine destroys the fields that were created to hold ecosys
!>  initial condition information
!
!-----------------------------------------------------------------------

    subroutine ocn_init_global_ocean_destroy_ecosys_fields()!{{{
        deallocate(tracerLat % array)
        deallocate(tracerLon % array)
    end subroutine ocn_init_global_ocean_destroy_ecosys_fields!}}}

!***********************************************************************
!
!  routine ocn_init_validate_global_ocean
!
!> \brief   Validation for global ocean test case
!> \author  Doug Jacobsen
!> \date    03/04/2014
!> \details
!>  This routine validates the configuration options for the global ocean test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_global_ocean(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout), target :: iocontext
      integer, intent(out) :: iErr

      type (mpas_io_context_type), pointer :: iocontext_ptr
      type (MPAS_IO_Handle_type) :: inputFile
      character (len=StrKIND), pointer :: config_init_configuration, &
                                          config_global_ocean_depth_file, &
                                          config_global_ocean_depth_dimname, &
                                          config_global_ocean_temperature_file, &
                                          config_global_ocean_salinity_file, &
                                          config_global_ocean_tracer_nlat_dimname, &
                                          config_global_ocean_tracer_nlon_dimname, &
                                          config_global_ocean_tracer_ndepth_dimname, &
                                          config_global_ocean_topography_file, &
                                          config_global_ocean_topography_nlat_dimname, &
                                          config_global_ocean_topography_nlon_dimname, &
                                          config_global_ocean_windstress_file, &
                                          config_global_ocean_windstress_nlat_dimname, &
                                          config_global_ocean_windstress_nlon_dimname, &
                                          config_global_ocean_land_ice_topo_file, &
                                          config_global_ocean_land_ice_topo_nlat_dimname, &
                                          config_global_ocean_land_ice_topo_nlon_dimname, &
                                          config_global_ocean_swData_file, &
                                          config_global_ocean_swData_nlon_dimname, &
                                          config_global_ocean_swData_nlat_dimname, &
                                          config_global_ocean_ecosys_file, &
                                          config_global_ocean_ecosys_nlat_dimname, &
                                          config_global_ocean_ecosys_nlon_dimname, &
                                          config_global_ocean_ecosys_ndepth_dimname

      integer, pointer :: config_vert_levels, config_global_ocean_tracer_vert_levels, &
                                          config_global_ocean_ecosys_vert_levels
      logical, pointer :: config_use_ecosysTracers
      logical, pointer :: landIceInitActive, config_global_ocean_depress_by_land_ice
      logical, pointer :: criticalPassagesActive, config_global_ocean_deepen_critical_passages

      iocontext_ptr => iocontext
      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)

      if(config_init_configuration .ne. trim('global_ocean')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)

      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_file', &
                                config_global_ocean_depth_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depth_dimname', &
                                config_global_ocean_depth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_temperature_file', &
                                config_global_ocean_temperature_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_salinity_file', &
                                config_global_ocean_salinity_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlat_dimname', &
                                config_global_ocean_tracer_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_nlon_dimname', &
                                config_global_ocean_tracer_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_ndepth_dimname', &
                                config_global_ocean_tracer_ndepth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_vert_levels', &
                                config_global_ocean_tracer_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_file', &
                                config_global_ocean_topography_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlat_dimname', &
                                config_global_ocean_topography_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_topography_nlon_dimname', &
                                config_global_ocean_topography_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_depress_by_land_ice', &
                                config_global_ocean_depress_by_land_ice)
      call mpas_pool_get_config(configPool, 'config_use_ecosysTracers', &
                                config_use_ecosysTracers)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_file', &
                                config_global_ocean_land_ice_topo_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_nlat_dimname', &
                                config_global_ocean_land_ice_topo_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_land_ice_topo_nlon_dimname', &
                                config_global_ocean_land_ice_topo_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_deepen_critical_passages', &
                                config_global_ocean_deepen_critical_passages)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_file', &
                                config_global_ocean_windstress_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlat_dimname', &
                                config_global_ocean_windstress_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_windstress_nlon_dimname', &
                                config_global_ocean_windstress_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_file', &
                                config_global_ocean_swData_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlat_dimname', &
                                config_global_ocean_swData_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_swData_nlon_dimname', &
                                config_global_ocean_swData_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_file', &
                                config_global_ocean_ecosys_file)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_nlat_dimname', &
                                config_global_ocean_ecosys_nlat_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_nlon_dimname', &
                                config_global_ocean_ecosys_nlon_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_ndepth_dimname', &
                                config_global_ocean_ecosys_ndepth_dimname)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_vert_levels', &
                                config_global_ocean_ecosys_vert_levels)

      call mpas_pool_get_package(packagePool, 'landIceInitActive', landIceInitActive)
      if ( config_global_ocean_depress_by_land_ice) then
         landIceInitActive = .true.
      end if

      call mpas_pool_get_package(packagePool, 'criticalPassagesActive', criticalPassagesActive)
      if ( config_global_ocean_deepen_critical_passages) then
         criticalPassagesActive = .true.
      end if

      if (trim(config_global_ocean_depth_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_depth_file'
         iErr = 1
         return
      end if

      inputFile = MPAS_io_open(config_global_ocean_depth_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
      if (iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: could not open file ', trim(config_global_ocean_depth_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, config_global_ocean_depth_dimname, nDepthOutput, iErr)

      call MPAS_io_close(inputFile, iErr)

      if (trim(config_global_ocean_temperature_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_temperature_file'
         iErr = 1
         return
      end if

      if (trim(config_global_ocean_salinity_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_salinity_file'
         iErr = 1
         return
      end if

      inputFile = MPAS_io_open(config_global_ocean_temperature_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
      if (iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: could not open file ', trim(config_global_ocean_temperature_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlat_dimname, nLatTracer, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_nlon_dimname, nLonTracer, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_tracer_ndepth_dimname, nDepthTracer, iErr)

      call MPAS_io_close(inputFile, iErr)

      if (config_global_ocean_tracer_vert_levels <= 0 .and. nDepthTracer > 0) then
         config_global_ocean_tracer_vert_levels = nDepthTracer
      else if(config_global_ocean_tracer_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Value of config_global_ocean_tracer_vert_levels=-1, ',&
              'but nDepthTracer was not correctly read from input file.'
         iErr = 1
      end if

      if (trim(config_global_ocean_windstress_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_windstress_file'
         iErr = 1
         return
      end if

      inputFile = MPAS_io_open(config_global_ocean_swData_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

      call MPAS_io_inq_dim(inputFile, config_global_ocean_swData_nlat_dimname, nLatSW, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_swData_nlon_dimname, nLonSW, iErr)

      call MPAS_io_close(inputFile, iErr)

      if (trim(config_global_ocean_topography_file) == 'none') then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_topography_file'
         iErr = 1
         return
      end if

      inputFile = MPAS_io_open(config_global_ocean_topography_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
      if (iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: could not open file ', trim(config_global_ocean_topography_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, config_global_ocean_topography_nlat_dimname, nLatTopo, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_topography_nlon_dimname, nLonTopo, iErr)

      call MPAS_io_close(inputFile, iErr)

      inputFile = MPAS_io_open(config_global_ocean_windstress_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
      if (iErr .ne. 0) then
         write(stderrUnit,*) 'ERROR: could not open file ', trim(config_global_ocean_windstress_file)
         return
      end if

      call MPAS_io_inq_dim(inputFile, config_global_ocean_windstress_nlat_dimname, nLatWind, iErr)
      call MPAS_io_inq_dim(inputFile, config_global_ocean_windstress_nlon_dimname, nLonWind, iErr)

      call MPAS_io_close(inputFile, iErr)

      if (config_vert_levels <= 0 .and. nDepthOutput > 0) then
         config_vert_levels = nDepthOutput
      else if(config_vert_levels <= 0) then
         write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Not given a usable value for vertical levels.'
         iErr = 1
      end if

      if ( config_use_ecosysTracers ) then
         if (trim(config_global_ocean_ecosys_file) == 'none') then
            write(stderrUnit,*) &
              'ERROR: Validation failed for global ocean. Invalid filename for config_global_ocean_windstress_file'
            iErr = 1
            return
         end if

         inputFile = MPAS_io_open(config_global_ocean_ecosys_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)

         call MPAS_io_inq_dim(inputFile, config_global_ocean_ecosys_nlat_dimname, nLatEcosys, iErr)
         call MPAS_io_inq_dim(inputFile, config_global_ocean_ecosys_nlon_dimname, nLonEcosys, iErr)
         call MPAS_io_inq_dim(inputFile, config_global_ocean_ecosys_ndepth_dimname, nDepthEcosys, iErr)

         call MPAS_io_close(inputFile, iErr)

         if (config_global_ocean_ecosys_vert_levels <= 0 .and. nDepthEcosys > 0) then
            config_global_ocean_ecosys_vert_levels = nDepthEcosys
         else if(config_global_ocean_ecosys_vert_levels <= 0) then
            write(stderrUnit,*) 'ERROR: Validation failed for global ocean. Value of config_global_ocean_ecosys_vert_levels=-1, ',&
                 'but nDepthEcosys was not correctly read from input file.'
            iErr = 1
         end if

      end if

      if ( config_global_ocean_depress_by_land_ice) then
         if (trim(config_global_ocean_land_ice_topo_file) == 'none') then
            write(stderrUnit,*) 'ERROR: Validation failed for global ocean. ', &
               'Invalid filename for config_global_ocean_land_ice_topo_file'
            iErr = 1
            return
         end if

         inputFile = MPAS_io_open(config_global_ocean_land_ice_topo_file, MPAS_IO_READ, MPAS_IO_NETCDF, iocontext_ptr, ierr=iErr)
         if (iErr .ne. 0) then
            write(stderrUnit,*) 'ERROR: could not open file ', trim(config_global_ocean_land_ice_topo_file)
            return
         end if

         call MPAS_io_inq_dim(inputFile, config_global_ocean_land_ice_topo_nlat_dimname, nLatLandIceThk, iErr)
         call MPAS_io_inq_dim(inputFile, config_global_ocean_land_ice_topo_nlon_dimname, nLonLandIceThk, iErr)

         call MPAS_io_close(inputFile, iErr)
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_global_ocean!}}}

!***********************************************************************

end module ocn_init_global_ocean

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
