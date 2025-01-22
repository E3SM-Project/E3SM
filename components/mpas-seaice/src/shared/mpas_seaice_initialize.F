!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_initialize
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_initialize

  use mpas_dmpar
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_stream_manager
  use mpas_io_units
  use mpas_abort
  use mpas_log, only: mpas_log_write
  use mpas_timer, only: mpas_timer_start, mpas_timer_stop

  implicit none

  private
  save

  public :: &
       seaice_init, &
       seaice_init_post_clock_advance

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init(&
       domain, &
       clock, &
       dt)!{{{

    use seaice_mesh, only: seaice_init_mesh
    use seaice_advection, only: seaice_init_advection
    use seaice_velocity_solver, only: seaice_init_velocity_solver
    use seaice_icepack, only: &
         seaice_init_icepack_constants, &
         seaice_init_icepack_physics_package_parameters, &
         seaice_init_icepack_physics_package_variables
    use seaice_column, only: &
         seaice_init_column_constants, &
         seaice_init_column_physics_package_parameters, &
         seaice_init_column_physics_package_variables
    use seaice_forcing, only: seaice_forcing_init, seaice_reset_coupler_fluxes
    use seaice_diagnostics, only: &
         seaice_set_testing_system_test_arrays
    use seaice_mesh_pool, only: &
         seaice_mesh_pool_create
    use seaice_prescribed, only: &
         seaice_init_prescribed_ice
    use seaice_special_boundaries, only: &
         seaice_init_special_boundaries, &
         seaice_set_special_boundaries_tracers

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_column_physics_type

    call mpas_log_write("seaice_init: Initialize sea-ice model")

    ! check configs
    call seaice_check_configs_coupled(domain)

    call mpas_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)

    ! set column constants
    if (trim(config_column_physics_type) == "icepack") then
       call seaice_init_icepack_constants()
    else if (trim(config_column_physics_type) == "column_package") then
       call seaice_init_column_constants()
    else
       call MPAS_log_write("Unknown config_column_physics_type: "//trim(config_column_physics_type), MPAS_LOG_CRIT)
    endif ! config_column_physics_type

    ! init testing system test arrays
    call seaice_set_testing_system_test_arrays(domain)

    ! initialize junk values
    call init_junk_values(domain)

    ! initialize special boundaries
    call seaice_init_special_boundaries(domain)

    ! initialize prescribed ice
    call seaice_init_prescribed_ice(domain)

    ! initialize landice mask if needed for testing
    call init_test_ice_shelf_mask(domain)

    ! mesh initialization
    call mpas_log_write(" Initialize mesh...")
    call seaice_init_mesh(domain)

    call mpas_log_write(" Initialize mesh pool...")
    call seaice_mesh_pool_create(domain)

    ! init the basic column physics package
    call mpas_log_write(" Initialize column parameters...")
    if (trim(config_column_physics_type) == "icepack") then
       call seaice_init_icepack_physics_package_parameters(domain)
    else if (trim(config_column_physics_type) == "column_package") then
       call seaice_init_column_physics_package_parameters(domain)
    endif ! config_column_physics_type

    ! init coupler fluxes
    call mpas_log_write(" Initialize coupler fields...")
    call initialize_coupler_fields(domain)

    ! initialize forcing
    call mpas_log_write(" Initialize forcing...")
    call seaice_forcing_init(domain, clock)

    ! init dynamics
    call mpas_log_write(" Initialize velocity solver...")
    call mpas_timer_start("Velocity solver init")
    call seaice_init_velocity_solver(domain)
    call mpas_timer_stop("Velocity solver init")

    ! init advection
    call mpas_log_write(" Initialize advection...")
    call seaice_init_advection(domain)

    ! column physics initialization
    call mpas_log_write(" Initialize column variables...")
    if (trim(config_column_physics_type) == "icepack") then
       call seaice_init_icepack_physics_package_variables(domain, clock)
    else if (trim(config_column_physics_type) == "column_package") then
       call seaice_init_column_physics_package_variables(domain, clock)
    endif ! config_column_physics_type

    ! init ice state
    call mpas_log_write(" Initialize ice state...")
    call init_ice_state(domain, clock)

    ! initial halo exchange
    call mpas_log_write(" Initial halo exchanges...")
    call initial_halo_exchanges(domain)

    ! coupler flux reset
    call mpas_log_write(" Initialize coupler fluxes...")
    call seaice_reset_coupler_fluxes(domain)

    ! special boundaries tracers
    call seaice_set_special_boundaries_tracers(domain)

    ! check constants
    call seaice_check_constants_coupled

  end subroutine seaice_init!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_post_clock_advance
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_post_clock_advance(&
       domain, &
       clock)!{{{

    use seaice_icepack, only: &
         seaice_init_icepack_shortwave
    use seaice_column, only: &
         seaice_init_column_shortwave

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_shortwave, &
         config_do_restart

    character(len=strKIND), pointer :: &
         config_column_physics_type

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)
    call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)
    call MPAS_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)

    if (config_use_column_physics .and. config_use_column_shortwave .and. .not. config_do_restart) then
       if (trim(config_column_physics_type) == "icepack") then
          call seaice_init_icepack_shortwave(domain, clock)
       else if (trim(config_column_physics_type) == "column_package") then
          call seaice_init_column_shortwave(domain, clock)
       endif ! config_column_physics_type
    endif

  end subroutine seaice_init_post_clock_advance

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_junk_values
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_junk_values(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34_RKIND

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)

       call MPAS_pool_get_array(mesh, "areaCell", areaCell)

       njunk = nCells + 1

       areaCell(njunk) = value_junk

       block => block % next
      end do

  end subroutine init_junk_values!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state(&
       domain, &
       clock)!{{{

    use seaice_testing, only: &
         seaice_init_square_test_case_hex, &
         seaice_init_square_point_test_case_hex
    use seaice_icepack, only: &
         seaice_icepack_aggregate
    use seaice_column, only: &
         seaice_column_aggregate

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_clock_type), intent(in) :: &
         clock

    type(block_type), pointer :: &
         block

    character(len=strKIND), pointer :: &
         config_initial_condition_type, &
         config_initial_velocity_type, &
         config_column_physics_type

    integer, dimension(:), pointer :: &
         interiorVertex

    type(MPAS_pool_type), pointer :: &
         icestate, &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCellInitial, &
         iceAreaCell

    logical, pointer :: &
         config_do_restart, &
         config_use_column_physics

    call MPAS_pool_get_config(domain % configs, "config_do_restart", config_do_restart)

    if (.not. config_do_restart) then

       call MPAS_pool_get_config(domain % configs, "config_initial_condition_type", config_initial_condition_type)
       call MPAS_pool_get_config(domain % configs, "config_initial_velocity_type", config_initial_velocity_type)

       ! set volumes/areas
       block => domain % blocklist
       do while (associated(block))

          if (trim(config_initial_condition_type) == "uniform") then

             call init_ice_state_uniform_ice(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "circle") then

             call init_ice_state_circle_of_ice(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "uniform_interior") then

             call init_ice_state_uniform_interior(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "uniform_1D") then

             call init_ice_state_uniform_1D(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "special") then

             call init_ice_state_special(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "random_coverage") then

             call init_ice_state_random_coverage(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "square") then

             call seaice_init_square_test_case_hex(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "square_point") then

             call seaice_init_square_point_test_case_hex(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "cice_default") then

             call init_ice_cice_default(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "single_cell") then

             call init_ice_single_cell(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "ridging") then

             call init_ice_ridging(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_condition_type) == "no_ice") then

             ! no ice initial condition
             continue

          else if (trim(config_initial_condition_type) == "none") then

             ! for test cases
             continue

          else

             call mpas_log_write(&
                  "init_ice_state: config_initial_condition_type unknown:"//trim(config_initial_condition_type), &
                  MPAS_LOG_CRIT)

          endif

          block => block % next
       end do

       ! halo exchanges
       call MPAS_dmpar_field_halo_exch(domain, 'iceAreaCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'iceVolumeCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'snowVolumeCategory', timeLevel=1)
       call MPAS_dmpar_field_halo_exch(domain, 'surfaceTemperature', timeLevel=1)

       ! other initialize variables
       block => domain % blocklist
       do while (associated(block))
          call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
          call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
          call MPAS_pool_get_array(icestate, "iceAreaCellInitial", iceAreaCellInitial)
          call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
          iceAreaCellInitial = iceAreaCell
          block => block % next
       enddo

       ! set velocities
       block => domain % blocklist
       do while (associated(block))

          if (trim(config_initial_velocity_type) == "uniform") then

             call init_ice_velocity_uniform(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "ocean") then

             call init_ice_velocity_ocean(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "random") then

             call init_ice_velocity_random(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "special") then

             call init_ice_velocity_special(&
                  block, &
                  domain % configs)

          else if (trim(config_initial_velocity_type) == "none") then

             continue

          else

             call mpas_log_write(&
                  "init_ice_state: config_initial_velocity_type unknown: "//trim(config_initial_velocity_type), &
                  MPAS_LOG_CRIT)

          endif

          block => block % next
       end do

       ! halo exchanges
       call MPAS_dmpar_field_halo_exch(domain, 'uVelocity')
       call MPAS_dmpar_field_halo_exch(domain, 'vVelocity')

    endif

    ! aggregate tracers even if restart
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_physics", config_use_column_physics)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_column_physics_type", config_column_physics_type)
    if (config_use_column_physics) then
       if (trim(config_column_physics_type) == "icepack") then
          call seaice_icepack_aggregate(domain)
       else if (trim(config_column_physics_type) == "column_package") then
          call seaice_column_aggregate(domain)
       endif ! config_column_physics_type
    endif

  end subroutine init_ice_state!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_uniform_ice(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceDegreesToRadians

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate

    integer, pointer :: &
         nCellsSolve

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         surfaceTemperatureCell

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)

    call MPAS_pool_get_dimension(mesh, "nCells", nCellsSolve)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)

    do iCell = 1, nCellsSolve

       if (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
           latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) then

          ! has ice
          iceAreaCategory(1,:,iCell)    = config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
          surfaceTemperature(1,:,iCell) = -1.0_RKIND

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND
          surfaceTemperature(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

    ! integrated quantities
    do iCell = 1, nCellsSolve

       iceAreaCell(iCell) = sum(iceAreaCategory(1,:,iCell))
       surfaceTemperatureCell(iCell) = -20.15_RKIND

    enddo ! iCell

  end subroutine init_ice_state_uniform_ice!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_1D
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date November, 2022
!> \details
!> Initializes sea ice in the first category only with constant
!> values and latitude bounds defined in the namelist.
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_uniform_1D(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceDegreesToRadians

    use ice_colpkg, only: &
         colpkg_init_trcr, &
         colpkg_enthalpy_snow

    use icepack_intfc, only: &
         icepack_init_enthalpy, &
         icepack_enthalpy_snow, &
         icepack_warnings_aborted

    use seaice_icepack, only: &
         seaice_icepack_write_warnings

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nSnowLayers, &
         nIceLayers, &
         nCategories

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         airTemperature, &
         seaSurfaceTemperature, &
         seaFreezingTemperature

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature, &
         snowEnthalpy, &
         iceEnthalpy, &
         iceSalinity

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         surfaceTemperatureCell

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    character(len=strKIND), pointer :: &
         config_column_physics_type

    integer :: &
         iCell, &
         iSnowLayer, &
         iIceLayer, &
         iCategory

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_config(configs, "config_column_physics_type", config_column_physics_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "initial", initial)

    call MPAS_pool_get_dimension(mesh, "nCells", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)
    call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

    call MPAS_pool_get_array(tracers_aggregate, "iceAreaCell", iceAreaCell)
    call MPAS_pool_get_array(tracers_aggregate, "surfaceTemperatureCell", surfaceTemperatureCell)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
    call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

    do iCell = 1, nCellsSolve

       iceAreaCategory(1,:,iCell)    = 0.0_RKIND
       iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
       snowVolumeCategory(1,:,iCell) = 0.0_RKIND
       surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)
       iceEnthalpy(:,:,iCell)        = 0.0_RKIND

       if (trim(config_column_physics_type) == "icepack") then

          do iSnowLayer = 1, nSnowLayers
             snowEnthalpy(iSnowLayer,:,iCell) = icepack_enthalpy_snow(0.0_RKIND)
          end do

          if (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
              latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) then

             ! has ice
             iceAreaCategory(1,1,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,1,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,1,iCell) = config_initial_snow_volume
             surfaceTemperature(1,1,iCell) = -1.0_RKIND

             call icepack_init_enthalpy(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,1,iCell), &
                  iceEnthalpy(:,1,iCell), &
                  snowEnthalpy(:,1,iCell))
             call seaice_icepack_write_warnings(icepack_warnings_aborted())
          endif

       else if (trim(config_column_physics_type) == "column_package") then

          do iSnowLayer = 1, nSnowLayers
             snowEnthalpy(iSnowLayer,:,iCell) = colpkg_enthalpy_snow(0.0_RKIND)
          end do

          if (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
              latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) then

             ! has ice
             iceAreaCategory(1,1,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,1,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,1,iCell) = config_initial_snow_volume
             surfaceTemperature(1,1,iCell) = -1.0_RKIND

             call colpkg_init_trcr(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,1,iCell), &
                  nIceLayers, &
                  nSnowLayers, &
                  iceEnthalpy(:,1,iCell), &
                  snowEnthalpy(:,1,iCell))
          endif

       endif ! config_column_physics_type

       iceAreaCell(iCell) = sum(iceAreaCategory(1,:,iCell))
       surfaceTemperatureCell(iCell) = -20.15_RKIND

       do iCategory = 1, nCategories
          do iIceLayer = 1, nIceLayers
             iceSalinity(iIceLayer,iCategory,iCell) = initialSalinityProfile(iIceLayer,iCell)
          enddo ! iIceLayer
       enddo ! iCategory
    enddo ! iCell

  end subroutine init_ice_state_uniform_1D!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_cice_default
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date March 2nd 1015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_cice_default(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceDegreesToRadians

    use icepack_intfc, only: &
         icepack_init_enthalpy, &
         icepack_enthalpy_snow, &
         icepack_warnings_aborted

    use seaice_icepack, only: &
         seaice_icepack_write_warnings

    use seaice_column, only: &
         seaice_column_init_trcr, &
         seaice_column_enthalpy_snow

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         airTemperature, &
         seaSurfaceTemperature, &
         seaFreezingTemperature

    integer, dimension(:), pointer :: &
         landIceMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceSalinity

    real(kind=RKIND), pointer :: &
         config_initial_latitude_north, &
         config_initial_latitude_south

    real(kind=RKIND), dimension(:), allocatable :: &
         initialCategoryIceArea, &
         initialCategoryIceThickness

    integer :: &
         iCell, &
         iCategory, &
         iIceLayer, &
         iSnowLayer

    real(kind=RKIND), parameter :: &
         initialCategorySnowThickness = 0.2_RKIND

    character(len=strKIND), pointer :: &
         config_column_physics_type

    call MPAS_pool_get_config(configs, "config_column_physics_type", config_column_physics_type)

    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "initial", initial)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)
    call MPAS_pool_get_array(ocean_coupling, "landIceMask", landIceMask)

    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
    call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

    ! initial volumes and areas
    allocate(initialCategoryIceArea(nCategories))
    allocate(initialCategoryIceThickness(nCategories))

    call initial_category_areas_and_volumes(&
         block, &
         initialCategoryIceArea, &
         initialCategoryIceThickness)

    do iCell = 1, nCellsSolve

       if (seaSurfaceTemperature(iCell) <= seaFreezingTemperature(iCell) + 0.2_RKIND .and. &
          (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
           latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) .and. &
          landIceMask(iCell) == 0) then

          ! has ice
          if (trim(config_column_physics_type) == "icepack") then
             do iCategory = 1, nCategories

                iceAreaCategory(1,iCategory,iCell)    = initialCategoryIceArea(iCategory)
                iceVolumeCategory(1,iCategory,iCell)  = initialCategoryIceArea(iCategory) * initialCategoryIceThickness(iCategory)
                snowVolumeCategory(1,iCategory,iCell) = min(iceAreaCategory(1,iCategory,iCell) * initialCategorySnowThickness, &
                     0.2_RKIND * iceVolumeCategory(1,iCategory,iCell))

                call icepack_init_enthalpy(&
                     airTemperature(iCell), &
                     seaFreezingTemperature(iCell), &
                     initialSalinityProfile(:,iCell), &
                     initialMeltingTemperatureProfile(:,iCell), &
                     surfaceTemperature(1,iCategory,iCell), &
                     iceEnthalpy(:,iCategory,iCell), &
                     snowEnthalpy(:,iCategory,iCell))
                call seaice_icepack_write_warnings(icepack_warnings_aborted())

             enddo ! iCategory
          else if (trim(config_column_physics_type) == "column_package") then
             do iCategory = 1, nCategories

                iceAreaCategory(1,iCategory,iCell)    = initialCategoryIceArea(iCategory)
                iceVolumeCategory(1,iCategory,iCell)  = initialCategoryIceArea(iCategory) * initialCategoryIceThickness(iCategory)
                snowVolumeCategory(1,iCategory,iCell) = min(iceAreaCategory(1,iCategory,iCell) * initialCategorySnowThickness, &
                     0.2_RKIND * iceVolumeCategory(1,iCategory,iCell))

                call seaice_column_init_trcr(&
                     airTemperature(iCell), &
                     seaFreezingTemperature(iCell), &
                     initialSalinityProfile(:,iCell), &
                     initialMeltingTemperatureProfile(:,iCell), &
                     surfaceTemperature(1,iCategory,iCell), &
                     nIceLayers, &
                     nSnowLayers, &
                     iceEnthalpy(:,iCategory,iCell), &
                     snowEnthalpy(:,iCategory,iCell))

             enddo ! iCategory
          endif ! config_column_physics_type

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND

          surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)
          if (trim(config_column_physics_type) == "icepack") then
             do iSnowLayer = 1, nSnowLayers
                snowEnthalpy(iSnowLayer,:,iCell) = icepack_enthalpy_snow(0.0_RKIND)
             end do
          else if (trim(config_column_physics_type) == "column_package") then
             do iSnowLayer = 1, nSnowLayers
                snowEnthalpy(iSnowLayer,:,iCell) = seaice_column_enthalpy_snow(0.0_RKIND)
             end do
          endif ! config_column_physics_type

       endif

    enddo ! iCell

    ! clean up
    deallocate(initialCategoryIceArea)
    deallocate(initialCategoryIceThickness)

    ! get vertical salinity profiles
    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iIceLayer = 1, nIceLayers
             iceSalinity(iIceLayer,iCategory,iCell) = initialSalinityProfile(iIceLayer,iCell)
          enddo ! iIceLayer
       enddo ! iCategory
    enddo ! iCell

  end subroutine init_ice_cice_default!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initial_category_areas_and_volumes
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date March 2nd 1015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initial_category_areas_and_volumes(&
       block, &
       initialCategoryIceArea, &
       initialCategoryIceThickness)

    use seaice_constants, only: &
         seaicePuny

    use icepack_intfc, only: &
         icepack_init_itd, &
         icepack_warnings_aborted

    use seaice_icepack, only: &
         seaice_icepack_write_warnings

    use seaice_column, only: &
         seaice_column_init_itd

    ! Note: the resulting average ice thickness
    ! tends to be less than hbar due to the
    ! nonlinear distribution of ice thicknesses

    type(block_type), intent(in) :: block

    real(kind=RKIND), dimension(:), intent(out) :: &
         initialCategoryIceArea, &
         initialCategoryIceThickness

    type(MPAS_pool_type), pointer :: &
         mesh, &
         initial

    real(kind=RKIND), dimension(:), pointer :: &
         categoryThicknessLimits

    real(kind=RKIND) :: &
         areaCategorySum

    integer, pointer :: &
         nCategories

    integer :: &
         iCategory

    logical, pointer :: &
         config_use_column_physics

    real(kind=RKIND), parameter :: &
         thicknessWithLargestArea = 3.0_RKIND ! initial ice thickness with greatest area

    logical :: &
         abortFlag

    character(len=strKIND) :: &
         abortMessage

    character(len=strKIND), pointer :: &
         config_column_physics_type

    call MPAS_pool_get_config(block % configs, "config_column_physics_type", config_column_physics_type)

    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)

    !--------------------------------------------------------
    ! init caregory limits if needed
    !--------------------------------------------------------

    call MPAS_pool_get_config(block % configs, "config_use_column_physics", config_use_column_physics)

    if (.not. config_use_column_physics) then

       if (trim(config_column_physics_type) == "icepack") then
          call icepack_init_itd(&
               categoryThicknessLimits)
          call seaice_icepack_write_warnings(icepack_warnings_aborted())
       else if (trim(config_column_physics_type) == "column_package") then
          call seaice_column_init_itd(&
               nCategories, &
               categoryThicknessLimits, &
               abortFlag, &
               abortMessage)
          if (abortFlag) then
             call mpas_log_write(&
                  "initial_category_areas_and_volumes: "//trim(abortMessage), &
                  MPAS_LOG_CRIT)
          endif
       endif ! config_column_physics_type

    endif

    !--------------------------------------------------------
    ! volumes
    !--------------------------------------------------------

    ! loop over categories to set initial volumes
    do iCategory = 1, nCategories-1

       ! middle of the category limits
       initialCategoryIceThickness(iCategory) = &
            0.5_RKIND * (categoryThicknessLimits(iCategory) + categoryThicknessLimits(iCategory+1))

    enddo ! iCategory

    ! thicknest category
    initialCategoryIceThickness(nCategories) = &
         categoryThicknessLimits(nCategories) + 1.0_RKIND

    !--------------------------------------------------------
    ! areas
    !--------------------------------------------------------

    ! initialize sum of areas in categories
    areaCategorySum = 0.0_RKIND

    do iCategory = 1, nCategories

       ! parabola, max at h=thicknessWithLargestArea, zero at h=0, 2*thicknessWithLargestArea
       initialCategoryIceArea(iCategory) = &
            max(0.0_RKIND, &
            (2.0_RKIND * thicknessWithLargestArea * initialCategoryIceThickness(iCategory) - &
             initialCategoryIceThickness(iCategory)**2))

       areaCategorySum = areaCategorySum + initialCategoryIceArea(iCategory)

    enddo ! iCategory

    ! normalize
    do iCategory = 1, nCategories

       initialCategoryIceArea(iCategory) = initialCategoryIceArea(iCategory) / &
            (areaCategorySum + seaicePuny / nCategories)

    enddo ! iCategory

  end subroutine initial_category_areas_and_volumes

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_single_cell
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 12th July 2022
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_single_cell(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceDegreesToRadians

    use ice_colpkg, only: &
         colpkg_init_trcr, &
         colpkg_enthalpy_snow

    use icepack_intfc, only: &
         icepack_init_enthalpy, &
         icepack_enthalpy_snow, &
         icepack_warnings_aborted

    use seaice_icepack, only: &
         seaice_icepack_write_warnings

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    real(kind=RKIND), dimension(:), pointer :: &
         airTemperature, &
         seaFreezingTemperature

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceSalinity

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume

    character(len=strKIND), pointer :: &
         config_column_physics_type

    integer :: &
         iCell, &
         iCategory, &
         iIceLayer, &
         iSnowLayer

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_column_physics_type", config_column_physics_type)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "initial", initial)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)

    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
    call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

    do iCell = 1, nCellsSolve

       if (trim(config_column_physics_type) == "icepack") then

          ! has ice
          do iCategory = 1, nCategories

             iceAreaCategory(1,iCategory,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,iCategory,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,iCategory,iCell) = config_initial_snow_volume

             call icepack_init_enthalpy(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  iceEnthalpy(:,iCategory,iCell), &
                  snowEnthalpy(:,iCategory,iCell))
             call seaice_icepack_write_warnings(icepack_warnings_aborted())

          enddo ! iCategory

       else if (trim(config_column_physics_type) == "column_package") then

          ! has ice
          do iCategory = 1, nCategories

             iceAreaCategory(1,iCategory,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,iCategory,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,iCategory,iCell) = config_initial_snow_volume

             call colpkg_init_trcr(&
                  airTemperature(iCell), &
                  seaFreezingTemperature(iCell), &
                  initialSalinityProfile(:,iCell), &
                  initialMeltingTemperatureProfile(:,iCell), &
                  surfaceTemperature(1,iCategory,iCell), &
                  nIceLayers, &
                  nSnowLayers, &
                  iceEnthalpy(:,iCategory,iCell), &
                  snowEnthalpy(:,iCategory,iCell))

          enddo ! iCategory

       endif ! config_column_physics_type

    enddo ! iCell

    ! get vertical salinity profiles
    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iIceLayer = 1, nIceLayers
             iceSalinity(iIceLayer,iCategory,iCell) = initialSalinityProfile(iIceLayer,iCell)
          enddo ! iIceLayer
       enddo ! iCategory
    enddo ! iCell

  end subroutine init_ice_single_cell

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_ridging
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 28th July 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_ridging(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceDegreesToRadians

    use icepack_intfc, only: &
         icepack_init_enthalpy, &
         icepack_enthalpy_snow, &
         icepack_warnings_aborted

    use seaice_icepack, only: &
         seaice_icepack_write_warnings

    use seaice_column, only: &
         seaice_column_init_trcr, &
         seaice_column_enthalpy_snow

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         atmos_coupling, &
         ocean_coupling, &
         initial

    integer, pointer :: &
         nCellsSolve, &
         nCategories, &
         nIceLayers, &
         nSnowLayers

    real(kind=RKIND), dimension(:), pointer :: &
         latCell, &
         airTemperature, &
         seaSurfaceTemperature, &
         seaFreezingTemperature, &
         categoryThicknessLimits

    integer, dimension(:), pointer :: &
         landIceMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         initialSalinityProfile, &
         initialMeltingTemperatureProfile

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         surfaceTemperature, &
         iceEnthalpy, &
         snowEnthalpy, &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         iceSalinity

    real(kind=RKIND), pointer :: &
         config_initial_latitude_north, &
         config_initial_latitude_south

    real(kind=RKIND), dimension(:), allocatable :: &
         initialCategoryIceArea, &
         initialCategoryIceThickness

    integer :: &
         iCell, &
         iCategory, &
         iIceLayer, &
         iSnowLayer

    real(kind=RKIND), parameter :: &
         initialCategorySnowThickness = 0.0_RKIND

    character(len=strKIND), pointer :: &
         config_column_physics_type

    call MPAS_pool_get_config(configs, "config_column_physics_type", config_column_physics_type)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)
    call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmos_coupling)
    call MPAS_pool_get_subpool(block % structs, "initial", initial)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nCategories", nCategories)
    call MPAS_pool_get_dimension(mesh, "nIceLayers", nIceLayers)
    call MPAS_pool_get_dimension(mesh, "nSnowLayers", nSnowLayers)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)
    call MPAS_pool_get_array(tracers, "iceEnthalpy", iceEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "snowEnthalpy", snowEnthalpy, 1)
    call MPAS_pool_get_array(tracers, "iceSalinity", iceSalinity, 1)

    call MPAS_pool_get_array(atmos_coupling, "airTemperature", airTemperature)

    call MPAS_pool_get_array(ocean_coupling, "seaSurfaceTemperature", seaSurfaceTemperature)
    call MPAS_pool_get_array(ocean_coupling, "seaFreezingTemperature", seaFreezingTemperature)
    call MPAS_pool_get_array(ocean_coupling, "landIceMask", landIceMask)

    call MPAS_pool_get_array(initial, "categoryThicknessLimits", categoryThicknessLimits)
    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)
    call MPAS_pool_get_array(initial, "initialMeltingTemperatureProfile", initialMeltingTemperatureProfile)

    ! initial thickness categories
    categoryThicknessLimits(1) = 0.0_RKIND
    categoryThicknessLimits(2) = 0.6_RKIND
    categoryThicknessLimits(3) = 1.4_RKIND
    categoryThicknessLimits(4) = 2.4_RKIND
    categoryThicknessLimits(5) = 3.6_RKIND

    ! initial volumes and areas
    allocate(initialCategoryIceArea(nCategories))
    allocate(initialCategoryIceThickness(nCategories))

    initialCategoryIceArea(1) = 0.05_RKIND
    initialCategoryIceArea(2) = 0.1_RKIND
    initialCategoryIceArea(3) = 0.3_RKIND
    initialCategoryIceArea(4) = 0.35_RKIND
    initialCategoryIceArea(5) = 0.2_RKIND

    initialCategoryIceThickness(1) = (categoryThicknessLimits(1) + categoryThicknessLimits(1)) * 0.5_RKIND
    initialCategoryIceThickness(2) = (categoryThicknessLimits(2) + categoryThicknessLimits(2)) * 0.5_RKIND
    initialCategoryIceThickness(3) = (categoryThicknessLimits(3) + categoryThicknessLimits(3)) * 0.5_RKIND
    initialCategoryIceThickness(4) = (categoryThicknessLimits(4) + categoryThicknessLimits(4)) * 0.5_RKIND
    initialCategoryIceThickness(5) = categoryThicknessLimits(5) + 1.0_RKIND

    if (trim(config_column_physics_type) == "icepack") then
       do iCell = 1, nCellsSolve

          if (seaSurfaceTemperature(iCell) <= seaFreezingTemperature(iCell) + 0.2_RKIND .and. &
               (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
               latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) .and. &
               landIceMask(iCell) == 0) then

             ! has ice
             do iCategory = 1, nCategories

                iceAreaCategory(1,iCategory,iCell)    = initialCategoryIceArea(iCategory)
                iceVolumeCategory(1,iCategory,iCell)  = initialCategoryIceArea(iCategory) * initialCategoryIceThickness(iCategory)
                snowVolumeCategory(1,iCategory,iCell) = min(iceAreaCategory(1,iCategory,iCell) * initialCategorySnowThickness, &
                     0.2_RKIND * iceVolumeCategory(1,iCategory,iCell))

                call icepack_init_enthalpy(&
                     airTemperature(iCell), &
                     seaFreezingTemperature(iCell), &
                     initialSalinityProfile(:,iCell), &
                     initialMeltingTemperatureProfile(:,iCell), &
                     surfaceTemperature(1,iCategory,iCell), &
                     iceEnthalpy(:,iCategory,iCell), &
                     snowEnthalpy(:,iCategory,iCell))
                call seaice_icepack_write_warnings(icepack_warnings_aborted())

             enddo ! iCategory

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND

             surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)
             do iSnowLayer = 1, nSnowLayers
                snowEnthalpy(iSnowLayer,:,iCell) = icepack_enthalpy_snow(0.0_RKIND)
             end do

          endif

       enddo ! iCell
    else if (trim(config_column_physics_type) == "column_package") then
       do iCell = 1, nCellsSolve

          if (seaSurfaceTemperature(iCell) <= seaFreezingTemperature(iCell) + 0.2_RKIND .and. &
               (latCell(iCell) > config_initial_latitude_north * seaiceDegreesToRadians .or. &
               latCell(iCell) < config_initial_latitude_south * seaiceDegreesToRadians) .and. &
               landIceMask(iCell) == 0) then

             ! has ice
             do iCategory = 1, nCategories

                iceAreaCategory(1,iCategory,iCell)    = initialCategoryIceArea(iCategory)
                iceVolumeCategory(1,iCategory,iCell)  = initialCategoryIceArea(iCategory) * initialCategoryIceThickness(iCategory)
                snowVolumeCategory(1,iCategory,iCell) = min(iceAreaCategory(1,iCategory,iCell) * initialCategorySnowThickness, &
                     0.2_RKIND * iceVolumeCategory(1,iCategory,iCell))

                call seaice_column_init_trcr(&
                     airTemperature(iCell), &
                     seaFreezingTemperature(iCell), &
                     initialSalinityProfile(:,iCell), &
                     initialMeltingTemperatureProfile(:,iCell), &
                     surfaceTemperature(1,iCategory,iCell), &
                     nIceLayers, &
                     nSnowLayers, &
                     iceEnthalpy(:,iCategory,iCell), &
                     snowEnthalpy(:,iCategory,iCell))

             enddo ! iCategory

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND

             surfaceTemperature(1,:,iCell) = seaFreezingTemperature(iCell)
             do iSnowLayer = 1, nSnowLayers
                snowEnthalpy(iSnowLayer,:,iCell) = seaice_column_enthalpy_snow(0.0_RKIND)
             end do

          endif

       enddo ! iCell
    endif ! config_column_physics_type

    ! clean up
    deallocate(initialCategoryIceArea)
    deallocate(initialCategoryIceThickness)

    ! get vertical salinity profiles
    call MPAS_pool_get_subpool(block % structs, "initial", initial)
    call MPAS_pool_get_array(initial, "initialSalinityProfile", initialSalinityProfile)

    do iCell = 1, nCellsSolve
       do iCategory = 1, nCategories
          do iIceLayer = 1, nIceLayers
             iceSalinity(iIceLayer,iCategory,iCell) = initialSalinityProfile(iIceLayer,iCell)
          enddo ! iIceLayer
       enddo ! iCategory
    enddo ! iCell

  end subroutine init_ice_ridging!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_interior_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_uniform_interior(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         boundary

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    integer, dimension(:), pointer :: &
         interiorCell

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    call MPAS_pool_get_array(boundary, "interiorCell", interiorCell)

    do iCell = 1, nCells

       if (interiorCell(iCell) == 1) then

          ! has ice
          iceAreaCategory(1,:,iCell)    = config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = config_initial_snow_volume

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine init_ice_state_uniform_interior!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_circle_of_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_circle_of_ice(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), pointer :: &
         earthRadius, &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume

    logical, pointer :: &
         config_rotate_cartesian_grid

    integer :: &
         iCell

    real(kind=RKIND) :: &
         circle_radius

    call MPAS_pool_get_config(configs, "config_earth_radius", earthRadius)

    circle_radius = 0.1_RKIND*earthRadius

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, 'config_rotate_cartesian_grid', config_rotate_cartesian_grid)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    if (config_rotate_cartesian_grid) then

       do iCell = 1, nCells

          ! Put a circle of ice at the North Pole (which will be rotated to the grid equator)
          if (zCell(iCell) > 0.0_RKIND .and. sqrt(xCell(iCell)**2 + yCell(iCell)**2) < circle_radius) then

             ! has ice
             iceAreaCategory(1,:,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
             surfaceTemperature(1,:,iCell) = 1.0_RKIND

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND
             surfaceTemperature(1,:,iCell) = 0.0_RKIND

          endif

       enddo ! iCell

    else   ! no rotation

       do iCell = 1, nCells

          ! Put a circle of ice on the equator
          if (xCell(iCell) > 0.0_RKIND .and. sqrt(yCell(iCell)**2 + zCell(iCell)**2) < circle_radius) then   ! Greenwich meridian
!          if (xCell(iCell) < 0.0_RKIND .and. sqrt(yCell(iCell)**2 + zCell(iCell)**2) < circle_radius) then  ! Date Line

             ! has ice
             iceAreaCategory(1,:,iCell)    = config_initial_ice_area
             iceVolumeCategory(1,:,iCell)  = config_initial_ice_volume
             snowVolumeCategory(1,:,iCell) = config_initial_snow_volume
             surfaceTemperature(1,:,iCell) = 1.0_RKIND

          else

             ! no ice
             iceAreaCategory(1,:,iCell)    = 0.0_RKIND
             iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
             snowVolumeCategory(1,:,iCell) = 0.0_RKIND
             surfaceTemperature(1,:,iCell) = 0.0_RKIND

          endif

       enddo ! iCell

    endif ! config_rotate_cartesian_grid

  end subroutine init_ice_state_circle_of_ice!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th July 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_special(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceRadiansToDegrees

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell1, &
         iCell2, &
         iCell0, &
         iEdge1, &
         iEdge2, &
         iEdgeOnCell0, &
         iEdgeOnCell1, &
         iEdgeOnCell2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)

    iceAreaCategory = 0.0_RKIND
    iceVolumeCategory = 0.0_RKIND
    snowVolumeCategory = 0.0_RKIND

    !iceAreaCategory = 1.0e-20_RKIND
    !iceVolumeCategory = 1.0_RKIND
    !snowVolumeCategory = 1.0_RKIND


    iCell1 = 2050
    iCell0 = 2051
    iCell2 = 2052

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell1 = 1, nEdgesOnCell(iCell1)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell1,iCell1)) then
             iEdge1 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell2 = 1, nEdgesOnCell(iCell2)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell2,iCell2)) then
             iEdge2 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    iceAreaCategory(:,:,iCell1) = 1.0_RKIND
    iceVolumeCategory(:,:,iCell1) = 1.0_RKIND

    iceAreaCategory(:,:,iCell2) = 1.0_RKIND
    iceVolumeCategory(:,:,iCell2) = 1.0_RKIND

  end subroutine init_ice_state_special!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_state_uniform_ice
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th July 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_state_random_coverage(&
       block, &
       configs)!{{{

    use seaice_constants, only: &
         seaiceRadiansToDegrees

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         tracers

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory, &
         snowVolumeCategory, &
         surfaceTemperature

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south

    integer :: &
         iCell

    real(kind=RKIND) :: &
         random, &
         random2

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory, 1)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature, 1)

    do iCell = 1, nCells

       call random_number(random)

       if (random > 0.5_RKIND) then

          call random_number(random2)

          ! has ice
          iceAreaCategory(1,:,iCell)    = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_ice_area
          iceVolumeCategory(1,:,iCell)  = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_ice_volume
          snowVolumeCategory(1,:,iCell) = 10.0_RKIND**(-11.0_RKIND * random2)!config_initial_snow_volume
          surfaceTemperature(1,:,iCell) = 1.0_RKIND

       else

          ! no ice
          iceAreaCategory(1,:,iCell)    = 0.0_RKIND
          iceVolumeCategory(1,:,iCell)  = 0.0_RKIND
          snowVolumeCategory(1,:,iCell) = 0.0_RKIND
          surfaceTemperature(1,:,iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine init_ice_state_random_coverage!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_uniform
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_uniform(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary, &
         ocean_coupling

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_coupling, "landIceMaskVertex", landIceMaskVertex)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1 .and. &
           landIceMaskVertex(iVertex) == 0) then

          uVelocity(iVertex) = config_initial_uvelocity
          vVelocity(iVertex) = config_initial_vvelocity

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_velocity_uniform!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_ocean
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_ocean(&
       block, &
       configs)!{{{

    use seaice_mesh, only: &
         seaice_interpolate_cell_to_vertex

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary, &
         ocean_coupling

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uOceanVelocity, &
         vOceanVelocity

    integer :: &
         iVertex

    real(kind=RKIND), dimension(:), allocatable :: &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "ocean_coupling", ocean_coupling)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_coupling, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "vOceanVelocity", vOceanVelocity)
    call MPAS_pool_get_array(ocean_coupling, "landIceMaskVertex", landIceMaskVertex)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    allocate(uOceanVelocityVertex(nVertices))
    allocate(vOceanVelocityVertex(nVertices))

    call seaice_interpolate_cell_to_vertex(mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call seaice_interpolate_cell_to_vertex(mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1 .and. &
           landIceMaskVertex(iVertex) == 0) then

          uVelocity(iVertex) = uOceanVelocityVertex(iVertex)
          vVelocity(iVertex) = vOceanVelocityVertex(iVertex)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

    deallocate(uOceanVelocityVertex)
    deallocate(vOceanVelocityVertex)

  end subroutine init_ice_velocity_ocean!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_random
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_random(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call random_seed()

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then

          call random_number(uVelocity(iVertex))
          call random_number(vVelocity(iVertex))

          uVelocity(iVertex) = config_initial_uvelocity * (uVelocity(iVertex) * 2.0_RKIND - 1.0_RKIND)
          vVelocity(iVertex) = config_initial_vvelocity * (vVelocity(iVertex) * 2.0_RKIND - 1.0_RKIND)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_velocity_random!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_ice_velocity_special
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_ice_velocity_special(&
       block, &
       configs)!{{{

    type(block_type), intent(inout) :: &
         block !< Input/Output:

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs !< Input:

    type (MPAS_pool_type), pointer :: &
         mesh, &
         velocity_solver, &
         boundary

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), pointer :: &
         config_initial_uvelocity, &
         config_initial_vvelocity

    integer :: &
         iVertex

    integer :: &
         iCell1, &
         iCell2, &
         iCell0, &
         iEdge1, &
         iEdge2, &
         iEdgeOnCell0, &
         iEdgeOnCell1, &
         iEdgeOnCell2

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         edgesOnCell

    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "edgesOnCell", edgesOnCell)

    uVelocity = 0.0_RKIND
    vVelocity = 0.0_RKIND

    iCell1 = 2050
    iCell0 = 2051
    iCell2 = 2052

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell1 = 1, nEdgesOnCell(iCell1)
          write(*,*) iEdgeOnCell0, iEdgeOnCell1, edgesOnCell(iEdgeOnCell0,iCell0), edgesOnCell(iEdgeOnCell1,iCell1)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell1,iCell1)) then
             iEdge1 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    do iEdgeOnCell0 = 1, nEdgesOnCell(iCell0)
       do iEdgeOnCell2 = 1, nEdgesOnCell(iCell2)
          if (edgesOnCell(iEdgeOnCell0,iCell0) == edgesOnCell(iEdgeOnCell2,iCell2)) then
             iEdge2 = edgesOnCell(iEdgeOnCell0,iCell0)
          endif
       enddo ! iEdgeOnCell1
    enddo ! iEdgeOnCell0

    uVelocity(:) = 1.0_RKIND

  end subroutine init_ice_velocity_special!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  initial_halo_exchanges
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date Januray 13th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initial_halo_exchanges(domain)

    type(domain_type), intent(inout) :: &
         domain

    logical, pointer :: &
         pkgWeakActive

    call MPAS_pool_get_package(domain % blocklist % packages, "pkgWeakActive", pkgWeakActive)

    if (pkgWeakActive) then

       call MPAS_dmpar_field_halo_exch(domain, 'normalVectorTriangle')

    endif

  end subroutine initial_halo_exchanges

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! initialize_coupler_fields
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 14th March 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine initialize_coupler_fields(domain)

    use icepack_intfc, only: &
         icepack_liquidus_temperature
    use seaice_icepack, only: &             !echmod - use icepack_intfc directly
         seaice_icepack_init_ocean_conc, &
         seaice_icepack_initial_air_drag_coefficient
    use seaice_column, only: &
         seaice_column_liquidus_temperature, &
         seaice_column_init_ocean_conc, &
         seaice_column_initial_air_drag_coefficient

    use seaice_constants, only: &
         seaiceStefanBoltzmann, &
         seaiceFreshWaterFreezingPoint

    type(domain_type) :: domain

    type(block_type), pointer :: block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         atmosCoupling, &
         atmosForcing, &
         oceanCoupling, &
         atmosFluxes, &
         drag, &
         biogeochemistry, &
         aerosols

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocity, &
         vAirVelocity, &
         seaSurfaceSalinity, &
         seaFreezingTemperature, &
         seaSurfaceTemperature, &
         windSpeed, &
         airDragCoefficient, &
         longwaveUp

    real(kind=RKIND), dimension(:), pointer :: &
         oceanNitrateConc, &
         oceanSilicateConc, &
         oceanAmmoniumConc, &
         oceanDMSConc, &
         oceanDMSPConc, &
         oceanHumicsConc, &
         carbonToNitrogenRatioAlgae, &
         carbonToNitrogenRatioDON

    real(kind=RKIND), dimension(:,:), pointer :: &
         atmosBioFluxes, &
         atmosBlackCarbonFlux, &
         atmosDustFlux, &
         oceanAlgaeConc, &
         oceanDOCConc, &
         oceanDICConc, &
         oceanDONConc, &
         oceanParticulateIronConc, &
         oceanDissolvedIronConc, &
         oceanZAerosolConc, &
         atmosAerosolFlux

    integer, pointer :: &
         maxAlgaeType, &
         maxDOCType, &
         maxDICType, &
         maxDONType, &
         maxIronType, &
         maxBCType, &
         maxDustType, &
         maxAerosolType, &
         nZBGCTracers

         logical, pointer :: &
         config_do_restart, &
         config_use_column_biogeochemistry, &
         config_use_column_physics, &
         config_use_aerosols, &
         config_use_zaerosols

    integer, pointer :: &
         nCells

    integer :: &
         iCell, &
         iBio, &
         indexj

    character(len=strKIND), pointer :: &
         config_column_physics_type

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_column_physics_type", config_column_physics_type)

       !-------------------------------------------------------------
       ! Physics fluxes received from atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCoupling)
       call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmosForcing)

       call MPAS_pool_get_array(atmosCoupling, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCoupling, "vAirVelocity", vAirVelocity)

       call MPAS_pool_get_array(atmosForcing, "windSpeed", windSpeed)

       windSpeed = sqrt(uAirVelocity**2 + vAirVelocity**2)

       !-------------------------------------------------------------
       ! Physics fluxes received from ocean
       !-------------------------------------------------------------

       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCoupling)
       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

       call MPAS_pool_get_dimension(mesh, "nCells", nCells)

       call MPAS_pool_get_array(oceanCoupling, "seaSurfaceSalinity", seaSurfaceSalinity)
       call MPAS_pool_get_array(oceanCoupling, "seaFreezingTemperature", seaFreezingTemperature)
       call MPAS_pool_get_array(oceanCoupling, "seaSurfaceTemperature", seaSurfaceTemperature)

       call MPAS_pool_get_config(block % configs, "config_do_restart", config_do_restart)

       if (trim(config_column_physics_type) == "icepack") then
          do iCell = 1, nCells
             seaFreezingTemperature(iCell) = icepack_liquidus_temperature(seaSurfaceSalinity(iCell))
          enddo ! iCell
       else if (trim(config_column_physics_type) == "column_package") then
          do iCell = 1, nCells
             seaFreezingTemperature(iCell) = seaice_column_liquidus_temperature(seaSurfaceSalinity(iCell))
          enddo ! iCell
       endif ! config_column_physics_type

       ! sea surface temperature is not initialized if we're restarting
       if (.not. config_do_restart) then
          seaSurfaceTemperature = seaFreezingTemperature
       endif

       !-------------------------------------------------------------
       ! fluxes sent to atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_config(block % configs, "config_use_column_physics", config_use_column_physics)

       if (config_use_column_physics) then

          call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxes)
          call MPAS_pool_get_array(atmosFluxes, "longwaveUp", longwaveUp)
          longwaveUp = -seaiceStefanBoltzmann * seaiceFreshWaterFreezingPoint**4

          call MPAS_pool_get_subpool(block % structs, "drag", drag)
          call MPAS_pool_get_array(drag, "airDragCoefficient", airDragCoefficient)
          if (trim(config_column_physics_type) == "icepack") then
             airDragCoefficient = seaice_icepack_initial_air_drag_coefficient()
          else if (trim(config_column_physics_type) == "column_package") then
             airDragCoefficient = seaice_column_initial_air_drag_coefficient()
          endif ! config_column_physics_type

       endif

       !-------------------------------------------------------------
       ! Bio fluxes received from atmosphere
       !-------------------------------------------------------------

       call MPAS_pool_get_config(block % configs, "config_use_aerosols", config_use_aerosols)

       if (config_use_aerosols) then
          call MPAS_pool_get_subpool(block % structs, "aerosols", aerosols)
          call MPAS_pool_get_array(aerosols, "atmosAerosolFlux", atmosAerosolFlux)

          atmosAerosolFlux = 1.e-12_RKIND
       endif

       call MPAS_pool_get_config(block % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
       call MPAS_pool_get_config(block % configs, "config_use_zaerosols", config_use_zaerosols)

       if (config_use_column_biogeochemistry .or. config_use_zaerosols) then

          call MPAS_pool_get_subpool(block % structs, "biogeochemistry", biogeochemistry)

          call MPAS_pool_get_array(biogeochemistry, "atmosBioFluxes", atmosBioFluxes)
          call MPAS_pool_get_array(biogeochemistry, "atmosBlackCarbonFlux", atmosBlackCarbonFlux)
          call MPAS_pool_get_array(biogeochemistry, "atmosDustFlux", atmosDustFlux)

          call MPAS_pool_get_dimension(block % dimensions, "maxAlgaeType", maxAlgaeType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDOCType", maxDOCType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDICType", maxDICType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDONType", maxDONType)
          call MPAS_pool_get_dimension(block % dimensions, "maxIronType", maxIronType)
          call MPAS_pool_get_dimension(block % dimensions, "maxBCType", maxBCType)
          call MPAS_pool_get_dimension(block % dimensions, "maxDustType", maxDustType)
          call MPAS_pool_get_dimension(block % dimensions, "maxAerosolType", maxAerosolType)
          call MPAS_pool_get_dimension(block % dimensions, "nZBGCTracers", nZBGCTracers)

          atmosBioFluxes = 0.0_RKIND
          do iCell = 1, nCells
             do iBio = 1, maxBCType
                atmosBlackCarbonFlux(iBio,iCell) =  1.e-12_RKIND
             enddo
             do iBio = 1, maxDustType
                atmosDustFlux(iBio,iCell) =  1.e-13_RKIND
             enddo
          enddo ! iCell

          !-------------------------------------------------------------
          ! Bio Concentrations received from ocean
          !-------------------------------------------------------------

          call MPAS_pool_get_array(biogeochemistry, "oceanAlgaeConc", oceanAlgaeConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDOCConc", oceanDOCConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDICConc", oceanDICConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDONConc", oceanDONConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanParticulateIronConc", oceanParticulateIronConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDissolvedIronConc", oceanDissolvedIronConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanNitrateConc", oceanNitrateConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanSilicateConc", oceanSilicateConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanAmmoniumConc", oceanAmmoniumConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDMSConc", oceanDMSConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanDMSPConc", oceanDMSPConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanHumicsConc", oceanHumicsConc)
          call MPAS_pool_get_array(biogeochemistry, "oceanZAerosolConc", oceanZAerosolConc)
          call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioAlgae", carbonToNitrogenRatioAlgae)
          call MPAS_pool_get_array(biogeochemistry, "carbonToNitrogenRatioDON", carbonToNitrogenRatioDON)

          if (trim(config_column_physics_type) == "icepack") then
             do iCell = 1, nCells

                call seaice_icepack_init_ocean_conc(&
                     oceanAmmoniumConc(iCell), &
                     oceanDMSPConc(iCell), &
                     oceanDMSConc(iCell), &
                     oceanAlgaeConc(:,iCell), &
                     oceanDOCConc(:,iCell), &
                     oceanDICConc(:,iCell), &
                     oceanDONConc(:,iCell), &
                     oceanDissolvedIronConc(:,iCell), &
                     oceanParticulateIronConc(:,iCell), &
                     oceanHumicsConc(iCell), &
                     oceanNitrateConc(iCell), &
                     oceanSilicateConc(iCell),&
                     oceanZAerosolConc(:,iCell), &
                     carbonToNitrogenRatioAlgae(:), &
                     carbonToNitrogenRatioDON(:))

             enddo ! iCell
          else if (trim(config_column_physics_type) == "column_package") then
             do iCell = 1, nCells

                call seaice_column_init_ocean_conc(&
                     oceanAmmoniumConc(iCell), &
                     oceanDMSPConc(iCell), &
                     oceanDMSConc(iCell), &
                     oceanAlgaeConc(:,iCell), &
                     oceanDOCConc(:,iCell), &
                     oceanDICConc(:,iCell), &
                     oceanDONConc(:,iCell), &
                     oceanDissolvedIronConc(:,iCell), &
                     oceanParticulateIronConc(:,iCell), &
                     oceanHumicsConc(iCell), &
                     oceanNitrateConc(iCell), &
                     oceanSilicateConc(iCell),&
                     oceanZAerosolConc(:,iCell), &
                     maxDICType, &
                     maxDONType, &
                     maxIronType, &
                     maxAerosolType, &
                     carbonToNitrogenRatioAlgae(:), &
                     carbonToNitrogenRatioDON(:))

             enddo ! iCell

          endif ! config_column_physics_type

       endif ! config_use_column_biogeochemistry

       block => block % next
    enddo

  end subroutine initialize_coupler_fields

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! init_test_ice_shelf_mask
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th April 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine init_test_ice_shelf_mask(domain)

    use seaice_constants, only: &
         seaiceDegreesToRadians, &
         seaiceRadiansToDegrees

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         oceanCouplingPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         landIceMask, &
         latCell, &
         lonCell

    integer, pointer :: &
         nCells

    integer :: &
         iCell

    real(kind=RKIND) :: &
         lonLimitWest = -66.10_RKIND, &
         lonLimitEast = -11.26_RKIND , &
         latLimit     = -71.94_RKIND

    logical, pointer :: &
         config_use_test_ice_shelf

    call MPAS_pool_get_config(domain % configs, "config_use_test_ice_shelf", config_use_test_ice_shelf)

    if (config_use_test_ice_shelf) then

       lonLimitWest = lonLimitWest + 360.0_RKIND
       lonLimitEast = lonLimitEast + 360.0_RKIND

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

          call MPAS_pool_get_array(oceanCouplingPool, "landIceMask", landIceMask)

          call MPAS_pool_get_array(meshPool, "latCell", latCell)
          call MPAS_pool_get_array(meshPool, "lonCell", lonCell)

          call MPAS_pool_get_dimension(block % dimensions, "nCells", nCells)

          do iCell = 1, nCells

             if (lonCell(iCell) > lonLimitWest * seaiceDegreesToRadians .and. &
                 lonCell(iCell) < lonLimitEast * seaiceDegreesToRadians .and. &
                 latCell(iCell) < latLimit     * seaiceDegreesToRadians) then

                landIceMask(iCell) = 1

             else

                landIceMask(iCell) = 0

             endif

          enddo ! iCell

          block => block % next
       enddo

    endif

  end subroutine init_test_ice_shelf_mask

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_check_constants_coupled
!
!> \brief
!> \author Elizabeth Hunke, LANL
!> \date 29 September 2023
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_check_constants_coupled

    use seaice_constants, only: &
         seaiceIceSnowEmissivity

#ifdef CCSMCOUPLED
    if (seaiceIceSnowEmissivity /= 1._RKIND) then
      call mpas_log_write(&
           'seaice_check_constants_coupled: seaiceIceSnowEmissivity out of bounds', &
           messageType=MPAS_LOG_CRIT)
    endif
#endif

  end subroutine seaice_check_constants_coupled

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_check_configs_coupled
!
!> \brief
!> \author Elizabeth Hunke, LANL
!> \date 3 November 2023
!> \details
!> Checks critical configuration settings for coupled simulations.
!> Aborts if they would cause conservation or other errors.
!> Warns for some non-critical, non-default settings.
!> Config consistency checks are in subroutine check_column_package_configs
!
!-----------------------------------------------------------------------

  subroutine seaice_check_configs_coupled(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! abort
    character(len=strKIND), pointer :: &
         config_frazil_coupling_type, &
         config_thermodynamics_type, &
         config_sea_freezing_temperature_type, &
         config_itd_conversion_type

    logical, pointer :: &
         config_use_column_shortwave, &
         config_use_column_vertical_thermodynamics, &
         config_use_column_ridging, &
         config_use_form_drag, &
         config_use_ocean_mixed_layer, &
         config_calc_surface_stresses, &
         config_calc_surface_temperature, &
         config_include_pond_freshwater_feedback, &
         config_update_ocean_fluxes

    ! write warnings
    character(len=strKIND), pointer :: &
         config_column_physics_type, &
         config_shortwave_type, &
         config_snow_redistribution_scheme, &
         config_category_bounds_type, &
         config_ridging_participation_function, &
         config_ridging_redistribution_function, &
         config_ocean_heat_transfer_type, &
         config_salt_flux_coupling_type, &       ! only available in icepack
         config_ocean_surface_type               ! only used in ice_comp_mct.F

    logical, pointer :: &
         config_use_column_physics, &
         config_use_column_itd_thermodynamics, & ! must be false for F case
         config_use_prescribed_ice, &
         config_use_snicar_ad, &
         config_use_shortwave_redistribution, &  ! only available in icepack
         config_use_column_snow_tracers, &
         config_use_snow_liquid_ponds, &
         config_use_high_frequency_coupling, &
         config_limit_air_temperatures, &        ! only used for CORE forcing
         config_use_congelation_basal_melt, &    ! deprecate with colpkg
         config_use_lateral_melt, &              ! deprecate with colpkg
         config_use_latent_processes             ! deprecate with colpkg

#ifdef CCSMCOUPLED
    ! abort
    call MPAS_pool_get_config(domain % configs, "config_category_bounds_type", config_category_bounds_type)
    call MPAS_pool_get_config(domain % configs, "config_frazil_coupling_type", config_frazil_coupling_type)
    call MPAS_pool_get_config(domain % configs, "config_thermodynamics_type", config_thermodynamics_type)
    call MPAS_pool_get_config(domain % configs, "config_sea_freezing_temperature_type", config_sea_freezing_temperature_type)
    call MPAS_pool_get_config(domain % configs, "config_itd_conversion_type", config_itd_conversion_type)
    call MPAS_pool_get_config(domain % configs, "config_use_column_shortwave", config_use_column_shortwave)
    call MPAS_pool_get_config(domain % configs, "config_use_column_vertical_thermodynamics", &
                                                 config_use_column_vertical_thermodynamics)
    call MPAS_pool_get_config(domain % configs, "config_use_column_ridging", config_use_column_ridging)
    call MPAS_pool_get_config(domain % configs, "config_use_form_drag", config_use_form_drag)
    call MPAS_pool_get_config(domain % configs, "config_use_ocean_mixed_layer", config_use_ocean_mixed_layer)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_stresses", config_calc_surface_stresses)
    call MPAS_pool_get_config(domain % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
    call MPAS_pool_get_config(domain % configs, "config_include_pond_freshwater_feedback", config_include_pond_freshwater_feedback)
    call MPAS_pool_get_config(domain % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)

    if (trim(config_category_bounds_type) == "single category") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_category_bounds_type must not be single category', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (trim(config_frazil_coupling_type) /= "external") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_frazil_coupling_type must be external', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (trim(config_thermodynamics_type) /= "mushy") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_thermodynamics_type must be mushy', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (trim(config_sea_freezing_temperature_type) /= "mushy") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_sea_freezing_temperature_type must be mushy', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (trim(config_itd_conversion_type) /= "linear remap") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_itd_conversion_type must be linear remap', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_use_column_shortwave) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_shortwave must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_use_column_vertical_thermodynamics) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_vertical_thermodynamics must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_use_column_ridging) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_ridging must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (config_use_form_drag) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_form_drag = true is not supported', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (config_use_ocean_mixed_layer) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_ocean_mixed_layer must be false in coupled runs', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_calc_surface_stresses) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_calc_surface_stresses must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_calc_surface_temperature) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_calc_surface_temperature must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (config_include_pond_freshwater_feedback) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_include_pond_freshwater_feedback must be false', &
           messageType=MPAS_LOG_CRIT)
    endif

    if (.not. config_update_ocean_fluxes) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_update_ocean_fluxes must be true', &
           messageType=MPAS_LOG_CRIT)
    endif

    ! write warnings
    call MPAS_pool_get_config(domain % configs, "config_column_physics_type", config_column_physics_type)
    call MPAS_pool_get_config(domain % configs, "config_shortwave_type", config_shortwave_type)
    call MPAS_pool_get_config(domain % configs, "config_snow_redistribution_scheme", config_snow_redistribution_scheme)
    call MPAS_pool_get_config(domain % configs, "config_category_bounds_type", config_category_bounds_type)
    call MPAS_pool_get_config(domain % configs, "config_ridging_participation_function", config_ridging_participation_function)
    call MPAS_pool_get_config(domain % configs, "config_ridging_redistribution_function", config_ridging_redistribution_function)
    call MPAS_pool_get_config(domain % configs, "config_ocean_heat_transfer_type", config_ocean_heat_transfer_type)
    call MPAS_pool_get_config(domain % configs, "config_salt_flux_coupling_type", config_salt_flux_coupling_type)
    call MPAS_pool_get_config(domain % configs, "config_ocean_surface_type", config_ocean_surface_type)

    call MPAS_pool_get_config(domain % configs, "config_use_column_physics", config_use_column_physics)
    call MPAS_pool_get_config(domain % configs, "config_use_column_itd_thermodynamics", config_use_column_itd_thermodynamics)
    call MPAS_pool_get_config(domain % configs, "config_use_prescribed_ice", config_use_prescribed_ice)
    call MPAS_pool_get_config(domain % configs, "config_use_column_snow_tracers", config_use_column_snow_tracers)
    call MPAS_pool_get_config(domain % configs, "config_use_snow_liquid_ponds", config_use_snow_liquid_ponds)
    call MPAS_pool_get_config(domain % configs, "config_use_snicar_ad", config_use_snicar_ad)
    call MPAS_pool_get_config(domain % configs, "config_use_shortwave_redistribution", config_use_shortwave_redistribution)
    call MPAS_pool_get_config(domain % configs, "config_use_high_frequency_coupling", config_use_high_frequency_coupling)
    call MPAS_pool_get_config(domain % configs, "config_limit_air_temperatures", config_limit_air_temperatures)
    call MPAS_pool_get_config(domain % configs, "config_use_latent_processes", config_use_latent_processes)             ! colpkg
    call MPAS_pool_get_config(domain % configs, "config_use_lateral_melt", config_use_lateral_melt)                     ! colpkg
    call MPAS_pool_get_config(domain % configs, "config_use_congelation_basal_melt", config_use_congelation_basal_melt) ! colpkg

    if (trim(config_column_physics_type) /= "icepack") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_column_physics_type /= icepack', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_shortwave_type) /= "dEdd") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_shortwave_type /= dEdd', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_snow_redistribution_scheme) /= "ITDrdg") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_snow_redistribution_scheme /= ITDrdg', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_category_bounds_type) /= "original") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_category_bounds_type /= original', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_ridging_participation_function) /= "exponential") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_ridging_participation_function /= exponential', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_ridging_redistribution_function) /= "exponential") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_ridging_redistribution_function /= exponential', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_ocean_heat_transfer_type) /= "constant") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_ocean_heat_transfer_type /= constant', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_salt_flux_coupling_type) /= "constant") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_salt_flux_coupling_type /= constant', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_ocean_surface_type) /= "free") then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_ocean_surface_type /= free', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_column_physics) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_physics = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_column_itd_thermodynamics) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_itd_thermodynamics = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (config_use_prescribed_ice) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_prescribed_ice = true', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_column_snow_tracers) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_column_snow_tracers = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_snow_liquid_ponds) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_snow_liquid_ponds = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_snicar_ad) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_snicar_ad = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (config_use_shortwave_redistribution) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_shortwave_redistribution = true', &
           messageType=MPAS_LOG_WARN)
    endif

    if (.not. config_use_high_frequency_coupling) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_high_frequency_coupling = false', &
           messageType=MPAS_LOG_WARN)
    endif

    if (config_limit_air_temperatures) then
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_limit_air_temperatures = true', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_column_physics_type) == 'icepack' .and. .not. config_use_latent_processes) then ! colpkg
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_latent_processes is not available in icepack', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_column_physics_type) == 'icepack' .and. .not. config_use_lateral_melt) then ! colpkg
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_lateral_melt is not available in icepack', &
           messageType=MPAS_LOG_WARN)
    endif

    if (trim(config_column_physics_type) == 'icepack' .and. .not. config_use_congelation_basal_melt) then ! colpkg
      call mpas_log_write(&
           'seaice_check_configs_coupled: config_use_congelation_basal_melt is not available in icepack', &
           messageType=MPAS_LOG_WARN)
    endif

#endif

  end subroutine seaice_check_configs_coupled

!-----------------------------------------------------------------------

end module seaice_initialize
