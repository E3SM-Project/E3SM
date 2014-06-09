module cice_initialize

  use mpas_grid_types
  use mpas_timekeeping

  implicit none

  private
  public :: &
       cice_init

contains

  !--------------------------------------------------------------------------

  subroutine cice_init(block, clock, dt)

    use mpas_configure  

    use cice_mesh, only: init_boundary
    use cice_advection, only: init_advection
    use cice_dynamics, only: init_dynamics
    use cice_testing, only: init_spherical_test_case
    use cice_diagnostics, only: test_case_output, prepare_plot_variables
    use cice_forcing, only: init_forcing, get_forcing

    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock

    real(kind=RKIND), intent(in) :: &
         dt
    
    type (MPAS_pool_type), pointer :: &
         configs

    logical, pointer :: &
         config_use_velocity_solver, &
         config_use_advection, &
         config_use_forcing

    configs => block % configs

    call MPAS_pool_get_config(configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_use_forcing", config_use_forcing)

    call init_junk_values(block)

    ! calculate boundary variables
    call init_boundary(block)

    ! init dynamics
    if (config_use_velocity_solver) &
         call init_dynamics(block, dt)
    
    ! init advection
    if (config_use_advection) &
         call init_advection(block)
    
    ! init forcing
    if (config_use_forcing) then

       call init_forcing(&
            block)
       
       call get_forcing(&
            block, &
            clock)
       
    endif

    ! init ice state
    call init_ice_state(&
       block, &
       configs)
    
    call test_case_output(block, 0)

    call prepare_plot_variables(block)

  end subroutine cice_init

  !--------------------------------------------------------------------------

  subroutine init_junk_values(block)

    type(block_type), intent(inout) :: block

    type(MPAS_pool_type), pointer :: &
         mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34

    integer, pointer :: &
         nCells
    
    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    njunk = nCells + 1
    
    areaCell(njunk) = value_junk

  end subroutine init_junk_values

  !--------------------------------------------------------------------------
  
  subroutine init_ice_state(&
       block, &
       configs)

    use mpas_configure
    use cice_testing, only: &
         init_square_test_case_hex, &
         init_square_point_test_case_hex

    type(block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer :: configs

    character(len=strKIND), pointer :: &
         config_initial_condition_type, &
         config_initial_velocity_type

    integer, dimension(:), pointer :: &
         interiorVertex

    call MPAS_pool_get_config(configs, "config_initial_condition_type", config_initial_condition_type)
    call MPAS_pool_get_config(configs, "config_initial_velocity_type", config_initial_velocity_type)

    ! set volumes/areas
    if (trim(config_initial_condition_type) == "uniform") then

       call init_ice_state_uniform_ice(&
            block, &
            configs)

    else if (trim(config_initial_condition_type) == "circle") then

       call init_ice_state_circle_of_ice(&
            block, &
            configs)

    else if (trim(config_initial_condition_type) == "square") then

       call init_square_test_case_hex(&
            block, &
            configs)

    else if (trim(config_initial_condition_type) == "square_point") then

       call init_square_point_test_case_hex(&
            block, &
            configs)

    else
       
       write(*,*) "Ice initial condition type unknown!"
       write(*,*) "config_initial_condition_type: ", trim(config_initial_condition_type)
       stop

    endif

    ! set velocities
    if (trim(config_initial_velocity_type) == "uniform") then

       call init_ice_velocity_uniform(&
            block, &
            configs)

    else if (trim(config_initial_velocity_type) == "ocean") then

       call init_ice_velocity_ocean(&
            block, &
            configs)

    else

       write(*,*) "Ice initial velocity type unknown!"
       write(*,*) "config_initial_velocity_type: ", trim(config_initial_velocity_type)
       stop

    endif

  end subroutine init_ice_state

  !--------------------------------------------------------------------------

  subroutine init_ice_state_uniform_ice(&
       block, &
       configs)

    use mpas_constants, only: &
         pii

    type(block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

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
         iCell

    real(kind=RKIND), parameter :: &
         rad_to_deg = pii / 180.0_RKIND

    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "latCell", latCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory)

    do iCell = 1, nCells

       if (latCell(iCell) >= config_initial_latitude_north * rad_to_deg .or. &
           latCell(iCell) <= config_initial_latitude_south * rad_to_deg) then

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

  end subroutine init_ice_state_uniform_ice

  !--------------------------------------------------------------------------

  subroutine init_ice_state_circle_of_ice(&
       block, &
       configs)

    type(block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

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
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume

    integer :: &
         iCell
    
    real(kind=RKIND), parameter :: &
         circle_radius = 1.0e6_RKIND
   
    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)

    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory)
    call MPAS_pool_get_array(tracers, "snowVolumeCategory", snowVolumeCategory)
    call MPAS_pool_get_array(tracers, "surfaceTemperature", surfaceTemperature)

    do iCell = 1, nCells
       
       if (sqrt(xCell(iCell)**2 + zCell(iCell)**2) < circle_radius .and. yCell(iCell) > 0.0_RKIND) then
          
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
    
  end subroutine init_ice_state_circle_of_ice

  !--------------------------------------------------------------------------

  subroutine init_ice_velocity_uniform(&
       block, &
       configs)
    
    type(block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

    type (MPAS_pool_type), pointer :: &
         mesh, &
         hexdyn, &
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
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    do iVertex = 1, nVertices

       if (interiorVertex(iVertex) == 1) then
          
          uVelocity(iVertex) = config_initial_uvelocity
          vVelocity(iVertex) = config_initial_vvelocity
          
       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

  end subroutine init_ice_velocity_uniform

  !--------------------------------------------------------------------------

  subroutine init_ice_velocity_ocean(&
       block, &
       configs)

    use cice_mesh, only: &
         interpolate_cell_to_vertex
 
    type(block_type), intent(inout) :: block

    type (MPAS_pool_type), pointer, intent(in) :: &
         configs

    type (MPAS_pool_type), pointer :: &
         mesh, &
         hexdyn, &
         boundary, &
         ocean_forcing

    integer, pointer :: &
         nVertices

    integer, dimension(:), pointer :: &
         interiorVertex

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
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "ocean_forcing", ocean_forcing)

    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_array(ocean_forcing, "uOceanVelocity", uOceanVelocity)
    call MPAS_pool_get_array(ocean_forcing, "vOceanVelocity", vOceanVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)

    allocate(uOceanVelocityVertex(nVertices))
    allocate(vOceanVelocityVertex(nVertices))

    call interpolate_cell_to_vertex(mesh, &
         uOceanVelocityVertex, &
         uOceanVelocity)

    call interpolate_cell_to_vertex(mesh, &
         vOceanVelocityVertex, &
         vOceanVelocity)

    do iVertex = 1, nVertices
    
       if (interiorVertex(iVertex) == 1) then

          uVelocity(iVertex) = uOceanVelocityVertex(iVertex)
          vVelocity(iVertex) = vOceanVelocityVertex(iVertex)

       else

          uVelocity(iVertex) = 0.0_RKIND
          vVelocity(iVertex) = 0.0_RKIND

       endif

    enddo ! iVertex

    deallocate(uOceanVelocityVertex)
    deallocate(vOceanVelocityVertex)

  end subroutine init_ice_velocity_ocean

  !--------------------------------------------------------------------------

end module cice_initialize
