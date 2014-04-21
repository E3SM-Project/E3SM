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
    !use cice_advection, only: init_advection
    use cice_dynamics, only: init_dynamics
    use cice_testing, only: init_spherical_test_case
    use cice_diagnostics, only: writeout_state
    use cice_forcing, only: init_forcing, get_forcing, forcing_file_comparison
    use cice_therm_zerolayer, only: init_semtner_forcing, init_semtner_test_case  
    use cice_unit_test, only: perform_unit_test

    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock

    real(kind=RKIND), intent(in) :: &
         dt
    
    type (MPAS_pool_type), pointer :: &
         mesh, &
         atmos_forcing, &
         icestate, &
         hexdyn, &
         boundary, &
         configs

    logical, pointer :: &
         config_use_velocity_solver, &
         config_use_advection, &
         config_use_forcing, &
         config_perform_unit_test

    character(len=strKIND), pointer :: &
         config_unit_test_type, &
         config_unit_test_subtype

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)
    call MPAS_pool_get_subpool(block % structs, "icestate", icestate)
    call MPAS_pool_get_subpool(block % structs, "hexdyn", hexdyn)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)

    configs => block % configs

    call MPAS_pool_get_config(configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(configs, "config_perform_unit_test", config_perform_unit_test)
    call MPAS_pool_get_config(configs, "config_unit_test_type", config_unit_test_type)
    call MPAS_pool_get_config(configs, "config_unit_test_subtype", config_unit_test_subtype)

    call init_junk_values(mesh)

    ! calculate boundary variables
    call init_boundary(block)

    ! init dynamics
    if (config_use_velocity_solver) &
         call init_dynamics(block, dt)
    
    ! init advection
    !if (config_use_advection) &
    !     call init_advection(mesh)
    
    ! init forcing
    if (config_use_forcing) then

       call init_forcing(&
            mesh, &
            atmos_forcing)
       
       call get_forcing(&
            mesh, &
            clock, &
            atmos_forcing)
       
    endif

    ! init ice state
    call init_ice_state(&
         configs, &
         mesh, &
         icestate, &
         hexdyn, &
         boundary, &
         atmos_forcing)
    
    if (config_perform_unit_test) &
         call perform_unit_test(block, config_unit_test_type, config_unit_test_subtype)
    
    ! write out the initial condition
    !call writeout_state(block, 0)

  end subroutine cice_init

  !--------------------------------------------------------------------------

  subroutine init_junk_values(mesh)

    type(MPAS_pool_type), intent(inout) :: mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34

    integer, pointer :: &
         nCells
    
    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    njunk = nCells + 1
    
    areaCell(njunk) = value_junk

  end subroutine init_junk_values

  !--------------------------------------------------------------------------
  
  subroutine init_ice_state(&
       configs, &
       mesh, &
       icestate, &
       hexdyn, &
       boundary, &
       atmos_forcing)

    use mpas_configure

    type (MPAS_pool_type), intent(in) :: mesh
    type (MPAS_pool_type), pointer :: icestate
    type (MPAS_pool_type), pointer :: hexdyn
    type (MPAS_pool_type), pointer :: boundary
    type (MPAS_pool_type), pointer :: atmos_forcing
    type (MPAS_pool_type), pointer :: configs

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         tracers                  ! (nTracers, nVertLevels, nCells)

    integer, pointer :: &
         index_iceAreaCell, &
         index_iceVolumeCell, &
         index_snowVolumeCell

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    character(len=strKIND), pointer :: &
         config_initial_condition_type, &
         config_initial_velocity_type

    real(kind=RKIND), pointer :: &
         config_initial_ice_area, &
         config_initial_ice_volume, &
         config_initial_snow_volume, &
         config_initial_latitude_north, &
         config_initial_latitude_south, &
         config_initial_uvelocity, &
         config_initial_vvelocity

    call MPAS_pool_get_array(icestate, "tracers", tracers)
    call MPAS_pool_get_array(icestate, "index_iceAreaCell", index_iceAreaCell)
    call MPAS_pool_get_array(icestate, "index_iceVolumeCell", index_iceVolumeCell)
    call MPAS_pool_get_array(icestate, "index_snowVolumeCell", index_snowVolumeCell)

    call MPAS_pool_get_array(hexdyn, "uVelocity", uVelocity)
    call MPAS_pool_get_array(hexdyn, "vVelocity", vVelocity)

    call MPAS_pool_get_config(configs, "config_initial_condition_type", config_initial_condition_type)
    call MPAS_pool_get_config(configs, "config_initial_velocity_type", config_initial_velocity_type)
    call MPAS_pool_get_config(configs, "config_initial_ice_area", config_initial_ice_area)
    call MPAS_pool_get_config(configs, "config_initial_ice_volume", config_initial_ice_volume)
    call MPAS_pool_get_config(configs, "config_initial_snow_volume", config_initial_snow_volume)
    call MPAS_pool_get_config(configs, "config_initial_latitude_north", config_initial_latitude_north)
    call MPAS_pool_get_config(configs, "config_initial_latitude_south", config_initial_latitude_south)
    call MPAS_pool_get_config(configs, "config_initial_uvelocity", config_initial_uvelocity)
    call MPAS_pool_get_config(configs, "config_initial_vvelocity", config_initial_vvelocity)

    ! set volumes/areas
    if (trim(config_initial_condition_type) == "uniform") then

       call init_ice_state_uniform_ice(&
            mesh, &
            tracers(index_iceAreaCell,1,:), &
            tracers(index_iceVolumeCell,1,:), &
            tracers(index_snowVolumeCell,1,:), &
            config_initial_ice_area, &
            config_initial_ice_volume, &
            config_initial_snow_volume, &
            config_initial_latitude_north, &
            config_initial_latitude_south)

    else if (trim(config_initial_condition_type) == "circle") then

       call init_ice_state_circle_of_ice(&
            mesh, &
            tracers(index_iceAreaCell,1,:), &
            tracers(index_iceVolumeCell,1,:), &
            tracers(index_snowVolumeCell,1,:), &
            config_initial_ice_area, &
            config_initial_ice_volume, &
            config_initial_snow_volume)

    else
       
       write(*,*) "Ice initial condition type unknown!"
       write(*,*) "config_initial_condition_type: ", trim(config_initial_condition_type)
       stop

    endif

    ! set velocities
    if (trim(config_initial_velocity_type) == "uniform") then

       call init_ice_velocity_uniform(&
            mesh, &
            uVelocity, &
            vVelocity, &
            config_initial_uvelocity, &
            config_initial_vvelocity)

    else

       write(*,*) "Ice initial velocity type unknown!"
       write(*,*) "config_initial_velocity_type: ", trim(config_initial_velocity_type)
       stop

    endif

  end subroutine init_ice_state

  !--------------------------------------------------------------------------

  subroutine init_ice_state_uniform_ice(&
       mesh, &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       iceAreaInitial, &
       iceVolumeInitial, &
       snowVolumeInitial, &
       latitudeLimitNorth, &
       latitudeLimitSouth)

    type (MPAS_pool_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    real(kind=RKIND), intent(in) :: &
         iceAreaInitial, &
         iceVolumeInitial, &
         snowVolumeInitial, &
         latitudeLimitNorth, &
         latitudeLimitSouth

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         pi = 3.14_RKIND, &
         rad_to_deg = pi / 180.0_RKIND

    integer, pointer :: &
         nCells

    real(kind=RKIND), dimension(:), pointer :: &
         latCell

    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "latCell", latCell)

    do iCell = 1, nCells

       if (latCell(iCell) >= latitudeLimitNorth * rad_to_deg .or. &
           latCell(iCell) <= latitudeLimitSouth * rad_to_deg) then

          ! has ice
          iceAreaCell(iCell)    = iceAreaInitial
          iceVolumeCell(iCell)  = iceVolumeInitial
          snowVolumeCell(iCell) = snowVolumeInitial
          
       else

          ! no ice
          iceAreaCell(iCell)    = 0.0_RKIND
          iceVolumeCell(iCell)  = 0.0_RKIND
          snowVolumeCell(iCell) = 0.0_RKIND

       endif

    enddo ! iCell

  end subroutine init_ice_state_uniform_ice

  !--------------------------------------------------------------------------

  subroutine init_ice_state_circle_of_ice(&
       mesh, &
       iceAreaCell, &
       iceVolumeCell, &
       snowVolumeCell, &
       iceAreaInitial, &
       iceVolumeInitial, &
       snowVolumeInitial)

    type (MPAS_pool_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    real(kind=RKIND), intent(in) :: &
         iceAreaInitial, &
         iceVolumeInitial, &
         snowVolumeInitial

    integer :: &
         iCell

    real(kind=RKIND) :: &
         x, y, z
    
    real(kind=RKIND), parameter :: &
         circle_radius = 1.0e6_RKIND

    integer, pointer :: &
         nCells
    
    real(kind=RKIND), dimension(:), pointer :: &
         xCell, &
         yCell, &
         zCell
   
    ! init variables
    call MPAS_pool_get_array(mesh, "nCells", nCells)
    call MPAS_pool_get_array(mesh, "xCell", xCell)
    call MPAS_pool_get_array(mesh, "yCell", yCell)
    call MPAS_pool_get_array(mesh, "zCell", zCell)

    do iCell = 1, nCells
       
       x = xCell(iCell)
       y = yCell(iCell)
       z = zCell(iCell)
       
       if (sqrt(x**2+z**2) < circle_radius .and. y > 0.0_RKIND) then
          
          ! has ice
          iceAreaCell(iCell)    = iceAreaInitial
          iceVolumeCell(iCell)  = iceVolumeInitial
          snowVolumeCell(iCell) = snowVolumeInitial
          
       else
          
          ! no ice
          iceAreaCell(iCell)    = 0.0_RKIND
          iceVolumeCell(iCell)  = 0.0_RKIND
          snowVolumeCell(iCell) = 0.0_RKIND
          
       endif
       
    enddo ! iCell
    
  end subroutine init_ice_state_circle_of_ice

  !--------------------------------------------------------------------------

  subroutine init_ice_velocity_uniform(&
       mesh, &
       uVelocity, &
       vVelocity, &
       InitialUVelocity, &
       InitialVVelocity)
    
    type (MPAS_pool_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity
    
    real(kind=RKIND), intent(in) :: &
         InitialUVelocity, &
         InitialVVelocity
    
    integer :: &
         iVertex

    integer, pointer :: &
         nVertices

    call MPAS_pool_get_array(mesh, "nVertices", nVertices)

    do iVertex = 1, nVertices

       uVelocity(iVertex) = InitialUVelocity
       vVelocity(iVertex) = InitialVVelocity

    enddo ! iVertex

  end subroutine init_ice_velocity_uniform

  !--------------------------------------------------------------------------

end module cice_initialize
