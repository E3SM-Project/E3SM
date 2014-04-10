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

    call init_junk_values(block % mesh)

    ! calculate boundary variables
    call init_boundary(block)

    ! init dynamics
    if (config_use_velocity_solver) &
         call init_dynamics(block, dt)
    
    ! init advection
    !if (config_use_advection) &
    !     call init_advection(block % mesh)
    
    ! init forcing
    if (config_use_forcing) then

       call init_forcing(&
            block % mesh, &
            block % atmos_forcing)
       
       call get_forcing(&
            block % mesh, &
            clock, &
            block % atmos_forcing)
       
    endif

    ! init ice state
    call init_ice_state(&
         block % mesh, &
         block % icestate, &
         block % hexdyn, &
         block % boundary, &
         block % atmos_forcing)
    
    if (config_perform_unit_test) &
         call perform_unit_test(block, config_unit_test_type, config_unit_test_subtype)
    
    ! write out the initial condition
    !call writeout_state(block, 0)

  end subroutine cice_init

  !--------------------------------------------------------------------------

  subroutine init_junk_values(mesh)

    type(mesh_type), intent(inout) :: mesh

    integer :: &
         njunk

    real(kind=RKIND), parameter :: &
         value_junk = -1.0e34

    njunk = mesh % nCells + 1
    
    mesh % areaCell % array(njunk) = value_junk

  end subroutine init_junk_values

  !--------------------------------------------------------------------------
  
  subroutine init_ice_state(&
       mesh, &
       icestate, &
       hexdyn, &
       boundary, &
       atmos_forcing)

    use mpas_configure

    type (mesh_type),          intent(in) :: mesh
    type (icestate_type),      pointer :: icestate
    type (hexdyn_type),        pointer :: hexdyn
    type (boundary_type),      pointer :: boundary
    type (atmos_forcing_type), pointer :: atmos_forcing

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         tracers                  ! (nTracers, nVertLevels, nCells)

    integer :: &
         index_iceAreaCell, &
         index_iceVolumeCell, &
         index_snowVolumeCell

    tracers => icestate % tracers % array

    index_iceAreaCell    = icestate % index_iceAreaCell
    index_iceVolumeCell  = icestate % index_iceVolumeCell
    index_snowVolumeCell = icestate % index_snowVolumeCell

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
            hexdyn % uVelocity % array, &
            hexdyn % vVelocity % array, &
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

    type (mesh_type), intent(in) :: mesh
    
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

    do iCell = 1, mesh % nCells

       if (mesh % latCell % array(iCell) >= latitudeLimitNorth * rad_to_deg .or. &
           mesh % latCell % array(iCell) <= latitudeLimitSouth * rad_to_deg) then

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

    type (mesh_type), intent(in) :: mesh
    
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

    do iCell = 1, mesh % nCells
       
       x = mesh % xCell % array(iCell)
       y = mesh % yCell % array(iCell)
       z = mesh % zCell % array(iCell)
       
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
    
    type (mesh_type), intent(in) :: mesh
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity
    
    real(kind=RKIND), intent(in) :: &
         InitialUVelocity, &
         InitialVVelocity
    
    integer :: &
         iVertex

    do iVertex = 1, mesh % nVertices

       uVelocity(iVertex) = InitialUVelocity
       vVelocity(iVertex) = InitialVVelocity

    enddo ! iVertex

  end subroutine init_ice_velocity_uniform

  !--------------------------------------------------------------------------

end module cice_initialize
