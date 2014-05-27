module cice_time_integration

  use mpas_grid_types
  use mpas_timekeeping

  implicit none

  integer :: nstep = 0

contains

  !--------------------------------------------------------------------------

  subroutine cice_timestep(block, clock, dt)

    use mpas_configure

    use cice_unit_test, only: &
         perform_unit_test

    use cice_dynamics, only: &
         run_dynamics

    use cice_advection, only: &
         run_advection

    use cice_therm_zerolayer, only: &
         therm_zerolayer

    use cice_forcing, only: &
         get_forcing

    use cice_diagnostics, only:&
         writeout_state, &
         prepare_plot_variables

    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock
    real(kind=RKIND), intent(in) :: dt

    real(kind=RKIND) :: &
         time

    type (MPAS_pool_type), pointer :: &
         configs, &
         mesh, &
         atmos_forcing

    logical, pointer :: &
         config_perform_unit_test, &
         config_use_forcing, &
         config_use_velocity_solver, &
         config_use_advection, &
         config_use_vertical_thermodynamics
         
    nstep = nstep + 1
    time = nstep * dt

    configs => block % configs

    call MPAS_pool_get_config(configs, "config_use_forcing", config_use_forcing)
    call MPAS_pool_get_config(configs, "config_use_velocity_solver", config_use_velocity_solver)
    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_use_vertical_thermodynamics", config_use_vertical_thermodynamics)
    call MPAS_pool_get_config(configs, "config_perform_unit_test", config_perform_unit_test)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "atmos_forcing", atmos_forcing)

    if (config_perform_unit_test) then
       call perform_unit_test(block, configs)
       return
    endif

    ! get forcing
    if (config_use_forcing) &
         call get_forcing( &
            mesh, &
            clock, &
            atmos_forcing)

    ! velocity solve
    if (config_use_velocity_solver) &
         call run_dynamics(block, dt, nstep)

    ! advection
    if (config_use_advection) &    
         call run_advection(block, dt)

    ! zero layer thermodynamics
    if (config_use_vertical_thermodynamics) &
         call therm_zerolayer(block)

    ! output
    !call writeout_state(block, nstep)

    call prepare_plot_variables(block)

  end subroutine cice_timestep

  !--------------------------------------------------------------------------

end module cice_time_integration
