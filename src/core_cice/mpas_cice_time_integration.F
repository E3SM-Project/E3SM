module cice_time_integration

  use mpas_grid_types
  use mpas_timekeeping

  integer :: nstep = 0

contains

  !--------------------------------------------------------------------------

  subroutine cice_timestep(block, clock, dt)

    use cice_dynamics, only: run_dynamics
    use cice_forcing, only: get_forcing

    use cice_testing, only: plot_cells_spherical, plot_pop_grid_cells, plot_state

    use cice_therm_zerolayer, only: therm_zerolayer

    implicit none

    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock
    real(kind=RKIND), intent(in) :: dt

    real(kind=RKIND) :: &
         time

    call get_forcing( &
         block % mesh, &
         clock, &
         block % atmos_forcing)

    nstep = nstep + 1
    time = nstep * dt

    call run_dynamics(block, dt, nstep)

  end subroutine cice_timestep

  !--------------------------------------------------------------------------

end module cice_time_integration
