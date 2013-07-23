module cice_time_integration

  use mpas_grid_types
  use mpas_timekeeping

contains

  !--------------------------------------------------------------------------

  subroutine cice_timestep(block, clock, dt)

    use cice_dynamics, only: run_dynamics
    use cice_forcing, only: get_forcing

    use cice_testing, only: plot_cells_spherical, plot_pop_grid

    use cice_dyn_hex_weak, only: spherical_test_strain_weak, &
                                 spherical_test_divergence_stress_weak, &
                                 spherical_test_divergence_stress_weak_polygon

    implicit none

    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock
    real(kind=RKIND), intent(in) :: dt

    integer :: i, j

    call spherical_test_strain_weak(&
         block % mesh,    &
         block % hexdyn,  &
         block % hexweak, &
         block % normal)

    !call spherical_test_divergence_stress_weak(&
    !     block % mesh,    &
    !     block % hexdyn,  &
    !     block % hexweak, &
    !     block % normal)

    stop

    !call run_dynamics(block, dt)

    !call get_forcing(block % mesh, clock, block % atmos_forcing)

    !call plot_cells_spherical(block % mesh, block % atmos_forcing % airTemperature % array(1:block % mesh % nCells))

    !stop

    !call plot_pop_grid(block % mesh, block % atmos_forcing % airTemperature % array(1:block % mesh % nCells)) ; stop

  end subroutine cice_timestep

  !--------------------------------------------------------------------------

end module cice_time_integration
