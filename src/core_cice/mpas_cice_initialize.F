module cice_initialize

  use mpas_grid_types
  use mpas_timekeeping

  private
  public :: &
       cice_init

contains

  !--------------------------------------------------------------------------

  subroutine cice_init(block, clock, dt)

    use cice_dynamics_shared, only: interior_vertices
    use cice_advection, only: init_advection
    use cice_dynamics, only: init_dynamics
    use cice_testing, only: init_spherical_test_case
    use cice_diagnostics, only: plot_pop_grid_cells, output_meridonal_strip
    use cice_forcing, only: init_forcing, get_forcing, forcing_file_comparison
    use cice_therm_zerolayer, only: init_semtner_forcing, init_semtner_test_case
    
    type(block_type), intent(inout) :: block
    type (MPAS_Clock_type), intent(in) :: clock

    real(kind=RKIND), intent(in) :: &
         dt

    !call forcing_file_comparison(block % mesh)
    !stop

    ! calculate which vertices are interior
    call interior_vertices(block % mesh, &
         block % boundary % interiorVertex % array)

    ! init dynamics
    call init_dynamics(block, dt)
    
    ! init advection
    call init_advection(block % mesh)
    
    ! init forcing
    !call init_forcing(block % mesh, &
    !     block % atmos_forcing)
    
    !call get_forcing(block % mesh, clock, block % atmos_forcing)

    !call plot_pop_grid_cells(block % mesh, block % atmos_forcing % airTemperature % array, "airTemperature.txt")
    
    
    call init_spherical_test_case(block % mesh, &
         block % icestate % tracers % array(block % icestate % index_iceAreaCell, 1, :), &
         block % icestate % tracers % array(block % icestate % index_iceVolumeCell, 1, :), &
         block % icestate % totalMassCell % array, &
         block % icestate % surfaceTemperature % array, &
         block % atmos_forcing % airTemperature % array, &
         block % hexdyn % uVelocity % array, &
         block % hexdyn % vVelocity % array, &
         block % boundary % interiorVertex % array)
    
    
    call output_meridonal_strip(block % mesh, block % icestate % tracers % array(block % icestate % index_iceAreaCell, 1, :), &
         "area_strip.txt", 10, 0, 1)
    
    !call plot_pop_grid_cells(block % mesh, block % icestate % tracers % array(block % icestate % index_iceAreaCell, 1, :), "iceAreaCell_ic.txt")
    !call plot_pop_grid_cells(block % mesh, block % icestate % tracers % array(block % icestate % index_iceVolumeCell, 1, :), "iceVolumeCell_ic.txt")


  end subroutine cice_init

  !--------------------------------------------------------------------------

end module cice_initialize
