!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_time_integration
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_time_integration

  use mpas_grid_types
  use mpas_timekeeping

  implicit none

  private
  save

  public :: &
       cice_timestep

  integer, private :: nstep = 0

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_timestep
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_timestep(&
       block, &
       clock, &
       dt, &
       itimestep)!{{{

    use mpas_configure

    use cice_unit_test, only: &
         cice_perform_unit_test

    use cice_velocity_solver, only: &
         cice_run_velocity_solver

    use cice_advection, only: &
         cice_run_advection

    use cice_therm_zerolayer, only: &
         cice_zerolayer_thermodynamics

    use cice_forcing, only: &
         cice_get_forcing

    use cice_diagnostics, only: &
         cice_prepare_output_variables, &
         cice_test_case_output

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    integer, intent(in) :: &
         itimestep !< Input: 

    real(kind=RKIND) :: &
         time

    type (MPAS_pool_type), pointer :: &
         configs

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

    if (config_perform_unit_test) then
       call cice_perform_unit_test(block, configs)
       return
    endif

    ! get forcing
    if (config_use_forcing) &
         call cice_get_forcing(block, clock)

    ! zero layer thermodynamics
    if (config_use_vertical_thermodynamics) &
         call cice_zerolayer_thermodynamics(block)

    ! velocity solve
    if (config_use_velocity_solver) &
         call cice_run_velocity_solver(block, dt, nstep)

    ! advection
    if (config_use_advection) &    
         call cice_run_advection(block, dt)

    ! prepare output variables
    call cice_prepare_output_variables(block)

    ! test case output
    call cice_test_case_output(block, itimestep)

  end subroutine cice_timestep!}}}

!-----------------------------------------------------------------------

end module cice_time_integration
