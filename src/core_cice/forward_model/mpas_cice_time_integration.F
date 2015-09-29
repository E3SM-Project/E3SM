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

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use cice_debug

  implicit none

  private
  save

  public :: &
       cice_timestep

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
       domain, &
       clock, &
       dt, &
       itimestep)!{{{

    use mpas_timer

    use cice_unit_test, only: &
         cice_perform_unit_test

    use cice_velocity_solver, only: &
         cice_run_velocity_solver

    use cice_advection, only: &
         cice_run_advection

    use cice_therm_zerolayer, only: &
         cice_zerolayer_thermodynamics

    use cice_diagnostics, only: &
         cice_run_diagnostics

    use cice_debug, only: &
         cice_prepare_output_variables, &
         cice_test_case_output, &
         cice_debug_write_out_fields

    use cice_column, only: &
         cice_column_predynamics_time_integration, &
         cice_column_dynamics_time_integration, &
         cice_column_postdynamics_time_integration, &
         cice_column_reinitialize_fluxes

    use cice_history, only: &
         cice_init_history_thermodynamics, &
!nj-bgc+
         cice_init_history_biogeochemistry
!nj-bgc-

    type(domain_type), intent(inout) :: &
         domain !< Input/Output: 

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input: 

    real(kind=RKIND), intent(in) :: &
         dt !< Input: 

    integer, intent(in) :: &
         itimestep !< Input: 

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         configs

    logical, pointer :: &
         config_perform_unit_test, &
         config_use_advection, &
         config_use_vertical_thermodynamics

    call cice_debug_write_out_fields(domain, "cice_timestep", 1)

    call store_timestep_counter(domain, itimestep)

    ! init history variables
    call cice_init_history_thermodynamics(domain)

!nj-bgc+
    ! init biogeochemistry history variables 
    call cice_init_history_biogeochemistry(domain)
!nj-bgc-

    configs => domain % configs

    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_use_vertical_thermodynamics", config_use_vertical_thermodynamics)
    call MPAS_pool_get_config(configs, "config_perform_unit_test", config_perform_unit_test)

    if (config_perform_unit_test) then
       call cice_perform_unit_test(domain)
       return
    endif

    ! pre dynamics column physics
    call mpas_timer_start("Column")
    call cice_column_predynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column")

    ! zero layer thermodynamics
    !call mpas_timer_start("thermodynamics")
    !if (config_use_vertical_thermodynamics) then

    !   block => domain % blocklist
    !   do while (associated(block))
    !      call cice_zerolayer_thermodynamics(block, clock)
    !      block => block % next
    !   end do
       
    !endif
    !call mpas_timer_stop("thermodynamics")

    ! velocity solve
    call mpas_timer_start("Velocity solve")
    call cice_run_velocity_solver(domain, clock, dt)
    call mpas_timer_stop("Velocity solve")

    ! advection
    call mpas_timer_start("Advection")
    if (config_use_advection) &
       call cice_run_advection(domain, clock, dt)
    call mpas_timer_stop("Advection")

    ! ridging
    call mpas_timer_start("Column")
    call cice_column_dynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column")

    ! shortwave
    call mpas_timer_start("Column")
    call cice_column_postdynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column")

    ! prepare output variables
    call mpas_timer_start("prepare out")

    block => domain % blocklist
    do while (associated(block))
       call cice_prepare_output_variables(block)
       block => block % next
    end do

    call mpas_timer_stop("prepare out")

    ! test case output
    block => domain % blocklist
    do while (associated(block))
       call cice_test_case_output(block, itimestep)
       block => block % next
    end do

    ! diagnostics
    block => domain % blocklist
    do while (associated(block))
       call cice_run_diagnostics(block, clock)
       block => block % next
    end do

    call flush (stderrUnit)

    call cice_debug_write_out_fields(domain, "cice_timestep", 2)

    ! reinitialize fluxes
    call cice_column_reinitialize_fluxes(domain)


  end subroutine cice_timestep!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  store_timestep_counter
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine store_timestep_counter(domain, itimestep)

    type(domain_type), intent(in) :: &
         domain

    integer, intent(in) :: &
         itimestep

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh

    integer, dimension(:), pointer :: &
         mesh_itimestep

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       
       call MPAS_pool_get_array(mesh, "itimestep", mesh_itimestep)    
       
       mesh_itimestep(1) = itimestep

       block => block % next
    end do
    
  end subroutine store_timestep_counter

!-----------------------------------------------------------------------

end module cice_time_integration
