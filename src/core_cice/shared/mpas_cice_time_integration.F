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

    use cice_column, only: &
         cice_column_predynamics_time_integration, &
         cice_column_dynamics_time_integration, &
         cice_column_postdynamics_time_integration, &
         cice_column_reinitialize_fluxes, &
         cice_column_reinitialize_diagnostics

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
         config_use_advection

    integer, pointer :: &
         config_dynamics_subcycle_number

    integer :: &
         iDynamicsSubcycle

    ! set halo timer
    call mpas_timer_start("Halo")
    call mpas_timer_stop("Halo")

    call store_timestep_counter(domain, itimestep)

    ! reinitialize diagnostics
    call cice_column_reinitialize_diagnostics(domain)

    configs => domain % configs

    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_perform_unit_test", config_perform_unit_test)

    if (config_perform_unit_test) then
       call cice_perform_unit_test(domain)
       return
    endif

    ! pre dynamics column physics
    call mpas_timer_start("Column")
    call cice_column_predynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column")

    ! loop of dynamcis subcycle
    call MPAS_pool_get_config(configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
    do iDynamicsSubcycle = 1, config_dynamics_subcycle_number

       ! velocity solve
       call mpas_timer_start("Velocity solver")
       call cice_run_velocity_solver(domain, clock)
       call mpas_timer_stop("Velocity solver")

       ! advection
       call mpas_timer_start("Advection")
       if (config_use_advection) &
            call cice_run_advection(domain, clock, dt)
       call mpas_timer_stop("Advection")

       ! ridging
       call mpas_timer_start("Column")
       call cice_column_dynamics_time_integration(domain, clock)
       call mpas_timer_stop("Column")

    enddo ! iDynamicsSubcycle

    ! shortwave
    call mpas_timer_start("Column")
    call cice_column_postdynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column")

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
