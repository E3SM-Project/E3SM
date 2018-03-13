!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_time_integration
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_time_integration

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping

  implicit none

  private
  save

  public :: &
       seaice_timestep, &
       seaice_timestep_finalize

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_timestep
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_timestep(&
       domain, &
       clock, &
       itimestep, &
       ierr)!{{{

    use mpas_timer

    use seaice_diagnostics, only: &
         seaice_set_time_diagnostics, &
         seaice_check_state, &
         seaice_load_balance_timers

    use seaice_unit_test, only: &
         seaice_perform_unit_test

    use seaice_velocity_solver, only: &
         seaice_run_velocity_solver

    use seaice_advection, only: &
         seaice_run_advection

    use seaice_column, only: &
         seaice_column_predynamics_time_integration, &
         seaice_column_dynamics_time_integration, &
         seaice_column_postdynamics_time_integration, &
         seaice_column_reinitialize_fluxes, &
         seaice_column_reinitialize_diagnostics_thermodynamics, &
         seaice_column_reinitialize_diagnostics_bgc, &
         seaice_column_reinitialize_diagnostics_dynamics

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type (MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, intent(in) :: &
         itimestep !< Input:

    integer, intent(inout) :: &
         ierr !< Input

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

#ifndef MPAS_PERF_MOD_TIMERS
    ! set halo timer
    call mpas_timer_start("Halo")
    call mpas_timer_stop("Halo")
#endif

    ! halo barrier diagnostics for coupling/forcing
    call seaice_load_balance_timers(domain, "coupling")

    ! set time diagnostics
    call seaice_set_time_diagnostics(domain)

    call store_timestep_counter(domain, itimestep)

    ! reinitialize diagnostics
    call mpas_timer_start("Reinitialize diagnostics thermodynamics/bgc")
    call seaice_column_reinitialize_diagnostics_thermodynamics(domain)
    call seaice_column_reinitialize_diagnostics_bgc(domain)
    call mpas_timer_stop("Reinitialize diagnostics thermodynamics/bgc")

    configs => domain % configs

    call MPAS_pool_get_config(configs, "config_use_advection", config_use_advection)
    call MPAS_pool_get_config(configs, "config_perform_unit_test", config_perform_unit_test)

    if (config_perform_unit_test) then
       call seaice_perform_unit_test(domain)
       return
    endif

    ! pre dynamics column physics
    call mpas_timer_start("Column pre-dynamics")
    call seaice_column_predynamics_time_integration(domain, clock, ierr)
    if (ierr > 0) return
    call mpas_timer_stop("Column pre-dynamics")

    ! loop of dynamcis subcycle
    call MPAS_pool_get_config(configs, "config_dynamics_subcycle_number", config_dynamics_subcycle_number)
    do iDynamicsSubcycle = 1, config_dynamics_subcycle_number

       ! reinitialize dynamics diagnostics
       call mpas_timer_start("Reinitialize diagnostics dynamics")
       call seaice_column_reinitialize_diagnostics_dynamics(domain)
       call mpas_timer_stop("Reinitialize diagnostics dynamics")

       ! velocity solve
       call mpas_timer_start("Velocity solver")
       call seaice_run_velocity_solver(domain, clock)
       call mpas_timer_stop("Velocity solver")

       ! advection
       call mpas_timer_start("Advection")
       if (config_use_advection) &
            call seaice_run_advection(domain, clock, ierr)
       call mpas_timer_stop("Advection")
       if (ierr > 0) return

       ! ridging
       call mpas_timer_start("Column")
       call seaice_column_dynamics_time_integration(domain, clock, ierr)
       call mpas_timer_stop("Column")
       if (ierr > 0) return

    enddo ! iDynamicsSubcycle

    ! shortwave
    call mpas_timer_start("Column post-dynamics")
    call seaice_column_postdynamics_time_integration(domain, clock)
    call mpas_timer_stop("Column post-dynamics")

    ! check the physical state of the model
    call seaice_check_state(domain)

  end subroutine seaice_timestep!}}}

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

    use seaice_constants, only: &
         seaiceSecondsPerDay

    type(domain_type), intent(in) :: &
         domain

    integer, intent(in) :: &
         itimestep

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         diagnosticsPool

    integer, dimension(:), pointer :: &
         mesh_itimestep

    real(kind=RKIND), pointer :: &
         timeAverageTestVariable, &
         config_dt

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)
       call MPAS_pool_get_array(mesh, "itimestep", mesh_itimestep)
       call MPAS_pool_get_array(diagnosticsPool, "timeAverageTestVariable", timeAverageTestVariable)

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)

       mesh_itimestep(1) = itimestep

       timeAverageTestVariable = timeAverageTestVariable + config_dt / seaiceSecondsPerDay

       block => block % next
    end do

  end subroutine store_timestep_counter

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_timestep_finalize
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 27th January 2016
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_timestep_finalize(&
       domain)!{{{

    use seaice_column, only: &
         seaice_column_reinitialize_fluxes

    type(domain_type), intent(in) :: &
         domain

    call seaice_column_reinitialize_fluxes(domain)

  end subroutine seaice_timestep_finalize

!-----------------------------------------------------------------------

end module seaice_time_integration
