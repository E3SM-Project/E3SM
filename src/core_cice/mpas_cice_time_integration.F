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

    use mpas_timer

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
         cice_run_diagnostics

    use cice_debug, only: &
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

    type (MPAS_pool_type), pointer :: &
         configs

    logical, pointer :: &
         config_perform_unit_test, &
         config_use_forcing, &
         config_use_velocity_solver, &
         config_use_advection, &
         config_use_vertical_thermodynamics

#define debug 0

    call store_timestep_counter(block, itimestep)

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
    call mpas_timer_start("forcing")
    if (config_use_forcing) &
         call cice_get_forcing(block, clock)
#if debug
    call check_solution(block, "post_force")
#endif
    call mpas_timer_stop("forcing")

    ! zero layer thermodynamics
    call mpas_timer_start("thermodynamics")
    if (config_use_vertical_thermodynamics) &
         call cice_zerolayer_thermodynamics(block, clock)
#if debug
    call check_solution(block, "post_therm")
#endif
    call mpas_timer_stop("thermodynamics")

    ! velocity solve
    call mpas_timer_start("velocity solve")
    if (config_use_velocity_solver) &
         call cice_run_velocity_solver(block, clock, dt)
#if debug
    call check_solution(block, "post_veloc")
#endif
    call mpas_timer_stop("velocity solve")

    ! advection
    call mpas_timer_start("advection")
    if (config_use_advection) &    
         call cice_run_advection(block, clock, dt)
#if debug
    call check_solution(block, "post_advec")
#endif
    call mpas_timer_stop("advection")

    ! prepare output variables
    call mpas_timer_start("prepare out")
    call cice_prepare_output_variables(block)
    call mpas_timer_stop("prepare out")

    ! test case output
    call cice_test_case_output(block, itimestep)

    ! diagnostics
    call cice_run_diagnostics(block, clock)

  end subroutine cice_timestep!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  check_solution
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine check_solution(block, message)

    use cice_constants, only: &
         iceAreaMinimum

    use cice_debug, only: &
         cice_plot_pop_grid_cells, &
         cice_plot_pop_grid_vertices

    type(block_type), intent(inout) :: &
         block !< Input/Output: 

    character(len=*), intent(in) :: &
         message

    type(MPAS_pool_type), pointer :: &
         mesh, &
         boundary, &
         tracers, &
         velocity_solver

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), allocatable :: &
         iceThickness

    integer, dimension(:), pointer :: &
         interiorVertex

    integer, dimension(:), pointer :: &
         itimestep

    character(len=200) :: filename

    integer :: &
         iCell

    integer, pointer :: &
         nCells

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "boundary", boundary)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocity_solver)

    call MPAS_pool_get_array(tracers, "iceAreaCategory", iceAreaCategory, 1)
    call MPAS_pool_get_array(tracers, "iceVolumeCategory", iceVolumeCategory, 1)

    call MPAS_pool_get_array(velocity_solver, "uVelocity", uVelocity)
    call MPAS_pool_get_array(velocity_solver, "vVelocity", vVelocity)

    call MPAS_pool_get_array(boundary, "interiorVertex", interiorVertex)     

    call MPAS_pool_get_array(mesh, "itimestep", itimestep)     
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)    

    write(*,*) itimestep(1), trim(message), &
         minval(iceAreaCategory), maxval(iceAreaCategory), &
         minval(iceVolumeCategory), maxval(iceVolumeCategory), &
         minval(uVelocity), maxval(uVelocity), &
         minval(vVelocity), maxval(vVelocity)
    write(*,*) itimestep(1), trim(message), &
         minloc(iceAreaCategory), maxloc(iceAreaCategory), &
         minloc(iceVolumeCategory), maxloc(iceVolumeCategory), &
         minloc(uVelocity), maxloc(uVelocity), &
         minloc(vVelocity), maxloc(vVelocity)

    write(filename,fmt='(a,a)') trim(message), "_iceAreaCategory.plt"
    call cice_plot_pop_grid_cells(mesh, iceAreaCategory(1,1,:), trim(filename), itimestep(1))

    write(filename,fmt='(a,a)') trim(message), "_iceVolumeCategory.plt"
    call cice_plot_pop_grid_cells(mesh, iceVolumeCategory(1,1,:), trim(filename), itimestep(1))

    write(filename,fmt='(a,a)') trim(message), "_uVelocity.plt"
    call cice_plot_pop_grid_vertices(mesh, uVelocity(:), interiorVertex, trim(filename), itimestep(1))

    write(filename,fmt='(a,a)') trim(message), "_vVelocity.plt"
    call cice_plot_pop_grid_vertices(mesh, vVelocity(:), interiorVertex, trim(filename), itimestep(1))

    write(filename, fmt='(a,a,i5.5,a)') "history/"//trim(message), "_iceAreaCategory_",itimestep(1),".txt"
    open(11,file=trim(filename))
    do iCell = 1, nCells
       write(11,*) iCell, iceAreaCategory(1,1,iCell)
    enddo ! iCell
    close(11)

    write(filename, fmt='(a,a,i5.5,a)') "history/"//trim(message), "_iceVolumeCategory_",itimestep(1),".txt"
    open(11,file=trim(filename))
    do iCell = 1, nCells
       write(11,*) iCell, iceVolumeCategory(1,1,iCell)
    enddo ! iCell
    close(11)

    open(11,file="history/concerns.txt",position='append')
    write(11,*) trim(message), itimestep(1)
    do iCell = 1, nCells
       if (iceVolumeCategory(1,1,iCell) < 0.0_RKIND .or. iceAreaCategory(1,1,iCell) < 0.0_RKIND) then
          write(11,*) iCell, iceAreaCategory(1,1,iCell), iceVolumeCategory(1,1,iCell)
       endif
       if (iceAreaCategory(1,1,iCell) > iceAreaMinimum .and. iceVolumeCategory(1,1,iCell) == 0.0_RKIND) then
          write(11,*) iCell, iceAreaCategory(1,1,iCell), iceVolumeCategory(1,1,iCell)
       endif
    enddo ! iCell
    close(11)

    allocate(iceThickness(nCells))

    write(filename, fmt='(a,a,i5.5,a)') "history/"//trim(message), "_iceThickness_",itimestep(1),".txt"
    open(11,file=trim(filename))

    do iCell = 1, nCells

       if (iceAreaCategory(1,1,iCell) > iceAreaMinimum) then
          iceThickness(iCell) = iceVolumeCategory(1,1,iCell) / iceAreaCategory(1,1,iCell)
       else
          iceThickness(iCell) = 0.0_RKIND
       endif

       write(11,*) iCell, iceThickness(iCell)

    enddo ! iCell

    close(11)

    write(filename,fmt='(a,a)') trim(message), "_iceThickness.plt"
    call cice_plot_pop_grid_cells(mesh, iceThickness(:), trim(filename), itimestep(1))

    deallocate(iceThickness)

  end subroutine check_solution

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

  subroutine store_timestep_counter(block, itimestep)

    type(block_type), intent(in) :: &
         block

    integer, intent(in) :: &
         itimestep

    type(MPAS_pool_type), pointer :: &
         mesh

    integer, dimension(:), pointer :: &
         mesh_itimestep

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)

    call MPAS_pool_get_array(mesh, "itimestep", mesh_itimestep)    

    mesh_itimestep(1) = itimestep

  end subroutine store_timestep_counter

!-----------------------------------------------------------------------

end module cice_time_integration
