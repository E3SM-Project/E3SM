module Simulation_Subsurface_class
  
  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Subsurface_class
  use PMC_Third_Party_class
  use PMC_Base_class
  use Realization_Subsurface_class
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"
  
  private

  type, public, extends(simulation_base_type) :: simulation_subsurface_type
    ! pointer to flow process model coupler
    class(pmc_subsurface_type), pointer :: flow_process_model_coupler
    ! pointer to reactive transport process model coupler
    class(pmc_subsurface_type), pointer :: rt_process_model_coupler
    ! pointer to realization object shared by flow and reactive transport
    class(realization_subsurface_type), pointer :: realization 
    ! regression object
    type(regression_type), pointer :: regression
    type(waypoint_list_type), pointer :: waypoint_list_subsurface
  contains
    procedure, public :: Init => SubsurfaceSimulationInit
    procedure, public :: JumpStart => SubsurfaceSimulationJumpStart
    procedure, public :: FinalizeRun => SubsurfaceFinalizeRun
    procedure, public :: Strip => SubsurfaceSimulationStrip
  end type simulation_subsurface_type
  
  public :: SubsurfaceSimulationCreate, &
            SubsurfaceSimulationInit, &
            SubsurfaceFinalizeRun, &
            SubsurfaceSimulationStrip, &
            SubsurfaceSimulationDestroy
  
contains

! ************************************************************************** !

function SubsurfaceSimulationCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(simulation_subsurface_type), pointer :: SubsurfaceSimulationCreate
  
#ifdef DEBUG
  print *, 'SimulationCreate'
#endif
  
  allocate(SubsurfaceSimulationCreate)
  call SubsurfaceSimulationCreate%Init(option)
  
end function SubsurfaceSimulationCreate

! ************************************************************************** !

subroutine SubsurfaceSimulationInit(this,option)
  ! 
  ! Initializes simulation values
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/22/13
  ! 
  use Waypoint_module
  use Option_module
  
  implicit none
  
  class(simulation_subsurface_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%flow_process_model_coupler)
  nullify(this%rt_process_model_coupler)
  nullify(this%realization)
  nullify(this%regression)
  this%waypoint_list_subsurface => WaypointListCreate()
  
end subroutine SubsurfaceSimulationInit

! ************************************************************************** !

subroutine SubsurfaceSimulationJumpStart(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  ! 
  use Logging_module
  use Output_module
  use Option_module
  use Output_Aux_module
  use Timestepper_Base_class

  implicit none
  
  class(simulation_subsurface_type) :: this

  class(timestepper_base_type), pointer :: master_timestepper
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_base_type), pointer :: tran_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  
#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceSimulationJumpStart()')
#endif

  nullify(master_timestepper)
  nullify(flow_timestepper)
  nullify(tran_timestepper)
  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE

  option => this%option
  output_option => this%output_option

  ! first time stepper is master
  master_timestepper => this%process_model_coupler_list%timestepper
  if (associated(this%flow_process_model_coupler)) then
    flow_timestepper => this%flow_process_model_coupler%timestepper
  endif
  if (associated(this%rt_process_model_coupler)) then
    tran_timestepper => this%rt_process_model_coupler%timestepper
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_timestepper%max_time_step < 0) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'')
    option%status = DONE
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(option,master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0) then
    if (output_option%print_initial_snap) snapshot_plot_flag = PETSC_TRUE
    if (output_option%print_initial_obs) observation_plot_flag = PETSC_TRUE
    if (output_option%print_initial_massbal) massbal_plot_flag = PETSC_FALSE
    call Output(this%realization,snapshot_plot_flag,observation_plot_flag, &
                massbal_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call printMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call printMsg(option)
    call printMsg(option,'') 
    option%status = DONE
    return
  endif

  ! increment plot number so that 000 is always the initial condition, and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_timestepper)) then
    if (.not.associated(flow_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      flow_timestepper%dt_max = flow_timestepper%cur_waypoint%dt_max
    endif
  endif  
  if (associated(tran_timestepper)) then
    if (.not.associated(tran_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start ' // &
        'time or simulation time needs to be extended on a restart.'
      call printMsg(option)
      option%status = FAIL
      return
    else
      tran_timestepper%dt_max = tran_timestepper%cur_waypoint%dt_max
    endif
  endif
           
  if (associated(flow_timestepper)) &
    flow_timestepper%start_time_step = flow_timestepper%steps + 1
  if (associated(tran_timestepper)) &
    tran_timestepper%start_time_step = tran_timestepper%steps + 1
  
end subroutine SubsurfaceSimulationJumpStart

! ************************************************************************** !

subroutine SubsurfaceFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_BE_class
  use Reaction_Sandbox_module, only : RSandboxDestroy
  use SrcSink_Sandbox_module, only : SSSandboxDestroyList
  use CLM_Rxn_module, only : RCLMRxnDestroy

  implicit none
  
  class(simulation_subsurface_type) :: this
  
  PetscErrorCode :: ierr
  
  class(timestepper_BE_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper

#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceFinalizeRun()')
#endif
  
  call SimulationBaseFinalizeRun(this)
  
  nullify(flow_timestepper)
  nullify(tran_timestepper)
  if (associated(this%flow_process_model_coupler)) then
    select type(ts => this%flow_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        flow_timestepper => ts
    end select
    call SSSandboxDestroyList()
  endif
  if (associated(this%rt_process_model_coupler)) then
    select type(ts => this%rt_process_model_coupler%timestepper)
      class is(timestepper_BE_type)
        tran_timestepper => ts
    end select
    call RSandboxDestroy()
    call RCLMRxnDestroy()
  endif
  
  call RegressionOutput(this%regression,this%realization, &
                        flow_timestepper,tran_timestepper)  
  
end subroutine SubsurfaceFinalizeRun

! ************************************************************************** !

subroutine SubsurfaceSimulationStrip(this)
  ! 
  ! Deallocates members of subsurface simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  implicit none
  
  class(simulation_subsurface_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SubsurfaceSimulationStrip()')
#endif
  
  call SimulationBaseStrip(this)
  call RealizationStrip(this%realization)
  deallocate(this%realization)
  nullify(this%realization)
  call RegressionDestroy(this%regression)
  call WaypointListDestroy(this%waypoint_list_subsurface)
  
end subroutine SubsurfaceSimulationStrip

! ************************************************************************** !

subroutine SubsurfaceSimulationDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  class(simulation_subsurface_type), pointer :: simulation
  
#ifdef DEBUG
  call printMsg(simulation%option,'SimulationDestroy()')
#endif
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SubsurfaceSimulationDestroy
  
end module Simulation_Subsurface_class
