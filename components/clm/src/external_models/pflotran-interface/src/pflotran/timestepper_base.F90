module Timestepper_Base_class
 
#include "petsc/finclude/petscsys.h"
  use Waypoint_module 
 
  use PFLOTRAN_Constants_module

  implicit none

  private
  
 
  PetscInt, parameter, public :: TS_CONTINUE = 0
  PetscInt, parameter, public :: TS_STOP_END_SIMULATION = 1
  PetscInt, parameter, public :: TS_STOP_MAX_TIME_STEP = 2
  PetscInt, parameter, public :: TS_STOP_WALLCLOCK_EXCEEDED = 3
  PetscInt, parameter, public :: TS_STOP_FAILURE = 4

  type, public :: timestepper_base_type
  
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: steps         ! The number of time steps taken by the code.
    PetscInt :: num_constant_time_steps   ! number of contiguous time_steps of constant size

    PetscInt :: max_time_step                ! Maximum number of time steps to be taken by the code.
    PetscInt :: max_time_step_cuts           ! Maximum number of timestep cuts within one time step.
    PetscReal :: time_step_reduction_factor  ! Scaling factor by which timestep is reduced.
    PetscReal :: time_step_max_growth_factor ! Maximum scaling factor by which timestep is increasd.
    PetscInt :: constant_time_step_threshold ! Steps needed after cutting to increase time step

    PetscInt :: cumulative_time_step_cuts    ! Total number of cuts in the timestep taken.    
    PetscReal :: cumulative_solver_time
    
    PetscReal :: dt
    PetscReal :: prev_dt
    PetscReal :: dt_init
    PetscReal :: dt_min
    PetscReal :: dt_max
    PetscBool :: revert_dt
    PetscInt :: num_contig_revert_due_to_sync
    PetscInt :: max_num_contig_revert
    
    PetscBool :: init_to_steady_state
    PetscBool :: run_as_steady_state
    PetscReal :: steady_state_rel_tol
    
    PetscBool :: time_step_cut_flag  ! flag toggled if timestep is cut

    PetscLogDouble :: start_time    
    PetscInt :: start_time_step ! the first time step of a given run
    PetscReal :: time_step_tolerance ! scalar used in determining time step size
    PetscReal :: target_time    ! time at end of "synchronized" time step 
    PetscBool :: print_ekg

    type(waypoint_type), pointer :: cur_waypoint
    type(waypoint_type), pointer :: prev_waypoint

  contains
    
    procedure, public :: ReadInput => TimestepperBaseRead
    procedure, public :: Init => TimestepperBaseInit
    procedure, public :: InitializeRun => TimestepperBaseInitializeRun
    procedure, public :: SetTargetTime => TimestepperBaseSetTargetTime
    procedure, public :: StepDT => TimestepperBaseStepDT
    procedure, public :: UpdateDT => TimestepperBaseUpdateDT
    procedure, public :: CheckpointBinary => TimestepperBaseCheckpointBinary
    procedure, public :: CheckpointHDF5 => TimestepperBaseCheckpointHDF5
    procedure, public :: RestartBinary => TimestepperBaseRestartBinary
    procedure, public :: RestartHDF5 => TimestepperBaseRestartHDF5
    procedure, public :: Reset => TimestepperBaseReset
    procedure, public :: WallClockStop => TimestepperBaseWallClockStop
    procedure, public :: PrintInfo => TimestepperBasePrintInfo
    procedure, public :: InputRecord => TimestepperBaseInputRecord
    procedure, public :: FinalizeRun => TimestepperBaseFinalizeRun
    procedure, public :: Strip => TimestepperBaseStrip
    procedure, public :: Destroy => TimestepperBaseDestroy
    
  end type timestepper_base_type
  
  type, public :: stepper_base_header_type
    PetscReal :: time
    PetscReal :: dt
    PetscReal :: prev_dt
    PetscInt :: num_steps
    PetscInt :: cumulative_time_step_cuts
    PetscInt :: num_constant_time_steps
    PetscInt :: num_contig_revert_due_to_sync
    PetscInt :: revert_dt
  end type stepper_base_header_type
  
  public :: TimestepperBaseCreate, &
            TimestepperBaseProcessKeyword, &
            TimestepperBaseStrip, &
            TimestepperBaseInit, &
            TimestepperBaseSetHeader, &
            TimestepperBaseGetHeader, &
            TimestepperBaseReset, &
            TimestepperBaseRegisterHeader, &
            TimestepperBasePrintInfo, &
            TimestepperBaseInputRecord

contains

! ************************************************************************** !

function TimestepperBaseCreate()
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  class(timestepper_base_type), pointer :: TimestepperBaseCreate
  
  class(timestepper_base_type), pointer :: this
  
  allocate(this)
  call this%Init()
  
  TimestepperBaseCreate => this
  
end function TimestepperBaseCreate

! ************************************************************************** !

subroutine TimestepperBaseInit(this)
  ! 
  ! Allocates and initializes a new Timestepper object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/01/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  class(timestepper_base_type) :: this
  
  this%name = ''
  this%steps = 0
  this%num_constant_time_steps = 0

  this%max_time_step = 999999
  this%max_time_step_cuts = 16
  this%time_step_reduction_factor = 0.5d0
  this%time_step_max_growth_factor = 2.d0
  this%constant_time_step_threshold = 5

  this%cumulative_time_step_cuts = 0    
  this%cumulative_solver_time = 0.d0

  this%start_time = 0.d0  
  this%start_time_step = 0
  this%time_step_tolerance = 0.1d0
  this%target_time = 0.d0
  
  this%prev_dt = 0.d0
  this%dt = 1.d0
  this%dt_init = 1.d0
  this%dt_min = 1.d-20   ! Ten zeptoseconds.
  this%dt_max = 3.1536d6 ! One-tenth of a year.  
  
  this%time_step_cut_flag = PETSC_FALSE
  
  this%init_to_steady_state = PETSC_FALSE
  this%steady_state_rel_tol = 1.d-8
  this%run_as_steady_state = PETSC_FALSE
  
  nullify(this%cur_waypoint)
  nullify(this%prev_waypoint)
  this%revert_dt = PETSC_FALSE
  this%num_contig_revert_due_to_sync = 0
  this%max_num_contig_revert = 2
  this%print_ekg = PETSC_FALSE
  
end subroutine TimestepperBaseInit

! ************************************************************************** !

subroutine TimestepperBaseInitializeRun(this,option)
  ! 
  ! Initializes the timestepper for the simulation.  This is more than just
  ! initializing parameters.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/21/14
  ! 

  use Option_module
  
  implicit none

  class(timestepper_base_type) :: this
  type(option_type) :: option
  
  call this%PrintInfo(option)
  option%time = this%target_time
  ! For the case where the second waypoint is a printout after the first time 
  ! step, we must increment the waypoint beyond the first (time=0.) waypoint.  
  ! Otherwise the second time step will be zero. - geh
  if (this%cur_waypoint%time < 1.d-40) then
    this%cur_waypoint => this%cur_waypoint%next
  endif

end subroutine TimestepperBaseInitializeRun

! ************************************************************************** !

subroutine TimestepperBaseRead(this,input,option)
  ! 
  ! Reads parameters associated with time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none

  class(timestepper_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseRead not supported.  Requires extension.'
  call printErrMsg(option)

end subroutine TimestepperBaseRead

! ************************************************************************** !

subroutine TimestepperBaseProcessKeyword(this,input,option,keyword)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  
  implicit none
  
  class(timestepper_base_type) :: this
  character(len=MAXWORDLENGTH) :: keyword
  type(input_type) :: input
  type(option_type) :: option

  select case(trim(keyword))

    case('NUM_STEPS_AFTER_TS_CUT')
      call InputReadInt(input,option,this%constant_time_step_threshold)
      call InputErrorMsg(input,option,'num_constant_time_steps_after_ts_cut', &
                         'TIMESTEPPER')
    case('MAX_STEPS')
      call InputReadInt(input,option,this%max_time_step)
      call InputErrorMsg(input,option,'max_time_step','TIMESTEPPER')
    case('MAX_TS_CUTS')
      call InputReadInt(input,option,this%max_time_step_cuts)
      call InputErrorMsg(input,option,'max_time_step_cuts','TIMESTEPPER')
    case('MAX_NUM_CONTIGUOUS_REVERTS')
      call InputReadInt(input,option,this%max_num_contig_revert)
      call InputErrorMsg(input,option,'max_num_contig_reverts','TIMESTEPPER')
    case('TIMESTEP_REDUCTION_FACTOR')
      call InputReadDouble(input,option,this%time_step_reduction_factor)
      call InputErrorMsg(input,option,'timestep reduction factor','TIMESTEPPER')
    case('TIMESTEP_MAXIMUM_GROWTH_FACTOR')
      call InputReadDouble(input,option,this%time_step_max_growth_factor)
      call InputErrorMsg(input,option,'timestep maximum growth factor', &
                         'TIMESTEPPER')
    case('TIMESTEP_OVERSTEP_TOLERANCE')
      call InputReadDouble(input,option,this%time_step_tolerance)
      call InputErrorMsg(input,option,'timestep overstep tolerance', &
                         'TIMESTEPPER')
    case('INITIALIZE_TO_STEADY_STATE')
      this%init_to_steady_state = PETSC_TRUE
      call InputReadDouble(input,option,this%steady_state_rel_tol)
      call InputErrorMsg(input,option,'steady state convergence relative &
                         &tolerance','TIMESTEPPER')
    case('RUN_AS_STEADY_STATE')
      this%run_as_steady_state = PETSC_TRUE
    case('PRINT_EKG')
      this%print_ekg = PETSC_TRUE
      option%print_ekg = PETSC_TRUE
    case('MAX_PRESSURE_CHANGE','MAX_TEMPERATURE_CHANGE', &
         'MAX_CONCENTRATION_CHANGE','MAX_SATURATION_CHANGE', &
         'PRESSURE_DAMPENING_FACTOR','SATURATION_CHANGE_LIMIT', &
         'PRESSURE_CHANGE_LIMIT','TEMPERATURE_CHANGE_LIMIT')
      option%io_buffer = 'Keyword "' // trim(keyword) // '" has been &
        &deprecated in TIMESTEPPER and moved to the FLOW PM OPTIONS block.'
      call printErrMsg(option)
    case default
      call InputKeywordUnrecognized(keyword,'TIMESTEPPER',option)
  end select

end subroutine TimestepperBaseProcessKeyword

! ************************************************************************** !

subroutine TimestepperBaseUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  ! 

  use PM_Base_class
  use Option_module
  
  implicit none

  class(timestepper_base_type) :: this
  class(pm_base_type) :: process_model
  
  process_model%option%io_buffer = 'TimestepperBaseStepDT must be extended.'
  call printErrMsg(process_model%option)

end subroutine TimestepperBaseUpdateDT

! ************************************************************************** !

subroutine TimestepperBaseSetTargetTime(this,sync_time,option,stop_flag, &
                                        snapshot_plot_flag, &
                                        observation_plot_flag, &
                                        massbal_plot_flag,checkpoint_flag)
  ! 
  ! Sets target time for timestepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none

  class(timestepper_base_type) :: this
  PetscReal :: sync_time
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag
  
  PetscReal :: target_time
  PetscReal :: dt
  PetscReal :: dt_max
  PetscInt :: cumulative_time_steps
  PetscInt :: max_time_step
  PetscReal :: max_time
  PetscReal :: tolerance
  PetscBool :: force_to_match_waypoint
  PetscBool :: equal_to_or_exceeds_waypoint
  PetscBool :: equal_to_or_exceeds_sync_time
  PetscBool :: revert_due_to_waypoint
  PetscBool :: revert_due_to_sync_time
  PetscBool :: truncated_due_to_next_dt_max
  PetscReal :: temp_time
  type(waypoint_type), pointer :: cur_waypoint, next_waypoint, prev_waypoint

!geh: for debugging
#ifdef DEBUG
  option%io_buffer = 'TimestepperBaseSetTargetTime()'
  call printMsg(option)
#endif
  
  if (this%time_step_cut_flag) then
    this%time_step_cut_flag = PETSC_FALSE
    !geh: pointing the cur_waypoint back may cause problems in the checkpoint
    !     file.  There is no way of knowing whether prev_waypoint is different
    !     from cur_waypoint as most of the time it will be identical.  I believe
    !     the only way around this is to check associated(cur,prev) to see if
    !     they differ and set a flag in the checkpoint file.  But even this will
    !     not work if more than one waypoint previous.
    this%cur_waypoint => this%prev_waypoint
  else
    ! If the maximum time step size decreased in the past step, need to set
    ! the time step size to the minimum of the this%prev_dt and 
    ! this%dt_max.  However, if we have to revert "max_num_contig_revert" 
    ! times in a row, throw away the old time step and move on.
    if (this%revert_dt) then
      if (this%num_contig_revert_due_to_sync < this%max_num_contig_revert) then
        this%dt = min(this%prev_dt,this%dt_max)
      endif
    endif
  endif
  this%revert_dt = PETSC_FALSE ! reset back to false
  revert_due_to_waypoint = PETSC_FALSE
  revert_due_to_sync_time = PETSC_FALSE
  
  dt = this%dt
  this%prev_dt = dt
  cur_waypoint => this%cur_waypoint
  ! need previous waypoint for reverting back on time step cut
  this%prev_waypoint => this%cur_waypoint
  ! dt_max must be lagged.  it can be updated below, but it must lag a waypoint.
  cumulative_time_steps = this%steps
  max_time_step = this%max_time_step
  tolerance = this%time_step_tolerance
!  target_time = this%target_time + dt

  do ! we cycle just in case the next waypoint is beyond the target_time
    dt_max = cur_waypoint%dt_max
    dt = min(dt,dt_max)
    ! ensure that the time step does not overstep the next waypoint time + 
    ! dtmax combination.
    target_time = this%target_time + dt

!---! This section of code ensures that no time step over steps the next 
    ! maximum time step (dt_max) if a waypoint is surpassed.
    force_to_match_waypoint = PETSC_FALSE
    if (associated(cur_waypoint%next)) then
      if (dt_max > cur_waypoint%next%dt_max .and. &
          dt > cur_waypoint%next%dt_max .and. &
          target_time > cur_waypoint%time) then
        if (this%target_time + cur_waypoint%next%dt_max < &
            cur_waypoint%time) then
          force_to_match_waypoint = PETSC_TRUE 
        else
          dt = cur_waypoint%next%dt_max
          target_time = this%target_time + dt
        endif
      endif
    endif
!---
    ! If a waypoint calls for a plot or change in src/sinks, adjust time step
    ! to match waypoint.
    force_to_match_waypoint = WaypointForceMatchToTime(cur_waypoint) .or. &
                              force_to_match_waypoint
    equal_to_or_exceeds_waypoint = target_time + tolerance*dt >= &
                                   cur_waypoint%time
    equal_to_or_exceeds_sync_time = target_time + tolerance*dt >= sync_time
    if (equal_to_or_exceeds_sync_time .and. sync_time < cur_waypoint%time) then
      ! flip back if the sync time arrives before the waypoint time.
      equal_to_or_exceeds_waypoint = PETSC_FALSE
    endif
    if (equal_to_or_exceeds_sync_time .or. &
        (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint)) then
      if (force_to_match_waypoint) then
        max_time = min(sync_time,cur_waypoint%time)
      else
        max_time = sync_time
      endif
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
      if (dt > dt_max .and. &
                                   ! 1 sec tolerance to avoid cancellation
          dabs(dt-dt_max) > 1.d0) then 
        dt = dt_max         ! error from waypoint%time - time
        target_time = target_time + dt
      else
        target_time = max_time
        if (equal_to_or_exceeds_waypoint) then
          ! Since the time step was cut to match the waypoint, we want to set 
          ! the time step back to its prior value after the waypoint is met.
          ! %revert_dt is a flag that does so above.
          if (force_to_match_waypoint) revert_due_to_waypoint = PETSC_TRUE
          if (cur_waypoint%print_snap_output) snapshot_plot_flag = PETSC_TRUE
          if (cur_waypoint%print_obs_output) observation_plot_flag = PETSC_TRUE
          if (cur_waypoint%print_msbl_output) massbal_plot_flag = PETSC_TRUE
          if (cur_waypoint%print_checkpoint) checkpoint_flag = PETSC_TRUE
        endif
        if (equal_to_or_exceeds_sync_time) then
          ! If the time step was cut to match the sync time, we want to set
          ! the time step back to its prior value.  However, if the time step
          ! is close to its full previous value, this constraint is unnecessary
          ! and limits the ability of process model couplers "below" to catch up
          ! with those above.  Thus the conditional (dt <= .5 prev_dt) below.
          !-Also note that if this timestepper is at a depth in the process 
          ! model coupler greater than 1 (not the top process model coupler)
          ! the timestepper will constantly be reverting to sync due to the
          ! tolerance applied above without the underlying conditional.
!          if (dt < 0.99d0 * this%prev_dt) then
          if (dt <= 0.5d0 * this%prev_dt) then
            revert_due_to_sync_time = PETSC_TRUE
          endif
        endif        
        if (max_time >= cur_waypoint%time) then
          cur_waypoint => cur_waypoint%next
        endif
      endif
      exit
    else if (target_time > cur_waypoint%time) then
      cur_waypoint => cur_waypoint%next
    else
      exit
    endif
  enddo
  ! subtract 1 from max_time_steps since we still have to complete the current
  ! time step

  if (revert_due_to_sync_time .or. revert_due_to_waypoint) then
    this%revert_dt = PETSC_TRUE
    if (revert_due_to_sync_time) then
      this%num_contig_revert_due_to_sync = &
        this%num_contig_revert_due_to_sync + 1
    endif
  else
    this%num_contig_revert_due_to_sync = 0
  endif

! if coupled with CLM, max_time_step IS unlimited
! because 'waypoint' is controled by the interface
! otherwise, PF will stop at some point but CLM not-yet done
#ifndef CLM_PFLOTRAN
  if (cumulative_time_steps >= max_time_step-1) then
    nullify(cur_waypoint)
    stop_flag = TS_STOP_MAX_TIME_STEP
  endif
#endif

  ! update maximum time step size to current waypoint value
  if (associated(cur_waypoint)) then
    dt_max = cur_waypoint%dt_max
  else if (stop_flag /= TS_STOP_MAX_TIME_STEP) then
    stop_flag = TS_STOP_END_SIMULATION ! stop after end of time step
  endif
  
  option%refactor_dt = dt
  this%dt = dt
  this%dt_max = dt_max
  this%target_time = target_time
  this%cur_waypoint => cur_waypoint

 end subroutine TimestepperBaseSetTargetTime

! ************************************************************************** !

subroutine TimestepperBaseStepDT(this,process_model,stop_flag)
  ! 
  ! Steps forward one step in time
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  ! 

  use PM_Base_class
  use Option_module
  use Output_module, only : Output
  
  implicit none

  class(timestepper_base_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag
  
  type(option_type), pointer :: option

  option => process_model%option
  
  option%io_buffer = 'TimestepperBaseStepDT must be extended.'
  call printErrMsg(option)
  
end subroutine TimestepperBaseStepDT

! ************************************************************************** !

subroutine TimestepperBasePrintInfo(this,option)
  ! 
  ! Prints information about time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  
  implicit none
  
  class(timestepper_base_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(this%name) // ' Time Stepper'
    write(string,*) this%max_time_step
    write(*,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(*,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(*,'("max cuts:",x,a)') trim(adjustl(string))
    write(string,*) this%time_step_reduction_factor
    write(*,'("ts reduction factor:",x,a)') trim(adjustl(string))
    write(string,*) this%time_step_max_growth_factor
    write(*,'("ts maximum growth factor:",x,a)') trim(adjustl(string))
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,*) 
    write(option%fid_out,'(a)') trim(this%name) // ' Time Stepper'
    write(string,*) this%max_time_step
    write(option%fid_out,'("max steps:",x,a)') trim(adjustl(string))
    write(string,*) this%constant_time_step_threshold
    write(option%fid_out,'("max constant cumulative time steps:",x,a)') &
      trim(adjustl(string))
    write(string,*) this%max_time_step_cuts
    write(option%fid_out,'("max cuts:",x,a)') trim(adjustl(string))
    write(string,*) this%time_step_reduction_factor
    write(option%fid_out,'("ts reduction factor:",x,a)') trim(adjustl(string))
    write(string,*) this%time_step_max_growth_factor
    write(option%fid_out,'("ts maximum growth factor:",x,a)') &
      trim(adjustl(string))
  endif    

end subroutine TimestepperBasePrintInfo

! ************************************************************************** !

subroutine TimestepperBaseInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  
  implicit none
  
  class(timestepper_base_type) :: this

#ifdef DEBUG
  write(*,*) 'TimestepperBaseInputRecord()'
#endif

  write(*,*) 'TimestepperBaseInputRecord must be extended for &
             &each timestepper mode.'
  stop

end subroutine TimestepperBaseInputRecord

! ************************************************************************** !

subroutine TimestepperBaseCheckpointBinary(this,viewer,option)
  ! 
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none



  class(timestepper_base_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseCheckpointBinary must be extended.'
  call printErrMsg(option)  
    
end subroutine TimestepperBaseCheckpointBinary

! ************************************************************************** !

subroutine TimestepperBaseCheckpointHDF5(this, chk_grp_id, option)
  ! 
  ! Checkpoints parameters/variables associated with a time stepper to a HDF5.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15
  ! 

  use Option_module

  implicit none
  
  class(timestepper_base_type) :: this
  PetscInt :: chk_grp_id
  type(option_type) :: option

  option%io_buffer = 'TimestepperBaseCheckpointHDF5 must be extended.'
  call printErrMsg(option)

end subroutine TimestepperBaseCheckpointHDF5

! ************************************************************************** !

subroutine TimestepperBaseRestartHDF5(this, chk_grp_id, option)
  ! 
  ! Restart parameters/variables associated with a time stepper to a HDF5.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/16/15
  ! 

  use Option_module

  implicit none

  class(timestepper_base_type) :: this
  PetscInt :: chk_grp_id
  type(option_type) :: option

  option%io_buffer = 'TimestepperBaseRestartHDF5 must be extended.'
  call printErrMsg(option)

end subroutine TimestepperBaseRestartHDF5

! ************************************************************************** !

subroutine TimestepperBaseRegisterHeader(this,bag,header)
  ! 
  ! Register header entries.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  ! 
#include "petsc/finclude/petscsys.h"  
  use petscsys
  use Option_module

  implicit none
  

  class(timestepper_base_type) :: this
  class(stepper_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  ! bagsize = 8 * 8 bytes = 64 bytes
  call PetscBagRegisterReal(bag,header%time,0,"time","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%dt,0,"dt","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%prev_dt,0,"prev_dt","", &
                            ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%num_steps,0,"num_steps","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%cumulative_time_step_cuts,0, &
                           "cumulative_time_step_cuts","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%num_constant_time_steps,0, &
                           "num_constant_time_steps","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%num_contig_revert_due_to_sync,0, &
                           "num_contig_revert_due_to_sync","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%revert_dt,0, &
                           "revert_dt","",ierr);CHKERRQ(ierr)
    
end subroutine TimestepperBaseRegisterHeader

! ************************************************************************** !

subroutine TimestepperBaseSetHeader(this,bag,header)
  ! 
  ! Sets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none
  

  class(timestepper_base_type) :: this
  class(stepper_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr

  header%time = this%target_time
  header%dt = this%dt
  header%prev_dt = this%prev_dt
  header%num_steps = this%steps
  header%cumulative_time_step_cuts = this%cumulative_time_step_cuts
  header%num_constant_time_steps = this%num_constant_time_steps
  header%num_contig_revert_due_to_sync = this%num_contig_revert_due_to_sync
  header%revert_dt = ZERO_INTEGER
  if (this%revert_dt) then
    header%revert_dt = ONE_INTEGER
  endif
    
end subroutine TimestepperBaseSetHeader

! ************************************************************************** !

subroutine TimestepperBaseRestartBinary(this,viewer,option)
  ! 
  ! Restarts parameters/variables associated with
  ! a time stepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(timestepper_base_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option
  
  option%io_buffer = 'TimestepperBaseRestartBinary must be extended.'
  call printErrMsg(option)  
    
end subroutine TimestepperBaseRestartBinary

! ************************************************************************** !

subroutine TimestepperBaseGetHeader(this,header)
  ! 
  ! Gets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none
  

  class(timestepper_base_type) :: this
  class(stepper_base_header_type) :: header
  
  this%target_time = header%time
  this%dt = header%dt
  this%prev_dt = header%prev_dt
  this%steps = header%num_steps
  this%cumulative_time_step_cuts = header%cumulative_time_step_cuts
  this%num_constant_time_steps = header%num_constant_time_steps
  this%num_contig_revert_due_to_sync = header%num_contig_revert_due_to_sync
  this%revert_dt = (header%revert_dt == ONE_INTEGER)
    
end subroutine TimestepperBaseGetHeader

! ************************************************************************** !

subroutine TimestepperBaseReset(this)
  ! 
  ! Zeros timestepper object members.
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/20/14
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  

  class(timestepper_base_type) :: this
  
  this%target_time = 0.d0
  this%dt = this%dt_init
  this%prev_dt = 0.d0
  this%steps = 0
  this%cumulative_time_step_cuts = 0
  this%num_constant_time_steps = 0
  this%num_contig_revert_due_to_sync = 0
  this%revert_dt = PETSC_FALSE
    
end subroutine TimestepperBaseReset

! ************************************************************************** !

function TimestepperBaseWallClockStop(this,option)
  ! 
  ! Stops time stepping when a prescribed wall clock time has been exceeded.
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/08/14
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(timestepper_base_type) :: this
  type(option_type) :: option
  
  PetscBool :: TimestepperBaseWallclockStop
  PetscLogDouble :: current_time, average_step_time
  PetscErrorCode :: ierr
  
  ! if a simulation wallclock duration time is set, check to see that the
  ! next time step will not exceed that value.  If it does, print the
  ! checkpoint and exit
  TimestepperBaseWallclockStop = PETSC_FALSE
  if (option%wallclock_stop_flag) then
    call PetscTime(current_time, ierr)
    average_step_time = (current_time-option%start_time)/ &
                        dble(this%steps-this%start_time_step+1) &
                        *2.d0  ! just to be safe, double it
    if (average_step_time + current_time > option%wallclock_stop_time) then
      TimestepperBaseWallclockStop = PETSC_TRUE
    endif
  endif
  
end function TimestepperBaseWallClockStop


! ************************************************************************** !

subroutine TimestepperBasePrintEKG(this)
  ! 
  ! Deallocates members of a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_base_type) :: this
  
  
  
end subroutine TimestepperBasePrintEKG

! ************************************************************************** !

recursive subroutine TimestepperBaseFinalizeRun(this,option)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none
  
  class(timestepper_base_type) :: this
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
#ifdef DEBUG
  call printMsg(option,'TimestepperBaseFinalizeRun()')
#endif
  
  if (OptionPrintToScreen(option)) then
    write(*,'(/,a," TS Base"," steps = ",i6," cuts = ",i6)') &
            trim(this%name), &
            this%steps, &
            this%cumulative_time_step_cuts
    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,'(a)') trim(this%name) // ' TS Base solver time = ' // &
      trim(adjustl(string)) // ' seconds'
  endif
  
end subroutine TimestepperBaseFinalizeRun

! ************************************************************************** !

subroutine TimestepperBaseStrip(this)
  ! 
  ! Deallocates members of a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  ! 

  implicit none
  
  class(timestepper_base_type) :: this
  
end subroutine TimestepperBaseStrip

! ************************************************************************** !

subroutine TimestepperBaseDestroy(this)
  ! 
  ! Deallocates a time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  class(timestepper_base_type) :: this
  
  call TimestepperBaseStrip(this)
    
!  deallocate(this)
!  nullify(this)
  
end subroutine TimestepperBaseDestroy

end module Timestepper_Base_class
