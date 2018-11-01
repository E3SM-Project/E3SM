module Timestepper_Surface_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Timestepper_Base_class
  use Solver_module
  use Waypoint_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  type, public, extends(timestepper_base_type) :: timestepper_surface_type
    PetscReal :: dt_max_allowable
    PetscReal :: surf_subsurf_coupling_flow_dt
    type(solver_type), pointer :: solver
  contains
    procedure, public :: CheckpointBinary => TimestepperSurfaceCheckpointBinary
    procedure, public :: Init => TimestepperSurfaceInit
    procedure, public :: RestartBinary => TimestepperSurfaceRestartBinary
    procedure, public :: Reset => TimestepperSurfaceReset
    procedure, public :: SetTargetTime => TimestepperSurfaceSetTargetTime
    procedure, public :: InputRecord => TimestepperSurfInputRecord
    procedure, public :: Strip => TimestepperSurfaceStrip
    procedure, public :: StepDT => TimestepperSurfaceStepDT
  end type timestepper_surface_type

  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: timestepper_surface_header_type
    real*8 :: dt_max_allowable
    real*8 :: surf_subsurf_coupling_flow_dt
  end type timestepper_surface_header_type
  PetscSizeT, parameter, private :: bagsize = 80 ! 64 (base) + 16 (BE)

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
#include "petsc/finclude/petscsys.h"
      use petscsys
      import :: timestepper_surface_header_type
      implicit none
      PetscBag :: bag
      class(timestepper_surface_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public TimestepperSurfaceSetTargetTime, &
         TimestepperSurfaceCreate

contains

! ************************************************************************** !

function TimestepperSurfaceCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/03/13
  ! 

  implicit none
  
  class(timestepper_surface_type), pointer :: TimestepperSurfaceCreate
  
  class(timestepper_surface_type), pointer :: surf_timestepper
  
  allocate(surf_timestepper)
  call surf_timestepper%Init()
  
  surf_timestepper%solver => SolverCreate()
  
  TimestepperSurfaceCreate => surf_timestepper
  
end function TimestepperSurfaceCreate

! ************************************************************************** !

subroutine TimestepperSurfaceInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/03/13
  ! 

  implicit none
  
  class (timestepper_surface_type) :: this

  call TimestepperBaseInit(this)

  this%dt_max_allowable = 0.d0
  this%surf_subsurf_coupling_flow_dt = 0.d0
  
end subroutine TimestepperSurfaceInit

! ************************************************************************** !

subroutine TimestepperSurfaceSetTargetTime(this,sync_time,option,stop_flag, &
                                           snapshot_plot_flag, &
                                           observation_plot_flag, &
                                           massbal_plot_flag,checkpoint_flag)
  ! 
  ! This routine ?
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/02/13
  ! 

  use Option_module

  implicit none

  class(timestepper_surface_type) :: this
  PetscReal :: sync_time
  type(option_type) :: option
  PetscInt :: stop_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag

  PetscReal :: dt
  PetscReal :: target_time
  PetscReal :: max_time
  PetscReal :: tolerance
  PetscBool :: equal_to_or_exceeds_waypoint
  PetscBool :: equal_to_or_exceeds_sync_time
  PetscBool :: force_to_match_waypoint
  type(waypoint_type), pointer :: cur_waypoint
  
  cur_waypoint => this%cur_waypoint

  dt = min(this%dt_max_allowable,this%dt_max)
  target_time = this%target_time + dt
  tolerance = this%time_step_tolerance

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  
  if (cur_waypoint%time < 1.d-40) then
    cur_waypoint => cur_waypoint%next
  endif

  force_to_match_waypoint = WaypointForceMatchToTime(cur_waypoint)
  equal_to_or_exceeds_waypoint = target_time + tolerance*dt >= cur_waypoint%time
  equal_to_or_exceeds_sync_time = target_time + tolerance*dt >= sync_time

  if (equal_to_or_exceeds_sync_time .or. &
      (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint)) then

    max_time = min(sync_time,cur_waypoint%time)
    ! decrement by time step size
    target_time = target_time - dt
    ! set new time step size based on max time
    dt = max_time - target_time
    target_time = target_time + dt

    if (Equal(max_time,cur_waypoint%time)) then
      if (cur_waypoint%print_snap_output) snapshot_plot_flag = PETSC_TRUE
      if (cur_waypoint%print_checkpoint) checkpoint_flag = PETSC_TRUE
    endif

    if (equal_to_or_exceeds_sync_time) then
      max_time = sync_time
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time
      target_time = target_time + dt

    endif

    if (equal_to_or_exceeds_waypoint .and. force_to_match_waypoint) then
      max_time = cur_waypoint%time
      ! decrement by time step size
      target_time = target_time - dt
      ! set new time step size based on max time
      dt = max_time - target_time

      target_time = target_time + dt

      if (cur_waypoint%print_snap_output) snapshot_plot_flag = PETSC_TRUE
      if (cur_waypoint%print_checkpoint) checkpoint_flag = PETSC_TRUE
    endif
  
  endif

  if (target_time >= cur_waypoint%time) then
    cur_waypoint => cur_waypoint%next
  endif
  this%dt = dt
  this%target_time = target_time
  this%cur_waypoint => cur_waypoint
  if (.not.associated(cur_waypoint)) stop_flag = TS_STOP_END_SIMULATION

end subroutine TimestepperSurfaceSetTargetTime

! ************************************************************************** !

subroutine TimestepperSurfaceStepDT(this,process_model,stop_flag)
  ! 
  ! This is a dummy routine added to be extended in timestepper_surface_type
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/03/13
  ! 

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Surface_Flow_class
  use Option_module
  use Output_module, only : Output
  use Surface_Flow_module
  
  implicit none

  class(timestepper_surface_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  PetscReal :: time
  PetscReal :: dtime
  PetscReal :: tmp
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  solver => this%solver
  option => process_model%option

  call process_model%PreSolve()

  call TSSetTimeStep(solver%ts,option%surf_flow_dt,ierr);CHKERRQ(ierr)
#if (PETSC_VERSION_MINOR >= 8)
  call TSSetStepNumber(solver%ts,ZERO_INTEGER,ierr);CHKERRQ(ierr)
#endif
  call TSSetExactFinalTime(solver%ts,TS_EXACTFINALTIME_MATCHSTEP, &
                           ierr);CHKERRQ(ierr)
  call TSSolve(solver%ts,process_model%solution_vec,ierr);CHKERRQ(ierr)
  call TSGetTime(solver%ts,time,ierr);CHKERRQ(ierr)
  call TSGetTimeStep(solver%ts,dtime,ierr);CHKERRQ(ierr)

  call process_model%PostSolve()

  this%steps = this%steps + 1

  if (option%print_screen_flag) then
    write(*, '(" SURFACE FLOW ",i6," Time= ",1pe12.5," Dt= ",1pe12.5," [",a1,"]")') &
      this%steps, &
      time/process_model%output_option%tconv, &
      dtime/process_model%output_option%tconv, &
      process_model%output_option%tunit
  endif

end subroutine TimestepperSurfaceStepDT

! ************************************************************************** !

subroutine TimestepperSurfaceCheckpointBinary(this,viewer,option)
  ! 
  ! This checkpoints parameters/variables associated with surface-timestepper
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/18/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_surface_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(timestepper_surface_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperSurfaceRegisterHeader(this,bag,header)
  call TimestepperSurfaceSetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperSurfaceCheckpointBinary

! ************************************************************************** !

subroutine TimestepperSurfaceRestartBinary(this,viewer,option)
  ! 
  ! This checkpoints parameters/variables associated with surface-timestepper
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/18/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_surface_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(timestepper_surface_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperSurfaceRegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call TimestepperSurfaceGetHeader(this,header)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperSurfaceRestartBinary

! ************************************************************************** !

subroutine TimestepperSurfaceRegisterHeader(this,bag,header)
  ! 
  ! This subroutine register header entries for surface-flow.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(timestepper_surface_type) :: this
  class(timestepper_surface_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  ! bagsize = 2 * 8 bytes = 16 bytes
  call PetscBagRegisterReal(bag,header%dt_max_allowable,0.d0, &
                           "dt_max_allowable","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterReal(bag,header%surf_subsurf_coupling_flow_dt,0.d0, &
                           "surf_subsurf_coupling_flow_dt","", &
                            ierr);CHKERRQ(ierr)

  call TimestepperBaseRegisterHeader(this,bag,header)

end subroutine TimestepperSurfaceRegisterHeader

! ************************************************************************** !

subroutine TimestepperSurfaceSetHeader(this,bag,header)
  ! 
  ! This subroutine sets values in checkpoint header.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_surface_type) :: this
  class(timestepper_surface_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  header%dt_max_allowable = this%dt_max_allowable
  header%surf_subsurf_coupling_flow_dt = this%surf_subsurf_coupling_flow_dt

  call TimestepperBaseSetHeader(this,bag,header)

end subroutine TimestepperSurfaceSetHeader

! ************************************************************************** !

subroutine TimestepperSurfaceGetHeader(this,header)
  ! 
  ! This subroutine gets values in checkpoint header.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 09/19/13
  ! 

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"

  class(timestepper_surface_type) :: this
  class(timestepper_surface_header_type) :: header

  PetscErrorCode :: ierr

  this%dt_max_allowable = header%dt_max_allowable
  this%surf_subsurf_coupling_flow_dt = header%surf_subsurf_coupling_flow_dt

  call TimestepperBaseGetHeader(this,header)

  call TSSetTime(this%solver%ts,this%target_time,ierr);CHKERRQ(ierr)

end subroutine TimestepperSurfaceGetHeader

! ************************************************************************** !

subroutine TimestepperSurfaceReset(this)

  implicit none

  class(timestepper_surface_type) :: this

  PetscErrorCode :: ierr

#if 0
  !TODO(Gautam): set these back to their initial values as if a simulation
  !              were initialized, but not yet run
  this%dt_max_allowable = header%dt_max_allowable
  this%surf_subsurf_coupling_flow_dt = header%surf_subsurf_coupling_flow_dt

  call TimestepperBaseReset(this)

  !TODO(Gautam): this%target_time is set to 0.d0 in TimestepperBaseReset(). Is
  !              that OK? - Glenn
  call TSSetTime(this%solver%ts,this%target_time,ierr);CHKERRQ(ierr)
#endif

end subroutine TimestepperSurfaceReset

! ************************************************************************** !

subroutine TimestepperSurfacePrintInfo(this,option)
  ! 
  ! Prints settings for base timestepper.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  ! 
  use Option_module

  implicit none
  
#include "petsc/finclude/petscts.h"  

  class(timestepper_surface_type) :: this
  type(option_type) :: option
  
  PetscErrorCode :: ierr
  
  if (OptionPrintToScreen(option)) then
    write(*,*) ' '
    write(*,*) 'Surface Flow TS Solver:'
    call TSView(this%solver%ts,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  call TimestepperBasePrintInfo(this,option)
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)
  
end subroutine TimestepperSurfacePrintInfo

! ************************************************************************** !

subroutine TimestepperSurfInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  
  implicit none
  
  class(timestepper_surface_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word
   
  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'max timestep size: '
  write(word,*) this%dt_max_allowable
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperSurfInputRecord

! ************************************************************************** !

subroutine TimestepperSurfaceStrip(this)
  ! 
  ! Deallocates members of a surface time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 

  implicit none
  
  class(timestepper_surface_type) :: this
  
  call TimestepperBaseStrip(this)
  call SolverDestroy(this%solver)

end subroutine TimestepperSurfaceStrip

! ************************************************************************** !

subroutine TimestepperSurfaceDestroy(this)
  ! 
  ! Deallocates a surface time stepper
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/14
  ! 

  implicit none
  
  class(timestepper_surface_type) :: this
  
  call TimestepperSurfaceStrip(this)
  
end subroutine TimestepperSurfaceDestroy

end module Timestepper_Surface_class
