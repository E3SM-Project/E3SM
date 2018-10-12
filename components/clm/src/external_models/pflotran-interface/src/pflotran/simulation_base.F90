module Simulation_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PMC_Base_class
  use PM_Base_class
  use Option_module
  use Output_Aux_module
  use Output_module
  use Simulation_Aux_module
  use Waypoint_module
  
  use PFLOTRAN_Constants_module

  implicit none

  
  private

  type, public :: simulation_base_type
    type(option_type), pointer :: option
    type(waypoint_list_type), pointer :: waypoint_list_outer ! for outer sync loop
    type(checkpoint_option_type), pointer :: checkpoint_option
    type(output_option_type), pointer :: output_option
    PetscInt :: stop_flag
    class(pmc_base_type), pointer :: process_model_coupler_list
    class(pm_base_type), pointer :: process_model_list
    type(simulation_aux_type), pointer :: sim_aux
  contains
    procedure, public :: Init => SimulationBaseInit
    procedure, public :: InitializeRun => SimulationBaseInitializeRun
    procedure, public :: InputRecord => SimulationInputRecord
    procedure, public :: JumpStart => SimulationBaseJumpStart
    procedure, public :: ExecuteRun
    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SimulationBaseFinalizeRun
    procedure, public :: Strip => SimulationBaseStrip
  end type simulation_base_type
  
  public :: SimulationBaseCreate, &
            SimulationBaseInit, &
            SimulationBaseInitializeRun, &
            SimulationInputRecordPrint, &
            SimulationInputRecord, &
            SimulationGetFinalWaypointTime, &
            SimulationBaseFinalizeRun, &
            SimulationBaseStrip, &
            SimulationBaseDestroy
  
contains

! ************************************************************************** !

function SimulationBaseCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Option_module

  implicit none
  
  class(simulation_base_type), pointer :: SimulationBaseCreate

  type(option_type), pointer :: option
  
  allocate(SimulationBaseCreate)
  call SimulationBaseCreate%Init(option)

end function SimulationBaseCreate

! ************************************************************************** !

subroutine SimulationBaseInit(this,option)
  ! 
  ! Initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 
  use Timestepper_Base_class, only : TS_CONTINUE
  use Option_module
  use Output_Aux_module
  use Waypoint_module

  implicit none
  
  class(simulation_base_type) :: this
  type(option_type), pointer :: option

  this%option => option
  this%waypoint_list_outer => WaypointListCreate()
  this%output_option => OutputOptionCreate()
  nullify(this%checkpoint_option)
  nullify(this%process_model_coupler_list)
  nullify(this%process_model_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = TS_CONTINUE

end subroutine SimulationBaseInit

! ************************************************************************** !

subroutine SimulationBaseInitializeRun(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module
  use Option_module
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  implicit none
  

  class(simulation_base_type) :: this

  PetscInt :: chk_grp_id
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseInitializeRun()')
#endif
  
  ! the user may request output of variable that do not exist for the 
  ! the requested process models; this routine should catch such issues.
  call OutputEnsureVariablesExist(this%output_option,this%option)
  if (associated(this%process_model_coupler_list)) then
    if (this%option%restart_flag) then
      if (index(this%option%restart_filename,'.chk') > 0) then
        call this%process_model_coupler_list%RestartBinary(viewer)
      elseif (index(this%option%restart_filename,'.h5') > 0) then
#if !defined(PETSC_HAVE_HDF5)
         this%option%io_buffer = 'HDF5 formatted restart not supported &
              &unless PFLOTRAN is compiled with HDF5 libraries enabled.'
         call printErrMsg(this%option)
#else
        call this%process_model_coupler_list%RestartHDF5(chk_grp_id)
#endif
      else
        this%option%io_buffer = 'Unknown restart filename format. ' // &
        'Only *.chk and *.h5 supported.'
        call printErrMsg(this%option)
      endif
    endif
  
    ! initialize performs overwrite of restart, if applicable
    call this%process_model_coupler_list%InitializeRun()  
    call this%JumpStart()
  endif
  
  call SimulationInputRecordPrint(this)
  call printMsg(this%option," ")
  call printMsg(this%option,"  Finished Initialization")
  call PetscLogEventEnd(logging%event_init,ierr);CHKERRQ(ierr)
  ! pushed in PFLOTRANInitializePostPetsc()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in FinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)
  
end subroutine SimulationBaseInitializeRun

! ************************************************************************** !

subroutine SimulationInputRecordPrint(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  use Checkpoint_module

  implicit none
  
  class(simulation_base_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  PetscBool :: is_open

  inquire(id, OPENED=is_open)
  if (is_open .and. OptionPrintToFile(this%option)) then
  !----------------------------------------------------------------------------
    ! print checkpoint information
    call CheckpointInputRecord(this%checkpoint_option,this%waypoint_list_outer)
  
    write(id,'(a)') ' '
    ! print process model coupler and process model information
    call this%process_model_coupler_list%InputRecord()
    
    ! print simulation-specific information
    call this%InputRecord()
  !----------------------------------------------------------------------------
  endif

end subroutine SimulationInputRecordPrint

! ************************************************************************** !

subroutine SimulationInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! This subroutine must be extended in the extended simulation objects.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 

  implicit none
  
  class(simulation_base_type) :: this

#ifdef DEBUG
  call printMsg(this%option,'SimulationInputRecord()')
#endif

  this%option%io_buffer = 'SimulationInputRecord must be extended for ' // &
    'each simulation mode.'
  call printErrMsg(this%option)

end subroutine SimulationInputRecord

! ************************************************************************** !

subroutine SimulationBaseJumpStart(this)
  ! 
  ! Gets the time stepping, etc. up and running
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  ! 
  use Option_module
  
  implicit none
  
  class(simulation_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseJumpStart()')
#endif

  this%option%io_buffer = 'SimulationBaseJumpStart must be extended for ' // &
    'each simulation mode.'
  call printErrMsg(this%option)
  
end subroutine SimulationBaseJumpStart

! ************************************************************************** !

subroutine ExecuteRun(this)
  ! 
  ! Initializes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Waypoint_module
  use Timestepper_Base_class, only : TS_CONTINUE
  use Checkpoint_module

  implicit none
  
  class(simulation_base_type) :: this
  
  PetscReal :: final_time
  PetscReal :: sync_time
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXSTRINGLENGTH) :: append_name

#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseExecuteRun()')
#endif

  if (.not.associated(this%process_model_coupler_list)) then
    return
  endif

  append_name = '-restart'

  final_time = SimulationGetFinalWaypointTime(this)
  cur_waypoint => this%waypoint_list_outer%first
  call WaypointSkipToTime(cur_waypoint,this%option%time)
  do
    if (this%stop_flag /= TS_CONTINUE) exit ! end simulation
    if (.not.associated(cur_waypoint)) exit
    call this%RunToTime(min(final_time,cur_waypoint%time))
    cur_waypoint => cur_waypoint%next
  enddo
  if (associated(this%process_model_coupler_list%checkpoint_option)) then
    call this%process_model_coupler_list%Checkpoint(append_name)
  endif

end subroutine ExecuteRun

! ************************************************************************** !

subroutine RunToTime(this,target_time)
  ! 
  ! Executes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Simulation_Aux_module

  implicit none

  class(simulation_base_type) :: this
  PetscReal :: target_time
  
  class(pmc_base_type), pointer :: cur_process_model_coupler
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseRunToTime()')
#endif
  
  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine RunToTime

! ************************************************************************** !

subroutine SimulationBaseFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  use Logging_module
  use Timestepper_Base_class, only : TS_STOP_WALLCLOCK_EXCEEDED
  
  implicit none
  
  class(simulation_base_type) :: this
  
  PetscErrorCode :: ierr
  
  class(pmc_base_type), pointer :: cur_process_model_coupler

#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseFinalizeRun()')
#endif
  
  if (this%stop_flag == TS_STOP_WALLCLOCK_EXCEEDED) then
    call printMsg(this%option,"Wallclock stop time exceeded.  Exiting!!!")
    call printMsg(this%option,"")
  endif
  
  if (associated(this%process_model_coupler_list)) then
    call this%process_model_coupler_list%FinalizeRun()
  endif
  
  ! pushed in InitializeRun()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)
  ! popped in OptionFinalize()
  call PetscLogStagePush(logging%stage(FINAL_STAGE),ierr);CHKERRQ(ierr)
  
end subroutine SimulationBaseFinalizeRun

! ************************************************************************** !

function SimulationGetFinalWaypointTime(this)
  ! 
  ! Returns the earliest final waypoint time
  ! from the top layer of process model
  ! couplers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/12/13
  ! 

  use Waypoint_module

  implicit none
  
  class(simulation_base_type) :: this
  
  PetscReal :: SimulationGetFinalWaypointTime

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscReal :: final_time
  
  SimulationGetFinalWaypointTime = 0.d0
  
  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    final_time = WaypointListGetFinalTime(cur_process_model_coupler% &
                                            waypoint_list)
    if (SimulationGetFinalWaypointTime < 1.d-40 .or. &
        final_time < SimulationGetFinalWaypointTime) then
      SimulationGetFinalWaypointTime = final_time
    endif
    cur_process_model_coupler => cur_process_model_coupler%peer
  enddo

end function SimulationGetFinalWaypointTime

! ************************************************************************** !

subroutine SimulationBaseStrip(this)
  ! 
  ! Deallocates members of simulation base
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 
  use Input_Aux_module
  use Waypoint_module
  use EOS_module
  
  implicit none
  
  class(simulation_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'SimulationBaseStrip()')
#endif
  call WaypointListDestroy(this%waypoint_list_outer)
  call SimAuxDestroy(this%sim_aux)
  call CheckpointOptionDestroy(this%checkpoint_option)
  call OutputOptionDestroy(this%output_option)
  if (associated(this%process_model_coupler_list)) then
    call this%process_model_coupler_list%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%process_model_coupler_list)
    nullify(this%process_model_coupler_list)
  endif
  call InputDbaseDestroy()

  call AllEOSDBaseDestroy()
  
end subroutine SimulationBaseStrip

! ************************************************************************** !

subroutine SimulationBaseDestroy(simulation)
  ! 
  ! Deallocates a simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  implicit none
  
  class(simulation_base_type), pointer :: simulation
  
#ifdef DEBUG
  call printMsg(simulation%option,'SimulationDestroy()')
#endif
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SimulationBaseDestroy
  
end module Simulation_Base_class
