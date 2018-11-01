! Process Model Coupler Base class
module PMC_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys  
  use PM_Base_class
  use Timestepper_Base_class
  use Option_module
  use Output_Aux_module
  use Waypoint_module
  use PM_Base_Pointer_module
  use Output_module, only : Output
  use Simulation_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
  PetscInt, parameter, public :: PM_CHILD = 0
  PetscInt, parameter, public :: PM_PEER = 1
  PetscInt, parameter, public :: PM_APPEND = 0
  PetscInt, parameter, public :: PM_INSERT = 1
  
  ! process model coupler type
  type, public :: pmc_base_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: stage
    PetscBool :: is_master
    type(option_type), pointer :: option
    type(checkpoint_option_type), pointer :: checkpoint_option
    class(timestepper_base_type), pointer :: timestepper
    class(pm_base_type), pointer :: pm_list
    type(waypoint_list_type), pointer :: waypoint_list
    class(pmc_base_type), pointer :: child
    class(pmc_base_type), pointer :: peer
    type(pm_base_pointer_type), pointer :: pm_ptr
    type(simulation_aux_type),pointer :: sim_aux
    procedure(Output), nopass, pointer :: Output
  contains
    procedure, public :: Init => PMCBaseInit
    procedure, public :: InitializeRun
    procedure, public :: InputRecord => PMCBaseInputRecord
    procedure, public :: CastToBase => PMCCastToBase
    procedure, public :: SetTimestepper => PMCBaseSetTimestepper
    procedure, public :: SetupSolvers => PMCBaseSetupSolvers
    procedure, public :: RunToTime => PMCBaseRunToTime
    procedure, public :: Checkpoint => PMCBaseCheckpoint
    procedure, public :: CheckpointBinary => PMCBaseCheckpointBinary
    procedure, public :: RestartBinary => PMCBaseRestartBinary
#if defined(PETSC_HAVE_HDF5)
    procedure, public :: CheckpointHDF5 => PMCBaseCheckpointHDF5
    procedure, public :: RestartHDF5 => PMCBaseRestartHDF5
#endif
    procedure, public :: FinalizeRun
    procedure, public :: OutputLocal
    procedure, public :: UpdateSolution => PMCBaseUpdateSolution
    procedure, public :: Destroy => PMCBaseDestroy
    procedure, public :: AccumulateAuxData
    procedure, public :: GetAuxData
    procedure, public :: SetAuxData
    procedure, public :: CheckNullPM => PMCBaseCheckNullPM
    !procedure, public :: SetChildPeerPtr => PMCBaseSetChildPeerPtr
  end type pmc_base_type
  
  abstract interface
    subroutine Synchronize(pmc)
      import pmc_base_type
      implicit none
        class(pmc_base_type) :: pmc
    end subroutine Synchronize
  end interface

  ! For checkpointing  
  type, public :: pmc_base_header_type
    PetscInt :: plot_number      ! in the checkpoint file format
    PetscInt :: times_per_h5_file! in the checkpoint file format
  end type pmc_base_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
#include "petsc/finclude/petscsys.h"
      use petscsys
      import :: pmc_base_header_type
      implicit none
      PetscBag :: bag
      class(pmc_base_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData   
    
  public :: PMCBaseCreate, &
            PMCBaseInit, &
            PMCBaseInputRecord, &
            PMCBaseSetChildPeerPtr, &
            PMCBaseStrip, &
            SetOutputFlags, &
            PMCCastToBase
  
contains

! ************************************************************************** !

function PMCBaseCreate()
  ! 
  ! Allocates and initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type), pointer :: PMCBaseCreate
  
  class(pmc_base_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCBase%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()

  PMCBaseCreate => pmc  
  
end function PMCBaseCreate

! ************************************************************************** !

subroutine PMCBaseInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
#ifdef DEBUG
  print *, 'PMCBase%Init()'
#endif
  
  this%name = 'PMCBase'
  this%stage = 0
  this%is_master = PETSC_FALSE
  nullify(this%option)
  nullify(this%checkpoint_option)
  nullify(this%timestepper)
  nullify(this%pm_list)
  nullify(this%waypoint_list)
  nullify(this%child)
  nullify(this%peer)
  nullify(this%sim_aux)
  this%Output => Null()
  
  allocate(this%pm_ptr)
  nullify(this%pm_ptr%pm)
  
end subroutine PMCBaseInit

! ************************************************************************** !

recursive subroutine PMCBaseInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 

  use PM_Base_class
  
  implicit none
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  ! print information about self
  write(id,'(a)') ' '
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') ' '
  write(id,'(a29)',advance='no') 'pmc: '
  write(id,'(a)') this%name
  if (associated(this%timestepper)) then
    call this%timestepper%inputrecord
  endif
  cur_pm => this%pm_list
  do ! loop through this pmc's process models
    if (.not.associated(cur_pm)) exit
    call cur_pm%inputrecord
    cur_pm => cur_pm%next
  enddo

  ! print information about child's pmc
  if (associated(this%child)) then
    call this%child%inputrecord
  endif
  ! print information about peer's pmc
  if (associated(this%peer)) then
    call this%peer%inputrecord
  endif
  
end subroutine PMCBaseInputRecord

! ************************************************************************** !

recursive subroutine PMCBaseSetChildPeerPtr(pmcA,relationship_to,pmcB, &
                                            pmcB_parent,position_instruction)
  ! 
  ! Orders pmcA under pmcB relative to the specified relationship (child or 
  ! peer) and insert/append instruction. If pmcB's parent is not relevant, a
  ! null dummy pointer of pmc_base_type should be passed in.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 02/15/2017
  ! 

  use PM_Base_class
  use Option_module
  
  implicit none
  
  class(pmc_base_type), pointer :: pmcA
  PetscInt :: relationship_to
  class(pmc_base_type), pointer :: pmcB
  class(pmc_base_type), pointer :: pmcB_parent
  class(pmc_base_type), pointer :: pmcB_dummy
  PetscInt :: position_instruction
  
  PetscInt :: new_relationship

  nullify(pmcB_dummy)
  
  if (.not.associated(pmcB)) then
    pmcA%option%io_buffer = 'pmcB passed PMCBaseSetChildPeerPtr is not &
                            &associated.'
    call printErrMsg(pmcA%option)
  endif

  select case(relationship_to)
  !--------------------------------------------------------
    case(PM_CHILD)
      if (associated(pmcB%child)) then
        new_relationship = PM_PEER
        call PMCBaseSetChildPeerPtr(pmcA,new_relationship,pmcB%child,pmcB, &
                                    position_instruction)
      else
        pmcB%child => pmcA
#ifdef DEBUG
        pmcA%option%io_buffer = trim(pmcA%name)// ' assigned as first child&
                                & of ' // trim(pmcB%name) // '.'
        call printMsg(pmcA%option)
#endif
      endif
  !--------------------------------------------------------
    case(PM_PEER)
      select case(position_instruction)
      !----------------------------------------------------
        case(PM_APPEND)
          if (associated(pmcB%peer)) then
            new_relationship = PM_PEER
            call PMCBaseSetChildPeerPtr(pmcA,new_relationship,pmcB%peer, &
                                        pmcB_dummy,position_instruction)
          else
            pmcB%peer => pmcA
#ifdef DEBUG
        pmcA%option%io_buffer = trim(pmcA%name) // ' assigned as peer of ' // &
                                trim(pmcB%name) // ' via "append".'
        call printMsg(option)
#endif
          endif
      !----------------------------------------------------
        case(PM_INSERT)
          pmcA%peer => pmcB
          if (associated(pmcB_parent)) then
            pmcB_parent%child => pmcA
#ifdef DEBUG
        pmcA%option%io_buffer = trim(pmcA%name)// ' assigned as first child&
                                & of ' // trim(pmcB%name) // ' via "insert".'
        call printMsg(option)
#endif
          else
            pmcA%option%io_buffer = 'Null pointer for pmcB_parent passed into &
                                     &PMCBaseSetChildPeerPtr.'
            call printErrMsg(pmcA%option)
          endif
      !----------------------------------------------------
      end select
  !--------------------------------------------------------
    case default
      pmcA%option%io_buffer = 'PMC relationship_to not understood in &
                              &PMCBaseSetChildPeerPtr.'
      call printErrMsg(pmcA%option)
  !--------------------------------------------------------
  end select
  
end subroutine PMCBaseSetChildPeerPtr

! ************************************************************************** !

function PMCCastToBase(this)
  ! 
  ! PMCBaseCastToBase: Casts an extended PMC to a pointer to the base class
  !                    in order to avoid a 'select type' statement when
  !                    pointing a pmc_base_type pointer to an extended class.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_base_type), target :: this
  
  class(pmc_base_type), pointer :: PMCCastToBase

  PMCCastToBase => this
  
end function PMCCastToBase

! ************************************************************************** !

subroutine PMCBaseSetupSolvers(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

#ifdef DEBUG
  call printMsg(this%option,'PMCBase%SetupSolvers()')
#endif

  ! For now there is nothing to be done here.  Most everything is done in
  ! PMCSubsurface
  
end subroutine PMCBaseSetupSolvers

! ************************************************************************** !

subroutine PMCBaseSetTimestepper(this,timestepper)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
  
  implicit none
  
  class(pmc_base_type) :: this
  class(timestepper_base_type), pointer :: timestepper

#ifdef DEBUG
  call printMsg(this%option,'PMCBase%SetTimestepper()')
#endif
  
  this%timestepper => timestepper
  
end subroutine PMCBaseSetTimestepper

! ************************************************************************** !

recursive subroutine InitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%InitializeRun()')
#endif
  
  if (associated(this%timestepper)) then
    call this%timestepper%InitializeRun(this%option)
  endif
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%InitializeRun()
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%child)) then
    call this%child%InitializeRun()
  endif
  
  if (associated(this%peer)) then
    call this%peer%InitializeRun()
  endif

end subroutine InitializeRun

! ************************************************************************** !

recursive subroutine PMCBaseRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Timestepper_Base_class
  use Checkpoint_module

  implicit none
  
  class(pmc_base_type), target :: this
  character(len=MAXSTRINGLENGTH) :: filename_append
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscBool :: failure
  PetscBool :: checkpoint_at_this_time_flag
  PetscBool :: snapshot_plot_at_this_time_flag
  PetscBool :: observation_plot_at_this_time_flag
  PetscBool :: massbal_plot_at_this_time_flag
  PetscBool :: checkpoint_at_this_timestep_flag
  PetscBool :: snapshot_plot_at_this_timestep_flag
  PetscBool :: observation_plot_at_this_timestep_flag
  PetscBool :: massbal_plot_at_this_timestep_flag
  PetscBool :: peer_already_run_to_time
  class(pm_base_type), pointer :: cur_pm
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return

  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)  
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE
  do
    if (local_stop_flag /= TS_CONTINUE) exit ! end simulation
    if (this%timestepper%target_time >= sync_time) exit
    
    call SetOutputFlags(this)
    checkpoint_at_this_time_flag = PETSC_FALSE
    snapshot_plot_at_this_time_flag = PETSC_FALSE
    observation_plot_at_this_time_flag = PETSC_FALSE
    massbal_plot_at_this_time_flag = PETSC_FALSE
    checkpoint_at_this_timestep_flag = PETSC_FALSE
    snapshot_plot_at_this_timestep_flag = PETSC_FALSE
    observation_plot_at_this_timestep_flag = PETSC_FALSE
    massbal_plot_at_this_timestep_flag = PETSC_FALSE
    
    call this%timestepper%SetTargetTime(sync_time,this%option, &
                                        local_stop_flag, &
                                        snapshot_plot_at_this_time_flag, &
                                        observation_plot_at_this_time_flag, &
                                        massbal_plot_at_this_time_flag, &
                                        checkpoint_at_this_time_flag)
    call this%timestepper%StepDT(this%pm_list,local_stop_flag)
    if (this%timestepper%time_step_cut_flag) then
      ! if timestep has been cut, all the I/O flags set above in 
      ! %SetTargetTime, which are based on waypoints times, not time step,
      ! should be turned off
      snapshot_plot_at_this_time_flag = PETSC_FALSE
      observation_plot_at_this_time_flag = PETSC_FALSE
      massbal_plot_at_this_time_flag = PETSC_FALSE
      checkpoint_at_this_time_flag = PETSC_FALSE
    endif
    if (local_stop_flag == TS_STOP_FAILURE) exit ! failure
    ! Have to loop over all process models coupled in this object and update
    ! the time step size.  Still need code to force all process models to
    ! use the same time step size if tightly or iteratively coupled.
    cur_pm => this%pm_list
    do
      if (.not.associated(cur_pm)) exit
      ! have to update option%time for conditions
      this%option%time = this%timestepper%target_time
      call cur_pm%UpdateSolution()
      call this%timestepper%UpdateDT(cur_pm)
      cur_pm => cur_pm%next
    enddo

    ! Accumulate data needed by process-model
    call this%AccumulateAuxData()

    ! Run underlying process model couplers
    if (associated(this%child)) then
      ! Set data needed by process-models
      call this%SetAuxData()
      call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
      ! Get data from other process-models
      call this%GetAuxData()
    endif

    ! only print output for process models of depth 0
    if (associated(this%Output)) then
      ! however, if we are using the modulus of the output_option%imod, we may
      ! still print
      snapshot_plot_at_this_timestep_flag = &
        (mod(this%timestepper%steps,this%pm_list% &
             output_option%periodic_snap_output_ts_imod) == 0)
      observation_plot_at_this_timestep_flag = &
        (mod(this%timestepper%steps,this%pm_list% &
             output_option%periodic_obs_output_ts_imod) == 0)
      massbal_plot_at_this_timestep_flag = &
        (mod(this%timestepper%steps,this%pm_list% &
             output_option%periodic_msbl_output_ts_imod) == 0)
      
      if (this%option%steady_state) &
        snapshot_plot_at_this_timestep_flag = PETSC_TRUE
      
      call this%Output(this%pm_list%realization_base, &
                       (snapshot_plot_at_this_time_flag .or. &
                        snapshot_plot_at_this_timestep_flag), &
                       (observation_plot_at_this_time_flag .or. &
                        observation_plot_at_this_timestep_flag), &
                       (massbal_plot_at_this_time_flag .or. &
                        massbal_plot_at_this_timestep_flag))
    endif
    
    if (this%is_master .and. associated(this%checkpoint_option)) then
      if (this%checkpoint_option%periodic_ts_incr > 0 .and. &
          mod(this%timestepper%steps, &
              this%checkpoint_option%periodic_ts_incr) == 0) then
        checkpoint_at_this_timestep_flag = PETSC_TRUE
      endif
    endif

    ! Checkpointing forces peers to be executed prior to the checkpoing.  If
    ! so, we need to skip the peer RunToTime outside the loop
    peer_already_run_to_time = PETSC_FALSE
    if (this%is_master .and. &
        (checkpoint_at_this_time_flag .or. &
         checkpoint_at_this_timestep_flag)) then
      ! if checkpointing, need to sync all other PMCs.  Those "below" are
      ! already in sync, but not those "next".
      ! Set data needed by process-model
      call this%SetAuxData()
      ! Run neighboring process model couplers
      if (associated(this%peer)) then
        call this%peer%RunToTime(this%timestepper%target_time,local_stop_flag)
        peer_already_run_to_time = PETSC_TRUE
      endif
      call this%GetAuxData()
      ! it is possible that two identical checkpoint files will be created,
      ! one at the time and another at the time step, but this is fine.
      if (checkpoint_at_this_time_flag) then
        filename_append = &
          CheckpointAppendNameAtTime(this%checkpoint_option, &
                                     this%option%time,this%option)
        call this%Checkpoint(filename_append)
      endif
      if (checkpoint_at_this_timestep_flag) then
        filename_append = &
          CheckpointAppendNameAtTimestep(this%checkpoint_option, &
                                         this%timestepper%steps, &
                                         this%option)
        call this%Checkpoint(filename_append)
      endif
    endif
    
    if (this%is_master) then
      if (this%timestepper%WallClockStop(this%option)) then
         local_stop_flag = TS_STOP_WALLCLOCK_EXCEEDED
      endif
    endif

  enddo
  
  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer) .and. .not.peer_already_run_to_time) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif
  
  stop_flag = max(stop_flag,local_stop_flag)
  
  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMCBaseRunToTime

! ************************************************************************** !

recursive subroutine PMCBaseUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%UpdateSolution()')
#endif
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    cur_pm => cur_pm%next
  enddo  
  
end subroutine PMCBaseUpdateSolution

! ************************************************************************** !

recursive subroutine FinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  character(len=MAXSTRINGLENGTH) :: string
  
#ifdef DEBUG
  call printMsg(this%option,'PMCBase%FinalizeRun()')
#endif
  
  if (associated(this%timestepper)) then
    call this%timestepper%FinalizeRun(this%option)
  endif

  if (associated(this%child)) then
    call this%child%FinalizeRun()
  endif
  
  if (associated(this%peer)) then
    call this%peer%FinalizeRun()
  endif
  
end subroutine FinalizeRun

! ************************************************************************** !

subroutine SetOutputFlags(this)
  ! 
  ! Toggles flags that determine whether output is printed
  ! to the screen and output file during a time step.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/29/13
  ! 

  use Option_module
  use Output_Aux_module
  
  implicit none
  
  class(pmc_base_type) :: this
  
  type(output_option_type), pointer :: output_option
  
  output_option => this%pm_list%output_option

  if (OptionPrintToScreen(this%option) .and. &
      mod(this%timestepper%steps,output_option%screen_imod) == 0) then
    this%option%print_screen_flag = PETSC_TRUE
  else
    this%option%print_screen_flag = PETSC_FALSE
  endif

  if (OptionPrintToFile(this%option) .and. &
      mod(this%timestepper%steps,output_option%output_file_imod) == 0) then
    this%option%print_file_flag = PETSC_TRUE
  else
    this%option%print_file_flag = PETSC_FALSE
      
  endif
  
end subroutine SetOutputFlags

! ************************************************************************** !

recursive subroutine OutputLocal(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this
  
  class(pm_base_type), pointer :: cur_pm
  
#ifdef DEBUG
  call printMsg(this%option,'PMC%Output()')
#endif
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
!    call Output(cur_pm%realization,snapshot_plot_flag,observation_plot_flag, &
!                massbal_plot_flag)
    cur_pm => cur_pm%next
  enddo
    
end subroutine OutputLocal

! ************************************************************************** !

recursive subroutine PMCBaseCheckpoint(this,filename_append)
  ! 
  ! Checkpoints PMC timestepper and state variables.
  ! 
  ! Author: Glenn Hammond
  ! Date: 2/2/16
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none

  class(pmc_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: filename_append
  
  PetscInt :: chk_grp_id
  PetscViewer :: viewer
  
  if (this%checkpoint_option%format == CHECKPOINT_BINARY .or. &
      this%checkpoint_option%format == CHECKPOINT_BOTH) then
    call this%CheckpointBinary(viewer,filename_append)
  endif
  if (this%checkpoint_option%format == CHECKPOINT_HDF5 .or. &
      this%checkpoint_option%format == CHECKPOINT_BOTH) then
#if !defined(PETSC_HAVE_HDF5)
    this%option%io_buffer = 'HDF5 formatted checkpointing not supported &
      &unless PFLOTRAN is compiled with HDF5 libraries enabled.'
    call printErrMsg(this%option)
#else
    call this%CheckpointHDF5(chk_grp_id,filename_append)
#endif
  endif

end subroutine PMCBaseCheckpoint

! ************************************************************************** !

recursive subroutine PMCBaseCheckpointBinary(this,viewer,append_name)
  ! 
  ! Checkpoints PMC timestepper and state variables.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module
  use Checkpoint_module, only : CheckpointOpenFileForWriteBinary, &
                                CheckPointWriteCompatibilityBinary

  implicit none

  class(pmc_base_type) :: this
  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: append_name
  
  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  type(pmc_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  ! if the top PMC
  if (this%is_master) then
    call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
    call PetscLogEventBegin(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)
    call CheckpointOpenFileForWriteBinary(viewer,append_name,this%option)
    call CheckPointWriteCompatibilityBinary(viewer,this%option)
    ! create header for storing local information specific to PMc
    call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
    call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
    call PMCBaseRegisterHeader(this,bag,header)
    call PMCBaseSetHeader(this,bag,header)
    call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
    call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  endif
  
  if (associated(this%timestepper)) then
    call this%timestepper%CheckpointBinary(viewer,this%option)
  endif
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%CheckpointBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%child)) then
    call this%child%CheckpointBinary(viewer,append_name)
  endif
  
  if (associated(this%peer)) then
    call this%peer%CheckpointBinary(viewer,append_name)
  endif
  
  if (this%is_master) then
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    viewer = PETSC_NULL_VIEWER
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to write to checkpoint file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
    
end subroutine PMCBaseCheckpointBinary

! ************************************************************************** !

subroutine PMCBaseRegisterHeader(this,bag,header)
  ! 
  ! Register header entries.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  ! bagsize = 2 * 8 bytes = 16 bytes
  call PetscBagRegisterInt(bag,header%plot_number,0, &
                           "plot number","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%times_per_h5_file,0, &
                           "times_per_h5_file","",ierr);CHKERRQ(ierr)

end subroutine PMCBaseRegisterHeader

! ************************************************************************** !

subroutine PMCBaseSetHeader(this,bag,header)
  ! 
  ! Sets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  PetscBag :: bag
  
  PetscErrorCode :: ierr
  
  header%plot_number = &
    this%pm_list%realization_base%output_option%plot_number

  header%times_per_h5_file = &
    this%pm_list%realization_base%output_option%times_per_h5_file

end subroutine PMCBaseSetHeader

! ************************************************************************** !

recursive subroutine PMCBaseRestartBinary(this,viewer)
  ! 
  ! Restarts PMC timestepper and state variables.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/26/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Logging_module
  use Checkpoint_module, only : CheckPointReadCompatibilityBinary

  implicit none
  

  class(pmc_base_type) :: this
  PetscViewer :: viewer

  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  type(pmc_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  ! if the top PMC, 
  if (this%is_master) then
    this%option%io_buffer = 'Restarting with checkpoint file "' // &
      trim(this%option%restart_filename) // '".'
    call printMsg(this%option)
    call PetscLogEventBegin(logging%event_restart,ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)
    call PetscViewerBinaryOpen(this%option%mycomm, &
                               this%option%restart_filename, &
                               FILE_MODE_READ,viewer,ierr);CHKERRQ(ierr)
    ! skip reading info file when loading, but not working
    call PetscViewerBinarySetSkipOptions(viewer,PETSC_TRUE,ierr);CHKERRQ(ierr)
    call CheckPointReadCompatibilityBinary(viewer,this%option)
    ! read pmc header
    call PetscBagCreate(this%option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
    call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
    call PMCBaseRegisterHeader(this,bag,header)
    call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
    call PMCBaseGetHeader(this,header)
    if (Initialized(this%option%restart_time)) then
      this%pm_list%realization_base%output_option%plot_number = 0
    endif
    call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  endif
  
  if (associated(this%timestepper)) then
    call this%timestepper%RestartBinary(viewer,this%option)
    if (Initialized(this%option%restart_time)) then
      ! simply a flag to set time back to zero, no matter what the restart
      ! time is set to.
      call this%timestepper%Reset()
      ! note that this sets the target time back to zero.
    endif
  
    ! Point cur_waypoint to the correct waypoint.
    !geh: there is a problem here in that the timesteppers "prev_waypoint"
    !     may not be set correctly if the time step does not converge. See
    !     top of TimestepperBaseSetTargetTime().
    call WaypointSkipToTime(this%timestepper%cur_waypoint, &
                            this%timestepper%target_time)
    !geh: this is a bit of a kludge.  Need to use the timestepper target time
    !     directly.  Time is needed to update boundary conditions within 
    !     this%UpdateSolution
    ! check to ensure that simulation is not restarted beyond the end of the
    ! prescribed final simulation time.
    if (.not.associated(this%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%timestepper%target_time* &
        this%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%pm_list%realization_base%output_option%tunit) // ').'
      call printErrMsg(this%option)
    endif
    this%option%time = this%timestepper%target_time
  endif
  
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call cur_pm%RestartBinary(viewer)
    cur_pm => cur_pm%next
  enddo
  
  if (associated(this%child)) then
    call this%child%RestartBinary(viewer)
  endif
  
  if (associated(this%peer)) then
    call this%peer%RestartBinary(viewer)
  endif
  
  if (this%is_master) then
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to read from restart file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_restart,ierr);CHKERRQ(ierr)
  endif
    
end subroutine PMCBaseRestartBinary

! ************************************************************************** !

subroutine PMCBaseGetHeader(this,header)
  ! 
  ! Gets values in checkpoint header.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/02/13
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module

  implicit none

  class(pmc_base_type) :: this
  class(pmc_base_header_type) :: header
  
  character(len=MAXSTRINGLENGTH) :: string

  this%pm_list%realization_base%output_option%plot_number = &
    header%plot_number

  ! Check the value of 'times_per_h5_file'
  if (header%times_per_h5_file /= &
      this%pm_list%realization_base%output_option%times_per_h5_file) then
    write(string,*) header%times_per_h5_file
    this%option%io_buffer = 'From checkpoint file: times_per_h5_file ' // trim(string)
    call printMsg(this%option)
    write(string,*) this%pm_list%realization_base%output_option%times_per_h5_file
    this%option%io_buffer = 'From inputdeck      : times_per_h5_file ' // trim(string)
    call printMsg(this%option)
    this%option%io_buffer = 'times_per_h5_file specified in inputdeck does not ' // &
      'match that stored in checkpoint file. Correct the inputdeck.'
    call printErrMsg(this%option)
  endif

  this%pm_list%realization_base%output_option%times_per_h5_file = &
    header%times_per_h5_file

end subroutine PMCBaseGetHeader

! ************************************************************************** !

#if defined(PETSC_HAVE_HDF5)
recursive subroutine PMCBaseCheckpointHDF5(this,chk_grp_id,append_name)
  !
  ! Checkpoints PMC timestepper and state variables in HDF5 format.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/29/15
  !
  use Logging_module
  use Checkpoint_module, only : CheckpointOpenFileForWriteHDF5, &
                                CheckPointWriteCompatibilityHDF5
  use hdf5

  implicit none

  class(pmc_base_type) :: this
  PetscInt :: chk_grp_id
  character(len=MAXSTRINGLENGTH) :: append_name

#if defined(SCORPIO_WRITE)
  integer :: h5_chk_grp_id
  integer :: h5_file_id
  integer :: pmc_grp_id
  integer :: pm_grp_id
  integer :: temp_id
#else
  integer(HID_T) :: h5_chk_grp_id
  integer(HID_T) :: h5_file_id
  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id
#endif

  class(pm_base_type), pointer :: cur_pm
  class(pmc_base_header_type), pointer :: header
  type(pmc_base_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscInt :: pmc_grp_id
  PetscSizeT :: bagsize
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err

  bagsize = size(transfer(dummy_header,dummy_char))

  ! if the top PMC
  if (this%is_master) then
    call PetscLogStagePush(logging%stage(OUTPUT_STAGE),ierr);CHKERRQ(ierr)
    call PetscLogEventBegin(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)
    call CheckpointOpenFileForWriteHDF5(h5_file_id, &
                                        h5_chk_grp_id, &
                                        append_name, this%option)
    call CheckPointWriteCompatibilityHDF5(h5_chk_grp_id, &
                                          this%option)
    call h5gcreate_f(h5_chk_grp_id, trim(this%name), &
                     h5_pmc_grp_id,hdf5_err, OBJECT_NAMELEN_DEFAULT_F)
    call PMCBaseSetHeaderHDF5(this, h5_pmc_grp_id, this%option)
    chk_grp_id = h5_chk_grp_id
  else
    h5_chk_grp_id = chk_grp_id
    call h5gcreate_f(h5_chk_grp_id, trim(this%name), &
                     h5_pmc_grp_id, hdf5_err, OBJECT_NAMELEN_DEFAULT_F)
  endif

  if (associated(this%timestepper)) then
    pmc_grp_id = h5_pmc_grp_id
    call this%timestepper%CheckpointHDF5(pmc_grp_id, this%option)
  endif

  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit

    call h5gcreate_f(h5_pmc_grp_id, trim(cur_pm%name), h5_pm_grp_id, &
         hdf5_err, OBJECT_NAMELEN_DEFAULT_F)
    call cur_pm%CheckpointHDF5(h5_pm_grp_id)
    call h5gclose_f(h5_pm_grp_id, hdf5_err)

    cur_pm => cur_pm%next
  enddo

  call h5gclose_f(h5_pmc_grp_id, hdf5_err)

  if (associated(this%child)) then
    call this%child%CheckpointHDF5(chk_grp_id,append_name)
  endif

  if (associated(this%peer)) then
    call this%peer%CheckpointHDF5(chk_grp_id,append_name)
  endif

  if (this%is_master) then
    call h5gclose_f(h5_chk_grp_id, hdf5_err)
    call h5fclose_f(h5_file_id,hdf5_err)
    call h5close_f(hdf5_err)
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to write to checkpoint file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_checkpoint,ierr);CHKERRQ(ierr)
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif

end subroutine PMCBaseCheckpointHDF5

! ************************************************************************** !

recursive subroutine PMCBaseRestartHDF5(this,chk_grp_id)
  ! 
  ! Restarts PMC timestepper and state variables from a HDF5
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/09/15
  ! 
  use Logging_module
  use hdf5
  use Checkpoint_module, only : CheckPointReadCompatibilityHDF5, &
                                CheckpointOpenFileForReadHDF5

  implicit none

  class(pmc_base_type) :: this
  PetscInt :: chk_grp_id

  class(pm_base_type), pointer :: cur_pm
  PetscLogDouble :: tstart, tend
  PetscErrorCode :: ierr
  PetscMPIInt :: hdf5_err

#if defined(SCORPIO_WRITE)
  integer :: h5_chk_grp_id
  integer :: h5_file_id
  integer :: h5_pmc_grp_id
  integer :: h5_pm_grp_id
#else
  integer(HID_T) :: h5_chk_grp_id
  integer(HID_T) :: h5_file_id
  integer(HID_T) :: h5_pmc_grp_id
  integer(HID_T) :: h5_pm_grp_id
#endif
  PetscInt :: pmc_grp_id

  ! if the top PMC
  if (this%is_master) then

    this%option%io_buffer = 'Restarting with checkpoint file "' // &
      trim(this%option%restart_filename) // '".'
    call printMsg(this%option)
    call PetscLogEventBegin(logging%event_restart, ierr);CHKERRQ(ierr)
    call PetscTime(tstart,ierr);CHKERRQ(ierr)

    call CheckpointOpenFileForReadHDF5(this%option%restart_filename, &
                                       h5_file_id, &
                                       h5_chk_grp_id, &
                                       this%option)

    call CheckPointReadCompatibilityHDF5(h5_chk_grp_id, &
                                         this%option)

    call h5gopen_f(h5_chk_grp_id, trim(this%name), &
                   h5_pmc_grp_id, hdf5_err)

    ! read pmc header
    call PMCBaseGetHeaderHDF5(this, h5_pmc_grp_id, this%option)

    if (Initialized(this%option%restart_time)) then
      this%pm_list%realization_base%output_option%plot_number = 0
    endif
    chk_grp_id = h5_chk_grp_id 
  else
    h5_chk_grp_id = chk_grp_id
    call h5gopen_f(h5_chk_grp_id, trim(this%name), &
                   h5_pmc_grp_id, hdf5_err)

  endif

  if (associated(this%timestepper)) then
    pmc_grp_id = h5_pmc_grp_id
    call this%timestepper%RestartHDF5(pmc_grp_id, this%option)

    if (Initialized(this%option%restart_time)) then
      ! simply a flag to set time back to zero, no matter what the restart
      ! time is set to.
      call this%timestepper%Reset()
      ! note that this sets the target time back to zero.
    endif

    ! Point cur_waypoint to the correct waypoint.
    !geh: there is a problem here in that the timesteppers "prev_waypoint"
    !     may not be set correctly if the time step does not converge. See
    !     top of TimestepperBaseSetTargetTime().
    call WaypointSkipToTime(this%timestepper%cur_waypoint, &
                            this%timestepper%target_time)
    !geh: this is a bit of a kludge.  Need to use the timestepper target time
    !     directly.  Time is needed to update boundary conditions within
    !     this%UpdateSolution
    ! check to ensure that simulation is not restarted beyond the end of the
    ! prescribed final simulation time.
    if (.not.associated(this%timestepper%cur_waypoint)) then
      write(this%option%io_buffer,*) this%timestepper%target_time/ &
        this%pm_list%realization_base%output_option%tconv
      this%option%io_buffer = 'Simulation is being restarted at a time that &
        &is at or beyond the end of checkpointed simulation (' // &
        trim(adjustl(this%option%io_buffer)) // &
        trim(this%pm_list%realization_base%output_option%tunit) // ').'
      call printErrMsg(this%option)
    endif
    this%option%time = this%timestepper%target_time
  endif

  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    call h5gopen_f(h5_pmc_grp_id, trim(cur_pm%name), h5_pm_grp_id, &
                   hdf5_err)
    call cur_pm%RestartHDF5(h5_pm_grp_id)
    call h5gclose_f(h5_pm_grp_id, hdf5_err)
    cur_pm => cur_pm%next
  enddo

  call h5gclose_f(h5_pmc_grp_id, hdf5_err)

  if (associated(this%child)) then
    call this%child%RestartHDF5(chk_grp_id)
  endif

  if (associated(this%peer)) then
    call this%peer%RestartHDF5(chk_grp_id)
  endif

  if (this%is_master) then
    call h5gclose_f(h5_chk_grp_id, hdf5_err)
    call h5fclose_f(h5_file_id, hdf5_err)
    call h5close_f(hdf5_err)
    call PetscTime(tend,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer, &
          '("      Seconds to read from restart file: ", f10.2)') &
      tend-tstart
    call printMsg(this%option)
    call PetscLogEventEnd(logging%event_restart,ierr);CHKERRQ(ierr)
  endif

end subroutine PMCBaseRestartHDF5

! ************************************************************************** !

subroutine PMCBaseSetHeaderHDF5(this, chk_grp_id, option)
  ! 
  ! Similar to PMCBaseSetHeader(), except this subroutine writes values in
  ! a HDF5.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/30/15
  ! 

  use Option_module
  use Checkpoint_module, only: CheckPointWriteIntDatasetHDF5
  use hdf5

  implicit none

  class(pmc_base_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
#else
  integer(HID_T) :: chk_grp_id
#endif
  type(option_type) :: option
  
#if defined(SCORPIO_WRITE)
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Output_plot_number" // CHAR(0)
  int_array(1) = this%pm_list%realization_base%output_option%plot_number
  call CheckPointWriteIntDatasetHDF5(chk_grp_id, &
                                     dataset_name, dataset_rank, &
                                     dims, start, length, stride, &
                                     int_array, option)

  dataset_name = "Output_times_per_h5_file" // CHAR(0)
  int_array(1) = this%pm_list%realization_base%output_option%times_per_h5_file
  call CheckPointWriteIntDatasetHDF5(chk_grp_id, &
                                     dataset_name, dataset_rank, &
                                     dims, start, length, stride, &
                                     int_array, option)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

end subroutine PMCBaseSetHeaderHDF5

! ************************************************************************** !

subroutine PMCBaseGetHeaderHDF5(this, chk_grp_id, option)
  ! 
  ! Similar to PMCBaseGetHeader(), except this subroutine reads values from
  ! a HDF5.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/16/15
  ! 

  use Option_module
  use Checkpoint_module, only: CheckPointReadIntDatasetHDF5
  use hdf5

  implicit none

  class(pmc_base_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: chk_grp_id
#else
  integer(HID_T) :: chk_grp_id
#endif
  type(option_type) :: option

#if defined(SCORPIO_WRITE)
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Output_plot_number" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
                                    dims, start, length, stride, &
                                    int_array, option)
  this%pm_list%realization_base%output_option%plot_number = int_array(1)

  dataset_name = "Output_times_per_h5_file" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(chk_grp_id, dataset_name, dataset_rank, &
                                    dims, start, length, stride, &
                                    int_array, option)
  this%pm_list%realization_base%output_option%times_per_h5_file = int_array(1)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

end subroutine PMCBaseGetHeaderHDF5
#endif

! ************************************************************************** !

subroutine AccumulateAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine AccumulateAuxData

! ************************************************************************** !

subroutine GetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine GetAuxData

! ************************************************************************** !

subroutine SetAuxData(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht,LBNL
  ! Date: 08/21/13
  ! 

  implicit none
  
  class(pmc_base_type) :: this

end subroutine SetAuxData

! ************************************************************************** !

subroutine PMCBaseUpdateMaterialProperties(this)
  !
  ! At a prescribed time, updates material properties based on instructions
  ! provided by a material update waypoint.
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/14
  
  implicit none
  
  class(pmc_base_type) :: this

end subroutine PMCBaseUpdateMaterialProperties

! ************************************************************************** !

recursive subroutine PMCBaseCheckNullPM(this,option)
  ! 
  ! This routine
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/10/14
  ! 
  use Option_module
  
  implicit none
  
  class(pmc_base_type) :: this
  type(option_type) :: option
  
  if (.not.associated(this%pm_list)) then
    option%io_buffer = 'Null PM in PMC "' // trim(this%name) // '".'
    call printErrMsg(option)
  endif
  
  if (associated(this%peer)) then
    call this%peer%CheckNullPM(option)
  endif
  
  if (associated(this%child)) then
    call this%child%CheckNullPM(option)
  endif

end subroutine PMCBaseCheckNullPM

! ************************************************************************** !

subroutine PMCBaseStrip(this)
  !
  ! Deallocates members of PMC Base.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_base_type) :: this

  ! these are destoyed elsewhere
  nullify(this%option)
  nullify(this%checkpoint_option)
  nullify(this%waypoint_list)

  if (associated(this%timestepper)) then
    call this%timestepper%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%timestepper)
    nullify(this%timestepper)
  endif
  if (associated(this%pm_list)) then
    ! destroy does not currently destroy; it strips
    call this%pm_list%Destroy()
    deallocate(this%pm_list)
    nullify(this%pm_list)
  endif
  if (associated(this%pm_ptr)) then
    nullify(this%pm_ptr%pm) ! solely a pointer
    deallocate(this%pm_ptr)
    nullify(this%pm_ptr)
  endif
  nullify(this%sim_aux)

end subroutine PMCBaseStrip

! ************************************************************************** !

recursive subroutine PMCBaseDestroy(this)
  ! 
  ! Deallocates a pmc object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_base_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMC%Destroy()')
#endif
  
  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif 
  
!  deallocate(pmc)
!  nullify(pmc)
  
end subroutine PMCBaseDestroy
  
end module PMC_Base_class
