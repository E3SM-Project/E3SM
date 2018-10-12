module Simulation_Surf_Subsurf_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Simulation_Surface_class
  use Simulation_Subsurface_class
  use Regression_module
  use Option_module
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Surface_class
  use Realization_Subsurface_class
  use Realization_Surface_class
  use Waypoint_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none

  private

  type, public, extends(simulation_subsurface_type) :: &
    simulation_surfsubsurface_type
    class(pmc_surface_type), pointer :: surf_flow_process_model_coupler
    class(realization_surface_type), pointer :: surf_realization
    type(waypoint_list_type), pointer :: waypoint_list_surfsubsurface
  contains
    procedure, public :: Init => SurfSubsurfaceSimulationInit
    procedure, public :: InitializeRun => SurfSubsurfaceInitializeRun
    procedure, public :: InputRecord => SurfSubsurfaceInputRecord
    procedure, public :: FinalizeRun => SurfSubsurfaceFinalizeRun
    procedure, public :: Strip => SurfSubsurfaceSimulationStrip
    procedure, public :: ExecuteRun => SurfSubsurfaceExecuteRun
    procedure, public :: RunToTime => SurfSubsurfaceSimulationRunToTime
  end type simulation_surfsubsurface_type

  public :: SurfSubsurfaceSimulationCreate, &
            SurfSubsurfaceSimulationInit, &
            SurfSubsurfaceFinalizeRun, &
            SurfSubsurfaceSimulationStrip, &
            SurfSubsurfaceSimulationDestroy

contains

! ************************************************************************** !

function SurfSubsurfaceSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(simulation_surfsubsurface_type), pointer :: SurfSubsurfaceSimulationCreate
  
  print *, 'SurfSubsurfaceSimulationCreate'
  
  allocate(SurfSubsurfaceSimulationCreate)
  call SurfSubsurfaceSimulationCreate%Init(option)
  
end function SurfSubsurfaceSimulationCreate

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationInit(this,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 
  use Waypoint_module
  use Option_module
  
  implicit none
  
  class(simulation_surfsubsurface_type) :: this
  type(option_type), pointer :: option
  
  call SubsurfaceSimulationInit(this,option)
  nullify(this%surf_realization)
  this%waypoint_list_surfsubsurface => WaypointListCreate()
  
end subroutine SurfSubsurfaceSimulationInit

! ************************************************************************** !

subroutine SurfSubsurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 
#include "petsc/finclude/petscviewer.h"
  use petscsys
  use Logging_module
  use Output_module
  use PMC_Surface_class

  implicit none
  
  class(simulation_surfsubsurface_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  PetscViewer :: viewer
  
  call printMsg(this%option,'Simulation%InitializeRun()')

  call this%process_model_coupler_list%InitializeRun()

  if (this%option%restart_flag) then
    call this%process_model_coupler_list%RestartBinary(viewer)
    cur_process_model_coupler => this%process_model_coupler_list
    select type(pmc => cur_process_model_coupler)
      class is(pmc_surface_type)
        select case(this%option%iflowmode)
          case (TH_MODE)
            call pmc%PMCSurfaceGetAuxDataAfterRestart()
          case default
            call printErrMsg(this%option,'SurfSubsurfaceInitializeRun ' // &
                  'not supported in current flow mode.')
        end select
    end select

  endif

end subroutine SurfSubsurfaceInitializeRun

! ************************************************************************** !

subroutine SurfSubsurfaceInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  use Output_module

  implicit none
  
  class(simulation_surfsubsurface_type) :: this
  
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
 
  write(id,'(a29)',advance='no') 'simulation type: '
  write(id,'(a)') 'surface-subsurface'

  ! print output file information
  call OutputInputRecord(this%output_option,this%waypoint_list_surfsubsurface)

end subroutine SurfSubsurfaceInputRecord

! ************************************************************************** !

subroutine SurfSubsurfaceExecuteRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Simulation_Base_class
  use Timestepper_Base_class, only : TS_CONTINUE
  use Checkpoint_module

  implicit none
  
  class(simulation_surfsubsurface_type) :: this

  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt
  character(len=MAXSTRINGLENGTH) :: append_name

  time = 0.d0
  time = this%option%time

  final_time = SimulationGetFinalWaypointTime(this)
  append_name = '-restart'

  call printMsg(this%option,'SurfSubsurfaceExecuteRun()')

  if (.not.associated(this%surf_realization)) then
    call this%RunToTime(final_time)

  else

    ! If simulation is decoupled surface-subsurface simulation, set
    ! dt_coupling to be dt_max
    if (Equal(this%surf_realization%dt_coupling,0.d0)) &
      this%surf_realization%dt_coupling = this%surf_realization%dt_max

    do
      if (time + this%surf_realization%dt_coupling > final_time) then
        dt = final_time-time
      else
        dt = this%surf_realization%dt_coupling
      endif

      time = time + dt
      call this%RunToTime(time)

      if (this%stop_flag /= TS_CONTINUE) exit ! end simulation

      if (time >= final_time) exit
    enddo

  endif
  if (associated(this%process_model_coupler_list%checkpoint_option)) then
    append_name = CheckpointFilename(append_name,this%option)
    call this%process_model_coupler_list%Checkpoint(append_name)
  endif

end subroutine SurfSubsurfaceExecuteRun

! ************************************************************************** !

subroutine SurfSubsurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  use Simulation_Base_class
  use Timestepper_Base_class

  implicit none
  
  class(simulation_surfsubsurface_type) :: this
  
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'SurfSubsurfaceFinalizeRun()')
  
  call SubsurfaceFinalizeRun(this)
  !call SurfaceFinalizeRun(this)
  
end subroutine SurfSubsurfaceFinalizeRun

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 
  use Waypoint_module
  use Simulation_Base_class

  implicit none
  
  class(simulation_surfsubsurface_type) :: this
  
  call printMsg(this%option,'SurfSubsurfaceSimulationStrip()')
  
  call SubsurfaceSimulationStrip(this)
  call RealizSurfStrip(this%surf_realization)
  deallocate(this%surf_realization)
  nullify(this%surf_realization)
  call WaypointListDestroy(this%waypoint_list_surfsubsurface)
 
end subroutine SurfSubsurfaceSimulationStrip

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationRunToTime(this,target_time)
  ! 
  ! This routine executes surface-subsurface simualation
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

#include "petsc/finclude/petscviewer.h"
  use petscsys
  use Option_module
  use Simulation_Aux_module

  implicit none

  class(simulation_surfsubsurface_type) :: this
  PetscReal :: target_time

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscViewer :: viewer

#ifdef DEBUG
  call printMsg(this%option,'RunToTime()')
#endif
  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine SurfSubsurfaceSimulationRunToTime

! ************************************************************************** !

subroutine SurfSubsurfaceSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/28/13
  ! 

  implicit none
  
  class(simulation_surfsubsurface_type), pointer :: simulation
  
  call printMsg(simulation%option,'SimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SurfSubsurfaceSimulationDestroy

end module Simulation_Surf_Subsurf_class
