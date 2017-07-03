module Simulation_Surface_class

  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Surface_class
  use PMC_Base_class
  use Realization_Surface_class
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public, extends(simulation_base_type) :: simulation_surface_type
    class(pmc_surface_type), pointer :: surf_flow_process_model_coupler
    class(realization_surface_type), pointer :: surf_realization
    type(regression_type), pointer :: regression
    type(waypoint_list_type), pointer :: waypoint_list_surface
  contains
    procedure, public :: Init => SurfaceSimulationInit
    procedure, public :: InputRecord => SurfaceSimInputRecord
    procedure, public :: InitializeRun => SurfaceInitializeRun
    procedure, public :: FinalizeRun => SurfaceFinalizeRun
    procedure, public :: Strip => SurfaceSimulationStrip
  end type simulation_surface_type

  public :: SurfaceSimulationCreate, &
            SurfaceSimulationInit, &
            SurfaceFinalizeRun, &
            SurfaceSimulationStrip, &
            SurfaceSimulationDestroy

contains

! ************************************************************************** !

function SurfaceSimulationCreate(option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(simulation_surface_type), pointer :: SurfaceSimulationCreate
  
  print *, 'SurfaceSimulationCreate'
  
  allocate(SurfaceSimulationCreate)
  call SurfaceSimulationCreate%Init(option)
  
end function SurfaceSimulationCreate

! ************************************************************************** !

subroutine SurfaceSimulationInit(this,option)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 
  use Waypoint_module
  use Option_module
  
  implicit none
  
  class(simulation_surface_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%regression)
  this%waypoint_list_surface => WaypointListCreate()
  
end subroutine SurfaceSimulationInit

! ************************************************************************** !

subroutine SurfaceSimInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  use Output_module
  
  implicit none
  
  class(simulation_surface_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
 
  write(id,'(a29)',advance='no') 'simulation type: '
  write(id,'(a)') 'surface'

  ! print output file information
  call OutputInputRecord(this%output_option,this%waypoint_list_surface)

end subroutine SurfaceSimInputRecord

! ************************************************************************** !

subroutine SurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Logging_module
  use Output_module

  implicit none
  
  class(simulation_surface_type) :: this

  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: cur_process_model_coupler_top
  class(pmc_base_type), pointer :: cur_process_model_coupler_below
  PetscInt :: depth
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'SurfaceInitializeRun: Simulation%InitializeRun()')

  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    depth = 0
    call cur_process_model_coupler%InitializeRun()
    cur_process_model_coupler => cur_process_model_coupler%peer
  enddo

  ! set depth in tree
  cur_process_model_coupler_top => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler_top)) exit
    depth = 0
    cur_process_model_coupler_below => cur_process_model_coupler_top%child
    do
      if (.not.associated(cur_process_model_coupler_below)) exit
      depth = depth + 1
      cur_process_model_coupler_below => cur_process_model_coupler_below%child
    enddo
    cur_process_model_coupler_top => cur_process_model_coupler_top%peer
  enddo

end subroutine SurfaceInitializeRun

! ************************************************************************** !

subroutine SurfaceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  use Timestepper_Base_class

  implicit none
  
  class(simulation_surface_type) :: this
  
  PetscErrorCode :: ierr

  class(timestepper_base_type), pointer :: surf_flow_timestepper

  call printMsg(this%option,'SurfaceFinalizeRun()')
  
  call SimulationBaseFinalizeRun(this)
  
  nullify(surf_flow_timestepper)
  surf_flow_timestepper => this%surf_flow_process_model_coupler%timestepper

end subroutine SurfaceFinalizeRun

! ************************************************************************** !

subroutine SurfaceSimulationStrip(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(simulation_surface_type) :: this
  
  call printMsg(this%option,'SurfaceSimulationStrip()')
  
  call SimulationBaseStrip(this)
  call RealizSurfStrip(this%surf_realization)
  deallocate(this%surf_realization)
  nullify(this%surf_realization)
  call RegressionDestroy(this%regression)
  call WaypointListDestroy(this%waypoint_list_surface)
  
end subroutine SurfaceSimulationStrip

! ************************************************************************** !

subroutine SurfaceSimulationDestroy(simulation)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/27/13
  ! 

  implicit none
  
  class(simulation_surface_type), pointer :: simulation
  
  call printMsg(simulation%option,'SurfaceSimulationDestroy()')
  
  if (.not.associated(simulation)) return
  
  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)
  
end subroutine SurfaceSimulationDestroy

end module Simulation_Surface_class
