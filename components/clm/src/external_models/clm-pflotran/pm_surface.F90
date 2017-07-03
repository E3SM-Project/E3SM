module PM_Surface_class

  use PM_Base_class
  use Realization_Surface_class
  use Communicator_Base_module
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscts.h"

  type, public, extends(pm_base_type) :: pm_surface_type
    class(realization_surface_type), pointer :: surf_realization
    class(communicator_type), pointer :: comm1
    PetscReal :: pressure_change_governor
    PetscReal :: temperature_change_governor
    PetscReal :: pressure_dampening_factor
    PetscReal :: pressure_change_limit
    PetscReal :: temperature_change_limit
  contains
    procedure, public :: Setup => PMSurfaceSetup
    procedure, public :: PMSurfaceSetRealization
    procedure, public :: InitializeRun => PMSurfaceInitializeRun
    procedure, public :: PreSolve => PMSurfacePreSolve
    procedure, public :: PostSolve => PMSurfacePostSolve
    procedure, public :: CheckpointBinary => PMSurfaceCheckpointBinary
    procedure, public :: RestartBinary => PMSurfaceRestartBinary
    procedure, public :: UpdateAuxVars => PMSurfaceUpdateAuxVars
    procedure, public :: InputRecord => PMSurfaceInputRecord
  end type pm_surface_type

  public :: PMSurfaceCreate, &
            PMSurfaceSetup, &
            PMSurfaceUpdateSolution, &
            PMSurfaceReadSelectCase, &
            PMSurfaceDestroy
  
contains

! ************************************************************************** !

subroutine PMSurfaceCreate(this)
  ! 
  ! Intializes shared members of surface process models
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  implicit none
  
  class(pm_surface_type) :: this
  
  this%pressure_change_governor = 5.d5
  this%temperature_change_governor = 5.d0
  this%pressure_dampening_factor = UNINITIALIZED_DOUBLE
  this%pressure_change_limit = UNINITIALIZED_DOUBLE
  this%temperature_change_limit = UNINITIALIZED_DOUBLE

  nullify(this%surf_realization)
  nullify(this%comm1)
  
  call PMBaseInit(this)

end subroutine PMSurfaceCreate

! ************************************************************************** !

subroutine PMSurfaceReadSelectCase(this,input,keyword,found,option)
  ! 
  ! Reads input file parameters associated with the subsurface flow process 
  !       model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/05/16

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_surface_type) :: this
  type(input_type) :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  type(option_type) :: option

  found = PETSC_TRUE
  select case(trim(keyword))

    case('MAX_PRESSURE_CHANGE')
      call InputReadDouble(input,option,this%pressure_change_governor)
      call InputDefaultMsg(input,option,'dpmxe')

    case('MAX_TEMPERATURE_CHANGE')
      call InputReadDouble(input,option,this%temperature_change_governor)
      call InputDefaultMsg(input,option,'dtmpmxe')

    case('PRESSURE_DAMPENING_FACTOR')
      call InputReadDouble(input,option,this%pressure_dampening_factor)
      call InputErrorMsg(input,option,'PRESSURE_DAMPENING_FACTOR', &
                          'TIMESTEPPER')

    case('PRESSURE_CHANGE_LIMIT')
      call InputReadDouble(input,option,this%pressure_change_limit)
      call InputErrorMsg(input,option,'PRESSURE_CHANGE_LIMIT', &
                          'TIMESTEPPER')

    case('TEMPERATURE_CHANGE_LIMIT')
      call InputReadDouble(input,option,this%temperature_change_limit)
      call InputErrorMsg(input,option,'TEMPERATURE_CHANGE_LIMIT', &
                          'TIMESTEPPER')
    case default
      found = PETSC_FALSE
  end select

end subroutine PMSurfaceReadSelectCase

! ************************************************************************** !

subroutine PMSurfaceSetup(this)
  ! 
  ! Initializes variables associated with subsurface process models
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Discretization_module
  use Communicator_Unstructured_class
  use Grid_module

  implicit none

  class(pm_surface_type) :: this

  ! set up communicator
  select case(this%surf_realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%option%io_buffer='Surface flow not supported on structured grids'
      call printErrMsg(this%option)
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  ! set the communicator
  call this%comm1%SetDM(this%surf_realization%discretization%dm_1dof)

end subroutine PMSurfaceSetup

! ************************************************************************** !

subroutine PMSurfaceSetRealization(this, surf_realization)
  ! 
  ! Initializes relization and PETSc vectors for solution and residual.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Realization_Surface_class
  use Grid_module

  implicit none

  class(pm_surface_type) :: this
  class(realization_surface_type), pointer :: surf_realization

  this%surf_realization => surf_realization
  this%realization_base => surf_realization

  this%solution_vec = surf_realization%surf_field%flow_xx
  this%residual_vec = surf_realization%surf_field%flow_r

end subroutine PMSurfaceSetRealization

! ************************************************************************** !

recursive subroutine PMSurfaceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  !

  implicit none

  class(pm_surface_type) :: this

end subroutine PMSurfaceInitializeRun

! ************************************************************************** !
subroutine PMSurfacePreSolve(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  use Global_module

  implicit none
  
  class(pm_surface_type) :: this
  
  this%option%io_buffer = 'PMSurfacePreSolve() must be extended.'
  call printErrMsg(this%option)  

end subroutine PMSurfacePreSolve

! ************************************************************************** !

subroutine PMSurfacePostSolve(this)
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Global_module

  implicit none
  
  class(pm_surface_type) :: this
  
  this%option%io_buffer = 'PMSurfacePostSolve() must be extended.'
  call printErrMsg(this%option)  
  
end subroutine PMSurfacePostSolve

! ************************************************************************** !

subroutine PMSurfaceUpdateSolution(this)
  !
  ! As a first step in updating the solution, update all flow-conditions.
  ! The solution will be updated by each child class of pm_surface_type.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Condition_module

  implicit none

  class(pm_surface_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE


  ! begin from RealizationUpdate()
  call FlowConditionUpdate(this%surf_realization%surf_flow_conditions, &
                           this%surf_realization%option, &
                           this%surf_realization%option%time)

  call RealizSurfAllCouplerAuxVars(this%surf_realization,force_update_flag)

end subroutine PMSurfaceUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceUpdateAuxVars(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14

  implicit none
  
  class(pm_surface_type) :: this

  this%option%io_buffer = 'PMSurfaceUpdateAuxVars() must be extended.'
  call printErrMsg(this%option)

end subroutine PMSurfaceUpdateAuxVars

! ************************************************************************** !

subroutine PMSurfaceCheckpointBinary(this,viewer)
  ! 
  ! This routine checkpoints data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Checkpoint_Surface_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_surface_type) :: this
  PetscViewer :: viewer

  call SurfaceCheckpointProcessModelBinary(viewer,this%surf_realization)

end subroutine PMSurfaceCheckpointBinary

! ************************************************************************** !

subroutine PMSurfaceRestartBinary(this,viewer)
  ! 
  ! This routine reads checkpoint data associated with surface-flow PM
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  use Checkpoint_Surface_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_surface_type) :: this
  PetscViewer :: viewer

  call SurfaceRestartProcessModelBinary(viewer,this%surf_realization)
  call this%UpdateAuxVars()
  call this%UpdateSolution()

end subroutine PMSurfaceRestartBinary

! ************************************************************************** !

recursive subroutine PMSurfaceFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  implicit none

  class(pm_surface_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMSurfaceFinalizeRun

! ************************************************************************** !

subroutine PMSurfaceInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_surface_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMSurfaceInputRecord

! ************************************************************************** !

subroutine PMSurfaceDestroy(this)
  ! 
  ! Destroys Surface process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/22/14
  ! 

  implicit none

  class(pm_surface_type) :: this

  call this%comm1%Destroy()

end subroutine PMSurfaceDestroy

end module PM_Surface_class
