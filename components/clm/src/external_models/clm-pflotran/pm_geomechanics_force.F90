module PM_Geomechanics_Force_class

  use PM_Base_class
  use Geomechanics_Realization_class
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

  type, public, extends(pm_base_type) :: pm_geomech_force_type
    class(realization_geomech_type), pointer :: geomech_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: Setup => PMGeomechForceSetup
    procedure, public :: PMGeomechForceSetRealization
    procedure, public :: InitializeRun => PMGeomechForceInitializeRun
    procedure, public :: FinalizeRun => PMGeomechForceFinalizeRun
    procedure, public :: InitializeTimestep => PMGeomechForceInitializeTimestep
    procedure, public :: Residual => PMGeomechForceResidual
    procedure, public :: Jacobian => PMGeomechForceJacobian
    procedure, public :: PreSolve => PMGeomechForcePreSolve
    procedure, public :: UpdateSolution => PMGeomechForceUpdateSolution
    procedure, public :: CheckpointBinary => PMGeomechForceCheckpointBinary
    procedure, public :: RestartBinary => PMGeomechForceRestartBinary
    procedure, public :: InputRecord => PMGeomechForceInputRecord
    procedure, public :: Destroy => PMGeomechForceDestroy
    procedure, public :: FinalizeTimestep => PMGeomechForceFinalizeTimestep
  end type pm_geomech_force_type

  public :: PMGeomechForceCreate

contains

! ************************************************************************** !

function PMGeomechForceCreate()
  ! 
  ! This routine creates
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  implicit none

  class(pm_geomech_force_type), pointer :: PMGeomechForceCreate

  class(pm_geomech_force_type), pointer :: geomech_force_pm

  allocate(geomech_force_pm)
  nullify(geomech_force_pm%option)
  nullify(geomech_force_pm%output_option)
  nullify(geomech_force_pm%geomech_realization)
  nullify(geomech_force_pm%comm1)

  call PMBaseInit(geomech_force_pm)

  PMGeomechForceCreate => geomech_force_pm

end function PMGeomechForceCreate

! ************************************************************************** !

subroutine PMGeomechForceSetup(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this

  ! set up communicator
  select case(this%geomech_realization%geomech_discretization%itype)
    case(STRUCTURED_GRID)
      this%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  !call this%comm1%SetDM(this%geomech_realization%geomech_discretization%dm_1dof)

end subroutine PMGeomechForceSetup

! ************************************************************************** !

recursive subroutine PMGeomechForceInitializeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Force_module, only : GeomechUpdateSolution

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForceInitializeRun

! ************************************************************************** !

recursive subroutine PMGeomechForceFinalizeRun(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  implicit none

  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%FinalizeRun()')
#endif

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMGeomechForceFinalizeRun

! ************************************************************************** !

subroutine PMGeomechForceSetRealization(this, geomech_realization)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this
  class(realization_geomech_type), pointer :: geomech_realization

  this%geomech_realization => geomech_realization
  this%realization_base => geomech_realization

  this%solution_vec = geomech_realization%geomech_field%disp_xx
  this%residual_vec = geomech_realization%geomech_field%disp_r

end subroutine PMGeomechForceSetRealization

! ************************************************************************** !

subroutine PMGeomechForceInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Force_module, only : GeomechanicsForceInitialGuess
  use Global_module
  
  implicit none
  
  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%InitializeTimestep()')
#endif

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GEOMECHANICS ",62("="))')
  endif
  
  call GeomechanicsForceInitialGuess(this%geomech_realization)
  
end subroutine PMGeomechForceInitializeTimestep

! ************************************************************************** !

subroutine PMGeomechForceResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Force_module, only : GeomechForceResidual

  implicit none
  
  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%Residual()')
#endif
  
  call GeomechForceResidual(snes,xx,r,this%geomech_realization,ierr)

end subroutine PMGeomechForceResidual

! ************************************************************************** !

subroutine PMGeomechForceJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Force_module, only : GeomechForceJacobian

  implicit none
  
  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%Jacobian()')
#endif
  
  call GeomechForceJacobian(snes,xx,A,B,this%geomech_realization,ierr)

end subroutine PMGeomechForceJacobian

! ************************************************************************** !

subroutine PMGeomechForcePreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForcePreSolve

! ************************************************************************** !

subroutine PMGeomechForceUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Force_module, only : GeomechUpdateSolution, &
                                        GeomechStoreInitialDisp, &
                                        GeomechForceUpdateAuxVars
  use Condition_module

  implicit none

  class(pm_geomech_force_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call GeomechUpdateSolution(this%geomech_realization)
  if (this%option%geomech_initial) then
    call GeomechStoreInitialDisp(this%geomech_realization)
    this%option%geomech_initial = PETSC_FALSE
  endif
  call GeomechForceUpdateAuxVars(this%geomech_realization)

end subroutine PMGeomechForceUpdateSolution

! ************************************************************************** !

subroutine PMGeomechForceFinalizeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Global_module

  implicit none
  
  class(pm_geomech_force_type) :: this
  
#ifdef PM_GEOMECH_FORCE_DEBUG  
  call printMsg(this%option,'PMGeomechForce%FinalizeTimestep()')
#endif

end subroutine PMGeomechForceFinalizeTimestep

! ************************************************************************** !

subroutine PMGeomechForceCheckpointBinary(this,viewer)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer
  
  call printErrMsg(this%option,'add code for checkpointing Geomech in PM approach')
  
end subroutine PMGeomechForceCheckpointBinary

! ************************************************************************** !

subroutine PMGeomechForceRestartBinary(this,viewer)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer
  
  call printErrMsg(this%option,'add code for restarting Geomech in PM approach')
  
end subroutine PMGeomechForceRestartBinary

! ************************************************************************** !

subroutine PMGeomechForceInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_geomech_force_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMGeomechForceInputRecord

! ************************************************************************** !

subroutine PMGeomechForceDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  ! 

  use Geomechanics_Realization_class, only : GeomechRealizDestroy

  implicit none
  
  class(pm_geomech_force_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_GEOMECH_FORCE_DEBUG
  call printMsg(this%option,'PMGeomechForce%Destroy()')
#endif

  call GeomechRealizDestroy(this%geomech_realization)

  call this%comm1%Destroy()
  
end subroutine PMGeomechForceDestroy

end module PM_Geomechanics_Force_class
