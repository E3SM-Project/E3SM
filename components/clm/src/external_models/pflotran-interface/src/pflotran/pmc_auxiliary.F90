module PMC_Auxiliary_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PMC_Base_class
  use PM_Auxiliary_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

  implicit none

  
  private
  type, public, extends(pmc_base_type) :: pmc_auxiliary_type
    class(pm_auxiliary_type), pointer :: pm_aux
  contains
    procedure, public :: Init => PMCAuxiliaryInit
    procedure, public :: RunToTime => PMCAuxiliaryRunToTime
    procedure, public :: Destroy => PMCAuxiliaryDestroy
  end type pmc_auxiliary_type
  
  public :: PMCAuxiliaryCreate
  
contains

! ************************************************************************** !

function PMCAuxiliaryCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_auxiliary_type), pointer :: PMCAuxiliaryCreate
  
  class(pmc_auxiliary_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCAuxiliary%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()
  
  PMCAuxiliaryCreate => pmc  
  
end function PMCAuxiliaryCreate

! ************************************************************************** !

subroutine PMCAuxiliaryInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 

  implicit none
  
  class(pmc_auxiliary_type) :: this
  
#ifdef DEBUG
  print *, 'PMCAuxiliary%Init()'
#endif
  
  call PMCBaseInit(this)
  this%name = 'PMCAuxiliary'

end subroutine PMCAuxiliaryInit

! ************************************************************************** !

recursive subroutine PMCAuxiliaryRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
!  use Init_Subsurface_module
  use Option_module
  
  implicit none
  
#include "petsc/finclude/petscviewer.h"  

  class(pmc_auxiliary_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  PetscInt :: local_stop_flag
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return
  
  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE
  ! if at end of simulation, skip update of material properties
  if (stop_flag /= TS_STOP_END_SIMULATION) then
    ! must use ierr here due to 32-/64-bit integer issues
    call this%pm_aux%Evaluate(sync_time,ierr)
    local_stop_flag = ierr
  endif
  
  ! Run underlying process model couplers
  if (associated(this%child)) then
    ! Set data needed by process-models
    call this%SetAuxData()
    call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
    ! Get data from other process-models
    call this%GetAuxData()
  endif

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer)) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif
  
  stop_flag = max(stop_flag,local_stop_flag)
  
  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMCAuxiliaryRunToTime

! ************************************************************************** !
!
! PMCAuxiliaryFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCAuxiliaryFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_auxiliary_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCAuxiliary%FinalizeRun()')
#endif
  
end subroutine PMCAuxiliaryFinalizeRun

! ************************************************************************** !

subroutine PMCAuxiliaryStrip(this)
  !
  ! Deallocates members of PMC Auxiliary.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_auxiliary_type) :: this

  call PMCBaseStrip(this)
  nullify(this%pm_aux)

end subroutine PMCAuxiliaryStrip

! ************************************************************************** !

recursive subroutine PMCAuxiliaryDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module

  implicit none
  
  class(pmc_auxiliary_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCAuxiliary%Destroy()')
#endif

  call PMCAuxiliaryStrip(this)
  
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
  
end subroutine PMCAuxiliaryDestroy
  
end module PMC_Auxiliary_class
