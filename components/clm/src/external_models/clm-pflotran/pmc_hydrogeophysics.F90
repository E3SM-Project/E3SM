module PMC_Hydrogeophysics_class

  use PMC_Base_class
  use Realization_Subsurface_class
  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
  
  type, public, extends(pmc_base_type) :: pmc_hydrogeophysics_type
    class(realization_subsurface_type), pointer :: realization
    Vec :: tracer_seq
    Vec :: saturation_seq
    Vec :: temperature_seq
    ! a pointer to xxx_mpi in simulation_hydrogeophysics_type
    Vec :: tracer_mpi 
    Vec :: saturation_mpi 
    Vec :: temperature_mpi 
    VecScatter :: pf_to_e4d_scatter
    PetscMPIInt :: pf_to_e4d_master_comm
  contains
    procedure, public :: Init => PMCHydrogeophysicsInit
    procedure, public :: InitializeRun => PMCHydrogeophysicsInitializeRun
    procedure, public :: RunToTime => PMCHydrogeophysicsRunToTime
    procedure, public :: FinalizeRun => PMCHydrogeophysicsFinalizeRun
    procedure, public :: Destroy => PMCHydrogeophysicsDestroy
    procedure, public :: GetAuxData => PMCHydrogeophysicsSynchronize
  end type pmc_hydrogeophysics_type
  
  public :: PMCHydrogeophysicsCreate
  
contains

! ************************************************************************** !

function PMCHydrogeophysicsCreate()
  ! 
  ! Allocates and initializes a new
  ! process_model_coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_hydrogeophysics_type), pointer :: PMCHydrogeophysicsCreate
  
  class(pmc_hydrogeophysics_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCHydrogeophysicsCreate => pmc  
  
end function PMCHydrogeophysicsCreate

! ************************************************************************** !

subroutine PMCHydrogeophysicsInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call PMCBaseInit(this)
  this%name = 'PMCHydrogeophysics'
  nullify(this%realization) 
  this%tracer_mpi = 0
  this%tracer_seq = 0
  this%saturation_mpi = 0
  this%saturation_seq = 0
  this%temperature_mpi = 0
  this%temperature_seq = 0
  this%pf_to_e4d_scatter = 0
  this%pf_to_e4d_master_comm = 0
  
end subroutine PMCHydrogeophysicsInit

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 
  use Realization_Base_class, only : RealizationGetVariable
  use Variables_module, only : POROSITY

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this

  PetscReal, pointer :: vec1_ptr(:)
  PetscReal, pointer :: vec2_ptr(:)
  PetscErrorCode :: ierr
  
  call printMsg(this%option,'PMCHydrogeophysics%InitializeRun()')
  
  if (associated(this%child)) then
    call this%child%InitializeRun()
  endif
  
  if (associated(this%peer)) then
    call this%peer%InitializeRun()
  endif

  ! send porosity once to E4D through tracer Vecs
  ! send as late as possible.
  call RealizationGetVariable(this%realization,this%realization%field%work, &
                              POROSITY,ZERO_INTEGER)
  call VecGetArrayF90(this%realization%field%work,vec1_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%tracer_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  vec2_ptr(:) = vec1_ptr(:)
  call VecRestoreArrayF90(this%realization%field%work,vec1_ptr, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%tracer_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  call VecScatterBegin(this%pf_to_e4d_scatter,this%tracer_mpi,this%tracer_seq, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%pf_to_e4d_scatter,this%tracer_mpi,this%tracer_seq, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

end subroutine PMCHydrogeophysicsInitializeRun

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Wrapper_Hydrogeophysics_module, only : HydrogeophysicsWrapperStep
  use Timestepper_Base_class, only : TS_CONTINUE

  implicit none
  
  class(pmc_hydrogeophysics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  class(pmc_base_type), pointer :: pmc_base
  PetscInt :: local_stop_flag
  
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE
  
  call HydrogeophysicsWrapperStep(sync_time, &
                                  this%tracer_mpi, &
                                  this%tracer_seq, &
                                  this%saturation_mpi, &
                                  this%saturation_seq, &
                                  this%temperature_mpi, &
                                  this%temperature_seq, &
                                  this%pf_to_e4d_scatter, &
                                  this%pf_to_e4d_master_comm,this%option)

  ! Run neighboring process model couplers
!  if (associated(this%child)) then
!    call this%child%RunToTime(sync_time,local_stop_flag)
!  endif

  ! Run neighboring process model couplers
!  if (associated(this%peer)) then
!    call this%peer%RunToTime(sync_time,local_stop_flag)
!  endif

  stop_flag = max(stop_flag,local_stop_flag)  
  
end subroutine PMCHydrogeophysicsRunToTime

! ************************************************************************** !

subroutine PMCHydrogeophysicsSynchronize(this)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Variables_module, only : PRIMARY_MOLALITY, LIQUID_SATURATION, TEMPERATURE
  use String_module
!  use Discretization_module

  implicit none
  
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscviewer.h"

  class(pmc_hydrogeophysics_type) :: this

  PetscReal, pointer :: vec1_ptr(:), vec2_ptr(:)
  PetscErrorCode :: ierr
  PetscInt, save :: num_calls = 0
  character(len=MAXSTRINGLENGTH) :: filename
  PetscViewer :: viewer
  Vec :: natural_vec
  PetscInt :: i

#if 1
  ! tracer
  call RealizationGetVariable(this%realization,this%realization%field%work, &
                              PRIMARY_MOLALITY,ONE_INTEGER,ZERO_INTEGER)
#else
  call DiscretizationCreateVector(this%realization%discretization,ONEDOF, &
                                  natural_vec,NATURAL,this%option)
  if (this%option%myrank == 0) then
    do i = 1, this%realization%patch%grid%nmax
      call VecSetValues(natural_vec,1,i-1,i*1.d0, &
                        INSERT_VALUES,ierr);CHKERRQ(ierr)
    enddo
  endif
  call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
  call DiscretizationNaturalToGlobal(this%realization%discretization, &
                                      natural_vec, &
                                      this%realization%field%work,ONEDOF)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
#endif
  call VecGetArrayF90(this%realization%field%work,vec1_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%tracer_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
!      vec1_ptr(:) = vec1_ptr(:) + num_calls
  vec2_ptr(:) = vec1_ptr(:)
!  print *, 'PMC update to tracer', vec2_ptr(16)
  call VecRestoreArrayF90(this%realization%field%work,vec1_ptr, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%tracer_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  
  ! liquid saturation
  call RealizationGetVariable(this%realization,this%realization%field%work, &
                              LIQUID_SATURATION,ZERO_INTEGER)
  call VecGetArrayF90(this%realization%field%work,vec1_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%saturation_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  vec2_ptr(:) = vec1_ptr(:)
  call VecRestoreArrayF90(this%realization%field%work,vec1_ptr, &
                          ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%saturation_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  
  ! temperature
  if (this%temperature_mpi /= 0) then
    call RealizationGetVariable(this%realization,this%realization%field%work, &
                                TEMPERATURE,ZERO_INTEGER)
    call VecGetArrayF90(this%realization%field%work,vec1_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%temperature_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
    vec2_ptr(:) = vec1_ptr(:)
    call VecRestoreArrayF90(this%realization%field%work,vec1_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(this%temperature_mpi,vec2_ptr,ierr);CHKERRQ(ierr)
  endif
  
#if 0
  filename = 'pf_tracer' // trim(StringFormatInt(num_calls)) // '.txt'
  call PetscViewerASCIIOpen(this%option%mycomm,filename,viewer, &
                            ierr);CHKERRQ(ierr)
  call VecView(this%realization%field%work,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  num_calls = num_calls + 1
  
end subroutine PMCHydrogeophysicsSynchronize

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Wrapper_Hydrogeophysics_module, only : HydrogeophysicsWrapperStop
  
  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  call printMsg(this%option,'PMCHydrogeophysics%FinalizeRun()')
  
  if (this%pf_to_e4d_master_comm /= MPI_COMM_NULL) then
    call HydrogeophysicsWrapperStop(this%option,this%pf_to_e4d_master_comm)
  endif
  
end subroutine PMCHydrogeophysicsFinalizeRun

! ************************************************************************** !

subroutine PMCHydrogeophysicsStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_hydrogeophysics_type) :: this
  
  PetscErrorCode :: ierr

  call PMCBaseStrip(this)
  nullify(this%realization)
  ! created in HydrogeophysicsInitialize()
  if (this%tracer_seq /= 0) then
    call VecDestroy(this%tracer_seq,ierr);CHKERRQ(ierr)
  endif
  this%tracer_seq = 0
  if (this%saturation_seq /= 0) then
    call VecDestroy(this%saturation_seq,ierr);CHKERRQ(ierr)
  endif
  this%saturation_seq = 0
  if (this%temperature_seq /= 0) then
    call VecDestroy(this%temperature_seq,ierr);CHKERRQ(ierr)
  endif
  this%temperature_seq = 0
  ! created in HydrogeophysicsInitialize()
  if (this%pf_to_e4d_scatter /= 0) then
    call VecScatterDestroy(this%pf_to_e4d_scatter, ierr);CHKERRQ(ierr)
  endif
  this%pf_to_e4d_scatter = 0
  ! these are solely pointers set in HydrogeophysicsInitialize()
  this%pf_to_e4d_master_comm = 0
  this%tracer_mpi = 0
  this%saturation_mpi = 0
  this%temperature_mpi = 0
  
end subroutine PMCHydrogeophysicsStrip

! ************************************************************************** !

recursive subroutine PMCHydrogeophysicsDestroy(this)
  ! 
  ! Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_hydrogeophysics_type) :: this

  PetscErrorCode :: ierr
  
  call printMsg(this%option,'PMCHydrogeophysics%Destroy()')
  
  call PMCHydrogeophysicsStrip(this)
  
  if (associated(this%child)) then
    call this%child%Destroy()
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
  endif  

end subroutine PMCHydrogeophysicsDestroy
  
end module PMC_Hydrogeophysics_class
