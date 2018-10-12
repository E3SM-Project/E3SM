module Communicator_Unstructured_class
#include "petsc/finclude/petscdm.h"
  use petscdm
  use Communicator_Base_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Explicit_module  
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(communicator_type) :: unstructured_communicator_type
    DM :: dm
    type(ugdm_type), pointer :: ugdm
  contains
    procedure, public :: SetDM => UnstructuredSetDM
    procedure, public :: GlobalToLocal => UnstructuredGlobalToLocal
    procedure, public :: LocalToGlobal => UnstructuredLocalToGlobal
    procedure, public :: LocalToLocal => UnstructuredLocalToLocal
    procedure, public :: GlobalToNatural => UnstructuredGlobalToNatural
    procedure, public :: NaturalToGlobal => UnstructuredNaturalToGlobal
    procedure, public :: AONaturalToPetsc => UnstructuredAONaturalToPetsc
!geh: finalization not yet supported by gfortran.
!    final :: UnstructuredCommunicatorDestroy
    procedure, public :: Destroy => UnstructuredCommunicatorDestroy
  end type unstructured_communicator_type

  public :: UnstructuredCommunicatorCreate
  
contains

! ************************************************************************** !

function UnstructuredCommunicatorCreate()
  ! 
  ! Allocates and initializes a new communicator
  ! object for unstructured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type), pointer :: &
    UnstructuredCommunicatorCreate
  
  class(unstructured_communicator_type), pointer :: communicator
  
  allocate(communicator)
  nullify(communicator%ugdm)
  communicator%dm = PETSC_NULL_DM

  UnstructuredCommunicatorCreate => communicator  
  
end function UnstructuredCommunicatorCreate

! ************************************************************************** !

subroutine UnstructuredSetDM(this,dm_ptr)
  ! 
  ! Sets pointer to DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use DM_Kludge_module

  implicit none
  
  class(unstructured_communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm = dm_ptr%dm
  this%ugdm => dm_ptr%ugdm
  
end subroutine UnstructuredSetDM

! ************************************************************************** !

subroutine UnstructuredGlobalToLocal(this,source,destination)
  ! 
  ! Performs global to local communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

!TODO(geh): move to communicator_base.F90

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMGlobalToLocalBegin(this%dm,source,INSERT_VALUES,destination, &
                            ierr);CHKERRQ(ierr)
  call DMGlobalToLocalEnd(this%dm,source,INSERT_VALUES,destination, &
                          ierr);CHKERRQ(ierr)
  
end subroutine UnstructuredGlobalToLocal

! ************************************************************************** !

subroutine UnstructuredLocalToGlobal(this,source,destination)
  ! 
  ! Performs local to global communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_ltog,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%ugdm%scatter_ltog,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      
end subroutine UnstructuredLocalToGlobal

! ************************************************************************** !

subroutine UnstructuredLocalToLocal(this,source,destination)
  ! 
  ! Performs local to local communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call VecScatterBegin(this%ugdm%scatter_ltol,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%ugdm%scatter_ltol,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine UnstructuredLocalToLocal

! ************************************************************************** !

subroutine UnstructuredGlobalToNatural(this,source,destination)
  ! 
  ! Performs global to natural communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_gton,source,destination, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%ugdm%scatter_gton,source,destination, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  
end subroutine UnstructuredGlobalToNatural

! ************************************************************************** !

subroutine UnstructuredNaturalToGlobal(this,source,destination)
  ! 
  ! Performs natural to global communication
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr

  call VecScatterBegin(this%ugdm%scatter_gton,source,destination, &
                       INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%ugdm%scatter_gton,source,destination, &
                     INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  
end subroutine UnstructuredNaturalToGlobal

! ************************************************************************** !

subroutine UnstructuredAONaturalToPetsc(this,array)
  ! 
  ! Maps indices in natural numbering to petsc
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/19/15
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  PetscInt :: array(:)

  AO :: ao
  PetscInt :: n
  PetscErrorCode :: ierr
  
  n = size(array)
  call AOApplicationToPetsc(this%ugdm%ao_natural_to_petsc,n,array, &
                            ierr);CHKERRQ(ierr)
  
end subroutine UnstructuredAONaturalToPetsc

! ************************************************************************** !

subroutine UnstructuredCommunicatorDestroy(this)
  ! 
  ! Deallocates a communicator object for
  ! unstructured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(unstructured_communicator_type) :: this
  
  PetscErrorCode :: ierr
  
  !geh: all DMs are currently destroyed in realization.  This DM is solely
  !     a pointer.  This will need to change, but skip for now.
  if (associated(this%ugdm)) then
    !call UGridDMDestroy(this%ugdm)
  endif
  nullify(this%ugdm)
  if (this%dm /= PETSC_NULL_DM) then
    !call DMDestroy(this%dm,ierr)
  endif
  this%dm = PETSC_NULL_DM 
  
end subroutine UnstructuredCommunicatorDestroy

end module Communicator_Unstructured_class
