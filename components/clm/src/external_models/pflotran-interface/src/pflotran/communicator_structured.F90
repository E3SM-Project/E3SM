module Communicator_Structured_class
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Communicator_Base_module
  use Grid_Structured_module  
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(communicator_type) :: structured_communicator_type
    DM :: dm
  contains
    procedure, public :: SetDM => StructuredSetDM
    procedure, public :: GlobalToLocal => StructuredGlobalToLocal
    procedure, public :: LocalToGlobal => StructuredLocalToGlobal
    procedure, public :: LocalToLocal => StructuredLocalToLocal
    procedure, public :: GlobalToNatural => StructuredGlobalToNatural
    procedure, public :: NaturalToGlobal => StructuredNaturalToGlobal
    procedure, public :: AONaturalToPetsc => StructuredAONaturalToPetsc
!geh: finalization not yet supported by gfortran.
!    final :: StructuredCommunicatorDestroy
    procedure, public :: Destroy => StructuredCommunicatorDestroy

  end type structured_communicator_type
  
  public :: StructuredCommunicatorCreate
  
contains

! ************************************************************************** !

function StructuredCommunicatorCreate()
  ! 
  ! Allocates and initializes a new communicator
  ! object for structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type), pointer :: StructuredCommunicatorCreate
  
  class(structured_communicator_type), pointer :: communicator
  
  allocate(communicator)
  communicator%dm = PETSC_NULL_DM

  StructuredCommunicatorCreate => communicator  
  
end function StructuredCommunicatorCreate

! ************************************************************************** !

subroutine StructuredSetDM(this,dm_ptr)
  ! 
  ! Sets pointer to DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use DM_Kludge_module

  implicit none
  
  class(structured_communicator_type) :: this
  type(dm_ptr_type) :: dm_ptr

  this%dm = dm_ptr%dm
  
end subroutine StructuredSetDM

! ************************************************************************** !

subroutine StructuredGlobalToLocal(this,source,destination)
  ! 
  ! Performs global to local communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMGlobalToLocalBegin(this%dm,source,INSERT_VALUES,destination, &
                            ierr);CHKERRQ(ierr)
  call DMGlobalToLocalEnd(this%dm,source,INSERT_VALUES,destination, &
                          ierr);CHKERRQ(ierr)
  
end subroutine StructuredGlobalToLocal

! ************************************************************************** !

subroutine StructuredLocalToGlobal(this,source,destination)
  ! 
  ! Performs local to global communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMLocalToGlobalBegin(this%dm,source,INSERT_VALUES,destination, &
                            ierr);CHKERRQ(ierr)
  call DMLocalToGlobalEnd(this%dm,source,INSERT_VALUES,destination, &
                          ierr);CHKERRQ(ierr)
  
end subroutine StructuredLocalToGlobal

! ************************************************************************** !

subroutine StructuredLocalToLocal(this,source,destination)
  ! 
  ! Performs local to local communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMLocalToLocalBegin(this%dm,source,INSERT_VALUES,destination, &
                           ierr);CHKERRQ(ierr)
  call DMLocalToLocalEnd(this%dm,source,INSERT_VALUES,destination, &
                         ierr);CHKERRQ(ierr)
  
end subroutine StructuredLocalToLocal

! ************************************************************************** !

subroutine StructuredGlobalToNatural(this,source,destination)
  ! 
  ! Performs global to natural communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMDAGlobalToNaturalBegin(this%dm,source,INSERT_VALUES,destination, &
                                ierr);CHKERRQ(ierr)
  call DMDAGlobalToNaturalEnd(this%dm,source,INSERT_VALUES,destination, &
                              ierr);CHKERRQ(ierr)
  
end subroutine StructuredGlobalToNatural

! ************************************************************************** !

subroutine StructuredNaturalToGlobal(this,source,destination)
  ! 
  ! Performs natural to global communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  Vec :: source
  Vec :: destination

  PetscErrorCode :: ierr
  
  call DMDANaturalToGlobalBegin(this%dm,source,INSERT_VALUES,destination, &
                                ierr);CHKERRQ(ierr)
  call DMDANaturalToGlobalEnd(this%dm,source,INSERT_VALUES,destination, &
                              ierr);CHKERRQ(ierr)
  
end subroutine StructuredNaturalToGlobal

! ************************************************************************** !

subroutine StructuredAONaturalToPetsc(this,array)
  ! 
  ! Maps indices in natural numbering to petsc using a DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/19/15
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  PetscInt :: array(:)

  AO :: ao
  PetscInt :: n
  PetscErrorCode :: ierr
  
  call DMDAGetAO(this%dm,ao,ierr);CHKERRQ(ierr)
  n = size(array)
  call AOApplicationToPetsc(ao,n,array,ierr);CHKERRQ(ierr)
  
end subroutine StructuredAONaturalToPetsc

! ************************************************************************** !

subroutine StructuredCommunicatorDestroy(this)
  ! 
  ! Deallocates a communicator object for
  ! structured grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/15/13
  ! 

  implicit none
  
  class(structured_communicator_type) :: this
  
  PetscErrorCode :: ierr
  
  if (this%dm /= PETSC_NULL_DM) then
    !geh: all DMs are currently destroyed in realization.  This DM is solely
    !     a pointer.  This will need to change, but skip for now.
    !call DMDestroy(this%dm,ierr)
  endif
  this%dm = PETSC_NULL_DM
  
end subroutine StructuredCommunicatorDestroy

end module Communicator_Structured_class
