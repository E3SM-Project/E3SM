module Data_Mediator_Vec_class
 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Data_Mediator_Base_class

  implicit none

  private
 
  type, public, extends(data_mediator_base_type) :: data_mediator_vec_type
    VecScatter :: scatter_ctx ! scatter context from vec to residual_vec
    Vec :: vec
  contains
    procedure, public :: Update => DataMediatorVecUpdate
    procedure, public :: Strip => DataMediatorVecStrip
  end type data_mediator_vec_type
  
  public :: DataMediatorVecCreate

contains

! ************************************************************************** !

function DataMediatorVecCreate()
  ! 
  ! Creates a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/24/15
  ! 
  
  implicit none

  class(data_mediator_vec_type), pointer :: DataMediatorVecCreate
  
  class(data_mediator_vec_type), pointer :: data_mediator
  
  allocate(data_mediator)
  call DataMediatorBaseCreate(data_mediator)
  data_mediator%vec = PETSC_NULL_VEC
  data_mediator%scatter_ctx = PETSC_NULL_VECSCATTER
  DataMediatorVecCreate => data_mediator

end function DataMediatorVecCreate

! ************************************************************************** !

recursive subroutine DataMediatorVecUpdate(this,data_mediator_vec,option)
  ! 
  ! Updates a data mediator object transfering data from
  ! the buffer into the PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/24/15
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none

  class(data_mediator_vec_type) :: this
  Vec :: data_mediator_vec
  type(option_type) :: option  
  
  PetscErrorCode :: ierr
  
  call VecScatterBegin(this%scatter_ctx,this%vec, &
                       data_mediator_vec,ADD_VALUES, &
                       SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(this%scatter_ctx,this%vec, &
                     data_mediator_vec,ADD_VALUES, &
                     SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                         
end subroutine DataMediatorVecUpdate

! ************************************************************************** !

recursive subroutine DataMediatorVecStrip(this)
  ! 
  ! Destroys a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/24/15
  ! 

  implicit none
  
  class(data_mediator_vec_type) :: this
  
  PetscErrorCode :: ierr
  
  ! update the next one
  if (associated(this%next)) then
    call this%next%Strip()
    deallocate(this%next)
    nullify(this%next)
  endif 
  
  ! Simply nullify the pointer as the dataset resides in a list to be
  ! destroyed separately.
  !nullify(data_mediator%dataset)
  if (this%scatter_ctx /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_ctx,ierr);CHKERRQ(ierr)
  endif
  if (this%vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine DataMediatorVecStrip

end module Data_Mediator_Vec_class
