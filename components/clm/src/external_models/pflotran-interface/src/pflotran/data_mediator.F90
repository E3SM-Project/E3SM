module Data_Mediator_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys

  use Data_Mediator_Base_class
  
  implicit none

  private

  public :: DataMediatorInit, &
            DataMediatorUpdate, &
            DataMediatorDestroy

contains

! ************************************************************************** !

subroutine DataMediatorInit(data_mediator_list, option)
  ! 
  ! Initializes data mediator object 
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/25/15
  ! 
  use Option_module

  implicit none
  
  class(data_mediator_base_type) :: data_mediator_list
  type(option_type) :: option
  
  
end subroutine DataMediatorInit

! ************************************************************************** !

subroutine DataMediatorUpdate(data_mediator_list,vec,option)
  ! 
  ! Updates a data mediator object transfering data from
  ! the buffer into the PETSc Vec
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/25/15
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none
  
  class(data_mediator_base_type), pointer :: data_mediator_list
  Vec :: vec
  type(option_type) :: option
  
  class(data_mediator_base_type), pointer :: cur_data_mediator
  PetscErrorCode :: ierr
  
  if (associated(data_mediator_list)) then
    call VecZeroEntries(vec,ierr);CHKERRQ(ierr)
    cur_data_mediator => data_mediator_list
    do
      if (.not.associated(cur_data_mediator)) exit
      call cur_data_mediator%Update(vec,option)
      cur_data_mediator => cur_data_mediator%next
    enddo
  endif
  
end subroutine DataMediatorUpdate

! ************************************************************************** !

subroutine DataMediatorDestroy(data_mediator_list)
  ! 
  ! Destroys a data mediator object list
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/25/15
  ! 
  implicit none
  
  class(data_mediator_base_type), pointer :: data_mediator_list
  
  if (.not.associated(data_mediator_list)) return
  
  call data_mediator_list%Strip()
  deallocate(data_mediator_list)
  nullify(data_mediator_list)  
  
end subroutine DataMediatorDestroy

end module Data_Mediator_module
