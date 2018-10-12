module Data_Mediator_Base_class
#include "petsc/finclude/petscvec.h"
  use petscvec

  use Dataset_Global_HDF5_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private
 
  type, public :: data_mediator_base_type
    character(len=MAXWORDLENGTH) :: name
    class(data_mediator_base_type), pointer :: next
  contains
    procedure, public :: AddToList => DataMediatorBaseAddToList
    procedure, public :: Update => DataMediatorBaseUpdate
    procedure, public :: Strip => DataMediatorBaseStrip
  end type data_mediator_base_type
  
  public :: DataMediatorBaseCreate

contains

! ************************************************************************** !

subroutine DataMediatorBaseCreate(this)
  ! 
  ! Creates a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/23/15
  ! 
  
  implicit none

  class(data_mediator_base_type) :: this
  
  ! Cannot allocate here.  Allocation takes place in daughter class  
  this%name = ''
  nullify(this%next)

end subroutine DataMediatorBaseCreate

! ************************************************************************** !

recursive subroutine DataMediatorBaseAddToList(this,list)
  ! 
  ! Adds a data mediator object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/23/15
  ! 

  implicit none
  
  class(data_mediator_base_type), target :: this
  class(data_mediator_base_type), pointer :: list

  if (associated(list)) then
    call this%AddToList(list%next)
  else
    list => this
  endif
  
end subroutine DataMediatorBaseAddToList

! ************************************************************************** !

recursive subroutine DataMediatorBaseUpdate(this,data_mediator_vec,option)
  ! 
  ! Adds contribution of data mediator object to vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/23/15
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none

  class(data_mediator_base_type) :: this
  Vec :: data_mediator_vec
  type(option_type) :: option
  print *, 'Must extend DataMediatorBaseUpdate.'
  stop
  
end subroutine DataMediatorBaseUpdate

! ************************************************************************** !

recursive subroutine DataMediatorBaseStrip(this)
  ! 
  ! Destroys a data mediator object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/24/15
  ! 

  implicit none
  
  class(data_mediator_base_type) :: this
  
  PetscErrorCode :: ierr
  
end subroutine DataMediatorBaseStrip

end module Data_Mediator_Base_class
