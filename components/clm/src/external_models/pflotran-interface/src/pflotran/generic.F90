module Generic_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  interface GenericParameterCreate
    module procedure GenericParameterCreate1
    module procedure GenericParameterCreate2
  end interface GenericParameterCreate

  type, public :: generic_parameter_type
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: rvalue
    PetscInt :: ivalue
    type(generic_parameter_type), pointer :: next
  end type generic_parameter_type

  public :: GenericParameterCreate, &
            GenericParameterAddToList, &
            GenericParameterGetNamedValue, &
            GenericParameterDestroy

contains

! ************************************************************************** !

function GenericParameterCreate1()
  ! 
  ! Allocates and initializes a new generic parameter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  implicit none
  
  type(generic_parameter_type), pointer :: GenericParameterCreate1
  
  GenericParameterCreate1 => &
    GenericParameterCreate2('',UNINITIALIZED_INTEGER,UNINITIALIZED_DOUBLE)
  
end function GenericParameterCreate1

! ************************************************************************** !

function GenericParameterCreate2(string,ivalue,rvalue)
  ! 
  ! Allocates and initializes a new generic parameter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  implicit none
  
  type(generic_parameter_type), pointer :: GenericParameterCreate2
  
  type(generic_parameter_type), pointer :: generic_parameter

  character(len=*) :: string
  PetscInt :: ivalue
  PetscReal :: rvalue
  
  allocate(generic_parameter)
  generic_parameter%name = trim(adjustl(string))
  generic_parameter%ivalue = ivalue
  generic_parameter%rvalue = rvalue
  nullify(generic_parameter%next)

  GenericParameterCreate2 => generic_parameter
  
end function GenericParameterCreate2

! ************************************************************************** !

subroutine GenericParameterAddToList(list,new_generic_parameter)
  ! 
  ! Deallocates generic parameter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  implicit none
  
  type(generic_parameter_type), pointer :: list
  type(generic_parameter_type), pointer :: new_generic_parameter
  
  type(generic_parameter_type), pointer :: cur_generic_parameter

  if (associated(list)) then
    cur_generic_parameter => list
    do
      if (.not.associated(cur_generic_parameter%next)) exit
      cur_generic_parameter => cur_generic_parameter%next
    enddo
    cur_generic_parameter%next => new_generic_parameter
  else
    list => new_generic_parameter
  endif
  
end subroutine GenericParameterAddToList

! ************************************************************************** !

subroutine GenericParameterGetNamedValue(list,string,ivalue,rvalue)
  ! 
  ! Returns the integer and/or double precision value associated with a string
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  use String_module

  implicit none
  
  type(generic_parameter_type), pointer :: list
  character(len=*) :: string
  PetscInt :: ivalue
  PetscReal :: rvalue
  
  type(generic_parameter_type), pointer :: cur_generic_parameter

  ivalue = UNINITIALIZED_INTEGER
  rvalue = UNINITIALIZED_DOUBLE

  cur_generic_parameter => list
  do
    if (.not.associated(cur_generic_parameter)) exit
    if (StringCompare(string,cur_generic_parameter%name)) then
      ivalue = cur_generic_parameter%ivalue
      rvalue = cur_generic_parameter%rvalue
    endif
    cur_generic_parameter => cur_generic_parameter%next
  enddo
  
end subroutine GenericParameterGetNamedValue

! ************************************************************************** !

recursive subroutine GenericParameterDestroy(generic_parameter)
  ! 
  ! Deallocates generic parameter object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/18
  ! 
  implicit none
  
  type(generic_parameter_type), pointer :: generic_parameter
  
  if (.not.associated(generic_parameter)) return

  deallocate(generic_parameter)
  nullify(generic_parameter)
  
end subroutine GenericParameterDestroy

end module Generic_module
