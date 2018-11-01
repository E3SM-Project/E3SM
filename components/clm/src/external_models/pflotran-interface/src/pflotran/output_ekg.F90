module Output_EKG_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Output_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: IUNIT_EKG = 87

  public :: OutputEKGInit, &
            OutputEKGFinalize
            
contains

! ************************************************************************** !

subroutine OutputEKGInit(option,num_steps)
  ! 
  ! Open EKG file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/15
  ! 
  use Option_module

  implicit none
  
  type(option_type) :: option
  PetscInt :: num_steps

  character(len=MAXSTRINGLENGTH) :: filename
  PetscBool :: lexists
  
  if (.not.option%print_ekg) return

  filename = trim(option%global_prefix) // trim(option%group_prefix) // '.ekg'
  if (OptionPrintToFile(option)) then
    inquire(file=filename,exist=lexists)
    if (num_steps == 0 .or. .not. lexists) then
      open(unit=IUNIT_EKG,file=filename,action="write",status="replace")
    else
      open(unit=IUNIT_EKG,file=filename,action="write",status="old", &
           position="append")
    endif  
  endif

end subroutine OutputEKGInit

! ************************************************************************** !

subroutine OutputEKGFinalize()
  ! 
  ! Closes the EKG file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/15
  !
  PetscBool :: lopened
  
  inquire(IUNIT_EKG,opened=lopened)
  if (lopened) then
    close(IUNIT_EKG)
  endif
  
end subroutine OutputEKGFinalize

end module Output_EKG_module
