module AuxVars_Base_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  type, public :: auxvar_base_type
    PetscReal :: effective_porosity ! factors in compressibility - common to all modes??
    PetscReal :: pert ! common to all modes to (perturbation for numerical jacobian)
  contains
    procedure, public :: Init => AuxVarBaseInit
  end type auxvar_base_type     

contains

! ************************************************************************** !

subroutine AuxVarBaseInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_base_type) :: this
  type(option_type) :: option

  !currently does nothing - could init the base members
  print *, 'Must extend InitBaseAuxVars '
  stop    

end subroutine AuxVarBaseInit


! ************************************************************************** !

end module AuxVars_Base_module

