#ifdef USE_PETSC_LIB


module ArrayDimThree

  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"


  type, public :: array_dim3_type
     PetscReal, dimension(3) :: arr
  end type array_dim3_type

end module ArrayDimThree
#endif
