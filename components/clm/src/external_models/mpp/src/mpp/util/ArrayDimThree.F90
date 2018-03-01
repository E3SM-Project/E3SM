module ArrayDimThree

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: array_dim3_type
     PetscReal, dimension(3) :: arr
  end type array_dim3_type

#endif

end module ArrayDimThree
