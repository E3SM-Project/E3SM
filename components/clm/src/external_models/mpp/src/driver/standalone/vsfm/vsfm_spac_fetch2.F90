program vsfm_spac_fetch2

#include <petsc/finclude/petsc.h>

  use vsfm_spac_fetch2_problem, only : run_vsfm_spac_fetch2_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  call run_vsfm_spac_fetch2_problem()

  call PetscFinalize(ierr)

end program vsfm_spac_fetch2
