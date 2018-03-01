program vsfm_sy1991

#include <petsc/finclude/petsc.h>

  use vsfm_sy1991_problem, only : run_vsfm_sy1991_problem
  use petscsys

  implicit none

  PetscErrorCode :: ierr

  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_vsfm_sy1991_problem()

  call PetscFinalize(ierr)

end program vsfm_sy1991
