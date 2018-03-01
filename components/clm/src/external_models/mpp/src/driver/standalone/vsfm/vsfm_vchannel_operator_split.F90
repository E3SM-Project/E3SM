program vsfm_vchannel_operator_split
  !
#include <petsc/finclude/petsc.h>
  !
  use vsfm_vchannel_problem_operator_split, only : run_vsfm_vchannel_problem_operator_split
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_vsfm_vchannel_problem_operator_split()

  call PetscFinalize(ierr)

end program vsfm_vchannel_operator_split
