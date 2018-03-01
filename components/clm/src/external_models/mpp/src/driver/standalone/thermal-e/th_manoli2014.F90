program th_manoli2014
  !
#include <petsc/finclude/petsc.h>
  !
  use th_manoli2014_problem , only : run_th_manoli2014_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  call run_th_manoli2014_problem()

  call PetscFinalize(ierr)

end program th_manoli2014
