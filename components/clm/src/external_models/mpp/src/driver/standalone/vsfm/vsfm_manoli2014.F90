program vsfm_manoli2014
  !
#include <petsc/finclude/petsc.h>
  !
  use vsfm_manoli2014_problem      , only : run_vsfm_manoli2014_problem
  use petscsys
  !
  implicit none
  !
  PetscErrorCode     :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD
  PETSC_COMM_SELF  = MPI_COMM_SELF

  write(*,*)'call run_vsfm_manoli2014_problem()'
  call run_vsfm_manoli2014_problem()

  call PetscFinalize(ierr)

end program vsfm_manoli2014
