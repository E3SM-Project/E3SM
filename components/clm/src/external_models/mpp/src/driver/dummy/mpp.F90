program main

#include <petsc/finclude/petsc.h>

  use MultiPhysicsProbVSFM , only : vsfm_mpp
  use mpp_varpar           , only : mpp_varpar_init
  use petscsys
  !
  implicit none
  !
  !
  PetscInt       :: nlevsoi, nlevgrnd, nlevsno
  PetscErrorCode :: ierr
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  PETSC_COMM_WORLD = MPI_COMM_WORLD

  call mpp_varpar_init(10, 15, 5, 20)

  call PetscFinalize(ierr)

end program main
