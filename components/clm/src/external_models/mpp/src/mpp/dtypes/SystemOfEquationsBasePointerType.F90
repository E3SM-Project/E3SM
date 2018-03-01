
module SystemOfEquationsBasePointerType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  use petscvec
  use petscmat
  use petscts
  use petscksp

  ! !USES:
  use SystemOfEquationsBaseType    , only : sysofeqns_base_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private


  ! PETSc requires the context (ctx) for procedures be declared as a "type"
  ! instead of a "class". Thus, this object is a workaround and subroutines
  ! defined here essentially call corresponding subroutines in
  ! SystemOfEquationsBaseType.F90
  !
  ! The approach implemented here is similar to the one taken in the
  ! PFLOTRAN (https://bitbucket.org/pflotran/pflotran-dev).

  type, public :: sysofeqns_base_pointer_type
     class(sysofeqns_base_type), pointer :: ptr
  end type sysofeqns_base_pointer_type

  public :: SOEResidual,  &
       SOEJacobian,  &
       SOEComputeRHS, &
       SOEComputeOperators
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine SOEResidual(snes, X, R, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc SNES
    !
    implicit none

    SNES                              :: snes
    Vec                               :: X
    Vec                               :: R
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%Residual(snes, X, R, ierr)

  end subroutine SOEResidual

  !------------------------------------------------------------------------
  subroutine SOEJacobian(snes, X, A, B, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc SNES
    !
    implicit none

    SNES                              :: snes
    Vec                               :: X
    Mat                               :: A
    Mat                               :: B
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%Jacobian(snes, X, A, B, ierr)

  end subroutine SOEJacobian

  !------------------------------------------------------------------------
  subroutine SOEComputeRHS(ksp, R, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc TS
    !
    implicit none

    KSP                               :: ksp
    Vec                               :: R
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%ComputeRHS(ksp, R, ierr)

  end subroutine SOEComputeRHS

  !------------------------------------------------------------------------
  subroutine SOEComputeOperators(ksp, A, B, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc KSP
    !
    implicit none

    KSP                               :: ksp
    Mat                               :: A
    Mat                               :: B
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%ComputeOperators(ksp, A, B, ierr)

  end subroutine SOEComputeOperators

#endif

end module SystemOfEquationsBasePointerType
