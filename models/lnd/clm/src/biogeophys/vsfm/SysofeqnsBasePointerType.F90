#ifdef USE_PETSC_LIB


module SystemOfEquationsBasePointerType

  ! !USES:
  use SystemOfEquationsBaseType    , only : sysofeqns_base_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscts.h"
#include "finclude/petscts.h90"
#include "finclude/petscksp.h"
#include "finclude/petscksp.h90"

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

  public :: SOEIFunction, &
       SOEIJacobian, &
       SOEResidual,  &
       SOEJacobian,  &
       SOEComputeRHS, &
       SOEComputeOperators
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine SOEIFunction(ts, t, U, Udot, F, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc TS
    !
    implicit none

    TS                                :: ts
    PetscReal                         :: t
    Vec                               :: U
    Vec                               :: Udot
    Vec                               :: F
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%IFunction(ts, t, U, Udot, F, ierr)

  end subroutine SOEIFunction

  !------------------------------------------------------------------------
  subroutine SOEIJacobian(ts, t, U, Udot, shift, A, B, this, ierr)
    !
    ! !DESCRIPTION:
    ! Wrapper subroutine when SoE uss PETSc TS
    !
    implicit none

    TS                                :: ts
    PetscReal                         :: t
    Vec                               :: U
    Vec                               :: Udot
    PetscReal                         :: shift
    Mat                               :: A
    Mat                               :: B
    type(sysofeqns_base_pointer_type) :: this
    PetscErrorCode                    :: ierr

    call this%ptr%IJacobian(ts, t, U, Udot, shift, A, B, ierr)

  end subroutine SOEIJacobian

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

end module SystemOfEquationsBasePointerType
#endif
