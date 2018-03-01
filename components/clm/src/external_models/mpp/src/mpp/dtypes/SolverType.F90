module SolverType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>
  use petscvec
  use petscmat
  use petscts
  use petscsnes
  use petscdm
  use petscdmda
  use petscsys
  use petscksp
  !
  ! !USES:
  use mpp_abortutils            , only : endrun
  use mpp_varctl                , only : iulog
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: solver_type
     PetscInt  :: petsc_solver_type            ! type of PETSc equation being solved (KSP, SNES, TS)

     DM        :: dm                           ! PETSc DM
     TS        :: ts                           ! PETSc TS
     SNES      :: snes                         ! PETSc SNES
     KSP       :: ksp                          ! PETSc KSP

     Vec       :: soln                         ! solution at current iteration + time step
     Vec       :: soln_prev                    ! solution vector at previous time step
     Vec       :: soln_prev_clm                ! solution vector at previous CLM time step
     Vec       :: res                          ! residual vector
     Vec       :: rhs                          ! used if SoE is a PETSc TS
     Mat       :: jac                          ! used if SoE is a PETSc TS/SNES
     Mat       :: Amat                         ! used if SoE is a PETSc KSP

     PetscInt  :: cumulative_newton_iterations ! Total number of Newton iterations
     PetscInt  :: cumulative_linear_iterations ! Total number of Linear iterations

     PetscBool :: use_dynamic_linesearch       ! Try another linesearch before cutting timestep

   contains
     procedure, public :: Init
     procedure, public :: GetSolverType
     procedure, public :: SetUseDynLineSearch
     procedure, public :: GetUseDynLineSearch
     procedure, public :: IncrementNewtonIterCount
     procedure, public :: IncrementLinearIterCount
  end type solver_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Init(this)
    !
    ! !DESCRIPTION:
    ! Initialize a solver_type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type) :: this

    this%cumulative_newton_iterations = 0
    this%cumulative_linear_iterations = 0

    this%use_dynamic_linesearch       = PETSC_FALSE

    this%petsc_solver_type            = 0
    this%dm                           = PETSC_NULL_DM
    this%ts                           = PETSC_NULL_TS
    this%snes                         = PETSC_NULL_SNES
    this%ksp                          = PETSC_NULL_KSP

    this%soln                         = PETSC_NULL_VEC
    this%soln_prev                    = PETSC_NULL_VEC
    this%soln_prev_clm                = PETSC_NULL_VEC
    this%res                          = PETSC_NULL_VEC
    this%rhs                          = PETSC_NULL_VEC
    this%jac                          = PETSC_NULL_MAT
    this%Amat                         = PETSC_NULL_MAT

  end subroutine Init

  !------------------------------------------------------------------------
  function GetSolverType(this)
    !
    ! !DESCRIPTION:
    ! Return solver type
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type) :: this
    !
    PetscInt :: GetSolverType

    GetSolverType = this%petsc_solver_type

  end function GetSolverType

  !------------------------------------------------------------------------
  function GetUseDynLineSearch(this)
    !
    ! !DESCRIPTION:
    ! Return setting of using dynamic linesearch
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type) :: this
    !
    PetscBool :: GetUseDynLineSearch

    GetUseDynLineSearch = this%use_dynamic_linesearch

  end function GetUseDynLineSearch

  !------------------------------------------------------------------------
  subroutine SetUseDynLineSearch(this, value)
    !
    ! !DESCRIPTION:
    ! Return setting of using dynamic linesearch
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type)    :: this
    PetscBool, intent(in) :: value

    this%use_dynamic_linesearch = value

  end subroutine SetUseDynLineSearch

  !------------------------------------------------------------------------
  subroutine IncrementNewtonIterCount(this, count)
    !
    ! !DESCRIPTION:
    ! Increase the cumulative newton iteration count
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type)   :: this
    PetscInt, intent(in) :: count

    this%cumulative_newton_iterations = this%cumulative_newton_iterations + count

  end subroutine IncrementNewtonIterCount

  !------------------------------------------------------------------------
  subroutine IncrementLinearIterCount(this, count)
    !
    ! !DESCRIPTION:
    ! Increase the cumulative linear iteration count
    !
    implicit none
    !
    ! !ARGUMENTS
    class(solver_type)   :: this
    PetscInt, intent(in) :: count

    this%cumulative_linear_iterations = this%cumulative_linear_iterations + count

  end subroutine IncrementLinearIterCount
#endif

end module SolverType
