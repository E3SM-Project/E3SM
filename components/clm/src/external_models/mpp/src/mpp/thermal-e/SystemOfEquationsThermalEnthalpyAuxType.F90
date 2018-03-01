
module SystemOfEquationsThermalEnthalpyAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_abortutils         , only : endrun
  use mpp_shr_log_mod        , only : errMsg => shr_log_errMsg
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: sysofeqns_thermal_enthalpy_auxvar_type

     PetscReal :: temperature     ! [K]

     PetscReal :: pressure   ! [Pa]

     PetscBool :: is_active       ! [-]
     
     ! Stores a value if auxvar corresponds to boundary condition or source sink
     PetscReal :: condition_value ! [depends on the type of condition]

     PetscInt  :: goveqn_id       ! [-]
     PetscInt  :: condition_id    ! [-]

     PetscBool :: is_in           ! [True/False] T = auxvar is for an internal control volume
     PetscBool :: is_bc           ! [True/False] T = auxvar is for a boundary condition
     PetscBool :: is_ss           ! [True/False] T = auxvar is for a source/sink condition

   contains
     procedure, public :: Init      => ThermalSOEAuxVarInit
  end type sysofeqns_thermal_enthalpy_auxvar_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermalSOEAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_thermal_enthalpy_auxvar_type) :: this

    this%temperature     = 273.15d0
    this%pressure        = 101325.d0

    this%is_active       = PETSC_FALSE

    this%condition_value = 0.d0

    this%goveqn_id       = 0
    this%condition_id    = 0

    this%is_in           = PETSC_FALSE
    this%is_bc           = PETSC_FALSE
    this%is_ss           = PETSC_FALSE

  end subroutine ThermalSOEAuxVarInit

#endif

end module SystemOfEquationsThermalEnthalpyAuxType
