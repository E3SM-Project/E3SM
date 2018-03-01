module ThermalKSPTemperatureBaseAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  !
  ! !PUBLIC TYPES:
  use petscsys

  implicit none

  private

  type, public :: therm_ksp_temp_base_auxvar_type

     ! primary unknown independent variable
     PetscReal :: temperature     ! [K]

     PetscBool :: is_active       !

     ! independent variable available from:
     !  - another governing equation, or
     !  - another system-of-equation
     PetscReal :: frac            ! [-]
     PetscReal :: dhsdT           ! [W/m^2/K]
     PetscReal :: dist_up         ! [m]
     PetscReal :: dist_dn         ! [m]

     ! derived quantities = f(state_variables, parameters)
     PetscReal :: therm_cond      ! [W/(m K)]
     PetscReal :: heat_cap_pva    ! [J/(m3 K)]

     ! If the auxvar corresponds to boundary condition
     ! or source sink, the value is stored in this
     ! variable
     PetscReal :: condition_value ! Depends

   contains
     procedure, public :: BaseInit               => ThermKSPTempBaseAuxVarInit
  end type therm_ksp_temp_base_auxvar_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempBaseAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_base_auxvar_type)   :: this

    this%temperature     = 273.15d0
    this%is_active       = PETSC_FALSE

    this%dhsdT           = 0.d0
    this%frac            = 0.d0
    this%dist_up         = 0.d0
    this%dist_dn         = 0.d0

    this%therm_cond      = 0.d0
    this%heat_cap_pva    = 0.d0
    this%condition_value = 0.d0

  end subroutine ThermKSPTempBaseAuxVarInit

#endif

end module ThermalKSPTemperatureBaseAuxType

