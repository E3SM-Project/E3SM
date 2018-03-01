
module SystemOfEquationsThermalAuxType

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

  type, public :: sysofeqns_thermal_auxvar_type

     PetscReal :: temperature     ! [K]

     PetscReal :: liq_areal_den   ! [kg/m^2] = h2osoi_liq
     PetscReal :: ice_areal_den   ! [kg/m^2] = h2osoi_ice
     PetscReal :: frac            ! [-]
     PetscReal :: snow_water      ! snow water (mm H2O)
     PetscInt  :: num_snow_layer  ! number of snow layer
     PetscReal :: dhsdT           ! [W/m^2/K]
     PetscReal :: dz              ! [m]
     PetscReal :: dist_up         ! [m]
     PetscReal :: dist_dn         ! [m]
     PetscReal :: tuning_factor

     PetscBool :: is_active       !
     
     PetscReal :: therm_cond      ! [W/(m K)]
     PetscReal :: heat_cap_pva    ! [J/(m3 K)]

     ! Stores a value if auxvar corresponds to boundary condition or source sink
     PetscReal :: condition_value ! [depends on the type of condition]

     PetscInt  :: goveqn_id       ! [-]
     PetscInt  :: condition_id    ! [-]

     PetscBool :: is_in           ! [True/False] T = auxvar is for an internal control volume
     PetscBool :: is_bc           ! [True/False] T = auxvar is for a boundary condition
     PetscBool :: is_ss           ! [True/False] T = auxvar is for a source/sink condition

   contains
     procedure, public :: Init      => ThermalSOEAuxVarInit
  end type sysofeqns_thermal_auxvar_type

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
    class(sysofeqns_thermal_auxvar_type) :: this

    this%temperature     = 273.15d0
    this%liq_areal_den   = 0.d0
    this%ice_areal_den   = 0.d0
    this%frac            = 0.d0
    this%snow_water      = 0.d0
    this%num_snow_layer  = 0

    this%dz              = 0.d0
    this%dist_up         = 0.d0
    this%dist_dn         = 0.d0
    this%dhsdT           = 0.d0
    this%tuning_factor   = 0.d0

    this%is_active       = PETSC_FALSE
    
    this%therm_cond      = 0.d0
    this%heat_cap_pva    = 0.d0

    this%condition_value = 0.d0
    this%goveqn_id       = 0
    this%condition_id    = 0

    this%is_in           = PETSC_FALSE
    this%is_bc           = PETSC_FALSE
    this%is_ss           = PETSC_FALSE

  end subroutine ThermalSOEAuxVarInit

#endif

end module SystemOfEquationsThermalAuxType
