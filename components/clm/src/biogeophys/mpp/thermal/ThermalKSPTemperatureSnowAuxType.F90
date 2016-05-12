#ifdef USE_PETSC_LIB

module  ThermalKSPTemperatureSnowAuxType
  !
  ! !USES:
  use ThermalKSPTemperatureBaseAuxType
  !
  ! !PUBLIC TYPES:
  implicit none

  private

#include "finclude/petscsys.h"

  type, public, extends(therm_ksp_temp_base_auxvar_type)  :: therm_ksp_temp_snow_auxvar_type

     PetscReal :: liq_areal_den       ! [kg/m^2] = h2osoi_liq
     PetscReal :: ice_areal_den       ! [kg/m^2] = h2osoi_ice
     PetscInt  :: num_snow_layer      ! number of snow layer
     PetscReal :: dz                  ! [m]
     PetscReal :: tuning_factor       !

   contains
     procedure, public :: Init                  => ThermKSPTempSnowAuxVarInit
     procedure, public :: AuxVarCompute         => ThermKSPTempSnowAuxVarCompute
  end type therm_ksp_temp_snow_auxvar_type

contains

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_snow_auxvar_type)   :: this

    call this%BaseInit()

    this%liq_areal_den            = 0.d0
    this%ice_areal_den            = 0.d0
    this%num_snow_layer           = 0
    this%dz = 0.d0
    this%tuning_factor      = 1.d0

  end subroutine ThermKSPTempSnowAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarCompute(this, dz, vol)
    !
    ! !DESCRIPTION:
    !
    use mpp_varcon      , only : cpice,  cpliq, tkair, tkice
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_ksp_temp_snow_auxvar_type)   :: this
    PetscReal                                :: dz
    PetscReal                                :: vol
    !
    ! LOCAL VARIABLES
    PetscReal :: bw

    if (.not.this%is_active) then
       return
    else
       bw = (this%ice_areal_den + this%liq_areal_den)/(this%frac * dz)
       this%therm_cond   = tkair + (7.75d-5*bw + 1.105d-6*bw*bw)*(tkice-tkair)
       this%heat_cap_pva = cpliq*this%liq_areal_den + cpice*this%ice_areal_den
    endif

    this%heat_cap_pva = this%heat_cap_pva/dz
  end subroutine ThermKSPTempSnowAuxVarCompute

end module ThermalKSPTemperatureSnowAuxType

#endif
