module ThermalEnthalpySoilAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use PorosityFunctionMod , only : porosity_params_type
  use SaturationFunction  , only : saturation_params_type
  use RichardsODEPressureAuxType, only : rich_ode_pres_auxvar_type
  !
  implicit none
  !
  private

  type, public, extends(rich_ode_pres_auxvar_type) :: therm_enthalpy_soil_auxvar_type

     PetscReal                    :: ul                       ! internal energy liquid [J kmol^{-1}]
     PetscReal                    :: hl                       ! enthalpy liquid [J kmol^{-3}]

     PetscReal                    :: dul_dP                   ! [J kmol^{-1} Pa^{-1}]
     PetscReal                    :: dhl_dP                   ! [J kmol^{-3} Pa^{-1}]

     PetscReal                    :: dul_dT                   ! [J kmol^{-1} K^{-1}]
     PetscReal                    :: dhl_dT                   ! [J kmol^{-3} K^{-1}]
     PetscReal                    :: dkr_dT                   ! [K^{-1}]
     PetscReal                    :: dsat_dT                  ! [K^{-1}]

     PetscReal                    :: Kel                      ! Kersten number liquid [-]
     PetscReal                    :: therm_cond_wet           ! wet thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: therm_cond_dry           ! dry thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: therm_cond               ! thermal conductivity [J s^{-1} m^{-3} K^{-1}]
     PetscReal                    :: therm_alpha

     PetscReal                    :: den_soil                 ! [kg m^{-3}]
     PetscReal                    :: heat_cap_soil            ! [J kg^{-1} K^{-1}]

     PetscInt                     :: int_energy_enthalpy_type ! [-]

   contains
     procedure, public :: Init    => ThermEnthalpyAuxVarInit
     procedure, public :: Set     => ThermEnthalpyAuxVarSet
     procedure, public :: Get     => ThermEnthalpyAuxVarGet
     procedure, public :: Compute => ThermEnthalpyAuxVarCompute     
  end type therm_enthalpy_soil_auxvar_type

  public :: ThermEnthalpyAuxVarCopy

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarInit(this)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use EOSWaterMod               , only : INT_ENERGY_ENTHALPY_CONSTANT
    use RichardsODEPressureAuxType, only : RichODEPressureAuxVarInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this

    call RichODEPressureAuxVarInit(this)
    
    this%ul                       = 0.d0
    this%hl                       = 0.d0

    this%dul_dT                   = 0.d0
    this%dhl_dT                   = 0.d0
    this%dkr_dT                   = 0.d0
    this%dsat_dT                  = 0.d0

    this%dul_dP                   = 0.d0
    this%dhl_dP                   = 0.d0

    this%Kel                      = 0.d0
    this%therm_cond_wet           = 0.d0
    this%therm_cond_dry           = 0.d0
    this%therm_cond               = 0.d0
    this%therm_alpha              = 0.d0

    this%perm(:)                  = 8.3913d-12

    this%den_soil                 = 0.d0
    this%heat_cap_soil            = 0.d0

    this%int_energy_enthalpy_type = INT_ENERGY_ENTHALPY_CONSTANT

  end subroutine ThermEnthalpyAuxVarInit

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarCopy(this, auxvar)
    !
    ! !DESCRIPTION:
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this
    class(therm_enthalpy_soil_auxvar_type) :: auxvar

    this%temperature              = auxvar%temperature
    this%pressure                 = auxvar%pressure

    this%ul                       = auxvar%ul
    this%hl                       = auxvar%hl
    this%den                      = auxvar%den
    this%vis                      = auxvar%vis
    this%sat                      = auxvar%sat
    this%kr                       = auxvar%kr

    this%dul_dT                   = auxvar%dul_dT
    this%dhl_dT                   = auxvar%dhl_dT
    this%dden_dT                  = auxvar%dden_dT
    this%dvis_dT                  = auxvar%dvis_dT
    this%dsat_dT                  = auxvar%dsat_dT
    this%dkr_dT                   = auxvar%dkr_dT

    this%dul_dP                   = auxvar%dul_dP
    this%dhl_dP                   = auxvar%dhl_dP
    this%dden_dP                  = auxvar%dden_dP
    this%dvis_dP                  = auxvar%dvis_dP
    this%dsat_dP                  = auxvar%dsat_dP
    this%dkr_dP                   = auxvar%dkr_dP
    this%dpor_dP                  = auxvar%dpor_dP

    this%Kel                      = auxvar%Kel
    this%therm_cond_wet           = auxvar%therm_cond_wet
    this%therm_cond_dry           = auxvar%therm_cond_dry
    this%therm_cond               = auxvar%therm_cond
    this%por                      = auxvar%por
    this%perm(:)                  = auxvar%perm(:)
    this%therm_alpha              = auxvar%therm_alpha

    this%den_soil                 = auxvar%den_soil
    this%heat_cap_soil            = auxvar%heat_cap_soil

    this%density_type             = auxvar%density_type
    this%int_energy_enthalpy_type = auxvar%int_energy_enthalpy_type

    this%condition_value          = auxvar%condition_value

    call this%porParams%Copy(auxvar%porParams)
    call this%satParams%Copy(auxvar%satParams)

  end subroutine ThermEnthalpyAuxVarCopy

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarSet(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) , intent(inout) :: this
    PetscInt                               , intent(in)    :: var_type
    PetscReal                              , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_TEMPERATURE)
       this%temperature     = variable_value

    case (VAR_BC_SS_CONDITION)
       this%condition_value = variable_value

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermEnthalpyAuxVarSet

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarGet(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) , intent(in)  :: this
    PetscInt                               , intent(in)  :: var_type
    PetscReal                              , intent(out) :: variable_value

    select case(var_type)
    case (VAR_TEMPERATURE)
       variable_value = this%temperature

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine ThermEnthalpyAuxVarGet

  !------------------------------------------------------------------------
  subroutine ThermEnthalpyAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    !
    use EOSWaterMod               , only : Density
    use EOSWaterMod               , only : Viscosity
    use PorosityFunctionMod       , only : PorosityFunctionComputation
    use SaturationFunction        , only : SatFunc_PressToSat
    use SaturationFunction        , only : SatFunc_PressToRelPerm
    use EOSWaterMod               , only : InternalEnergyAndEnthalpy
    use MultiPhysicsProbConstants , only : FMWH2O, PRESSURE_REF
    !
    implicit none
    !
    ! !ARGUMENTS
    class(therm_enthalpy_soil_auxvar_type) :: this
    !
    PetscReal            :: pressure
    PetscReal, parameter :: frac_liq_sat = 1.d0

    ! Compute saturation
    call SatFunc_PressToSat(this%satParams , this%pressure, this%sat, &
         this%dsat_dP)

    ! Compute relative permeability
    call SatFunc_PressToRelPerm(this%satParams, this%pressure, frac_liq_sat, &
         this%kr, this%dkr_dP)

    ! Compute porosity
    call PorosityFunctionComputation(this%porParams, this%pressure, this%por, &
         this%dpor_dP)

    pressure = this%pressure
    if (this%pressure < PRESSURE_REF)  pressure = PRESSURE_REF

    ! Compute density
    call Density(pressure, this%temperature, this%density_type, &
         this%den, this%dden_dP, this%dden_dT)

    ! Compute viscosity
    call Viscosity(pressure, this%temperature, this%vis, this%dvis_dP, &
         this%dvis_dT)

    ! Compute internal energy and enthalpy
    call InternalEnergyAndEnthalpy(pressure, this%temperature, &
         this%int_energy_enthalpy_type, this%den*FMWH2O, &
         this%dden_dT*FMWH2O, this%dden_dP*FMWH2O, &
         this%ul, this%hl, this%dul_dT, this%dhl_dT, &
         this%dul_dP, this%dhl_dP)

    this%Kel        = (this%sat + 1.d-6 )**(this%therm_alpha)
    this%therm_cond = this%therm_cond_wet*this%Kel + &
                      this%therm_cond_dry*(1.d0 - this%Kel)

  end subroutine ThermEnthalpyAuxVarCompute

#endif

end module ThermalEnthalpySoilAuxType
