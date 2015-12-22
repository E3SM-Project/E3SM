#ifdef USE_PETSC_LIB


module RichardsODEPressureAuxType

  ! !USES:
  use clm_varctl          , only : iulog
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use PorosityFunctionMod , only : porosity_params_type
  use SaturationFunction  , only : saturation_params_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  type, public :: rich_ode_pres_auxvar_type

     ! primary unknown independent variable
     PetscReal :: pressure        ! [Pa]

     PetscReal :: pressure_prev   ! [Pa]

     ! independent variable available from:
     !  - another governing equation, or
     !  - another system-of-equation
     PetscReal :: temperature     ! [K]
     PetscReal :: frac_liq_sat    ! [-]

     ! If the auxvar corresponds to boundary condition
     ! or source sink, the value is stored in this
     ! variable
     PetscReal :: condition_value ! Depends

     ! parameters
     PetscReal :: perm(3)         ! [m^2]
     PetscReal :: por             ! [-]
     PetscInt  :: density_type    ! [-]
     PetscReal :: conductance     ! [s^{-1}]

     ! derived quantities = f(state_variables, parameters)
     PetscReal :: vis             ! [Pa s]
     PetscReal :: kr              ! [-]
     PetscReal :: sat             ! [-]
     PetscReal :: den             ! [kg m^{-3}]

     PetscReal :: dvis_dP         ! [s]
     PetscReal :: dpor_dP         ! [Pa^{-1}]
     PetscReal :: dkr_dP          ! [Pa^{-1}]
     PetscReal :: dsat_dP         ! [Pa^{-1}]
     PetscReal :: dden_dP         ! [kg m^{-3} Pa^{-1}]

     type(porosity_params_type)   :: porParams
     type(saturation_params_type) :: satParams

   contains
     procedure, public :: Init          => RichODEPressureAuxVarInit
     procedure, public :: Copy          => RichODEPressureAuxVarCopy
     procedure, public :: SetValue      => RichODEPressureAuxVarSetValue
     procedure, public :: GetValue      => RichODEPressureAuxVarGetValue
     procedure, public :: AuxVarCompute => RichODEPressureAuxVarCompute
  end type rich_ode_pres_auxvar_type

  !------------------------------------------------------------------------
contains


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    ! !USES:
    use EOSWaterMod         , only : DENSITY_CONSTANT
    use PorosityFunctionMod , only : PorosityFunctionInit
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this

    this%pressure                = 0.d0
    this%temperature             = 273.15d0 + 25.d0
    this%frac_liq_sat            = 1.d0
    this%condition_value         = 0.d0
    this%pressure_prev           = 3.5355d3
    this%perm(:)                 = 0.d0
    this%por                     = 0.d0
    this%conductance             = 0.d0

    this%satParams%sat_func_type = 0
    this%satParams%sat_res       = 0.d0
    this%satParams%alpha         = 0.d0
    this%satParams%bc_lambda     = 0.d0
    this%satParams%vg_m          = 0.d0
    this%satParams%vg_n          = 0.d0

    call PorosityFunctionInit(this%porParams)

    this%vis                     = 0.d0
    this%dvis_dP                 = 0.d0
    this%dpor_dP                 = 0.d0
    this%dkr_dP                  = 0.d0
    this%dsat_dP                 = 0.d0
    this%dden_dP                 = 0.d0
    this%density_type            = DENSITY_CONSTANT

  end subroutine RichODEPressureAuxVarInit

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarCopy(this, auxvar)
    !
    ! !DESCRIPTION:
    ! Copies content of `auxvar`
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this
    class(rich_ode_pres_auxvar_type) :: auxvar

    this%pressure                = auxvar%pressure
    this%temperature             = auxvar%temperature
    this%frac_liq_sat            = auxvar%frac_liq_sat
    this%condition_value         = auxvar%condition_value
    this%pressure_prev           = auxvar%pressure_prev
    this%perm(:)                 = auxvar%perm(:)
    this%por                     = auxvar%por
    this%conductance             = auxvar%conductance

    this%satParams%sat_func_type = auxvar%satParams%sat_func_type
    this%satParams%sat_res       = auxvar%satParams%sat_res
    this%satParams%alpha         = auxvar%satParams%alpha
    this%satParams%bc_lambda     = auxvar%satParams%bc_lambda
    this%satParams%vg_m          = auxvar%satParams%vg_m
    this%satParams%vg_n          = auxvar%satParams%vg_n

    call this%porParams%Copy(auxvar%porParams)
    call this%satParams%Copy(auxvar%satParams)

    this%vis                     = auxvar%vis
    this%dvis_dP                 = auxvar%dvis_dP
    this%dpor_dP                 = auxvar%dpor_dP
    this%dkr_dP                  = auxvar%dkr_dP
    this%dsat_dP                 = auxvar%dsat_dP
    this%dden_dP                 = auxvar%dden_dP
    this%density_type            = auxvar%density_type

  end subroutine RichODEPressureAuxVarCopy

  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarSetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_PREV
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type), intent(inout) :: this
    PetscInt, intent(in)                            :: var_type
    PetscReal, intent(in)                           :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       this%pressure        = variable_value
    case (VAR_TEMPERATURE)
       this%temperature     = variable_value
    case (VAR_PRESSURE_PREV)
       this%pressure_prev   = variable_value
    case (VAR_BC_SS_CONDITION)
       this%condition_value = variable_value
    case default
       write(iulog,*) 'In RichODEPressureAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichODEPressureAuxVarSetValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarGetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants, only : VAR_LIQ_SAT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type) :: this
    PetscInt, intent(in)             :: var_type
    PetscReal, intent(out)           :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       variable_value = this%pressure
    case (VAR_TEMPERATURE)
       variable_value = this%temperature
    case (VAR_BC_SS_CONDITION)
       variable_value = this%condition_value
    case (VAR_LIQ_SAT)
       variable_value = this%sat
    case default
       write(iulog,*) 'In RichODEPressureAuxVarGetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichODEPressureAuxVarGetValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    ! Computes various secondary quantities (sat, den, etc) based on
    ! the primary quantity (pressure).
    !
    ! !USES:
    use EOSWaterMod           , only : Density
    use EOSWaterMod           , only : Viscosity
    use PorosityFunctionMod   , only : PorosityFunctionComputation
    use SaturationFunction    , only : SatFunc_PressToSat
    use SaturationFunction    , only : SatFunc_PressToRelPerm
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_auxvar_type)   :: this
    !
    ! LOCAL VARIABLES

    ! Compute saturation
    call SatFunc_PressToSat(this%satParams , &
         this%pressure                     , &
         this%sat                          , &
         this%dsat_dP                        &
         )

    ! Compute relative permeability
    call SatFunc_PressToRelPerm(this%satParams , &
         this%pressure                         , &
         this%frac_liq_sat                     , &
         this%kr                               , &
         this%dkr_dP                             &
         )

    ! Compute density
    call Density(this%pressure                      , &
         this%temperature                           , &
         this%density_type                          , &
         this%den                                   , &
         this%dden_dP                                 &
         )

    ! Compute viscosity
    call Viscosity(this%pressure                    , &
         this%temperature                           , &
         this%vis                                   , &
         this%dvis_dP                                 &
         )

    ! Compute porosity
    call PorosityFunctionComputation(this%porParams , &
         this%pressure                              , &
         this%por                                   , &
         this%dpor_dP                                 &
         )

  end subroutine RichODEPressureAuxVarCompute


end module RichardsODEPressureAuxType

#endif
