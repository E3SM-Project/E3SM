module RichardsODEPressureConnAuxType

#ifdef USE_PETSC_LIB

  use mpp_varctl          , only : iulog
  use mpp_abortutils      , only : endrun
  use mpp_shr_log_mod     , only : errMsg => shr_log_errMsg
  use SaturationFunction  , only : saturation_params_type
  !
  implicit none
  private

#include <petsc/finclude/petsc.h>

  type, public :: rich_ode_pres_conn_auxvar_type

     PetscReal :: pressure_up ! [Pa]
     PetscReal :: pressure_dn ! [Pa]
     PetscReal :: kr

     PetscInt  :: flux_type
     PetscInt  :: conductance_type
     PetscReal :: conductance
     PetscReal :: conductance_up
     PetscReal :: conductance_dn
     PetscReal :: dkr_dP_up
     PetscReal :: dkr_dP_dn
     PetscReal :: krg
     PetscReal :: dkrg_dP_up
     PetscReal :: dkrg_dP_dn
     PetscReal :: conductance_upwind_weight

     type(saturation_params_type) :: satParams_up
     type(saturation_params_type) :: satParams_dn

   contains
     procedure, public :: Init          => RichODEPressureConnAuxInit
     procedure, public :: SetRealValue  => RichODEPressureConnAuxSetRealValue
     procedure, public :: SetIntValue   => RichODEPressureConnAuxSetIntValue
     procedure, public :: AuxVarCompute => RichODEPressureConnAuxVarCompute
  end type rich_ode_pres_conn_auxvar_type

contains
  
  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize the object
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use petscsys
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) :: this

    this%pressure_up               = 0.d0
    this%pressure_dn               = 0.d0
    this%conductance               = 0.d0
    this%conductance_up            = 0.d0
    this%conductance_dn            = 0.d0

    this%flux_type                 = DARCY_FLUX_TYPE
    this%conductance_type          = CONDUCTANCE_CAMPBELL_TYPE
    this%dkr_dP_up                 = 0.d0
    this%dkr_dP_dn                 = 0.d0
    this%kr                        = 0.d0
    this%krg                       = 0.d0
    this%dkrg_dP_up                = 0.d0
    this%dkrg_dP_dn                = 0.d0
    this%conductance_upwind_weight = 0.d0

    call this%satParams_up%Init()
    call this%satParams_dn%Init()
    
  end subroutine RichODEPressureConnAuxInit

  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetRealValue(this, var_type, variable_value )
    !
    ! !DESCRIPTION:
    ! Set a real value
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_UP
    use MultiPhysicsProbConstants, only : VAR_PRESSURE_DN
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_UP
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_DN
    use MultiPhysicsProbConstants, only : VAR_CAMPBELL_HE
    use MultiPhysicsProbConstants, only : VAR_CAMPBELL_N
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) , intent(inout) :: this
    PetscInt                              , intent(in)    :: var_type
    PetscReal                             , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_CONDUCTANCE)
       this%conductance = variable_value

    case (VAR_PRESSURE_UP)
       this%pressure_up = variable_value

    case (VAR_PRESSURE_DN)
       this%pressure_dn = variable_value

    case (VAR_CONDUCTANCE_UP)
       this%conductance_up = variable_value

    case (VAR_CONDUCTANCE_DN)
       this%conductance_dn = variable_value

    case default
       write(iulog,*) 'In RichODEPressureConnAuxSetRealValue: unknown type ',var_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select
    
     end subroutine RichODEPressureConnAuxSetRealValue

  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxSetIntValue(this, var_type, variable_value )
    !
    ! !DESCRIPTION:
    ! Set an integer value
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_FLUX_TYPE
    use MultiPhysicsProbConstants, only : VAR_CONDUCTANCE_TYPE
    use MultiPhysicsProbConstants, only : DARCY_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_FLUX_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_MANOLI_TYPE
    !
    implicit none
    !
    class(rich_ode_pres_conn_auxvar_type) , intent(inout) :: this
    PetscInt                              , intent(in)    :: var_type
    PetscInt                              , intent(in)    :: variable_value

    select case(var_type)
    case (VAR_FLUX_TYPE)
       select case(variable_value)
       case (DARCY_FLUX_TYPE)
          this%flux_type = variable_value
          
       case (CONDUCTANCE_FLUX_TYPE)
          this%flux_type = variable_value
          
       case default
          write(iulog,*)'Unknown value for VAR_FLUX_TYPE = ',variable_value
          call endrun(msg=errMsg(__FILE__,__LINE__))
       end select

    case (VAR_CONDUCTANCE_TYPE)
       select case(variable_value)
       case (CONDUCTANCE_CAMPBELL_TYPE)
          this%conductance_type = CONDUCTANCE_CAMPBELL_TYPE

       case (CONDUCTANCE_MANOLI_TYPE)
          this%conductance_type = CONDUCTANCE_MANOLI_TYPE

       case (0)
          ! Do nothing

       case default
          write(iulog,*)'Unknown value for VAR_CONDUCTANCE_TYPE = ',variable_value
          call endrun(msg=errMsg(__FILE__,__LINE__))
       end select

    case default
       write(iulog,*) 'In RichODEPressureConnAuxSetIntValue: unknown type ',var_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichODEPressureConnAuxSetIntValue


  !------------------------------------------------------------------------
  subroutine RichODEPressureConnAuxVarCompute(this)
    !
    ! !DESCRIPTION:
    ! Computes various secondary quantities (sat, den, etc) based on
    ! the primary quantity (pressure).
    !
    ! !USES:
    use SaturationFunction       , only : SatFunc_PressToRelPerm
    use MultiPhysicsProbConstants, only : CONDUCTANCE_CAMPBELL_TYPE
    use MultiPhysicsProbConstants, only : CONDUCTANCE_MANOLI_TYPE
    !
    implicit none
    !
    ! !ARGUMENTS
    class(rich_ode_pres_conn_auxvar_type) :: this
    PetscReal                             :: frac_liq_sat
    PetscReal                             :: kr_up, kr_dn
    PetscReal                             :: krg_up, krg_dn
    PetscReal                             :: denom
    PetscReal                             :: dkr_up_dP_up, dkr_dn_dP_dn

    select case(this%conductance_type)
    case(CONDUCTANCE_CAMPBELL_TYPE)
       if (this%satParams_up%relperm_func_type == 0 .and. &
            this%satParams_dn%relperm_func_type == 0  ) then

          this%kr         = 1.d0
          this%dkr_dP_up  = 0.d0
          this%dkr_dP_dn  = 0.d0
          this%krg        = this%kr        * this%conductance
          this%dkrg_dP_up = this%dkr_dP_up * this%conductance
          this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

       else

          if (this%satParams_up%relperm_func_type == 0) then

             frac_liq_sat = 1.d0

             call SatFunc_PressToRelPerm( &
                  this%satParams_dn, this%pressure_dn, &
                  frac_liq_sat, kr_dn, dkr_dn_dP_dn)

             this%kr         = kr_dn
             this%dkr_dP_up  = 0.d0
             this%dkr_dP_dn  = dkr_dn_dP_dn
             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

          elseif (this%satParams_dn%relperm_func_type == 0) then

             call SatFunc_PressToRelPerm( &
                  this%satParams_up, this%pressure_up, &
                  frac_liq_sat, kr_up, dkr_up_dP_up)

             this%kr         = kr_up
             this%dkr_dP_up  = dkr_up_dP_up
             this%dkr_dP_dn  = 0.d0
             this%krg        = this%kr        * this%conductance
             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance

          else

             call SatFunc_PressToRelPerm( &
                  this%satParams_up, this%pressure_up, &
                  frac_liq_sat, kr_up, dkr_up_dP_up)

             call SatFunc_PressToRelPerm( &
                  this%satParams_dn, this%pressure_dn, &
                  frac_liq_sat, kr_dn, dkr_dn_dP_dn)

             this%kr = &
                  (       this%conductance_upwind_weight  * kr_up + &
                  (1.d0 - this%conductance_upwind_weight) * kr_dn)

             this%dkr_dP_up =  &
                  this%conductance_upwind_weight          * dkr_up_dP_up

             this%dkr_dP_dn = &
                  (1.d0 - this%conductance_upwind_weight) * dkr_dn_dP_dn

             this%krg        = this%kr        * this%conductance
             this%dkrg_dP_up = this%dkr_dP_up * this%conductance
             this%dkrg_dP_dn = this%dkr_dP_dn * this%conductance
          end if

       endif

    case (CONDUCTANCE_MANOLI_TYPE)

       call SatFunc_PressToRelPerm(              &
            this%satParams_up, this%pressure_up, &
            frac_liq_sat, kr_up, dkr_up_dP_up)

       call SatFunc_PressToRelPerm(              &
            this%satParams_dn, this%pressure_dn, &
            frac_liq_sat, kr_dn, dkr_dn_dP_dn)

       krg_up          = kr_up * this%conductance_up
       krg_dn          = kr_dn * this%conductance_dn
       denom           = krg_up + krg_dn

       this%krg        = krg_up * krg_dn / denom
       this%dkrg_dP_up = (krg_dn/denom)**2.d0 * dkr_up_dP_up * this%conductance_up
       this%dkrg_dP_dn = (krg_up/denom)**2.d0 * dkr_dn_dP_dn * this%conductance_dn

    case default
       write(iulog,*) 'Unknown conductance_type ', this%conductance_type
       call endrun(msg=errMsg(__FILE__,__LINE__))
    end select

  end subroutine RichODEPressureConnAuxVarCompute

#endif

end module RichardsODEPressureConnAuxType
