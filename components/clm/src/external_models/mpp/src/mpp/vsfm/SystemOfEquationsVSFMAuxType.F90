
module SystemOfEquationsVSFMAuxType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl         , only : iulog
  use mpp_abortutils         , only : endrun
  use mpp_shr_log_mod        , only : errMsg => shr_log_errMsg
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: sysofeqns_vsfm_auxvar_type

     PetscReal :: pressure               ! [Pa]
     PetscReal :: dxi_dtime              ! [kg m^{-3} s^{-1}]

     PetscReal :: temperature            ! [K]
     PetscReal :: liq_sat                ! [-]
     PetscReal :: frac_liq_sat           ! [-]
     PetscReal :: mass                   ! [kg]
     PetscReal :: soil_matrix_pot        ! [m]

     ! Stores a value if auxvar corresponds to boundary condition or source sink
     PetscReal :: condition_value        ! [depends on the type of condition]
     PetscReal :: lateral_mass_exchanged ! [kg]
     PetscReal :: boundary_mass_exchanged! [kg]
     PetscReal :: mass_flux              ! [kg/s]

     PetscInt  :: goveqn_id              ! [-]
     PetscInt  :: condition_id           ! [-]

     PetscBool :: is_in                  ! [True/False] T = auxvar is for an internal control volume
     PetscBool :: is_bc                  ! [True/False] T = auxvar is for a boundary condition
     PetscBool :: is_ss                  ! [True/False] T = auxvar is for a source/sink condition
     PetscBool :: is_active              ! [True/False] T = auxvar is for a grid cell that is active

   contains
     procedure, public :: Init     => VSFMSOEAuxVarInit
     procedure, public :: SetValue => VSFMSOEAuxVarSetValue
     procedure, public :: GetValue => VSFMSOEAuxVarGetValue
  end type sysofeqns_vsfm_auxvar_type

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine VSFMSOEAuxVarInit(this)
    !
    ! !DESCRIPTION:
    ! Initialize an auxiliary variable
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_auxvar_type) :: this

    this%pressure                = 0.d0
    this%dxi_dtime               = 0.d0
    this%temperature             = 273.15d0 + 25.d0
    this%liq_sat                 = 0.d0
    this%frac_liq_sat            = 1.d0
    this%mass                    = 0.d0
    this%soil_matrix_pot         = 0.d0
    this%condition_value         = 0.d0
    this%lateral_mass_exchanged  = 0.d0
    this%boundary_mass_exchanged = 0.d0
    this%mass_flux               = 0.d0
    this%goveqn_id               = 0
    this%condition_id            = 0

    this%is_in                   = PETSC_FALSE
    this%is_bc                   = PETSC_FALSE
    this%is_ss                   = PETSC_FALSE
    this%is_active               = PETSC_TRUE

  end subroutine VSFMSOEAuxVarInit

  !------------------------------------------------------------------------
  subroutine VSFMSOEAuxVarSetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    ! Sets value for a particular variable type
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_DXI_DTIME
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants, only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants, only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants, only : VAR_MASS
    use MultiPhysicsProbConstants, only : VAR_SOIL_MATRIX_POT
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_auxvar_type) :: this
    PetscInt, intent(in)              :: var_type
    PetscReal                         :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       this%pressure        = variable_value
    case (VAR_TEMPERATURE)
       this%temperature     = variable_value
    case (VAR_DXI_DTIME)
       this%dxi_dtime       = variable_value
    case (VAR_BC_SS_CONDITION)
       this%condition_value = variable_value
    case (VAR_LIQ_SAT)
       this%liq_sat         = variable_value
    case (VAR_FRAC_LIQ_SAT)
       this%frac_liq_sat    = variable_value
    case (VAR_MASS)
       this%mass            = variable_value
    case (VAR_SOIL_MATRIX_POT)
       this%soil_matrix_pot = variable_value
    case default
       write(iulog,*) 'In VSFMSOEAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEAuxVarSetValue

  !------------------------------------------------------------------------
  subroutine VSFMSOEAuxVarGetValue(this, var_type, variable_value)
    !
    ! !DESCRIPTION:
    ! Returns value for a particular variable type
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_DXI_DTIME
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants, only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants, only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants, only : VAR_MASS
    use MultiPhysicsProbConstants, only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants, only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants, only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants, only : VAR_MASS_FLUX
    !
    implicit none
    !
    ! !ARGUMENTS
    class(sysofeqns_vsfm_auxvar_type) :: this
    PetscInt, intent(in)              :: var_type
    PetscReal, intent(out)            :: variable_value

    select case(var_type)
    case (VAR_PRESSURE)
       variable_value = this%pressure
    case (VAR_TEMPERATURE)
       variable_value = this%temperature
    case (VAR_DXI_DTIME)
       variable_value = this%dxi_dtime
    case (VAR_BC_SS_CONDITION)
       variable_value = this%condition_value
    case (VAR_LIQ_SAT)
       variable_value = this%liq_sat
    case (VAR_FRAC_LIQ_SAT)
       variable_value = this%frac_liq_sat
    case (VAR_MASS)
       variable_value = this%mass
    case (VAR_SOIL_MATRIX_POT)
       variable_value = this%soil_matrix_pot
    case (VAR_LATERAL_MASS_EXCHANGED)
       variable_value = this%lateral_mass_exchanged
    case (VAR_BC_MASS_EXCHANGED)
       variable_value = this%boundary_mass_exchanged
    case (VAR_MASS_FLUX)
       variable_value = this%mass_flux
    case default
       write(iulog,*) 'In VSFMSOEAuxVarGetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine VSFMSOEAuxVarGetValue

#endif

end module SystemOfEquationsVSFMAuxType
