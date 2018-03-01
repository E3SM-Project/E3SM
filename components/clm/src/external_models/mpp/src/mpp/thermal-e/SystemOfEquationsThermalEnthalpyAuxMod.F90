
module SystemOfEquationsThermalEnthalpyAuxMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                              , only : iulog
  use mpp_abortutils                          , only : endrun
  use mpp_shr_log_mod                         , only : errMsg => shr_log_errMsg
  use SystemOfEquationsThermalEnthalpyAuxType , only : sysofeqns_thermal_enthalpy_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  public :: SOEThermalEnthalpyAuxSetRData
  public :: SOEThermalEnthalpyAuxGetRData
  public :: SOEThermalEnthalpyAuxSetBData

contains
  
  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyAuxSetRData(auxvars, var_type, nauxvar, auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                          :: var_type
    PetscInt                                      :: nauxvar
    PetscInt                                      :: auxvar_offset
    PetscReal                                     :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                      :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%temperature = data_1d(iauxvar)
       enddo

    case (VAR_PRESSURE)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%pressure = data_1d(iauxvar)
       enddo

    case (VAR_BC_SS_CONDITION)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%condition_value = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'ERROR: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine SOEThermalEnthalpyAuxSetRData

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyAuxGetRData(auxvars, var_type, nauxvar, &
       auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    ! Extracts data from the SoE auxvar
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                                   :: var_type
    PetscInt                                               :: nauxvar
    PetscInt                                               :: auxvar_offset
    PetscReal                                              :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                               :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(iauxvar + auxvar_offset)%temperature
       enddo

    case (VAR_PRESSURE)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(iauxvar + auxvar_offset)%pressure
       enddo

    case (VAR_BC_SS_CONDITION)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(iauxvar + auxvar_offset)%condition_value
       enddo

    case default
       write(iulog,*) 'ERROR: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine SOEThermalEnthalpyAuxGetRData

  !------------------------------------------------------------------------
  subroutine SOEThermalEnthalpyAuxSetBData(auxvars, var_type, nauxvar, auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_ACTIVE
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_enthalpy_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                          :: var_type
    PetscInt                                      :: nauxvar
    PetscInt                                      :: auxvar_offset
    PetscBool, pointer                            :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                      :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_ACTIVE)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%is_active = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'ERROR: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine SOEThermalEnthalpyAuxSetBData

#endif
end module SystemOfEquationsThermalEnthalpyAuxMod
