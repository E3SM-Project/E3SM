
module ThermalEnthalpySoilAuxMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  use mpp_varctl                 , only : iulog
  use mpp_abortutils             , only : endrun
  use mpp_shr_log_mod            , only : errMsg => shr_log_errMsg
  use ThermalEnthalpySoilAuxType , only : therm_enthalpy_soil_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  public :: ThermalEnthalpySoilAuxVarSetRValues
  public :: ThermalEnthalpySoilAuxVarGetRValues

contains
  
  !------------------------------------------------------------------------
  subroutine ThermalEnthalpySoilAuxVarSetRValues(auxvars, var_type, nauxvar, auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(therm_enthalpy_soil_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                           :: var_type
    PetscInt                                       :: nauxvar
    PetscInt, pointer                              :: auxvar_ids(:)
    PetscReal, pointer                             :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                       :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, size(data_1d)
          auxvars(auxvar_ids(iauxvar))%temperature = data_1d(iauxvar)
       enddo

    case (VAR_PRESSURE)
       do iauxvar = 1, size(data_1d)
          auxvars(auxvar_ids(iauxvar))%pressure = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalEnthalpySoilAuxVarSetRValues

  
  !------------------------------------------------------------------------
  subroutine ThermalEnthalpySoilAuxVarGetRValues(auxvars, var_type, nauxvar, auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_DZ
    use MultiPhysicsProbConstants, only : VAR_THERMAL_COND
    !
    implicit none
    !
    ! !ARGUMENTS
    type(therm_enthalpy_soil_auxvar_type) , pointer    :: auxvars(:)
    PetscInt                              , intent(in) :: var_type
    PetscInt                                           :: nauxvar
    PetscInt                              , pointer    :: auxvar_ids(:)
    PetscReal                             , pointer    :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                           :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_TEMPERATURE)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%temperature
       enddo

    case default
       write(iulog,*) 'In ThermKSPTempSoilAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermalEnthalpySoilAuxVarGetRValues

#endif

end module ThermalEnthalpySoilAuxMod
