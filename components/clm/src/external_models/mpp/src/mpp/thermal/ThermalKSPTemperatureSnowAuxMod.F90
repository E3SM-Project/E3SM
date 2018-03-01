module ThermalKSPTemperatureSnowAuxMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                      , only : iulog
  use mpp_abortutils                      , only : endrun
  use mpp_shr_log_mod                     , only : errMsg => shr_log_errMsg
  use ThermalKSPTemperatureSnowAuxType, only : therm_ksp_temp_snow_auxvar_type
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  public :: ThermKSPTempSnowAuxVarSetRValues
  public :: ThermKSPTempSnowAuxVarGetRValues

contains
  
  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarSetRValues(auxvars, var_type, nauxvar, auxvar_ids, data_1d)
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
    type(therm_ksp_temp_snow_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                   :: var_type
    PetscInt                               :: nauxvar
    PetscInt, pointer                               :: auxvar_ids(:)
    PetscReal, pointer                              :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: iauxvar

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

    case (VAR_DZ)
       do iauxvar = 1, size(data_1d)
          auxvars(auxvar_ids(iauxvar))%dz = data_1d(iauxvar)
       enddo

    case (VAR_THERMAL_COND)
       do iauxvar = 1, size(data_1d)
          auxvars(auxvar_ids(iauxvar))%therm_cond = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'In ThermKSPTempSnowAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSnowAuxVarSetRValues

  
  !------------------------------------------------------------------------
  subroutine ThermKSPTempSnowAuxVarGetRValues(auxvars, var_type, nauxvar, auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_FRAC
    use MultiPhysicsProbConstants, only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants, only : VAR_ACTIVE
    use MultiPhysicsProbConstants, only : VAR_DZ
    use MultiPhysicsProbConstants, only : VAR_DIST_UP
    !
    implicit none
    !
    ! !ARGUMENTS
    type(therm_ksp_temp_snow_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                   :: var_type
    PetscInt                               :: nauxvar
    PetscInt, pointer                      :: auxvar_ids(:)
    PetscReal, pointer                         :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: iauxvar

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

    case (VAR_FRAC)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%frac
       enddo

    case (VAR_THERMAL_COND)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%therm_cond
       enddo

    case (VAR_ACTIVE)
       do iauxvar = 1, size(data_1d)
          if (auxvars(auxvar_ids(iauxvar))%is_active) then
             data_1d(iauxvar) = 1.d0
          else
             data_1d(iauxvar) = 0.d0
          endif
       enddo

    case (VAR_DZ)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%dz
       enddo

    case (VAR_DIST_UP)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%dist_up
       enddo

    case default
       write(iulog,*) 'In ThermKSPTempSnowAuxVarSetValue: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine ThermKSPTempSnowAuxVarGetRValues

#endif
end module ThermalKSPTemperatureSnowAuxMod
