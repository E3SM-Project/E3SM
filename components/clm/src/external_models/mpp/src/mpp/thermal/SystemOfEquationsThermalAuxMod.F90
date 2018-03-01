
module SystemOfEquationsThermalAuxMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl                      , only : iulog
  use mpp_abortutils                      , only : endrun
  use mpp_shr_log_mod                     , only : errMsg => shr_log_errMsg
  use SystemOfEquationsThermalAuxType , only : sysofeqns_thermal_auxvar_type
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  public :: SOEThermalAuxSetRData
  public :: SOEThermalAuxSetIData
  public :: SOEThermalAuxSetBData

contains
  
  !------------------------------------------------------------------------
  subroutine SOEThermalAuxSetRData(auxvars, var_type, nauxvar, auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    ! Used by CLM to set values of boundary conditions and source-sink
    ! terms for the Thermal solver.
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants, only : VAR_LIQ_AREAL_DEN
    use MultiPhysicsProbConstants, only : VAR_ICE_AREAL_DEN
    use MultiPhysicsProbConstants, only : VAR_FRAC
    use MultiPhysicsProbConstants, only : VAR_SNOW_WATER
    use MultiPhysicsProbConstants, only : VAR_NUM_SNOW_LYR
    use MultiPhysicsProbConstants, only : VAR_DHS_DT
    use MultiPhysicsProbConstants, only : VAR_DZ
    use MultiPhysicsProbConstants, only : VAR_DIST_UP
    use MultiPhysicsProbConstants, only : VAR_DIST_DN
    use MultiPhysicsProbConstants, only : VAR_TUNING_FACTOR
    use MultiPhysicsProbConstants, only : VAR_THERMAL_COND
    use MultiPhysicsProbConstants, only : VAR_HEAT_CAP
    use MultiPhysicsProbConstants, only : VAR_BC_SS_CONDITION
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                          :: var_type
    PetscInt                                      :: nauxvar
    PetscInt                                      :: auxvar_offset
    PetscReal, pointer                            :: data_1d(:)
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

    case (VAR_LIQ_AREAL_DEN)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%liq_areal_den = data_1d(iauxvar)
       enddo

    case (VAR_ICE_AREAL_DEN)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%ice_areal_den = data_1d(iauxvar)
       enddo

    case (VAR_FRAC)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%frac = data_1d(iauxvar)
       enddo

    case (VAR_SNOW_WATER)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%snow_water = data_1d(iauxvar)
       enddo

    case (VAR_DHS_DT)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%dhsdT = data_1d(iauxvar)
       enddo

    case (VAR_DZ)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%dz = data_1d(iauxvar)
       enddo

    case (VAR_DIST_UP)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%dist_up = data_1d(iauxvar)
       enddo

    case (VAR_DIST_DN)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%dist_dn = data_1d(iauxvar)
       enddo

    case (VAR_THERMAL_COND)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%therm_cond = data_1d(iauxvar)
       enddo

    case (VAR_HEAT_CAP)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%heat_cap_pva = data_1d(iauxvar)
       enddo

    case (VAR_BC_SS_CONDITION)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%condition_value = data_1d(iauxvar)
       enddo

    case (VAR_TUNING_FACTOR)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%tuning_factor = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'ERROR: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine SOEThermalAuxSetRData

  !------------------------------------------------------------------------
  subroutine SOEThermalAuxSetIData(auxvars, var_type, nauxvar, auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_NUM_SNOW_LYR
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_auxvar_type), pointer :: auxvars(:)
    PetscInt, intent(in)                          :: var_type
    PetscInt                                      :: nauxvar
    PetscInt                                      :: auxvar_offset
    PetscInt, pointer                             :: data_1d(:)
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
    case (VAR_NUM_SNOW_LYR)
       do iauxvar = 1, size(data_1d)
          auxvars(iauxvar + auxvar_offset)%num_snow_layer = data_1d(iauxvar)
       enddo

    case default
       write(iulog,*) 'ERROR: unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))

    end select

  end subroutine SOEThermalAuxSetIData

  !------------------------------------------------------------------------
  subroutine SOEThermalAuxSetBData(auxvars, var_type, nauxvar, auxvar_offset, data_1d)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : VAR_ACTIVE
    !
    implicit none
    !
    ! !ARGUMENTS
    type (sysofeqns_thermal_auxvar_type), pointer :: auxvars(:)
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

  end subroutine SOEThermalAuxSetBData

#endif
end module SystemOfEquationsThermalAuxMod
