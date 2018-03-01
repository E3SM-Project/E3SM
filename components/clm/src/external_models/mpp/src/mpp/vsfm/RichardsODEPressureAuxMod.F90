module RichardsODEPressureAuxMod

#ifdef USE_PETSC_LIB

  use mpp_varctl                 , only : iulog
  use mpp_abortutils             , only : endrun
  use mpp_shr_log_mod            , only : errMsg => shr_log_errMsg
  use RichardsODEPressureAuxType , only : rich_ode_pres_auxvar_type
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include <petsc/finclude/petsc.h>

  public :: RichardsODEPressureAuxVarSetRValues
  public :: RichardsODEPressureAuxVarGetRValues

contains
  
  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarSetRValues(auxvars, var_type, nauxvar, &
       auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_TEMPERATURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(rich_ode_pres_auxvar_type) , pointer    :: auxvars(:)
    PetscInt                        , intent(in) :: var_type
    PetscInt                        , intent(in) :: nauxvar
    PetscInt                        , pointer    :: auxvar_ids(:)
    PetscReal                       , pointer    :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                     :: iauxvar

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


    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureAuxVarSetRValues

  
  !------------------------------------------------------------------------
  subroutine RichardsODEPressureAuxVarGetRValues(auxvars, var_type, nauxvar, &
       auxvar_ids, data_1d)
    !
    ! !DESCRIPTION:
    !
    use MultiPhysicsProbConstants, only : VAR_PRESSURE
    !
    implicit none
    !
    ! !ARGUMENTS
    type(rich_ode_pres_auxvar_type) , pointer    :: auxvars(:)
    PetscInt                        , intent(in) :: var_type
    PetscInt                        , intent(in) :: nauxvar
    PetscInt                        , pointer    :: auxvar_ids(:)
    PetscReal                       , pointer    :: data_1d(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                     :: iauxvar

    if (size(data_1d) > nauxvar) then
       write(iulog,*) 'ERROR: size(data_1d) > nauxvar'
       write(iulog,*) 'size(data_1d) = ',size(data_1d)
       write(iulog,*) 'nauxvar       = ', nauxvar
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    select case(var_type)
    case (VAR_PRESSURE)
       do iauxvar = 1, size(data_1d)
          data_1d(iauxvar) = auxvars(auxvar_ids(iauxvar))%pressure
       enddo

    case default
       write(iulog,*) 'Unknown var_type'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine RichardsODEPressureAuxVarGetRValues

#endif

end module RichardsODEPressureAuxMod
