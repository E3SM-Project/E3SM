#ifdef USE_PETSC_LIB


module PorosityFunctionMod

  ! !USES:
  use clm_varctl         , only : iulog
  use abortutils         , only : endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  PetscInt, parameter, public :: POROSITY_CONSTANT    = 1
  PetscInt, parameter, public :: POROSITY_LINEAR      = 2

  type, public :: porosity_params_type

     PetscInt  :: porosityId         ! Identify porosity model using, e.g., `POROSITY_LINEAR`
     PetscReal :: porosity_base      ! Base Porosity [-]

     PetscReal :: pressure_reference ! Reference pressure used in prosity models [Pa]
     PetscReal :: lin_mod_slope      ! Parameter sepecific to linear porosity model [Pa^{-1}]

   contains
     procedure, public :: Copy     => PorosityFunctionCopy
  end type porosity_params_type

  public :: PorosityFunctionInit
  public :: PorosityFunctionSetConstantModel
  public :: PorosityFunctionSetLinearModel
  public :: PorosityFunctionComputation

contains


  !------------------------------------------------------------------------
  subroutine PorosityFunctionInit (this)
    !
    ! !DESCRIPTION:
    ! Initializes
    !
    implicit none
    !
    ! !ARGUMENTS
    type(porosity_params_type), intent (out) :: this

    this%porosityId         = -1
    this%porosity_base      = 0.d0
    this%pressure_reference = 0.d0
    this%lin_mod_slope      = 0.d0

  end subroutine PorosityFunctionInit

  !------------------------------------------------------------------------
  subroutine PorosityFunctionSetConstantModel (porParams, por_base)
    !
    ! !DESCRIPTION:
    ! Sets values corresponding to a constant porosity model.
    !
    implicit none
    !
    ! !ARGUMENTS
    type(porosity_params_type), intent (out) :: porParams
    PetscReal, intent (in)                   :: por_base

    porParams%porosityId    = POROSITY_CONSTANT
    porParams%porosity_base = por_base

  end subroutine PorosityFunctionSetConstantModel

  !------------------------------------------------------------------------
  subroutine PorosityFunctionSetLinearModel (porParams, por_base, &
       press_base, slope)
    !
    ! !DESCRIPTION:
    ! Sets values corresponding to a linear porosity model.
    !
    implicit none
    !
    ! !ARGUMENTS
    type(porosity_params_type), intent (out) :: porParams
    PetscReal, intent (in)                   :: por_base
    PetscReal, intent (in)                   :: press_base
    PetscReal, intent (in)                   :: slope

    porParams%porosityId         = POROSITY_LINEAR
    porParams%porosity_base      = por_base
    porParams%pressure_reference = press_base
    porParams%lin_mod_slope      = slope

  end subroutine PorosityFunctionSetLinearModel

  !------------------------------------------------------------------------
  subroutine PorosityFunctionComputation (porParams, P, por, dpor_dP)
    !
    ! !DESCRIPTION:
    ! Given a pressure value, returns  porosity and derivative of prosity
    ! w.r.t pressure
    !
    implicit none
    !
    ! !ARGUMENTS
    type(porosity_params_type), intent (in) :: porParams
    PetscReal, intent (in)                  :: P
    PetscReal, intent (out)                 :: por
    PetscReal, intent (out)                 :: dpor_dP

    select case(porParams%porosityId)
    case (POROSITY_CONSTANT)
       call PorosityFunctionComputationConstantModel(porParams, por, dpor_dP)
    case (POROSITY_LINEAR)
       call PorosityFunctionComputationLinearModel(porParams, P, por, dpor_dP)
    case default
       write(iulog,*)'PorosityFunctionComputation: Unknown porosityId'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine PorosityFunctionComputation

  !------------------------------------------------------------------------
  subroutine PorosityFunctionComputationConstantModel (porParams, por, dpor_dP)
    !
    ! !DESCRIPTION:
    ! Given a pressure value, returns  porosity and derivative of prosity
    ! w.r.t pressure for constant porosity model
    !
    implicit none

    ! !ARGUMENTS
    type(porosity_params_type), intent (in) :: porParams
    PetscReal, intent (out)                 :: por
    PetscReal, intent (out)                 :: dpor_dP

    por      = porParams%porosity_base
    dpor_dP  = 0.d0

  end subroutine PorosityFunctionComputationConstantModel

  !------------------------------------------------------------------------
  subroutine PorosityFunctionComputationLinearModel (porParams, P, por, dpor_dP)
    !
    ! !DESCRIPTION:
    ! Given a pressure value, returns  porosity and derivative of prosity
    ! w.r.t pressure for linear porosity model
    !
    implicit none
    !
    ! !ARGUMENTS
    type(porosity_params_type), intent (in) :: porParams
    PetscReal, intent (in)                  :: P
    PetscReal, intent (out)                 :: por
    PetscReal, intent (out)                 :: dpor_dP

    por      = porParams%porosity_base  + &
         (P - porParams%pressure_reference)*porParams%lin_mod_slope
    dpor_dP  = porParams%lin_mod_slope

  end subroutine PorosityFunctionComputationLinearModel

  !------------------------------------------------------------------------
  subroutine PorosityFunctionCopy (this, porParams)
    !
    ! !DESCRIPTION:
    ! Makes a copy of values from another porosity parameters (porParams).
    !
    implicit none
    !
    ! !ARGUMENTS
    class(porosity_params_type), intent (out) :: this
    class(porosity_params_type), intent (in)  :: porParams

    select case(porParams%porosityId)
    case (POROSITY_CONSTANT)

       call PorosityFunctionSetConstantModel (this, porParams%porosity_base)

    case (POROSITY_LINEAR)

       call PorosityFunctionSetLinearModel (this, &
            porParams%porosity_base,              &
            porParams%pressure_reference,         &
            porParams%lin_mod_slope)

    case default
       write(iulog,*)'PorosityFunctionCopy: Unknown porosityId'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select


  end subroutine PorosityFunctionCopy

end module PorosityFunctionMod

#endif
