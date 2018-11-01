module AuxVars_Flow_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module

  implicit none

  private

  type, public, extends(auxvar_base_type) :: auxvar_flow_type
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: pc(:)     ! capillary pressure (iphase-1)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal, pointer :: mobility(:) ! relative perm / dynamic viscosity
    PetscReal, pointer :: viscosity(:) ! dynamic viscosity
    PetscInt, pointer :: table_idx(:)
    !PetscReal, pointer :: dsat_dp(:,:)
    !PetscReal, pointer :: dden_dp(:,:)
    !PetscReal, pointer :: dsat_dt(:)
    !PetscReal, pointer :: dden_dt(:)
    !PetscReal, pointer :: dmobility_dp(:)
  contains
    !procedure, public :: Init => InitAuxVarFlow
  end type auxvar_flow_type

  public :: AuxVarFlowInit, AuxVarFlowStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowInit(this,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  !

  use Option_module

  implicit none

  class(auxvar_flow_type) :: this
  type(option_type) :: option

  allocate(this%pres(option%nphase))
  this%pres = 0.d0
  allocate(this%pc(option%nphase - ONE_INTEGER))
  this%pc = 0.0d0
  allocate(this%sat(option%nphase))
  this%sat = 0.d0
  allocate(this%den(option%nphase))
  this%den = 0.d0
  allocate(this%den_kg(option%nphase))
  this%den_kg = 0.d0
  allocate(this%mobility(option%nphase))
  this%mobility = 0.d0
  allocate(this%viscosity(option%nphase))
  this%viscosity = 0.d0
  if (option%neos_table_indices > 0) then
    allocate(this%table_idx(option%neos_table_indices))
    this%table_idx = 1
  end if


end subroutine AuxVarFlowInit

! ************************************************************************** !

subroutine AuxVarFlowStrip(this)
  !
  ! AuxVarFlowDestroy: Deallocates a flow auxiliary object
  !
  ! Author: Paolo Orsini
  ! Date: 8/5/16
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_type) :: this

  call DeallocateArray(this%pres)
  call DeallocateArray(this%pc)
  call DeallocateArray(this%sat)
  call DeallocateArray(this%den)
  call DeallocateArray(this%den_kg)
  call DeallocateArray(this%mobility)
  call DeallocateArray(this%viscosity)
  call DeallocateArray(this%table_idx)

end subroutine AuxVarFlowStrip
! ************************************************************************** !

end module AuxVars_Flow_module
