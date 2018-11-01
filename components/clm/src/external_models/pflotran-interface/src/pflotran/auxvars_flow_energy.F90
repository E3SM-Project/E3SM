module AuxVars_FlowEnergy_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module

  implicit none
  
  private 

  type, public, extends(auxvar_flow_type) :: auxvar_flow_energy_type
    PetscReal :: temp
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
  contains
   !..............
  end type auxvar_flow_energy_type

  public :: AuxVarFlowEnergyInit, AuxVarFlowEnergyStrip

contains

! ************************************************************************** !

subroutine AuxVarFlowEnergyInit(this,option)
  ! 
  ! Initialise energy auxiliary variables
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/16
  ! 

  use Option_module

  implicit none

  class(auxvar_flow_energy_type) :: this
  type(option_type) :: option

  this%temp = 0.d0
  allocate(this%H(option%nphase))
  this%H = 0.d0
  allocate(this%U(option%nphase))
  this%U = 0.d0

end subroutine AuxVarFlowEnergyInit

! ************************************************************************** !

subroutine AuxVarFlowEnergyStrip(this)
  ! 
  ! AuxVarFlowDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 8/5/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_energy_type) :: this

  call DeallocateArray(this%H)  
  call DeallocateArray(this%U)  

end subroutine AuxVarFlowEnergyStrip

! ************************************************************************** !

end module AuxVars_FlowEnergy_module

