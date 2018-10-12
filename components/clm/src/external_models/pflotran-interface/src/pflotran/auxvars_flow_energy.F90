module AuxVars_FlowEnergy_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_flow_type) :: auxvar_flow_energy_type
    PetscReal :: temp
    PetscReal, pointer :: H(:) ! MJ/kmol
    PetscReal, pointer :: U(:) ! MJ/kmol
  contains
   !..............
  end type auxvar_flow_energy_type

  public :: AuxVarFlowEnergyStrip

contains

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

