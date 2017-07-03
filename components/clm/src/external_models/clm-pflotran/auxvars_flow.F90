module AuxVars_Flow_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_base_type) :: auxvar_flow_type
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal, pointer :: mobility(:) ! relative perm / dynamic viscosity
    !PetscReal, pointer :: dsat_dp(:,:)
    !PetscReal, pointer :: dden_dp(:,:)
    !PetscReal, pointer :: dsat_dt(:)
    !PetscReal, pointer :: dden_dt(:)
    !PetscReal, pointer :: dmobility_dp(:)
  contains
    !procedure, public :: Init => InitAuxVarFlow
  end type auxvar_flow_type

  public :: AuxVarFlowStrip

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
 
  !this could be used to initialize flow part of auxvar

  !currently does nothing - could init the base members
  !print *, 'Must extend InitAuxVarFlow '
  !stop    

end subroutine AuxVarFlowInit

! ************************************************************************** !

subroutine AuxVarFlowStrip(this)
  ! 
  ! AuxVarFlowDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 8/5/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_type) :: this

  call DeallocateArray(this%pres)  
  call DeallocateArray(this%sat)  
  call DeallocateArray(this%den)  
  call DeallocateArray(this%den_kg)  
  call DeallocateArray(this%mobility)  

end subroutine AuxVarFlowStrip
! ************************************************************************** !

end module AuxVars_Flow_module

