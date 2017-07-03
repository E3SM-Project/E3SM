module AuxVars_TOilIms_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public, extends(auxvar_flow_energy_type) :: auxvar_toil_ims_type
    ! no data at the moment
  contains
    procedure, public :: Init => AuxVarTOilImsInit
    procedure, public :: Strip => AuxVarTOilImsStrip
  end type auxvar_toil_ims_type

  public :: AuxVarTOilImsStrip

contains

! ************************************************************************** !

subroutine AuxVarTOilImsInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_toil_ims_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0
  
  !part of this initialization could be common to all modes
  ! need to agree a standard filling for the pressure:
  ! nphases and capillaty pressure afterwards? (can be 2 for three phases) 

  !the indended part below should go on AuxVarsFlowInit, however we need a
  !standard for the pressure array entry, e.g. 
  !nphases and capillaty pressure afterwards? (pcs be 2 for three phases)
  !otherwise define a pc array which is mode dependent
    !two phases (water,oil) and capillary pressure
    allocate(this%pres(option%nphase+ONE_INTEGER))
    this%pres = 0.d0
    allocate(this%sat(option%nphase))
    this%sat = 0.d0
    allocate(this%den(option%nphase))
    this%den = 0.d0
    allocate(this%den_kg(option%nphase))
    this%den_kg = 0.d0
    allocate(this%mobility(option%nphase))
    this%mobility = 0.d0


  !this could go on AuxVarsFlowEnergyInit
    this%temp = 0.d0
    allocate(this%H(option%nphase))
    this%H = 0.d0
    allocate(this%U(option%nphase))
    this%U = 0.d0
  !  

end subroutine AuxVarTOilImsInit

! ************************************************************************** !

subroutine AuxVarTOilImsStrip(this)
  ! 
  ! TOilImsAuxVarDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_toil_ims_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

end subroutine AuxVarTOilImsStrip
! ************************************************************************** !

end module AuxVars_TOilIms_module

