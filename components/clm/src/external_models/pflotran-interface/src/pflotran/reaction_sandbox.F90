module Reaction_Sandbox_module

  use Reaction_Sandbox_Base_class
  
  ! Add new reacton sandbox classes here.
  use Reaction_Sandbox_SomDec_class
  use Reaction_Sandbox_PlantN_class
  use Reaction_Sandbox_Langmuir_class
  use Reaction_Sandbox_Microbial_class
  use Reaction_Sandbox_Nitrif_class
  use Reaction_Sandbox_Denitr_class
  use Reaction_Sandbox_Degas_class
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  class(reaction_sandbox_base_type), pointer, public :: rxn_sandbox_list

  interface RSandboxRead
    module procedure RSandboxRead1
    module procedure RSandboxRead2
  end interface
  
  interface RSandboxDestroy
    module procedure RSandboxDestroy1
    module procedure RSandboxDestroy2
  end interface
  
  public :: RSandboxInit, &
            RSandboxRead, &
            RSandboxSkipInput, &
            RSandboxSetup, &
            RSandbox, &
            RSandboxUpdateKineticState, &
            RSandboxAuxiliaryPlotVariables, &
            RSandboxDestroy

contains

! ************************************************************************** !

subroutine RSandboxInit(option)
  ! 
  ! Initializes the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  ! 
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(rxn_sandbox_list)) then
    call RSandboxDestroy()
  endif
  nullify(rxn_sandbox_list)

end subroutine RSandboxInit

! ************************************************************************** !

subroutine RSandboxSetup(reaction,option)
  ! 
  ! Calls all the initialization routines for all reactions in
  ! the sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  ! 

  use Option_module
  use Reaction_Aux_module, only : reaction_type 
  
  implicit none
  
  type(reaction_type) :: reaction
  type(option_type) :: option
  
  class(reaction_sandbox_base_type), pointer :: cur_sandbox  

  ! sandbox reactions
  cur_sandbox => rxn_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Setup(reaction,option)
    cur_sandbox => cur_sandbox%next
  enddo 

end subroutine RSandboxSetup

! ************************************************************************** !

subroutine RSandboxRead1(input,option)
  ! 
  ! Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/16/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option

  call RSandboxRead(rxn_sandbox_list,input,option)

end subroutine RSandboxRead1

! ************************************************************************** !

subroutine RSandboxRead2(local_sandbox_list,input,option)
  ! 
  ! RSandboxRead: Reads input deck for reaction sandbox parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  class(reaction_sandbox_base_type), pointer :: local_sandbox_list  
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  class(reaction_sandbox_base_type), pointer :: new_sandbox, cur_sandbox
  
  nullify(new_sandbox)
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,REACTION_SANDBOX')
    call StringToUpper(word)   

    select case(trim(word))
      ! Add new cases statements for new reacton sandbox classes here.
      case('SOMDECOMP')
        new_sandbox => SomDecCreate()
      case('PLANTN')
        new_sandbox => PlantNCreate()
      case('CLM-MICROBIAL')
        new_sandbox => MicrobialCreate()
      case('NITRIFICATION')
        new_sandbox => NitrifCreate()
      case('DENITRIFICATION')
        new_sandbox => DenitrCreate()
      case('DEGAS')
        new_sandbox => degasCreate()
      case('LANGMUIR')
        new_sandbox => LangmuirCreate()
      case default
        call InputKeywordUnrecognized(word,'CHEMISTRY,REACTION_SANDBOX',option)
    end select
    
    call new_sandbox%ReadInput(input,option)
    
    if (.not.associated(local_sandbox_list)) then
      local_sandbox_list => new_sandbox
    else
      cur_sandbox => local_sandbox_list
      do
        if (.not.associated(cur_sandbox%next)) exit
        cur_sandbox => cur_sandbox%next
      enddo
      cur_sandbox%next => new_sandbox
    endif
  enddo
  
end subroutine RSandboxRead2

! ************************************************************************** !

subroutine RSandboxAuxiliaryPlotVariables(list,reaction,option)
  ! 
  ! Adds auxilairy plot variables to the list to be printed
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/16
  ! 

  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  
  implicit none

  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  type(reaction_type) :: reaction
  
  class(reaction_sandbox_base_type), pointer :: cur_reaction
  
  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
    call cur_reaction%AuxiliaryPlotVariables(list,reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandboxAuxiliaryPlotVariables

! ************************************************************************** !

subroutine RSandboxSkipInput(input,option)
  ! 
  ! Intelligently skips over REACTION_SANDBOX block
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  
  class(reaction_sandbox_base_type), pointer :: dummy_list
  
  nullify(dummy_list)
  call RSandboxRead(dummy_list,input,option)
  call RSandboxDestroy(dummy_list)
  
end subroutine RSandboxSkipInput

! ************************************************************************** !

subroutine RSandbox(Residual,Jacobian,compute_derivative,rt_auxvar, &
                    global_auxvar,material_auxvar,reaction,option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class, only: material_auxvar_type
  
  implicit none

  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  class(reaction_sandbox_base_type), pointer :: cur_reaction
  
  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
      call cur_reaction%Evaluate(Residual,Jacobian,compute_derivative, &
                                 rt_auxvar,global_auxvar,material_auxvar, &
                                 reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandbox

! ************************************************************************** !

subroutine RSandboxUpdateKineticState(rt_auxvar,global_auxvar, &
                                      material_auxvar,reaction,option)
  ! 
  ! Updates volume fractions, etc. at the end of a time step.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/06/16
  ! 

  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class, only: material_auxvar_type
  
  implicit none

  type(option_type) :: option
  type(reaction_type) :: reaction
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  class(reaction_sandbox_base_type), pointer :: cur_reaction
  
  cur_reaction => rxn_sandbox_list
  do
    if (.not.associated(cur_reaction)) exit
      call cur_reaction%UpdateKineticState(rt_auxvar,global_auxvar, &
                                           material_auxvar,reaction,option)
    cur_reaction => cur_reaction%next
  enddo

end subroutine RSandboxUpdateKineticState

! ************************************************************************** !

subroutine RSandboxDestroy1()
  ! 
  ! Destroys master sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/16/13
  ! 

  implicit none

  call RSandboxDestroy(rxn_sandbox_list)
  
end subroutine RSandboxDestroy1

! ************************************************************************** !

subroutine RSandboxDestroy2(local_sandbox_list)
  ! 
  ! Destroys arbitrary sandbox list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/08/12
  ! 

  implicit none

  class(reaction_sandbox_base_type), pointer :: local_sandbox_list

  class(reaction_sandbox_base_type), pointer :: cur_sandbox, prev_sandbox
  
  ! sandbox reactions
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    prev_sandbox => cur_sandbox%next
    call cur_sandbox%Destroy()
    deallocate(cur_sandbox)
    cur_sandbox => prev_sandbox
  enddo  

end subroutine RSandboxDestroy2

end module Reaction_Sandbox_module
