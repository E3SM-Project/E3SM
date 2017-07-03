module Reaction_Solid_Soln_Aux_module
  
  use Reaction_Mineral_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public :: solid_solution_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: num_stoich_solid
    character(len=MAXWORDLENGTH), pointer :: stoich_solid_names(:)
    PetscInt, pointer :: stoich_solid_ids(:)
#if 0
    PetscInt :: num_end_member
    type(stoichiometric_solid_type), pointer :: stoich_solid
#endif    
    type(solid_solution_type), pointer :: next
  end type solid_solution_type

#if 0
  type, public :: stoichiometric_solid_type
    type(mineral_rxn_type), pointer :: mineral ! stoichiometric solid
    type(mineral_rxn_type), pointer :: end_members
    type(stoichiometric_solid_type), pointer :: next
  end type stoichiometric_solid_type
    
  type, public :: solid_solution_rxn_type
    character(len=MAXSTRINGLENGTH) :: database_filename
    PetscInt :: num_dbase_temperatures
    PetscReal, pointer :: dbase_temperatures(:)
    type(solid_solution_type), pointer :: list
    type(mineral_type), pointer :: mineral
  end type solid_solution_rxn_type
#endif
  
  public :: SolidSolutionCreate, &
            SolidSolutionDestroy
             
contains

#if 0

! ************************************************************************** !

function SolidSolutionReactionCreate()
  ! 
  ! Allocate and initialize solid solution reaction
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  ! 

  implicit none
  
  type(solid_solution_rxn_type), pointer :: SolidSolutionReactionCreate
  
  type(solid_solution_rxn_type), pointer :: solid_solution_rxn

  allocate(solid_solution_rxn)
  
  solid_solution_rxn%num_dbase_temperatures = 0
  nullify(solid_solution_rxn%dbase_temperatures)

  nullify(solid_solution_rxn%list)
  
  solid_solution_rxn%mineral => MineralReactionCreate()

  SolidSolutionReactionCreate => solid_solution_rxn
  
end function SolidSolutionReactionCreate
#endif

! ************************************************************************** !

function SolidSolutionCreate()
  ! 
  ! Allocate and initialize solid solution object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/17/12
  ! 

  implicit none
  
  type(solid_solution_type), pointer :: SolidSolutionCreate
  
  type(solid_solution_type), pointer :: solid_solution

  allocate(solid_solution)
  
  solid_solution%name = ''
  solid_solution%num_stoich_solid = 0
#if 0  
  solid_solution%num_end_member = 0
#endif

  nullify(solid_solution%stoich_solid_names)
  nullify(solid_solution%stoich_solid_ids)
  nullify(solid_solution%next)
  
  SolidSolutionCreate => solid_solution
  
end function SolidSolutionCreate

#if 0

! ************************************************************************** !

function StoichiometricSolidCreate()
  ! 
  ! Allocate and initialize stoichiometric solid
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/17/12
  ! 

  implicit none
  
  type(stoichiometric_solid_type), pointer :: StoichiometricSolidCreate
  
  type(stoichiometric_solid_type), pointer :: stoich_solid

  allocate(stoich_solid)
  
  nullify(stoich_solid%mineral)
  nullify(stoich_solid%end_members) ! nullify the list for now
  nullify(stoich_solid%next)
  
  StoichiometricSolidCreate => stoich_solid
  
end function StoichiometricSolidCreate

! ************************************************************************** !

subroutine StoichiometricSolidDestroy(stoich_solid)
  ! 
  ! Deallocates solid solution object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/17/12
  ! 

  implicit none
  
  type(stoichiometric_solid_type), pointer :: stoich_solid
  
  type(mineral_rxn_type), pointer :: cur_mineral, prev_mineral

  if (.not.associated(stoich_solid)) return
  
  ! Do not use recursion here as it may result in a large drain on the stack

  cur_mineral => stoich_solid%end_members
  do
    if (.not.associated(cur_mineral)) exit
    prev_mineral => cur_mineral
    cur_mineral => cur_mineral%next
    call MineralDestroy(prev_mineral)
  enddo    
  call MineralDestroy(stoich_solid%mineral)
  
  deallocate(stoich_solid)
  nullify(stoich_solid)
  
end subroutine StoichiometricSolidDestroy
#endif

! ************************************************************************** !

recursive subroutine SolidSolutionDestroy(solid_solution)
  ! 
  ! Deallocates solid solution object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/17/12
  ! 

  implicit none
  
  type(solid_solution_type), pointer :: solid_solution
  
#if 0  
  type(stoichiometric_solid_type), pointer :: cur_stoich_solid, &
                                              prev_stoich_solid
#endif

  if (.not.associated(solid_solution)) return
  
  ! recursive
  call SolidSolutionDestroy(solid_solution%next)

#if 0  
  ! I don't want to destroy recursively here as the memory use may
  ! be to large for large solid solutions
  cur_stoich_solid => solid_solution%stoich_solid
  do 
    if (.not.associated(cur_stoich_solid)) exit
    prev_stoich_solid => cur_stoich_solid
    cur_stoich_solid => cur_stoich_solid%next
    call StoichiometricSolidDestroy(prev_stoich_solid)
  enddo
#endif
  deallocate(solid_solution%stoich_solid_names)
  nullify(solid_solution%stoich_solid_names)
  deallocate(solid_solution%stoich_solid_ids)
  nullify(solid_solution%stoich_solid_ids)
  
  deallocate(solid_solution)
  nullify(solid_solution)
  
end subroutine SolidSolutionDestroy

#if 0

! ************************************************************************** !

subroutine SolidSolutionReactionDestroy(solid_solution)
  ! 
  ! Deallocates a solid solution object
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  ! 

  implicit none

  type(solid_solution_rxn_type), pointer :: solid_solution
  
  ! recursive
  call SolidSolutionDestroy(solid_solution%list)
  
  deallocate(solid_solution)
  nullify(solid_solution)

end subroutine SolidSolutionReactionDestroy
#endif

end module Reaction_Solid_Soln_Aux_module
