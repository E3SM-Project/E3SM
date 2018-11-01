module Reaction_Immobile_Aux_module
  
  use Reaction_Database_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  type, public :: immobile_species_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: molar_weight
    PetscBool :: print_me
    type(immobile_species_type), pointer :: next    
  end type immobile_species_type
  
  type, public :: immobile_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    character(len=MAXWORDLENGTH), pointer :: constraint_aux_string(:)
    PetscBool, pointer :: external_dataset(:)
  end type immobile_constraint_type
  
  type, public :: immobile_decay_rxn_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: rate_constant
    PetscReal :: half_life
    PetscBool :: print_me
    type(immobile_decay_rxn_type), pointer :: next
  end type immobile_decay_rxn_type
  
  type, public :: immobile_type

    PetscInt :: nimmobile
    PetscBool :: print_all
    
    type(immobile_species_type), pointer :: list
    type(immobile_decay_rxn_type), pointer :: decay_rxn_list

    ! immobile species
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscBool, pointer :: print_me(:)    
    
    ! decay rxn
    PetscInt :: ndecay_rxn
    PetscInt, pointer :: decayspecid(:)
    PetscReal, pointer :: decay_rate_constant(:)    

  end type immobile_type
  
  interface GetImmobileSpeciesIDFromName
    module procedure GetImmobileSpeciesIDFromName1
    module procedure GetImmobileSpeciesIDFromName2
  end interface  

  public :: ImmobileCreate, &
            ImmobileSpeciesCreate, &
            ImmobileConstraintCreate, &
            ImmobileDecayRxnCreate, &
            ImmobileGetCount, &
            ImmobileConstraintDestroy, &
            GetImmobileSpeciesIDFromName, &
            ImmobileDestroy
             
contains

! ************************************************************************** !

function ImmobileCreate()
  ! 
  ! Allocate and initialize immobile object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/11/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  type(immobile_type), pointer :: ImmobileCreate
  
  type(immobile_type), pointer :: immobile

  allocate(immobile)  
  nullify(immobile%list)
  nullify(immobile%decay_rxn_list)
  immobile%nimmobile = 0
  immobile%print_all = PETSC_FALSE
  nullify(immobile%names)
  nullify(immobile%print_me)
  
  immobile%ndecay_rxn = 0
  nullify(immobile%decayspecid)
  nullify(immobile%decay_rate_constant)

  ImmobileCreate => immobile
  
end function ImmobileCreate

! ************************************************************************** !

function ImmobileSpeciesCreate()
  ! 
  ! Allocate and initialize a immobile species object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  type(immobile_species_type), pointer :: ImmobileSpeciesCreate
  
  type(immobile_species_type), pointer :: species

  allocate(species)  
  species%id = 0
  species%name = ''
  species%molar_weight = 0.d0
  species%print_me = PETSC_FALSE
  nullify(species%next)

  ImmobileSpeciesCreate => species
  
end function ImmobileSpeciesCreate

! ************************************************************************** !

function ImmobileConstraintCreate(immobile,option)
  ! 
  ! Creates a immobile constraint object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none
  
  type(immobile_type) :: immobile
  type(option_type) :: option
  type(immobile_constraint_type), pointer :: ImmobileConstraintCreate

  type(immobile_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(immobile%nimmobile))
  constraint%names = ''
  allocate(constraint%constraint_conc(immobile%nimmobile))
  constraint%constraint_conc = 0.d0
  allocate(constraint%constraint_aux_string(immobile%nimmobile))
  constraint%constraint_aux_string = ''
  allocate(constraint%external_dataset(immobile%nimmobile))
  constraint%external_dataset = PETSC_FALSE

  ImmobileConstraintCreate => constraint

end function ImmobileConstraintCreate

! ************************************************************************** !

function ImmobileDecayRxnCreate()
  ! 
  ! Allocate and initialize a immobile decay reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/31/15
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
    
  type(immobile_decay_rxn_type), pointer :: ImmobileDecayRxnCreate

  type(immobile_decay_rxn_type), pointer :: rxn
  
  allocate(rxn)
  rxn%id = 0
  rxn%species_name = ''
  rxn%rate_constant = 0.d0
  rxn%half_life = 0.d0
  rxn%print_me = PETSC_FALSE
  nullify(rxn%next)
  
  ImmobileDecayRxnCreate => rxn
  
end function ImmobileDecayRxnCreate

! ************************************************************************** !

function ImmobileGetCount(immobile)
  ! 
  ! Returns the number of immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  implicit none
  
  PetscInt :: ImmobileGetCount
  type(immobile_type) :: immobile

  type(immobile_species_type), pointer :: immobile_species

  ImmobileGetCount = 0
  immobile_species => immobile%list
  do
    if (.not.associated(immobile_species)) exit
    ImmobileGetCount = ImmobileGetCount + 1
    immobile_species => immobile_species%next
  enddo

end function ImmobileGetCount

! ************************************************************************** !

function GetImmobileSpeciesIDFromName1(name,immobile,option)
  ! 
  ! Returns the id of named immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(immobile_type) :: immobile
  type(option_type) :: option
  
  PetscInt :: GetImmobileSpeciesIDFromName1

  GetImmobileSpeciesIDFromName1 = &
    GetImmobileSpeciesIDFromName2(name,immobile,PETSC_TRUE,option)
  
end function GetImmobileSpeciesIDFromName1

! ************************************************************************** !

function GetImmobileSpeciesIDFromName2(name,immobile,return_error,option)
  ! 
  ! Returns the id of named immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/13
  ! 

  use Option_module
  use String_module
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: name
  type(immobile_type) :: immobile
  PetscBool :: return_error
  type(option_type) :: option

  PetscInt :: GetImmobileSpeciesIDFromName2

  type(immobile_species_type), pointer :: species
  PetscInt :: i

  GetImmobileSpeciesIDFromName2 = UNINITIALIZED_INTEGER
  
  ! if the primary species name list exists
  if (associated(immobile%names)) then
    do i = 1, size(immobile%names)
      if (StringCompare(name,immobile%names(i), &
                        MAXWORDLENGTH)) then
        GetImmobileSpeciesIDFromName2 = i
        exit
      endif
    enddo
  else
    species => immobile%list
    i = 0
    do
      if (.not.associated(species)) exit
      i = i + 1
      if (StringCompare(name,species%name,MAXWORDLENGTH)) then
        GetImmobileSpeciesIDFromName2 = i
        exit
      endif
      species => species%next
    enddo
  endif

  if (return_error .and. GetImmobileSpeciesIDFromName2 <= 0) then
    option%io_buffer = 'Species "' // trim(name) // &
      '" not founds among immobile species in GetImmobileSpeciesIDFromName().'
    call printErrMsg(option)
  endif
  
end function GetImmobileSpeciesIDFromName2

! ************************************************************************** !

subroutine ImmobileSpeciesDestroy(species)
  ! 
  ! Deallocates a immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  implicit none
    
  type(immobile_species_type), pointer :: species

  if (.not.associated(species)) return
  
  deallocate(species)  
  nullify(species)

end subroutine ImmobileSpeciesDestroy

! ************************************************************************** !

recursive subroutine ImmobileDecayRxnDestroy(rxn)
  ! 
  ! Deallocates a general reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/31/15
  ! 

  implicit none
    
  type(immobile_decay_rxn_type), pointer :: rxn

  if (.not.associated(rxn)) return
  
  call ImmobileDecayRxnDestroy(rxn%next)
  nullify(rxn%next)
  deallocate(rxn)  
  nullify(rxn)

end subroutine ImmobileDecayRxnDestroy

! ************************************************************************** !

subroutine ImmobileConstraintDestroy(constraint)
  ! 
  ! Destroys a colloid constraint object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/12/10
  ! 

  use Utility_module, only: DeallocateArray

  implicit none
  
  type(immobile_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%constraint_aux_string)
  call DeallocateArray(constraint%external_dataset)
  
  deallocate(constraint)
  nullify(constraint)

end subroutine ImmobileConstraintDestroy

! ************************************************************************** !

subroutine ImmobileDestroy(immobile)
  ! 
  ! Deallocates a immobile object
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  ! 

  use Utility_module, only: DeallocateArray
  
  implicit none

  type(immobile_type), pointer :: immobile
  
  type(immobile_species_type), pointer :: cur_immobile_species, &
                                         prev_immobile_species

  if (.not.associated(immobile)) return

  ! immobile species
  cur_immobile_species => immobile%list
  do
    if (.not.associated(cur_immobile_species)) exit
    prev_immobile_species => cur_immobile_species
    cur_immobile_species => cur_immobile_species%next
    call ImmobileSpeciesDestroy(prev_immobile_species)
  enddo    
  nullify(immobile%list)
  
  call DeallocateArray(immobile%names)
  call DeallocateArray(immobile%print_me)
  
  call ImmobileDecayRxnDestroy(immobile%decay_rxn_list)
  call DeallocateArray(immobile%decayspecid)
  call DeallocateArray(immobile%decay_rate_constant)
  
  deallocate(immobile)
  nullify(immobile)

end subroutine ImmobileDestroy

end module Reaction_Immobile_Aux_module
