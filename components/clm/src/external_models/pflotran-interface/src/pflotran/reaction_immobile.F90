module Reaction_Immobile_module

  use Reaction_Immobile_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  public :: ImmobileRead, &
            ImmobileDecayRxnRead, &
            ImmobileProcessConstraint, &
            RImmobileDecay

contains

! ************************************************************************** !

subroutine ImmobileRead(immobile,input,option)
  ! 
  ! Reads immobile species
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/02/13
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(immobile_type) :: immobile
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(immobile_species_type), pointer :: new_immobile_species, &
                                          prev_immobile_species
           
  ! find end of list if it exists
  if (associated(immobile%list)) then
    new_immobile_species => immobile%list
    do
      if (.not.associated(new_immobile_species%next)) exit
      new_immobile_species => new_immobile_species%next
    enddo
    prev_immobile_species => new_immobile_species
    nullify(new_immobile_species)
  else
    nullify(prev_immobile_species)
  endif
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    ! this count is required for comparisons prior to BasisInit()
    immobile%nimmobile = immobile%nimmobile + 1          
    new_immobile_species => ImmobileSpeciesCreate()
    call InputReadWord(input,option,new_immobile_species%name,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                        'CHEMISTRY,IMMOBILE_SPECIES')
    if (.not.associated(prev_immobile_species)) then
      immobile%list => new_immobile_species
      new_immobile_species%id = 1
    else
      prev_immobile_species%next => new_immobile_species
      new_immobile_species%id = prev_immobile_species%id + 1
    endif
    prev_immobile_species => new_immobile_species
    nullify(new_immobile_species)
  enddo   

end subroutine ImmobileRead

! ************************************************************************** !

subroutine ImmobileDecayRxnRead(immobile,input,option)
  ! 
  ! Reads chemical species
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Units_module
  
  implicit none
  
  type(immobile_type) :: immobile
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string
  type(immobile_decay_rxn_type), pointer :: immobile_decay_rxn
  type(immobile_decay_rxn_type), pointer :: cur_immobile_decay_rxn

  error_string = 'CHEMISTRY,IMMOBILE_DECAY_REACTION'
  
  immobile%ndecay_rxn = immobile%ndecay_rxn + 1
        
  immobile_decay_rxn => ImmobileDecayRxnCreate()
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)   

    select case(trim(word))
      case('SPECIES_NAME')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'species name',error_string)
        immobile_decay_rxn%species_name = word
      case('RATE_CONSTANT')
        internal_units = '1/sec'
        call InputReadDouble(input,option,immobile_decay_rxn%rate_constant)  
        call InputErrorMsg(input,option,'rate cosntant', &
                             'CHEMISTRY,IMMOBILE_DECAY_REACTION') 
        call InputReadAndConvertUnits(input,immobile_decay_rxn%rate_constant, &
                                      internal_units, &
                                      trim(error_string)//',rate constant', &
                                      option)
      case('HALF_LIFE')
        internal_units = 'sec'
        call InputReadDouble(input,option,immobile_decay_rxn%half_life)
        call InputErrorMsg(input,option,'half life',error_string)
        call InputReadAndConvertUnits(input,immobile_decay_rxn%half_life, &
                                      internal_units, &
                                      trim(error_string)//',half life',option)
        ! convert half life to rate constant
        immobile_decay_rxn%rate_constant = &
          -1.d0*log(0.5d0)/immobile_decay_rxn%half_life
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  if (Uninitialized(immobile_decay_rxn%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT or HALF_LIFE must be set in ' // &
      'IMMOBILE_DECAY_REACTION.'
    call printErrMsg(option)
  endif
  if (.not.associated(immobile%decay_rxn_list)) then
    immobile%decay_rxn_list => immobile_decay_rxn
    immobile_decay_rxn%id = 1
  else
    cur_immobile_decay_rxn => immobile%decay_rxn_list
    do
      if (.not.associated(cur_immobile_decay_rxn%next)) then
        cur_immobile_decay_rxn%next => immobile_decay_rxn
        immobile_decay_rxn%id = cur_immobile_decay_rxn%id + 1
        exit
      endif
      cur_immobile_decay_rxn => cur_immobile_decay_rxn%next
    enddo
  endif
  nullify(immobile_decay_rxn)

end subroutine ImmobileDecayRxnRead

! ************************************************************************** !

subroutine ImmobileProcessConstraint(immobile,constraint_name, &
                                    constraint,option)
  ! 
  ! Initializes constraints based on immobile
  ! species in system
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module  
  
  implicit none
  
  type(immobile_type), pointer :: immobile
  character(len=MAXWORDLENGTH) :: constraint_name
  type(immobile_constraint_type), pointer :: constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: iimmobile, jimmobile
  
  character(len=MAXWORDLENGTH) :: immobile_name(immobile%nimmobile)
  character(len=MAXWORDLENGTH) :: constraint_aux_string(immobile%nimmobile)
  PetscReal :: constraint_conc(immobile%nimmobile)
  PetscBool :: external_dataset(immobile%nimmobile)
  
  if (.not.associated(constraint)) return
  
  immobile_name = ''
  constraint_aux_string = ''
  external_dataset = PETSC_FALSE
  do iimmobile = 1, immobile%nimmobile
    found = PETSC_FALSE
    do jimmobile = 1, immobile%nimmobile
      if (StringCompare(constraint%names(iimmobile), &
                        immobile%names(jimmobile), &
                        MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
                'Immobile species "' // trim(constraint%names(iimmobile)) // &
                '" from CONSTRAINT "' // trim(constraint_name) // &
                '" not found among immobile species.'
      call printErrMsg(option)
    else
      immobile_name(iimmobile) = constraint%names(iimmobile)
      constraint_conc(iimmobile) = &
        constraint%constraint_conc(iimmobile)
      constraint_aux_string(iimmobile) = &
        constraint%constraint_aux_string(iimmobile)
      external_dataset(iimmobile) = constraint%external_dataset(iimmobile)
    endif  
  enddo
  constraint%names = immobile_name
  constraint%constraint_conc = constraint_conc
  constraint%constraint_aux_string = constraint_aux_string
  constraint%external_dataset = external_dataset

end subroutine ImmobileProcessConstraint

! ************************************************************************** !

subroutine RImmobileDecay(Res,Jac,compute_derivative,rt_auxvar, &
                          global_auxvar,material_auxvar,reaction, &
                          option)
  ! 
  ! Computes decay of biomass species 
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/31/15
  ! 
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Res(reaction%ncomp)
  PetscReal :: Jac(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  
  PetscInt :: icomp, irxn, immobile_id
  PetscReal :: rate_constant, rate, volume

  PetscInt, parameter :: iphase = 1

  volume = material_auxvar%volume

  do irxn = 1, reaction%immobile%ndecay_rxn ! for each reaction
    
    ! we assume only one chemical component involved in decay reaction
    icomp = reaction%immobile%decayspecid(irxn)
    ! units = m^3 bulk/sec = [1/sec] * [m^3 bulk]
    rate_constant = reaction%immobile%decay_rate_constant(irxn)*volume
    ! rate [mol/sec] = [m^3 bulk/sec] * [mol/m^3 bulk]
    rate = rate_constant*rt_auxvar%immobile(icomp)
    immobile_id = reaction%offset_immobile + icomp
    
    ! units = mol/sec              ! implicit stoichiometry of -1.d0 (- -1.d0*)
    Res(immobile_id) = Res(immobile_id) + rate

    if (.not. compute_derivative) cycle
    ! units = (mol/sec)*(m^3/mol) = m^3/sec
    Jac(immobile_id,immobile_id) = Jac(immobile_id,immobile_id) + rate_constant
    
  enddo  ! loop over reactions
    
end subroutine RImmobileDecay

end module Reaction_Immobile_module
