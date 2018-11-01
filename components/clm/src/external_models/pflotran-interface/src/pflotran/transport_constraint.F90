module Transport_Constraint_module
 
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  
  use Reaction_Mineral_Aux_module
  use Reaction_Immobile_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  ! concentration subcondition types
  PetscInt, parameter, public :: CONSTRAINT_NULL = 0
  PetscInt, parameter, public :: CONSTRAINT_FREE = 1
  PetscInt, parameter, public :: CONSTRAINT_TOTAL = 2
  PetscInt, parameter, public :: CONSTRAINT_LOG = 3
  PetscInt, parameter, public :: CONSTRAINT_PH = 4
  PetscInt, parameter, public :: CONSTRAINT_MINERAL = 5
  PetscInt, parameter, public :: CONSTRAINT_GAS = 6
  PetscInt, parameter, public :: CONSTRAINT_CHARGE_BAL = 7
  PetscInt, parameter, public :: CONSTRAINT_TOTAL_SORB = 9
  PetscInt, parameter, public :: CONSTRAINT_SUPERCRIT_CO2 = 10

  type, public :: tran_constraint_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name         
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(guess_constraint_type), pointer :: free_ion_guess
    type(mineral_constraint_type), pointer :: minerals

    type(colloid_constraint_type), pointer :: colloids
    type(immobile_constraint_type), pointer :: immobile_species
    PetscBool :: requires_equilibration
    type(tran_constraint_type), pointer :: next    
  end type tran_constraint_type
  
  type, public :: tran_constraint_ptr_type
    type(tran_constraint_type), pointer :: ptr
  end type tran_constraint_ptr_type
  
  type, public :: tran_constraint_list_type
    PetscInt :: num_constraints
    type(tran_constraint_type), pointer :: first
    type(tran_constraint_type), pointer :: last
    type(tran_constraint_ptr_type), pointer :: array(:)    
  end type tran_constraint_list_type
  
  type, public :: tran_constraint_coupler_type
    character(len=MAXWORDLENGTH) :: constraint_name   
    PetscReal :: time
    PetscInt :: num_iterations
    character(len=MAXWORDLENGTH) :: time_units
    type(aq_species_constraint_type), pointer :: aqueous_species
    type(guess_constraint_type), pointer :: free_ion_guess
    type(mineral_constraint_type), pointer :: minerals

    type(colloid_constraint_type), pointer :: colloids
    type(immobile_constraint_type), pointer :: immobile_species
    type(global_auxvar_type), pointer :: global_auxvar
    type(reactive_transport_auxvar_type), pointer :: rt_auxvar
    type(tran_constraint_coupler_type), pointer :: next   
  end type tran_constraint_coupler_type
      
  public :: TranConstraintAddToList, &
            TranConstraintInitList, &
            TranConstraintDestroyList, &
            TranConstraintGetPtrFromList, &
            TranConstraintCreate, &
            TranConstraintRead, &
            TranConstraintDestroy, &
            TranConstraintCouplerCreate, &
            TranConstraintCouplerDestroy
    
contains

! ************************************************************************** !

function TranConstraintCreate(option)
  ! 
  ! Creates a transport constraint (set of concentrations
  ! and constraints for setting boundary or initial
  ! condition).
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_constraint_type), pointer :: TranConstraintCreate
  
  type(tran_constraint_type), pointer :: constraint
  
  allocate(constraint)
  nullify(constraint%aqueous_species)
  nullify(constraint%free_ion_guess)
  nullify(constraint%minerals)

  nullify(constraint%colloids)
  nullify(constraint%immobile_species)
  nullify(constraint%next)
  constraint%id = 0
  constraint%name = ''
  constraint%requires_equilibration = PETSC_FALSE
  
  TranConstraintCreate => constraint

end function TranConstraintCreate

! ************************************************************************** !

function TranConstraintCouplerCreate(option)
  ! 
  ! Creates a coupler that ties a constraint to a
  ! transport condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_constraint_coupler_type), pointer :: TranConstraintCouplerCreate
  
  type(tran_constraint_coupler_type), pointer :: coupler
  
  allocate(coupler)
  nullify(coupler%aqueous_species)
  nullify(coupler%free_ion_guess)
  nullify(coupler%minerals)

  nullify(coupler%colloids)
  nullify(coupler%immobile_species)
  
  coupler%num_iterations = 0
  nullify(coupler%rt_auxvar)
  nullify(coupler%global_auxvar)
  
  nullify(coupler%next)
  coupler%constraint_name = ''
  coupler%time = 0.d0
  coupler%time_units = ''
  
  TranConstraintCouplerCreate => coupler

end function TranConstraintCouplerCreate

! ************************************************************************** !

subroutine TranConstraintRead(constraint,reaction,input,option)
  ! 
  ! Reads a transport constraint from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Option_module
  use Input_Aux_module
  use Units_module
  use String_module
  use Logging_module

  implicit none
  
  type(tran_constraint_type) :: constraint
  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: block_string
  PetscInt :: icomp, imnrl, iimmobile

  PetscInt :: length
  type(aq_species_constraint_type), pointer :: aq_species_constraint
  type(guess_constraint_type), pointer :: free_ion_guess_constraint
  type(mineral_constraint_type), pointer :: mineral_constraint

  type(colloid_constraint_type), pointer :: colloid_constraint
  type(immobile_constraint_type), pointer :: immobile_constraint
  PetscErrorCode :: ierr
  PetscReal :: tempreal

  call PetscLogEventBegin(logging%event_tran_constraint_read, &
                          ierr);CHKERRQ(ierr)

  ! read the constraint
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONSTRAINT')
        
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONSTRAINT')   
      
    select case(trim(word))

      case('CONC','CONCENTRATIONS')

        aq_species_constraint => &
          AqueousSpeciesConstraintCreate(reaction,option)

        block_string = 'CONSTRAINT, CONCENTRATIONS'
        icomp = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit  
          
          icomp = icomp + 1        
          
          if (icomp > reaction%naqcomp) then
            option%io_buffer = 'Number of concentration constraints ' // &
                               'exceeds number of primary chemical ' // &
                               'components in constraint: ' // &
                                trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,aq_species_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'aqueous species name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(aq_species_constraint%names(icomp))
          call printMsg(option)
          
          call InputReadDouble(input,option, &
                               aq_species_constraint%constraint_conc(icomp))
          call InputErrorMsg(input,option,'concentration',block_string)
          
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputDefaultMsg(input,option, &
                               trim(block_string) // ' constraint_type')
          length = len_trim(word)
          if (length > 0) then
            call StringToUpper(word)
            select case(word)
              case('F','FREE')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_FREE
              case('T','TOTAL')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
              case('TOTAL_SORB')
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_TOTAL_SORB
              case('S')
                option%io_buffer = '"S" constraint type no longer ' // &
                  'supported as of March 4, 2013.'
                call printErrMsg(option)
              case('P','PH')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_PH
              case('L','LOG')
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_LOG
              case('M','MINERAL','MNRL') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_MINERAL
              case('G','GAS') 
                aq_species_constraint%constraint_type(icomp) = CONSTRAINT_GAS
              case('SC','CONSTRAINT_SUPERCRIT_CO2') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_SUPERCRIT_CO2
              case('Z','CHG') 
                aq_species_constraint%constraint_type(icomp) = &
                  CONSTRAINT_CHARGE_BAL
              case default
                call InputKeywordUnrecognized(word, &
                       'CONSTRAINT,CONCENTRATION,TYPE',option)
            end select 
            
            if (aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_MINERAL .or. &
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_GAS .or.&
                aq_species_constraint%constraint_type(icomp) == &
                  CONSTRAINT_SUPERCRIT_CO2) then
              call InputReadWord(input,option,aq_species_constraint% &
                                 constraint_aux_string(icomp), &
                                 PETSC_TRUE)
              call InputErrorMsg(input,option,'constraining species name', &
                                 block_string)
            else
              call InputReadWord(input,option,word,PETSC_FALSE)
              if (input%ierr == 0) then
                call StringToUpper(word)
                select case(word)
                  case('DATASET')
                    call InputReadWord(input,option,aq_species_constraint% &
                                       constraint_aux_string(icomp),PETSC_TRUE)
                    call InputErrorMsg(input,option,'dataset name', &
                                       block_string)
                    aq_species_constraint%external_dataset(icomp) = PETSC_TRUE
                end select
              endif
            endif
          else
            aq_species_constraint%constraint_type(icomp) = CONSTRAINT_TOTAL
          endif  
        
        enddo  
        
        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is less than ' // &
                   'number of primary species in aqueous constraint.'
          call printErrMsg(option)        
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of concentration constraints is greater than ' // &
                   'number of primary species in aqueous constraint.'
          call printWrnMsg(option)        
        endif
        
        if (associated(constraint%aqueous_species)) &
          call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
        constraint%aqueous_species => aq_species_constraint 
        
      case('FREE_ION_GUESS')

        free_ion_guess_constraint => GuessConstraintCreate(reaction,option)

        block_string = 'CONSTRAINT, FREE_ION_GUESS'
        icomp = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit  
          
          icomp = icomp + 1        
          
          if (icomp > reaction%naqcomp) then
            option%io_buffer = 'Number of free ion guess constraints ' // &
                               'exceeds number of primary chemical ' // &
                               'components in constraint: ' // &
                                trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option, &
                             free_ion_guess_constraint%names(icomp), &
                             PETSC_TRUE)
          call InputErrorMsg(input,option,'free ion guess name',block_string)
          option%io_buffer = 'Constraint Species: ' // &
                             trim(free_ion_guess_constraint%names(icomp))
          call printMsg(option)
          
          call InputReadDouble(input,option,free_ion_guess_constraint%conc(icomp))
          call InputErrorMsg(input,option,'free ion guess',block_string)
        enddo

        if (icomp < reaction%naqcomp) then
          option%io_buffer = &
                   'Number of free ion guess constraints is less than ' // &
                   'number of primary species in aqueous constraint.'
          call printErrMsg(option)        
        endif
        if (icomp > reaction%naqcomp) then
          option%io_buffer = &
                   'Number of free ion guess constraints is greater than ' // &
                   'number of primary species in aqueous constraint.'
          call printWrnMsg(option)        
        endif
        
        if (associated(constraint%free_ion_guess)) &
          call GuessConstraintDestroy(constraint%free_ion_guess)
        constraint%free_ion_guess => free_ion_guess_constraint
        nullify(free_ion_guess_constraint)

      case('MNRL','MINERALS')

        mineral_constraint => MineralConstraintCreate(reaction%mineral,option)

        block_string = 'CONSTRAINT, MINERALS'
        imnrl = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit          
          
          imnrl = imnrl + 1

          if (imnrl > reaction%mineral%nkinmnrl) then
            option%io_buffer = &
                     'Number of mineral constraints exceeds number of ' // &
                     'kinetic minerals in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,mineral_constraint%names(imnrl), &
                             PETSC_TRUE)
          call InputErrorMsg(input,option,'mineral name',block_string)
          option%io_buffer = 'Constraint Minerals: ' // &
                             trim(mineral_constraint%names(imnrl))
          call printMsg(option)

          ! volume fraction
          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            input%buf = trim(string)
            call InputReadWord(input,option,mineral_constraint% &
                                constraint_vol_frac_string(imnrl),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' VOL FRAC')
            mineral_constraint%external_voL_frac_dataset(imnrl) = PETSC_TRUE
          else
            call InputReadDouble(input,option, &
                                 mineral_constraint%constraint_vol_frac(imnrl))
            call InputErrorMsg(input,option,'volume fraction',block_string)
          endif

          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            input%buf = trim(string)
            call InputReadWord(input,option,mineral_constraint% &
                                constraint_area_string(imnrl),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' SURF AREA')
            mineral_constraint%external_area_dataset(imnrl) = PETSC_TRUE
          else
            ! specific surface area
            call InputReadDouble(input,option, &
                                 mineral_constraint%constraint_area(imnrl))
            call InputErrorMsg(input,option,'surface area',block_string)
          endif
          ! read units if they exist. conversion takes place later
          call InputReadWord(input,option,mineral_constraint% &
                             constraint_area_units(imnrl),PETSC_TRUE)
          if (InputError(input)) then
            mineral_constraint%constraint_area_units(imnrl) = 'm^2/m^3'
            input%err_buf = trim(mineral_constraint%names(imnrl)) // &
                             ' SPECIFIC SURFACE_AREA UNITS'
            call InputDefaultMsg(input,option)
          endif
        enddo  
        
        if (imnrl < reaction%mineral%nkinmnrl) then
          option%io_buffer = &
                   'Mineral lists in constraints must provide a volume ' // &
                   'fraction and surface area for all kinetic minerals ' // &
                   '(listed under MINERAL_KINETICS card in CHEMISTRY), ' // &
                   'regardless of whether or not they are present (just ' // &
                   'assign a zero volume fraction if not present).'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%minerals)) then
          call MineralConstraintDestroy(constraint%minerals)
        endif
        constraint%minerals => mineral_constraint 

      case('COLL','COLLOIDS')

        colloid_constraint => ColloidConstraintCreate(reaction,option)

        block_string = 'CONSTRAINT, COLLOIDS'
        icomp = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit          
          
          icomp = icomp + 1

          if (icomp > reaction%ncoll) then
            option%io_buffer = &
                     'Number of colloid constraints exceeds number of ' // &
                     'colloids in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option,colloid_constraint%names(icomp), &
                          PETSC_TRUE)
          call InputErrorMsg(input,option,'colloid name',block_string)
          option%io_buffer = 'Constraint Colloids: ' // &
                             trim(colloid_constraint%names(icomp))
          call printMsg(option)
          call InputReadDouble(input,option, &
                               colloid_constraint%constraint_conc_mob(icomp))
          call InputErrorMsg(input,option,'mobile concentration',block_string)
          call InputReadDouble(input,option, &
                               colloid_constraint%constraint_conc_imb(icomp))
          call InputErrorMsg(input,option,'immobile concentration', &
                             block_string)
        
        enddo  
        
        if (icomp < reaction%ncoll) then
          option%io_buffer = &
                   'Colloid lists in constraints must provide mobile ' // &
                   'and immobile concentrations for all colloids ' // &
                   '(listed under the COLLOIDS card in CHEMISTRY), ' // &
                   'regardless of whether or not they are present (just ' // &
                   'assign a small value (e.g. 1.d-40) if not present).'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%colloids)) then
          call ColloidConstraintDestroy(constraint%colloids)
        endif
        constraint%colloids => colloid_constraint 

        
        
      case('IMMOBILE')

        immobile_constraint => &
          ImmobileConstraintCreate(reaction%immobile,option)

        block_string = 'CONSTRAINT, IMMOBILE'
        iimmobile = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,block_string)
          
          if (InputCheckExit(input,option)) exit          
          
          iimmobile = iimmobile + 1

          if (iimmobile > reaction%immobile%nimmobile) then
            option%io_buffer = &
                     'Number of immobile constraints exceeds number of ' // &
                     'immobile species in constraint: ' // &
                      trim(constraint%name)
            call printErrMsg(option)
          endif
          
          call InputReadWord(input,option, &
                             immobile_constraint%names(iimmobile),PETSC_TRUE)
          call InputErrorMsg(input,option,'immobile name',block_string)
          option%io_buffer = 'Constraint Immobile: ' // &
                             trim(immobile_constraint%names(iimmobile))
          call printMsg(option)

          ! concentration
          string = trim(input%buf)
          call InputReadWord(string,word,PETSC_TRUE,ierr)
          ! if a dataset
          if (StringCompareIgnoreCase(word,'DATASET')) then
            input%buf = trim(string)
            call InputReadWord(input,option,immobile_constraint% &
                                constraint_aux_string(iimmobile),PETSC_TRUE)
            call InputErrorMsg(input,option,'dataset name', &
                               trim(block_string) // ' concentration')
            immobile_constraint%external_dataset(iimmobile) = PETSC_TRUE
            ! set vol frac to NaN to catch bugs
            tempreal = -1.d0
            immobile_constraint%constraint_conc(iimmobile) = sqrt(tempreal)
          else
            call InputReadDouble(input,option, &
                                 immobile_constraint%constraint_conc(iimmobile))
            call InputErrorMsg(input,option,'concentration',block_string)
          endif

          ! read units if they exist
          internal_units = 'mol/m^3'
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) then
            input%err_buf = trim(immobile_constraint%names(iimmobile)) // &
                             ' IMMOBILE CONCENTRATION UNITS'
            call InputDefaultMsg(input,option)
          else
            immobile_constraint%constraint_conc(iimmobile) = &
              immobile_constraint%constraint_conc(iimmobile) * &
              UnitsConvertToInternal(word,internal_units,option)
          endif
        enddo  
        
        if (iimmobile < reaction%immobile%nimmobile) then
          option%io_buffer = &
                   'Immobile lists in constraints must provide a ' // &
                   'concentration for all immobile species ' // &
                   '(listed under IMMOBILE card in CHEMISTRY), ' // &
                   'regardless of whether or not they are present.'
          call printErrMsg(option)        
        endif
        
        if (associated(constraint%immobile_species)) then
          call ImmobileConstraintDestroy(constraint%immobile_species)
        endif
        constraint%immobile_species => immobile_constraint 
        
      case default
        call InputKeywordUnrecognized(word,'CONSTRAINT',option)
    end select 
  
  enddo  
  
  call PetscLogEventEnd(logging%event_tran_constraint_read,ierr);CHKERRQ(ierr)

end subroutine TranConstraintRead

! ************************************************************************** !

subroutine TranConstraintInitList(list)
  ! 
  ! Initializes a transport constraint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none

  type(tran_constraint_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_constraints = 0

end subroutine TranConstraintInitList

! ************************************************************************** !

subroutine TranConstraintAddToList(new_constraint,list)
  ! 
  ! Adds a new constraint to a transport constraint
  ! list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none
  
  type(tran_constraint_type), pointer :: new_constraint
  type(tran_constraint_list_type) :: list
  
  list%num_constraints = list%num_constraints + 1
  new_constraint%id = list%num_constraints
  if (.not.associated(list%first)) list%first => new_constraint
  if (associated(list%last)) list%last%next => new_constraint
  list%last => new_constraint
  
end subroutine TranConstraintAddToList

! ************************************************************************** !

function TranConstraintGetPtrFromList(constraint_name,constraint_list)
  ! 
  ! Returns a pointer to the constraint matching
  ! constraint_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  ! 

  use String_module

  implicit none
  
  type(tran_constraint_type), pointer :: TranConstraintGetPtrFromList
  character(len=MAXWORDLENGTH) :: constraint_name
  type(tran_constraint_list_type) :: constraint_list
 
  PetscInt :: length
  type(tran_constraint_type), pointer :: constraint
    
  nullify(TranConstraintGetPtrFromList)
  constraint => constraint_list%first
  
  do 
    if (.not.associated(constraint)) exit
    length = len_trim(constraint_name)
    if (length == len_trim(constraint%name) .and. &
        StringCompare(constraint%name,constraint_name, &
                        length)) then
      TranConstraintGetPtrFromList => constraint
      return
    endif
    constraint => constraint%next
  enddo
  
end function TranConstraintGetPtrFromList

! ************************************************************************** !

subroutine TranConstraintDestroy(constraint)
  ! 
  ! Deallocates a constraint
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none
  
  type(tran_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return

  if (associated(constraint%aqueous_species)) &
    call AqueousSpeciesConstraintDestroy(constraint%aqueous_species)
  nullify(constraint%aqueous_species)
  if (associated(constraint%free_ion_guess)) &
    call GuessConstraintDestroy(constraint%free_ion_guess)
  nullify(constraint%free_ion_guess)
  if (associated(constraint%minerals)) &
    call MineralConstraintDestroy(constraint%minerals)
  nullify(constraint%minerals)
  if (associated(constraint%colloids)) &
    call ColloidConstraintDestroy(constraint%colloids)
  nullify(constraint%colloids)
  if (associated(constraint%immobile_species)) &
    call ImmobileConstraintDestroy(constraint%immobile_species)
  nullify(constraint%immobile_species)

  deallocate(constraint)
  nullify(constraint)

end subroutine TranConstraintDestroy

! ************************************************************************** !

subroutine TranConstraintDestroyList(constraint_list)
  ! 
  ! Deallocates a list of constraints
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  implicit none
  
  type(tran_constraint_list_type), pointer :: constraint_list
  
  type(tran_constraint_type), pointer :: constraint, prev_constraint
  
  if (.not.associated(constraint_list)) return
  
  constraint => constraint_list%first
  do 
    if (.not.associated(constraint)) exit
    prev_constraint => constraint
    constraint => constraint%next
    call TranConstraintDestroy(prev_constraint)
  enddo
  
  constraint_list%num_constraints = 0
  nullify(constraint_list%first)
  nullify(constraint_list%last)
  if (associated(constraint_list%array)) deallocate(constraint_list%array)
  nullify(constraint_list%array)
  
  deallocate(constraint_list)
  nullify(constraint_list)

end subroutine TranConstraintDestroyList

! ************************************************************************** !

subroutine TranConstraintCouplerDestroy(coupler_list)
  ! 
  ! Destroys a constraint coupler linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use Option_module
  
  implicit none
  
  type(tran_constraint_coupler_type), pointer :: coupler_list
  
  type(tran_constraint_coupler_type), pointer :: cur_coupler, prev_coupler
  
  cur_coupler => coupler_list
  
  do
    if (.not.associated(cur_coupler)) exit
    prev_coupler => cur_coupler
    cur_coupler => cur_coupler%next
    if (associated(prev_coupler%rt_auxvar)) then
      call RTAuxVarDestroy(prev_coupler%rt_auxvar)
    endif
    nullify(prev_coupler%rt_auxvar)
    if (associated(prev_coupler%global_auxvar)) then
      call GlobalAuxVarDestroy(prev_coupler%global_auxvar)
    endif
    nullify(prev_coupler%global_auxvar)
    nullify(prev_coupler%aqueous_species)
    nullify(prev_coupler%minerals)

    nullify(prev_coupler%colloids)
    nullify(prev_coupler%immobile_species)
    nullify(prev_coupler%next)
    deallocate(prev_coupler)
    nullify(prev_coupler)
  enddo
  
  nullify(coupler_list)
  
end subroutine TranConstraintCouplerDestroy

end module Transport_Constraint_module
