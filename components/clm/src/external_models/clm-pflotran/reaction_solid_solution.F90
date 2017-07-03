module Reaction_Solid_Solution_module

  use Reaction_Mineral_Aux_module
  use Reaction_Aux_module
  use Reaction_Solid_Soln_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  public :: SolidSolutionReadFromInputFile, &
            SolidSolutionLinkNamesToIDs
            
contains

! ************************************************************************** !

subroutine SolidSolutionReadFromInputFile(solid_solution_list,input, &
                                          option)
  ! 
  ! Reads solid solution from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/16/12
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  
  implicit none
  
  type(solid_solution_type), pointer :: solid_solution_list
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name  
  character(len=MAXWORDLENGTH) :: card
  PetscInt, parameter :: max_stoich_solid_names = 200
  PetscInt :: stoich_solid_count
  character(len=MAXWORDLENGTH) :: stoich_solid_names(max_stoich_solid_names)
  type(solid_solution_type), pointer :: solid_solution, prev_solid_solution
  
  nullify(prev_solid_solution)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    ! name of solid solution
    solid_solution => SolidSolutionCreate()
    if (associated(prev_solid_solution)) then
      prev_solid_solution%next => solid_solution
      prev_solid_solution => solid_solution
    else
      solid_solution_list => solid_solution
    endif
    
    call InputReadWord(input,option,solid_solution%name,PETSC_TRUE)  
    call InputErrorMsg(input,option,'Solid Solution Name', &
                       'CHEMISTRY,SOLID_SOLUTIONS')

    stoich_solid_count = 0
    stoich_solid_names = ''
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      if (InputCheckExit(input,option)) exit
      
      stoich_solid_count = stoich_solid_count + 1
      
      if (stoich_solid_count > max_stoich_solid_names) then
        write(string,*) max_stoich_solid_names
        option%io_buffer = '# stoichmetric solids exceeds max (' // &
          trim(adjustl(string)) // ').  ' // &
          'Increase variable "max_stoich_solid_names" in ' // &
          'SolidSolutionReadFromInputFile.'
        call printErrMsg(option)
      endif
      call InputReadWord(input,option, &
                         stoich_solid_names(stoich_solid_count), &
                         PETSC_TRUE)  
      call InputErrorMsg(input,option,'Stoichiometric Solid Name', &
                         'CHEMISTRY,SOLID_SOLUTIONS')
    enddo
    
    allocate(solid_solution%stoich_solid_ids(stoich_solid_count))
    solid_solution%stoich_solid_ids = 0
    allocate(solid_solution%stoich_solid_names(stoich_solid_count))
    solid_solution%stoich_solid_names(1:stoich_solid_count) = &
      stoich_solid_names(1:stoich_solid_count)
    solid_solution%num_stoich_solid = stoich_solid_count

#if 0
    string = input%buf
    call InputReadWord(input,option,card,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SOLID_SOLUTIONS')
    call StringToUpper(card)
    select case(card)
      case('DATABASE')
        call InputReadNChars(string, &
                             solid_solution_rxn%database_filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE,input%ierr)  
        call InputErrorMsg(input,option,'keyword', &
                           'CHEMISTRY,SOLID_SOLUTIONS,DATABASE FILENAME')        
      case default
        solid_solution => SolidSolutionCreate()
        call InputReadWord(input,option,solid_solution%name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'keyword','CHEMISTRY,SOLID_SOLUTIONS')   
        if (.not.associated(solid_solution_rxn%list)) then
          solid_solution_rxn%list => solid_solution
        endif
        if (associated(prev_solid_solution)) then
          prev_solid_solution%next => solid_solution
        endif
        prev_solid_solution => solid_solution
        nullify(solid_solution)
    end select
#endif  
  enddo
  
end subroutine SolidSolutionReadFromInputFile

! ************************************************************************** !

subroutine SolidSolutionLinkNamesToIDs(solid_solution_list, &
                                       mineral_reaction, &
                                       option)
  ! 
  ! SolidSolutionReadFromDatabase: Reads solid solution from the database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/20/12
  ! 
  use Option_module
  use String_module
  use Reaction_Mineral_Aux_module
  
  implicit none
  
  type(solid_solution_type), pointer :: solid_solution_list
  type(mineral_type), pointer :: mineral_reaction
  type(option_type) :: option
  
  type(solid_solution_type), pointer :: cur_solid_soln
  PetscInt :: istoich_solid
  PetscInt :: ikinmnrl
  
  cur_solid_soln => solid_solution_list
  do
    if (.not.associated(cur_solid_soln)) exit
    do istoich_solid = 1, cur_solid_soln%num_stoich_solid
      do ikinmnrl = 1, mineral_reaction%nkinmnrl
        if (StringCompareIgnoreCase( &
              mineral_reaction%kinmnrl_names(ikinmnrl), &
              cur_solid_soln%stoich_solid_names(istoich_solid))) then
          cur_solid_soln%stoich_solid_ids(istoich_solid) = ikinmnrl
          exit
        endif
      enddo
    enddo
    cur_solid_soln => cur_solid_soln%next
  enddo

end subroutine SolidSolutionLinkNamesToIDs

#if 0

! ************************************************************************** !

subroutine SolidSolutionReadFromDatabase(solid_solution_rxn,option)
  ! 
  ! Reads solid solution from the database
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/20/12
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module
  use Reaction_Mineral_module
  
  implicit none
  
  type(solid_solution_rxn_type) :: solid_solution_rxn
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: name  
  character(len=MAXWORDLENGTH) :: card
  type(input_type), pointer :: input
  type(solid_solution_type), pointer :: solid_solution, prev_solid_solution
  type(stoichiometric_solid_type), pointer :: stoich_solid, prev_stoich_solid
  type(mineral_rxn_type), pointer :: mineral, prev_end_member
  PetscInt :: itemp
  PetscBool :: found
           
  if (len_trim(solid_solution_rxn%database_filename) < 2) then
    option%io_buffer = 'Database filename not included in input deck.'
    call printErrMsg(option)
  endif
  input => InputCreate(IUNIT_TEMP,solid_solution_rxn%database_filename,option)

  do ! loop over every entry in the database
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'SolidSolutionReadFromDatabase')

    call InputReadWord(input,option,card,PETSC_TRUE)  
    call InputErrorMsg(input,option,'keyword','CHEMISTRY,SOLID_SOLUTIONS')
    call StringToUpper(card)

    select case(card)
      case('SOLID_SOLUTION')
      
        if (solid_solution_rxn%num_dbase_temperatures == 0) then
          option%io_buffer = 'Temperatures must be defined prior to ' // &
                             'reading solid solution.'
          call printErrMsg(option)
        endif
      
        call InputReadWord(input,option,name,PETSC_TRUE)
        call InputErrorMsg(input,option,'SOLID_SOLUTION Name', &
                           'CHEMISTRY,SOLID_SOLUTIONS')
        solid_solution => solid_solution_rxn%list
        found = PETSC_FALSE
        do
          if (.not.associated(solid_solution)) exit
          if (StringCompare(name,solid_solution%name,MAXWORDLENGTH)) then
            found = PETSC_TRUE
            exit
          endif
          solid_solution => solid_solution%next
        enddo
        ! if solid solution not in list, skip to end of solid solution
        if (.not.found) then 
          call InputSkipToEND(input,option,card)
        endif
        nullify(prev_stoich_solid)
        nullify(prev_end_member)
      case('STOICHIOMETRIC_SOLID','END_MEMBER')
        mineral => MineralRxnCreate()
        call InputReadWord(input,option,mineral%name,PETSC_TRUE)  
        call InputErrorMsg(input,option,'keyword','CHEMISTRY,MINERALS')    
        call MineralReadFromDatabase(mineral, &
                                   solid_solution_rxn%num_dbase_temperatures, &
                                   input,option)
        ! assign mineral to stoichiometric solid or end member list
        select case(card)
          case('STOICHIOMETRIC_SOLID')
            solid_solution%num_stoich_solid = &
              solid_solution%num_stoich_solid + 1
            stoich_solid => StoichiometricSolidCreate()
            stoich_solid%mineral => mineral
            if (associated(prev_stoich_solid)) then
              prev_stoich_solid%next => stoich_solid
            else
              solid_solution%stoich_solid => stoich_solid
              prev_stoich_solid => stoich_solid
            endif
            nullify(prev_end_member)
          case('END_MEMBER')
            solid_solution%num_end_member = &
              solid_solution%num_end_member + 1
            if (associated(prev_end_member)) then
              prev_end_member%next => mineral
            else
              stoich_solid%end_members => mineral
              prev_end_member => mineral
            endif
        end select
        nullify(mineral)  
      case('TEMPERATURES')
        string = 'Temperatures in SolidSolutionReadFromDatabase'
        call UtilityReadRealArray(solid_solution_rxn%dbase_temperatures, &
                                  ZERO_INTEGER,string,input,option)
        solid_solution_rxn%num_dbase_temperatures = &
          size(solid_solution_rxn%dbase_temperatures)
      case default
        call InputKeywordUnrecognized(word, &
                     'SOLID SOLUTION,DATABASE',option)
    end select
  enddo
    
end subroutine SolidSolutionReadFromDatabase
#endif

end module Reaction_Solid_Solution_module
