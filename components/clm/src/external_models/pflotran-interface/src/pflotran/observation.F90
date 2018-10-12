module Observation_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Region_module
  use Connection_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  

  PetscInt, parameter, public :: OBSERVATION_SCALAR = 1
  PetscInt, parameter, public :: OBSERVATION_FLUX = 2
  PetscInt, parameter, public :: OBSERVATION_AT_CELL_CENTER = 1
  PetscInt, parameter, public :: OBSERVATION_AT_COORDINATE = 2

  type, public :: observation_type
    ! all added variables must be included in ObservationCreateFromObservation
    PetscInt :: id
    PetscInt :: itype
    PetscBool :: print_velocities
    PetscBool :: print_secondary_data(5)          ! first entry is for temp., second is for conc. and third is for mineral vol frac., fourth is rate, fifth is SI
    PetscBool :: at_cell_center
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: linkage_name
    type(connection_set_type), pointer :: connection_set
    type(region_type), pointer :: region
    type(observation_type), pointer :: next
  end type observation_type
  
  type, public :: observation_list_type
    PetscInt :: num_observations
    type(observation_type), pointer :: first
    type(observation_type), pointer :: last
    type(observation_type), pointer :: array(:)
  end type observation_list_type

  public :: ObservationCreate, ObservationDestroy, ObservationRead, &
            ObservationAddToList, ObservationInitList, ObservationDestroyList, &
            ObservationGetPtrFromList, ObservationRemoveFromList

  interface ObservationCreate
    module procedure ObservationCreate1
    module procedure ObservationCreateFromObservation
  end interface
    
contains

! ************************************************************************** !

function ObservationCreate1()
  ! 
  ! Create object that stores observation regions
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none
  
  type(observation_type), pointer :: ObservationCreate1
  
  type(observation_type), pointer :: observation
  
  allocate(observation)
  
  observation%name = ""
  observation%linkage_name = ""
  observation%id = 0
  observation%itype = OBSERVATION_SCALAR
  observation%print_velocities = PETSC_FALSE
  observation%at_cell_center = PETSC_TRUE
  observation%print_secondary_data = PETSC_FALSE
  nullify(observation%region)
  nullify(observation%next)
  
  ObservationCreate1 => observation

end function ObservationCreate1

! ************************************************************************** !

function ObservationCreateFromObservation(observation)
  ! 
  ! Create object that stores observation regions
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none
  
  type(observation_type), pointer :: ObservationCreateFromObservation
  type(observation_type), pointer :: observation

  type(observation_type), pointer :: new_observation
  
  new_observation => ObservationCreate1()
  
  new_observation%name = observation%name
  new_observation%linkage_name = observation%linkage_name
  new_observation%id = observation%id
  new_observation%itype = observation%itype
  new_observation%print_velocities = observation%print_velocities
  new_observation%at_cell_center = observation%at_cell_center
  new_observation%print_secondary_data = &
       observation%print_secondary_data
  ! keep these null for now to catch bugs
  nullify(new_observation%region)
  nullify(new_observation%next)
  
  ObservationCreateFromObservation => new_observation

end function ObservationCreateFromObservation

! ************************************************************************** !

subroutine ObservationRead(observation,input,option)
  ! 
  ! Reads observation data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(observation_type) :: observation
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','OBSERVATION')   
      
    select case(trim(keyword))
    
      case('BOUNDARY_CONDITION')
        call InputReadWord(input,option,observation%linkage_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'boundary condition name','OBSERVATION')
        option%flow%store_fluxes = PETSC_TRUE
        option%transport%store_fluxes = PETSC_TRUE
        observation%itype = OBSERVATION_FLUX
      case('REGION')
        call InputReadWord(input,option,observation%linkage_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'region name','OBSERVATION')
        observation%itype = OBSERVATION_SCALAR
      case('VELOCITY')
        observation%print_velocities = PETSC_TRUE
      case('SECONDARY_TEMPERATURE')
        if (option%use_mc) then
          observation%print_secondary_data(1) = PETSC_TRUE
        else
          option%io_buffer = 'Keyword SECONDARY_TEMPERATURE can only be used' // &
                             ' MULTIPLE_CONTINUUM keyword'
          call printErrMsg(option)
        endif
      case('SECONDARY_CONCENTRATION')
        if (option%use_mc) then
          observation%print_secondary_data(2) = PETSC_TRUE
        else
          option%io_buffer = 'Keyword SECONDARY_CONCENTRATION can only be used' // &
                             ' MULTIPLE_CONTINUUM keyword'
          call printErrMsg(option)
        endif
      case('SECONDARY_MINERAL_VOLFRAC')
        if (option%use_mc) then
          observation%print_secondary_data(3) = PETSC_TRUE
        else
          option%io_buffer = 'Keyword SECONDARY_MINERAL_VOLFRAC can ' // &
                             'only be used MULTIPLE_CONTINUUM keyword'
          call printErrMsg(option)
        endif
      case('SECONDARY_MINERAL_RATE')
        if (option%use_mc) then
          observation%print_secondary_data(4) = PETSC_TRUE
        else
          option%io_buffer = 'Keyword SECONDARY_MINERAL_RATE can ' // &
                             'only be used MULTIPLE_CONTINUUM keyword'
          call printErrMsg(option)
        endif
      case('SECONDARY_MINERAL_SI')
        if (option%use_mc) then
          observation%print_secondary_data(5) = PETSC_TRUE
        else
          option%io_buffer = 'Keyword SECONDARY_MINERAL_SI can ' // &
                             'only be used MULTIPLE_CONTINUUM keyword'
          call printErrMsg(option)
        endif
      case('AT_CELL_CENTER')
        observation%at_cell_center = PETSC_TRUE
      case('AT_COORDINATE')
        observation%at_cell_center = PETSC_FALSE
      case default
        call InputKeywordUnrecognized(keyword,'OBSERVATION',option)
    end select 
  
  enddo  

end subroutine ObservationRead

! ************************************************************************** !

subroutine ObservationInitList(list)
  ! 
  ! Initializes a observation list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none

  type(observation_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_observations = 0

end subroutine ObservationInitList

! ************************************************************************** !

subroutine ObservationAddToList(new_observation,list)
  ! 
  ! Adds a new observation to a observation list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none
  
  type(observation_type), pointer :: new_observation
  type(observation_list_type) :: list
  
  list%num_observations = list%num_observations + 1
  new_observation%id = list%num_observations
  if (.not.associated(list%first)) list%first => new_observation
  if (associated(list%last)) list%last%next => new_observation
  list%last => new_observation
  
end subroutine ObservationAddToList

! ************************************************************************** !

subroutine ObservationRemoveFromList(observation,list)
  ! 
  ! Removes a observation from a observation list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none
  
  type(observation_type), pointer :: observation
  type(observation_list_type) :: list
  
  type(observation_type), pointer :: cur_observation, prev_observation
  
  cur_observation => list%first
  nullify(prev_observation)
  
  do
    if (.not.associated(cur_observation)) exit
    if (associated(cur_observation,observation)) then
      if (associated(prev_observation)) then
        prev_observation%next => cur_observation%next
      else
        list%first => cur_observation%next
      endif
      if (.not.associated(cur_observation%next)) then
        list%last => prev_observation
      endif
      list%num_observations = list%num_observations-1
      call ObservationDestroy(cur_observation)
      return
    endif
    prev_observation => cur_observation
    cur_observation => cur_observation%next
  enddo
  
end subroutine ObservationRemoveFromList

! ************************************************************************** !

function ObservationGetPtrFromList(observation_name,observation_list)
  ! 
  ! Returns a pointer to the observation matching &
  ! observation_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  use String_module

  implicit none
  
  type(observation_type), pointer :: ObservationGetPtrFromList
  character(len=MAXWORDLENGTH) :: observation_name
  type(observation_list_type) :: observation_list
 
  PetscInt :: length
  type(observation_type), pointer :: observation
    
  nullify(ObservationGetPtrFromList)
  observation => observation_list%first
  
  do 
    if (.not.associated(observation)) exit
    length = len_trim(observation_name)
    if (length == len_trim(observation%name) .and. &
        StringCompare(observation%name,observation_name, &
                        length)) then
      ObservationGetPtrFromList => observation
      return
    endif
    observation => observation%next
  enddo
  
end function ObservationGetPtrFromList

! ************************************************************************** !

subroutine ObservationDestroyList(observation_list)
  ! 
  ! Deallocates a list of observations
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  ! 

  implicit none
  
  type(observation_list_type), pointer :: observation_list
  
  type(observation_type), pointer :: observation, prev_observation
  
  if (.not.associated(observation_list)) return
  
  observation => observation_list%first
  do 
    if (.not.associated(observation)) exit
    prev_observation => observation
    observation => observation%next
    call ObservationDestroy(prev_observation)
  enddo
  
  observation_list%num_observations = 0
  nullify(observation_list%first)
  nullify(observation_list%last)
  if (associated(observation_list%array)) deallocate(observation_list%array)
  nullify(observation_list%array)
  
  deallocate(observation_list)
  nullify(observation_list)

end subroutine ObservationDestroyList

! ************************************************************************** !

subroutine ObservationDestroy(observation)
  ! 
  ! Deallocates a observation
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(observation_type), pointer :: observation
  
  PetscInt :: i
  
  if (.not.associated(observation)) return
  
  nullify(observation%region)
  nullify(observation%connection_set)
  deallocate(observation)
  nullify(observation)

end subroutine ObservationDestroy

end module Observation_module
