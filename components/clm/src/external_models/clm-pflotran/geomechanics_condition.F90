module Geomechanics_Condition_module
 
!  use Global_Aux_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module
  
  use PFLOTRAN_Constants_module
  
  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

#if 0
!geh: no longer needed
  PetscInt, parameter :: NULL = 0
  PetscInt, parameter :: STEP = 1
  PetscInt, parameter :: LINEAR = 2
  
  type, public :: geomech_condition_dataset_type
    type(time_series_type), pointer :: time_series
    class(dataset_base_type), pointer :: dataset
  end type geomech_condition_dataset_type  
#endif

  type, public :: geomech_condition_type
    PetscInt :: id    ! id from which condition can be referenced
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name    ! name of condition (e.g. boundary)
    PetscInt :: num_sub_conditions
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(time_storage_type), pointer :: default_time_storage
    type(geomech_sub_condition_type), pointer :: displacement_x
    type(geomech_sub_condition_type), pointer :: displacement_y
    type(geomech_sub_condition_type), pointer :: displacement_z
    type(geomech_sub_condition_type), pointer :: force_x ! Added force conditions 09/19/2013, SK
    type(geomech_sub_condition_type), pointer :: force_y 
    type(geomech_sub_condition_type), pointer :: force_z
    type(geomech_sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(geomech_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type geomech_condition_type
    
  type, public :: geomech_sub_condition_type
    PetscInt :: itype ! integer describing type of condition
    PetscInt :: isubtype
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    class(dataset_base_type), pointer :: dataset
  end type geomech_sub_condition_type
  
  type, public :: geomech_sub_condition_ptr_type
    type(geomech_sub_condition_type), pointer :: ptr
  end type geomech_sub_condition_ptr_type
    
  type, public :: geomech_condition_ptr_type
    type(geomech_condition_type), pointer :: ptr
  end type geomech_condition_ptr_type
  
  type, public :: geomech_condition_list_type
    PetscInt :: num_conditions
    type(geomech_condition_type), pointer :: first
    type(geomech_condition_type), pointer :: last
    type(geomech_condition_type), pointer :: array(:)    
  end type geomech_condition_list_type

  public :: GeomechConditionCreate, &
            GeomechConditionDestroy, &
            GeomechConditionRead, &
            GeomechConditionAddToList, &
            GeomechConditionInitList, &
            GeomechConditionDestroyList, &
            GeomechConditionGetPtrFromList, &
            GeomechConditionUpdate, &
            GeomechConditionPrint, &
            GeomechConditionIsTransient
            
  
contains

! ************************************************************************** !

function GeomechConditionCreate(option)
  ! 
  ! Creates a condition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/07/13
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(geomech_condition_type), pointer :: GeomechConditionCreate
  
  type(geomech_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)
  nullify(condition%force_x)
  nullify(condition%force_y)
  nullify(condition%force_z)
  nullify(condition%sub_condition_ptr)
  nullify(condition%itype)
  nullify(condition%next)
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%num_sub_conditions = 0
  condition%name = ''
  
  GeomechConditionCreate => condition

end function GeomechConditionCreate

! ************************************************************************** !

function GeomechSubConditionCreate(ndof)
  ! 
  ! Creates a sub_condition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module
  
  implicit none
  
  type(geomech_sub_condition_type), pointer :: GeomechSubConditionCreate
  
  PetscInt :: ndof
  
  type(geomech_sub_condition_type), pointer :: sub_condition
  class(dataset_ascii_type), pointer :: dataset_ascii
  
  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = 0
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''
  nullify(sub_condition%dataset)  

  ! by default, all dataset are of type dataset_ascii_type, unless overwritten
  dataset_ascii => DatasetAsciiCreate()
  call DatasetAsciiInit(dataset_ascii)
  dataset_ascii%array_width = ndof
  dataset_ascii%data_type = DATASET_REAL
  sub_condition%dataset => dataset_ascii
  nullify(dataset_ascii)

  GeomechSubConditionCreate => sub_condition

end function GeomechSubConditionCreate

! ************************************************************************** !

subroutine GeomechSubConditionVerify(option, condition, sub_condition_name, &
                                     sub_condition, default_time_storage, &
                                     destroy_if_null)
  ! 
  ! Verifies the data in a subcondition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module
  use Dataset_module

  implicit none

  type(option_type) :: option
  type(geomech_condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(geomech_sub_condition_type), pointer :: sub_condition
  type(time_storage_type), pointer :: default_time_storage
  PetscBool :: destroy_if_null

  character(len=MAXSTRINGLENGTH) :: header
  
  if (.not.associated(sub_condition)) return
  
 ! dataset is not optional
  if (.not.(associated(sub_condition%dataset%rarray) .or. &
            associated(sub_condition%dataset%rbuffer) .or. &
            ! if a dataset name is read, instead of data at this point
            len_trim(sub_condition%dataset%name) > 0)) then
    if (destroy_if_null) call GeomechSubConditionDestroy(sub_condition)
    return
  endif
  
  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call printErrMsg(option)
  endif
  
  header = 'GEOMECHANICS_CONDITION' // '/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Value(s)'
  call DatasetVerify(sub_condition%dataset,default_time_storage, &
                     header,option)

end subroutine GeomechSubConditionVerify

! ************************************************************************** !

subroutine GeomechConditionRead(condition,input,option)
  ! 
  ! Reads a condition from the input file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Condition_module
  
  implicit none
  
  type(geomech_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, internal_units
  type(geomech_sub_condition_type), pointer :: sub_condition_ptr,  &
                                       displacement_x, displacement_y, &
                                       displacement_z
  type(geomech_sub_condition_type), pointer :: force_x, force_y, force_z 
  PetscReal :: default_time
  PetscInt :: default_iphase
  character(len=MAXWORDLENGTH) :: default_ctype
  PetscInt :: default_itype
  PetscInt :: array_size, idof
  PetscBool :: found
  PetscBool :: destroy_if_null
  PetscErrorCode :: ierr
  PetscInt :: num_sub_conditions
  PetscInt :: count
  !geh: may not need default_time_storage
  type(time_storage_type), pointer :: default_time_storage

  default_time = 0.d0
  default_iphase = 0
  
  !geh: may not need default_time_storage
  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP  
  
#if 0
!geh: no longer needed
  call GeomechConditionDatasetInit(default_geomech_dataset)
  default_geomech_dataset%time_series => TimeSeriesCreate()
  default_geomech_dataset%time_series%rank = 1
  default_geomech_dataset%time_series%interpolation_method = STEP
  default_geomech_dataset%time_series%is_cyclic = PETSC_FALSE
#endif

  displacement_x => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_y => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_z => GeomechSubConditionCreate(ONE_INTEGER)
  force_x => GeomechSubConditionCreate(ONE_INTEGER)
  force_y => GeomechSubConditionCreate(ONE_INTEGER)
  force_z => GeomechSubConditionCreate(ONE_INTEGER)
  displacement_x%name = 'displacement_x'
  displacement_y%name = 'displacement_y'
  displacement_z%name = 'displacement_z'
  force_x%name = 'force_x'
  force_y%name = 'force_y'
  force_z%name = 'force_z'
  
  condition%time_units = 'yr'
  condition%length_units = 'm'
  
  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  displacement_x%units = 'm'
  displacement_y%units = 'm'
  displacement_z%units = 'm'
  force_x%units = 'N'
  force_y%units = 'N'
  force_z%units = 'N'

  default_ctype = 'dirichlet'
  default_itype = DIRICHLET_BC

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('UNITS') ! read default units for condition arguments
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
          end select
        enddo
      case('CYCLIC')
        ! by default, is_cyclic is set to PETSC_FALSE
        default_time_storage%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')   
        call StringToLower(word)
        select case(word)
          case('step')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_STEP
          case('linear') 
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(trim(word))
            case('PRESSURE')
            case('DISPLACEMENT_X')
              sub_condition_ptr => displacement_x
            case('DISPLACEMENT_Y')
              sub_condition_ptr => displacement_y
            case('DISPLACEMENT_Z')
              sub_condition_ptr => displacement_z
            case('FORCE_X')
              sub_condition_ptr => force_x 
            case('FORCE_Y')
              sub_condition_ptr => force_y 
            case('FORCE_Z')
              sub_condition_ptr => force_z 
            case default
              call InputKeywordUnrecognized(word, &
                     'geomechanics condition type',option)
          end select
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'TYPE','CONDITION')   
          call StringToLower(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case default
              call InputKeywordUnrecognized(word, &
                     'geomechanics condition bc type',option)
          end select
        enddo
      case('TIME','TIMES')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION')   
      case('DISPLACEMENT_X')
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 displacement_x%dataset, &
                                 displacement_x%units, &
                                 internal_units)
      case('DISPLACEMENT_Y')
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 displacement_y%dataset, &
                                 displacement_y%units, &
                                 internal_units) 
      case('DISPLACEMENT_Z')
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 displacement_z%dataset, &
                                 displacement_z%units, &
                                 internal_units)
      case('FORCE_X')
      internal_units = 'N'
        call ConditionReadValues(input,option,word, &
                                 force_x%dataset, &
                                 force_x%units, &
                                 internal_units)
      case('FORCE_Y')
        internal_units = 'N'
        call ConditionReadValues(input,option,word, &
                                 force_y%dataset, &
                                 force_y%units, &
                                 internal_units)
      case('FORCE_Z')
        internal_units = 'N'
        call ConditionReadValues(input,option,word, &
                                 force_z%dataset, &
                                 force_z%units, &
                                 internal_units)
      case default
        call InputKeywordUnrecognized(word, &
                     'geomechanics condition',option)
    end select 
  
  enddo  
  
  word = 'displacement_x'
  call GeomechSubConditionVerify(option,condition,word,displacement_x, &
                                 default_time_storage, &
                                 PETSC_TRUE)
  word = 'displacement_y'
  call GeomechSubConditionVerify(option,condition,word,displacement_y, &
                                 default_time_storage, &
                                 PETSC_TRUE)
  word = 'displacement_z'
  call GeomechSubConditionVerify(option,condition,word,displacement_z, &
                                 default_time_storage, &
                                 PETSC_TRUE)

  word = 'force_x'
  call GeomechSubConditionVerify(option,condition,word,force_x, &
                                 default_time_storage, &
                                 PETSC_TRUE)

  word = 'force_y'
  call GeomechSubConditionVerify(option,condition,word,force_y, &
                                 default_time_storage, &
                                 PETSC_TRUE)

  word = 'force_z'
  call GeomechSubConditionVerify(option,condition,word,force_z, &
                                 default_time_storage, &
                                 PETSC_TRUE)



  num_sub_conditions = 0
  if (associated(displacement_x)) then
    condition%displacement_x => displacement_x
    num_sub_conditions = num_sub_conditions + 1
    condition%displacement_x%isubtype = ONE_INTEGER
  endif                         

  if (associated(displacement_y)) then
    condition%displacement_y => displacement_y
    num_sub_conditions = num_sub_conditions + 1
    condition%displacement_y%isubtype = TWO_INTEGER    
  endif                         

  if (associated(displacement_z)) then
    condition%displacement_z => displacement_z
    num_sub_conditions = num_sub_conditions + 1
    condition%displacement_z%isubtype = THREE_INTEGER
  endif                         

  if (associated(force_x)) then
    condition%force_x => force_x 
    num_sub_conditions = num_sub_conditions + 1
    condition%force_x%isubtype = FOUR_INTEGER 
  endif                         

  if (associated(force_y)) then
    condition%force_y => force_y 
    num_sub_conditions = num_sub_conditions + 1
    condition%force_y%isubtype = FIVE_INTEGER 
  endif                         

  if (associated(force_z)) then
    condition%force_z => force_z 
    num_sub_conditions = num_sub_conditions + 1
    condition%force_z%isubtype = THREE_INTEGER
  endif                         

  if (num_sub_conditions == 0) then
    option%io_buffer = 'displacement/force condition null in condition: ' // &
                        trim(condition%name)
    call printErrMsg(option)   
  endif

  condition%num_sub_conditions = num_sub_conditions
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo

  ! SK: I am using isubtype to differentiate between x, y, z in sub_condition_ptr
  ! since all of the displacements need not be specified.
  count = 0
  if (associated(displacement_x)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => displacement_x
  endif
  if (associated(displacement_y)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => displacement_y
  endif
  if (associated(displacement_z)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => displacement_z
  endif
  if (associated(force_x)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => force_x 
  endif
  if (associated(force_y)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => force_y 
  endif
  if (associated(force_z)) then
    count = count + 1
    condition%sub_condition_ptr(count)%ptr => force_z 
  endif    
  
  condition%default_time_storage => default_time_storage
    
end subroutine GeomechConditionRead

! ************************************************************************** !

subroutine GeomechConditionPrint(condition,option)
  ! 
  ! Prints Geomech condition info
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module

  implicit none
  
  type(geomech_condition_type) :: condition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Geomech Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') &
    trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') &
    trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') &
    trim(condition%length_units)
  
  do i=1, condition%num_sub_conditions
    call GeomechConditionPrintSubCondition(&
                                        condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)
  
end subroutine GeomechConditionPrint

! ************************************************************************** !

subroutine GeomechConditionPrintSubCondition(subcondition,option)
  ! 
  ! Prints Geomech subcondition info
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module

  implicit none
  
  type(geomech_sub_condition_type) :: subcondition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  select case(subcondition%itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
  end select
  100 format(6x,'Type: ',a)  
  write(option%fid_out,100) trim(string)
  
  110 format(6x,a)  

  write(option%fid_out,110) 'Geomech Dataset:'
  if (associated(subcondition%dataset)) then
!geh    call DatasetPrint(subcondition%dataset,option)
    option%io_buffer = 'TODO(geh): add DatasetPrint()'
    call printMsg(option)
  endif
            
end subroutine GeomechConditionPrintSubCondition

! ************************************************************************** !

subroutine GeomechConditionUpdate(condition_list,option,time)
  ! 
  ! Updates a transient condition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Option_module
  use Dataset_module
  
  implicit none
  
  type(geomech_condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  
  type(geomech_condition_type), pointer :: condition
  type(geomech_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition   
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr
      
      if (associated(sub_condition)) then
        call DatasetUpdate(sub_condition%dataset,time,option)
      endif
      
    enddo
      
    condition => condition%next
    
  enddo
  
end subroutine GeomechConditionUpdate

! ************************************************************************** !

subroutine GeomechConditionInitList(list)
  ! 
  ! Initializes a condition list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  implicit none

  type(geomech_condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine GeomechConditionInitList

! ************************************************************************** !

subroutine GeomechConditionAddToList(new_condition,list)
  ! 
  ! Adds a new condition to a condition list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  implicit none
  
  type(geomech_condition_type), pointer :: new_condition
  type(geomech_condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine GeomechConditionAddToList

! ************************************************************************** !

function GeomechConditionGetPtrFromList(condition_name,condition_list)
  ! 
  ! Returns a pointer to the condition matching &
  ! condition_name
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use String_module
  
  implicit none
  
  type(geomech_condition_type), pointer :: GeomechConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(geomech_condition_list_type) :: condition_list
 
  PetscInt :: length
  type(geomech_condition_type), pointer :: condition
    
  nullify(GeomechConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      GeomechConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function GeomechConditionGetPtrFromList

! ************************************************************************** !

function GeomechConditionIsTransient(condition)
  ! 
  ! Returns PETSC_TRUE for geomech condition if
  ! it is transient
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  implicit none
  
  type(geomech_condition_type) :: condition
 
  PetscBool :: GeomechConditionIsTransient
  
  GeomechConditionIsTransient = PETSC_FALSE

  if (GeomechSubConditionIsTransient(condition%displacement_x) .or. &
      GeomechSubConditionIsTransient(condition%displacement_y) .or. &
      GeomechSubConditionIsTransient(condition%displacement_z)) then
    GeomechConditionIsTransient = PETSC_TRUE
  endif
 
  if (GeomechSubConditionIsTransient(condition%force_x) .or. &
      GeomechSubConditionIsTransient(condition%force_y) .or. &
      GeomechSubConditionIsTransient(condition%force_z)) then
    GeomechConditionIsTransient = PETSC_TRUE
  endif
  
 
end function GeomechConditionIsTransient

! ************************************************************************** !

function GeomechSubConditionIsTransient(sub_condition)
  ! 
  ! Returns PETSC_TRUE for geomech sub condition
  ! if it is transient
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/12/13
  ! 

  use Dataset_module

  implicit none
  
  type(geomech_sub_condition_type), pointer :: sub_condition
  
  PetscBool :: GeomechSubConditionIsTransient
  
  GeomechSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (DatasetIsTransient(sub_condition%dataset)) then
      GeomechSubConditionIsTransient = PETSC_TRUE
    endif
  endif  
  
end function GeomechSubConditionIsTransient

! ************************************************************************** !

subroutine GeomechConditionDestroyList(condition_list)
  ! 
  ! Deallocates a list of conditions
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(geomech_condition_list_type), pointer :: condition_list
  
  type(geomech_condition_type), pointer :: condition, &
                                                       prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call GeomechConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine GeomechConditionDestroyList

! ************************************************************************** !

subroutine GeomechConditionDestroy(condition)
  ! 
  ! Deallocates a condition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(geomech_condition_type), pointer :: condition
  
  PetscInt :: i
  
  if (.not.associated(condition)) return
  
  if (associated(condition%sub_condition_ptr)) then
    do i=1,condition%num_sub_conditions
      call GeomechSubConditionDestroy(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  if (associated(condition%itype)) deallocate(condition%itype)
  nullify(condition%itype)
  
  nullify(condition%displacement_x)
  nullify(condition%displacement_y)
  nullify(condition%displacement_z)
  nullify(condition%force_x)
  nullify(condition%force_y)
  nullify(condition%force_z)
  
  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine GeomechConditionDestroy

! ************************************************************************** !

subroutine GeomechSubConditionDestroy(sub_condition)
  ! 
  ! Destroys a sub_condition
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/04/08
  ! 

  use Dataset_module
  use Dataset_Ascii_class
  
  implicit none
  
  type(geomech_sub_condition_type), pointer :: sub_condition
  class(dataset_ascii_type), pointer :: dataset_ascii
  
  if (.not.associated(sub_condition)) return
  
  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(sub_condition%dataset)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine GeomechSubConditionDestroy

end module Geomechanics_Condition_module
