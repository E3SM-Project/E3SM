module Condition_module
 
!  use Reaction_Aux_module
!  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module
  
  use Transport_Constraint_module
!  use Reaction_Surface_Complexation_Aux_module  
!  use Reaction_Mineral_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

  type, public :: flow_condition_type
    PetscInt :: id                          ! id from which condition can be referenced
    PetscBool :: is_transient
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name    ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: iphase
    PetscInt, pointer :: itype(:)
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(time_storage_type), pointer :: default_time_storage
    class(dataset_base_type), pointer :: datum
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: saturation
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: well
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: energy_rate
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_general_condition_type), pointer :: general
    type(flow_toil_ims_condition_type), pointer :: toil_ims  
    type(flow_well_condition_type), pointer :: flow_well  ! flow_well to avoid conflict with well
    ! any new sub conditions must be added to FlowConditionIsTransient
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type flow_condition_type

  ! data structure for general phase
  type, public :: flow_general_condition_type
    type(flow_sub_condition_type), pointer :: liquid_pressure
    type(flow_sub_condition_type), pointer :: gas_pressure
    type(flow_sub_condition_type), pointer :: gas_saturation
    type(flow_sub_condition_type), pointer :: mole_fraction
    type(flow_sub_condition_type), pointer :: relative_humidity
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: gas_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_general_condition_type

  ! data structure for toil_ims
  type, public :: flow_toil_ims_condition_type
    !type(flow_sub_condition_type), pointer :: liquid_pressure
    !type(flow_sub_condition_type), pointer :: oil_pressure
    !type(flow_sub_condition_type), pointer :: oil_saturation
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: saturation
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: rate
    !PO: to add well when adding well capabilities.
    type(flow_sub_condition_type), pointer :: liquid_flux
    type(flow_sub_condition_type), pointer :: oil_flux
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_sub_condition_type), pointer :: owc   ! oil water contact 
    type(flow_sub_condition_type), pointer :: liq_press_grad ! water piezometric head gradient
    ! any new sub conditions must be added to FlowConditionIsTransient
  end type flow_toil_ims_condition_type

  type, public :: flow_well_condition_type
    type(flow_sub_condition_type), pointer :: pressure
    type(flow_sub_condition_type), pointer :: rate
    type(flow_sub_condition_type), pointer :: temperature
    !when needed add here other variables such as WOR, WGR, etc 
    !any new sub conditions must be added to FlowConditionIsTransient
  end type flow_well_condition_type
    
  type, public :: flow_sub_condition_type
    PetscInt :: itype                  ! integer describing type of condition
    PetscInt :: isubtype
    character(len=MAXWORDLENGTH) :: ctype ! character string describing type of condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    PetscReal :: aux_real(2)
    class(dataset_base_type), pointer :: gradient
    class(dataset_base_type), pointer :: dataset
  end type flow_sub_condition_type
  
  type, public :: sub_condition_ptr_type
    type(flow_sub_condition_type), pointer :: ptr
  end type sub_condition_ptr_type
    
  type, public :: condition_ptr_type
    type(flow_condition_type), pointer :: ptr
  end type condition_ptr_type
  
  type, public :: condition_list_type
    PetscInt :: num_conditions
    type(flow_condition_type), pointer :: first
    type(flow_condition_type), pointer :: last
    type(flow_condition_type), pointer :: array(:)    
  end type condition_list_type
  
  type, public :: tran_condition_type
    PetscInt :: id                     ! id from which condition can be referenced
    PetscInt :: itype                  ! integer describing type of condition
    PetscBool :: is_transient
    character(len=MAXWORDLENGTH) :: name  ! name of condition (e.g. initial, recharge)
    type(tran_constraint_coupler_type), pointer :: constraint_coupler_list
    type(tran_constraint_coupler_type), pointer :: cur_constraint_coupler
    type(tran_condition_type), pointer :: next
  end type tran_condition_type
  
  type, public :: tran_condition_ptr_type
    type(tran_condition_type), pointer :: ptr
  end type tran_condition_ptr_type
  
  type, public :: tran_condition_list_type
    PetscInt :: num_conditions
    type(tran_condition_type), pointer :: first
    type(tran_condition_type), pointer :: last
    type(tran_condition_ptr_type), pointer :: array(:)    
  end type tran_condition_list_type
  
  public :: FlowConditionCreate, FlowConditionDestroy, FlowConditionRead, &
            FlowConditionGeneralRead, FlowConditionTOilImsRead, &
            FlowConditionAddToList, FlowConditionInitList, &
            FlowConditionDestroyList, &
            FlowConditionGetPtrFromList, FlowConditionUpdate, &
            FlowConditionPrint, &
            TranConditionCreate, &
            TranConditionAddToList, TranConditionInitList, &
            TranConditionDestroyList, TranConditionGetPtrFromList, &
            TranConstraintAddToList, TranConstraintInitList, &
            TranConstraintDestroyList, TranConstraintGetPtrFromList, &
            TranConditionRead, TranConstraintRead, &
            TranConditionUpdate, &
            FlowConditionIsTransient, &
            ConditionReadValues, &
            GetSubConditionName, &
            FlowConditionUnknownItype, &
            FlowCondInputRecord, &
            TranCondInputRecord
    
contains

! ************************************************************************** !

function FlowConditionCreate(option)
  ! 
  ! Creates a condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_condition_type), pointer :: FlowConditionCreate
  
  type(flow_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%saturation)
  nullify(condition%rate)
  nullify(condition%energy_rate)
  nullify(condition%energy_flux)
  nullify(condition%well)
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%sub_condition_ptr)
  nullify(condition%general)
  nullify(condition%toil_ims)
  nullify(condition%flow_well)
  nullify(condition%itype)
  nullify(condition%next)
  nullify(condition%datum)
  nullify(condition%default_time_storage)
  condition%is_transient = PETSC_FALSE
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%iphase = 0
  condition%num_sub_conditions = 0
  condition%name = ''
  
  FlowConditionCreate => condition

end function FlowConditionCreate

! ************************************************************************** !

function TranConditionCreate(option)
  ! 
  ! Creates a transport condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(tran_condition_type), pointer :: TranConditionCreate
  
  type(tran_condition_type), pointer :: condition
  
  allocate(condition)
  nullify(condition%constraint_coupler_list)
  nullify(condition%cur_constraint_coupler)
  nullify(condition%next)
  condition%id = 0
  condition%itype = 0
  condition%name = ''

  TranConditionCreate => condition

end function TranConditionCreate

! ************************************************************************** !

function FlowGeneralConditionCreate(option)
  ! 
  ! Creates a condition for general mode
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/26/11
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_general_condition_type), pointer :: FlowGeneralConditionCreate
  
  type(flow_general_condition_type), pointer :: general_condition
  
  allocate(general_condition)
  nullify(general_condition%liquid_pressure)
  nullify(general_condition%gas_pressure)
  nullify(general_condition%gas_saturation)
  nullify(general_condition%relative_humidity)
  nullify(general_condition%mole_fraction)
  nullify(general_condition%temperature)
  nullify(general_condition%liquid_flux)
  nullify(general_condition%gas_flux)
  nullify(general_condition%energy_flux)
  nullify(general_condition%rate)

  FlowGeneralConditionCreate => general_condition

end function FlowGeneralConditionCreate

! ************************************************************************** !

function FlowTOilImsConditionCreate(option)
  ! 
  ! Creates a condition for toil_ims mode
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 09/9/2015
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_toil_ims_condition_type), pointer :: FlowTOilImsConditionCreate
  
  type(flow_toil_ims_condition_type), pointer :: toil_ims_condition

  allocate(toil_ims_condition)
  nullify(toil_ims_condition%pressure)
  nullify(toil_ims_condition%saturation)
  nullify(toil_ims_condition%temperature)
  nullify(toil_ims_condition%enthalpy)
  nullify(toil_ims_condition%rate)
  nullify(toil_ims_condition%liquid_flux)
  nullify(toil_ims_condition%oil_flux)
  nullify(toil_ims_condition%energy_flux)
  nullify(toil_ims_condition%owc)
  nullify(toil_ims_condition%liq_press_grad) 

  FlowTOilImsConditionCreate => toil_ims_condition

end function FlowTOilImsConditionCreate

! ************************************************************************** !

function FlowWellConditionCreate(option)
  ! 
  ! Creates a condition for toil_ims mode
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(flow_well_condition_type), pointer :: FlowWellConditionCreate
  
  type(flow_well_condition_type), pointer :: flow_well_condition

  allocate(flow_well_condition)
  nullify(flow_well_condition%pressure)
  nullify(flow_well_condition%rate)
  nullify(flow_well_condition%temperature)

  FlowWellConditionCreate => flow_well_condition

end function FlowWellConditionCreate


! ************************************************************************** !

function FlowGeneralSubConditionPtr(sub_condition_name,general, &
                                    option)
  ! 
  ! Returns a pointer to a subcondition, creating
  ! them if necessary
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/09/11
  ! 

  use Option_module
  use Input_Aux_module, only : InputKeywordUnrecognized

  implicit none

  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_general_condition_type) :: general
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowGeneralSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('LIQUID_PRESSURE')
      if (associated(general%liquid_pressure)) then
        sub_condition_ptr => general%liquid_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%liquid_pressure => sub_condition_ptr
      endif
    case('GAS_PRESSURE')
      if (associated(general%gas_pressure)) then
        sub_condition_ptr => general%gas_pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','GAS_SATURATION')
      if (associated(general%gas_saturation)) then
        sub_condition_ptr => general%gas_saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_saturation => sub_condition_ptr
      endif
    case('TEMPERATURE')
      if (associated(general%temperature)) then
        sub_condition_ptr => general%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%temperature => sub_condition_ptr
      endif
    case('RELATIVE_HUMIDITY')
      if (associated(general%relative_humidity)) then
        sub_condition_ptr => general%relative_humidity
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%relative_humidity => sub_condition_ptr
      endif
    case('MOLE_FRACTION')
      if (associated(general%mole_fraction)) then
        sub_condition_ptr => general%mole_fraction
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%mole_fraction => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(general%liquid_flux)) then
        sub_condition_ptr => general%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%liquid_flux => sub_condition_ptr
      endif
    case('GAS_FLUX')
      if (associated(general%gas_flux)) then
        sub_condition_ptr => general%gas_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%gas_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(general%energy_flux)) then
        sub_condition_ptr => general%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        general%energy_flux => sub_condition_ptr
      endif
    case('RATE')
      if (associated(general%rate)) then
        sub_condition_ptr => general%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        general%rate => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(sub_condition_name, &
                                    'general condition,type',option)
  end select

  FlowGeneralSubConditionPtr => sub_condition_ptr

end function FlowGeneralSubConditionPtr

! ************************************************************************** !

function FlowTOilImsSubConditionPtr(sub_condition_name,toil_ims, &
                                    option)
  ! 
  ! Returns a pointer to a subcondition, creating
  ! them if necessary for toil_ims flow mode
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 09/09/2015
  ! 

  use Option_module
  use Input_Aux_module, only : InputKeywordUnrecognized

  implicit none

  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_toil_ims_condition_type) :: toil_ims
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowTOilImsSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE')
      if (associated(toil_ims%pressure)) then
        sub_condition_ptr => toil_ims%pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%pressure => sub_condition_ptr
      endif
    case('LIQUID_SATURATION','OIL_SATURATION')
      if (associated(toil_ims%saturation)) then
        sub_condition_ptr => toil_ims%saturation
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%saturation => sub_condition_ptr
      endif
    case('TEMPERATURE')
      if (associated(toil_ims%temperature)) then
        sub_condition_ptr => toil_ims%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%temperature => sub_condition_ptr
      endif
    case('ENTHALPY')
      if (associated(toil_ims%enthalpy)) then
        sub_condition_ptr => toil_ims%enthalpy
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%enthalpy => sub_condition_ptr
      endif
    case('LIQUID_FLUX')
      if (associated(toil_ims%liquid_flux)) then
        sub_condition_ptr => toil_ims%liquid_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%liquid_flux => sub_condition_ptr
      endif
    case('OIL_FLUX')
      if (associated(toil_ims%oil_flux)) then
        sub_condition_ptr => toil_ims%oil_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%oil_flux => sub_condition_ptr
      endif
    case('ENERGY_FLUX')
      if (associated(toil_ims%energy_flux)) then
        sub_condition_ptr => toil_ims%energy_flux
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%energy_flux => sub_condition_ptr
      endif
    case('OWC')
      if (associated(toil_ims%owc)) then
        sub_condition_ptr => toil_ims%owc
      else
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%owc => sub_condition_ptr
      endif
    case('WATER_PRESSURE_GRAD')
      if (associated(toil_ims%liq_press_grad)) then
        sub_condition_ptr => toil_ims%liq_press_grad
      else
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%liq_press_grad => sub_condition_ptr
      endif
    !a condition can have either RATE or WELL_RATE (target/limits)
    case('RATE')
      if (associated(toil_ims%rate)) then
        sub_condition_ptr => toil_ims%rate
      else
        ! energy rate is loaded in the third record
        sub_condition_ptr => FlowSubConditionCreate(THREE_INTEGER)
        toil_ims%rate => sub_condition_ptr
      endif
    !a condition can have either RATE or WELL_RATE (target/limits)  
    case('WELL_RATE')
      if (associated(toil_ims%rate)) then
        sub_condition_ptr => toil_ims%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        toil_ims%rate => sub_condition_ptr
      endif
    case default
      call InputKeywordUnrecognized(sub_condition_name, &
                                    'toil_ims condition,type',option)
  end select

  FlowTOilImsSubConditionPtr => sub_condition_ptr

end function FlowTOilImsSubConditionPtr

! ************************************************************************** !

function FlowWellSubConditionPtr(sub_condition_name,flow_well, &
                                    option)
  ! 
  ! Returns a pointer to a subcondition, creating
  ! them if necessary for flow_well 
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 6/03/2016
  ! 

  use Option_module
  use Input_Aux_module, only : InputKeywordUnrecognized

  implicit none

  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_well_condition_type) :: flow_well
  type(option_type) :: option

  type(flow_sub_condition_type), pointer :: FlowWellSubConditionPtr
  type(flow_sub_condition_type), pointer :: sub_condition_ptr

  select case(sub_condition_name)
    case('WELL_PRESSURE')
      if (associated(flow_well%pressure)) then
        sub_condition_ptr => flow_well%pressure
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        flow_well%pressure => sub_condition_ptr
      endif
    !a condition can have either RATE or WELL_RATE (target/limits)  
    case('WELL_RATE')
      if (associated(flow_well%rate)) then
        sub_condition_ptr => flow_well%rate
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        flow_well%rate => sub_condition_ptr
      endif
    case('WELL_TEMPERATURE')
      if (associated(flow_well%temperature)) then
        sub_condition_ptr => flow_well%temperature
      else
        sub_condition_ptr => FlowSubConditionCreate(ONE_INTEGER)
        flow_well%temperature => sub_condition_ptr
      endif

  end select

  FlowWellSubConditionPtr => sub_condition_ptr

end function FlowWellSubConditionPtr

! ************************************************************************** !

function FlowSubConditionCreate(ndof)
  ! 
  ! Creates a sub_condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  ! 

  use Dataset_Ascii_class
  use Option_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: FlowSubConditionCreate
  
  PetscInt :: ndof
  
  type(flow_sub_condition_type), pointer :: sub_condition
  class(dataset_ascii_type), pointer :: dataset_ascii
  
  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = NULL_CONDITION
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''
  sub_condition%aux_real = UNINITIALIZED_DOUBLE
  nullify(sub_condition%gradient)
  nullify(sub_condition%dataset)

  ! by default, all dataset are of type dataset_ascii_type, unless overwritten
  dataset_ascii => DatasetAsciiCreate()
  call DatasetAsciiInit(dataset_ascii)
  dataset_ascii%array_width = ndof
  dataset_ascii%data_type = DATASET_REAL
  sub_condition%dataset => dataset_ascii
  nullify(dataset_ascii)
  
  FlowSubConditionCreate => sub_condition

end function FlowSubConditionCreate

! ************************************************************************** !

function GetFlowSubCondFromArrayByName(sub_condition_ptr_list,name)
  ! 
  ! returns a pointer to a subcondition with
  ! matching name
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  ! 

  use Input_Aux_module
  use String_module
  
  implicit none
  
  type(flow_sub_condition_type), pointer :: GetFlowSubCondFromArrayByName
  type(sub_condition_ptr_type), pointer :: sub_condition_ptr_list(:)
  character(len=MAXWORDLENGTH) :: name
  
  PetscInt :: idof
  PetscInt :: length
  
  nullify(GetFlowSubCondFromArrayByName)
  length = len_trim(name)
  do idof = 1, size(sub_condition_ptr_list)
    if (length == len_trim(sub_condition_ptr_list(idof)%ptr%name) .and. &
        StringCompare(name,sub_condition_ptr_list(idof)%ptr%name,length)) then
      GetFlowSubCondFromArrayByName => sub_condition_ptr_list(idof)%ptr
      return
    endif
  enddo

  print *, 'GetFlowSubCondFromArrayByName() needs to be updated to include' // &
           'the general_condition_type.'
  stop
  
end function GetFlowSubCondFromArrayByName

! ************************************************************************** !

subroutine FlowSubConditionVerify(option, condition, sub_condition_name, &
                                  sub_condition, default_time_storage, &
                                  destroy_if_null)
  ! 
  ! Verifies the data in a subcondition
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  ! 
  use Time_Storage_module
  use Option_module
  use Dataset_module

  implicit none
  
  type(option_type) :: option
  type(flow_condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_sub_condition_type), pointer :: sub_condition
  type(time_storage_type), pointer :: default_time_storage
  PetscBool :: destroy_if_null

  character(len=MAXSTRINGLENGTH) :: header

  if (.not.associated(sub_condition)) return

  ! dataset is not optional
  if (.not.(associated(sub_condition%dataset%rarray) .or. &
            associated(sub_condition%dataset%rbuffer) .or. &
            ! if a dataset name is read, instead of data at this point
            len_trim(sub_condition%dataset%name) > 0 .or. &
            sub_condition%itype /= NULL_CONDITION)) then
    if (destroy_if_null) call FlowSubConditionDestroy(sub_condition)
    return
  endif
  
!  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
  if (sub_condition%itype == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call printErrMsg(option)
  endif
  
  header = 'SUBSURFACE/FLOW_CONDITION' // '/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Value(s)'
  call DatasetVerify(sub_condition%dataset,default_time_storage, &
                     header,option)
  header = 'SUBSURFACE/FLOW_CONDITION' // '/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Gradient'
  call DatasetVerify(sub_condition%gradient,default_time_storage, &
                     header,option)

end subroutine FlowSubConditionVerify

! ************************************************************************** !

subroutine FlowConditionRead(condition,input,option)
  ! 
  ! Reads a condition from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: rate_string
  character(len=MAXWORDLENGTH) :: internal_units
  type(flow_sub_condition_type), pointer :: pressure, flux, temperature, &
                                       concentration, enthalpy, rate, well,&
                                       sub_condition_ptr, saturation, &
                                       energy_rate, energy_flux
  PetscReal :: default_time
  PetscInt :: default_iphase
  PetscInt :: idof
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii  
  
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  default_time = 0.d0
  default_iphase = 0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'
  
  pressure => FlowSubConditionCreate(option%nphase)
  pressure%name = 'pressure'
  flux => pressure
  rate => FlowSubConditionCreate(option%nflowspec)
  rate%name = 'rate'
  energy_rate => FlowSubConditionCreate(ONE_INTEGER)
  energy_rate%name = 'energy_rate'
  energy_flux => FlowSubConditionCreate(ONE_INTEGER)
  energy_flux%name = 'energy_flux'
  well => FlowSubConditionCreate(7 + option%nflowspec)
  well%name = 'well'
  saturation => FlowSubConditionCreate(option%nphase)
  saturation%name = 'saturation'
  temperature => FlowSubConditionCreate(ONE_INTEGER)
  temperature%name = 'temperature'
  concentration => FlowSubConditionCreate(ONE_INTEGER)
  concentration%name = 'concentration'
  enthalpy => FlowSubConditionCreate(option%nphase)
  enthalpy%name = 'enthalpy'

  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  rate%units = 'kg/s'
  energy_rate%units = 'W'
  energy_flux%units = 'W/m^2'
  well%units = 'Pa'
  saturation%units = ' '
  temperature%units = 'C'
  concentration%units = 'M'
  enthalpy%units = 'kJ/mol'

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
            case('Pa','KPa')
              pressure%units = trim(word)
            case('kg/s','kg/yr')
              rate%units = trim(word)
            case('W','J/yr')
              energy_rate%units = trim(word)
            case('W/m^2','J/m^2/yr')
              energy_flux%units = trim(word)
            case('m/s','m/yr')
              flux%units = trim(word)
            case('C','K')
              temperature%units = trim(word)
            case('M','mol/L')
              concentration%units = trim(word)
            case('kJ/mol')
              enthalpy%units = trim(word)
            case default
              call InputKeywordUnrecognized(word,'condition,units',option)
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
          case default
            call InputKeywordUnrecognized(word,'condition,interpolation', &
                                          option)
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
              sub_condition_ptr => pressure
              internal_units = 'Pa'
            case('RATE')
              sub_condition_ptr => rate
              internal_units = 'unitless/sec'
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
              internal_units = 'MJ/sec|MW'
            case('WELL')
              sub_condition_ptr => well
              internal_units = 'Pa'
            case('FLUX')
              sub_condition_ptr => flux
              internal_units = 'meter/sec'
            case('ENERGY_FLUX')
              sub_condition_ptr => energy_flux
              internal_units = 'MW/m^2|MJ/sec-m^2'
            case('SATURATION')
              sub_condition_ptr => saturation
              internal_units = 'unitless'
            case('TEMPERATURE')
              sub_condition_ptr => temperature
              internal_units = 'C'
            case('CONCENTRATION')
              sub_condition_ptr => concentration
              internal_units = 'unitless'
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
              internal_units = 'MJ/mol'
            case default
              call InputKeywordUnrecognized(word,'condition,type',option)
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
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'
            case('energy_rate')
              sub_condition_ptr%itype = ENERGY_RATE_SS
              rate_string = 'MJ/sec|MW'
            case('heterogeneous_energy_rate')
              sub_condition_ptr%itype = HET_ENERGY_RATE_SS
              rate_string = 'MJ/sec|MW'
            case('scaled_mass_rate','scaled_volumetric_rate', &
                 'scaled_energy_rate')
              select case(word)
                case('scaled_mass_rate')
                  sub_condition_ptr%itype = SCALED_MASS_RATE_SS
                  rate_string = 'kg/sec'
                case('scaled_volumetric_rate')
                  sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
                  rate_string = 'm^3/sec'
                case('scaled_energy_rate')
                  sub_condition_ptr%itype = SCALED_ENERGY_RATE_SS
                  rate_string = 'MW|MJ/sec'
              end select
              ! store name of type for error messaging below.
              string = word
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" ' // trim(string)
                    call InputKeywordUnrecognized(word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, ' // &
                  'VOLUME, PERM subtypes in '// &
                  'flow condition "' // trim(condition%name) // &
                  '" ' // trim(string)
                call printErrMsg(option)              
                endif
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('well','production_well', 'injection_well')
              sub_condition_ptr%itype = WELL_SS
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('unit_gradient')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may ' // &
                  'only be associated with a PRESSURE flow condition.'
                call printErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC 
            case('heterogeneous_volumetric_rate')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'
            case('heterogeneous_mass_rate')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'
            case('heterogeneous_dirichlet')
              sub_condition_ptr%itype = HET_DIRICHLET
            case('heterogeneous_surface_seepage')
              sub_condition_ptr%itype = HET_SURF_SEEPAGE_BC
            case('spillover')
              sub_condition_ptr%itype = SPILLOVER_BC
            case('surface_dirichlet')
              sub_condition_ptr%itype = SURFACE_DIRICHLET
            case('surface_zero_gradheight')
              sub_condition_ptr%itype = SURFACE_ZERO_GRADHEIGHT
            case('surface_spillover')
              sub_condition_ptr%itype = SURFACE_SPILLOVER
            case default
              call InputKeywordUnrecognized(word,'condition bc type',option)
          end select
        enddo
      case('TIME','TIMES')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION')   
      case('IPHASE')
        call InputReadInt(input,option,default_iphase)
        call InputErrorMsg(input,option,'IPHASE','CONDITION')   
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii) 
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 condition%datum,word,internal_units)
      case('GRADIENT','GRAD')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
              internal_units = 'Pa/meter'
            case('RATE')
              sub_condition_ptr => rate
              internal_units = 'kg/sec-meter'
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
              internal_units = 'MW/meter|MJ/sec-meter'
            case('WELL')
              sub_condition_ptr => well
              internal_units = 'Pa/meter'
            case('FLUX')
              sub_condition_ptr => flux
              internal_units = 'm/sec-m|unitless/sec'
            case('SATURATION')
              sub_condition_ptr => saturation
              internal_units = 'unitless/meter'
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
              internal_units = 'temperature/m'
            case('CONC','CONCENTRATION')
              sub_condition_ptr => concentration
              internal_units = 'unitless'
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
              internal_units = 'kJ/mol-meter'
            case default
              call InputKeywordUnrecognized(word, &
                     'FLOW CONDITION,GRADIENT,TYPE',option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
      case('TEMPERATURE','TEMP')
        internal_units = 'C'
        call ConditionReadValues(input,option,word, &
                                 temperature%dataset, &
                                 temperature%units,internal_units)
      case('ENTHALPY','H')
        internal_units = 'kJ/mol'
        call ConditionReadValues(input,option,word, &
                                 enthalpy%dataset, &
                                 enthalpy%units,internal_units)
      case('PRESSURE','PRES','PRESS')
        internal_units = 'Pa'
        call ConditionReadValues(input,option,word, &
                                 pressure%dataset, &
                                 pressure%units,internal_units)
      case('RATE')
        internal_units = rate_string
        call ConditionReadValues(input,option,word, &
                                 rate%dataset, &
                                 rate%units,internal_units)
      case('ENERGY_FLUX')
        input%force_units = PETSC_TRUE
        internal_units = 'MW/m^2|MJ/m^2-sec'
        call ConditionReadValues(input,option,word, &
                                 energy_flux%dataset, &
                                 energy_flux%units,internal_units)
        input%force_units = PETSC_FALSE
      case('ENERGY_RATE')
        input%force_units = PETSC_TRUE
        internal_units = 'MJ/sec|MW'
        input%err_buf = word
        call ConditionReadValues(input,option,word, &
                                 energy_rate%dataset, &
                                 energy_rate%units,internal_units)
        input%force_units = PETSC_FALSE
      case('WELL')
        internal_units = 'Pa'
        call ConditionReadValues(input,option,word, &
                                 well%dataset, &
                                 well%units,internal_units)
      case('FLUX','VELOCITY','VEL')
        internal_units = 'meter/sec'
        call ConditionReadValues(input,option,word, &
                                 pressure%dataset, &
                                 pressure%units,internal_units)
      case('CONC','CONCENTRATION')
        internal_units = 'unitless'
        call ConditionReadValues(input,option,word, &
                                 concentration%dataset, &
                                 concentration%units,internal_units)
      case('SAT','SATURATION')
        internal_units = 'unitless'
        call ConditionReadValues(input,option,word, &
                                 saturation%dataset, &
                                 saturation%units,internal_units)
      case('CONDUCTANCE')
        call InputReadDouble(input,option,pressure%aux_real(1))
        call InputErrorMsg(input,option,'CONDUCTANCE','CONDITION')   
      case default
        call InputKeywordUnrecognized(word,'flow condition',option)
    end select 
  
  enddo  
  
  ! check whether
  if (default_iphase == 0) then
    option%io_buffer = '"iphase" not set in condition; set to 1'
    call printWrnMsg(option)
    condition%iphase = 1
  else
    condition%iphase = default_iphase    
  endif
  
  ! datum is not required
  string = trim(condition%name) // '/' // 'Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! check to ensure that a rate condition is not of type pressure   
  if (associated(rate)) then
    select case(rate%itype)
      case(DIRICHLET_BC,NEUMANN_BC,HYDROSTATIC_BC,UNIT_GRADIENT_BC, &
           CONDUCTANCE_BC,ZERO_GRADIENT_BC,SEEPAGE_BC,SURFACE_DIRICHLET, &
           SURFACE_SPILLOVER)
        option%io_buffer = 'RATE condition must not be of type: dirichlet, ' // &
          'neumann, zero_gradient, dirichlet_zero_gradient, hydrostatic, ' // &
          'seepage, or conductance".'
        call printErrMsg(option)
    end select
  endif
  ! check to ensure that a pressure condition is not of type rate   
  if (associated(pressure)) then                          
    select case(pressure%itype)
      case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS)
        option%io_buffer = 'PRESSURE or FLUX condition must not be of type: ' // &
          'mass_rate, scaled_mass_rate, volumetric_rate, ' // &
          'scaled_volumetric_rate, equilibrium, or production_well.'
        call printErrMsg(option)
    end select
  endif

  ! verify the datasets
  word = 'pressure/flux'
  call FlowSubConditionVerify(option,condition,word,pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_flux'
  call FlowSubConditionVerify(option,condition,word,energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_rate'
  call FlowSubConditionVerify(option,condition,word,energy_rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'well'
  call FlowSubConditionVerify(option,condition,word,well, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'saturation'
  call FlowSubConditionVerify(option,condition,word,saturation, &
                              default_time_storage, &
                              PETSC_TRUE)

  word = 'concentration'
  call FlowSubConditionVerify(option,condition,word,concentration, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)

  select case(option%iflowmode)
    case(G_MODE)
      option%io_buffer = 'General mode not supported in original FlowConditionRead.'
      call printMsg(option)
    case(TOIL_IMS_MODE)
      option%io_buffer = 'TOilIms mode not supported in original FlowConditionRead.'
      call printMsg(option)
    case(MPH_MODE,IMS_MODE,FLASH2_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
     
      
      if (.not.associated(temperature) .and. .not.associated(energy_rate)) then
        option%io_buffer = 'temperature and energy rate condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      if (associated(temperature)) then
        condition%temperature => temperature
      endif
      if (associated(energy_flux)) then
        condition%energy_flux => energy_flux
      endif
      if (associated(energy_rate)) then
        condition%energy_rate => energy_rate
      endif
      
      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%concentration => concentration
      
      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%enthalpy => enthalpy
      
      condition%num_sub_conditions = 4
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 4
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(THREE_INTEGER)%ptr => concentration
      if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy
      if (associated(energy_rate)) &
        condition%sub_condition_ptr(FOUR_INTEGER)%ptr => energy_rate
        
      allocate(condition%itype(FIVE_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      condition%itype(TWO_INTEGER) = temperature%itype
      condition%itype(THREE_INTEGER) = concentration%itype
      if (associated(enthalpy)) condition%itype(FOUR_INTEGER) = concentration%itype
      if (associated(energy_rate)) condition%itype(FOUR_INTEGER) = energy_rate%itype

    case(TH_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif

      if (.not.associated(temperature) .and. .not.associated(energy_rate) &
          .and. .not.associated(energy_flux)) then
        option%io_buffer = 'temperature, energy_flux, and energy_rate ' // &
          'condition null in condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      if (associated(temperature) .and. associated(energy_rate) ) then
        option%io_buffer = 'Both, temperature and energy_rate cannot be ' // &
                            'specified in condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      if (associated(temperature)) condition%temperature => temperature
      if (associated(energy_flux)) condition%energy_flux => energy_flux
      if (associated(energy_rate)) condition%energy_rate => energy_rate

      if (associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition not supported in TH mode: ' // &
                            trim(condition%name)
        call printErrMsg(option)
      endif
      if (associated(enthalpy)) condition%enthalpy => enthalpy
      
      condition%num_sub_conditions = TWO_INTEGER
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 2
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      if (associated(saturation)) condition%sub_condition_ptr(ONE_INTEGER)%ptr &
                                  => saturation
      if ( associated(temperature)) &
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      if (associated(energy_flux)) condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_flux
      if (associated(energy_rate)) condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_rate

      allocate(condition%itype(TWO_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      if (associated(saturation)) condition%itype(ONE_INTEGER) = &
                                    saturation%itype
      if (associated(temperature)) condition%itype(TWO_INTEGER) = temperature%itype
      if (associated(energy_flux)) condition%itype(TWO_INTEGER) = energy_flux%itype
      if (associated(energy_rate)) condition%itype(TWO_INTEGER) = energy_rate%itype

!#if 0
    case(MIS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well)) then
        option%io_buffer = 'pressure and rate condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
      
      if (.not.associated(concentration)) then
        option%io_buffer = 'concentration condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%concentration => concentration

#if 0
      if (.not.associated(temperature)) then
        option%io_buffer = 'temperature condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%temperature => temperature
      
      if (.not.associated(enthalpy)) then
        option%io_buffer = 'enthalpy condition null in condition: ' // &
                            trim(condition%name)      
        call printErrMsg(option)
      endif                         
      condition%enthalpy => enthalpy
#endif

      condition%num_sub_conditions = 2
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, 2
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs in problem
      if (associated(pressure)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(rate)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(well)) condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
!     condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      condition%sub_condition_ptr(TWO_INTEGER)%ptr => concentration
!     if (associated(enthalpy)) condition%sub_condition_ptr(FOUR_INTEGER)%ptr => enthalpy

      allocate(condition%itype(TWO_INTEGER))
      condition%itype = 0
      if (associated(pressure)) condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(rate)) condition%itype(ONE_INTEGER) = rate%itype
      if (associated(well)) condition%itype(ONE_INTEGER) = well%itype
      condition%itype(TWO_INTEGER) = concentration%itype
!#endif
    
    case(RICHARDS_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate) .and. &
          .not.associated(saturation) .and. .not.associated(well)) then
        option%io_buffer = 'pressure, rate and saturation condition null in ' // &
                           'condition: ' // trim(condition%name)
        call printErrMsg(option)      
      endif
      
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(well)) then
        condition%well => well
      endif
            
      condition%num_sub_conditions = 1
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      if (associated(pressure)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      elseif (associated(saturation)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => saturation
      elseif (associated(rate)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      elseif (associated(well)) then
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => well
      endif                         

      allocate(condition%itype(ONE_INTEGER))
      if (associated(pressure)) then 
        condition%itype(ONE_INTEGER) = pressure%itype
      else if (associated(saturation)) then
        condition%itype(ONE_INTEGER) = saturation%itype
      else if (associated(rate)) then
        condition%itype(ONE_INTEGER) = rate%itype
      else if (associated(well)) then
        condition%itype(ONE_INTEGER) = well%itype
      endif
      
      ! these are not used with richards
      if (associated(temperature)) call FlowSubConditionDestroy(temperature)
      if (associated(enthalpy)) call FlowSubConditionDestroy(enthalpy)

  end select

  condition%default_time_storage => default_time_storage
  
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionRead

! ************************************************************************** !

subroutine FlowConditionGeneralRead(condition,input,option)
  ! 
  ! Reads a condition from the input file for
  ! general mode
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/14/11
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module
  
  use General_Aux_module
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string, internal_units
  character(len=MAXWORDLENGTH) :: word
  type(flow_general_condition_type), pointer :: general
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscBool :: default_is_cyclic
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii  
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'

  default_time = 0.d0
  default_iphase = 0
  
  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP
  
  select case(option%iflowmode)
    case(G_MODE)
      general => FlowGeneralConditionCreate(option)
      condition%general => general
  end select
  
  ! read the condition
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
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
          select case(option%iflowmode)
            case(G_MODE)
              sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                              option)
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
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_string = 'kg/sec'                                
            case('scaled_mass_rate')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'                                
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_mass_rate type'
                    call InputKeywordUnrecognized(word,string,option)
                end select
              else      
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, ' // &
                  'VOLUME, PERM subtypes in '// &
                  'flow condition "' // trim(condition%name) // &
                  '" scaled_mass_rate type'
                call printErrMsg(option)
              endif
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'                                  
            case('scaled_volumetric_rate')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'                                   
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, ' // &
                  'VOLUME, PERM subtypes in '// &
                  'flow condition "' // trim(condition%name) // &
                  '" scaled_volumetric_rate type'
                call printErrMsg(option)
              endif
            case('heterogeneous_volumetric_rate')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'                                 
            case('heterogeneous_mass_rate')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'                                 
            case('heterogeneous_dirichlet')
              sub_condition_ptr%itype = HET_DIRICHLET
            case('heterogeneous_surface_seepage')
              sub_condition_ptr%itype = HET_SURF_SEEPAGE_BC
            case default
              call InputKeywordUnrecognized(word,'flow condition,type',option)
          end select
        enddo
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)        
        internal_units = 'meter'
        call ConditionReadValues(input,option,word,condition%datum, &
                                 word,internal_units)
      case('GRADIENT')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(option%iflowmode)
            case(G_MODE)
              sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                              option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units = 'unitless/meter'
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
      case('CONDUCTANCE')
        word = 'LIQUID_PRESSURE'
        select case(option%iflowmode)
          case(G_MODE)
            sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                            option)
        end select
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')   
      case('LIQUID_PRESSURE','GAS_PRESSURE','LIQUID_SATURATION', &
           'GAS_SATURATION','TEMPERATURE','MOLE_FRACTION','RATE', &
           'LIQUID_FLUX','GAS_FLUX','ENERGY_FLUX','RELATIVE_HUMIDITY')
        select case(option%iflowmode)
          case(G_MODE)
            sub_condition_ptr => FlowGeneralSubConditionPtr(word,general, &
                                                            option)
        end select
        select case(trim(word))
          case('LIQUID_PRESSURE','GAS_PRESSURE')
            internal_units = 'Pa'
          case('LIQUID_SATURATION','GAS_SATURATION','MOLE_FRACTION', &
               'RELATIVE_HUMIDITY')
            internal_units = 'unitless'
          case('TEMPERATURE')
            internal_units = 'C'
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = trim(rate_string) // ',' // trim(rate_string) //&
                             ',MJ/sec|MW'
          case('LIQUID_FLUX','GAS_FLUX')
            internal_units = 'meter/sec'
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = 'MW/m^2|MJ/m^2-sec'
        end select
        call ConditionReadValues(input,option,word, &
                                 sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units,internal_units)
        input%force_units = PETSC_FALSE
        select case(word)
          case('LIQUID_SATURATION') ! convert to gas saturation
            if (associated(sub_condition_ptr%dataset%rbuffer)) then
              sub_condition_ptr%dataset%rbuffer(:) = 1.d0 - &
                sub_condition_ptr%dataset%rbuffer(:)
            endif
            sub_condition_ptr%dataset%rarray(:) = 1.d0 - &
              sub_condition_ptr%dataset%rarray(:)
        end select
      case default
        call InputKeywordUnrecognized(word,'flow condition',option)
    end select 
  
  enddo  
  
  ! datum is not required
  string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name) // '/Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! need mole fraction and some sort of saturation
  if (.false.) then
    ! neumann or mass/volumetric flux
    ! need temperature
  !  condition%sub_condition_ptr(GENERAL_FLUX_DOF)%ptr => general%flux
    if (.not.associated(general%mole_fraction) .and. &
        .not.associated(general%gas_saturation)) then
      option%io_buffer = 'General Phase flux condition must include ' // &
        'a MOLE_FRACTION or GAS/LIQUID_SATURATION.'
      call printErrMsg(option)
    endif
    if (associated(general%mole_fraction) .and. &
        associated(general%gas_saturation)) then
      option%io_buffer = 'General Phase flux condition must include ' // &
        'only a MOLE_FRACTION or GAS/LIQUID_SATURATION, not both.'
      call printErrMsg(option)
    endif
    if (.not.associated(general%temperature)) then
      option%io_buffer = 'General Phase flux condition must include ' // &
        'a temperature'
      call printErrMsg(option)
    endif
  else
    if (associated(general%rate)) then
      condition%iphase = ANY_STATE
    elseif (associated(general%liquid_flux) .and. &
            associated(general%gas_flux) .and. &
            (associated(general%energy_flux) .or. &
             associated(general%temperature))) then
      condition%iphase = ANY_STATE
    else
      ! some sort of dirichlet-based pressure, temperature, etc.
      if (.not.associated(general%liquid_pressure) .and. &
          .not.associated(general%gas_pressure)) then
        option%io_buffer = 'General Phase non-rate condition must ' // &
          'include a liquid or gas pressure'
        call printErrMsg(option)
      endif
      if (.not.associated(general%mole_fraction) .and. &
          .not.associated(general%relative_humidity) .and. &
          .not.associated(general%gas_saturation)) then
        option%io_buffer = 'General Phase non-rate condition must ' // &
          'include mole fraction, relative humidity, or gas/liquid saturation'
        call printErrMsg(option)
      endif
      if (.not.associated(general%temperature)) then
        option%io_buffer = 'General Phase non-rate condition must ' // &
          'include temperature'
        call printErrMsg(option)
      endif
      if (associated(general%gas_pressure) .and. &
          associated(general%gas_saturation)) then
        ! two phase condition
        condition%iphase = TWO_PHASE_STATE
      else if (associated(general%liquid_pressure) .and. &
               associated(general%mole_fraction)) then
        ! liquid phase condition
        condition%iphase = LIQUID_STATE
      else if (associated(general%gas_pressure) .and. &
               (associated(general%mole_fraction) .or. &
                associated(general%relative_humidity))) then
        ! gas phase condition
        condition%iphase = GAS_STATE
      endif
    endif
    if (condition%iphase == NULL_STATE) then
      option%io_buffer = 'General Phase non-rate/flux condition contains ' // &
        'an unsupported combination of primary dependent variables.'
      call printErrMsg(option)
    endif
  endif
    
  ! verify the datasets
  word = 'liquid pressure'
  call FlowSubConditionVerify(option,condition,word,general%liquid_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas pressure'
  call FlowSubConditionVerify(option,condition,word,general%gas_pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas saturation'
  call FlowSubConditionVerify(option,condition,word,general%gas_saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'relative humidity'
  call FlowSubConditionVerify(option,condition,word,general%relative_humidity, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'mole fraction'
  call FlowSubConditionVerify(option,condition,word,general%mole_fraction, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,general%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,general%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'gas flux'
  call FlowSubConditionVerify(option,condition,word,general%gas_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,general%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,general%rate, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(general%liquid_pressure)) &
    i = i + 1
  if (associated(general%gas_pressure)) &
    i = i + 1
  if (associated(general%gas_saturation)) &
    i = i + 1
  if (associated(general%relative_humidity)) &
    i = i + 1
  if (associated(general%mole_fraction)) &
    i = i + 1
  if (associated(general%temperature)) &
    i = i + 1
  if (associated(general%liquid_flux)) &
    i = i + 1
  if (associated(general%gas_flux)) &
    i = i + 1
  if (associated(general%energy_flux)) &
    i = i + 1
  if (associated(general%rate)) &
    i = i + 1
  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(general%liquid_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%liquid_pressure
  endif
  if (associated(general%gas_pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_pressure
  endif  
  if (associated(general%gas_saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_saturation
  endif  
  if (associated(general%relative_humidity)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%relative_humidity
  endif  
  if (associated(general%mole_fraction)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%mole_fraction
  endif  
  if (associated(general%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%temperature
  endif  
  if (associated(general%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%liquid_flux
  endif  
  if (associated(general%gas_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%gas_flux
  endif  
  if (associated(general%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%energy_flux
  endif  
  if (associated(general%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => general%rate
  endif

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo
  
  condition%default_time_storage => default_time_storage
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionGeneralRead

! ************************************************************************** !

subroutine FlowConditionTOilImsRead(condition,input,option)
  ! 
  ! Reads a condition from the input file for
  ! toil_ims mode
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 9/9/2015
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module
  
  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module
  
  implicit none
  
  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: rate_string, internal_units
  character(len=MAXWORDLENGTH) :: word
  type(flow_toil_ims_condition_type), pointer :: toil_ims
  type(flow_well_condition_type), pointer :: flow_well
  type(flow_sub_condition_type), pointer :: sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: default_iphase
  class(dataset_base_type), pointer :: default_flow_dataset
  class(dataset_base_type), pointer :: default_gradient
  PetscInt :: idof, i
  PetscBool :: default_is_cyclic
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii  
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  rate_string = 'not_assigned'
  internal_units = 'not_assigned'

  default_time = 0.d0
  default_iphase = 0
  
  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP
  
  !select case(option%iflowmode)
  !  ! do we really need this select case??
  !  case(TOIL_IMS_MODE)
  !    toil_ims => FlowTOilImsConditionCreate(option)
  !    condition%toil_ims => toil_ims
  !end select

  toil_ims => FlowTOilImsConditionCreate(option)
  condition%toil_ims => toil_ims

  flow_well => FlowWellConditionCreate(option)  
  condition%flow_well => flow_well 

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
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
          !select case(option%iflowmode)
          !  ! do we need this select case?
          !  case(TOIL_IMS_MODE)
          !    sub_condition_ptr => FlowTOilImsSubConditionPtr(word,toil_ims, &
          !                                                    option)
          !end select
          select case(word)
            case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE', &
                 'LIQUID_SATURATION', 'OIL_SATURATION','TEMPERATURE','RATE', &
                 'LIQUID_FLUX','OIL_FLUX', 'ENERGY_FLUX','ENTHALPY','OWC', &
                 'WATER_PRESSURE_GRAD')

              sub_condition_ptr => FlowTOilImsSubConditionPtr(word,toil_ims, &
                                                              option)
            case('WELL_PRESSURE','WELL_RATE','WELL_TEMPERATURE')
              sub_condition_ptr => FlowWellSubConditionPtr(word,flow_well, &
                                                           option)
            case default
              call InputKeywordUnrecognized(word,'flow condition',option)                                          
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
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS 
              rate_string = 'kg/sec'                               
            case('scaled_mass_rate')
              sub_condition_ptr%itype = SCALED_MASS_RATE_SS
              rate_string = 'kg/sec'                                  
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = &
                      trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                case('perm')
                  sub_condition_ptr%isubtype = SCALE_BY_PERM
                case default
                  string = 'flow condition "' // trim(condition%name) // &
                    '" scaled_mass_rate type'
                  call InputKeywordUnrecognized(word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, ' // &
                  'VOLUME, PERM subtypes in '// &
                  'flow condition "' // trim(condition%name) // &
                  '" scaled_mass_rate type'
                call printErrMsg(option)
              endif
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'                                
            case('scaled_volumetric_rate')
              sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
              rate_string = 'm^3/sec'                                    
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" scaled_volumetric_rate type'
                    call InputKeywordUnrecognized(word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, ' // &
                  'VOLUME, PERM subtypes in '// &
                  'flow condition "' // trim(condition%name) // &
                  '" scaled_volumetric_rate type'
                call printErrMsg(option)
              endif
            case('heterogeneous_volumetric_rate')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_string = 'm^3/sec'                                      
            case('heterogeneous_mass_rate')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_string = 'kg/sec'                                         
            case('heterogeneous_dirichlet')
              sub_condition_ptr%itype = HET_DIRICHLET
            case('heterogeneous_surface_seepage')
              sub_condition_ptr%itype = HET_SURF_SEEPAGE_BC
            case('mass_rate_target')
              sub_condition_ptr%itype = WELL_MASS_RATE_TARGET 
              rate_string = 'kg/sec'                               
            case('mass_rate_max')
              sub_condition_ptr%itype = WELL_MASS_RATE_MAX
              rate_string = 'kg/sec'               
            case('mass_rate_min')
              sub_condition_ptr%itype = WELL_MASS_RATE_MIN
              rate_string = 'kg/sec'               
            case('vol_rate_target')
              sub_condition_ptr%itype = WELL_VOL_RATE_TARGET
              rate_string = 'm^3/sec'               
            case('vol_rate_max')
              sub_condition_ptr%itype = WELL_VOL_RATE_MAX
              rate_string = 'm^3/sec'               
            case('vol_rate_min')
              sub_condition_ptr%itype = WELL_VOL_RATE_MIN
              rate_string = 'm^3/sec'               
            case('bhp')
              sub_condition_ptr%itype = WELL_BHP
            case('bhp_min')
              sub_condition_ptr%itype = WELL_BHP_MIN
            case('bhp_max')
              sub_condition_ptr%itype = WELL_BHP_MAX                
            case default
              call InputKeywordUnrecognized(word,'flow condition,type',option)
          end select
        enddo

      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)        
        internal_units = 'meter'
        call ConditionReadValues(input,option,word,condition%datum, &
                                 word,internal_units)
      case('GRADIENT')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')
          
          if (InputCheckExit(input,option)) exit          
          
          if (InputError(input)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')   
          call StringToUpper(word)
          select case(option%iflowmode)
            case(TOIL_IMS_MODE)
              sub_condition_ptr => &
                 FlowTOilImsSubConditionPtr(word,toil_ims,option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          internal_units = 'unitless/meter'
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word, internal_units)
          nullify(sub_condition_ptr)
        enddo
      case('CONDUCTANCE')
        word = 'PRESSURE'
        select case(option%iflowmode)
          case(TOIL_IMS_MODE)
            sub_condition_ptr => FlowTOilImsSubConditionPtr(word,toil_ims, &
                                                            option)
        end select
        call InputReadDouble(input,option,sub_condition_ptr%aux_real(1))
        call InputErrorMsg(input,option,'LIQUID_CONDUCTANCE','CONDITION')   
      case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE','LIQUID_SATURATION', &
           'OIL_SATURATION','TEMPERATURE','RATE', 'LIQUID_FLUX','OIL_FLUX', &
           'ENERGY_FLUX','ENTHALPY','OWC','WATER_PRESSURE_GRAD','WELL_RATE', &
           'WELL_PRESSURE','WELL_TEMPERATURE')
        !select case(option%iflowmode)
        !  case(TOIL_IMS_MODE)
        !    sub_condition_ptr => FlowTOilImsSubConditionPtr(word,toil_ims, &
        !                                                    option)
        !end select
        select case(word)
          case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE', &
               'LIQUID_SATURATION', 'OIL_SATURATION','TEMPERATURE','RATE', &
               'LIQUID_FLUX','OIL_FLUX', 'ENERGY_FLUX','ENTHALPY','OWC', &
                'WATER_PRESSURE_GRAD')
            sub_condition_ptr => FlowTOilImsSubConditionPtr(word,toil_ims, &
                                                            option)
          case('WELL_RATE','WELL_PRESSURE','WELL_TEMPERATURE')
            sub_condition_ptr => FlowWellSubConditionPtr(word,flow_well, &
                                                         option)
        end select
        
        select case(trim(word))
        !give a type to pass FlowSubConditionVerify.
          case('OWC','WATER_PRESSURE_GRAD')
            sub_condition_ptr%itype = DIRICHLET_BC 
            sub_condition_ptr%ctype = 'dirichlet'
        end select

        select case(trim(word))
          case('PRESSURE','OIL_PRESSURE','WATER_PRESSURE')
            internal_units = 'Pa'
          case('LIQUID_SATURATION','OIL_SATURATION')
            internal_units = 'unitless'
          case('TEMPERATURE')
            internal_units = 'C'
          case('OWC')
            internal_units = 'meter'
          case('WATER_PRESSURE_GRAD')
            internal_units = 'Pa/meter'
          case('RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = trim(rate_string) // ',' // trim(rate_string) //&
                             ',MJ/sec|MW'
          case('LIQUID_FLUX','OIL_FLUX')
            internal_units = 'meter/sec'
          case('ENERGY_FLUX')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = 'MW/m^2|MJ/sec-m^2'
          case('ENTHALPY')
            internal_units = 'MJ/mol'
          case('WELL_RATE')
            input%force_units = PETSC_TRUE
            input%err_buf = word
            internal_units = trim(rate_string)
          case('WELL_PRESSURE')
            internal_units = 'Pa'
          case('WELL_TEMPERATURE')
            internal_units = 'C'
        end select
        call ConditionReadValues(input,option,word, &
                                 sub_condition_ptr%dataset, &
                                 sub_condition_ptr%units,internal_units)
        input%force_units = PETSC_FALSE
        select case(word)
          case('LIQUID_SATURATION') ! convert to oil saturation
            if (associated(sub_condition_ptr%dataset%rbuffer)) then
              sub_condition_ptr%dataset%rbuffer(:) = 1.d0 - &
                sub_condition_ptr%dataset%rbuffer(:)
            endif
            sub_condition_ptr%dataset%rarray(:) = 1.d0 - &
              sub_condition_ptr%dataset%rarray(:)
        end select
      case default
        call InputKeywordUnrecognized(word,'flow condition',option)
    end select 
  
  enddo  
  
  ! datum, owc, and liq_press_grad are not required
  string = 'SUBSURFACE/FLOW_CONDITION' // trim(condition%name) // '/Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! phase condition should never be used in TOilIms
  condition%iphase = ZERO_INTEGER

  ! unless the coondtion is a rate or pressure bhp (i.e. a bhp controlled well)
  ! - pressure is required
  ! - water or oil saturation is required
  ! - temperature is required
  !if (.not.associated(toil_ims%rate)) then
  if ( (.not.associated(toil_ims%rate)) .and. &
       ( .not.( associated(flow_well%rate).or. &
                associated(flow_well%pressure) ) ) & 
     ) then
    ! this branch is executed for sub_conditions that are not a rate or well
    ! some sort of dirichlet-based pressure, temperature, etc.
    if (.not.associated(toil_ims%pressure)) then
      option%io_buffer = 'TOilIms Phase non-rate condition must ' // &
        'include a pressure'
      call printErrMsg(option)
    endif

    if (.not.associated(toil_ims%saturation) ) then
      option%io_buffer = 'TOilIms Phase non-rate condition must ' // &
        'include liquid or oil saturation'
      call printErrMsg(option)
    endif
    if (.not.associated(toil_ims%temperature)) then
      option%io_buffer = 'TOilIms Phase non-rate condition must ' // &
        'include temperature'
      call printErrMsg(option)
    endif
  endif

  ! control that enthalpy is used for src/sink only
  if ( (.not.associated(toil_ims%rate)) .and. &
        associated(toil_ims%enthalpy)  ) then
      option%io_buffer = 'TOilIms Enthlapy condition is not' // &
       'currently supported for boundary & initial conditions'
      call printErrMsg(option)          
  end if
  ! within a src/sink either temp or enthalpy can be defined   
  if (associated(toil_ims%rate)) then
    if ( associated(toil_ims%temperature).and. &
        associated(toil_ims%enthalpy) &
       ) then
      option%io_buffer = 'TOilIms Rate condition can ' // &
       'have either temp or enthalpy'
      call printErrMsg(option)      
    end if
    ! only dirich condition supported for src/sink temp or enthalpy
    if ( ( associated(toil_ims%temperature).and. &
          (toil_ims%temperature%itype /= DIRICHLET_BC) &
         ) .or. &
         ( associated(toil_ims%enthalpy).and. &
          (toil_ims%enthalpy%itype /= DIRICHLET_BC ) &
         ) &
       ) then
      option%io_buffer = 'TOilIms Src/Sink; only dirichlet type ' // &
       'is supported for temperature and enthalpy conditions'
      call printErrMsg(option)      
    end if

    ! in the casew below enthalpy or temperature overwrite energy rate
    !if (  ( associated(toil_ims%temperature).or. &
    !       associated(toil_ims%enthalpy) &
    !     ) .and. &
    !     ( size(toil_ims%rate%dataset%rarray) == THREE_INTEGER ) &
    !   ) then 
    !  option%io_buffer = 'TOilIms Src/Sink error: ' // &
    !   'either define enery rate or temperature/enthalpy value'
    !  call printErrMsg(option)      
    !end if
  end if ! end if rate

  !TODO- in case of well_rate and/or well_pressure, well_temp must be present  
  !     


  ! verify the datasets
  word = 'pressure'
  call FlowSubConditionVerify(option,condition,word,toil_ims%pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil saturation'
  call FlowSubConditionVerify(option,condition,word,toil_ims%saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,toil_ims%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,toil_ims%enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'liquid flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%liquid_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%oil_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy flux'
  call FlowSubConditionVerify(option,condition,word,toil_ims%energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'oil water contact'
  call FlowSubConditionVerify(option,condition,word,toil_ims%owc, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'water pressure gradient'
  call FlowSubConditionVerify(option,condition,word,toil_ims%liq_press_grad, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,toil_ims%rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'well pressure'
  call FlowSubConditionVerify(option,condition,word,flow_well%pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'well rate'
  call FlowSubConditionVerify(option,condition,word,flow_well%rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'well temperature'
  call FlowSubConditionVerify(option,condition,word,flow_well%temperature, &
                              default_time_storage, &
                              PETSC_TRUE)

  condition%num_sub_conditions = 0
  i = 0
  if (associated(toil_ims%pressure)) &
    i = i + 1
  if (associated(toil_ims%saturation)) &
    i = i + 1
  if (associated(toil_ims%temperature)) &
    i = i + 1
  if (associated(toil_ims%enthalpy)) &
    i = i + 1
  if (associated(toil_ims%liquid_flux)) &
    i = i + 1
  if (associated(toil_ims%oil_flux)) &
    i = i + 1
  if (associated(toil_ims%energy_flux)) &
    i = i + 1
  if (associated(toil_ims%owc)) &
    i = i + 1
  if (associated(toil_ims%liq_press_grad)) &
    i = i + 1
  if (associated(toil_ims%rate)) &
    i = i + 1
  if (associated(flow_well%pressure)) &
    i = i + 1
  if (associated(flow_well%rate)) &
    i = i + 1
  if (associated(flow_well%temperature)) &
    i = i + 1
  ! assing number of sub_condition
  condition%num_sub_conditions = i
  allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    nullify(condition%sub_condition_ptr(idof)%ptr)
  enddo
  i = 0
  if (associated(toil_ims%pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%pressure
  endif
  if (associated(toil_ims%saturation)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%saturation
  endif  
  if (associated(toil_ims%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%temperature
  endif 
  if (associated(toil_ims%enthalpy)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%enthalpy
  endif  
  if (associated(toil_ims%liquid_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%liquid_flux
  endif  
  if (associated(toil_ims%oil_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%oil_flux
  endif  
  if (associated(toil_ims%energy_flux)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%energy_flux
  endif
  if (associated(toil_ims%owc)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%owc
  endif
  if (associated(toil_ims%liq_press_grad)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%liq_press_grad
  endif  
  if (associated(toil_ims%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => toil_ims%rate
  endif
  if (associated(flow_well%pressure)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => flow_well%pressure
  endif
  if (associated(flow_well%rate)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => flow_well%rate
  endif
  if (associated(flow_well%temperature)) then
    i = i + 1
    condition%sub_condition_ptr(i)%ptr => flow_well%temperature
  endif

  ! set condition types
  allocate(condition%itype(condition%num_sub_conditions))
  do idof = 1, condition%num_sub_conditions
    condition%itype(idof) = condition%sub_condition_ptr(idof)%ptr%itype
  enddo
  
  condition%default_time_storage => default_time_storage
    
  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionTOilImsRead

! ************************************************************************** !

subroutine TranConditionRead(condition,constraint_list,reaction,input,option)
  ! 
  ! Reads a transport condition from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module  
  use Units_module
  use Reaction_Aux_module
  
  implicit none
  
  type(tran_condition_type) :: condition
  type(tran_constraint_list_type) :: constraint_list
  type(reaction_type) :: reaction
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(tran_constraint_type), pointer :: constraint
  type(tran_constraint_coupler_type), pointer :: constraint_coupler, cur_coupler
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, internal_units
  PetscReal :: default_time
  character(len=MAXWORDLENGTH) :: default_time_units
  PetscInt :: default_iphase
  PetscInt :: default_itype
  PetscBool :: found
  PetscInt :: icomp
  PetscBool :: minerals_exist
  PetscErrorCode :: ierr
  PetscReal :: conversion

  call PetscLogEventBegin(logging%event_tran_condition_read, &
                          ierr);CHKERRQ(ierr)

  default_time = 0.d0
  default_iphase = 0
  default_time_units = ''

  ! read the condition
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')
          
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','CONDITION')   
      
    select case(trim(word))
    
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'TYPE','CONDITION')   
        call StringToLower(word)
        select case(word)
            case('dirichlet')
              condition%itype = DIRICHLET_BC
            case('dirichlet_zero_gradient')
              condition%itype = DIRICHLET_ZERO_GRADIENT_BC
            case('equilibrium')
              condition%itype = EQUILIBRIUM_SS
            case('neumann')
              condition%itype = NEUMANN_BC
            case('mole','mole_rate')
              condition%itype = MASS_RATE_SS
            case('zero_gradient')
              condition%itype = ZERO_GRADIENT_BC
            case default
              call InputKeywordUnrecognized(word,'transport condition type', &
                                            option)
        end select
      case('TIME')
        call InputReadDouble(input,option,default_time)
        call InputErrorMsg(input,option,'TIME','CONDITION') 
      case('TIME_UNITS') 
        call InputReadWord(input,option,word,PETSC_TRUE) 
        call InputErrorMsg(input,option,'UNITS','CONDITION')   
        select case(trim(word))     
          case('s','sec','min','m','hr','h','d','day','y','yr')
            default_time_units = trim(word)         
          case default
            option%io_buffer = 'Units "' // trim(word) // '" not recognized.'
            call printErrMsg(option)
        end select          
      case('CONSTRAINT_LIST')
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONSTRAINT')
              
          if (InputCheckExit(input,option)) exit  

          constraint_coupler => TranConstraintCouplerCreate(option)
          call InputReadDouble(input,option,constraint_coupler%time)
          call InputErrorMsg(input,option,'time','CONSTRAINT_LIST') 
          ! time units are optional  
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'constraint name','CONSTRAINT_LIST') 
          ! read constraint name
          call InputReadWord(input,option,constraint_coupler%constraint_name, &
                             PETSC_TRUE)
          if (InputError(input)) then
            constraint_coupler%time_units = default_time_units
            constraint_coupler%constraint_name = trim(word)
          else
            constraint_coupler%time_units = word
          endif
          ! convert time units
          if (len_trim(constraint_coupler%time_units) > 0) then
            internal_units = 'sec'
            constraint_coupler%time = constraint_coupler%time* &
              UnitsConvertToInternal(constraint_coupler%time_units, &
                                     internal_units,option)
          endif
          ! add to end of list
          if (.not.associated(condition%constraint_coupler_list)) then
            condition%constraint_coupler_list => constraint_coupler
          else
            cur_coupler => condition%constraint_coupler_list
            do
              if (.not.associated(cur_coupler%next)) exit
              cur_coupler => cur_coupler%next
            enddo
            cur_coupler%next => constraint_coupler
          endif
        enddo
      case('CONSTRAINT')
        constraint => TranConstraintCreate(option)
        constraint_coupler => TranConstraintCouplerCreate(option)
        call InputReadWord(input,option,constraint%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'constraint','name') 
        option%io_buffer = 'Constraint: ' // trim(constraint%name)
        call printMsg(option)
        call TranConstraintRead(constraint,reaction,input,option)
        call TranConstraintAddToList(constraint,constraint_list)
        constraint_coupler%aqueous_species => constraint%aqueous_species
        constraint_coupler%minerals => constraint%minerals
        constraint_coupler%surface_complexes => constraint%surface_complexes
        constraint_coupler%colloids => constraint%colloids
        constraint_coupler%immobile_species => constraint%immobile_species
        constraint_coupler%time = default_time
        ! add to end of coupler list
        if (.not.associated(condition%constraint_coupler_list)) then
          condition%constraint_coupler_list => constraint_coupler
        else
          cur_coupler => condition%constraint_coupler_list
          do
            if (.not.associated(cur_coupler%next)) exit
            cur_coupler => cur_coupler%next
          enddo
          cur_coupler%next => constraint_coupler
        endif        
      case default
        call InputKeywordUnrecognized(word,'transport condition',option)
    end select 
  
  enddo  

  if (.not.associated(condition%constraint_coupler_list)) then
    option%io_buffer = 'No CONSTRAINT or CONSTRAINT_LIST defined in &
                       &Transport Condition "' // trim(condition%name) // '".'
    call printErrMsg(option)
  endif

  if (len_trim(default_time_units) > 0) then
    internal_units = 'sec'
    conversion = UnitsConvertToInternal(default_time_units,internal_units, &
                                        option)
    cur_coupler => condition%constraint_coupler_list
    do
      if (.not.associated(cur_coupler)) exit
      if (len_trim(cur_coupler%time_units) == 0) then
        cur_coupler%time = cur_coupler%time*conversion
      endif
      cur_coupler => cur_coupler%next
    enddo
  endif

  call PetscLogEventEnd(logging%event_tran_condition_read,ierr);CHKERRQ(ierr)

end subroutine TranConditionRead

! ************************************************************************** !

subroutine ConditionReadValues(input,option,keyword,dataset_base,units, &
                               data_internal_units)
  ! 
  ! Read the value(s) of a condition variable
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  use Logging_module
  use HDF5_Aux_module
  use Units_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  class(dataset_base_type), pointer :: dataset_base
  character(len=MAXWORDLENGTH) :: units
  character(len=MAXSTRINGLENGTH), pointer :: internal_unit_strings(:)
  character(len=MAXWORDLENGTH) :: data_internal_units
  
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i, icount
  PetscInt :: icol
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_buffer(:)
  type(input_type), pointer :: input2
  PetscErrorCode :: ierr

#if defined(PETSC_HAVE_HDF5)  
  integer(HID_T) :: file_id
  integer(HID_T) :: prop_id
  PetscMPIInt :: hdf5_err
#endif

  call PetscLogEventBegin(logging%event_flow_condition_read_values, &
                          ierr);CHKERRQ(ierr)
  
  ! dataset_base, though of type dataset_base_type, should always be created
  ! as dataset_ascii_type.
  dataset_ascii => DatasetAsciiCast(dataset_base)
  if (.not.associated(dataset_ascii)) then
    ! The dataset was not of type dataset_asci and was likely set to a different
    ! type.  There is a bug in the input file.
    option%io_buffer = 'Dataset associated with ' // trim(keyword) // &
      ' in the input file is already associated with a different dataset ' // &
      'type.  Check for duplicate definitions of ' // trim(keyword) // '.'
    call printErrMsg(option)
  endif

  nullify(input2)
  filename = ''
  realization_word = ''
  hdf5_path = ''

  internal_unit_strings => StringSplit(data_internal_units,',')
  
  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToLower(word)
  length = len_trim(word)
  if (StringStartsWithAlpha(word)) then
    if (length == FOUR_INTEGER .and. StringCompare(word,'file',FOUR_INTEGER)) then 
      input%err_buf2 = trim(keyword) // ', FILE'
      input%err_buf = 'keyword'
      call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
      if (input%ierr == 0) then
        filename = string2
      else
        option%io_buffer = 'The ability to read realization dependent ' // &
          'datasets outside the DATASET block is no longer supported'
        call printErrMsg(option)
      endif
    
      if (len_trim(filename) < 2) then
        option%io_buffer = 'No filename listed under Flow_Condition: ' // &
                           trim(keyword)
        call printErrMsg(option)
      endif

      if (index(filename,'.h5') > 0) then
        write(option%io_buffer,'("Reading of HDF5 datasets for flow ", &
                                 &"conditions not currently supported.")')
        call printErrMsg(option)
#if 0      
#if !defined(PETSC_HAVE_HDF5)
        write(option%io_buffer,'("PFLOTRAN must be compiled with HDF5 to ", &
                                 &"read HDF5 formatted flow conditions.")')
        call printErrMsg(option)
#else   
        if (len_trim(hdf5_path) < 1) then
          option%io_buffer = 'No hdf5 path listed under Flow_Condition: ' // &
                             trim(keyword)
          call printErrMsg(option)
        endif

        call h5open_f(hdf5_err)
        option%io_buffer = 'Opening hdf5 file: ' // trim(filename)
        call printMsg(option)
        call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
        call h5pset_fapl_mpio_f(prop_id,option%mycomm,MPI_INFO_NULL,hdf5_err)
#endif
        call HDF5OpenFileReadOnly(filename,file_id,prop_id,option)
        call h5pclose_f(prop_id,hdf5_err)

        hdf5_path = trim(hdf5_path) // trim(realization_word)
        call HDF5ReadNDimRealArray(option,file_id,hdf5_path,ndims,dims, &
                                   real_buffer)
        option%io_buffer = 'Closing hdf5 file: ' // trim(filename)
        call printMsg(option)  
        call h5fclose_f(file_id,hdf5_err)
        call h5close_f(hdf5_err)
      
        ! dims(1) = size of array
        ! dims(2) = number of data point in time
        if (dims(1)-1 == flow_dataset%time_series%rank) then
          ! alright, the 2d data is layed out in C-style.  now place it in
          ! the appropriate arrays
          allocate(flow_dataset%time_series%times(dims(2)))
          flow_dataset%time_series%times = UNINITIALIZED_DOUBLE
          allocate(flow_dataset%time_series%values(flow_dataset%time_series%rank,dims(2))) 
          flow_dataset%time_series%values = UNINITIALIZED_DOUBLE
          icount = 1
          do i = 1, dims(2)
            flow_dataset%time_series%times(i) = real_buffer(icount)
            icount = icount + 1
            do icol = 1, flow_dataset%time_series%rank
              flow_dataset%time_series%values(icol,i) = real_buffer(icount)
              icount = icount + 1
            enddo
          enddo  
        else
          option%io_buffer = 'HDF condition data set rank does not match' // &
            'rank of internal data set.  Email Glenn for additions'
          call printErrMsg(option)
        endif
        if (associated(dims)) deallocate(dims)
        nullify(dims)
        if (associated(real_buffer)) deallocate(real_buffer)
        nullify(real_buffer)
#endif    
#endif    
! if 0
      else
        i = index(filename,'.',PETSC_TRUE)
        if (i > 2) then
          filename = filename(1:i-1) // trim(realization_word) // filename(i:)
        else
          filename = trim(filename) // trim(realization_word)
        endif
        input2 => InputCreate(IUNIT_TEMP,filename,option)
        input2%force_units = input%force_units
        call DatasetAsciiRead(dataset_ascii,input2,data_internal_units,option)
        dataset_ascii%filename = filename
        call InputDestroy(input2)
      endif
    else if (StringCompare(word,'dataset')) then
      call InputReadWord(input,option,word,PETSC_TRUE)    
      input%err_buf2 = trim(keyword) // ', DATASET'
      input%err_buf = 'dataset name'
      call InputErrorMsg(input,option)
      call DatasetDestroy(dataset_base)
      dataset_base => DatasetBaseCreate()
      dataset_base%name = word
    else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then 
      call DatasetAsciiRead(dataset_ascii,input,data_internal_units,option)
    else
      option%io_buffer = 'Keyword "' // trim(word) // &
        '" not recognized in when reading condition values for "' // &
        trim(keyword) // '".'
      call printErrMsg(option)
    endif
  else
    input%buf = trim(string2)
    allocate(dataset_ascii%rarray(dataset_ascii%array_width))
    do icol=1,dataset_ascii%array_width
      call InputReadDouble(input,option,dataset_ascii%rarray(icol))
      write(input%err_buf,'(a,i2)') trim(keyword) // &
                                    ' dataset_values, icol = ', icol
      input%err_buf2 = 'CONDITION'
      call InputErrorMsg(input,option) 
    enddo
    string2 = input%buf
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (InputError(input)) then
      call InputCheckMandatoryUnits(input,option)
      word = trim(keyword) // ' UNITS'
      call InputDefaultMsg(input,option,word)
    else
      input%buf = string2
      units = ''
      do icol=1,dataset_ascii%array_width
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,'CONDITION')   
        dataset_ascii%rarray(icol) = UnitsConvertToInternal(word, &
                                     internal_unit_strings(icol),option) * &
                                     dataset_ascii%rarray(icol)
        units = trim(units) // ' ' // trim(word)
      enddo
    endif
  endif
  
  deallocate(internal_unit_strings)
  nullify(internal_unit_strings)  

  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
                        ierr);CHKERRQ(ierr)

end subroutine ConditionReadValues

! ************************************************************************** !

subroutine FlowConditionPrint(condition,option)
  ! 
  ! Prints flow condition info
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  ! 

  use Option_module
  use Dataset_module

  implicit none
  
  type(flow_condition_type) :: condition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Flow Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') trim(condition%length_units)
  
100 format(6x,a)  
  write(option%fid_out,100) 'Datum:'
  if (associated(condition%datum)) then
    call DatasetPrint(condition%datum,option)
  endif
  
  do i=1, condition%num_sub_conditions
    call FlowConditionPrintSubCondition(condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)
  
end subroutine FlowConditionPrint

! ************************************************************************** !

subroutine FlowConditionPrintSubCondition(subcondition,option)
  ! 
  ! Prints flow subcondition info
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  ! 

  use Option_module
  use Dataset_module

  implicit none
  
  type(flow_sub_condition_type) :: subcondition
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  string = GetSubConditionName(subcondition%itype)

  105 format(6x,'Type: ',a)  
  write(option%fid_out,105) trim(string)
  
  110 format(6x,a)  
  
  write(option%fid_out,110) 'Gradient:'
  if (associated(subcondition%gradient)) then
    call DatasetPrint(subcondition%gradient,option)
  endif

  write(option%fid_out,110) 'Data:'
  if (associated(subcondition%dataset)) then
    call DatasetPrint(subcondition%dataset,option)
  endif
            
end subroutine FlowConditionPrintSubCondition

! ************************************************************************** !

function GetSubConditionName(subcon_itype)
  ! 
  ! SubConditionName: Return name of subcondition
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/16/13
  ! 

  implicit none

  PetscInt :: subcon_itype

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: GetSubConditionName

  select case(subcon_itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(DIRICHLET_ZERO_GRADIENT_BC)
      string = 'dirichlet-zero gradient'
    case(MASS_RATE_SS)
      string = 'mass_rate'
    case(WELL_SS)
      string = 'well'
    case(HYDROSTATIC_BC)
      string = 'hydrostatic'
    case(CONDUCTANCE_BC)
      string = 'conductance'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
    case(SEEPAGE_BC)
      string = 'seepage'
    case(VOLUMETRIC_RATE_SS)
      string = 'volumetric rate'
    case(EQUILIBRIUM_SS)
      string = 'equilibrium'
    case(UNIT_GRADIENT_BC)
      string = 'unit gradient'
    case(SCALED_MASS_RATE_SS)
      string = 'scaled mass rate'
    case(SCALED_VOLUMETRIC_RATE_SS)
      string = 'scaled volumetric rate'
    case(HET_VOL_RATE_SS)
      string = 'heterogeneous volumetric rate'
    case(HET_MASS_RATE_SS)
      string = 'heterogeneous mass rate'
    case(HET_DIRICHLET)
      string = 'heterogeneous dirichlet'
    case(ENERGY_RATE_SS)
      string = 'energy rate'
    case(SCALED_ENERGY_RATE_SS)
      string = 'scaled energy rate'
    case(HET_ENERGY_RATE_SS)
      string = 'heterogeneous energy rate'
    case(HET_SURF_SEEPAGE_BC)
      string = 'heterogeneous surface seepage'
    case(SPILLOVER_BC)
      string = 'spillover'
    case(WELL_MASS_RATE_TARGET)
      string = 'well target mass rate'
    case(WELL_MASS_RATE_MAX)
      string = 'well maximum mass rate'
    case(WELL_MASS_RATE_MIN)
      string = 'well minimum mass rate'
    case(WELL_VOL_RATE_TARGET)
      string = 'well target volumetric rate'
    case(WELL_VOL_RATE_MAX)
      string = 'well maximum volumetric rate'
    case(WELL_VOL_RATE_MIN)
      string = 'well minimum volumetric rate'
    case(WELL_BHP)
      string = 'well bottom hole pressure'
    case(WELL_BHP_MIN)
      string = 'well minimum bottom hole pressure'
    case(WELL_BHP_MAX)
      string = 'well maximum bottom hole pressure'
    case(SURFACE_DIRICHLET)
      string = 'surface_dirichlet'
    case(SURFACE_ZERO_GRADHEIGHT)
      string = 'surface_zero_gradheight'
    case(SURFACE_SPILLOVER)
      string = 'surface_spillover'
    end select

  GetSubConditionName = trim(string)

end function GetSubConditionName

! ************************************************************************** !

subroutine FlowConditionUpdate(condition_list,option,time)
  ! 
  ! Updates a transient condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  use Option_module
  use Dataset_module
  
  implicit none
  
  type(condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  
  type(flow_condition_type), pointer :: condition
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition  
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    call DatasetUpdate(condition%datum,time,option)
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr
      
      if (associated(sub_condition)) then
        call DatasetUpdate(sub_condition%dataset,time,option)
        call DatasetUpdate(sub_condition%gradient,time,option)
      endif
      
    enddo
      
    condition => condition%next
    
  enddo
  
end subroutine FlowConditionUpdate

! ************************************************************************** !

subroutine TranConditionUpdate(condition_list,option,time)
  ! 
  ! Updates a transient transport condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  use Option_module
  
  implicit none
  
  type(tran_condition_list_type) :: condition_list
  type(option_type) :: option
  PetscReal :: time
  
  type(tran_condition_type), pointer :: condition
  
  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    
    do
      if (associated(condition%cur_constraint_coupler%next)) then
        if (time >= condition%cur_constraint_coupler%next%time) then
          condition%cur_constraint_coupler => &
            condition%cur_constraint_coupler%next
        else
          exit
        endif
      else 
        exit
      endif
    enddo
    condition => condition%next
    
  enddo
  
end subroutine TranConditionUpdate

! ************************************************************************** !

subroutine FlowConditionInitList(list)
  ! 
  ! Initializes a condition list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none

  type(condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine FlowConditionInitList

! ************************************************************************** !

subroutine FlowConditionAddToList(new_condition,list)
  ! 
  ! Adds a new condition to a condition list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(flow_condition_type), pointer :: new_condition
  type(condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine FlowConditionAddToList

! ************************************************************************** !

function FlowConditionGetPtrFromList(condition_name,condition_list)
  ! 
  ! Returns a pointer to the condition matching &
  ! condition_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use String_module
  
  implicit none
  
  type(flow_condition_type), pointer :: FlowConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list
 
  PetscInt :: length
  type(flow_condition_type), pointer :: condition
    
  nullify(FlowConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      FlowConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function FlowConditionGetPtrFromList

! ************************************************************************** !

subroutine TranConditionInitList(list)
  ! 
  ! Initializes a transport condition list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  ! 

  implicit none

  type(tran_condition_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine TranConditionInitList

! ************************************************************************** !

subroutine TranConditionAddToList(new_condition,list)
  ! 
  ! Adds a new condition to a transport condition list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  ! 

  implicit none
  
  type(tran_condition_type), pointer :: new_condition
  type(tran_condition_list_type) :: list
  
  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition
  
end subroutine TranConditionAddToList

! ************************************************************************** !

function TranConditionGetPtrFromList(condition_name,condition_list)
  ! 
  ! Returns a pointer to the condition matching
  ! condition_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/08
  ! 

  use String_module

  implicit none
  
  type(tran_condition_type), pointer :: TranConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(tran_condition_list_type) :: condition_list
 
  PetscInt :: length
  type(tran_condition_type), pointer :: condition
    
  nullify(TranConditionGetPtrFromList)
  condition => condition_list%first
  
  do 
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                        length)) then
      TranConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo
  
end function TranConditionGetPtrFromList

! ************************************************************************** !

function FlowConditionIsTransient(condition)
  ! 
  ! Returns PETSC_TRUE
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Dataset_module

  implicit none
  
  type(flow_condition_type) :: condition
  
  PetscBool :: FlowConditionIsTransient
  
  FlowConditionIsTransient = PETSC_FALSE

  if (DatasetIsTransient(condition%datum) .or. &
      FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%concentration) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%well) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%energy_rate) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. &
      FlowConditionTOilImsIsTransient(condition%toil_ims) .or. &
      FlowWellConditionIsTransient(condition%flow_well) .or. &
      FlowConditionGeneralIsTransient(condition%general)) then
    FlowConditionIsTransient = PETSC_TRUE
  endif

end function FlowConditionIsTransient

! ************************************************************************** !

function FlowConditionGeneralIsTransient(condition)
  ! 
  ! Returns PETSC_TRUE
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Dataset_module

  implicit none
  
  type(flow_general_condition_type), pointer :: condition
  
  PetscBool :: FlowConditionGeneralIsTransient
  
  FlowConditionGeneralIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return
  
  if (FlowSubConditionIsTransient(condition%liquid_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_pressure) .or. &
      FlowSubConditionIsTransient(condition%gas_saturation) .or. &
      FlowSubConditionIsTransient(condition%relative_humidity) .or. &
      FlowSubConditionIsTransient(condition%mole_fraction) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%gas_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux)) then
    FlowConditionGeneralIsTransient = PETSC_TRUE
  endif
  
end function FlowConditionGeneralIsTransient

! ************************************************************************** !

function FlowConditionTOilImsIsTransient(condition)
  ! 
  ! Returns PETSC_TRUE
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/20/15
  ! 

  use Dataset_module

  implicit none
  
  type(flow_toil_ims_condition_type), pointer :: condition
  
  PetscBool :: FlowConditionTOilImsIsTransient
  
  FlowConditionTOilImsIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return
  
  if (FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%liquid_flux) .or. &
      FlowSubConditionIsTransient(condition%oil_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. & 
      FlowSubConditionIsTransient(condition%owc) .or. &
      FlowSubConditionIsTransient(condition%liq_press_grad)) then
    FlowConditionTOilImsIsTransient = PETSC_TRUE
  endif
  
end function FlowConditionTOilImsIsTransient

! ************************************************************************** !

function FlowWellConditionIsTransient(condition)
  ! 
  ! Returns PETSC_TRUE
  ! 
  ! Author: Paolo Orsini
  ! Date: 8/05/16
  ! 

  use Dataset_module

  implicit none
  
  type(flow_well_condition_type), pointer :: condition
  
  PetscBool :: FlowWellConditionIsTransient
  
  FlowWellConditionIsTransient = PETSC_FALSE

  if (.not.associated(condition)) return
  
  if (FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%temperature)) then
    FlowWellConditionIsTransient = PETSC_TRUE
  endif

end function FlowWellConditionIsTransient

! ************************************************************************** !

function FlowSubConditionIsTransient(sub_condition)
  ! 
  ! Returns PETSC_TRUE
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Dataset_module

  implicit none
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  PetscBool :: FlowSubConditionIsTransient
  
  FlowSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (DatasetIsTransient(sub_condition%dataset) .or. &
        DatasetIsTransient(sub_condition%gradient)) then
      FlowSubConditionIsTransient = PETSC_TRUE
    endif
  endif  
  
end function FlowSubConditionIsTransient

! ************************************************************************** !

function FlowConditionUnknownItype(condition,message,type_name)
  ! 
  ! Returns a string indicating which flow condition has a wrong type.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/21/16
  ! 
  implicit none
  
  type(flow_condition_type) :: condition
  character(len=*) :: message
  character(len=*) :: type_name
  
  character(len=MAXSTRINGLENGTH) :: FlowConditionUnknownItype

  FlowConditionUnknownItype = 'Unknown TYPE (' // trim(type_name) // &
    ') for ' // trim(message) // ' within FLOW_CONDITION "' // &
    trim(condition%name) // '".'
  
end function FlowConditionUnknownItype

! **************************************************************************** !

subroutine FlowCondInputRecord(flow_condition_list,option)
  ! 
  ! Prints ingested flow condition information to 
  ! the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/19/2016
  ! 
  use Option_module
  use Dataset_Base_class

  implicit none
  
  type(condition_list_type), pointer :: flow_condition_list
  type(option_type), pointer :: option

  type(flow_condition_type), pointer :: cur_fc
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT
  
  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'FLOW CONDITIONS'
  
  cur_fc => flow_condition_list%first
  do
    if (.not.associated(cur_fc)) exit
    write(id,'(a29)',advance='no') 'flow condition name: '
    write(id,'(a)') adjustl(trim(cur_fc%name))
    if (cur_fc%num_sub_conditions > 0) then
      do k = 1,cur_fc%num_sub_conditions
        write(id,'(a29)',advance='no') 'sub condition name: '
        write(id,'(a)') adjustl(trim(cur_fc%sub_condition_ptr(k)%ptr%name))
        write(id,'(a29)',advance='no') 'sub condition type: '
        write(id,'(a)') adjustl(trim(cur_fc%sub_condition_ptr(k)%ptr%ctype))
        if (associated(cur_fc%sub_condition_ptr(k)%ptr%dataset)) then
          call DatasetBasePrint(cur_fc%sub_condition_ptr(k)%ptr%dataset,option)
          ! DatasetBasePrint doesn't seem to do anything?
        endif
        if (associated(cur_fc%sub_condition_ptr(k)%ptr%gradient)) then
          call DatasetBasePrint(cur_fc%sub_condition_ptr(k)%ptr%gradient,option)
        endif
      enddo
    endif
    
    write(id,'(a29)') '---------------------------: '
    cur_fc => cur_fc%next
  enddo
  
end subroutine FlowCondInputRecord

! **************************************************************************** !

subroutine TranCondInputRecord(tran_condition_list,option)
  ! 
  ! Prints ingested transport condition information to 
  ! the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/19/2016
  ! 
  use Option_module
  use Dataset_Base_class
  Use Transport_Constraint_module

  implicit none
  
  type(tran_condition_list_type), pointer :: tran_condition_list
  type(option_type), pointer :: option

  type(tran_constraint_coupler_type), pointer :: cur_tcon_coupler
  type(tran_condition_type), pointer :: cur_tc
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT
  
  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'TRANSPORT CONDITIONS'
  
  cur_tc => tran_condition_list%first
  do
    if (.not.associated(cur_tc)) exit
    write(id,'(a29)',advance='no') 'transport condition name: '
    write(id,'(a)') adjustl(trim(cur_tc%name))
    write(id,'(a29)',advance='no') 'transport condition type: '
    select case (cur_tc%itype)
      case(DIRICHLET_BC)
        write(id,'(a)') 'dirichlet'
      case(DIRICHLET_ZERO_GRADIENT_BC)
        write(id,'(a)') 'dirichlet_zero_gradient'
      case(EQUILIBRIUM_SS)
        write(id,'(a)') 'equilibrium'
      case(NEUMANN_BC)
        write(id,'(a)') 'neumann'
      case(MASS_RATE_SS)
        write(id,'(a)') 'mole_rate'
      case(ZERO_GRADIENT_BC)
        write(id,'(a)') 'zero_gradient'
    end select
    write(id,'(a29)',advance='no') 'is transient?: '
    if (cur_tc%is_transient) then
      write(id,'(a)') 'YES'
    else
      write(id,'(a)') 'NO'
    endif
    cur_tcon_coupler => cur_tc%constraint_coupler_list
    do
      if (.not.associated(cur_tcon_coupler)) exit
      write(id,'(a29)',advance='no') 'transport constraint name: '
      write(id,'(a)') adjustl(trim(cur_tcon_coupler%constraint_name))
      
      ! aqueous species concentraion constraint
      if (associated(cur_tcon_coupler%aqueous_species)) then
        do k = 1,size(cur_tcon_coupler%aqueous_species%names)
          write(id,'(a29)',advance='no') 'aqueous species constraint: '
          write(string,*) trim(cur_tcon_coupler%aqueous_species%names(k))
          select case (cur_tcon_coupler%aqueous_species%constraint_type(k))
            case(CONSTRAINT_FREE)
              string = trim(string) // ', free'
            case(CONSTRAINT_TOTAL)
              string = trim(string) // ', total'
            case(CONSTRAINT_TOTAL_SORB)
              string = trim(string) // ', total_sorb'
            case(CONSTRAINT_PH)
              string = trim(string) // ', ph'
            case(CONSTRAINT_LOG)
              string = trim(string) // ', log'
            case(CONSTRAINT_MINERAL)
              string = trim(string) // ', mineral'
            case(CONSTRAINT_GAS)
              string = trim(string) // ', gas'
            case(CONSTRAINT_SUPERCRIT_CO2)
              string = trim(string) // ', super critical CO2'
            case(CONSTRAINT_CHARGE_BAL)
              string = trim(string) // ', charge balance'
          end select
          write(word1,*) cur_tcon_coupler%aqueous_species%constraint_conc(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' mol'
        enddo
      endif
      
      ! free-ion guess constraint
      if (associated(cur_tcon_coupler%free_ion_guess)) then
        do k = 1,size(cur_tcon_coupler%free_ion_guess%names)
          write(id,'(a29)',advance='no') 'free ion guess constraint: '
          write(string,*) trim(cur_tcon_coupler%free_ion_guess%names(k))
          write(word1,*) cur_tcon_coupler%free_ion_guess%conc(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' mol'
        enddo
      endif
      
      ! mineral constraint
      if (associated(cur_tcon_coupler%minerals)) then
        do k = 1,size(cur_tcon_coupler%minerals%names)
          write(id,'(a29)',advance='no') 'mineral vol. frac. constraint: '
          write(string,*) trim(cur_tcon_coupler%minerals%names(k))
          write(word1,*) cur_tcon_coupler%minerals%constraint_vol_frac(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' m^3/m^3'
          write(id,'(a29)',advance='no') 'mineral area constraint: '
          write(string,*) trim(cur_tcon_coupler%minerals%names(k))
          write(word1,*) cur_tcon_coupler%minerals%constraint_area(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' m^2/m^3'
        enddo
      endif
      
      ! surface complexes constraint
      if (associated(cur_tcon_coupler%surface_complexes)) then
        do k = 1,size(cur_tcon_coupler%surface_complexes%names)
          write(id,'(a29)',advance='no') 'surface complex constraint: '
          write(string,*) trim(cur_tcon_coupler%surface_complexes%names(k))
          write(word1,*) cur_tcon_coupler%surface_complexes%constraint_conc(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' mol/m^3'
        enddo
      endif
      
      ! colloids constraint
      if (associated(cur_tcon_coupler%colloids)) then
        do k = 1,size(cur_tcon_coupler%colloids%names)
          write(id,'(a29)',advance='no') 'colloid mobile constraint: '
          write(string,*) trim(cur_tcon_coupler%colloids%names(k))
          write(word1,*) cur_tcon_coupler%colloids%constraint_conc_mob(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' mol'
          write(id,'(a29)',advance='no') 'colloid immobile constraint: '
          write(string,*) trim(cur_tcon_coupler%colloids%names(k))
          write(word1,*) cur_tcon_coupler%colloids%constraint_conc_imb(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) & 
                          // ' mol'
        enddo
      endif
      
      ! immobile species constraint
      if (associated(cur_tcon_coupler%immobile_species)) then
        do k = 1,size(cur_tcon_coupler%immobile_species%names)
          write(id,'(a29)',advance='no') 'immobile species constraint: '
          write(string,*) trim(cur_tcon_coupler%immobile_species%names(k))
          write(word1,*) cur_tcon_coupler%immobile_species%constraint_conc(k)
          write(id,'(a)') trim(string) // ', ' // adjustl(trim(word1)) &
                          // ' mol/m^3'
        enddo
      endif
      
      cur_tcon_coupler => cur_tcon_coupler%next
    enddo
    
    write(id,'(a29)') '---------------------------: '
    cur_tc => cur_tc%next
  enddo
  
end subroutine TranCondInputRecord
  
! ************************************************************************** !

subroutine FlowConditionDestroyList(condition_list)
  ! 
  ! Deallocates a list of conditions
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(condition_list_type), pointer :: condition_list
  
  type(flow_condition_type), pointer :: condition, prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call FlowConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine FlowConditionDestroyList

! ************************************************************************** !

subroutine FlowConditionDestroy(condition)
  ! 
  ! Deallocates a condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Dataset_module
  use Dataset_Ascii_class
  use Utility_module
  
  implicit none
  
  type(flow_condition_type), pointer :: condition
  
  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscInt :: i
  
  if (.not.associated(condition)) return
  
  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(condition%datum)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)

  if (associated(condition%sub_condition_ptr)) then
    ! only nullify the pointers; don't destroy as they are destroyed below
    do i=1,condition%num_sub_conditions
      nullify(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  call DeallocateArray(condition%itype)
  
  call FlowSubConditionDestroy(condition%pressure)
  call FlowSubConditionDestroy(condition%saturation)
  call FlowSubConditionDestroy(condition%rate)
  call FlowSubConditionDestroy(condition%well)
  call FlowSubConditionDestroy(condition%temperature)
  call FlowSubConditionDestroy(condition%concentration)
  call FlowSubConditionDestroy(condition%enthalpy)
  call FlowSubConditionDestroy(condition%energy_rate)

  call TimeStorageDestroy(condition%default_time_storage)
  call FlowGeneralConditionDestroy(condition%general)
  call FlowToilConditionDestroy(condition%toil_ims)
  call FlowWellConditionDestroy(condition%flow_well)  

  nullify(condition%next)  
  
  deallocate(condition)
  nullify(condition)

end subroutine FlowConditionDestroy

! ************************************************************************** !

subroutine FlowGeneralConditionDestroy(general_condition)
  ! 
  ! Destroys a general mode condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/26/11
  ! 

  use Option_module
  
  implicit none
  
  type(flow_general_condition_type), pointer :: general_condition
  
  if (.not.associated(general_condition)) return

  call FlowSubConditionDestroy(general_condition%liquid_pressure)
  call FlowSubConditionDestroy(general_condition%gas_pressure)
  call FlowSubConditionDestroy(general_condition%gas_saturation)
  call FlowSubConditionDestroy(general_condition%relative_humidity)
  call FlowSubConditionDestroy(general_condition%mole_fraction)
  call FlowSubConditionDestroy(general_condition%temperature)
  call FlowSubConditionDestroy(general_condition%liquid_flux)
  call FlowSubConditionDestroy(general_condition%gas_flux)
  call FlowSubConditionDestroy(general_condition%energy_flux)
  call FlowSubConditionDestroy(general_condition%rate)

  deallocate(general_condition)
  nullify(general_condition)

end subroutine FlowGeneralConditionDestroy

! ************************************************************************** !

subroutine FlowToilConditionDestroy(toil_ims_condition)
  ! 
  ! Destroys a toil_ims mode condition
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/06/15
  ! 

  use Option_module
  
  implicit none
  
  type(flow_toil_ims_condition_type), pointer :: toil_ims_condition

  if (.not.associated(toil_ims_condition)) return

  call FlowSubConditionDestroy(toil_ims_condition%pressure)
  call FlowSubConditionDestroy(toil_ims_condition%saturation)
  call FlowSubConditionDestroy(toil_ims_condition%temperature)
  call FlowSubConditionDestroy(toil_ims_condition%enthalpy)
  call FlowSubConditionDestroy(toil_ims_condition%liquid_flux)
  call FlowSubConditionDestroy(toil_ims_condition%oil_flux)
  call FlowSubConditionDestroy(toil_ims_condition%energy_flux)
  call FlowSubConditionDestroy(toil_ims_condition%rate)
  call FlowSubConditionDestroy(toil_ims_condition%owc)
  call FlowSubConditionDestroy(toil_ims_condition%liq_press_grad)

  deallocate(toil_ims_condition)
  nullify(toil_ims_condition)

end subroutine FlowToilConditionDestroy

! ************************************************************************** !

subroutine FlowWellConditionDestroy(flow_well_condition)
  ! 
  ! Destroys a toil_ims mode condition
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/06/15
  ! 

  use Option_module
  
  implicit none
  
  type(flow_well_condition_type), pointer :: flow_well_condition

  if (.not.associated(flow_well_condition)) return

  call FlowSubConditionDestroy(flow_well_condition%pressure)
  call FlowSubConditionDestroy(flow_well_condition%rate)
  call FlowSubConditionDestroy(flow_well_condition%temperature)

  deallocate(flow_well_condition)
  nullify(flow_well_condition)

end subroutine FlowWellConditionDestroy

! ************************************************************************** !

subroutine FlowSubConditionDestroy(sub_condition)
  ! 
  ! Destroys a sub_condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  ! 

  use Dataset_module
  use Dataset_Ascii_class

  implicit none
  
  type(flow_sub_condition_type), pointer :: sub_condition
  
  class(dataset_ascii_type), pointer :: dataset_ascii
  
  if (.not.associated(sub_condition)) return
  
  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(sub_condition%dataset)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)
  dataset_ascii => DatasetAsciiCast(sub_condition%gradient)
  call DatasetAsciiDestroy(dataset_ascii)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine FlowSubConditionDestroy

! ************************************************************************** !

subroutine TranConditionDestroyList(condition_list)
  ! 
  ! Deallocates a list of conditions
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(tran_condition_list_type), pointer :: condition_list
  
  type(tran_condition_type), pointer :: condition, prev_condition
  
  if (.not.associated(condition_list)) return
  
  condition => condition_list%first
  do 
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call TranConditionDestroy(prev_condition)
  enddo
  
  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)
  
  deallocate(condition_list)
  nullify(condition_list)

end subroutine TranConditionDestroyList

! ************************************************************************** !

subroutine TranConditionDestroy(condition)
  ! 
  ! Deallocates a condition
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(tran_condition_type), pointer :: condition
  
  if (.not.associated(condition)) return
  
  if (associated(condition%constraint_coupler_list)) &
    call TranConstraintCouplerDestroy(condition%constraint_coupler_list)

  deallocate(condition)
  nullify(condition)

end subroutine TranConditionDestroy

end module Condition_module
