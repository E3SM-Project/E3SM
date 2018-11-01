module Condition_module

#include "petsc/finclude/petscsys.h"
  use petscsys

!  use Reaction_Aux_module
!  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module

  use Transport_Constraint_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

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
    type(flow_sub_condition_type), pointer :: temperature
    type(flow_sub_condition_type), pointer :: concentration
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: energy_rate
    type(flow_sub_condition_type), pointer :: energy_flux
    ! any new sub conditions must be added to FlowConditionIsTransient
    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)
    type(flow_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type flow_condition_type

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
  nullify(condition%temperature)
  nullify(condition%concentration)
  nullify(condition%enthalpy)
  nullify(condition%sub_condition_ptr)
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

  print *, 'GetFlowSubCondFromArrayByName() needs to be updated to include &
           &the general_condition_type.'
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

  header = 'SUBSURFACE/FLOW_CONDITION/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Value(s)'
  call DatasetVerify(sub_condition%dataset,default_time_storage, &
                     header,option)
  header = 'SUBSURFACE/FLOW_CONDITION/' // &
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
  character(len=MAXWORDLENGTH) :: rate_unit_string
  character(len=MAXWORDLENGTH) :: energy_rate_unit_string
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

  rate_unit_string = 'not_assigned'
  energy_rate_unit_string = 'not_assigned'
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
  
    internal_units = 'not_assigned'

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
            case('RATE')
              sub_condition_ptr => rate
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
            case('WELL')
              sub_condition_ptr => well
            case('FLUX')
              sub_condition_ptr => flux
            case('ENERGY_FLUX')
              sub_condition_ptr => energy_flux
            case('SATURATION')
              sub_condition_ptr => saturation
            case('TEMPERATURE')
              sub_condition_ptr => temperature
            case('CONCENTRATION')
              sub_condition_ptr => concentration
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
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
              rate_unit_string = 'kg/sec'
            case('energy_rate')
              sub_condition_ptr%itype = ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('heterogeneous_energy_rate')
              sub_condition_ptr%itype = HET_ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('scaled_mass_rate','scaled_volumetric_rate', &
                 'scaled_energy_rate')
              select case(word)
                case('scaled_mass_rate')
                  sub_condition_ptr%itype = SCALED_MASS_RATE_SS
                  rate_unit_string = 'kg/sec'
                case('scaled_volumetric_rate')
                  sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
                  rate_unit_string = 'm^3/sec'
                case('scaled_energy_rate')
                  sub_condition_ptr%itype = SCALED_ENERGY_RATE_SS
                  energy_rate_unit_string = 'MW|MJ/sec'
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
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // '" ' // trim(string)
                call printErrMsg(option)
                endif
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = CONDUCTANCE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('seepage')
              sub_condition_ptr%itype = SEEPAGE_BC
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('unit_gradient')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may &
                  &only be associated with a PRESSURE flow condition.'
                call printErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC
            case('heterogeneous_volumetric_rate')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('heterogeneous_mass_rate')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('heterogeneous_dirichlet')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('heterogeneous_seepage')
              sub_condition_ptr%itype = HET_SEEPAGE_BC
            case('heterogeneous_conductance')
              sub_condition_ptr%itype = HET_CONDUCTANCE_BC
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
          internal_units = 'not_assigned'
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
        internal_units = rate_unit_string
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
        internal_units = energy_rate_unit_string
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
           SURFACE_SPILLOVER,HET_DIRICHLET_BC,HET_SEEPAGE_BC,&
           HET_CONDUCTANCE_BC)
        option%io_buffer = 'RATE condition must not be of type: dirichlet, &
          &neumann, zero_gradient, dirichlet_zero_gradient, hydrostatic, &
          &seepage, or conductance".'
        call printErrMsg(option)
    end select
  endif
  ! check to ensure that a pressure condition is not of type rate
  if (associated(pressure)) then
    select case(pressure%itype)
      case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS)
        option%io_buffer = 'PRESSURE or FLUX condition must not be of type: &
          &mass_rate, scaled_mass_rate, volumetric_rate, &
          &scaled_volumetric_rate, equilibrium, or production_well.'
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

    case(TH_MODE)
      if (.not.associated(pressure) .and. .not.associated(rate)&
           .and. .not.associated(well) .and. .not.associated(saturation)) then
        option%io_buffer = 'pressure, rate and saturation condition null in &
                           &condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif

      if (.not.associated(temperature) .and. .not.associated(energy_rate) &
          .and. .not.associated(energy_flux)) then
        option%io_buffer = 'temperature, energy_flux, and energy_rate &
          &condition null in condition: ' // trim(condition%name)
        call printErrMsg(option)
      endif
      if (associated(temperature) .and. associated(energy_rate) ) then
        option%io_buffer = 'Both, temperature and energy_rate cannot be &
                           &specified in condition: ' // trim(condition%name)
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
    

  end select

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionRead

! ************************************************************************** !

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

subroutine ConditionReadValues(input,option,keyword,dataset_base, &
                               data_external_units,data_internal_units)
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
  character(len=*) :: data_external_units
  character(len=*) :: data_internal_units

  character(len=MAXSTRINGLENGTH), pointer :: internal_unit_strings(:)
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i, icount
  PetscInt :: icol
  PetscInt :: ndims
  PetscInt, pointer :: dims(:)
  PetscReal, pointer :: real_buffer(:)
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
      ' in the input file is already associated with a different dataset &
      &type.  Check for duplicate definitions of ' // trim(keyword) // '.'
    call printErrMsg(option)
  endif

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
    if (length == FOUR_INTEGER .and. &
        StringCompare(word,'file',FOUR_INTEGER)) then 
      input%err_buf2 = trim(keyword) // ', FILE'
      input%err_buf = 'keyword'
      call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
      if (input%ierr == 0) then
        filename = string2
      else
        option%io_buffer = 'The ability to read realization dependent &
          &datasets outside the DATASET block is no longer supported'
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
          allocate(flow_dataset%time_series%values(flow_dataset% &
                                                     time_series%rank,dims(2))) 
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
          option%io_buffer = 'HDF condition data set rank does not match &
            &rank of internal data set.  Email Glenn for additions'
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
        error_string = 'CONDITION,' // trim(keyword) // ',FILE'
        call DatasetAsciiReadFile(dataset_ascii,filename,data_external_units, &
                                  data_internal_units,error_string,option)
        dataset_ascii%filename = filename
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
      error_string = 'CONDITION,' // trim(keyword) // ',LIST'
      call DatasetAsciiReadList(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
    else
      option%io_buffer = 'Keyword "' // trim(word) // &
        '" not recognized in when reading condition values for "' // &
        trim(keyword) // '".'
      call printErrMsg(option)
    endif
  else
    input%buf = trim(string2)
    error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
    call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
#if 0
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
#endif
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
    case(HET_DIRICHLET_BC)
      string = 'heterogeneous dirichlet'
    case(HET_SEEPAGE_BC)
      string = 'heterogeneous seepage'
    case(HET_CONDUCTANCE_BC)
      string = 'heterogeneous conductance'
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

subroutine FlowConditionUpdate(condition_list,option)
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

  type(flow_condition_type), pointer :: condition
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    call DatasetUpdate(condition%datum,option)
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr

      if (associated(sub_condition)) then
        call DatasetUpdate(sub_condition%dataset,option)
        call DatasetUpdate(sub_condition%gradient,option)
      endif

    enddo

    condition => condition%next

  enddo

end subroutine FlowConditionUpdate

! ************************************************************************** !

subroutine TranConditionUpdate(condition_list,option)
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

  type(tran_condition_type), pointer :: condition

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit

    do
      if (associated(condition%cur_constraint_coupler%next)) then
        if (option%time >= condition%cur_constraint_coupler%next%time) then
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
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%energy_rate) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) &
      ) then
    FlowConditionIsTransient = PETSC_TRUE
  endif

end function FlowConditionIsTransient

! ************************************************************************** !


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
  call FlowSubConditionDestroy(condition%temperature)
  call FlowSubConditionDestroy(condition%concentration)
  call FlowSubConditionDestroy(condition%enthalpy)
  call FlowSubConditionDestroy(condition%energy_rate)

  call TimeStorageDestroy(condition%default_time_storage)

  nullify(condition%next)

  deallocate(condition)
  nullify(condition)

end subroutine FlowConditionDestroy


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
