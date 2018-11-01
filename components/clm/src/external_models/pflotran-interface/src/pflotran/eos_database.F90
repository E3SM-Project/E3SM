module EOSData_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none

  private

  PetscInt, parameter, public :: EOS_DENSITY = 1
  PetscInt, parameter, public :: EOS_ENTHALPY = 2
  PetscInt, parameter, public :: EOS_VISCOSITY = 3
  PetscInt, parameter, public :: EOS_INTERNAL_ENERGY = 4
  PetscInt, parameter, public :: EOS_FVF = 5
  PetscInt, parameter, public :: EOS_RS = 6
  PetscInt, parameter, public :: EOS_COMPRESSIBILITY = 7
  PetscInt, parameter, public :: EOS_VISCOSIBILITY = 8
  !add here other properties and derivatives, then increase EOS_MAX_PROP_NUM
  !and add the default user and internal units in EOSDataBaseInit
  PetscInt, parameter, public :: EOS_MAX_PROP_NUM = 8

  !type, public :: eos_data_base_type
  type, abstract, public :: eos_data_base_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscInt :: num_p   ! number of pressure intervals
    PetscInt :: num_t   ! number of temperature points
    PetscInt :: num_prop ! number of properties in the database
    PetscReal :: press_unit_conv_factor
    character(len=MAXWORDLENGTH) :: press_internal_units
    character(len=MAXWORDLENGTH) :: press_user_units
    PetscReal :: temp_unit_conv_factor
    character(len=MAXWORDLENGTH) :: temp_internal_units
    character(len=MAXWORDLENGTH) :: temp_user_units
  contains
    procedure :: EOSDataBaseInit
    procedure, public :: EOSPropPresent
    procedure :: ReadUserUnits
    procedure :: UnitConversionFactors
    procedure :: ConvertFVFtoMolarDensity
    procedure :: EOSDataBaseStrip
    procedure :: EOSSetConstGradExtrap
    procedure, public :: SetDefaultInternalUnits
    procedure, public :: SetMetricUnits
    procedure, public :: AddEOSProp
    procedure, public :: SetupVarLinLogInterp
  end type

  type, public, extends(eos_data_base_type) :: eos_database_type
    character(len=MAXWORDLENGTH) :: file_name
    PetscReal :: dp      ! uniform pressure interval
    PetscReal :: dt      ! uniform temperature interval
    class(lookup_table_uniform_type), pointer :: lookup_table_uni
  contains
    procedure, public :: Read => EOSDatabaseRead
    procedure, public :: EOSProp => EOSPropLinearInterp
    procedure, public :: EOSPropGrad => EOSPropGradDatabase
  end type

  type, public, extends(eos_data_base_type) :: eos_table_type
    PetscReal :: temperature !temperature value for isothermal data
    PetscInt :: n_indices !number of indices to be saved for lookup
    PetscInt :: first_index !location of first index in auxvars
    class(lookup_table_general_type), pointer :: lookup_table_gen
    class(eos_table_type), pointer :: next => null()
  contains
    procedure, public :: Read => EOSTableRead
    procedure, public :: EOSProp => EOSPropTable
    procedure, public :: EOSPropGrad => EOSPropGradTable
  end type

  type, public :: eos_table_list_type
    PetscInt :: num_eos_tables
    class(eos_table_type), pointer :: first
    class(eos_table_type), pointer :: last
    class(eos_table_type), pointer :: array(:)
  end type eos_table_list_type

  type(eos_table_list_type), pointer, public :: eos_table_list => null()

  public :: EOSDatabaseCreate, &
            EOSDatabaseDestroy, &
            EOSTableCreate, &
            EOSTableDestroy, &
            EOSTableInitList, &
            EOSTableAddToList, &
            EOSTableProcessList, &
            EOSTableDestroyList



contains

! ************************************************************************** !

subroutine EOSDataBaseInit(this)

  implicit none

  class(eos_data_base_type) :: this


  PetscInt :: i_prop
  this%name = ''
  this%id = UNINITIALIZED_INTEGER
  this%num_p = UNINITIALIZED_INTEGER
  this%num_t = UNINITIALIZED_INTEGER
  this%num_prop = UNINITIALIZED_INTEGER
  this%press_unit_conv_factor = 1.0d0
  this%temp_unit_conv_factor = 1.0d0
  this%press_internal_units = 'Pa'
  this%temp_internal_units = 'C'
  this%press_user_units = 'Pa'
  this%temp_user_units = 'C'

end subroutine EOSDataBaseInit

! ************************************************************************** !

function EOSPropPresent(this,prop_iname)
  !
  ! Author: Paolo Orsini
  ! Date: 12/18/15
  !
  ! Checks if a property is defined in the database

  implicit none

  class(eos_data_base_type) :: this
  PetscInt, intent(in) :: prop_iname
  PetscBool :: EOSPropPresent

  select type(this)
    class is(eos_database_type)
      EOSPropPresent = &
              associated(this%lookup_table_uni%var_array(prop_iname)%ptr)
    class is(eos_table_type)
      EOSPropPresent = &
            associated(this%lookup_table_gen%var_array(prop_iname)%ptr) 
  end select


end function EOSPropPresent

! ************************************************************************** !

subroutine ReadUserUnits(this,input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Read unit EOS data unit and compute convert them to internal

  use Input_Aux_module
  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, internal_units, user_units
  character(len=MAXSTRINGLENGTH) :: error_string
  type(lookup_table_var_ptr_type), pointer :: prop_array(:)
  PetscInt :: prop_idx
  prop_idx = 0
  
  select type(this)
    class is(eos_database_type)
      prop_array =>  this%lookup_table_uni%var_array
    class is(eos_table_type)
      prop_array =>  this%lookup_table_gen%var_array
  end select

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    select case(keyword)
      case('PRESSURE')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        this%press_user_units = trim(user_units)
      case('TEMPERATURE')
        !only Celsius currently supported - currently unit must be defined
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          if (trim(word) /= "C") then
            option%io_buffer = "EOS data only Temperatures in Celcius " // &
                                "are supported"
            call printErrMsg(option)
          end if
        end if
        this%temp_unit_conv_factor = 1.0
        this%temp_user_units = 'C'
      case('DENSITY')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        prop_array(EOS_DENSITY)%ptr%user_units = trim(user_units)
      case('ENTHALPY')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        prop_array(EOS_ENTHALPY)%ptr%user_units = trim(user_units)        
      case('INTERNAL_ENERGY')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        prop_array(EOS_INTERNAL_ENERGY)%ptr%user_units = trim(user_units)
      case('VISCOSITY')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        prop_array(EOS_VISCOSITY)%ptr%user_units = trim(user_units)
      case('FVF')
        call InputReadWord(input,option,user_units,PETSC_TRUE)
        prop_array(EOS_FVF)%ptr%user_units = trim(user_units)
      case('RS')
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         prop_array(EOS_RS)%ptr%user_units = trim(user_units)
      case('COMPRESSIBILITY')
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         prop_array(EOS_COMPRESSIBILITY)%ptr%user_units = trim(user_units)
      case('VISCOSIBILITY')
         call InputReadWord(input,option,user_units,PETSC_TRUE)
         prop_array(EOS_VISCOSIBILITY)%ptr%user_units = trim(user_units)
      case default
        error_string = trim(error_string) // ': ' // this%name // &
        ': EOS DATA units'
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  end do

  call this%UnitConversionFactors(option)
  nullify(prop_array)

end subroutine ReadUserUnits

! ************************************************************************** !
subroutine UnitConversionFactors(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17 - modified 04/18/2018
  !
  ! Process user input units and compute properties conversion fator

  use Option_module
  use Units_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option

  PetscInt :: i_prop, data_idx

  !covert pressure
  this%press_unit_conv_factor = UnitsConvertToInternal(this%press_user_units, &
                                  this%press_internal_units, option)

  !no temperature conversion - currently allowed only C
  
  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      call this%lookup_table_uni%LookupTableVarConvFactors(option)
    class is(eos_table_type)
      call this%lookup_table_gen%LookupTableVarConvFactors(option)
  end select    

end subroutine UnitConversionFactors

! ************************************************************************** !

subroutine SetDefaultInternalUnits(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17
  !
  ! Set default internal units - to read pvt tables

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option

  type(lookup_table_var_ptr_type), pointer :: prop_array(:) => null()
  PetscInt :: prop_idx

  this%press_internal_units = 'Pa'
  this%temp_internal_units = 'C'
 
  nullify(prop_array)

  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      prop_array => this%lookup_table_uni%var_array
    class is(eos_table_type)
      prop_array => this%lookup_table_gen%var_array
  end select
 
  do prop_idx = 1,size(prop_array(:))
    if ( associated(prop_array(prop_idx)%ptr) ) then
      select case (prop_array(prop_idx)%ptr%iname)
        case(EOS_DENSITY)
          prop_array(prop_idx)%ptr%internal_units = 'kg/m^3'
        case(EOS_ENTHALPY)
          prop_array(prop_idx)%ptr%internal_units = 'J/kg'
        case(EOS_VISCOSITY)
          prop_array(prop_idx)%ptr%internal_units = 'Pa.s'
        case(EOS_INTERNAL_ENERGY)
          prop_array(prop_idx)%ptr%internal_units = 'J/kg'
        case(EOS_FVF)
          prop_array(prop_idx)%ptr%internal_units = 'm^3/m^3'
        case(EOS_RS)
          prop_array(prop_idx)%ptr%internal_units = 'm^3/m^3'
        case(EOS_COMPRESSIBILITY)
          prop_array(prop_idx)%ptr%internal_units = '1/Pa'
        case(EOS_VISCOSIBILITY)
          prop_array(prop_idx)%ptr%internal_units = '1/Pa'
      end select
    end if
  end do

  nullify(prop_array)

end subroutine SetDefaultInternalUnits

! ************************************************************************** !

subroutine SetMetricUnits(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17
  !
  ! Set METRIC user units of a pvt tables

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option

  type(lookup_table_var_ptr_type), pointer :: prop_array(:) => null()
  PetscInt :: prop_idx

  this%press_user_units = 'Bar'
  this%temp_user_units = 'C'

  nullify(prop_array)
  
  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      prop_array => this%lookup_table_uni%var_array
    class is(eos_table_type)
      prop_array => this%lookup_table_gen%var_array
  end select
  
  do prop_idx = 1,size(prop_array(:))
    if ( associated(prop_array(prop_idx)%ptr) ) then
      select case (prop_array(prop_idx)%ptr%iname)
        case(EOS_DENSITY)
          prop_array(prop_idx)%ptr%user_units = 'kg/m^3'
        case(EOS_ENTHALPY)
          prop_array(prop_idx)%ptr%user_units = 'kJ/kg'
        case(EOS_VISCOSITY)
          prop_array(prop_idx)%ptr%user_units = 'cP'
        case(EOS_INTERNAL_ENERGY)
          prop_array(prop_idx)%ptr%user_units = 'kJ/kg'
        case(EOS_FVF)
          prop_array(prop_idx)%ptr%user_units = 'm^3/m^3'
        case(EOS_RS)
          prop_array(prop_idx)%ptr%user_units = 'm^3/m^3'
        case(EOS_COMPRESSIBILITY)
          prop_array(prop_idx)%ptr%user_units = '1/Bar'
        case(EOS_VISCOSIBILITY)
          prop_array(prop_idx)%ptr%user_units = '1/Bar'
      end select
    end if
  end do  

  nullify(prop_array)

  call this%UnitConversionFactors(option)

end subroutine SetMetricUnits

! ************************************************************************** !

subroutine AddEOSProp(this,new_var,option)
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/18
  !
  ! Add new var eos properties to eos data (to variable list) 

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(lookup_table_var_type), pointer :: new_var
  type(option_type) :: option

  type(lookup_table_var_ptr_type), pointer :: var_array(:) => null()
  type(lookup_table_var_list_type), pointer :: vars => null()
  character(len=MAXSTRINGLENGTH) :: string 
  character(len=MAXWORDLENGTH) :: word_error

  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      vars => this%lookup_table_uni%vars
      var_array => this%lookup_table_uni%var_array
    class is(eos_table_type)
      vars => this%lookup_table_gen%vars
      var_array => this%lookup_table_gen%var_array
  end select

  string = 'Adding Property to eos database or table'

  !if( Uninitialized(new_var%extrapolation_itype) ) then
  !  option%io_buffer = UninitializedMessage('Extrapolation type',string)
  !  call printErrMsg(option)
  !end if

  if( Uninitialized(new_var%iname) ) then
    option%io_buffer = UninitializedMessage('EOS property iname',string)
    call printErrMsg(option)
  end if

  if (.not.EOSPropExistInDictionary(new_var%iname)) then
    write(word_error,*) new_var%iname
    option%io_buffer = 'EOS property with iname = ' // word_error
    call printErrMsg(option)
  end if

  call LookupTableVarAddToList(new_var,vars)
    
  var_array(new_var%iname)%ptr => new_var
  
  nullify(vars)
  nullify(var_array)

end subroutine AddEOSProp

! ************************************************************************** !

subroutine EOSSetConstGradExtrap(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 06/02/18
  !
  ! Set constant gradient extrapolation method for all lookup variables

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option
  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      call this%lookup_table_uni%SetupConstGradExtrap(option)
    class is(eos_table_type)
      call this%lookup_table_gen%SetupConstGradExtrap(option)
  end select    

end subroutine EOSSetConstGradExtrap

! ************************************************************************** !

function EOSPropExistInDictionary(property_iname)
  !
  ! Author: Paolo Orsini
  ! Date: 05/05/18
  !
  ! Check if the EOS property integer name exist in the dictionary
  ! defined at the beginning of the module 

  implicit none

  PetscBool :: EOSPropExistInDictionary
  PetscInt, intent(in) :: property_iname
  
  EOSPropExistInDictionary = .false.
  
  select case(property_iname)
    case(EOS_DENSITY)
      EOSPropExistInDictionary = .true.
    case(EOS_ENTHALPY)
      EOSPropExistInDictionary = .true.
    case(EOS_VISCOSITY)
      EOSPropExistInDictionary = .true.
    case(EOS_INTERNAL_ENERGY)
      EOSPropExistInDictionary = .true.
    case(EOS_FVF)
      EOSPropExistInDictionary = .true.
    case(EOS_RS)
      EOSPropExistInDictionary = .true.      
    case(EOS_COMPRESSIBILITY)
      EOSPropExistInDictionary = .true.
    case(EOS_VISCOSIBILITY)
      EOSPropExistInDictionary = .true.      
  end select

end function EOSPropExistInDictionary

! ************************************************************************** !

subroutine SetupVarLinLogInterp(this,var_iname,option)
  !
  ! Author: Paolo Orsini
  ! Date: 05/04/18
  !
  ! Add new var eos properties to eos data (to variable list) 

  use Option_module

  implicit none

  class(eos_data_base_type) :: this
  type(option_type) :: option
  PetscInt, intent(in) :: var_iname

select type(this)
  class is(eos_database_type)
    call this%lookup_table_uni%SetupVarLinLogInterp(var_iname,option)
  class is(eos_table_type)
    call this%lookup_table_gen%SetupVarLinLogInterp(var_iname,option)
end select

end subroutine SetupVarLinLogInterp

! ************************************************************************** !

subroutine ConvertFVFtoMolarDensity(this,FMW,reference_density_kg)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Replaces FVF with molar density

  implicit none

  class(eos_data_base_type) :: this
  PetscReal, intent(in) :: FMW
  PetscReal, intent(in) :: reference_density_kg
 
  PetscInt :: data_idx
  PetscReal, pointer :: var_data(:,:) => null()
  type(lookup_table_var_ptr_type), pointer :: var_array(:) => null()

  !EOS properties conversion factors 
  select type(this)
    class is(eos_database_type)
      var_array => this%lookup_table_uni%var_array
      var_data => this%lookup_table_uni%var_data
    class is(eos_table_type)
      var_array => this%lookup_table_gen%var_array
      var_data => this%lookup_table_gen%var_data
  end select  
 
  data_idx = var_array(EOS_FVF)%ptr%data_idx
  !          kmol/rm^3 = kg/sm^3 * kmol/kg * sm^3/rm^3
  var_data(data_idx,:) = reference_density_kg / FMW / var_data(data_idx,:)
  !from this point on in the data_idx there is not FVF but EOS_DENSITY
  !modify variable pointig this
  var_array(EOS_FVF)%ptr%iname = EOS_DENSITY
  var_array(EOS_FVF)%ptr%internal_units = 'kmol/m^3'
  var_array(EOS_FVF)%ptr%user_units = 'kmol/m^3'
  var_array(EOS_FVF)%ptr%conversion_factor = 1.0
  !change pointer
  var_array(EOS_DENSITY)%ptr => var_array(EOS_FVF)%ptr
  !FVF not available any longer
  nullify(var_array(EOS_FVF)%ptr)
  
  nullify(var_array)
  nullify(var_data)

end subroutine ConvertFVFtoMolarDensity

! ************************************************************************** !

subroutine EOSDataBaseStrip(this)

  use Utility_module

  implicit none

  class(eos_data_base_type) :: this

end subroutine EOSDataBaseStrip

! ************************************************************************** !

function EOSDatabaseCreate(filename,dbase_name)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !

  implicit none

  class(eos_database_type), pointer :: EOSDatabaseCreate
  character(len=MAXWORDLENGTH) :: filename
  character(len=*) :: dbase_name
  
  PetscInt :: i_var

  allocate(EOSDatabaseCreate)
  call EOSDatabaseCreate%EOSDataBaseInit()
  EOSDatabaseCreate%name = trim(dbase_name)
  EOSDatabaseCreate%file_name = trim(filename)
  nullify(EOSDatabaseCreate%lookup_table_uni)

  !initialize new data structure
  EOSDatabaseCreate%lookup_table_uni => LookupTableCreateUniform(TWO_INTEGER)
  call EOSDatabaseCreate%lookup_table_uni%LookupTableVarsInit(EOS_MAX_PROP_NUM)
  

end function EOSDatabaseCreate

! ************************************************************************** !

subroutine EOSDatabaseRead(this,option)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15 - modified 04/18/2018
  !
  ! Reads the the an EOS database from a text file for one or more
  ! phase properties.
  !
  ! Database format (for look up table)
  ! Header: NUM_P, NUM_T and DATA_LIST_ORDER
  ! DATA
  ! column 1 = pressure
  ! column 2 = temperature
  ! column 2+1, column 2+n: properties in the order given in DATA_LIST_ORDER
  !
  ! The dataset must be followed by NUM_P * NUM_T lines, not interrupted
  ! by a commented line.
  ! Each line must include P,T, and all properties listed in DATA_LIST_ORDER
  !
  ! each property listed depends on P and T, prop(P,T).
  !
  ! temperature values must be equispaced (same dt in the entire table)
  ! pressure values must be equispaced (same dp in the entire table)
  ! The data must be ordered so for growing P and T, with T looping faster.
  !
  ! This format repeates unnecessary P and T values in favours o readibility
  ! and creation of small dataset by hand. Inefficient for for large datasets.
  !

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(eos_database_type) :: this
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word, internal_units, user_units
  character(len=MAXSTRINGLENGTH) :: error_string = 'EOS_DATABASE'
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: prop_idx, prop_count, i_idx, j_idx
  PetscInt :: data_size
  PetscInt :: data_idx
  PetscReal :: tempreal
  PetscBool :: pres_present, temp_present
  type(lookup_table_var_type), pointer :: db_var => null()

  type(input_type), pointer :: input_table


  if (len_trim(this%file_name) < 1) then
    option%io_buffer = 'FILENAME must be specified for EOS_DATABASE.'
    call printErrMsg(option)
  endif

  input_table => InputCreate(IUNIT_TEMP,this%file_name,option)
  input_table%ierr = 0
  input_table%force_units = PETSC_FALSE

  !if ( option%myrank == 0 ) then
    option%io_buffer = 'Reading database = ' // this%file_name
    call printMsg(option)
  !end if

  !set deafult user units - identical to internal units
  this%press_user_units = 'Pa'
  this%temp_user_units = 'C'  

  !reading the database file header
  do

    call InputReadPflotranString(input_table,option)
    !if (InputCheckExit(input_table,option)) exit
    if (InputError(input_table)) exit

    call InputReadWord(input_table,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input_table,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(keyword)
      case('NUM_P')
        call InputReadInt(input_table,option,this%num_p)
        call InputErrorMsg(input_table,option,'number of dp',error_string)
      case('NUM_T')
        call InputReadInt(input_table,option,this%num_t)
        call InputErrorMsg(input_table,option,'number of dt',error_string)
      case('DATA_LIST_ORDER')
        pres_present = PETSC_FALSE; 
        temp_present = PETSC_FALSE;
        prop_count = 0
        do
          call InputReadPflotranString(input_table,option)
          if (InputCheckExit(input_table,option)) exit
          call InputReadWord(input_table,option,word,PETSC_TRUE)
          select case(word)
            case('PRESSURE')
              pres_present = PETSC_TRUE
            case('TEMPERATURE')
              temp_present = PETSC_TRUE
            case('DENSITY')
              prop_count = prop_count + 1
              internal_units = 'kg/m^3'
              user_units = 'kg/m^3'
              db_var => CreateLookupTableVar(EOS_DENSITY,internal_units, &
                                  user_units,prop_count)
              call LookupTableVarAddToList(db_var,this%lookup_table_uni%vars)
              this%lookup_table_uni%var_array(EOS_DENSITY)%ptr => db_var              
            case('ENTHALPY')
              prop_count = prop_count + 1
              internal_units = 'J/kg'
              user_units = 'J/kg'
              db_var => CreateLookupTableVar(EOS_ENTHALPY,internal_units, &
                                  user_units,prop_count)
              call LookupTableVarAddToList(db_var,this%lookup_table_uni%vars)
              this%lookup_table_uni%var_array(EOS_ENTHALPY)%ptr => db_var
            case('INTERNAL_ENERGY')
              prop_count = prop_count + 1
              internal_units = 'J/kg'
              user_units = 'J/kg'
              db_var => CreateLookupTableVar(EOS_INTERNAL_ENERGY, & 
                                        internal_units,user_units,prop_count)
              call LookupTableVarAddToList(db_var,this%lookup_table_uni%vars)
              this%lookup_table_uni%var_array(EOS_INTERNAL_ENERGY)%ptr => &
                                                                   db_var
            case('VISCOSITY')
              prop_count = prop_count + 1
              internal_units = 'Pa-s'
              user_units = 'Pa-s'
              db_var => CreateLookupTableVar(EOS_VISCOSITY,internal_units, &
                                             user_units,prop_count)
              call LookupTableVarAddToList(db_var,this%lookup_table_uni%vars)
              this%lookup_table_uni%var_array(EOS_VISCOSITY)%ptr => db_var
            case default
              error_string = trim(error_string) // ': ' // this%file_name // &
              ': DATA_LIST_ORDER'
              call InputKeywordUnrecognized(keyword,error_string,option)
          end select
        end do
        this%num_prop = prop_count
        ! go back to DATA_LIST_ORDER - read variable again to comput
        ! unit covnertion factors
        string = "DATA_LIST_ORDER"
        call InputFindStringInFile(input_table,option,string)
        call this%ReadUserUnits(input_table,option)

        if (.not.pres_present) then
          option%io_buffer = 'PRESSURE must be present in any EOS_DATABASE.'
          call printErrMsg(option)
        end if
        if (.not.temp_present) then
          option%io_buffer = 'TEMPERATURE must be present in any EOS_DATABASE.'
          call printErrMsg(option)
        end if
      case('DATA')
        exit
      case default
        error_string = trim(error_string) // ': ' // this%file_name
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select

  end do

  nullify(db_var)

  data_size = this%num_p * this%num_t

  this%num_prop = prop_count
  allocate(this%lookup_table_uni%var_data(prop_count,data_size))

  !two-dims lookup table created in EOSDatabaseCreate
  this%lookup_table_uni%dims(1) = this%num_t
  this%lookup_table_uni%dims(2) = this%num_p
  allocate(this%lookup_table_uni%axis1%values(this%num_t))
  this%lookup_table_uni%axis1%values(1:this%num_t) = UNINITIALIZED_DOUBLE
  allocate(this%lookup_table_uni%axis2%values(this%num_p))
  this%lookup_table_uni%axis2%values(1:this%num_p) = UNINITIALIZED_DOUBLE

  !TODO
  ! start loading data - at the moment using Input facility, however this file
  ! can be large. TODO. Implement more efficient solutions:
  ! - using read(,) without InputReadDouble filter
  ! - adding the option of reading a .h5 file where the database is defined

  ! go to data - first time to load axis1 and 2 values
  string = "DATA"
  call InputFindStringInFile(input_table,option,string)

  do j_idx = 1,this%num_p

    do i_idx = 1, this%num_t
      call InputReadPflotranString(input_table,option)

      call InputReadDouble(input_table,option,tempreal)
      call InputErrorMsg(input_table,option, &
                           'VALUE', 'EOS_DATABASE PRESS_VALUE')
      ! convert pressure to internal units - Pa
      this%lookup_table_uni%axis2%values(j_idx) = &
          tempreal * this%press_unit_conv_factor
      !reads temperature values - repeated this%num_p times - not efficient
      call InputReadDouble(input_table,option, &
                           this%lookup_table_uni%axis1%values(i_idx))
      call InputErrorMsg(input_table,option, &
                         'VALUE', 'EOS_DATABASE TEMP_VALUE')

      prop_count = i_idx + (j_idx-1) * this%num_t
      do prop_idx = 1,this%num_prop
        call InputReadDouble(input_table,option, &
                        this%lookup_table_uni%var_data(prop_idx,prop_count))        
        call InputErrorMsg(input_table,option,&
                           'VALUE','EOS_DATABASE PROP_VALUE')
      end do

    end do

  end do

  !unit covnersions
  call this%lookup_table_uni%VarPointAndUnitConv(option)
  !setup constant gradient extrapolation
  call this%EOSSetConstGradExtrap(option)
  !allocate lookup var gradients 
  call this%lookup_table_uni%LookupTableVarInitGradients(option)

  call InputDestroy(input_table)

end subroutine EOSDatabaseRead

! ************************************************************************** !

subroutine EOSPropLinearInterp(this,T,P,prop_iname,prop_value,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 12/12/15
  !
  ! interpolates a single EOS property from the EOS database
  ! Note: when more properties must be extracted from the same EOS database,
  !       i.e. properties listed for the same values and range of P and T,
  !       all propoertis should be extracted at the same time to perform
  !       only once the look up operations, which are:
  !       (a) LookupTableIndexUniform, i.e. axis look up
  !       (b) P,T location checks within a single 2D domain, (p,p+dp; t,t+dt)
  !           currently done within LookupTableInterpolate2DUniform
  !
  !       TODO: add a method lookup_table_uni to extract mulitdimensional data
  !             at the same time (data plus function)

  implicit none

  class(eos_database_type) :: this
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: prop_value ! database units (SI)
  PetscErrorCode, intent(out) :: ierr

  !PetscReal :: EOSEOSProp !property given in units specified in the database

  ierr = 0
  !insert check P and/or T out of range
  ! T smaller then minimum T in the database
  ! errors here a passed back to the EOSxxxPhase module
  ! and to ModeXXXComputeAuxVar where cells ids are available to check where
  ! P and T are out of range.
  ! TODO: (i) define error id identifiers for each EOSPhase, or for all EOSs.
  !       (ii) define a function that handles errors ids and print error message.
  !       (iii) can then remove the print statment below
  !       note: need to minimise the number if-checks for efficiency.
  if ( T < this%lookup_table_uni%axis1%values(1) ) then
    ierr = 101
    print*, "EOSEOSProp - T smaller than min val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop
  end if
  if ( T > this%lookup_table_uni%axis1%values(this%num_t) ) then
    ierr = 102
    print*, "EOSEOSProp - T larger than max val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop
  end if
  if ( P < this%lookup_table_uni%axis2%values(1) ) then
    ierr = 103
    print*, "EOSEOSProp - P smaller than min val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop
  end if
  if ( P > this%lookup_table_uni%axis2%values(this%num_p) ) then
    ierr = 104
    print*, "EOSEOSProp - P larger than max val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop
  end if

  this%lookup_table_uni%data => &
                    this%lookup_table_uni%var_array(prop_iname)%ptr%data
                             !T       P      optional
  !this%lookup_table_uni%Sample(lookup1,lookup2,lookup3)
  prop_value = this%lookup_table_uni%Sample(T,P)

  nullify(this%lookup_table_uni%data)

end subroutine EOSPropLinearInterp

! ************************************************************************** !

subroutine EOSPropGradDatabase(this,T,P,prop_iname,sample,dSampledT, &
                               dSampledP,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 06/04/18
  !
  ! interpolates a single EOS property from the EOS database 
  ! and compute its gradient with respect to the P and T of the database 
  ! (i.e. raw table gradients)

  implicit none

  class(eos_database_type) :: this
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: sample ! PFLOTRAN intenral units (SI)
  PetscReal, intent(out) :: dSampledT ! internal PFLOTRAN units ([Sample]/[C])
  PetscReal, intent(out) :: dSampledP ! internal PFLOTRAN units ([Sample]/[Pa])  
  PetscErrorCode, intent(out) :: ierr

  ierr = 0

  !bound checks to be included later as part of reporting
  !go trgough the table list and save the outbound in a report data strucutre?
  !could be done in auxvra compute using a reusable table list specific func

  call this%lookup_table_uni%SampleAndGradient(prop_iname,T,P)
  
  sample = this%lookup_table_uni%var_array(prop_iname)%ptr%sample
  dSampledT = this%lookup_table_uni%var_array(prop_iname)%ptr%sample_grad(1)
  dSampledP = this%lookup_table_uni%var_array(prop_iname)%ptr%sample_grad(2)

end subroutine EOSPropGradDatabase

! ************************************************************************** !

subroutine EOSDatabaseDestroy(eos_database)
  !
  ! Author: Paolo Orsini
  ! Date: 12/14/15
  !
  ! destroys EOS database

  use Utility_module

  implicit none

  class(eos_database_type), pointer :: eos_database

  if (.not.associated(eos_database)) return

  call eos_database%EOSDataBaseStrip()
  call LookupTableDestroy(eos_database%lookup_table_uni)

  deallocate(eos_database)
  nullify(eos_database)

end subroutine EOSDatabaseDestroy

! ************************************************************************** !

function EOSTableCreate(table_name,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Create EOS table

  use Option_module

  implicit none

  class(eos_table_type), pointer :: EOSTableCreate
  character(len=*) :: table_name
  type(option_type) :: option

  allocate(EOSTableCreate)
  call EOSTableCreate%EOSDataBaseInit()
  EOSTableCreate%name = trim(table_name)
  EOSTableCreate%n_indices = 0
  EOSTableCreate%first_index = 0
  nullify(EOSTableCreate%lookup_table_gen)

  !create table without dimensions 
  EOSTableCreate%lookup_table_gen => LookupTableCreateGeneral()
  call EOSTableCreate%lookup_table_gen%LookupTableVarsInit(EOS_MAX_PROP_NUM)

end function EOSTableCreate

! ************************************************************************** !

subroutine EOSTableRead(this,input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Reads the PVT table data

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none

  class(eos_table_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string = 'EOS_PVT_TABLE'
  class(lookup_table_axis_type), pointer :: axis1
  PetscReal, pointer :: var_data(:,:)
  PetscReal, pointer :: press_data_array(:,:)
  PetscInt :: n_temp_count
  PetscInt :: i_temp, i_press, i_prop
  PetscInt :: n_press_count, n_press_count_first, n_press_count_cur
  PetscInt :: temp_array_size, press_array_size
  PetscReal, pointer :: temp_array(:)

  n_press_count_first = 0
  n_temp_count = 0
  n_press_count = 0
  n_press_count_cur = 0
  temp_array_size = 100 !estimate for temperature points
  allocate(temp_array(temp_array_size))
  temp_array = UNINITIALIZED_DOUBLE
  press_array_size = 1000 !estimate of pressure points from all tables
  ! num_prop+1 = + 1 make space for the pressure
  allocate(press_data_array(this%num_prop+1,press_array_size))
  press_data_array = UNINITIALIZED_DOUBLE

  !reading pvt table
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    if (InputError(input)) then
      option%io_buffer = 'Error found as reading PVT table'
      call printErrMsg(option)
    end if
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(keyword)
      case('DATA_UNITS')
        call this%ReadUserUnits(input,option)
      case('DATA')
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit
          if (InputError(input)) then
            option%io_buffer = 'Error found as reading PVT table - DATA block'
            call printErrMsg(option)
          end if
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'word',error_string)
          call StringToUpper(keyword)
          select case(word)
            case('TEMPERATURE')
              n_temp_count = n_temp_count + 1
              if (n_temp_count > temp_array_size) then
                call reallocateRealArray(temp_array,temp_array_size)
              end if
              call InputReadDouble(input,option,temp_array(n_temp_count))
              call InputErrorMsg(input,option,&
                                 'VALUE','EOS_PVT_TABLE - TEMPERATURE')
              !read pvt data for a given temperature
              !call ReadPressureTable(input,option,this%num_prop,&
              !                       n_press_count_cur,data_array)
              call ReadPressureTable(input,option,n_press_count, &
                                     n_press_count_cur,press_data_array)
              if (n_temp_count == 1) then
                n_press_count_first = n_press_count_cur
              else
                if ( n_press_count_cur /= n_press_count_first ) then
                  option%io_buffer =  'PVT Table = ' // this%name // &
                                      "PVT tables with pressure points that" // &
                                      "varies with T not currently supported"
                  call printErrMsg(option)
                end if
              end if
            case default
              option%io_buffer = 'PVT Table DATA block must contain ' // &
                                 'TEMPERATURE blocks'
              call printErrMsg(option)
          end select
        end do
      case default
        option%io_buffer = 'PVT Table = ' // 'this%name'
    end select
  end do

  this%num_t = n_temp_count
  this%num_p = n_press_count

  if (this%num_t ==1) then !isothermal press_data_array
    this%temperature = temp_array(1)
    ! create general 1D lookup table
    !this%lookup_table_gen => LookupTableCreateGeneral(ONE_INTEGER) 
    allocate(axis1)   
    allocate(axis1%values(this%num_p))  
    do i_press = 1,this%num_p
       axis1%values(i_press) = press_data_array(1,i_press) * &
                               this%press_unit_conv_factor
    end do
      this%lookup_table_gen%dim = 1
      this%lookup_table_gen%dims(1) = this%num_p 
      this%lookup_table_gen%axis1 => axis1  
    nullify(axis1)
    this%n_indices = 1

  else if (this%num_t > 1) then
    !creat general 2D lookip table and load the temperature values in axis1
    this%temperature = UNINITIALIZED_DOUBLE
    allocate(axis1)
    allocate(axis1%values(this%num_t))
    do i_temp = 1,this%num_t
      axis1%values(i_temp) = temp_array(i_temp) * this%temp_unit_conv_factor
    end do
    allocate(this%lookup_table_gen%axis2)
    call LookupTableAxisInit(this%lookup_table_gen%axis2)
    this%lookup_table_gen%axis2%saved_index2 = 1    
    allocate(this%lookup_table_gen%axis2%values(this%num_p))      
    do i_press = 1,this%num_p
      this%lookup_table_gen%axis2%values(i_press) = &
          press_data_array(1,i_press) * this%press_unit_conv_factor
    end do
    this%lookup_table_gen%axis1 => axis1
    this%lookup_table_gen%dim = 2
    this%lookup_table_gen%dims(1) = this%num_t
    this%lookup_table_gen%dims(2) = n_press_count_first      
    this%n_indices = 3
    nullify(axis1) 
  end if

  nullify(var_data)
  allocate(var_data(this%num_prop,n_press_count))
  do i_press = 1,n_press_count
    do i_prop=1,this%num_prop      
      !the conversion is done later - 
      var_data(i_prop,i_press) = press_data_array(i_prop+1,i_press)  
    end do
  end do
  
  this%lookup_table_gen%var_data => var_data
  !convert units
  call this%lookup_table_gen%VarPointAndUnitConv(option)
  !setup constant gradient extrapolation method
  call this%EOSSetConstGradExtrap(option)
  !allocate lookup var gradients 
  call this%lookup_table_gen%LookupTableVarInitGradients(option)

  nullify(var_data)

  call DeallocateArray(press_data_array)
  call DeallocateArray(temp_array) 


end subroutine EOSTableRead

! ************************************************************************** !

subroutine ReadPressureTable(input,option,press_idx,n_press,press_data_array)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Reads the PVT fields vs pressure
  !
  ! "*" entries not yet supported
  !PO TODO: when ecountering "*" as entry in the data field interpolates
  ! using previous and following data field.
  ! read the field data in each line using InputReadDouble,
  ! if getting an error, read as "word" and check if "*"
  ! compare each word vs "*" using function "StringCompare"
  ! Wehn finning "*", map location into an array to poste post-process later
  ! for the interpolation. "*" should not be entered in the first and last
  ! record
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none

  !class(eos_database_type) :: this
  !PetscInt, intent(in) :: num_fields
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscInt, intent(inout) :: press_idx
  PetscInt, intent(out) :: n_press
  PetscReal, pointer :: press_data_array(:,:)

  PetscInt :: i_data
  PetscInt :: size_rank2
  PetscInt :: num_fields !this is size_rank1 - include the pressure field

  num_fields = size(press_data_array,1)
  size_rank2 = size(press_data_array,2)

  n_press = 0
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    if (InputError(input)) then
      option%io_buffer = 'Error PVT table - reading Pressure Table'
      call printErrMsg(option)
    end if
    n_press = n_press + 1
    press_idx = press_idx + 1
    if ( press_idx > size_rank2 ) then
      !each time doubles the size of rank 2
      !tmp_array_size overwritten by new size
      call reallocateRealArray(press_data_array,size_rank2)
    end if
    do i_data = 1, num_fields
      call InputReadDouble(input,option,press_data_array(i_data,press_idx))
      call InputErrorMsg(input,option,&
                         'VALUE','EOS_PVT_TABLE - PRESSURE TABLE VALUE')
    end do
  end do

end subroutine ReadPressureTable

! ************************************************************************** !

!subroutine EOSPropTable(this,T,P,prop_iname,prop_value,ierr,iP1,iT,iP2)
subroutine EOSPropTable(this,T,P,prop_iname,prop_value,indices,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 10/20/17
  !
  ! interpolates a single EOS property from a PVT Table
  !
  ! if table%dim == 1, only pressure lookup
  !  indices(eos_table%first_index) = Pressure index
  ! else if table%dim == 2, one temperature lookup and two pressure lookups
  !  indices(eos_table%first_index+1) = iP1, i.e. Pressure_index_1
  !  indices(eos_table%first_index+2) = iP2, i.e. Pressure_index_2
  ! the case table%dim == 2, with one temperature lookup and
  !                       one pressure lookups not supported
  !
  ! Note: when more properties must be extracted from the same PVT table,
  !       i.e. properties listed for the same values and range of P and T,
  !       all propoertis should be extracted at the same time to perform
  !       only once the look up operations, which are:
  !       (a) LookupTableIndexGeneral, i.e. axis look up
  !       (b) P,T location checks within a single 2D domain, (p,p+dp; t,t+dt)
  !           currently done within LookupTableInterpolate2DGeneral
  !
  !    PO TODO: add a method lookup_table_gen to extract mulitdimensional data
  !             at the same time (data plus function)

  implicit none

  class(eos_table_type) :: this
  PetscReal, intent(in) :: T              ! temperature [C]
  PetscReal, intent(in) :: P              ! pressure [Pa]
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: prop_value ! database units (SI)
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, intent(inout) :: indices(:)

  ierr = 0
  !no errors are assigned if P,T go out of range, since the values
  !are extrapolated using the edge values - however warnning are printed
  !on the screen to signal extrapolation

  this%lookup_table_gen%data => &
                    this%lookup_table_gen%var_array(prop_iname)%ptr%data

  if (this%lookup_table_gen%dim == ONE_INTEGER) then
    if ( P < this%lookup_table_gen%axis1%values(1) ) then
      print*, "EOSPropTable - P smaller than min val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    end if
    if ( P > this%lookup_table_gen%axis1%values(this%num_p) ) then
      print*, "EOSPropTable - P larger than max val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    end if
    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)
    !T       P      optional
    !this%lookup_table_gen%Sample(lookup1,lookup2,lookup3)
    prop_value = this%lookup_table_gen%Sample(P)
    !must copy back iP1 because
    ! saved_index = indices(eos_table%first_index) is a copying operation not pointing
    ! saved_index should be a pointer to save an assigment operation
    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
  else if (this%lookup_table_gen%dim == TWO_INTEGER) then
    if ( T < this%lookup_table_gen%axis1%values(1) ) then
      print*, "EOSPropTable - T smaller than min val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated for Temp val [°C] = ", T
    end if
    if ( T > this%lookup_table_gen%axis1%values(this%num_t) ) then
      print*, "EOSPropTable - T larger than max val in table"
      print*, "table name = ", this%name
      print*, "Property extrapolated Temp val [°C] = ", T
    end if
    ! Due to the possible irregular shape of Pmax(T) and Pmin(T)
    ! a out of bound check requires the computation of Pmin & Pmax in the table
    ! the check below is valid if:
    ! Pmin = axis1%values(1) and Pmax=axis2%values(this%num_p)
    !if ( P < this%lookup_table_gen%axis1%values(1) ) then
    !  print*, "EOSEOSProp - P smaller than min val in EOSdatabase"
    !  print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    !end if
    !if ( P > this%lookup_table_gen%axis2%values(this%num_p) ) then
    !  print*, "EOSEOSProp - P larger than max val in EOSdatabase"
    !  print*, "Property extrapolated for Press val [Mpa] = ", P*1.d-6
    !end if
    !this%lookup_table_gen%axis1%saved_index = iT
    !this%lookup_table_gen%axis2%saved_index = iP1
    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)
    this%lookup_table_gen%axis2%saved_index = indices(this%first_index+1)
    this%lookup_table_gen%axis2%saved_index2 = indices(this%first_index+2)
    !if (present(iP2)) then
    !  this%lookup_table_gen%axis2%saved_index2 = iP2
    !else
    !  this%lookup_table_gen%axis2%saved_index2 = iP1
    !end if
    !T       P      optional
    !this%lookup_table_gen%Sample(lookup1,lookup2,lookup3)
    prop_value = this%lookup_table_gen%Sample(T,P)
    !must copy back iT, iP1 and iP2 because
    ! saved_index = iT is a copying operation not a pointer
    ! saved_index should be a pointer to save an assigment operation
    !iT = this%lookup_table_gen%axis1%saved_index
    !iP1 = this%lookup_table_gen%axis2%saved_index
    !if (present(iP2)) then
    !  iP2 = this%lookup_table_gen%axis2%saved_index2
    !else
    !  iP2 = iP1
    !end if
    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
    indices(this%first_index+1)= this%lookup_table_gen%axis2%saved_index
    indices(this%first_index+2) = this%lookup_table_gen%axis2%saved_index2
  end if

  nullify(this%lookup_table_gen%data)

end subroutine EOSPropTable

! ************************************************************************** !

subroutine EOSPropGradTable(this,T,P,prop_iname,sample,dSampledT,dSampledP, &
                            indices,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 10/20/17
  !
  ! interpolates a single EOS property from a PVT Table and computes its 
  ! gradient with respect to T and P as given in the PVT table
  !
  ! if table%dim == 1, only pressure lookup
  !  indices(eos_table%first_index) = Pressure index
  ! else if table%dim == 2, one temperature lookup and two pressure lookups
  !  indices(eos_table%first_index+1) = iP1, i.e. Pressure_index_1
  !  indices(eos_table%first_index+2) = iP2, i.e. Pressure_index_2
  ! the case table%dim == 2, with one temperature lookup and
  !                       one pressure lookups not supported
  !

  implicit none

  class(eos_table_type) :: this
  PetscReal, intent(in) :: T              ! temperature [C]
  PetscReal, intent(in) :: P              ! pressure [Pa]
  PetscInt, intent(in) :: prop_iname
  PetscReal, intent(out) :: sample ! PFLOTRAN intenral units (SI)
  PetscReal, intent(out) :: dSampledT ! internal PFLOTRAN units ([Sample]/[C])
  PetscReal, intent(out) :: dSampledP ! internal PFLOTRAN units ([Sample]/[Pa])
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, intent(inout) :: indices(:)

  ierr = 0

  !bound checks to be included later as part of reporting
  !go trgough the table list and save the outbound in a report data strucutre?
  !could be done in auxvar compute using a reusable table list specific func

  if (this%lookup_table_gen%dim == ONE_INTEGER) then

    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)

    call this%lookup_table_gen%SampleAndGradient(prop_iname,P)
    sample = this%lookup_table_gen%var_array(prop_iname)%ptr%sample
    dSampledP = this%lookup_table_gen%var_array(prop_iname)%ptr%sample_grad(1)
    dSampledT = 0.0d0
    !must copy back the index because
    ! saved_index = indices(eos_table%first_index) is a copying operation not pointing
    ! saved_index should be a pointer to save two assigment operations (in/out)
    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
  else if (this%lookup_table_gen%dim == TWO_INTEGER) then

    this%lookup_table_gen%axis1%saved_index = indices(this%first_index)
    this%lookup_table_gen%axis2%saved_index = indices(this%first_index+1)
    this%lookup_table_gen%axis2%saved_index2 = indices(this%first_index+2)
    
    call this%lookup_table_gen%SampleAndGradient(prop_iname,T,P)
    sample = this%lookup_table_gen%var_array(prop_iname)%ptr%sample
    dSampledT = this%lookup_table_gen%var_array(prop_iname)%ptr%sample_grad(1)
    dSampledP = this%lookup_table_gen%var_array(prop_iname)%ptr%sample_grad(2)

    indices(this%first_index) = this%lookup_table_gen%axis1%saved_index
    indices(this%first_index+1)= this%lookup_table_gen%axis2%saved_index
    indices(this%first_index+2) = this%lookup_table_gen%axis2%saved_index2
  end if

end subroutine EOSPropGradTable

! ************************************************************************** !

subroutine EOSTableDestroy(eos_table)
  !
  ! Author: Paolo Orsini
  ! Date: 10/21/17
  !
  ! destroys EOS Table

  use Utility_module

  implicit none

  class(eos_table_type), pointer :: eos_table

  if (.not.associated(eos_table)) return

  call eos_table%EOSDataBaseStrip()
  call LookupTableDestroy(eos_table%lookup_table_gen)

  deallocate(eos_table)
  nullify(eos_table)

end subroutine EOSTableDestroy

! ************************************************************************** !

subroutine EOSTableInitList()
  !
  ! Initializes eos_table_list
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !

  implicit none

  !type(eos_table_list_type), pointer :: list

  allocate(eos_table_list)

  nullify(eos_table_list%first)
  nullify(eos_table_list%last)
  nullify(eos_table_list%array)
  eos_table_list%num_eos_tables = 0

end subroutine EOSTableInitList

! ************************************************************************** !

subroutine EOSTableAddToList(new_eos_table,list)
  !
  ! Adds a new eos_table to an eos table list
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !
  implicit none

  class(eos_table_type), pointer :: new_eos_table
  type(eos_table_list_type) :: list

  list%num_eos_tables = list%num_eos_tables + 1
  new_eos_table%id = list%num_eos_tables
  if (.not.associated(list%first)) list%first => new_eos_table
  if (associated(list%last)) list%last%next => new_eos_table
  list%last => new_eos_table

end subroutine EOSTableAddToList

! ************************************************************************** !

subroutine EOSTableProcessList(option)
  !
  ! Loop through EOS table list and build indices
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17

  use Option_module

  implicit none

  type(option_type) :: option

  type(eos_table_list_type), pointer :: list
  class(eos_table_type), pointer :: eos_table

  list => eos_table_list

  if (.not.associated(list)) return

  option%neos_table_indices = 0
  !loop over  EOS tables and create a map to store save indices in auxvars
  eos_table => list%first
  do
    if (.not.(associated(eos_table))) exit
    eos_table%first_index = option%neos_table_indices + 1
    option%neos_table_indices = option%neos_table_indices + eos_table%n_indices
    ! add here other operation that must be performed on all eos tables
    ! move to next table if it is associated
    if( .not.(associated(eos_table%next))) exit
    eos_table => eos_table%next
  enddo

end subroutine EOSTableProcessList

! ************************************************************************** !

subroutine EOSTableDestroyList()
  !
  ! Destroy a list of EOS tables
  !
  ! Author: Paolo Orsini
  ! Date: 10/24/17
  !

  implicit none

  type(eos_table_list_type), pointer :: list

  class(eos_table_type), pointer :: eos_table, prev_eos_table

  list => eos_table_list
  
  if (.not.associated(list)) return

  eos_table => list%first
  do
    if (.not.associated(eos_table)) exit
    prev_eos_table => eos_table
    eos_table => eos_table%next
    call EOSTableDestroy(prev_eos_table)
    nullify(prev_eos_table)
  enddo

  list%num_eos_tables = 0
  nullify(list%first)
  nullify(list%last)
  if (associated(list%array)) deallocate(list%array)
  nullify(list%array)

  deallocate(list)
  nullify(list)

  nullify(eos_table_list)

end subroutine EOSTableDestroyList

! ************************************************************************** !

end module EOSData_module
