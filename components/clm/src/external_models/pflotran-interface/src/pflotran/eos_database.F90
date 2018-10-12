module EOSDatabase_module
 
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none
  
  private
 
#include "petsc/finclude/petscsys.h"
 
  PetscInt, parameter, public :: EOS_DENSITY = 1
  PetscInt, parameter, public :: EOS_ENTHALPY = 2
  PetscInt, parameter, public :: EOS_VISCOSITY = 3 
  !add here other properties and derivatives and increase MAX_PROP_NUM
  PetscInt, parameter, public :: MAX_PROP_NUM = 3

  type, public :: eos_database_type 
    character(len=MAXWORDLENGTH) :: dbase_name
    character(len=MAXWORDLENGTH) :: file_name
    PetscInt :: num_dp   ! number of pressure intervals
    PetscInt :: num_dt   ! number of temperature intervals
    PetscInt :: num_prop ! number of properties in the database
    PetscReal :: dp      ! uniform pressure interval
    PetscReal :: dt      ! uniform temperature interval 
    PetscInt :: data_to_prop_map(MAX_PROP_NUM) !map the data idx to the prop. 
    PetscInt :: prop_to_data_map(MAX_PROP_NUM) !map the prop to the dat idx
    PetscReal, pointer :: data(:,:)
    class(lookup_table_uniform_type), pointer :: lookup_table
  contains
    procedure, public :: Read => EOSDatabaseRead 
    procedure, public :: EOSProp => EOSPropLinearInterp
    procedure, public :: EOSPropPresent
  end type

  public :: EOSDatabaseCreate, &
            EOSDatabaseDestroy

contains

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

  allocate(EOSDatabaseCreate)
  EOSDatabaseCreate%dbase_name = dbase_name 
  EOSDatabaseCreate%file_name = filename
  EOSDatabaseCreate%num_dp = UNINITIALIZED_INTEGER
  EOSDatabaseCreate%num_dt = UNINITIALIZED_INTEGER
  EOSDatabaseCreate%dp = UNINITIALIZED_DOUBLE
  EOSDatabaseCreate%dt = UNINITIALIZED_DOUBLE 
  EOSDatabaseCreate%data_to_prop_map(1:MAX_PROP_NUM) = UNINITIALIZED_INTEGER
  EOSDatabaseCreate%prop_to_data_map(1:MAX_PROP_NUM) = UNINITIALIZED_INTEGER

  nullify(EOSDatabaseCreate%lookup_table)
  nullify(EOSDatabaseCreate%data)

end function EOSDatabaseCreate

! ************************************************************************** !

subroutine EOSDatabaseRead(this,option)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  ! 
  ! Reads the the an EOS database from a text file for one or more 
  ! phase properties.  
  ! 
  ! Database format (for look up table)
  ! Header: NUM_DP, NUM_DT and DATA_LIST_ORDER
  ! DATA
  ! column 1 = pressure
  ! column 2 = temperature
  ! column 2+1, column 2+n: properties in the order given in DATA_LIST_ORDER
  !
  ! The dataset must be followed by NUM_DP * NUM_DT lines, not interrupted
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

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string = 'EOS_DATABASE'
  character(len=MAXSTRINGLENGTH) :: string 
  PetscInt :: prop_idx, prop_count, i_idx, j_idx 
  PetscInt :: data_size
  PetscReal :: tempreal
  PetscBool :: pres_present, temp_present

  type(input_type), pointer :: input_table


  if (len_trim(this%file_name) < 1) then
    option%io_buffer = 'FILENAME must be specified for EOS_DATABASE.'
    call printErrMsg(option)
  endif

  input_table => InputCreate(IUNIT_TEMP,this%file_name,option)
  input_table%ierr = 0

  !if ( option%myrank == 0 ) then
    option%io_buffer = 'Reading database = ' // this%file_name
    call printMsg(option) 
  !end if

  !reading the database file header
  do

    call InputReadPflotranString(input_table,option)
    !if (InputCheckExit(input_table,option)) exit
    if (InputError(input_table)) exit

    call InputReadWord(input_table,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input_table,option,'keyword',error_string)
    call StringToUpper(keyword)   
    select case(keyword)
      case('NUM_DP')
        call InputReadInt(input_table,option,this%num_dp)
        call InputErrorMsg(input_table,option,'number of dp',error_string)
      case('NUM_DT')
        call InputReadInt(input_table,option,this%num_dt)
        call InputErrorMsg(input_table,option,'number of dt',error_string)        
      case('DATA_LIST_ORDER')  
        pres_present = PETSC_FALSE; temp_present = PETSC_FALSE;
        prop_idx = 0
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
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_DENSITY 
              this%prop_to_data_map(EOS_DENSITY) = prop_idx 
            case('ENTHALPY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_ENTHALPY 
              this%prop_to_data_map(EOS_ENTHALPY) = prop_idx 
            case('VISCOSITY')
              prop_idx = prop_idx + 1
              this%data_to_prop_map(prop_idx) = EOS_VISCOSITY 
              this%prop_to_data_map(EOS_VISCOSITY) = prop_idx 
            case default
              error_string = trim(error_string) // ': ' // this%file_name // &
              ': DATA_LIST_ORDER'
              call InputKeywordUnrecognized(keyword,error_string,option)
          end select
        end do
        this%num_prop = prop_idx
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

  data_size = this%num_dp * this%num_dt

  this%num_prop = prop_idx
  allocate(this%data(prop_idx,data_size))

  !create lookup table  
  this%lookup_table => LookupTableCreateUniform(TWO_INTEGER)
  this%lookup_table%dims(1) = this%num_dt
  this%lookup_table%dims(2) = this%num_dp
  allocate(this%lookup_table%axis1%values(this%num_dt))
  this%lookup_table%axis1%values(1:this%num_dt) = UNINITIALIZED_DOUBLE
  allocate(this%lookup_table%axis2%values(this%num_dp))
  this%lookup_table%axis2%values(1:this%num_dp) = UNINITIALIZED_DOUBLE

  !TODO
  ! start loading data - at the moment using Input facility, however this file
  ! can be large. TODO. Implement more efficient solutions:
  ! - using read(,) without InputReadDouble filter
  ! - adding the option of reading a .h5 file where the database is defined  

  ! go to data - first time to load axis1 and 2 values
  string = "DATA"
  call InputFindStringInFile(input_table,option,string)

  do j_idx = 1,this%num_dp 
       
    do i_idx = 1, this%num_dt
      call InputReadPflotranString(input_table,option)

      call InputReadDouble(input_table,option,tempreal) 
      call InputErrorMsg(input_table,option, &
                           'VALUE', 'EOS_DATABASE PRESS_VALUE') 
      ! convert MPa in Pa
      this%lookup_table%axis2%values(j_idx) = tempreal * 1.0d6
   
      ! this is repeated this%num_dp times - not efficient
      call InputReadDouble(input_table,option, &
                           this%lookup_table%axis1%values(i_idx))
      call InputErrorMsg(input_table,option, &
                         'VALUE', 'EOS_DATABASE TEMP_VALUE') 
    
      prop_count = i_idx + (j_idx-1) * this%num_dt
      do prop_idx = 1,this%num_prop
        call InputReadDouble(input_table,option,this%data(prop_idx,prop_count))
        call InputErrorMsg(input_table,option,&
                           'VALUE','EOS_DATABASE PROP_VALUE')
      end do      

    end do

  end do

  call InputDestroy(input_table)

end subroutine EOSDatabaseRead

! ************************************************************************** !
function EOSPropPresent(this,prop_iname)
  ! 
  ! Author: Paolo Orsini
  ! Date: 12/18/15
  ! 
  ! Checks if a property is defined in the database

  implicit none

  class(eos_database_type) :: this
  PetscInt, intent(in) :: prop_iname
  PetscBool :: EOSPropPresent
                                        
  EOSPropPresent = Initialized(this%prop_to_data_map(prop_iname))

end function EOSPropPresent

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
  !       TODO: add a method lookup_table to extract mulitdimensional data 
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
  if ( T < this%lookup_table%axis1%values(1) ) then 
    ierr = 101
    print*, "EOSEOSProp - T smaller than min val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop  
  end if 
  if ( T > this%lookup_table%axis1%values(this%num_dt) ) then
    ierr = 102
    print*, "EOSEOSProp - T larger than max val in EOSdatabase"
    print*, "Temp val [°C] = ", T
    stop  
  end if
  if ( P < this%lookup_table%axis2%values(1) ) then
    ierr = 103
    print*, "EOSEOSProp - P smaller than min val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop  
  end if
  if ( P > this%lookup_table%axis2%values(this%num_dp) ) then
    ierr = 104
    print*, "EOSEOSProp - P larger than max val in EOSdatabase"
    print*, "Press val [Mpa] = ", P*1.d-6
    stop  
  end if

  this%lookup_table%data => this%data(this%prop_to_data_map(prop_iname),:)

                             !T       P      optional
  !this%lookup_table%Sample(lookup1,lookup2,lookup3) 
  prop_value = this%lookup_table%Sample(T,P) 

  nullify(this%lookup_table%data)

end subroutine EOSPropLinearInterp

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

  !deallocate(eos_database%data)
  !nullify(eos_database%data) 
  call DeallocateArray(eos_database%data)

  call LookupTableDestroy(eos_database%lookup_table) 

  deallocate(eos_database)
  nullify(eos_database)

end subroutine EOSDatabaseDestroy

! ************************************************************************** !

end module EOSDatabase_module
