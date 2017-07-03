module Fracture_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: frac_init_pres_index = 1
  PetscInt, parameter, public :: frac_alt_pres_index = 2
  PetscInt, parameter, public :: frac_max_poro_index = 3
  PetscInt, parameter, public :: frac_poro_exp_index = 4
  PetscInt, parameter, public :: frac_change_perm_x_index = 1
  PetscInt, parameter, public :: frac_change_perm_y_index = 2
  PetscInt, parameter, public :: frac_change_perm_z_index = 3
  
  type, public :: fracture_type
    PetscReal :: init_pressure
    PetscReal :: altered_pressure
    PetscReal :: maximum_porosity
    PetscReal :: porosity_exponent
    PetscReal :: change_perm_x
    PetscReal :: change_perm_y
    PetscReal :: change_perm_z
  contains
    procedure, public :: Read => FractureRead
  end type fracture_type
  
!  class(fracture_type), pointer, public :: fracture

  public :: FractureInit, &
            FractureCreate, &
            FractureSetInitialPressure, &
            FractureAuxVarInit, &
            FracturePropertytoAux, &
            FractureDestroy, &
            FracturePoroEvaluate, &
            FracturePermEvaluate
  
  contains

! ************************************************************************** !

function FractureCreate()
  !
  ! Author: Heeho Park
  ! Date: 4/7/15
  !

  implicit none
  
  class(fracture_type), pointer :: FractureCreate
  class(fracture_type), pointer :: fracture
  
  allocate(fracture)
  call FractureInit(fracture)
  
  FractureCreate => fracture
  
end function FractureCreate

! ************************************************************************** !

subroutine FractureInit(this)
  !
  ! Author: Heeho Park
  ! Date: 4/7/2015
  !

  implicit none
  
  class(fracture_type), pointer :: this
  
  this%init_pressure = UNINITIALIZED_DOUBLE
  this%altered_pressure = UNINITIALIZED_DOUBLE
  this%maximum_porosity = UNINITIALIZED_DOUBLE
  this%porosity_exponent = UNINITIALIZED_DOUBLE
  this%change_perm_x = 0.d0
  this%change_perm_y = 0.d0
  this%change_perm_z = 0.d0

end subroutine FractureInit

! ************************************************************************** !

subroutine FractureAuxVarInit(fracture_material,auxvar)
  !
  ! Author: Heeho Park
  ! Date: 7/8/2015
  !

  use Material_Aux_class
  
  implicit none
  
  class(fracture_type), pointer :: fracture_material
  class(material_auxvar_type), intent(inout) :: auxvar

  if (associated(fracture_material)) then
    allocate(auxvar%fracture)
    allocate(auxvar%fracture%properties(4))
    allocate(auxvar%fracture%vector(3))
    auxvar%fracture%properties = 0.d0
    auxvar%fracture%vector = 0.d0
  endif

end subroutine FractureAuxVarInit

! ************************************************************************** !

subroutine FracturePropertytoAux(auxvar,fracture_property)
  !
  ! Author: Heeho Park
  ! Date: 7/8/2015
  !

  use Material_Aux_class
  
  implicit none

  class(material_auxvar_type), intent(inout) :: auxvar
  class(fracture_type), pointer :: fracture_property

  
  auxvar%fracture%properties(frac_init_pres_index) = &
    fracture_property%init_pressure
  auxvar%fracture%properties(frac_alt_pres_index) = &
    fracture_property%altered_pressure
  auxvar%fracture%properties(frac_max_poro_index) = &
    fracture_property%maximum_porosity
  auxvar%fracture%properties(frac_poro_exp_index) = &
    fracture_property%porosity_exponent
  auxvar%fracture%vector(frac_change_perm_x_index) = &
    fracture_property%change_perm_x
  auxvar%fracture%vector(frac_change_perm_y_index) = &
    fracture_property%change_perm_y
  auxvar%fracture%vector(frac_change_perm_z_index) = &
    fracture_property%change_perm_z

end subroutine FracturePropertytoAux

! ************************************************************************** !

subroutine FractureRead(this,input,option)
  ! 
  ! Author: Heeho Park
  ! Date: 4/7/15
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  
  implicit none
  
  class(fracture_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: word
  
  do
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option, &
                                    'MATERIAL_PROPERTY,WIPP-FRACTURE')
          
      if (InputCheckExit(input,option)) exit
          
      if (InputError(input)) exit
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,'keyword', &
                          'MATERIAL_PROPERTY,WIPP-FRACTURE')   
      select case(trim(word))
        case('INITIATING_PRESSURE')
          call InputReadDouble(input,option, &
                                this%init_pressure)
          call InputErrorMsg(input,option, &
                              'initiating pressure of fracturing', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('ALTERED_PRESSURE')
          call InputReadDouble(input,option, &
                                this%altered_pressure)
          call InputErrorMsg(input,option, &
                              'altered pressure of fracturing', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('MAXIMUM_FRACTURE_POROSITY')
          call InputReadDouble(input,option, &
                                this%maximum_porosity)
          call InputErrorMsg(input,option, &
                              'maximum fracture porosity', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('FRACTURE_EXPONENT')
          call InputReadDouble(input,option, &
                              this%porosity_exponent)
          call InputErrorMsg(input,option, &
                          'dimensionless fracture exponent for porosity', &
                              'MATERIAL_PROPERTY,WIPP-FRACTURE')
        case('ALTER_PERM_X')
          this%change_perm_x = 1.d0
        case('ALTER_PERM_Y')
          this%change_perm_y = 1.d0
        case('ALTER_PERM_Z')
          this%change_perm_z = 1.d0
        case default
          call InputKeywordUnrecognized(word, &
                  'MATERIAL_PROPERTY,WIPP-FRACTURE',option)
      end select
    enddo

end subroutine FractureRead

! ************************************************************************** !

subroutine FractureSetInitialPressure(fracture,initial_cell_pressure)
  !
  ! Sets the pressure referenced in fracture
  !
  use Material_Aux_class

  implicit none
  
  type(fracture_auxvar_type) :: fracture
  PetscReal, intent(in) :: initial_cell_pressure
  
  fracture%properties(frac_init_pres_index) = &
    fracture%properties(frac_init_pres_index) + initial_cell_pressure
  fracture%properties(frac_alt_pres_index) = &
    fracture%properties(frac_alt_pres_index) + &
    fracture%properties(frac_init_pres_index)

end subroutine FractureSetInitialPressure

! ************************************************************************** !

subroutine FracturePoroEvaluate(auxvar,pressure,compressed_porosity, &
                                dcompressed_porosity_dp)
  !
  ! Calculates porosity induced by fracture BRAGFLO_6.02_UM Eq. (136)
  ! 4.10 Pressure-Induced Fracture Treatment
  !
  ! Author: Heeho Park
  ! Date: 03/12/15
  !

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  
!  class(fracture_type) :: this
  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: pressure
  PetscReal, intent(out) :: compressed_porosity
  PetscReal, intent(out) :: dcompressed_porosity_dp
  
  PetscReal :: Ci, Ca
  PetscReal :: P0, Pa, Pi
  PetscReal :: phia, phi0

       ! convert bulk compressibility to pore compressibility
  Ci = auxvar%soil_properties(soil_compressibility_index) / &
       auxvar%porosity_base
  P0 = auxvar%soil_properties(soil_reference_pressure_index)
  Pa = auxvar%fracture%properties(frac_alt_pres_index)
  Pi = auxvar%fracture%properties(frac_init_pres_index)
  phia = auxvar%fracture%properties(frac_max_poro_index)
  phi0 = auxvar%porosity_base
  
  if (.not.associated(MaterialCompressSoilPtr, &
                      MaterialCompressSoilBRAGFLO)) then
    option%io_buffer = 'WIPP Fracture Function must be used with ' // &
      'BRAGFLO soil compressibility function.'
    call printErrMsg(option)
  endif
  
  if (pressure < Pi) then
    call MaterialCompressSoil(auxvar,pressure, compressed_porosity, &
                              dcompressed_porosity_dp)
  else if (pressure > Pi .and. pressure < Pa) then
    Ca = Ci*(1.d0 - 2.d0 * (Pa-P0)/(Pa-Pi)) + &
      2.d0/(Pa-Pi)*log(phia/phi0)
    compressed_porosity = phi0 * exp(Ci*(pressure-P0) + &
      ((Ca-Ci)*(pressure-Pi)**2.d0)/(2.d0*(Pa-Pi)))
    !mathematica solution
    dcompressed_porosity_dp = exp(Ci*(pressure-P0) + &
      ((Ca-Ci)*(pressure-Pi)**2.d0) / (2.d0*(Pa-Pi))) * &
      phi0 * (Ci + ((Ca-Ci)*(pressure-Pi)) / (Pa-Pi))
  else if (pressure >= Pa) then
    compressed_porosity = phia
    dcompressed_porosity_dp = 0.d0
  endif

end subroutine FracturePoroEvaluate

! ************************************************************************** !
                                
subroutine FracturePermEvaluate(auxvar,permeability,altered_perm, &
                                    daltered_perm_dp,dist)
  !
  ! Calculates permeability induced by fracture BRAGFLO_6.02_UM Eq. (136)
  ! 4.10 Pressure-Induced Fracture Treatment
  !
  ! Author: Heeho Park
  ! Date: 03/12/15
  !

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(option_type) :: option
  
!  class(fracture_type) :: this
  class(material_auxvar_type), intent(in) :: auxvar
  PetscReal, intent(in) :: permeability
  PetscReal, intent(out) :: altered_perm
  PetscReal, intent(out) :: daltered_perm_dp
  PetscReal :: dist(-1:3)

  PetscReal :: phii, n
  PetscReal :: phi

  if (dot_product(dist(1:3),auxvar%fracture%vector) < 1.d-40) return
  
  phi = auxvar%porosity
  phii = auxvar%porosity_base
  n = auxvar%fracture%properties(frac_poro_exp_index)

  if (.not.associated(MaterialCompressSoilPtr, &
                      MaterialCompressSoilBRAGFLO)) then
    option%io_buffer = 'WIPP Fracture Function must be used with ' // &
      'BRAGFLO soil compressibility function.'
    call printErrMsg(option)
  endif
  
  ! phi = altered porosity
  ! phii = porosity at initiating pressure
  altered_perm = permeability * (phi/phii)**n
  !the derivative ignored at this time since it requires additional parameters
  !and calculations. we'll will need cell_pressure as an input and other 
  !parameters used in MaterialFracturePorosityWIPP
  daltered_perm_dp = 1.d-10
  ! Mathematica solution
  !Ca = Ci*(1.d0 - 2.d0 * (Pa-P0)/(Pa-Pi)) + &
  !   2.d0/(Pa-Pi)*log(phia/phi0)
  !a = exp(Ci*(pressure-P0) + &
  !    ((Ca-Ci)*(pressure-Pi)**2.d0) / (2.d0*(Pa-Pi))) * &
  !    phi0 * (Ci + ((Ca-Ci)*(pressure-Pi)) / (Pa-Pi))
  !b = exp(Ci*(pressure-P0) + &
  !    ((Ca-Ci)*(pressure-Pi)**2.d0) / (2.d0*(Pa-Pi)))
  !daltered_perm_dp = a * permeability * n * (b*phi0/phii)**(n-1)
  
end subroutine FracturePermEvaluate

! ************************************************************************** !

subroutine FractureDestroy(this)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(fracture_type), pointer :: this
  
  if (.not.associated(this)) return

  deallocate(this)
  nullify(this)

end subroutine FractureDestroy

end module Fracture_module

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
  
module Creep_Closure_module
  
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type, public :: creep_closure_type
    character(len=MAXWORDLENGTH) :: material_name
    PetscInt :: imat
    PetscInt :: num_times
    PetscInt :: num_values_per_time
    class(lookup_table_general_type), pointer :: lookup_table
  contains
    procedure, public :: Read => CreepClosureRead
    procedure, public :: Evaluate => CreepClosureEvaluate
    procedure, public :: Test => CreepClosureTest
  end type creep_closure_type
  
  class(creep_closure_type), pointer, public :: creep_closure
  
  interface CreepClosureDestroy
    module procedure CreepClosureDestroy1
    module procedure CreepClosureDestroy2
  end interface

  public :: CreepClosureInit, &
            CreepClosureCreate, &
            CreepClosureDestroy
  
  contains
  
! ************************************************************************** !

subroutine CreepClosureInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  if (associated(creep_closure)) then
    call CreepClosureDestroy(creep_closure)
  endif
  nullify(creep_closure)
  
end subroutine CreepClosureInit

! ************************************************************************** !

function CreepClosureCreate()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: CreepClosureCreate
  
  allocate(CreepClosureCreate)
  CreepClosureCreate%material_name = ''
  CreepClosureCreate%imat = UNINITIALIZED_INTEGER
  CreepClosureCreate%num_times = UNINITIALIZED_INTEGER
  CreepClosureCreate%num_values_per_time = UNINITIALIZED_INTEGER
  nullify(CreepClosureCreate%lookup_table)
  
end function CreepClosureCreate

! ************************************************************************** !

subroutine CreepClosureRead(this,input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  
  implicit none
  
  class(creep_closure_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: keyword, word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_string = 'CREEP_CLOSURE'
  type(input_type), pointer :: input2
  PetscInt :: temp_int
  PetscReal :: time_units_conversion

  time_units_conversion = 1.d0
  filename = ''
  input%ierr = 0

  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('FILENAME') 
        call InputReadWord(input,option,filename,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename',error_string)
      case('MATERIAL') 
        call InputReadWord(input,option,this%material_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'material',error_string)
     case default
        call InputKeywordUnrecognized(keyword,'CREEP_CLOSURE',option)
    end select
  enddo
  
  if (len_trim(filename) < 1) then
    option%io_buffer = 'FILENAME must be specified for CREEP_CLOSURE.'
    call printErrMsg(option)
  endif
  
  this%lookup_table => LookupTableCreateGeneral(TWO_INTEGER)
  error_string = 'CREEP_CLOSURE file'
  input2 => InputCreate(IUNIT_TEMP,filename,option)
  input2%ierr = 0
  do
  
    call InputReadPflotranString(input2,option)

    if (InputError(input2)) exit

    call InputReadWord(input2,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input2,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('NUM_TIMES') 
        call InputReadInt(input2,option,this%num_times)
        call InputErrorMsg(input2,option,'number of times',error_string)
      case('NUM_VALUES_PER_TIME') 
        call InputReadInt(input2,option,this%num_values_per_time)
        call InputErrorMsg(input2,option,'number of pressure',error_string)
      case('TIME_UNITS') 
        internal_units = 'sec'
        call InputReadWord(input2,option,word,PETSC_TRUE) 
        call InputErrorMsg(input2,option,'UNITS','CONDITION')   
        call StringToLower(word)
        time_units_conversion = UnitsConvertToInternal(word, &
                                internal_units,option)
      case('TIME')
        if (Uninitialized(this%num_times) .or. &
            Uninitialized(this%num_values_per_time)) then
          option%io_buffer = 'NUM_TIMES and NUM_VALUES_PER_TIME must be ' // &
            'specified prior to reading the corresponding arrays.'
          call printErrMsg(option)
        endif
        this%lookup_table%dims(1) = this%num_times
        this%lookup_table%dims(2) = this%num_values_per_time
        temp_int = this%num_times*this%num_values_per_time
        allocate(this%lookup_table%axis1%values(this%num_times))
        allocate(this%lookup_table%axis2%values(temp_int))
        allocate(this%lookup_table%data(temp_int))
        string = 'TIME in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%axis1%values, &
                              NEG_ONE_INTEGER,string, &
                              input2,option)
        this%lookup_table%axis1%values = this%lookup_table%axis1%values * &
          time_units_conversion
      case('PRESSURE') 
        string = 'PRESSURE in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%axis2%values, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
      case('POROSITY') 
        string = 'POROSITY in CREEP_CLOSURE'
        call UtilityReadArray(this%lookup_table%data, &
                              NEG_ONE_INTEGER, &
                              string,input2,option)
     case default
        error_string = trim(error_string) // ': ' // filename
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  call InputDestroy(input2)
  
  if (size(this%lookup_table%axis1%values) /= this%num_times) then
    option%io_buffer = 'Number of times does not match NUM_TIMES.'
    call printErrMsg(option)
  endif  
  if (size(this%lookup_table%axis2%values) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of pressures does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif
  if (size(this%lookup_table%data) /= &
    this%num_times*this%num_values_per_time) then
    option%io_buffer = 'Number of porosities does not match NUM_TIMES * ' // &
                       'NUM_VALUES_PER_TIME.'
    call printErrMsg(option)
  endif
  
end subroutine CreepClosureRead

! ************************************************************************** !

function CreepClosureEvaluate(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  PetscReal :: CreepClosureEvaluate
  
  CreepClosureEvaluate = this%lookup_table%Sample(time,pressure)
  
end function CreepClosureEvaluate


! ************************************************************************** !

subroutine CreepClosureTest(this,time,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(creep_closure_type) :: this
  PetscReal :: time
  PetscReal :: pressure
  
  print *, time, pressure, this%Evaluate(time,pressure)
  
end subroutine CreepClosureTest

! ************************************************************************** !

subroutine CreepClosureDestroy1()
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  call CreepClosureDestroy(creep_closure)

end subroutine CreepClosureDestroy1

! ************************************************************************** !

subroutine CreepClosureDestroy2(creep_closure)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(creep_closure_type), pointer :: creep_closure
  
  if (.not.associated(creep_closure)) return

  call LookupTableDestroy(creep_closure%lookup_table)
  deallocate(creep_closure)
  nullify(creep_closure)

end subroutine CreepClosureDestroy2

end module Creep_Closure_module

! ************************************************************************** !
! ************************************************************************** !
! ************************************************************************** !
  
module Klinkenberg_module
  
  use PFLOTRAN_Constants_module
  use Lookup_Table_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type, public :: klinkenberg_type
    PetscReal :: a
    PetscReal :: b
  contains
    procedure, public :: Read => KlinkenbergRead
    procedure, public :: Evaluate => KlinkenbergEvaluate
    procedure, public :: Test => KlinkenbergTest
  end type klinkenberg_type
  
  class(klinkenberg_type), pointer, public :: klinkenberg
  
  interface KlinkenbergDestroy
    module procedure KlinkenbergDestroy1
    module procedure KlinkenbergDestroy2
  end interface

  public :: KlinkenbergInit, &
            KlinkenbergCreate, &
            KlinkenbergDestroy
  
  contains
  
! ************************************************************************** !

subroutine KlinkenbergInit()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  if (associated(klinkenberg)) then
    call KlinkenbergDestroy(klinkenberg)
  endif
  nullify(klinkenberg)
  
end subroutine KlinkenbergInit

! ************************************************************************** !

function KlinkenbergCreate()
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(klinkenberg_type), pointer :: KlinkenbergCreate
  
  allocate(KlinkenbergCreate)
  KlinkenbergCreate%a = UNINITIALIZED_DOUBLE
  KlinkenbergCreate%b = UNINITIALIZED_DOUBLE
  
end function KlinkenbergCreate

! ************************************************************************** !

subroutine KlinkenbergRead(this,input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Units_module
  
  implicit none
  
  class(klinkenberg_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string = 'KLINKENBERG_EFFECT'

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('A') 
        call InputReadDouble(input,option,this%a)
        call InputErrorMsg(input,option,'a',error_string)
      case('B') 
        call InputReadDouble(input,option,this%b)
        call InputErrorMsg(input,option,'b',error_string)
     case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
  if (Uninitialized(this%a)) then
    option%io_buffer = &
      'Parameter "a" must be included to model the Klinkenberg Effect.'
    call printErrMsg(option)
  endif
  if (Uninitialized(this%b)) then
    option%io_buffer = &
      'Parameter "b" must be included to model the Klinkenberg Effect.'
    call printErrMsg(option)
  endif
  
end subroutine KlinkenbergRead

! ************************************************************************** !

function KlinkenbergEvaluate(this,liquid_permeability,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(klinkenberg_type) :: this
  PetscReal :: liquid_permeability
  PetscReal :: pressure
  
  PetscReal :: KlinkenbergEvaluate
  PetscReal :: gas_permeability
  
  gas_permeability = liquid_permeability * &
                        (1.d0 + &
                         this%b * (liquid_permeability**this%a) / &
                           pressure)
  KlinkenbergEvaluate = gas_permeability
  
end function KlinkenbergEvaluate

! ************************************************************************** !

subroutine KlinkenbergTest(this,liquid_permeability,pressure)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  implicit none
  
  class(klinkenberg_type) :: this
  PetscReal :: liquid_permeability
  PetscReal :: pressure
  
  print *, liquid_permeability, pressure, &
           this%Evaluate(liquid_permeability,pressure)
  
end subroutine Klinkenbergtest

! ************************************************************************** !

subroutine KlinkenbergDestroy1()
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  call KlinkenbergDestroy(klinkenberg)

end subroutine KlinkenbergDestroy1

! ************************************************************************** !

subroutine KlinkenbergDestroy2(klinkenberg)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 

  implicit none
  
  class(klinkenberg_type), pointer :: klinkenberg
  
  if (.not.associated(klinkenberg)) return

  deallocate(klinkenberg)
  nullify(klinkenberg)

end subroutine KlinkenbergDestroy2

end module Klinkenberg_module

! ************************************************************************** !

module WIPP_module
  
  use PFLOTRAN_Constants_module
  use Creep_Closure_module

  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  type :: wipp_type
    class(creep_closure_type), pointer :: creep_closure
  end type wipp_type
  
  type(wipp_type), pointer, public :: wipp
  
  interface WIPPDestroy
    module procedure WIPPDestroy1
    module procedure WIPPDestroy2
  end interface
  
  public :: WIPPInit, &
            WIPPGetPtr, &
            WIPPRead, &
            WIPPDestroy

contains


! ************************************************************************** !

subroutine WIPPInit()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  if (associated(wipp)) then
    call WIPPDestroy(wipp)
  endif
  nullify(wipp)  
  
end subroutine WIPPInit

! ************************************************************************** !

function WIPPGetPtr()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  type(wipp_type), pointer :: WIPPGetPtr

  if (.not.associated(wipp)) then
    allocate(wipp)
    nullify(wipp%creep_closure)
  endif
  
  WIPPGetPtr => wipp
  
end function WIPPGetPtr

! ************************************************************************** !

subroutine WIPPRead(input,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/13/14
  ! 
  use Option_module
  use Input_Aux_module
  use String_module
  use Creep_Closure_module
  
  implicit none
  
  type(input_type), pointer :: input
  type(option_type) :: option
  
  type(wipp_type), pointer :: wipp
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string = 'WIPP'

  wipp => WIPPGetPtr()
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('CREEP_CLOSURE')
        call CreepClosureInit()
        creep_closure => CreepClosureCreate()
        call creep_closure%Read(input,option)
        option%flow%transient_porosity = PETSC_TRUE
        wipp%creep_closure => creep_closure      
     case default
        call InputKeywordUnrecognized(keyword,error_string,option)
    end select
  enddo
  
end subroutine WIPPRead

! ************************************************************************** !

subroutine WIPPDestroy1()
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  call WIPPDestroy(wipp)

end subroutine WIPPDestroy1

! ************************************************************************** !

subroutine WippDestroy2(wipp)
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/15
  !

  implicit none
  
  type(wipp_type), pointer :: wipp
  
  if (.not.associated(wipp)) return

  call CreepClosureDestroy(wipp%creep_closure)
  deallocate(wipp)
  nullify(wipp)

end subroutine WippDestroy2

end module WIPP_module
