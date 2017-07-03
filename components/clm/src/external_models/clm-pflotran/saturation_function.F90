module Saturation_Function_module
 

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
 
  type, public :: saturation_function_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: saturation_function_ctype
    PetscInt :: saturation_function_itype
    character(len=MAXWORDLENGTH) :: permeability_function_ctype
    PetscInt :: permeability_function_itype
    PetscBool :: print_me
    PetscReal, pointer :: Sr(:)
    PetscReal, pointer :: Kr0(:)
    PetscReal :: m
    PetscReal :: lambda
    PetscReal :: alpha
    PetscReal :: pcwmax
    PetscReal :: betac
    PetscReal :: power
    PetscInt :: hysteresis_id
    PetscInt :: hysteresis_params(6)
    PetscReal :: sat_spline_low
    PetscReal :: sat_spline_high
    PetscReal :: sat_spline_coefficients(3)
    PetscReal :: pres_spline_low
    PetscReal :: pres_spline_high
    PetscReal :: pres_spline_coefficients(4)
    PetscReal :: ani_A       ! parameters for anisotropic relative permeability model
    PetscReal :: ani_B       ! parameters for anisotropic relative permeability model
    PetscReal :: ani_C       ! parameters for anisotropic relative permeability model

    type(saturation_function_type), pointer :: next
  end type saturation_function_type
  
  type, public :: saturation_function_ptr_type
    type(saturation_function_type), pointer :: ptr
  end type saturation_function_ptr_type

  interface SaturationFunctionCompute
    module procedure SaturationFunctionCompute1
    module procedure SaturationFunctionCompute2
    module procedure SaturationFunctionCompute3
  end interface
  
  public :: SaturationFunctionCreate, &
            SaturationFunctionDestroy, &
            SaturationFunctionAddToList, &
            SaturationFunctionCompute, &
            SaturatFuncConvertListToArray, &
            SatFunctionComputePolynomial, &
            PermFunctionComputePolynomial, &
            SaturationFunctionRead, &
            SatFuncGetLiqRelPermFromSat, &
            SatFuncGetGasRelPermFromSat, &
            SatFuncGetCapillaryPressure, &
            SaturationFunctionGetID, &
            SatFuncComputeIcePExplicit, &
            CapillaryPressureThreshold, &
            SatFuncComputeIcePKImplicit, &
            SatFuncComputeIcePKExplicit, &
            SatFuncComputeIceDallAmico, &
            SatFuncComputeIcePKExplicitNoCryo
            
  ! Saturation function 
  PetscInt, parameter, public :: VAN_GENUCHTEN = 1
  PetscInt, parameter, public :: BROOKS_COREY = 2
  PetscInt, parameter, public :: THOMEER_COREY = 3
  PetscInt, parameter, public :: NMT_EXP = 4
  PetscInt, parameter, public :: PRUESS_1 = 5
  PetscInt, parameter, public :: LINEAR_MODEL = 6
  PetscInt, parameter, public :: VAN_GENUCHTEN_PARKER = 7
  PetscInt, parameter, public :: VAN_GENUCHTEN_DOUGHTY = 8
  PetscInt, parameter, public :: LEVERETT = 9

  ! Permeability function
  PetscInt, parameter :: DEFAULT = 0
  PetscInt, parameter, public :: BURDINE = 1
  PetscInt, parameter, public :: MUALEM = 2
  PetscInt, parameter, public :: FATT_KLIKOFF = 3

contains

! ************************************************************************** !

function SaturationFunctionCreate(option)
  ! 
  ! Creates a saturation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 
  
  use Option_module
  
  implicit none

  type(saturation_function_type), pointer :: SaturationFunctionCreate
  type(option_type) :: option
  
  type(saturation_function_type), pointer :: saturation_function
  
  allocate(saturation_function)
  saturation_function%id = 0
  saturation_function%name = ''
  saturation_function%saturation_function_ctype = 'VAN_GENUCHTEN'
  saturation_function%saturation_function_itype = VAN_GENUCHTEN
  saturation_function%permeability_function_ctype = 'MUALEM'
  saturation_function%permeability_function_itype = MUALEM
  saturation_function%print_me = PETSC_FALSE
  allocate(saturation_function%Sr(option%nphase))
  saturation_function%Sr = 0.d0
  allocate(saturation_function%Kr0(option%nphase))
  saturation_function%Kr0 = 1.0d0
  saturation_function%m = 0.d0
  saturation_function%lambda = 0.d0
  saturation_function%alpha = 0.d0
  saturation_function%pcwmax = 1.d9
  saturation_function%betac = 0.d0
  saturation_function%power = 0.d0
  saturation_function%hysteresis_id = 0
  saturation_function%hysteresis_params = 0
  saturation_function%sat_spline_low = 0.d0
  saturation_function%sat_spline_high = 0.d0
  saturation_function%sat_spline_coefficients = 0.d0
  saturation_function%pres_spline_low = 0.d0
  saturation_function%pres_spline_high = 0.d0
  saturation_function%pres_spline_coefficients = 0.d0
  saturation_function%ani_A = 0.d0
  saturation_function%ani_B = 0.d0
  saturation_function%ani_C = 0.d0
  nullify(saturation_function%next)
  SaturationFunctionCreate => saturation_function

end function SaturationFunctionCreate

! ************************************************************************** !

subroutine SaturationFunctionRead(saturation_function,input,option)
  ! 
  ! Reads in contents of a saturation_function card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  type(saturation_function_type) :: saturation_function
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscInt :: iphase
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal

  select case(option%iflowmode)
    case(G_MODE)
      option%io_buffer = 'SATURATION_FUNCTION card is no longer ' // &
        'supported for GENERAL mode.  Please use CHARACTERISTIC_' // &
        'CURVES card defined on the PFLOTRAN wiki.'
      call printErrMsg(option)
  end select
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('PERMEABILITY_FUNCTION_TYPE') 
        call InputReadWord(input,option, &
                           saturation_function%permeability_function_ctype, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'permeability function type', &
                           'SATURATION_FUNCTION')
      case('SATURATION_FUNCTION_TYPE') 
        call InputReadWord(input,option, &
                           saturation_function%saturation_function_ctype, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation function type', &
                           'SATURATION_FUNCTION')
      case('PERMEABILITY_END_POINT')
        select case(option%iflowmode)
          case(FLASH2_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','PERMEABILITY_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Kr0(iphase))
            word = trim(keyword) // 'permeabiliy end point'
            call InputErrorMsg(input,option,word,'PERMEABILITY_FUNCTION')
          case(MPH_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','PERMEABILITY_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Kr0(iphase))
            word = trim(keyword) // 'permeabiliy end point'
            call InputErrorMsg(input,option,word,'PERMEABILITY_FUNCTION')
          case(IMS_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','PERMEABILITY_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Kr0(iphase))
            word = trim(keyword) // 'permeabiliy end point'
            call InputErrorMsg(input,option,word,'PERMEABILITY_FUNCTION')
        end select
      case('RESIDUAL_SATURATION') 
        select case(option%iflowmode)
          case(FLASH2_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(MPH_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(IMS_MODE)
            call InputReadWord(input,option,keyword,PETSC_TRUE)
            call InputErrorMsg(input,option,'keyword','SATURATION_FUNCTION')
            call StringToUpper(keyword)   
            select case(trim(keyword))
              case('WATER','WATER_PHASE','LIQUID','LIQUID_PHASE')
                iphase = 1
              case('CO2','CO2_PHASE','GAS','GAS_PHASE')
                iphase = 2
              case('OIL','OIL_PHASE','NAPL','NAPL_PHASE')
                iphase = 3
            end select
            call InputReadDouble(input,option,saturation_function%Sr(iphase))
            word = trim(keyword) // ' residual saturation'
            call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
          case(G_MODE)
            iphase = 0
            string = input%buf
            call InputReadDouble(input,option,tempreal)
!            call InputErrorMsg(input,option,'residual saturation','SATURATION_FUNCTION')
            if (input%ierr /= 0) then
              input%ierr = 0
              input%buf = string
              call InputReadWord(input,option,keyword,PETSC_TRUE)
              call InputErrorMsg(input,option,'phase', &
                                 'SATURATION_FUNCTION,RESIDUAL_SATURATION')
              call StringToUpper(keyword)   
              select case(trim(keyword))
                case('LIQUID','LIQUID_PHASE')
                  iphase = 1
                case('GAS','GAS_PHASE')
                  iphase = 2
                case default
                  call InputKeywordUnrecognized(keyword, &
                    'SATURATION_FUNCTION,RESIDUAL_SATURATION',option)
              end select
              call InputReadDouble(input,option,tempreal)
              word = trim(keyword) // ' residual saturation'
              call InputErrorMsg(input,option,word,'SATURATION_FUNCTION')
            else
              ! if missing phase keyword, assume for all phases and set
              ! buffer to value
            endif
            if (iphase > 0) then
              saturation_function%Sr(iphase) = tempreal
            else
              saturation_function%Sr(:) = tempreal
            endif
          case(RICHARDS_MODE,TH_MODE)
            call InputReadDouble(input,option,saturation_function%Sr(1))
            call InputErrorMsg(input,option,'residual saturation','SATURATION_FUNCTION')
        end select
      case('LAMBDA') 
        call InputReadDouble(input,option,saturation_function%lambda)
        call InputErrorMsg(input,option,'lambda','SATURATION_FUNCTION')
        saturation_function%m = saturation_function%lambda
      case('ALPHA') 
        call InputReadDouble(input,option,saturation_function%alpha)
        call InputErrorMsg(input,option,'alpha','SATURATION_FUNCTION')
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcwmax)
        call InputErrorMsg(input,option,'maximum capillary pressure','SATURATION_FUNCTION')
      case('BETAC') 
        call InputReadDouble(input,option,saturation_function%betac)
        call InputErrorMsg(input,option,'betac','SATURATION_FUNCTION')
      case('POWER') 
        call InputReadDouble(input,option,saturation_function%power)
        call InputErrorMsg(input,option,'power','SATURATION_FUNCTION')
      case('ANI_A') 
        call InputReadDouble(input,option,saturation_function%ani_A)
        call InputErrorMsg(input,option,'ani_A','SATURATION_FUNCTION')
      case('ANI_B') 
        call InputReadDouble(input,option,saturation_function%ani_B)
        call InputErrorMsg(input,option,'ani_B','SATURATION_FUNCTION')
      case('ANI_C') 
        call InputReadDouble(input,option,saturation_function%ani_C)
        call InputErrorMsg(input,option,'ani_C','SATURATION_FUNCTION')
      case('VERIFY') 
        saturation_function%print_me = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'SATURATION_FUNCTION',option)
    end select 
  
  enddo 
  
  if (saturation_function%m < 1.d-40 .and. &
      .not. &
      (StringCompare(saturation_function%name,'default',SEVEN_INTEGER) .or. &
       StringCompareIgnoreCase(saturation_function%saturation_function_ctype, &
                               'leverett',EIGHT_INTEGER))) then
    option%io_buffer = 'Saturation function parameter "m" not set ' // &
                       'properly in saturation function "' // &
                       trim(saturation_function%name) // '".'
    call printErrMsg(option)
  endif
  
  call SaturationFunctionSetTypes(saturation_function,option)

end subroutine SaturationFunctionRead

! ************************************************************************** !

subroutine SaturationFunctionSetTypes(saturation_function,option)
  ! 
  ! Sets the integer values that map the saturation and permeability  
  ! functions to a specific type
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/23/14

  use Option_module
  use String_module
  
  implicit none

  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  ! set permeability function integer type
  call StringToUpper(saturation_function%permeability_function_ctype)
  select case(trim(saturation_function%permeability_function_ctype))
    case('DEFAULT')
      saturation_function%permeability_function_itype = DEFAULT
    case('BURDINE')
      saturation_function%permeability_function_itype = BURDINE
    case('MUALEM')
      saturation_function%permeability_function_itype = MUALEM
    case('NMT_EXP')
      saturation_function%permeability_function_itype = NMT_EXP
    case('PRUESS_1')
      saturation_function%permeability_function_itype = PRUESS_1
    case('VAN_GENUCHTEN_PARKER')
      saturation_function%permeability_function_itype = VAN_GENUCHTEN_PARKER
    case('FATT_KLIKOFF')
      saturation_function%permeability_function_itype = FATT_KLIKOFF
    case default
      option%io_buffer = 'Permeability function type "' // &
                          trim(saturation_function%permeability_function_ctype) // &
                          '" not recognized ' // &
                          ' in saturation function ' // &
                          trim(saturation_function%name)
      call printErrMsg(option)
  end select
    
  ! set saturation function integer type
  call StringToUpper(saturation_function%saturation_function_ctype)
  select case(trim(saturation_function%saturation_function_ctype))
    case('VAN_GENUCHTEN')
      saturation_function%saturation_function_itype = VAN_GENUCHTEN
    case('BROOKS_COREY')
      saturation_function%saturation_function_itype = BROOKS_COREY
    case('LINEAR_MODEL')
      saturation_function%saturation_function_itype = LINEAR_MODEL
    case('THOMEER_COREY')
      saturation_function%saturation_function_itype = THOMEER_COREY
    case('NMT_EXP')
      saturation_function%saturation_function_itype = NMT_EXP
    case('PRUESS_1')
      saturation_function%saturation_function_itype = PRUESS_1
    case('VAN_GENUCHTEN_PARKER')
      saturation_function%saturation_function_itype = VAN_GENUCHTEN_PARKER
    case('VAN_GENUCHTEN_DOUGHTY')
      saturation_function%saturation_function_itype = VAN_GENUCHTEN_DOUGHTY
    case('LEVERETT')
      saturation_function%saturation_function_itype = LEVERETT
    case default
      option%io_buffer = 'Saturation function type "' // &
                          trim(saturation_function%saturation_function_ctype) // &
                          '" not recognized ' // &
                          ' in saturation function ' // &
                          trim(saturation_function%name)
      call printErrMsg(option)
  end select  

end subroutine SaturationFunctionSetTypes

! ************************************************************************** !

subroutine SatFunctionComputePolynomial(option,saturation_function)
  ! 
  ! Computes a spline spanning the
  ! discontinuity in Brooks Corey
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/27/08
  ! 
  
  use Option_module
  use Utility_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  
  PetscReal :: b(4)
  PetscReal :: pressure_1, pressure_2
  PetscReal :: saturation_1, saturation_2
  
  PetscReal :: n

  select case(saturation_function%saturation_function_itype) 
    case(BROOKS_COREY)

      ! polynomial fitting pc as a function of saturation
      ! 1.05 is essentially pc*alpha (i.e. pc = 1.05/alpha)
      saturation_1 = 1.05d0**(-saturation_function%lambda)
      saturation_2 = 1.d0
  
      saturation_function%sat_spline_low = saturation_1
      saturation_function%sat_spline_high = saturation_2

      b = 0.d0
      ! fill right hand side
      ! capillary pressure at 1
      b(1) = 1.05d0/saturation_function%alpha 
      ! capillary pressure at 2
      b(2) = 0.d0
      ! derivative of pressure at saturation_1
      ! pc = Se**(-1/lambda)/alpha
      ! dpc_dSe = -1/lambda*Se**(-1/lambda-1)/alpha
      b(3) = -1.d0/saturation_function%lambda* &
             saturation_1**(-1.d0/saturation_function%lambda-1.d0)/ &
             saturation_function%alpha

      call QuadraticPolynomialSetup(saturation_1,saturation_2,b(1:3), &
                                    ! indicates derivative given at 1
                                    PETSC_TRUE) 
      
      saturation_function%sat_spline_coefficients(1:3) = b(1:3)

      ! polynomial fitting saturation as a function of pc
      !geh: cannot invert the pressure/saturation relationship above
      !     since it can result in saturations > 1 with both
      !     quadratic and cubic polynomials
      ! fill matix with values
      pressure_1 = 0.95/saturation_function%alpha
      pressure_2 = 1.05/saturation_function%alpha
      
      saturation_function%pres_spline_low = pressure_1
      saturation_function%pres_spline_high = pressure_2
  
      b = 0.d0
      ! Se at 1
      b(1) = 1.d0
      ! Se at 2
      b(2) = (pressure_2*saturation_function%alpha)** &
             (-saturation_function%lambda)
      ! derivative of Se at 1
      b(3) = 0.d0 
      ! derivative of Se at 2
      b(4) = -saturation_function%lambda/pressure_2* &
               (pressure_2*saturation_function%alpha)** &
                 (-saturation_function%lambda)

      call CubicPolynomialSetup(pressure_1,pressure_2,b)

      saturation_function%pres_spline_coefficients(1:4) = b(1:4)
      
  case(VAN_GENUCHTEN)
 
      ! return for now
      return

      !geh: keep for now
#if 0
      ! geh: needs to be refactored to consider capillary pressure
      ! fill matix with values
      ! these are capillary pressures
      pressure_low = 0  ! saturated
      pressure_high = 0.01d0*option%reference_pressure  ! just below saturated
  
      saturation_function%spline_low = pressure_low
      saturation_function%spline_high = pressure_high
    
      n = 1.d0/(1.d0 - saturation_function%m)
      b(1) = (1.d0+(pressure_high*saturation_function%alpha)**n)** &
               (-saturation_function%m)
      b(2) = 1.d0
      b(3) = -saturation_function%m*n*saturation_function%alpha* &
             (saturation_function%alpha*pressure_high)**(n-1.d0)* &
             (1.d0+(saturation_function%alpha*pressure_high)**n)** &
               (-saturation_function%m-1.d0)
      b(4) = 0.d0
  
      call CubicPolynomialSetup(pressure_high,pressure_low,b)

      saturation_function%spline_coefficients(1:4) = b(1:4)
#endif

  end select
  
end subroutine SatFunctionComputePolynomial

! ************************************************************************** !

subroutine PermFunctionComputePolynomial(option,saturation_function)
  ! 
  ! Computes a spline spanning the
  ! discontinuity in Brooks Corey
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/27/12
  ! 
  
  use Option_module
  use Utility_module
  
  implicit none

  type(option_type) :: option
  type(saturation_function_type) :: saturation_function
  
  PetscReal :: b(4)
  PetscReal :: Se_high, Se_low, one_over_m, Se_one_over_m, m
  
  select case(saturation_function%saturation_function_itype) 

    case(BROOKS_COREY)

    case(VAN_GENUCHTEN)
 
#ifdef MUALEM_SPLINE
      ! fill matix with values
      Se_low = 0.99d0  ! saturated
      Se_high = 1.d0  ! just below saturated
  
      saturation_function%spline_low = Se_low
      saturation_function%spline_high = Se_high
    
      m = saturation_function%m
      one_over_m = 1.d0/m
      Se_one_over_m = Se_low**one_over_m
      b(1) = 1.d0
      b(2) = sqrt(Se_low)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
      b(3) = 0.d0
      b(4) = 0.5d0*b(2)/Se_low+ &
             2.d0*Se_low**(one_over_m-0.5d0)* &
             (1.d0-Se_one_over_m)**(m-1.d0)* &
             (1.d0-(1.d0-Se_one_over_m)**m)
  
      call CubicPolynomialSetup(Se_high,Se_low,b)
  
      saturation_function%spline_coefficients(1:4) = b(1:4)
#endif

  end select
  
end subroutine PermFunctionComputePolynomial

! ************************************************************************** !

subroutine SaturationFunctionAddToList(saturation_function,list)
  ! 
  ! Adds a saturation function to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  type(saturation_function_type), pointer :: list

  type(saturation_function_type), pointer :: cur_saturation_function
  
  if (associated(list)) then
    cur_saturation_function => list
    ! loop to end of list
    do
      if (.not.associated(cur_saturation_function%next)) exit
      cur_saturation_function => cur_saturation_function%next
    enddo
    cur_saturation_function%next => saturation_function
  else
    list => saturation_function
  endif
  
end subroutine SaturationFunctionAddToList

! ************************************************************************** !

subroutine SaturatFuncConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the
  ! saturation functions in the list
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/11/07
  ! 

  use String_module
  use Option_module
  
  implicit none
  
  type(saturation_function_type), pointer :: list
  type(saturation_function_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  type(saturation_function_type), pointer :: cur_saturation_function
  PetscInt :: count

  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    cur_saturation_function => cur_saturation_function%next
  enddo
  
  if (associated(array)) deallocate(array)
  allocate(array(count))
  
  count = 0
  cur_saturation_function => list
  do 
    if (.not.associated(cur_saturation_function)) exit
    count = count + 1
    cur_saturation_function%id = count
    array(count)%ptr => cur_saturation_function
    if (cur_saturation_function%print_me .and. &
        option%myrank == option%io_rank) then
      call SaturationFunctionVerify(cur_saturation_function,option)
    endif
    cur_saturation_function => cur_saturation_function%next
  enddo

end subroutine SaturatFuncConvertListToArray

! ************************************************************************** !

subroutine SaturationFunctionCompute1(capillary_pressure,saturation, &
                                      relative_perm, &
                                      dsat_dpres,dkr_dpres, &
                                      saturation_function, &
                                      auxvar1,auxvar2, &
                                      option)
  ! 
  ! Computes the saturation and relative permeability
  ! (and associated derivatives) as a function of
  ! capillary pressure
  ! 
  ! Author: Glenn Hammond
  ! Date: 2/9/12
  ! 
  use Option_module
  
  implicit none

  PetscReal :: capillary_pressure, saturation, relative_perm
  PetscReal :: dsat_dpres, dkr_dpres
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1,auxvar2
  type(option_type) :: option

  PetscBool :: switch_to_saturated
  
  call SaturationFunctionCompute2(capillary_pressure,saturation, &
                                  relative_perm, &
                                  dsat_dpres,dkr_dpres, &
                                  saturation_function, &
                                  auxvar1,auxvar2, &
                                  switch_to_saturated,option)

end subroutine SaturationFunctionCompute1

! ************************************************************************** !

subroutine SaturationFunctionCompute2(capillary_pressure,saturation, &
                                      relative_perm, &
                                      dsat_dpres,dkr_dpres, &
                                      saturation_function, &
                                      auxvar1,auxvar2, &
                                      switch_to_saturated,option)
  ! 
  ! Computes the saturation and relative permeability
  ! (and associated derivatives) as a function of
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !                                    
  ! Author: Glenn Hammond
  ! Date: 12/11/07
  ! 
  use Option_module
  use Utility_module, only:CubicPolynomialEvaluate
  
  implicit none

  PetscReal :: capillary_pressure, saturation, relative_perm 
  PetscReal :: dsat_dpres, dkr_dpres
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1,auxvar2
  PetscBool :: switch_to_saturated
  type(option_type) :: option

  PetscInt :: iphase
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha, pcmax
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_dpc, dsat_dpc, dkr_dpc
  PetscReal :: dkr_dSe, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: por, perm
  PetscReal :: Fg, a, Pd, PHg, pct_over_pcmax, pc_over_pcmax, pc_log_ratio

  PetscReal, parameter :: pc_alpha_n_epsilon = 1.d-15
  
  iphase = 1
  dsat_dpres = 0.d0
  dkr_dpres = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      ! reference #1
#if 0
      if (pc < saturation_function%spline_low) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      else if (pc < saturation_function%spline_high) then
        Sr = saturation_function%Sr(iphase)
        call CubicPolynomialEvaluate(saturation_function%spline_coefficients, &
                                     pc,Se,dSe_dpc)
        saturation = Sr + (1.d0-Sr)*Se
        dsat_dpc = (1.d0-Sr)*dSe_dpc
#else
      if (capillary_pressure <= 0.d0) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
#endif        
      else
        alpha = saturation_function%alpha
        pc = capillary_pressure
        m = saturation_function%m
        n = 1.d0/(1.d0-m)
        pc_alpha = pc*alpha
        pc_alpha_n = pc_alpha**n
        !geh:  This conditional does not catch potential cancelation in 
        !      the dkr_dSe deriviative calculation.  Therefore, I am setting
        !      an epsilon here
!        if (1.d0 + pc_alpha_n == 1.d0) then ! check for zero perturbation
        if (pc_alpha_n < pc_alpha_n_epsilon) then 
          saturation = 1.d0
          relative_perm = 1.d0
          switch_to_saturated = PETSC_TRUE
          return
        endif
        one_plus_pc_alpha_n = 1.d0+pc_alpha_n
        Sr = saturation_function%Sr(iphase)
        Se = one_plus_pc_alpha_n**(-m)
        dSe_dpc = -m*n*alpha*pc_alpha_n/(pc_alpha*one_plus_pc_alpha_n**(m+1))
        saturation = Sr + (1.d0-Sr)*Se
        dsat_dpc = (1.d0-Sr)*dSe_dpc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'vG Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'vG Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = Se*Se*(1.d0-(1.d0-Se_one_over_m)**m)
          dkr_dSe = 2.d0*relative_perm/Se + &
                   Se*Se_one_over_m*(1.d0-Se_one_over_m)**(m-1.d0)
          dkr_dpc = dkr_dSe*dSe_dpc
        case(MUALEM)
          ! reference #1
#ifdef MUALEM_SPLINE
          if (Se > saturation_function%spline_low) then
            call CubicPolynomialEvaluate( &
              saturation_function%spline_coefficients, &
              Se,relative_perm,dkr_dSe)
          else
#endif          
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
          dkr_dSe = 0.5d0*relative_perm/Se+ &
                   2.d0*Se**(one_over_m-0.5d0)* &
                        (1.d0-Se_one_over_m)**(m-1.d0)* &
                        (1.d0-(1.d0-Se_one_over_m)**m)
#ifdef MUALEM_SPLINE
          endif
#endif          
          dkr_dpc = dkr_dSe*dSe_dpc
        case default
          option%io_buffer = 'Unknown relative permeabilty function' 
          call printErrMsg(option)
      end select
    case(BROOKS_COREY)
      ! reference #1
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      pc = capillary_pressure
      if (pc < saturation_function%pres_spline_low) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      else if (pc < saturation_function%pres_spline_high) then
        lambda = saturation_function%lambda
        Sr = saturation_function%Sr(iphase)
        call CubicPolynomialEvaluate(saturation_function% &
                                     pres_spline_coefficients, &
                                     pc,Se,dSe_dpc)
        saturation = Sr + (1.d0-Sr)*Se
        dsat_dpc = (1.d0-Sr)*dSe_dpc
      else
        lambda = saturation_function%lambda
        Sr = saturation_function%Sr(iphase)
        pc_alpha_neg_lambda = (pc*alpha)**(-lambda)
        Se = pc_alpha_neg_lambda
        dSe_dpc = -lambda/pc*pc_alpha_neg_lambda
        saturation = Sr + (1.d0-Sr)*Se
!        dsat_dpc = -lambda*(1.d0-Sr)/pc*pc_alpha_neg_lambda
        dsat_dpc = (1.d0-Sr)*dSe_dpc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'BC Saturation > 1:', saturation
      else if (saturation > 1.d0 .or. saturation < Sr) then
        print *, option%myrank, 'BC Saturation < Sr:', saturation, Sr
      endif
      ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          power = 3.d0+2.d0/lambda
          relative_perm = Se**power
          dkr_dSe = power*relative_perm/Se
          dkr_dpc = dkr_dSe*dSe_dpc
        case(MUALEM)
          ! reference #1
          power = 2.5d0+2.d0/lambda
          relative_perm = Se**power
          dkr_dSe = power*relative_perm/Se
          dkr_dpc = dkr_dSe*dSe_dpc
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(LINEAR_MODEL)
      ! Added by Bwalya Malama 01/30/2014
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      lambda = saturation_function%lambda
      pcmax = saturation_function%pcwmax
      pc = capillary_pressure

      if (capillary_pressure <= 0.d0) then
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      else
        Sr = saturation_function%Sr(iphase)
        Se = (pcmax-pc)/(pcmax-one_over_alpha)
        dSe_dpc = -1.d0/(pcmax-one_over_alpha)
        saturation = Sr + (1.d0-Sr)*Se
        dsat_dpc = (1.d0-Sr)*dSe_dpc
      endif
      if (saturation > 1.d0) then
        print *, option%myrank, 'BC Saturation > 1:', saturation
      else if (saturation < Sr) then
        print *, option%myrank, 'BC Saturation < Sr:', saturation, Sr
      endif
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          relative_perm = Se
          dkr_dSe = 1.d0
        case(MUALEM)
          power = 5.d-1
          pct_over_pcmax = one_over_alpha/pcmax
          pc_over_pcmax = 1.d0-(1.d0-pct_over_pcmax)*Se
          pc_log_ratio = log(pc_over_pcmax)/log(pct_over_pcmax)
          relative_perm = (Se**power)*(pc_log_ratio**2.d0)
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(THOMEER_COREY)
      pc = capillary_pressure
      por = auxvar1
      perm = auxvar2*1.013202d15 ! convert from m^2 to mD
      Fg = saturation_function%alpha
      a = saturation_function%m
      Pd = 100.d0*por/sqrt(perm/(3.8068d0*(Fg**(-1.334d0)))) ! psi
      PHg = 9.63051d-4*pc
      if (PHg > Pd) then
        saturation = 1.d0-exp(-Fg/log10(PHg/Pd))
#if 0
        alpha = pc*(1.d0+1.d-8)
        m = 9.63051d-4*alpha
        n = 1.d0-exp(-Fg/log10(m/Pd))
        n = (n-saturation)/(alpha-pc)
#endif        
        dsat_dpc = (saturation-1.d0)*Fg/(log10(PHg/Pd)**2.d0)/(pc*2.30258509d0)
        ! Sr assumed to be zero
        relative_perm = saturation**a
        dkr_dpc = a*saturation**(a-1.d0)*dsat_dpc
      else
        saturation = 1.d0
        relative_perm = 1.d0
        switch_to_saturated = PETSC_TRUE
        return
      endif
    case default
      option%io_buffer = 'Unknown saturation function'
      call printErrMsg(option)
  end select

  dsat_dpres = -dsat_dpc 
  dkr_dpres = -dkr_dpc

end subroutine SaturationFunctionCompute2

! ************************************************************************** !

subroutine SaturationFunctionCompute3(capillary_pressure,saturation, &
                                      saturation_function, &
                                      option)
  ! 
  ! Computes just saturation as a function of
  ! capillary pressure
  ! 
  ! Author: Glenn Hammond
  ! Date: 2/9/12
  ! 
  use Option_module
  
  implicit none

  PetscReal :: capillary_pressure, saturation
  PetscReal :: relative_perm_dummy, dsat_dpres_dummy, dkr_dpres_dummy
  type(saturation_function_type) :: saturation_function
  PetscReal :: auxvar1_dummy, auxvar2_dummy
  type(option_type) :: option

  PetscBool :: switch_to_saturated_dummy
  
  call SaturationFunctionCompute2(capillary_pressure,saturation, &
                                  relative_perm_dummy, &
                                  dsat_dpres_dummy,dkr_dpres_dummy, &
                                  saturation_function, &
                                  auxvar1_dummy,auxvar2_dummy, &
                                  switch_to_saturated_dummy,option)

end subroutine SaturationFunctionCompute3

! ************************************************************************** !

subroutine SatFuncComputeIcePExplicit(liquid_pressure, temperature, &
                                      ice_saturation, &
                                      liquid_saturation, gas_saturation, &
                                      liquid_relative_perm, dsl_pl, & 
                                      dsl_temp, dsg_pl, dsg_temp, dsi_pl, &
                                      dsi_temp, dkr_pl, dkr_temp, &
                                      saturation_function, pth, option)
  ! 
  ! Computes the saturation of ice, water vapor
  ! and liquid water given the saturation function
  ! temperature, water vapor pressure and liquid
  ! water pressure
  !
  ! Model based on Painter, Comp. Geosci, 2011
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 11/14/11
  ! 

  use Option_module
  use PFLOTRAN_Constants_module
 
implicit none

  PetscReal :: liquid_pressure, temperature
  PetscReal :: ice_saturation, liquid_saturation, gas_saturation
  PetscReal :: liquid_relative_perm
  PetscReal :: dsl_pl, dsl_temp
  PetscReal :: dsg_pl, dsg_temp
  PetscReal :: dsi_pl, dsi_temp
  PetscReal :: dkr_pl
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscReal :: alpha, lambda, m, n
  PetscReal :: pc, Se, one_over_m, Se_one_over_m, dSe_dpc, dkr_dpc
  PetscReal :: dkr_dSe, power
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: function_A, function_B
  PetscReal :: pc_il, gamma, pc_il_alpha, pc_il_alpha_n, Se_temp
  PetscReal :: one_plus_pc_il_alpha_n
  PetscReal :: dfunc_A_temp
  PetscReal :: dfunc_B_pl
  PetscReal :: liq_sat_one_over_m, dkr_ds_liq, dkr_temp
  PetscReal :: pth, dSe_dpc_at_pth
  PetscReal, parameter :: den_ice = 9.167d2 !in kg/m3 at 273.15K
  PetscReal, parameter :: interfacial_tensions_ratio = 2.33
  PetscReal, parameter :: T_0 = 273.15d0 !in K
  
  PetscReal :: dsi_dpl, dsg_dpl, dsl_dpl
  PetscReal :: dsi_dT, dsg_dT, dsl_dT
            
  dsl_pl = 0.d0
  dsl_temp = 0.d0
  dsg_pl = 0.d0
  dsg_temp = 0.d0
  dsi_pl = 0.d0
  dsi_temp = 0.d0
  dkr_pl = 0.d0
  dkr_temp = 0.d0
  dkr_ds_liq = 0.d0
  
  ! compute saturation
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      if (liquid_pressure >= option%reference_pressure) then
        function_B = 1.d0
        dfunc_B_pl = 0.d0
      else
        alpha = saturation_function%alpha
        pc = option%reference_pressure - liquid_pressure
        m = saturation_function%m
        n = 1.d0/(1.d0 - m)
        pc_alpha = pc*alpha
        pc_alpha_n = pc_alpha**n
        one_plus_pc_alpha_n = 1.d0 + pc_alpha_n
        Se = one_plus_pc_alpha_n**(-m)
        dSe_dpc = -m*n*alpha*pc_alpha_n/(pc_alpha*one_plus_pc_alpha_n**(m+1))
        if (pc >= pth) then
          dSe_dpc_at_pth = -m*n*(1.d0 + (alpha*pth)**n)**(-1.d0-m)*(alpha**n*pth**(n-1.d0))
          Se = (pc - 1.d8)*dSe_dpc_at_pth
          dSe_dpc = dSe_dpc_at_pth
        ! write (*,*) option%myrank, 'pc:', pc, 'Se:', Se, 'dSe_dpc', dSe_dpc 
        endif 
        function_B = 1.d0/Se
        dfunc_B_pl = 1.d0/(Se**(2.d0))*dSe_dpc        
      endif
      if (temperature >= 0.d0) then
        function_A = 1.d0
        dfunc_A_temp = 0.d0
      else
        gamma = den_ice*HEAT_OF_FUSION*interfacial_tensions_ratio
        pc_il = gamma*(-(temperature))/T_0
        alpha = saturation_function%alpha
        m = saturation_function%m
        n = 1.d0/(1.d0 - m)
        pc_il_alpha = pc_il*alpha
        pc_il_alpha_n = pc_il_alpha**n
        one_plus_pc_il_alpha_n = 1.d0 + pc_il_alpha_n
        Se_temp = one_plus_pc_il_alpha_n**(-m)
        function_A = 1.d0/Se_temp
        dfunc_A_temp = (gamma/T_0)*1.d0/(Se_temp**(2.d0))*(-m)* &
                       ((one_plus_pc_il_alpha_n)**(-m - 1.d0))*n* &
                       (pc_il**(n - 1.d0))*(alpha**n)
      endif           
    case default
      option%io_buffer = 'Ice module only supports Van Genuchten'
      call printErrMsg(option)
  end select
  
  liquid_saturation = 1.d0/(function_A + function_B - 1.d0)
  gas_saturation = liquid_saturation*(function_B - 1.d0)
  ice_saturation = liquid_saturation*(function_A - 1.d0)

  dsl_pl = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_B_pl)
  dsl_temp = - 1.d0/(function_A + function_B - 1.d0)**(2.d0)*(dfunc_A_temp)
  
  dsg_pl = dsl_pl*(function_B - 1.d0) + liquid_saturation*dfunc_B_pl
  dsg_temp = dsl_temp*(function_B - 1.d0)
  
  dsi_pl = dsl_pl*(function_A - 1.d0)
  dsi_temp = dsl_temp*(function_A - 1.d0) + liquid_saturation*dfunc_A_temp
  
  if (liquid_saturation > 1.d0) then
    print *, option%myrank, 'vG Liquid Saturation > 1:', liquid_saturation
  else if (liquid_saturation < 0.d0) then
    print *, option%myrank, 'vG Liquid Saturation < 0:', liquid_saturation
  endif

  if (gas_saturation > 1.d0) then
    print *, option%myrank, 'vG Gas Saturation > 1:', gas_saturation
  else if (gas_saturation < 0.d0) then
    print *, option%myrank, 'vG Gas Saturation < 0:', gas_saturation
  endif
 
  if (ice_saturation > 1.d0) then
    print *, option%myrank, 'vG Ice Saturation > 1:', ice_saturation
  else if (ice_saturation < 0.d0) then
    print *, option%myrank, 'vG Ice Saturation < 0:', ice_saturation
  endif

  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      if (liquid_saturation == 1.d0) then
        liquid_relative_perm = 1.d0
        dkr_ds_liq = 0.d0
      else
        m = saturation_function%m
        one_over_m = 1.d0/m
        liq_sat_one_over_m = liquid_saturation**one_over_m
        liquid_relative_perm = sqrt(liquid_saturation)* &
                               (1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
        dkr_ds_liq = 0.5d0*liquid_relative_perm/liquid_saturation + &
                     2.d0*liquid_saturation**(one_over_m - 0.5d0)* &
                     (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                     (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
      endif
        dkr_pl = dkr_ds_liq*dsl_pl
        dkr_temp = dkr_ds_liq*dsl_temp
    case default
      option%io_buffer = 'Ice module only supports Mualem' 
      call printErrMsg(option)
  end select
   
end subroutine SatFuncComputeIcePExplicit

! ************************************************************************** !

subroutine ComputeVGAndDerivative(alpha,lambda,Pc,S,dS_dPc)
  ! 
  ! Evaluates van Genunchten saturation function and
  ! its derivative with respect to capillary pressure
  ! for given capillary pressure
  !
  ! Note: Derivative wrt capillary pressure and not liquid pressure 
  ! is evaluated
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  implicit none

  PetscReal :: alpha, lambda, gamma  
  PetscReal :: Pc, S, dS_dPc
  
  gamma = 1.d0/(1.d0 - lambda)
  if (Pc > 0.d0) then
    S =  (1.d0 + (alpha*Pc)**gamma)**(-lambda)
    dS_dPc = (-lambda)*((1.d0 + (alpha*Pc)**gamma)**(-lambda - 1.d0))* &
         (gamma*alpha*(alpha*Pc)**(gamma - 1.d0))
  else
    S = 1.d0
    dS_dPc = 0.d0
  endif
 
end subroutine ComputeVGAndDerivative

! ************************************************************************** !

subroutine ComputeInvVGAndDerivative(alpha,lambda,sat,Sinv,dSinv_dsat)
  ! 
  ! Evaluates inverse of van Genunchten saturation function
  ! and its derivative with respect to saturation at a given saturation
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  implicit none
  
  PetscReal :: alpha, lambda, gamma
  PetscReal :: sat, Sinv, dSinv_dsat
  
  gamma = 1.d0/(1.d0 - lambda)
  if (sat == 1.d0) then
    Sinv = 0.d0
    dSinv_dsat = 0.d0
  else
    Sinv = 1.d0/alpha*((sat)**(-1.d0/lambda) - 1.d0)**(1.d0/gamma)
    dSinv_dsat = 1.d0/alpha*1.d0/gamma*((sat**(-1.d0/lambda) - 1.d0)**(1.d0/gamma - &
      1.d0))*(-1.d0/lambda)*(sat**(-1.d0/lambda - 1.d0))
  endif
  
end subroutine ComputeInvVGAndDerivative

! ************************************************************************** !

subroutine CalculateImplicitIceFunc(alpha,lambda,Pcgl,T,s_i,func_val,dfunc_val)
  ! 
  ! Evaluates the value of implicit equation whose
  ! solution is ice saturation
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  implicit none
  
  PetscReal :: alpha, lambda
  PetscReal :: Pcgl, T, s_i, func_val
  PetscReal :: temp_term, PC 
  PetscReal :: sat, dsat, sat_term
  PetscReal :: sat_inv, dsat_inv
  PetscReal :: sat_PC, dsat_dpc, dfunc_val
  PetscReal, parameter :: beta = 2.33          ! dimensionless -- ratio of surf. tens
  PetscReal, parameter :: rho_i = 9.167d2      ! in kg/m^3
  PetscReal, parameter :: T_0 = 273.15         ! in K
  
  if (T >= 0.d0) then   ! T is in C
    temp_term = 0.d0
  else
    temp_term = -beta*rho_i*HEAT_OF_FUSION*T/T_0
  endif
  call ComputeVGAndDerivative(alpha,lambda,Pcgl,sat,dsat)
  sat_term = (s_i + (1.d0 - s_i)*sat)
  call ComputeInvVGAndDerivative(alpha,lambda,sat_term,sat_inv,dsat_inv)
  PC = temp_term + sat_inv
  call ComputeVGAndDerivative(alpha,lambda,PC,sat_PC,dsat_dPC)
  func_val = (1.d0 - s_i)*sat - sat_PC
  dfunc_val = -sat - dsat_dpc*dsat_inv*(1.d0 - sat)
  
end subroutine CalculateImplicitIceFunc

! ************************************************************************** !

subroutine CalcPhasePartitionIceNewt(alpha,lambda,Pcgl,T,s_g,s_i,s_l)
  ! 
  ! Solves the implicit constitutive relation
  ! to calculate saturations of ice, liquid
  ! and vapor phases
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  implicit none
  PetscReal :: alpha, lambda
  PetscReal :: Pcgl, T, s_g, s_i, s_l
  PetscReal :: func_val, dfunc_val 
  PetscReal :: x, x_new, sat, dsat
  PetscInt :: iter
  PetscReal, parameter :: eps = 1.d-12
  PetscInt, parameter :: maxit = 10
  
  
  if (T >= 0.d0) then
    x = 0.d0
  else
    x = 5.d-1          ! Initial guess
    do iter = 1,maxit
      call CalculateImplicitIceFunc(alpha,lambda,Pcgl,T,x,func_val,dfunc_val)
!       print *, 'iteration:', iter, 'value:', x, 'inormr:', abs(func_val), &
!      'dfunc_val:', dfunc_val, 'dfunc_val_num:', dfunc_val_num
      if (abs(func_val) < eps) exit
      x_new = x - func_val/dfunc_val
      if (x_new >= 1.d0) then
        x_new = 1.d0 - 1.d-8
      endif
      if (x_new <= 0.d0) then
        x_new = 1.d-8
      endif
      x = x_new
    enddo
  endif
  
  call ComputeVGAndDerivative(alpha,lambda,Pcgl,sat,dsat)
  s_i = x
  s_l = (1.d0 - x)*sat
  s_g = 1.d0 - s_l - s_i

end subroutine CalcPhasePartitionIceNewt

! ************************************************************************** !

subroutine CalcPhasePartitionIceBis(alpha,lambda,Pcgl,T,s_g,s_i,s_l)
  ! 
  ! Solves the implicit constitutive relation
  ! to calculate saturations of ice, liquid
  ! and vapor phases using bisection method
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  implicit none
  
  PetscReal :: alpha, lambda
  PetscReal :: Pcgl, T, s_g, s_i, s_l
  PetscReal :: max_s, min_s
  PetscReal :: x, F_max_s, F_min_s, sol, F_x
  PetscReal :: sat, dsat
  PetscInt :: i
  PetscReal :: dF_min_s, dF_max_s, dF_x
  PetscInt, parameter :: maxit = 100
  PetscReal, parameter :: tol = 1.d-12
  
  max_s = 1.d0
  min_s = 0.d0
  
  if (T >= 0.d0) then
    sol = 0.d0
  else
    do i = 1, maxit
      call CalculateImplicitIceFunc(alpha,lambda,Pcgl,T,min_s,F_min_s,dF_min_s)
      call CalculateImplicitIceFunc(alpha,lambda,Pcgl,T,max_s,F_max_s,dF_max_s)
      if (abs(F_min_s) < tol) then
        sol = min_s
        exit
      endif
      if (abs(F_max_s) < tol) then
        sol = max_s
        exit
      endif
      x = min_s + (max_s - min_s)/2.d0
      call CalculateImplicitIceFunc(alpha,lambda,Pcgl,T,x,F_x,dF_x)
      ! print *, i, x, F_x, max_s, min_s
      if (abs(F_x) < tol) then
        sol = x
        exit
      endif
      if (F_min_s*F_x > 0.d0) then
        min_s = x
      else
        max_s = x
      endif
    enddo
  endif
    
  call ComputeVGAndDerivative(alpha,lambda,Pcgl,sat,dsat)
  s_i = sol
  s_l = (1.d0 - sol)*sat
  s_g = 1.d0 - s_l - s_i

end subroutine  CalcPhasePartitionIceBis

! ************************************************************************** !

subroutine CalcPhasePartitionIceDeriv(alpha,lambda,Pcgl,T,s_g,s_i,s_l, &
                                      dsg_dpl,dsg_dT,dsi_dpl,dsi_dT, &
                                      dsl_dpl,dsl_dT)
  ! 
  ! Solves the implicit constitutive relation
  ! to calculate saturations of ice, liquid
  ! and vapor phases
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 
                                          
  implicit none
  
  PetscReal :: alpha, lambda
  PetscReal :: Pcgl, T
  PetscReal :: dsg_dpl, dsg_dT
  PetscReal :: dsi_dpl, dsi_dT
  PetscReal :: dsl_dpl, dsl_dT
  PetscReal :: s_g, s_i, s_l
  PetscReal :: sat, dsat
  PetscReal :: sat_inv, dsat_inv
  PetscReal :: PC, sat_PC, dsat_dpc
  PetscReal :: G, dS_dpl, temp_term
  PetscReal :: L, M, N
  PetscReal, parameter :: beta = 2.33          ! dimensionless -- ratio of surf. tens
  PetscReal, parameter :: rho_i = 9.167d2      ! in kg/m^3
  PetscReal, parameter :: T_0 = 273.15         ! in K

#if 0
  PetscReal :: dsi_dpl_num, dsi_dT_num
  PetscReal :: dsg_dpl_num, dsg_dT_num
  PetscReal :: dsl_dpl_num, dsl_dT_num
  PetscReal :: s_g_pinc, s_i_pinc, s_l_pinc
  PetscReal :: s_g_Tinc, s_i_Tinc, s_l_Tinc
  PetscReal, parameter :: delta = 1.d-8
#endif

  dsi_dpl = 0.d0
  dsi_dT = 0.d0
  dsg_dpl = 0.d0
  dsg_dT = 0.d0
  dsl_dpl = 0.d0
  dsl_dT = 0.d0

#if 0
  dsi_dpl_num = 0.d0
  dsg_dpl_num = 0.d0
  dsl_dpl_num = 0.d0
  dsi_dT_num = 0.d0
  dsg_dT_num = 0.d0
  dsl_dT_num = 0.d0
#endif
  
  ! Calculate the derivatives of saturation with respect to pl
  call CalcPhasePartitionIceNewt(alpha,lambda,Pcgl,T,s_g,s_i,s_l)
  if (T >= 0.d0) then
    temp_term = 0.d0
  else
    temp_term = -beta*rho_i*HEAT_OF_FUSION*T/T_0
  endif
  call ComputeInvVGAndDerivative(alpha,lambda,(s_i + s_l),sat_inv,dsat_inv)
  PC = temp_term + sat_inv 
  call ComputeVGAndDerivative(alpha,lambda,PC,sat_PC,dsat_dPC)
  call ComputeVGAndDerivative(alpha,lambda,Pcgl,sat,dsat)
  G = dsat_dpc*dsat_inv
  dS_dpl = dsat*(-1.d0)
  if (G == 1.d0) then
    dsi_dpl = 0.d0
    dsl_dpl = (1.d0 - s_i)*dS_dpl
  else
    dsi_dpl = (1.d0 - s_i)/(G/(1.d0 - G) + sat)*dS_dpl
    dsl_dpl = dsi_dpl*G/(1.d0 - G)
  endif
  dsg_dpl = -dsi_dpl - dsl_dpl

  
  ! Calculate the derivatives of saturation with respect to temp
  L = dsat_dpc
  if (T >= 0.d0) then
    M = 0.d0
  else
    M = temp_term/T
  endif
  N = dsat_inv
  dsi_dT = -L*M/(L*N + (1.d0 - L*N)*sat)
  dsl_dT = -dsi_dT*sat
  dsg_dT = -dsi_dT - dsl_dT

#if 0     
  ! Numerical derivatives      
  
  call CalcPhasePartitionIceNewt(alpha,lambda,Pcgl*(1.d0 + delta),T,s_g_pinc, &
                                 s_i_pinc,s_l_pinc)
                                 
  if (T >= 0.d0) then
    dsi_dT_num = 0.d0
    dsg_dT_num = 0.d0
    dsl_dT_num = 0.d0 
  else                                 
    call CalcPhasePartitionIceNewt(alpha,lambda,Pcgl,T*(1.d0 + delta),s_g_Tinc, &
                                   s_i_Tinc,s_l_Tinc)
    dsi_dT_num = (s_i_Tinc - s_i)/(T*delta)
    dsg_dT_num = (s_g_Tinc - s_g)/(T*delta)
    dsl_dT_num = (s_l_Tinc - s_l)/(T*delta)
  endif
                                   
  dsi_dpl_num = (s_i_pinc - s_i)/(delta*Pcgl)*(-1.d0) ! -1.d0 factor for dPcgl/dpl
  dsg_dpl_num = (s_g_pinc - s_g)/(delta*Pcgl)*(-1.d0)
  dsl_dpl_num = (s_l_pinc - s_l)/(delta*Pcgl)*(-1.d0) 
  
  print *, 'analytical-press:', 'dsg_dpl:', dsg_dpl, &
           'dsi_dpl:', dsi_dpl, 'dsl_dpl:', dsl_dpl 
  print *, 'numerical-press:' , 'dsg_dpl:', dsg_dpl_num, &
           'dsi_dpl:', dsi_dpl_num, 'dsl_dpl:', dsl_dpl_num
  
  print *, 'analytical-temp:', 'dsg_dT:', dsg_dT, &
           'dsi_dT:', dsi_dT, 'dsl_dT:', dsl_dT
  print *, 'numerical-temp:', 'dsg_dT:', dsg_dT_num, &
           'dsi_dT:', dsi_dT_num, 'dsl_dT:', dsl_dT_num
 
  dsi_dpl = dsi_dpl_num
  dsg_dpl = dsg_dpl_num
  dsl_dpl = dsl_dpl_num
  
  dsi_dT = dsi_dT_num
  dsg_dT = dsg_dT_num
  dsl_dT = dsl_dT_num
#endif
   

end subroutine CalcPhasePartitionIceDeriv 

! ************************************************************************** !

subroutine SatFuncComputeIcePKImplicit(pl,T,s_i,s_l,s_g,kr,dsl_dpl, & 
                                     dsl_dT,dsg_dpl,dsg_dT,dsi_dpl, &
                                     dsi_dT,dkr_dpl,dkr_dT, &
                                     saturation_function,pth,option)
  ! 
  ! Calculates the saturations of water phases
  ! and their derivative with respect to liquid
  ! pressure
  
  ! Model based on Painter & Karra implicit model VZJ (2014)
  ! 
  ! Author: Satish Karra
  ! Date: 10/16/12
  ! 

  use Option_module
 
implicit none

  PetscReal :: pl, T
  PetscReal :: s_i, s_g, s_l, kr
  PetscReal :: dkr_dpl, dkr_dT
  PetscReal :: dkr_dsl
  PetscReal :: alpha, m, Pcgl
  PetscReal :: one_over_m
  PetscReal :: liq_sat_one_over_m
  PetscReal :: dsl_dpl, dsl_dT
  PetscReal :: dsi_dpl, dsi_dT
  PetscReal :: dsg_dpl, dsg_dT
  PetscReal :: pth
  
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option


  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      alpha = saturation_function%alpha
      Pcgl = option%reference_pressure - pl
      m = saturation_function%m      
      call CalcPhasePartitionIceDeriv(alpha,m,Pcgl,T,s_g,s_i,s_l,dsg_dpl, &
                                      dsg_dT,dsi_dpl,dsi_dT,dsl_dpl, &
                                      dsl_dT)
    case default  
      option%io_buffer = 'Only van Genuchten supported with ice'
      call printErrMsg(option)
  end select
 
  ! Check for bounds on saturations         
  if (s_l > 1.d0) then
    print *, option%myrank, 'vG Liquid Saturation > 1:', s_l
  else if (s_l < 0.d0) then
    print *, option%myrank, 'vG Liquid Saturation < 0:', s_l
  endif

  if (s_g > 1.d0) then
    print *, option%myrank, 'vG Gas Saturation > 1:', s_g
  else if (s_g < 0.d0) then
    print *, option%myrank, 'vG Gas Saturation < 0:', s_g
  endif
 
  if (s_i > 1.d0) then
    print *, option%myrank, 'vG Ice Saturation > 1:', s_i
  else if (s_i < 0.d0) then
    print *, option%myrank, 'vG Ice Saturation < 0:', s_i
  endif
 
  ! Calculate relative permeability
  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      if (s_l == 1.d0) then
        kr = 1.d0
        dkr_dsl = 0.d0
      else
        m = saturation_function%m
        one_over_m = 1.d0/m
        liq_sat_one_over_m = s_l**one_over_m
        kr = sqrt(s_l)*(1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
        dkr_dsl = 0.5d0*kr/s_l + &
                  2.d0*s_l**(one_over_m - 0.5d0)* &
                  (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                  (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
      endif
        dkr_dpl = dkr_dsl*dsl_dpl
        dkr_dT = dkr_dsl*dsl_dT
    case default
      option%io_buffer = 'Ice module only supports Mualem' 
      call printErrMsg(option)
  end select 

#if 0  
  write(*,*) 'rank:', option%myrank, 'sl:', s_l, &
  'sg:', s_g, 'si:', s_i, 'dsl_pl:', dsl_dpl, &
  'dsl_temp:', dsl_dT, 'dsg_pl:', dsg_dpl, 'dsg_temp:', dsg_dT, &
  'dsi_pl:', dsi_dpl, 'dsi_temp:', dsi_dT, 'kr:', kr, &
  'dkr_pl:', dkr_dpl, 'dkr_temp:', dkr_dT, 'pl:', pl, 'T:', T  
#endif

end subroutine SatFuncComputeIcePKImplicit

! ************************************************************************** !

subroutine SatFuncComputeIcePKExplicit(pl,T,s_i,s_l,s_g,kr,dsl_dpl, & 
                                       dsl_dT,dsg_dpl,dsg_dT,dsi_dpl, &
                                       dsi_dT,dkr_dpl,dkr_dT, &
                                       saturation_function,pth,option)
  ! 
  ! Calculates the saturations of water phases
  ! and their derivative with respect to liquid
  ! pressure, temperature
  !
  ! Model used: explicit model from Painter & Karra, VJZ (2014).
  ! 
  ! Author: Satish Karra
  ! Date: 01/13/13
  ! 

  use Option_module
 
implicit none

  PetscReal :: pl, T
  PetscReal :: s_i, s_g, s_l, kr
  PetscReal :: dkr_dpl, dkr_dT
  PetscReal :: dkr_dsl
  PetscReal :: alpha, m, Pcgl
  PetscReal :: one_over_m
  PetscReal :: liq_sat_one_over_m
  PetscReal :: dsl_dpl, dsl_dT
  PetscReal :: dsi_dpl, dsi_dT
  PetscReal :: dsg_dpl, dsg_dT
  PetscReal :: pth
  PetscReal :: S, dS, Sinv, dSinv
  PetscReal, parameter :: beta = 2.2           ! dimensionless -- ratio of surf. tension
  PetscReal, parameter :: rho_l = 9.998d2      ! in kg/m^3
  PetscReal, parameter :: T_0 = 273.15         ! in K
  PetscReal, parameter :: L_f = 3.34d5         ! in J/kg
  PetscReal :: T_f, theta, X, Y, dS_dX

  
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option


  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      T = T + T_0  ! convert to K 
      alpha = saturation_function%alpha
      Pcgl = option%reference_pressure - pl
      m = saturation_function%m      
      call ComputeVGAndDerivative(alpha,m,Pcgl,S,dS)
      call ComputeInvVGAndDerivative(alpha,m,S,Sinv,dSinv)
      T_f = T_0 - 1.d0/beta*T_0/L_f/rho_l*Sinv
      theta = (T-T_0)/T_0
      if (T < T_f) then
        X = -beta*theta*L_f*rho_l
        call ComputeVGAndDerivative(alpha,m,X,s_l,dS_dX)
        s_i = 1.d0 - s_l/S
        dsl_dT = 1.d0/T_0*dS_dX*(-beta*rho_l*L_f)
        dsl_dpl = 0.d0
        dsi_dT = -1.d0/S*dsl_dT
        dsi_dpl = -1.d0*s_l/S**2*dS 
      else
        s_l = S
        s_i = 0.d0
        dsl_dpl = -dS
        dsl_dT = 0.d0
        dsi_dpl = 0.d0
        dsi_dT = 0.d0
      endif
      s_g = 1.d0 - s_l - s_i
      dsg_dpl = -dsl_dpl - dsi_dpl
      dsg_dT = -dsl_dT - dsi_dT  
      T = T - T_0 ! change back to C 
    case default  
      option%io_buffer = 'Only van Genuchten supported with ice'
      call printErrMsg(option)
  end select
 
  ! Check for bounds on saturations         
  if (s_l > 1.d0) then
    print *, option%myrank, 'vG Liquid Saturation > 1:', s_l
  else if (s_l < 0.d0) then
    print *, option%myrank, 'vG Liquid Saturation < 0:', s_l
  endif

  if (s_g > 1.d0) then
    print *, option%myrank, 'vG Gas Saturation > 1:', s_g
  else if (s_g < 0.d0) then
    print *, option%myrank, 'vG Gas Saturation < 0:', s_g
  endif
 
  if (s_i > 1.d0) then
    print *, option%myrank, 'vG Ice Saturation > 1:', s_i
  else if (s_i < 0.d0) then
    print *, option%myrank, 'vG Ice Saturation < 0:', s_i
  endif
 
  ! Calculate relative permeability
  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      if (s_l == 1.d0) then
        kr = 1.d0
        dkr_dsl = 0.d0
      else
        m = saturation_function%m
        one_over_m = 1.d0/m
        liq_sat_one_over_m = s_l**one_over_m
        kr = sqrt(s_l)*(1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
        dkr_dsl = 0.5d0*kr/s_l + &
                  2.d0*s_l**(one_over_m - 0.5d0)* &
                  (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                  (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
      endif
        dkr_dpl = dkr_dsl*dsl_dpl
        dkr_dT = dkr_dsl*dsl_dT
    case default
      option%io_buffer = 'Ice module only supports Mualem' 
      call printErrMsg(option)
  end select 

#if 0
  print *, '======================================' 
  write(*,*) 'rank:', option%myrank, 'sl:', s_l, &
  'sg:', s_g, 'si:', s_i, 'dsl_pl:', dsl_dpl, &
  'dsl_temp:', dsl_dT, 'dsg_pl:', dsg_dpl, 'dsg_temp:', dsg_dT, &
  'dsi_pl:', dsi_dpl, 'dsi_temp:', dsi_dT, 'kr:', kr, &
  'dkr_pl:', dkr_dpl, 'dkr_temp:', dkr_dT, 'pl:', pl, 'T:', T  
#endif

end subroutine SatFuncComputeIcePKExplicit

! ************************************************************************** !

subroutine SatFuncComputeIcePKExplicitNoCryo(pl,T,s_i,s_l,s_g,kr,dsl_dpl, & 
                                             dsl_dT,dsg_dpl,dsg_dT,dsi_dpl, &
                                             dsi_dT,dkr_dpl,dkr_dT, &
                                             saturation_function,pth,option)
  ! 
  ! Calculates the saturations of water phases
  ! and their derivative with respect to liquid
  ! pressure, temperature
  !
  ! Model used: explicit model from Painter & Karra, VJZ (2014).
  ! Removed Cryosuction. For this, we assume that s_i + s_l = S*(pcgl), that is
  ! the total water content does not change
  ! 
  ! Author: Satish Karra
  ! Date: 01/13/13
  ! 

  use Option_module
 
implicit none

  PetscReal :: pl, T
  PetscReal :: s_i, s_g, s_l, kr
  PetscReal :: dkr_dpl, dkr_dT
  PetscReal :: dkr_dsl
  PetscReal :: alpha, m, Pcgl
  PetscReal :: one_over_m
  PetscReal :: liq_sat_one_over_m
  PetscReal :: dsl_dpl, dsl_dT
  PetscReal :: dsi_dpl, dsi_dT
  PetscReal :: dsg_dpl, dsg_dT
  PetscReal :: pth
  PetscReal :: S, dS, Sinv, dSinv
  PetscReal, parameter :: beta = 2.2           ! dimensionless -- ratio of surf. tension
  PetscReal, parameter :: rho_l = 9.998d2      ! in kg/m^3
  PetscReal, parameter :: T_0 = 273.15         ! in K
  PetscReal, parameter :: L_f = 3.34d5         ! in J/kg
  PetscReal :: T_f, theta, X, Y, dS_dX

  
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option


  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      T = T + T_0  ! convert to K 
      alpha = saturation_function%alpha
      Pcgl = option%reference_pressure - pl
      m = saturation_function%m      
      call ComputeVGAndDerivative(alpha,m,Pcgl,S,dS)
      call ComputeInvVGAndDerivative(alpha,m,S,Sinv,dSinv)
      T_f = T_0 - 1.d0/beta*T_0/L_f/rho_l*Sinv
      theta = (T-T_0)/T_0
      if (T < T_f) then
        X = -beta*theta*L_f*rho_l
        call ComputeVGAndDerivative(alpha,m,X,s_l,dS_dX)
        s_i = S - s_l
        dsl_dT = 1.d0/T_0*dS_dX*(-beta*rho_l*L_f)
        dsl_dpl = 0.d0
        dsi_dT = -dsl_dT
        dsi_dpl = -dS - dsl_dpl 
      else
        s_l = S
        s_i = 0.d0
        dsl_dpl = -dS
        dsl_dT = 0.d0
        dsi_dpl = 0.d0
        dsi_dT = 0.d0
      endif
      s_g = 1.d0 - s_l - s_i
      dsg_dpl = -dsl_dpl - dsi_dpl
      dsg_dT = -dsl_dT - dsi_dT  
      T = T - T_0 ! change back to C 
    case default  
      option%io_buffer = 'Only van Genuchten supported with ice'
      call printErrMsg(option)
  end select
 
  ! Check for bounds on saturations         
  if (s_l > 1.d0) then
    print *, option%myrank, 'vG Liquid Saturation > 1:', s_l
  else if (s_l < 0.d0) then
    print *, option%myrank, 'vG Liquid Saturation < 0:', s_l
  endif

  if (s_g > 1.d0) then
    print *, option%myrank, 'vG Gas Saturation > 1:', s_g
  else if (s_g < 0.d0) then
    print *, option%myrank, 'vG Gas Saturation < 0:', s_g
  endif
 
  if (s_i > 1.d0) then
    print *, option%myrank, 'vG Ice Saturation > 1:', s_i
  else if (s_i < 0.d0) then
    print *, option%myrank, 'vG Ice Saturation < 0:', s_i
  endif
 
  ! Calculate relative permeability
  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      if (s_l == 1.d0) then
        kr = 1.d0
        dkr_dsl = 0.d0
      else
        m = saturation_function%m
        one_over_m = 1.d0/m
        liq_sat_one_over_m = s_l**one_over_m
        kr = sqrt(s_l)*(1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
        dkr_dsl = 0.5d0*kr/s_l + &
                  2.d0*s_l**(one_over_m - 0.5d0)* &
                  (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                  (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
      endif
        dkr_dpl = dkr_dsl*dsl_dpl
        dkr_dT = dkr_dsl*dsl_dT
    case default
      option%io_buffer = 'Ice module only supports Mualem' 
      call printErrMsg(option)
  end select 

#if 0
  print *, '======================================' 
  write(*,*) 'rank:', option%myrank, 'sl:', s_l, &
  'sg:', s_g, 'si:', s_i, 'dsl_pl:', dsl_dpl, &
  'dsl_temp:', dsl_dT, 'dsg_pl:', dsg_dpl, 'dsg_temp:', dsg_dT, &
  'dsi_pl:', dsi_dpl, 'dsi_temp:', dsi_dT, 'kr:', kr, &
  'dkr_pl:', dkr_dpl, 'dkr_temp:', dkr_dT, 'pl:', pl, 'T:', T  
#endif

end subroutine SatFuncComputeIcePKExplicitNoCryo

! ************************************************************************** !

subroutine SatFuncComputeIceDallAmico(pl, T, &
                                      p_fh2o, &
                                      dp_fh2o_dP, &
                                      dp_fh2o_dT, &
                                      s_i, s_l, s_g, &
                                      kr, &
                                      dsl_dpl, dsl_dT, &
                                      dsg_dpl, dsg_dT, &
                                      dsi_dpl, dsi_dT, &
                                      dkr_dpl, dkr_dT, &
                                      saturation_function, &
                                      option)
  !
  ! Calculates the saturations of water phases  and their derivative with
  ! respect to liquid pressure and temperature.
  !
  ! Model used: Dall'Amico (2010) and Dall' Amico et al. (2011)
  !
  ! Author: Gautam Bisht
  ! Date: 02/24/14
  !

  use Option_module

  implicit none

  PetscReal :: pl, T
  PetscReal :: p_fh2o, dp_fh2o_dP, dp_fh2o_dT
  PetscReal :: s_i, s_g, s_l
  PetscReal :: kr
  PetscReal :: dsl_dpl, dsl_dT
  PetscReal :: dsi_dpl, dsi_dT
  PetscReal :: dsg_dpl, dsg_dT
  PetscReal :: dkr_dpl, dkr_dT
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscReal :: Se,Sr
  PetscReal :: dkr_dsl, dkr_dSe
  PetscReal :: alpha
  PetscReal :: m
  PetscReal :: Pc0, Pc1
  PetscReal :: one_over_m
  PetscReal :: liq_sat_one_over_m
  PetscReal :: S0,S1
  PetscReal :: dS0,dS1
  PetscReal :: H, dH_dT
  PetscReal :: T_star
  PetscReal :: theta
  PetscReal :: x
  PetscReal :: dummy
  PetscBool :: switch
  PetscReal :: numer
  PetscReal :: denom
  PetscReal :: fct
  PetscReal :: T_star_th
  PetscReal :: T_star_min
  PetscReal :: T_star_max

  !PetscReal, parameter :: beta = 2.2           ! dimensionless -- ratio of surf. tension
  PetscReal, parameter :: beta = 1             ! dimensionless [assumed as 1.d0]
  PetscReal, parameter :: rho_l = 9.998d2      ! in kg/m^3
  PetscReal, parameter :: T_0 = 273.15         ! in K
  PetscReal, parameter :: L_f = 3.34d5         ! in J/kg
  PetscReal, parameter :: k = 1.d6

  T_star_th = 5.d-1 ! [K]

  s_g = 0.d0
  dsg_dpl = 0.d0
  dsg_dT = 0.d0

  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)

      T = T + T_0  ! convert to K
      alpha = saturation_function%alpha
      m = saturation_function%m

      Pc0 = option%reference_pressure - pl

      T_star = T_0 - 1.d0/beta*T_0/L_f/rho_l*Pc0

      if (T<T_star) then
        H = 1.d0
        dH_dT = 0.d0
      else
        H = 0.d0
        dH_dT = 0.d0
      endif

      !GB: Add an option to swich between step-function and smooth approximation
      !    of step function.
      !x = (T - T_star)*k
      !H = 0.5d0 - atan(x)/PI
      T_star_min = T_star - T_star_th
      T_star_max = T_star

      if (T < T_star_min) then
        H = 1.d0
        dH_dT = 0.d0
      else if (T > T_star_max) then
        H = 0.d0
        dH_dT = 0.d0
      else
        numer = T          - T_star_min
        denom = T_star_max - T_star_min
        fct   = 1.d0 - (numer/denom)**2.d0
        H     = fct**2.d0
        dH_dT = -4.d0*numer*fct/(denom*denom)
      endif

      theta = (T - T_star)/T_star
      Pc1 = Pc0 - beta*theta*L_f*rho_l*H

      p_fh2o = option%reference_pressure - Pc1
      dp_fh2o_dT = -beta*L_f*rho_l*H/T_star
      dp_fh2o_dP = 1.d0 - T*T_0/T_star/T_star*H
      p_fh2o = pl
      dp_fh2o_dT = 0.d0
      dp_fh2o_dP = 1.d0

      ! dummy and switch are not used here
      call SaturationFunctionCompute2(Pc0,S0,dummy,dS0,dummy,saturation_function,dummy,dummy,switch,option)
      call SaturationFunctionCompute2(Pc1,S1,dummy,dS1,dummy,saturation_function,dummy,dummy,switch,option)

      ! convert dS/dpsi to dS/dpl
      dS0 = -dS0
      dS1 = -dS1

      s_l = S1
      s_i = S0 - S1

      dsl_dpl = -dS1*(1.0d0 - T*T_0/T_star/T_star*H)
      dsi_dpl = -dS0 - dsl_dpl

      dsl_dT = dS1*(-beta*L_f*rho_l*H/T_star - beta*theta*L_f*rho_l*dH_dT)
      dsi_dT = -dsl_dT

      T = T - T_0 ! change back to C

    case default
      option%io_buffer = 'Only van Genuchten supported with ice'
      call printErrMsg(option)
  end select

  ! Calculate relative permeability
  select case(saturation_function%permeability_function_itype)
    case(MUALEM)
      Sr = saturation_function%Sr(1)
      Se = (s_l-Sr)/(1.0d0-Sr)
      if ( abs(Se-1.d0) < 1.0d-12 ) then
        kr = 1.d0
        dkr_dsl = 0.d0
      else
        m = saturation_function%m
        one_over_m = 1.d0/m
        liq_sat_one_over_m = Se**one_over_m
        kr = sqrt(Se)*(1.d0 - (1.d0 - liq_sat_one_over_m)**m)**2.d0
        dkr_dSe = 0.5d0*kr/Se + &
                  2.d0*Se**(one_over_m - 0.5d0)* &
                 (1.d0 - liq_sat_one_over_m)**(m - 1.d0)* &
                 (1.d0 - (1.d0 - liq_sat_one_over_m)**m)
        dkr_dsl = dkr_dSe / ( 1.0d0 - Sr )
      endif
      dkr_dpl = dkr_dsl*dsl_dpl
      dkr_dT = dkr_dsl*dsl_dT
    case default
      option%io_buffer = 'Ice module only supports Mualem'
      call printErrMsg(option)
  end select

end subroutine SatFuncComputeIceDallAmico

! ************************************************************************** !

subroutine SatFuncGetLiqRelPermFromSat(saturation,relative_perm,dkr_dSe, &
                                       saturation_function,iphase, &
                                       derivative,option)
  ! 
  ! Calculates relative permeability from
  ! phase saturation
  !                                    
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/05/11
  ! 

  use Option_module
  
  implicit none

  PetscReal :: saturation, relative_perm, dkr_dSe
  PetscReal :: power, pct_over_pcmax, pc_over_pcmax, pc_log_ratio
  PetscReal :: pcmax, one_over_alpha, alpha, lambda
  PetscInt :: iphase
  type(saturation_function_type) :: saturation_function
  PetscBool :: derivative
  type(option_type) :: option

  PetscReal :: m, Sr
  PetscReal :: Se, one_over_m, Se_one_over_m
  
  Sr = saturation_function%Sr(iphase)
  Se = (saturation-Sr)/(1.d0-Sr)
  
  ! if saturation is below residual saturation (this is possible with 
  ! two-phase with dry-out), need to bail out.
  if (Se <= 0.d0) then
    relative_perm = 0.d0
    dkr_dSe = 0.d0
    return
  else if (Se > 1.d0) then
    Se = 1.d0
  endif
    
  ! compute relative permeability
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
    ! compute relative permeability
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          m = saturation_function%m
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = Se*Se*(1.d0-(1.d0-Se_one_over_m)**m)
          if (derivative) then
            dkr_dSe = 2.d0*relative_perm/Se + &
                 Se*Se_one_over_m*(1.d0-Se_one_over_m)**(m-1.d0)
          endif
        case(MUALEM)
          ! reference #1
          m = saturation_function%m
          one_over_m = 1.d0/m
          Se_one_over_m = Se**one_over_m
          relative_perm = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
          if (derivative) then
            dkr_dSe = 0.5d0*relative_perm/Se+ &
                 2.d0*Se**(one_over_m-0.5d0)* &
                      (1.d0-Se_one_over_m)**(m-1.d0)* &
                      (1.d0-(1.d0-Se_one_over_m)**m)
          endif
        case default
          option%io_buffer = 'Unknown relative permeabilty function' 
          call printErrMsg(option)
      end select
    case(BROOKS_COREY)
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          lambda = saturation_function%lambda
          power = 3.d0+2.d0/lambda
          relative_perm = Se**power
          dkr_dSe = power*relative_perm/Se
        case(MUALEM)
          ! reference #1
          lambda = saturation_function%lambda
          power = 2.5d0+2.d0/lambda
          relative_perm = Se**power
          dkr_dSe = power*relative_perm/Se
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(LINEAR_MODEL)
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          relative_perm = Se
          if (derivative) then
            dkr_dSe = 1.d0
          endif
        case(MUALEM)
          power = 5.d-1
          alpha = saturation_function%alpha
          one_over_alpha = 1.d0/alpha
          pcmax = saturation_function%pcwmax
          
          pct_over_pcmax = one_over_alpha/pcmax
          pc_over_pcmax = 1.d0-(1.d0-pct_over_pcmax)*Se
          pc_log_ratio = log(pc_over_pcmax)/log(pct_over_pcmax)
          relative_perm = (Se**power)*(pc_log_ratio**2.d0)
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(LEVERETT)
      select case(saturation_function%permeability_function_itype)
        case(FATT_KLIKOFF)
          relative_perm = Se**3
          if (derivative) then
            dkr_dSe = 3.d0 * Se**2
          endif
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
  end select
  
end subroutine SatFuncGetLiqRelPermFromSat

! ************************************************************************** !

subroutine SatFuncGetGasRelPermFromSat(liquid_saturation, &
                                       gas_relative_perm, &
                                       saturation_function,option)
  ! 
  ! Calculates gas phase relative permeability from liquid saturation
  !
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/05/11
  ! 

  use Option_module
  
  implicit none

  PetscReal :: liquid_saturation
  PetscReal :: gas_relative_perm
  PetscReal :: power, pct_over_pcmax, pc_over_pcmax, pc_log_ratio
  PetscReal :: pcmax, one_over_alpha, alpha, liq_relative_perm
  type(saturation_function_type) :: saturation_function
  PetscBool :: derivative
  type(option_type) :: option

  PetscReal :: Srl, Srg
  PetscReal :: S_star, S_hat
  PetscReal :: Sg
  PetscReal :: lambda
  PetscReal :: m
  
  Srl = saturation_function%Sr(LIQUID_PHASE)
  Srg = saturation_function%Sr(GAS_PHASE)
  S_star = (liquid_saturation-Srl)/(1.d0-Srl)
  S_hat = (liquid_saturation-Srl)/(1.d0-Srl-Srg)

  gas_relative_perm = 0.d0

  if (S_hat >= 1.d0) then
    return
  else if (S_hat <= 0.d0) then
    gas_relative_perm = 1.d0
    return
  endif
  Sg = 1.d0-S_hat
  
  ! compute relative permeability
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
    ! compute relative permeability
      m = saturation_function%m
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          gas_relative_perm = Sg*Sg*(1.d0-S_hat**(1.d0/m))**m
        case(MUALEM)
          ! reference #1
          gas_relative_perm = sqrt(Sg)*(1.d0-S_hat**(1.d0/m))**(2.d0*m)
        case default
          option%io_buffer = 'Unknown relative permeabilty function' 
          call printErrMsg(option)
      end select
    case(BROOKS_COREY)
      lambda = saturation_function%lambda
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          ! reference #1
          gas_relative_perm = Sg*Sg* &
                              (1.d0-S_hat**(1.d0+2.d0/lambda))
        case(MUALEM)
          ! reference #1
          gas_relative_perm = sqrt(Sg)* &
                              (1.d0-S_hat**(1.d0+1.d0/lambda))**2.d0
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(LINEAR_MODEL)
      option%io_buffer = &
        'Linear model not yet supported in SatFuncGetGasRelPermFromSat.'
      call printErrMsg(option)
      select case(saturation_function%permeability_function_itype)
        case(BURDINE)
          gas_relative_perm = Sg  
        case(MUALEM)
          power = 5.d-1
          alpha = saturation_function%alpha
          one_over_alpha = 1.d0/alpha
          pcmax = saturation_function%pcwmax
          
          pct_over_pcmax = one_over_alpha/pcmax
          pc_over_pcmax = 1.d0-(1.d0-pct_over_pcmax)*S_star
          pc_log_ratio = log(pc_over_pcmax)/log(pct_over_pcmax)
          liq_relative_perm = (S_star**power)*(pc_log_ratio**2.d0)
          
          gas_relative_perm = Sg**power * liq_relative_perm * S_hat**(-power)
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
    case(LEVERETT)
      select case(saturation_function%permeability_function_itype)
        case(FATT_KLIKOFF)
          gas_relative_perm = Sg**3
        case default
          option%io_buffer = 'Unknown relative permeabilty function'
          call printErrMsg(option)
      end select
  end select

end subroutine SatFuncGetGasRelPermFromSat

! ************************************************************************** !

subroutine CapillaryPressureThreshold(saturation_function,cap_threshold,option)
  ! 
  ! Computes the capillary pressure threshold
  ! after which instead of van Genuchten a linear function is used upto 100 Mpa
  ! capillary pressure. The saturation at 100 Mpa is set to zero
  ! This threshold value depends only on van Genuchten parameters alpha and lambda
  ! This is used mainly for ice problem, so that the pressure doesnt go to large
  ! negative values
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/12/12
  ! 

  use Option_module
  
  implicit none
  
  PetscReal :: alpha, lambda, cap_threshold
  type(option_type) :: option
  type(saturation_function_type) :: saturation_function

  
  PetscReal :: gamma, p_new, res_value, jac_value, p_old
  PetscReal, parameter :: eps = 1.d-8
  PetscInt, parameter :: maxit = 100
  PetscInt :: i 
  
  alpha = saturation_function%alpha
  lambda = saturation_function%m
  alpha = alpha*1.d6
  gamma = 1.d0/(1.d0 - lambda)
  
  p_old = 99.d0
  
  
  do i = 1, maxit
    call ResidualCapPressThre(p_old,alpha,lambda,gamma,res_value)
    call JacobianCapPressThre(p_old,alpha,lambda,gamma,jac_value)
    p_new = p_old - res_value/jac_value
    !write (*,*) 'rank:', option%myrank, 'iter:', i, 'p_new:', p_new, 'p_old:', p_old, &
    !  'residual:', res_value, 'jacobian:', jac_value
    if ((abs(p_new - p_old) < eps)) exit
    p_old = p_new
  enddo
  
  cap_threshold = p_old*1.d6 ! convert to Pa
  
end subroutine CapillaryPressureThreshold

! ************************************************************************** !

subroutine ResidualCapPressThre(p,alpha,lambda,gamma,res)
  ! 
  ! Computes the residual to calculate capillary pressure
  ! thresold in the subroutine CapillaryPressureThreshold
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/12/12
  ! 

  implicit none
  
  PetscReal :: p, alpha, lambda, gamma, res
  
  res = lambda*gamma*alpha**gamma*p**(gamma-1.0)*1.d2 - &
        (alpha*p)**gamma*(1.d0 + gamma*lambda) - 1.d0

end subroutine ResidualCapPressThre

! ************************************************************************** !

subroutine JacobianCapPressThre(p,alpha,lambda,gamma,jac)
  ! 
  ! Computes the jacobian to calculate capillary pressure
  ! thresold in the subroutine CapillaryPressureThreshold
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 09/12/12
  ! 

  implicit none
  
  PetscReal :: p, alpha, lambda, gamma, jac
  
  jac = lambda*gamma*alpha**gamma*(gamma - 1.d0)*p**(gamma - 2.d0)*1.d2 - &
        alpha**gamma*gamma*p**(gamma - 1.d0)*(1.d0 + gamma*lambda)
  
 
end subroutine JacobianCapPressThre

! ************************************************************************** !

subroutine SatFuncGetCapillaryPressure(capillary_pressure,saturation, &
                                       temp,saturation_function,option)
  ! 
  ! Computes the capillary pressure as a function of
  ! pressure
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/03/09
  ! 

  use Option_module
  use Utility_module, only : QuadraticPolynomialEvaluate
  
  implicit none

  PetscReal :: capillary_pressure, saturation, temp
  type(saturation_function_type) :: saturation_function
  type(option_type) :: option

  PetscInt :: iphase
  PetscReal :: alpha, lambda, m, n, Sr, one_over_alpha
  PetscReal :: Se, derivative
  PetscReal :: pc_alpha, pc_alpha_n, one_plus_pc_alpha_n
  PetscReal :: pc_alpha_neg_lambda, pcmax
  PetscReal :: f, sigma, tk, os

  iphase = 1

  Sr = saturation_function%Sr(iphase)
  if (saturation <= Sr) then
    capillary_pressure = saturation_function%pcwmax
    return
  else if (saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
    
  select case(saturation_function%saturation_function_itype)
    case(VAN_GENUCHTEN)
      alpha = saturation_function%alpha
      m = saturation_function%m
      n = 1.d0/(1.d0-m)
      Se = (saturation-Sr)/(1.d0-Sr)
      one_plus_pc_alpha_n = Se**(-1.d0/m)
      pc_alpha_n = one_plus_pc_alpha_n - 1.d0
      pc_alpha = pc_alpha_n**(1.d0/n)
      capillary_pressure = pc_alpha/alpha
    case(BROOKS_COREY)
      alpha = saturation_function%alpha
      lambda = saturation_function%lambda
      Sr = saturation_function%Sr(iphase)
      Se = (saturation-Sr)/(1.d0-Sr)
      if (Se > saturation_function%sat_spline_low) then
        call QuadraticPolynomialEvaluate(saturation_function% &
                                         sat_spline_coefficients(1:3), &
                                         Se,capillary_pressure,derivative)
      else
        pc_alpha_neg_lambda = Se
        capillary_pressure = (pc_alpha_neg_lambda**(-1.d0/lambda))/alpha
      endif
    case(LEVERETT)
      tk = 273.15d0 + temp
      Se = (saturation-Sr)/(1.d0-Sr)
      os = 1.d0-Se
      f = os*(1.417d0 + os*(-2.120d0 + 1.263d0*os))
      sigma = 1.d0 - 0.625d0 * (374.15d0 - tk)/H2O_CRITICAL_TEMPERATURE
      sigma = sigma * 0.2358d0 * &
              ((374.15d0 - tk)/H2O_CRITICAL_TEMPERATURE)**1.256d0
      capillary_pressure = 632455.53d0 * sigma * f
    case(LINEAR_MODEL)
      alpha = saturation_function%alpha
      one_over_alpha = 1.d0/alpha
      pcmax = saturation_function%pcwmax
      Sr = saturation_function%Sr(iphase)
      Se = (saturation-Sr)/(1.d0-Sr)
      capillary_pressure = (one_over_alpha-pcmax)*Se + pcmax
#if 0
    case(THOMEER_COREY)
      pc = option%reference_pressure-pressure
      por = auxvar1
      perm = auxvar2*1.013202d15 ! convert from m^2 to mD
      Fg = saturation_function%alpha
      a = saturation_function%m
      Pd = 100.d0*por/sqrt(perm/(3.8068d0*(Fg**(-1.334d0)))) ! psi
      PHg = 9.63051d-4*pc
      if (PHg > Pd) then
        saturation = 1.d0-exp(-Fg/log10(PHg/Pd))
#if 0
        alpha = pc*(1.d0+1.d-8)
        m = 9.63051d-4*alpha
        n = 1.d0-exp(-Fg/log10(m/Pd))
        n = (n-saturation)/(alpha-pc)
#endif        
        dsat_dpc = (saturation-1.d0)*Fg/(log10(PHg/Pd)**2.d0)/(pc*2.30258509d0)
        ! Sr assumed to be zero
        relative_perm = saturation**a
        dkr_dpc = a*saturation**(a-1.d0)*dsat_dpc
      else
        saturation = 1.d0
        relative_perm = 1.d0
        return
      endif
#endif
    case default
      option%io_buffer = 'Unknown saturation function'
      call printErrMsg(option)
  end select

  capillary_pressure = min(capillary_pressure,saturation_function%pcwmax)

end subroutine SatFuncGetCapillaryPressure

! ************************************************************************** !

function SaturationFunctionGetID(saturation_function_list, &
                                 saturation_function_name, &
                                 material_property_name, option)
  ! 
  ! Returns the ID of the saturation function named
  ! "saturation_function_name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  type(saturation_function_type), pointer :: saturation_function_list
  character(len=MAXWORDLENGTH) :: saturation_function_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: SaturationFunctionGetID
  PetscBool :: found
  type(saturation_function_type), pointer :: cur_saturation_function

  found = PETSC_FALSE
  cur_saturation_function => saturation_function_list
  do 
    if (.not.associated(cur_saturation_function)) exit
    if (StringCompare(saturation_function_name, &
                      cur_saturation_function%name,MAXWORDLENGTH)) then
      found = PETSC_TRUE
      SaturationFunctionGetID = cur_saturation_function%id
      return
    endif
    cur_saturation_function => cur_saturation_function%next
  enddo
  if (.not.found) then
    option%io_buffer = 'Saturation function "' // &
             trim(saturation_function_name) // &
             '" in material property "' // &
             trim(material_property_name) // &
             '" not found among available saturation functions.'
    call printErrMsg(option)    
  endif

end function SaturationFunctionGetID

! ************************************************************************** !

subroutine SaturationFunctionVerify(saturation_function,option)
  ! 
  ! Evaluates saturation function curves for plotting external to PFLOTRAN
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/28/14
  ! 

  use Option_module
  
  implicit none

  type(saturation_function_type) :: saturation_function
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: pc, pc_increment, pc_max
  PetscReal :: sat, dummy_real
  PetscInt :: count, i
  PetscReal :: x(101), y(101), krl(101), krg(101)

  if (.not.(saturation_function%saturation_function_itype == VAN_GENUCHTEN .or.&
            saturation_function%saturation_function_itype == BROOKS_COREY)) then
    return
  endif

  ! calculate saturation as a function of capillary pressure
  ! start at 1 Pa up to maximum capillary pressure
  pc_max = saturation_function%pcwmax
  pc = 1.d0
  pc_increment = 1.d0
  count = 0
  do
    if (pc > pc_max) exit
    count = count + 1
    call SaturationFunctionCompute(pc,sat,saturation_function,option)
    x(count) = pc
    y(count) = sat
    if (pc > 0.99*pc_increment*10.d0) pc_increment = pc_increment*10.d0
    pc = pc + pc_increment
  enddo
  
  write(string,*) saturation_function%name
  string = trim(saturation_function%name) // '_pc_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation"' 
  do i = 1, count
    write(86,'(2es14.6)') x(i), y(i)
  enddo
  close(86)
  
  ! calculate capillary pressure as a function of saturation
  do i = 1, 101
    sat = dble(i-1)*0.01d0
    call SatFuncGetCapillaryPressure(pc,sat,option%reference_temperature, &
                                     saturation_function,option)
    x(i) = sat
    y(i) = pc
    call SatFuncGetLiqRelPermFromSat(sat,krl(i),dummy_real, &
                                     saturation_function,ONE_INTEGER, &
                                     PETSC_FALSE,option)
    call SatFuncGetGasRelPermFromSat(sat,krg(i),saturation_function,option)
  enddo  
  count = 101
  
  write(string,*) saturation_function%name
  string = trim(saturation_function%name) // '_sat_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure"' 
  do i = 1, count
    write(86,'(2es14.6)') x(i), y(i)
  enddo
  close(86)
  
  write(string,*) saturation_function%name
  string = trim(saturation_function%name) // '_krl.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "liquid relative permeability"' 
  do i = 1, count
    write(86,'(2es14.6)') x(i), krl(i)
  enddo
  close(86)
  
  write(string,*) saturation_function%name
  string = trim(saturation_function%name) // '_krg.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "gas relative permeability"' 
  do i = 1, count
    write(86,'(2es14.6)') x(i), krg(i)
  enddo
  close(86)
  
end subroutine SaturationFunctionVerify  
  
! ************************************************************************** !

recursive subroutine SaturationFunctionDestroy(saturation_function)
  ! 
  ! Destroys a saturation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  type(saturation_function_type), pointer :: saturation_function
  
  if (.not.associated(saturation_function)) return
  
  call SaturationFunctionDestroy(saturation_function%next)
    
  if (associated(saturation_function%Sr)) deallocate(saturation_function%Sr)
  nullify(saturation_function%Sr)

  if (associated(saturation_function%Kr0)) deallocate(saturation_function%Kr0)
  nullify(saturation_function%Kr0)
    
  deallocate(saturation_function)
  nullify(saturation_function)
  
end subroutine SaturationFunctionDestroy

end module Saturation_Function_module
