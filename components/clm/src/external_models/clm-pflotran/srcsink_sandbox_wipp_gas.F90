module SrcSink_Sandbox_WIPP_Gas_class

! Sandbox srcsink for WIPP gas generation source terms
  
  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: WIPP_GAS_WATER_SATURATION_INDEX = 1
  PetscInt, parameter, public :: WIPP_GAS_TEMPERATURE_INDEX = 2
  
  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_wipp_gas_type
    PetscReal :: inundated_corrosion_rate
    PetscReal :: inundated_degradation_rate
    PetscReal :: humid_corrosion_factor
    PetscReal :: humid_degradation_factor
    PetscReal :: h2_fe_ratio
    PetscReal :: h2_ch2o_ratio  
    PetscReal :: alpharxn
    PetscReal :: satwick
    PetscReal :: smin
  contains
    procedure, public :: ReadInput => WIPPGasGenerationRead
    procedure, public :: Setup => WIPPGasGenerationSetup
    procedure, public :: Evaluate => WIPPGasGenerationSrcSink
    procedure, public :: Destroy => WIPPGasGenerationDestroy
  end type srcsink_sandbox_wipp_gas_type

  public :: WIPPGasGenerationCreate

contains

! ************************************************************************** !

function WIPPGasGenerationCreate()
  ! 
  ! Allocates WIPP gas generation src/sink object.
  ! 
  ! Author: Glenn Hammond, Edit: Heeho Park, Jenn Frederick
  ! Date: 04/11/14, 05/15/14, 01/09/2017
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type), pointer :: WIPPGasGenerationCreate

  allocate(WIPPGasGenerationCreate)
  call SSSandboxBaseInit(WIPPGasGenerationCreate)
  WIPPGasGenerationCreate%inundated_corrosion_rate = 3.0d-8   ! [1/sec]
  WIPPGasGenerationCreate%inundated_degradation_rate = 1.5d-7 ! [1/sec]
  WIPPGasGenerationCreate%humid_corrosion_factor = 1.0d-3
  WIPPGasGenerationCreate%humid_degradation_factor = 0.2d0
  WIPPGasGenerationCreate%h2_fe_ratio = 1.3081d0
  WIPPGasGenerationCreate%h2_ch2o_ratio = 1.1100d0
  WIPPGasGenerationCreate%alpharxn = -1.d3
  WIPPGasGenerationCreate%satwick = 0.25d0
  WIPPGasGenerationCreate%smin = 1.5d-2
  nullify(WIPPGasGenerationCreate%next)  
      
end function WIPPGasGenerationCreate

! ************************************************************************** !

subroutine WIPPGasGenerationRead(this,input,option)
  ! 
  ! Reads input deck for WIPP gas generation src/sink parameters
  ! 
  ! Author: Glenn Hammond, Edit: Heeho Park, Jenn Frederick
  ! Date: 04/11/14, 05/15/14, 01/09/2017
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXSTRINGLENGTH) :: error_strg
  PetscBool :: found
  
  error_strg = 'SRCSINK_SANDBOX,WIPP-GAS_GENERATION'
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_strg)
    call StringToUpper(word)   

    call SSSandboxBaseSelectCase(this,input,option,word,found)
    if (found) cycle
    
    select case(trim(word))
      case('INUNDATED_CORROSION_RATE')
        internal_units = 'unitless/sec'
        call InputReadDouble(input,option,this%inundated_corrosion_rate)
        call InputErrorMsg(input,option,'inundated_corrosion_rate',error_strg)
      case('INUNDATED_DEGRADATION_RATE')
        internal_units = 'unitless/sec'
        call InputReadDouble(input,option,this%inundated_degradation_rate)
        call InputErrorMsg(input,option,'inundated_degradation_rate',error_strg)
      case('HUMID_CORROSION_FACTOR')
        call InputReadDouble(input,option,this%humid_corrosion_factor)
        call InputErrorMsg(input,option,'humid_corrosion_factor',error_strg)
      case('HUMID_DEGRADATION_FACTOR')
        call InputReadDouble(input,option,this%humid_degradation_factor)
        call InputErrorMsg(input,option,'humid_degradation_factor',error_strg)
      case('H2_FE_RATIO')
        call InputReadDouble(input,option,this%h2_fe_ratio)
        call InputErrorMsg(input,option,'h2_fe_ratio',error_strg)
      case('H2_CH2O_RATIO')
        call InputReadDouble(input,option,this%h2_ch2o_ratio)
        call InputErrorMsg(input,option,'h2_ch2o_ratio',error_strg)
      case('ALPHARXN')
        call InputReadDouble(input,option,this%alpharxn)
        call InputErrorMsg(input,option,'alpharxn',error_strg)
      case('SOCMIN')
        call InputReadDouble(input,option,this%smin)
        call InputErrorMsg(input,option,'socmin',error_strg)
      case('SATWICK')
        call InputReadDouble(input,option,this%satwick)
        call InputErrorMsg(input,option,'satwick',error_strg)
      case default
        call InputKeywordUnrecognized(word, &
          'SRCSINK_SANDBOX,WIPP-GAS_GENERATION',option)
    end select
  enddo

end subroutine WIPPGasGenerationRead

! ************************************************************************** !

subroutine WIPPGasGenerationSetup(this,grid,option)
  ! 
  ! Sets up the WIPP gas generation src/sink
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  use Option_module
  use Grid_module
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  type(grid_type) :: grid
  type(option_type) :: option
  
  call SSSandboxBaseSetup(this,grid,option)

end subroutine WIPPGasGenerationSetup 

! ************************************************************************** !

subroutine WIPPGasGenerationSrcSink(this,Residual,Jacobian, &
                                    compute_derivative, &
                                    material_auxvar,aux_real,option)
  ! 
  ! Evaluates src/sink storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond, Edited: Heeho Park, Jenn Frederick
  ! Date: 04/11/14, 05/15/14, 01/09/2017
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  use EOS_Gas_module
  
  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this  
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
  
  ! gas generation calculation variables
  PetscReal :: water_saturation
  PetscReal :: s_eff
  PetscReal :: humid_corrosion_rate
  PetscReal :: humid_degradation_rate
  PetscReal :: corrosion_gas_rate
  PetscReal :: degradation_gas_rate
  PetscReal :: gas_generation_rate
  
  ! enthalpy calculation variables
  PetscReal :: T        ! temperature [C]
  PetscReal :: dummy_P  ! pressure [Pa]
  PetscReal :: H        ! enthalpy [J/kmol]
  PetscReal :: dH_dT    ! derivative enthalpy wrt temperature
  PetscReal :: dH_dP    ! derivative enthalpy wrt pressure
  PetscReal :: U        ! internal energy [J/kmol]
  PetscReal :: dU_dT    ! deriv. internal energy wrt temperature
  PetscReal :: dU_dP    ! deriv. internal energy wrt pressure
  PetscErrorCode :: ierr
  
!  PetscReal :: test_value
  
  ! get water saturation from each grid cell
  water_saturation = aux_real(WIPP_GAS_WATER_SATURATION_INDEX)
  ! calculate the effective saturation which takes into account wicking sat,
  ! alpharxn, and s_min (see BRAGFLO v6.02 User Manual, section 4.13.1)
  s_eff = water_saturation - this%smin + this%satwick*(1.d0 - &
          exp(200.d0*this%alpharxn* &
          (max((water_saturation-this%smin),0.d0))**2.d0))
  ! calculate gas generation parameters using s_eff instead of water_saturation
  humid_corrosion_rate = this%humid_corrosion_factor * &
                         this%inundated_corrosion_rate
  humid_degradation_rate = this%humid_degradation_factor * & 
                           this%inundated_degradation_rate
  ! Eq. 143 BRAGFLO v6.02 User's Manual, section 4.13.3, pg 55
  corrosion_gas_rate = this%inundated_corrosion_rate * s_eff + &
                       humid_corrosion_rate * (1.d0 - s_eff)
  ! Eq. 148 BRAGFLO v6.02 User's Manual, section 4.13.4, pg 55
  degradation_gas_rate = this%inundated_degradation_rate * s_eff + &
                         humid_degradation_rate * (1.d0 - s_eff)
  ! gas generation unit: mol of H2 / (m^3 * s)
  gas_generation_rate = this%h2_fe_ratio * corrosion_gas_rate + &
                         this%h2_ch2o_ratio * degradation_gas_rate
  ! convert to kmol/s -> mol/(m^3 * s) * (volume of the cell) * 
  !                       (1 kmol/1000 mol)
  gas_generation_rate = gas_generation_rate * material_auxvar%volume * 1.d-3
  !  test_value = gas_generation_rate * 2.01588d0 / material_auxvar%volume

  ! positive is inflow 
  ! kmol/s
  Residual(TWO_INTEGER) = gas_generation_rate 
  
  T = aux_real(WIPP_GAS_TEMPERATURE_INDEX)
  call EOSGasEnergy(T,dummy_P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
  ! energy equation
  ! units = MJ/s -> enthalpy(J/kmol) * (MJ/1E+6J) * gas_generation_rate (kmol/s)
  Residual(THREE_INTEGER) = H * 1.d-6 * gas_generation_rate
  
  if (compute_derivative) then
    
    ! jacobian something

  endif
  
end subroutine WIPPGasGenerationSrcSink

! ************************************************************************** !

subroutine WIPPGasGenerationDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_gas_type) :: this
  
  call SSSandboxBaseDestroy(this)

end subroutine WIPPGasGenerationDestroy

end module SrcSink_Sandbox_WIPP_Gas_class
