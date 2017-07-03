module SrcSink_Sandbox_WIPP_Well_class

! Sandbox srcsink for WIPP well source terms
  
  use PFLOTRAN_Constants_module
  use SrcSink_Sandbox_Base_class
  
  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: WIPP_WELL_LIQUID_MOBILITY = 1
  PetscInt, parameter, public :: WIPP_WELL_GAS_MOBILITY = 2
  PetscInt, parameter, public :: WIPP_WELL_LIQUID_PRESSURE = 3
  PetscInt, parameter, public :: WIPP_WELL_GAS_PRESSURE = 4
  PetscInt, parameter, public :: WIPP_WELL_LIQUID_ENTHALPY = 5
  PetscInt, parameter, public :: WIPP_WELL_GAS_ENTHALPY = 6
  PetscInt, parameter, public :: WIPP_WELL_XMOL_AIR_IN_LIQUID = 7
  PetscInt, parameter, public :: WIPP_WELL_XMOL_WATER_IN_GAS = 8
  PetscInt, parameter, public :: WIPP_WELL_LIQUID_DENSITY = 9
  PetscInt, parameter, public :: WIPP_WELL_GAS_DENSITY = 10

  type, public, &
    extends(srcsink_sandbox_base_type) :: srcsink_sandbox_wipp_well_type
    PetscReal :: well_pressure      ! Pa
    PetscReal :: productivity_index ! m^3
  contains
    procedure, public :: ReadInput => WIPPWellRead
    procedure, public :: Setup => WIPPWellSetup
    procedure, public :: Evaluate => WIPPWellSrcSink
    procedure, public :: Destroy => WIPPWellDestroy
  end type srcsink_sandbox_wipp_well_type

  public :: WIPPWellCreate

contains

! ************************************************************************** !

function WIPPWellCreate()
  ! 
  ! Allocates WIPP well src/sink object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_well_type), pointer :: WIPPWellCreate

  allocate(WIPPWellCreate)
  call SSSandboxBaseInit(WIPPWellCreate)  
  nullify(WIPPWellCreate%next)  
      
end function WIPPWellCreate

! ************************************************************************** !

subroutine WIPPWellRead(this,input,option)
  ! 
  ! Reads input deck for WIPP well src/sink parameters
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(srcsink_sandbox_wipp_well_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units
  PetscBool :: found
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'SRCSINK_SANDBOX,WIPP')
    call StringToUpper(word)   

    call SSSandboxBaseSelectCase(this,input,option,word,found)
    if (found) cycle
    
    select case(trim(word))
      case('WELL_PRESSURE')
        call InputReadDouble(input,option,this%well_pressure)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,WIPP,WELL')
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr == 0) then
          internal_units = 'Pa'
          this%well_pressure = this%well_pressure * &
               UnitsConvertToInternal(word,internal_units,option)
        endif
      case('WELL_PRODUCTIVITY_INDEX')
        call InputReadDouble(input,option,this%productivity_index)
        call InputErrorMsg(input,option,word,'SOURCE_SINK_SANDBOX,WIPP,WELL')
      case default
        call InputKeywordUnrecognized(word,'SRCSINK_SANDBOX,WIPP',option)
    end select
  enddo

end subroutine WIPPWellRead

! ************************************************************************** !

subroutine WIPPWellSetup(this,grid,option)
    
  use Option_module
  use Grid_module
  
  implicit none
  
  class(srcsink_sandbox_wipp_well_type) :: this
  type(grid_type) :: grid
  type(option_type) :: option
    
  call SSSandboxBaseSetup(this,grid,option)

end subroutine WIPPWellSetup 

! ************************************************************************** !

subroutine WIPPWellSrcSink(this,Residual,Jacobian,compute_derivative, &
                           material_auxvar,aux_real,option)
  ! 
  ! Evaluates src/sink storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_wipp_well_type) :: this  
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
  PetscReal :: q_liquid, q_gas
  
  ! q_i = rho_i*I*(kr_i/mu_i)*(p_i-p_well)
  ! units: m^3/s
  
  q_liquid = -1.d0 * &
            this%productivity_index * &
            aux_real(WIPP_WELL_LIQUID_MOBILITY) * &
            (aux_real(WIPP_WELL_LIQUID_PRESSURE) - this%well_pressure)
  q_gas = -1.d0 * &
          this%productivity_index * &
          aux_real(WIPP_WELL_GAS_MOBILITY) * &
          (aux_real(WIPP_WELL_GAS_PRESSURE) - this%well_pressure)
  ! positive is inflow
  ! units = kmol/s
  ! water equation
  Residual(ONE_INTEGER) = q_liquid * aux_real(WIPP_WELL_LIQUID_DENSITY) * &
                          (1.d0-aux_real(WIPP_WELL_XMOL_AIR_IN_LIQUID)) + &
                          q_gas * aux_real(WIPP_WELL_XMOL_WATER_IN_GAS) * &
                          aux_real(WIPP_WELL_GAS_DENSITY)
  ! air equation
  Residual(TWO_INTEGER) = q_liquid * aux_real(WIPP_WELL_XMOL_AIR_IN_LIQUID) * & 
                          aux_real(WIPP_WELL_LIQUID_DENSITY) + &
                          q_gas * aux_real(WIPP_WELL_GAS_DENSITY) * &
                          (1.d0-aux_real(WIPP_WELL_XMOL_WATER_IN_GAS))
  ! energy equation
  ! units = MJ/s
  Residual(THREE_INTEGER) = q_liquid * aux_real(WIPP_WELL_LIQUID_DENSITY) * &
                            aux_real(WIPP_WELL_LIQUID_ENTHALPY) + &
                            q_gas * aux_real(WIPP_WELL_GAS_DENSITY) * &
                            aux_real(WIPP_WELL_GAS_ENTHALPY)
  
  if (associated(this%instantaneous_mass_rate)) then
    this%instantaneous_mass_rate(:) = -1.d0*Residual(:)
  endif
                            
  if (compute_derivative) then
    
    ! jacobian something

  endif
  
end subroutine WIPPWellSrcSink

! ************************************************************************** !

subroutine WIPPWellDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/14
  ! 

  implicit none
  
  class(srcsink_sandbox_wipp_well_type) :: this
  
  call SSSandboxBaseDestroy(this)

end subroutine WIPPWellDestroy

end module SrcSink_Sandbox_WIPP_Well_class
