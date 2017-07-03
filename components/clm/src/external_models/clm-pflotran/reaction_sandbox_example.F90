module Reaction_Sandbox_Example_class

! 1. Change all references to "Example" as desired to rename the module and
!    and subroutines within the module. 

  use Reaction_Sandbox_Base_class
  
  use Global_Aux_module
  use Reactive_Transport_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

! 2. Add module variables here.  Note that one must use the PETSc data types 
!    PetscInt, PetscReal, PetscBool to declare variables of type integer
!    float/real*8, and logical respectively.  E.g.
!
! PetscReal, parameter :: formula_weight_of_water = 18.01534d0
  
  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_example_type
! 3. Add variables/arrays associated with new reaction.  
    character(len=MAXWORDLENGTH) :: species_name
    PetscInt :: species_id
    PetscReal :: rate_constant
  contains
    procedure, public :: ReadInput => ExampleRead
    procedure, public :: Setup => ExampleSetup
    procedure, public :: Evaluate => ExampleReact
    procedure, public :: Destroy => ExampleDestroy
  end type reaction_sandbox_example_type

  public :: ExampleCreate

contains

! ************************************************************************** !

function ExampleCreate()
  ! 
  ! Allocates example reaction object.
  ! 
  ! Author: John Doe (replace in all subroutine headers with name of developer)
  ! Date: 00/00/00 (replace in all subroutine headers with current date)
  ! 

  implicit none
  
  class(reaction_sandbox_example_type), pointer :: ExampleCreate

! 4. Add code to allocate object and initialized all variables to zero and
!    nullify all pointers. E.g.
  allocate(ExampleCreate)
  ExampleCreate%species_name = ''
  ExampleCreate%species_id = 0
  ExampleCreate%rate_constant = 0.d0
  nullify(ExampleCreate%next)  
      
end function ExampleCreate

! ************************************************************************** !

subroutine ExampleRead(this,input,option)
  ! 
  ! Reads input deck for example reaction parameters (if any)
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal
  
  implicit none
  
  class(reaction_sandbox_example_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units
  
  do 
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,TEMPLATE')
    call StringToUpper(word)   

    select case(trim(word))

      ! Example Input:

      ! CHEMISTRY
      !   ...
      !   REACTION_SANDBOX
      !   : begin user-defined input
      !     TEMPLATE
      !       EXAMPLE_INTEGER 1
      !       EXAMPLE_INTEGER_ARRAY 2 3 4
      !     END
      !   : end user defined input
      !   END
      !   ...
      ! END

! 5. Add case statement for reading variables.
      case('SPECIES_NAME')
! 6. Read the variable
        ! Read the character string indicating which of the primary species
        ! is being decayed.
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)  
! 7. Inform the user of any errors if not read correctly.
        call InputErrorMsg(input,option,'species_name', &
                           'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')    
! 8. Repeat for other variables
      case('RATE_CONSTANT')
        ! Read the double precision rate constant
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,'rate_constant', &
                           'CHEMISTRY,REACTION_SANDBOX,EXAMPLE')
        ! Read the units
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (InputError(input)) then
          ! If units do not exist, assume default units of 1/s which are the
          ! standard internal PFLOTRAN units for this rate constant.
          input%err_buf = 'REACTION_SANDBOX,EXAMPLE,RATE CONSTANT UNITS'
          call InputDefaultMsg(input,option)
        else              
          ! If units exist, convert to internal units of 1/s
          internal_units = 'unitless/sec'
          this%rate_constant = this%rate_constant * &
            UnitsConvertToInternal(word,internal_units,option)
        endif
      case default
        call InputKeywordUnrecognized(word, &
                     'CHEMISTRY,REACTION_SANDBOX,TEMPLATE',option)
    end select
  enddo
  
end subroutine ExampleRead

! ************************************************************************** !

subroutine ExampleSetup(this,reaction,option)
  ! 
  ! Sets up the example reaction either with parameters either
  ! read from the input deck or hardwired.
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Reaction_Aux_module, only : reaction_type, GetPrimarySpeciesIDFromName
  use Option_module

  implicit none
  
  class(reaction_sandbox_example_type) :: this
  type(reaction_type) :: reaction
  type(option_type) :: option

! 9. Add code to initialize 
  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)
      
end subroutine ExampleSetup

! ************************************************************************** !

subroutine ExampleReact(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar,reaction, &
                         option)
  ! 
  ! Evaluates reaction storing residual and/or Jacobian
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_class
  
  implicit none
  
  class(reaction_sandbox_example_type) :: this  
  type(option_type) :: option
  type(reaction_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water
  
  ! Description of subroutine arguments:

  ! Residual - 1D array storing residual entries in units mol/sec
  ! Jacobian - 2D array storing Jacobian entires in units kg water/sec
  !
  !  Jacobian [kg water/sec] * dc [mol/kg water] = -Res [mol/sec]
  !
  ! compute_derivative - Flag indicating whether analtical derivative should
  !   be calculated.  The user must provide either the analytical derivatives 
  !   or a numerical approximation unless always running with 
  !   NUMERICAL_JACOBIAN_RXN defined in input deck.  If the use of 
  !   NUMERICAL_JACOBIAN_RXN is assumed, the user should provide an error 
  !   message when compute_derivative is true.  E.g.
  !
  !   option%io_buffer = 'NUMERICAL_JACOBIAN_RXN must always be used ' // &
  !                      'due to assumptions in Example'
  !   call printErrMsg(option)
  !
  ! rt_auxvar - Object holding chemistry information (e.g. concentrations,
  !   activity coefficients, mineral volume fractions, etc.).  See
  !   reactive_transport_aux.F90.  
  !
  !   Useful variables:
  !     rt_auxvar%total(:,iphase) - total component concentrations 
  !                                 [mol/L water] for phase
  !     rt_auxvar%pri_molal(:) - free ion concentrations [mol/kg water]
  !     rt_auxvar%pri_act_coef(:) - activity coefficients for primary species
  !     rt_auxvar%aqueous%dtotal(:,iphase) - derivative of total component
  !                 concentration with respect to free ion [kg water/L water]
  !
  ! global_auxvar - Object holding information on flow (e.g. saturation,
  !   density, viscosity, temperature, etc)
  !
  !   Useful variables:
  !     global_auxvar%den(iphase) - fluid density [mol/m^3] 
  !     global_auxvar%den_kg(iphase) - fluid density [kg/m^3] 
  !     global_auxvar%sat(iphase) - saturation 
  !     global_auxvar%temp - temperature [C]
  !
  ! porosity - effective porosity of grid cell [m^3 pore/m^3 bulk]                     
  ! volume - volume of grid cell [m^3]
  ! reaction - Provides access to variable describing chemistry.  E.g.
  !   reaction%ncomp - # chemical degrees of freedom (mobile and immobile)
  !   reaction%naqcomp - # chemical degrees of freedom on water
  !   reaction%primary_species_names(:) - names of primary species
  !
  ! option - Provides handle for controlling simulation, catching and
  !          reporting errors.
  
  ! Unit of the residual must be in moles/second
  ! global_auxvar%sat(iphase) = saturation of cell
  ! 1.d3 converts m^3 water -> L water
  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3
  ! always subtract contribution from residual
  Residual(this%species_id) = Residual(this%species_id) - &
    this%rate_constant * &  ! 1/sec
    L_water * & ! L water
    ! rt_auxvar%total(this%species_id,iphase) = species total component 
    !   concentration
    rt_auxvar%total(this%species_id,iphase) ! mol/L water
    
  
  
  if (compute_derivative) then

! 11. If using an analytical Jacobian, add code for Jacobian evaluation

    ! always add contribution to Jacobian
    ! units = (mol/sec)*(kg water/mol) = kg water/sec
    Jacobian(this%species_id,this%species_id) = &
    Jacobian(this%species_id,this%species_id) + &
      this%rate_constant * & ! 1/sec
      L_water * & ! L water
                  ! kg water/L water
      ! rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) = 
      !   derivative of total component concentration with respect to the
      !   free ion concentration of the same species.
      rt_auxvar%aqueous%dtotal(this%species_id,this%species_id,iphase) 

  endif
  
end subroutine ExampleReact

! ************************************************************************** !

subroutine ExampleDestroy(this)
  ! 
  ! Destroys allocatable or pointer objects created in this
  ! module
  ! 
  ! Author: John Doe
  ! Date: 00/00/00
  ! 

  implicit none
  
  class(reaction_sandbox_example_type) :: this  

! 12. Add code to deallocate contents of the example object

end subroutine ExampleDestroy

end module Reaction_Sandbox_Example_class
