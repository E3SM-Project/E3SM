module Fluid_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys 
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: fluid_property_type
    PetscReal :: tort_bin_diff
    PetscReal :: vap_air_diff_coef
    PetscReal :: exp_binary_diff
    PetscReal :: enh_binary_diff_coef
    PetscReal :: diff_base
    PetscReal :: diff_exp
    character(len=MAXWORDLENGTH) :: phase_name
    PetscInt :: phase_id
    PetscReal :: diffusion_coefficient
    PetscReal :: gas_diffusion_coefficient
    PetscReal :: diffusion_activation_energy
    PetscReal :: nacl_concentration
    type(fluid_property_type), pointer :: next
  end type fluid_property_type
  
  public :: FluidPropertyCreate, FluidPropertyDestroy, &
            FluidPropertyRead, FluidPropertyAddToList

contains

! ************************************************************************** !

function FluidPropertyCreate()
  ! 
  ! Creates a fluid property object
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 
  
  implicit none

  type(fluid_property_type), pointer :: FluidPropertyCreate
  
  type(fluid_property_type), pointer :: fluid_property
  
  allocate(fluid_property)
  fluid_property%tort_bin_diff = 0.d0
  fluid_property%vap_air_diff_coef = 2.13d-5
  fluid_property%exp_binary_diff = 0.d0
  fluid_property%enh_binary_diff_coef = 0.d0
  fluid_property%diff_base = 0.d0
  fluid_property%diff_exp = 0.d0
  fluid_property%phase_name = 'LIQUID'
  fluid_property%phase_id = 0
  fluid_property%diffusion_coefficient = 1.d-9
  fluid_property%gas_diffusion_coefficient = 2.13D-5
  ! for liquid, one can use 12.6 kJ/mol as an activation energy
  fluid_property%diffusion_activation_energy = 0.d0
  fluid_property%nacl_concentration = 0.d0
  nullify(fluid_property%next)
  FluidPropertyCreate => fluid_property

end function FluidPropertyCreate

! ************************************************************************** !

subroutine FluidPropertyRead(fluid_property,input,option)
  ! 
  ! Reads in contents of a fluid property card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Units_module

  implicit none
  
  type(fluid_property_type) :: fluid_property
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXWORDLENGTH) :: internal_units

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','FLUID_PROPERTY')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('PHASE') 
        call InputReadWord(input,option,fluid_property%phase_name,PETSC_TRUE)
        call InputErrorMsg(input,option,'phase','FLUID_PROPERTY')
      case('DIFFUSION_COEFFICIENT','LIQUID_DIFFUSION_COEFFICIENT') 
        call InputReadDouble(input,option,fluid_property%diffusion_coefficient)
        call InputErrorMsg(input,option,'diffusion coefficient', &
                           'FLUID_PROPERTY')
        call InputReadAndConvertUnits(input, &
                                      fluid_property%diffusion_coefficient, &
                         'm^2/sec','FLUID_PROPERTY,diffusion_coeffient',option)
      case('DIFFUSION_ACTIVATION_ENERGY') 
        call InputReadDouble(input,option, &
                             fluid_property%diffusion_activation_energy)
        call InputErrorMsg(input,option,'diffusion activation energy', &
                           'FLUID_PROPERTY')
        call InputReadAndConvertUnits(input, &
                fluid_property%diffusion_activation_energy, &
                'J/mol','FLUID_PROPERTY,diffusion activation energy',option)
      case('GAS_DIFFUSION_COEFFICIENT') 
        call InputReadDouble(input,option, &
                             fluid_property%gas_diffusion_coefficient)
        call InputErrorMsg(input,option,'gas diffusion coefficient', &
                           'FLUID_PROPERTY')
      case default
        call InputKeywordUnrecognized(keyword,'FLUID_PROPERTY',option)
    end select
    
  enddo  

  if (.not.(StringCompareIgnoreCase(fluid_property%phase_name,'LIQUID') .or. &
            StringCompareIgnoreCase(fluid_property%phase_name,'GAS'))) then
    option%io_buffer = 'PHASE in FLUID_PROPERTY should be LIQUID or GAS.'
    call printErrMsg(option)
  endif


end subroutine FluidPropertyRead

! ************************************************************************** !

subroutine FluidPropertyAddToList(fluid_property,list)
  ! 
  ! Adds a thermal property to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  implicit none
  
  type(fluid_property_type), pointer :: fluid_property
  type(fluid_property_type), pointer :: list

  type(fluid_property_type), pointer :: cur_fluid_property
  
  if (associated(list)) then
    cur_fluid_property => list
    ! loop to end of list
    do
      if (.not.associated(cur_fluid_property%next)) exit
      cur_fluid_property => cur_fluid_property%next
    enddo
    cur_fluid_property%next => fluid_property
  else
    list => fluid_property
  endif
  
end subroutine FluidPropertyAddToList

! ************************************************************************** !

recursive subroutine FluidPropertyDestroy(fluid_property)
  ! 
  ! Destroys a fluid property
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  implicit none
  
  type(fluid_property_type), pointer :: fluid_property
  
  if (.not.associated(fluid_property)) return
  
  call FluidPropertyDestroy(fluid_property%next)

  deallocate(fluid_property)
  nullify(fluid_property)
  
end subroutine FluidPropertyDestroy

end module Fluid_module
