module Integral_Flux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  

  PetscInt, parameter, public :: INTEGRATE_FLOW = 1
  PetscInt, parameter, public :: INTEGRATE_TRANSPORT = 2

  type, public :: integral_flux_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    type(point3d_type), pointer :: coordinates(:)
    PetscBool :: invert_direction
    PetscInt, pointer :: connections(:)
    PetscReal, pointer :: integral_value(:)
    type(integral_flux_type), pointer :: next
  end type integral_flux_type
  
  type, public :: integral_flux_list_type
    PetscInt :: num_integral_fluxes
    type(integral_flux_type), pointer :: first
    type(integral_flux_type), pointer :: last
    type(integral_flux_type), pointer :: array(:)
  end type integral_flux_list_type

  public :: IntegralFluxCreate, &
            IntegralFluxDestroy, &
            IntegralFluxRead, &
            IntegralFluxAddToList, &
            IntegralFluxInitList, &
            IntegralFluxDestroyList, &
            IntegralFluxGetPtrFromList, &
            IntegralFluxSizeStorage, &
            IntegralFluxUpdate, &
            IntegralFluxGetInstantaneous

contains

! ************************************************************************** !

function IntegralFluxCreate()
  ! 
  ! Create object that stores integral flux information
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxCreate
  
  type(integral_flux_type), pointer :: integral_flux
  
  allocate(integral_flux)
  
  integral_flux%name = ''
  integral_flux%id = 0
  integral_flux%invert_direction = PETSC_FALSE
  nullify(integral_flux%connections)
  nullify(integral_flux%coordinates)
  nullify(integral_flux%integral_value)
  nullify(integral_flux%next)
  
  IntegralFluxCreate => integral_flux

end function IntegralFluxCreate

! ************************************************************************** !

subroutine IntegralFluxRead(integral_flux,input,option)
  ! 
  ! Reads integral flux data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','integral_flux')   
      
    select case(trim(keyword))
      case('NAME')
        call InputReadWord(input,option,integral_flux%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','INTEGRAL_FLUX')    
      case('INVERT_DIRECTION')
        integral_flux%invert_direction = PETSC_TRUE
      case('COORDINATES')
        call GeometryReadCoordinates(input,option,integral_flux%name, &
                                     integral_flux%coordinates)
      case default
        call InputKeywordUnrecognized(keyword,'INTEGRAL_FLUX',option)
    end select 
  
  enddo  

  if (len_trim(integral_flux%name) < 1) then
    option%io_buffer = 'All INTEGRAL_FLUXes must have a name.'
    call printErrMsg(option)
  endif

end subroutine IntegralFluxRead

! ************************************************************************** !

subroutine IntegralFluxSizeStorage(integral_flux,option)
  ! 
  ! Sizes the arrays that store the integrated flux
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  type(option_type) :: option
  
  allocate(integral_flux%integral_value(option%nflowdof+option%ntrandof))
  integral_flux%integral_value = 0.d0

end subroutine IntegralFluxSizeStorage

! ************************************************************************** !

subroutine IntegralFluxUpdate(integral_flux_list,internal_fluxes, &
                              boundary_fluxes,iflag,option)
  ! 
  ! Updates the stored integrated value of each integral flux measurement
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_list_type) :: integral_flux_list
  PetscReal :: internal_fluxes(:,:)
  PetscReal :: boundary_fluxes(:,:)
  PetscInt :: iflag
  type(option_type) :: option
  
  type(integral_flux_type), pointer :: integral_flux
  PetscReal, allocatable :: sum_array(:)
  PetscInt :: offset
  PetscInt :: num_values
  PetscReal :: dt
 
  if (.not.associated(integral_flux_list%first)) return

  select case(iflag)
    case(INTEGRATE_FLOW)
      offset = 0
      num_values = option%nflowdof
      dt = option%flow_dt
    case(INTEGRATE_TRANSPORT)
      offset = option%nflowdof
      num_values = option%ntrandof
      dt = option%tran_dt
    case default
      offset = -1 ! to catch bugs
  end select
  
  allocate(sum_array(num_values))
  integral_flux => integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    sum_array = 0.d0
    call IntegralFluxGetInstantaneous(integral_flux, internal_fluxes, &
                                      boundary_fluxes,num_values, &
                                      sum_array,option)
    integral_flux%integral_value(offset+1:offset+num_values) = &
      integral_flux%integral_value(offset+1:offset+num_values) + &
      sum_array(1:num_values)*dt
    integral_flux => integral_flux%next
  enddo
  deallocate(sum_array)

end subroutine IntegralFluxUpdate


! ************************************************************************** !

subroutine IntegralFluxGetInstantaneous(integral_flux, internal_fluxes, &
                                        boundary_fluxes,num_values, &
                                        sum_array,option)
  ! 
  ! Returns the instantaneous mole flux for an integral flux object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use Option_module
  
  implicit none
  
  type(integral_flux_type) :: integral_flux
  PetscReal :: internal_fluxes(:,:)
  PetscReal :: boundary_fluxes(:,:)
  PetscInt :: num_values
  PetscReal :: sum_array(:)
  type(option_type) :: option
  
  PetscInt :: i
  PetscInt :: iconn
  
  sum_array = 0.d0
  if (.not.associated(integral_flux%connections)) return
  
  do i = 1, size(integral_flux%connections)
    iconn = integral_flux%connections(i)
    if (iconn > 0) then ! internal
      sum_array(1:num_values) = sum_array(1:num_values) + &
                                internal_fluxes(1:num_values,iconn)
    else ! boundary
      sum_array(1:num_values) = sum_array(1:num_values) + &
                                boundary_fluxes(1:num_values,-iconn)
    endif
  enddo
  if (integral_flux%invert_direction) then
    sum_array = -1.d0 * sum_array
  endif

end subroutine IntegralFluxGetInstantaneous

! ************************************************************************** !

subroutine IntegralFluxInitList(list)
  ! 
  ! Initializes a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none

  type(integral_flux_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_integral_fluxes = 0

end subroutine IntegralFluxInitList

! ************************************************************************** !

subroutine IntegralFluxAddToList(new_integral_flux,list)
  ! 
  ! Adds a new integral_flux to a integral flux list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_type), pointer :: new_integral_flux
  type(integral_flux_list_type) :: list
  
  list%num_integral_fluxes = list%num_integral_fluxes + 1
  new_integral_flux%id = list%num_integral_fluxes
  if (.not.associated(list%first)) list%first => new_integral_flux
  if (associated(list%last)) list%last%next => new_integral_flux
  list%last => new_integral_flux
  
end subroutine IntegralFluxAddToList

! ************************************************************************** !

function IntegralFluxGetPtrFromList(integral_flux_name,integral_flux_list)
  ! 
  ! Returns a pointer to the integral flux matching integral_flux_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  use String_module

  implicit none
  
  type(integral_flux_type), pointer :: IntegralFluxGetPtrFromList
  character(len=MAXWORDLENGTH) :: integral_flux_name
  type(integral_flux_list_type) :: integral_flux_list
 
  PetscInt :: length
  type(integral_flux_type), pointer :: integral_flux
    
  nullify(IntegralFluxGetPtrFromList)
  integral_flux => integral_flux_list%first
  
  do 
    if (.not.associated(integral_flux)) exit
    length = len_trim(integral_flux_name)
    if (length == len_trim(integral_flux%name) .and. &
        StringCompare(integral_flux%name,integral_flux_name, &
                        length)) then
      IntegralFluxGetPtrFromList => integral_flux
      return
    endif
    integral_flux => integral_flux%next
  enddo
  
end function IntegralFluxGetPtrFromList

! ************************************************************************** !

subroutine IntegralFluxDestroyList(integral_flux_list)
  ! 
  ! Deallocates a list of integral fluxes
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 

  implicit none
  
  type(integral_flux_list_type), pointer :: integral_flux_list
  
  type(integral_flux_type), pointer :: integral_flux, prev_integral_flux
  
  if (.not.associated(integral_flux_list)) return
  
  integral_flux => integral_flux_list%first
  do 
    if (.not.associated(integral_flux)) exit
    prev_integral_flux => integral_flux
    integral_flux => integral_flux%next
    call IntegralFluxDestroy(prev_integral_flux)
  enddo
  
  integral_flux_list%num_integral_fluxes = 0
  nullify(integral_flux_list%first)
  nullify(integral_flux_list%last)
  if (associated(integral_flux_list%array)) deallocate(integral_flux_list%array)
  nullify(integral_flux_list%array)
  
  deallocate(integral_flux_list)
  nullify(integral_flux_list)

end subroutine IntegralFluxDestroyList

! ************************************************************************** !

subroutine IntegralFluxDestroy(integral_flux)
  ! 
  ! Deallocates a integral flux
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/20/14
  ! 
  use Utility_module
  
  implicit none
  
  type(integral_flux_type), pointer :: integral_flux
  
  PetscInt :: i
  
  if (.not.associated(integral_flux)) return
  
  if (associated(integral_flux%coordinates)) &
    deallocate(integral_flux%coordinates)
  nullify(integral_flux%coordinates)
  call DeallocateArray(integral_flux%connections)
  call DeallocateArray(integral_flux%integral_value)
  deallocate(integral_flux)
  nullify(integral_flux)

end subroutine IntegralFluxDestroy

end module Integral_Flux_module
