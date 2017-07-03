module Geomechanics_Debug_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, public :: geomech_debug_type
    PetscBool :: vecview_residual
    PetscBool :: vecview_solution
    PetscBool :: matview_Jacobian
    PetscBool :: matview_Jacobian_detailed
    PetscBool :: norm_Jacobian
    PetscBool :: print_numerical_derivatives
    PetscBool :: print_couplers
    character(len=MAXSTRINGLENGTH) :: coupler_string
    PetscBool :: print_waypoints
  end type geomech_debug_type
   
  public :: GeomechDebugCreate, GeomechDebugRead
  
contains

! ************************************************************************** !

function GeomechDebugCreate()
  ! 
  ! Create object that stores debugging options
  ! for geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(geomech_debug_type), pointer :: GeomechDebugCreate
  type(geomech_debug_type), pointer :: debug
  
  allocate(debug)
  
  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Jacobian = PETSC_FALSE
  debug%matview_Jacobian_detailed = PETSC_FALSE
  debug%norm_Jacobian = PETSC_FALSE
  debug%print_numerical_derivatives = PETSC_FALSE
  debug%print_couplers = PETSC_FALSE
  debug%coupler_string = ''
  debug%print_waypoints = PETSC_FALSE

  GeomechDebugCreate => debug

end function GeomechDebugCreate

! ************************************************************************** !

subroutine GeomechDebugRead(debug,input,option)
  ! 
  ! Reads debugging data from the input file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
    
  type(geomech_debug_type) :: debug
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS_DEBUG')   
      
    select case(trim(keyword))
    
      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = PETSC_TRUE
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = PETSC_TRUE
      case('PRINT_JACOBIAN','MATVIEW_JACOBIAN','VIEW_JACOBIAN')
        debug%matview_Jacobian = PETSC_TRUE
      case('PRINT_JACOBIAN_NORM','NORM_JACOBIAN')
        debug%norm_Jacobian = PETSC_TRUE
      case('PRINT_COUPLERS','PRINT_COUPLER')
        debug%print_couplers = PETSC_TRUE
        debug%coupler_string = trim(adjustl(input%buf))
      case('PRINT_JACOBIAN_DETAILED','MATVIEW_JACOBIAN_DETAILED', &
           'VIEW_JACOBIAN_DETAILED')
        debug%matview_Jacobian_detailed = PETSC_TRUE
      case('PRINT_NUMERICAL_DERIVATIVES','VIEW_NUMERICAL_DERIVATIVES')
        debug%print_numerical_derivatives = PETSC_TRUE
      case('WAYPOINTS')
        debug%print_waypoints = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'GEOMECHANICS_DEBUG',option)
    end select 
  
  enddo  

end subroutine GeomechDebugRead

end module Geomechanics_Debug_module
