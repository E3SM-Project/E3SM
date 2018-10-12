module Debug_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private
  

  type, public :: debug_type
    PetscBool :: vecview_residual
    PetscBool :: vecview_solution
    PetscBool :: matview_Jacobian
    PetscBool :: matview_Jacobian_detailed
    PetscBool :: norm_Jacobian

    PetscBool :: binary_format

    PetscBool :: print_numerical_derivatives

    PetscBool :: print_couplers
    PetscBool :: print_regions
    character(len=MAXSTRINGLENGTH) :: coupler_string
    PetscBool :: print_waypoints
  end type debug_type

  public :: DebugCreate, DebugRead, DebugCreateViewer, DebugDestroy
  
contains

! ************************************************************************** !

function DebugCreate()
  ! 
  ! Create object that stores debugging options for PFLOW
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 

  implicit none
  
  type(debug_type), pointer :: DebugCreate
  
  type(debug_type), pointer :: debug
  
  allocate(debug)
  
  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Jacobian = PETSC_FALSE
  debug%matview_Jacobian_detailed = PETSC_FALSE
  debug%norm_Jacobian = PETSC_FALSE

  debug%binary_format = PETSC_FALSE
  
  debug%print_numerical_derivatives = PETSC_FALSE
  
  debug%print_couplers = PETSC_FALSE
  debug%print_regions = PETSC_FALSE
  debug%coupler_string = ''
  debug%print_waypoints = PETSC_FALSE

  DebugCreate => debug

end function DebugCreate

! ************************************************************************** !

subroutine DebugRead(debug,input,option)
  ! 
  ! Reads debugging data from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 

  use Option_module
  use Input_Aux_module
  
  implicit none
    
  type(debug_type) :: debug
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','DEBUG')   
      
    select case(trim(keyword))
    
      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = PETSC_TRUE
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = PETSC_TRUE
      case('PRINT_JACOBIAN','MATVIEW_JACOBIAN','VIEW_JACOBIAN')
        debug%matview_Jacobian = PETSC_TRUE
      case('PRINT_JACOBIAN_NORM','NORM_JACOBIAN')
        debug%norm_Jacobian = PETSC_TRUE
      case('PRINT_REGIONS')
        debug%print_regions = PETSC_TRUE
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
      case('BINARY_FORMAT')
        debug%binary_format = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'DEBUG',option)
    end select 
  
  enddo  

end subroutine DebugRead

! ************************************************************************** !

subroutine DebugCreateViewer(debug,viewer_name_prefix,option,viewer)
  !
  ! Creates a PETSc viewer for saving PETSc vector or matrix in ASCII or
  ! binary format
  !
  ! Author: Gautam Bisht
  ! Date: 09/23/14
  !

  use Option_module
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none

  type(debug_type), pointer :: debug
  character(len=MAXSTRINGLENGTH), intent(in) :: viewer_name_prefix
  type(option_type) :: option
  PetscViewer, intent (inout) :: viewer

  character(len=MAXWORDLENGTH) :: viewer_name
  PetscErrorCode :: ierr

  if (debug%binary_format) then
    viewer_name = trim(adjustl(viewer_name_prefix)) // '.bin'
    call PetscViewerBinaryOpen(option%mycomm,viewer_name, &
                               FILE_MODE_WRITE,viewer,ierr);CHKERRQ(ierr)
  else
    viewer_name = trim(viewer_name_prefix) // '.out'
    call PetscViewerASCIIOpen(option%mycomm,viewer_name,viewer, &
                              ierr);CHKERRQ(ierr)
  endif


end subroutine DebugCreateViewer

! ************************************************************************** !

subroutine DebugDestroy(debug)
  ! 
  ! Deallocates memory associated with debug object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 
  implicit none
  
  type(debug_type), pointer :: debug

  if (.not.associated(debug)) return
  
  deallocate(debug)
  nullify(debug)
  
end subroutine DebugDestroy

end module Debug_module
