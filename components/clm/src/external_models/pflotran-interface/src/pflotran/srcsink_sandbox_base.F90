module SrcSink_Sandbox_Base_class
  
  use PFLOTRAN_Constants_module
  use Geometry_module

  implicit none
  
  private
  
#include "petsc/finclude/petscsys.h"

  type, abstract, public :: srcsink_sandbox_base_type
    PetscInt :: local_cell_id
    PetscInt :: natural_cell_id
    type(point3d_type) :: coordinate    
    PetscReal, pointer :: instantaneous_mass_rate(:)
    PetscReal, pointer :: cumulative_mass(:)
    class(srcsink_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => SSSandboxBaseRead
    procedure, public :: Setup => SSSandboxBaseSetup
    procedure, public :: Update => SSSandboxBaseUpdate
    procedure, public :: Evaluate => SSSandboxBaseEvaluate
    procedure, public :: Destroy => SSSandboxBaseDestroy    
  end type srcsink_sandbox_base_type
  
  public :: SSSandboxBaseInit, &
            SSSandboxBaseSetup, &
            SSSandboxBaseRead, &
            SSSandboxBaseSelectCase, &
            SSSandboxBaseDestroy
  
contains

! ************************************************************************** !

subroutine SSSandboxBaseInit(this)
    
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
    
  this%coordinate%x = UNINITIALIZED_DOUBLE
  this%coordinate%y = UNINITIALIZED_DOUBLE
  this%coordinate%z = UNINITIALIZED_DOUBLE
  this%local_cell_id = UNINITIALIZED_INTEGER
  this%natural_cell_id = UNINITIALIZED_INTEGER
  nullify(this%instantaneous_mass_rate)
  nullify(this%cumulative_mass)
  nullify(this%next)
  
end subroutine SSSandboxBaseInit

! ************************************************************************** !

subroutine SSSandboxBaseSetup(this,grid,option)
    
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Grid_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscInt :: local_id, ghosted_id
  PetscInt :: i, iflag
  PetscErrorCode :: ierr

  local_id = 0
  if (Initialized(this%coordinate%x)) then
    call GridGetLocalIDFromCoordinate(grid,this%coordinate,option,local_id)
  else if (Initialized(this%natural_cell_id)) then
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2A(ghosted_id) == this%natural_cell_id) then
        local_id = grid%nG2L(ghosted_id)
        if (local_id > 0) exit
      endif
    enddo
  else
    option%io_buffer = 'Source/sink in SSSandbox not associate with the &
      &domain through either a CELL_ID or COORDINATE.'
    call printErrMsg(option)
  endif  

  ! check to ensure that only one grid cell is mapped
  iflag = 0
  if (local_id > 0) then
    iflag = 1
    this%local_cell_id = local_id
  endif
  call MPI_Allreduce(MPI_IN_PLACE,iflag,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr)
  if (iflag > 1) then
    option%io_buffer = 'More than one grid cell mapped in SSSandboxBaseSetup.'
    call printErrMsg(option)
  else if (iflag == 0) then
    option%io_buffer = 'No grid cells mapped in SSSandboxBaseSetup.'
    call printErrMsg(option)
  endif
  
end subroutine SSSandboxBaseSetup 

! ************************************************************************** !

subroutine SSSandboxBaseRead(this,input,option)
    
  use Option_module
  use Input_Aux_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
end subroutine SSSandboxBaseRead  

! ************************************************************************** !

subroutine SSSandboxBaseSelectCase(this,input,option,keyword,found)
    
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use Input_Aux_module
  use Geometry_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  
  character(len=MAXSTRINGLENGTH) :: error_string
  
  error_string = 'SOURCE_SINK_SANDBOX'
  
  found = PETSC_TRUE
  select case(trim(keyword))
    case('REGION')
      option%io_buffer = 'The REGION card has been deprecated in &
        &Source/Sink Sandbox.  Please switch to using a COORDINATE and &
        &defining one Src/Sink block for each coordinate.'
      call printErrMsg(option)
    case('COORDINATE')
      call GeometryReadCoordinate(input,option,this%coordinate,error_string)
    case('CELL_ID')
      call InputReadInt(input,option,this%natural_cell_id)
      call InputErrorMsg(input,option,'cell id',error_string)
    case default
      found = PETSC_FALSE
  end select   
  
end subroutine SSSandboxBaseSelectCase

! ************************************************************************** !

subroutine SSSandboxBaseUpdate(this,option)
    
  use Option_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(option_type) :: option
  
  if (associated(this%cumulative_mass)) then
    this%cumulative_mass(:) = this%cumulative_mass(:) + &
      option%flow_dt*this%instantaneous_mass_rate(:)
  endif
  
end subroutine SSSandboxBaseUpdate   

! ************************************************************************** !

subroutine SSSandboxBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
                                 material_auxvar,aux_real,option)
  
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  type(option_type) :: option
  PetscBool :: compute_derivative
  PetscReal :: Residual(option%nflowdof)
  PetscReal :: Jacobian(option%nflowdof,option%nflowdof)
  class(material_auxvar_type) :: material_auxvar
  PetscReal :: aux_real(:)
      
end subroutine SSSandboxBaseEvaluate

! ************************************************************************** !

subroutine SSSandboxBaseDestroy(this)

  use Utility_module
  
  implicit none
  
  class(srcsink_sandbox_base_type) :: this
  
  call DeallocateArray(this%instantaneous_mass_rate)
  call DeallocateArray(this%cumulative_mass)

end subroutine SSSandboxBaseDestroy  

end module SrcSink_Sandbox_Base_class
