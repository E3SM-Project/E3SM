module PM_Surface_Flow_class

  use PM_Base_class
  use PM_Surface_class
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscts.h"

  type, public, extends(pm_surface_type) :: pm_surface_flow_type
  contains
    procedure, public :: Read => PMSurfaceFlowRead
    procedure, public :: UpdateTimestep => PMSurfaceFlowUpdateTimestep
    procedure, public :: PreSolve => PMSurfaceFlowPreSolve
    procedure, public :: PostSolve => PMSurfaceFlowPostSolve
    procedure, public :: UpdateSolution => PMSurfaceFlowUpdateSolution
    procedure, public :: Destroy => PMSurfaceFlowDestroy
    procedure, public :: RHSFunction => PMSurfaceFlowRHSFunction
    procedure, public :: UpdateAuxVars => PMSurfaceFlowUpdateAuxVars
    procedure, public :: InputRecord => PMSurfaceFlowInputRecord
  end type pm_surface_flow_type

  public :: PMSurfaceFlowCreate, &
            PMSurfaceFlowDTExplicit

contains

! ************************************************************************** !

function PMSurfaceFlowCreate()
  ! 
  ! Creates Surface-flow process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type), pointer :: PMSurfaceFlowCreate

  class(pm_surface_flow_type), pointer :: surface_flow_pm

  allocate(surface_flow_pm)
  call PMSurfaceCreate(surface_flow_pm)
  surface_flow_pm%name = 'Surface Flow'

  PMSurfaceFlowCreate => surface_flow_pm

end function PMSurfaceFlowCreate

! ************************************************************************** !

subroutine PMSurfaceFlowRead(this,input)
  ! 
  ! Reads input file parameters associated with the Surface process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module
  use Option_module

  implicit none

  class(pm_surface_flow_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option

  error_string = 'Surface Flow Options'

  input%ierr = 0
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMSurfaceReadSelectCase(this,input,word,found,option)
    if (found) cycle

    select case(trim(word))
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo

end subroutine PMSurfaceFlowRead

! ************************************************************************** !

subroutine PMSurfaceFlowPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE FLOW ",64("="))')
  endif

end subroutine PMSurfaceFlowPreSolve

! ************************************************************************** !

subroutine PMSurfaceFlowUpdateSolution(this)
  !
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowUpdateSolution
  use Condition_module

  implicit none

  class(pm_surface_flow_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

  call PMSurfaceUpdateSolution(this)
  call SurfaceFlowUpdateSolution(this%surf_realization)

end subroutine PMSurfaceFlowUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceFlowRHSFunction(this,ts,time,xx,ff,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowRHSFunction

  implicit none

  class(pm_surface_flow_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr

  call SurfaceFlowRHSFunction(ts,time,xx,ff,this%surf_realization,ierr)

end subroutine PMSurfaceFlowRHSFunction

! ************************************************************************** !

subroutine PMSurfaceFlowUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowComputeMaxDt

  implicit none

  class(pm_surface_flow_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceFlowComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceFlowUpdateTimestep

! ************************************************************************** !

subroutine PMSurfaceFlowDTExplicit(this,dt_max)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Surface_Flow_module, only : SurfaceFlowComputeMaxDt

  implicit none

  class(pm_surface_flow_type) :: this
  PetscReal :: dt_max

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceFlowComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt_max = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceFlowDTExplicit

! ************************************************************************** !

subroutine PMSurfaceFlowPostSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  use Discretization_module
  use Surface_Field_module
  use Surface_Flow_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pm_surface_flow_type) :: this

  PetscReal, pointer :: xx_p(:)
  PetscInt :: local_id
  type(surface_field_type), pointer :: surf_field 
  PetscErrorCode :: ierr

  surf_field => this%surf_realization%surf_field

  ! Ensure evolved solution is +ve
  call VecGetArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  do local_id = 1,this%surf_realization%discretization%grid%nlmax
    if (xx_p(local_id)<1.d-15) xx_p(local_id) = 0.d0
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Update aux vars
  call SurfaceFlowUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceFlowPostSolve

! ************************************************************************** !
subroutine PMSurfaceFlowUpdateAuxVars(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/30/14

  use Surface_Flow_module

  implicit none

  class(pm_surface_flow_type) :: this

  call SurfaceFlowUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceFlowUpdateAuxVars

! ************************************************************************** !

subroutine PMSurfaceFlowInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_surface_flow_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMSurfaceFlowInputRecord

! ************************************************************************** !

subroutine PMSurfaceFlowDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/11/13
  ! 

  implicit none

  class(pm_surface_flow_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceFlowDestroy()')
#endif

#ifndef SIMPLIFY
!  call SurfaceFlowDestroy(this%surf_realization)
#endif
  call this%comm1%Destroy()

end subroutine PMSurfaceFlowDestroy

end module PM_Surface_Flow_class
