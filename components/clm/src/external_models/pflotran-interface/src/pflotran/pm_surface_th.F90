module PM_Surface_TH_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Surface_class
  use Realization_Surface_class
  use Realization_Subsurface_class
  use Communicator_Base_module
  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_surface_type) :: pm_surface_th_type
  contains
    procedure, public :: Read => PMSurfaceTHRead
    procedure, public :: UpdateTimestep => PMSurfaceTHUpdateTimestep
    procedure, public :: PreSolve => PMSurfaceTHPreSolve
    procedure, public :: PostSolve => PMSurfaceTHPostSolve
    procedure, public :: UpdateSolution => PMSurfaceTHUpdateSolution
    procedure, public :: Destroy => PMSurfaceTHDestroy
    procedure, public :: RHSFunction => PMSurfaceTHRHSFunction
    procedure, public :: UpdateAuxVars => PMSurfaceTHUpdateAuxVars
    procedure, public :: InputRecord => PMSurfaceTHInputRecord
  end type pm_surface_th_type

  public :: PMSurfaceTHCreate, &
            PMSurfaceTHDTExplicit

contains

! ************************************************************************** !

function PMSurfaceTHCreate()
  ! 
  ! Creates Surface TH process model
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  implicit none

  class(pm_surface_th_type), pointer :: PMSurfaceTHCreate

  class(pm_surface_th_type), pointer :: surface_th_pm

  allocate(surface_th_pm)
  call PMSurfaceCreate(surface_th_pm)
  surface_th_pm%name = 'Surface+TH Flow'

  PMSurfaceTHCreate => surface_th_pm

end function PMSurfaceTHCreate

! ************************************************************************** !

subroutine PMSurfaceTHRead(this,input)
  ! 
  ! Reads input file parameters associated with the Surface TH process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use Surface_TH_Aux_module
 
  implicit none
  
  class(pm_surface_th_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'Surface TH Options'
  
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
  
end subroutine PMSurfaceTHRead

! ************************************************************************** !

subroutine PMSurfaceTHPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  implicit none

  PetscErrorCode :: ierr
  class(pm_surface_th_type) :: this

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," SURFACE TH FLOW ",61("="))')
  endif

end subroutine PMSurfaceTHPreSolve

! ************************************************************************** !

subroutine PMSurfaceTHUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHUpdateSolution
  use Condition_module

  implicit none

  class(pm_surface_th_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

  call PMSurfaceUpdateSolution(this)
  call SurfaceTHUpdateSolution(this%surf_realization)

end subroutine PMSurfaceTHUpdateSolution

! ************************************************************************** !

subroutine PMSurfaceTHRHSFunction(this,ts,time,xx,ff,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHRHSFunction

  implicit none

  class(pm_surface_th_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr

  call SurfaceTHRHSFunction(ts,time,xx,ff,this%surf_realization,ierr)

end subroutine PMSurfaceTHRHSFunction

! ************************************************************************** !

subroutine PMSurfaceTHUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                     num_newton_iterations,tfac, &
                                     time_step_max_growth_factor)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHComputeMaxDt

  implicit none

  class(pm_surface_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt = min(0.9d0*dt_max_glb,this%surf_realization%dt_max, &
           time_step_max_growth_factor*dt)

end subroutine PMSurfaceTHUpdateTimestep

! ************************************************************************** !

subroutine PMSurfaceTHDTExplicit(this,dt_max)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHComputeMaxDt

  implicit none

  class(pm_surface_th_type) :: this
  PetscReal :: dt_max

  PetscReal :: dt_max_glb
  PetscErrorCode :: ierr
  PetscReal :: dt_max_loc

  call SurfaceTHComputeMaxDt(this%surf_realization,dt_max_loc)
  call MPI_Allreduce(dt_max_loc,dt_max_glb,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                     MPI_MIN,this%option%mycomm,ierr)
  dt_max = min(0.9d0*dt_max_glb,this%surf_realization%dt_max)

end subroutine PMSurfaceTHDTExplicit

! ************************************************************************** !

subroutine PMSurfaceTHPostSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use Discretization_module
  use Surface_Field_module
  use Surface_TH_module

  implicit none

  class(pm_surface_th_type) :: this

  PetscReal, pointer :: xx_p(:)
  PetscInt :: local_id
  PetscInt :: istart, iend
  type(surface_field_type), pointer :: surf_field 
  type(grid_type),pointer :: surf_grid
  PetscErrorCode :: ierr

  surf_grid => this%surf_realization%discretization%grid
  surf_field => this%surf_realization%surf_field

  ! Ensure evolved solution is +ve
  call VecGetArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  do local_id = 1,this%surf_realization%discretization%grid%nlmax
    iend = local_id*this%option%nflowdof
    istart = iend - this%option%nflowdof + 1
    if (xx_p(istart) < 1.d-15) then
      xx_p(istart) = 0.d0
      xx_p(iend) = 0.d0
    endif
  enddo
  call VecRestoreArrayF90(surf_field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  ! First, update the solution vector
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

  ! Update aux vars
  call SurfaceTHUpdateTemperature(this%surf_realization)
  call SurfaceTHUpdateAuxVars(this%surf_realization)

  ! Update the temperature due to atmospheric forcing using an implicit
  ! time integration.
  call SurfaceTHImplicitAtmForcing(this%surf_realization)
  call SurfaceTHUpdateAuxVars(this%surf_realization)
  call DiscretizationGlobalToLocal(this%surf_realization%discretization, &
          surf_field%flow_xx,surf_field%flow_xx_loc,NFLOWDOF)

end subroutine PMSurfaceTHPostSolve

! ************************************************************************** !
subroutine PMSurfaceTHUpdateAuxVars(this)
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 04/30/14

  use Surface_TH_module

  implicit none

  class(pm_surface_th_type) :: this

  call SurfaceTHUpdateAuxVars(this%surf_realization)

end subroutine PMSurfaceTHUpdateAuxVars

! ************************************************************************** !

subroutine PMSurfaceTHInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_surface_th_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMSurfaceTHInputRecord

! ************************************************************************** !

subroutine PMSurfaceTHDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 07/23/13
  ! 

  use Surface_TH_module, only : SurfaceTHDestroy

  implicit none

  class(pm_surface_th_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

#ifdef PM_SURFACE_FLOW_DEBUG
  call printMsg(this%option,'PMSurfaceTHDestroy()')
#endif

  call SurfaceTHDestroy(this%surf_realization)
  call this%comm1%Destroy()

end subroutine PMSurfaceTHDestroy

end module PM_Surface_TH_class
