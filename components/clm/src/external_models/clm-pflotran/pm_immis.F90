module PM_Immis_class

  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"

  type, public, extends(pm_subsurface_flow_type) :: pm_immis_type
  contains
    procedure, public :: Read => PMImmisRead
    procedure, public :: InitializeTimestep => PMImmisInitializeTimestep
    procedure, public :: Residual => PMImmisResidual
    procedure, public :: Jacobian => PMImmisJacobian
    procedure, public :: UpdateTimestep => PMImmisUpdateTimestep
    procedure, public :: PreSolve => PMImmisPreSolve
    procedure, public :: PostSolve => PMImmisPostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMImmisCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMImmisCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMImmisTimeCut
    procedure, public :: UpdateSolution => PMImmisUpdateSolution
    procedure, public :: UpdateAuxVars => PMImmisUpdateAuxVars
    procedure, public :: MaxChange => PMImmisMaxChange
    procedure, public :: ComputeMassBalance => PMImmisComputeMassBalance
    procedure, public :: InputRecord => PMImmisInputRecord
    procedure, public :: Destroy => PMImmisDestroy
  end type pm_immis_type
  
  public :: PMImmisCreate
  
contains

! ************************************************************************** !

function PMImmisCreate()
  ! 
  ! Creates Immiscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type), pointer :: PMImmisCreate

  class(pm_immis_type), pointer :: immis_pm
  
  allocate(immis_pm)

  call PMSubsurfaceFlowCreate(immis_pm)
  immis_pm%name = 'Immisible Flow'

  PMImmisCreate => immis_pm
  
end function PMImmisCreate

! ************************************************************************** !

subroutine PMImmisRead(this,input)
  ! 
  ! Reads input file parameters associated with the Immis process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use Immis_Aux_module
 
  implicit none
  
  class(pm_immis_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'Immiscible Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,word,found,option)
    if (found) cycle
    
    select case(trim(word))
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMImmisRead

! ************************************************************************** !

subroutine PMImmisInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisInitializeTimestep
  
  implicit none
  
  class(pm_immis_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," IMMISCIBLE FLOW ",61("="))')
  endif
  
  call ImmisInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)         
  
end subroutine PMImmisInitializeTimestep

! ************************************************************************** !

subroutine PMImmisPreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
end subroutine PMImmisPreSolve

! ************************************************************************** !

subroutine PMImmisPostSolve(this)
  ! 
  ! PMImmisUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
end subroutine PMImmisPostSolve

! ************************************************************************** !

subroutine PMImmisUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_immis_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: uus
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      utmp = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
      uus= this%saturation_change_governor/(this%max_saturation_change+1.d-6)  
      ut = min(up,utmp,uus)
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    dt_p = fac * dt * (1.d0 + up)

    dtt = min(dt_tfac,dt_p)
  endif
  
  if (dtt > 2.d0 * dt) dtt = 2.d0 * dt
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  
end subroutine PMImmisUpdateTimestep

! ************************************************************************** !

subroutine PMImmisResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisResidual

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call ImmisResidual(snes,xx,r,this%realization,ierr)

end subroutine PMImmisResidual

! ************************************************************************** !

subroutine PMImmisJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisJacobian

  implicit none
  
  class(pm_immis_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call ImmisJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMImmisJacobian
    
#if 0

! ************************************************************************** !

subroutine PMImmisCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePre

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call ImmisCheckUpdatePre(line_search,X,dX,changed,this%realization,ierr)

end subroutine PMImmisCheckUpdatePre

! ************************************************************************** !

subroutine PMImmisCheckUpdatePost(this,line_search,P0,dP,P1,dX_changed, &
                                  X1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisCheckUpdatePost

  implicit none
  
  class(pm_immis_type) :: this
  SNESLineSearch :: line_search
  Vec :: P0
  Vec :: dP
  Vec :: P1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  call ImmisCheckUpdatePost(line_search,P0,dP,P1,dX_changed, &
                               X1_changed,this%realization,ierr)

end subroutine PMImmisCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMImmisTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisTimeCut

  implicit none
  
  class(pm_immis_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call ImmisTimeCut(this%realization)

end subroutine PMImmisTimeCut

! ************************************************************************** !

subroutine PMImmisUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisUpdateSolution

  implicit none
  
  class(pm_immis_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call ImmisUpdateSolution(this%realization)

end subroutine PMImmisUpdateSolution     

! ************************************************************************** !

subroutine PMImmisUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Immis_module, only : ImmisUpdateAuxVars
    
  implicit none
  
  class(pm_immis_type) :: this

  call ImmisUpdateAuxVars(this%realization)

end subroutine PMImmisUpdateAuxVars   

! ************************************************************************** !

subroutine PMImmisMaxChange(this)
  ! 
  ! Not needed given PMImmisMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisMaxChange

  implicit none
  
  class(pm_immis_type) :: this
  
  call ImmisMaxChange(this%realization,this%max_pressure_change, &
                      this%max_temperature_change,this%max_saturation_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%max_pressure_change,this%max_temperature_change, &
          this%max_saturation_change
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%max_pressure_change,this%max_temperature_change, &
          this%max_saturation_change
  endif  

end subroutine PMImmisMaxChange

! ************************************************************************** !

subroutine PMImmisComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisComputeMassBalance

  implicit none
  
  class(pm_immis_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call ImmisComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMImmisComputeMassBalance

! ************************************************************************** !

recursive subroutine PMImmisFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_immis_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMImmisFinalizeRun

! ************************************************************************** !

subroutine PMImmisInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_immis_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'immiscible'

end subroutine PMImmisInputRecord

! ************************************************************************** !

subroutine PMImmisDestroy(this)
  ! 
  ! Destroys Immiscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Immis_module, only : ImmisDestroy

  implicit none
  
  class(pm_immis_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call ImmisDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMImmisDestroy
  
end module PM_Immis_class
