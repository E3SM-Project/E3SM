module PM_Miscible_class

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

  type, public, extends(pm_subsurface_flow_type) :: pm_miscible_type
  contains
    procedure, public :: Read => PMMiscibleRead
    procedure, public :: InitializeTimestep => PMMiscibleInitializeTimestep
    procedure, public :: Residual => PMMiscibleResidual
    procedure, public :: Jacobian => PMMiscibleJacobian
    procedure, public :: UpdateTimestep => PMMiscibleUpdateTimestep
    procedure, public :: PreSolve => PMMisciblePreSolve
    procedure, public :: PostSolve => PMMisciblePostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMMiscibleCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMiscibleCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMiscibleTimeCut
    procedure, public :: UpdateSolution => PMMiscibleUpdateSolution
    procedure, public :: UpdateAuxVars => PMMiscibleUpdateAuxVars
    procedure, public :: MaxChange => PMMiscibleMaxChange
    procedure, public :: ComputeMassBalance => PMMiscibleComputeMassBalance
    procedure, public :: InputRecord => PMMiscibleInputRecord
    procedure, public :: Destroy => PMMiscibleDestroy
  end type pm_miscible_type
  
  public :: PMMiscibleCreate
  
contains

! ************************************************************************** !

function PMMiscibleCreate()
  ! 
  ! Creates Miscible process models shell
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type), pointer :: PMMiscibleCreate

  class(pm_miscible_type), pointer :: miscible_pm
  
#ifdef PM__DEBUG  
  print *, 'PMMiscibleCreate()'
#endif  

  allocate(miscible_pm)

  call PMSubsurfaceFlowCreate(miscible_pm)
  miscible_pm%name = 'Miscible Flow'

  PMMiscibleCreate => miscible_pm
  
end function PMMiscibleCreate

! ************************************************************************** !

subroutine PMMiscibleRead(this,input)
  ! 
  ! Reads input file parameters associated with the Miscible process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use Miscible_Aux_module
 
  implicit none
  
  class(pm_miscible_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'Miscible Options'
  
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
  
end subroutine PMMiscibleRead

! ************************************************************************** !

subroutine PMMiscibleInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleInitializeTimestep
  
  implicit none
  
  class(pm_miscible_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MISCIBLE FLOW ",63("="))')
  endif
  
  call MiscibleInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)         
  
end subroutine PMMiscibleInitializeTimestep

! ************************************************************************** !

subroutine PMMisciblePreSolve(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
end subroutine PMMisciblePreSolve

! ************************************************************************** !

subroutine PMMisciblePostSolve(this)
  ! 
  ! PMMiscibleUpdatePostSolve:
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  implicit none
  
  class(pm_miscible_type) :: this
  
end subroutine PMMisciblePostSolve

! ************************************************************************** !

subroutine PMMiscibleUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_miscible_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: uc
  PetscReal :: uus
  PetscReal :: dtt
  PetscReal :: dt_p
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_MISCIBLE_DEBUG  
  call printMsg(this%option,'PMMiscible%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      utmp = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
      uc = this%xmol_change_governor/(this%max_xmol_change+1.d-6)
      uus= this%saturation_change_governor/(this%max_saturation_change+1.d-6)      
      ut = min(up,utmp,uc,uus)
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
  
end subroutine PMMiscibleUpdateTimestep

! ************************************************************************** !

subroutine PMMiscibleResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleResidual

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call MiscibleResidual(snes,xx,r,this%realization,ierr)

end subroutine PMMiscibleResidual

! ************************************************************************** !

subroutine PMMiscibleJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleJacobian

  implicit none
  
  class(pm_miscible_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call MiscibleJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMMiscibleJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePre

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call MiscibleCheckUpdatePre(line_search,X,dX,changed,this%realization,ierr)

end subroutine PMMiscibleCheckUpdatePre

! ************************************************************************** !

subroutine PMMiscibleCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleCheckUpdatePost

  implicit none
  
  class(pm_miscible_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  call MiscibleCheckUpdatePost(line_search,X0,dX,X1,dX_changed, &
                               X1_changed,this%realization,ierr)

end subroutine PMMiscibleCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMiscibleTimeCut(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleTimeCut

  implicit none
  
  class(pm_miscible_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call MiscibleTimeCut(this%realization)

end subroutine PMMiscibleTimeCut

! ************************************************************************** !

subroutine PMMiscibleUpdateSolution(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleUpdateSolution

  implicit none
  
  class(pm_miscible_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call MiscibleUpdateSolution(this%realization)

end subroutine PMMiscibleUpdateSolution     

! ************************************************************************** !

subroutine PMMiscibleUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Miscible_module, only : MiscibleUpdateAuxVars
  
  implicit none
  
  class(pm_miscible_type) :: this

  call MiscibleUpdateAuxVars(this%realization)

end subroutine PMMiscibleUpdateAuxVars   

! ************************************************************************** !

subroutine PMMiscibleMaxChange(this)
  ! 
  ! Not needed given PMMiscibleMaxChange is called in PostSolve
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Mphase_module, only : MphaseMaxChange

  implicit none
  
  class(pm_miscible_type) :: this
  
  !geh: yes, call Mphase.  No need to replicate code  
  call MphaseMaxChange(this%realization,this%max_pressure_change, &
                       this%max_temperature_change, &
                       this%max_saturation_change, &
                       this%max_xmol_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%max_pressure_change,this%max_temperature_change, &
          this%max_xmol_change,this%max_saturation_change
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4," dcmx= ",1pe12.4," dsmx= ",1pe12.4)') &
          this%max_pressure_change,this%max_temperature_change, &
          this%max_xmol_change,this%max_saturation_change
  endif     

end subroutine PMMiscibleMaxChange

! ************************************************************************** !

subroutine PMMiscibleComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleComputeMassBalance

  implicit none
  
  class(pm_miscible_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call MiscibleComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMMiscibleComputeMassBalance

! ************************************************************************** !

subroutine PMMiscibleInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_miscible_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'miscible'

end subroutine PMMiscibleInputRecord

! ************************************************************************** !

subroutine PMMiscibleDestroy(this)
  ! 
  ! Destroys Miscible process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 11/27/13
  ! 

  use Miscible_module, only : MiscibleDestroy

  implicit none
  
  class(pm_miscible_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call MiscibleDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMMiscibleDestroy
  
end module PM_Miscible_class
