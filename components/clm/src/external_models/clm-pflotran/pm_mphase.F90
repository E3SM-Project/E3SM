module PM_Mphase_class

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

  type, public, extends(pm_subsurface_flow_type) :: pm_mphase_type
  contains
    procedure, public :: Read => PMMphaseRead
    procedure, public :: InitializeTimestep => PMMphaseInitializeTimestep
    procedure, public :: Residual => PMMphaseResidual
    procedure, public :: Jacobian => PMMphaseJacobian
    procedure, public :: UpdateTimestep => PMMphaseUpdateTimestep
    procedure, public :: PreSolve => PMMphasePreSolve
    procedure, public :: PostSolve => PMMphasePostSolve
#if 0
    procedure, public :: CheckUpdatePre => PMMphaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMphaseCheckUpdatePost
#endif
    procedure, public :: TimeCut => PMMphaseTimeCut
    procedure, public :: UpdateSolution => PMMphaseUpdateSolution
    procedure, public :: UpdateAuxVars => PMMphaseUpdateAuxVars
    procedure, public :: MaxChange => PMMphaseMaxChange
    procedure, public :: ComputeMassBalance => PMMphaseComputeMassBalance
    procedure, public :: InputRecord => PMMphaseInputRecord
    procedure, public :: Destroy => PMMphaseDestroy
  end type pm_mphase_type
  
  public :: PMMphaseCreate
  
contains

! ************************************************************************** !

function PMMphaseCreate()
  ! 
  ! Creates Mphase process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_mphase_type), pointer :: PMMphaseCreate

  class(pm_mphase_type), pointer :: mphase_pm

  allocate(mphase_pm)

  call PMSubsurfaceFlowCreate(mphase_pm)
  mphase_pm%name = 'Mphase CO2 Flow'

  PMMphaseCreate => mphase_pm
  
end function PMMphaseCreate

! ************************************************************************** !

subroutine PMMphaseRead(this,input)
  ! 
  ! Reads input file parameters associated with the Mphase process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use Mphase_Aux_module
 
  implicit none
  
  class(pm_mphase_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'Mphase Options'
  
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
  
end subroutine PMMphaseRead

! ************************************************************************** !

subroutine PMMphaseInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseInitializeTimestep
  
  implicit none
  
  class(pm_mphase_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)         

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," MPHASE FLOW ",65("="))')
  endif
  
  call MphaseInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)         
  
end subroutine PMMphaseInitializeTimestep

! ************************************************************************** !

subroutine PMMphasePreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Grid_module
  use Patch_module
  use Global_Aux_module
  use Coupler_module
  use Connection_module  

  implicit none
  
  class(pm_mphase_type) :: this

  type(reaction_type), pointer :: reaction
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscInt :: na_id, cl_id
  
  reaction => this%realization%reaction
  option => this%realization%option
#if 1
  if (associated(reaction)) then
    if (associated(reaction%species_idx)) then
      patch => this%realization%patch
      global_auxvars => patch%aux%Global%auxvars
      if (associated(global_auxvars(1)%m_nacl)) then
        na_id = reaction%species_idx%na_ion_id
        cl_id = reaction%species_idx%cl_ion_id
 
        grid => patch%grid
        rt_auxvars => patch%aux%RT%auxvars
        global_auxvars => patch%aux%Global%auxvars

        if (na_id > 0 .and. cl_id > 0) then
          do ghosted_id = 1, grid%ngmax
            if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
            !geh - Ignore inactive cells with inactive materials
            if (patch%imat(ghosted_id) <= 0) cycle
            global_auxvars(ghosted_id)%m_nacl(1) = &
              rt_auxvars(ghosted_id)%pri_molal(na_id)
            global_auxvars(ghosted_id)%m_nacl(2) = &
              rt_auxvars(ghosted_id)%pri_molal(cl_id)
          enddo
        else    
          do ghosted_id = 1, grid%ngmax
            if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
            !geh - Ignore inactive cells with inactive materials
            if (patch%imat(ghosted_id) <= 0) cycle
            global_auxvars(ghosted_id)%m_nacl = option%m_nacl
          enddo
        endif
      
        rt_auxvars => this%realization%patch%aux%RT%auxvars_bc
        global_auxvars => this%realization%patch%aux%Global%auxvars_bc
    
        boundary_condition => patch%boundary_condition_list%first
        sum_connection = 0 
        do 
          if (.not.associated(boundary_condition)) exit
          cur_connection_set => boundary_condition%connection_set
          if (na_id > 0 .and. cl_id > 0) then
            do iconn = 1, cur_connection_set%num_connections
              sum_connection = sum_connection + 1
              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle
              global_auxvars(sum_connection)%m_nacl(1) = &
                rt_auxvars(sum_connection)%pri_molal(na_id)
              global_auxvars(sum_connection)%m_nacl(2) = &
                rt_auxvars(sum_connection)%pri_molal(cl_id)
            enddo
          else    
            do iconn = 1, cur_connection_set%num_connections
              sum_connection = sum_connection + 1
              local_id = cur_connection_set%id_dn(iconn)
              ghosted_id = grid%nL2G(local_id)
              if (patch%imat(ghosted_id) <= 0) cycle
              global_auxvars(sum_connection)%m_nacl = option%m_nacl
            enddo
          endif        
          boundary_condition => boundary_condition%next
        enddo
      endif
    endif
  endif   
#endif

end subroutine PMMphasePreSolve

! ************************************************************************** !

subroutine PMMphasePostSolve(this)
  ! 
  ! PMMphaseUpdatePostSolve:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_mphase_type) :: this
  
end subroutine PMMphasePostSolve

! ************************************************************************** !

subroutine PMMphaseUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_mphase_type) :: this
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
  
end subroutine PMMphaseUpdateTimestep

! ************************************************************************** !

subroutine PMMphaseResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseResidual

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call MphaseResidual(snes,xx,r,this%realization,ierr)

end subroutine PMMphaseResidual

! ************************************************************************** !

subroutine PMMphaseJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseJacobian

  implicit none
  
  class(pm_mphase_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call MphaseJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMMphaseJacobian
    
#if 0

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePre

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  call MphaseCheckUpdatePre(line_search,X,dX,changed,this%realization,ierr)

end subroutine PMMphaseCheckUpdatePre

! ************************************************************************** !

subroutine PMMphaseCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                   X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseCheckUpdatePost

  implicit none
  
  class(pm_mphase_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  call MphaseCheckUpdatePost(line_search,X0,dX,X1,dX_changed, &
                               X1_changed,this%realization,ierr)

end subroutine PMMphaseCheckUpdatePost
#endif

! ************************************************************************** !

subroutine PMMphaseTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseTimeCut

  implicit none
  
  class(pm_mphase_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call MphaseTimeCut(this%realization)

end subroutine PMMphaseTimeCut

! ************************************************************************** !

subroutine PMMphaseUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseUpdateSolution

  implicit none
  
  class(pm_mphase_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call MphaseUpdateSolution(this%realization)

end subroutine PMMphaseUpdateSolution     

! ************************************************************************** !

subroutine PMMphaseUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Mphase_module, only : MphaseUpdateAuxVars
  
  implicit none
  
  class(pm_mphase_type) :: this

  call MphaseUpdateAuxVars(this%realization)

end subroutine PMMphaseUpdateAuxVars   

! ************************************************************************** !

subroutine PMMphaseMaxChange(this)
  ! 
  ! Not needed given MphaseMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseMaxChange

  implicit none
  
  class(pm_mphase_type) :: this
  
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

end subroutine PMMphaseMaxChange

! ************************************************************************** !

subroutine PMMphaseComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseComputeMassBalance

  implicit none
  
  class(pm_mphase_type) :: this
  PetscReal :: mass_balance_array(:)
  
  !geh: currently does not include "trapped" mass
  !call MphaseComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMMphaseComputeMassBalance

! ************************************************************************** !

subroutine PMMphaseInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_mphase_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'mphase'

end subroutine PMMphaseInputRecord

! ************************************************************************** !

subroutine PMMphaseDestroy(this)
  ! 
  ! Destroys Mphase process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Mphase_module, only : MphaseDestroy

  implicit none
  
  class(pm_mphase_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call MphaseDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMMphaseDestroy
  
end module PM_Mphase_class
