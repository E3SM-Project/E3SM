module PM_TH_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
!geh: using TH_module here fails with gfortran (internal compiler error)
!  use TH_module
  use Realization_Subsurface_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_th_type
    class(communicator_type), pointer :: commN
  contains
    procedure, public :: Setup => PMTHSetup
    procedure, public :: Read => PMTHRead
    procedure, public :: InitializeTimestep => PMTHInitializeTimestep
    procedure, public :: Residual => PMTHResidual
    procedure, public :: Jacobian => PMTHJacobian
    procedure, public :: UpdateTimestep => PMTHUpdateTimestep
    procedure, public :: PreSolve => PMTHPreSolve
    procedure, public :: PostSolve => PMTHPostSolve
    procedure, public :: CheckUpdatePre => PMTHCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHCheckUpdatePost
    procedure, public :: TimeCut => PMTHTimeCut
    procedure, public :: UpdateSolution => PMTHUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHUpdateAuxVars
    procedure, public :: MaxChange => PMTHMaxChange
    procedure, public :: ComputeMassBalance => PMTHComputeMassBalance
    procedure, public :: InputRecord => PMTHInputRecord
    procedure, public :: Destroy => PMTHDestroy
  end type pm_th_type
  
  public :: PMTHCreate
  
contains

! ************************************************************************** !

function PMTHCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_th_type), pointer :: PMTHCreate

  class(pm_th_type), pointer :: th_pm
  
#ifdef PM_TH_DEBUG
  print *, 'PMTHCreate()'
#endif  

  allocate(th_pm)

  nullify(th_pm%commN)

  call PMSubsurfaceFlowCreate(th_pm)
  th_pm%name = 'TH Flow'

  PMTHCreate => th_pm
  
end function PMTHCreate

! ************************************************************************** !

subroutine PMTHRead(this,input)
  ! 
  ! Reads input file parameters associated with the TH process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use TH_Aux_module
 
  implicit none
  
  class(pm_th_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option
  
  error_string = 'TH Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,word,found, &
                                        error_string,option)
    if (found) cycle
    
    select case(trim(word))
      case('ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,th_itol_scaled_res)
        call InputDefaultMsg(input,option,'th_itol_scaled_res')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,th_itol_rel_update)
        call InputDefaultMsg(input,option,'th_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE        
      case('FREEZING')
        option%use_th_freezing = PETSC_TRUE
        option%io_buffer = ' TH: using FREEZING submode!'
        call printMsg(option)
        ! Override the default setting for TH-mode with freezing
        call EOSWaterSetDensity('PAINTER')
        call EOSWaterSetEnthalpy('PAINTER')
      case('ICE_MODEL')
        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('PAINTER_EXPLICIT')
            option%ice_model = PAINTER_EXPLICIT
          case ('PAINTER_KARRA_IMPLICIT')
            option%ice_model = PAINTER_KARRA_IMPLICIT
          case ('PAINTER_KARRA_EXPLICIT')
            option%ice_model = PAINTER_KARRA_EXPLICIT
          case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
            option%ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
          case ('PAINTER_KARRA_EXPLICIT_SMOOTH')
            option%ice_model = PAINTER_KARRA_EXPLICIT_SMOOTH
            call InputReadDouble(input,option,tempreal)
            call InputDefaultMsg(input,option,'freezing-thawing smoothing')
            if(tempreal > 1.d-10) option%frzthw_halfwidth = tempreal
          case ('DALL_AMICO')
            option%ice_model = DALL_AMICO
          case default
            option%io_buffer = 'Cannot identify the specificed ice model.' // &
             'Specify PAINTER_EXPLICIT or PAINTER_KARRA_IMPLICIT' // &
             ' or PAINTER_KARRA_EXPLICIT or PAINTER_KARRA_EXPLICIT_NOCRYO ' // &
             ' or PAINTER_KARRA_EXPLICIT_SMOOTH or DALL_AMICO.'
            call printErrMsg(option)
          end select
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMTHRead

! ************************************************************************** !

subroutine PMTHSetup(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowSetup(this)
  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%commN%SetDM(this%realization%discretization%dm_nflowdof)

end subroutine PMTHSetup

! ************************************************************************** !

subroutine PMTHInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THInitializeTimestep

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%icap_loc, &
                               this%realization%field%icap_loc)
  call this%comm1%LocalToLocal(this%realization%field%ithrm_loc, &
                               this%realization%field%ithrm_loc)
  call this%comm1%LocalToLocal(this%realization%field%iphas_loc, &
                               this%realization%field%iphas_loc)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," TH FLOW ",69("="))')
  endif
  
  call THInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)
  
end subroutine PMTHInitializeTimestep

! ************************************************************************** !

subroutine PMTHPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMTHPreSolve

! ************************************************************************** !

subroutine PMTHPostSolve(this)
  ! 
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this
  
end subroutine PMTHPostSolve

! ************************************************************************** !

subroutine PMTHUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac, &
                              time_step_max_growth_factor)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: dtt
  PetscReal :: dt_u
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_TH_DEBUG
  call printMsg(this%option,'PMTH%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      utmp = this%temperature_change_governor/ &
             (this%max_temperature_change+1.d-5)
      ut = min(up,utmp)
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    utmp = this%temperature_change_governor/ &
           (this%max_temperature_change+1.d-5)
    ut = min(up,utmp)
    dt_u = fac * dt * (1.d0 + ut)

    dtt = min(dt_tfac,dt_u)
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  
end subroutine PMTHUpdateTimestep

! ************************************************************************** !

subroutine PMTHResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THResidual

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call THResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHResidual

! ************************************************************************** !

subroutine PMTHJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THJacobian

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call THJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTHJacobian

! ************************************************************************** !

subroutine PMTHCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module
  use TH_Aux_module
  use Global_Aux_module

  implicit none
  
  class(pm_th_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P0, P1, P_R, delP, delP_old
  PetscReal :: scale, press_limit, temp_limit
  PetscInt :: iend, istart
  
  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars

  if (Initialized(this%pressure_change_limit)) then

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    press_limit = dabs(this%pressure_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = X_p(istart)
      delP = dX_p(istart)
      if (press_limit < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &         
          grid%nG2A(grid%nL2G(local_id)),press_limit,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),press_limit),delP)
      dX_p(istart) = delP
    enddo
    
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)

  endif
  
  if (dabs(this%temperature_change_limit) > 0.d0) then
      
    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    temp_limit = dabs(this%temperature_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = X_p(iend)
      delP = dX_p(iend)
      if (abs(delP) > abs(temp_limit)) then
        write(option%io_buffer,'("dT_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),temp_limit,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),temp_limit),delP)
      dX_p(iend) = delP
    enddo
    
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    
  endif


  if (Initialized(this%pressure_dampening_factor)) then
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = this%pressure_dampening_factor

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      delP = dX_p(istart)
      P0 = X_p(istart)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1 
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call printMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dX_p(istart) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine PMTHCheckUpdatePre

! ************************************************************************** !

subroutine PMTHCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Secondary_Continuum_Aux_module
  use TH_module
  use TH_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Patch_module

  implicit none
  
  class(pm_th_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: ithrm_loc_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(TH_parameter_type), pointer :: TH_parameter
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscReal :: Res(2)
  PetscReal :: inf_norm, global_inf_norm
  PetscReal :: vol_frac_prim
  PetscInt :: istart, iend
  
  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  TH_auxvars => patch%aux%TH%auxvars
  TH_parameter => patch%aux%TH%TH_parameter
  global_auxvars => patch%aux%Global%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  material_auxvars => patch%aux%Material%auxvars
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  if (this%check_post_convergence) then
    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%ithrm_loc,ithrm_loc_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    
    inf_norm = 0.d0
    vol_frac_prim = 1.d0
    
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
    
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      
      if (option%use_mc) then
        vol_frac_prim = TH_sec_heat_vars(local_id)%epsilon
      endif

      call THAccumulation(TH_auxvars(ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           TH_parameter%dencpr(int(ithrm_loc_p(ghosted_id))), &
                           option,vol_frac_prim,Res)
                                                        
      inf_norm = max(inf_norm,min(dabs(dX_p(istart)/X1_p(istart)), &
                                  dabs(r_p(istart)/Res(1)), &
                                  dabs(dX_p(iend)/X1_p(iend)), &
                                  dabs(r_p(iend)/Res(2))))
    enddo
    call MPI_Allreduce(inf_norm,global_inf_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    option%converged = PETSC_TRUE
    if (global_inf_norm > th_itol_scaled_res) &
      option%converged = PETSC_FALSE
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMTHCheckUpdatePost

! ************************************************************************** !

subroutine PMTHTimeCut(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THTimeCut

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call THTimeCut(this%realization)

end subroutine PMTHTimeCut

! ************************************************************************** !

subroutine PMTHUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THUpdateSolution, THUpdateSurfaceBC

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call THUpdateSolution(this%realization)
  if (this%option%surf_flow_on) &
    call THUpdateSurfaceBC(this%realization)

end subroutine PMTHUpdateSolution     

! ************************************************************************** !

subroutine PMTHUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use TH_module, only : THUpdateAuxVars
  
  implicit none
  
  class(pm_th_type) :: this

  call THUpdateAuxVars(this%realization)

end subroutine PMTHUpdateAuxVars   

! ************************************************************************** !

subroutine PMTHMaxChange(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THMaxChange

  implicit none
  
  class(pm_th_type) :: this
  
  call THMaxChange(this%realization,this%max_pressure_change, &
                   this%max_temperature_change)
#ifndef CLM_PFLOTRAN
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4," dtmpmx= ",1pe12.4)') &
      this%max_pressure_change,this%max_temperature_change
  endif
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4, &
      & " dtmpmx= ",1pe12.4)') &
      this%max_pressure_change,this%max_temperature_change
  endif 
#endif
end subroutine PMTHMaxChange

! ************************************************************************** !

subroutine PMTHComputeMassBalance(this,mass_balance_array)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THComputeMassBalance

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call THComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHComputeMassBalance

! ************************************************************************** !

subroutine PMTHInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_th_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'thermo-hydro'

end subroutine PMTHInputRecord

! ************************************************************************** !

subroutine PMTHDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call this%commN%Destroy()

  ! preserve this ordering
  call THDestroy(this%realization%patch)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHDestroy

end module PM_TH_class
