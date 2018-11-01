module PM_RT_class
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
!geh: using Reactive_Transport_module here fails with gfortran (internal 
!     compiler error)
!  use Reactive_Transport_module
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_rt_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    class(communicator_type), pointer :: commN
    ! local variables
    PetscBool :: steady_flow
    PetscReal :: tran_weight_t0
    PetscReal :: tran_weight_t1
    PetscBool :: check_post_convergence
    ! these govern the size of subsequent time steps
    PetscReal :: max_concentration_change
    PetscReal :: max_volfrac_change
    PetscReal :: volfrac_change_governor
    PetscReal :: cfl_governor
    PetscBool :: temperature_dependent_diffusion
    ! for transport only
    PetscBool :: transient_porosity
    PetscBool :: include_gas_phase
  contains
    procedure, public :: Setup => PMRTSetup
    procedure, public :: Read => PMRTRead
    procedure, public :: PMRTSetRealization
    procedure, public :: InitializeRun => PMRTInitializeRun
    procedure, public :: FinalizeRun => PMRTFinalizeRun
    procedure, public :: InitializeTimestep => PMRTInitializeTimestep
    procedure, public :: FinalizeTimestep => PMRTFinalizeTimestep
    procedure, public :: Residual => PMRTResidual
    procedure, public :: Jacobian => PMRTJacobian
    procedure, public :: UpdateTimestep => PMRTUpdateTimestep
    procedure, public :: PreSolve => PMRTPreSolve
    procedure, public :: PostSolve => PMRTPostSolve
    procedure, public :: AcceptSolution => PMRTAcceptSolution
    procedure, public :: CheckUpdatePre => PMRTCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRTCheckUpdatePost
    procedure, public :: CheckConvergence => PMRTCheckConvergence
    procedure, public :: TimeCut => PMRTTimeCut
    procedure, public :: UpdateSolution => PMRTUpdateSolution1
    procedure, public :: UpdateAuxVars => PMRTUpdateAuxVars
    procedure, public :: MaxChange => PMRTMaxChange
    procedure, public :: ComputeMassBalance => PMRTComputeMassBalance
    procedure, public :: SetTranWeights => SetTranWeights
    procedure, public :: CheckpointBinary => PMRTCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMRTCheckpointHDF5
    procedure, public :: RestartBinary => PMRTRestartBinary
    procedure, public :: RestartHDF5 => PMRTRestartHDF5
    procedure, public :: InputRecord => PMRTInputRecord
    procedure, public :: Destroy => PMRTDestroy
  end type pm_rt_type
  
  type, public, extends(pm_base_header_type) :: pm_rt_header_type
    PetscInt :: checkpoint_activity_coefs
  end type pm_rt_header_type  
  
  public :: PMRTCreate

contains

! ************************************************************************** !

function PMRTCreate()
  ! 
  ! Creates reactive transport process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_rt_type), pointer :: PMRTCreate

  class(pm_rt_type), pointer :: rt_pm
  
#ifdef PM_RT_DEBUG  
  print *, 'PMRTCreate()'
#endif
  
  allocate(rt_pm)
  nullify(rt_pm%option)
  nullify(rt_pm%output_option)
  nullify(rt_pm%realization)
  nullify(rt_pm%comm1)
  nullify(rt_pm%commN)
  
  ! local variables
  rt_pm%steady_flow = PETSC_FALSE
  rt_pm%tran_weight_t0 = 0.d0
  rt_pm%tran_weight_t1 = 0.d0
  rt_pm%check_post_convergence = PETSC_FALSE
  rt_pm%max_concentration_change = 0.d0
  rt_pm%max_volfrac_change = 0.d0
  rt_pm%volfrac_change_governor = 1.d0
  rt_pm%cfl_governor = UNINITIALIZED_DOUBLE
  rt_pm%temperature_dependent_diffusion = PETSC_FALSE
  ! these flags can only be true for transport only
  rt_pm%transient_porosity = PETSC_FALSE
  rt_pm%include_gas_phase = PETSC_FALSE

  call PMBaseInit(rt_pm)
  rt_pm%name = 'Reactive Transport'
  
  PMRTCreate => rt_pm
  
end function PMRTCreate

! ************************************************************************** !

subroutine PMRTRead(this,input)
  ! 
  ! Reads input file parameters associated with the reactive transport 
  ! process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/25/16
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Reactive_Transport_Aux_module
 
  implicit none
  
  class(pm_rt_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  option => this%option
  
  error_string = 'Reactive Transport Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    
    select case(trim(word))
      case('GLOBAL_IMPLICIT','OPERATOR_SPLIT','OPERATOR_SPLITTING')
      case('MAX_VOLUME_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%volfrac_change_governor)
        call InputDefaultMsg(input,option,'maximum volume fraction change')
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,rt_itol_rel_update)
        call InputDefaultMsg(input,option,'rt_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE
      case('NUMERICAL_JACOBIAN')
        option%transport%numerical_derivatives = PETSC_TRUE
      case('INCLUDE_GAS_PHASE')
        this%include_gas_phase = PETSC_TRUE
        option%transport%nphase = 2
      case('TEMPERATURE_DEPENDENT_DIFFUSION')
        this%temperature_dependent_diffusion = PETSC_TRUE
      case('MAX_CFL')
        call InputReadDouble(input,option,this%cfl_governor)
        call InputErrorMsg(input,option,'MAX_CFL', &
                           'SUBSURFACE_TRANSPORT OPTIONS')
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMRTRead

! ************************************************************************** !

subroutine PMRTSetup(this)
  ! 
  ! Initializes variables associated with reactive transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

#ifndef SIMPLIFY
  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 
#endif  
  use Reactive_Transport_Aux_module, only : reactive_transport_param_type
  
  implicit none
  
  class(pm_rt_type) :: this

  type(reactive_transport_param_type), pointer :: rt_parameter

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Setup()')
#endif

  rt_parameter => this%realization%patch%aux%RT%rt_parameter
  
  ! pass down flags from PMRT class
  rt_parameter%temperature_dependent_diffusion = &
    this%temperature_dependent_diffusion

#ifndef SIMPLIFY  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%commN%SetDM(this%realization%discretization%dm_ntrandof)
#endif

  ! set the communicator
  this%comm1 => this%realization%comm1

  ! only set these flags if transport only
  if (this%option%nflowdof == 0) then
    if (associated(this%realization%reaction)) then
      if (this%realization%reaction%update_porosity & !.or. &
!          this%realization%reaction%update_tortuosity .or. &
!          this%realization%reaction%update_mnrl_surf_with_porosity &
          ) then
        this%transient_porosity = PETSC_TRUE
      endif
    endif
  endif
  
end subroutine PMRTSetup

! ************************************************************************** !

subroutine PMRTSetRealization(this,realization)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Subsurface_class  

  implicit none
  
  class(pm_rt_type) :: this
  class(realization_subsurface_type), pointer :: realization

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%SetRealization()')
#endif
  
  this%realization => realization
  this%realization_base => realization
  
  if (realization%reaction%use_log_formulation) then
    this%solution_vec = realization%field%tran_log_xx
  else
    this%solution_vec = realization%field%tran_xx
  endif
  this%residual_vec = realization%field%tran_r
  
end subroutine PMRTSetRealization

! ************************************************************************** !

recursive subroutine PMRTInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Reactive_Transport_module, only : RTUpdateEquilibriumState
  use Condition_Control_module
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Reactive_Transport_module, only : RTUpdateAuxVars, &
                                        RTClearActivityCoefficients
  use Variables_module, only : POROSITY
  use Material_Aux_class, only : POROSITY_MINERAL 
  use Material_module, only : MaterialGetAuxVarVecLoc

  implicit none
  
  class(pm_rt_type) :: this
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%InitializeRun()')
#endif

  ! check for uninitialized flow variables
  call RealizUnInitializedVarsTran(this%realization)

  if (this%transient_porosity) then
    call RealizationCalcMineralPorosity(this%realization)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_MINERAL)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity0)
    call VecCopy(this%realization%field%porosity0, &
                 this%realization%field%porosity_t,ierr);CHKERRQ(ierr)
    call VecCopy(this%realization%field%porosity0, &
                 this%realization%field%porosity_tpdt,ierr);CHKERRQ(ierr)
  endif
  
  ! restart
  if (this%option%restart_flag .and. &
      this%option%overwrite_restart_transport) then
    call RTClearActivityCoefficients(this%realization)
    call CondControlAssignTranInitCond(this%realization)  
  endif
  
  ! update boundary concentrations so that activity coefficients can be 
  ! calculated at first time step
  call RTUpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  ! pass PETSC_FALSE to turn off update of kinetic state variables
  call PMRTUpdateSolution2(this,PETSC_FALSE)
  
#if 0
  if (this%option%jumpstart_kinetic_sorption .and. &
      this%option%time < 1.d-40) then
    ! only user jumpstart for a restarted simulation
    if (.not. this%option%restart_flag) then
      this%option%io_buffer = 'Only use JUMPSTART_KINETIC_SORPTION on a ' // &
        'restarted simulation.  ReactionEquilibrateConstraint() will ' // &
        'appropriately set sorbed initial concentrations for a normal ' // &
        '(non-restarted) simulation.'
      call printErrMsg(this%option)
    endif
  endif
#endif  
    
end subroutine PMRTInitializeRun

! ************************************************************************** !

subroutine PMRTInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTInitializeTimestep, &
                                        RTUpdateActivityCoefficients
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_TIMESTEP
  use Global_module
  use Material_module

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: time
 
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%InitializeTimestep()')
#endif
  
  this%option%tran_dt = this%option%dt

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," REACTIVE TRANSPORT ",58("="))')
  endif
  
  ! interpolate flow parameters/data
  ! this must remain here as these weighted values are used by both
  ! RTInitializeTimestep and RTTimeCut (which calls RTInitializeTimestep)
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call this%SetTranWeights()
    if (this%option%flow%transient_porosity) then
      ! weight material properties (e.g. porosity)
      call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                                 this%tran_weight_t0, &
                                 this%realization%field,this%comm1)
    endif
    ! set densities and saturations to t
    call GlobalWeightAuxVars(this%realization,this%tran_weight_t0)
  else if (this%transient_porosity) then
    this%tran_weight_t0 = 0.d0
    call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                               this%tran_weight_t0, &
                               this%realization%field,this%comm1)
  endif

  call RTInitializeTimestep(this%realization)

  if (this%realization%reaction%act_coef_update_frequency == &
      ACT_COEF_FREQUENCY_TIMESTEP) then
    call RTUpdateActivityCoefficients(this%realization,PETSC_TRUE,PETSC_TRUE)
  endif

end subroutine PMRTInitializeTimestep

! ************************************************************************** !

subroutine PMRTPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTUpdateTransportCoefs
  use Global_module  
  use Material_module
  use Data_Mediator_module

  implicit none
  
  class(pm_rt_type) :: this
  
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdatePreSolve()')
#endif
  
  ! set densities and saturations to t+dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    if (this%option%flow%transient_porosity) then
      ! weight material properties (e.g. porosity)
      call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                                 this%tran_weight_t1, &
                                 this%realization%field,this%comm1)
    endif
    call GlobalWeightAuxVars(this%realization,this%tran_weight_t1)
  else if (this%transient_porosity) then
    this%tran_weight_t1 = 1.d0
    call MaterialWeightAuxVars(this%realization%patch%aux%Material, &
                               this%tran_weight_t1, &
                               this%realization%field,this%comm1)
  endif

  call RTUpdateTransportCoefs(this%realization)
  
#if 0
  ! the problem here is that activity coefficients will be updated every time
  ! presolve is called, regardless of TS vs NI.  We need to split this out.
  if (this%realization%reaction%act_coef_update_frequency /= &
      ACT_COEF_FREQUENCY_OFF) then
    call RTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_TRUE,PETSC_TRUE)
!       The below is set within RTUpdateAuxVarsPatch() when 
!         PETSC_TRUE,PETSC_TRUE,* are passed
!       patch%aux%RT%auxvars_up_to_date = PETSC_TRUE 
  endif
#endif

  if (this%realization%reaction%use_log_formulation) then
    call VecCopy(this%realization%field%tran_xx, &
                 this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
    call VecLog(this%realization%field%tran_log_xx,ierr);CHKERRQ(ierr)
  endif
  
  call DataMediatorUpdate(this%realization%tran_data_mediator_list, &
                          this%realization%field%tran_mass_transfer, &
                          this%realization%option)
  
end subroutine PMRTPreSolve

! ************************************************************************** !

subroutine PMRTPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%PostSolve()')
#endif
  
end subroutine PMRTPostSolve

! ************************************************************************** !

subroutine PMRTFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/03/13
  ! 

  use Reactive_Transport_module, only : RTMaxChange
  use Variables_module, only : POROSITY
  use Material_module, only : MaterialGetAuxVarVecLoc
  use Material_Aux_class, only : POROSITY_MINERAL 
  use Global_module

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: time  
  PetscErrorCode :: ierr

  if (this%transient_porosity) then
    call VecCopy(this%realization%field%porosity_tpdt, &
                 this%realization%field%porosity_t,ierr);CHKERRQ(ierr)
    call RealizationUpdatePropertiesTS(this%realization)
    call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 POROSITY,POROSITY_MINERAL)
    call this%comm1%LocalToGlobal(this%realization%field%work_loc, &
                                  this%realization%field%porosity_tpdt)
  endif
  
  call RTMaxChange(this%realization,this%max_concentration_change, &
                   this%max_volfrac_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      this%max_concentration_change, &
      this%max_concentration_change/this%option%tran_dt
    if (this%realization%reaction%mineral%nkinmnrl > 0) then
      write(*,'("               dvfmx= ",1pe12.4," dvf/dt= ",1pe12.4, &
            &" [1/s]")') &
        this%max_volfrac_change, this%max_volfrac_change/this%option%tran_dt
    endif
  endif

#ifndef CLM_PFLOTRAN
!the following generates a large ascii file if coupled with CLM duo to long run steps
  if (this%option%print_file_flag) then  
    write(this%option%fid_out,&
            '("  --> max chng: dcmx= ",1pe12.4,"  dc/dt= ",1pe12.4, &
            &" [mol/s]")') &
      this%max_concentration_change, &
      this%max_concentration_change/this%option%tran_dt
    if (this%realization%reaction%mineral%nkinmnrl > 0) then
      write(this%option%fid_out, &
        '("               dvfmx= ",1pe12.4," dvf/dt= ",1pe12.4," [1/s]")') &
        this%max_volfrac_change, this%max_volfrac_change/this%option%tran_dt
    endif
  endif
#endif

end subroutine PMRTFinalizeTimestep

! ************************************************************************** !

function PMRTAcceptSolution(this)
  ! 
  ! PMRichardsAcceptSolution:
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_rt_type) :: this
  
  PetscBool :: PMRTAcceptSolution
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%AcceptSolution()')
#endif
  ! do nothing
  PMRTAcceptSolution = PETSC_TRUE
  
end function PMRTAcceptSolution

! ************************************************************************** !

subroutine PMRTUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac, &
                              time_step_max_growth_factor)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: dtt, uvf, dt_vf, dt_tfac, fac
  PetscInt :: ifac
  PetscReal, parameter :: pert = 1.d-20
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdateTimestep()')  
#endif
  
  if (this%volfrac_change_governor < 1.d0) then
    ! with volume fraction potentially scaling the time step.
    if (iacceleration > 0) then
      fac = 0.5d0
      if (num_newton_iterations >= iacceleration) then
        fac = 0.33d0
        uvf = 0.d0
      else
        uvf = this%volfrac_change_governor/(this%max_volfrac_change+pert)
      endif
      dtt = fac * dt * (1.d0 + uvf)
    else
      ifac = max(min(num_newton_iterations,size(tfac)),1)
      dt_tfac = tfac(ifac) * dt

      fac = 0.5d0
      uvf= this%volfrac_change_governor/(this%max_volfrac_change+pert)
      dt_vf = fac * dt * (1.d0 + uvf)

      dtt = min(dt_tfac,dt_vf)
    endif
  else
    ! original implementation
    dtt = dt
    if (num_newton_iterations <= iacceleration) then
      if (num_newton_iterations <= size(tfac)) then
        dtt = tfac(num_newton_iterations) * dt
      else
        dtt = 0.5d0 * dt
      endif
    else
      dtt = 0.5d0 * dt
    endif
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: see comment above under flow stepper
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  
end subroutine PMRTUpdateTimestep

! ************************************************************************** !

recursive subroutine PMRTFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%PMRTFinalizeRun()')
#endif
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMRTFinalizeRun

! ************************************************************************** !

subroutine PMRTResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTResidual

  implicit none
  
  class(pm_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Residual()')  
#endif
  
  call RTResidual(snes,xx,r,this%realization,ierr)

end subroutine PMRTResidual

! ************************************************************************** !

subroutine PMRTJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTJacobian

  implicit none
  
  class(pm_rt_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%Jacobian()')  
#endif

  call RTJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMRTJacobian

! ************************************************************************** !

subroutine PMRTCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! In the case of the log formulation, ensures that the update
  ! vector does not exceed a prescribed tolerance
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/16/09
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Option_module
  use Reaction_Aux_module

  implicit none
  
  class(pm_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: C_p(:)
  PetscReal, pointer :: dC_p(:)
  type(grid_type), pointer :: grid
  type(reaction_type), pointer :: reaction
  PetscReal :: ratio, min_ratio
  PetscReal, parameter :: min_allowable_scale = 1.d-10
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i, n
  
  grid => this%realization%patch%grid
  reaction => this%realization%reaction
  
  call VecGetArrayF90(dX,dC_p,ierr);CHKERRQ(ierr)

  if (reaction%use_log_formulation) then
    ! C and dC are actually lnC and dlnC
    dC_p = dsign(1.d0,dC_p)*min(dabs(dC_p),reaction%max_dlnC)
    ! at this point, it does not matter whether "changed" is set to true, 
    ! since it is not checkied in PETSc.  Thus, I don't want to spend 
    ! time checking for changes and performing an allreduce for log 
    ! formulation.
    if (Initialized(reaction%truncated_concentration)) then
      call VecGetArrayReadF90(X,C_p,ierr);CHKERRQ(ierr)
      dC_p = min(C_p-log(reaction%truncated_concentration),dC_p)
      call VecRestoreArrayReadF90(X,C_p,ierr);CHKERRQ(ierr)
    endif
  else
    call VecGetLocalSize(X,n,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X,C_p,ierr);CHKERRQ(ierr)
    
    if (Initialized(reaction%truncated_concentration)) then
      dC_p = min(dC_p,C_p-reaction%truncated_concentration)
    else
      ! C^p+1 = C^p - dC^p
      ! if dC is positive and abs(dC) larger than C
      ! we need to scale the update
      
      ! compute smallest ratio of C to dC
#if 0
      min_ratio = 1.d0/maxval(dC_p/C_p)
#else
      min_ratio = 1.d20 ! large number
      do i = 1, n
        if (C_p(i) <= dC_p(i)) then
          ratio = abs(C_p(i)/dC_p(i))
          if (ratio < min_ratio) min_ratio = ratio
        endif
      enddo
#endif
      ratio = min_ratio
    
      ! get global minimum
      call MPI_Allreduce(ratio,min_ratio,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                         MPI_MIN,this%realization%option%mycomm,ierr)
                       
      ! scale if necessary
      if (min_ratio < 1.d0) then
        if (min_ratio < this%realization%option%min_allowable_scale) then
          write(string,'(es10.3)') min_ratio
          string = 'The update of primary species concentration is being ' // &
            'scaled by a very small value (i.e. ' // &
            trim(adjustl(string)) // &
            ') to prevent negative concentrations.  This value is too ' // &
            'small and will likely cause the solver to mistakenly ' // &
            'converge based on the infinity norm of the update vector. ' // &
            'In this case, it is recommended that you use the ' // &
            'LOG_FORMULATION for chemistry or truncate concentrations ' // &
            '(TRUNCATE_CONCENTRATION <float> in CHEMISTRY block). ' // &
            'If that does not work, please send your input deck to ' // &
            'pflotran-dev@googlegroups.com.'
          this%realization%option%io_buffer = string
          call printErrMsg(this%realization%option)
        endif
        ! scale by 0.99 to make the update slightly smaller than the min_ratio
        dC_p = dC_p*min_ratio*0.99d0
        changed = PETSC_TRUE
      endif
    endif
    call VecRestoreArrayReadF90(X,C_p,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(dX,dC_p,ierr);CHKERRQ(ierr)

end subroutine PMRTCheckUpdatePre

! ************************************************************************** !

subroutine PMRTCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                               X1_changed,ierr)
  ! 
  ! Checks convergence after to update
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/14
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Secondary_Continuum_module, only : SecondaryRTUpdateIterate
  use Output_EKG_module
  use Reactive_Transport_Aux_module

  implicit none
  
  class(pm_rt_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch  
  PetscReal, pointer :: C0_p(:)
  PetscReal, pointer :: dC_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:)  
  PetscBool :: converged_due_to_rel_update
  PetscBool :: converged_due_to_residual
  PetscReal :: max_relative_change
  PetscReal :: max_scaled_residual
  PetscInt :: converged_flag
  PetscInt :: temp_int
  PetscReal :: max_relative_change_by_dof(this%option%ntrandof)
  PetscReal :: global_max_rel_change_by_dof(this%option%ntrandof)
  PetscMPIInt :: mpi_int
  PetscInt :: local_id, offset, idof, index
  PetscReal :: tempreal
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  converged_flag = 0
  if (this%check_post_convergence) then
    converged_due_to_rel_update = PETSC_FALSE
    converged_due_to_residual = PETSC_FALSE
    call VecGetArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    max_relative_change = maxval(dabs(dC_p(:)/C0_p(:)))
    call VecRestoreArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)
    max_scaled_residual = maxval(dabs(r_p(:)/accum_p(:)))
    call VecRestoreArrayReadF90(field%tran_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%tran_accum,accum_p,ierr);CHKERRQ(ierr)
    converged_due_to_rel_update = (Initialized(rt_itol_rel_update) .and. &
                                   max_relative_change < rt_itol_rel_update)
    converged_due_to_residual = (Initialized(rt_itol_scaled_res) .and. &
                                max_scaled_residual < rt_itol_scaled_res)
    if (converged_due_to_rel_update .or. converged_due_to_residual) then
      converged_flag = 1
    endif
  endif
  
  ! get global minimum
  call MPI_Allreduce(converged_flag,temp_int,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MIN,this%realization%option%mycomm,ierr)

  option%converged = PETSC_FALSE
  if (temp_int == 1) then
    option%converged = PETSC_TRUE
  endif
  
  if (option%use_mc) then  
    call SecondaryRTUpdateIterate(line_search,X0,dX,X1,dX_changed, &
                                  X1_changed,this%realization,ierr)
  endif
  
  if (this%print_ekg) then
    call VecGetArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    max_relative_change_by_dof = -1.d20
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%ntrandof
      do idof = 1, option%ntrandof
        index = idof + offset
        tempreal = dabs(dC_p(index)/C0_p(index))
        max_relative_change_by_dof(idof) = &
          max(max_relative_change_by_dof(idof),tempreal)
      enddo
    enddo
    call VecRestoreArrayReadF90(dX,dC_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,C0_p,ierr);CHKERRQ(ierr)
    mpi_int = option%ntrandof
    call MPI_Allreduce(MPI_IN_PLACE,max_relative_change_by_dof,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,this%option%mycomm,ierr)
    if (OptionPrintToFile(option)) then
100 format("REACTIVE TRANSPORT  NEWTON_ITERATION ",30es16.8)
      write(IUNIT_EKG,100) max_relative_change_by_dof(:)
    endif    
  endif

end subroutine PMRTCheckUpdatePost

! ************************************************************************** !

subroutine PMRTCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/15/17
  ! 
  use Convergence_module

  implicit none

  class(pm_rt_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  call ConvergenceTest(snes,it,xnorm,unorm,fnorm,reason, &
                       this%realization%patch%grid, &
                       this%option,this%solver,ierr)

end subroutine PMRTCheckConvergence

! ************************************************************************** !

subroutine PMRTTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTTimeCut

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%TimeCut()')
#endif
  
  this%option%tran_dt = this%option%dt
  if (this%option%nflowdof > 0 .and. .not. this%steady_flow) then
    call this%SetTranWeights()
  endif
  call RTTimeCut(this%realization)

end subroutine PMRTTimeCut

! ************************************************************************** !

subroutine PMRTUpdateSolution1(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module
  use Condition_module

  implicit none
  
  class(pm_rt_type) :: this
                                ! update kinetics
  call PMRTUpdateSolution2(this,PETSC_TRUE)
  
end subroutine PMRTUpdateSolution1

! ************************************************************************** !

subroutine PMRTUpdateSolution2(this, update_kinetics)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module
  use Condition_module
  use Integral_Flux_module

  implicit none
  
  class(pm_rt_type) :: this
  PetscBool :: update_kinetics
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%UpdateSolution()')
#endif
  
  ! begin from RealizationUpdate()
  call TranConditionUpdate(this%realization%transport_conditions, &
                           this%realization%option)
  if (associated(this%realization%uniform_velocity_dataset)) then
    call RealizUpdateUniformVelocity(this%realization)
  endif  
  ! end from RealizationUpdate()
  ! The update of status must be in this order!
  call RTUpdateEquilibriumState(this%realization)
  if (update_kinetics) &
    call RTUpdateKineticState(this%realization)
  
!TODO(geh): MassTransfer
!geh - moved to RTPreSolve()
!  call MassTransferUpdate(this%realization%rt_data_mediator_list, &
!                          this%realization%patch%grid, &
!                          this%realization%option)
  
  if (this%realization%option%compute_mass_balance_new) then
    call RTUpdateMassBalance(this%realization)
  endif
  if (this%option%transport%store_fluxes) then
    call IntegralFluxUpdate(this%realization%patch%integral_flux_list, &
                            this%realization%patch%internal_tran_fluxes, &
                            this%realization%patch%boundary_tran_fluxes, &
                            INTEGRATE_TRANSPORT,this%option)
  endif

end subroutine PMRTUpdateSolution2     

! ************************************************************************** !

subroutine PMRTUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Reactive_Transport_module, only : RTUpdateAuxVars
  
  implicit none
  
  class(pm_rt_type) :: this
                                      ! cells      bcs         act coefs
  call RTUpdateAuxVars(this%realization,PETSC_TRUE,PETSC_FALSE,PETSC_FALSE)

end subroutine PMRTUpdateAuxVars  

! ************************************************************************** !

subroutine PMRTMaxChange(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTMaxChange

  implicit none
  
  class(pm_rt_type) :: this
  
#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%MaxChange()')
#endif

  print *, 'PMRTMaxChange not implemented'
  stop
!  call RTMaxChange(this%realization)

end subroutine PMRTMaxChange

! ************************************************************************** !

subroutine PMRTComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTComputeMassBalance

  implicit none
  
  class(pm_rt_type) :: this
  PetscReal :: mass_balance_array(:)

#ifdef PM_RT_DEBUG  
  call printMsg(this%option,'PMRT%MassBalance()')
#endif

#ifndef SIMPLIFY 
  call RTComputeMassBalance(this%realization,mass_balance_array)
#endif

end subroutine PMRTComputeMassBalance

! ************************************************************************** !

subroutine SetTranWeights(this)
  ! 
  ! Sets the weights at t0 or t1 for transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/17/11; 04/03/13
  ! 

  use Option_module

  implicit none
  
  class(pm_rt_type) :: this

  PetscReal :: flow_dt
  PetscReal :: flow_t0
  PetscReal :: flow_t1

  ! option%tran_time is the time at beginning of transport step
  flow_t0 = this%realization%patch%aux%Global%time_t
  flow_t1 = this%realization%patch%aux%Global%time_tpdt
  flow_dt = flow_t1-flow_t0
  this%tran_weight_t0 = max(0.d0,(this%option%time-flow_t0)/flow_dt)
  this%tran_weight_t1 = min(1.d0, &
                            (this%option%time+this%option%tran_dt-flow_t0)/ &
                            flow_dt)

end subroutine SetTranWeights

! ************************************************************************** !

subroutine PMRTCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints flow reactive transport process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/29/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module

  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION, &
                               REACTION_AUXILIARY
  
  implicit none

  interface PetscBagGetData

! ************************************************************************** !

    subroutine PetscBagGetData(bag,header,ierr)
#include "petsc/finclude/petscsys.h"
      use petscsys
      import :: pm_rt_header_type
      implicit none
      PetscBag :: bag
      class(pm_rt_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData 

  PetscViewer :: viewer
  class(pm_rt_type) :: this
  PetscErrorCode :: ierr

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  PetscInt :: i

  class(pm_rt_header_type), pointer :: header
  type(pm_rt_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  
  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = PETSC_NULL_VEC

  bagsize = size(transfer(dummy_header,dummy_char))
  
  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%checkpoint_activity_coefs,0, &
                           "checkpoint_activity_coefs","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%ndof,0, &
                           "ndof","",ierr);CHKERRQ(ierr)
  if (associated(realization%reaction)) then
    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      header%checkpoint_activity_coefs = ONE_INTEGER
    else
      header%checkpoint_activity_coefs = ZERO_INTEGER
    endif
  else
    header%checkpoint_activity_coefs = ZERO_INTEGER
  endif
  !geh: %ndof should be pushed down to the base class, but this is not possible
  !     as long as option%ntrandof is used.
  header%ndof = option%ntrandof
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  
  if (option%ntrandof > 0) then
    call VecView(field%tran_xx, viewer, ierr);CHKERRQ(ierr)
    ! create a global vec for writing below 
    if (global_vec == PETSC_NULL_VEC) then
      call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
    endif
    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      ! allocated vector
      do i = 1, realization%reaction%naqcomp
        call RealizationGetVariable(realization,global_vec, &
                                   PRIMARY_ACTIVITY_COEF,i)
        call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
      enddo
      do i = 1, realization%reaction%neqcplx
        call RealizationGetVariable(realization,global_vec, &
                                   SECONDARY_ACTIVITY_COEF,i)
        call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
      enddo
    endif
    ! mineral volume fractions for kinetic minerals
    if (realization%reaction%mineral%nkinmnrl > 0) then
      do i = 1, realization%reaction%mineral%nkinmnrl
        call RealizationGetVariable(realization,global_vec, &
                                   MINERAL_VOLUME_FRACTION,i)
        call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
      enddo
    endif
    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction%nauxiliary> 0) then
      do i = 1, realization%reaction%nauxiliary
        call RealizationGetVariable(realization,global_vec, &
                                    REACTION_AUXILIARY,i)
        call VecView(global_vec,viewer,ierr);CHKERRQ(ierr)
      enddo
    endif
  endif

  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine PMRTCheckpointBinary

! ************************************************************************** !

subroutine PMRTRestartBinary(this,viewer)
  ! 
  ! Restarts flow reactive transport process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/29/13
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION, &
                               REACTION_AUXILIARY
  
  implicit none

  interface PetscBagGetData

! ************************************************************************** !

    subroutine PetscBagGetData(bag,header,ierr)
#include "petsc/finclude/petscsys.h"
      use petscsys
      import :: pm_rt_header_type
      implicit none
      PetscBag :: bag
      class(pm_rt_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData 

  PetscViewer :: viewer
  class(pm_rt_type) :: this
  PetscErrorCode :: ierr

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec, local_vec
  PetscInt :: i

  class(pm_rt_header_type), pointer :: header
  type(pm_rt_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  
  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid
  
  global_vec = PETSC_NULL_VEC
  local_vec = PETSC_NULL_VEC
  
  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(this%option%mycomm, bagsize, bag, ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag, header, ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%checkpoint_activity_coefs,0, &
                           "checkpoint_activity_coefs","",ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%ndof,0, &
                           "ndof","",ierr);CHKERRQ(ierr)
  call PetscBagLoad(viewer, bag, ierr);CHKERRQ(ierr)
  option%ntrandof = header%ndof
  
  call VecLoad(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
  call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                    field%tran_xx_loc,NTRANDOF)
  call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)

  if (global_vec == PETSC_NULL_VEC) then
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
  endif    
  if (header%checkpoint_activity_coefs == ONE_INTEGER) then
    call DiscretizationCreateVector(discretization,ONEDOF,local_vec, &
                                    LOCAL,option)
    do i = 1, realization%reaction%naqcomp
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                        local_vec,ONEDOF)
      call RealizationSetVariable(realization,local_vec,LOCAL, &
                                  PRIMARY_ACTIVITY_COEF,i)
    enddo
    do i = 1, realization%reaction%neqcplx
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call DiscretizationGlobalToLocal(discretization,global_vec, &
                                        local_vec,ONEDOF)
      call RealizationSetVariable(realization,local_vec,LOCAL, &
                                  SECONDARY_ACTIVITY_COEF,i)
    enddo
  endif
  ! mineral volume fractions for kinetic minerals
  if (realization%reaction%mineral%nkinmnrl > 0) then
    do i = 1, realization%reaction%mineral%nkinmnrl
      ! have to load the vecs no matter what
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      if (.not.option%transport%no_restart_mineral_vol_frac) then
        call RealizationSetVariable(realization,global_vec,GLOBAL, &
                                    MINERAL_VOLUME_FRACTION,i)
      endif
    enddo
  endif
  ! auxiliary data for reactions (e.g. cumulative mass)
  if (realization%reaction%nauxiliary> 0) then
    do i = 1, realization%reaction%nauxiliary
      call VecLoad(global_vec,viewer,ierr);CHKERRQ(ierr)
      call RealizationSetVariable(realization,global_vec,GLOBAL, &
                                  REACTION_AUXILIARY,i)
    enddo
  endif
    
  ! We are finished, so clean up.
  if (global_vec /= PETSC_NULL_VEC) then
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  endif
  if (local_vec /= PETSC_NULL_VEC) then
    call VecDestroy(local_vec,ierr);CHKERRQ(ierr)
  endif
  
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)
  
  if (realization%reaction%use_full_geochemistry) then
                                     ! cells     bcs        act coefs.
    call RTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  endif
  ! do not update kinetics.
  call PMRTUpdateSolution2(this,PETSC_FALSE)
  
end subroutine PMRTRestartBinary

! ************************************************************************** !

subroutine PMRTCheckpointHDF5(this, pm_grp_id)
  ! 
  ! Checkpoints flow reactive transport process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  ! 

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  class(pm_rt_type) :: this
  integer :: pm_grp_id
  type(option_type) :: option
  print *, 'PFLOTRAN must be compiled with HDF5 to ' // &
        'write HDF5 formatted checkpoint file. Darn.'
  stop
#else

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module

  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION, &
                               REACTION_AUXILIARY
  use hdf5
  use Checkpoint_module, only: CheckPointWriteIntDatasetHDF5
  use HDF5_module, only : HDF5WriteDataSetFromVec

  implicit none

  class(pm_rt_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

#if defined(SCORPIO_WRITE)
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: i
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  if (associated(realization%reaction)) then
    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then
      int_array(1) = ONE_INTEGER
    else
      int_array(1) = ZERO_INTEGER
    endif
  else
    int_array(1) = ZERO_INTEGER
  endif

  dataset_name = "Checkpoint_Activity_Coefs" // CHAR(0)
  call CheckPointWriteIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                     dims, start, length, stride, &
                                     int_array, option)

  dataset_name = "NDOF" // CHAR(0)
  int_array(1) = option%ntrandof
  call CheckPointWriteIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                     dims, start, length, stride, &
                                     int_array, option)

  !geh: %ndof should be pushed down to the base class, but this is not possible
  !     as long as option%ntrandof is used.

  if (option%ntrandof > 0) then

    call DiscretizationCreateVector(realization%discretization, NTRANDOF, &
                                     natural_vec, NATURAL, option)
    call DiscretizationGlobalToNatural(realization%discretization, &
                                       field%tran_xx, &
                                       natural_vec, NTRANDOF)
    dataset_name = "Primary_Variable" // CHAR(0)
    call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)
    call VecDestroy(natural_vec, ierr); CHKERRQ(ierr)

    ! create a global vec for writing below
    call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                      global_vec,GLOBAL,option)
    call DiscretizationCreateVector(realization%discretization, ONEDOF, &
                                     natural_vec, NATURAL, option)

    if (realization%reaction%checkpoint_activity_coefs .and. &
        realization%reaction%act_coef_update_frequency /= &
        ACT_COEF_FREQUENCY_OFF) then

      do i = 1, realization%reaction%naqcomp
        call RealizationGetVariable(realization,global_vec, &
                                    PRIMARY_ACTIVITY_COEF,i)
        call DiscretizationGlobalToNatural(realization%discretization, &
                                           global_vec, natural_vec, ONEDOF)
        write(dataset_name,*) i
        dataset_name = 'Aq_comp_' // trim(adjustl(dataset_name))
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)
      enddo

      do i = 1, realization%reaction%neqcplx
        call RealizationGetVariable(realization,global_vec, &
                                   SECONDARY_ACTIVITY_COEF,i)
        call DiscretizationGlobalToNatural(realization%discretization, &
                                           global_vec, natural_vec, ONEDOF)
        write(dataset_name,*) i
        dataset_name = 'Eq_cplx_' // trim(adjustl(dataset_name))
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)
      enddo
    endif

    ! mineral volume fractions for kinetic minerals
    if (realization%reaction%mineral%nkinmnrl > 0) then
      do i = 1, realization%reaction%mineral%nkinmnrl
        call RealizationGetVariable(realization,global_vec, &
                                   MINERAL_VOLUME_FRACTION,i)
        call DiscretizationGlobalToNatural(realization%discretization, &
                                           global_vec,natural_vec,ONEDOF)
        write(dataset_name,*) i
        dataset_name = 'Kinetic_mineral_' // trim(adjustl(dataset_name))
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)
      enddo
    endif

    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction%nauxiliary> 0) then
      do i = 1, realization%reaction%nauxiliary
        call RealizationGetVariable(realization,global_vec, &
                                    REACTION_AUXILIARY,i)
        call DiscretizationGlobalToNatural(realization%discretization, &
                                           global_vec, natural_vec, ONEDOF)
        write(dataset_name,*) i
        dataset_name = 'Reaction_auxiliary_' // trim(adjustl(dataset_name))
        call HDF5WriteDataSetFromVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)
      enddo
    endif
    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

   endif
#endif

end subroutine PMRTCheckpointHDF5

! ************************************************************************** !

subroutine PMRTRestartHDF5(this, pm_grp_id)
  ! 
  ! Checkpoints flow reactive transport process model
  ! 
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  ! 

#if  !defined(PETSC_HAVE_HDF5)
  implicit none
  class(pm_rt_type) :: this
  integer :: pm_grp_id
  type(option_type) :: option
  print *, 'PFLOTRAN must be compiled with HDF5 to ' // &
        'write HDF5 formatted checkpoint file. Darn.'
  stop
#else

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Field_module
  use Discretization_module
  use Grid_module
  use Reactive_Transport_module, only : RTUpdateAuxVars
  use Reaction_Aux_module, only : ACT_COEF_FREQUENCY_OFF
  use Variables_module, only : PRIMARY_ACTIVITY_COEF, &
                               SECONDARY_ACTIVITY_COEF, &
                               MINERAL_VOLUME_FRACTION, &
                               REACTION_AUXILIARY
  use hdf5
  use Checkpoint_module, only: CheckPointReadIntDatasetHDF5
  use HDF5_module, only : HDF5ReadDataSetInVec

  implicit none

  class(pm_rt_type) :: this
#if defined(SCORPIO_WRITE)
  integer :: pm_grp_id
#else
  integer(HID_T) :: pm_grp_id
#endif

#if defined(SCORPIO_WRITE)
  integer, pointer :: dims(:)
  integer, pointer :: start(:)
  integer, pointer :: stride(:)
  integer, pointer :: length(:)
#else
  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
#endif

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers 
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  Vec :: local_vec
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: i
  PetscInt :: checkpoint_activity_coefs
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  grid => realization%patch%grid

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Checkpoint_Activity_Coefs" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                    dims, start, length, stride, &
                                    int_array, option)
  checkpoint_activity_coefs = int_array(1)
  
  dataset_name = "NDOF" // CHAR(0)
  int_array(1) = option%ntrandof
  call CheckPointReadIntDatasetHDF5(pm_grp_id, dataset_name, dataset_rank, &
                                    dims, start, length, stride, &
                                    int_array, option)
  option%ntrandof = int_array(1)
  
  !geh: %ndof should be pushed down to the base class, but this is not possible
  !     as long as option%ntrandof is used.

  if (option%ntrandof > 0) then

    call DiscretizationCreateVector(discretization, NTRANDOF, &
                                     natural_vec, NATURAL, option)
    dataset_name = "Primary_Variable" // CHAR(0)
    call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
                             pm_grp_id, H5T_NATIVE_DOUBLE)
    call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                       field%tran_xx, &
                                       NTRANDOF)
    call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                    field%tran_xx_loc,NTRANDOF)
    call VecCopy(field%tran_xx,field%tran_yy,ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec, ierr); CHKERRQ(ierr)

    ! create a global vec for reading
    call DiscretizationCreateVector(discretization,ONEDOF, &
                                    global_vec,GLOBAL,option)
    call DiscretizationCreateVector(discretization, ONEDOF, &
                                    natural_vec, NATURAL, option)
    call DiscretizationCreateVector(discretization,ONEDOF,local_vec, &
                                    LOCAL,option)

    if (checkpoint_activity_coefs == ONE_INTEGER) then

      do i = 1, realization%reaction%naqcomp
        write(dataset_name,*) i
        dataset_name = 'Aq_comp_' // trim(adjustl(dataset_name))
        call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)

        call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                           global_vec, ONEDOF)
        call DiscretizationGlobalToLocal(discretization, global_vec, &
                                         local_vec, ONEDOF)
        call RealizationSetVariable(realization, local_vec, LOCAL, &
                                    PRIMARY_ACTIVITY_COEF,i)
      enddo

      do i = 1, realization%reaction%neqcplx
        write(dataset_name,*) i
        dataset_name = 'Eq_cplx_' // trim(adjustl(dataset_name))
        call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)

        call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                           global_vec, ONEDOF)
        call DiscretizationGlobalToLocal(discretization, global_vec, &
                                         local_vec, ONEDOF)
        call RealizationSetVariable(realization, local_vec, LOCAL, &
                                   SECONDARY_ACTIVITY_COEF, i)
      enddo
    endif

    ! mineral volume fractions for kinetic minerals
    if (realization%reaction%mineral%nkinmnrl > 0) then
      do i = 1, realization%reaction%mineral%nkinmnrl
        write(dataset_name,*) i
        dataset_name = 'Kinetic_mineral_' // trim(adjustl(dataset_name))
        call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)

        call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                           global_vec, ONEDOF)
        call DiscretizationGlobalToLocal(discretization, global_vec, &
                                         local_vec, ONEDOF)
        call RealizationSetVariable(realization, local_vec, LOCAL, &
                                   MINERAL_VOLUME_FRACTION,i)
      enddo
    endif

    ! auxiliary data for reactions (e.g. cumulative mass)
    if (realization%reaction%nauxiliary> 0) then
      do i = 1, realization%reaction%nauxiliary
        write(dataset_name,*) i
        dataset_name = 'Reaction_auxiliary_' // trim(adjustl(dataset_name))
        call HDF5ReadDataSetInVec(dataset_name, option, natural_vec, &
           pm_grp_id, H5T_NATIVE_DOUBLE)

        call DiscretizationNaturalToGlobal(discretization, natural_vec, &
                                           global_vec, ONEDOF)
        call DiscretizationGlobalToLocal(discretization, global_vec, &
                                         local_vec, ONEDOF)
        call RealizationSetVariable(realization,local_vec,LOCAL, &
                                    REACTION_AUXILIARY,i)
      enddo
    endif

    call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

  endif

  if (realization%reaction%use_full_geochemistry) then
                                     ! cells     bcs        act coefs.
    call RTUpdateAuxVars(realization,PETSC_FALSE,PETSC_TRUE,PETSC_FALSE)
  endif
  ! do not update kinetics.
  call PMRTUpdateSolution2(this,PETSC_FALSE)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)

#endif

end subroutine PMRTRestartHDF5

! ************************************************************************** !

subroutine PMRTInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_rt_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMRTInputRecord

! ************************************************************************** !

subroutine PMRTDestroy(this)
  ! 
  ! Destroys RT process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Reactive_Transport_module, only : RTDestroy

  implicit none
  
  class(pm_rt_type) :: this

  call RTDestroy(this%realization)
  ! destroyed in realization
  nullify(this%comm1)
  nullify(this%option)
  nullify(this%output_option)
  call this%commN%Destroy()
  if (associated(this%commN)) deallocate(this%commN)
  nullify(this%commN)  

end subroutine PMRTDestroy
  
end module PM_RT_class
