module PM_General_class

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

  type, public, extends(pm_subsurface_flow_type) :: pm_general_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
  contains
    procedure, public :: Read => PMGeneralRead
    procedure, public :: InitializeRun => PMGeneralInitializeRun
    procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    procedure, public :: Residual => PMGeneralResidual
    procedure, public :: Jacobian => PMGeneralJacobian
    procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMGeneralPreSolve
    procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    procedure, public :: UpdateAuxVars => PMGeneralUpdateAuxVars
    procedure, public :: MaxChange => PMGeneralMaxChange
    procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: InputRecord => PMGeneralInputRecord
    procedure, public :: CheckpointBinary => PMGeneralCheckpointBinary
    procedure, public :: RestartBinary => PMGeneralRestartBinary
    procedure, public :: Destroy => PMGeneralDestroy
  end type pm_general_type
  
  public :: PMGeneralCreate
  
contains

! ************************************************************************** !

function PMGeneralCreate()
  ! 
  ! Creates General process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_general_type), pointer :: PMGeneralCreate

  class(pm_general_type), pointer :: general_pm
  
#ifdef PM_GENERAL_DEBUG  
  print *, 'PMGeneralCreate()'
#endif  

  allocate(general_pm)
  allocate(general_pm%max_change_ivar(6))
  general_pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                                LIQUID_MOLE_FRACTION, TEMPERATURE, &
                                GAS_SATURATION]
  allocate(general_pm%max_change_isubvar(6))
                                   ! UNINITIALIZED_INTEGER avoids zeroing of 
                                   ! pressures not represented in phase
                                       ! 2 = air in xmol(air,liquid)
  general_pm%max_change_isubvar = [0,0,0,2,0,0]
  
  call PMSubsurfaceFlowCreate(general_pm)
  general_pm%name = 'General Multiphase Flow'

  PMGeneralCreate => general_pm
  
end function PMGeneralCreate

! ************************************************************************** !

subroutine PMGeneralRead(this,input)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/15
  !
  use General_module
  use General_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_general_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'General Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found,option)    
    if (found) cycle
    
    select case(trim(keyword))
      case('ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,general_itol_scaled_res)
        call InputDefaultMsg(input,option,'general_itol_scaled_res')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,general_itol_rel_update)
        call InputDefaultMsg(input,option,'general_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE        
      case('TOUGH2_ITOL_SCALED_RESIDUAL')
        ! since general_tough2_itol_scaled_res_e1 is an array, we must read
        ! the tolerance into a double and copy it to the array.
        tempreal = UNINITIALIZED_DOUBLE
        call InputReadDouble(input,option,tempreal)
        ! tempreal will remain uninitialized if the read fails.
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e1')
        if (Initialized(tempreal)) then
          general_tough2_itol_scaled_res_e1 = tempreal
        endif
        call InputReadDouble(input,option,general_tough2_itol_scaled_res_e2)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e2')
        general_tough2_conv_criteria = PETSC_TRUE
        this%check_post_convergence = PETSC_TRUE
      case('T2_ITOL_SCALED_RESIDUAL_TEMP')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option, &
                           'tough_itol_scaled_residual_e1 for temperature', &
                           error_string)
        general_tough2_itol_scaled_res_e1(3,:) = tempreal
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,'window epsilon',error_string)
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,'gas component formula wt.', &
                           error_string)
      case('TWO_PHASE_ENERGY_DOF')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'two_phase_energy_dof',error_string)
        call GeneralAuxSetEnergyDOF(word,option)
      case('ISOTHERMAL')
        general_isothermal = PETSC_TRUE
      case('NO_AIR')
        general_no_air = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,general_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,general_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,general_damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
#if 0        
      case('GOVERN_MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%dPmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable pressure change', &
                           error_string)
      case('GOVERN_MAXIMUM_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,this%dTmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable temperature change', &
                           error_string)
      case('GOVERN_MAXIMUM_SATURATION_CHANGE')
        call InputReadDouble(input,option,this%dSmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable saturation change', &
                           error_string)
      case('GOVERN_MAXIMUM_MOLE_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%dXmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable mole fraction change', &
                           error_string)
#endif
      case('DEBUG_CELL')
        call InputReadInt(input,option,general_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id',error_string)
      case('NO_TEMP_DEPENDENT_DIFFUSION')
        general_temp_dep_gas_air_diff = PETSC_FALSE
      case('DIFFUSE_XMASS')
        general_diffuse_xmol = PETSC_FALSE
      case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_TRUE
      case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_FALSE
      case('ANALYTICAL_DERIVATIVES')
        general_analytical_derivatives = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'GENERAL Mode',option)
    end select
    
  enddo  
  
  if (general_isothermal .and. &
      general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
    option%io_buffer = 'Isothermal GENERAL mode may only be run with ' // &
                       'temperature as the two phase energy dof.'
    call printErrMsg(option)
  endif

end subroutine PMGeneralRead

! ************************************************************************** !

recursive subroutine PMGeneralInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 

  use Realization_Base_class
  
  implicit none
  
  class(pm_general_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 6
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMGeneralInitializeRun

! ************************************************************************** !

subroutine PMGeneralInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_general_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)
                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," GENERAL FLOW ",64("="))')
  endif
  
  call GeneralInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)                                 
  
end subroutine PMGeneralInitializeTimestep

! ************************************************************************** !

subroutine PMGeneralPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPreSolve

! ************************************************************************** !

subroutine PMGeneralPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPostSolve

! ************************************************************************** !

subroutine PMGeneralUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  type(field_type), pointer :: field
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
    ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
    us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
    umin = min(up,ut,ux,us)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

  if (Initialized(this%cfl_governor)) then
    ! Since saturations are not stored in global_auxvar for general mode, we
    ! must copy them over for the CFL check
    ! liquid saturation
    field => this%realization%field
    call RealizationGetVariable(this%realization,field%work, &
                                LIQUID_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               LIQUID_SATURATION,TIME_NULL)
    call RealizationGetVariable(this%realization,field%work, &
                                GAS_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               GAS_SATURATION,TIME_NULL)
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  endif

end subroutine PMGeneralUpdateTimestep

! ************************************************************************** !

subroutine PMGeneralResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralResidual

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call GeneralResidual(snes,xx,r,this%realization,ierr)

end subroutine PMGeneralResidual

! ************************************************************************** !

subroutine PMGeneralJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralJacobian

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call GeneralJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMGeneralJacobian

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use General_Aux_module
  use Global_Aux_module
  
  implicit none
  
  class(pm_general_type) :: this
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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
#ifdef DEBUG_GENERAL_INFO
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  character(len=MAXWORDLENGTH) :: cell_id_word
  PetscInt, parameter :: max_cell_id = 10
  PetscInt :: cell_id, cell_locator(0:max_cell_id)
  PetscInt :: i, ii
#endif
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: liquid_pressure_index, gas_pressure_index, air_pressure_index
  PetscInt :: saturation_index, temperature_index, xmol_index
  PetscInt :: lid, gid, apid, cpid, vpid, spid
  PetscInt :: pgas_index
  PetscReal :: liquid_pressure0, liquid_pressure1, del_liquid_pressure
  PetscReal :: gas_pressure0, gas_pressure1, del_gas_pressure
  PetscReal :: air_pressure0, air_pressure1, del_air_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: xmol0, xmol1, del_xmol
  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: min_pressure
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  gen_auxvars => this%realization%patch%aux%General%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  patch => this%realization%patch

  spid = option%saturation_pressure_id
  apid = option%air_pressure_id

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  ! truncation
  ! air mole fraction in liquid phase must be truncated.  we do not use scaling
  ! here because of the very small values.  just truncation.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    select case(global_auxvars(ghosted_id)%istate)
      case(LIQUID_STATE)
        xmol_index = offset + GENERAL_LIQUID_STATE_X_MOLE_DOF
        if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
          ! we use 1.d-10 since cancelation can occur with smaller values
          dX_p(xmol_index) = X_p(xmol_index)
        endif
      case(TWO_PHASE_STATE)
        pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
        if (X_p(pgas_index) - dX_p(pgas_index) < &
            gen_auxvars(ZERO_INTEGER,ghosted_id)% &
              pres(option%saturation_pressure_id)) then
          dX_p(pgas_index) = X_p(pgas_index) - &
            gen_auxvars(ZERO_INTEGER,ghosted_id)% &
              pres(option%saturation_pressure_id)
        endif
    end select
  enddo

  scale = initial_scale
  if (general_max_it_before_damping > 0 .and. &
      newton_iteration > general_max_it_before_damping) then
    scale = general_damping_factor
  endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!#define TRUNCATE_LIQUID_PRESSURE
!! TRUNCATE_GAS/AIR_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_GAS_PRESSURE
!#define TRUNCATE_AIR_PRESSURE

#ifdef DEBUG_GENERAL_INFO
  cell_locator = 0
#endif

  ! scaling
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof
    temp_scale = 1.d0
#ifdef DEBUG_GENERAL_INFO
    cell_id = grid%nG2A(ghosted_id)
    write(cell_id_word,*) cell_id
    cell_id_word = '(Cell ' // trim(adjustl(cell_id_word)) // '): '
#endif
    select case(global_auxvars(ghosted_id)%istate)
      case(LIQUID_STATE)
        liquid_pressure_index  = offset + GENERAL_LIQUID_PRESSURE_DOF
        temperature_index  = offset + GENERAL_ENERGY_DOF
        dX_p(liquid_pressure_index) = dX_p(liquid_pressure_index) * &
                                      general_pressure_scale
        temp_scale = 1.d0
        del_liquid_pressure = dX_p(liquid_pressure_index)
        liquid_pressure0 = X_p(liquid_pressure_index)
        liquid_pressure1 = liquid_pressure0 - del_liquid_pressure
        del_temperature = dX_p(temperature_index)
        temperature0 = X_p(temperature_index)
        temperature1 = temperature0 - del_temperature
#ifdef LIMIT_MAX_PRESSURE_CHANGE
        if (dabs(del_liquid_pressure) > general_max_pressure_change) then
          temp_real = dabs(general_max_pressure_change/del_liquid_pressure)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Liquid pressure change scaled to truncate at max_pressure_change: '
          call printMsg(option,string)
          write(string2,*) liquid_pressure0
          string = '  Liquid Pressure 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) liquid_pressure1
          string = '  Liquid Pressure 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_liquid_pressure
          string = 'Liquid Pressure change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif 
!LIMIT_MAX_PRESSURE_CHANGE
#ifdef TRUNCATE_LIQUID_PRESSURE
        ! truncate liquid pressure change to prevent liquid pressure from 
        ! dropping below the air pressure while in the liquid state
        min_pressure = gen_auxvars(ZERO_INTEGER,ghosted_id)%pres(apid) + &
                       gen_auxvars(ZERO_INTEGER,ghosted_id)%pres(spid)
        if (liquid_pressure1 < min_pressure) then
          temp_real = tolerance * (liquid_pressure0 - min_pressure)
          if (dabs(del_liquid_pressure) > 1.d-40) then
            temp_real = dabs(temp_real / del_liquid_pressure)
          else
            option%io_buffer = 'Something is seriously wrong in ' // &
              'GeneralCheckUpdatePre(liquid<min).  Contact Glenn through ' // &
              'pflotran-dev@googlegroups.com.'
            call printErrMsg(option)
          endif
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Liquid pressure change scaled to prevent liquid ' // &
            'pressure from dropping below air pressure: '
          call printMsg(option,string)
          write(string2,*) liquid_pressure0
          string = '  Liquid pressure 0: ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) liquid_pressure1
          string = '  Liquid pressure 1: ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_liquid_pressure
          string = '  pressure change  : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif 
!TRUNCATE_LIQUID_PRESSURE  
#ifdef LIMIT_MAX_TEMPERATURE_CHANGE
        if (dabs(del_temperature) > max_temperature_change) then
          temp_real = dabs(max_temperature_change/del_temperature)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Temperature change scaled to truncate at max_temperature_change: '
          call printMsg(option,string)
          write(string2,*) temperature0
          string = '  Temperature 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temperature1
          string = '  Temperature 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_temperature
          string = 'Temperature change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif 
!LIMIT_MAX_TEMPERATURE_CHANGE        
      case(TWO_PHASE_STATE)
        gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
!        air_pressure_index = offset + 2
        saturation_index = offset + GENERAL_GAS_SATURATION_DOF
        temperature_index  = offset + GENERAL_ENERGY_DOF
        dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                   general_pressure_scale
        if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then                                   
          air_pressure_index = offset + GENERAL_ENERGY_DOF
          dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                     general_pressure_scale
          del_air_pressure = dX_p(air_pressure_index)
          air_pressure0 = X_p(air_pressure_index)
          air_pressure1 = air_pressure0 - del_air_pressure
        endif
        temp_scale = 1.d0
        del_gas_pressure = dX_p(gas_pressure_index)
        gas_pressure0 = X_p(gas_pressure_index)
        gas_pressure1 = gas_pressure0 - del_gas_pressure
        del_saturation = dX_p(saturation_index)
        saturation0 = X_p(saturation_index)
        saturation1 = saturation0 - del_saturation
#ifdef LIMIT_MAX_PRESSURE_CHANGE
        if (dabs(del_gas_pressure) > general_max_pressure_change) then
          temp_real = dabs(general_max_pressure_change/del_gas_pressure)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas pressure change scaled to truncate at max_pressure_change: '
          call printMsg(option,string)
          write(string2,*) gas_pressure0
          string = '  Gas Pressure 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) gas_pressure1
          string = '  Gas Pressure 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_gas_pressure
          string = 'Gas Pressure change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif
#ifdef TRUNCATE_GAS_PRESSURE
        if (gas_pressure1 <= 0.d0) then
          if (dabs(del_gas_pressure) > 1.d-40) then
            temp_real = tolerance * dabs(gas_pressure0 / del_gas_pressure)
#ifdef DEBUG_GENERAL_INFO
            if (cell_locator(0) < max_cell_id) then
              cell_locator(0) = cell_locator(0) + 1
              cell_locator(cell_locator(0)) = ghosted_id
            endif
            string = trim(cell_id_word) // &
              'Gas pressure change scaled to prevent gas ' // &
              'pressure from dropping below zero: '
            call printMsg(option,string)
            write(string2,*) gas_pressure0
            string = '  Gas pressure 0   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) gas_pressure1
            string = '  Gas pressure 1   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) -1.d0*del_gas_pressure
            string = '  pressure change  : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) temp_real
            string = '          scaling  : ' // adjustl(string2)
            call printMsg(option,string)
#endif
            temp_scale = min(temp_scale,temp_real)
          endif
        endif
#endif 
!TRUNCATE_GAS_PRESSURE
#ifdef TRUNCATE_AIR_PRESSURE
        if (air_pressure1 <= 0.d0) then
          if (dabs(del_air_pressure) > 1.d-40) then
            temp_real = tolerance * dabs(air_pressure0 / del_air_pressure)
#ifdef DEBUG_GENERAL_INFO
            if (cell_locator(0) < max_cell_id) then
              cell_locator(0) = cell_locator(0) + 1
              cell_locator(cell_locator(0)) = ghosted_id
            endif
            string = trim(cell_id_word) // &
              'Air pressure change scaled to prevent air ' // &
              'pressure from dropping below zero: '
            call printMsg(option,string)
            write(string2,*) air_pressure0
            string = '  Air pressure 0   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) air_pressure1
            string = '  Air pressure 1   : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) -1.d0*del_air_pressure
            string = '  pressure change  : ' // adjustl(string2)
            call printMsg(option,string)
            write(string2,*) temp_real
            string = '          scaling  : ' // adjustl(string2)
            call printMsg(option,string)
#endif
            temp_scale = min(temp_scale,temp_real)
          endif
        endif
#endif 
!TRUNCATE_AIR_PRESSURE
#if defined(TRUNCATE_GAS_PRESSURE) && defined(TRUNCATE_AIR_PRESSURE)
        ! have to factor in scaled update from previous conditionals
        gas_pressure1 = gas_pressure0 - temp_scale * del_gas_pressure
        air_pressure1 = air_pressure0 - temp_scale * del_air_pressure
        if (gas_pressure1 <= air_pressure1) then
!          temp_real = (air_pressure0 - gas_pressure0) / &
!                      (temp_scale * (del_air_pressure - del_gas_pressure))
!          temp_real = temp_real * tolerance * temp_scale
          ! refactored to prevent divide by 0
          temp_real = temp_scale * (del_air_pressure - del_gas_pressure)
          if (temp_real > 0.d0) then
            temp_real = (air_pressure0 - gas_pressure0) / temp_real
            temp_real = temp_real * tolerance * temp_scale
          endif
          if (temp_real <= 0.d0) then
            ! add info on pressures here
            option%io_buffer = 'Something is seriously wrong in ' // &
              'GeneralCheckUpdatePre(gas<=air).  Contact Glenn through ' // &
              'pflotran-dev@googlegroups.com.'
            call printErrMsg(option)
          endif
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas/Air pressure change scaled again to prevent gas ' // &
            'pressure from dropping below air pressure: '
          call printMsg(option,string)
          write(string2,*) gas_pressure0
          string = '  Gas pressure 0       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) gas_pressure1
          string = '  Gas pressure 1       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*temp_real*del_gas_pressure
          string = '  Gas pressure change  : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) air_pressure0
          string = '  Air pressure 0       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) air_pressure1
          string = '  Air pressure 1       : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*temp_real*del_air_pressure
          string = '  Air pressure change  : ' // adjustl(string2)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif 
!TRUNCATE_GAS_PRESSURE && TRUNCATE_AIR_PRESSURE
#ifdef LIMIT_MAX_SATURATION_CHANGE
        if (dabs(del_saturation) > max_saturation_change) then
          temp_real = dabs(max_saturation_change/del_saturation)
#ifdef DEBUG_GENERAL_INFO
          if (cell_locator(0) < max_cell_id) then
            cell_locator(0) = cell_locator(0) + 1
            cell_locator(cell_locator(0)) = ghosted_id
          endif
          string = trim(cell_id_word) // &
            'Gas saturation change scaled to truncate at ' // &
            'max_saturation_change: '
          call printMsg(option,string)
          write(string2,*) saturation0
          string = '  Saturation 0    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) saturation1
          string = '  Saturation 1    : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) -1.d0*del_saturation
          string = 'Saturation change : ' // adjustl(string2)
          call printMsg(option,string)
          write(string2,*) temp_real
          string = '          scaling  : ' // adjustl(string2)
          call printMsg(option,string)
#endif
          temp_scale = min(temp_scale,temp_real)
        endif
#endif 
!LIMIT_MAX_SATURATION_CHANGE        
      case(GAS_STATE) 
        gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
        air_pressure_index = offset + GENERAL_GAS_STATE_AIR_PRESSURE_DOF
        dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                   general_pressure_scale
        dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                   general_pressure_scale
    end select
    scale = min(scale,temp_scale) 
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  if (scale < 0.9999d0) then
#ifdef DEBUG_GENERAL_INFO
    string  = '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    call printMsg(option,string)
    write(string2,*) scale, (grid%nG2A(cell_locator(i)),i=1,cell_locator(0))
    string = 'Final scaling: : ' // adjustl(string2)
    call printMsg(option,string)
    do i = 1, cell_locator(0)
      ghosted_id = cell_locator(i)
      offset = (ghosted_id-1)*option%nflowdof
      write(string2,*) grid%nG2A(ghosted_id)
      string = 'Cell ' // trim(adjustl(string2))
      write(string2,*) global_auxvars(ghosted_id)%istate
      string = trim(string) // ' (State = ' // trim(adjustl(string2)) // ') '
      call printMsg(option,string)
      ! for some reason cannot perform array operation on dX_p(:)
      write(string2,*) (X_p(offset+ii),ii=1,3)
      string = '   Orig. Solution: ' // trim(adjustl(string2))
      call printMsg(option,string)
      write(string2,*) (X_p(offset+ii)-dX_p(offset+ii),ii=1,3)
      string = '  Solution before: ' // trim(adjustl(string2))
      call printMsg(option,string)
      write(string2,*) (X_p(offset+ii)-scale*dX_p(offset+ii),ii=1,3)
      string = '   Solution after: ' // trim(adjustl(string2))
      call printMsg(option,string)
    enddo
    string  = '++++++++++++++++++++++++++++++++++++++++++++++++++++++'
    call printMsg(option,string)
#endif
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMGeneralCheckUpdatePre

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  use Output_EKG_module
  
  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(general_auxvar_type), pointer :: general_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(material_parameter_type), pointer :: material_parameter
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset , ival, idof
#ifdef DEBUG_GENERAL_INFO
  PetscInt :: icell_max_rel_update(3), icell_max_scaled_residual(3)
  PetscInt :: istate_max_rel_update(3), istate_max_scaled_residual(3)
  character(len=2) :: state_char
#endif
  PetscReal :: dX_X0, R_A, R
  PetscReal :: inf_norm_rel_update(3,3), global_inf_norm_rel_update(3,3)
  PetscReal :: inf_norm_scaled_residual(3,3), global_inf_norm_scaled_residual(3,3)
  PetscReal :: inf_norm_update(3,3), global_inf_norm_update(3,3)
  PetscReal :: inf_norm_residual(3,3), global_inf_norm_residual(3,3)
  PetscReal :: two_norm_residual(3,3), global_two_norm_residual(3,3)
  PetscReal, parameter :: inf_pres_tol = 1.d-1
  PetscReal, parameter :: inf_temp_tol = 1.d-5
  PetscReal, parameter :: inf_sat_tol = 1.d-6
  PetscReal, parameter :: inf_xmol_tol = 1.d-6
  !geh: note the scaling by 0.d0 several lines down which prevent false 
  !     convergence
  PetscReal, parameter :: inf_norm_update_tol(3,3) = &
    reshape([inf_pres_tol,inf_xmol_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_sat_tol], &
            shape(inf_norm_update_tol)) * &
            0.d0
  PetscReal :: temp(3,12), global_temp(3,12)
  PetscMPIInt :: mpi_int
  PetscBool :: converged_abs_update
  PetscBool :: converged_rel_update
  PetscBool :: converged_scaled_residual
  PetscInt :: istate
  PetscReal :: t_over_v
  PetscReal :: two_norm
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  general_auxvars => patch%aux%General%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  option%converged = PETSC_FALSE
  if (this%check_post_convergence) then
    call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
#ifdef DEBUG_GENERAL_INFO
    icell_max_rel_update = 0
    istate_max_rel_update = 0
    icell_max_scaled_residual = 0
    istate_max_scaled_residual = 0
#endif
    inf_norm_update(:,:) = -1.d20
    inf_norm_rel_update(:,:) = -1.d20
    inf_norm_scaled_residual(:,:) = -1.d20
    inf_norm_residual(:,:) = -1.d20
    two_norm_residual(:,:) = 0.d0
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      istate = global_auxvars(ghosted_id)%istate
      do idof = 1, option%nflowdof
        ival = offset+idof
        R = r_p(ival)
#ifdef DEBUG_GENERAL_INFO
        two_norm_residual(idof,istate) = two_norm_residual(idof,istate) + R*R
#endif
        inf_norm_residual(idof,istate) = max(inf_norm_residual(idof,istate), &
                                             dabs(R))
        if (general_tough2_conv_criteria) then
          !geh: scale by t_over_v to match TOUGH2 residual units. see equation
          !     B.5 of TOUGH2 user manual (LBNL-43134)
          t_over_v = option%flow_dt/material_auxvars(ghosted_id)%volume
          if (accum_p2(ival)*t_over_v < general_tough2_itol_scaled_res_e2) then
            R_A = dabs(R*t_over_v)
          else
            R_A = dabs(R/accum_p2(ival))
          endif
        else
          R_A = dabs(R/accum_p(ival))
        endif
        dX_X0 = dabs(dX_p(ival)/X0_p(ival))
        inf_norm_update(idof,istate) = max(inf_norm_update(idof,istate), &
                                           dabs(dX_p(ival)))
        if (inf_norm_rel_update(idof,istate) < dX_X0) then
#ifdef DEBUG_GENERAL_INFO
          if (maxval(inf_norm_rel_update(idof,:)) < dX_X0) then
            icell_max_rel_update(idof) = grid%nG2A(ghosted_id)
            istate_max_rel_update(idof) = global_auxvars(ghosted_id)%istate
          endif
#endif
          inf_norm_rel_update(idof,istate) = dX_X0
        endif
        if (inf_norm_scaled_residual(idof,istate) < R_A) then
#ifdef DEBUG_GENERAL_INFO
          if (maxval(inf_norm_scaled_residual(idof,:)) < R_A) then
            icell_max_scaled_residual(idof) = grid%nG2A(ghosted_id)
            istate_max_scaled_residual(idof) = global_auxvars(ghosted_id)%istate
          endif
#endif
          inf_norm_scaled_residual(idof,istate) = R_A
        endif
      enddo
    enddo
    temp(1:3,1:3) = inf_norm_update(:,:)
    temp(1:3,4:6) = inf_norm_rel_update(:,:)
    temp(1:3,7:9) = inf_norm_scaled_residual(:,:)
    temp(1:3,10:12) = inf_norm_residual(:,:)
    mpi_int = 36
    call MPI_Allreduce(temp,global_temp,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    global_inf_norm_update(:,:) = global_temp(1:3,1:3)
    global_inf_norm_rel_update(:,:) = global_temp(1:3,4:6)
    global_inf_norm_scaled_residual(:,:) = global_temp(1:3,7:9)
    global_inf_norm_residual(:,:) = global_temp(1:3,10:12)

#ifdef DEBUG_GENERAL_INFO
    mpi_int = 9
    call MPI_Allreduce(two_norm_residual,global_two_norm_residual,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)
    global_two_norm_residual = sqrt(global_two_norm_residual)
#endif

    converged_abs_update = PETSC_TRUE
    converged_scaled_residual = PETSC_TRUE
    do istate = 1, 3
      do idof = 1, option%nflowdof
        if (global_inf_norm_update(idof,istate) > &
          inf_norm_update_tol(idof,istate)) then
          converged_abs_update = PETSC_FALSE
        endif
        if (general_tough2_conv_criteria) then
          if (global_inf_norm_scaled_residual(idof,istate) > &
            general_tough2_itol_scaled_res_e1(idof,istate)) then
            converged_scaled_residual = PETSC_FALSE
          endif
        endif
      enddo  
    enddo  
    converged_rel_update = maxval(global_inf_norm_rel_update) < &
                                  general_itol_rel_update
    if (.not.general_tough2_conv_criteria) then
      converged_scaled_residual = maxval(global_inf_norm_scaled_residual) < &
                                  general_itol_scaled_res
    endif
#if 0
    do idof = 1, option%nflowdof
      if (global_inf_norm(idof) > option%flow%post_convergence_tol) then
        converged_rel_update = PETSC_FALSE
#ifdef DEBUG_GENERAL_INFO
        select case(istate_max(idof))
          case(1)
            state_char = 'L'
          case(2)
            state_char = 'G'
          case(3)
            state_char = '2P'
        end select
        write(*,'(''-+ '',a3,i2,''('',i5,''):'',es12.4, &
                 &'' dX_X/dX/X:'',3es12.4, &
                 &'' R_A/R/A:'',3es12.4)') state_char,idof, &
           icell_max(idof),global_inf_norm(idof), &
           dX_X1_max(idof), dX_max(idof),  X1_max(idof), &
           R_A_max(idof), R_max(idof), A_max(idof)
#endif
      endif
    enddo
#endif
#ifdef DEBUG_GENERAL_INFO
    write(*,'(4x,''-+  dpl:'',es12.4,''  dxa:'',es12.4,''  dt:'',es12.4)') &
      (max(global_inf_norm_update(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+  dpg:'',es12.4,''  dpa:'',es12.4,''  dt:'',es12.4)') &
      (max(global_inf_norm_update(idof,2),0.d0),idof=1,3)
    if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
      write(*,'(4x,''-+  dpg:'',es12.4,''  dsg:'',es12.4,''  dt:'',es12.4)') &
        (max(global_inf_norm_update(idof,3),0.d0),idof=1,3)
    else
      write(*,'(4x,''-+  dpg:'',es12.4,''  dsg:'',es12.4,'' dpa:'',es12.4)') &
        (max(global_inf_norm_update(idof,3),0.d0),idof=1,3)
    endif
    write(*,'(4x,''-+ rupl:'',es12.4,'' ruxa:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+ rupg:'',es12.4,'' rupa:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,2),0.d0),idof=1,3)
    write(*,'(4x,''-+ rupg:'',es12.4,'' rusg:'',es12.4,'' rut:'',es12.4)') &
      (max(global_inf_norm_rel_update(idof,3),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,1),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,2),0.d0),idof=1,3)
    write(*,'(4x,''-+  srl:'',es12.4,''  srg:'',es12.4,'' sre:'',es12.4)') &
      (max(global_inf_norm_scaled_residual(idof,3),0.d0),idof=1,3)
    write(*,'(4x,''-+ ru1 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(1), istate_max_rel_update(1), &
      X0_p((icell_max_rel_update(1)-1)*3+1), &
      -1.d0*dX_p((icell_max_rel_update(1)-1)*3+1), &
      r_p((icell_max_rel_update(1)-1)*3+1)
    write(*,'(4x,''-+ ru2 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(2), istate_max_rel_update(2), &
      X0_p((icell_max_rel_update(2)-1)*3+2), &
      -1.d0*dX_p((icell_max_rel_update(2)-1)*3+2), &
      r_p((icell_max_rel_update(2)-1)*3+2)
    write(*,'(4x,''-+ ru3 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_rel_update(3), istate_max_rel_update(3), &
      X0_p((icell_max_rel_update(3)-1)*3+3), &
      -1.d0*dX_p((icell_max_rel_update(3)-1)*3+3), &
      r_p((icell_max_rel_update(3)-1)*3+3)
    write(*,'(4x,''-+ sr1 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(1), istate_max_scaled_residual(1), &
      X0_p((icell_max_scaled_residual(1)-1)*3+1), &
      -1.d0*dX_p((icell_max_scaled_residual(1)-1)*3+1), &
      r_p((icell_max_scaled_residual(1)-1)*3+1)
    write(*,'(4x,''-+ sr2 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(2), istate_max_scaled_residual(2), &
      X0_p((icell_max_scaled_residual(2)-1)*3+2), &
      -1.d0*dX_p((icell_max_scaled_residual(2)-1)*3+2), &
      r_p((icell_max_scaled_residual(2)-1)*3+2)
    write(*,'(4x,''-+ sr3 icell:'',i7,''  st:'',i3,''  X:'',es11.3, &
              &''  dX:'',es11.3,''  R:'',es11.3)') &
      icell_max_scaled_residual(3), istate_max_scaled_residual(3), &
      X0_p((icell_max_scaled_residual(3)-1)*3+3), &
      -1.d0*dX_p((icell_max_scaled_residual(3)-1)*3+3), &
      r_p((icell_max_scaled_residual(3)-1)*3+3)
#endif
    option%converged = PETSC_FALSE
    if (converged_abs_update .or. converged_rel_update .or. &
        converged_scaled_residual) then
      option%converged = PETSC_TRUE
    endif
    call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
                               
    if (this%print_ekg) then
      call VecNorm(field%flow_r,NORM_2,two_norm,ierr);CHKERRQ(ierr)
      if (OptionPrintToFile(option)) then
  100 format("GENERAL NEWTON_ITERATION ",100es11.3)
        write(IUNIT_EKG,100) &
          global_inf_norm_update(:,:), &
          global_inf_norm_rel_update(:,:), &
          global_inf_norm_scaled_residual(:,:), &
          global_inf_norm_residual(:,:), &
          two_norm
      endif    
    endif                               
  endif                               

end subroutine PMGeneralCheckUpdatePost

! ************************************************************************** !

subroutine PMGeneralTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralTimeCut

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call GeneralTimeCut(this%realization)

end subroutine PMGeneralTimeCut

! ************************************************************************** !

subroutine PMGeneralUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralUpdateSolution, &
                             GeneralMapBCAuxVarsToGlobal

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call GeneralUpdateSolution(this%realization)
  call GeneralMapBCAuxVarsToGlobal(this%realization)

end subroutine PMGeneralUpdateSolution     

! ************************************************************************** !

subroutine PMGeneralUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use General_module, only : GeneralUpdateAuxVars

  implicit none
  
  class(pm_general_type) :: this

  call GeneralUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMGeneralUpdateAuxVars   

! ************************************************************************** !

subroutine PMGeneralMaxChange(this)
  ! 
  ! Not needed given GeneralMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use General_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_general_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  
  PetscErrorCode :: ierr
  
  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  do i = 1, 6
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4, &
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif

  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  ! ignore air pressure as it jumps during phase change
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_xmol_change = max_change_global(4)
  this%max_temperature_change = max_change_global(5)
  this%max_saturation_change = max_change_global(6)
  
end subroutine PMGeneralMaxChange

! ************************************************************************** !

subroutine PMGeneralComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralComputeMassBalance

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call GeneralComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMGeneralComputeMassBalance

! ************************************************************************** !

subroutine PMGeneralInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_general_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'general'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_RESIDUAL: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMGeneralInputRecord

! ************************************************************************** !

subroutine PMGeneralCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMGeneralCheckpointBinary

! ************************************************************************** !

subroutine PMGeneralRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  
end subroutine PMGeneralRestartBinary
! ************************************************************************** !

subroutine PMGeneralDestroy(this)
  ! 
  ! Destroys General process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralDestroy

  implicit none
  
  class(pm_general_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call GeneralDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMGeneralDestroy
  
end module PM_General_class
