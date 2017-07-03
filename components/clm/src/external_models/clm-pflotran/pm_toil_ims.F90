module PM_TOilIms_class

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

  type, public, extends(pm_subsurface_flow_type) :: pm_toil_ims_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
  contains
    ! all the routines below needs to be replaced, uncomment as I develop them
    procedure, public :: Read => PMTOilImsRead
    procedure, public :: InitializeRun => PMTOilImsInitializeRun
    procedure, public :: InitializeTimestep => PMTOilImsInitializeTimestep
    procedure, public :: Residual => PMTOilImsResidual
    procedure, public :: Jacobian => PMTOilImsJacobian
    procedure, public :: UpdateTimestep => PMTOilImsUpdateTimestep
    procedure, public :: PreSolve => PMTOilImsPreSolve
    !procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMTOilImsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTOilImsCheckUpdatePost
    procedure, public :: TimeCut => PMTOilImsTimeCut
    procedure, public :: UpdateSolution => PMTOilImsUpdateSolution
    procedure, public :: UpdateAuxVars => PMTOilImsUpdateAuxVars
    procedure, public :: MaxChange => PMTOilImsMaxChange
    !procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: CheckpointBinary => PMTOilImsCheckpointBinary
    procedure, public :: RestartBinary => PMTOilImsRestartBinary
    procedure, public :: InputRecord => PMTOilImsInputRecord
    procedure, public :: Destroy => PMTOilImsDestroy
  end type pm_toil_ims_type
  
  public :: PMToilImsCreate
  
contains

! ************************************************************************** !

function PMTOilImsCreate()
  ! 
  ! Creates TOilIms process models shell
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 9/8/15
  ! 
  use Variables_module, only : LIQUID_PRESSURE, OIL_PRESSURE, OIL_SATURATION, &
                               TEMPERATURE
  use TOilIms_module, only : TOilImsDefaultSetup
  implicit none
  
  class(pm_toil_ims_type), pointer :: PMToilImsCreate

  class(pm_toil_ims_type), pointer :: toil_ims_pm
  
!#ifdef PM_TOIL_IMS_DEBUG  
  print *, 'PMTOilImsCreate()'
!#endif  

  allocate(toil_ims_pm)

  allocate(toil_ims_pm%max_change_ivar(4))
  toil_ims_pm%max_change_ivar = [LIQUID_PRESSURE, OIL_PRESSURE, &
                                OIL_SATURATION, TEMPERATURE]
  allocate(toil_ims_pm%max_change_isubvar(4))
  toil_ims_pm%max_change_isubvar = [0,0,0,0]
  
  call PMSubsurfaceFlowCreate(toil_ims_pm)
  toil_ims_pm%name = 'TOilIms Flow'

  call TOilImsDefaultSetup()

  PMTOilImsCreate => toil_ims_pm
  
end function PMTOilImsCreate

! ************************************************************************** !

subroutine PMTOilImsRead(this,input)
  ! 
  ! Reads input specific to pm_toil_Ims.
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: Date: 9/9/15
  !
  use PM_TOilIms_Aux_module 
  use Input_Aux_module
  use String_module
  use Option_module
  use TOilIms_module, only : TOilImsFluxDipcSetup

  implicit none

  class(pm_toil_ims_type) :: this  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'TOilIms Options'  

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
        call InputReadDouble(input,option,toil_ims_itol_scaled_res)
        call InputDefaultMsg(input,option,'toil_ims_itol_scaled_res')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,toil_ims_itol_rel_update)
        call InputDefaultMsg(input,option,'toil_ims_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE        
      case('TOUGH2_ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,tempreal)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e1')
        toil_ims_tgh2_itol_scld_res_e1 = tempreal
        call InputReadDouble(input,option,toil_ims_tgh2_itol_scld_res_e2)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e2')
        toil_ims_tough2_conv_criteria = PETSC_TRUE
        this%check_post_convergence = PETSC_TRUE
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,toil_ims_window_epsilon)
        call InputErrorMsg(input,option,'window epsilon',error_string)
      case('ISOTHERMAL')
        toil_ims_isothermal = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,toil_ims_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,toil_ims_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,toil_ims_damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
#if 0
      case('GOVERN_MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%pressure_change_governor)
        call InputErrorMsg(input,option,'maximum allowable pressure change', &
                           error_string)
      case('GOVERN_MAXIMUM_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,this%temperature_change_governor)
        call InputErrorMsg(input,option, &
                           'maximum allowable temperature change', &
                           error_string)
      case('GOVERN_MAXIMUM_SATURATION_CHANGE')
        call InputReadDouble(input,option,this%saturation_change_governor)
        call InputErrorMsg(input,option,'maximum allowable saturation change', &
                           error_string)
#endif
      case('DEBUG_CELL')
        call InputReadInt(input,option,toil_ims_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id',error_string)
      ! might need some input here for the thermal diffusion model
      !case('NO_TEMP_DEPENDENT_DIFFUSION')
      !  general_temp_dep_gas_air_diff = PETSC_FALSE
      !case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_TRUE
      !case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
      !  general_harmonic_diff_density = PETSC_FALSE
      case('FLUX_DIPC')
        call TOilImsFluxDipcSetup()
      case default
        call InputKeywordUnrecognized(keyword,'TOIL_IMS Mode',option)
    end select
    
  enddo  
  
end subroutine PMTOilImsRead

! ************************************************************************** !

recursive subroutine PMTOilImsInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  use Realization_Base_class
  
  implicit none
  
  class(pm_toil_ims_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,FOUR_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 4
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo


  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMTOilImsInitializeRun

! ************************************************************************** !

subroutine PMTOilImsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use TOilIms_module, only : TOilImsInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_toil_ims_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)
                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," TOIL_IMS FLOW ",64("="))')
  endif
  
  call TOilImsInitializeTimestep(this%realization)

  call PMSubsurfaceFlowInitializeTimestepB(this)                                 
  
end subroutine PMTOilImsInitializeTimestep

! ************************************************************************** !

subroutine PMTOilImsPreSolve(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  implicit none

  class(pm_toil_ims_type) :: this

  ! currently does nothing - could add here explicit iitialization
  ! for highly het. problems

end subroutine PMTOilImsPreSolve

! ************************************************************************** !

subroutine PMTOilImsUpdateAuxVars(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15

  use TOilIms_module, only : TOilImsUpdateAuxVars

  implicit none
  
  class(pm_toil_ims_type) :: this

  call TOilImsUpdateAuxVars(this%realization)

end subroutine PMTOilImsUpdateAuxVars   

! ************************************************************************** !

subroutine PMTOilImsUpdateSolution(this)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/23/15
  ! 

  use TOilIms_module, only : TOilImsUpdateSolution, &
                             TOilImsMapBCAuxVarsToGlobal 

  implicit none
  
  class(pm_toil_ims_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call TOilImsUpdateSolution(this%realization)
  call TOilImsMapBCAuxVarsToGlobal(this%realization)

end subroutine PMTOilImsUpdateSolution     


! ************************************************************************** !

subroutine PMTOilImsUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! Date modified : 08/08/16 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, OIL_SATURATION

  implicit none
  
  class(pm_toil_ims_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, us, umin
  PetscReal :: dtt
  type(field_type), pointer :: field
    
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    !up = this%dPmax_allowable/(this%dPmax+0.1)
    !ut = this%dTmax_allowable/(this%dTmax+1.d-5)
    !us = this%dSmax_allowable/(this%dSmax+1.d-5)
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
    us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
    umin = min(up,ut,us)
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
                                OIL_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               OIL_SATURATION,TIME_NULL)
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  endif

end subroutine PMTOilImsUpdateTimestep

! ************************************************************************** !

subroutine PMTOilImsResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/15
  ! 

  use TOilIms_module, only : TOilImsResidual

  implicit none
  
  class(pm_toil_ims_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  ! in theroy call for material properties update - currently does nothing 
  call PMSubsurfaceFlowUpdatePropertiesNI(this) 
  call TOilImsResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTOilImsResidual

! ************************************************************************** !

subroutine PMTOilImsJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/07/15
  ! 

  use TOilIms_module, only : ToilImsJacobian

  implicit none
  
  class(pm_toil_ims_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call TOilImsJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTOilImsJacobian


! ************************************************************************** !

!subroutine PMTOilImsCheckUpdatePre(this,line_search,P,dP,changed,ierr)
!  ! 
!  ! Author: Paolo Orsini (OGS)
!  ! Date: 10/22/15
!  ! 
!
!  use TOilIms_module, only : TOilImsCheckUpdatePre
!
!  implicit none
!  
!  class(pm_toil_ims_type) :: this
!  SNESLineSearch :: line_search
!  Vec :: P
!  Vec :: dP
!  PetscBool :: changed
!  PetscErrorCode :: ierr
!  
!  call TOilImsCheckUpdatePre(line_search,P,dP,changed,this%realization,ierr)
!
!end subroutine PMTOilImsCheckUpdatePre

! ************************************************************************** !
! use function below when switching to latest code
subroutine PMTOilImsCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Paolo Orsini (OGS)
  ! Date: 11/09/15
  ! 
  !use Realization_Subsurface_class
  use Grid_module
  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module
  !use Global_Aux_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none
  
  class(pm_toil_ims_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X_p(:), dX_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  !type(toil_ims_auxvar_type), pointer :: toil_auxvars(:,:)
  !type(global_auxvar_type), pointer :: global_auxvars(:)  

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset

  PetscInt :: pressure_index, saturation_index, temperature_index

  PetscReal :: pressure0, pressure1, del_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation

  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration

  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  !toil_auxvars => this%realization%patch%aux%TOil_ims%auxvars
  !global_auxvars => this%realization%patch%aux%Global%auxvars

  patch => this%realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  ! truncation
  ! Oil Saturation must be truncated.  we do not use scaling
  ! here because of the very small values.  just truncation.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    saturation_index = offset + TOIL_IMS_SATURATION_DOF
    if ( (X_p(saturation_index) - dX_p(saturation_index)) < 0.d0 ) then
      ! we use 1.d-6 since cancelation can occur with smaller values
      ! this threshold is imposed in the initial condition
      dX_p(saturation_index) = X_p(saturation_index)
    end if
  enddo

  scale = initial_scale
  if (toil_ims_max_it_before_damping > 0 .and. &
      newton_iteration > toil_ims_max_it_before_damping) then
    scale = toil_ims_damping_factor
  endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
!!#define LIMIT_MAX_TEMPERATURE_CHANGE
!! TRUNCATE_PRESSURE is needed for times when the solve wants
!! to pull them negative.
!#define TRUNCATE_PRESSURE

  ! scaling
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    offset = (local_id-1)*option%nflowdof
    temp_scale = 1.d0
    pressure_index = offset + TOIL_IMS_PRESSURE_DOF
    saturation_index = offset + TOIL_IMS_SATURATION_DOF
    temperature_index  = offset + TOIL_IMS_ENERGY_DOF
    dX_p(pressure_index) = dX_p(pressure_index) * toil_ims_pressure_scale
    temp_scale = 1.d0
    del_pressure = dX_p(pressure_index)
    pressure0 = X_p(pressure_index)
    pressure1 = pressure0 - del_pressure
    del_saturation = dX_p(saturation_index)
    saturation0 = X_p(saturation_index)
    saturation1 = saturation0 - del_saturation
#ifdef LIMIT_MAX_PRESSURE_CHANGE
    if (dabs(del_pressure) > toil_ims_max_pressure_change) then
      temp_real = dabs(toil_ims_max_pressure_change/del_pressure)
      temp_scale = min(temp_scale,temp_real)
     endif
#endif
#ifdef TRUNCATE_PRESSURE
    if (pressure1 <= 0.d0) then
      if (dabs(del_pressure) > 1.d-40) then
        temp_real = tolerance * dabs(pressure0 / del_pressure)
        temp_scale = min(temp_scale,temp_real)
      endif
    endif
#endif 
!TRUNCATE_PRESSURE

#ifdef LIMIT_MAX_SATURATION_CHANGE
    if (dabs(del_saturation) > max_saturation_change) then
       temp_real = dabs(max_saturation_change/del_saturation)
       temp_scale = min(temp_scale,temp_real)
    endif
#endif 
!LIMIT_MAX_SATURATION_CHANGE        
    scale = min(scale,temp_scale) 
  enddo

  temp_scale = scale
  call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION, &
                     MPI_MIN,option%mycomm,ierr)

  ! it performs an homogenous scaling using the smallest scaling factor
  ! over all subdomains domains
  if (scale < 0.9999d0) then
    dX_p = scale*dX_p
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMTOilImsCheckUpdatePre

! ************************************************************************** !

! ************************************************************************** !

!subroutine PMTOilImsCheckUpdatePost(this,line_search,P0,dP,P1,dP_changed, &
!                                    P1_changed,ierr)
!  ! 
!  ! Author: Paolo Orsini
!  ! Date: 11/09/15
!  ! 
!
!  use TOilIms_module, only : TOilImsCheckUpdatePost
!
!  implicit none
!  
!  class(pm_toil_ims_type) :: this
!  SNESLineSearch :: line_search
!  Vec :: P0
!  Vec :: dP
!  Vec :: P1
!  PetscBool :: dP_changed
!  PetscBool :: P1_changed
!  PetscErrorCode :: ierr
!  
!  call TOilImsCheckUpdatePost(line_search,P0,dP,P1,dP_changed, &
!                                   P1_changed,this%realization,ierr)
!
!end subroutine PMTOilImsCheckUpdatePost

! ************************************************************************** !
! use function below when switching to the latest code
subroutine PMTOilImsCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! 
  !use Global_Aux_module
  !use TOilIms_Aux_module
  use PM_TOilIms_Aux_module
  use Grid_module
  use Option_module
  !use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  !use Output_EKG_module
  
  implicit none
  
  class(pm_toil_ims_type) :: this
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
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset , ival, idof
  PetscReal :: dX_X0, R_A, R

  PetscReal :: inf_norm_rel_update(3), global_inf_norm_rel_update(3)
  PetscReal :: inf_norm_scaled_residual(3), global_inf_norm_scaled_residual(3)
  PetscReal :: inf_norm_update(3), global_inf_norm_update(3)
  PetscReal :: inf_norm_residual(3), global_inf_norm_residual(3)
  PetscReal :: two_norm_residual(3), global_two_norm_residual(3)
  PetscReal, parameter :: inf_pres_tol = 1.d-1
  PetscReal, parameter :: inf_temp_tol = 1.d-5
  PetscReal, parameter :: inf_sat_tol = 1.d-6
  !geh: note the scaling by 0.d0 several lines down which prevent false 
  !     convergence 
  ! PO scaling by 0 kill the inf_norm_update convergence criteria
  PetscReal, parameter :: inf_norm_update_tol(3) = &
    reshape([inf_pres_tol,inf_sat_tol,inf_temp_tol], &
            shape(inf_norm_update_tol)) * &
            0.d0
  PetscReal :: temp(12), global_temp(12)
  PetscMPIInt :: mpi_int
  PetscBool :: converged_abs_update
  PetscBool :: converged_rel_update
  PetscBool :: converged_scaled_residual
  PetscReal :: t_over_v
 
  grid => this%realization%patch%grid 
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch ! in patch imat for active/inactive cells
  material_auxvars => patch%aux%Material%auxvars 
 
  ! it indicates that neither dX of the updated solution are modified 
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  option%converged = PETSC_FALSE
  if (this%check_post_convergence) then
    call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)

    inf_norm_update(:) = -1.d20
    inf_norm_rel_update(:) = -1.d20
    inf_norm_scaled_residual(:) = -1.d20
    inf_norm_residual(:) = -1.d20
    two_norm_residual(:) = 0.d0
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      do idof = 1, option%nflowdof
        ival = offset+idof
        R = r_p(ival)
        inf_norm_residual(idof) = max(inf_norm_residual(idof),dabs(R))
        if (toil_ims_tough2_conv_criteria) then
          !geh: scale by t_over_v to match TOUGH2 residual units. see equation
          !     B.5 of TOUGH2 user manual (LBNL-43134)
          t_over_v = option%flow_dt/material_auxvars(ghosted_id)%volume
          if (accum_p2(ival)*t_over_v < toil_ims_tgh2_itol_scld_res_e2) then
            R_A = dabs(R*t_over_v)
          else
            R_A = dabs(R/accum_p2(ival))
          endif
        else
          R_A = dabs(R/accum_p(ival))
        endif
        dX_X0 = dabs(dX_p(ival)/X0_p(ival))
        inf_norm_update(idof) = max(inf_norm_update(idof),dabs(dX_p(ival)))
        if (inf_norm_rel_update(idof) < dX_X0) then
          inf_norm_rel_update(idof) = dX_X0
        endif
        if (inf_norm_scaled_residual(idof) < R_A) then
          inf_norm_scaled_residual(idof) = R_A
        endif
      enddo
    enddo
    temp(1:3) = inf_norm_update(:)
    temp(4:6) = inf_norm_rel_update(:)
    temp(7:9) = inf_norm_scaled_residual(:)
    temp(10:12) = inf_norm_residual(:)
    mpi_int = 12
    call MPI_Allreduce(temp,global_temp,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
    global_inf_norm_update(:) = global_temp(1:3)
    global_inf_norm_rel_update(:) = global_temp(4:6)
    global_inf_norm_scaled_residual(:) = global_temp(7:9)
    global_inf_norm_residual(:) = global_temp(10:12)

    converged_abs_update = PETSC_TRUE
    converged_scaled_residual = PETSC_TRUE
    do idof = 1, option%nflowdof
      ! imposing inf_norm_update <= inf_norm_update_tol for convergence
      if (global_inf_norm_update(idof) > inf_norm_update_tol(idof)) then
        converged_abs_update = PETSC_FALSE
      endif
      if (toil_ims_tough2_conv_criteria) then
        if (global_inf_norm_scaled_residual(idof) > &
            toil_ims_tgh2_itol_scld_res_e1(idof)) then
          converged_scaled_residual = PETSC_FALSE
        endif
      endif
    enddo  

    if (.not.toil_ims_tough2_conv_criteria) then
      converged_scaled_residual = maxval(global_inf_norm_scaled_residual) < &
                                  toil_ims_itol_scaled_res
    endif

    ! global_inf_norm_rel_update alway >0 because dabs values 
    ! when not inu, toil_ims_itol_rel_update < 0  because assigned uninitialized value (-999)
    converged_rel_update = maxval(global_inf_norm_rel_update) < &
                                  toil_ims_itol_rel_update  

   ! converged_rel_update = maxval(global_inf_norm_rel_update) < &
   !                        option%flow%inf_rel_update_tol
   ! if (toil_ims_tough2_conv_criteria) then
   !   converged_scaled_residual = maxval(global_inf_norm_scaled_residual) < &
   !                               toil_ims_tgh2_itol_scld_res_e1
   ! else
   !   converged_scaled_residual = maxval(global_inf_norm_scaled_residual) < &
   !                               option%flow%inf_scaled_res_tol
   ! endif
#if 0
    do idof = 1, option%nflowdof
      if (global_inf_norm(idof) > option%flow%post_convergence_tol) then
        converged_rel_update = PETSC_FALSE
      endif
    enddo
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

  endif

end subroutine PMTOilImsCheckUpdatePost

! ************************************************************************** !

subroutine PMTOilImsTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use TOilIms_module, only : TOilImsTimeCut

  implicit none
  
  class(pm_toil_ims_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call TOilImsTimeCut(this%realization)

end subroutine PMTOilImsTimeCut

! ************************************************************************** !

! ************************************************************************** !

subroutine PMTOilImsMaxChange(this)
  ! 
  ! Not needed given ToilImsMaxChange is called in PostSolve
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  !use General_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, OIL_PRESSURE, OIL_SATURATION, &
                               TEMPERATURE
  implicit none
  
  class(pm_toil_ims_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(4)
  PetscReal :: max_change_global(4)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  PetscErrorCode :: ierr
  
  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  ! max changes loaded in this%max_change_ivar(i) with the following order:
  ! 1. LIQUID_PRESSURE, 2. OIL_PRESSURE, 3.OIL_SATURATION, 4.TEMPERATURE

  max_change_global = 0.d0
  max_change_local = 0.d0
  do i = 1,4
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
  call MPI_Allreduce(max_change_local,max_change_global,FOUR_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpo= ",1pe12.4,&
      & "  dso= ",1pe12.4,&
      & " dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpo= ",1pe12.4,&
      & "  dso= ",1pe12.4, &
      & " dt= ",1pe12.4)') &
      max_change_global(1:4)
  endif
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_saturation_change = max_change_global(3)
  this%max_temperature_change = max_change_global(4)
  
end subroutine PMTOilImsMaxChange

! ************************************************************************** !

subroutine PMTOilImsCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with General PM
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_toil_ims_type) :: this
  PetscViewer :: viewer
 
  ! currently doing this but it is not needed 
  call GlobalGetAuxVarVecLoc(this%realization, &
                             this%realization%field%iphas_loc, &
                             STATE,ZERO_INTEGER)
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMTOilImsCheckpointBinary

! ************************************************************************** !

subroutine PMTOilImsRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with General PM
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_toil_ims_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  ! currently doing this but it is not needed for TOIL_IMS
  call GlobalSetAuxVarVecLoc(this%realization, &
                             this%realization%field%iphas_loc, &
                             STATE,ZERO_INTEGER)
  
end subroutine PMTOilImsRestartBinary

! ************************************************************************** !

subroutine PMTOilImsInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_toil_ims_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'thermal oil immiscible'

end subroutine PMTOilImsInputRecord

! ************************************************************************** !

subroutine PMTOilImsDestroy(this)
  ! 
  ! Destroys General process model
  ! 
  ! Author: Paolo Orsini
  ! Date: 11/09/15
  ! 

  use TOilIms_module, only : TOilImsDestroy

  implicit none
  
  class(pm_toil_ims_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call TOilImsDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMTOilImsDestroy

! ************************************************************************** !

end module PM_TOilIms_class
