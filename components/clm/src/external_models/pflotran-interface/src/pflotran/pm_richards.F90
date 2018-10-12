module PM_Richards_class

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

  type, public, extends(pm_subsurface_flow_type) :: pm_richards_type
  contains
    procedure, public :: Read => PMRichardsRead
    procedure, public :: InitializeTimestep => PMRichardsInitializeTimestep
    procedure, public :: Residual => PMRichardsResidual
    procedure, public :: Jacobian => PMRichardsJacobian
    procedure, public :: UpdateTimestep => PMRichardsUpdateTimestep
    procedure, public :: PreSolve => PMRichardsPreSolve
    procedure, public :: PostSolve => PMRichardsPostSolve
    procedure, public :: CheckUpdatePre => PMRichardsCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMRichardsCheckUpdatePost
    procedure, public :: TimeCut => PMRichardsTimeCut
    procedure, public :: UpdateSolution => PMRichardsUpdateSolution
    procedure, public :: UpdateAuxVars => PMRichardsUpdateAuxVars
    procedure, public :: MaxChange => PMRichardsMaxChange
    procedure, public :: ComputeMassBalance => PMRichardsComputeMassBalance
    procedure, public :: Destroy => PMRichardsDestroy
  end type pm_richards_type
  
  public :: PMRichardsCreate
  
contains

! ************************************************************************** !

function PMRichardsCreate()
  ! 
  ! Creates Richards process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_richards_type), pointer :: PMRichardsCreate

  class(pm_richards_type), pointer :: richards_pm
  
  allocate(richards_pm)
  call PMSubsurfaceFlowCreate(richards_pm)
  richards_pm%name = 'Richards Flow'

  PMRichardsCreate => richards_pm
  
end function PMRichardsCreate

! ************************************************************************** !

subroutine PMRichardsRead(this,input)
  ! 
  ! Reads input file parameters associated with the Richards process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use Richards_Aux_module
 
  implicit none
  
  class(pm_richards_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscReal :: Mannings_coeff
  PetscBool :: found

  option => this%option
  
  error_string = 'Richards Options'
  
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
      case('ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,richards_itol_scaled_res)
        call InputDefaultMsg(input,option,'itol_scaled_residual')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,richards_itol_rel_update)
        call InputDefaultMsg(input,option,'richards_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE
      case('INLINE_SURFACE_REGION')
        option%inline_surface_flow = PETSC_TRUE
        call InputReadWord(input,option,word,PETSC_FALSE)
        option%inline_surface_region_name = word
      case('INLINE_SURFACE_MANNINGS_COEFF')
        call InputReadDouble(input,option,Mannings_coeff)
        option%inline_surface_Mannings_coeff = Mannings_coeff
      case default
        call InputKeywordUnrecognized(word,error_string,option)
    end select
  enddo
  
end subroutine PMRichardsRead

! ************************************************************************** !

subroutine PMRichardsInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsInitializeTimestep
  
  implicit none
  
  class(pm_richards_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," RICHARDS FLOW ",63("="))')
  endif
  
  call RichardsInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)
  
end subroutine PMRichardsInitializeTimestep

! ************************************************************************** !

subroutine PMRichardsPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_richards_type) :: this

end subroutine PMRichardsPreSolve

! ************************************************************************** !

subroutine PMRichardsPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none
  
  class(pm_richards_type) :: this
  
end subroutine PMRichardsPostSolve

! ************************************************************************** !

subroutine PMRichardsUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
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
      ut = up
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
  
end subroutine PMRichardsUpdateTimestep

! ************************************************************************** !

subroutine PMRichardsResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsResidual

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call RichardsResidual(snes,xx,r,this%realization,ierr)

end subroutine PMRichardsResidual

! ************************************************************************** !

subroutine PMRichardsJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsJacobian

  implicit none
  
  class(pm_richards_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call RichardsJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMRichardsJacobian

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Characteristic_Curves_module
  use Patch_module
  use Richards_Aux_module
  use Global_Aux_module
  use Patch_module
  
  implicit none
  
  class(pm_richards_type) :: this
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
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P_R, P0, P1, delP
  PetscReal :: scale, sat, sat_pert, pert, pc_pert, press_pert, delP_pert
  
  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars

  if (Initialized(this%saturation_change_limit)) then

    changed = PETSC_TRUE

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    pert = dabs(this%saturation_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      sat = global_auxvars(ghosted_id)%sat(1)
      sat_pert = sat - sign(1.d0,sat-0.5d0)*pert
      call patch%characteristic_curves_array( &
             patch%sat_func_id(ghosted_id))%ptr% &
             saturation_function%CapillaryPressure(sat_pert,pc_pert,option)
      press_pert = option%reference_pressure - pc_pert
      P0 = X_p(local_id)
      delP = dX_p(local_id)
      delP_pert = dabs(P0 - press_pert)
      if (delP_pert < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),delP_pert,dabs(delP)
        call printMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),delP_pert),delP)
      dX_p(local_id) = delP
    enddo
    
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)

  endif

  if (Initialized(this%pressure_dampening_factor)) then
    changed = PETSC_TRUE
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = this%pressure_dampening_factor

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      delP = dX_p(local_id)
      P0 = X_p(local_id)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1 
        call printMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
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
        write(option%io_buffer,'("Residual:",es15.7)') r_p(local_id)
        call printMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dX_p(local_id) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine PMRichardsCheckUpdatePre

! ************************************************************************** !

subroutine PMRichardsCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                     X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Richards_Aux_module
  use Global_Aux_module
  use Material_Aux_class
  use Patch_module
  use Richards_Common_module

  implicit none
  
  class(pm_richards_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(richards_auxvar_type), pointer :: rich_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart
  PetscReal :: Res(1)
  PetscReal :: inf_norm, global_inf_norm
  
  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  rich_auxvars => patch%aux%Richards%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  option%converged = PETSC_FALSE
  if (this%check_post_convergence) then
    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    
    inf_norm = 0.d0
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      istart = (local_id-1)*option%nflowdof + 1

      if (patch%imat(ghosted_id) <= 0) cycle
    
      call RichardsAccumulation(rich_auxvars(ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                option,Res)
      inf_norm = max(inf_norm,min(dabs(dX_p(local_id)/X0_p(local_id)), &
                                  dabs(r_p(istart)/Res(1))))
    enddo
    call MPI_Allreduce(inf_norm,global_inf_norm,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MAX,option%mycomm,ierr)
    option%converged = PETSC_TRUE
    if (global_inf_norm > richards_itol_scaled_res) &
      option%converged = PETSC_FALSE
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X0,X0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine PMRichardsCheckUpdatePost

! ************************************************************************** !

subroutine PMRichardsTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsTimeCut

  implicit none
  
  class(pm_richards_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call RichardsTimeCut(this%realization)

end subroutine PMRichardsTimeCut

! ************************************************************************** !

subroutine PMRichardsUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsUpdateSolution, &
                              RichardsUpdateSurfacePress

  implicit none
  
  class(pm_richards_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call RichardsUpdateSolution(this%realization)
  if (this%option%surf_flow_on) &
    call RichardsUpdateSurfacePress(this%realization)

end subroutine PMRichardsUpdateSolution     

! ************************************************************************** !

subroutine PMRichardsUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Richards_module, only : RichardsUpdateAuxVars
  
  implicit none
  
  class(pm_richards_type) :: this

  call RichardsUpdateAuxVars(this%realization)

end subroutine PMRichardsUpdateAuxVars   

! ************************************************************************** !

subroutine PMRichardsMaxChange(this)
  ! 
  ! Not needed given RichardsMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsMaxChange

  implicit none
  
  class(pm_richards_type) :: this
  
  call RichardsMaxChange(this%realization,this%max_pressure_change)
  if (this%option%print_screen_flag) then
    write(*,'("  --> max chng: dpmx= ",1pe12.4)') this%max_pressure_change
  endif
#ifndef CLM_PFLOTRAN
  if (this%option%print_file_flag) then
    write(this%option%fid_out,'("  --> max chng: dpmx= ",1pe12.4)') &
      this%max_pressure_change
  endif    
#endif
end subroutine PMRichardsMaxChange

! ************************************************************************** !

subroutine PMRichardsComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsComputeMassBalance

  implicit none
  
  class(pm_richards_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call RichardsComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMRichardsComputeMassBalance

! ************************************************************************** !

subroutine PMRichardsDestroy(this)
  ! 
  ! Destroys Richards process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Richards_module, only : RichardsDestroy

  implicit none
  
  class(pm_richards_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  ! preserve this ordering
  call RichardsDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMRichardsDestroy
  
end module PM_Richards_class
