module PMC_Geomechanics_class

  use PMC_Base_class
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  private

  type, public, extends(pmc_base_type) :: pmc_geomechanics_type
    class(realization_subsurface_type), pointer :: subsurf_realization
    class(realization_geomech_type), pointer :: geomech_realization
  contains
    procedure, public :: Init => PMCGeomechanicsInit
    procedure, public :: SetupSolvers => PMCGeomechanicsSetupSolvers
    procedure, public :: RunToTime => PMCGeomechanicsRunToTime
    procedure, public :: GetAuxData => PMCGeomechanicsGetAuxData
    procedure, public :: SetAuxData => PMCGeomechanicsSetAuxData
    procedure, public :: Destroy => PMCGeomechanicsDestroy
  end type pmc_geomechanics_type

  public :: PMCGeomechanicsCreate

contains

! ************************************************************************** !

function PMCGeomechanicsCreate()
  ! 
  ! This routine allocates and initializes a new object.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(pmc_geomechanics_type), pointer :: PMCGeomechanicsCreate
  
  class(pmc_geomechanics_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCGeomechanicsCreate%Create()'
#endif

  allocate(pmc)
  call pmc%Init()
  
  PMCGeomechanicsCreate => pmc  
  
end function PMCGeomechanicsCreate

! ************************************************************************** !

subroutine PMCGeomechanicsInit(this)
  ! 
  ! This routine initializes a new process model coupler object.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  implicit none
  
  class(pmc_geomechanics_type) :: this
  
#ifdef DEBUG
  print *, 'PMCGeomechanics%Init()'
#endif

  call PMCBaseInit(this)
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsInit

! ************************************************************************** !

subroutine PMCGeomechanicsSetupSolvers(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 
  use Convergence_module
  use Geomechanics_Discretization_module
  use Timestepper_Base_class
  use Timestepper_Steady_class
  use PM_Base_class
  use PM_Base_Pointer_module
  use Option_module
  use Solver_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscpc.h"

  class(pmc_geomechanics_type) :: this

  class(realization_geomech_type), pointer :: geomech_realization
  class(geomech_discretization_type), pointer :: geomech_discretization
  class(timestepper_steady_type), pointer :: ts_steady
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

#ifdef DEBUG
  call printMsg(this%option,'PMCGeomechanicsSetupSolvers')
#endif

  option => this%option
  geomech_realization => this%geomech_realization
  geomech_discretization => geomech_realization%geomech_discretization
  select type(ts => this%timestepper)
    class is (timestepper_steady_type)
      ts_steady => ts
      solver => ts%solver
  end select 

  call printMsg(option,"  Beginning setup of GEOMECH SNES ")

  if (solver%J_mat_type == MATAIJ) then
    option%io_buffer = 'AIJ matrix not supported for geomechanics.'
    call printErrMsg(option)
  endif

  call SolverCreateSNES(solver,option%mycomm)
  call SNESSetOptionsPrefix(solver%snes, "geomech_", &
                            ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  if (solver%Jpre_mat_type == '') then
    solver%Jpre_mat_type = solver%J_mat_type
  endif
  call GeomechDiscretizationCreateJacobian(geomech_realization% &
                                           geomech_discretization,NGEODOF, &
                                           solver%Jpre_mat_type, &
                                           solver%Jpre,option)

  solver%J = solver%Jpre
  call MatSetOptionsPrefix(solver%Jpre,"geomech_", &
                            ierr);CHKERRQ(ierr)

  ! by default turn off line search
  call SNESGetLineSearch(solver%snes,linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC, &
                              ierr);CHKERRQ(ierr)


  ! Have PETSc do a SNES_View() at the end of each solve if verbosity > 0.
  if (option%verbosity >= 1) then
    string = '-geomech_snes_view'
    call PetscOptionsInsertString(PETSC_NULL_OBJECT, &
                                   string, ierr);CHKERRQ(ierr)
  endif

  option%io_buffer = 'Solver: ' // trim(solver%ksp_type)
  call printMsg(option)
  option%io_buffer = 'Preconditioner: ' // trim(solver%pc_type)
  call printMsg(option)

  ! shell for custom convergence test.  The default SNES convergence test
  ! is call within this function.
  ts_steady%convergence_context => &
             ConvergenceContextCreate(solver,option, &
                                      this%subsurf_realization%patch%grid)
  call SNESSetConvergenceTest(solver%snes,ConvergenceTest, &
                              ts_steady%convergence_context, &
                              PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)            
  call SNESSetFunction(solver%snes, &
                       this%pm_ptr%pm%residual_vec, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                       PMResidual, &
                       this%pm_ptr%pm, &
#else
                       PMResidualPtr, &
                       this%pm_ptr, &
#endif
                       ierr);CHKERRQ(ierr)
  call SNESSetJacobian(solver%snes, &
                       solver%J, &
                       solver%Jpre, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                       PMJacobian, &
                       this%pm_ptr%pm, &
#else
                       PMJacobianPtr, &
                       this%pm_ptr, &
#endif
                       ierr);CHKERRQ(ierr)

  call SolverSetSNESOptions(solver,option)

  call printMsg(option,"  Finished setting up GEOMECH SNES ")


end subroutine PMCGeomechanicsSetupSolvers

! ************************************************************************** !

recursive subroutine PMCGeomechanicsRunToTime(this,sync_time,stop_flag)
  ! 
  ! This routine runs the geomechanics simulation.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Timestepper_Base_class
  use Option_module
  use PM_Base_class
  use Output_Geomechanics_module

  implicit none

  class(pmc_geomechanics_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  PetscInt :: local_stop_flag
  PetscBool :: snapshot_plot_flag
  PetscBool :: observation_plot_flag
  PetscBool :: massbal_plot_flag
  PetscBool :: checkpoint_flag  
    
  class(pm_base_type), pointer :: cur_pm

  this%option%io_buffer = trim(this%name) // ':' // trim(this%pm_list%name)
  call printVerboseMsg(this%option)
  
  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = 0

  call SetOutputFlags(this)
  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE

  call this%timestepper%SetTargetTime(sync_time,this%option,local_stop_flag, &
                                      snapshot_plot_flag, &
                                      observation_plot_flag, &
                                      massbal_plot_flag,checkpoint_flag)
  call this%timestepper%StepDT(this%pm_list,local_stop_flag)
  
  ! Check if it is initial solve
  if (this%timestepper%steps == 1) then
    this%option%geomech_initial = PETSC_TRUE
  endif

  ! Have to loop over all process models coupled in this object and update
  ! the time step size.  Still need code to force all process models to
  ! use the same time step size if tightly or iteratively coupled.
  cur_pm => this%pm_list
  do
    if (.not.associated(cur_pm)) exit
    ! have to update option%time for conditions
    this%option%time = this%timestepper%target_time
    call cur_pm%UpdateSolution()
    ! Geomechanics PM does not have an associate time 
    !call this%timestepper%UpdateDT(cur_pm)
    cur_pm => cur_pm%next
  enddo

  ! Run underlying process model couplers
  if (associated(this%child)) then
    ! Set data needed by process-model
    call this%SetAuxData()
    call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
  endif
  
  if (this%timestepper%time_step_cut_flag) then
    snapshot_plot_flag = PETSC_FALSE
  endif
  ! however, if we are using the modulus of the output_option%imod, we may
  ! still print
  if (mod(this%timestepper%steps,this%pm_list% &
          output_option%periodic_snap_output_ts_imod) == 0) then
    snapshot_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_obs_output_ts_imod) == 0) then
    observation_plot_flag = PETSC_TRUE
  endif
  if (mod(this%timestepper%steps,this%pm_list%output_option% &
          periodic_msbl_output_ts_imod) == 0) then
    massbal_plot_flag = PETSC_TRUE
  endif
    
  call OutputGeomechanics(this%geomech_realization,snapshot_plot_flag, &
                          observation_plot_flag,massbal_plot_flag)
  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer)) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)

end subroutine PMCGeomechanicsRunToTime

! ************************************************************************** !

subroutine PMCGeomechanicsSetAuxData(this)
  ! 
  ! This routine updates data in simulation_aux that is required by other
  ! process models.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module
  use Grid_module

  implicit none

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"

  class(pmc_geomechanics_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(grid_type), pointer :: grid
  PetscInt :: local_id
  PetscScalar, pointer :: por0_p(:)
  PetscScalar, pointer :: por_p(:)
  PetscScalar, pointer :: strain_p(:)
  PetscErrorCode :: ierr
  PetscReal :: trace_epsilon
  PetscReal :: por_new

  ! If at initialization stage, do nothing
  if (this%timestepper%steps == 0) return

  select type(pmc => this)
    class is(pmc_geomechanics_type)
      if (this%option%geomech_subsurf_coupling == GEOMECH_TWO_WAY_COUPLED) then

        grid => pmc%subsurf_realization%patch%grid

        ! Save strain dataset in sim_aux%subsurf_strain
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%strain, &
                             pmc%sim_aux%subsurf_strain, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%strain, &
                           pmc%sim_aux%subsurf_strain, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
                             
        ! Save stress dataset in sim_aux%subsurf_stress
        call VecScatterBegin(pmc%sim_aux%geomechanics_to_subsurf, &
                             pmc%geomech_realization%geomech_field%stress, &
                             pmc%sim_aux%subsurf_stress, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(pmc%sim_aux%geomechanics_to_subsurf, &
                           pmc%geomech_realization%geomech_field%stress, &
                           pmc%sim_aux%subsurf_stress, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

        ! Update porosity dataset in sim_aux%subsurf_por
        call VecGetArrayF90(pmc%sim_aux%subsurf_por0, por0_p,  &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_por, por_p,  &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_strain, strain_p,  &
                            ierr);CHKERRQ(ierr)

        do local_id = 1, grid%nlmax
          trace_epsilon = strain_p((local_id - 1)*SIX_INTEGER + ONE_INTEGER) + &
                          strain_p((local_id - 1)*SIX_INTEGER + TWO_INTEGER) + &
                          strain_p((local_id - 1)*SIX_INTEGER + THREE_INTEGER)
          por_new = por0_p(local_id)/(1.d0 + (1.d0 - por0_p(local_id))*trace_epsilon)
          por_p(local_id) = por_new
        enddo

        call VecRestoreArrayF90(pmc%sim_aux%subsurf_por0, por0_p,  &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_strain, strain_p,  &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_por, por_p,  &
                                ierr);CHKERRQ(ierr)

      endif

  end select

end subroutine PMCGeomechanicsSetAuxData

! ************************************************************************** !

subroutine PMCGeomechanicsGetAuxData(this)
  ! 
  ! This routine updates data for geomechanics simulation from other process
  ! models.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! 

  use Option_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module

  implicit none

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"

  class(pmc_geomechanics_type) :: this

  PetscErrorCode :: ierr

  select type(pmc => this)
    class is(pmc_geomechanics_type)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_pres, &
                           pmc%geomech_realization%geomech_field%press, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                         pmc%sim_aux%subsurf_pres, &
                         pmc%geomech_realization%geomech_field%press, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call VecScatterBegin(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(pmc%sim_aux%subsurf_to_geomechanics, &
                           pmc%sim_aux%subsurf_temp, &
                           pmc%geomech_realization%geomech_field%temp, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%press, &
                            pmc%geomech_realization%geomech_field%press_loc, &
                            ONEDOF)

      call GeomechDiscretizationGlobalToLocal( &
                            pmc%geomech_realization%geomech_discretization, &
                            pmc%geomech_realization%geomech_field%temp, &
                            pmc%geomech_realization%geomech_field%temp_loc, &
                            ONEDOF)

  end select

end subroutine PMCGeomechanicsGetAuxData

! ************************************************************************** !

subroutine PMCGeomechanicsStrip(this)
  !
  ! Deallocates members of PMC Geomechanics.
  !
  ! Author: Satish Karra
  ! Date: 06/01/16
  
  implicit none
  
  class(pmc_geomechanics_type) :: this

  call PMCBaseStrip(this)
  ! realizations destroyed elsewhere
  nullify(this%subsurf_realization)
  nullify(this%geomech_realization)

end subroutine PMCGeomechanicsStrip

! ************************************************************************** !

recursive subroutine PMCGeomechanicsDestroy(this)
  ! 
  ! Author: Satish Karra
  ! Date: 06/01/16
  ! 
  use Option_module

  implicit none
  
  class(pmc_geomechanics_type) :: this
  
#ifdef DEBUG
  call printMsg(this%option,'PMCGeomechanics%Destroy()')
#endif

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif
  
  call PMCGeomechanicsStrip(this)

end subroutine PMCGeomechanicsDestroy

end module PMC_Geomechanics_class
