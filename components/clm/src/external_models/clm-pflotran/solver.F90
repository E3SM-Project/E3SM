module Solver_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
 
#include "petsc/finclude/petscsys.h"

#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petscts.h"
! If the PETSc release is 3.3 or lower, then include petscpcmg.h.
! If using an older version of petsc-dev and petscpcmg.h is required, 
! it can be used by having the makefile turn on HAVE_PETSCPCMG_H.
#if (((PETSC_VERSION_RELEASE) && ((PETSC_VERSION_MAJOR<3) || ((PETSC_VERSION_MAJOR==3) && (PETSC_VERSION_MINOR<=3)))) || (HAVE_PETSCPCMG_H))
#include "petsc/finclude/petscpcmg.h"
#endif

!#include "petsc/finclude/petscpcmg.h"

  type, public :: solver_type
    PetscInt :: itype            ! type: flow or transport
    PetscReal :: linear_atol       ! absolute tolerance
    PetscReal :: linear_rtol       ! relative tolerance
    PetscReal :: linear_dtol       ! divergence tolerance
    PetscInt :: linear_max_iterations     ! maximum number of iterations
    PetscReal :: linear_zero_pivot_tol  ! zero pivot tolerance for LU
    PetscBool :: linear_stop_on_failure ! flag determines whether the code is
                                        ! killed when the solver fails, as
                                        ! opposed ot cutting the step.
    PetscReal :: newton_atol       ! absolute tolerance
    PetscReal :: newton_rtol       ! relative tolerance
    PetscReal :: newton_stol       ! relative tolerance (relative to previous iteration)
    PetscReal :: newton_dtol       ! divergence tolerance
    PetscReal :: newton_inf_res_tol    ! infinity tolerance for residual
    PetscReal :: newton_inf_upd_tol    ! infinity tolerance for update
    PetscReal :: newton_inf_rel_update_tol ! infinity norm on relative update (c(i)-c(i-1))/c(i-1)
    PetscReal :: newton_inf_scaled_res_tol ! infinity norm on scale residual (r(i)/accum(i))
    PetscReal :: newton_inf_res_tol_sec  ! infinity tolerance for secondary continuum residual
    PetscInt :: newton_max_iterations     ! maximum number of iterations
    PetscInt :: newton_min_iterations     ! minimum number of iterations
    PetscInt :: newton_maxf      ! maximum number of function evaluations
    PetscReal :: max_norm          ! maximum norm for divergence
    PetscBool :: use_galerkin_mg  ! If true, precondition linear systems with 
                                   ! Galerkin-type geometric multigrid.
    PetscInt :: galerkin_mg_levels  ! Number of discretization levels for 
                                    ! the Galerkin MG (includes finest level).
    PetscInt :: galerkin_mg_levels_x
    PetscInt :: galerkin_mg_levels_y
    PetscInt :: galerkin_mg_levels_z

    ! Jacobian matrix
    Mat :: J    ! Jacobian
    Mat :: Jpre ! Jacobian to be used in preconditioner
    MatType :: J_mat_type
    MatType :: Jpre_mat_type

    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

    Mat, pointer :: interpolation(:)
      ! Hierarchy of interpolation operators for Galerkin multigrid.

    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType :: pc_type
    KSP ::  ksp
    PC ::  pc
    TS :: ts
    
    PetscBool :: inexact_newton

    PetscBool :: print_convergence
    PetscBool :: print_detailed_convergence
    PetscBool :: print_linear_iterations
    PetscBool :: check_infinity_norm
    PetscBool :: print_ekg
            
  end type solver_type
  
  public :: SolverCreate, &
            SolverDestroy, &
            SolverReadLinear, &
            SolverReadNewton, &
            SolverCreateSNES, &
            SolverSetSNESOptions, &
            SolverCreateTS, &
            SolverPrintNewtonInfo, &
            SolverPrintLinearInfo, &
            SolverCheckCommandLine, &
            SolverLinearPrintFailedReason
  
contains

! ************************************************************************** !

function SolverCreate()
  ! 
  ! Allocates and initializes a new (empty) Solver object
  ! Note that this does not create the PETSc solver contexts associated
  ! with the Solver.  These contexts are created via a subsequent call to
  ! SolverCreateSNES().
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(solver_type), pointer :: SolverCreate
  
  type(solver_type), pointer :: solver
  
  allocate(solver)
  
  ! initialize to default values
  solver%itype = NULL_CLASS
  solver%linear_atol = PETSC_DEFAULT_REAL
  solver%linear_rtol = PETSC_DEFAULT_REAL
  solver%linear_dtol = PETSC_DEFAULT_REAL
  solver%linear_max_iterations = PETSC_DEFAULT_INTEGER
  solver%linear_zero_pivot_tol = UNINITIALIZED_DOUBLE
  solver%linear_stop_on_failure = PETSC_FALSE
  
  solver%newton_atol = PETSC_DEFAULT_REAL
  solver%newton_rtol = PETSC_DEFAULT_REAL
  solver%newton_stol = PETSC_DEFAULT_REAL
  solver%newton_dtol = PETSC_DEFAULT_REAL 
  solver%max_norm = 1.d20     ! set to a large value
  solver%newton_inf_res_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_upd_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_rel_update_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_scaled_res_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_res_tol_sec = 1.d-10
  solver%newton_max_iterations = PETSC_DEFAULT_INTEGER
  solver%newton_min_iterations = 1
  solver%newton_maxf = PETSC_DEFAULT_INTEGER

  solver%use_galerkin_mg = PETSC_FALSE
  solver%galerkin_mg_levels = 1
  solver%galerkin_mg_levels_x = 1
  solver%galerkin_mg_levels_y = 1
  solver%galerkin_mg_levels_z = 1
  
  solver%J = 0
  solver%Jpre = 0
  solver%J_mat_type = MATBAIJ
  solver%Jpre_mat_type = ''
!  solver%interpolation = 0
  nullify(solver%interpolation)
  solver%matfdcoloring = 0
  solver%snes = 0
  solver%ksp_type = KSPBCGS
  solver%pc_type = ""
  solver%ksp = 0
  solver%pc = 0
  solver%ts = 0
  
  solver%inexact_newton = PETSC_FALSE
  
  solver%print_convergence = PETSC_TRUE
  solver%print_detailed_convergence = PETSC_FALSE
  solver%print_linear_iterations = PETSC_FALSE
  solver%check_infinity_norm = PETSC_TRUE
  solver%print_ekg = PETSC_FALSE
    
  SolverCreate => solver
  
end function SolverCreate

! ************************************************************************** !

subroutine SolverCreateSNES(solver,comm)
  ! 
  ! Create PETSc SNES object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call SNESCreate(comm,solver%snes,ierr);CHKERRQ(ierr)
  call SNESSetFromOptions(solver%snes,ierr);CHKERRQ(ierr)

  ! grab handles for ksp and pc
  call SNESGetKSP(solver%snes,solver%ksp,ierr);CHKERRQ(ierr)
  call KSPGetPC(solver%ksp,solver%pc,ierr);CHKERRQ(ierr)

end subroutine SolverCreateSNES

! ************************************************************************** !

subroutine SolverSetSNESOptions(solver, option)
  ! 
  ! Sets options for SNES
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  SNESLineSearch :: linesearch
  KSP, pointer :: sub_ksps(:)
  PC :: pc
  PetscInt :: nsub_ksp
  PetscInt :: first_sub_ksp
  PetscErrorCode :: ierr
  PetscInt :: i
  
  ! if ksp_type or pc_type specified in input file, set them here
  if (len_trim(solver%ksp_type) > 1) then
    call KSPSetType(solver%ksp,solver%ksp_type,ierr);CHKERRQ(ierr)
  endif
  if (len_trim(solver%pc_type) > 1) then
    call PCSetType(solver%pc,solver%pc_type,ierr);CHKERRQ(ierr)
  endif

  call KSPSetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                        solver%linear_dtol,solver%linear_max_iterations, &
                        ierr);CHKERRQ(ierr)
  ! as of PETSc 3.7, we need to turn on error reporting due to zero pivots
  ! as PETSc no longer reports zero pivots for very small concentrations
  !geh: this gets overwritten by ksp->errorifnotconverted
  if (solver%linear_stop_on_failure) then
    call KSPSetErrorIfNotConverged(solver%ksp,PETSC_TRUE,ierr);CHKERRQ(ierr)
  endif

  ! allow override from command line
  call KSPSetFromOptions(solver%ksp,ierr);CHKERRQ(ierr)
  call PCSetFromOptions(solver%pc,ierr);CHKERRQ(ierr)
    
  ! get the ksp_type and pc_type incase of command line override.
  call KSPGetType(solver%ksp,solver%ksp_type,ierr);CHKERRQ(ierr)
  call PCGetType(solver%pc,solver%pc_type,ierr);CHKERRQ(ierr)

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(solver%snes, solver%newton_atol, solver%newton_rtol, &
                         solver%newton_stol,solver%newton_max_iterations, &
                         solver%newton_maxf,ierr);CHKERRQ(ierr)

  ! set inexact newton, currently applies default settings
  if (solver%inexact_newton) then
    call SNESKSPSetUseEW(solver%snes,PETSC_TRUE,ierr);CHKERRQ(ierr)
  endif

!  call SNESLineSearchSet(solver%snes,SNESLineSearchNo,PETSC_NULL)

  ! Setup for n-level Galerkin multigrid.
  if (solver%use_galerkin_mg) then
    call PCSetType(solver%pc, PCMG,ierr);CHKERRQ(ierr)
    call PCMGSetLevels(solver%pc, solver%galerkin_mg_levels, &
                       PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
    do i=1,solver%galerkin_mg_levels-1
      call PCMGSetInterpolation(solver%pc, i, solver%interpolation(i), &
                                ierr);CHKERRQ(ierr)
      call PCMGSetGalerkin(solver%pc,ierr);CHKERRQ(ierr)
    enddo
  endif
  
  ! allow override from command line; for some reason must come before
  ! LineSearchParams, or they crash
  call SNESSetFromOptions(solver%snes,ierr);CHKERRQ(ierr)

  ! the below must come after SNESSetFromOptions
  ! PETSc no longer performs a shift on matrix diagonals by default.  We 
  ! force the shift since it helps alleviate zero pivots.
  call PCFactorSetShiftType(solver%pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
  if (solver%pc_type == PCBJACOBI) then
    call KSPSetup(solver%ksp,ierr);CHKERRQ(ierr)
    call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                            PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
    allocate(sub_ksps(nsub_ksp))
    sub_ksps = 0
    call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                            sub_ksps,ierr);CHKERRQ(ierr)
    do i = 1, nsub_ksp
      call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
      call PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
    enddo
    deallocate(sub_ksps)
    nullify(sub_ksps)
  elseif (.not.(solver%pc_type == PCLU .or. solver%pc_type == PCILU)) then
    option%io_buffer = 'PCFactorShiftType for PC ' // &
      trim(solver%pc_type) // ' is not supported at this time.'
    call printErrMsg(option)
  endif
  
  if (Initialized(solver%linear_zero_pivot_tol)) then
    call PCFactorSetZeroPivot(solver%pc,solver%linear_zero_pivot_tol, &
                              ierr);CHKERRQ(ierr)
    if (solver%pc_type == PCBJACOBI) then
      call KSPSetup(solver%ksp,ierr);CHKERRQ(ierr)
      call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
      allocate(sub_ksps(nsub_ksp))
      sub_ksps = 0
      call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              sub_ksps,ierr);CHKERRQ(ierr)
      do i = 1, nsub_ksp
        call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
        call PCFactorSetZeroPivot(pc,solver%linear_zero_pivot_tol, &
                                  ierr);CHKERRQ(ierr)
      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    elseif (.not.(solver%pc_type == PCLU .or. solver%pc_type == PCILU)) then
      option%io_buffer = 'PCFactorSetZeroPivot for PC ' // &
        trim(solver%pc_type) // ' is not supported at this time.'
      call printErrMsg(option)
    endif
  endif

  call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetTolerances(linesearch, solver%newton_stol,       &
          PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_INTEGER, ierr);CHKERRQ(ierr)

  call SNESGetTolerances(solver%snes,solver%newton_atol,solver%newton_rtol, &
                         solver%newton_stol,solver%newton_max_iterations, &
                         solver%newton_maxf,ierr);CHKERRQ(ierr)

  call KSPGetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                         solver%linear_dtol,solver%linear_max_iterations, &
                        ierr);CHKERRQ(ierr)

end subroutine SolverSetSNESOptions

! ************************************************************************** !

subroutine SolverCreateTS(solver,comm)
  ! 
  ! This routine creates PETSc TS object.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 01/18/13
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call TSCreate(comm,solver%ts,ierr);CHKERRQ(ierr)
  call TSSetFromOptions(solver%ts,ierr);CHKERRQ(ierr)

end subroutine SolverCreateTS

! ************************************************************************** !

subroutine SolverReadLinear(solver,input,option)
  ! 
  ! Reads parameters associated with linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscErrorCode :: ierr
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2, prefix
  character(len=MAXSTRINGLENGTH) :: string

  select case(solver%itype)
    case(FLOW_CLASS)
      prefix = '-flow_'
    case(TRANSPORT_CLASS)
      prefix = '-tran_'
  end select

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','LINEAR SOLVER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('SOLVER_TYPE','SOLVER','KRYLOV_TYPE','KRYLOV','KSP','KSP_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'ksp_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PREONLY')
            solver%ksp_type = KSPPREONLY
          case('GMRES')
            solver%ksp_type = KSPGMRES
          case('FGMRES')
            solver%ksp_type = KSPFGMRES
          case('BCGS','BICGSTAB','BI-CGSTAB')
            solver%ksp_type = KSPBCGS
          case('IBCGS','IBICGSTAB','IBI-CGSTAB')
            solver%ksp_type = KSPIBCGS
          case('RICHARDSON')
            solver%ksp_type = KSPRICHARDSON
          case('CG')
            solver%ksp_type = KSPCG
          case('DIRECT')
            solver%ksp_type = KSPPREONLY
            solver%pc_type = PCLU
          case('ITERATIVE','KRYLOV')
            solver%ksp_type = KSPBCGS
            solver%pc_type = PCBJACOBI
          case default
            option%io_buffer  = 'Krylov solver type: ' // trim(word) // &
                                ' unknown.'
            call printErrMsg(option)
        end select

      case('PRECONDITIONER_TYPE','PRECONDITIONER','PC','PC_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'pc_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PCNONE')
            solver%pc_type = PCNONE
          case('ILU','PCILU')
            solver%pc_type = PCILU
          case('LU','PCLU')
            solver%pc_type = PCLU
          case('BJACOBI','BLOCK_JACOBI')
            solver%pc_type = PCBJACOBI
          case('JACOBI')
            solver%pc_type = PCJACOBI
          case('ASM','ADDITIVE_SCHWARZ')
            solver%pc_type = PCASM
          case('HYPRE')
            solver%pc_type = PCHYPRE
          case('SHELL')
            solver%pc_type = PCSHELL
          case default
            option%io_buffer  = 'Preconditioner type: ' // trim(word) // &
                                ' unknown.'
            call printErrMsg(option)
        end select

      case('HYPRE_OPTIONS')
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          call InputReadWord(input,option,keyword,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword', &
                             'LINEAR SOLVER, HYPRE options')   
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'type', &
                                 'LINEAR SOLVER, HYPRE options')  
              call StringToLower(word)
              select case(trim(word))
                case('pilut','parasails','boomeramg','euclid')
                  string = trim(prefix) // 'pc_hypre_type'
                  call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                            trim(string),trim(word), &
                                            ierr);CHKERRQ(ierr)
                case default
                  option%io_buffer  = 'HYPRE preconditioner type: ' // &
                                      trim(word) // ' unknown.'
                  call printErrMsg(option)
              end select
            case('BOOMERAMG_CYCLE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG cycle type', &
                                 'LINEAR SOLVER, HYPRE options')  
              call StringToLower(word)
              string = trim(prefix) // 'pc_hypre_boomeramg_cycle_type'
              select case(trim(word))
                case('V')
                  call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                            trim(string),'1', &
                                            ierr);CHKERRQ(ierr)
                case('W')
                  call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                            trim(string),'2', &
                                            ierr);CHKERRQ(ierr)
                case default
                  option%io_buffer  = 'HYPRE BoomerAMG cycle type: ' &
                                      // trim(word) // ' unknown.'
                  call printErrMsg(option)
              end select
            case('BOOMERAMG_MAX_LEVELS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_levels'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MAX_ITER')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum iterations', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_iter'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TOL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG convergence tolerance', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_tol'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TRUNCFACTOR')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG interpolation truncation factor', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_truncfactor'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG # levels aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_nl'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NUM_PATHS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG # paths for aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_num_paths'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_STRONG_THRESHOLD')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG threshold for strong connectivity', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_strong_threshold'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                         'BoomerAMG number of grid sweeps up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_all'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG number of grid sweeps down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_down'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG number of grid sweeps up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_up'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG number of grid sweeps for coarse level', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG relaxation type for up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_all'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG relaxation type for down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_down'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_up'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for coarse grids', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG outer relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) //  &
                       'pc_hypre_boomeramg_outer_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                              'BoomerAMG outer relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // &
                       'pc_hypre_boomeramg_outer_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NO_CF')
              string = trim(prefix) // 'pc_hypre_boomeramg_no_CF'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MEASURE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG measure type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_measure_type'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_COARSEN_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG coarsen type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_coarsen_type'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_INTERPOLATION_TYPE','BOOMERAMG_INTERP_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG interpolation type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_interp_type'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_COARSEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG set nodal coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_coarsen'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_RELAXATION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG nodal relaxation via Schwarz', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_relaxation'
              call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case default
              option%io_buffer  = 'HYPRE option: ' // trim(keyword) // &
                                  ' unknown.'
              call printErrMsg(option)
          end select
        enddo

      case('ATOL')
        call InputReadDouble(input,option,solver%linear_atol)
        call InputErrorMsg(input,option,'linear_atol','LINEAR_SOLVER')

      case('RTOL')
        call InputReadDouble(input,option,solver%linear_rtol)
        call InputErrorMsg(input,option,'linear_rtol','LINEAR_SOLVER')

      case('DTOL')
        call InputReadDouble(input,option,solver%linear_dtol)
        call InputErrorMsg(input,option,'linear_dtol','LINEAR_SOLVER')
   
      case('MAXIT')
        call InputReadInt(input,option,solver%linear_max_iterations)
        call InputErrorMsg(input,option,'linear_max_iterations','LINEAR_SOLVER')

      case('ZERO_PIVOT_TOL','LU_ZERO_PIVOT_TOL')
        call InputReadDouble(input,option,solver%linear_zero_pivot_tol)
        call InputErrorMsg(input,option,'linear_zero_pivot_tol', &
                           'LINEAR_SOLVER')

      case('STOP_ON_FAILURE')
        solver%linear_stop_on_failure = PETSC_TRUE

      case('MUMPS')
        string = trim(prefix) // 'pc_factor_mat_solver_package'
        word = 'mumps'
        call PetscOptionsSetValue(PETSC_NULL_OBJECT, &
                                  trim(string),trim(word),ierr);CHKERRQ(ierr)
   
      case default
        call InputKeywordUnrecognized(keyword,'LINEAR_SOLVER',option)
    end select 
  
  enddo  

end subroutine SolverReadLinear

! ************************************************************************** !

subroutine SolverReadNewton(solver,input,option)
  ! 
  ! Reads parameters associated with linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2

  input%ierr = 0
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','NEWTON SOLVER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case ('INEXACT_NEWTON')
        solver%inexact_newton = PETSC_TRUE

      case ('NO_PRINT_CONVERGENCE')
        solver%print_convergence = PETSC_FALSE

      case ('NO_INF_NORM','NO_INFINITY_NORM')
        solver%check_infinity_norm = PETSC_FALSE

      case('MAXIMUM_NEWTON_ITERATIONS')
        call InputReadInt(input,option,solver%newton_max_iterations)
        call InputErrorMsg(input,option,'maximum newton iterations', &
                           'NEWTON_SOLVER')

      case('MINIMUM_NEWTON_ITERATIONS')
        call InputReadInt(input,option,solver%newton_min_iterations)
        call InputErrorMsg(input,option,'minimum newton iterations', &
                           'NEWTON_SOLVER')

      case ('PRINT_DETAILED_CONVERGENCE')
        solver%print_detailed_convergence = PETSC_TRUE

      case ('PRINT_LINEAR_ITERATIONS')
        solver%print_linear_iterations = PETSC_TRUE

      case('ATOL')
        call InputReadDouble(input,option,solver%newton_atol)
        call InputErrorMsg(input,option,'newton_atol','NEWTON_SOLVER')

      case('RTOL')
        call InputReadDouble(input,option,solver%newton_rtol)
        call InputErrorMsg(input,option,'newton_rtol','NEWTON_SOLVER')

      case('STOL')
        call InputReadDouble(input,option,solver%newton_stol)
        call InputErrorMsg(input,option,'newton_stol','NEWTON_SOLVER')
      
      case('DTOL')
        call InputReadDouble(input,option,solver%newton_dtol)
        call InputErrorMsg(input,option,'newton_dtol','NEWTON_SOLVER')

      case('MAX_NORM')
        call InputReadDouble(input,option,solver%max_norm)
        call InputErrorMsg(input,option,'max_norm','NEWTON_SOLVER')
   
      case('ITOL', 'INF_TOL', 'ITOL_RES', 'INF_TOL_RES')
        call InputReadDouble(input,option,solver%newton_inf_res_tol)
        call InputErrorMsg(input,option,'newton_inf_res_tol','NEWTON_SOLVER')
   
      case('ITOL_UPDATE', 'INF_TOL_UPDATE')
        call InputReadDouble(input,option,solver%newton_inf_upd_tol)
        call InputErrorMsg(input,option,'newton_inf_upd_tol','NEWTON_SOLVER')

      case('ITOL_SCALED_RESIDUAL')
        option%io_buffer = 'Flow NEWTON_SOLVER ITOL_SCALED_RESIDUAL is ' // &
          'now specific to each process model and must be defined in ' // &
          'the SIMULATION/PROCESS_MODELS/SUBSURFACE_FLOW/OPTIONS block.'
        call printErrMsg(option)
          
      case('ITOL_RELATIVE_UPDATE')
        option%io_buffer = 'Flow NEWTON_SOLVER ITOL_RELATIVE_UPDATE is ' // &
          'now specific to each process model and must be defined in ' // &
          'the SIMULATION/PROCESS_MODELS/SUBSURFACE_FLOW/OPTIONS block.'
        call printErrMsg(option)

      case('ITOL_SEC','ITOL_RES_SEC','INF_TOL_SEC')
        if (.not.option%use_mc) then
          option%io_buffer = 'NEWTON ITOL_SEC not supported without ' // &
            'MULTIPLE_CONTINUUM keyword.'
          call printErrMsg(option)
        endif
        if (.not.solver%itype == TRANSPORT_CLASS) then
          option%io_buffer = 'NEWTON ITOL_SEC supported in ' // &
            'TRANSPORT only.'
          call printErrMsg(option)        
        endif         
        call InputReadDouble(input,option,solver%newton_inf_res_tol_sec)
        call InputErrorMsg(input,option,'newton_inf_res_tol_sec', &
                           'NEWTON_SOLVER')
   
      case('MAXIT')
        call InputReadInt(input,option,solver%newton_max_iterations)
        call InputErrorMsg(input,option,'maximum newton iterations', &
                           'NEWTON_SOLVER')

      case('MAXF')
        call InputReadInt(input,option,solver%newton_maxf)
        call InputErrorMsg(input,option,'newton_maxf','NEWTON_SOLVER')

      case('MATRIX_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('BAIJ')
            solver%J_mat_type = MATBAIJ
          case('AIJ')
!           solver%J_mat_type = MATBAIJ
            solver%J_mat_type = MATAIJ
          case('MFFD','MATRIX_FREE')
            solver%J_mat_type = MATMFFD
          case('HYPRESTRUCT')
            solver%J_mat_type = MATHYPRESTRUCT
          case default
            option%io_buffer = 'Matrix type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select
        
      case('PRECONDITIONER_MATRIX_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('BAIJ')
            solver%Jpre_mat_type = MATBAIJ
          case('AIJ')
!           solver%Jpre_mat_type = MATBAIJ
            solver%Jpre_mat_type = MATAIJ
          case('MFFD','MATRIX_FREE')
            solver%Jpre_mat_type = MATMFFD
          case('HYPRESTRUCT')
             solver%Jpre_mat_type = MATHYPRESTRUCT
          case('SHELL')
             solver%Jpre_mat_type = MATSHELL
          case default
            option%io_buffer  = 'Preconditioner Matrix type: ' // trim(word) // ' unknown.'
            call printErrMsg(option)
        end select
        
      case default
        call InputKeywordUnrecognized(keyword,'NEWTON_SOLVER',option)
    end select 
  
  enddo  

end subroutine SolverReadNewton

! ************************************************************************** !

subroutine SolverPrintLinearInfo(solver,header,option)
  ! 
  ! Prints information about linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  
  implicit none
  
  type(solver_type) :: solver
  character(len=*) :: header  
  type(option_type) :: option

  PetscInt :: fid

#if !defined(PETSC_HAVE_MUMPS)
  if (option%mycommsize > 1) then
    if (solver%ksp_type == KSPPREONLY .and. solver%pc_type == PCLU) then
      option%io_buffer = 'Direct solver (KSPPREONLY + PCLU) not ' // &
        ' supported when running in parallel.  Switch to SOLVER ITERATIVE.'
      call printErrMsg(option)
    endif
  endif
#endif
  
  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header) // ' Linear Solver'
    write(*,'("   solver:  ",a)') trim(solver%ksp_type)
    write(*,'("  precond:  ",a)') trim(solver%pc_type)
    write(*,'("     atol:",1pe12.4)') solver%linear_atol
    write(*,'("     rtol:",1pe12.4)') solver%linear_rtol
    write(*,'("     dtol:",1pe12.4)') solver%linear_dtol
    write(*,'(" max iter:",i7)') solver%linear_max_iterations
    if (Initialized(solver%linear_zero_pivot_tol)) then
      write(*,'("pivot tol:",1pe12.4)') solver%linear_zero_pivot_tol
    endif
  endif
  
  if (OptionPrintToFile(option)) then
    fid = option%fid_out
    write(fid,*) 
    write(fid,'(a)') trim(header) // ' Linear Solver'
    write(fid,'("   solver:  ",a)') trim(solver%ksp_type)
    write(fid,'("  precond:  ",a)') trim(solver%pc_type)
    write(fid,'("     atol:",1pe12.4)') solver%linear_atol
    write(fid,'("     rtol:",1pe12.4)') solver%linear_rtol
    write(fid,'("     dtol:",1pe12.4)') solver%linear_dtol
    write(fid,'(" max iter:",i7)') solver%linear_max_iterations
    if (Initialized(solver%linear_zero_pivot_tol)) then
      write(fid,'("pivot tol:",1pe12.4)') solver%linear_zero_pivot_tol
    endif
  endif

end subroutine SolverPrintLinearInfo

! ************************************************************************** !

subroutine SolverPrintNewtonInfo(solver,header,option)    
  ! 
  ! Prints information about Newton solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  character(len=*) :: header  
  type(option_type) :: option
  PetscInt :: fid

  if (OptionPrintToScreen(option)) then
    write(*,*) 
    write(*,'(a)') trim(header) // ' Newton Solver'
    write(*,'("        atol:",1pe12.4)') solver%newton_atol
    write(*,'("        rtol:",1pe12.4)') solver%newton_rtol
    write(*,'("        stol:",1pe12.4)') solver%newton_stol
    write(*,'("        dtol:",1pe12.4)') solver%newton_dtol
    write(*,'("     maxnorm:",1pe12.4)') solver%max_norm
    write(*,'("   inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(*,'("   inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(*,'("inftolrelupd:",1pe12.4)') solver%newton_inf_rel_update_tol
    write(*,'("inftolsclres:",1pe12.4)') solver%newton_inf_scaled_res_tol
    write(*,'("    max iter:",i6)') solver%newton_max_iterations
    write(*,'("    min iter:",i6)') solver%newton_min_iterations
    write(*,'("        maxf:",i6)') solver%newton_maxf
    write(*,*) 
    if (len_trim(solver%J_mat_type) > 2) then
      write(*,'("matrix type:",a20)') solver%J_mat_type
    endif
    if (len_trim(solver%Jpre_mat_type) > 2) then
      write(*,'("precond. matrix type:",a20)') solver%Jpre_mat_type
    endif
    if (solver%inexact_newton) then
      write(*,'("inexact newton: on")')
    else
      write(*,'("inexact newton: off")')
    endif
        
    if (solver%print_convergence) then
      write(*,'("print convergence: on")')
    else
      write(*,'("print convergence: off")')
    endif
        
    if (solver%print_detailed_convergence) then
      write(*,'("print detailed convergence: on")')
    else
      write(*,'("print detailed convergence: off")')
    endif
        
    if (solver%check_infinity_norm) then
      write(*,'("check infinity norm: on")')
    else
      write(*,'("check infinity norm: off")')
    endif
  endif

  if (OptionPrintToFile(option)) then
    fid = option%fid_out
    write(fid,*) 
    write(fid,'(a)') trim(header) // ' Newton Solver'
    write(fid,'("        atol:",1pe12.4)') solver%newton_atol
    write(fid,'("        rtol:",1pe12.4)') solver%newton_rtol
    write(fid,'("        stol:",1pe12.4)') solver%newton_stol
    write(fid,'("        dtol:",1pe12.4)') solver%newton_dtol
    write(fid,'("     maxnorm:",1pe12.4)') solver%max_norm
    write(fid,'("   inftolres:",1pe12.4)') solver%newton_inf_res_tol
    write(fid,'("   inftolupd:",1pe12.4)') solver%newton_inf_upd_tol
    write(fid,'("inftolrelupd:",1pe12.4)') solver%newton_inf_rel_update_tol
    write(fid,'("inftolsclres:",1pe12.4)') solver%newton_inf_scaled_res_tol
    write(fid,'("    max iter:",i6)') solver%newton_max_iterations
    write(fid,'("    min iter:",i6)') solver%newton_min_iterations
    write(fid,'("        maxf:",i6)') solver%newton_maxf
    write(fid,*) 
    if (len_trim(solver%J_mat_type) > 2) then
      write(fid,'("matrix type:",a20)') solver%J_mat_type
    endif
    if (len_trim(solver%Jpre_mat_type) > 2) then
      write(fid,'("precond. matrix type:",a20)') solver%Jpre_mat_type
    endif
    if (solver%inexact_newton) then
      write(fid,'("inexact newton: on")')
    else
      write(fid,'("inexact newton: off")')
    endif
        
    if (solver%print_convergence) then
      write(fid,'("print convergence: on")')
    else
      write(fid,'("print convergence: off")')
    endif
        
    if (solver%print_detailed_convergence) then
      write(fid,'("print detailed convergence: on")')
    else
      write(fid,'("print detailed convergence: off")')
    endif
        
    if (solver%check_infinity_norm) then
      write(fid,'("check infinity norm: on")')
    else
      write(fid,'("check infinity norm: off")')
    endif
  endif

end subroutine SolverPrintNewtonInfo

! ************************************************************************** !

subroutine SolverCheckCommandLine(solver)
  ! 
  ! Parses the command line for various solver
  ! options.
  ! Note: In order to use the PETSc OptionsPrefix associated with
  ! solver%snes in parsing the options, the call to SolverCheckCommandLine()
  ! should come after the SNESSetOptionsPrefix(solver%snes,...) call.
  ! 
  ! Author: Richard Tran Mills
  ! Date: 05/09/2008
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: prefix
  character(len=MAXSTRINGLENGTH) :: mat_type
  PetscBool :: is_present

  if (solver%snes /= 0) then
    call SNESGetOptionsPrefix(solver%snes, prefix, ierr);CHKERRQ(ierr)
  else
    prefix = PETSC_NULL_CHARACTER
  endif

  ! Parse the options to determine if the matrix type has been specified.
  call PetscOptionsGetString(PETSC_NULL_OBJECT,prefix, '-mat_type', mat_type, &
                             is_present,ierr);CHKERRQ(ierr)
  if (is_present) solver%J_mat_type = trim(mat_type)
  
  call PetscOptionsGetString(PETSC_NULL_OBJECT,prefix, '-pre_mat_type', &
                             mat_type, is_present,ierr);CHKERRQ(ierr)
  if (is_present) solver%Jpre_mat_type = trim(mat_type)

  ! Parse the options for the Galerkin multigrid solver.
  ! Users can specify the number of levels of coarsening via the
  ! 'galerkin_mg N' option, which will set the number of levels in the 
  ! x, y, and z directions all to N.  For semi-coarsening, however, 
  ! it is possible to set the number of levels in each direction 
  ! individually via options such as '-galerkin_mg_x N', which would 
  ! override the number of levels in the x direction set by '-galerkin_mg'.
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,prefix, '-galerkin_mg', &
                          solver%galerkin_mg_levels, solver%use_galerkin_mg, &
                          ierr);CHKERRQ(ierr)
  if (solver%use_galerkin_mg) then
    solver%galerkin_mg_levels_x = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_y = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_z = solver%galerkin_mg_levels
  endif

  call PetscOptionsGetInt(PETSC_NULL_OBJECT,prefix, '-galerkin_mg_x', &
                          solver%galerkin_mg_levels_x, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,prefix, '-galerkin_mg_y', &
                          solver%galerkin_mg_levels_y, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(PETSC_NULL_OBJECT,prefix, '-galerkin_mg_z', &
                          solver%galerkin_mg_levels_z, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE

  if (solver%use_galerkin_mg) then
    solver%J_mat_type = MATAIJ
      ! Must use AIJ above, as BAIJ is not supported for Galerkin MG solver.
    solver%galerkin_mg_levels = max(solver%galerkin_mg_levels_x, &
                                    solver%galerkin_mg_levels_y, &
                                    solver%galerkin_mg_levels_z)
  endif
                             

end subroutine SolverCheckCommandLine

! ************************************************************************** !

subroutine SolverLinearPrintFailedReason(solver,option)    
  ! 
  ! Prints the reason for the solver failing
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/16
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  KSP, pointer :: sub_ksps(:)
  PC :: pc
  Mat :: mat
  PCType :: pc_type
  PetscInt :: i
  PetscInt :: nsub_ksp
  PetscInt :: first_sub_ksp
  KSPConvergedReason :: ksp_reason
  PCFailedReason :: pc_failed_reason, global_pc_failed_reason
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: zero_pivot_tol, zero_pivot
  character(len=MAXWORDLENGTH) :: word, word2
  PetscInt :: irow, temp_int
  PetscErrorCode :: ierr

  call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr);CHKERRQ(ierr)
  select case(ksp_reason)
    case(KSP_DIVERGED_ITS)
      option%io_buffer = ' -> KSPReason: Diverged due to iterations'
    case(KSP_DIVERGED_DTOL)
      option%io_buffer = ' -> KSPReason: Diverged due to dtol'
    case(KSP_DIVERGED_BREAKDOWN)
      option%io_buffer = ' -> KSPReason: Diverged due to breakdown'
    case(KSP_DIVERGED_BREAKDOWN_BICG)
      option%io_buffer = ' -> KSPReason: Diverged due to breakdown bicg'
    case(KSP_DIVERGED_NONSYMMETRIC)
      option%io_buffer = ' -> KSPReason: Diverged due to nonsymmetric'
    case(KSP_DIVERGED_INDEFINITE_PC)
      option%io_buffer = ' -> KSPReason: Diverged due to indefinite PC'
    case(KSP_DIVERGED_NANORINF)
      option%io_buffer = ' -> KSPReason: Diverged due to NaN or Inf PC'
    case(KSP_DIVERGED_INDEFINITE_MAT)
      option%io_buffer = ' -> KSPReason: Diverged due to indefinite matix'
    case(KSP_DIVERGED_PCSETUP_FAILED)
      option%io_buffer = ' -> KSPReason: Diverged due to PC setup failed'
      pc = solver%pc
      call PCGetType(pc,pc_type,ierr);CHKERRQ(ierr)
      call PCGetSetUpFailedReason(pc,pc_failed_reason, &
                                  ierr);CHKERRQ(ierr)
      ! have to perform global reduction on pc_failed_reason
      temp_int = pc_failed_reason
      call MPI_Allreduce(MPI_IN_PLACE,temp_int,ONE_INTEGER_MPI, &
                         MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
      global_pc_failed_reason = temp_int
      if (global_pc_failed_reason == PC_SUBPC_ERROR) then
        if (pc_type == PCBJACOBI) then
          call PCBJacobiGetSubKSP(pc,nsub_ksp,first_sub_ksp, &
                                  PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
          allocate(sub_ksps(nsub_ksp))
          sub_ksps = 0
          call PCBJacobiGetSubKSP(pc,nsub_ksp,first_sub_ksp, &
                                  sub_ksps,ierr);CHKERRQ(ierr)
          if (nsub_ksp > 1) then
            option%io_buffer = 'NSUB_KSP > 1.  What to do?  Email pflotran&
              &-dev@googlegroups.com.'
            call printErrMsg(option)
          endif
          do i = 1, nsub_ksp
            call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
            call PCGetSetUpFailedReason(pc,pc_failed_reason, &
                                        ierr);CHKERRQ(ierr)
          enddo
          deallocate(sub_ksps)
          nullify(sub_ksps)
        else
          option%io_buffer = 'Error in SUB PC of unknown type "' // &
            trim(pc_type) // '".'
          call printErrMsg(option)
        endif
      endif
      ! have to perform global reduction (again) on pc_failed_reason
      temp_int = pc_failed_reason
      call MPI_Allreduce(MPI_IN_PLACE,temp_int,ONE_INTEGER_MPI, &
                         MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
      global_pc_failed_reason = temp_int
      select case(global_pc_failed_reason)
        case(PC_FACTOR_STRUCT_ZEROPIVOT,PC_FACTOR_NUMERIC_ZEROPIVOT)
          select case(solver%itype)
            case(FLOW_CLASS)
              string = 'Flow'
            case(TRANSPORT_CLASS)
              string = 'Transport'
          end select
          call PCFactorGetZeroPivot(pc,zero_pivot_tol, &
                                    ierr);CHKERRQ(ierr)
          write(word,*) zero_pivot_tol
#if PETSC_VERSION_GT(3,7,3)
          ! In parallel, some processes will not have a zero pivot and
          ! will report zero as the error.  We must skip these processes.
          zero_pivot = 1.d20
          ! note that this is not the global pc reason
          select case(pc_failed_reason)
            case(PC_FACTOR_STRUCT_ZEROPIVOT,PC_FACTOR_NUMERIC_ZEROPIVOT)
            call PCFactorGetMatrix(pc,mat,ierr);CHKERRQ(ierr)
            call MatFactorGetErrorZeroPivot(mat,zero_pivot,irow, &
                                            ierr);CHKERRQ(ierr)
          end select
          call MPI_Allreduce(MPI_IN_PLACE,zero_pivot,ONE_INTEGER_MPI, &
                             MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
          write(word2,*) zero_pivot
#endif
          option%io_buffer = 'PC Setup failed for ' // trim(string) // &
            '. The ' // trim(string) // ' preconditioner zero pivot &
            &tolerance (' // trim(adjustl(word)) // &
#if PETSC_VERSION_GT(3,7,3)
            ') is too large due to a zero pivot of ' // &
            trim(adjustl(word2)) // '. Please set a ZERO_PIVOT_TOL smaller &
            &than that value or email pflotran-dev@googlegroups.com with &
            &this information for guildance.'
#else
            ') is too large. Please run PFLOTRAN with STOP_ON_FAILURE &
            &added to the respective SOLVER block to determine the &
            &needed ZERO_PIVOT_TOL in that solver block.'
#endif
          call printErrMsg(option)
      end select
    case default
      write(option%io_buffer,'('' -> KSPReason: Unknown: '',i2)') &
        ksp_reason
  end select
  call printMsg(option)

end subroutine SolverLinearPrintFailedReason

! ************************************************************************** !

subroutine SolverDestroy(solver)
  ! 
  ! Deallocates a solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(solver_type), pointer :: solver
  
  PetscErrorCode :: ierr
  PetscInt :: i

  if (.not.associated(solver)) return

  if (solver%Jpre == solver%J) then
    solver%Jpre = 0
  else if (solver%Jpre /= 0) then
    call MatDestroy(solver%Jpre,ierr);CHKERRQ(ierr)
  endif
  if (solver%J /= 0) then
    call MatDestroy(solver%J,ierr);CHKERRQ(ierr)
  endif
  if (associated(solver%interpolation)) then
    do i=1,solver%galerkin_mg_levels-1
      call MatDestroy(solver%interpolation(i),ierr);CHKERRQ(ierr)
    enddo
    deallocate(solver%interpolation)
  endif
  if (solver%matfdcoloring /= 0) then
    call MatFDColoringDestroy(solver%matfdcoloring,ierr);CHKERRQ(ierr)
  endif

  if (solver%snes /= 0) then
    call SNESDestroy(solver%snes,ierr);CHKERRQ(ierr)
  endif
  if (solver%ts /= 0) then
    call TSDestroy(solver%ts,ierr);CHKERRQ(ierr)
  endif

  solver%ksp = 0
  solver%pc = 0
    
  deallocate(solver)
  nullify(solver)
  
end subroutine SolverDestroy
  
end module Solver_module
