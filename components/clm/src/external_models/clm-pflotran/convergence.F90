module Convergence_module

  use Solver_module
  use Option_module
  use Grid_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscsnes.h"
#include "petsc/finclude/petsclog.h"

  type, public :: convergence_context_type
    type(solver_type), pointer :: solver
    type(option_type), pointer :: option
    type(grid_type), pointer :: grid
  end type convergence_context_type


  public :: ConvergenceContextCreate, ConvergenceTest, &
            ConvergenceContextDestroy
  
contains

! ************************************************************************** !

function ConvergenceContextCreate(solver,option,grid)
  ! 
  ! Creates a context containing pointer
  ! for convergence subroutines
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 

  implicit none
  
  type(convergence_context_type), pointer :: ConvergenceContextCreate
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  
  type(convergence_context_type), pointer :: context
  
  allocate(context)
  context%solver => solver
  context%option => option
  context%grid => grid

  ConvergenceContextCreate => context

end function ConvergenceContextCreate

! ************************************************************************** !

subroutine ConvergenceTest(snes_,i_iteration,xnorm,unorm,fnorm,reason,context, &
                           ierr)
  ! 
  ! User defined convergence test
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 

  implicit none
  
  SNES :: snes_
  PetscInt :: i_iteration
  PetscReal :: xnorm ! 2-norm of updated solution
  PetscReal :: unorm ! 2-norm of update. PETSc refers to this as snorm
  PetscReal :: fnorm ! 2-norm of updated residual
  SNESConvergedReason :: reason
  type(convergence_context_type) :: context
  PetscErrorCode :: ierr
  
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  
  Vec :: solution_vec
  Vec :: update_vec
  Vec :: residual_vec
  PetscReal :: inorm_solution  !infinity norms
  PetscReal :: inorm_update  
  PetscReal :: inorm_residual  
  
  PetscInt :: i, ndof, max_index, min_index
  PetscReal, allocatable :: fnorm_solution_stride(:)
  PetscReal, allocatable :: fnorm_update_stride(:)
  PetscReal, allocatable :: fnorm_residual_stride(:)
  PetscReal, allocatable :: inorm_solution_stride(:)
  PetscReal, allocatable :: inorm_update_stride(:)
  PetscReal, allocatable :: inorm_residual_stride(:)
  
  PetscReal :: norm1_solution
  PetscReal :: norm1_update
  PetscReal :: norm1_residual
  PetscReal, allocatable :: norm1_solution_stride(:)
  PetscReal, allocatable :: norm1_update_stride(:)
  PetscReal, allocatable :: norm1_residual_stride(:)

  KSP :: ksp
  
  PetscInt, allocatable :: imax_solution(:)
  PetscInt, allocatable :: imax_update(:)
  PetscInt, allocatable :: imax_residual(:)
  PetscReal, allocatable :: max_solution_val(:)
  PetscReal, allocatable :: max_update_val(:)
  PetscReal, allocatable :: max_residual_val(:)
  
  PetscInt, allocatable :: imin_solution(:)
  PetscInt, allocatable :: imin_update(:)
  PetscInt, allocatable :: imin_residual(:)
  PetscReal, allocatable :: min_solution_val(:)
  PetscReal, allocatable :: min_update_val(:)
  PetscReal, allocatable :: min_residual_val(:)
  PetscReal, pointer :: vec_ptr(:)
  
  character(len=MAXSTRINGLENGTH) :: string, string2, string3, sec_string
  PetscBool :: print_sol_norm_info = PETSC_FALSE
  PetscBool :: print_upd_norm_info = PETSC_FALSE
  PetscBool :: print_res_norm_info = PETSC_FALSE
  PetscBool :: print_norm_by_dof_info = PETSC_FALSE
  PetscBool :: print_max_val_and_loc_info = PETSC_FALSE
  PetscBool :: print_1_norm_info = PETSC_FALSE
  PetscBool :: print_2_norm_info = PETSC_FALSE
  PetscBool :: print_inf_norm_info = PETSC_FALSE

  PetscInt :: sec_reason

!typedef enum {/* converged */
!              SNES_CONVERGED_FNORM_ABS         =  2, /* ||F|| < atol */
!              SNES_CONVERGED_FNORM_RELATIVE    =  3, /* ||F|| < rtol*||F_initial|| */
!              SNES_CONVERGED_SNORM_RELATIVE    =  4, /* Newton computed step size small; || delta x || < stol || x ||*/
!              SNES_CONVERGED_ITS               =  5, /* maximum iterations reached */
!              SNES_CONVERGED_TR_DELTA          =  7,
!              /* diverged */
!              SNES_DIVERGED_FUNCTION_DOMAIN     = -1, /* the new x location passed the function is not in the domain of F */
!              SNES_DIVERGED_FUNCTION_COUNT      = -2,
!              SNES_DIVERGED_LINEAR_SOLVE        = -3, /* the linear solve failed */
!              SNES_DIVERGED_FNORM_NAN           = -4,
!              SNES_DIVERGED_MAX_IT              = -5,
!              SNES_DIVERGED_LINE_SEARCH         = -6, /* the line search failed */
!              SNES_DIVERGED_INNER               = -7, /* inner solve failed */
!              SNES_DIVERGED_LOCAL_MIN           = -8, /* || J^T b || is small, implies converged to local minimum of F() */
!              SNES_CONVERGED_ITERATING          =  0} SNESConvergedReason;

  solver => context%solver
  option => context%option
  grid => context%grid

  if (option%use_touch_options) then
    string = 'detailed_convergence'
    if (OptionCheckTouch(option,string)) then
      if (solver%print_detailed_convergence) then
        solver%print_detailed_convergence = PETSC_FALSE
      else
        solver%print_detailed_convergence = PETSC_TRUE
      endif
    endif
  endif

  !geh: We must check the convergence here as i_iteration initializes
  !     snes->ttol for subsequent iterations.
  call SNESConvergedDefault(snes_,i_iteration,xnorm,unorm,fnorm,reason, &
                            PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
#if 0
  if (i_iteration == 0 .and. &
      option%print_screen_flag .and. solver%print_convergence) then
    write(*,'(i3," 2r:",es9.2)') i_iteration, fnorm
  endif
#endif

  ! for some reason (e.g. negative saturation/mole fraction in multiphase),
  ! we are forcing extra newton iterations
  if (option%force_newton_iteration) then
    reason = 0
!   reason = -1
    return
  endif
  
! Checking if norm exceeds divergence tolerance
!geh: inorm_residual is being used without being calculated.
!      if (fnorm > solver%max_norm .or. unorm > solver%max_norm .or. &
!        inorm_residual > solver%max_norm) then
  
  if (option%out_of_table) then
    reason = -9
  endif
   
  if (option%converged) then
    reason = 12
    ! set back to false
    option%converged = PETSC_FALSE
  endif
    
!  if (reason <= 0 .and. solver%check_infinity_norm) then
  if (solver%check_infinity_norm) then
  
    call SNESGetFunction(snes_,residual_vec,PETSC_NULL_OBJECT, &
                         PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)

    call VecNorm(residual_vec,NORM_INFINITY,inorm_residual,ierr);CHKERRQ(ierr)

    if (i_iteration > 0) then
      call SNESGetSolutionUpdate(snes_,update_vec,ierr);CHKERRQ(ierr)
      call VecNorm(update_vec,NORM_INFINITY,inorm_update,ierr);CHKERRQ(ierr)
    else
      inorm_update = 0.d0
    endif

    if (inorm_residual < solver%newton_inf_res_tol) then
      reason = 10
    else
!      if (reason > 0 .and. inorm_residual > 100.d0*solver%newton_inf_res_tol) &
!        reason = 0
    endif

    if (inorm_update < solver%newton_inf_upd_tol .and. i_iteration > 0) then
      reason = 11
    endif
    
    if (inorm_residual > solver%max_norm) then
      reason = -10
    endif
 
    ! This is to check if the secondary continuum residual convergences
    ! for nonlinear problems specifically transport
    if (solver%itype == TRANSPORT_CLASS .and. option%use_mc .and. &
       reason > 0 .and. i_iteration > 0) then
      if (option%infnorm_res_sec < solver%newton_inf_res_tol_sec) then
        sec_reason = 1
      else
        reason = 0
      endif
    endif
    
    ! force the minimum number of iterations
    if (i_iteration < solver%newton_min_iterations) then
      reason = 0
    endif

    if (option%print_screen_flag .and. solver%print_convergence) then
      i = int(reason)
      select case(i)
        case(-10)
          string = 'max_norm'
        case(-9)
          string = 'out_of_EOS_table'
        case(2)
          string = 'atol'
        case(3)
          string = 'rtol'
        case(4)
          string = 'stol'
        case(10)
          string = 'itol_res'
        case(11)
          string = 'itol_upd'
        case(12)
          string = 'itol_post_check'
        case default
          write(string,'(i3)') reason
      end select
      if (option%use_mc .and. option%ntrandof > 0 .and. solver%itype == &
          TRANSPORT_CLASS) then
        i = int(sec_reason) 
        select case(i)
          case(1)
            sec_string = 'itol_res_sec'
          case default
            write(sec_string,'(i3)') sec_reason
        end select
        write(*,'(i3," 2r:",es9.2, &
                & " 2x:",es9.2, &
                & " 2u:",es9.2, &
                & " ir:",es9.2, &
                & " iu:",es9.2, &
                & " irsec:",es9.2, &
                & " rsn: ",a, ", ",a)') &
                i_iteration, fnorm, xnorm, unorm, inorm_residual, &
                inorm_update, option%infnorm_res_sec, &
                trim(string), trim(sec_string)
      else
        write(*,'(i3," 2r:",es9.2, &
                & " 2x:",es9.2, &
                & " 2u:",es9.2, &
                & " ir:",es9.2, &
                & " iu:",es9.2, &
                & " rsn: ",a)') &
                i_iteration, fnorm, xnorm, unorm, inorm_residual, &
                inorm_update, trim(string)        
      endif
    endif
  else
  
    ! This is to check if the secondary continuum residual convergences
    ! for nonlinear problems specifically transport
    if (solver%itype == TRANSPORT_CLASS .and. option%use_mc .and. &
       reason > 0 .and. i_iteration > 0) then
      if (option%infnorm_res_sec < solver%newton_inf_res_tol_sec) then
        reason = 13
      else
        reason = 0
      endif
    endif
    
    ! force the minimum number of iterations
    if (i_iteration < solver%newton_min_iterations) then
      reason = 0
    endif

    if (option%print_screen_flag .and. solver%print_convergence) then
      i = int(reason)
      select case(i)
        case(-9)
          string = 'out_of_EOS_table'
        case(2)
          string = 'atol'
        case(3)
          string = 'rtol'
        case(4)
          string = 'stol'
        case(10)
          string = 'itol_res'
        case(11)
          string = 'itol_upd'
        case(12)
          string = 'itol_post_check'
        case(13)
          string = 'itol_res_sec'
        case default
          write(string,'(i3)') reason
      end select
      write(*,'(i3," 2r:",es10.2, &
              & " 2u:",es10.2, &
              & 32x, &
              & " rsn: ",a)') i_iteration, fnorm, unorm, trim(string)
      if (solver%print_linear_iterations) then
        call KSPGetIterationNumber(solver%ksp,i,ierr);CHKERRQ(ierr)
        write(option%io_buffer,'("   Linear Solver Iterations: ",i6)') i
        call printMsg(option)
      endif
    endif
  endif    

  if (solver%print_detailed_convergence) then

    call SNESGetSolution(snes_,solution_vec,ierr);CHKERRQ(ierr)
    ! the ctx object should really be PETSC_NULL_OBJECT.  A bug in petsc
    call SNESGetFunction(snes_,residual_vec,PETSC_NULL_OBJECT, &
                         PETSC_NULL_INTEGER, &
                         ierr);CHKERRQ(ierr)
    call SNESGetSolutionUpdate(snes_,update_vec,ierr);CHKERRQ(ierr)
    
    ! infinity norms
    call VecNorm(solution_vec,NORM_INFINITY,inorm_solution,ierr);CHKERRQ(ierr)
    call VecNorm(update_vec,NORM_INFINITY,inorm_update,ierr);CHKERRQ(ierr)
    call VecNorm(residual_vec,NORM_INFINITY,inorm_residual,ierr);CHKERRQ(ierr)

    call VecNorm(solution_vec,NORM_1,norm1_solution,ierr);CHKERRQ(ierr)
    call VecNorm(update_vec,NORM_1,norm1_update,ierr);CHKERRQ(ierr)
    call VecNorm(residual_vec,NORM_1,norm1_residual,ierr);CHKERRQ(ierr)
    
    call VecGetBlockSize(solution_vec,ndof,ierr);CHKERRQ(ierr)
    
    allocate(fnorm_solution_stride(ndof))
    allocate(fnorm_update_stride(ndof))
    allocate(fnorm_residual_stride(ndof))
    allocate(inorm_solution_stride(ndof))
    allocate(inorm_update_stride(ndof))
    allocate(inorm_residual_stride(ndof))
    allocate(norm1_solution_stride(ndof))
    allocate(norm1_update_stride(ndof))
    allocate(norm1_residual_stride(ndof))
    
    allocate(imax_solution(ndof))
    allocate(imax_update(ndof))
    allocate(imax_residual(ndof))
    allocate(max_solution_val(ndof))
    allocate(max_update_val(ndof))
    allocate(max_residual_val(ndof))

    allocate(imin_solution(ndof))
    allocate(imin_update(ndof))
    allocate(imin_residual(ndof))
    allocate(min_solution_val(ndof))
    allocate(min_update_val(ndof))
    allocate(min_residual_val(ndof))

    call VecStrideNormAll(solution_vec,NORM_1,norm1_solution_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(update_vec,NORM_1,norm1_update_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(residual_vec,NORM_1,norm1_residual_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(solution_vec,NORM_2,fnorm_solution_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(update_vec,NORM_2,fnorm_update_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(residual_vec,NORM_2,fnorm_residual_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(solution_vec,NORM_INFINITY,inorm_solution_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(update_vec,NORM_INFINITY,inorm_update_stride, &
                          ierr);CHKERRQ(ierr)
    call VecStrideNormAll(residual_vec,NORM_INFINITY,inorm_residual_stride, &
                          ierr);CHKERRQ(ierr)
    
    ! can't use VecStrideMaxAll since the index location is not currently supported.
    do i=1,ndof
      call VecStrideMax(solution_vec,i-1,imax_solution(i),max_solution_val(i), &
                        ierr);CHKERRQ(ierr)
      call VecStrideMax(update_vec,i-1,imax_update(i),max_update_val(i), &
                        ierr);CHKERRQ(ierr)
      call VecStrideMax(residual_vec,i-1,imax_residual(i),max_residual_val(i), &
                        ierr);CHKERRQ(ierr)
      ! tweak the index to get the cell id from the mdof vector
      imax_solution(i) = GridIndexToCellID(solution_vec,imax_solution(i),grid,GLOBAL)
      imax_update(i) = GridIndexToCellID(update_vec,imax_update(i),grid,GLOBAL)
      imax_residual(i) = GridIndexToCellID(residual_vec,imax_residual(i),grid,GLOBAL)
!      imax_solution(i) = imax_solution(i)/ndof
!      imax_update(i) = imax_update(i)/ndof
!      imax_residual(i) = imax_residual(i)/ndof
    enddo

    do i=1,ndof
      call VecStrideMin(solution_vec,i-1,imin_solution(i),min_solution_val(i), &
                        ierr);CHKERRQ(ierr)
      call VecStrideMin(update_vec,i-1,imin_update(i),min_update_val(i), &
                        ierr);CHKERRQ(ierr)
      call VecStrideMin(residual_vec,i-1,imin_residual(i),min_residual_val(i), &
                        ierr);CHKERRQ(ierr)
      ! tweak the index to get the cell id from the mdof vector
      imin_solution(i) = GridIndexToCellID(solution_vec,imin_solution(i),grid,GLOBAL)
      imin_update(i) = GridIndexToCellID(update_vec,imax_update(i),grid,GLOBAL)
      imin_residual(i) = GridIndexToCellID(residual_vec,imin_residual(i),grid,GLOBAL)
!      imin_solution(i) = imin_solution(i)/ndof
!      imin_update(i) = imin_update(i)/ndof
!      imin_residual(i) = imin_residual(i)/ndof
    enddo

    if (option%print_screen_flag) then
      select case(reason)
        case (10)
          string = "CONVERGED_USER_NORM_INF_REL"
        case (11)
          string = "CONVERGED_USER_NORM_INF_UPD"
        case(SNES_CONVERGED_FNORM_ABS)
          string = "SNES_CONVERGED_FNORM_ABS"
        case(SNES_CONVERGED_FNORM_RELATIVE)
          string = "SNES_CONVERGED_FNORM_RELATIVE"
        case(SNES_CONVERGED_SNORM_RELATIVE)
          string = "SNES_CONVERGED_SNORM_RELATIVE"
        case(SNES_CONVERGED_ITS)
          string = "SNES_CONVERGED_ITS"
        case(SNES_CONVERGED_TR_DELTA)
          string = "SNES_CONVERGED_TR_DELTA"
  !      case(SNES_DIVERGED_FUNCTION_DOMAIN)
  !        string = "SNES_DIVERGED_FUNCTION_DOMAIN"
        case(SNES_DIVERGED_FUNCTION_COUNT)
          string = "SNES_DIVERGED_FUNCTION_COUNT"
        case(SNES_DIVERGED_LINEAR_SOLVE)
          string = "SNES_DIVERGED_LINEAR_SOLVE"
        case(SNES_DIVERGED_FNORM_NAN)
          string = "SNES_DIVERGED_FNORM_NAN"
        case(SNES_DIVERGED_MAX_IT)
          string = "SNES_DIVERGED_MAX_IT"
        case(SNES_DIVERGED_LINE_SEARCH)
          string = "SNES_DIVERGED_LINE_SEARCH"
        case(SNES_DIVERGED_LOCAL_MIN)
          string = "SNES_DIVERGED_LOCAL_MIN"
        case(SNES_CONVERGED_ITERATING)
          string = "SNES_CONVERGED_ITERATING"
        case default
          string = "UNKNOWN"
      end select

      ! uncomment the lines below to determine data printed
      
      print_sol_norm_info = PETSC_TRUE  ! solution_vec norm information
      print_upd_norm_info = PETSC_TRUE  ! update_vec norm information
      print_res_norm_info = PETSC_TRUE  ! residual_vec norm information
    
      !print_norm_by_dof_info = PETSC_TRUE
      print_max_val_and_loc_info = PETSC_TRUE

      !print_1_norm_info = PETSC_TRUE
      print_2_norm_info = PETSC_TRUE
      print_inf_norm_info = PETSC_TRUE

      print *
      print *, 'reason: ', reason, ' - ', trim(string)
      print *, 'SNES iteration :', i_iteration
      call SNESGetKSP(snes_,ksp,ierr);CHKERRQ(ierr)
      call KSPGetIterationNumber(ksp,i,ierr);CHKERRQ(ierr)
      print *, 'KSP iterations :', i
      if (print_1_norm_info) then
        if (print_sol_norm_info) print *, 'norm_1_solution:   ', norm1_solution
        if (print_upd_norm_info) print *, 'norm_1_update:     ', norm1_update
        if (print_res_norm_info) print *, 'norm_1_residual:   ', norm1_residual
      endif
      if (print_2_norm_info) then
        if (print_sol_norm_info) print *, 'norm_2_solution:   ', fnorm_solution_stride
        if (print_upd_norm_info) print *, 'norm_2_update:     ', fnorm_update_stride
        if (print_res_norm_info) print *, 'norm_2_residual:   ', fnorm_residual_stride
      endif
      if (print_inf_norm_info) then
        if (print_sol_norm_info) print *, 'norm_inf_solution: ', inorm_solution
        if (print_upd_norm_info) print *, 'norm_inf_update:   ', inorm_update
        if (print_res_norm_info) print *, 'norm_inf_residual: ', inorm_residual
      endif
      if (print_max_val_and_loc_info) then
        print *, 'max/min locations (zero-based index) by dof:'
        do i=1,ndof
          print *, '  dof: ', i
          if (print_sol_norm_info) then
            print *, '    solution_vec max: ', imax_solution(i), max_solution_val(i)
            print *, '    solution_vec min: ', imin_solution(i), min_solution_val(i)
          endif
          if (print_upd_norm_info) then ! since update is -dx, need to invert
            print *, '    update_vec max:   ', imin_update(i), -min_update_val(i)
            print *, '    update_vec min:   ', imax_update(i), -max_update_val(i)
          endif
          if (print_res_norm_info) then
            print *, '    residual_vec max: ', imax_residual(i), max_residual_val(i)
            print *, '    residual_vec min: ', imin_residual(i), min_residual_val(i)
          endif
        enddo
      endif
      if (print_norm_by_dof_info) then
        print *, 'norm by dof:'
        do i=1,ndof
          print *, '  dof: ', i
          if (print_sol_norm_info) then
            if (print_1_norm_info) &
              print *, '    norm_1_solution:   ', norm1_solution_stride(i)
            if (print_2_norm_info) &
              print *, '    norm_2_solution:   ', fnorm_solution_stride(i)
            if (print_inf_norm_info) &
              print *, '    norm_inf_solution: ', inorm_solution_stride(i)
            if (print_1_norm_info .or. print_2_norm_info .or. &
                print_inf_norm_info) print *, '    -'
          endif
          if (print_upd_norm_info) then
            if (print_1_norm_info) &
              print *, '    norm_1_update:   ', norm1_update_stride(i)
            if (print_2_norm_info) &
              print *, '    norm_2_update:   ', fnorm_update_stride(i)
            if (print_inf_norm_info) &
              print *, '    norm_inf_update: ', inorm_update_stride(i)
            if (print_1_norm_info .or. print_2_norm_info .or. &
                print_inf_norm_info) print *, '    -'
          endif
          if (print_res_norm_info) then
            if (print_1_norm_info) &
              print *, '    norm_1_residual:   ', norm1_residual_stride(i)
            if (print_2_norm_info) &
              print *, '    norm_2_residual:   ', fnorm_residual_stride(i)
            if (print_inf_norm_info) &
              print *, '    norm_inf_residual: ', inorm_residual_stride(i)
          endif
        enddo
      endif
      print *
    endif
    
    deallocate(fnorm_solution_stride)
    deallocate(fnorm_update_stride)
    deallocate(fnorm_residual_stride)
    deallocate(inorm_solution_stride)
    deallocate(inorm_update_stride)
    deallocate(inorm_residual_stride)
    deallocate(norm1_solution_stride)
    deallocate(norm1_update_stride)
    deallocate(norm1_residual_stride)
    
    deallocate(imax_solution)
    deallocate(imax_update)
    deallocate(imax_residual)
    deallocate(max_solution_val)
    deallocate(max_update_val)
    deallocate(max_residual_val)

    deallocate(imin_solution)
    deallocate(imin_update)
    deallocate(imin_residual)
    deallocate(min_solution_val)
    deallocate(min_update_val)
    deallocate(min_residual_val)
    
  endif
  
end subroutine ConvergenceTest

! ************************************************************************** !

subroutine ConvergenceContextDestroy(context)
  ! 
  ! Destroy context
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 

  implicit none
  
  type(convergence_context_type), pointer :: context
  
  if (.not.associated(context)) return
  
  nullify(context%solver)
  nullify(context%option)
  nullify(context%grid)
  
  deallocate(context)
  nullify(context)

end subroutine ConvergenceContextDestroy

end module Convergence_module
