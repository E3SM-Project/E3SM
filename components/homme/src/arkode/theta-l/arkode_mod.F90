#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module arkode_mod

  use element_state,                 only: timelevels
  use derivative_mod,                only: derivative_t
  use HommeNVector,                  only: NVec_t
  use hybrid_mod,                    only: hybrid_t
  use hybvcoord_mod,                 only: hvcoord_t
  use kinds,                         only: real_kind
  use fsundials_matrix_mod,          only: SUNMatrix
  use fsundials_linearsolver_mod,    only: SUNLinearSolver
  use fsundials_nonlinearsolver_mod, only: SUNNonlinearSolver
  use fsundials_nvector_mod,         only: N_Vector
  use iso_c_binding

  implicit none

  private

  integer, parameter :: freelevels = timelevels-35
  ! Note that 35 is an estimate, but needs to be at least > 30
  ! If a larger Krylov subspace is desired, timelevels should be
  ! increased.
  integer, parameter :: max_niters = 10

  type(c_ptr) :: arkode_mem ! ARKode memory object

  type(SUNMatrix), pointer          :: sunmat ! sundials matrix object (dummy)
  type(SUNLinearSolver), pointer    :: sunls  ! sundials linear solver object
  type(SUNNonlinearSolver), pointer :: sunnls ! sundials nonlinear solver object

  ! Nonlinear solver convergence test variables
  real(c_double) :: crate ! convergence rate estimate
  real(c_double) :: delp  ! previous correction norm

  ! ARKode namelist variables
  real(real_kind), public :: rel_tol = 1.d-8
  real(real_kind), public :: abs_tol = 1.d-8  ! val < 0 indicates array atol
  logical, public         :: calc_nonlinear_stats = .true.
  logical, public         :: use_column_solver = .true.

  ! data type for passing ARKode parameters
  type :: parameter_list
    ! Linear Solver Info (flag for columnwise/GMRES, GMRES parameters)
    integer         :: precLR ! preconditioning: 0=none, 1=left, 2=right, 3=left+right
    integer         :: gstype ! Gram-Schmidt orthogonalization: 1=modified, 2=classical
    integer         :: maxl = freelevels ! max size of Krylov subspace (# of iterations/vectors)
    real(real_kind) :: lintol ! linear convergence tolerance factor (0 indicates default)
    ! General Iteration Info
    real(real_kind) :: rtol ! relative tolerance for iteration convergence
    real(real_kind) :: atol(6) ! absolute tolerances (u,v,w,phinh,vtheta_dp,dp3d)
  end type parameter_list

  ! dummy data type that enables an array of N_Vector pointers
  type N_VectorPtr
    type(N_Vector), pointer :: ptr
  end type N_VectorPtr

  public :: parameter_list, evolve_solution, get_hvcoord_ptr
  public :: update_nonlinear_stats, finalize_nonlinear_stats

  save

  type(hvcoord_t), pointer      :: hvcoord_ptr
  type(hybrid_t), pointer       :: hybrid_ptr
  type(derivative_t), pointer   :: deriv_ptr
  type(parameter_list), pointer :: param_ptr
  type(N_VectorPtr)             :: atol
  type(N_VectorPtr)             :: y(3)
  type(c_funptr)                :: efun_ptr, ifun_ptr
  real(real_kind)               :: dt_save, eta_ave_w_save, rout(40)
  integer                       :: imex_save
  logical                       :: initialized = .false.
  ! variables used for nonlinear solver stats, not necessary for use of ARKode
  integer :: total_nonlinear_iterations = 0
  integer :: max_nonlinear_iterations = 0
  integer :: num_timesteps = 0

contains

  subroutine update_nonlinear_stats(timesteps, nonlinear_iters)
    !-----------------------------------------------------------------
    ! Description: update ARKode performance statistics
    !   Arguments:
    !          timesteps - (int, input, optional) specify # of timesteps
    !    nonlinear_iters - (int, input, optional) specify # of nonlinear iters
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use farkode_arkstep_mod, only: FARKStepGetNumNonlinSolvIters
    use parallel_mod,        only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    integer, intent(in), optional :: timesteps
    integer, intent(in), optional :: nonlinear_iters

    ! local variables
    integer(c_int) :: ierr
    integer(c_long) :: num_iters(1)

    !======= Internals ============


    if (present(timesteps)) then
      num_timesteps = num_timesteps + timesteps
    else
      num_timesteps = num_timesteps + 1
    end if
    if (present(nonlinear_iters)) then
      max_nonlinear_iterations = max(max_nonlinear_iterations, nonlinear_iters)
      total_nonlinear_iterations = total_nonlinear_iterations + nonlinear_iters
    else
      ierr = FARKStepGetNumNonlinSolvIters(arkode_mem, num_iters)
      if (ierr /= 0) then
        call abortmp('FARKStepGetNumNonlinSolvIters failed')
      end if
      max_nonlinear_iterations = max(max_nonlinear_iterations, num_iters(1))
      total_nonlinear_iterations = total_nonlinear_iterations + num_iters(1)
    end if
    return
  end subroutine update_nonlinear_stats

  !=================================================================

  subroutine finalize_nonlinear_stats(comm, my_rank, master_rank, comm_size)
    !-----------------------------------------------------------------
    ! Description: print ARKode performance statistics
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use mpi

    !======= Declarations =========
    implicit none

    ! calling variables
    integer, intent(in) :: comm
    integer, intent(in) :: my_rank
    integer, intent(in) :: master_rank
    integer, intent(in) :: comm_size

    ! local variables
    integer :: max_result, sum_result, ierr

    !======= Internals ============
    call MPI_Reduce(max_nonlinear_iterations, max_result, 1, MPI_INTEGER, &
                    MPI_MAX, master_rank, comm, ierr)
    call MPI_Reduce(total_nonlinear_iterations, sum_result, 1, MPI_INTEGER, &
                    MPI_SUM, master_rank, comm, ierr)
    if (my_rank == master_rank) then
      print *, 'ARKode Nonlinear Solver Statistics:'
      print '(2x,A,i9)','Max num nonlin iters   =', max_result
      print '(2x,A,i9)','Total num nonlin iters =', sum_result
      print '(2x,A,i9)','Total num timesteps    =', num_timesteps
      print '(2x,A,f9.2)','Avg num nonlin iters   =', sum_result/float(comm_size*num_timesteps)
    end if

    return
  end subroutine finalize_nonlinear_stats

  !=================================================================

  subroutine get_hvcoord_ptr(hvcoord)
    !-----------------------------------------------------------------
    ! Description: sets pointer to current hvcoord object
    !   Arguments:
    !     hvcoord - (obj*, output) hvcoord object pointer
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use hybvcoord_mod,  only: hvcoord_t

    !======= Declarations =========
    implicit none

    ! calling variables
    type(hvcoord_t),    intent(out) :: hvcoord

    !======= Internals ============
    hvcoord = hvcoord_ptr

    return
  end subroutine get_hvcoord_ptr

  !=================================================================


  function evolve_solution(elem, nets, nete, deriv, hvcoord, hybrid, &
                           dt, eta_ave_w, n0, np1, arkode_parameters, &
                           table_set) result(ierr)
    !-----------------------------------------------------------------
    ! Description: resets internal ARKode solution without memory allocation
    !   Arguments:
    !                 elem - (obj*, input) element objects
    !                 nets - (int, input) starting index for elem array
    !                 nete - (int, input) ending index for elem array
    !                   tl - (obj*, input) timelevel object
    !                deriv - (obj*, input) deriv object
    !              hvcoord - (obj*, input) hvcoord object
    !               hybrid - (obj*, input) hybrid object
    !                   dt - (real, input) current timestep
    !            eta_ave_w - (real, input) average flux value
    !                   n0 - (int, input) timelevel holding current solution
    !                  np1 - (int, input) timelevel for evolved solution
    !    arkode_parameters - (parameter, input) object for arkode parameters
    !            table_set - (butcher_table_set, input) Butcher table(s)
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use arkode_tables,       only: butcher_table_set
    use derivative_mod,      only: derivative_t
    use element_mod,         only: element_t
    use farkode_mod,         only: ARK_ONE_STEP
    use farkode_arkstep_mod, only: FARKStepEvolve, FARKStepSetFixedStep
    use HommeNVector,        only: MakeHommeNVector, SetHommeNVectorPar
    use hybrid_mod,          only: hybrid_t
    use hybvcoord_mod,       only: hvcoord_t
    use kinds,               only: real_kind
    use parallel_mod,        only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(derivative_t), target, intent(in) :: deriv
    type(hvcoord_t), target,    intent(in) :: hvcoord
    type(hybrid_t), target,     intent(in) :: hybrid
    type(parameter_list),       intent(in) :: arkode_parameters
    type(butcher_table_set),    intent(in) :: table_set
    real(real_kind),            intent(in) :: dt, eta_ave_w
    integer,                    intent(in) :: nets, nete, n0, np1
    type(element_t),            intent(inout) :: elem(:)

    ! local variables
    real(real_kind)              :: tstart
    real(real_kind)              :: tstop(1)
    integer(C_INT)               :: iflag, ierr
    integer                      :: i
    !======= Internals ============
    ! Set NVector communicator
    call SetHommeNVectorPar(hybrid%par)

    ! specify start time to be 0.0 so stage time available in efun & ifun
    tstart = 0.d0

    ! store variables for efun and ifun
    dt_save = dt
    eta_ave_w_save = eta_ave_w
    imex_save = table_set%imex
    hybrid_ptr => hybrid
    deriv_ptr => deriv
    hvcoord_ptr => hvcoord

    ! Initialize or reinitialize ARKode
    if (.not.initialized) then
      call initialize(elem, nets, nete, hybrid%par, n0, tstart, &
                      arkode_parameters, table_set)
      initialized = .true.
    else
      call reinitialize(n0, tstart, arkode_parameters)
    end if

    ! Set ARKode time step
    ierr = FARKStepSetFixedStep(arkode_mem, dt)
    if (ierr /= 0) then
      call abortmp('FARKStepSetFixedStep failed')
    endif

    ! Take single step
    tstop(1) = tstart + dt
    ierr = FARKStepEvolve(arkode_mem, tstart + dt, y(np1)%ptr, tstop, ARK_ONE_STEP)

  end function evolve_solution

  !=================================================================

  subroutine reinitialize(n0, tstart, arkode_parameters)
    !-----------------------------------------------------------------
    ! Description: resets internal ARKode solution without memory allocation
    !   Arguments:
    !                   n0 - (int, input) timelevel holding current solution
    !               tstart - (real, input) time to start ARKode solve at
    !    arkode_parameters - (parameter, input) arkode parameter object
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use farkode_arkstep_mod, only: FARKStepSetMaxNonlinIters, FARKStepReInit, &
                                   FARKStepSVtolerances
    use HommeNVector,        only: FN_VGetContent
    use kinds,               only: real_kind
    use parallel_mod,        only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parameter_list), target, intent(in)  :: arkode_parameters
    integer,                      intent(in)  :: n0
    real(real_kind),              intent(in)  :: tstart

    ! local variables
    type(NVec_t), pointer         :: atol_content => NULL()
    type(parameter_list), pointer :: ap
    integer(C_INT)                :: ierr
    integer                       :: ie, tl

    !======= Internals ============
    ap => arkode_parameters

    atol_content => FN_VGetContent(atol%ptr)

    ! update rtol and atol
    tl = atol_content%tl_idx
    do ie=atol_content%nets,atol_content%nete
      atol_content%elem(ie)%state%v(:,:,1,:,tl) = ap%atol(1)
      atol_content%elem(ie)%state%v(:,:,2,:,tl) = ap%atol(2)
      atol_content%elem(ie)%state%w_i(:,:,:,tl) = ap%atol(3)
      atol_content%elem(ie)%state%phinh_i(:,:,:,tl) = ap%atol(4)
      atol_content%elem(ie)%state%vtheta_dp(:,:,:,tl) = ap%atol(5)
      atol_content%elem(ie)%state%dp3d(:,:,:,tl) = ap%atol(6)
    end do
    ierr = FARKStepSVtolerances(arkode_mem, ap%rtol, atol%ptr)
    if (ierr /= 0) then
      call abortmp('FARKStepSVtolerances failed')
    endif

    ! reinitialize ARKode with current solution
    ierr = FARKStepReInit(arkode_mem, efun_ptr, ifun_ptr, tstart, y(n0)%ptr)
    if (ierr /= 0) then
      call abortmp('FARKStepReInit failed')
    endif

    if (c_associated(ifun_ptr)) then
      ! reset max number of nonlinear iterations per stage
      ierr = FARKStepSetMaxNonlinIters(arkode_mem, max_niters)
      if (ierr /= 0) then
        call abortmp('FARKStepSetMaxNonlinIters failed')
      endif
    end if

    return
  end subroutine reinitialize

  !=================================================================

  subroutine initialize(elem, nets, nete, par, n0, tstart, &
                        arkode_parameters, table_set)
    !-----------------------------------------------------------------
    ! Description: allocates memory for and initializes ARKode
    !   Arguments:
    !                 elem - (obj*, input) element objects
    !                 nets - (int, input) starting index for elem array
    !                 nete - (int, input) ending index for elem array
    !                  par - (obj*, input) parallel object
    !                   n0 - (int, input) timelevel holding current solution
    !               tstart - (real, input) time to start ARKode solve at
    !    arkode_parameters - (parameter, input) object for arkode parameters
    !            table_set - (butcher_table_set, input) object for Butcher table(s)
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use arkode_tables,            only: butcher_table_set
    use element_mod,              only: element_t
    use farkode_mod,              only: FARKodeButcherTable_Create
    use farkode_arkstep_mod,      only: FARKStepCreate, FARKStepSetTables, &
                                        FARKStepSetMaxNonlinIters, FARKStepSVtolerances, &
                                        FARKStepSetDiagnostics, FARKStepWriteParameters, &
                                        FARKStepSetLinearSolver, &
                                        FARKStepSetLinearSolutionScaling, &
                                        FARKStepSetLinSysFn, FARKStepSetEpsLin, &
                                        FARKStepSetNonlinearSolver
    use HommeNVector,             only: NVec_t, MakeHommeNVector, SetHommeNVectorPar
    use HommeSUNLinSol,           only: FSUNMatrix_HOMME, FSUNLinSol_HOMME, FARKodeLinSysFn
    use fsunnonlinsol_newton_mod, only: FSUNNonlinSol_Newton, FSUNNonlinSolSetConvTestFn_Newton
    use fsunlinsol_spgmr_mod,     only: FSUNLinSol_SPGMR, FSUNLinSol_SPGMRSetGSType
    use kinds,                    only: real_kind
    use parallel_mod,             only: parallel_t, abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parallel_t),                   intent(in) :: par
    type(parameter_list), target,       intent(in) :: arkode_parameters
    type(butcher_table_set), target,    intent(in) :: table_set
    real(real_kind),                    intent(in) :: tstart
    integer,                            intent(in) :: nets, nete, n0
    type(element_t),                    intent(inout) :: elem(:)

    ! local variables
    type(parameter_list), pointer          :: ap
    type(butcher_table_set), pointer       :: ts
    real(real_kind)                        :: rpar(1)
    type(c_ptr)                            :: Te
    type(c_ptr)                            :: Ti
    type(c_ptr)                            :: file_ptr
    real(real_kind)                        :: a(table_set%s*table_set%s)
    integer(C_INT)                         :: ierr
    integer                                :: i, j

    !======= Internals ============
    ap => arkode_parameters
    ts => table_set

    if (par%masterproc) print *,"Initializing ARKode"
    ! 'create' NVec_t objects that will correspond to the original 3 HOMME
    ! timelevels, assuming that tl%nm1, tl%n0, and tl%np1 are taken from the
    ! set {1,2,3}
    do i=1,3
      y(i)%ptr => MakeHommeNVector(elem, nets, nete, i, ierr)
      if (ierr /= 0) then
        call abortmp('Error in MakeHommeNVector')
      end if
    end do

    ! save rtol, set data in 4th timelevel to atol values,
    ! and 'create' NVec_t object
    do i=nets,nete
      elem(i)%state%v(:,:,1,:,4) = ap%atol(1)
      elem(i)%state%v(:,:,2,:,4) = ap%atol(2)
      elem(i)%state%w_i(:,:,:,4) = ap%atol(3)
      elem(i)%state%phinh_i(:,:,:,4) = ap%atol(4)
      elem(i)%state%vtheta_dp(:,:,:,4) = ap%atol(5)
      elem(i)%state%dp3d(:,:,:,4) = ap%atol(6)
    end do
    atol%ptr => MakeHommeNVector(elem, nets, nete, 4, ierr)
    if (ierr /= 0) then
      call abortmp('Error in MakeHommeNVector')
    end if

    ! set explicit and/or implicit RHS pointers, allocate ARKode memory,
    ! and set tolerances
    if (ts%imex == 0) then
      efun_ptr = c_null_funptr
      ifun_ptr = c_funloc(ifun)
    else if (ts%imex == 1) then
      efun_ptr = c_funloc(efun)
      ifun_ptr = c_null_funptr
    else if (ts%imex == 2) then
      efun_ptr = c_funloc(efun)
      ifun_ptr = c_funloc(ifun)
    else
      call abortmp('Invalid value of table_set%imex in arkode_mod.F90')
    end if
    arkode_mem = FARKStepCreate(efun_ptr, ifun_ptr, tstart, y(n0)%ptr)
    ierr = FARKStepSVtolerances(arkode_mem, ap%rtol, atol%ptr)
    if (ierr /= 0) then
      call abortmp('FARKStepSVtolerances failed')
    end if

    ! Set IRK Butcher table for implicit or ARK problems
    if (c_associated(ifun_ptr)) then

      ! flatten Ai to C array (row-major)
      do i=1,ts%s
        do j = 1,ts%s
          a(ts%s*(i-1)+j) = ts%Ai(i,j)
        end do
      end do
      ! create implicit table
      Ti = FARKodeButcherTable_Create(ts%s, ts%q, ts%p, ts%ci, a, ts%bi, ts%bi2)
    end if

    ! Set ERK Butcher table for explicit or ARK problems
    if (c_associated(efun_ptr)) then
      ! flatten Ae to C array (row-major)
      do i=1,ts%s
        do j = 1,ts%s
          a(ts%s*(i-1)+j) = ts%Ae(i,j)
        end do
      end do
      ! create explicit table
      Te = FARKodeButcherTable_Create(ts%s, ts%q, ts%p, ts%ce, a, ts%be, ts%be2)

    else
      call abortmp('arkode_init: invalid imex parameter value')
    end if

    ierr = FARKStepSetTables(arkode_mem, ts%q, ts%p, Ti, Te)
    if (ierr /= 0) then
      call abortmp('FARKStepSetTables failed')
    end if


    ! Set solver if implicit or imex problem
    if (c_associated(ifun_ptr)) then

       nullify(sunmat)
       nullify(sunls)
       nullify(sunnls)

       if (use_column_solver) then

          ! create a dummy matrix
          sunmat => FSUNMatrix_HOMME()
          if (.not.associated(sunmat)) then
             call abortmp('arkode_init: FSUNMatrix_HOMME failed')
          end if

          ! create the HOMME columnwise direct solver wrapper
          sunls => FSUNLinSol_HOMME(arkode_mem)
          if (.not.associated(sunls)) then
             call abortmp('arkode_init: FSUNLinSol_HOMME failed')
          end if

       else

          ! use the SUNDIALS GMRES linear solver (and set options)
          sunls => FSUNLinSol_SPGMR(y(n0)%ptr, ap%precLR, ap%maxl)
          if (.not.associated(sunls)) then
             call abortmp('arkode_init: FSUNLinSol_SPGMR failed')
          end if

          ! set the orthogonalization type
          ierr = FSUNLinSol_SPGMRSetGSType(sunls, ap%gstype)
          if (ierr /= 0) then
             call abortmp('arkode_init: FSUNLinSol_SPGMRSetGSType failed')
          end if

       endif

       ! attach the linear solver and matrix to ARKode
       ierr = FARKStepSetLinearSolver(arkode_mem, sunls, sunmat)
       if (ierr /= 0) then
          call abortmp('arkode_init: FARKStepSetLinearSolver failed')
       end if

       ! set linear solver related options
       if (use_column_solver) then

          ! set a dummy linear system function
          ierr = FARKStepSetLinSysFn(arkode_mem, c_funloc(FARKodeLinSysFn))
          if (ierr /= 0) then
             call abortmp('arkode_init: FARKStepSetLinearSolver failed')
          end if

          ! disable linear system solution scaling
          ierr = FARKStepSetLinearSolutionScaling(arkode_mem, 0)
          if (ierr /= 0) then
             call abortmp('arkode_init: FARKStepSetLinearSolver failed')
          end if

       else

          ! set the linear solve tolerance factor
          ierr = FARKStepSetEpsLin(arkode_mem, ap%lintol)
          if (ierr /= 0) then
             call abortmp('arkode_init: FARKSpilsSetEpsLin failed')
          end if

       end if

       ! create SUNDIALS Newton solver
       sunnls => FSUNNonlinSol_Newton(y(n0)%ptr)
       if (.not.associated(sunnls)) then
          call abortmp('arkode_init: FSUNNonlinSol_Newton failed')
       end if

       ! attach the nonlinear solver
       ierr = FARKStepSetNonlinearSolver(arkode_mem, sunnls)
       if (ierr /= 0) then
          call abortmp('arkode_init: FARKStepSetNonlinearSolver failed')
       endif

       ! set custom convergence test function
       ierr = FSUNNonlinSolSetConvTestFn_Newton(sunnls, c_funloc(convtest), c_null_ptr)
       if (ierr /= 0) then
          call abortmp('arkode_init: FSUNNonlinSolSetConvTestFn_Newton failed')
       endif

       ! set max number of nonlinear iterations per stage
       ierr = FARKStepSetMaxNonlinIters(arkode_mem, max_niters)
       if (ierr /= 0) then
          call abortmp('arkode_init: FARKStepSetMaxNonlinIters failed')
       endif

    endif

    return
  end subroutine initialize

  integer(c_int) function convtest(sunnls, sunvec_y, sunvec_del, tol, sunvec_ewt, mem) result(ierr) bind(C)
    !-----------------------------------------------------------------
    ! Description: convtest provides a custom conergence test function
    ! to the default Newton nonlinear solver
    !
    ! Arguments:
    !      sunnls - (SUNNonlinearSolver, input) nonlinear solver object
    !    sunvec_y - (N_Vector, input) NLS correction vector
    !  sunvec_del - (N_Vector, input) NLS correction update vector
    !         tol - (dbl, input) solve tolerance in WRMS norm
    !  sunvec_ewt - (N_Vector, input) weights for WRMS norm
    !         mem - (void* ptr, input) user supplied memory (unused)
    !        ierr - (cint, output) return flag: 0 = success,
    !                    0 > not yet converged, 0 < fail/divergent
    !-----------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use fsundials_nvector_mod,         only: FN_VWrmsNorm
    use fsundials_nonlinearsolver_mod, only: FSUNNonlinSolGetCurIter, &
         SUN_NLS_SUCCESS, &
         SUN_NLS_CONTINUE

    implicit none

    type(SUNNonlinearSolver) :: sunnls
    type(N_Vector)           :: sunvec_y
    type(N_Vector)           :: sunvec_del
    real(c_double), value    :: tol
    type(N_Vector)           :: sunvec_ewt
    type(c_ptr), value       :: mem

    real(c_double) :: delnrm
    real(c_double) :: dcon
    integer(c_int) :: citer(1)

    ! compute the norm of the correction
    delnrm = FN_VWrmsNorm(sunvec_del, sunvec_ewt)

    ! get the current nonlinear solver iteration
    ierr = FSUNNonlinSolGetCurIter(sunnls, citer)
    if (ierr /= 0) then
       return
    end if

    ! update the convergence rate
    if (citer(1) == 0) then
       crate = 1.0d0
    else
       crate = max(0.3d0 * crate, delnrm / delp)
    end if

    ! save norm for next iteration
    delp = delnrm

    ! compute scaled norm to test convergence
    dcon = min(crate, 1.0d0) * delnrm / tol

    ! check for convergence
    if (dcon < 1.0d0) then
       ierr = SUN_NLS_SUCCESS
    else
       ierr = SUN_NLS_CONTINUE
    end if

  end function convtest

  function efun(t, sunvec_y, sunvec_ydot, data) result(ierr) bind(C)
    !-----------------------------------------------------------------
    ! Description: efun provides the explicit portion of the right
    !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
    !
    ! Arguments:
    !           t - (dbl, input) current time
    !    sunvec_y - (N_Vector, input) N_Vector containing current solution
    ! sunvec_ydot - (N_Vector, input) N_Vector to hold right-hand side function
    !        data - (void(*), input) user parameter data (unused here)
    !        ierr - (cint, output) return flag: 0=>success,
    !                    1=>recoverable error, -1=>non-recoverable error
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use HommeNVector, only: NVec_t, FN_VGetContent
    use fsundials_nvector_mod, only: N_Vector

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t
    type(N_Vector)        :: sunvec_y
    type(N_Vector)        :: sunvec_ydot
    type(c_ptr)           :: data
    integer(c_int)        :: ierr

    ! local variables
    type(NVec_t), pointer :: y => NULL()
    type(NVec_t), pointer :: ydot => NULL()
    real(real_kind)       :: dt, eta_ave_w, bval, cval, scale1, scale2
    integer               :: imex, ie
    integer               :: farkefun

    !======= Internals ============

    y => FN_VGetContent(sunvec_y)
    ydot => FN_VGetContent(sunvec_ydot)
    ierr = farkefun(t, y, ydot, hvcoord_ptr, hybrid_ptr, deriv_ptr, &
                    imex_save, dt_save, eta_ave_w_save)

  end function efun

  function ifun(t, sunvec_y, sunvec_ydot, data) result(ierr) bind(C)
    !-----------------------------------------------------------------
    ! Description: ifun provides the implicit portion of the right
    !     hand side function for the ODE:   dy/dt = fi(t,y) + fe(t,y)
    !
    ! Arguments:
    !           t - (dbl, input) current time
    !    sunvec_y - (N_Vector, input) N_Vector containing current solution
    ! sunvec_ydot - (N_Vector, input) N_Vector to hold right-hand side function
    !        data - (void(*), input) user parameter data (unused here)
    !        ierr - (cint, output) return flag: 0=>success,
    !                    1=>recoverable error, -1=>non-recoverable error
    !-----------------------------------------------------------------
    !======= Inclusions ===========
    use HommeNVector, only: NVec_t, FN_VGetContent
    use fsundials_nvector_mod, only: N_Vector

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value  :: t
    type(N_Vector)         :: sunvec_y
    type(N_Vector)         :: sunvec_ydot
    type(c_ptr)            :: data
    integer(c_int)         :: ierr

    ! local variables
    type(NVec_t), pointer :: y => NULL()
    type(NVec_t), pointer :: ydot => NULL()
    real(real_kind)       :: dt, eta_ave_w, bval, cval, scale1, scale2
    integer               :: imex, ie
    integer               :: farkifun

    !======= Internals ============

    y => FN_VGetContent(sunvec_y)
    ydot => FN_VGetContent(sunvec_ydot)
    ierr = farkifun(t, y, ydot, hvcoord_ptr, hybrid_ptr, deriv_ptr,  &
                    imex_save, dt_save, eta_ave_w_save)

  end function ifun

end module arkode_mod
