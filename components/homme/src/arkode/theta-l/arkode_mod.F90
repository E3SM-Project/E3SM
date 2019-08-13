#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define ARK324_ARK 1
#define ARK436_ARK 2
#define ARK453_ARK 3
#define ARS222_ARK 4
#define ARS232_ARK 5
#define ARS233_ARK 6
#define ARS343_ARK 7
#define ARS443_ARK 8
#define SSP3333B_ARK 9
#define SSP3333C_ARK 10
#define RK2_ARK 11
#define KGU35_ARK 12
#define IMKG232_ARK 13
#define IMKG242_ARK 14
#define IMKG243_ARK 15
#define IMKG252_ARK 16
#define IMKG253_ARK 17
#define IMKG254_ARK 18
#define IMKG342_ARK 19
#define IMKG343_ARK 20
#define IMKG353_ARK 21
#define IMKG354_ARK 22




module arkode_mod

  use element_state,  only: timelevels
  use derivative_mod, only: derivative_t
  use HommeNVector,   only: NVec_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind
  use iso_c_binding

  implicit none

  private

  integer, parameter :: max_stage_num = 10
  integer, parameter :: freelevels = timelevels-35
  ! Note that 35 is an estimate, but needs to be at least > 30
  ! If a larger Krylov subspace is desired, timelevels should be
  ! increased.
  integer, parameter :: max_niters = 10

  ! ARKode namelist variables
  real(real_kind), public :: rel_tol = 1.d-8
  real(real_kind), public :: abs_tol = 1.d-8  ! val < 0 indicates array atol
  logical, public         :: calc_nonlinear_stats = .true.
  logical, public         :: use_column_solver = .true.

  ! data type for passing ARKode Butcher table names
  type :: table_list
    integer :: ARK324   = ARK324_ARK
    integer :: ARK436   = ARK436_ARK
    integer :: ARK453   = ARK453_ARK
    integer :: ARS222   = ARS222_ARK
    integer :: ARS232   = ARS232_ARK
    integer :: ARS233   = ARS233_ARK
    integer :: ARS343   = ARS343_ARK
    integer :: ARS443   = ARS443_ARK
    integer :: SSP3333B = SSP3333B_ARK
    integer :: SSP3333C = SSP3333C_ARK
    integer :: RK2      = RK2_ARK
    integer :: KGU35    = KGU35_ARK
    integer :: IMKG232  = IMKG232_ARK
    integer :: IMKG242  = IMKG242_ARK
    integer :: IMKG243  = IMKG243_ARK
    integer :: IMKG252  = IMKG252_ARK
    integer :: IMKG253  = IMKG253_ARK
    integer :: IMKG254  = IMKG254_ARK
    integer :: IMKG342  = IMKG342_ARK
    integer :: IMKG343  = IMKG343_ARK
    integer :: IMKG353  = IMKG353_ARK
    integer :: IMKG354  = IMKG353_ARK
  end type table_list

  ! data type for passing ARKode parameters
  type :: parameter_list
    ! *RK Method Information
    integer         :: imex ! 0=implicit, 1=explicit, 2=imex
    integer         :: s ! number of stages
    integer         :: q ! method order
    integer         :: p ! embedded method order
    ! Explicit Butcher Table
    real(real_kind) :: Ae(max_stage_num,max_stage_num)
    real(real_kind) :: be(max_stage_num)
    real(real_kind) :: be2(max_stage_num)
    real(real_kind) :: ce(max_stage_num)
    ! Implicit Butcher Table
    real(real_kind) :: Ai(max_stage_num,max_stage_num)
    real(real_kind) :: bi(max_stage_num)
    real(real_kind) :: bi2(max_stage_num)
    real(real_kind) :: ci(max_stage_num)
    ! Linear Solver Info (flag for columnwise/GMRES, GMRES parameters)
    integer         :: precLR ! preconditioning: 0=none, 1=left, 2=right, 3=left+right
    integer         :: gstype ! Gram-Schmidt orthogonalization: 1=modified, 2=classical
    integer         :: maxl = freelevels ! max size of Krylov subspace (# of iterations/vectors)
    real(real_kind) :: lintol ! linear convergence tolerance factor (0 indicates default)
    ! General Iteration Info
    real(real_kind) :: rtol ! relative tolerance for iteration convergence
    real(real_kind) :: atol(6) ! absolute tolerances (u,v,w,phinh,vtheta_dp,dp3d)
  end type parameter_list

  public :: parameter_list, update_arkode, get_solution_ptr, get_hvcoord_ptr
  public :: get_qn0, get_RHS_vars
  public :: max_stage_num, table_list, set_Butcher_tables
  public :: update_nonlinear_stats, finalize_nonlinear_stats

  save

  type(hvcoord_t), pointer      :: hvcoord_ptr
  type(hybrid_t), pointer       :: hybrid_ptr
  type(derivative_t), pointer   :: deriv_ptr
  type(parameter_list), pointer :: param_ptr
  type(NVec_t), target          :: y_F(3), atol_F
  type(c_ptr)                   :: y_C(3), atol_C
  real(real_kind)               :: dt_save, eta_ave_w_save, rout(40)
  integer                       :: imex_save, qn0_save
  logical                       :: initialized = .false.
  integer(C_LONG)               :: iout(40)
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

    !======= Declarations =========
    implicit none

    ! calling variables
    integer, intent(in), optional :: timesteps
    integer, intent(in), optional :: nonlinear_iters

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
      max_nonlinear_iterations = max(max_nonlinear_iterations, iout(11))
      total_nonlinear_iterations = total_nonlinear_iterations + iout(11)
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

  subroutine get_solution_ptr(np1, ynp1)
    !-----------------------------------------------------------------
    ! Description: sets pointer to NVector for specified timelevel
    !   Arguments:
    !     np1 - (int, input) timelevel
    !    ynp1 - (cptr, output) NVector pointer
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    integer,      intent(in) :: np1
    type(c_ptr), intent(out) :: ynp1

    !======= Internals ============
    ynp1 = y_C(np1)

    return
  end subroutine get_solution_ptr

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

  subroutine get_qn0(qn0)
    !-----------------------------------------------------------------
    ! Description: obtains current qn0 value
    !   Arguments:
    !     qn0 - (int, output) qn0 value
    !-----------------------------------------------------------------

    !======= Declarations =========
    implicit none

    ! calling variables
    integer, intent(out) :: qn0

    !======= Internals ============
    qn0 = qn0_save

    return
  end subroutine get_qn0

  !=================================================================

  subroutine get_RHS_vars(imex, dt, eta_ave_w, hybrid, deriv)
    !-----------------------------------------------------------------
    ! Description: sets variables and objects needed to compute RHS
    !   Arguments:
    !        imex - (int, output) variable for imex specification
    !          dt - (real, output) variable for timestep size
    !   eta_ave_w - (real, output) variable for eta_ave_w value
    !      hybrid - (obj*, output) hybrid object pointer
    !       deriv - (obj*, output) deriv object pointer
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use derivative_mod, only: derivative_t
    use hybrid_mod,     only: hybrid_t
    use kinds,          only: real_kind

    !======= Declarations =========
    implicit none

    ! calling variables
    type(hybrid_t),     intent(out) :: hybrid
    type(derivative_t), intent(out) :: deriv
    integer,            intent(out) :: imex
    real(real_kind),    intent(out) :: dt, eta_ave_w

    !======= Internals ============
    imex = imex_save
    dt = dt_save
    eta_ave_w = eta_ave_w_save
    hybrid = hybrid_ptr
    deriv = deriv_ptr

    return
  end subroutine get_RHS_vars

  !=================================================================

  subroutine update_arkode(elem, nets, nete, deriv, hvcoord, hybrid, &
                           dt, eta_ave_w, n0, qn0, arkode_parameters)
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
    !                  qn0 - (int, input) timelevel for tracer mass
    !    arkode_parameters - (parameter, input) object for arkode parameters
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use derivative_mod,   only: derivative_t
    use element_mod,      only: element_t
    use HommeNVector,     only: NVec_t, MakeHommeNVector, SetHommeNVectorPar
    use hybrid_mod,       only: hybrid_t
    use hybvcoord_mod,    only: hvcoord_t
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp
    use iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(derivative_t), target, intent(in) :: deriv
    type(hvcoord_t), target,    intent(in) :: hvcoord
    type(hybrid_t), target,     intent(in) :: hybrid
    type(parameter_list),       intent(in) :: arkode_parameters
    real(real_kind),            intent(in) :: dt, eta_ave_w
    integer,                    intent(in) :: nets, nete, n0, qn0
    type(element_t),            intent(inout) :: elem(:)

    ! local variables
    real(real_kind) :: tstart
    integer(C_INT)  :: iflag, ierr
    integer         :: i

    !======= Internals ============
    ! Set NVector communicator
    call SetHommeNVectorPar(hybrid%par)

    ! specify start time to be 0.0 so stage time available in farkefun & farkifun
    tstart = 0.d0

    ! store variables for farkefun and farkifun
    dt_save = dt
    eta_ave_w_save = eta_ave_w
    imex_save = arkode_parameters%imex
    qn0_save = qn0
    hybrid_ptr => hybrid
    deriv_ptr => deriv
    hvcoord_ptr => hvcoord

    ! Initialize or reinitialize ARKode
    if (.not.initialized) then
      call initialize(elem, nets, nete, hybrid%par, n0, qn0, tstart, &
                      arkode_parameters)
      initialized = .true.
    else
      call reinitialize(n0, tstart, arkode_parameters)
    end if

    ! Set ARKode time step
    call farksetrin('FIXED_STEP', dt, ierr)
    if (ierr /= 0) then
      call abortmp('farksetrin failed')
    endif

    return
  end subroutine update_arkode

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
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp
    use iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parameter_list), target, intent(in)  :: arkode_parameters
    integer,                      intent(in)  :: n0
    real(real_kind),              intent(in)  :: tstart

    ! local variables
    type(parameter_list), pointer :: ap
    integer(C_INT)                :: iatol, ierr
    integer(C_LONG)               :: ival
    integer                       :: ie

    !======= Internals ============
    ap => arkode_parameters

    ! update rtol and atol
    do ie=atol_F%nets,atol_F%nete
      atol_F%elem(ie)%state%v(:,:,1,:,atol_F%tl_idx) = ap%atol(1)
      atol_F%elem(ie)%state%v(:,:,2,:,atol_F%tl_idx) = ap%atol(2)
      atol_F%elem(ie)%state%w_i(:,:,:,atol_F%tl_idx) = ap%atol(3)
      atol_F%elem(ie)%state%phinh_i(:,:,:,atol_F%tl_idx) = ap%atol(4)
      atol_F%elem(ie)%state%vtheta_dp(:,:,:,atol_F%tl_idx) = ap%atol(5)
      atol_F%elem(ie)%state%dp3d(:,:,:,atol_F%tl_idx) = ap%atol(6)
    end do

    ! reinitialize ARKode with current solution and tolerances
    iatol = 2
    call farkreinit(tstart, y_C(n0), ap%imex, iatol, ap%rtol, atol_C, ierr)
    if (ierr /= 0) then
      call abortmp('arkode_reinit: farkreinit failed')
    endif

    ! reset max number of nonlinear iterations per stage
    ival = max_niters
    call farksetiin('MAX_NITERS', ival, ierr)
    if (ierr /= 0) then
      call abortmp('arkode_reinit: farksetinn(MAX_NITERS) failed')
    endif

    return
  end subroutine reinitialize

  !=================================================================

  subroutine initialize(elem, nets, nete, par, n0, qn0, tstart, arkode_parameters)
    !-----------------------------------------------------------------
    ! Description: allocates memory for and initializes ARKode
    !   Arguments:
    !                 elem - (obj*, input) element objects
    !                 nets - (int, input) starting index for elem array
    !                 nete - (int, input) ending index for elem array
    !                  par - (obj*, input) parallel object
    !                   n0 - (int, input) timelevel holding current solution
    !                  qn0 - (int, input) timelevel for tracer mass
    !               tstart - (real, input) time to start ARKode solve at
    !    arkode_parameters - (parameter, input) object for arkode parameters
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use element_mod,      only: element_t
    use HommeNVector,     only: NVec_t, MakeHommeNVector, SetHommeNVectorPar
    use kinds,            only: real_kind
    use parallel_mod,     only: parallel_t, abortmp
    use iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parallel_t),             intent(in) :: par
    type(parameter_list), target, intent(in) :: arkode_parameters
    real(real_kind),              intent(in) :: tstart
    integer,                      intent(in) :: nets, nete, n0, qn0
    type(element_t),              intent(inout) :: elem(:)

    ! local variables
    type(parameter_list), pointer :: ap
    real(real_kind)               :: rpar(1)
    real(real_kind)               :: A_C1(arkode_parameters%s*arkode_parameters%s)
    real(real_kind)               :: A_C2(arkode_parameters%s*arkode_parameters%s)
    integer(C_INT)                :: idef, iatol, ierr
    integer(C_LONG)               :: ipar(1), ival
    integer                       :: i, j

    !======= Internals ============
    ap => arkode_parameters

    ! initialize error flag and placeholders for optional farmalloc outputs
    ierr = 0
    iout = 0
    rout = 0.d0

    if (par%masterproc) print *,"Initializing ARKode"
    ! 'create' NVec_t objects that will correspond to the original 3 HOMME
    ! timelevels, assuming that tl%nm1, tl%n0, and tl%np1 are taken from the
    ! set {1,2,3}
    do i=1,3
      call MakeHommeNVector(elem, nets, nete, i, y_F(i), ierr)
      if (ierr /= 0) then
        call abortmp('Error in MakeHommeNVector')
      end if
      ! get C pointer
      y_C(i) = c_loc(y_F(i))
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
    call MakeHommeNVector(elem, nets, nete, 4, atol_F, ierr)
    if (ierr /= 0) then
      call abortmp('Error in MakeHommeNVector')
    end if
    ! get C pointer
    atol_C = c_loc(atol_F)

    ! initialize ARKode data & operators
    idef = 4  ! flag specifying which SUNDIALS solver will be used (4=ARKode)
    call fnvextinit(idef, ierr)
    if (ierr /= 0) then
       call abortmp('arkode_init: fnvextinit failed')
    end if

    ! ARKode dataspace
    iatol = 2
    call farkmalloc(tstart, y_C(n0), ap%imex, iatol, ap%rtol, atol_C, &
                    iout, rout, ipar, rpar, ierr)
    if (ierr /= 0) then
       call abortmp('arkode_init: farkmalloc failed')
    end if

    ! Set IRK Butcher table for implicit problems
    if (ap%imex == 0) then
      ! flatten Ai to C array (row-major)
      do i=1,ap%s
        do j = 1,ap%s
          A_C1(ap%s*(i-1)+j) = ap%Ai(i,j)
        end do
      end do
      ! set table
      call farksetirktable(ap%s, ap%q, ap%p, ap%ci, A_C1, ap%bi, ap%bi2, ierr)
      if (ierr /= 0) then
        call abortmp('arkode_init: farksetirktable failed')
      end if

    ! Set ERK Butcher table for explicit problems
    else if (ap%imex == 1) then
      ! flatten Ae to C array (row-major)
      do i=1,ap%s
        do j = 1,ap%s
          A_C1(ap%s*(i-1)+j) = ap%Ae(i,j)
        end do
      end do
      ! set table
      call farkseterktable(ap%s, ap%q, ap%p, ap%ce, A_C1, ap%be, ap%be2, ierr)
      if (ierr /= 0) then
        call abortmp('arkode_init: farkseterktable failed')
      end if

    ! Set ARK Butcher table for IMEX problems
    else if (ap%imex == 2) then
      ! flatten Ae and Ai to C arrays (row-major)
      do i=1,ap%s
        do j = 1,ap%s
          A_C1(ap%s*(i-1)+j) = ap%Ai(i,j)
          A_C2(ap%s*(i-1)+j) = ap%Ae(i,j)
        end do
      end do
      ! set tables
      call farksetarktables(ap%s, ap%q, ap%p, ap%ci, ap%ce, A_C1, A_C2, &
                          ap%bi, ap%be, ap%bi2, ap%be2, ierr)
      if (ierr /= 0) then
        call abortmp('arkode_init: farksetarktable failed')
      end if
    else
      call abortmp('arkode_init: invalid imex parameter value')
    end if

    ! Set linear solve if implicit or imex problem
    if (ap%imex == 0 .or. ap%imex == 2) then
    !      To indicate that the implicit problem is linear, make the following
    !      call.  The argument specifies whether the linearly implicit problem
    !      changes as the problem evolves (1) or not (0)
    !  lidef = 0
  !  call farksetiin('LINEAR', lidef, ierr)
  !  if (ierr /= 0) then
  !     write(0,*) ' arkode_init: farksetiin failed'
  !  endif


      if (use_column_solver) then
      ! use the HOMME columnwise direct solver
        call FColumnSolInit(ierr)
        if (ierr /= 0) then
          call abortmp('arkode_init: FColumnSolInit failed')
        end if

      else
      ! use the GMRES linear solver (and set options)
        idef = 4 ! flag specifying which SUNDIALS solver will be used (4=ARKode)
        call FSunSPGMRInit(idef, ap%precLR, ap%maxl, ierr)
        if (ierr /= 0) then
          call abortmp('arkode_init: FSunSPGMRInit failed')
        end if
        call FSunSPGMRSetGSType(idef, ap%gstype, ierr)
        if (ierr /= 0) then
          call abortmp('arkode_init: FSunSPGMRSetGSType failed')
        end if
        call FARKSpilsInit(ierr)
        if (ierr /= 0) then
          call abortmp('arkode_init: FARKSpilsInit failed')
        end if
        call FARKSpilsSetEpsLin(ap%lintol, ierr)
        if (ierr /= 0) then
          call abortmp('arkode_init: FARKSpilsSetEpsLin failed')
        end if

        !      Indicate to use our own preconditioner setup/solve routines (otherwise
        !      preconditioning is disabled)
        if (ap%precLR /= 0) then
          idef = 1
          call farkspilssetprec(idef, ierr)
          if (ierr /= 0) then
            write(0,*) ' arkode_init: farkspilssetprec failed'
          endif
        endif

        !      Indicate to use our own Jacobian-vector product routine (otherwise it
        !      uses a finite-difference approximation)
        !idef = 1
        !call farkspilssetjac(idef, ierr)
        !if (ierr /= 0) then
        !   write(0,*) ' arkode_init: farkspilssetjac failed'
        !endif

      endif
    endif

    ! reset max number of nonlinear iterations per stage
    ival = max_niters
    call farksetiin('MAX_NITERS', ival, ierr)
    if (ierr /= 0) then
      call abortmp('arkode_reinit: farksetinn(MAX_NITERS) failed')
    endif

    if (par%masterproc) call farkwriteparameters(ierr)



    return
  end subroutine initialize

  !=================================================================

  subroutine set_Butcher_tables(arkode_parameters, table_name)
    !-----------------------------------------------------------------
    ! Description: sets Butcher tables for ARKode solver
    !   Arguments:
    !    arkode_parameters - (parameter, in/output) object for arkode parameters
    !           table_name - (integer, input) constant identifying table name
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parameter_list), target, intent(inout) :: arkode_parameters
    integer,                      intent(in)    :: table_name

    ! local variables
    type(parameter_list), pointer :: ap
    real(real_kind) :: beta, gamma, delta, b1, b2
    real(real_kind) :: a(max_stage_num), ahat(max_stage_num)
    real(real_kind) :: b(max_stage_num), dhat(max_stage_num)

    !======= Internals ============
    ap => arkode_parameters

    select case (table_name)

    case (RK2_ARK)
        ap%imex = 1 ! explicit
        ap%s = 2 ! 2 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:2,1:2) = 0.d0
        ap%Ae(2,1) = 0.5d0
        ! Explicit Butcher Table (vectors)
        ap%ce(1:2) = (/ 0.d0, 0.5d0 /)
        ap%be(1:2) = (/ 0.d0, 1.d0 /)

      case (KGU35_ARK)
        ap%imex = 1 ! explicit
        ap%s = 5 ! 5 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:5,1:5) = 0.d0
        ap%Ae(2,1) = 0.2d0
        ap%Ae(3,1:2) = (/  0.d0, 0.2d0 /)
        ap%Ae(4,1:3) = (/  0.d0,  0.d0, 1.d0/3.d0 /)
        ap%Ae(5,1:4) = (/  0.d0,  0.d0,      0.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:5) = (/   0.d0, 0.2d0, 0.2d0, 1.d0/3.d0, 2.d0/3.d0 /)
        ap%be(1:5) = (/ 0.25d0,  0.d0,  0.d0,      0.d0,    0.75d0 /)

      case (ARS232_ARK)
        ap%imex = 2 ! imex
        ap%s = 3 ! 3 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        delta = -2.d0*sqrt(2.d0)/3.d0
        gamma = 1.d0 - 1.d0/sqrt(2.d0)
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:3,1:3) = 0.d0
        ap%Ai(2,1:2) = (/ 0.d0,      gamma /)
        ap%Ai(3,1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:3) = (/ 0.d0, gamma, 1.d0 /)
        ap%bi(1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:3,1:3) = 0.d0
        ap%Ae(2,1) =  gamma
        ap%Ae(3,1:2) = (/ delta, 1.d0-delta /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:3) = ap%ci(1:3)
        ap%be(1:3) = ap%bi(1:3)

      case (ARK453_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:5,1:5) = 0.d0
        ap%Ai(2,1:2) = (/ -0.22284985318525410d0, 0.32591194130117247d0 /)
        ap%Ai(3,1:3) = (/ -0.46801347074080545d0, 0.86349284225716961d0, &
                            0.32591194130117247d0 /)
        ap%Ai(4,1:4) = (/ -0.46509906651927421d0, 0.81063103116959553d0, &
                            0.61036726756832357d0, 0.32591194130117247d0 /)
        ap%Ai(5,1:5) = (/ 0.87795339639076675d0, -0.72692641526151547d0, &
                            0.75204137157372720d0, -0.22898029400415088d0, &
                            0.32591194130117247d0 /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:5) = (/ 0.d0, 0.1030620881159184d0, &
                          0.72139131281753662d0, 1.28181117351981733d0, &
                          1.d0 /)
        ap%bi(1:5) = (/ 0.87795339639076672d0, -0.72692641526151549d0, &
                          0.7520413715737272d0, -0.22898029400415090d0, &
                          0.32591194130117246d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:5,1:5) = 0.d0
        ap%Ae(2,1) = 0.10306208811591838d0
        ap%Ae(3,1:2) = (/ -0.94124866143519894d0, 1.6626399742527356d0 /)
        ap%Ae(4,1:3) = (/ -1.3670975201437765d0, 1.3815852911016873d0, &
                            1.2673234025619065d0 /)
        ap%Ae(5,1:4) = (/ -0.81287582068772448d0, 0.81223739060505738d0, &
                            0.90644429603699305d0, 0.094194134045674111d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:5) = ap%ci(1:5)
        ap%be(1:5) = ap%bi(1:5)

      case (ARS233_ARK)
        ap%imex = 2 ! imex
        ap%s = 3 ! 3 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        gamma = 1.d0/6.d0*(3.d0 + sqrt(3.d0))
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:3,1:3) = 0.d0
        ap%Ai(2,1:2) = (/ 0.d0,          gamma /)
        ap%Ai(3,1:3) = (/ 0.d0, 1.d0-2.d0*gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:3) = (/ 0.d0, gamma, 1.d0-gamma /)
        ap%bi(1:3) = (/ 0.d0, 0.5d0,      0.5d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:3,1:3) = 0.d0
        ap%Ae(2,1) = gamma
        ap%Ae(3,1:2) = (/ gamma-1.d0, 2.d0*(1.d0-gamma) /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:3) = ap%ci(1:3)
        ap%be(1:3) = ap%bi(1:3)

      case (ARS222_ARK)
        ap%imex = 2 ! imex
        ap%s = 3 ! 3 stage
        ap%q = 2 ! 2rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        gamma = 1.d0 - 1.d0/sqrt(2.d0)
        delta = 1.d0 - 1.d0/(2.d0*gamma)
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:3,1:3) = 0.d0
        ap%Ai(2,1:2) = (/ 0.d0,      gamma /)
        ap%Ai(3,1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:3) = (/ 0.d0,      gamma,  1.d0 /)
        ap%bi(1:3) = (/ 0.d0, 1.d0-gamma, gamma /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:3,1:3) = 0.d0
        ap%Ae(2,1) = gamma
        ap%Ae(3,1:2) = (/ delta, 1.d0-delta /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:3) = ap%ci(1:3)
        ap%be(1:3) = ap%bi(1:3)

      case (ARS343_ARK)
        ap%imex = 2 ! imex
        ap%s = 4 ! 4 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        gamma = 0.4358665215084590d0
        b1 = -3.d0*gamma*gamma/2.d0 + 4.d0*gamma - 1.d0/4.d0
        b2 = 3.d0*gamma*gamma/2.d0 - 5.d0*gamma + 5.d0/4.d0
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:4,1:4) = 0.d0
        ap%Ai(2,1:2) = (/ 0.d0, gamma /)
        ap%Ai(3,1:3) = (/ 0.d0, (1.d0-gamma)/2.d0, gamma /)
        ap%Ai(4,1:4) = (/ 0.d0,                b1,    b2, gamma /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:4) = (/ 0.d0, gamma, (1.d0+gamma)/2.d0,  1.d0 /)
        ap%bi(1:4) = (/ 0.d0,    b1,                b2, gamma /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:4,1:4) = 0.d0
        ap%Ae(2,1) = gamma
        ap%Ae(4,2) = 0.5529291480359398d0
        ap%Ae(4,3) = 0.5529291480359398d0
        ap%Ae(4,1) = 1.d0 - ap%Ae(4,2) - ap%Ae(4,3)
        ap%Ae(3,1) = ap%Ae(4,2)*(2.d0 - 9.d0*gamma + 3.d0*gamma*gamma)/2.d0 &
                    +ap%Ae(4,3)*(11.d0 - 42.d0*gamma + 15.d0*gamma*gamma)/4.d0 &
                    -7.d0/2.d0 + 13.d0*gamma - 9.d0*gamma*gamma/2.d0
        ap%Ae(3,2) = ap%Ae(4,2)*(-2.d0 + 9.d0*gamma - 3.d0*gamma*gamma)/2.d0 &
                    +ap%Ae(4,3)*(-11.d0 + 42.d0*gamma - 15.d0*gamma*gamma)/4.d0 &
                    +4.d0 - 25.d0*gamma/2.d0 + 9.d0*gamma*gamma/2.d0
        ! Explicit Butcher Table (vectors)
        ap%ce(1:4) = ap%ci(1:4)
        ap%be(1:4) = ap%bi(1:4)

      case (ARS443_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:5,1:5) = 0.d0
        ap%Ai(2,1:2) = (/ 0.d0,  1.d0/2.d0 /)
        ap%Ai(3,1:3) = (/ 0.d0,  1.d0/6.d0,  1.d0/2.d0 /)
        ap%Ai(4,1:4) = (/ 0.d0, -1.d0/2.d0,  1.d0/2.d0, 1.d0/2.d0 /)
        ap%Ai(5,1:5) = (/ 0.d0,  3.d0/2.d0, -3.d0/2.d0, 1.d0/2.d0, 1.d0/2.d0 /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:5) = (/ 0.d0, 1.d0/2.d0,  2.d0/3.d0, 1.d0/2.d0,      1.d0 /)
        ap%bi(1:5) = (/ 0.d0, 3.d0/2.d0, -3.d0/2.d0, 1.d0/2.d0, 1.d0/2.d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:5,1:5) = 0.d0
        ap%Ae(2,1) = 1.d0/2.d0
        ap%Ae(3,1:2) = (/ 11.d0/18.d0, 1.d0/18.d0 /)
        ap%Ae(4,1:3) = (/   5.d0/6.d0, -5.d0/6.d0, 1.d0/2.d0 /)
        ap%Ae(5,1:4) = (/   1.d0/4.d0,  7.d0/4.d0, 3.d0/4.d0, -7.d0/4.d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:5) = ap%ci(1:5)
        ap%be(1:5) = (/ 1.d0/4.d0, 7.d0/4.d0, 3.d0/4.d0, -7.d0/4.d0, 0.d0 /)

      case (ARK324_ARK)
        ap%imex = 2 ! imex
        ap%s = 4 ! 4 stage
        ap%q = 3 ! 3rd order
        ap%p = 2 ! 2nd order embedding
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:4,1:4) = 0.d0
        ap%Ai(2,1:2) = (/ 1767732205903.d0/4055673282236.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ap%Ai(3,1:3) = (/ 2746238789719.d0/10658868560708.d0, &
                         -640167445237.d0/6845629431997.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ap%Ai(4,1:4) = (/ 1471266399579.d0/7840856788654.d0, &
                         -4482444167858.d0/7529755066697.d0, &
                          11266239266428.d0/11593286722821.d0, &
                          1767732205903.d0/4055673282236.d0 /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:4) = (/ 0.d0, 1767732205903.d0/2027836641118.d0, 3.d0/5.d0, 1.d0 /)
        ap%bi(1) = 1471266399579.d0/7840856788654.d0
        ap%bi(2) = -4482444167858.d0/7529755066697.d0
        ap%bi(3) =  11266239266428.d0/11593286722821.d0
        ap%bi(4) = 1767732205903.d0/4055673282236.d0
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:4,1:4) = 0.d0
        ap%Ae(2,1) = 1767732205903.d0/2027836641118.d0
        ap%Ae(3,1:2) = (/ 5535828885825.d0/10492691773637.d0, &
                          788022342437.d0/10882634858940.d0 /)
        ap%Ae(4,1:3) = (/ 6485989280629.d0/16251701735622.d0, &
                         -4246266847089.d0/9704473918619.d0, &
                          10755448449292.d0/10357097424841.d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:4) = ap%ci(1:4)
        ap%be(1:4) = ap%bi(1:4)
        ! Embedding
        ap%be2(1) = 2756255671327.d0/12835298489170.d0
        ap%be2(2) = -10771552573575.d0/22201958757719.d0
        ap%be2(3) = 9247589265047.d0/10645013368117.d0
        ap%be2(4) = 2193209047091.d0/5459859503100.d0

      case (ARK436_ARK)
        ap%imex = 2 ! imex
        ap%s = 6 ! 6 stage
        ap%q = 4 ! 4th order
        ap%p = 3 ! 3rd order embedding
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:6,1:6) = 0.d0
        ap%Ai(2,1:2) = (/ 1.d0/4.d0, 1.d0/4.d0 /)
        ap%Ai(3,1:3) = (/ 8611.d0/62500.d0, -1743.d0/31250.d0, 1.d0/4.d0 /)
        ap%Ai(4,1:4) = (/ 5012029.d0/34652500.d0, -654441.d0/2922500.d0, &
                          174375.d0/388108.d0, 1.d0/4.d0 /)
        ap%Ai(5,1:5) = (/ 15267082809.d0/155376265600.d0, &
                          -71443401.d0/120774400.d0, 730878875.d0/902184768.d0, &
                          2285395.d0/8070912.d0, 1.d0/4.d0 /)
        ap%Ai(6,1:6) = (/ 82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, &
                          69875.d0/102672.d0, -2260.d0/8211.d0, 1.d0/4.d0 /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:6) = (/ 0.d0, 1.d0/2.d0, 83.d0/250.d0, 31.d0/50.d0, &
                      17.d0/20.d0, 1.d0 /)
        ap%bi(1:6) = (/ 82889.d0/524892.d0, 0.d0, 15625.d0/83664.d0, &
                      69875.d0/102672.d0, -2260.d0/8211.d0, 1.d0/4.d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:6,1:6) = 0.d0
        ap%Ae(2,1) = 1.d0/2.d0
        ap%Ae(3,1:2) = (/ 13861.d0/62500.d0, 6889.d0/62500.d0 /)
        ap%Ae(4,1:3) = (/ -116923316275.d0/2393684061468.d0, &
                          -2731218467317.d0/15368042101831.d0, &
                          9408046702089.d0/11113171139209.d0 /)
        ap%Ae(5,1:4) = (/ -451086348788.d0/2902428689909.d0, &
                          -2682348792572.d0/7519795681897.d0, &
                          12662868775082.d0/11960479115383.d0, &
                          3355817975965.d0/11060851509271.d0 /)
        ap%Ae(6,1:5) = (/ 647845179188.d0/3216320057751.d0, &
                          73281519250.d0/8382639484533.d0, &
                          552539513391.d0/3454668386233.d0, &
                          3354512671639.d0/8306763924573.d0, 4040.d0/17871.d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:6) = ap%ci(1:6)
        ap%be(1:6) = ap%bi(1:6)
        ! Embedding
        ap%be2(1:6) = (/ 4586570599.d0/29645900160.d0, 0.d0, 178811875.d0/945068544.d0, &
                        814220225.d0/1159782912.d0, -3700637.d0/11593932.d0, &
                        61727.d0/225920.d0 /)

      case (SSP3333B_ARK)
        ap%imex = 2 ! imex
        ap%s = 3 ! 3 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:3,1:3) = 0.d0
        ap%Ai(2,1:2) = (/      0.d0,       1.d0 /)
        ap%Ai(3,1:3) = (/ 1.d0/6.d0, -1.d0/3.d0, 2.d0/3.d0 /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:3) = (/ 0.d0, 1.d0, 0.5d0 /)
        ap%bi(1:3) = (/ 1.d0/6.d0, 1.d0/6.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:3,1:3) = 0.d0
        ap%Ae(2,1) = 1.d0
        ap%Ae(3,1:2) = (/ 0.25d0, 0.25d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:3) = ap%ci(1:3)
        ap%be(1:3) = ap%bi(1:3)

      case (SSP3333C_ARK)
        ap%imex = 2 ! imex
        ap%s = 3 ! 3 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        beta = sqrt(3.d0)/5.d0 + 0.5d0
        gamma = -1.d0/8.d0*(sqrt(3.d0) + 1.d0)
        ! Implicit Butcher Table (matrix)
        ap%Ai(1:3,1:3) = 0.d0
        ap%Ai(2,1:2) = (/ 4.d0*gamma + 2.d0*beta, 1.d0 - 4.d0*gamma - 2.d0*beta /)
        ap%Ai(3,1:3) = (/   0.5d0 - beta - gamma,                         gamma, beta /)
        ! Implicit Butcher Table (vectors)
        ap%ci(1:3) = (/ 0.d0, 1.d0, 0.5d0 /)
        ap%bi(1:3) = (/ 1.d0/6.d0, 1.d0/6.d0, 2.d0/3.d0 /)
        ! Explicit Butcher Table (matrix)
        ap%Ae(1:3,1:3) = 0.d0
        ap%Ae(2,1) = 1.d0
        ap%Ae(3,1:2) = (/ 0.25d0, 0.25d0 /)
        ! Explicit Butcher Table (vectors)
        ap%ce(1:3) = ap%ci(1:3)
        ap%be(1:3) = ap%bi(1:3)

      case (IMKG232_ARK)
        ap%imex = 2 ! imex
        ap%s = 4 ! 4 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.d0))
        a(1:3) = (/ 0.5d0, 0.5d0, 1.d0 /)
        ahat(1:3) = (/ 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:2) = (/ delta, delta /)
        b(1:2) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG242_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.d0))
        a(1:4) = (/ 0.25d0, 1.d0/3.d0, 0.5d0, 1.d0 /)
        ahat(1:4) = (/ 0.d0, 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:3) = (/ 0.d0, delta, delta /)
        b(1:3) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG243_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0+sqrt(3.d0)/6.d0
        a(1:4) = (/ 0.25d0, 1.d0/3.d0, 0.5d0, 1.d0 /)
        ahat(1:4) = (/ 0.d0, 1.d0/6.d0, -sqrt(3.d0)/6.d0, 1.d0 /)
        dhat(1:3) = (/ delta, delta, delta /)
        b(1:3) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG252_ARK)
        ap%imex = 2 ! imex
        ap%s = 6 ! 6 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(2.d0+sqrt(2.0))
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, 0.d0, 0.d0, -0.5d0*(sqrt(2.d0)+1.d0), 1.d0 /)
        dhat(1:4) = (/ 0.d0, 0.d0, delta, delta /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG253_ARK)
        ap%imex = 2 ! imex
        ap%s = 6 ! 6 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0+sqrt(3.d0)/6.d0
        gamma = 0.25d0*sqrt(3.d0)*(1.d0+sqrt(3.d0)/3.d0)*((sqrt(3.d0)/3.d0-1.d0)**2-2.d0)
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, 0.d0, gamma, -sqrt(3.d0)/6.d0, 1.d0 /)
        dhat(1:4) = (/ 0.d0, delta, delta, delta /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG254_ARK)
        ap%imex = 2 ! imex
        ap%s = 6 ! 6 stage
        ap%q = 2 ! 2nd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:5) = (/ 0.25d0, 1.d0/6.d0, 3.d0/8.d0, 0.5d0, 1.d0 /)
        ahat(1:5) = (/ 0.d0, -1.d0/20.d0, 1.25d0, -0.5d0, 1.d0 /)
        dhat(1:4) = (/ -0.5d0, 1.d0, 1.d0, 1.d0 /)
        b(1:4) = 0.d0
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG342_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        delta = 0.5d0*(1.d0+sqrt(3.d0)/3.d0)
        a(1:4) = (/ 0.25d0, 2.d0/3.d0, 1.d0/3.d0, 0.75d0 /)
        ahat(1:4) = (/ 0.d0, (1.d0-sqrt(3.d0))/6.d0, -(1.d0+sqrt(3.d0))/6.d0, 0.75d0 /)
        dhat(1:3) = (/ 0.d0, delta, delta /)
        b(1:3) = (/ 0.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)


      case (IMKG343_ARK)
        ap%imex = 2 ! imex
        ap%s = 5 ! 5 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:4) = (/ 0.25d0, 2.d0/3.d0, 1.d0/3.d0, 0.75d0 /)
        ahat(1:4) = (/ 0.d0, -1.d0/3.d0, -2.d0/3.d0, 0.75d0 /)
        dhat(1:3) = (/ -1.d0/3.d0, 1.d0, 1.d0 /)
        b(1:3) = (/ 0.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)

      case (IMKG353_ARK)
        ap%imex = 2 ! imex
        ap%s = 6 ! 6 stage
        ap%q = 3 ! 3rd order
        ap%p = 0 ! no embedded order
        ap%be2 = 0.d0 ! no embedded explicit method
        ! IMEX-KG vectors
        a(1:5) = (/ -0.017391304347826087d0, -.92d0, 5.d0/3.d0, 1.d0/3.d0, 3.d0/4.d0 /)
        ahat(1:2) = (/ 0.3075640504095504d0, -1.2990164859879263d0 /)
        ahat(3:5) = (/ 1.2516666666666665d0, -0.8166666666666668d0, 3d0/4d0/)
        dhat(1:4) = (/ -0.2981612530370581d0, .415d0, .415d0, 1.15d0 /)
        b(1:4) = (/ 1.d0, -1.d0, 1.d0/3.d0, 0.25d0 /)
        ! set IMEX-KG Butcher table
        call set_IMKG_Butcher_tables(arkode_parameters, ap%s, a, ahat, dhat, b)


      case default
        call abortmp('Unknown ARKode Butcher table name')
    end select

  end subroutine set_Butcher_tables
  !=================================================================

  subroutine set_IMKG_Butcher_tables(arkode_parameters, s, a, ahat, dhat, b)
    !-----------------------------------------------------------------
    ! Description: sets Butcher tables in IMEX-KG format:
    !   Ae:   0 |    0             Ai:     0 |       0
    !        c1 |   a1    0            chat1 |   ahat1 dhat1
    !        c2 |   b1   a2  0         chat2 |    b1 ahat2 dhat2
    !        c3 |   b2    0 a3 0       chat3 |    b2     0 ahat3 dhat3
    !           |                            |
    !        cq | bq-1 0 ... 0 aq 0    chatq |  bq-1     0    ...    0 ahatq 0
    !        ----------------------    -----------------------------------------
    !           | bq-1 0 ... 0 aq 0          |  bq-1     0    ...    0 ahatq 0
    !
    !        ci = sum_j Ae(i,j)        chati = sum_j Ai(i,j)
    !
    !   Arguments:
    !    arkode_parameters - (parameter, in/output) object for arkode parameters
    !                    s - (integer, input) number of stages
    !                    a - (real(max_stage_num), input) alpha vector
    !                 ahat - (real(max_stage_num), input) alpha_hat vector
    !                 dhat - (real(max_stage_num), in/output) delta_hat vector
    !                    b - (real(max_stage_num), input) beta vector
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use kinds,            only: real_kind
    use parallel_mod,     only: abortmp

    !======= Declarations =========
    implicit none

    ! calling variables
    type(parameter_list), target, intent(inout) :: arkode_parameters
    integer,                      intent(in)    :: s
    real(real_kind),              intent(in)    :: a(max_stage_num)
    real(real_kind),              intent(in)    :: ahat(max_stage_num)
    real(real_kind),              intent(inout) :: dhat(max_stage_num)
    real(real_kind),              intent(in)    :: b(max_stage_num)

    ! local variables
    type(parameter_list), pointer :: ap

    !======= Internals ============
    ap => arkode_parameters

    ! for easier implementation, extend dhat vector by setting s-1 index to zero
    dhat(s-1) = 0.d0

    ! set matrices and c vectors according to IMEX-KG format
    ap%Ai(1:s,1:s) = 0.d0
    ap%Ae(1:s,1:s) = 0.d0
    ap%ci(1) = 0.d0
    ap%ce(1) = 0.d0
    if (s > 1) then
      ap%Ai(2,1:2) = (/ ahat(1), dhat(1) /)
      ap%Ae(2,1) = a(1)
      ap%ci(2) = ahat(1)+dhat(1)
      ap%ce(2) = a(1)
    end if
    if (s > 2) then
      ap%Ai(3,1:3) = (/ b(1), ahat(2), dhat(2) /)
      ap%Ae(3,1:2) = (/ b(1), a(2) /)
      ap%ci(3) = b(1)+ahat(2)+dhat(2)
      ap%ce(3) = b(1)+a(2)
    end if
    if (s > 3) then
      ap%Ai(4,1:4) = (/ b(2), 0.d0, ahat(3), dhat(3) /)
      ap%Ae(4,1:3) = (/ b(2), 0.d0, a(3) /)
      ap%ci(4) = b(2)+ahat(3)+dhat(3)
      ap%ce(4) = b(2)+a(3)
    end if
    if (s > 4) then
      ap%Ai(5,1:5) = (/ b(3), 0.d0, 0.d0, ahat(4), dhat(4) /)
      ap%Ae(5,1:4) = (/ b(3), 0.d0, 0.d0, a(4) /)
      ap%ci(5) = b(3)+ahat(4)+dhat(4)
      ap%ce(5) = b(3)+a(4)
    end if
    if (s > 5) then
      ap%Ai(6,1:6) = (/ b(4), 0.d0, 0.d0, 0.d0, ahat(5), dhat(5) /)
      ap%Ae(6,1:5) = (/ b(4), 0.d0, 0.d0, 0.d0, a(5) /)
      ap%ci(6) = b(4)+ahat(5)+dhat(5)
      ap%ce(6) = b(4)+a(5)
    end if
    if (s > 6) then
      call abortmp('ARKode IMEX-KG only implemented for 6 stages or less')
    end if

    ! set b vectors
    ap%bi(1:s) = ap%Ai(s,1:s)
    ap%be(1:s) = ap%Ae(s,1:s)

  end subroutine set_IMKG_Butcher_tables

end module arkode_mod
