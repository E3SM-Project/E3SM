#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module arkode_mod

  use element_state,  only: timelevels
  use derivative_mod, only: derivative_t
  use hybrid_mod,     only: hybrid_t
  use hybvcoord_mod,  only: hvcoord_t
  use kinds,          only: real_kind
  use iso_c_binding

  implicit none

  private

  integer, parameter :: max_stage_num = 10
  integer, parameter :: freelevels = timelevels-30
  ! Note that 30 is an estimate, but needs to be at least > 25
  ! If a larger Krylov subspace is desired, timelevels should be
  ! increased.

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
    ! GMRES Linear Solver Info
    integer         :: precLR ! preconditioning: 0=none, 1=left, 2=right, 3=left+right
    integer         :: gstype ! Gram-Schmidt orthogonalization: 1=modified, 2=classical
    integer         :: maxl = freelevels ! max size of Krylov subspace (# of iterations/vectors)
    real(real_kind) :: lintol ! linear convergence tolerance factor (0 indicates default)
    ! General Iteration Info
    integer         :: iatol=1 ! indices of atol to use: 1=1, 2=all
    real(real_kind) :: rtol ! relative tolerance for iteration convergence
    real(real_kind) :: atol(6) ! absolute tolerances (u,v,w,phinh,theta_dp_cp,dp3d)
  end type parameter_list

  public :: parameter_list, update_arkode, get_solution_ptr, get_RHS_vars

  save

  type(hvcoord_t), pointer    :: hvcoord_ptr
  type(hybrid_t), pointer     :: hybrid_ptr
  type(derivative_t), pointer :: deriv_ptr
  type(c_ptr)                 :: y_C(3), atol_C
  real(real_kind)             :: dt_save
  real(real_kind)             :: eta_ave_w_save
  integer                     :: imex_save, qn0_save

  logical :: initialized = .false.

contains

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

  subroutine get_RHS_vars(imex, qn0, dt, eta_ave_w, hvcoord, hybrid, deriv)
    !-----------------------------------------------------------------
    ! Description: sets variables and objects needed to compute RHS
    !   Arguments:
    !     hvcoord - (obj*, output) hvcoord object pointer
    !      hybrid - (obj*, output) hybrid object pointer
    !       deriv - (obj*, output) deriv object pointer
    !-----------------------------------------------------------------

    !======= Inclusions ===========
    use derivative_mod, only: derivative_t
    use hybrid_mod,     only: hybrid_t
    use hybvcoord_mod,  only: hvcoord_t
    use kinds,          only: real_kind

    !======= Declarations =========
    implicit none

    ! calling variables
    type(hvcoord_t),    intent(out) :: hvcoord
    type(hybrid_t),     intent(out) :: hybrid
    type(derivative_t), intent(out) :: deriv
    integer,            intent(out) :: imex, qn0
    real(real_kind),    intent(out) :: dt, eta_ave_w

    !======= Internals ============
    imex = imex_save
    qn0 = qn0_save
    dt = dt_save
    eta_ave_w = eta_ave_w_save
    hvcoord = hvcoord_ptr
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
    integer(C_INT)  :: ierr
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
    type(parameter_list),  intent(in)  :: arkode_parameters
    integer,               intent(in)  :: n0
    real(real_kind),       intent(in)  :: tstart

    ! local variables
    integer(C_INT)  :: ierr

    !======= Internals ============
    if (arkode_parameters%iatol == 1) then
      call farkreinit(tstart, y_C(n0), arkode_parameters%imex, arkode_parameters%iatol, &
                      arkode_parameters%rtol, arkode_parameters%atol(1), ierr)
    else if (arkode_parameters%iatol == 2) then
      call farkreinit(tstart, y_C(n0), arkode_parameters%imex, arkode_parameters%iatol, &
                      arkode_parameters%rtol, atol_C, ierr)
    end if
    if (ierr /= 0) then
      call abortmp('arkode_init: farkreinit failed')
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
    type(NVec_t), target          :: y(3), z
    type(parameter_list), pointer :: ap
    real(real_kind)               :: rout(40), rpar(1)
    real(real_kind)               :: A_C1(arkode_parameters%s*arkode_parameters%s)
    real(real_kind)               :: A_C2(arkode_parameters%s*arkode_parameters%s)
    integer(C_INT)                :: idef, ierr
    integer(C_LONG)               :: iout(40), ipar(1)
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
      call MakeHommeNVector(elem, nets, nete, i, y(i), ierr)
      if (ierr /= 0) then
        call abortmp('Error in MakeHommeNVector')
      end if
      ! get C pointer
      y_C(i) = c_loc(y(i))
    end do

    if (ap%iatol == 2) then
      ! set data in 4th timelevel to atol values and 'create' NVec_t object
      ! NOTE: this is where one could implement spatially targetted convergence criteria
      do i=nets,nete
        elem(i)%state%v(:,:,1,:,4) = ap%atol(1)
        elem(i)%state%v(:,:,2,:,4) = ap%atol(2)
        elem(i)%state%w(:,:,:,4) = ap%atol(3)
        elem(i)%state%phinh(:,:,:,4) = ap%atol(4)
        elem(i)%state%theta_dp_cp(:,:,:,4) = ap%atol(5)
        elem(i)%state%dp3d(:,:,:,4) = ap%atol(6)
      end do
      call MakeHommeNVector(elem, nets, nete, 4, z, ierr)
      if (ierr /= 0) then
        call abortmp('Error in MakeHommeNVector')
      end if
      ! get C pointer
      atol_C = c_loc(z)
    end if

    ! initialize ARKode data & operators
    idef = 4  ! flag specifying which SUNDIALS solver will be used (4=ARKode)
    call fnvextinit(idef, ierr)
    if (ierr /= 0) then
       call abortmp('arkode_init: fnvextinit failed')
    end if

    ! ARKode dataspace
    if (ap%iatol == 1) then
      call farkmalloc(tstart, y_C(n0), ap%imex, ap%iatol, ap%rtol, ap%atol(1), &
                      iout, rout, ipar, rpar, ierr)
    else if (ap%iatol == 2) then
      call farkmalloc(tstart, y_C(n0), ap%imex, ap%iatol, ap%rtol, atol_C, &
                      iout, rout, ipar, rpar, ierr)
    end if
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

      ! indicate use of GMRES linear solver
      call farkspgmr(ap%precLR, ap%gstype, ap%maxl, ap%lintol, ierr)
      if (ierr /= 0) then
        call abortmp('arkode_init: farkspgmr failed')
      end if

      !      Indicate to use our own Jacobian-vector product routine (otherwise it
      !      uses a finite-difference approximation)
      !idef = 1
      !call farkspilssetjac(idef, ierr)
      !if (ierr /= 0) then
      !   write(0,*) ' arkode_init: farkspilssetjac failed'
      !endif

      !      Indicate to use our own preconditioner setup/solve routines (otherwise
      !      preconditioning is disabled)
      !idef = 1
      !call farkspilssetprec(idef, ierr)
      !if (ierr /= 0) then
      !   write(0,*) ' arkode_init: farkspilssetprec failed'
      !endif
    end if

    if (par%masterproc) call farkwriteparameters(ierr)



    return
  end subroutine initialize

  !=================================================================

end module arkode_mod
