!-------------------------------------------------------------------------------
! Daniel R. Reynolds, SMU Mathematics
! David J. Gardner, LLNL
! Copyright 2017; all rights reserved
!-------------------------------------------------------------------------------
! This is the implementation file for a Fortran columnwise linear
! solver package for HOMME+ARKode.  This is specifically created
! to pair with the NVECTOR_EXTERNAL vector implementation.  It is
! assumed that the linear solver 'setup' is empty, and that the
! 'solve' routine internally sets up the linear system matrix
! on-the-fly and uses it in the solve process.
!-------------------------------------------------------------------------------

module HommeSUNLinSol
  !-----------------------------------------------------------------------------
  ! Fortran user-defined SUNLinearSolver interface
  !-----------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fsundials_linearsolver_mod

  use HommeNVector, only: NVec_t, FN_VGetContent

  save

  type(c_ptr), pointer :: arkmem ! arkode memory
  type(c_ptr)          :: ycur_ptr

contains

  function FSUNMatrix_HOMME() result(sunmat)
    !---------------------------------------------------------------------------
    ! Create a SUNDIALS matrix obejct. This is a dummy object that does nothing
    ! but needs to be provided to ARKode since it expects a direct linear solver
    ! object to have a corresponding matrix object.
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNMatrix), pointer     :: sunmat ! matrix object
    type(SUNMatrix_Ops), pointer :: matops ! matrix operations

    ! allocate SUNMatrix structure
    sunmat => FSUNMatNewEmpty()

    ! access the SUNMatrix ops structure
    call c_f_pointer(sunmat%ops, matops)

    ! set matrix operations
    matops%getid   = c_funloc(FSUNMatGetID_HOMME)
    matops%clone   = c_funloc(FSUNMatClone_HOMME)
    matops%destroy = c_funloc(FSUNMatDestroy_HOMME)

  end function FSUNMatrix_HOMME

  !=============================================================================

  integer(SUNMatrix_ID) function FSUNMatGetID_HOMME(sunmat_a) result(id) bind(C)
    !---------------------------------------------------------------------------
    ! Return the matrix ID
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNMatrix) :: sunmat_a

    id = SUNMATRIX_CUSTOM

    return

  end function FSUNMatGetID_HOMME

  !=============================================================================

  function FSUNMatClone_HOMME(sunmat_a) result(b_ptr) bind(C)
    !---------------------------------------------------------------------------
    ! Clone the matrix (dummy operation, makes an empty clone)
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNMatrix) :: sunmat_a
    type(c_ptr)     :: b_ptr

    type(SUNMatrix), pointer :: sunmat_b
    integer(c_int),  pointer :: retval

    ! allocate SUNMatrix structure
    sunmat_b => FSUNMatNewEmpty()

    ! copy operations from a into b
    retval = FSUNMatCopyOps(sunmat_a, sunmat_b)

    ! set the c_ptr output
    b_ptr = c_loc(sunmat_b)

    return

  end function FSUNMatClone_HOMME

  !=============================================================================

  subroutine FSUNMatDestroy_HOMME(sunmat)
    !---------------------------------------------------------------------------
    ! Free the matrix
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNMatrix) :: sunmat

    ! deallocate overall SUNMatrix structure
    call FSUNMatFreeEmpty(sunmat)

  end subroutine FSUNMatDestroy_HOMME

  !=============================================================================

  function FSUNLinSol_HOMME(arkode_mem) result(sunls)
    !---------------------------------------------------------------------------
    ! Create a SUNDIALS linear solver object that wraps the native HOMME column
    ! solver
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use farkode_arkstep_mod

    implicit none

    type(c_ptr), target, intent(in)    :: arkode_mem
    type(SUNLinearSolver), pointer     :: sunls
    type(SUNLinearSolver_Ops), pointer :: lsops

    !======= Internals ============

    ! save arkode mem to make ARKStep calls in the solve function
    arkmem => arkode_mem

    ! allocate N_Vector* to obtain current solution in solve function
    ycur_ptr = FN_VNewVectorArray(1)

    ! allocate SUNLinearSolver structure
    sunls => FSUNLinSolNewEmpty()

    ! access the SUNLinearSolver ops structure
    call c_f_pointer(sunls%ops, lsops)

    ! set linear solver operations
    lsops%gettype = c_funloc(FSUNLinSolGetType_HOMME)
    lsops%solve   = c_funloc(FSUNLinSolSolve_HOMME)
    lsops%free    = c_funloc(FSUNLinSolFree_HOMME)

  end function FSUNLinSol_HOMME

  !=============================================================================

  integer(SUNLinearSolver_Type) function FSUNLinSolGetType_HOMME(sunls) &
       result(type) bind(C)
    !---------------------------------------------------------------------------
    ! Return the linear solver type
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNLinearSolver) :: sunls

    type = SUNLINEARSOLVER_DIRECT

    return

  end function FSUNLinSolGetType_HOMME

  !=============================================================================

  integer(c_int) function FSUNLinSolSolve_HOMME(sunls_s, sunmat_a, sunvec_x, sunvec_b, tol) &
       result(ier) bind(C)
    !---------------------------------------------------------------------------
    ! SUNLinSolSolve_HOMME is the C interface routine to call the
    ! Fortran-supplied FCOLUMNSOL_SOLVE routine.
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use farkode_arkstep_mod
    use HommeNVector, only: NVec_t, FN_VGetContent

    implicit none

    type(SUNLinearSolver) :: sunls_s
    type(SUNMatrix)       :: sunmat_a
    type(N_Vector)        :: sunvec_x
    type(N_Vector)        :: sunvec_b
    real(c_double), value :: tol

    type(N_Vector), pointer :: sunvec_ycur
    real(c_double)          :: tcur(1)
    real(c_double)          :: gamma(1)
    type(NVec_t), pointer   :: ycur
    type(Nvec_t), pointer   :: b
    type(NVec_t), pointer   :: x
    integer                 :: FColumnSolSolve

    ! initialize return flag to success
    ier = 0

    ! access current time
    ier = FARKStepGetCurrentTime(arkmem, tcur)

    ! access current state
    ier = FARKStepGetCurrentState(arkmem, ycur_ptr)
    sunvec_ycur => FN_VGetVecAtIndexVectorArray(ycur_ptr, 0)

    ! access current gamma value
    ier = FARKStepGetCurrentGamma(arkmem, gamma)

    ! call Fortran routine to do operation
    ycur => FN_VGetContent(sunvec_ycur)
    b => FN_VGetContent(sunvec_b)
    x => FN_VGetContent(sunvec_x)
    ier = FColumnSolSolve(b, tcur, ycur, gamma, x)

    return

  end function FSUNLinSolSolve_HOMME

  !=============================================================================

  subroutine FSUNLinSolFree_HOMME()
    !---------------------------------------------------------------------------
    ! Free interface to Fortran column solver
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none

    type(SUNLinearSolver) :: sunls

    ! detach arkode mem pointer
    nullify(arkmem)

    ! free vector pointer
    call FN_VDestroyVectorArray(ycur_ptr, 1)

    ! deallocate overall SUNLinearSolver structure
    call FSUNLinSolFreeEmpty(sunls)

  end subroutine FSUNLinSolFree_HOMME

  !=============================================================================

  integer(c_int) function FARKodeLinSysFn(t, y, fy, A, M, jok, jcur, gamma, user_data, tmp1, tmp2, tmp3) &
       result(ier) bind(C)
    !---------------------------------------------------------------------------
    ! dummy linear system function
    !---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    implicit none
    real(c_double)  :: t
    type(N_Vector)  :: y
    type(N_Vector)  :: fy
    type(SUNMatrix) :: A
    type(SUNMatrix) :: M
    integer(c_int)  :: jok
    integer(c_int)  :: jcur
    real(c_double)  :: gamma
    type(c_ptr)     :: user_data
    type(N_Vector)  :: tmp1
    type(N_Vector)  :: tmp2
    type(N_Vector)  :: tmp3

    ier = 0

    return

  end function FARKodeLinSysFn

  !=============================================================================

end module HommeSUNLinSol
