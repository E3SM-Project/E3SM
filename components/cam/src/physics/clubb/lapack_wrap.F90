!-----------------------------------------------------------------------
! $Id: lapack_wrap.F90 6849 2014-04-22 21:52:30Z charlass@uwm.edu $
!===============================================================================
module lapack_wrap

! Description:
!   Wrappers for the band diagonal and tridiagonal direct matrix
!   solvers contained in the LAPACK library.

! References:
!   LAPACK--Linear Algebra PACKage
!   URL: <http://www.netlib.org/lapack/>
!-----------------------------------------------------------------------
  use constants_clubb, only:  & 
    fstderr ! Variable(s)

  use error_code, only:  & 
    clubb_singular_matrix,  & ! Variable(s)
    clubb_bad_lapack_arg, & 
    clubb_var_equals_NaN, &
    clubb_no_error

  use clubb_precision, only: &
    core_rknd, & ! Variable(s)
    dp

  implicit none

  ! Simple routines
  public :: tridag_solve, band_solve

  ! Expert routines
  public :: tridag_solvex, band_solvex

  private :: lapack_isnan

  ! A best guess for what the precision of a single precision and double
  ! precision float is in LAPACK.  Hopefully this will work more portably on
  ! architectures like Itanium than the old code -dschanen 11 Aug 2011
  integer, parameter, private :: &
    sp = kind ( 0.0 )

  private ! Set Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine tridag_solvex( solve_type, ndim, nrhs, &
                            supd, diag, subd, rhs, &
                            solution, rcond, err_code )

! Description:
!   Solves a tridiagonal system of equations (expert routine).

! References:
!   <http://www.netlib.org/lapack/single/sgtsvx.f>
!   <http://www.netlib.org/lapack/double/dgtsvx.f>

! Notes:
!   More expensive than the simple routine, but tridiagonal
!   decomposition is still relatively cheap.
!-----------------------------------------------------------------------
    use error_code, only: &
      clubb_at_least_debug_level ! Logical function

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    external ::  & 
      sgtsvx,  & ! Single-prec. General Tridiagonal Solver eXpert
      dgtsvx     ! Double-prec. General Tridiagonal Solver eXpert

    intrinsic :: kind

    ! Input variables
    character(len=*), intent(in) ::  & 
      solve_type ! Used to write a message if this fails

    integer, intent(in) ::  & 
      ndim,  & ! N-dimension of matrix
      nrhs     ! # of right hand sides to back subst. after LU-decomp.

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout), dimension(ndim) ::  & 
      diag,       & ! Main diagonal
      subd, supd    ! Sub and super diagonal

    real( kind = core_rknd ), intent(inout), dimension(ndim,nrhs) ::  & 
      rhs ! RHS input

    ! The estimate of the reciprocal of the condition number on the LHS matrix.
    ! If rcond is < machine precision the matrix is singular to working
    ! precision, and info == ndim+1.  If rcond == 0, then the LHS matrix
    ! is singular.  This condition is indicated by a return code of info > 0.
    real( kind = core_rknd ), intent(out) :: rcond

    integer, intent(out) ::  & 
      err_code ! Used to determine when a decomp. failed

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(ndim,nrhs) ::  & 
      solution ! Solution

    ! Local Variables
    ! These contain the decomposition of the matrix
    real( kind = core_rknd ), dimension(ndim-1) :: dlf, duf
    real( kind = core_rknd ), dimension(ndim)   :: df
    real( kind = core_rknd ), dimension(ndim-2) :: du2

    integer, dimension(ndim) ::  & 
      ipivot  ! Index of pivots done during decomposition

    integer, dimension(ndim) ::  & 
      iwork   ! `scrap' array


    real( kind = core_rknd ), dimension(nrhs) ::  & 
      ferr,  & ! Forward error estimate
      berr     ! Backward error estimate

    real( kind = core_rknd ), dimension(3*ndim) ::  & 
      work  ! `Scrap' array

    integer :: info ! Diagnostic output

    integer :: i  ! Array index

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
!    $                   DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
!    $                   WORK, IWORK, INFO )
!-----------------------------------------------------------------------

    if ( kind( diag(1) ) == dp ) then
      call dgtsvx( "Not Factored", "No Transpose lhs", ndim, nrhs,  & 
                   subd(2:ndim), diag, supd(1:ndim-1),  & 
                   dlf, df, duf, du2, ipivot,  & 
                   rhs, ndim, solution, ndim, rcond, & 
                   ferr, berr, work, iwork, info )

    else if ( kind( diag(1) ) == sp ) then
      call sgtsvx( "Not Factored", "No Transpose lhs", ndim, nrhs,  & 
                   subd(2:ndim), diag, supd(1:ndim-1),  & 
                   dlf, df, duf, du2, ipivot,  & 
                   rhs, ndim, solution, ndim, rcond, & 
                   ferr, berr, work, iwork, info )

    else
      stop "tridag_solvex: Cannot resolve the precision of real datatype"

    end if

    ! Print diagnostics for when ferr is large
    if ( clubb_at_least_debug_level( 2 ) .and. any( ferr > 1.e-3_core_rknd ) ) then

      write(fstderr,*) "Warning, large error est. for: " // trim( solve_type )

      do i = 1, nrhs, 1
        write(fstderr,*) "rhs # ", i, "tridag forward error est. =", ferr(i)
        write(fstderr,*) "rhs # ", i, "tridag backward error est. =", berr(i)
      end do

      write(fstderr,'(2(a20,e15.6))') "rcond est. = ", rcond, & 
        "machine epsilon = ", epsilon( diag(1) )
    end if

    select case( info )
    case( :-1 )
      write(fstderr,*) trim( solve_type )// & 
        "illegal value in argument", -info
      err_code = clubb_bad_lapack_arg

    case( 0 )
      ! Success!
      if ( lapack_isnan( ndim, nrhs, solution ) ) then
        err_code = clubb_var_equals_NaN 
      else
        err_code = clubb_no_error
      end if

    case( 1: )
      if ( info == ndim+1 ) then
        write(fstderr,*) trim( solve_type) // & 
          " Warning: matrix is singular to working precision."
        write(fstderr,'(a,e12.5)')  & 
          "Estimate of the reciprocal of the condition number: ", rcond
        err_code = clubb_no_error
      else
        write(fstderr,*) solve_type// & 
          " singular matrix."
        err_code = clubb_singular_matrix
      end if

    end select

    return
  end subroutine tridag_solvex

!-----------------------------------------------------------------------
  subroutine tridag_solve & 
             ( solve_type, ndim, nrhs, &
               supd, diag, subd, rhs, &
               solution, err_code )

! Description:
!   Solves a tridiagonal system of equations (simple routine)

! References:
!   <http://www.netlib.org/lapack/single/sgtsv.f>
!   <http://www.netlib.org/lapack/double/dgtsv.f>
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)
#ifndef NDEBUG
#ifdef ARCH_MIC_KNL
    use, intrinsic :: ieee_exceptions
    !use, intrinsic :: ieee_arithmetic
#endif
#endif
    implicit none

    ! External
    external ::  & 
      sgtsv,  & ! Single-prec. General Tridiagonal Solver eXpert
      dgtsv     ! Double-prec. General Tridiagonal Solver eXpert

    intrinsic :: kind

    ! Input variables
    character(len=*), intent(in) ::  & 
      solve_type ! Used to write a message if this fails

    integer, intent(in) ::  & 
      ndim,  & ! N-dimension of matrix
      nrhs     ! # of right hand sides to back subst. after LU-decomp.

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout), dimension(ndim) ::  & 
      diag,       & ! Main diagonal
      subd, supd ! Sub and super diagonal

    real( kind = core_rknd ), intent(inout), dimension(ndim,nrhs) ::  & 
      rhs ! RHS input

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(ndim,nrhs) ::  & 
      solution ! Solution


    integer, intent(out) ::  & 
      err_code ! Used to determine when a decomp. failed

    ! Local Variables

    real( kind = dp ), dimension(ndim) :: &
      subd_dp, supd_dp, diag_dp

    real( kind = dp ), dimension(ndim,nrhs) :: &
      rhs_dp

    integer :: info ! Diagnostic output

!-----------------------------------------------------------------------
!       *** The LAPACK Routine ***
!       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!-----------------------------------------------------------------------


    if ( kind( diag(1) ) == dp ) then
#ifndef NDEBUG
#ifdef ARCH_MIC_KNL
      ! when floating-point exceptions are turned on, this call was failing with a div-by-zero on KNL with Intel/MKL. Solution 
      ! was to turn off exceptions only here at this call (and only for machine with ARCH_MIC_KNL defined) (github 1183)
      call ieee_set_halting_mode(IEEE_DIVIDE_BY_ZERO, .false.) ! Turn off stopping on div-by-zero only
      !call ieee_set_halting_mode(ieee_usual, .false.)  ! Turn off stopping on all floating-point exceptions
#endif
#endif
      call dgtsv( ndim, nrhs, subd(2:ndim), diag, supd(1:ndim-1),  &
                  rhs, ndim, info )
#ifndef NDEBUG
#ifdef ARCH_MIC_KNL
      call ieee_set_halting_mode(IEEE_DIVIDE_BY_ZERO, .true.) ! Turn back on stopping on div-by-zero only
      !call ieee_set_halting_mode(ieee_usual, .true.)  ! Turn off stopping on all floating-point exceptions
#endif
#endif

    else if ( kind( diag(1) ) == sp ) then
      call sgtsv( ndim, nrhs, subd(2:ndim), diag, supd(1:ndim-1),  & 
                  rhs, ndim, info )

    else
      !stop "tridag_solve: Cannot resolve the precision of real datatype"
      ! Eric Raut Aug 2013: Force double precision
      subd_dp = real( subd, kind=dp )
      diag_dp = real( diag, kind=dp )
      supd_dp = real( supd, kind=dp )
      rhs_dp = real( rhs, kind=dp )
      call dgtsv( ndim, nrhs, subd_dp(2:ndim), diag_dp, supd_dp(1:ndim-1),  &
                  rhs_dp, ndim, info )
      subd = real( subd_dp, kind=core_rknd )
      diag = real( diag_dp, kind=core_rknd )
      supd = real( supd_dp, kind=core_rknd )
      rhs = real( rhs_dp, kind=core_rknd )
    end if

    select case( info )
    case( :-1 )
      write(fstderr,*) trim( solve_type )// & 
        " illegal value in argument", -info
      err_code = clubb_bad_lapack_arg

      solution = -999._core_rknd

    case( 0 )
      ! Success!
      if ( lapack_isnan( ndim, nrhs, rhs ) ) then
        err_code = clubb_var_equals_NaN 
      else
        err_code = clubb_no_error
      end if

      solution = rhs

    case( 1: )
      write(fstderr,*) trim( solve_type )//" singular matrix."
      err_code = clubb_singular_matrix

      solution = -999._core_rknd

    end select

    return
  end subroutine tridag_solve

!-----------------------------------------------------------------------
  subroutine band_solvex( solve_type, nsup, nsub, ndim, nrhs,  & 
                          lhs, rhs, solution, rcond, err_code )
! Description:
!   Restructure and then solve a band diagonal system, with
!   diagnostic output

! References:
!   <http://www.netlib.org/lapack/single/sgbsvx.f>
!   <http://www.netlib.org/lapack/double/dgbsvx.f>

! Notes:
!   I found that due to the use of sgbcon/dgbcon it is much
!   more expensive to use this on most systems than the simple
!   driver. Use this version only if you don't case about compute time.
!   Also note that this version equilibrates the lhs and does an iterative
!   refinement of the solutions, which results in a slightly different answer
!   than the simple driver does. -dschanen 24 Sep 2008
!-----------------------------------------------------------------------
    use error_code, only: &
      clubb_at_least_debug_level ! Logical function

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    external ::  & 
      sgbsvx,  & ! Single-prec. General Band Solver eXpert
      dgbsvx  ! Double-prec. General Band Solver eXpert

    intrinsic :: eoshift, kind, trim

    ! Input Variables
    character(len=*), intent(in) :: solve_type

    integer, intent(in) :: & 
      nsup,  & ! Number of superdiagonals
      nsub,  & ! Number of subdiagonals
      ndim,  & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs     ! Number of RHS's to back substitute for

    real( kind = core_rknd ), dimension(nsup+nsub+1,ndim), intent(inout) ::  & 
      lhs ! Left hand side
    real( kind = core_rknd ), dimension(ndim,nrhs), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim,nrhs), intent(out) :: &
      solution

    ! The estimate of the reciprocal condition number of matrix
    ! after equilibration (if done).
    real( kind = core_rknd ), intent(out) ::  & 
      rcond

    integer, intent(out) :: err_code ! Valid calculation?

    ! Local Variables

    ! Workspaces
    real( kind = core_rknd ), dimension(3*ndim)  :: work
    integer, dimension(ndim) :: iwork

    real( kind = core_rknd ), dimension(2*nsub+nsup+1,ndim) :: & 
      lulhs ! LU Decomposition of the LHS

    integer, dimension(ndim) ::  & 
      ipivot

    real( kind = core_rknd ), dimension(nrhs) ::  & 
      ferr, berr ! Forward and backward error estimate

    real( kind = core_rknd ), dimension(ndim) ::  & 
      rscale, cscale ! Row and column scale factors for the LHS

    integer ::  & 
      info,   & ! If this doesn't come back as 0, something went wrong
      offset, & ! Loop iterator
      imain,  & ! Main diagonal of the matrix
      i         ! Loop iterator

    character ::  & 
      equed ! Row equilibration status


!-----------------------------------------------------------------------
!       Reorder Matrix to use LAPACK band matrix format (5x6)

!       Shift example:

!       [    *        *     lhs(1,1) lhs(1,2) lhs(1,3) lhs(1,4) ] (2)=>
!       [    *     lhs(2,1) lhs(2,2) lhs(2,3) lhs(2,4) lhs(2,5) ] (1)=>
!       [ lhs(3,1) lhs(3,2) lhs(3,3) lhs(3,4) lhs(3,5) lhs(3,6) ]
! <=(1) [ lhs(4,2) lhs(4,3) lhs(4,4) lhs(4,5) lhs(4,6)    *     ]
! <=(2) [ lhs(5,3) lhs(5,4) lhs(5,5) lhs(5,6)    *        *     ]

!       The '*' indicates unreferenced elements.
!       For additional bands above and below the main diagonal, the
!       shifts to the left or right increases by the distance from the
!       main diagonal of the matrix.
!-----------------------------------------------------------------------

    imain = nsup + 1

    ! For the offset, (+) is left, and (-) is right

    ! Sub diagonals
    do offset = 1, nsub, 1
      lhs(imain+offset, 1:ndim) & 
      = eoshift( lhs(imain+offset, 1:ndim), offset )
    end do

    ! Super diagonals
    do offset = 1, nsup, 1
      lhs(imain-offset, 1:ndim) & 
      = eoshift( lhs(imain-offset, 1:ndim), -offset )
    end do

!-----------------------------------------------------------------------
!     *** The LAPACK Routine ***
!     SUBROUTINE SGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,
!    $                   LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX,
!    $                   RCOND, FERR, BERR, WORK, IWORK, INFO )
!-----------------------------------------------------------------------

    if ( kind( lhs(1,1) ) == dp ) then
      call dgbsvx( 'Equilibrate lhs', 'No Transpose lhs', & 
                   ndim, nsub, nsup, nrhs, & 
                   lhs, nsup+nsub+1, lulhs, 2*nsub+nsup+1,  & 
                   ipivot, equed, rscale, cscale, & 
                   rhs, ndim, solution, ndim, & 
                   rcond, ferr, berr, work, iwork, info )

    else if ( kind( lhs(1,1) ) == sp ) then
      call sgbsvx( 'Equilibrate lhs', 'No Transpose lhs', & 
                   ndim, nsub, nsup, nrhs, & 
                   lhs, nsup+nsub+1, lulhs, 2*nsub+nsup+1, & 
                   ipivot, equed, rscale, cscale, & 
                   rhs, ndim, solution, ndim, & 
                   rcond, ferr, berr, work, iwork, info )

    else
      stop "band_solvex: Cannot resolve the precision of real datatype"
      ! One implication of this is that CLUBB cannot be used with quad
      ! precision variables without a quad precision band diagonal solver
    end if

! %% debug
!       select case ( equed )
!       case ('N')
!         print *, "No equilib. was required for lhs."
!       case ('R')
!         print *, "Row equilib. was done on lhs."
!       case ('C')
!         print *, "Column equilib. was done on lhs."
!       case ('B')
!         print *, "Row and column equilib. was done on lhs."
!       end select

!       write(*,'(a,e12.5)') "Row scale : ", rscale
!       write(*,'(a,e12.5)') "Column scale: ", cscale
!       write(*,'(a,e12.5)') "Estimate of the reciprocal of the "//
!                            "condition number: ", rcond
!       write(*,'(a,e12.5)') "Forward Error Estimate: ", ferr
!       write(*,'(a,e12.5)') "Backward Error Estimate: ", berr
! %% end debug

    ! Diagnostic information
    if ( clubb_at_least_debug_level( 2 ) .and. any( ferr > 1.e-3_core_rknd ) ) then

      write(fstderr,*) "Warning, large error est. for: " // trim( solve_type )

      do i = 1, nrhs, 1
        write(fstderr,*) "rhs # ", i, "band_solvex forward error est. =", ferr(i)
        write(fstderr,*) "rhs # ", i, "band_solvex backward error est. =", berr(i)
      end do

      write(fstderr,'(2(a20,e15.6))') "rcond est. = ", rcond, & 
        "machine epsilon = ", epsilon( lhs(1,1) )
    end if

    select case( info )

    case( :-1 )
      write(fstderr,*) trim( solve_type )// & 
        " illegal value for argument", -info
      err_code = clubb_bad_lapack_arg

    case( 0 )
      ! Success!
      if ( lapack_isnan( ndim, nrhs, solution ) ) then
        err_code = clubb_var_equals_NaN 
      else
        err_code = clubb_no_error
      end if

    case( 1: )
      if ( info == ndim+1 ) then
        write(fstderr,*) trim( solve_type )// & 
          " Warning: matrix singular to working precision."
        write(fstderr,'(a,e12.5)')  & 
          "Estimate of the reciprocal of the"// & 
          " condition number: ", rcond
        err_code = clubb_no_error
      else
        write(fstderr,*) trim( solve_type )// & 
          " band solver: singular matrix"
        err_code = clubb_singular_matrix
      end if

    end select

    return
  end subroutine band_solvex

!-----------------------------------------------------------------------
  subroutine band_solve( solve_type, nsup, nsub, ndim, nrhs,  & 
                          lhs, rhs, solution, err_code )
! Description:
!   Restructure and then solve a band diagonal system

! References:
!   <http://www.netlib.org/lapack/single/sgbsv.f>
!   <http://www.netlib.org/lapack/double/dgbsv.f>
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    external ::  & 
      sgbsv,  & ! Single-prec. General Band Solver
      dgbsv  ! Double-prec. General Band Solver

    intrinsic :: eoshift, kind, trim

    ! Input Variables
    character(len=*), intent(in) :: solve_type

    integer, intent(in) :: & 
      nsup,  & ! Number of superdiagonals
      nsub,  & ! Number of subdiagonals
      ndim,  & ! The order of the LHS Matrix, i.e. the # of linear equations
      nrhs  ! Number of RHS's to solve for

    ! Note: matrix lhs is intent(in), not intent(inout)
    ! as in the subroutine band_solvex( )
    real( kind = core_rknd ), dimension(nsup+nsub+1,ndim), intent(in) ::  & 
      lhs ! Left hand side
    real( kind = core_rknd ), dimension(ndim,nrhs), intent(inout) ::  & 
      rhs ! Right hand side(s)

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim,nrhs), intent(out) :: solution

    integer, intent(out) :: err_code ! Valid calculation?

    ! Local Variables

    ! Workspaces
    real( kind = core_rknd ), dimension(2*nsub+nsup+1,ndim) :: & 
      lulhs ! LU Decomposition of the LHS

    real( kind = dp ), dimension(2*nsub+nsup+1,ndim) :: &
      lulhs_dp

    real( kind = dp ), dimension(ndim,nrhs) :: &
      rhs_dp

    integer, dimension(ndim) ::  & 
      ipivot

    integer ::  & 
      info,   & ! If this doesn't come back as 0, something went wrong
      offset, & ! Loop iterator
      imain  ! Main diagonal of the matrix

    ! Copy LHS into Decomposition scratch space
    lulhs = 0.0_core_rknd
    lulhs(nsub+1:2*nsub+nsup+1, 1:ndim) = lhs(1:nsub+nsup+1, 1:ndim)

!-----------------------------------------------------------------------
!       Reorder LU Matrix to use LAPACK band matrix format

!       Shift example for lulhs matrix (note the extra bands):

!       [    +        +        +        +        +        +     ]
!       [    +        +        +        +        +        +     ]
!       [    *        *     lhs(1,1) lhs(1,2) lhs(1,3) lhs(1,4) ] (2)=>
!       [    *     lhs(2,1) lhs(2,2) lhs(2,3) lhs(2,4) lhs(2,5) ] (1)=>
!       [ lhs(3,1) lhs(3,2) lhs(3,3) lhs(3,4) lhs(3,5) lhs(3,6) ]
! <=(1) [ lhs(4,2) lhs(4,3) lhs(4,4) lhs(4,5) lhs(4,6)    *     ]
! <=(2) [ lhs(5,3) lhs(5,4) lhs(5,5) lhs(5,6)    *        *     ]
!       [    +        +        +        +        +        +     ]
!       [    +        +        +        +        +        +     ]

!       The '*' indicates unreferenced elements.
!       The '+' indicates an element overwritten during decomposition.
!       For additional bands above and below the main diagonal, the
!       shifts to the left or right increases by the distance from the
!       main diagonal of the matrix.
!-----------------------------------------------------------------------

    ! Reorder lulhs, omitting the additional 2*nsub bands
    ! that are used for the LU decomposition of the matrix.

    imain = nsub + nsup + 1

    ! For the offset, (+) is left, and (-) is right

    ! Sub diagonals
    do offset = 1, nsub, 1
      lulhs(imain+offset, 1:ndim) & 
      = eoshift( lulhs(imain+offset, 1:ndim), offset )
    end do

    ! Super diagonals
    do offset = 1, nsup, 1
      lulhs(imain-offset, 1:ndim) & 
      = eoshift( lulhs(imain-offset, 1:ndim), -offset )
    end do

!-----------------------------------------------------------------------
!       *** LAPACK routine ***
!       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
!-----------------------------------------------------------------------

    if ( kind( lhs(1,1) ) == dp ) then
      call dgbsv( ndim, nsub, nsup, nrhs, lulhs, nsub*2+nsup+1,  & 
                  ipivot, rhs, ndim, info )

    else if ( kind( lhs(1,1) ) == sp ) then
      call sgbsv( ndim, nsub, nsup, nrhs, lulhs, nsub*2+nsup+1,  & 
                  ipivot, rhs, ndim, info )

    else
      !stop "band_solve: Cannot resolve the precision of real datatype"
      ! One implication of this is that CLUBB cannot be used with quad
      ! precision variables without a quad precision band diagonal solver
      ! Eric Raut Aug 2013: force double precision
      lulhs_dp = real( lulhs, kind=dp )
      rhs_dp = real( rhs, kind=dp )
      call dgbsv( ndim, nsub, nsup, nrhs, lulhs_dp, nsub*2+nsup+1,  &
                  ipivot, rhs_dp, ndim, info )
      rhs = real( rhs_dp, kind=core_rknd )
    end if

    select case( info )

    case( :-1 )
      write(fstderr,*) trim( solve_type )// & 
        " illegal value for argument ", -info
      err_code = clubb_bad_lapack_arg

      solution = -999._core_rknd

    case( 0 )
      ! Success!
      if ( lapack_isnan( ndim, nrhs, rhs ) ) then
        err_code = clubb_var_equals_NaN 
      else
        err_code = clubb_no_error
      end if

      solution = rhs

    case( 1: )
      write(fstderr,*) trim( solve_type )//" band solver: singular matrix"
      err_code = clubb_singular_matrix

      solution = -999._core_rknd

    end select

    return
  end subroutine band_solve

!-----------------------------------------------------------------------
  logical function lapack_isnan( ndim, nrhs, variable )

! Description:
!   Check for NaN values in a variable using the LAPACK subroutines

! References:
!   <http://www.netlib.org/lapack/single/sisnan.f>
!   <http://www.netlib.org/lapack/double/disnan.f>
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none
#ifdef NO_LAPACK_ISNAN /* Used for older LAPACK libraries that don't have sisnan/disnan */

    intrinsic :: any

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    real( kind = core_rknd ), dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

    lapack_isnan = any( variable(:,1:nrhs) /= variable(:,1:nrhs) )
#else
    logical, external :: sisnan, disnan 

    integer, intent(in) :: &
      ndim, & ! Size of variable
      nrhs    ! Number of right hand sides

    real( kind = core_rknd ), dimension(ndim,nrhs), intent(in) :: &
      variable ! Variable to check

    integer :: k, j

    ! ---- Begin Code ----

    lapack_isnan = .false.

    if ( kind( variable ) == dp ) then
      do k = 1, ndim
        do j = 1, nrhs
          lapack_isnan = disnan( variable(k,j) )
          if ( lapack_isnan ) exit
        end do
        if ( lapack_isnan ) exit
      end do
    else if ( kind( variable ) == sp ) then
      do k = 1, ndim
        do j = 1, nrhs
          lapack_isnan = sisnan( variable(k,j) )
          if ( lapack_isnan ) exit
        end do
        if ( lapack_isnan ) exit
      end do
    else
      stop "lapack_isnan: Cannot resolve the precision of real datatype"
    end if
#endif /* NO_LAPACK_ISNAN */

    return
  end function lapack_isnan

end module lapack_wrap
