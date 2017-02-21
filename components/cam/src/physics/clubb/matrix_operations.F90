!-----------------------------------------------------------------------
! $Id: matrix_operations.F90 7655 2015-04-29 19:22:26Z raut@uwm.edu $
!===============================================================================
module matrix_operations

  implicit none


  public :: symm_covar_matrix_2_corr_matrix, Cholesky_factor, &
    row_mult_lower_tri_matrix, print_lower_triangular_matrix, &
    get_lower_triangular_matrix, set_lower_triangular_matrix, &
    mirror_lower_triangular_matrix

  private :: Symm_matrix_eigenvalues

  private ! Default scope

  contains
 
!-----------------------------------------------------------------------
  subroutine symm_covar_matrix_2_corr_matrix( ndim, covar, corr )

! Description:
!   Convert a matrix of covariances in to a matrix of correlations.
!   This only does the computation the lower triangular portion of the 
!   matrix.
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! double precision

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    integer, intent(in) :: ndim

    real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: &
      covar ! Covariance Matrix [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim,ndim), intent(out) :: & 
      corr ! Correlation Matrix [-]

    ! Local Variables
    integer :: i, j

    ! ---- Begin Code ----

    corr = 0._core_rknd ! Initialize to 0

    do i = 1, ndim
      do j = 1, i
        corr(i,j) = covar(i,j) / sqrt( covar(i,i) * covar(j,j) )
      end do
    end do

    return
  end subroutine symm_covar_matrix_2_corr_matrix
!-----------------------------------------------------------------------
  subroutine row_mult_lower_tri_matrix( ndim, xvector, tmatrix_in, tmatrix_out )

! Description:
!   Do a row-wise multiply of the elements of a lower triangular matrix.
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! double precision

    implicit none


    ! Input Variables
    integer, intent(in) :: ndim

    real( kind = core_rknd ), dimension(ndim), intent(in) :: & 
      xvector ! Factors to be multiplied across a row [units vary]

    ! Input Variables
    real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: &
      tmatrix_in ! nxn matrix (usually a correlation matrix) [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim,ndim), intent(inout) :: &
      tmatrix_out ! nxn matrix (usually a covariance matrix) [units vary]

    ! Local Variables
    integer :: i, j

    ! ---- Begin Code ----

    do i = 1, ndim
      do j = 1, i
        tmatrix_out(i,j) = tmatrix_in(i,j) * xvector(i)
      end do
    end do

    return
  end subroutine row_mult_lower_tri_matrix

!-------------------------------------------------------------------------------
  subroutine Cholesky_factor( ndim, a_input, a_scaling, a_Cholesky, l_scaled )
!  Description:
!    Create a Cholesky factorization of a_input.
!    If the factorization fails we use a modified a_input matrix and attempt
!    to factorize again.
!
!  References:
!    <http://www.netlib.org/lapack/explore-html/a00868.html> dpotrf
!    <http://www.netlib.org/lapack/explore-html/a00860.html> dpoequ
!    <http://www.netlib.org/lapack/explore-html/a00753.html> dlaqsy
!-------------------------------------------------------------------------------
    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use constants_clubb, only: &
      fstderr ! Constant

    use clubb_precision, only: & 
      core_rknd

    implicit none

    ! External
    external :: dpotrf, dpoequ, dlaqsy, & ! LAPACK subroutines
                spotrf, spoequ, slaqsy

    ! Constant Parameters
    integer, parameter :: itermax = 10 ! Max iterations of the modified method

    real( kind = core_rknd), parameter :: d_coef = 0.1_core_rknd 
       ! Coefficient applied if the decomposition doesn't work

    ! Input Variables
    integer, intent(in) :: ndim

    real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: a_input

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim), intent(out) :: a_scaling

    real( kind = core_rknd ), dimension(ndim,ndim), intent(out) :: a_Cholesky

    logical, intent(out) :: l_scaled

    ! Local Variables
    real( kind = core_rknd ), dimension(ndim) :: a_eigenvalues
    real( kind = core_rknd ), dimension(ndim,ndim) ::  a_corr, a_scaled

    real( kind = core_rknd ) :: tau, d_smallest

    real( kind = core_rknd ) :: amax, scond
    integer :: info
    integer :: i, j, iter

    character :: equed

    logical :: l_dp

    ! ---- Begin code ----

    a_scaled = a_input ! Copy input array into output array

!   do i = 1, n
!     do j = 1, n
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    equed = 'N'

    if ( kind( 0.0_core_rknd ) == kind( 0.0d0 ) ) then
      l_dp = .true.
    else if ( kind( 0.0_core_rknd ) == kind( 0.0 ) ) then
      l_dp = .false.
    else
      stop "Precision is not single or double precision in Cholesky_factor"
    end if

    ! Compute scaling for a_input
    if ( l_dp ) then
      call dpoequ( ndim, a_input, ndim, a_scaling, scond, amax, info )
    else
      call spoequ( ndim, a_input, ndim, a_scaling, scond, amax, info )
    end if

    if ( info == 0 ) then
      ! Apply scaling to a_input
      if ( l_dp ) then
        call dlaqsy( 'Lower', ndim, a_scaled, ndim, a_scaling, scond, amax, equed )
      else
        call slaqsy( 'Lower', ndim, a_scaled, ndim, a_scaling, scond, amax, equed )
      end if
    end if

    ! Determine if scaling was necessary
    if ( equed == 'Y' ) then
      l_scaled = .true.
      a_Cholesky = a_scaled
    else
      l_scaled = .false.
      a_Cholesky = a_input
    end if

    do iter = 1, itermax

      if ( l_dp ) then
        call dpotrf( 'Lower', ndim, a_Cholesky, ndim, info )
      else
        call spotrf( 'Lower', ndim, a_Cholesky, ndim, info )
      end if

      select case( info )
      case( :-1 )
        write(fstderr,*) "Cholesky_factor " // & 
          " illegal value for argument ", -info
        stop
      case( 0 )
        ! Success!
        if ( clubb_at_least_debug_level( 1 ) .and. iter > 1 ) then
          write(fstderr,*) "a_factored (worked)="
          do i = 1, ndim
            do j = 1, i
              write(fstderr,'(g10.3)',advance='no') a_Cholesky(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if
        exit
      case( 1: )
        if ( clubb_at_least_debug_level( 1 ) ) then
          ! This shouldn't happen now that the s and t Mellor(chi/eta) elements have been
          ! modified to never be perfectly correlated, but it's here just in case.
          ! -dschanen 10 Sept 2010
          write(fstderr,*) "Cholesky_factor: leading minor of order ", &
            info, " is not positive definite."
          write(fstderr,*) "factorization failed."
          write(fstderr,*) "a_input="
          do i = 1, ndim
            do j = 1, i
              write(fstderr,'(g10.3)',advance='no') a_input(i,j)
            end do
            write(fstderr,*) ""
          end do
          write(fstderr,*) "a_Cholesky="
          do i = 1, ndim
            do j = 1, i
              write(fstderr,'(g10.3)',advance='no') a_Cholesky(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if

        if ( clubb_at_least_debug_level( 2 ) ) then
          call Symm_matrix_eigenvalues( ndim, a_input, a_eigenvalues )
          write(fstderr,*) "a_eigenvalues="
          do i = 1, ndim
            write(fstderr,'(g10.3)',advance='no') a_eigenvalues(i)
          end do
          write(fstderr,*) ""

          call symm_covar_matrix_2_corr_matrix( ndim, a_input, a_corr )
          write(fstderr,*) "a_correlations="
          do i = 1, ndim
            do j = 1, i
              write(fstderr,'(g10.3)',advance='no') a_corr(i,j)
            end do
            write(fstderr,*) ""
          end do
        end if

        if ( iter == itermax ) then
          write(fstderr,*) "iteration =", iter, "itermax =", itermax
          stop "Fatal error in Cholesky_factor"
        else if ( clubb_at_least_debug_level( 1 ) ) then
          ! Adding a STOP statement to prevent this problem from slipping under
          ! the rug.
          stop "Fatal error in Cholesky_factor"
          write(fstderr,*) "Attempting to modify matrix to allow factorization."
        end if

        if ( l_scaled ) then
          a_Cholesky = a_scaled
        else
          a_Cholesky = a_input
        end if
        ! The number used for tau here is case specific to the Sigma covariance
        ! matrix in the latin hypercube code and is not at all general.
        ! Tau should be a number that is small relative to the other diagonal
        ! elements of the matrix to have keep the error caused by modifying 'a' low.
        ! -dschanen 30 Aug 2010
        d_smallest = a_Cholesky(1,1)
        do i = 2, ndim
          if ( d_smallest > a_Cholesky(i,i) ) d_smallest = a_Cholesky(i,i)
        end do
        ! Use the smallest element * d_coef * iteration
        tau = d_smallest * d_coef * real( iter, kind=core_rknd ) 

!       print *, "tau =", tau, "d_smallest = ", d_smallest

        do i = 1, ndim
          do j = 1, ndim
            if ( i == j ) then
              a_Cholesky(i,j) = a_Cholesky(i,j) + tau ! Add tau to the diagonal
            else
              a_Cholesky(i,j) = a_Cholesky(i,j)
            end if
          end do
        end do

        if ( clubb_at_least_debug_level( 2 ) ) then
          call Symm_matrix_eigenvalues( ndim, a_Cholesky, a_eigenvalues )
          write(fstderr,*) "a_modified eigenvalues="
          do i = 1, ndim
            write(fstderr,'(e10.3)',advance='no') a_eigenvalues(i)
          end do
          write(fstderr,*) ""
        end if

      end select ! info
    end do ! 1..itermax

    return
  end subroutine Cholesky_factor

!----------------------------------------------------------------------
  subroutine Symm_matrix_eigenvalues( ndim, a_input, a_eigenvalues )

!   Description:
!     Computes the eigevalues of a_input
!
!   References:
!     None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant

    use clubb_precision, only: &
      core_rknd ! double precision

    implicit none

    ! External
    external :: dsyev, ssyev ! LAPACK subroutine(s)

    ! Parameters
    integer, parameter :: &
      lwork = 180 ! This is the optimal value I obtained for an n of 5 -dschanen 31 Aug 2010

    ! Input Variables
    integer, intent(in) :: ndim

    real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: a_input

    ! Output Variables
    real( kind = core_rknd ), dimension(ndim), intent(out) :: a_eigenvalues

    ! Local Variables
    real( kind = core_rknd ), dimension(ndim,ndim) :: a_scratch

    real( kind = core_rknd ), dimension(lwork) :: work

    integer :: info
!   integer :: i, j
    ! ---- Begin code ----

    a_scratch = a_input

!   do i = 1, ndim
!     do j = 1, ndim
!       write(6,'(e10.3)',advance='no') a(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause

    if ( kind( 0.0_core_rknd ) == kind( 0.0d0 ) ) then
      call dsyev( 'No eigenvectors', 'Lower', ndim, a_scratch, ndim, &
                  a_eigenvalues, work, lwork, info )
    else if ( kind( 0.0_core_rknd ) == kind( 0.0 ) ) then
      call ssyev( 'No eigenvectors', 'Lower', ndim, a_scratch, ndim, &
                  a_eigenvalues, work, lwork, info )
    else
      stop "Precision is not single or double in Symm_matrix_eigenvalues"
    end if

    select case( info )
    case( :-1 )
      write(fstderr,*) "Symm_matrix_eigenvalues:" // & 
        " illegal value for argument ", -info
      stop
    case( 0 )
      ! Success!

    case( 1: )
      write(fstderr,*) "Symm_matrix_eigenvalues: Algorithm failed to converge."
      stop
    end select

    return
  end subroutine Symm_matrix_eigenvalues
!-------------------------------------------------------------------------------
  subroutine set_lower_triangular_matrix( d_variables, index1, index2, xpyp, &
                                          matrix )
! Description:
!   Set a value for the lower triangular portion of a matrix.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! user defined precision

    implicit none

    ! External
    intrinsic :: max, min

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of variates
      index1, index2 ! Indices for 2 variates (the order doesn't matter)

    real( kind = core_rknd ), intent(in) :: &
      xpyp ! Value for the matrix (usually a correlation or covariance) [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(inout) :: &
      matrix ! The lower triangular matrix

    integer :: i,j

    ! ---- Begin Code ----

    ! Reverse these to set the values of upper triangular matrix
    i = max( index1, index2 )
    j = min( index1, index2 )

    if( i > 0 .and. j > 0 ) then
      matrix(i,j) = xpyp
    end if

    return
  end subroutine set_lower_triangular_matrix
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine get_lower_triangular_matrix( d_variables, index1, index2, matrix, &
                                          xpyp )
! Description:
!   Returns a value from the lower triangular portion of a matrix.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    implicit none

    ! External
    intrinsic :: max, min

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of variates
      index1, index2 ! Indices for 2 variates (the order doesn't matter)

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      matrix ! The covariance matrix

    real( kind = core_rknd ), intent(out) :: &
      xpyp ! Value from the matrix (usually a correlation or covariance) [units vary]

    integer :: i,j

    ! ---- Begin Code ----

    ! Reverse these to set the values of upper triangular matrix
    i = max( index1, index2 )
    j = min( index1, index2 )

    xpyp = matrix(i,j)

    return
  end subroutine get_lower_triangular_matrix

!-----------------------------------------------------------------------
  subroutine print_lower_triangular_matrix( iunit, ndim, matrix )

! Description:
!   Print the values of lower triangular matrix to a file or console.

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iunit, & ! File I/O logical unit (usually 6 for stdout and 0 for stderr)
      ndim     ! Dimension of the matrix

    real( kind = core_rknd ), dimension(ndim,ndim), intent(in) :: &
      matrix ! Lower triangular matrix [units vary]

    ! Local Variables
    integer :: i, j

    ! ---- Begin Code ----

    do i = 1, ndim
      do j = 1, i
        write(iunit,fmt='(g15.6)',advance='no') matrix(i,j)
      end do
      write(iunit,fmt=*) "" ! newline
    end do

    return
  end subroutine print_lower_triangular_matrix

  !-----------------------------------------------------------------------
  subroutine mirror_lower_triangular_matrix( nvars, matrix )

  ! Description:
  !   Mirrors the elements of a lower triangular matrix to the upper
  !   triangle so that it is symmetric.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd  ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nvars ! Number of variables in each dimension of square matrix

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nvars,nvars), intent(inout) :: &
      matrix  ! Lower triangluar square matrix

    ! Local Variables
    integer :: row, col

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( nvars > 1 ) then

      do col=2, nvars
        do row=1, col-1

          matrix(row,col) = matrix(col,row)

        end do
      end do

    end if ! nvars > 1

    return

  end subroutine mirror_lower_triangular_matrix
  !-----------------------------------------------------------------------

end module matrix_operations
