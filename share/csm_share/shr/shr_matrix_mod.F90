module shr_matrix_mod

  ! This module contains routines for working with matrices.

  ! If possible, use routines from BLAS / LAPACK. This module should only contain
  ! routines that go beyond what is available in BLAS / LAPACK.

#include "shr_assert.h"
  use shr_kind_mod, only : r8 => SHR_KIND_R8
  use shr_log_mod,  only : errMsg => shr_log_errMsg

  implicit none
  private
  save

  public :: tridiagonal_inverse  ! invert a tridiagonal matrix

contains

  !-----------------------------------------------------------------------
  subroutine tridiagonal_inverse(a, b, c, Tinv)
    !
    ! !DESCRIPTION:
    ! Inverts a tridiagonal matrix.
    !
    ! All input / output arrays should have the same number of elements.
    !
    ! Author: Bill Lipscomb
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: a(:)  ! Center diagonal
    real(r8), intent(in)  :: b(:)  ! Upper diagonal (superdiagonal); b(n) is ignored
    real(r8), intent(in)  :: c(:)  ! Lower diagonal (subdiagonal); c(n) is ignored
    real(r8), intent(out) :: Tinv(:,:)  ! Inverse matrix
    !
    ! !LOCAL VARIABLES:
    integer :: n ! matrix dimension
    integer :: i, j, ii
    real(r8) :: theta(0:size(a))
    real(r8) :: phi(1:(size(a)+1))
    real(r8) :: detT       ! determinant of inverse matrix
    real(r8) :: b_product  ! cumulative product of b coefficients
    real(r8) :: c_product  ! cumulative product of c coefficients

    character(len=*), parameter :: subname = 'tridiagonal_inverse'
    !-----------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Here is the formula for coefficients of the inverse of a tridiagonal matrix:
    !
    !       | a_1   b_1                       |
    !       |                                 |
    !       | c_1   a_2   b_2                 |
    !       |                                 |
    !   T = |       c_2   a_3   ...           |
    !       |                                 |
    !       |             ...   ...     b_n-1 |
    !       |                                 |
    !       |                   c_n-1   a_n   |
    !
    ! Tinv(i,j) = (-1)^(i+j) * (b_i ... b_{j-1}) * theta_{i-1} * phi_{j+1} / theta_n  if i <= j
    !
    ! Tinv(i,j) = (-1)^(i+j) * (c_j ... b_{i-1}) * theta_{j-1} * phi_{i+1} / theta_n  if i > j
    !
    ! where theta_0 = 1
    !       theta_1 = a_1
    !       theta_i = a_i*theta_{i-1} - b_{i-1}*c_{i-1}*theta_{i-2} for i = 2, 3, ..., n
    !
    !       phi_{n+1} = 1
    !       phi_n     = a_n
    !       phi_i     = a_i*phi_{i+1} - b_i*c_i*phi_{i+2} for i = n-1, n-2, ..., 1
    !
    ! Note: For i = j, the b products are evaluated as b_i ... b_{j-1} = 1.
    !------------------------------------------------------------------------

    n = size(a)
    SHR_ASSERT((size(b) == n), errMsg(__FILE__, __LINE__))
    SHR_ASSERT((size(c) == n), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((shape(Tinv) == [n,n]), errMsg(__FILE__, __LINE__))

    ! Compute theta recursively
    ! Note: theta(n) is the determinant

    theta(0) = 1._r8
    theta(1) = a(1)

    do i = 2, n
       theta(i) = a(i)*theta(i-1) - b(i-1)*c(i-1)*theta(i-2)
    enddo

    detT = theta(n)

    ! Compute phi recursively

    phi(n+1) = 1._r8
    phi(n)   = a(n)

    do i = n-1, 1, -1
       phi(i) = a(i)*phi(i+1) - b(i)*c(i)*phi(i+2)
    enddo

    ! Compute coefficients of Tinv

    do j = 1, n
       do i = 1, n

          if (i <= j) then

             ! compute product of b terms from i to j-1

             b_product = 1  ! if i = j
             if (i < j) then
                do ii = i, j-1
                   b_product = b_product*b(ii)
                enddo
             endif

             ! compute coefficient

             Tinv(i,j) = (-1)**(i+j) * b_product * theta(i-1) * phi(j+1) / detT

          else  ! i > j

             ! compute product of c terms from j to i-1

             c_product = 1
             do ii = j, i-1
                c_product = c_product*c(ii)
             enddo

             ! compute coefficient

             Tinv(i,j) = (-1)**(i+j) * c_product * theta(j-1) * phi(i+1) / detT

          endif  ! i >= j

       enddo   ! i
    enddo      ! j


  end subroutine tridiagonal_inverse

end module shr_matrix_mod
