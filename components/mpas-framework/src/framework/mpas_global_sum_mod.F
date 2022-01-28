!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
! mpas_global_sum_mod
!
!! \file
!! \brief This file contains reproducible sum functions for common data
!!        data types (R8, R4, I8, I4) and array dimensions (scalar
!!        up to 6-d arrays).  A second interface allows for the sum of 
!!        a product of two arrays that is useful for either dot products
!!        or multiplicative masks. Finally, a multi-field option permits
!!        computing sums for multiple fields (of the same shape, type)
!!        at the same time to reduce the number of messages needed. For
!!        multi-field sums, only dimensions up to 5 are supported.
!!
!!        All sums should give identical answers independent of domain
!!        decomposition or ordering of elements.
!!
!!        The interfaces are overloaded, so to compute the sum of an
!!        array or scalar, the call should be as follows with mpiComm
!!        being the MPI communicator defined for the model or sub-model:
!!
!!        For sums of scalars across MPI tasks, use:
!!           sum = global_sum(scalar, mpiComm)
!!        For sums of arrays over the full range of array indices, use:
!!           sum = global_sum(array, mpiComm)
!!        For sums of arrays over a subset of array indices, use:
!!           sum = global_sum(array, mpiComm, indexRange)
!!            where indexRange is an integer array
!!            of size (2*numDimensions) that contain the start, end
!!            index of each dimension. For example, an MPAS
!!            array(nz, nCells) might have an indexRange defined as
!!                indexRange(1) = kmin, indexRange(2) = kmax,
!!                indexRange(3) = 1,    indexRange(4) = nCellsOwned
!!
!!        To compute a sum over products (eg for a dot product or for
!!        using a multiplicative mask), the interface is the same with 
!!        an extra argument for the second array/scalar, eg:
!!           sum = global_sum(scalar1, scalar2, mpiComm)
!!           sum = global_sum(array1, array2, mpiComm)
!!           sum = global_sum(array1, array2, mpiComm, indexRange)
!!
!!        To compute sums of multiple fields, the interfaces are the
!!        same as above with an added nFld in the name, so:
!!           sum(1:nFlds) = global_sum_nfld(scalar(:), mpiComm)
!!           sum(1:nFlds) = global_sum_nfld(array, mpiComm)
!!           sum(1:nFlds) = global_sum_nfld(array, mpiComm, indexRange)
!!           sum(1:nFlds) = global_sum_nfld(scalar1(:), scalar2(:), mpiComm)
!!           sum(1:nFlds) = global_sum_nfld(array1, array2, mpiComm)
!!           sum(1:nFlds) = global_sum_nfld(array1, array2, mpiComm, indexRange)
!!        The multiple fields must be in an array with the field index
!!        as the last (right-most) index in the array. If an indexRange
!!        is provided, the same range is used for all fields. For the
!!        sum of a product - the two product arrays must have the same
!!        shape and size (including number of fields)
!!
!!        For reproducibility, single precision sums are computed in
!!        double precision. For double precision sums, we use the
!!        double-double algorithm of Donald Knuth, further improved by
!!        David A Bailey and presented in:
!!        He, Yun and Chris Ding, 2001: Using Accurate Arithmetics to
!!          Improve Numerical Reproducibility and Stability in Parallel
!!          Applications, J. of Supercomputing, 18, 259-277.
!
!***********************************************************************

module mpas_global_sum_mod

   implicit none
   private

   include "mpif.h"

   !--------------------------------------------------------------------
   ! Public functions
   !--------------------------------------------------------------------

   public :: mpas_global_sum, mpas_global_sum_nfld

   interface mpas_global_sum
      module procedure globalSumR8d0
      module procedure globalSumR4d0
      module procedure globalSumI8d0
      module procedure globalSumI4d0
      module procedure globalSumR8d1
      module procedure globalSumR4d1
      module procedure globalSumI8d1
      module procedure globalSumI4d1
      module procedure globalSumR8d2
      module procedure globalSumR4d2
      module procedure globalSumI8d2
      module procedure globalSumI4d2
      module procedure globalSumR8d3
      module procedure globalSumR4d3
      module procedure globalSumI8d3
      module procedure globalSumI4d3
      module procedure globalSumR8d4
      module procedure globalSumR4d4
      module procedure globalSumI8d4
      module procedure globalSumI4d4
      module procedure globalSumR8d5
      module procedure globalSumR4d5
      module procedure globalSumI8d5
      module procedure globalSumI4d5
      module procedure globalSumR8d6
      module procedure globalSumR4d6
      module procedure globalSumI8d6
      module procedure globalSumI4d6
      module procedure globalSumProdR8d1
      module procedure globalSumProdR4d1
      module procedure globalSumProdI8d1
      module procedure globalSumProdI4d1
      module procedure globalSumProdR8d2
      module procedure globalSumProdR4d2
      module procedure globalSumProdI8d2
      module procedure globalSumProdI4d2
      module procedure globalSumProdR8d3
      module procedure globalSumProdR4d3
      module procedure globalSumProdI8d3
      module procedure globalSumProdI4d3
      module procedure globalSumProdR8d4
      module procedure globalSumProdR4d4
      module procedure globalSumProdI8d4
      module procedure globalSumProdI4d4
      module procedure globalSumProdR8d5
      module procedure globalSumProdR4d5
      module procedure globalSumProdI8d5
      module procedure globalSumProdI4d5
      module procedure globalSumProdR8d6
      module procedure globalSumProdR4d6
      module procedure globalSumProdI8d6
      module procedure globalSumProdI4d6
   end interface mpas_global_sum

   interface mpas_global_sum_nfld
      module procedure globalSumNfldR8d0
      module procedure globalSumNfldR4d0
      module procedure globalSumNfldI8d0
      module procedure globalSumNfldI4d0
      module procedure globalSumNfldR8d1
      module procedure globalSumNfldR4d1
      module procedure globalSumNfldI8d1
      module procedure globalSumNfldI4d1
      module procedure globalSumNfldR8d2
      module procedure globalSumNfldR4d2
      module procedure globalSumNfldI8d2
      module procedure globalSumNfldI4d2
      module procedure globalSumNfldR8d3
      module procedure globalSumNfldR4d3
      module procedure globalSumNfldI8d3
      module procedure globalSumNfldI4d3
      module procedure globalSumNfldR8d4
      module procedure globalSumNfldR4d4
      module procedure globalSumNfldI8d4
      module procedure globalSumNfldI4d4
      module procedure globalSumNfldR8d5
      module procedure globalSumNfldR4d5
      module procedure globalSumNfldI8d5
      module procedure globalSumNfldI4d5
      module procedure globalSumProdNfldR8d1
      module procedure globalSumProdNfldR4d1
      module procedure globalSumProdNfldI8d1
      module procedure globalSumProdNfldI4d1
      module procedure globalSumProdNfldR8d2
      module procedure globalSumProdNfldR4d2
      module procedure globalSumProdNfldI8d2
      module procedure globalSumProdNfldI4d2
      module procedure globalSumProdNfldR8d3
      module procedure globalSumProdNfldR4d3
      module procedure globalSumProdNfldI8d3
      module procedure globalSumProdNfldI4d3
      module procedure globalSumProdNfldR8d4
      module procedure globalSumProdNfldR4d4
      module procedure globalSumProdNfldI8d4
      module procedure globalSumProdNfldI4d4
      module procedure globalSumProdNfldR8d5
      module procedure globalSumProdNfldR4d5
      module procedure globalSumProdNfldI8d5
      module procedure globalSumProdNfldI4d5
   end interface mpas_global_sum_nfld

   !--------------------------------------------------------------------
   ! Private parameters
   !--------------------------------------------------------------------

   logical :: r8sumInitialized = .false.

   integer :: MPI_SUMDD ! special MPI operator for reproducible r8 sum

   integer, parameter :: & ! supported data types
      I4 = selected_int_kind(6),  &
      I8 = selected_int_kind(13), &
      R4 = selected_real_kind(6), &
      R8 = selected_real_kind(12)

!***********************************************************************

contains

!***********************************************************************

   function globalSumR8d0(scalar, mpiComm) result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), intent(in) :: scalar ! value to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Local variables
      integer :: ierr ! local error flag
      ! complex numbers are used to efficiently store values and
      ! residuals in the reproducible sum
      complex(R8) :: localTmp, globalTmp

      !---- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! Reproducible sums
      localTmp = cmplx(scalar,0.d0,R8)
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localTmp, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d0

!***********************************************************************

   function globalSumR4d0(scalar, mpiComm) result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), intent(in) :: scalar ! value to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Local variables
      integer :: ierr ! local error flag
      real(R8) :: &
         localSum, globalTmp ! accumulate in double for reproducibility

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      localSum = scalar
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision

   end function globalSumR4d0

!***********************************************************************

   function globalSumI8d0(scalar, mpiComm) result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), intent(in) :: scalar ! value to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Local variables
      integer :: ierr ! local error flag

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      call MPI_ALLREDUCE(scalar, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d0

!***********************************************************************

   function globalSumI4d0(scalar, mpiComm) result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), intent(in) :: scalar ! value to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Local variables
      integer :: ierr ! local error flag

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      call MPI_ALLREDUCE(scalar, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d0

!***********************************************************************

   function globalSumR8d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do i=imin,imax
         t1 = array(i) + real(localSum)
         e  = t1 - array(i)
         t2 = ((real(localSum) - e) &
               + (array(i) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d1

!***********************************************************************

   function globalSumR4d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do i=imin,imax
         localSum = localSum + array(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d1

!***********************************************************************

   function globalSumI8d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do i=imin,imax
         localSum = localSum + array(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d1

!***********************************************************************

   function globalSumI4d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum
      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do i=imin,imax
         localSum = localSum + array(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d1

!***********************************************************************

   function globalSumR8d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do j=jmin,jmax
      do i=imin,imax
         t1 = array(i,j) + real(localSum)
         e  = t1 - array(i,j)
         t2 = ((real(localSum) - e) &
               + (array(i,j) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d2

!***********************************************************************

   function globalSumR4d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d2

!***********************************************************************

   function globalSumI8d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d2

!***********************************************************************

   function globalSumI4d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d2

!***********************************************************************

   function globalSumR8d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         t1 = array(i,j,k) + real(localSum)
         e  = t1 - array(i,j,k)
         t2 = ((real(localSum) - e) &
               + (array(i,j,k) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d3

!***********************************************************************

   function globalSumR4d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d3

!***********************************************************************

   function globalSumI8d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d3

!***********************************************************************

   function globalSumI4d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d3

!***********************************************************************

   function globalSumR8d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         t1 = array(i,j,k,l) + real(localSum)
         e  = t1 - array(i,j,k,l)
         t2 = ((real(localSum) - e) &
               + (array(i,j,k,l) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d4

!***********************************************************************

   function globalSumR4d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d4

!***********************************************************************

   function globalSumI8d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d4

!***********************************************************************

   function globalSumI4d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d4

!***********************************************************************

   function globalSumR8d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         t1 = array(i,j,k,l,m) + real(localSum)
         e  = t1 - array(i,j,k,l,m)
         t2 = ((real(localSum) - e) &
               + (array(i,j,k,l,m) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d5

!***********************************************************************

   function globalSumR4d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d5

!***********************************************************************

   function globalSumI8d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d5

!***********************************************************************

   function globalSumI4d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d5

!***********************************************************************

   function globalSumR8d6(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         t1 = array(i,j,k,l,m,n) + real(localSum)
         e  = t1 - array(i,j,k,l,m,n)
         t2 = ((real(localSum) - e) &
               + (array(i,j,k,l,m,n) - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumR8d6

!***********************************************************************

   function globalSumR4d6(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumR4d6

!***********************************************************************

   function globalSumI8d6(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI8d6

!***********************************************************************

   function globalSumI4d6(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:,:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumI4d6

!***********************************************************************

   function globalSumProdR8d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do i=imin,imax
         prod = array(i)*array2(i)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d1

!***********************************************************************

   function globalSumProdR4d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do i=imin,imax
         localSum = localSum + array(i)*array2(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumProdR4d1

!***********************************************************************

   function globalSumProdI8d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do i=imin,imax
         localSum = localSum + array(i)*array2(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d1

!***********************************************************************

   function globalSumProdI4d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do i=imin,imax
         localSum = localSum + array(i)*array2(i)
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d1

!***********************************************************************

   function globalSumProdR8d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do j=jmin,jmax
      do i=imin,imax
         prod = array(i,j)*array2(i,j)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d2

!***********************************************************************

   function globalSumProdR4d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)*array2(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumProdR4d2

!***********************************************************************

   function globalSumProdI8d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)*array2(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d2

!***********************************************************************

   function globalSumProdI4d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j)*array2(i,j)
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d2

!***********************************************************************

   function globalSumProdR8d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         prod = array(i,j,k)*array2(i,j,k)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d3

!***********************************************************************

   function globalSumProdR4d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)*array2(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumProdR4d3

!***********************************************************************

   function globalSumProdI8d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)*array2(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d3

!***********************************************************************

   function globalSumProdI4d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k)*array2(i,j,k)
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d3

!***********************************************************************

   function globalSumProdR8d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         prod = array(i,j,k,l)*array2(i,j,k,l)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d4

!***********************************************************************

   function globalSumProdR4d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)*array2(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   end function globalSumProdR4d4

!***********************************************************************

   function globalSumProdI8d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)*array2(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d4

!***********************************************************************

   function globalSumProdI4d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l)*array2(i,j,k,l)
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d4

!***********************************************************************

   function globalSumProdR8d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         prod = array(i,j,k,l,m)*array2(i,j,k,l,m)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d5

!***********************************************************************

   function globalSumProdR4d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)*array2(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   endfunction globalSumProdR4d5

!***********************************************************************

   function globalSumProdI8d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)*array2(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d5

!***********************************************************************

   function globalSumProdI4d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m)*array2(i,j,k,l,m)
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d5

!***********************************************************************

   function globalSumProdR8d6(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R8) :: globalSum ! final sum as result

      !*** Input variables
      real(R8), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8) :: localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      localSum = cmplx(0.d0,0.d0)
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         prod = array(i,j,k,l,m,n)*array2(i,j,k,l,m,n)
         t1 = prod + real(localSum)
         e  = t1 - prod
         t2 = ((real(localSum) - e) + (prod - (t1 - e))) &
                 + aimag(localSum)
         ! The result is t1 + t2, after normalization.
         localSum = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum = real(globalTmp)

   end function globalSumProdR8d6

!***********************************************************************

   function globalSumProdR4d6(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      real(R4) :: globalSum ! final sum as result

      !*** Input variables
      real(R4), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum
      real(R8) :: &
         localSum,  &! sum of local domain accumulated in double prec
         globalTmp   ! global double precision sum

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0.d0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)*array2(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = globalTmp ! return to single precision for result

   endfunction globalSumProdR4d6

!***********************************************************************

   function globalSumProdI8d6(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(8) :: globalSum ! final sum as result

      !*** Input variables
      integer(8), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum
      integer(8) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)*array2(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI8d6

!***********************************************************************

   function globalSumProdI4d6(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Output (return) variable
      integer(4) :: globalSum ! final sum as result

      !*** Input variables
      integer(4), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         m, mmin, mmax,&! iterator and index range for sum
         n, nmin, nmax  ! iterator and index range for sum

      integer(4) :: localSum  ! sum of local domain

      !--- Begin code

      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
         nmin = indxRange(11)
         nmax = indxRange(12)
      else ! full array size
         imin = 1
         imax = size(array, dim=1)
         jmin = 1
         jmax = size(array, dim=2)
         kmin = 1
         kmax = size(array, dim=3)
         lmin = 1
         lmax = size(array, dim=4)
         mmin = 1
         mmax = size(array, dim=5)
         nmin = 1
         nmax = size(array, dim=6)
      endif
 
      ! Compute the local sum of the input array
      localSum = 0
      do n=nmin,nmax
      do m=mmin,mmax
      do l=lmin,lmax
      do k=kmin,kmax
      do j=jmin,jmax
      do i=imin,imax
         localSum = localSum + array(i,j,k,l,m,n)*array2(i,j,k,l,m,n)
      end do
      end do
      end do
      end do
      end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, 1, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumProdI4d6

!***********************************************************************

   function globalSumNfldR8d0(scalars, mpiComm) result(globalSum)

      !*** Input variables
      real(R8), dimension(:), intent(in) :: &
         scalars ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Output (return) variable
      real(R8), dimension(size(scalars)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer :: nFlds ! number of independent scalars to sum
      integer :: ierr ! local error flag
      ! complex numbers are used to efficiently store values and
      ! residuals in the reproducible sum
      complex(R8), dimension(size(scalars)) :: localTmp, globalTmp

      !---- Begin code
      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFlds = size(scalars)

      ! Reproducible sums
      localTmp(:) = cmplx(scalars(:),0.d0,R8)
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localTmp, globalTmp, nFlds, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d0

!***********************************************************************

   function globalSumNfldR4d0(scalars, mpiComm) result(globalSum)

      !*** Input variables
      real (R4), dimension(:), intent(in) :: &
         scalars ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Output (return) variable
      real (R4), dimension(size(scalars)) :: &
         globalSum ! final sum as result

      !*** Local variables
      real (R8), dimension(size(scalars)) :: &
         local, sumTmp ! perform sums in double precision
      integer :: nFlds ! number of independent scalars to sum
      integer :: ierr ! local error flag

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      nFlds = size(scalars)
      local(:) = scalars(:)
      call MPI_ALLREDUCE(local, sumTmp, nFlds, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, mpiComm, ierr)
      globalSum = sumTmp ! convert double back to real

   end function globalSumNfldR4d0

!***********************************************************************

   function globalSumNfldI4d0(scalars, mpiComm) result(globalSum)

      !*** Input variables
      integer, dimension(:), intent(in) :: &
         scalars ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Output (return) variable
      integer, dimension(size(scalars)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer :: nFlds ! number of independent scalars to sum
      integer :: ierr ! local error flag

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      nFlds = size(scalars)
      call MPI_ALLREDUCE(scalars, globalSum, nFlds, MPI_INTEGER, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d0

!***********************************************************************

   function globalSumNfldI8d0(scalars, mpiComm) result(globalSum)

      !*** Input variables
      integer(I8), dimension(:), intent(in) :: &
         scalars ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI

      !*** Output (return) variable
      integer(I8), dimension(size(scalars)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer :: nFlds ! number of independent scalars to sum
      integer :: ierr ! local error flag

      !---- Begin code
      ! Simply call MPI allreduce to sum across ranks and return to all
      nFlds = size(scalars)
      call MPI_ALLREDUCE(scalars, globalSum, nFlds, MPI_INTEGER8, &
                         MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d0

!***********************************************************************

   function globalSumNfldR8d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=2)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      do ifld=1,nFields
         localSum(ifld) = cmplx(0.d0,0.d0)
         do i=imin,imax
            t1 = array(i,ifld) + real(localSum(ifld))
            e  = t1 - array(i,ifld)
            t2 = ((real(localSum(ifld)) - e) &
                  + (array(i,ifld) - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d1

!***********************************************************************

   function globalSumNfldR4d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=2)) :: &
         localSum, &! sum of local domain, accum in double
         globalTmp  ! final sum in double precision

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      ! Accumulate in double precision for reproducibility
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      ! Accumulate in double precision for reproducibility
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:) ! convert double back to real

   end function globalSumNfldR4d1

!***********************************************************************

   function globalSumNfldI4d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:), intent(in) :: array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=2)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d1

!***********************************************************************

   function globalSumNfldI8d1(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=2)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d1

!***********************************************************************

   function globalSumNfldR8d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=3)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do j=jmin,jmax
         do i=imin,imax
            t1 = array(i,j,ifld) + real(localSum(ifld))
            e  = t1 - array(i,j,ifld)
            t2 = ((real(localSum(ifld)) - e) &
                  + (array(i,j,ifld) - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d2

!***********************************************************************

   function globalSumNfldR4d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=3)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp global sum in double

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumNfldR4d2

!***********************************************************************

   function globalSumNfldI4d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=3)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d2

!***********************************************************************

   function globalSumNfldI8d2(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=3)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d2

!***********************************************************************

   function globalSumNfldR8d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=4)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            t1 = array(i,j,k,ifld) + real(localSum(ifld))
            e  = t1 - array(i,j,k,ifld)
            t2 = ((real(localSum(ifld)) - e) &
                  + (array(i,j,k,ifld) - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d3

!***********************************************************************

   function globalSumNfldR4d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=4)) :: &
         localSum, &! sum of local domain
         globalTmp  ! global sum temp in double

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumNfldR4d3

!***********************************************************************

   function globalSumNfldI4d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=4)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d3

!***********************************************************************

   function globalSumNfldI8d3(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=4)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d3

!***********************************************************************

   function globalSumNfldR8d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=5)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            t1 = array(i,j,k,l,ifld) + real(localSum(ifld))
            e  = t1 - array(i,j,k,l,ifld)
            t2 = ((real(localSum(ifld)) - e) &
                  + (array(i,j,k,l,ifld) - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d4

!***********************************************************************

   function globalSumNfldR4d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=5)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp global sum in double

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumNfldR4d4

!***********************************************************************

   function globalSumNfldI4d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=5)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d4

!***********************************************************************

   function globalSumNfldI8d4(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=5)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d4

!***********************************************************************

   function globalSumNfldR8d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e      ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=6)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            t1 = array(i,j,k,l,m,ifld) + real(localSum(ifld))
            e  = t1 - array(i,j,k,l,m,ifld)
            t2 = ((real(localSum(ifld)) - e) &
                  + (array(i,j,k,l,m,ifld) - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0,R8)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumNfldR8d5

!***********************************************************************

   function globalSumNfldR4d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=6)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp global sum in double

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumNfldR4d5

!***********************************************************************

   function globalSumNfldI4d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=6)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI4d5

!***********************************************************************

   function globalSumNfldI8d5(array, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:,:,:), intent(in) :: &
         array ! values to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=6)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumNfldI8d5

!***********************************************************************

   function globalSumProdNfldR8d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=2)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Accumulate the local sum using Knuth's algorithm
      do ifld=1,nFields
         localSum(ifld) = cmplx(0.d0,0.d0)
         do i=imin,imax
            prod = array(i,ifld)*array2(i,ifld)
            t1 = prod + real(localSum(ifld))
            e  = t1 - prod
            t2 = ((real(localSum(ifld)) - e) + (prod - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, MPI_DOUBLE_COMPLEX, &
                         MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumProdNfldR8d1

!***********************************************************************

   function globalSumProdNfldR4d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=2)) :: &
         localSum, &! sum of local domain, accum in double
         globalTmp  ! final sum in double precision

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      ! Accumulate in double precision for reproducibility
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array (i,ifld)* &
                                              array2(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      ! Accumulate in double precision for reproducibility
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:) ! convert double back to real

   end function globalSumProdNfldR4d1

!***********************************************************************

   function globalSumProdNfldI4d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=2)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array (i,ifld)* &
                                              array2(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI4d1

!***********************************************************************

   function globalSumProdNfldI8d1(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:), intent(in) :: &
         array, array2 ! values to be multiplied, then summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=2)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=2)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=2) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + array (i,ifld)* &
                                              array2(i,ifld)
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI8d1

!***********************************************************************

   function globalSumProdNfldR8d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=3)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do j=jmin,jmax
         do i=imin,imax
            prod = array(i,j,ifld)*array2(i,j,ifld)
            t1 = prod + real(localSum(ifld))
            e  = t1 - prod
            t2 = ((real(localSum(ifld)) - e) + (prod - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumProdNfldR8d2

!***********************************************************************

   function globalSumProdNfldR4d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=3)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp global sum in double

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,ifld)* &
                             array2(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumProdNfldR4d2

!***********************************************************************

   function globalSumProdNfldI4d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      integer, dimension(size(array,dim=3)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,ifld)* &
                             array2(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI4d2

!***********************************************************************

   function globalSumProdNfldI8d2(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=3)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         i, imin, imax,&! iterator and index range for sum
         j, jmin, jmax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=3)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=3) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,ifld)* &
                             array2(i,j,ifld)
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI8d2

!***********************************************************************

   function globalSumProdNfldR8d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=4)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            prod = array(i,j,k,ifld)*array2(i,j,k,ifld)
            t1 = prod + real(localSum(ifld))
            e  = t1 - prod
            t2 = ((real(localSum(ifld)) - e) + (prod - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumProdNfldR8d3

!***********************************************************************

   function globalSumProdNfldR4d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=4)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp for global sum in double

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,ifld)* &
                             array2(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumProdNfldR4d3

!***********************************************************************

   function globalSumProdNfldI4d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=4)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,ifld)* &
                             array2(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI4d3

!***********************************************************************

   function globalSumProdNfldI8d3(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=4)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=4)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=4) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,ifld)* &
                             array2(i,j,k,ifld)
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI8d3

!***********************************************************************

   function globalSumProdNfldR8d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=5)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            prod = array(i,j,k,l,ifld)*array2(i,j,k,l,ifld)
            t1 = prod + real(localSum(ifld))
            e  = t1 - prod
            t2 = ((real(localSum(ifld)) - e) + (prod - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumProdNfldR8d4

!***********************************************************************

   function globalSumProdNfldR4d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=5)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp for global sum in double

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,ifld)* &
                             array2(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumProdNfldR4d4

!***********************************************************************

   function globalSumProdNfldI4d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=5)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,ifld)* &
                             array2(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI4d4

!***********************************************************************

   function globalSumProdNfldI8d4(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=5)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=5)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=5) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,ifld)* &
                             array2(i,j,k,l,ifld)
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI8d4

!***********************************************************************

   function globalSumProdNfldR8d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R8), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R8), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum

      real (R8) :: &
         t1, t2, e, prod ! local temporaries

      ! complex numbers are used to efficiently store values and
      ! residuals for the MPI operator
      complex(R8), dimension(size(array,dim=6)) :: &
         localSum, globalTmp

      !--- Begin code

      ! If this is the first R8 call, initialize the reproducible
      ! sum operator
      if (.not. r8sumInitialized) call globalSumInit()

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         ! Accumulate the local sum using Knuth's algorithm
         localSum(ifld) = cmplx(0.d0,0.d0)
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            prod = array(i,j,k,l,m,ifld)*array2(i,j,k,l,m,ifld)
            t1 = prod + real(localSum(ifld))
            e  = t1 - prod
            t2 = ((real(localSum(ifld)) - e) + (prod - (t1 - e))) &
                    + aimag(localSum(ifld))
            ! The result is t1 + t2, after normalization.
            localSum(ifld) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      globalTmp(:) = cmplx(0.d0,0.d0)
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_COMPLEX, MPI_SUMDD, mpiComm, ierr)
      globalSum(:) = real(globalTmp(:))

   end function globalSumProdNfldR8d5

!***********************************************************************

   function globalSumProdNfldR4d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      real(R4), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      real(R4), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      real(R8), dimension(size(array,dim=6)) :: &
         localSum, &! sum of local domain
         globalTmp  ! temp for global sum in double

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,m,ifld)* &
                             array2(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalTmp, nFields, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, mpiComm, ierr)
      globalSum(:) = globalTmp(:)

   end function globalSumProdNfldR4d5

!***********************************************************************

   function globalSumProdNfldI4d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer, dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer, dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer, dimension(size(array,dim=6)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,m,ifld)* &
                             array2(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI4d5

!***********************************************************************

   function globalSumProdNfldI8d5(array, array2, mpiComm, indxRange) &
      result(globalSum)

      !*** Input variables
      integer(I8), dimension(:,:,:,:,:,:), intent(in) :: &
         array, array2 ! values whose product is to be summed
      integer, intent(in) :: mpiComm ! communicator for MPI
      integer, intent(in), dimension(:), optional :: &
         indxRange ! optional restricted index range for sum

      !*** Output (return) variable
      integer(I8), dimension(size(array,dim=6)) :: &
         globalSum ! final sum as result

      !*** Local variables
      integer ::       &
         ierr,         &! local error flag
         ifld, nFields,&! number of fields to be summed and loop indx
         m, mmin, mmax,&! iterator and index range for sum
         l, lmin, lmax,&! iterator and index range for sum
         k, kmin, kmax,&! iterator and index range for sum
         j, jmin, jmax,&! iterator and index range for sum
         i, imin, imax  ! iterator and index range for sum
      integer(I8), dimension(size(array,dim=6)) :: &
         localSum  ! sum of local domain

      !--- Begin code

      nFields = size(array,dim=6) 
      ! If requested, restrict the range of indices over which to sum
      if (present(indxRange)) then
         imin = indxRange(1)
         imax = indxRange(2)
         jmin = indxRange(3)
         jmax = indxRange(4)
         kmin = indxRange(5)
         kmax = indxRange(6)
         lmin = indxRange(7)
         lmax = indxRange(8)
         mmin = indxRange(9)
         mmax = indxRange(10)
      else ! full array size
         imin = 1
         imax = size(array,dim=1)
         jmin = 1
         jmax = size(array,dim=2)
         kmin = 1
         kmax = size(array,dim=3)
         lmin = 1
         lmax = size(array,dim=4)
         mmin = 1
         mmax = size(array,dim=5)
      endif
 
      ! Compute the local sum of the input array
      do ifld=1,nFields
         localSum(ifld) = 0.d0
         do m=mmin,mmax
         do l=lmin,lmax
         do k=kmin,kmax
         do j=jmin,jmax
         do i=imin,imax
            localSum(ifld) = localSum(ifld) + &
                             array (i,j,k,l,m,ifld)* &
                             array2(i,j,k,l,m,ifld)
         end do
         end do
         end do
         end do
         end do
      end do

      ! Call MPI allreduce to compute the final global sum
      call MPI_ALLREDUCE(localSum, globalSum, nFields, &
                         MPI_INTEGER8, MPI_SUM, mpiComm, ierr)

   end function globalSumProdNfldI8d5

!***********************************************************************
!
!  This subroutine is provided to compute the reproducible sum of two
!  double precision scalars using the double-double algorithm of Knuth
!  referenced in module comments at top. This function is mainly for
!  defining the MPI operator for computing the sum across tasks. For
!  local sums, this algorithm is replicated in all the r8 interfaces
!  above. The routine computes the individual sum dda+ddb for len
!  pairs of values. The operands are defined as complex numbers in the
!  algorithm to track the value and its residual. The result is
!  returned in ddb. The itype argument is required for MPI operators
!  but is ignored here because it is only called for double precision.
!  
!-----------------------------------------------------------------------

   subroutine ddsum (dda, ddb, len, itype)

      !*** Inputs

      integer, intent(in) :: &
         len,        &! the number of individual sums to compute
         itype        ! id for data type (required for MPI operators)

      !*** Input/output

      complex(R8), intent(inout), dimension(len) :: &
         dda, ddb     ! the two numbers to be summed, the result of
                      ! dda+ddb is stored in ddb on return to
                      ! correspond to MPI reduction calls
                      ! The algorithms stores a number and its 
                      ! residual as a complex number

      !*** Local variables

      integer i       ! loop iterator

      real (R8) :: &
         e, t1, t2    ! local temps for residuals and manipulation

      !-------------
      ! Begin code

      do i = 1, len
         ! Compute dda + ddb using Knuth's trick.
         t1 = real(dda(i)) + real(ddb(i))
         e = t1 - real(dda(i))
         t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e))) &
            + aimag(dda(i)) + aimag(ddb(i))
         ! The result is t1 + t2, after normalization.
         ddb(i) = cmplx (t1 + t2, t2 - ((t1 + t2) - t1), R8)
      enddo

   !--------------------------------------------------------------------

   end subroutine ddsum

!***********************************************************************
!
!  This subroutine initialized the MPI operatore for reproducible sums
!  for double precision in MPI.
!  
!-----------------------------------------------------------------------

   subroutine globalSumInit()

      integer :: ierr  ! local error flag

      ! Define the MPI operator for reproducible sums
      call MPI_OP_CREATE(ddsum, .true., MPI_SUMDD, ierr)

      ! Set the initialized flag so this will not be called again
      r8sumInitialized = .true.

   end subroutine globalSumInit

!***********************************************************************

end module mpas_global_sum_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
