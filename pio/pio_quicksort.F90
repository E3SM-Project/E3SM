!
! this code borrowed from cam
!
! http://www.ccsm.ucar.edu/models/atm-cam/docs/cam3.0/cam30fv-browser/
!   html_code/utils/quicksort.F90.html


module pio_quicksort


! sort routine to arrange array elements from smallest to largest
!
! grabbed from A millers web site http://users.bigpond.net.au/amiller/
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
! pjr added module declaration
! mvr modified integer array to intent inout - may now be any integer 
!     array that gets sorted along with associated real array
!
! 2007/10/2 rml modified for integer input and initialize the order array
!


implicit none
save
private
public quick_sort
contains


RECURSIVE SUBROUTINE quick_sort(list, order)
  use pio_kinds, only : i4
implicit none

integer(i4), DIMENSION (:), INTENT(INOUT)  :: list
INTEGER(i4), DIMENSION (:), INTENT(OUT)  :: order

! Local variable
INTEGER :: i
integer :: count

  count=size(list)
  do i=1,count
    order(i)=i
  end do

CALL quick_sort_1(1, SIZE(list))

CONTAINS


RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
integer                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1



SUBROUTINE interchange_sort(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
integer                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort

end module pio_quicksort
