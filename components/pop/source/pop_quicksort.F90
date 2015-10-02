!
! this code borrowed from cam (with some minor modifications)
!
! http://www.ccsm.ucar.edu/models/atm-cam/docs/cam3.0/cam30fv-browser/
!   html_code/utils/quicksort.F90.html


module pop_quicksort


! sort routine to arrange array elements from smallest to largest
!
! just sort the array (no auxiallary array also sorted)
! also pass in the size since we may have extra stuff at the end



implicit none
private
public pop_quick_sort
contains


RECURSIVE SUBROUTINE pop_quick_sort(list, list_size)
  use kinds_mod, only : i4
implicit none

integer(i4), DIMENSION (:), INTENT(INOUT)  :: list
integer(i4), intent(in) :: list_size

CALL quick_sort_1(1, list_size)

CONTAINS


RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j
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
INTEGER             :: i, j
integer                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE pop_quick_sort

end module pop_quicksort
