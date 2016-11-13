program test_driver

use circle, only: circle_area, pi, r8

implicit none

! Roundoff level tolerance.
real(r8), parameter :: tol = 1.e-15_r8

call assert(abs(pi - circle_area(1.0_r8)) <= tol, "Circle has wrong area.")

contains

  subroutine assert(val, msg)
    logical, intent(in) :: val
    character(len=*), intent(in) :: msg

    if (.not. val) then
       print *, msg
       stop 1
    end if

  end subroutine assert

end program test_driver
