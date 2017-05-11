program test_driver

use interpolate_1d

implicit none

! Tolerance for rounding error.
real(r8), parameter :: tol = 1.e-15_r8

real(r8), parameter :: test_xs(3) = [0._r8, 1._r8, 2._r8]
real(r8), parameter :: test_ys(3) = [0._r8, 1._r8, 4._r8]

real(r8) :: y

! Check that interpolating at x = 0.5 gives y = ~0.5.
y = interpolate_point(test_xs, test_ys, 0.5_r8)
call assert(abs(y - 0.5_r8) <= tol, &
     "Interpolation within first interval is wrong.")

! Same with 1.5 giving ~2.5.
y = interpolate_point(test_xs, test_ys, 1.5_r8)
call assert(abs(y - 2.5_r8) <= tol, &
     "Interpolation within second interval is wrong.")

! Test output at lower boundary.
y = interpolate_point(test_xs, test_ys, -1._r8)
call assert(abs(y - 0._r8) <= tol, &
     "Interpolate_point extrapolates incorrectly below lower bound.")

! Upper boundary.
y = interpolate_point(test_xs, test_ys, 3._r8)
call assert(abs(y - 4._r8) <= tol, &
     "Interpolate_point extrapolates incorrectly above upper bound.")

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
