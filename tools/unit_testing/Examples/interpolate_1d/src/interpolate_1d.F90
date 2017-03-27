module interpolate_1d

implicit none

integer, parameter :: r8 = selected_real_kind(12)

contains

! Using values of an input function (xs, ys), and x at a particular point,
! interpolate to find y at that point.
! Out of bounds values get the value at the boundary.

! Precondition: xs must have non-zero length and be monotonically
! increasing.
! Precondition: ys must be at least as long as xs.
function interpolate_point(xs, ys, x_pt) result(y_pt)
  real(r8), intent(in) :: xs(:), ys(:)
  real(r8), intent(in) :: x_pt
  real(r8) :: y_pt

  integer :: i
  real(r8) :: weight

  if (x_pt < xs(1)) then
     y_pt = ys(1)
     return
  end if

  do i = 2, size(xs)
     if (x_pt < xs(i)) then
        weight = (x_pt-xs(i-1)) / (xs(i)-xs(i-1))
        y_pt = ys(i-1) + (ys(i)-ys(i-1)) * weight
        return
     end if
  end do

  y_pt = ys(size(xs))

end function interpolate_point

end module interpolate_1d
