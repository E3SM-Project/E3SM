module circle

integer, parameter :: r8 = selected_real_kind(12)
real(r8), parameter :: pi = 3.14159265358979323846_r8

contains

function circle_area(r)
  real(r8), intent(in) :: r
  real(r8) :: circle_area

  circle_area = pi*r*r

end function circle_area

end module circle
