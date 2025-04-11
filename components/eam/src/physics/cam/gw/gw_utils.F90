module gw_utils

!
! This module contains utility code for the gravity wave modules.
!

implicit none
private
save

! Real kind for gravity wave parameterization.
integer, public, parameter :: r8 = selected_real_kind(12)

! Public interface
public :: get_unit_vector
public :: dot_2d
public :: midpoint_interp

contains

! Take two components of a vector, and find the unit vector components and
! total magnitude.
elemental subroutine get_unit_vector(u, v, u_n, v_n, mag)
  real(r8), intent(in) :: u
  real(r8), intent(in) :: v
  real(r8), intent(out) :: u_n
  real(r8), intent(out) :: v_n
  real(r8), intent(out) :: mag

  mag = sqrt(u*u + v*v)

  if (mag > 0._r8) then
     u_n = u/mag
     v_n = v/mag
  else
     u_n = 0._r8
     v_n = 0._r8
  end if

end subroutine get_unit_vector

! Elemental version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
real(r8) elemental function dot_2d(u1, v1, u2, v2)
  real(r8), intent(in) :: u1, v1
  real(r8), intent(in) :: u2, v2

  dot_2d = u1*u2 + v1*v2

end function dot_2d

! Pure function that interpolates the values of the input array along
! dimension 2. This is obviously not a very generic routine, unlike, say,
! CAM's lininterp. But it's used often enough that it seems worth providing
! here.
pure function midpoint_interp(arr) result(interp)
  real(r8), intent(in) :: arr(:,:)
  real(r8) :: interp(size(arr,1),size(arr,2)-1)

  integer :: i

  do i = 1, size(interp,2)
     interp(:,i) = 0.5_r8 * (arr(:,i)+arr(:,i+1))
  end do

end function midpoint_interp

end module gw_utils
