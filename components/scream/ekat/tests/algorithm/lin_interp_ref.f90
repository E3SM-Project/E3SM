#include "ekat_config.f"

module scream_ut

  use iso_c_binding, only: c_int, c_double, c_float

  implicit none

#ifdef EKAT_DOUBLE_PRECISION
  integer, parameter :: c_real = c_double
#else
  integer, parameter :: c_real = c_float
#endif

contains

  subroutine linear_interp_c(x1,x2,y1,y2,km1,km2,ncol,minthresh) bind(c)

    integer(kind=c_int), value, intent(in) :: km1, km2, ncol
    real(kind=c_real), value, intent(in) :: minthresh
    real(kind=c_real), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(kind=c_real), intent(in) :: x2(ncol,km2)
    real(kind=c_real), intent(out) :: y2(ncol,km2)

    call linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)

  end subroutine linear_interp_c

  !==============================================================
  ! Linear interpolation to get values on various grids

  subroutine linear_interp(x1,x2,y1,y2,km1,km2,ncol,minthresh)
    use iso_c_binding
    implicit none

    integer(kind=c_int), intent(in) :: km1, km2
    integer(kind=c_int), intent(in) :: ncol
    real(kind=c_real), intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real(kind=c_real), intent(in) :: x2(ncol,km2)
    real(kind=c_real), intent(in) :: minthresh
    real(kind=c_real), intent(out) :: y2(ncol,km2)

    integer :: k1, k2, i

    do i=1,ncol
       do k2=1,km2
          if( x2(i,k2) <= x1(i,1) ) then
             y2(i,k2) = y1(i,1) + (y1(i,2)-y1(i,1))*(x2(i,k2)-x1(i,1))/(x1(i,2)-x1(i,1))
          elseif( x2(i,k2) >= x1(i,km1) ) then
             y2(i,k2) = y1(i,km1) + (y1(i,km1)-y1(i,km1-1))*(x2(i,k2)-x1(i,km1))/(x1(i,km1)-x1(i,km1-1))
          else
             do k1 = 2,km1
                if( (x2(i,k2)>=x1(i,k1-1)).and.(x2(i,k2)<x1(i,k1)) ) then
                   y2(i,k2) = y1(i,k1-1) + (y1(i,k1)-y1(i,k1-1))*(x2(i,k2)-x1(i,k1-1))/(x1(i,k1)-x1(i,k1-1))
                endif
             enddo ! end k1 loop
          endif

          if (y2(i,k2) .lt. minthresh) then
             y2(i,k2) = minthresh
          endif

       enddo ! end k2 loop
    enddo ! i loop

    return

end subroutine linear_interp

end module scream_ut
