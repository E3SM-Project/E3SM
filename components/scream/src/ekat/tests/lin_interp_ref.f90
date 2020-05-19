#include "scream_config.f"

#ifdef SCREAM_DOUBLE_PRECISION
#define c_real c_double
#else
#define c_real c_float
#endif

module scream_ut

  implicit none

contains

  subroutine linear_interp_c(x1,x2,y1,y2,km1,km2,ncol,minthresh) bind(c)

    use iso_c_binding

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
    implicit none

    integer, intent(in) :: km1, km2
    integer, intent(in) :: ncol
    real, intent(in) :: x1(ncol,km1), y1(ncol,km1)
    real, intent(in) :: x2(ncol,km2)
    real, intent(in) :: minthresh
    real, intent(out) :: y2(ncol,km2)

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
