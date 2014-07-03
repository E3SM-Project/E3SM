module time_utils

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use cam_logfile,  only : iulog

  private
  public  :: flt_date, moz_findplb

contains

  subroutine moz_findplb( x, nx, xval, index )
    !-----------------------------------------------------------------------
    ! 	... find periodic lower bound
    ! 	search the input array for the lower bound of the interval that
    ! 	contains the input value.  the returned index satifies:
    ! 	x(index) .le. xval .lt. x(index+1)
    ! 	assume the array represents values in one cycle of a periodic coordinate.
    ! 	so, if xval .lt. x(1), then index=nx.
    !-----------------------------------------------------------------------

    implicit none

    !-----------------------------------------------------------------------
    !	... dummy args
    !-----------------------------------------------------------------------
    integer, intent(in)  :: nx
    integer, intent(out) :: index
    real(r8), intent(in) :: x(nx)               ! strictly increasing array
    real(r8), intent(in) :: xval

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: i

    if( xval < x(1) .or. xval >= x(nx) ) then
       index = nx
       return
    end if

    do i = 2,nx
       if( xval < x(i) ) then
          index = i - 1
          exit
       end if
    end do

  end subroutine moz_findplb

  real(r8) function flt_date( ncdate, ncsec )
    !----------------------------------------------------------------------- 
    ! Purpose: Convert date and seconds of day to floating point days since
    !          0001/01/01
    !-----------------------------------------------------------------------
    use time_manager, only : timemgr_datediff
    implicit none

    !-----------------------------------------------------------------------
    !	... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)   :: ncdate      ! Current date as yyyymmdd
    integer, intent(in)   :: ncsec       ! Seconds of day for current date

    integer :: refymd = 00010101
    integer :: reftod = 0

    call timemgr_datediff(refymd, reftod, ncdate, ncsec, flt_date)
  end function flt_date

end module time_utils
