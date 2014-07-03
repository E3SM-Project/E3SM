
subroutine bnddyi (ncdate, ncsec, doy)
!----------------------------------------------------------------------- 
! 
! Purpose: Convert date and seconds of day to floating point calendar day, for
!          boundary dataset handling
! 
! Method: Use table of days per month to do conversion
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils,   only: endrun
   use cam_logfile,  only: iulog
!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments------------------------------------
!
! Arguments
!
   integer, intent(in) :: ncdate      ! Current date as yymmdd or yyyymmdd
   integer, intent(in) :: ncsec       ! Seconds of day for current date

   real(r8), intent(out) :: doy       ! Day of year
!
! Local Variables
!
   integer mnth        ! Month number
   integer mday        ! Day number of month
   integer jdcon(12)   ! Starting day number for each month
   save jdcon
   data jdcon/0,31,59,90,120,151,181,212,243,273,304,334/
!
! Decode month and day
!
   mnth = mod(ncdate,10000)/100
   if (mnth < 1 .or. mnth > 12) then
      write(iulog,*)'BNDDYI: Bad month index=', mnth
      call endrun
   end if
   mday = mod(ncdate,100)
   doy = jdcon(mnth) + mday + ncsec/86400._r8

   if (doy < 1._r8 .or. doy > 366._r8) then
      write(iulog,*)'BNDDYI: bad day of year = ',doy
      call endrun
   end if
!
   return
end subroutine bnddyi
