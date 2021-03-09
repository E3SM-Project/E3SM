subroutine output_dateinfo (ncid, dateid, datesecid, timeid, iyr1out, &
                            mon1out, iyrnout, monnout, monlen)
   use prec

   implicit none

   integer, intent(in) :: ncid        ! netcdf file id
   integer, intent(in) :: dateid      ! date id
   integer, intent(in) :: datesecid   ! id for seconds of current date
   integer, intent(in) :: timeid      ! time variable id
   integer, intent(in) :: iyr1out     ! start year written to output file (default iyr1rd)
   integer, intent(in) :: mon1out     ! start month written to output file (default mon1rd)
   integer, intent(in) :: iyrnout     ! end year written to output file (default iyrnrd)
   integer, intent(in) :: monnout     ! end month written to output file
   integer, intent(in) :: monlen(12)  ! length of months
   
   integer :: date(1)                 ! array to please lf95
   integer :: datesec(1)              ! array to please lf95
   integer :: yr                      ! year
   integer :: mo                      ! month
   integer :: start(1)                ! initial value
   integer, parameter :: count(1) = 1 ! number of values to write to netcdf file

   real(r8) :: time(1)                ! array to please lf95
   real(r8) :: prvrem                 ! remaining length of prv month
   real(r8) :: halfmo                 ! half the length of current month
   
   if (mon1out < 1 .or. mon1out > 12) then
      call err_exit ('output_dateinfo: mon1out must be between 1 and 12')
   end if

   if (monnout < 1 .or. monnout > 12) then
      call err_exit ('output_dateinfo: monnout must be between 1 and 12')
   end if

   yr = iyr1out
   mo = mon1out
   start(1) = 0
   time(1)  = 0.
   prvrem   = 0.
   do while (yr < iyrnout .or. yr == iyrnout .and. mo <= monnout)
      start(1) = start(1) + 1
      date(1) = yr*10000 + mo*100 + monlen(mo)/2 + 1

      datesec(1) = 0
      if (mod(monlen(mo), 2) /= 0) datesec(1) = 43200

      halfmo  = monlen(mo)/2 + datesec(1)/86400._r8       ! half of the current month length
      time(1) = time(1) + prvrem + halfmo                 ! plus remaining half of prv month
      prvrem  = monlen(mo) - halfmo                       ! set remainder for next iteration
      call wrap_nf_put_vara_int (ncid, dateid, start, count, date)
      call wrap_nf_put_vara_int (ncid, datesecid, start, count, datesec)
      call wrap_nf_put_vara_double (ncid, timeid, start, count, time)
      mo = mo + 1
      if (mo == 13) then
         mo = 1
         yr = yr + 1
      end if
   end do

   return
end subroutine output_dateinfo
