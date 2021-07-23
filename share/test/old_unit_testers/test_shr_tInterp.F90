program test_shr_tInterp
use shr_kind_mod
use test_mod
use shr_tInterp_mod
use shr_cal_mod,     only : shr_cal_noleap
use shr_const_mod,   only : SHR_CONST_CDAY

implicit none

integer :: date_lb, date_ub, date_in
integer :: sec_lb, sec_ub, sec_in
real(SHR_KIND_R8) :: f1, f2
character(SHR_KIND_CS) :: alogo
character(SHR_KIND_CS) :: calendar_name = shr_cal_noleap
real(SHR_KIND_R8) :: expected(2), values(2)
integer :: rc
integer, parameter :: LIN_TEST = 1, LOWER_TEST = 2, UPPER_TEST = 3, &
                      NEAREST_TEST = 4, num_tests = 4, num_times = 47
integer :: n, i

call test_init( num_tests*num_times+3 )
do n = 1, num_tests
   if (      n == LIN_TEST     )then
      alogo = 'linear'
   else if ( n == LOWER_TEST   )then
      alogo = 'lower'
   else if ( n == UPPER_TEST   )then
      alogo = 'upper'
   else if ( n == NEAREST_TEST )then
      alogo = 'nearest'
   end if

   write(*,*) "Test type: ", trim(alogo)

   date_lb = 20010101
   date_ub = 20010102
   sec_lb  = 0
   sec_ub  = 0
   date_in = 20010101
   sec_in  = 0
   do i = 1, num_times
      write(*,*) "seconds in ", sec_in
      if (      n == LIN_TEST     )then
         f1 = sec_in / SHR_CONST_CDAY
         expected = (/ 1.0_SHR_KIND_R8 - f1, f1 /)
      else if ( n == LOWER_TEST   )then
         expected = (/ 1.0_SHR_KIND_R8, 0.0_SHR_KIND_R8 /)
      else if ( n == UPPER_TEST   )then
         expected = (/ 0.0_SHR_KIND_R8, 1.0_SHR_KIND_R8 /)
      else if ( n == NEAREST_TEST )then
         if ( sec_in <= SHR_CONST_CDAY /2 )then
            expected = (/ 1.0_SHR_KIND_R8, 0.0_SHR_KIND_R8 /)
         else
            expected = (/ 0.0_SHR_KIND_R8, 1.0_SHR_KIND_R8 /)
         end if
      end if
      call shr_tInterp_getFactors( date_lb, sec_lb, date_ub, sec_ub, date_in, &
                                   sec_in, f1, f2, calendar_name, algo=alogo )
      values(1) = f1
      values(2) = f2
      if ( alogo == "linear" )then
         call test_close( values, expected, 1.e-10_SHR_KIND_R8, "Test if factors are as expected" )
      else
         call test_is( values, expected, "Test if factors are as expected" )
      end if
      sec_in  = sec_in + 1800
   end do
end do

!  Error tests
call shr_tInterp_setAbort( flag=.false. )

alogo = 'linear'

! lb and ub dates are the same
date_lb = 20010101
date_ub = 20010101
sec_lb  = 1457
sec_ub  = 1456
date_in = 20010101
sec_in  = 1456
call shr_tInterp_getFactors( date_lb, sec_lb, date_ub, sec_ub, date_in, &
                             sec_in, f1, f2, calendar_name, algo=alogo, rc=rc )
call test_is( rc, expected=1, description="Test that aborts if ub < lb date" )

! unrecognized alogorithm name

alogo = 'zztop'
call shr_tInterp_getFactors( date_lb, sec_lb, date_ub, sec_ub, date_in, &
                             sec_in, f1, f2, calendar_name, algo=alogo, rc=rc )
call test_is( rc, expected=1, description="Test that recognizes a bad alogo name" )

! Test that abort if input date is outside of interval of lb and ub

alogo = 'linear'
date_lb = 20010101
date_ub = 20010115
sec_lb  = 0
sec_ub  = 0
date_in = 20010205
sec_in  = 1456
call shr_tInterp_getFactors( date_lb, sec_lb, date_ub, sec_ub, date_in, &
                             sec_in, f1, f2, calendar_name, algo=alogo, rc=rc )
call test_is( rc, expected=1, description="Test that aborts for linear if input date is outside range of lb and ub dates" )

call test_final( )

end program test_shr_tInterp
