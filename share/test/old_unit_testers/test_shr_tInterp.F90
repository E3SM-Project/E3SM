program test_shr_tInterp
use shr_kind_mod
use test_mod
use shr_tInterp_mod
use shr_cal_mod,     only : shr_cal_noleap, shr_cal_gregorian, shr_cal_solarDecl, shr_cal_date2julian
use shr_const_mod,   only : SHR_CONST_CDAY
use shr_orb_mod,     only : shr_orb_context_type, shr_orb_cosz

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

call test_init( num_tests*num_times+10 )
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

call test_equation_of_time()
call test_final( )

contains

subroutine test_equation_of_time()
   type(shr_orb_context_type) :: context
   real(SHR_KIND_R8) :: lon(2), lat(2), legacy(2), actual(2), reference(2), average(2)
   real(SHR_KIND_R8) :: delta, eccf, eot_seconds, calday
   integer :: mode, tod, point
   character(SHR_KIND_CS) :: calendar

   context%eccen = 0.0167_SHR_KIND_R8
   context%obliqr = 0.409_SHR_KIND_R8
   context%lambm0 = -0.0327_SHR_KIND_R8
   context%mvelpp = 4.94_SHR_KIND_R8
   lon = (/0._SHR_KIND_R8, 0.2_SHR_KIND_R8/)
   lat = (/0.4_SHR_KIND_R8, -0.2_SHR_KIND_R8/)
   context%use_equation_of_time = .false.
   call shr_tInterp_getCosz(legacy,lon,lat,20000320,43200,context%eccen,context%mvelpp, &
      context%lambm0,context%obliqr,shr_cal_noleap)
   call shr_tInterp_getCosz(actual,lon,lat,20000320,43200,context%eccen,context%mvelpp, &
      context%lambm0,context%obliqr,shr_cal_noleap,context)
   call test_is(all(transfer(actual,(/0_SHR_KIND_I8,0_SHR_KIND_I8/)) == &
      transfer(legacy,(/0_SHR_KIND_I8,0_SHR_KIND_I8/))), .true., 'EOT-off forcing remains bitwise unchanged')

   context%use_equation_of_time = .true.
   do mode = 1, 2
    calendar = shr_cal_noleap
    if (mode == 2) calendar = shr_cal_gregorian
    call shr_cal_solarDecl(context,20000320,43200._SHR_KIND_R8,calendar,delta,eccf,eot_seconds)
    call shr_cal_date2julian(20000320,43200,calday,calendar)
    do point = 1, 2
       reference(point) = shr_orb_cosz(calday,lat(point),lon(point),delta,eot_seconds=eot_seconds)
    end do
    call shr_tInterp_getCosz(actual,lon,lat,20000320,43200,context%eccen,context%mvelpp, &
       context%lambm0,context%obliqr,calendar,context)
    call test_close(actual,reference,1.e-13_SHR_KIND_R8,'Forcing instantaneous EOT matches orbital context')

    reference = 0._SHR_KIND_R8
    do tod = 43200, 45600, 1200
       call shr_cal_solarDecl(context,20000320,real(tod,SHR_KIND_R8),calendar,delta,eccf,eot_seconds)
       call shr_cal_date2julian(20000320,tod,calday,calendar)
       do point = 1, 2
        reference(point) = reference(point) + &
           shr_orb_cosz(calday,lat(point),lon(point),delta,eot_seconds=eot_seconds) / 3._SHR_KIND_R8
       end do
    end do
    call shr_tInterp_getAvgCosz(average,lon,lat,20000320,43200,20000320,46800, &
       context%eccen,context%mvelpp,context%lambm0,context%obliqr,1200,calendar,context)
    call test_close(average,reference,1.e-13_SHR_KIND_R8,'Forcing average recomputes EOT at each sample')
    call shr_tInterp_getCosz(actual,lon,lat,20000320,0,context%eccen,context%mvelpp, &
       context%lambm0,context%obliqr,calendar,context)
    call test_is(all(actual == 0._SHR_KIND_R8), .true., 'Enabled forcing has zero nighttime normalization')
   end do
end subroutine test_equation_of_time

end program test_shr_tInterp
