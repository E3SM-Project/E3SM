 program test_shr_sys
!
! Simple unit-test program for the shr_sys_mod module.
!
! Erik Kluzek
!
! $Id: test_shr_sys.F90 7482 2007-11-07 20:54:58Z erik $
!
 use shr_kind_mod, only: SHR_KIND_I8, SHR_KIND_R8
 use shr_sys_mod, only: shr_sys_irtc, shr_sys_system, shr_sys_flush, &
                        shr_sys_getenv, shr_sys_chdir, shr_sys_sleep, &
                        shr_sys_abort
 implicit none
 real(SHR_KIND_R8) sum
 integer i
#if (defined IRIX64 || defined CRAY || defined AIX)
 integer(kind=8):: irtc0, irtcf
 integer(kind=8):: irtc
#endif
 integer(SHR_KIND_I8):: sirtc0, sirtcf, rate
 integer rcode
 character(len=90) val
 real(SHR_KIND_R8) :: sec

 print *, "Unit-tester for shr_sys_mod"
 print *, "First lets test the shr_sys_irtc function"
#if (defined IRIX64 || defined CRAY || defined AIX)
 irtc0 = irtc( )
#endif
 sirtc0 = shr_sys_irtc( )
 sum = 0.0_SHR_KIND_R8
 do i = 1, 10000000
   sum = sum + exp( (i*5.0_SHR_KIND_R8*3.14159265_SHR_KIND_R8) / (i + 10.0_SHR_KIND_R8) )
 end do
 sirtcf = shr_sys_irtc( )
#if (defined IRIX64 || defined CRAY || defined AIX)
 print *, 'irtc call:         ', irtcf - irtc0
#endif
#if (defined IRIX64 || defined CRAY || defined AIX)
 irtcf = irtc( )
#endif
 print *, 'shr_sys_irtc call: ', sirtcf - sirtc0
 print *, 'Test the getenv call'
 call shr_sys_getenv( "LOGNAME", val, rcode ) 
 print *, "value of LOGNAME = ", val
 print *, 'Test the chdir call (just do a chdir .)'
 call shr_sys_system( "pwd", rcode )
 call shr_sys_chdir( ".", rcode )
 call shr_sys_system( "pwd", rcode )
 sec = 55.0_SHR_KIND_R8
 print *, 'Test the shr_sys_sleep call for a ', sec, ' second sleep'
#if (defined IRIX64 || defined CRAY || defined AIX)
 irtc0 = irtc( )
#endif
 sirtc0 = shr_sys_irtc( )
 call shr_sys_sleep( sec )
 sirtcf = shr_sys_irtc( rate )
#if (defined IRIX64 || defined CRAY || defined AIX)
 irtcf = irtc( )
#endif
#if (defined IRIX64 || defined CRAY || defined AIX)
 print *, 'irtc call:         ', irtcf - irtc0
 print *, 'irtc call:         ', irtcf, irtc0
#endif
 print *, 'shr_sys_irtc call: ', sirtcf - sirtc0, ' seconds: ', (sirtcf - sirtc0)/rate
 print *, 'shr_sys_irtc call: ', sirtcf, sirtc0
 print *, 'Test the shr_sys_flush call'
 call shr_sys_flush( 6 )
 print *, 'PASS'
 print *, 'Next test should abort appropriatly -- if it does so -- tests PASS'
 print *, 'Finally test the shr_sys_abort call'
 call shr_sys_abort
 print *, 'abort call does NOT abort code -- something is wrong'
 print *, 'FAIL'
 end program test_shr_sys
