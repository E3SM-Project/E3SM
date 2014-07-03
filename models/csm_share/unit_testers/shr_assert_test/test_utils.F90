module test_utils

use shr_kind_mod, only: &
     shr_kind_in, &
     shr_kind_cs, &
     shr_kind_cl

use shr_sys_mod, only: &
     shr_sys_flush

implicit none
private
save

! List public members
public :: CESM_UT_Test
public :: CESM_UT_Suite
public :: start_suite
public :: start_test
public :: pass_if_true
public :: set_passed
public :: set_failed
public :: end_suite

! Test type
type :: CESM_UT_Test

   character(len=shr_kind_cl) :: name
   integer(shr_kind_in) :: num

end type CESM_UT_Test

! Suite type
type :: CESM_UT_Suite

   character(len=shr_kind_cl) :: name
   integer(shr_kind_in) :: log_unit
   integer(shr_kind_in) :: num_tests = 0
   integer(shr_kind_in) :: num_passed = 0
   integer(shr_kind_in) :: num_failed = 0

end type CESM_UT_Suite

! Internal test status codes
integer(shr_kind_in), parameter :: pass_code = 1
integer(shr_kind_in), parameter :: fail_code = 2

contains

function start_suite(name, log_unit) result(suite)
  character(len=*), intent(in) :: name
  integer(shr_kind_in), intent(in) :: log_unit
  type(CESM_UT_Suite) :: suite

! Should be able to do the following, but PGI 11 has bug.
!  suite = CESM_UT_Suite(trim(name), log_unit)
  suite%name = trim(name)
  suite%log_unit = log_unit

  write(suite%log_unit,*) "Starting suite: ",trim(suite%name)
  write(suite%log_unit,*) ""
  call shr_sys_flush(suite%log_unit)

end function start_suite

function start_test(suite, name) result(test)
  type(CESM_UT_suite), intent(inout) :: suite
  character(len=*), intent(in) :: name
  type(CESM_UT_Test) :: test

  suite%num_tests = suite%num_tests + 1
  test = CESM_UT_Test(trim(name), suite%num_tests)

  write(suite%log_unit,"(A,I5,A,A)") &
       "Starting test #",test%num,": ",trim(test%name)
  call shr_sys_flush(suite%log_unit)

end function start_test

subroutine pass_if_true(suite, test, flag)
  type(CESM_UT_suite), intent(inout) :: suite
  type(CESM_UT_Test), intent(inout) :: test
  logical, intent(in) :: flag

  if (flag) then
     call set_passed(suite, test)
  else
     call set_failed(suite, test)
  end if

end subroutine pass_if_true

subroutine set_passed(suite, test)
  type(CESM_UT_suite), intent(inout) :: suite
  type(CESM_UT_Test), intent(inout) :: test

  suite%num_passed = suite%num_passed + 1

  call print_test_status(suite, test, pass_code)

end subroutine set_passed

subroutine set_failed(suite, test)
  type(CESM_UT_suite), intent(inout) :: suite
  type(CESM_UT_Test), intent(inout) :: test

  suite%num_failed = suite%num_failed + 1

  call print_test_status(suite, test, fail_code)

end subroutine set_failed

subroutine print_test_status(suite, test, code)
  type(CESM_UT_suite), intent(inout) :: suite
  type(CESM_UT_Test), intent(inout) :: test
  integer(shr_kind_in), intent(in) :: code

  write(suite%log_unit,"(A,I5,4A)") &
       "Test #",test%num,", ",trim(test%name),", returned with status: ", &
       stat_c2s(code)
  write(suite%log_unit,*) ""
  call shr_sys_flush(suite%log_unit)

end subroutine print_test_status

subroutine end_suite(suite)
  type(CESM_UT_suite), intent(inout) :: suite

  write(suite%log_unit,*) "End of suite ",trim(suite%name)
  write(suite%log_unit,*) "Total tests run: ",suite%num_tests
  write(suite%log_unit,*) "Reported passed: ",suite%num_passed
  write(suite%log_unit,*) "Reported failed: ",suite%num_failed
  if (suite%num_tests == suite%num_passed .and. &
       suite%num_failed == 0) then
     write(suite%log_unit,*) "Overall status: PASSED"
  else
     write(suite%log_unit,*) "Overall status: FAILED"
  end if

  call shr_sys_flush(suite%log_unit)

end subroutine end_suite

! Internal code to string conversion utility.
pure function stat_c2s(code) result(str)
  integer(shr_kind_in), intent(in) :: code
  character(len=shr_kind_cs) :: str

  select case (code)
  case (pass_code)
     str = "PASS"
  case (fail_code)
     str = "FAIL"
  case default
     str = "UNKNOWN"
  end select

end function stat_c2s

end module test_utils
