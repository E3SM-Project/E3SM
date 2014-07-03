module logical_tests

use test_utils, only: &
     CESM_UT_TEST, CESM_UT_Suite, &
     start_test, pass_if_true

use shr_sys_mod, only: &
     pull_aborted

implicit none
private
save

public :: run_assert_tests

type(CESM_UT_Test) :: test

contains

subroutine run_assert_tests( suite )

  use shr_assert_mod, only: &
       shr_assert

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  logical :: test_array(5)

  test = start_test(suite, "vanilla assert true")
  call shr_assert(.true.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "vanilla assert false")
  call shr_assert(.false.)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "vanilla assert array true")
  test_array = .true.
  call shr_assert(test_array)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "vanilla assert array mixed")
  test_array = .true.
  test_array(4) = .false.
  call shr_assert(test_array)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "vanilla assert array false")
  test_array = .false.
  call shr_assert(test_array)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_assert_tests

end module logical_tests
