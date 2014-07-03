module domain_tests

use test_utils, only: &
     CESM_UT_TEST, CESM_UT_Suite, &
     start_test, pass_if_true

use shr_kind_mod, only: &
     r4 => shr_kind_r4, &
     r8 => shr_kind_r8, &
     i4 => shr_kind_i4, &
     i8 => shr_kind_i8

! Use quiet NaN so this works even if floating-point trapping is on.
use shr_infnan_mod, only: &
     nan => shr_infnan_qnan, &
     nan_type => shr_infnan_nan_type, &
     assignment(=)

use shr_sys_mod, only: &
     pull_aborted

use shr_assert_mod, only: &
     shr_assert_in_domain

implicit none
private
save

public :: run_null_tests
public :: run_nan_tests
public :: run_lt_tests
public :: run_gt_tests
public :: run_le_tests
public :: run_ge_tests
public :: run_eq_tests
public :: run_ne_tests
public :: run_multi_tests
public :: run_inf_tests

type(CESM_UT_Test) :: test

! To get a grip on the mind-boggling combinatorics here,
! limit unit tests to these four dimension/kind combinations.
integer(i4) :: var_1d_i4(3)
integer(i8) :: var_0d_i8
real(r4) :: var_0d_r4
real(r8) :: var_2d_r8(2,2)

contains

subroutine run_null_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Assert absolutely nothing about each variable.

  test = start_test(suite, "assert nothing i4 1d")
  var_1d_i4 = 0_i4
  call shr_assert_in_domain(var_1d_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nothing i8 0d")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nothing r4 0d")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nothing r8 2d")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8)
  call pass_if_true(suite, test, .not. pull_aborted())

end subroutine run_null_tests

subroutine run_nan_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that a variable is (not) NaN.
  ! Note that integers are never NaN.

  ! 1d r4

  test = start_test(suite, "assert not nan i4 1d (should succeed)")
  var_1d_i4 = 0_i4
  call shr_assert_in_domain(var_1d_i4, is_nan=.false.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nan i4 1d (should abort)")
  var_1d_i4 = 1_i4
  call shr_assert_in_domain(var_1d_i4, is_nan=.true.)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert not nan i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, is_nan=.false.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nan i8 0d (should abort)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, is_nan=.true.)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert not nan r4 0d (should succeed)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, is_nan=.false.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nan r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, is_nan=.true.)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert nan r4 0d (NaN, should succeed)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, is_nan=.true.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert not nan r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, is_nan=.false.)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert not nan r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  call shr_assert_in_domain(var_2d_r8, is_nan=.false.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert nan r8 2d (should abort)")
  var_2d_r8 = nan
  var_2d_r8(1,2) = 0._r8
  call shr_assert_in_domain(var_2d_r8, is_nan=.true.)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert nan r8 2d (NaN, should succeed)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, is_nan=.true.)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert not nan r8 2d (NaN, should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,2) = nan
  call shr_assert_in_domain(var_2d_r8, is_nan=.false.)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_nan_tests

subroutine run_lt_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is less than some other value.

  ! 1d i4

  test = start_test(suite, "assert lt i4 1d (should succeed)")
  var_1d_i4 = 0_i4
  call shr_assert_in_domain(var_1d_i4, lt=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert lt i4 1d (should abort)")
  var_1d_i4 = 0_i4
  var_1d_i4(2) = 1_i4
  call shr_assert_in_domain(var_1d_i4, lt=1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert lt i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, lt=1_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert lt i8 0d (should abort)")
  var_0d_i8 = 1_i8
  call shr_assert_in_domain(var_0d_i8, lt=1_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert lt r4 0d (should succeed)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, lt=1._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert lt r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, lt=-1._r4)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert lt r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, lt=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert lt r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  call shr_assert_in_domain(var_2d_r8, lt=1._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert lt r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = 1._r8
  call shr_assert_in_domain(var_2d_r8, lt=1._r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert lt r8 2d (NaN, should abort)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, lt=1._r8)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_lt_tests

subroutine run_gt_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is less than some other value.

  ! 1d i4

  test = start_test(suite, "assert gt i4 1d (should succeed)")
  var_1d_i4 = 0_i4
  call shr_assert_in_domain(var_1d_i4, gt=-1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert gt i4 1d (should abort)")
  var_1d_i4 = 0_i4
  var_1d_i4(2) = -1_i4
  call shr_assert_in_domain(var_1d_i4, gt=-1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert gt i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, gt=-1_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert gt i8 0d (should abort)")
  var_0d_i8 = -1_i8
  call shr_assert_in_domain(var_0d_i8, gt=-1_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert gt r4 0d (should succeed)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, gt=-1._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert gt r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, gt=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert gt r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, gt=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert gt r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  call shr_assert_in_domain(var_2d_r8, gt=-1._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert gt r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = -1._r8
  call shr_assert_in_domain(var_2d_r8, gt=-1._r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert gt r8 2d (NaN, should abort)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, gt=1._r8)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_gt_tests

subroutine run_le_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is less than or equal to some other
  ! value.

  ! 1d i4

  test = start_test(suite, "assert le i4 1d (should succeed)")
  var_1d_i4 = 1_i4
  call shr_assert_in_domain(var_1d_i4, le=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert le i4 1d (should abort)")
  var_1d_i4 = 0_i4
  var_1d_i4(2) = 2_i4
  call shr_assert_in_domain(var_1d_i4, le=1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert le i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, le=1_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert le i8 0d (should abort)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, le=-1_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert le r4 0d (should succeed)")
  var_0d_r4 = -1._r4
  call shr_assert_in_domain(var_0d_r4, le=-1._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert le r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, le=-1._r4)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert le r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, le=0._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert le r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = 1._r8
  call shr_assert_in_domain(var_2d_r8, le=1._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert le r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = 1._r8
  call shr_assert_in_domain(var_2d_r8, le=0._r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert le r8 2d (NaN, should abort)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, le=1._r8)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_le_tests

subroutine run_ge_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is greater than or equal to some
  ! other value.

  ! 1d i4

  test = start_test(suite, "assert ge i4 1d (should succeed)")
  var_1d_i4 = 1_i4
  call shr_assert_in_domain(var_1d_i4, ge=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ge i4 1d (should abort)")
  var_1d_i4 = 0_i4
  var_1d_i4(2) = 2_i4
  call shr_assert_in_domain(var_1d_i4, ge=1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert ge i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, ge=-1_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ge i8 0d (should abort)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, ge=1_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert ge r4 0d (should succeed)")
  var_0d_r4 = -1._r4
  call shr_assert_in_domain(var_0d_r4, ge=-1._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ge r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, ge=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert ge r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, ge=0._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert ge r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = 1._r8
  call shr_assert_in_domain(var_2d_r8, ge=-1._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ge r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = -1._r8
  call shr_assert_in_domain(var_2d_r8, ge=0._r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert ge r8 2d (NaN, should abort)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, ge=-1._r8)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_ge_tests

subroutine run_eq_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is equal to some other value.

  ! 1d i4

  test = start_test(suite, "assert eq i4 1d (should succeed)")
  var_1d_i4 = 1_i4
  call shr_assert_in_domain(var_1d_i4, eq=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert eq i4 1d (should abort)")
  var_1d_i4 = 0_i4
  var_1d_i4(2) = 2_i4
  call shr_assert_in_domain(var_1d_i4, eq=1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert eq i8 0d (should succeed)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, eq=0_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert eq i8 0d (should abort)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, eq=1_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert eq r4 0d (should succeed)")
  var_0d_r4 = -1._r4
  call shr_assert_in_domain(var_0d_r4, eq=-1._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert eq r4 0d (should abort)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, eq=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert eq r4 0d (NaN, should abort)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, eq=0._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert eq r8 2d (should succeed)")
  var_2d_r8 = 0._r8
  call shr_assert_in_domain(var_2d_r8, eq=0._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert eq r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(2,1) = -1._r8
  call shr_assert_in_domain(var_2d_r8, eq=0._r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert eq r8 2d (NaN, should abort)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, eq=var_2D_r8(1,1))
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_eq_tests

subroutine run_ne_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Check assertions that the variable is not equal to some other value.

  ! 1d i4

  test = start_test(suite, "assert ne i4 1d (should succeed)")
  var_1d_i4 = 0_i4
  call shr_assert_in_domain(var_1d_i4, ne=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne i4 1d (should abort)")
  var_1d_i4 = 1_i4
  var_1d_i4(2) = 2_i4
  call shr_assert_in_domain(var_1d_i4, ne=1_i4)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d i8

  test = start_test(suite, "assert ne i8 0d (should succeed)")
  var_0d_i8 = 1_i8
  call shr_assert_in_domain(var_0d_i8, ne=0_i8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne i8 0d (should abort)")
  var_0d_i8 = 0_i8
  call shr_assert_in_domain(var_0d_i8, ne=0_i8)
  call pass_if_true(suite, test, pull_aborted())

  ! 0d r4

  test = start_test(suite, "assert ne r4 0d (should succeed)")
  var_0d_r4 = -1._r4
  call shr_assert_in_domain(var_0d_r4, ne=0._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne r4 0d (NaN, should succeed)")
  var_0d_r4 = nan
  call shr_assert_in_domain(var_0d_r4, ne=0._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne r4 0d (should abort)")
  var_0d_r4 = 1._r4
  call shr_assert_in_domain(var_0d_r4, ne=1._r4)
  call pass_if_true(suite, test, pull_aborted())

  ! 2d r8

  test = start_test(suite, "assert ne r8 2d (should succeed)")
  var_2d_r8 = -1._r8
  var_2d_r8(2,1) = 2._r8
  call shr_assert_in_domain(var_2d_r8, ne=0._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne r8 2d (NaN, should succeed)")
  var_2d_r8 = nan
  call shr_assert_in_domain(var_2d_r8, ne=var_2D_r8(1,1))
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ne r8 2d (should abort)")
  var_2d_r8 = 0._r8
  call shr_assert_in_domain(var_2d_r8, ne=0._r8)
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_ne_tests

subroutine run_multi_tests( suite )

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Test multiple assertions at once.


  test = start_test(suite, "assert multi i4 1d (should succeed)")
  var_1d_i4 = 1_i4
  call shr_assert_in_domain(var_1d_i4, is_nan=.false., le=5_i4, eq=1_i4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert multi i8 0d (should abort)")
  var_0d_i8 = -2_i8
  call shr_assert_in_domain(var_0d_i8, is_nan=.true., ne=-2_i8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert multi r4 0d (should succeed)")
  var_0d_r4 = 0._r4
  call shr_assert_in_domain(var_0d_r4, lt=1._r4, ge=0._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert multi r8 2d (should abort)")
  var_2d_r8 = 0._r8
  var_2d_r8(1,1) = nan
  call shr_assert_in_domain(var_2d_r8, is_nan=.false., ne=var_2d_r8(1,1))
  call pass_if_true(suite, test, pull_aborted())

end subroutine run_multi_tests

subroutine run_inf_tests( suite )

  use shr_infnan_mod, only: &
       pinf => shr_infnan_posinf, &
       ninf => shr_infnan_neginf

  ! Suite to run in.
  type(CESM_UT_Suite), intent(inout) :: suite

  ! Local Inf/NaN
  real(r4) :: pinf_r4, ninf_r4, nan_r4
  real(r8) :: pinf_r8, ninf_r8, nan_r8

  pinf_r4 = pinf
  ninf_r4 = ninf
  nan_r4  = nan
  pinf_r8 = pinf
  ninf_r8 = ninf
  nan_r8  = nan

  ! Run some of the above tests with infinite values.

  test = start_test(suite, "assert pinf vs inf r4 0d (should succeed)")
  var_0d_r4 = pinf
  call shr_assert_in_domain(var_0d_r4, le=pinf_r4, gt=ninf_r4, &
       eq=pinf_r4, ne=ninf_r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert pinf vs nan r4 0d (should abort)")
  var_0d_r4 = pinf
  call shr_assert_in_domain(var_0d_r4, is_nan=.true.)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert ninf vs nan r4 0d (should succeed)")
  var_0d_r4 = ninf
  call shr_assert_in_domain(var_0d_r4, is_nan=.false., ne=nan_r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ninf vs finite r4 0d (should succeed)")
  var_0d_r4 = ninf
  call shr_assert_in_domain(var_0d_r4, lt=2._r4, le=-3._r4, ne=-10._r4)
  call pass_if_true(suite, test, .not. pull_aborted())

  test = start_test(suite, "assert ninf vs inf r8 2d (should abort)")
  var_2d_r8 = ninf
  call shr_assert_in_domain(var_2d_r8, eq=pinf_r8)
  call pass_if_true(suite, test, pull_aborted())

  test = start_test(suite, "assert mixed vs finite r8 2d (should succeed)")
  var_2d_r8 = pinf
  var_2d_r8(1,2) = 1._r8
  call shr_assert_in_domain(var_2d_r8, is_nan=.false., gt=ninf_r8, &
       ge=1._r8, ne=5._r8)
  call pass_if_true(suite, test, .not. pull_aborted())

end subroutine run_inf_tests

end module domain_tests
