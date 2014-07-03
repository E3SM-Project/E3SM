program test_shr_assert

use shr_log_mod, only: &
     shr_log_unit

use test_utils, only: &
     CESM_UT_Suite, start_suite, end_suite

use logical_tests, only: &
     run_assert_tests

use domain_tests, only: &
     run_null_tests, &
     run_nan_tests, &
     run_lt_tests, &
     run_gt_tests, &
     run_le_tests, &
     run_ge_tests, &
     run_eq_tests, &
     run_ne_tests, &
     run_multi_tests, &
     run_inf_tests

type(CESM_UT_Suite) :: suite

open(20,file="./test.log")
shr_log_unit = 20

suite = start_suite("shr_assert_mod",shr_log_unit)

call run_assert_tests(suite)

call run_null_tests(suite)

call run_nan_tests(suite)

call run_lt_tests(suite)

call run_gt_tests(suite)

call run_le_tests(suite)

call run_ge_tests(suite)

call run_eq_tests(suite)

call run_ne_tests(suite)

call run_multi_tests(suite)

call run_inf_tests(suite)

call end_suite(suite)

close(shr_log_unit)

end program test_shr_assert
