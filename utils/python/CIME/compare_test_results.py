from __future__ import print_function

import CIME.compare_namelists, CIME.simple_compare
from CIME.utils import expect
from CIME.test_status import *
from CIME.hist_utils import compare_baseline
from CIME.case import Case

import os, glob

###############################################################################
def compare_history(testcase_dir_for_test, baseline_name, baseline_root, log_id):
###############################################################################
    with Case(testcase_dir_for_test) as case:
        baseline_full_dir = os.path.join(baseline_root, baseline_name, case.get_value("CASEBASEID"))
        outfile_suffix = "%s.%s"%(baseline_name, log_id)
        result, comments = compare_baseline(case, baseline_dir=baseline_full_dir,
                                            outfile_suffix=outfile_suffix)

        return result, comments

###############################################################################
def compare_test_results(baseline_name, baseline_root, test_root, compiler, test_id=None, compare_tests=None):
    """Compares with baselines for all matching tests

    Outputs results for each test to stdout (one line per test); possible status
    codes are: PASS, FAIL, SKIP. (A SKIP denotes a test that did not make it to
    the run phase or a test for which the run phase did not pass: we skip
    baseline comparisons in this case.)

    In addition, creates files named compare.log.BASELINE_NAME.TIMESTAMP in each
    test directory, which contain more detailed output. Also creates
    *.cprnc.out.BASELINE_NAME.TIMESTAMP files in each run directory.

    Returns True if all tests generated either PASS or SKIP results, False if
    there was at least one FAIL result.

    """
###############################################################################
    test_id_glob = "*%s*%s*" % (compiler, baseline_name) if test_id is None else "*%s" % test_id
    test_status_files = glob.glob("%s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))
    expect(test_status_files, "No matching test cases found in for %s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))

    # ID to use in the log file names, to avoid file name collisions with
    # earlier files that may exist.
    log_id = CIME.utils.get_timestamp()
    logfile_name = "compare.log.%s.%s"%(baseline_name, log_id)

    all_pass_or_skip = True

    for test_status_file in test_status_files:
        test_dir = os.path.dirname(test_status_file)
        ts = TestStatus(test_dir=test_dir)
        test_name = ts.get_name()
        if (compare_tests in [[], None] or CIME.utils.match_any(test_name, compare_tests)):
            CIME.utils.append_status(
                msg = ("Comparing against baseline with compare_test_results:\n" +
                       "Baseline: %s\n"%(baseline_name) +
                       "In baseline_root: %s"%(baseline_root)),
                caseroot = test_dir,
                sfile = logfile_name)

            run_result = ts.get_status(RUN_PHASE)
            compare_comment = ""
            detailed_comments = ""
            if (run_result is None):
                compare_result = "SKIP"
                compare_comment = "Test did not make it to run phase"
                do_compare = False
            elif (run_result != TEST_PASS_STATUS):
                compare_result = "SKIP"
                compare_comment = "Run phase did not pass"
                do_compare = False
            else:
                do_compare = True

            if do_compare:
                # Compare hist files
                success, detailed_comments = compare_history(test_dir, baseline_name, baseline_root, log_id)
                if success:
                    compare_result = TEST_PASS_STATUS
                else:
                    compare_result = TEST_FAIL_STATUS
                    all_pass_or_skip = False

                # Following the logic in SystemTestsCommon._compare_baseline:
                # We'll print the comment if it's a brief one-liner; otherwise
                # the comment will only appear in the log file
                if "\n" not in detailed_comments:
                    compare_comment = detailed_comments

            brief_result = "%s %s %s"%(compare_result, test_name, BASELINE_PHASE)
            if compare_comment:
                brief_result += " %s"%(compare_comment)
            print(brief_result)

            CIME.utils.append_status(
                msg = brief_result,
                caseroot = test_dir,
                sfile = logfile_name)

            if detailed_comments:
                CIME.utils.append_status(
                    msg = "Detailed comments:\n" + detailed_comments,
                    caseroot = test_dir,
                    sfile = logfile_name)

    return all_pass_or_skip

