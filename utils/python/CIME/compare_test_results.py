import CIME.compare_namelists, CIME.simple_compare
from CIME.utils import expect
from CIME.test_status import *
from CIME.hist_utils import compare_baseline
from CIME.case import Case

import os, glob

###############################################################################
def compare_history(testcase_dir_for_test, baseline_name, baseline_root):
###############################################################################
    with Case(testcase_dir_for_test) as case:
        baseline_full_dir = os.path.join(baseline_root, baseline_name, case.get_value("CASEBASEID"))
        result, comments = compare_baseline(case, baseline_dir=baseline_full_dir)
        if result:
            return True, None
        else:
            logging.info(comments)
            return False, "Diff'd"

###############################################################################
def compare_test_results(baseline_name, baseline_root, test_root, compiler, test_id=None, compare_tests=None):
###############################################################################
    test_id_glob = "*%s*%s*" % (compiler, baseline_name) if test_id is None else "*%s" % test_id
    test_status_files = glob.glob("%s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))
    expect(test_status_files, "No matching test cases found in for %s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))

    broken_compares = []
    for test_status_file in test_status_files:
        test_dir = os.path.dirname(test_status_file)
        ts = TestStatus(test_dir=test_dir)
        test_name = ts.get_name()
        if (compare_tests in [[], None] or CIME.utils.match_any(test_name, compare_tests)):
            overall_result = ts.get_overall_test_status()

            # Compute hist status, False implies it diffed
            run_result = ts.get_status(RUN_PHASE)
            if (run_result is None):
                broken_compares.append((test_name, "no run phase"))
                logging.warning("Test '%s' did not make it to run phase" % test_name)
                hist_no_compare = True
            elif (run_result != TEST_PASS_STATUS):
                broken_compares.append((test_name, "test did not pass"))
                logging.warning("Test '%s' did not pass, not safe to compare" % test_name)
                hist_no_compare = True
            else:
                hist_no_compare = False

            # Now, do the compare
            if hist_no_compare:
                logging.info("Cannot compare test: %s, overall status: %s" % (test_name, overall_result))
            else:

                logging.info("###############################################################################")
                logging.info("Comparing results for test: %s, most recent result: %s" % (test_name, overall_result))
                logging.info("###############################################################################")

                # Compare hist files
                if (not hist_no_compare):
                    success, reason = compare_history(test_dir, baseline_name, baseline_root)
                    if (not success):
                        broken_compares.append((test_name, reason))

    # Make sure user knows that some tests were not compareed
    success = True
    for broken_compare, reason in broken_compares:
        logging.warning("COMPARE FAILED FOR TEST: %s, reason %s" % (broken_compare, reason))
        success = False

    return success
