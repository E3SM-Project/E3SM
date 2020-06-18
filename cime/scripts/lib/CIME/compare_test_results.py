import CIME.compare_namelists, CIME.simple_compare
from CIME.utils import append_status, EnvironmentContext
from CIME.test_status import *
from CIME.hist_utils import compare_baseline, get_ts_synopsis
from CIME.case import Case
from CIME.test_utils import get_test_status_files

import os, logging

###############################################################################
def append_status_cprnc_log(msg, logfile_name, test_dir):
###############################################################################
    try:
        append_status(msg, logfile_name, caseroot=test_dir)
    except IOError:
        pass

###############################################################################
def compare_namelists(case, baseline_name, baseline_root, logfile_name):
###############################################################################
    log_lvl = logging.getLogger().getEffectiveLevel()
    logging.disable(logging.CRITICAL)
    success = case.case_cmpgen_namelists(compare=True, compare_name=baseline_name, baseline_root=baseline_root, logfile_name=logfile_name)
    logging.getLogger().setLevel(log_lvl)
    return success

###############################################################################
def compare_history(case, baseline_name, baseline_root, log_id):
###############################################################################
    real_user = case.get_value("REALUSER")
    with EnvironmentContext(USER=real_user):
        baseline_full_dir = os.path.join(baseline_root, baseline_name, case.get_value("CASEBASEID"))

        outfile_suffix = "{}.{}".format(baseline_name, log_id)
        try:
            result, comments = compare_baseline(case, baseline_dir=baseline_full_dir,
                                                outfile_suffix=outfile_suffix)
        except IOError:
            result, comments = compare_baseline(case, baseline_dir=baseline_full_dir,
                                                outfile_suffix=None)

        return result, comments

###############################################################################
def compare_test_results(baseline_name, baseline_root, test_root, compiler, test_id=None, compare_tests=None, namelists_only=False, hist_only=False):
###############################################################################
    """
    Compares with baselines for all matching tests

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
    test_status_files = get_test_status_files(test_root, compiler, test_id=test_id)

    # ID to use in the log file names, to avoid file name collisions with
    # earlier files that may exist.
    log_id = CIME.utils.get_timestamp()

    all_pass_or_skip = True

    for test_status_file in test_status_files:
        test_dir = os.path.dirname(test_status_file)
        ts = TestStatus(test_dir=test_dir)
        test_name = ts.get_name()
        if (compare_tests in [[], None] or CIME.utils.match_any(test_name, compare_tests)):

            if (not hist_only):
                nl_compare_result = None
                nl_compare_comment = ""
                nl_result = ts.get_status(SETUP_PHASE)
                if (nl_result is None):
                    nl_compare_result = "SKIP"
                    nl_compare_comment = "Test did not make it to setup phase"
                    nl_do_compare = False
                else:
                    nl_do_compare = True
            else:
                nl_do_compare = False

            detailed_comments = ""
            if (not namelists_only):
                compare_result = None
                compare_comment = ""
                run_result = ts.get_status(RUN_PHASE)
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
            else:
                do_compare = False

            if nl_do_compare or do_compare:
                with Case(test_dir) as case:

                    if baseline_name is None:
                        baseline_name = case.get_value("BASELINE_NAME_CMP")
                        if not baseline_name:
                            baseline_name = CIME.utils.get_current_branch(repo=CIME.utils.get_cime_root())

                    if baseline_root is None:
                        baseline_root = case.get_value("BASELINE_ROOT")

                    logfile_name = "compare.log.{}.{}".format(baseline_name.replace("/", "_"), log_id)

                    append_status_cprnc_log(
                        "Comparing against baseline with compare_test_results:\n"
                        "Baseline: {}\n In baseline_root: {}".format(baseline_name, baseline_root),
                        logfile_name,
                        test_dir)

                    if nl_do_compare:
                        nl_success = compare_namelists(case, baseline_name, baseline_root, logfile_name)
                        if nl_success:
                            nl_compare_result = TEST_PASS_STATUS
                            nl_compare_comment = ""
                        else:
                            nl_compare_result = TEST_FAIL_STATUS
                            nl_compare_comment = "See {}/{}".format(test_dir, logfile_name)
                            all_pass_or_skip = False

                    if do_compare:
                        success, detailed_comments = compare_history(case, baseline_name, baseline_root, log_id)
                        if success:
                            compare_result = TEST_PASS_STATUS
                        else:
                            compare_result = TEST_FAIL_STATUS
                            all_pass_or_skip = False

                        compare_comment = get_ts_synopsis(detailed_comments)

            brief_result = ""
            if not hist_only:
                brief_result += "{} {} {} {}\n".format(nl_compare_result, test_name, NAMELIST_PHASE, nl_compare_comment)

            if not namelists_only:
                brief_result += "{} {} {}".format(compare_result, test_name, BASELINE_PHASE)
                if compare_comment:
                    brief_result += " {}".format(compare_comment)
                brief_result += "\n"

            print(brief_result)

            append_status_cprnc_log(brief_result, logfile_name, test_dir)

            if detailed_comments:
                append_status_cprnc_log("Detailed comments:\n" + detailed_comments, logfile_name, test_dir)

    return all_pass_or_skip
