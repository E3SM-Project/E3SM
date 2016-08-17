import CIME.compare_namelists, CIME.simple_compare
from CIME.test_scheduler import NAMELIST_PHASE
from CIME.utils import run_cmd, expect
from CIME.XML.machines import Machines
from CIME.test_status import *
from CIME.hist_utils import generate_baseline, compare_baseline
from CIME.case import Case

import os, glob, time

_MACHINE = Machines()

###############################################################################
def bless_namelists(test_name, report_only, force, baseline_name):
###############################################################################
    # Be aware that restart test will overwrite the original namelist files
    # with versions of the files that should not be blessed. This forces us to
    # re-run create_test.

    # Update namelist files
    print "Test '%s' had a namelist diff" % test_name
    if (not report_only and
        (force or raw_input("Update namelists (y/n)? ").upper() in ["Y", "YES"])):
        stat, _, err = run_cmd("create_test -n -g %s -b %s" % (test_name, baseline_name))
        if stat != 0:
            return False, "Namelist regen failed: '%s'" % err
        else:
            return True, None

###############################################################################
def bless_history(test_name, testcase_dir_for_test, baseline_name, report_only, force):
###############################################################################
    with Case(testcase_dir_for_test) as case:
        baseline_full_dir = os.path.join(case.get_value("BASELINE_ROOT"), case.get_value("COMPILER"), baseline_name, case.get_value("CASEBASEID"))
        result, comments = compare_baseline(case, baseline_dir=baseline_full_dir)
        if result:
            logging.warning("Test '%s' was marked as DIFF but compare_baseline did not find diff?" % test_name)
            return False, "No diff found"
        else:
            print comments
            if (not report_only and
                (force or raw_input("Update this diff (y/n)? ").upper() in ["Y", "YES"])):
                result, comments = generate_baseline(case, baseline_dir=baseline_full_dir)
                if not result:
                    logging.warning("Hist file bless FAILED for test %s" % test_name)
                    return False, "Generate baseline failed: %s" % comments
                else:
                    print comments
                    return True, None
            else:
                return True, None

###############################################################################
def bless_test_results(baseline_name, test_root, compiler, test_id=None, namelists_only=False, hist_only=False, report_only=False, force=False, bless_tests=None):
###############################################################################
    test_id_glob = "*%s*%s*" % (compiler, baseline_name) if test_id is None else "*%s" % test_id
    test_status_files = glob.glob("%s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))
    expect(test_status_files, "No matching test cases found in for %s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))

    broken_blesses = []
    for test_status_file in test_status_files:
        ts = TestStatus(test_dir=os.path.dirname(test_status_file))
        test_name = ts.get_name()
        if (bless_tests in [[], None] or CIME.utils.match_any(test_name, bless_tests)):
            overall_result = ts.get_overall_test_status()

            # Compute namelists status, False implies it diffed
            if (not hist_only):
                if (ts.get_status(NAMELIST_PHASE) is not None):
                    nl_no_bless = ts.get_status(NAMELIST_PHASE) == TEST_PASS_STATUS
                else:
                    logging.warning("Test '%s' did not make it to namelist phase" % test_name)
                    broken_blesses.append((test_name, "no namelist phase"))
                    nl_no_bless = True
            else:
                nl_no_bless = True

            # Compute hist status, False implies it diffed
            if (not namelists_only):
                run_result = ts.get_status(RUN_PHASE)
                baseline_comp_result = ts.get_status("%s_baseline" % COMPARE_PHASE)
                if (run_result is None):
                    broken_blesses.append((test_name, "no run phase"))
                    logging.warning("Test '%s' did not make it to run phase" % test_name)
                    hist_no_bless = True
                elif (run_result == TEST_PASS_STATUS):
                    if (baseline_comp_result is None):
                        broken_blesses.append((test_name, "no history compare performed"))
                        logging.warning("Test '%s' had no history compare phase" % test_name)
                        hist_no_bless = True
                    else:
                        hist_no_bless = baseline_comp_result == TEST_PASS_STATUS
                else:
                    broken_blesses.append((test_name, "test did not pass"))
                    logging.warning("Test '%s' did not pass, not safe to bless" % test_name)
                    hist_no_bless = True
            else:
                hist_no_bless = True

            # Now, do the bless
            if ( (nl_no_bless and hist_no_bless) or (nl_no_bless and namelists_only) or (hist_no_bless and hist_only) ):
                print "Nothing to bless for test:", test_name, " overall status:", overall_result
            else:

                print "###############################################################################"
                print "Blessing results for test:", test_name, "most recent result:", overall_result
                print "###############################################################################"
                if not force:
                    time.sleep(2)

                # Get testcase dir for this test
                if (test_id is None):
                    # The full name already contains the compiler, so we just need to glob for the branch name
                    globs = glob.glob("%s/%s*%s*" % (test_root, test_name, baseline_name))
                else:
                    globs = glob.glob("%s/%s%s" % (test_root, test_name, test_id_glob))

                if (len(globs) != 1):
                    logging.warning("Expected exactly one match for testcase area for test '%s', found '%s'" % (test_name, globs))
                    broken_blesses.append((test_name, "multiple matching testcase dirs"))
                    continue

                testcase_dir_for_test = globs[0]

                # Bless namelists
                if (not nl_no_bless):
                    success, reason = bless_namelists(test_name, report_only, force, baseline_name)
                    if not success:
                        broken_blesses.append(test_name, reason)

                # Bless hist files
                if (not hist_no_bless):
                    success, reason = bless_history(test_name, testcase_dir_for_test, baseline_name, report_only, force)
                    if (not success):
                        broken_blesses.append((test_name, reason))

    # Make sure user knows that some tests were not blessed
    success = True
    for broken_bless, reason in broken_blesses:
        logging.warning("FAILED TO BLESS TEST: %s, reason %s" % (broken_bless, reason))
        success = False

    return success
