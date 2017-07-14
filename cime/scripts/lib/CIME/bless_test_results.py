import CIME.compare_namelists, CIME.simple_compare
from CIME.test_scheduler import NAMELIST_PHASE
from CIME.utils import run_cmd, expect, get_scripts_root, get_model
from CIME.test_status import *
from CIME.hist_utils import generate_baseline, compare_baseline
from CIME.case import Case

import os, glob, time

###############################################################################
def bless_namelists(test_name, report_only, force, baseline_name, baseline_root):
###############################################################################
    # Be aware that restart test will overwrite the original namelist files
    # with versions of the files that should not be blessed. This forces us to
    # re-run create_test.

    # Update namelist files
    print("Test '{}' had namelist diff".format(test_name))
    if (not report_only and
        (force or raw_input("Update namelists (y/n)? ").upper() in ["Y", "YES"])):
        create_test_gen_args = " -g {} ".format(baseline_name if get_model() == "cesm" else " -g -b {} ".format(baseline_name))
        stat, _, err = run_cmd("{}/create_test {} -n {} --baseline-root {} -o".format(get_scripts_root(), test_name, create_test_gen_args, baseline_root))
        if stat != 0:
            return False, "Namelist regen failed: '{}'".format(err)
        else:
            return True, None
    else:
        return True, None

###############################################################################
def bless_history(test_name, testcase_dir_for_test, baseline_name, baseline_root, compiler, report_only, force):
###############################################################################
    with Case(testcase_dir_for_test) as case:
        if get_model() == "acme":
            baseline_full_dir = os.path.join(baseline_root, compiler, baseline_name, case.get_value("CASEBASEID"))
        else:
            baseline_full_dir = os.path.join(baseline_root, baseline_name, case.get_value("CASEBASEID"))

        result, comments = compare_baseline(case, baseline_dir=baseline_full_dir, outfile_suffix=None)
        if result:
            return True, None
        else:
            print(comments)
            if (not report_only and
                (force or raw_input("Update this diff (y/n)? ").upper() in ["Y", "YES"])):
                result, comments = generate_baseline(case, baseline_dir=baseline_full_dir)
                if not result:
                    logging.warning("Hist file bless FAILED for test {}".format(test_name))
                    return False, "Generate baseline failed: {}".format(comments)
                else:
                    print(comments)
                    return True, None
            else:
                return True, None

###############################################################################
def bless_test_results(baseline_name, baseline_root, test_root, compiler, test_id=None, namelists_only=False, hist_only=False, report_only=False, force=False, bless_tests=None, no_skip_pass=False):
###############################################################################
    test_id_glob = "*{}*{}*".format(compiler, baseline_name) if test_id is None else "*{}".format(test_id)
    test_status_files = glob.glob("{}/{}/{}".format(test_root, test_id_glob, TEST_STATUS_FILENAME))
    expect(test_status_files, "No matching test cases found in for {}/{}/{}".format(test_root, test_id_glob, TEST_STATUS_FILENAME))

    broken_blesses = []
    for test_status_file in test_status_files:
        test_dir = os.path.dirname(test_status_file)
        ts = TestStatus(test_dir=test_dir)
        test_name = ts.get_name()
        if (bless_tests in [[], None] or CIME.utils.match_any(test_name, bless_tests)):
            overall_result = ts.get_overall_test_status()

            # See if we need to bless namelist
            if (not hist_only):
                if no_skip_pass:
                    nl_bless = True
                else:
                    nl_bless = ts.get_status(NAMELIST_PHASE) != TEST_PASS_STATUS
            else:
                nl_bless = False

            # See if we need to bless baselines
            if (not namelists_only):
                run_result = ts.get_status(RUN_PHASE)
                if (run_result is None):
                    broken_blesses.append((test_name, "no run phase"))
                    logging.warning("Test '{}' did not make it to run phase".format(test_name))
                    hist_bless = False
                elif (run_result != TEST_PASS_STATUS):
                    broken_blesses.append((test_name, "test did not pass"))
                    logging.warning("Test '{}' did not pass, not safe to bless".format(test_name))
                    hist_bless = False
                elif no_skip_pass:
                    hist_bless = True
                else:
                    hist_bless = ts.get_status(BASELINE_PHASE) != TEST_PASS_STATUS
            else:
                hist_bless = False

            # Now, do the bless
            if not nl_bless and not hist_bless:
                print("Nothing to bless for test:", test_name, " overall status:", overall_result)
            else:

                print("###############################################################################")
                print("Blessing results for test:", test_name, "most recent result:", overall_result)
                print("###############################################################################")
                if not force:
                    time.sleep(2)

                # Bless namelists
                if nl_bless:
                    success, reason = bless_namelists(test_name, report_only, force, baseline_name, baseline_root)
                    if not success:
                        broken_blesses.append((test_name, reason))

                # Bless hist files
                if hist_bless:
                    if "HOMME" in test_name:
                        success = False
                        reason = "HOMME tests cannot be blessed with bless_for_tests"
                    else:
                        success, reason = bless_history(test_name, test_dir, baseline_name, baseline_root, compiler, report_only, force)

                    if (not success):
                        broken_blesses.append((test_name, reason))

    # Make sure user knows that some tests were not blessed
    success = True
    for broken_bless, reason in broken_blesses:
        logging.warning("FAILED TO BLESS TEST: {}, reason {}".format(broken_bless, reason))
        success = False

    return success
