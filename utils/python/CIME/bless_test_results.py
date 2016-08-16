import CIME.compare_namelists, CIME.simple_compare
from CIME.test_scheduler import NAMELIST_PHASE
from CIME.utils import run_cmd, run_cmd_no_fail, expect
from CIME.XML.machines import Machines
from CIME.test_status import *

import os, glob, time

_MACHINE = Machines()

###############################################################################
def bless_namelists(test_name, baseline_dir_for_test, testcase_dir_for_test, report_only, force):
###############################################################################
    namelist_files = []

    for root, _, files in os.walk(baseline_dir_for_test):
        if (root == baseline_dir_for_test):
            rel_root = ""
        else:
            rel_root = root.replace("%s/" % baseline_dir_for_test, "")

        for file_ in files:
            rel_file = os.path.join(rel_root, file_)

            baseline_file = os.path.join(baseline_dir_for_test, rel_file)
            testcase_file = os.path.join(testcase_dir_for_test, rel_file)

            if (os.path.isfile(testcase_file) and (rel_root == "CaseDocs" or file_.startswith("user_nl"))):
                if (CIME.compare_namelists.is_namelist_file(baseline_file)):
                    if (not CIME.compare_namelists.compare_namelist_files(baseline_file, testcase_file, test_name)):
                        print "Namelist files '%s' and '%s' did not match" % (baseline_file, testcase_file)
                        print
                        if (not report_only and
                            (force or raw_input("Update this file (y/n)? ").upper() in ["Y", "YES"])):
                            namelist_files.append((rel_file, rel_file))
                else:
                    if (not CIME.simple_compare.compare_files(baseline_file, testcase_file, test_name)):
                        print "Simple files '%s' and '%s' did not match" % (baseline_file, testcase_file)
                        print
                        if (not report_only and
                            (force or raw_input("Update this file (y/n)? ").upper() in ["Y", "YES"])):
                            namelist_files.append((rel_file, rel_file))

    # Update namelist files
    if (namelist_files):
        CIME.utils.safe_copy(testcase_dir_for_test, baseline_dir_for_test, namelist_files)

###############################################################################
def bless_history(test_name, baseline_tag, baseline_dir_for_test, testcase_dir_for_test, report_only, force):
###############################################################################
    # Get user that test was run as (affects loc of hist files)
    acme_root = run_cmd_no_fail("./xmlquery CESMSCRATCHROOT -value", from_dir=testcase_dir_for_test)

    case = os.path.basename(testcase_dir_for_test)
    cime_root = CIME.utils.get_cime_root()
    compgen = os.path.join(cime_root, "scripts", "Tools", "component_compgen_baseline.sh")
    machine_env = os.path.join(testcase_dir_for_test, ".env_mach_specific.sh")
    run_dir = os.path.join(acme_root, case, "run")
    cprnc_loc = _MACHINE.get_value("CCSM_CPRNC")
    compgen_cmd = "source %s && %s -baseline_dir %s -testcase %s -testcase_base %s -test_dir %s -cprnc_exe %s -generate_tag %s" % \
                  (machine_env, compgen, baseline_dir_for_test, case, test_name, run_dir, cprnc_loc, baseline_tag)

    check_compare = os.path.join(cime_root, "scripts", "Tools", "component_write_comparefail.pl")
    stat, out, _ = run_cmd("%s %s 2>&1" % (check_compare, run_dir))

    if (stat != 0):
        # found diff, offer rebless
        print out

        if (not report_only and
            (force or raw_input("Update this diff (y/n)? ").upper() in ["Y", "YES"])):
            stat = run_cmd(compgen_cmd, verbose=True)[0]
            if (stat != 0):
                logging.warning("Hist file bless FAILED for test %s" % test_name)
                return False, "Bless command failed"
            else:
                return True, None
        else:
            return True, None
    else:
        logging.warning("Test '%s' was marked as DIFF but cprnc did not find diff?" % test_name)
        return False, "No diff found or missing baseline"

###############################################################################
def bless_test_results(baseline_name, test_root, compiler, test_id=None, namelists_only=False, hist_only=False, report_only=False, force=False, bless_tests=None):
###############################################################################
    test_id_glob = "*%s*%s*" % (compiler, baseline_name) if test_id is None else "*%s" % test_id
    test_status_files = glob.glob("%s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))
    expect(test_status_files, "No matching test cases found in for %s/%s/%s" % (test_root, test_id_glob, TEST_STATUS_FILENAME))

    baseline_root = _MACHINE.get_value("CCSM_BASELINE")
    baseline_tag  = os.path.join(compiler, baseline_name)
    baseline_area = os.path.join(baseline_root, baseline_tag)

    # The env_mach_specific script may need these to be defined
    compiler = _MACHINE.get_default_compiler() # this MUST match compiler that cprnc was built with
    os.environ["COMPILER"]  = compiler
    os.environ["MPILIB"] = _MACHINE.get_default_MPIlib(attributes={"compiler":compiler})

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
                time.sleep(2)

                # Get baseline dir for this test
                baseline_dir_for_test = os.path.join(baseline_area, test_name)
                if (not os.path.isdir(baseline_dir_for_test)):
                    logging.warning("Problem, baseline dir '%s' does not exist" % baseline_dir_for_test)
                    broken_blesses.append((test_name, "missing baseline dir"))
                    continue

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
                    bless_namelists(test_name, baseline_dir_for_test, testcase_dir_for_test, report_only, force)

                # Bless hist files
                if (not hist_no_bless):
                    success, reason = bless_history(test_name, baseline_tag, baseline_dir_for_test, testcase_dir_for_test, report_only, force)
                    if (not success):
                        broken_blesses.append((test_name, reason))

    # Make sure user knows that some tests were not blessed
    for broken_bless, reason in broken_blesses:
        logging.warning("FAILED TO BLESS TEST: %s, reason %s" % (broken_bless, reason))
