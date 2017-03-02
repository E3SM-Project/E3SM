"""
Library for case.cmpgen_namelists.
"""

from CIME.XML.standard_module_setup import *

from CIME.preview_namelists import create_namelists
from CIME.compare_namelists import is_namelist_file, compare_namelist_files
from CIME.simple_compare import compare_files
from CIME.utils import append_status
from CIME.test_status import *

import os, shutil, traceback, stat, glob

logger = logging.getLogger(__name__)

def _do_full_nl_comp(case, test, compare_name, baseline_root=None):
    test_dir       = case.get_value("CASEROOT")
    casedoc_dir    = os.path.join(test_dir, "CaseDocs")
    baseline_root  = case.get_value("BASELINE_ROOT") if baseline_root is None else baseline_root

    all_match         = True
    baseline_dir      = os.path.join(baseline_root, compare_name, test)
    baseline_casedocs = os.path.join(baseline_dir, "CaseDocs")

    # Start off by comparing everything in CaseDocs except a few arbitrary files (ugh!)
    # TODO: Namelist files should have consistent suffix
    all_items_to_compare = [item for item in glob.glob("%s/*" % casedoc_dir)\
                            if "README" not in os.path.basename(item)\
                            and not item.endswith("doc")\
                            and not item.endswith("prescribed")\
                            and not os.path.basename(item).startswith(".")] + \
                            glob.glob("%s/*user_nl*" % test_dir)

    comments = "NLCOMP\n"
    for item in all_items_to_compare:
        baseline_counterpart = os.path.join(baseline_casedocs \
                                            if os.path.dirname(item).endswith("CaseDocs") \
                                            else baseline_dir,os.path.basename(item))
        if not os.path.exists(baseline_counterpart):
            comments += "Missing baseline namelist '%s'\n" % baseline_counterpart
            all_match = False
        else:
            if is_namelist_file(item):
                success, current_comments = compare_namelist_files(baseline_counterpart, item, test)
            else:
                success, current_comments = compare_files(baseline_counterpart, item, test)

            all_match &= success
            if not success:
                comments += "Comparison failed between '%s' with '%s'\n" % (item, baseline_counterpart)

            comments += current_comments

    logging.info(comments)
    return all_match, comments

def _do_full_nl_gen(case, test, generate_name, baseline_root=None):
    test_dir       = case.get_value("CASEROOT")
    casedoc_dir    = os.path.join(test_dir, "CaseDocs")
    baseline_root  = case.get_value("BASELINE_ROOT") if baseline_root is None else baseline_root

    baseline_dir      = os.path.join(baseline_root, generate_name, test)
    baseline_casedocs = os.path.join(baseline_dir, "CaseDocs")

    if not os.path.isdir(baseline_dir):
        os.makedirs(baseline_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)

    if os.path.isdir(baseline_casedocs):
        shutil.rmtree(baseline_casedocs)

    shutil.copytree(casedoc_dir, baseline_casedocs)
    os.chmod(baseline_casedocs, stat.S_IRWXU | stat.S_IRWXG | stat.S_IXOTH | stat.S_IROTH)
    for item in glob.glob("%s/*" % baseline_casedocs):
        os.chmod(item, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP)

    for item in glob.glob(os.path.join(test_dir, "user_nl*")):
        preexisting_baseline = os.path.join(baseline_dir, os.path.basename(item))
        if (os.path.exists(preexisting_baseline)):
            os.remove(preexisting_baseline)

        shutil.copy2(item, baseline_dir)
        os.chmod(preexisting_baseline, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP)

def case_cmpgen_namelists(case, compare=False, generate=False, compare_name=None, generate_name=None, baseline_root=None, logfile_name="TestStatus.log"):
    expect(case.get_value("TEST"), "Only makes sense to run this for a test case")

    caseroot, casebaseid = case.get_value("CASEROOT"), case.get_value("CASEBASEID")

    if not compare:
        compare = case.get_value("COMPARE_BASELINE")
    if not generate:
        generate = case.get_value("GENERATE_BASELINE")

    if not compare and not generate:
        logging.debug("No namelists compares requested")
        return True

    # create namelists for case if they haven't been already
    casedocs = os.path.join(caseroot, "CaseDocs")
    if not os.path.exists(os.path.join(casedocs, "drv_in")):
        create_namelists(case)

    test_name = casebaseid if casebaseid is not None else case.get_value("CASE")
    with TestStatus(test_dir=caseroot, test_name=test_name) as ts:
        try:
            # Inside this try are where we catch non-fatal errors, IE errors involving
            # baseline operations which may not directly impact the functioning of the viability of this case
            if compare and not compare_name:
                compare_name = case.get_value("BASELINE_NAME_CMP")
                expect(compare_name, "Was asked to do baseline compare but unable to determine baseline name")
                logging.info("Comparing namelists with baselines '%s'" % compare_name)
            if generate and not generate_name:
                generate_name = case.get_value("BASELINE_NAME_GEN")
                expect(generate_name, "Was asked to do baseline generation but unable to determine baseline name")
                logging.info("Generating namelists to baselines '%s'" % generate_name)

            success = True
            output = ""
            if compare:
                success, output = _do_full_nl_comp(case, test_name, compare_name, baseline_root)
            if generate:
                _do_full_nl_gen(case, test_name, generate_name, baseline_root)
        except:
            ts.set_status(NAMELIST_PHASE, TEST_FAIL_STATUS)
            success = False
            warn = "Exception during namelist operations:\n%s\n%s" % (sys.exc_info()[1], traceback.format_exc())
            output += warn
            logging.warning(warn)
        finally:
            ts.set_status(NAMELIST_PHASE, TEST_PASS_STATUS if success else TEST_FAIL_STATUS)
            append_status(output, logfile_name, caseroot=caseroot)

        return success

