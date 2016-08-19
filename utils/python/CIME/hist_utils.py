"""
Functions for actions pertaining to history files.
"""

from CIME.XML.standard_module_setup import *

import logging, glob, os, shutil, re
logger = logging.getLogger(__name__)

def _iter_model_file_substrs(case):

    models = case.get_compset_components()
    models.append('cpl')
    for model in models:
        yield model

def _get_all_hist_files(testcase, model, from_dir, suffix=""):
    suffix = (".%s" % suffix) if suffix else ""

    # Match hist files produced by run
    test_hists = glob.glob("%s/%s.%s*.h?.nc%s" % (from_dir, testcase, model, suffix))
    test_hists.extend(glob.glob("%s/%s.%s*.h.nc%s" % (from_dir, testcase, model, suffix)))

    # Match multi-instance files produced by run
    test_hists.extend(glob.glob("%s/%s.%s*.h?.*.nc%s" % (from_dir, testcase, model, suffix)))
    test_hists.extend(glob.glob("%s/%s.%s*.h.*.nc%s" % (from_dir, testcase, model, suffix)))

    # suffix == "" implies baseline comparison, baseline hist files have simpler names
    if suffix == "":
        test_hists.extend(glob.glob("%s/%s.h.nc" % (from_dir, model)))
        test_hists.extend(glob.glob("%s/%s.h?.nc" % (from_dir, model)))

    test_hists.sort()
    return test_hists

def _get_latest_hist_files(testcase, model, from_dir, suffix=""):
    test_hists = _get_all_hist_files(testcase, model, from_dir, suffix)
    latest_files = {}
    histlist = []
    date_match = re.compile(r"(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d\d\d\d).nc")
    for hist in test_hists:
        ext = get_extension(model, hist)
        if ext in latest_files:
            s1 = date_match.search(hist)
            if s1 is None:
                latest_files[ext] = hist
                continue
            (yr,month,day,time) = s1.group(1,2,3,4)
            s2 = date_match.search(latest_files[ext])
            (pyr,pmonth,pday,ptime) = s2.group(1,2,3,4)
            if yr > pyr or (yr == pyr and month > pmonth) or \
                    (yr == pyr and month == pmonth and day > pday) or \
                    (yr == pyr and month == pmonth and day == pday and time > ptime):
                latest_files[ext] = hist
                logger.debug("ext %s hist %s %s"%(ext,hist,latest_files))
        else:
            latest_files[ext] = hist

    for key in latest_files.keys():
        histlist.append(latest_files[key])
    return histlist


def move(case, suffix):
    """
    Change the suffix for the most recent batch of hist files in a case. This can
    allow you to temporarily "save" these files so they won't be blown away if you
    re-run the case.

    case - The case containing the files you want to save
    suffix - The string suffix you want to add to saved files, this can be used to find them later.
    """
    rundir   = case.get_value("RUNDIR")
    testcase = case.get_value("CASE")

    # Loop over models
    comments = "Moving hist files to suffix '%s'\n" % suffix
    num_moved = 0
    for model in _iter_model_file_substrs(case):
        comments += "  Moving hist files for model '%s'\n" % model
        test_hists = _get_latest_hist_files(testcase, model, rundir)
        num_moved += len(test_hists)
        for test_hist in test_hists:
            new_file = "%s.%s" % (test_hist, suffix)
            if os.path.exists(new_file):
                os.remove(new_file)

            comments += "    Copying '%s' to '%s'\n" % (test_hist, new_file)
            shutil.copy(test_hist, new_file)

    expect(num_moved > 0, "move failed: no hist files found in rundir '%s'" % rundir)

    return comments

def _compare_hists(case, from_dir1, from_dir2, suffix1="", suffix2=""):
    if from_dir1 == from_dir2:
        expect(suffix1 != suffix2, "Comparing files to themselves?")

    testcase = case.get_value("CASE")
    all_success = True
    num_compared = 0
    comments = "Comparing hists for case '%s' dir1='%s', suffix1='%s',  dir2='%s' suffix2='%s'\n" % \
        (testcase, from_dir1, suffix1, from_dir2, suffix2)

    for model in _iter_model_file_substrs(case):
        comments += "  comparing model '%s'\n" % model
        hists1 = _get_latest_hist_files(testcase, model, from_dir1, suffix1)
        hists2 = _get_latest_hist_files(testcase, model, from_dir2, suffix2)
        len_hist1 = len(hists1)
        len_hist2 = len(hists2)
        if len_hist1 == 0:
            comments += " no hist files found for model %s\n"%model
            continue
        if len_hist1 > len_hist2:
            comments += "    num hists does not match %d != %d\n" % (len_hist1, len_hist2)
            all_success = False
            hists2 += ['MISSING'] * (len_hist1 - len_hist2)
        if len_hist1 < len_hist2:
            comments += "    num hists does not match %d != %d\n" % (len_hist1, len_hist2)
            all_success = False
            hists1 += ['MISSING'] * (len_hist2 - len_hist1)

        num_compared += len(hists1)
        for hist1, hist2 in zip(hists1, hists2):
            if hist1 == "MISSING":
                comments += "     No match for file %s found in %s\n"%(hist2, from_dir1)
            elif hist2 == "MISSING":
                comments += "     No match for file %s found in %s\n"%(hist1, from_dir2)
            else:
                success, cprnc_comments = cprnc(hist1, hist2, case, from_dir1)
                if success:
                    comments += "    %s matched %s\n" % (hist1, hist2)
                else:
                    comments += "    %s did NOT match %s\n" % (hist1, hist2)
                    comments += cprnc_comments + "\n"
                    all_success = False

    expect(num_compared > 0, "Did not compare any hist files for suffix1='%s' suffix2='%s', dir1='%s', dir2='%s'\nComments=%s" %
           (suffix1, suffix2, from_dir1, from_dir2, comments))
    return all_success, comments

def compare_test(case, suffix1, suffix2):
    """
    Compares two sets of component history files in the testcase directory

    case - The case containing the hist files to compare
    suffix1 - The suffix that identifies the first batch of hist files
    suffix1 - The suffix that identifies the second batch of hist files

    returns (SUCCESS, comments)
    """
    rundir   = case.get_value("RUNDIR")

    return _compare_hists(case, rundir, rundir, suffix1, suffix2)

def cprnc(file1, file2, case, rundir):
    """
    Run cprnc to compare two individual nc files

    file1 - the full or relative path of the first file
    file2 - the full or relative path of the second file
    case - the case containing the files
    rundir - the rundir for the case

    returns True if the files matched
    """
    cprnc_exe = case.get_value("CCSM_CPRNC")
    basename = os.path.basename(file1)
    stat, out, _ = run_cmd("%s -m %s %s 2>&1 | tee %s/%s.cprnc.out" % (cprnc_exe, file1, file2, rundir, basename))
    return (stat == 0 and "files seem to be IDENTICAL" in out, out)

def compare_baseline(case, baseline_dir=None):
    """
    compare the current test output to a baseline result

    case - The case containing the hist files to be compared against baselines
    baseline_dir - Optionally, specify a specific baseline dir, otherwise it will be computed from case config

    returns (SUCCESS, comments)
    SUCCESS means all hist files matched their corresponding baseline
    """
    rundir   = case.get_value("RUNDIR")
    if baseline_dir is None:
        baselineroot = case.get_value("BASELINE_ROOT")
        basecmp_dir = os.path.join(baselineroot, case.get_value("BASECMP_CASE"))
        dirs_to_check = (baselineroot, basecmp_dir)
    else:
        basecmp_dir = baseline_dir
        dirs_to_check = (basecmp_dir,)

    for bdir in dirs_to_check:
        if not os.path.isdir(bdir):
            return False, "ERROR baseline directory '%s' does not exist" % bdir

    return _compare_hists(case, rundir, basecmp_dir)

def get_extension(model, filepath):
    """
    For a hist file for the given model, return what we call the "extension"

    model - The component model
    filepath - The path of the hist file

    >>> get_extension("cpl", "cpl.hi.nc")
    'hi'
    >>> get_extension("cpl", "cpl.h.nc")
    'h'
    >>> get_extension("cpl", "cpl.h1.nc")
    'h1'
    >>> get_extension("cpl", "TESTRUNDIFF_Mmpi-serial.f19_g16_rx1.A.melvin_gnu.C.fake_testing_only_20160816_164150-20160816_164240.cpl.hi.0.nc.base")
    'hi'
    >>> get_extension("cpl", "TESTRUNDIFF_Mmpi-serial.f19_g16_rx1.A.melvin_gnu.C.fake_testing_only_20160816_164150-20160816_164240.cpl.h.nc")
    'h'
    """
    basename = os.path.basename(filepath)
    ext_regex = re.compile(r'.*%s.*[.](h.?)([.][^.]+)?[.]nc' % model)
    m = ext_regex.match(basename)
    expect(m is not None, "Failed to get extension for file '%s'" % filepath)
    return m.groups()[0]

def generate_baseline(case, baseline_dir=None):
    """
    copy the current test output to baseline result

    case - The case containing the hist files to be copied into baselines
    baseline_dir - Optionally, specify a specific baseline dir, otherwise it will be computed from case config

    returns (SUCCESS, comments)
    """
    rundir   = case.get_value("RUNDIR")
    if baseline_dir is None:
        baselineroot = case.get_value("BASELINE_ROOT")
        basegen_dir = os.path.join(baselineroot, case.get_value("BASEGEN_CASE"))
    else:
        basegen_dir = baseline_dir
    testcase = case.get_value("CASE")

    if not os.path.isdir(basegen_dir):
        os.makedirs(basegen_dir)

    comments = "Generating baselines into '%s'\n" % basegen_dir
    num_gen = 0
    for model in _iter_model_file_substrs(case):
        comments += "  generating for model '%s'\n" % model
        hists =  _get_latest_hist_files(testcase, model, rundir)
        logger.debug("latest_files: %s"%hists)
        num_gen += len(hists)
        for hist in hists:
            ext = get_extension(model, hist)
            baseline = os.path.join(basegen_dir, os.path.basename(hist))
            if os.path.exists(baseline):
                os.remove(baseline)

            shutil.copy(hist, baseline)
            comments += "    generating baseline '%s' from file %s\n" % (baseline, hist)

    expect(num_gen > 0, "Could not generate any hist files for case '%s', something is seriously wrong" % testcase)

    return True, comments
