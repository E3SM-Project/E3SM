"""
Functions for actions pertaining to history files.
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.component import Component

import logging, glob, os, shutil, re

def _iter_model_file_substrs():
    drv_comp = Component()
    models = drv_comp.get_valid_model_components()

    for model in models:
        model = "cpl" if model == "DRV" else model.lower()
        yield model

def _get_all_hist_files(testcase, model, suffix="", from_dir=os.getcwd()):
    suffix = (".%s" % suffix) if suffix else ""
    test_hists = glob.glob("%s/%s.%s*.h?.nc%s" % (from_dir, testcase, model, suffix))
    test_hists.extend(glob.glob("%s/%s.%s*.h.nc%s" % (from_dir, testcase, model, suffix)))
    test_hists.extend(glob.glob("%s/%s.%s*.h?.*.nc%s" % (from_dir, testcase, model, suffix)))
    test_hists.extend(glob.glob("%s/%s.%s*.h.*.nc%s" % (from_dir, testcase, model, suffix)))
    if suffix == "":
        test_hists.extend(glob.glob("%s/%s.h.nc" % (from_dir, model)))
        test_hists.extend(glob.glob("%s/%s.h?.nc" % (from_dir, model)))
    test_hists.sort()
    return test_hists

def move(case, suffix):
    """
    Change the suffix for the most recent batch of hist files
    """
    rundir   = case.get_value("RUNDIR")
    testcase = case.get_value("CASE")

    # Loop over models
    comments = "Moving hist files to suffix '%s'\n" % suffix
    num_moved = 0
    for model in _iter_model_file_substrs():
        comments += "  Moving hist files for model '%s'\n" % model
        test_hists = _get_all_hist_files(testcase, model, from_dir=rundir)
        num_moved += len(test_hists)
        for test_hist in test_hists:
            new_file = "%s.%s" % (test_hist, suffix)
            if os.path.exists(new_file):
                os.remove(new_file)

            comments += "    Copying '%s' to '%s'\n" % (test_hist, new_file)
            shutil.copy(test_hist, new_file)

    expect(num_moved > 0, "move failed: no hist files found in rundir '%s'" % rundir)

    return comments

def _compare_hists(case, suffix1="", suffix2="", from_dir1=os.getcwd(), from_dir2=os.getcwd()):

    testcase = case.get_value("CASE")
    all_success = True
    num_compared = 0
    comments = "Comparing hists for case '%s' dir1='%s', suffix1='%s',  dir2='%s' suffix2='%s'\n" % \
        (testcase, from_dir1, suffix1, from_dir2, suffix2)

    for model in _iter_model_file_substrs():
        comments += "  comparing model '%s'\n" % model
        hists1 = _get_all_hist_files(testcase, model, suffix1, from_dir1)
        hists2 = _get_all_hist_files(testcase, model, suffix2, from_dir2)
        if len(hists1) != len(hists2):
            comments += "    num hists does not match %d != %d\n" % (len(hists1), len(hists2))
            all_success = False
            continue

        num_compared += len(hists1)
        for hist1, hist2 in zip(hists1, hists2):
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
    Compares two component history files in the testcase directory

    returns (SUCCESS, comments)
    """
    rundir   = case.get_value("RUNDIR")

    return _compare_hists(case, suffix1, suffix2, rundir, rundir)

def cprnc(file1, file2, case, rundir):
    cprnc_exe = case.get_value("CCSM_CPRNC")
    basename = os.path.basename(file1)
    stat, out, _ = run_cmd("%s %s %s 2>&1 | tee %s/%s.cprnc.out" % (cprnc_exe, file1, file2, rundir, basename))
    return (stat == 0 and "IDENTICAL" in out, out)

def compare_baseline(case, baseline_dir=None):
    """
    compare the current test output to a baseline result
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
            return False, "ERROR %s does not exist" % bdir

    return _compare_hists(case, from_dir1=rundir, from_dir2=basecmp_dir)

def get_extension(model, filepath):
    """
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
    for model in _iter_model_file_substrs():
        comments += "  generating for model '%s'\n" % model
        hists = _get_all_hist_files(testcase, model, from_dir=rundir)
        num_gen += len(hists)
        for hist in hists:
            ext = get_extension(model, hist)
            basename = "%s.%s.nc" % (model, ext)
            baseline = os.path.join(basegen_dir, basename)
            if os.path.exists(baseline):
                os.remove(baseline)

            shutil.copy(hist, baseline)
            comments += "    generating baseline '%s'\n" % baseline

    return True, comments
