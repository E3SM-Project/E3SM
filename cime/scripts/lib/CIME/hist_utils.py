"""
Functions for actions pertaining to history files.
"""
from CIME.XML.standard_module_setup import *
from CIME.test_status import TEST_NO_BASELINES_COMMENT, TEST_STATUS_FILENAME
from CIME.utils import get_current_commit, get_timestamp, get_model, safe_copy

import logging, os, re, stat, filecmp
logger = logging.getLogger(__name__)

BLESS_LOG_NAME = "bless_log"

NO_COMPARE     = "had no compare counterpart"
NO_ORIGINAL    = "had no original counterpart"
DIFF_COMMENT   = "did NOT match"
# COMPARISON_COMMENT_OPTIONS should include all of the above: these are any of the special
# comment strings that describe the reason for a comparison failure
COMPARISON_COMMENT_OPTIONS = set([NO_COMPARE,
                                  NO_ORIGINAL,
                                  DIFF_COMMENT])
# Comments that indicate a true baseline comparison problem - not simply a BFAIL
NON_BFAIL_COMMENT_OPTIONS = COMPARISON_COMMENT_OPTIONS - set([NO_COMPARE])

def _iter_model_file_substrs(case):
    models = case.get_compset_components()
    models.append('cpl')
    for model in models:
        yield model

def _get_all_hist_files(model, from_dir, file_extensions, suffix="", ref_case=None):
    suffix = (".{}".format(suffix) if suffix else "")

    test_hists = []
    # Match hist files produced by run
    for extension in file_extensions:
        if 'initial' in extension:
            continue
        if extension.endswith('$'):
            extension = extension[:-1]
        string = model+r'\d?_?(\d{4})?\.'+extension+suffix+'$'
        logger.debug ("Regex is {}".format(string))
        pfile = re.compile(string)
        test_hists.extend([os.path.join(from_dir,f) for f in os.listdir(from_dir) if pfile.search(f)])

    if ref_case:
        test_hists = [h for h in test_hists if not (ref_case in os.path.basename(h))]

    test_hists = list(set(test_hists))
    test_hists.sort()
    logger.debug("_get_all_hist_files returns {} for model {}".format(test_hists, model))
    return test_hists

def _get_latest_hist_files(model, from_dir, file_extensions, suffix="", ref_case=None):
    test_hists = _get_all_hist_files(model, from_dir, file_extensions, suffix=suffix, ref_case=ref_case)
    latest_files = {}
    histlist = []
    for hist in test_hists:
        ext = get_extension(model, hist)
        latest_files[ext] = hist

    for key in latest_files.keys():
        histlist.append(latest_files[key])
    return histlist

def copy(case, suffix):
    """Copy the most recent batch of hist files in a case, adding the given suffix.

    This can allow you to temporarily "save" these files so they won't be blown
    away if you re-run the case.

    case - The case containing the files you want to save
    suffix - The string suffix you want to add to saved files, this can be used to find them later.
    """
    rundir   = case.get_value("RUNDIR")
    ref_case = case.get_value("RUN_REFCASE")
    # Loop over models
    archive = case.get_env("archive")
    comments = "Copying hist files to suffix '{}'\n".format(suffix)
    num_copied = 0
    for model in _iter_model_file_substrs(case):
        comments += "  Copying hist files for model '{}'\n".format(model)
        if model == 'cpl':
            file_extensions = archive.get_hist_file_extensions(archive.get_entry('drv'))
        else:
            file_extensions = archive.get_hist_file_extensions(archive.get_entry(model))
        test_hists = _get_latest_hist_files(model, rundir, file_extensions, ref_case=ref_case)
        num_copied += len(test_hists)
        for test_hist in test_hists:
            new_file = "{}.{}".format(test_hist, suffix)
            if os.path.exists(new_file):
                os.remove(new_file)

            comments += "    Copying '{}' to '{}'\n".format(test_hist, new_file)

            # Need to copy rather than move in case there are some history files
            # that will need to continue to be filled on the next phase; this
            # can be the case for a restart run.
            #
            # (If it weren't for that possibility, a move/rename would be more
            # robust here: The problem with a copy is that there can be
            # confusion after the second run as to which files were created by
            # the first run and which by the second. For example, if the second
            # run fails to output any history files, the test will still pass,
            # because the test system will think that run1's files were output
            # by run2. But we live with that downside for the sake of the reason
            # noted above.)
            safe_copy(test_hist, new_file)

    expect(num_copied > 0, "copy failed: no hist files found in rundir '{}'".format(rundir))

    return comments

def rename_all_hist_files(case, suffix):
    """Renaming all hist files in a case, adding the given suffix.

    case - The case containing the files you want to save
    suffix - The string suffix you want to add to saved files, this can be used to find them later.
    """
    rundir   = case.get_value("RUNDIR")
    ref_case = case.get_value("RUN_REFCASE")
    # Loop over models
    archive = case.get_env("archive")
    comments = "Renaming hist files by adding suffix '{}'\n".format(suffix)
    num_renamed = 0
    for model in _iter_model_file_substrs(case):
        comments += "  Renaming hist files for model '{}'\n".format(model)
        if model == 'cpl':
            file_extensions = archive.get_hist_file_extensions(archive.get_entry('drv'))
        else:
            file_extensions = archive.get_hist_file_extensions(archive.get_entry(model))
        test_hists = _get_all_hist_files(model, rundir, file_extensions, ref_case=ref_case)
        num_renamed += len(test_hists)
        for test_hist in test_hists:
            new_file = "{}.{}".format(test_hist, suffix)
            if os.path.exists(new_file):
                os.remove(new_file)

            comments += "    Renaming '{}' to '{}'\n".format(test_hist, new_file)

           #safe_copy(test_hist, new_file)
           #os.remove(test_hist)
            os.rename(test_hist, new_file)

    expect(num_renamed > 0, "renaming failed: no hist files found in rundir '{}'".format(rundir))

    return comments

def _hists_match(model, hists1, hists2, suffix1="", suffix2=""):
    """
    return (num in set 1 but not 2 , num in set 2 but not 1, matchups)

    >>> hists1 = ['FOO.G.cpl.h1.nc', 'FOO.G.cpl.h2.nc', 'FOO.G.cpl.h3.nc']
    >>> hists2 = ['cpl.h2.nc', 'cpl.h3.nc', 'cpl.h4.nc']
    >>> _hists_match('cpl', hists1, hists2)
    (['FOO.G.cpl.h1.nc'], ['cpl.h4.nc'], [('FOO.G.cpl.h2.nc', 'cpl.h2.nc'), ('FOO.G.cpl.h3.nc', 'cpl.h3.nc')])
    >>> hists1 = ['FOO.G.cpl.h1.nc.SUF1', 'FOO.G.cpl.h2.nc.SUF1', 'FOO.G.cpl.h3.nc.SUF1']
    >>> hists2 = ['cpl.h2.nc.SUF2', 'cpl.h3.nc.SUF2', 'cpl.h4.nc.SUF2']
    >>> _hists_match('cpl', hists1, hists2, 'SUF1', 'SUF2')
    (['FOO.G.cpl.h1.nc.SUF1'], ['cpl.h4.nc.SUF2'], [('FOO.G.cpl.h2.nc.SUF1', 'cpl.h2.nc.SUF2'), ('FOO.G.cpl.h3.nc.SUF1', 'cpl.h3.nc.SUF2')])
    >>> hists1 = ['cam.h0.1850-01-08-00000.nc']
    >>> hists2 = ['cam_0001.h0.1850-01-08-00000.nc','cam_0002.h0.1850-01-08-00000.nc']
    >>> _hists_match('cam', hists1, hists2, '', '')
    ([], [], [('cam.h0.1850-01-08-00000.nc', 'cam_0001.h0.1850-01-08-00000.nc'), ('cam.h0.1850-01-08-00000.nc', 'cam_0002.h0.1850-01-08-00000.nc')])
    >>> hists1 = ['cam_0001.h0.1850-01-08-00000.nc.base','cam_0002.h0.1850-01-08-00000.nc.base']
    >>> hists2 = ['cam_0001.h0.1850-01-08-00000.nc.rest','cam_0002.h0.1850-01-08-00000.nc.rest']
    >>> _hists_match('cam', hists1, hists2, 'base', 'rest')
    ([], [], [('cam_0001.h0.1850-01-08-00000.nc.base', 'cam_0001.h0.1850-01-08-00000.nc.rest'), ('cam_0002.h0.1850-01-08-00000.nc.base', 'cam_0002.h0.1850-01-08-00000.nc.rest')])
    """
    normalized1, normalized2 = [], []
    multi_normalized1, multi_normalized2 = [], []
    multiinst = False

    for hists, suffix, normalized, multi_normalized in [(hists1, suffix1, normalized1, multi_normalized1), (hists2, suffix2, normalized2, multi_normalized2)]:
        for hist in hists:
            normalized_name = hist[hist.rfind(model):]
            if suffix != "":
                expect(normalized_name.endswith(suffix), "How did '{}' not have suffix '{}'".format(hist, suffix))
                normalized_name = normalized_name[:len(normalized_name) - len(suffix) - 1]

            m = re.search("(.+)_[0-9]{4}(.+.nc)",normalized_name)
            if m is not None:
                multiinst = True
                multi_normalized.append(m.group(1)+m.group(2))

            normalized.append(normalized_name)

    set_of_1_not_2 = set(normalized1) - set(normalized2)
    set_of_2_not_1 = set(normalized2) - set(normalized1)

    one_not_two = sorted([hists1[normalized1.index(item)] for item in set_of_1_not_2])
    two_not_one = sorted([hists2[normalized2.index(item)] for item in set_of_2_not_1])

    both = set(normalized1) & set(normalized2)

    match_ups = sorted([ (hists1[normalized1.index(item)], hists2[normalized2.index(item)]) for item in both])

    # Special case - comparing multiinstance to single instance files

    if multi_normalized1 != multi_normalized2:
        # in this case hists1 contains multiinstance hists2 does not
        if set(multi_normalized1) == set(normalized2):
            for idx, norm_hist1 in enumerate(multi_normalized1):
                for idx1, hist2 in enumerate(hists2):
                    norm_hist2 = normalized2[idx1]
                    if norm_hist1 == norm_hist2:
                        match_ups.append((hists1[idx], hist2))
                        if hist2 in two_not_one:
                            two_not_one.remove(hist2)
                        if hists1[idx] in one_not_two:
                            one_not_two.remove(hists1[idx])
        # in this case hists2 contains multiinstance hists1 does not
        if set(multi_normalized2) == set(normalized1):
            for idx, norm_hist2 in enumerate(multi_normalized2):
                for idx1, hist1 in enumerate(hists1):
                    norm_hist1 = normalized1[idx1]
                    if norm_hist2 == norm_hist1:
                        match_ups.append((hist1, hists2[idx]))
                        if hist1 in one_not_two:
                            one_not_two.remove(hist1)
                        if hists2[idx] in two_not_one:
                            two_not_one.remove(hists2[idx])

    if not multiinst:
        expect(len(match_ups) + len(set_of_1_not_2) == len(hists1), "Programming error1")
        expect(len(match_ups) + len(set_of_2_not_1) == len(hists2), "Programming error2")

    return one_not_two, two_not_one, match_ups

def _compare_hists(case, from_dir1, from_dir2, suffix1="", suffix2="", outfile_suffix=""):
    if from_dir1 == from_dir2:
        expect(suffix1 != suffix2, "Comparing files to themselves?")

    casename = case.get_value("CASE")
    testcase = case.get_value("TESTCASE")
    casedir = case.get_value("CASEROOT")
    all_success = True
    num_compared = 0
    comments = "Comparing hists for case '{}' dir1='{}', suffix1='{}',  dir2='{}' suffix2='{}'\n".format(casename, from_dir1, suffix1, from_dir2, suffix2)
    multiinst_driver_compare = False
    archive = case.get_env('archive')
    ref_case = case.get_value("RUN_REFCASE")
    for model in _iter_model_file_substrs(case):
        if model == 'cpl' and suffix2 == 'multiinst':
            multiinst_driver_compare = True
        comments += "  comparing model '{}'\n".format(model)
        if model == 'cpl':
            file_extensions = archive.get_hist_file_extensions(archive.get_entry('drv'))
        else:
            file_extensions = archive.get_hist_file_extensions(archive.get_entry(model))
        hists1 = _get_latest_hist_files(model, from_dir1, file_extensions, suffix=suffix1, ref_case=ref_case)
        hists2 = _get_latest_hist_files(model, from_dir2, file_extensions, suffix=suffix2, ref_case=ref_case)
        if len(hists1) == 0 and len(hists2) == 0:
            comments += "    no hist files found for model {}\n".format(model)
            continue

        one_not_two, two_not_one, match_ups = _hists_match(model, hists1, hists2, suffix1, suffix2)
        for item in one_not_two:
            comments += "    File '{}' {} in '{}' with suffix '{}'\n".format(item, NO_COMPARE, from_dir2, suffix2)
            all_success = False

        for item in two_not_one:
            comments += "    File '{}' {} in '{}' with suffix '{}'\n".format(item, NO_ORIGINAL, from_dir1, suffix1)
            all_success = False

        num_compared += len(match_ups)

        for hist1, hist2 in match_ups:
            success, cprnc_log_file = cprnc(model, hist1, hist2, case, from_dir1,
                                            multiinst_driver_compare=multiinst_driver_compare,
                                            outfile_suffix=outfile_suffix)
            if success:
                comments += "    {} matched {}\n".format(hist1, hist2)
            else:
                comments += "    {} {} {}\n".format(hist1, DIFF_COMMENT, hist2)
                comments += "    cat " + cprnc_log_file + "\n"
                expected_log_file = os.path.join(casedir, os.path.basename(cprnc_log_file))
                if not (os.path.exists(expected_log_file) and filecmp.cmp(cprnc_log_file, expected_log_file)):
                    try:
                        safe_copy(cprnc_log_file, casedir)
                    except (OSError, IOError) as _:
                        logger.warning("Could not copy {} to {}".format(cprnc_log_file, casedir))

                all_success = False
    # PFS test may not have any history files to compare.
    if num_compared == 0 and testcase != "PFS":
        all_success = False
        comments += "Did not compare any hist files! Missing baselines?\n"

    comments += "PASS" if all_success else "FAIL"

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

def cprnc(model, file1, file2, case, rundir, multiinst_driver_compare=False, outfile_suffix=""):
    """
    Run cprnc to compare two individual nc files

    file1 - the full or relative path of the first file
    file2 - the full or relative path of the second file
    case - the case containing the files
    rundir - the rundir for the case
    outfile_suffix - if non-blank, then the output file name ends with this
        suffix (with a '.' added before the given suffix).
        Use None to avoid permissions issues in the case dir.

    returns (True if the files matched, log_name)
    """
    cprnc_exe = case.get_value("CCSM_CPRNC")
    basename = os.path.basename(file1)
    multiinst_regex = re.compile(r'.*%s[^_]*(_[0-9]{4})[.]h.?[.][^.]+?[.]nc' % model)
    mstr = ''
    mstr1 = ''
    mstr2 = ''
    #  If one is a multiinstance file but the other is not add an instance string
    m1 = multiinst_regex.match(file1)
    m2 = multiinst_regex.match(file2)
    if m1 is not None:
        mstr1 = m1.group(1)
    if m2 is not None:
        mstr2 = m2.group(1)
    if mstr1 != mstr2:
        mstr = mstr1+mstr2

    output_filename = os.path.join(rundir, "{}{}.cprnc.out".format(basename, mstr))
    if outfile_suffix:
        output_filename += ".{}".format(outfile_suffix)

    if outfile_suffix is None:
        cpr_stat, out, _ = run_cmd("{} -m {} {}".format(cprnc_exe, file1, file2), combine_output=True)
    else:
        cpr_stat = run_cmd("{} -m {} {}".format(cprnc_exe, file1, file2), combine_output=True, arg_stdout=output_filename)[0]
        with open(output_filename, "r") as fd:
            out = fd.read()

    if multiinst_driver_compare:
        #  In a multiinstance test the cpl hist file will have a different number of
        # dimensions and so cprnc will indicate that the files seem to be DIFFERENT
        # in this case we only want to check that the fields we are able to compare
        # have no differences.
        return (cpr_stat == 0 and " 0 had non-zero differences" in out, output_filename)
    else:
        return (cpr_stat == 0 and "files seem to be IDENTICAL" in out, output_filename)

def compare_baseline(case, baseline_dir=None, outfile_suffix=""):
    """
    compare the current test output to a baseline result

    case - The case containing the hist files to be compared against baselines
    baseline_dir - Optionally, specify a specific baseline dir, otherwise it will be computed from case config
    outfile_suffix - if non-blank, then the cprnc output file name ends with
        this suffix (with a '.' added before the given suffix). if None, no output file saved.

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
            return False, "ERROR {} baseline directory '{}' does not exist".format(TEST_NO_BASELINES_COMMENT,bdir)

    success, comments = _compare_hists(case, rundir, basecmp_dir, outfile_suffix=outfile_suffix)
    if get_model() == "e3sm":
        bless_log = os.path.join(basecmp_dir, BLESS_LOG_NAME)
        if os.path.exists(bless_log):
            last_line = open(bless_log, "r").readlines()[-1]
            comments += "\n  Most recent bless: {}".format(last_line)

    return success, comments

def get_extension(model, filepath):
    """
    For a hist file for the given model, return what we call the "extension"

    model - The component model
    filepath - The path of the hist file

    >>> get_extension("cpl", "cpl.hi.nc")
    'hi'
    >>> get_extension("cpl", "cpl.h.nc")
    'h'
    >>> get_extension("cpl", "cpl.h1.nc.base")
    'h1'
    >>> get_extension("cpl", "TESTRUNDIFF.cpl.hi.0.nc.base")
    'hi'
    >>> get_extension("cpl", "TESTRUNDIFF_Mmpi-serial.f19_g16_rx1.A.melvin_gnu.C.fake_testing_only_20160816_164150-20160816_164240.cpl.h.nc")
    'h'
    >>> get_extension("clm","clm2_0002.h0.1850-01-06-00000.nc")
    '0002.h0'
    >>> get_extension("pop","PFS.f09_g16.B1850.cheyenne_intel.allactive-default.GC.c2_0_b1f2_int.pop.h.ecosys.nday1.0001-01-02.nc")
    'h'
    >>> get_extension("mom", "ga0xnw.mom6.frc._0001_001.nc")
    'frc'
    >>> get_extension("mom", "ga0xnw.mom6.sfc.day._0001_001.nc")
    'sfc.day'
    >>> get_extension("mom", "bixmc5.mom6.prog._0001_01_05_84600.nc")
    'prog'
    >>> get_extension("mom", "bixmc5.mom6.hm._0001_01_03_42300.nc")
    'hm'
    >>> get_extension("mom", "bixmc5.mom6.hmz._0001_01_03_42300.nc")
    'hmz'
    """
    basename = os.path.basename(filepath)
    m = None
    ext_regexes = []

    # First add any model-specific extension regexes; these will be checked before the
    # general regex
    if model == "mom":
        # Need to check 'sfc.day' specially: the embedded '.' messes up the
        # general-purpose regex
        ext_regexes.append(r'sfc\.day')

    # Now add the general-purpose extension regex
    ext_regexes.append(r'\w+')

    for ext_regex in ext_regexes:
        full_regex_str = model+r'\d?_?(\d{4})?\.('+ext_regex+r')[-\w\.]*\.nc\.?'
        full_regex = re.compile(full_regex_str)
        m = full_regex.search(basename)
        if m is not None:
            break

    expect(m is not None, "Failed to get extension for file '{}'".format(filepath))

    if m.group(1) is not None:
        result = m.group(1)+'.'+m.group(2)
    else:
        result = m.group(2)

    return result

def generate_teststatus(testdir, baseline_dir):
    """
    CESM stores it's TestStatus file in baselines. Do not let exceptions
    escape from this function.
    """
    if get_model() == "cesm":
        try:
            if not os.path.isdir(baseline_dir):
                os.makedirs(baseline_dir)

            safe_copy(os.path.join(testdir, TEST_STATUS_FILENAME), baseline_dir)
        except Exception as e:
            logger.warning("Could not copy {} to baselines, {}".format(os.path.join(testdir, TEST_STATUS_FILENAME), str(e)))

def generate_baseline(case, baseline_dir=None, allow_baseline_overwrite=False):
    """
    copy the current test output to baseline result

    case - The case containing the hist files to be copied into baselines
    baseline_dir - Optionally, specify a specific baseline dir, otherwise it will be computed from case config
    allow_baseline_overwrite must be true to generate baselines to an existing directory.

    returns (SUCCESS, comments)
    """
    rundir   = case.get_value("RUNDIR")
    ref_case = case.get_value("RUN_REFCASE")
    if baseline_dir is None:
        baselineroot = case.get_value("BASELINE_ROOT")
        basegen_dir = os.path.join(baselineroot, case.get_value("BASEGEN_CASE"))
    else:
        basegen_dir = baseline_dir
    testcase = case.get_value("CASE")
    archive = case.get_env('archive')

    if not os.path.isdir(basegen_dir):
        os.makedirs(basegen_dir)

    if (os.path.isdir(os.path.join(basegen_dir,testcase)) and
        not allow_baseline_overwrite):
        expect(False, " Cowardly refusing to overwrite existing baseline directory")

    comments = "Generating baselines into '{}'\n".format(basegen_dir)
    num_gen = 0
    for model in _iter_model_file_substrs(case):
        comments += "  generating for model '{}'\n".format(model)
        if model == 'cpl':
            file_extensions = archive.get_hist_file_extensions(archive.get_entry('drv'))
        else:
            file_extensions = archive.get_hist_file_extensions(archive.get_entry(model))
        hists =  _get_latest_hist_files(model, rundir, file_extensions, ref_case=ref_case)
        logger.debug("latest_files: {}".format(hists))
        num_gen += len(hists)
        for hist in hists:
            basename = hist[hist.rfind(model):]
            baseline = os.path.join(basegen_dir, basename)
            if os.path.exists(baseline):
                os.remove(baseline)

            safe_copy(hist, baseline)
            comments += "    generating baseline '{}' from file {}\n".format(baseline, hist)

    # copy latest cpl log to baseline
    # drop the date so that the name is generic
    if case.get_value("COMP_INTERFACE") == "nuopc":
        cplname = "med"
    else:
        cplname = "cpl"
    newestcpllogfile = case.get_latest_cpl_log(coupler_log_path=case.get_value("RUNDIR"), cplname=cplname)
    if newestcpllogfile is None:
        logger.warning("No {}.log file found in directory {}".format(cplname,case.get_value("RUNDIR")))
    else:
        safe_copy(newestcpllogfile, os.path.join(basegen_dir, "{}.log.gz".format(cplname)))
    testname = case.get_value("TESTCASE")
    expect(num_gen > 0 or testname == "PFS", "Could not generate any hist files for case '{}', something is seriously wrong".format(os.path.join(rundir, testcase)))
    #make sure permissions are open in baseline directory
    for root, _, files in os.walk(basegen_dir):
        for name in files:
            try:
                os.chmod(os.path.join(root,name), stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
            except OSError:
                # We tried. Not worth hard failure here.
                pass

    if get_model() == "e3sm":
        bless_log = os.path.join(basegen_dir, BLESS_LOG_NAME)
        with open(bless_log, "a") as fd:
            fd.write("sha:{} date:{}\n".format(get_current_commit(repo=case.get_value("CIMEROOT")),
                                               get_timestamp(timestamp_format="%Y-%m-%d_%H:%M:%S")))

    return True, comments

def get_ts_synopsis(comments):
    r"""
    Reduce case diff comments down to a single line synopsis so that we can put
    something in the TestStatus file. It's expected that the comments provided
    to this function came from compare_baseline, not compare_tests.

    >>> get_ts_synopsis('')
    ''
    >>> get_ts_synopsis('big error')
    'big error'
    >>> get_ts_synopsis('big error\n')
    'big error'
    >>> get_ts_synopsis('stuff\n    File foo had no compare counterpart in bar with suffix baz\nPass\n')
    'ERROR BFAIL some baseline files were missing'
    >>> get_ts_synopsis('stuff\n    File foo did NOT match bar with suffix baz\nPass\n')
    'DIFF'
    >>> get_ts_synopsis('stuff\n    File foo did NOT match bar with suffix baz\n    File foo had no compare counterpart in bar with suffix baz\nPass\n')
    'DIFF'
    >>> get_ts_synopsis('File foo had no compare counterpart in bar with suffix baz\n File foo had no original counterpart in bar with suffix baz\n')
    'DIFF'
    """
    if not comments:
        return ""
    elif "\n" not in comments.strip():
        return comments.strip()
    else:
        has_bfails = False
        has_real_fails = False
        for line in comments.splitlines():
            if NO_COMPARE in line:
                has_bfails = True
            for non_bfail_comment in NON_BFAIL_COMMENT_OPTIONS:
                if non_bfail_comment in line:
                    has_real_fails = True

        if has_real_fails:
            return "DIFF"
        elif has_bfails:
            return "ERROR {} some baseline files were missing".format(TEST_NO_BASELINES_COMMENT)
        else:
            return ""
