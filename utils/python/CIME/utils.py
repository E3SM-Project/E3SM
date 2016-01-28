"""
Common functions used by cime python scripts
"""
import logging
import sys
import os
import time

_MODEL = None
_CIMEROOT = None
# Return this error code if the scripts worked but tests failed
TESTS_FAILED_ERR_CODE = 165

def expect(condition, error_msg):
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: FAIL: error2
    """
    if (not condition):
        raise SystemExit('ERROR: '+error_msg)
    

def get_python_libs_location_within_cime():
    """
    From within CIME, return subdirectory of python libraries
    """
    return os.path.join("utils", "python")

def get_cime_root():
    """
    Return the absolute path to the root of CIME that contains this script

    >>> os.path.isdir(os.path.join(get_cime_root(), get_acme_scripts_location_within_cime()))
    True
    """
    global _CIMEROOT
    if (_CIMEROOT is None):
        try:
            _CIMEROOT = os.environ["CIMEROOT"]
        except KeyError:
            script_absdir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
            assert script_absdir.endswith(get_python_libs_location_within_cime()), script_absdir
            _CIMEROOT = os.path.abspath(os.path.join(script_absdir,"..",".."))
    logging.info( "CIMEROOT is " + _CIMEROOT)
    return _CIMEROOT


def set_model(model):
    global _MODEL
    _MODEL = model


def get_model():
    global _MODEL
    if (_MODEL is None):
        try:
            _MODEL = os.environ["CIME_MODEL"]
        except KeyError:
            modelroot = os.path.join(get_cime_root(), "cime_config")
            models = os.listdir(modelroot)
            msg = "Environment variable CIME_MODEL must be set to one of: "
            msg += ", ".join([model for model in models if os.path.isdir(os.path.join(modelroot,model)) and model != "xml_schemas"])
            expect(False, msg)

    return _MODEL

_hack=object()

def run_cmd(cmd, ok_to_fail=False, input_str=None, from_dir=None, verbose=None,
            arg_stdout=_hack, arg_stderr=_hack):
    """
    Wrapper around subprocess to make it much more convenient to run shell commands

    >>> run_cmd('echo foo')
    'foo'

    >>> run_cmd('ls file_i_hope_doesnt_exist')
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: Command: 'ls file_i_hope_doesnt_exist' failed with error 'ls: cannot access file_i_hope_doesnt_exist: No such file or directory'

    >>> run_cmd('ls file_i_hope_doesnt_exist', ok_to_fail=True)[0] != 0
    True

    >>> run_cmd('grep foo', input_str='foo')
    'foo'
    """
    import subprocess # Not safe to do globally, module not available in older pythons

    # Real defaults for these value should be subprocess.PIPE
    if (arg_stdout is _hack):
        arg_stdout = subprocess.PIPE
    if (arg_stderr is _hack):
        arg_stderr = subprocess.PIPE

    logging.info("RUN: %s" % cmd)

    if (input_str is not None):
        stdin = subprocess.PIPE
    else:
        stdin = None

    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=arg_stdout,
                            stderr=arg_stderr,
                            stdin=stdin,
                            cwd=from_dir)
    
    output, errput = proc.communicate(input_str)
    output = output.strip() if output is not None else output
    errput = errput.strip() if errput is not None else errput
    stat = proc.wait()

    logging.info("  stat: %d\n" % stat)
    logging.info("  output: %s\n" % output)
    logging.info("  errput: %s\n" % errput)

    if (ok_to_fail):
        return stat, output, errput
    else:
        if (arg_stderr is not None):
            errput = errput if errput is not None else open(arg_stderr.name, "r").read()
            expect(stat == 0, "Command: '%s' failed with error '%s'" % (cmd, errput))
        else:
            expect(stat == 0, "Command: '%s' failed. See terminal output" % cmd)
        return output

def check_minimum_python_version(major, minor):
    """
    Check your python version.

    >>> check_minimum_python_version(sys.version_info[0], sys.version_info[1])
    >>>
    """
    expect(sys.version_info[0] == major and sys.version_info[1] >= minor,
           "Python %d.%d+ is required, you have %d.%d" %
           (major, minor, sys.version_info[0], sys.version_info[1]))

def normalize_case_id(case_id):
    """
    Given a CIME case_id, return it in form TESTCASE.GRID.COMPSET.PLATFORM

    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel.G.20151121')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    """
    sep_count = case_id.count(".")
    expect(sep_count in [3, 5],
           "Case needs to be in form: TESTCASE.GRID.COMPSET.PLATFORM  or  TESTCASE.GRID.COMPSET.PLATFORM.GC.TESTID")
    if (sep_count == 5):
        return ".".join(case_id.split(".")[:-2])
    else:
        return case_id

def parse_test_name(test_name):
    """
    Given a CIME test name TESTCASE[_CASEOPTS].GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]],
    return each component of the testname with machine and compiler split

    >>> parse_test_name('ERS.fe12_123.JGF')
    ['ERS', None, 'fe12_123', 'JGF', None, None, None]
    >>> parse_test_name('ERS_D.fe12_123.JGF')
    ['ERS', ['D'], 'fe12_123', 'JGF', None, None, None]
    >>> parse_test_name('ERS_D_P1.fe12_123.JGF')
    ['ERS', ['D', 'P1'], 'fe12_123', 'JGF', None, None, None]
    >>> parse_test_name('ERS.fe12_123.JGF.machine_compiler')
    ['ERS', None, 'fe12_123', 'JGF', 'machine', 'compiler', None]
    >>> parse_test_name('ERS.fe12_123.JGF.machine_compiler.test-mods')
    ['ERS', None, 'fe12_123', 'JGF', 'machine', 'compiler', 'test/mods']
    """
    rv = [None] * 6
    num_dots = test_name.count(".")
    expect(num_dots >= 2 and num_dots <= 4,
           "'%s' does not look like a CIME test name, expect TESTCASE.GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]]" % test_name)

    rv[0:num_dots+1] = test_name.split(".")
    testcase_field_underscores = rv[0].count("_")
    rv.insert(1, None) # Make room for caseopts
    if (testcase_field_underscores > 0):
        full_str = rv[0]
        rv[0]    = full_str.split("_")[0]
        rv[1]    = full_str.split("_")[1:]

    if (num_dots >= 3):
        expect(rv[4].count("_") == 1,
               "Expected 4th item of '%s' ('%s') to be in form machine_compiler" % (test_name, rv[4]))
        rv[4:5] = rv[4].split("_")
        rv.pop()

    if (rv[-1] is not None):
        rv[-1] = rv[-1].replace("-", "/")

    return rv

def get_full_test_name(test, machine, compiler, testmod=None):
    """
    Given a CIME test name, return in form TESTCASE.GRID.COMPSET.MACHINE_COMPILER[.TESTMODS]
    Use the machine, compiler, and testmod provided to fill out the name if needed

    >>> get_full_test_name("ERS.ne16_fe16.JGF", "melvin", "gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu'
    >>> get_full_test_name("ERS.ne16_fe16.JGF.melvin_gnu.mods", "melvin", "gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu.mods'
    >>> get_full_test_name("ERS.ne16_fe16.JGF", "melvin", "gnu", "mods/test")
    'ERS.ne16_fe16.JGF.melvin_gnu.mods-test'
    """
    if (test.count(".") == 2):
        return "%s.%s_%s%s" % (test, machine, compiler, "" if testmod is None else ".%s" % testmod.replace("/", "-"))
    else:
        _, _, _, _, test_machine, test_compiler, test_testmod = parse_test_name(test)
        expect(machine == test_machine,
               "Found testname/machine mismatch, test is '%s', your current machine is '%s'" % (test, machine))
        expect(compiler == test_compiler,
               "Found testname/compiler mismatch, test is '%s', your current compiler is '%s'" % (test, compiler))
        if (test_testmod is None):
            return "%s%s" % (test, "" if testmod is None else ".%s" % testmod.replace("/", "-"))
        else:
            return test

def get_current_branch(repo=None):
    """
    Return the name of the current branch for a repository

    >>> get_current_branch() is not None
    True
    """
    if ("GIT_BRANCH" in os.environ):
        # This approach works better for Jenkins jobs because the Jenkins
        # git plugin does not use local tracking branches, it just checks out
        # to a commit
        branch = os.environ["GIT_BRANCH"]
        if (branch.startswith("origin/")):
            branch = branch.replace("origin/", "", 1)
        return branch
    else:
        stat, output, _ = run_cmd("git symbolic-ref HEAD", from_dir=repo, ok_to_fail=True)
        if (stat != 0):
            return None
        else:
            return output.replace("refs/heads/", "")

def get_current_commit(short=False, repo=None):
    """
    Return the sha1 of the current HEAD commit

    >>> get_current_commit() is not None
    True
    """
    output = run_cmd("git rev-parse %s HEAD" % ("--short" if short else ""), from_dir=repo)
    return output

def get_acme_scripts_location_within_cime():
    """
    From within CIME, return subdirectory where ACME scripts live.
    """
    return "scripts-python"

def get_cime_location_within_acme():
    """
    From within ACME, return subdirectory where CIME lives.
    """
    return "cime"

def get_model_config_location_within_cime(model=get_model()):
    return os.path.join("cime_config", model)

def get_acme_root():
    """
    Return the absolute path to the root of ACME that contains this script
    """
    cime_absdir = get_cime_root()
    assert cime_absdir.endswith(get_cime_location_within_acme()), cime_absdir
    return os.path.normpath(cime_absdir[:len(cime_absdir)-len(get_cime_location_within_acme())])

def get_acme_scripts_root():
    """
    Get absolute path to acme scripts

    >>> os.path.isdir(get_acme_scripts_root())
    True
    """
    return os.path.join(get_cime_root(), get_acme_scripts_location_within_cime())

def get_python_libs_root():
    """
    Get absolute path to acme scripts

    >>> os.path.isdir(get_python_libs_root())
    True
    """
    return os.path.join(get_cime_root(), get_python_libs_location_within_cime())

def get_model_config_root(model=get_model()):
    """
    Get absolute path to acme config area"

    >>> os.path.isdir(get_model_config_root())
    True
    """
    return os.path.join(get_cime_root(), get_model_config_location_within_cime(model))

def stop_buffering_output():
    """
    All stdout, stderr will not be buffered after this is called.
    """
    sys.stdout.flush()
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

def start_buffering_output():
    """
    All stdout, stderr will be buffered after this is called. This is python's
    default behavior.
    """
    sys.stdout.flush()
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w')

def match_any(item, re_list):
    """
    Return true if item matches any regex in re_list
    """
    for regex_str in re_list:
        regex = re.compile(regex_str)
        if (regex.match(item)):
            return True

    return False

def safe_copy(src_dir, tgt_dir, file_map):
    """
    Copies a set of files from one dir to another. Works even if overwriting a
    read-only file. Files can be relative paths and the relative path will be
    matched on the tgt side.
    """
    import shutil
    for src_file, tgt_file in file_map:
        full_tgt = os.path.join(tgt_dir, tgt_file)
        full_src = src_file if os.path.isabs(src_file) else os.path.join(src_dir, src_file)
        expect(os.path.isfile(full_src), "Source dir '%s' missing file '%s'" % (src_dir, src_file))
        if (os.path.isfile(full_tgt)):
            os.remove(full_tgt)
        shutil.copy2(full_src, full_tgt)

def find_proc_id(proc_name=None,
                 children_only=False,
                 of_parent=None):
    """
    Children implies recursive.
    """
    expect(proc_name is not None or children_only,
           "Must provide proc_name if not searching for children")
    expect(not (of_parent is not None and not children_only),
           "of_parent only used with children_only")

    parent = of_parent if of_parent is not None else os.getpid()

    pgrep_cmd = "pgrep %s %s" % (proc_name if proc_name is not None else "",
                                 "-P %d" % parent if children_only else "")
    stat, output, errput = run_cmd(pgrep_cmd, ok_to_fail=True)
    expect(stat in [0, 1], "pgrep failed with error: '%s'" % errput)

    rv = set([int(item.strip()) for item in output.splitlines()])
    if (children_only):
        pgrep_cmd = "pgrep -P %s" % parent
        stat, output, errput = run_cmd(pgrep_cmd, ok_to_fail=True)
        expect(stat in [0, 1], "pgrep failed with error: '%s'" % errput)

        for child in output.splitlines():
            rv = rv.union(set(find_proc_id(proc_name, children_only, int(child.strip()))))

    return list(rv)

def get_utc_timestamp(timestamp_format="%Y%m%d_%H%M%S"):
    """
    Get a string representing the current UTC time in format: YYMMDD_HHMMSS

    The format can be changed if needed.
    """
    utc_time_tuple = time.gmtime()
    return time.strftime(timestamp_format, utc_time_tuple)

