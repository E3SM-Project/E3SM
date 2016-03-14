"""
Common functions used by cime python scripts
Warning: you cannot use CIME Classes in this module as it causes circular dependancies
"""
import logging
import logging.config
import sys
import os
import time
import re
if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

# Return this error code if the scripts worked but tests failed
TESTS_FAILED_ERR_CODE = 165
logger = logging.getLogger(__name__)

def expect(condition, error_msg):
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: error2
    """
    if (not condition):
        # Uncomment these to bring up a debugger when an expect fails
        #import pdb
        #pdb.set_trace()
        raise SystemExit("ERROR: %s" % error_msg)

# Should only be called from get_cime_config()
def _read_cime_config_file():
    """
    READ the config file in ~/.cime, this file may contain
    [main]
    CIMEROOT=/path/to/cime
    CIME_MODEL=acme,cesm
    PROJECT=someprojectnumber
    """
    cime_config_file = os.path.abspath(os.path.join(os.path.expanduser("~"),
                                                  ".cime","config"))
    cime_config = config_parser()
    if(os.path.isfile(cime_config_file)):
        cime_config.read(cime_config_file)
    else:
        logger.debug("File %s not found" % cime_config_file)
        cime_config.add_section('main')

    return cime_config

_CIMECONFIG = None
def get_cime_config():
    global _CIMECONFIG
    if (not _CIMECONFIG):
        _CIMECONFIG = _read_cime_config_file()

    return _CIMECONFIG

def get_python_libs_location_within_cime():
    """
    From within CIME, return subdirectory of python libraries
    """
    return os.path.join("utils", "python")

def get_cime_root():
    """
    Return the absolute path to the root of CIME that contains this script

    >>> os.path.isdir(os.path.join(get_cime_root(), get_scripts_location_within_cime()))
    True
    """
    cime_config = get_cime_config()
    if(cime_config.has_option('main','CIMEROOT')):
        cimeroot = cime_config.get('main','CIMEROOT')
    else:
        try:
            cimeroot = os.environ["CIMEROOT"]
        except KeyError:
            script_absdir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
            assert script_absdir.endswith(get_python_libs_location_within_cime()), script_absdir
            cimeroot = os.path.abspath(os.path.join(script_absdir,"..",".."))
        cime_config.set('main','CIMEROOT',cimeroot)
    logger.debug( "CIMEROOT is " + cimeroot)
    return cimeroot

def set_model(model):
    """
    Set the model to be used in this session
    """
    cime_config = get_cime_config()
    cime_config.set('main','MODEL',model)

def get_model():
    """
    Get the currently configured model value

    >>> set_model('rocky')
    >>> print get_model()
    rocky
    """
    cime_config = get_cime_config()
    model = None
    if (cime_config.has_option('main','MODEL')):
        model = cime_config.get('main','MODEL')
    else:
        model = os.environ.get("CIME_MODEL")
        if (model is not None):
            set_model(model)
        else:
            modelroot = os.path.join(get_cime_root(), "cime_config")
            models = os.listdir(modelroot)
            msg = "Environment variable CIME_MODEL must be set to one of: "
            msg += ", ".join([model for model in models
                              if os.path.isdir(os.path.join(modelroot,model))
                              and model != "xml_schemas"])
            expect(False, msg)

    return model

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

    if (verbose or logger.level == logging.DEBUG):
        logger.info("RUN: %s" % cmd)

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

    if (verbose or logger.level == logging.DEBUG):
        if stat != 0:
            logger.info("  stat: %d\n" % stat)
        if output:
            logger.info("  output: %s\n" % output)
        if errput:
            logger.info("  errput: %s\n" % errput)

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
    Given a case_id, return it in form TESTCASE.GRID.COMPSET.PLATFORM

    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel.test-mod')
    'ERT.ne16_g37.B1850C5.skybridge_intel.test-mod'
    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel.G.20151121')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel.test-mod.G.20151121')
    'ERT.ne16_g37.B1850C5.skybridge_intel.test-mod'
    """
    sep_count = case_id.count(".")
    expect(sep_count >= 3 and sep_count <= 6,
           "Case '%s' needs to be in form: TESTCASE.GRID.COMPSET.PLATFORM[.TESTMOD]  or  TESTCASE.GRID.COMPSET.PLATFORM[.TESTMOD].GC.TESTID" % case_id)
    if (sep_count in [5, 6]):
        return ".".join(case_id.split(".")[:-2])
    else:
        return case_id

def parse_test_name(test_name):
    """
    Given a CIME test name TESTCASE[_CASEOPTS].GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]],
    return each component of the testname with machine and compiler split

    >>> parse_test_name('ERS')
    ['ERS', None, None, None, None, None, None]
    >>> parse_test_name('ERS.fe12_123')
    ['ERS', None, 'fe12_123', None, None, None, None]
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
    expect(num_dots <= 4,
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

def get_full_test_name(partial_test, grid=None, compset=None, machine=None, compiler=None, testmod=None):
    """
    Given a partial CIME test name, return in form TESTCASE.GRID.COMPSET.MACHINE_COMPILER[.TESTMODS]
    Use the additional args to fill out the name if needed

    >>> get_full_test_name("ERS", grid="ne16_fe16", compset="JGF", machine="melvin", compiler="gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu'
    >>> get_full_test_name("ERS.ne16_fe16", compset="JGF", machine="melvin", compiler="gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu'
    >>> get_full_test_name("ERS.ne16_fe16.JGF", machine="melvin", compiler="gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu'
    >>> get_full_test_name("ERS.ne16_fe16.JGF.melvin_gnu.mods", machine="melvin", compiler="gnu")
    'ERS.ne16_fe16.JGF.melvin_gnu.mods'
    >>> get_full_test_name("ERS.ne16_fe16.JGF", machine="melvin", compiler="gnu", testmod="mods/test")
    'ERS.ne16_fe16.JGF.melvin_gnu.mods-test'
    """
    _, _, partial_grid, partial_compset, partial_machine, partial_compiler, partial_testmod = parse_test_name(partial_test)

    required_fields = [
        (partial_grid, grid, "grid"),
        (partial_compset, compset, "compset"),
        (partial_machine, machine, "machine"),
        (partial_compiler, compiler, "compiler"),
        ]

    result = partial_test

    for partial_val, arg_val, name in required_fields:
        if (partial_val is None):
            # Add to result based on args
            expect(arg_val is not None,
                   "Could not fill-out test name, partial string '%s' had no %s information and you did not provide any" % (partial_test, name))
            result = "%s%s%s" % (result, "_" if name == "compiler" else ".", arg_val)
        elif (arg_val is not None):
            expect(arg_val == partial_val,
                   "Mismatch in field %s, partial string '%s' indicated it should be '%s' but you provided '%s'" % (name, partial_test, partial_val, arg_val))

    if (partial_testmod is None):
        if (testmod is None):
            # No testmod for this test and that's OK
            pass
        else:
            result += ".%s" % testmod.replace("/", "-")
    elif (testmod is not None):
        expect(arg_val == partial_val,
               "Mismatch in field testmod, partial string '%s' indicated it should be '%s' but you provided '%s'" % (partial_test, partial_testmod, testmod))

    return result

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

def get_scripts_location_within_cime():
    """
    From within CIME, return subdirectory where scripts live.
    """
    return "scripts"

def get_cime_location_within_acme():
    """
    From within acme, return subdirectory where CIME lives.
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

def get_scripts_root():
    """
    Get absolute path to scripts

    >>> os.path.isdir(get_scripts_root())
    True
    """
    return os.path.join(get_cime_root(), get_scripts_location_within_cime())

def get_python_libs_root():
    """
    Get absolute path to scripts

    >>> os.path.isdir(get_python_libs_root())
    True
    """
    return os.path.join(get_cime_root(), get_python_libs_location_within_cime())

def get_model_config_root(model=get_model()):
    """
    Get absolute path to model config area"

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

def get_project():
    """
    Hierarchy for choosing PROJECT:
    1. Environment variable PROJECT
    2 environment variable ACCOUNT   (this is for backward compatibility)
    3. File $HOME/.cime/config   (this is new)
    4 File $HOME/.cesm_proj  (again - backward compatibility)
    5 config_machines.xml
    """
    project = os.environ.get("PROJECT")
    if (project is not None):
        logger.info("Using project from env PROJECT "+project)
        return project
    project = os.environ.get("ACCOUNT")
    if (project is not None):
        logger.info("Using project from env ACCOUNT "+project)
        return project

    cime_config = get_cime_config()
    if (cime_config.has_option('main','PROJECT')):
        project = cime_config.get('main','PROJECT')
        if (project is not None):
            logger.info("Using project from .cime/config "+project)
            return project
        projectfile = os.path.abspath(os.path.join(os.path.expanduser("~"),
                                                   ".cesm_proj"))
        if (os.path.isfile(projectfile)):
            with open(projectfile,'r') as myfile:
                project = myfile.read()
                cime_config.set('main','PROJECT',project)

    return project

def setup_standard_logging_options(parser):
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print extra information")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="Print debug information (very verbose)")
    parser.add_argument("-s", "--silent", action="store_true",
                        help="Print only warnings and error messages")

def handle_standard_logging_options(args):
    LOGGING = {
        'version': 1,
        'disable_existing_loggers': False,
        'propagate': True,
        'formatters': {
            'verbose': {
                'format': '%(name)-12s %(levelname)-8s %(message)s',
                },
            'debug': {
                'format': '%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                },
            'default':{
                'format': '%(message)s',
                },
            },
        'handlers': {
            'console1': {
                'class': 'logging.StreamHandler',
                'formatter': 'default',
                'level' : logging.INFO,
                },
            'file1' : {
                'class': 'logging.FileHandler',
                'mode':'w',
                'formatter':'debug',
                'filename' : '%s.log'%sys.argv[0],
                'level' : logging.DEBUG,
                },
            },
        'root': {
            'handlers': ['console1'],
            }
        }
    root_logger=logging.getLogger()
    # DEBUG trumps INFO trumps Silent (WARN)
    if (args.debug == True):
        LOGGING['root']['handlers'].append('file1')
        LOGGING['handlers']['console1']['formatter'] = 'debug'
        root_logger.setLevel(logging.DEBUG)
        logger.warn("Log level set to DEBUG")
    elif (args.verbose == True):
        LOGGING['handlers']['console1']['formatter'] = 'verbose'
        logger.warn("Log level set to VERBOSE")
    elif (args.silent == True):
        root_logger.setLevel(logging.WARN)

    logging.config.dictConfig(LOGGING)

def get_logging_options():
    """
    Use to pass same logging options as was used for current
    executable to subprocesses.
    """
    root_logger = logging.getLogger()

    if (root_logger.level == logging.INFO):
        return "--verbose"
    elif (root_logger.level == logging.DEBUG):
        return "--debug"
    elif (root_logger.level == logging.WARN):
        return "--silent"
    else:
        return ""

def convert_to_type(value, type_str, vid=""):
    """
    Convert value from string to another type.
    vid is only for generating better error messages.
    """
    if value is not None:

        if type_str == "char":
            pass

        elif type_str == "integer":
            try:
                value = int(value)
            except ValueError:
                expect(False, "Entry %s was listed as type int but value '%s' is not valid int" % (vid, value))

        elif type_str == "logical":
            expect(value in ["TRUE", "FALSE","true","false"],
                   "Entry %s was listed as type logical but had val '%s' instead of TRUE or FALSE" % (vid, value))
            value = value == "TRUE" or value == "true"

        elif type_str == "real":
            try:
                value = float(value)
            except:
                expect(False, "Entry %s was listed as type real but value '%s' is not valid real" % (vid, value))

        else:
            expect(False, "Unknown type '%s'" % type_str)

    return value

def convert_to_string(value, type_str, vid=""):
    """
    Convert value back to string.
    vid is only for generating better error messages.
    """
    if value is not None:
        if type_str == "char":
            expect(type(value) is str, "Wrong type for entry id '%s'" % vid)
        elif type_str == "integer":
            expect(type(value) is int, "Wrong type for entry id '%s'" % vid)
            value = str(value)
        elif type_str == "logical":
            expect(type(value) is bool, "Wrong type for entry id '%s'" % vid)
            value = "TRUE" if value else "FALSE"
        elif type_str == "real":
            expect(type(value) is float, "Wrong type for entry id '%s'" % vid)
            value = str(value)
        else:
            expect(False, "Unknown type '%s'" % type_str)

    return value
