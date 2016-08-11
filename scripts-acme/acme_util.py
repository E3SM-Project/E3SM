"""
Common functions used by acme python scripts
"""

import sys, socket, re, os, time

_VERBOSE = False

_MACHINE_INFO = None

# batch-system-name -> ( cmd-to-list-all-jobs-for-user, cmd-to-delete-job )
# TODO -> This info should be derived from config_batch.xml
BATCH_INFO = \
{
    "slurm" : (
        "squeue -o '%i' -h -u USER",
        "scancel"
    ),
    "pbs" : (
        "qselect -u USER",
        "qdel"
    ),
    "cobalt" : (
        "qstat -u USER | tail -n+3 | awk '{print $1}'",
        "qdel"
    ),
}

# Don't know if this belongs here longterm
# machine -> default project for nightly process
_MACHINE_PROJECTS = {
    "redsky"    : "fy150001",
    "skybridge" : "fy150001",
    "edison"    : "acme",
    "corip1"    : "acme",
    "blues"     : "ACME",
    "titan"     : "cli115",
    "mira"      : "HiRes_EarthSys_2",
    "cetus"     : "HiRes_EarthSys_2",
}

# Return this error code if the scripts worked but tests failed
TESTS_FAILED_ERR_CODE = 100

###############################################################################
def expect(condition, error_msg):
###############################################################################
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
        raise SystemExit("FAIL: %s" % error_msg)

###############################################################################
def warning(msg):
###############################################################################
    """
    Print a warning to stderr
    """
    print >> sys.stderr, "WARNING:", msg

###############################################################################
def verbose_print(msg, override=None):
###############################################################################
    if ( (_VERBOSE and not override is False) or override):
        print msg

###############################################################################
def set_verbosity(verbose):
###############################################################################
    global _VERBOSE
    _VERBOSE = verbose

_hack=object()
###############################################################################
def run_cmd(cmd, ok_to_fail=False, input_str=None, from_dir=None, verbose=None,
            arg_stdout=_hack, arg_stderr=_hack):
###############################################################################
    """
    Wrapper around subprocess to make it much more convenient to run shell commands

    >>> run_cmd('echo foo')
    'foo'

    >>> run_cmd('ls file_i_hope_doesnt_exist')
    Traceback (most recent call last):
        ...
    SystemExit: FAIL: Command: 'ls file_i_hope_doesnt_exist' failed with error 'ls: cannot access file_i_hope_doesnt_exist: No such file or directory'

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

    verbose_print("RUN: %s" % cmd, verbose)

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

    verbose_print("  stat: %d\n" % stat, verbose)
    verbose_print("  output: %s\n" % output, verbose)
    verbose_print("  errput: %s\n" % errput, verbose)

    if (ok_to_fail):
        return stat, output, errput
    else:
        if (arg_stderr is not None):
            errput = errput if errput is not None else open(arg_stderr.name, "r").read()
            expect(stat == 0, "Command: '%s' failed with error '%s'" % (cmd, errput))
        else:
            expect(stat == 0, "Command: '%s' failed. See terminal output" % cmd)
        return output

###############################################################################
def check_minimum_python_version(major, minor):
###############################################################################
    """
    Check your python version.

    >>> check_minimum_python_version(sys.version_info[0], sys.version_info[1])
    >>>
    """
    expect(sys.version_info[0] == major and sys.version_info[1] >= minor,
           "Python %d.%d+ is required, you have %d.%d" %
           (major, minor, sys.version_info[0], sys.version_info[1]))

###############################################################################
def normalize_case_id(case_id):
###############################################################################
    """
    Given an ACME case_id, return it in form TESTCASE.GRID.COMPSET.PLATFORM

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
           "Case '%s' needs to be in form: TESTCASE.GRID.COMPSET.PLATFORM  or  TESTCASE.GRID.COMPSET.PLATFORM.GC.TESTID" % case_id)
    if (sep_count in [5, 6]):
        return ".".join(case_id.split(".")[:-2])
    else:
        return case_id

###############################################################################
def parse_test_name(test_name):
###############################################################################
    """
    Given an ACME test name TESTCASE[_CASEOPTS].GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]],
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
           "'%s' does not look like an ACME test name, expect TESTCASE.GRID.COMPSET[.MACHINE_COMPILER[.TESTMODS]]" % test_name)

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

###############################################################################
def get_full_test_name(test, machine, compiler, testmod=None):
###############################################################################
    """
    Given an ACME test name, return in form TESTCASE.GRID.COMPSET.MACHINE_COMPILER[.TESTMODS]
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

###############################################################################
def probe_machine_name():
###############################################################################
    """
    Use the hostname of your machine to probe for the ACME name for this
    machine.

    >>> probe_machine_name() is not None
    True
    """
    hostname = socket.getfqdn()

    machines = get_machines()
    for machine in machines:
        regex_str = get_machine_info("NODENAME_REGEX", machine=machine)
        if (regex_str):
            regex = re.compile(regex_str)
            if (regex.match(hostname)):
                return machine

    return None

###############################################################################
def get_current_branch(repo=None):
###############################################################################
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

###############################################################################
def get_current_commit(short=False, repo=None):
###############################################################################
    """
    Return the sha1 of the current HEAD commit

    >>> get_current_commit() is not None
    True
    """
    output = run_cmd("git rev-parse %s HEAD" % ("--short" if short else ""), from_dir=repo)
    return output

###############################################################################
def get_acme_scripts_location_within_cime():
###############################################################################
    """
    From within CIME, return subdirectory where ACME scripts live.
    """
    return "scripts-acme"

###############################################################################
def get_cime_location_within_acme():
###############################################################################
    """
    From within ACME, return subdirectory where CIME lives.
    """
    return "cime"

###############################################################################
def get_acme_scripts_location_within_acme():
###############################################################################
    """
    From within ACME, return subdirectory where ACME scripts live.
    """
    return os.path.join(get_cime_location_within_acme(), get_acme_scripts_location_within_cime())

###############################################################################
def get_cime_root():
###############################################################################
    """
    Return the absolute path to the root of CIME that contains this script

    >>> os.path.isdir(os.path.join(get_cime_root(), get_acme_scripts_location_within_cime()))
    True
    """
    acme_script_absdir = os.path.abspath(os.path.dirname(__file__))
    assert acme_script_absdir.endswith(get_acme_scripts_location_within_cime())
    return os.path.normpath(acme_script_absdir[:len(acme_script_absdir)-len(get_acme_scripts_location_within_cime())])

###############################################################################
def get_acme_root():
###############################################################################
    """
    Return the absolute path to the root of ACME that contains this script

    >>> os.path.isdir(os.path.join(get_acme_root(), '.git'))
    True
    """
    cime_absdir = get_cime_root()
    assert cime_absdir.endswith(get_cime_location_within_acme()), cime_absdir
    return os.path.normpath(cime_absdir[:len(cime_absdir)-len(get_cime_location_within_acme())])

###############################################################################
def get_acme_scripts_root():
###############################################################################
    """
    Get absolute path to acme scripts
    """
    return os.path.abspath(os.path.dirname(__file__))

###############################################################################
def stop_buffering_output():
###############################################################################
    """
    All stdout, stderr will not be buffered after this is called.
    """
    sys.stdout.flush()
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

###############################################################################
def start_buffering_output():
###############################################################################
    """
    All stdout, stderr will be buffered after this is called. This is python's
    default behavior.
    """
    sys.stdout.flush()
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w')

###############################################################################
def match_any(item, re_list):
###############################################################################
    """
    Return true if item matches any regex in re_list
    """
    for regex_str in re_list:
        regex = re.compile(regex_str)
        if (regex.match(item)):
            return True

    return False

###############################################################################
def safe_copy(src_dir, tgt_dir, file_map):
###############################################################################
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

###############################################################################
def find_proc_id(proc_name=None,
                 children_only=False,
                 of_parent=None):
###############################################################################
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

###############################################################################
def get_batch_system(machine=None):
###############################################################################
    return get_machine_info("BATCH_SYSTEM", machine=machine)

###############################################################################
def get_my_queued_jobs():
###############################################################################
    """
    Return a list of job ids for the current user
    """
    import getpass
    batch_system = get_batch_system()
    expect(batch_system is not None, "Failed to probe batch system")

    list_cmd = BATCH_INFO[batch_system][0].replace("USER", getpass.getuser())
    return run_cmd(list_cmd).split()

###############################################################################
def delete_jobs(jobs):
###############################################################################
    """
    Return a list of job ids for the current user.

    Returns (status, output, errput)
    """
    batch_system = get_batch_system()
    expect(batch_system is not None, "Failed to probe batch system")

    del_cmd = "%s %s" % (BATCH_INFO[batch_system][1], " ".join(jobs))
    return run_cmd(del_cmd, ok_to_fail=True, verbose=True)

###############################################################################
def parse_config_machines():
###############################################################################
    """
    """
    global _MACHINE_INFO
    if (_MACHINE_INFO is None):
        _MACHINE_INFO = {}
        import xml.etree.ElementTree as ET
        config_machines_xml = os.path.join(get_cime_root(), "machines-acme", "config_machines.xml")
        tree = ET.parse(config_machines_xml)
        root = tree.getroot()
        expect(root.tag == "config_machines",
               "The given XML file is not a valid list of machine configurations.")

        # Each child of this root is a machine entry.
        for machine in root:
            if (machine.tag == "machine"):
                expect("MACH" in machine.attrib, "Invalid machine entry found for machine '%s'" % machine)
                mach_name = machine.attrib["MACH"]
                expect(mach_name not in _MACHINE_INFO, "Duplicate machine entry '%s'" % mach_name)
                data = {}
                for item in machine:
                    item_data = "" if item.text is None else item.text.strip()
                    if (item.tag in ["COMPILERS", "MPILIBS"]):
                        item_data = [strip_item.strip() for strip_item in item_data.split(",")]
                    elif(item.tag == "batch_system"):
                        item_data = item.attrib["type"]

                    data[item.tag.upper()] = item_data

                _MACHINE_INFO[mach_name] = data
            else:
                warning("Ignoring unrecognized tag: '%s'" % machine.tag)

###############################################################################
def get_machine_info(items, machine=None, user=None, project=None, case=None, raw=False):
###############################################################################
    """
    Return information on machine. If no arg provided, probe for machine.

    If only asked for one thing, will just return the value. If asked for multiple
    things, will return list.

    Info returned as tuple:
    (compiler, test_suite, use_batch, project, testroot, baseline_root, proxy)

    >>> parse_config_machines()

    >>> get_machine_info(["NODENAME_REGEX", "TESTS"], machine="skybridge")
    ['skybridge-login', 'acme_integration']

    >>> get_machine_info("CESMSCRATCHROOT", machine="melvin", user="jenkins")
    '/home/jenkins/acme/scratch'

    >>> get_machine_info("EXEROOT", machine="melvin", user="jenkins", case="Foo")
    '/home/jenkins/acme/scratch/Foo/bld'
    """
    parse_config_machines()

    import getpass
    user = getpass.getuser() if user is None else user

    if (machine is None):
        machine = probe_machine_name()
    expect(machine is not None, "Failed to probe machine. Please provide machine to whatever script you just ran")
    expect(machine in _MACHINE_INFO, "No info for machine '%s'" % machine)

    if (isinstance(items, str)):
        items = [items]

    rv = []
    if (raw):
        for item in items:
            rv.append(_MACHINE_INFO[machine][item])
    else:
        reference_re = re.compile(r'\$(\w+)')
        env_ref_re   = re.compile(r'\$ENV\{(\w+)\}')
        for item in items:
            item_data = _MACHINE_INFO[machine][item] if item in _MACHINE_INFO[machine] else None
            if (isinstance(item_data, str)):
                for m in env_ref_re.finditer(item_data):
                    env_var = m.groups()[0]
                    expect(env_var in os.environ,
                           "Field '%s' for machine '%s' refers to undefined env var '%s'" % (item, machine, env_var))
                    item_data = item_data.replace(m.group(), os.environ[env_var])

                for m in reference_re.finditer(item_data):
                    ref = m.groups()[0]
                    if (ref in _MACHINE_INFO[machine]):
                        item_data = item_data.replace(m.group(), get_machine_info(ref, machine=machine, user=user, project=project))

                item_data = item_data.replace("$USER", user)
                # Need extra logic to handle case where user string was brought in from env ($HOME)
                if (user != getpass.getuser()):
                    item_data = item_data.replace(getpass.getuser(), user, 1)

                if ("$PROJECT" in item_data):
                    project = get_machine_project(machine=machine) if project is None else project
                    expect(project is not None, "Cannot evaluate '%s' without project information" % item)
                    item_data = item_data.replace("$PROJECT", project)

                # $CASE is another special case, it can only be provided by user
                # TODO: It actually comes from one of the env xml files
                if ("$CASE" in item_data):
                    expect(case is not None, "Data for '%s' required case information but none provided" % item)
                    item_data = item_data.replace("$CASE", case)

            rv.append(item_data)

    if (len(rv) == 1):
        return rv[0]
    else:
        return rv

###############################################################################
def get_machines():
###############################################################################
    """
    Return all machines defined by the config_machines.xml
    """
    parse_config_machines()
    return _MACHINE_INFO.keys()

###############################################################################
def get_machine_project(machine=None):
###############################################################################
    """
    Return default project account to use on this machine

    >>> get_machine_project("skybridge")
    'fy150001'
    """
    if ("PROJECT" in os.environ):
        return os.environ["PROJECT"]

    if (machine is None):
        machine = probe_machine_name()

    if (machine in _MACHINE_PROJECTS):
        return _MACHINE_PROJECTS[machine]
    else:
        return None

###############################################################################
def does_machine_have_batch(machine=None):
###############################################################################
    """
    Return if this machine has a batch system

    >>> does_machine_have_batch("melvin")
    False
    >>> does_machine_have_batch("skybridge")
    True
    """
    batch_system = get_machine_info("BATCH_SYSTEM", machine=machine)
    return not (batch_system is None or batch_system == "none")

###############################################################################
def get_utc_timestamp(timestamp_format="%Y%m%d_%H%M%S"):
###############################################################################
    """
    Get a string representing the current UTC time in format: YYMMDD_HHMMSS

    The format can be changed if needed.
    """
    utc_time_tuple = time.gmtime()
    return time.strftime(timestamp_format, utc_time_tuple)
