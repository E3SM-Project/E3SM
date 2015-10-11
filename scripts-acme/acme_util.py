"""
Common functions used by acme python scripts
"""

import sys, socket, re, os, time

_VERBOSE = False

MACHINE_NODENAMES = [
    ("redsky", re.compile(r"redsky-login")),
    ("skybridge", re.compile(r"skybridge-login")),
    ("melvin", re.compile(r"melvin")),
    ("edison", re.compile(r"edison")),
    ("blues", re.compile(r"blogin")),
    ("titan", re.compile(r"titan")),
    ("mira", re.compile(r"mira")),
    ("cetus", re.compile(r"cetus")),
]

# batch-system-name -> ( cmd-to-list-all-jobs-for-user, cmd-to-delete-job )
BATCH_INFO = \
{
    "slurm" : (
        "squeue -o '%i' -h -u",
        "scancel"
    ),
    "pbs" : (
        "qselect -u",
        "qdel"
    ),
}

# TODO: Get these values from ACME XML files instead of duplicating here.

# machine -> defaults (compiler, test_suite, use_batch, project, testroot, baseline_root, proxy, maintain_baselines)
MACHINE_INFO = {
    "redsky"    : (
        "intel",
        "acme_integration",
        True,
        "fy150001",
        "/gscratch/<USER>/acme_scratch",
        "/projects/ccsm/ccsm_baselines",
        "wwwproxy.sandia.gov:80",
        True,
    ),
    "skybridge" : (
        "intel",
        "acme_integration",
        True,
        "fy150001",
        "/gscratch/<USER>/acme_scratch/skybridge",
        "/projects/ccsm/ccsm_baselines",
        "wwwproxy.sandia.gov:80",
        True,
    ),
    "melvin"    : (
        "gnu",
        "acme_developer",
        False,
        "ignore",
        "/home/<USER>/acme/scratch",
        "/home/jgfouca/acme/baselines",
        "sonproxy.sandia.gov:80",
        True,
    ),
    "edison"    : (
        "intel",
        "acme_developer",
        True,
        "acme",
        "/scratch1/scratchdirs/<USER>/acme_scratch",
        "/project/projectdirs/acme/baselines",
        None,
        False,
    ),
    "blues"    : (
        "pgi",
        "acme_developer",
        True,
        "ACME",
        "/lcrc/project/ACME/<USER>/acme_scratch",
        "/lcrc/group/earthscience/acme_baselines",
        None,
        False,
    ),
    "titan"    : (
        "pgi",
        "acme_developer",
        True,
        "cli115",
        "/lustre/atlas/scratch/<USER>/<PROJECT>",
        "/lustre/atlas1/cli900/world-shared/cesm/acme/baselines",
        None,
        False,
    ),
    "mira"     : (
        "ibm",
        "acme_developer",
        True,
        "HiRes_EarthSys",
        "/projects/<PROJECT>/<USER>",
        "/projects/ccsm/ccsm_baselines",
        None,
        False,
    ),
    "cetus"     : (
        "ibm",
        "acme_developer",
        True,
        "HiRes_EarthSys",
        "/projects/<PROJECT>/<USER>",
        "/projects/ccsm/ccsm_baselines",
        None,
        False,
    ),
}

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
    Given an ACME case_id, return it in form TEST.GRID.COMPSET.PLATFORM

    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    >>> normalize_case_id('ERT.ne16_g37.B1850C5.skybridge_intel.G.20151121')
    'ERT.ne16_g37.B1850C5.skybridge_intel'
    """
    sep_count = case_id.count(".")
    expect(sep_count in [3, 5],
           "Case needs to be in form: TEST.GRID.COMPSET.PLATFORM  or  TEST.GRID.COMPSET.PLATFORM.GC.TESTID")
    if (sep_count == 5):
        return ".".join(case_id.split(".")[:-2])
    else:
        return case_id

###############################################################################
def probe_machine_name():
###############################################################################
    """
    Use the hostname of your machine to probe for the ACME name for this
    machine.

    >>> probe_machine_name() is not None
    True
    """
    hostname = socket.gethostname().split(".")[0]

    for machine_name, machine_re in MACHINE_NODENAMES:
        if (machine_re.match(hostname) is not None):
            return machine_name

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
        stat, output, errput = run_cmd("git symbolic-ref HEAD", from_dir=repo, ok_to_fail=True)
        if (stat != 0):
            warning("Couldn't get current git branch, error: '%s'" % errput)
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
    os.path.join(get_cime_location_within_acme(), get_acme_scripts_location_within_cime())

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
def probe_batch_system():
###############################################################################
    import distutils.spawn
    for batch_system, cmds in BATCH_INFO.iteritems():
        exe = cmds[0].split()[0]
        exe_path = distutils.spawn.find_executable(exe)
        if (exe_path is not None):
            return batch_system

    return None

###############################################################################
def get_my_queued_jobs(batch_system=None):
###############################################################################
    """
    Return a list of job ids for the current user
    """
    import getpass
    if (batch_system is None):
        batch_system = probe_batch_system()
    expect(batch_system is not None, "Failed to probe batch system")

    list_cmd = "%s %s" % (BATCH_INFO[batch_system][0], getpass.getuser())
    return run_cmd(list_cmd).split()

###############################################################################
def get_batch_system_info(batch_system=None):
###############################################################################
    """
    Return information on batch system. If no arg provided, probe for batch
    system.

    Info returned as tuple (SUBMIT CMD, DELETE CMD)
    """
    if (batch_system is None):
        batch_system = probe_batch_system()
    expect(batch_system is not None, "Failed to probe batch system")
    expect(batch_system in BATCH_INFO, "No info for batch system '%s'" % batch_system)

    return BATCH_INFO[batch_system]

###############################################################################
def get_machine_info(machine=None, user=None, project=None, raw=False):
###############################################################################
    """
    Return information on machine. If no arg provided, probe for machine.

    Info returned as tuple:
    (compiler, test_suite, use_batch, project, testroot, baseline_root, proxy)
    """
    import getpass
    user = getpass.getuser() if user is None else user

    if (machine is None):
        machine = probe_machine_name()
    expect(machine is not None, "Failed to probe machine")
    expect(machine in MACHINE_INFO, "No info for machine '%s'" % machine)

    machine_info_copy = list(MACHINE_INFO[machine])
    machine_info_copy[3] = project if project is not None else machine_info_copy[3]
    project = machine_info_copy[3]

    if (raw):
        return machine_info_copy
    else:
        return [item.replace("<USER>", user).replace("<PROJECT>", project) if type(item) is str else item for item in machine_info_copy]

###############################################################################
def get_utc_timestamp(timestamp_format="%Y%m%d_%H%M%S"):
###############################################################################
    """
    Get a string representing the current UTC time in format: YYMMDD_HHMMSS

    The format can be changed if needed.
    """
    utc_time_tuple = time.gmtime()
    return time.strftime(timestamp_format, utc_time_tuple)
