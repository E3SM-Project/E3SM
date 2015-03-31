"""
Common functions used by acme python scripts
"""

import sys, socket, re, os

_VERBOSE = False

MACHINE_NODENAMES = [
    ("redsky", re.compile(r"redsky-login")),
    ("skybridge", re.compile(r"skybridge-login")),
    ("melvin", re.compile(r"melvin"))
]

###############################################################################
def expect(condition, error_msg):
###############################################################################
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.
    """
    if (not condition):
        raise SystemExit(error_msg)

###############################################################################
def warning(msg):
###############################################################################
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
    output, errput = proc.communicate()
    stat = proc.wait()

    verbose_print("  stat: %d\n" % stat, verbose)
    verbose_print("  output: %s\n" % output, verbose)
    verbose_print("  errput: %s\n" % errput, verbose)

    if (ok_to_fail):
        return stat, output, errput
    else:
        expect(stat == 0, "Command: '%s' failed with error '%s'" % (cmd, errput))
        return output

###############################################################################
def check_minimum_python_version(major, minor):
###############################################################################
    expect(sys.version_info[0] == major and sys.version_info[1] >= minor,
           "Python %d.%d+ is required, you have %d.%d" %
           (major, minor, sys.version_info[0], sys.version_info[1]))

###############################################################################
def normalize_case_id(case_id):
###############################################################################
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
        output = run_cmd("git symbolic-ref HEAD", from_dir=repo)
        return output.replace("refs/heads/", "").strip()
