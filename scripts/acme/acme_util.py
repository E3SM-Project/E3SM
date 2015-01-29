"""
Common functions used by acme python scripts
"""

import sys

_VERBOSE = False

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

###############################################################################
def run_cmd(cmd, ok_to_fail=False, input_str=None, from_dir=None, verbose=None):
###############################################################################
    import subprocess # Not safe to do globally, module not available in older pythons

    verbose_print("RUN: %s" % cmd, verbose)

    if (input_str is not None):
        stdin = subprocess.PIPE
    else:
        stdin = None

    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
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
    expect(sys.version_info.major == major and sys.version_info.minor >= minor,
           "Python %d.%d+ is required, you have %d.%d" %
           (major, minor, sys.version_info.major, sys.version_info.minor))
