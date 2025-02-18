"""
Utilities
"""

import os, sys, re, signal, subprocess, site, time, argparse
from importlib import import_module
from pathlib import Path

###############################################################################
def expect(condition, error_msg, exc_type=SystemExit, error_prefix="ERROR:"):
###############################################################################
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: error2
    """
    if not condition:
        msg = error_prefix + " " + error_msg
        raise exc_type(msg)

###############################################################################
def run_cmd(cmd, input_str=None, from_dir=None, verbose=None, dry_run=False,
            arg_stdout=subprocess.PIPE, arg_stderr=subprocess.PIPE, env=None, combine_output=False):
###############################################################################
    """
    Wrapper around subprocess to make it much more convenient to run shell commands

    >>> run_cmd('ls file_i_hope_doesnt_exist')[0] != 0
    True
    """
    arg_stderr = subprocess.STDOUT if combine_output else arg_stderr
    from_dir = str(from_dir) if from_dir else from_dir

    if verbose:
        print("RUN: {}\nFROM: {}".format(cmd, os.getcwd() if from_dir is None else from_dir))

    if dry_run:
        return 0, "", ""

    if input_str is not None:
        stdin = subprocess.PIPE
        input_str = input_str.encode('utf-8')
    else:
        stdin = None

    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=arg_stdout,
                            stderr=arg_stderr,
                            stdin=stdin,
                            cwd=from_dir,
                            env=env)

    output, errput = proc.communicate(input_str)
    if output is not None:
        try:
            output = output.decode('utf-8', errors='ignore')
            output = output.strip()
        except AttributeError:
            pass
    if errput is not None:
        try:
            errput = errput.decode('utf-8', errors='ignore')
            errput = errput.strip()
        except AttributeError:
            pass

    stat = proc.wait()

    return stat, output, errput

###############################################################################
def run_cmd_no_fail(cmd, input_str=None, from_dir=None, verbose=None, dry_run=False,
                    arg_stdout=subprocess.PIPE, arg_stderr=subprocess.PIPE, env=None, combine_output=False, exc_type=SystemExit):
###############################################################################
    """
    Wrapper around subprocess to make it much more convenient to run shell commands.
    Expects command to work. Just returns output string.

    >>> run_cmd_no_fail('echo foo') == 'foo'
    True
    >>> run_cmd_no_fail('echo THE ERROR >&2; false') # doctest:+ELLIPSIS
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: Command: 'echo THE ERROR >&2; false' failed with error ...

    >>> run_cmd_no_fail('grep foo', input_str='foo') == 'foo'
    True
    >>> run_cmd_no_fail('echo THE ERROR >&2', combine_output=True) == 'THE ERROR'
    True
    """
    stat, output, errput = run_cmd(cmd, input_str=input_str, from_dir=from_dir, verbose=verbose, dry_run=dry_run,
                                   arg_stdout=arg_stdout, arg_stderr=arg_stderr, env=env, combine_output=combine_output)
    if stat != 0:
        # If command produced no errput, put output in the exception since we
        # have nothing else to go on.
        errput = output if not errput else errput
        if errput is None:
            errput = ""

        expect(False, "Command: '{}' failed with error '{}' from dir '{}'".format(cmd, errput, os.getcwd() if from_dir is None else from_dir), exc_type=exc_type)

    return output

###############################################################################
def run_cmd_assert_result(test_obj, cmd, from_dir=None, expect_works=True, env=None, verbose=False):
###############################################################################
    """
    Run a shell command from a unittest.
    """
    from_dir = Path() if from_dir is None else from_dir
    stat, output, errput = run_cmd(cmd, from_dir=from_dir, env=env, verbose=verbose)
    problem = None
    if expect_works and stat != 0:
        problem = "SHOULD HAVE WORKED"
    elif not expect_works and stat == 0:
        problem = "SHOULD NOT HAVE WORKED"

    if problem is not None:
        msg = \
"""
COMMAND: %s
FROM_DIR: %s
%s
OUTPUT: %s
ERRPUT: %s
""" % (cmd, from_dir, problem, output, errput)
        test_obj.assertTrue(False, msg=msg)

    return output

###############################################################################
def check_minimum_python_version(major, minor):
###############################################################################
    """
    Check your python version.

    >>> check_minimum_python_version(sys.version_info[0], sys.version_info[1])
    >>>
    """
    msg = "Python " + str(major) + ", minor version " + str(minor) + " is required, you have " + str(sys.version_info[0]) + "." + str(sys.version_info[1])
    expect(sys.version_info[0] > major or
           (sys.version_info[0] == major and sys.version_info[1] >= minor), msg)

###############################################################################
def convert_to_seconds(time_str):
###############################################################################
    """
    Convert time value in [[HH:]MM:]SS to seconds

    >>> convert_to_seconds("42")
    42
    >>> convert_to_seconds("01:01:01")
    3661
    """
    components = time_str.split(":")
    expect(len(components) < 4, "Unusual time string: '{}'".format(time_str))

    components.reverse()
    result = 0
    for idx, component in enumerate(components):
        result += int(component) * pow(60, idx)

    return result

###############################################################################
def convert_to_babylonian_time(seconds):
###############################################################################
    """
    Convert time value to seconds to HH:MM:SS

    >>> convert_to_babylonian_time(3661)
    '01:01:01'
    """
    hours = int(seconds / 3600)
    seconds %= 3600
    minutes = int(seconds / 60)
    seconds %= 60

    return "{:02d}:{:02d}:{:02d}".format(hours, minutes, seconds)

###############################################################################
def format_time(time_format, input_format, input_time):
###############################################################################
    """
    Converts the string input_time from input_format to time_format
    Valid format specifiers are "%H", "%M", and "%S"
    % signs must be followed by an H, M, or S and then a separator
    Separators can be any string without digits or a % sign
    Each specifier can occur more than once in the input_format,
    but only the first occurence will be used.
    An example of a valid format: "%H:%M:%S"
    Unlike strptime, this does support %H >= 24

    >>> format_time("%H:%M:%S", "%H", "43")
    '43:00:00'
    >>> format_time("%H  %M", "%M,%S", "59,59")
    '0  59'
    >>> format_time("%H, %S", "%H:%M:%S", "2:43:9")
    '2, 09'
    """
    input_fields = input_format.split("%")
    expect(input_fields[0] == input_time[:len(input_fields[0])],
           "Failed to parse the input time; does not match the header string")
    input_time = input_time[len(input_fields[0]):]
    timespec = {"H": None, "M": None, "S": None}
    maxvals = {"M": 60, "S": 60}
    DIGIT_CHECK = re.compile('[^0-9]*')
    # Loop invariants given input follows the specs:
    # field starts with H, M, or S
    # input_time starts with a number corresponding with the start of field
    for field in input_fields[1:]:
        # Find all of the digits at the start of the string
        spec = field[0]
        value_re = re.match(r'\d*', input_time)
        expect(value_re is not None,
               "Failed to parse the input time for the '{}' specifier, expected an integer".format(spec))
        value = value_re.group(0)
        expect(spec in timespec, "Unknown time specifier '" + spec + "'")
        # Don't do anything if the time field is already specified
        if timespec[spec] is None:
            # Verify we aren't exceeding the maximum value
            if spec in maxvals:
                expect(int(value) < maxvals[spec],
                       "Failed to parse the '{}' specifier: A value less than {:d} is expected".format(spec, maxvals[spec]))
            timespec[spec] = value
        input_time = input_time[len(value):]
        # Check for the separator string
        expect(len(re.match(DIGIT_CHECK, field).group(0)) == len(field),
               "Numbers are not permissible in separator strings")
        expect(input_time[:len(field) - 1] == field[1:],
               "The separator string ({}) doesn't match '{}'".format(field[1:], input_time))
        input_time = input_time[len(field) - 1:]
    output_fields = time_format.split("%")
    output_time = output_fields[0]
    # Used when a value isn't given
    min_len_spec = {"H": 1, "M": 2, "S": 2}
    # Loop invariants given input follows the specs:
    # field starts with H, M, or S
    # output_time
    for field in output_fields[1:]:
        expect(field == output_fields[-1] or len(field) > 1,
               "Separator strings are required to properly parse times")
        spec = field[0]
        expect(spec in timespec, "Unknown time specifier '" + spec + "'")
        if timespec[spec] is not None:
            output_time += "0" * (min_len_spec[spec] - len(timespec[spec]))
            output_time += timespec[spec]
        else:
            output_time += "0" * min_len_spec[spec]
        output_time += field[1:]
    return output_time

###############################################################################
def get_timestamp(timestamp_format="%Y%m%d_%H%M%S", utc_time=False):
###############################################################################
    """
    Get a string representing the current UTC time in format: YYYYMMDD_HHMMSS

    The format can be changed if needed.
    """
    if utc_time:
        time_tuple = time.gmtime()
    else:
        time_tuple = time.localtime()
    return time.strftime(timestamp_format, time_tuple)

###############################################################################
class SharedArea(object):
###############################################################################
    """
    Enable 0002 umask within this manager
    """

    def __init__(self, new_perms=0o002):
        self._orig_umask = None
        self._new_perms  = new_perms

    def __enter__(self):
        self._orig_umask = os.umask(self._new_perms)

    def __exit__(self, *_):
        os.umask(self._orig_umask)

###############################################################################
class Timeout(object):
###############################################################################
    """
    A context manager that implements a timeout. By default, it
    will raise exception, but a custon function call can be provided.
    Provided None as seconds makes this class a no-op
    """
    def __init__(self, seconds, action=None):
        self._seconds = seconds
        self._action  = action if action is not None else self._handle_timeout

    def _handle_timeout(self, *_):
        raise RuntimeError("Timeout expired")

    def __enter__(self):
        if self._seconds is not None:
            signal.signal(signal.SIGALRM, self._action)
            signal.alarm(self._seconds)

    def __exit__(self, *_):
        if self._seconds is not None:
            signal.alarm(0)

###############################################################################
def median(items):
###############################################################################
    """
    >>> items = [2.3]
    >>> median(items)
    2.3
    >>> items = [2.3, 8.1, 3.4, 1.5, 11, 3.42321]
    >>> median(items)
    3.4116049999999998
    >>> items = [2.3, 8.1, 3.4, 1.5, 11, 3.42321, -3.1]
    >>> median(items)
    3.4
    """
    if not items:
        return None
    else:
        quotient, remainder = divmod(len(items), 2)
        return sorted(items)[quotient] if remainder else sum(sorted(items)[quotient - 1:quotient + 1]) / 2.

###############################################################################
def ensure_pip():
###############################################################################
    """
    Ensures that pip is available. Notice that we cannot use the _ensure_pylib_impl
    function below, since it would cause circular dependencies. This one has to
    be done by hand.
    """
    try:
        import pip # pylint: disable=unused-import

    except ModuleNotFoundError:
        # Use ensurepip for installing pip
        import ensurepip
        ensurepip.bootstrap(user=True)

        # needed to "rehash" available libs
        site.main() # pylint: disable=no-member

        import pip # pylint: disable=unused-import

###############################################################################
def pip_install_lib(pip_libname):
###############################################################################
    """
    Ask pip to install a version of a package which is >= min_version
    """
    # Installs will use pip, so we need to ensure it is available
    ensure_pip()

    # Note: --trusted-host may not work for ancient versions of python
    #       --upgrade makes sure we get the latest version, even if one is already installed
    stat, _, err = run_cmd("{} -m pip install --upgrade {} --trusted-host files.pythonhosted.org --user".format(sys.executable, pip_libname))
    expect(stat == 0, "Failed to install {}, cannot continue:\n{}".format(pip_libname, err))

    # needed to "rehash" available libs
    site.main() # pylint: disable=no-member

###############################################################################
def package_version_ok(pkg, min_version=None):
###############################################################################
    """
    Checks that the loaded package's version is >= that the minimum required one.
    If no minimum version is passed, then return True
    """
    from pkg_resources import parse_version

    return True if min_version is None else parse_version(pkg.__version__) >= parse_version(min_version)

###############################################################################
def _ensure_pylib_impl(libname, min_version=None, pip_libname=None):
###############################################################################
    """
    Internal method, clients should not call this directly; please use of the
    public ensure_XXX methods. If one does not exist, we will need to evaluate
    if we want to add a new outside dependency.
    """

    install = False
    try:
        pkg = import_module(libname)

        if not package_version_ok(pkg,min_version):
            print("Detected version for package {} is too old: detected {}, required >= {}. Will attempt to upgrade the package locally".format(libname, pkg.__version__,min_version))
            install = True

    except ImportError:
        print("Detected missing package {}. Will attempt to install locally".format(libname))
        pip_libname = pip_libname if pip_libname else libname

        install = True

    if install:
        pip_install_lib(pip_libname)
        pkg = import_module(libname)

    expect(package_version_ok(pkg,min_version),
           "Error! Could not find version {} for package {}.".format(min_version,libname))

# We've accepted these outside dependencies
def ensure_yaml():   _ensure_pylib_impl("yaml", pip_libname="pyyaml",min_version='5.1')
def ensure_pylint(): _ensure_pylib_impl("pylint")
def ensure_psutil(): _ensure_pylib_impl("psutil")
def ensure_netcdf4(): _ensure_pylib_impl("netCDF4")

###############################################################################
class GoodFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawDescriptionHelpFormatter
):
###############################################################################
    """
    We want argument default info to be added but we also want to
    preserve formatting in the description string.
    """
