"""
Utilities
"""

import os, sys, re, signal

def expect(condition, error_msg, exc_type=SystemExit, error_prefix="ERROR:"):
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

_hack=object()
def run_cmd(cmd, input_str=None, from_dir=None, verbose=None,
            arg_stdout=_hack, arg_stderr=_hack, env=None, combine_output=False):
    """
    Wrapper around subprocess to make it much more convenient to run shell commands

    >>> run_cmd('ls file_i_hope_doesnt_exist')[0] != 0
    True
    """
    import subprocess # Not safe to do globally, module not available in older pythons

    # Real defaults for these value should be subprocess.PIPE
    if arg_stdout is _hack:
        arg_stdout = subprocess.PIPE

    if arg_stderr is _hack:
        arg_stderr = subprocess.STDOUT if combine_output else subprocess.PIPE

    if verbose:
        print("RUN: {}\nFROM: {}".format(cmd, os.getcwd() if from_dir is None else from_dir))

    if (input_str is not None):
        stdin = subprocess.PIPE
    else:
        stdin = None

    proc = subprocess.Popen(cmd,
                            shell=True,
                            stdout=arg_stdout,
                            stderr=arg_stderr,
                            stdin=stdin,
                            cwd=from_dir,
                            env=env,
                            universal_newlines=True)

    output, errput = proc.communicate(input_str)
    if output is not None:
        try:
            output = output.strip()
        except AttributeError:
            pass
    if errput is not None:
        try:
            errput = errput.strip()
        except AttributeError:
            pass

    stat = proc.wait()

    return stat, output, errput

def run_cmd_no_fail(cmd, input_str=None, from_dir=None, verbose=None,
                    arg_stdout=_hack, arg_stderr=_hack, env=None, combine_output=False, exc_type=SystemExit):
    """
    Wrapper around subprocess to make it much more convenient to run shell commands.
    Expects command to work. Just returns output string.

    >>> run_cmd_no_fail('echo foo') == 'foo'
    True
    >>> run_cmd_no_fail('echo THE ERROR >&2; false') # doctest:+ELLIPSIS
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: Command: 'echo THE ERROR >&2; false' failed with error ...

    >>> run_cmd_no_fail('grep foo', input_str=b'foo') == 'foo'
    True
    >>> run_cmd_no_fail('echo THE ERROR >&2', combine_output=True) == 'THE ERROR'
    True
    """
    stat, output, errput = run_cmd(cmd, input_str, from_dir, verbose, arg_stdout, arg_stderr, env, combine_output)
    if stat != 0:
        # If command produced no errput, put output in the exception since we
        # have nothing else to go on.
        errput = output if not errput else errput
        if errput is None:
            errput = ""

        expect(False, "Command: '{}' failed with error '{}' from dir '{}'".format(cmd, errput, os.getcwd() if from_dir is None else from_dir), exc_type=exc_type)

    return output

def check_minimum_python_version(major, minor):
    """
    Check your python version.

    >>> check_minimum_python_version(sys.version_info[0], sys.version_info[1])
    >>>
    """
    msg = "Python " + str(major) + ", minor version " + str(minor) + " is required, you have " + str(sys.version_info[0]) + "." + str(sys.version_info[1])
    expect(sys.version_info[0] > major or
           (sys.version_info[0] == major and sys.version_info[1] >= minor), msg)

def convert_to_seconds(time_str):
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

def convert_to_babylonian_time(seconds):
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

def format_time(time_format, input_format, input_time):
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

class SharedArea(object):
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

class Timeout(object):
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

def median(items):
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

def get_current_branch(repo=None):
    """
    Return the name of the current branch for a repository

    >>> if "GIT_BRANCH" in os.environ:
    ...     get_current_branch() is not None
    ... else:
    ...     os.environ["GIT_BRANCH"] = "foo"
    ...     get_current_branch() == "foo"
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
        stat, output, _ = run_cmd("git symbolic-ref HEAD", from_dir=repo)
        if (stat != 0):
            return None
        else:
            return output.replace("refs/heads/", "")

def get_current_commit(short=False, repo=None, tag=False, commit="HEAD"):
    """
    Return the sha1 of the current HEAD commit

    >>> get_current_commit() is not None
    True
    """
    if tag:
        rc, output, _ = run_cmd("git describe --tags $(git log -n1 --pretty='%h')", from_dir=repo)
    else:
        rc, output, _ = run_cmd("git rev-parse {} {}".format("--short" if short else "", commit), from_dir=repo)

    return output if rc == 0 else None

def get_current_head(repo=None):
    """
    Return current head, preferring branch name if possible
    """
    branch = get_current_branch(repo=repo)
    if not branch:
        return get_current_commit(repo=repo)
    else:
        return branch

def is_repo_clean(repo=None):
    rc, output, _ = run_cmd("git status --porcelain --untracked-files=no", combine_output=True, from_dir=repo)
    return rc == 0 and output == ""
