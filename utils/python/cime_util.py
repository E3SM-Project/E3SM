"""
Common functions used by cime python scripts
"""

import sys, socket, re, os, time, logging, getpass
from CIME.utils import expect, get_cime_root, get_model, set_model, get_python_libs_location_within_cime
from CIME.XML.Files import Files
from CIME.XML.Machines import Machines

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
    "mira"      : "HiRes_EarthSys",
    "cetus"     : "HiRes_EarthSys",
    "yellowstone" : "P93300606",
}

# Return this error code if the scripts worked but tests failed
TESTS_FAILED_ERR_CODE = 165
###############################################################################
<<<<<<< HEAD
=======
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
>>>>>>> master


###############################################################################

###############################################################################

###############################################################################
###############################################################################


###############################################################################
def probe_machine_name():
###############################################################################
    """
    Use the hostname of your machine to probe for the CIME name for this
    machine.

    >>> probe_machine_name() is not None
    True
    """
    parse_config_machines()
    hostname = socket.gethostname().split(".")[0]
    machine = _MACHINE_INFO.find_machine_from_regex(hostname)
    return machine


###############################################################################

###############################################################################

###############################################################################
def get_batch_system(machine=None):
###############################################################################
    return _MACHINE_INFO.get_value("BATCH_SYSTEM")

###############################################################################
def get_my_queued_jobs():
###############################################################################
    """
    Return a list of job ids for the current user
    """
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
    Moving toward an object oriented model, replace _MACHINE_INFO dict with an object of class Machines
    """
    global _MACHINE_INFO
    if (_MACHINE_INFO is None):
        files = Files()
        config_machines = files.get_resolved_value(files.get_value('MACHINES_SPEC_FILE'))
        _MACHINE_INFO = Machines(config_machines)

###############################################################################
def get_machine_info(items, machine=None):
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

    >>> get_machine_info("CESMSCRATCHROOT", machine="melvin").replace(getpass.getuser(), "USER")
    '/home/USER/acme/scratch'

    >>> get_machine_info("EXEROOT", machine="melvin").replace(getpass.getuser(), "USER")
    '/home/USER/acme/scratch/$CASE/bld'
    """
    parse_config_machines()

    result = []

    if (machine is None):
        machine = _MACHINE_INFO.name if _MACHINE_INFO.name is not None else probe_machine_name()

    expect(machine is not None, "Failed to probe machine. Please provide machine to whatever script you just ran")

    _MACHINE_INFO.set_machine(machine)

    if (type(items) == str):
        result = _MACHINE_INFO.get_value(items)
        if(result is not None):
            result = _MACHINE_INFO.get_resolved_value(result)
    else:
        for item in items:
            thisresult = _MACHINE_INFO.get_value(item)
            if (thisresult is not None):
                thisresult = _MACHINE_INFO.get_resolved_value(thisresult)
            result.append(thisresult)
    return result

###############################################################################
def get_machines():
###############################################################################
    """
    Return all machines defined by the config_machines.xml
    """
    parse_config_machines()

    return _MACHINE_INFO.list_available_machines()

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
    parse_config_machines()
    if (machine is not None):
        _MACHINE_INFO.set_machine(machine)
    batch_system = _MACHINE_INFO.get_node("batch_system")
    return not (batch_system is None or batch_system[0].get('type') == "none")

###############################################################################
<<<<<<< HEAD
=======
def get_utc_timestamp(timestamp_format="%Y%m%d_%H%M%S"):
###############################################################################
    """
    Get a string representing the current UTC time in format: YYMMDD_HHMMSS

    The format can be changed if needed.
    """
    utc_time_tuple = time.gmtime()
    return time.strftime(timestamp_format, utc_time_tuple)

###############################################################################
def setup_standard_logging_options(parser):
###############################################################################
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Print extra information")

    parser.add_argument("-d", "--debug", action="store_true",
                        help="Print debug information (very verbose)")

###############################################################################
def handle_standard_logging_options(args):
###############################################################################
    root_logger = logging.getLogger()

    if (args.verbose == True):
        root_logger.setLevel(logging.INFO)
    if (args.debug == True):
        root_logger.setLevel(logging.DEBUG)
>>>>>>> master
