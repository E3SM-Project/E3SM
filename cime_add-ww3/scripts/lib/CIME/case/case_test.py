"""
Run a testcase.
case_test is a member of class Case from case.py
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, find_system_test, append_testlog, find_proc_id
from CIME.SystemTests.system_tests_common import *

import sys, signal

def _iter_signal_names():
    for signame in [item for item in dir(signal) if item.startswith("SIG") and not item.startswith("SIG_")]:
        yield signame

def _signal_handler(signum, _):
    name = "Unknown"
    for signame in _iter_signal_names():
        if signum == getattr(signal, signame):
            name = signame

    # Terminate children
    proc_ids = find_proc_id(children_only=True)
    for proc_id in proc_ids:
        try:
            os.kill(proc_id, signal.SIGKILL)
        except OSError:
            # If the batch system killed the entire process group, these
            # processes might already be dying
            pass

    # Throw an exception so SystemTest infrastructure can handle this error
    expect(False, "Job killed due to receiving signal {:d} ({})".format(signum, name))

def _set_up_signal_handlers():
    """
    Add handles for all signals that might be used to abort a test

    We need to handle a wide variety due to different implementations of the
    timeout mechanism for different batch systems.
    """
    for signame in ["SIGINT", "SIGTERM", "SIGXCPU", "SIGUSR1", "SIGUSR2"]:
        signum = getattr(signal, signame)
        signal.signal(signum, _signal_handler)

def case_test(self, testname=None, reset=False, skip_pnl=False):
    if testname is None:
        testname = self.get_value('TESTCASE')

    expect(testname is not None, "testname argument not resolved")
    logging.warning("Running test for {}".format(testname))

    _set_up_signal_handlers()

    try:
        # The following line can throw exceptions if the testname is
        # not found or the test constructor throws. We need to be
        # sure to leave TestStatus in the appropriate state if that
        # happens.
        test = find_system_test(testname, self)(self)
    except BaseException:
        caseroot = self.get_value("CASEROOT")
        with TestStatus(test_dir=caseroot) as ts:
            ts.set_status(RUN_PHASE, TEST_FAIL_STATUS, comments="failed to initialize")
        append_testlog(str(sys.exc_info()[1]))
        raise

    if reset:
        logger.info("Reset test to initial conditions and exit")
        # pylint: disable=protected-access
        test._resetup_case(RUN_PHASE)
        return True
    success = test.run(skip_pnl=skip_pnl)

    return success
