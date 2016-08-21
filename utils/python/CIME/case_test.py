"""
Run a testcase.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, find_system_test, append_status
from CIME.SystemTests.system_tests_common import *

import sys

def case_test(case, testname=None):
    if testname is None:
        testname = case.get_value('TESTCASE')

    expect(testname is not None, "testname argument not resolved")
    logging.warn("Running test for %s" % testname)

    try:
        # The following line can throw exceptions if the testname is
        # not found or the test constructor throws. We need to be
        # sure to leave TestStatus in the appropriate state if that
        # happens.
        test = find_system_test(testname, case)(case)
    except:
        caseroot = case.get_value("CASEROOT")
        with TestStatus(test_dir=caseroot) as ts:
            ts.set_status(RUN_PHASE, TEST_FAIL_STATUS, comments="failed to initialize")
        append_status(sys.exc_info()[1], sfile="TestStatus.log")
        raise

    success = test.run()

    return success
