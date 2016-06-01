"""
Run a testcase.
"""

from CIME.utils import expect
from CIME.case import Case
from CIME.SystemTests.sms import SMS
from CIME.SystemTests.seq import SEQ
from CIME.SystemTests.ers import ERS
from CIME.SystemTests.nck import NCK
from CIME.SystemTests.err import ERR
from CIME.SystemTests.eri import ERI
from CIME.SystemTests.system_tests_common import *

def case_test(case, testname=None):
    if testname is None:
        testname = case.get_value('TESTCASE')

    expect(testname is not None, "testname argument not resolved")
    logging.warn("Running test for %s" % testname)

    test = globals()[testname](case)

    success = test.run()

    test.report()

    if case.get_value("GENERATE_BASELINE"):
        test.generate_baseline()

    if case.get_value("COMPARE_BASELINE"):
        test.compare_baseline()

    return success
