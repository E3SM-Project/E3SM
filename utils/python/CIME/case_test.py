"""
Run a testcase.
"""

from CIME.utils import expect, find_system_test
from CIME.SystemTests.system_tests_common import *

def case_test(case, testname=None):
    if testname is None:
        testname = case.get_value('TESTCASE')

    expect(testname is not None, "testname argument not resolved")
    logging.warn("Running test for %s" % testname)

    test = find_system_test(testname, case)(case)
    with test:
        success = test.run()

    return success
