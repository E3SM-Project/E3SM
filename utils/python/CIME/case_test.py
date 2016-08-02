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

    try:
        test = find_system_test(testname, case)(case)
        success = test.run()

        test.report()

        if case.get_value("GENERATE_BASELINE"):
            test.generate_baseline()

        if case.get_value("COMPARE_BASELINE"):
            test.compare_baseline()
    except:
        # An uncaught except MUST cause the test to report FAIL
        test.fail_test()
        raise

    return success
