"""
Run a testcase.
"""

from CIME.utils import expect

from CIME.SystemTests.eri import ERI # pylint: disable=unused-import
from CIME.SystemTests.err import ERR # pylint: disable=unused-import
from CIME.SystemTests.erp import ERP # pylint: disable=unused-import
from CIME.SystemTests.ers import ERS # pylint: disable=unused-import
from CIME.SystemTests.ert import ERT # pylint: disable=unused-import
from CIME.SystemTests.lii import LII # pylint: disable=unused-import
from CIME.SystemTests.nck import NCK # pylint: disable=unused-import
from CIME.SystemTests.pea import PEA # pylint: disable=unused-import
from CIME.SystemTests.pem import PEM # pylint: disable=unused-import
from CIME.SystemTests.pet import PET # pylint: disable=unused-import
from CIME.SystemTests.pfs import PFS # pylint: disable=unused-import
from CIME.SystemTests.sms import SMS # pylint: disable=unused-import
from CIME.SystemTests.seq import SEQ # pylint: disable=unused-import
from CIME.SystemTests.ssp import SSP # pylint: disable=unused-import

from CIME.SystemTests.system_tests_common import *

def case_test(case, testname=None):
    if testname is None:
        testname = case.get_value('TESTCASE')

    expect(testname is not None, "testname argument not resolved")
    logging.warn("Running test for %s" % testname)

    try:
        test = globals()[testname](case)

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
