"""
Run a testcase.
"""

from CIME.utils import expect
from CIME.case import Case

from CIME.SystemTests.eri import ERI
from CIME.SystemTests.err import ERR
from CIME.SystemTests.erp import ERP
from CIME.SystemTests.ers import ERS
from CIME.SystemTests.ert import ERT
from CIME.SystemTests.lii import LII
from CIME.SystemTests.nck import NCK
from CIME.SystemTests.pea import PEA
from CIME.SystemTests.pem import PEM
from CIME.SystemTests.pet import PET
from CIME.SystemTests.pfs import PFS
from CIME.SystemTests.sms import SMS
from CIME.SystemTests.seq import SEQ
from CIME.SystemTests.ssp import SSP

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
        test._runstatus = "FAIL"
        raise

    return success
