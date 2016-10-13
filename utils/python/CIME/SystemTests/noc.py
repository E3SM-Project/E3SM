"""
CIME NOC test  This class inherits from SystemTestsCompareTwo
The NOC test is a multi-instance validation for single instance ocean (default length)
  do an initial run test with NINST 2 (other than ocn), with mod to instance 1 (suffix: inst1_base, inst2_mod)
  do an initial run test with NINST 2 (other than ocn), with mod to instance 2 (suffix: inst1_base, inst2_mod)
  compare inst1_base with inst2_base
  compare inst1_mod  with inst2_mod
"""
from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

logger = logging.getLogger(__name__)

class NOC(SystemTestsCompareTwo):

    def __init__(self, case):
        """
        initialize an object interface to the NOC system test
        """
        #expectedrunvars = ["CONTINUE_RUN", "REST_OPTION", "HIST_OPTION", "HIST_N"]
        SystemTestsCompareTwo.__init__(self, case)
