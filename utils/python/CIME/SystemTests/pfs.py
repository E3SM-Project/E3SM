"""
CIME performance test  This class inherits from SystemTestsCommon
"""
from CIME.XML.standard_module_setup import *
from system_tests_common import SystemTestsCommon

class PFS(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the PFS system test
        """
        SystemTestsCommon.__init__(self, case)
