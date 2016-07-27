"""
Implementation of the CIME REP test

This test verifies that two identical runs give bit-for-bit results
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

class REP(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case)

        self.set_test_suffix('rep2')
        self.set_description_first_phase('Run 1')
        self.set_description_second_phase('Run 2')

    def _setup_first_phase(self):
        pass

    def _setup_second_phase(self):
        pass


