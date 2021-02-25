"""
Implementation of the CIME REP test

This test verifies that two identical runs give bit-for-bit results
"""

from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo

class REP(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds = False,
                                       run_two_suffix = 'rep2')

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        pass
