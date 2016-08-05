"""
Implementation of the CIME REP test

This test verifies that two identical runs give bit-for-bit results
"""

from CIME.SystemTests.system_tests_compare_two_clone import SystemTestsCompareTwoClone

class REP(SystemTestsCompareTwo):

    def __init__(self, case):
        SystemTestsCompareTwoClone.__init__(self, case,
                                            run_two_suffix = 'rep2')

    def _common_setup(self):
        # TODO(wjs, 2016-07-27) Many of these settings are made for most tests -
        # especially setting HIST_OPTION and HIST_N. Should these be moved to
        # some common place to avoid duplication by all tests (and in case some
        # other "standard" settings are added in the future)?
        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")

    def _case_one_setup(self):
        pass

    def _case_two_setup(self):
        pass


