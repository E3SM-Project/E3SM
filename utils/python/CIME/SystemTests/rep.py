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

    def _common_setup(self):
        # TODO(wjs, 2016-07-27) Many of these settings are made for most tests -
        # especially setting HIST_OPTION and HIST_N. Should these be moved to
        # some common place to avoid duplication by all tests (and in case some
        # other "standard" settings are added in the future)?
        self._case.set_value("CONTINUE_RUN",False)
        self._case.set_value("REST_OPTION","none")
        self._case.set_value("HIST_OPTION","$STOP_OPTION")
        self._case.set_value("HIST_N","$STOP_N")

    def _setup_first_phase(self):
        pass

    def _setup_second_phase(self):
        pass


