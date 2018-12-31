#!/usr/bin/env python

import unittest
import os
from CIME import test_status
from CIME.tests.custom_assertions_test_status import CustomAssertionsTestStatus

class TestTestStatus(CustomAssertionsTestStatus):

    _TESTNAME = 'fake_test'

    def setUp(self):
        self._ts = test_status.TestStatus(test_dir=os.path.join('nonexistent', 'path'),
                                          test_name=self._TESTNAME,
                                          no_io=True)
        self._set_core_phases_to_pass()

    def _set_core_phases_to_pass(self):
        """Set all core phases of self._ts to pass status"""
        with self._ts:
            for phase in test_status.CORE_PHASES:
                self._ts.set_status(phase, test_status.TEST_PASS_STATUS)

    # ------------------------------------------------------------------------
    # Tests of TestStatus.phase_statuses_dump
    # ------------------------------------------------------------------------

    def test_psdump_corePhasesPass(self):
        output = self._ts.phase_statuses_dump()
        self.assert_core_phases(output, self._TESTNAME, fails=[])


if __name__ == '__main__':
    unittest.main()
