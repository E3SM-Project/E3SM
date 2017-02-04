"""
Implementation of the CIME data assimilation test: Tests simple DA script
which does not change the CAM input. Compares answers to non-DA run.

"""

import os.path
import logging
import glob

import CIME.XML.standard_module_setup as sms
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.utils import expect
from CIME.case_setup import case_setup

###############################################################################
class DAE(SystemTestsCompareTwo):
###############################################################################
    """
    Implementation of the CIME data assimilation test: Tests simple DA script
    which does not change the CAM input. Compares answers to non-DA run.
    Refers to a faux data assimilation script in the
    cime/utils/data_assimilation directory
    """

    ###########################################################################
    def __init__(self, case):
    ###########################################################################
        SystemTestsCompareTwo.__init__(self, case,
                                       separate_builds=False,
                                       run_two_suffix='da',
                                       run_one_description='no data assimilation',
                                       run_two_description='data assimilation')

    ###########################################################################
    def _case_one_setup(self):
    ###########################################################################
        case_setup(self._case, test_mode=True, reset=True)

    ###########################################################################
    def _case_two_setup(self):
    ###########################################################################
        # Set up data assimilation in config_tests.xml once that's ready
        # We need to find the utils/data_assimilation directory
        # LIB_DIR should be our parent dir
        da_dir = os.path.join(os.path.dirname(sms.LIB_DIR), "data_assimilation")
        expect(os.path.isdir(da_dir), "ERROR: da_dir, '%s', does not exist"%da_dir)
        da_file = os.path.join(da_dir, "da_no_data_mod.sh")
        expect(os.path.isfile(da_file), "ERROR: da_file, '%s', does not exist"%da_file)

        # Set up two data assimilation cycles each half of the full run
        self._case.set_value("DATA_ASSIMILATION", "TRUE")
        self._case.set_value("DATA_ASSIMILATION_SCRIPT", da_file)
        self._case.set_value("DATA_ASSIMILATION_CYCLES", 2)
        stopn = self._case.get_value("STOP_N")
        expect((stopn % 2) == 0, "ERROR: DAE test requires that STOP_N be even")
        stopn = stopn / 2
        self._case.set_value("STOP_N", stopn)

        self._case.flush()

    ###########################################################################
    def run_phase(self):
    ###########################################################################
        SystemTestsCompareTwo.run_phase(self)
        # Do some checks on the data assimilation 'output' from case2
        self._activate_case2()
        case_root = self._get_caseroot2()
        da_files = glob.glob(os.path.join(case_root, 'da.log.*'))
        if da_files is None:
            logger = logging.getLogger(__name__)
            logger.warning("No DA files in %s", os.path.join(case_root, 'da.log.*'))

        da_cycles = self._case.get_value("DATA_ASSIMILATION_CYCLES")
        expect((da_files is not None) and (len(da_files) == da_cycles),
               "ERROR: There were %d DA cycles in run but %d DA files were found"%
               (da_cycles, len(da_files) if da_files is not None else 0))
        da_files.sort()
        cycle_num = 0
        for fname in da_files:
            found_caseroot = False
            found_cycle = False
            with open(fname) as dfile:
                for line in dfile:
                    expect(line[0:5] != 'ERROR', "ERROR, error line found in %s"%fname)
                    if line[0:8] == 'caseroot':
                        found_caseroot = True
                    elif line[0:5] == 'cycle':
                        found_cycle = True
                        expect(int(line[7:]) == cycle_num,
                               "ERROR: Wrong cycle (%d) found in %s (expected %d)"%
                               (int(line[7:]), fname, cycle_num))
                    else:
                        expect(False, "ERROR: Unrecognized line ('%s') found in %s"%(line, fname))

                # End of for loop
                expect(found_caseroot, "ERROR: No caseroot found in %s"%fname)
                expect(found_cycle, "ERROR: No cycle found in %s"%fname)
            # End of with
            cycle_num = cycle_num + 1
