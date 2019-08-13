"""
Implementation of the CIME data assimilation test:
Compares standard run with run broken into two data assimilation cycles.
Runs a simple DA script on each cycle which performs checks but does not
change any model state (restart files). Compares answers of two runs.

"""

import os.path
import logging
import glob
import gzip

import CIME.XML.standard_module_setup as sms
from CIME.SystemTests.system_tests_compare_two import SystemTestsCompareTwo
from CIME.utils import expect

###############################################################################
class DAE(SystemTestsCompareTwo):
###############################################################################
    """
    Implementation of the CIME data assimilation test:
    Compares standard run with a run broken into two data assimilation cycles.
    Runs a simple DA script on each cycle which performs checks but does not
    change any model state (restart files). Compares answers of two runs.
    Refers to a faux data assimilation script in the
    cime/scripts/data_assimilation directory
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
        # Even though there may be test mods turning on data assimilation,
        #   case1 is the control so turn it off
        self._case.set_value("DATA_ASSIMILATION_SCRIPT", "")
        self._case.set_value("DATA_ASSIMILATION_CYCLES", 1)

    ###########################################################################
    def _case_two_setup(self):
    ###########################################################################
        # Allow testmods to set an assimilation script
        if len(self._case.get_value("DATA_ASSIMILATION_SCRIPT")) == 0:
            # We need to find the scripts/data_assimilation directory
            # LIB_DIR should be our parent dir
            da_dir = os.path.join(os.path.dirname(sms.LIB_DIR), "data_assimilation")
            expect(os.path.isdir(da_dir), "ERROR: da_dir, '{}', does not exist".format(da_dir))
            da_file = os.path.join(da_dir, "da_no_data_mod.sh")
            expect(os.path.isfile(da_file), "ERROR: da_file, '{}', does not exist".format(da_file))
            # Set up two data assimilation cycles each half of the full run
            self._case.set_value("DATA_ASSIMILATION_SCRIPT", da_file)

        # We need at least 2 DA cycles
        da_cycles = self._case.get_value("DATA_ASSIMILATION_CYCLES")
        if da_cycles < 2:
            da_cycles = 2
            self._case.set_value("DATA_ASSIMILATION_CYCLES", da_cycles)
        stopn = self._case.get_value("STOP_N")
        expect((stopn % da_cycles) == 0, "ERROR: DAE test with {0} cycles requires that STOP_N be divisible by {0}".format(da_cycles))
        stopn = int(stopn / da_cycles)
        self._case.set_value("STOP_N", stopn)

        self._case.flush()

    ###########################################################################
    def run_phase(self): # pylint: disable=arguments-differ
    ###########################################################################
        # Clean up any da.log files in case this is a re-run.
        self._activate_case2()
        case_root = self._get_caseroot2()
        rundir2 = self._case.get_value("RUNDIR")
        da_files = glob.glob(os.path.join(rundir2, 'da.log.*'))
        for file_ in da_files:
            os.remove(file_)
        # End for

        # CONTINUE_RUN ends up TRUE, set it back in case this is a re-run.
        with self._case:
            self._case.set_value("CONTINUE_RUN", False)
            # Turn off post DA in case this is a re-run
            for comp in self._case.get_values("COMP_CLASSES"):
                if comp == "ESP":
                    continue
                else:
                    self._case.set_value("DATA_ASSIMILATION_{}".format(comp), False)

        # Start normal run here
        self._activate_case1()
        SystemTestsCompareTwo.run_phase(self)

        # Do some checks on the data assimilation 'output' from case2
        self._activate_case2()
        da_files = glob.glob(os.path.join(rundir2, 'da.log.*'))
        if da_files is None:
            logger = logging.getLogger(__name__)
            path = os.path.join(case_root, 'da.log.*')
            logger.warning("No DA files in {}".format(path))

        da_cycles = self._case.get_value("DATA_ASSIMILATION_CYCLES")
        expect((da_files is not None) and (len(da_files) == da_cycles),
               "ERROR: There were {:d} DA cycles in run but {:d} DA files were found".format(da_cycles, len(da_files) if da_files is not None else 0))
        da_files.sort()
        cycle_num = 0
        compset = self._case.get_value("COMPSET")
        # Special case for DWAV so we can make sure other variables are set
        is_dwav = '_DWAV' in compset
        for fname in da_files:
            found_caseroot = False
            found_cycle = False
            found_signal = 0
            found_init = 0
            if is_dwav:
                expected_init = self._case.get_value("NINST_WAV")
            else:
                # Expect a signal from every instance of every DA component
                expected_init = 0
                for comp in self._case.get_values("COMP_CLASSES"):
                    if comp == "ESP":
                        continue
                    elif self._case.get_value("DATA_ASSIMILATION_{}".format(comp)):
                        expected_init = expected_init + self._case.get_value("NINST_{}".format(comp))

            # Adjust expected initial run and post-DA numbers
            if cycle_num == 0:
                expected_signal = 0
            else:
                expected_signal = expected_init
                expected_init = 0

            with gzip.open(fname, "r") as dfile:
                for bline in dfile:
                    line = bline.decode("utf-8")
                    expect(not 'ERROR' in line, "ERROR, error line {} found in {}".format(line, fname))
                    if 'caseroot' in line[0:8]:
                        found_caseroot = True
                    elif 'cycle' in line[0:5]:
                        found_cycle = True
                        expect(int(line[7:]) == cycle_num,
                               "ERROR: Wrong cycle ({:d}) found in {} (expected {:d})".format(int(line[7:]), fname, cycle_num))
                    elif 'resume signal' in line:
                        found_signal = found_signal + 1
                        expect('Post-DA resume signal found' in line[0:27],
                               "ERROR: bad post-DA message found in {}".format(fname))
                    elif 'Initial run' in line:
                        found_init = found_init + 1
                        expect('Initial run signal found' in line[0:24],
                               "ERROR: bad Initial run message found in {}".format(fname))
                    else:
                        expect(False, "ERROR: Unrecognized line ('{}') found in {}".format(line, fname))

                # End for
                expect(found_caseroot, "ERROR: No caseroot found in {}".format(fname))
                expect(found_cycle, "ERROR: No cycle found in {}".format(fname))
                expect(found_signal == expected_signal,
                       "ERROR: Expected {} post-DA resume signal message(s), {} found in {}".format(expected_signal, found_signal, fname))
                expect(found_init == expected_init,
                       "ERROR: Expected {} Initial run message(s), {} found in {}".format(expected_init, found_init, fname))
            # End with
            cycle_num = cycle_num + 1
        # End for
