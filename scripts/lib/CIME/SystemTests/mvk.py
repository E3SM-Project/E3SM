"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

import subprocess
import importlib

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup

import CIME.test_status


logger = logging.getLogger(__name__)

try:
    import livvkit
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'livvkit'])
    livvkit = importlib.import_module('livvkit')
    logger.info('Installed LIVVkit v{} via pip.'.format(livvkit.__version__))

# FIXME: remove once eve is a package, or inc. eve into CIME/SystemTest/Test_utils
eve_lib_dir = '/ccs/home/kennedy/EVE/eve'
sys.path.append(eve_lib_dir)
import eve


# Build executable with multiple instances
ninst = 20


class MVK(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the MVK test
        """
        SystemTestsCommon.__init__(self, case)

        logger.info('MVK:INIT: {}'.format(self._compare_baseline()))


    def build_phase(self, sharedlib_only=False, model_only=False):
        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warn('Starting to build multi-instance exe')
            for comp in self._case.get_values("COMP_CLASSES"):
                self._case.set_value('NTHRDS_{}'.format(comp), 1)

                ntasks = self._case.get_value("NTASKS_{}".format(comp))
                logger.info('MVK:NTASKS_{}: {}'.format(comp, ntasks))

                self._case.set_value('NTASKS_{}'.format(comp), ntasks*ninst)
                if comp != 'CPL':
                    self._case.set_value('NINST_{}'.format(comp), ninst)

            self._case.set_value('ATM_NCPL', 18)

            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        # =================================================================
        # Run-time settings.
        # Do this already in build_phase so that we can check the xml and
        # namelist files before job starts running.
        # =================================================================

        # namelist specifications for each instance
        for iinst in range(1, ninst+1):
            with open('user_nl_cam_{:04d}'.format(iinst), 'w') as nl_atm_file:
                nl_atm_file.write('new_random = .true.\n')
                nl_atm_file.write('pertlim = 1.0e-10\n')
                nl_atm_file.write('seed_custom = {}\n'.format(iinst))


    def run_phase(self):

        self.run_indv()


    # def _generate_baseline(self):
    #     run_dir = self._case.get_value("RUNDIR")
    #     basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), self._case.get_value("BASEGEN_CASE"))
    #     logger.info('MVK:GENDIR: {}'.format(basegen_dir))
    #
    #     with self._test_status:
    #         self._test_status.set_status(CIME.test_status.GENERATE_PHASE, CIME.test_status.TEST_PASS_STATUS)


    # def _compare_baseline(self):
    #     logger.info('MVK:CMPR: {}'.format(self._case.get_value("COMPARE_BASELINE")))
    #     logger.info('MVK:CMPR: {}'.format(type(self._case.get_value("COMPARE_BASELINE"))))
    #
    #     run_dir = self._case.get_value("RUNDIR")
    #     case_name = self._case.get_value("CASE")
    #     basegen_case = self._case.get_value("BASEGEN_CASE")
    #     basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), self._case.get_value("BASEGEN_CASE"))
    #
    #     eve_config = {
    #         "MVK_{}".format(case_name): {
    #             "module": os.path.join(eve_lib_dir, "extensions/ks.py"),
    #             "case1": case_name,
    #             "dir1": run_dir,
    #             "case2": basegen_case,
    #             "dir2": basegen_dir,
    #             "ninst": ninst,
    #             "critical": 13
    #         }
    #     }
    #     logger.info('MVK:EVECONFIG: {}'.format(eve_config))
    #
    #     with self._test_status:
    #         self._test_status.set_status(CIME.test_status.BASELINE_PHASE, CIME.test_status.TEST_PASS_STATUS)