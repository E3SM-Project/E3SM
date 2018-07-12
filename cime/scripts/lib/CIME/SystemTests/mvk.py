"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

import subprocess
import importlib
import shutil
import stat
import json

from CIME.XML.standard_module_setup import *
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.hist_utils import _get_all_hist_files, BLESS_LOG_NAME
from CIME.utils import append_testlog, get_current_commit, get_timestamp, get_model
from CIME.test_status import *


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


    def _generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        with self._test_status:
            # generate baseline

            # BEGIN: modified CIME.hist_utils.generate_baseline
            rundir = self._case.get_value("RUNDIR")
            basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"),
                                       self._case.get_value("BASEGEN_CASE"))
            testcase = self._case.get_value("CASE")

            if not os.path.isdir(basegen_dir):
                os.makedirs(basegen_dir)

            if os.path.isdir(os.path.join(basegen_dir, testcase)):
                expect(False, " Cowardly refusing to overwrite existing baseline directory")

            comments = "Generating baselines into '{}'\n".format(basegen_dir)
            num_gen = 0

            model = 'cam'
            comments += "  generating for model '{}'\n".format(model)
            hists = _get_all_hist_files(testcase, model, rundir)
            logger.debug("mvk_hist_files: {}".format(hists))

            num_gen += len(hists)
            for hist in hists:
                basename = hist[hist.rfind(model):]
                baseline = os.path.join(basegen_dir, basename)
                if os.path.exists(baseline):
                    os.remove(baseline)

                shutil.copy(hist, baseline)
                comments += "    generating baseline '{}' from file {}\n".format(baseline, hist)

            newestcpllogfile = self._case.get_latest_cpl_log(coupler_log_path=self._case.get_value("LOGDIR"))
            if newestcpllogfile is None:
                logger.warning("No cpl.log file found in log directory {}".format(self._case.get_value("LOGDIR")))
            else:
                shutil.copyfile(newestcpllogfile,
                                os.path.join(basegen_dir, "cpl.log.gz"))

            expect(num_gen > 0, "Could not generate any hist files for case '{}', something is seriously wrong".format(
                os.path.join(rundir, testcase)))
            # make sure permissions are open in baseline directory
            for root, _, files in os.walk(basegen_dir):
                for name in files:
                    try:
                        os.chmod(os.path.join(root, name),
                                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
                    except OSError:
                        # We tried. Not worth hard failure here.
                        pass

            if get_model() == "e3sm":
                bless_log = os.path.join(basegen_dir, BLESS_LOG_NAME)
                with open(bless_log, "a") as fd:
                    fd.write("sha:{} date:{}\n".format(get_current_commit(repo=self._case.get_value("CIMEROOT")),
                                                       get_timestamp(timestamp_format="%Y-%m-%d_%H:%M:%S")))
            # END: modified CIME.hist_utils.generate_baseline

            append_testlog(comments)
            status = TEST_PASS_STATUS
            baseline_name = self._case.get_value("BASEGEN_CASE")
            self._test_status.set_status("{}".format(GENERATE_PHASE), status, comments=os.path.dirname(baseline_name))
            basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), self._case.get_value("BASEGEN_CASE"))
            # copy latest cpl log to baseline
            # drop the date so that the name is generic
            newestcpllogfiles = self._get_latest_cpl_logs()
            for cpllog in newestcpllogfiles:
                m = re.search(r"/(cpl.*.log).*.gz", cpllog)
                if m is not None:
                    baselog = os.path.join(basegen_dir, m.group(1)) + ".gz"
                    shutil.copyfile(cpllog,
                                    os.path.join(basegen_dir, baselog))


    def _compare_baseline(self):
        with self._test_status:
            resubmit = int(self._case.get_value("RESUBMIT"))
            logger.info('MVK:CMPR: Resubmits: {}'.format(resubmit))
            if resubmit > 0:
                # This is here because the comparison is run for each submission
                # and we only want to compare once the whole run is finished. We
                # need to return a pass here to continue the submission process.
                self._test_status.set_status(CIME.test_status.BASELINE_PHASE, CIME.test_status.TEST_PASS_STATUS)
                return

            # logger.info('MVK:CMPR: Compare basline? {}'.format(self._case.get_value("COMPARE_BASELINE")))
            #              MVK:CMPR: Compare basline? True
            # logger.info('MVK:CMPR: Baseline root: {}'.format(self._case.get_value("BASELINE_ROOT")))
            #              MVK:CMPR: Baseline root: /lustre/atlas/proj-shared/cli106/kennedy/baselines
            # logger.info('MVK:CMPR: Baseline comparison name: {}'.format(self._case.get_value("BASECMP_CASE")))
            #              MVK:CMPR: Baseline comparison name: jhkennedy/cime/mvk/MVK_PL.ne4_oQU240.FC5AV1C-04P2.titan_pgi

            run_dir = self._case.get_value("RUNDIR")
            case_name = self._case.get_value("CASE")
            basecmp_case = self._case.get_value("BASECMP_CASE")
            base_dir = os.path.join(self._case.get_value("BASELINE_ROOT"), basecmp_case)
            base_name = os.path.basename(basecmp_case)

            eve_config = {
                "{}".format(case_name.split('.')[-1]): {
                    "module": os.path.join(eve_lib_dir, "extensions/ks.py"),
                    "case1": case_name,
                    "dir1": run_dir,
                    "case2": base_name,
                    "dir2": base_dir,
                    "ninst": ninst,
                    "critical": 13
                }
            }

            with open(os.path.join(run_dir, '.'.join([case_name,'json'])), 'w') as config_file:
                json.dump(eve_config, config_file, indent=4)

        self._test_status.set_status(CIME.test_status.BASELINE_PHASE, CIME.test_status.TEST_PASS_STATUS)


