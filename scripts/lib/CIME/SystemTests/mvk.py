"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

import os
import json
import logging

from distutils import dir_util

import CIME.test_status
import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup
from CIME.XML.machines import Machines


import evv4esm  # pylint: disable=import-error
from evv4esm.__main__ import main as evv  # pylint: disable=import-error

evv_lib_dir = os.path.abspath(os.path.dirname(evv4esm.__file__))
logger = logging.getLogger(__name__)

NINST = 20


class MVK(SystemTestsCommon):

    def __init__(self, case):
        """
        initialize an object interface to the MVK test
        """
        SystemTestsCommon.__init__(self, case)

        if self._case.get_value("RESUBMIT") == 0 \
                and self._case.get_value("GENERATE_BASELINE") is False:
            self._case.set_value("COMPARE_BASELINE", True)
        else:
            self._case.set_value("COMPARE_BASELINE", False)

    def build_phase(self, sharedlib_only=False, model_only=False):
        # Only want this to happen once. It will impact the sharedlib build
        # so it has to happen there.
        if not model_only:
            logging.warning('Starting to build multi-instance exe')
            for comp in self._case.get_values("COMP_CLASSES"):
                self._case.set_value('NTHRDS_{}'.format(comp), 1)

                ntasks = self._case.get_value("NTASKS_{}".format(comp))

                self._case.set_value('NTASKS_{}'.format(comp), ntasks * NINST)
                if comp != 'CPL':
                    self._case.set_value('NINST_{}'.format(comp), NINST)

            self._case.set_value('ATM_NCPL', 18)

            self._case.flush()

            case_setup(self._case, test_mode=False, reset=True)

        self.build_indv(sharedlib_only=sharedlib_only, model_only=model_only)

        for iinst in range(1, NINST + 1):
            with open('user_nl_cam_{:04d}'.format(iinst), 'w') as nl_atm_file:
                nl_atm_file.write('new_random = .true.\n')
                nl_atm_file.write('pertlim = 1.0e-10\n')
                nl_atm_file.write('seed_custom = {}\n'.format(iinst))

    def _generate_baseline(self):
        """
        generate a new baseline case based on the current test
        """
        super(MVK, self)._generate_baseline()

        with CIME.utils.SharedArea():
            basegen_dir = os.path.join(self._case.get_value("BASELINE_ROOT"),
                                       self._case.get_value("BASEGEN_CASE"))

            rundir = self._case.get_value("RUNDIR")
            ref_case = self._case.get_value("RUN_REFCASE")

            model = 'cam'
            env_archive = self._case.get_env("archive")
            hists = env_archive.get_all_hist_files(self._case.get_value("CASE"), model, rundir, [r'h\d*.*\.nc'], ref_case=ref_case)
            logger.debug("MVK additional baseline files: {}".format(hists))
            hists = [os.path.join(rundir,hist) for hist in hists]
            for hist in hists:
                basename = hist[hist.rfind(model):]
                baseline = os.path.join(basegen_dir, basename)
                if os.path.exists(baseline):
                    os.remove(baseline)

                CIME.utils.safe_copy(hist, baseline, preserve_meta=False)

    def _compare_baseline(self):
        with self._test_status:
            if int(self._case.get_value("RESUBMIT")) > 0:
                # This is here because the comparison is run for each submission
                # and we only want to compare once the whole run is finished. We
                # need to return a pass here to continue the submission process.
                self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                             CIME.test_status.TEST_PASS_STATUS)
                return

            self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                         CIME.test_status.TEST_FAIL_STATUS)

            run_dir = self._case.get_value("RUNDIR")
            case_name = self._case.get_value("CASE")
            base_dir = os.path.join(self._case.get_value("BASELINE_ROOT"),
                                    self._case.get_value("BASECMP_CASE"))

            test_name = "{}".format(case_name.split('.')[-1])
            evv_config = {
                test_name: {
                    "module": os.path.join(evv_lib_dir, "extensions", "ks.py"),
                    "test-case": "Test",
                    "test-dir": run_dir,
                    "ref-case": "Baseline",
                    "ref-dir": base_dir,
                    "var-set": "default",
                    "ninst": NINST,
                    "critical": 13
                }
            }

            json_file = os.path.join(run_dir, '.'.join([case_name, 'json']))
            with open(json_file, 'w') as config_file:
                json.dump(evv_config, config_file, indent=4)

            evv_out_dir = os.path.join(run_dir, '.'.join([case_name, 'evv']))
            evv(['-e', json_file, '-o', evv_out_dir])

            with open(os.path.join(evv_out_dir, 'index.json')) as evv_f:
                evv_status = json.load(evv_f)

            comments = ""
            for evv_elem in evv_status['Data']['Elements']:
                if evv_elem['Type'] == 'ValSummary' \
                        and evv_elem['TableTitle'] == 'Kolmogorov-Smirnov test':
                    comments = "; ".join("{}: {}".format(key, val) for key, val
                                         in evv_elem['Data'][test_name][''].items())
                    if evv_elem['Data'][test_name]['']['Test status'].lower() == 'pass':
                        self._test_status.set_status(CIME.test_status.BASELINE_PHASE,
                                                     CIME.test_status.TEST_PASS_STATUS)
                    break

                status = self._test_status.get_status(CIME.test_status.BASELINE_PHASE)
                mach_name = self._case.get_value("MACH")
                mach_obj = Machines(machine=mach_name)
                htmlroot = CIME.utils.get_htmlroot(mach_obj)
                urlroot = CIME.utils.get_urlroot(mach_obj)
                if htmlroot is not None:
                    with CIME.utils.SharedArea():
                        dir_util.copy_tree(evv_out_dir, os.path.join(htmlroot, 'evv', case_name), preserve_mode=False)
                    if urlroot is None:
                        urlroot = "[{}_URL]".format(mach_name.capitalize())
                    viewing = "{}/evv/{}/index.html".format(urlroot, case_name)
                else:
                    viewing = "{}\n" \
                              "    EVV viewing instructions can be found at: " \
                              "        https://github.com/E3SM-Project/E3SM/blob/master/cime/scripts/" \
                              "climate_reproducibility/README.md#test-passfail-and-extended-output" \
                              "".format(evv_out_dir)

                comments = "{} {} for test '{}'.\n" \
                           "    {}\n" \
                           "    EVV results can be viewed at:\n" \
                           "        {}".format(CIME.test_status.BASELINE_PHASE, status, test_name, comments, viewing)

                CIME.utils.append_testlog(comments, self._orig_caseroot)
