"""
Reproducible Climate Statistics testing using multi-instance capability.

This test runs multiple EAMxx instances with different perturbation seeds
and uses statistical tests to verify that the climate state is identical
between different runs.

RCS inherits from SystemTestsCommon and only overrides:
- setup_phase (to setup multi-instance with per-instance perturbed seeds)
- _generate_baseline (move needed hist files to baseline directory)
- _compare_baseline methods (to run the statistical tests)

RCS relies on two util files:
- rcs_perts.py: functions to duplicate and modify yaml files
- rcs_stats.py: functions to conduct statistical testing
"""

import os
import glob
import logging
import sys

import CIME.test_status
import CIME.utils
from CIME.status import append_testlog
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup

logger = logging.getLogger(__name__)


# pylint: disable=too-few-public-methods
class RCS(SystemTestsCommon):
    """Reproducible Climate Statistics Test using multi-instance capability"""

    # File pattern for multi-instance ensemble output files
    # Use ???? as placeholder for 4-digit instance number
    ENSEMBLE_FILE_PATTERN = "*.scream_????.h.AVERAGE.*.nc"

    # pylint: disable=too-many-arguments, too-many-positional-arguments
    def setup_phase(
        self,
        clean=False,
        test_mode=False,
        reset=False,
        keep=False,
        disable_git=False,
    ):
        """setup phase implementation"""
        # first call the parent method and flush
        self.setup_indv(
            clean=clean,
            test_mode=test_mode,
            reset=reset,
            keep=keep,
            disable_git=disable_git,
        )
        self._case.flush()
        # and again...?
        case_setup(self._case, test_mode=False, reset=True)

        # get run directory
        run_dir = self._case.get_value("RUNDIR")
        # get n_inst
        n_inst = int(self._case.get_value("NINST_ATM"))
        # return early if n_inst <= 1
        # we really don't want people to run this test with n_inst=1
        if n_inst <= 1:
            msg = (
                f"NINST_ATM = {n_inst}. This test requires NINST_ATM > 1. "
                "Consider setting NINST_ATM > 1 in your env_run.xml "
                "or use _C# specifier in test name for a multi-driver "
                "multi-instance setup (producing # pelayout copies), "
                "or _N# for a single-driver multi-instance setup "
                "(dividing specified pelayout among # instances)."
            )
            raise ValueError(msg)

        # get rcs_perts functions
        # but first add the directory to sys.path if not already there
        rcs_perts_path = os.path.join(
            os.path.dirname(__file__), 'rcs_perts.py'
        )
        if not os.path.exists(rcs_perts_path):
            raise ImportError(
                f"Cannot find rcs_perts.py at {rcs_perts_path}"
            )
        if os.path.dirname(__file__) not in sys.path:
            sys.path.insert(0, os.path.dirname(__file__))
        # pylint: disable=import-outside-toplevel
        from rcs_perts import duplicate_yaml_file, update_yaml_file

        # duplicate the yaml files n_inst times
        duplicate_yaml_file(f"{run_dir}/data/scream_input.yaml", n_inst)
        duplicate_yaml_file(f"{run_dir}/data/monthly_average.yaml", n_inst)
        # Let's update the perturbation properties inside the yaml files
        # this handles unique seeds and unique output files manually
        for i in range(1, n_inst + 1):
            yaml_file = f"{run_dir}/data/scream_input.yaml_{i:04d}"
            out_file = f"{run_dir}/data/monthly_average.yaml_{i:04d}"
            if not os.path.isfile(yaml_file):
                raise FileNotFoundError(
                    f"File {yaml_file} does not exist.")
            if not os.path.isfile(out_file):
                raise FileNotFoundError(f"File {out_file} does not exist.")
            update_yaml_file(yaml_file, i, "pert")
            update_yaml_file(out_file, i, "out")

    def _generate_baseline(self):
        """generate a new baseline case based on the current test"""
        # might as well call the parent method first
        super()._generate_baseline()

        with CIME.utils.SharedArea():
            # get the baseline and run directories
            base_gen_dir = os.path.join(
                self._case.get_value("BASELINE_ROOT"),
                self._case.get_value("BASEGEN_CASE"),
            )
            run_dir = self._case.get_value("RUNDIR")

            # Get all files that match the ensemble pattern
            hists = glob.glob(
                os.path.join(run_dir, self.ENSEMBLE_FILE_PATTERN)
            )
            hist_files = [os.path.basename(h) for h in hists]

            for hist in hist_files:
                src = os.path.join(run_dir, hist)
                tgt = os.path.join(base_gen_dir, hist)
                # remove baselines if they exist
                # this is safe because cime forces users to use -o
                if os.path.exists(tgt):
                    os.remove(tgt)

                # log and copy
                logger.info(
                    "Copying ... \n \t %s \n ... to ... \n \t %s \n\n",
                    src, tgt
                )
                CIME.utils.safe_copy(src, tgt, preserve_meta=False)

    def _compare_baseline(self):
        """compare phase implementation"""
        with self._test_status as ts:
            # if we are resubmitting, then we don't do the comparison
            if int(self._case.get_value("RESUBMIT")) > 0:
                ts.set_status(
                    CIME.test_status.BASELINE_PHASE,
                    CIME.test_status.TEST_PASS_STATUS
                )
                return

            # set to FAIL to start with, will update later
            ts.set_status(
                CIME.test_status.BASELINE_PHASE,
                CIME.test_status.TEST_FAIL_STATUS
            )

            # get the run and baseline directories
            run_dir = self._case.get_value("RUNDIR")
            base_dir = os.path.join(
                self._case.get_value("BASELINE_ROOT"),
                self._case.get_value("BASECMP_CASE"),
            )

            # launch the statistics tests
            # first, import rcs_stats funcs from the other file
            rcs_stats_path = os.path.join(
                os.path.dirname(__file__), 'rcs_stats.py'
            )
            if not os.path.exists(rcs_stats_path):
                raise ImportError(
                    f"Cannot find rcs_stats.py at {rcs_stats_path}"
                )
            # Add the directory to sys.path if not already there
            if os.path.dirname(__file__) not in sys.path:
                sys.path.insert(0, os.path.dirname(__file__))
            # note be extra safe and import whole file
            # because we want to avoid import errors of needed pkgs
            # pylint: disable=import-outside-toplevel
            import rcs_stats as rcss
            # now, launch
            comments, new_ts = rcss.run_stats_comparison(
                run_dir,
                base_dir,
                analysis_type="spatiotemporal",
                test_type="ks",
                alpha=0.01,
                run_file_pattern=self.ENSEMBLE_FILE_PATTERN,
                base_file_pattern=self.ENSEMBLE_FILE_PATTERN,
            )

            if new_ts == "PASS":
                out_ts = CIME.test_status.TEST_PASS_STATUS
            else:
                out_ts = CIME.test_status.TEST_FAIL_STATUS

            # log the results and set the test status
            append_testlog(comments, self._orig_caseroot)
            ts.set_status(
                CIME.test_status.BASELINE_PHASE, out_ts
            )
