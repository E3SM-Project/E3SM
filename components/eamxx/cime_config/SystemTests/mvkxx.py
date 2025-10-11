"""
Multivariate test for climate reproducibility using the Kolmogrov-Smirnov (K-S)
test and based on The CESM/E3SM model's multi-instance capability is used to
conduct an ensemble of simulations starting from different initial conditions.

This class inherits from SystemTestsCommon.
"""

import os
import glob
import shutil
import logging

import CIME.utils
from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup


logger = logging.getLogger(__name__)

num_instances = os.getenv("MVKXX_NINST", "6")
n_inst = int(num_instances)

logger.info("MVKxx using n_inst=%d", n_inst)


def duplicate_yaml_files(yaml_file, num_copies):
    """Duplicate a YAML file into multiple copies with four-digit suffixes."""

    if not os.path.isfile(yaml_file):
        raise FileNotFoundError(f"The file {yaml_file} does not exist.")

    for i in range(1, num_copies + 1):
        new_file = f"{yaml_file}_{i:04d}"
        shutil.copyfile(yaml_file, new_file)


def update_yaml_perturbation_seed(yaml_file, seed, pert_out):
    """Update the perturbation seed in a YAML file using basic text manipulation."""

    # Read the file content
    with open(yaml_file, "r", encoding="utf-8") as file:
        lines = file.readlines()

    if pert_out == "pert":
        found_seed = False
        found_output = False
        new_lines = []

        # Process each line
        for line in lines:
            if line.strip().startswith("perturbation_random_seed:"):
                # replace perturbation_random_seed: 0
                #    with perturbation_random_seed: <seed>
                new_lines.append(
                    line.replace(
                        "perturbation_random_seed: 0",
                        f"perturbation_random_seed: {seed}",
                    )
                )
                found_seed = True
            elif "monthly_average.yaml" in line.strip():
                # replace "monthly_average.yaml"
                #    with "monthly_average.yaml_{seed:04d}"
                new_lines.append(
                    line.replace(
                        "monthly_average.yaml",
                        f"monthly_average.yaml_{seed:04d}"
                    )
                )
                found_output = True
            else:
                new_lines.append(line)

        if not found_seed:
            raise ValueError(f"'perturbation_random_seed' NOT in {yaml_file}")
        if not found_output:
            raise ValueError(f"'monthly_average.yaml' NOT in {yaml_file}")

        # Write back to file
        with open(yaml_file, "w", encoding="utf-8") as file:
            file.writelines(new_lines)

    elif pert_out == "out":
        # Track if we found and updated required fields
        new_lines = []

        # Process each line
        for line in lines:
            if line.strip().startswith("filename_prefix:"):
                # replace ".scream" with ".scream_{seed:04d}"
                new_lines.append(line.replace(
                    ".scream", f".scream_{seed:04d}"))
                found_seed = True
            else:
                new_lines.append(line)

        # Add missing sections if needed
        if not found_seed:
            raise ValueError(f"Couldn't find 'filename_prefix' in {yaml_file}")

        # Write back to file
        with open(yaml_file, "w", encoding="utf-8") as file:
            file.writelines(new_lines)


class MVKxx(SystemTestsCommon):
    """Multivariate K-S test using multi-instance capability"""

    def __init__(self, case, **kwargs):
        """
        initialize an object interface to the MVKxx test
        """
        SystemTestsCommon.__init__(self, case, **kwargs)

        if self._case.get_value("MODEL") == "e3sm":
            self.component = "scream"
        else:
            self.component = "cam"

    def setup_phase(self, clean=False, test_mode=False, reset=False, keep=False, disable_git=False):
        """setup phase implementation"""

        self.setup_indv(
            clean=clean,
            test_mode=test_mode,
            reset=reset,
            keep=keep,
            disable_git=disable_git,
        )

        logging.info("Starting to set up multi-instance exe")
        for comp in self._case.get_values("COMP_CLASSES"):
            n_tasks = self._case.get_value("NTASKS_{}".format(comp))
            self._case.set_value(f"NTASKS_{comp}", n_tasks * n_inst)
            if comp != "CPL":
                self._case.set_value(f"NINST_{comp}", n_inst)

        self._case.flush()

        case_setup(self._case, test_mode=False, reset=True)
        run_dir = self._case.get_value("RUNDIR")
        duplicate_yaml_files(run_dir + "/data/scream_input.yaml", n_inst)
        duplicate_yaml_files(run_dir + "/data/monthly_average.yaml", n_inst)

        # Let's update the perturbation seed in the YAML files
        for i in range(1, n_inst + 1):
            yaml_file = f"{run_dir}/data/scream_input.yaml_{i:04d}"
            out_file = f"{run_dir}/data/monthly_average.yaml_{i:04d}"
            if not os.path.isfile(yaml_file):
                raise FileNotFoundError(f"File {yaml_file} does not exist.")
            if not os.path.isfile(out_file):
                raise FileNotFoundError(f"File {out_file} does not exist.")
            update_yaml_perturbation_seed(yaml_file, i, "pert")
            update_yaml_perturbation_seed(out_file, i, "out")

    def run_phase(self):
        """override run phase ?"""
        self.run_indv()

    def _generate_baseline(self):
        """generate a new baseline case based on the current test"""
        super()._generate_baseline()

        with CIME.utils.SharedArea():
            base_gen_dir = os.path.join(
                self._case.get_value("BASELINE_ROOT"),
                self._case.get_value("BASEGEN_CASE"),
            )

            run_dir = self._case.get_value("RUNDIR")
            ref_case = self._case.get_value("RUN_REFCASE")

            env_archive = self._case.get_env("archive")
            hists = env_archive.get_all_hist_files(
                self._case.get_value("CASE"), self.component, run_dir, ref_case=ref_case
            )
            # for eamxx, we need to get all files that have
            # *scream_????.h.*.nc added to this list
            more_hists = glob.glob(os.path.join(
                run_dir, "*scream_????.h.AVERAGE.*.nc"))
            # before copying, let's also rename some files
            # current pattern is scream_????.h.AVERAGE.nmonths_x1.????-??.nc
            # desired pattern is scream_????.h.????-??.nc
            for hist in more_hists:
                if "scream" in hist:
                    # get rid of the AVERAGE.nmonths_x1.
                    new_hist = hist.replace("AVERAGE.nmonths_x1.", "")
                    os.rename(hist, new_hist)
                    # add it to hists
                    hists.append(new_hist)
            logger.debug("MVKxx additional baseline files: %s", hists)
            hists = [os.path.join(run_dir, hist) for hist in hists]
            for hist in hists:
                base_name = hist[hist.rfind(self.component):]
                baseline = os.path.join(base_gen_dir, base_name)
                if os.path.exists(baseline):
                    os.remove(baseline)

                CIME.utils.safe_copy(hist, baseline, preserve_meta=False)
