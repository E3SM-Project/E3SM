"""
TODO: Add description here.

This class inherits from SystemTestsCommon.
"""

import os
import shutil

from CIME.SystemTests.system_tests_common import SystemTestsCommon
from CIME.case.case_setup import case_setup


class MVKxx(SystemTestsCommon):
    """Multivariate K-S test using multi-instance capability"""

    def setup_phase(
        self, clean=False, test_mode=False, reset=False, keep=False, disable_git=False
    ):
        """setup phase implementation"""

        self.setup_indv(
            clean=clean,
            test_mode=test_mode,
            reset=reset,
            keep=keep,
            disable_git=disable_git,
        )

        self._case.flush()

        case_setup(self._case, test_mode=False, reset=True)
        run_dir = self._case.get_value("RUNDIR")

        n_inst = int(self._case.get_value("NINST_ATM"))
        if n_inst > 1:
            duplicate_yaml_files(run_dir + "/data/scream_input.yaml", n_inst)
            duplicate_yaml_files(
                run_dir + "/data/monthly_average.yaml", n_inst)
            # Let's update the perturbation seed in the YAML files
            for i in range(1, n_inst + 1):
                yaml_file = f"{run_dir}/data/scream_input.yaml_{i:04d}"
                out_file = f"{run_dir}/data/monthly_average.yaml_{i:04d}"
                if not os.path.isfile(yaml_file):
                    raise FileNotFoundError(
                        f"File {yaml_file} does not exist.")
                if not os.path.isfile(out_file):
                    raise FileNotFoundError(f"File {out_file} does not exist.")
                update_yaml_perturbation_seed(yaml_file, i, "pert")
                update_yaml_perturbation_seed(out_file, i, "out")
        else:
            msg = (
                f"NINST_ATM = {n_inst}. This test requires NINST_ATM > 1. "
                "Consider setting NINST_ATM > 1 in your env_run.xml "
                "or use _C# specifier in test name for a multi-driver "
                "multi-instance setup (producing # pelayout copies), "
                "or _N# for a single-driver multi-instance setup "
                "(dividing specified pelayout among # instances)."
            )
            raise ValueError(msg)


# TEST UTILITY FUNCTIONS BELOW
# - duplicate_yaml_files is used in setup
# - update_yaml_perturbation_seed is used in setup


def duplicate_yaml_files(yaml_file, num_copies):
    """Duplicate a YAML file into multiple copies with four-digit suffixes."""

    if not os.path.isfile(yaml_file):
        raise FileNotFoundError(f"The file {yaml_file} does not exist.")

    for i in range(1, num_copies + 1):
        new_file = f"{yaml_file}_{i:04d}"
        shutil.copyfile(yaml_file, new_file)


def update_yaml_perturbation_seed(yaml_file, seed, pert_out):
    """Update the perturbation seed."""

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
                        "monthly_average.yaml", f"monthly_average.yaml_{seed:04d}"
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
