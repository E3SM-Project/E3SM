"""
Perturbation functions for RCS system test.
"""

import os
import shutil


def duplicate_yaml_file(yaml_file, num_copies):
    """Duplicate a YAML file into multiple copies with four-digit suffixes."""

    if not os.path.isfile(yaml_file):
        raise FileNotFoundError(f"The file {yaml_file} does not exist.")

    for i in range(1, num_copies + 1):
        new_file = f"{yaml_file}_{i:04d}"
        shutil.copyfile(yaml_file, new_file)


def update_yaml_file(yaml_file, seed, pert_out):
    """Update YAML input and output files with perturbation details."""

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
                #    with perturbation_random_seed: {seed}
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
        found_prefix = False

        new_lines = []
        # Process each line
        for line in lines:
            if line.strip().startswith("filename_prefix:"):
                # replace ".scream" with ".scream_{seed:04d}"
                new_lines.append(line.replace(
                    ".scream", f".scream_{seed:04d}"))
                found_prefix = True
            else:
                new_lines.append(line)

        if not found_prefix:
            raise ValueError(f"Couldn't find 'filename_prefix' in {yaml_file}")

        # Write the new lines back to file
        with open(yaml_file, "w", encoding="utf-8") as file:
            file.writelines(new_lines)
