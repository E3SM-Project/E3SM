#!usr/bin/env python
from __future__ import division, print_function

import argparse
import os

"""
Usage: metrics_checker.py [options]
Example: python metrics_checker.py -t /lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/previous_output/all_sets_v2_6_1_20220328_d62f554/ -r /lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/all_sets/

Options:
  -t FILE, Path to test e3sm_diags results directory
  -r FILE, Path to reference e3sm_diags results directory

About:
This script is used to compare seasonal mean tables between a rest and reference e3sm_diags run,
and to print out lines of variables being changed in test.
"""


parser = argparse.ArgumentParser()
parser.add_argument(
    "--ref_path",
    "-r",
    dest="ref",
    help="Path to reference e3sm_diags output",
    metavar="FILE",
)
parser.add_argument(
    "--test_path",
    "-t",
    dest="test",
    help="Path to test e3sm_diags output",
    metavar="FILE",
)

args = parser.parse_args()

seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

ref_path = args.ref
test_path = args.test


def compare_metrics(ref_path, test_path, season):
    fref = os.path.join(ref_path, "viewer/table-data", f"{season}_metrics_table.csv")
    ftest = os.path.join(test_path, "viewer/table-data", f"{season}_metrics_table.csv")
    try:
        with open(fref, "r") as ref, open(ftest, "r") as test:
            file_ref = ref.readlines()
            header = file_ref[0]
            file_test = test.readlines()
            # print(file_test)
            num_matching = -1
            num_missing = 0
            num_ref = -1  # header lines are same therefor to -1
            num_addition = len(file_test) - len(file_ref)
            for line in file_ref:
                num_ref = num_ref + 1
                varid = line.split(",")[0]
                if line not in file_test:
                    print(f"Found difference in {season}", line)
                    matching_varid = [s for s in file_test if varid in s]
                    if len(matching_varid):
                        print(header)
                        print("ref :", line)
                        print("test:", matching_varid[0])
                    else:
                        num_missing = num_missing + 1
                        print(f"{varid} is missing in test dataset")
                else:
                    num_matching = num_matching + 1

            for line in file_test:
                varid = line.split(",")[0]
                if line not in file_ref:
                    matching_varid = [s for s in file_ref if varid in s]
                    if not len(matching_varid):
                        print(f"{varid} is added in test dataset")
                else:
                    num_matching = num_matching + 1
            print(
                f"\nSUMMARY for {season}: {num_matching} out of {num_ref} have matching metrics with ref files,\n           {num_missing} variables are missing in test datasets;\n           {num_addition} more variables are present in test data."
            )

    except Exception as e:
        print("Failed to open file:" + str(e))


for season in seasons:
    compare_metrics(ref_path, test_path, season)
