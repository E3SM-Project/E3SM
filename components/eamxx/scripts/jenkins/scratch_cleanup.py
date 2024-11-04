#!/usr/bin/env python3

"""
Clean up old files in the scratch area for mappy.
"""

from pathlib import Path
import re, time, shutil, sys, argparse

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""\n{0} [-c HOURS]
OR
{0} --help

\033[1mEXAMPLES:\033[0m
    \033[1;32m# Purge files older than 20 hours \033[0m
    > {0} -c 20
""".format(Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-c", "--cutoff", type=int, default=30, help="The cutoff age for purging in hours")

    parser.add_argument("-d", "--dry-run", action="store_true", help="Do a dry run, don't actually remove files")

    args = parser.parse_args(args[1:])

    return args

###############################################################################
def scratch_cleanup(cutoff, dry_run):
###############################################################################
    scratch = Path('/ascldap/users/e3sm-jenkins/acme/scratch')

    timestamp_re = re.compile(r'.*(20[0-9]{6}_[0-9]{6}).*')
    timestamps = set()
    for item in scratch.iterdir():
        basename = item.name
        re_match = timestamp_re.match(basename)
        if re_match:
            timestamps.add(re_match.groups()[0])

    tformat = "%Y%m%d_%H%M%S"
    curr_time = time.time()

    for timestamp in timestamps:
        timestamp_time = time.mktime(time.strptime(timestamp, tformat))
        age_in_hours = (curr_time - timestamp_time) / 3600
        if age_in_hours > cutoff:
            print(f"Timestamp {timestamp} is {age_in_hours} hours old and corresponding files will be removed")
            files_to_remove = scratch.glob(f"*{timestamp}*")
            for file_to_remove in files_to_remove:
                print(f"  Removing {file_to_remove}")
                if not dry_run:
                    if file_to_remove.is_dir():
                        shutil.rmtree(file_to_remove)
                    else:
                        file_to_remove.unlink()

    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = scratch_cleanup(**vars(parse_command_line(sys.argv, description)))

    sys.exit(0 if success else 1)

###############################################################################

if __name__ == "__main__":
    _main_func(__doc__)
