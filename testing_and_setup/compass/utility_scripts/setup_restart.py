#!/usr/bin/env python
"""
Modifies a namelist file provided with the -f flag.  If the
file named by config_restart_timestamp_name exists, set
config_do_restart = .true. and config_start_time to the
value supplied through the -s flag (default is 'file')
to enable a restart run.

This script is intended to be called at the end of a
job script that will be run multiple times, with
the first beginning from an initial-condition file
and subsequent runs continuing from a restart file.
"""

import os
import argparse
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--fileName", dest="fileName", help="A namelist file to be changed to restart mode", metavar="FILE", required=True)
parser.add_argument("-s", "--startTime", dest="startTime", help="The new value to assign to config_start_time (default is 'file')", metavar="STARTTIME")

args = parser.parse_args()


if args.startTime is None:
    args.startTime = "'file'"

lines = []
restartPointer = 'Restart_timestamp'
inFile = open(args.fileName, 'r')
for line in inFile:
    if 'config_restart_timestamp_name' in line:
        restartPointer = line.split("=")[-1]
        restartPointer = restartPointer.lstrip(" \t'")
        restartPointer = restartPointer.rstrip(" \t\n'")
    lines.append(line)
inFile.close()

if not os.path.exists(restartPointer):
    # nothing to do
    sys.exit(0)

outFile = open(args.fileName, 'w')
for line in lines:
    if 'config_do_restart' in line:
        line = "    config_do_restart = .true.\n"
    if 'config_start_time' in line:
        line = "    config_start_time = %s\n"%args.startTime
    outFile.write(line)
outFile.close()

sys.exit(0)
