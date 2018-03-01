#!/usr/bin/env python
"""
Reads the restart pointer named in the namelist file
provided through the -f flag.  If the file named by
config_restart_timestamp_name does not exists, exits 0
(meaning the run is not done).  If the file does exist,
checks that the current date does not exceed the end
date specified by the -e flag with format
YYYY-MM-DD_hh:mm:ss. (Note: only years between 0001
and 9999 are supported.)  If this date is exceeded, 
exits 1, indicating that the run should stop.  
Otherwise, exits 0, indicating the run has not yet
completed.

Note that this script will correctly detect that
a simulation has completed only if a restart file is
written at the desired end date.  Otherwise, the
script will detect that the simulation should be
finished only after the first restart file is written
out that exceeds the end date.

An example snippet from a bash job script using this
script to perform 12 month-long runs in sequence
might look like the following:

...
runsPerJob=12
endDate=0002-01-01_00:00:00

for iter in (`seq 1 $runsPerJob`)
do
   ./check_progress.py -f namelist.ocean -e $endDate
    if ($? != 0) then
        echo 'run finished'
        exit 1
    endif

    ./run.py

    if ($? != 0) then
        echo 'run failed'
        exit 1
    endif

    ./setup_restart.py -f namelist.ocean -s "'file'"
done

In this example, the ./run.py script could perform a number
of steps following the model run such as modifying forcing,
updating the model state (e.g. for "offline" coupling),
or performing sanity checks.

A number of jobs could be chained with this script. If one of
these jobs exits with exit status 1, this could be made to
cancel dependent jobs (depending on the capabilities of the 
job scheduler).  Even if not, the next job would immediately
exit once it runs, thus using negligible computing time.
"""

import os
import argparse
from datetime import datetime
import sys

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--fileName", dest="fileName", help="A namelist file in which to find the name of the restart pointer", metavar="FILE", required=True)
parser.add_argument("-e", "--endDate", dest="endDate", help="End date of the run", metavar="DATE", required=True)

args = parser.parse_args()

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

inFile = open(restartPointer, 'r')
for line in inFile:
    line = line.strip(" \t\n")
    restartDate = datetime.strptime(line, '%Y-%m-%d_%H:%M:%S')
    break
inFile.close()
endDate = datetime.strptime(args.endDate, '%Y-%m-%d_%H:%M:%S')
if restartDate >= endDate:
    print 'Run has completed.'
    sys.exit(1)
sys.exit(0)
