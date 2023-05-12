#!/bin/bash -e
{{ batchdirectives }}

# template to create a case run shell script. This should only ever be called
# by case.submit when on batch. Use case.submit from the command line to run your case.

# cd to case
caseroot={{ caseroot }}
cd $caseroot

# Set PYTHONPATH so we can make cime calls if needed
LIBDIR={{ cimeroot }}
export PYTHONPATH=$LIBDIR:$PYTHONPATH

# setup environment
source .env_mach_specific.sh

# get new lid
lid=$(python3 -c 'import CIME.utils; print(CIME.utils.new_lid())')
export LID=$lid

# Clean/make timing dirs
RUNDIR=$(./xmlquery RUNDIR --value)
if [ -e $RUNDIR/timing ]; then
    /bin/rm -rf $RUNDIR/timing
fi
mkdir -p $RUNDIR/timing/checkpoints

# minimum namelist action
./preview_namelists --component cpl
#./preview_namelists # uncomment for full namelist generation

# uncomment for lockfile checking
# ./check_lockedfiles

# setup OMP_NUM_THREADS
export OMP_NUM_THREADS=$(./xmlquery THREAD_COUNT --value)

# save prerun provenance?

# MPIRUN!
cd $(./xmlquery RUNDIR --value)
{{ mpirun }}

# get timing
python3 -c '
import os
from CIME.case import Case
from CIME.get_timing import get_timing

with Case("'$caseroot'", read_only=True) as case:
    get_timing(case, "'$lid'")
'

# save logs?

# save postrun provenance?

# resubmit ?
