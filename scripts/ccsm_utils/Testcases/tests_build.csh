#!/bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv     || exit -2

./$CASE.build

