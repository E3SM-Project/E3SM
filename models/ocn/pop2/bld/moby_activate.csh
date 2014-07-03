#!/bin/csh -f
set verbose

source ./Tools/ccsm_getenv || exit -2

rm -f LockedFiles/env_build*  >& /dev/null

./xmlchange -file env_build.xml -id OCN_TRACER_MODULES     -val "iage moby"  || exit -900
./xmlchange -file env_build.xml -id OCN_TRACER_MODULES_OPT -val darwin       || exit -910
./xmlchange -file env_build.xml -id OCN_SUBMODEL           -val moby         || exit -920
./xmlchange -file env_build.xml -id BUILD_COMPLETE -val FALSE                || exit -930
./xmlchange -file env_build.xml -id BUILD_STATUS -val 0                      || exit -940

echo "CCSM $CASE mobysetup FINISHED SUCCESSFULLY"       

exit 0
