#!/bin/csh -f
set verbose

rm -f LockedFiles/env_build*  >& /dev/null

./xmlchange -file env_build.xml -id OCN_TRACER_MODULES     -val "iage"  || exit -900
./xmlchange -file env_build.xml -id OCN_TRACER_MODULES_OPT -val ""      || exit -910
./xmlchange -file env_build.xml -id OCN_SUBMODEL           -val ""      || exit -920
./xmlchange -file env_build.xml -id BUILD_COMPLETE         -val FALSE   || exit -930
./xmlchange -file env_build.xml -id BUILD_STATUS           -val 0       || exit -940

set $CASE = `./xmlquery CASE -value`
echo "CCSM $CASE mobysetup FINISHED SUCCESSFULLY"       

exit 0
