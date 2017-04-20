#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 3 on a mac running Darwin
#
# usage: ./jobscript-...

# launch the simulation
EXEC=../../../test_execs/theta-nlev40/theta-nlev40
#openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
openmpiexec -n 1 $EXEC < ./namelist-lowres.nl                           # launch simulation
