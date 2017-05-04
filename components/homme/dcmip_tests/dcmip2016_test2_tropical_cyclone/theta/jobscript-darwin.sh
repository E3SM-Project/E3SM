#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 2 on a mac running Darwin
#
# usage: ./jobscript-...

# launch the simulation
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                         # launch simulation
