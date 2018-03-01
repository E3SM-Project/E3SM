#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-1 on a mac running Darwin
#
# usage: ./jobscript-...

# launch the simulation
EXEC=../../../test_execs/theta-nlev60/theta-nlev60
openmpiexec -n 6 $EXEC < ./namelist-h-lowres.nl                         # launch simulation
