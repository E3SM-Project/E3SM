#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 3 on a mac running Darwin
#
# usage: ./jobscript-...

# launch the simulation
EXEC=../../../test_execs/preqx-nlev40-interp/preqx-nlev40-interp        # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
