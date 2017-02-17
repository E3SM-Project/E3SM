#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-0 on a mac running Darwin
#
# usage: ./jobscript-...

#EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
