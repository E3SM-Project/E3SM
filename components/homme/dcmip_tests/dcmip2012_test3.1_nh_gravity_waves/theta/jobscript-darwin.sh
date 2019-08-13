#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on a mac running Darwin
#

EXEC=../../../test_execs/theta-nlev20/theta-nlev20
openmpiexec -n 6 $EXEC < ./namelist-h-lowres.nl                         # launch simulation
