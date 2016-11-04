#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on a mac running Darwin
#

EXEC=../../test_execs/preqx-nlev20-interp/preqx-nlev20-interp           # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
