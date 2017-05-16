#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-1 on a mac running Darwin
#

# launch the simulation
EXEC=../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp           # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                          # launch simulation
