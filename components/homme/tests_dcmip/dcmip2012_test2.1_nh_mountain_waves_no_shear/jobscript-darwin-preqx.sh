#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-1 on a mac running Darwin
#

EXEC=../../test_execs/preqx-nlev30/preqx-nlev30                         # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-default.nl                          # launch simulation
