#!/bin/bash
#
#   Jobscript for launching dcmip2012 test3-1 on a mac running Darwin
#

EXEC=../../../test_execs/pese-nlev20/pese-nlev20                        # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
