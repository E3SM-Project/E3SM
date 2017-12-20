#!/bin/bash
#
#   Jobscript for launching dcmip2012 test1-3 on a mac running Darwin
#
# usage: ./jobscript-...

EXEC=../../../test_execs/pese-nlev60/pese-nlev60                        # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
