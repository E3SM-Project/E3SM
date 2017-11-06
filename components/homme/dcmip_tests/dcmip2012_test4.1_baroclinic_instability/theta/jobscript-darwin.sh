#!/bin/bash
#
#   Jobscript for launching dcmip2012 test4-1 on a mac running Darwin
#
# usage: ./jobscript-...

EXEC=../../../test_execs/theta-nlev30/theta-nlev30
openmpiexec -n 6 $EXEC < ./nh-x1000.nl                           # launch simulation
