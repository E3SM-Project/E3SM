#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 2 on a mac running Darwin
#
# usage: ./jobscript-...

# 4dg resolution
EXEC=../../../test_execs/theta-nlev30/theta-nlev30
cp ./namelist-r400.nl input.nl
openmpiexec -n 6 $EXEC < input.nl                          
