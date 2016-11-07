#!/bin/bash
#
#   Jobscript for launching dcmip2012 test2-0 on a mac running Darwin
#
# usage: ./jobscript-...

EXEC=../../test_execs/preqx-nlev30/preqx-nlev30-native   # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl            # launch simulation
