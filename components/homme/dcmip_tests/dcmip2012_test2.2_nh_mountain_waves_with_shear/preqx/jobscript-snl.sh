#!/bin/tcsh 


set NCPU = 40 

# NH model
set EXEC = ../../../test_execs/theta-nlev60/theta-nlev60    
set namelist = namelist-default.nl


# hydrostatic theta
#set EXEC = ../../../test_execs/theta-nlev60/theta-nlev60    
#set namelist = namelist-default.nl


# hydrostatic preqx
#set EXEC =../../../test_execs/preqx-nlev60-interp/preqx-nlev60-interp        # set name of executable
#set namelist = namelist-default.nl

\cp -f $namelist input.nl
date
mpirun -np $NCPU $EXEC < input.nl
date
