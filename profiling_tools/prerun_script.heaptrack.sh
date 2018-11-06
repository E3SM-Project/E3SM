#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


if [ $1 = "CPPFLAGS" ]; then
  echo "-g"
  exit
elif [ $1 = "LDFLAGS" ]; then
  exit
elif [ $1 = "README" ]; then
  echo "This script runs heaptrack, a memory profiling tool."
  exit
elif [ $1 = "POSTPROCESS" ]; then
  exit
elif [ $1 = "PROFILE" ]; then
  tool=/project/projectdirs/acme/software/heaptrack/knl/bin/heaptrack
  echo "0 $tool `./xmlquery --value EXEROOT`/e3sm.exe" > run/mpmd.${LID}.conf;
  echo "1-$( echo `./xmlquery --value TOTALPES`-1 | bc) `./xmlquery --value EXEROOT`/e3sm.exe" >> run/mpmd.${LID}.conf  
fi

