#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ $1 = "CPPFLAGS" ]; then
  exit
elif [ $1 = "LDFLAGS" ]; then
  exit
elif [ $1 = "README" ]; then
  echo "This is a dummy script that runs no profiling tool."
  exit
elif [ $1 = "POSTPROCESS" ]; then
  exit
elif [ $1 = "PROFILE" ]; then
  echo "0-$( echo `./xmlquery --value TOTALPES`-1 | bc) `./xmlquery --value EXEROOT`/e3sm.exe" >> run/mpmd.${LID}.conf  
fi
