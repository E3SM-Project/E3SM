#!/usr/bin/env bash
caseroot=$1
file=$2
cd $caseroot
if hash ncks 2>/dev/null; then
 ncks --glb compset=`./xmlquery --value COMPSET` --glb grid=`./xmlquery --value GRID` --glb testtype="uf" --glb compiler=`./xmlquery --value COMPILER` --glb machineid=`./xmlquery --value MACH`  $file $file.tmp
  mv $file.tmp $file
else
  echo "This script requires the ncks tool"
fi
