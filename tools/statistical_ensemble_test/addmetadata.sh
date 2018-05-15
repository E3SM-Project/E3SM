#!/usr/bin/env bash
#
# Adds metadata to netcdf statistical ensemble test files.
#
Args=("$@")
i=0
while [ $i -le ${#Args[@]} ]; do
  case ${Args[$i]} in
    --caseroot )
      i=$((i+1))
      caseroot=${Args[$i]}
      if [ ! -d ${caseroot} ]; then
        echo "ERROR: caseroot not found: $caseroot"
        exit 2
      fi
      if [ ! -f ${caseroot}/xmlquery ]; then
	echo "ERROR: Directory $caseroot does not appear to be a cesm case directory"
	exit 3
      fi
    ;;
    --histfile )
      i=$((i+1))
      histfile=${Args[$i]}
      if [ ! -f ${histfile} ]; then
	echo "ERROR: file not found $histfile"
        exit 4
      fi
    ;;
  esac
  i=$((i+1))
done

cd $caseroot

if hash ncks 2>/dev/null; then
 ncks --glb compset=`./xmlquery --value COMPSET` --glb grid=`./xmlquery --value GRID` --glb testtype="uf" --glb compiler=`./xmlquery --value COMPILER` --glb machineid=`./xmlquery --value MACH`  --glb model_version=`./xmlquery --value MODEL_VERSION`   $histfile $histfile.tmp
  mv $histfile.tmp $histfile
else
  echo "This script requires the ncks tool"
fi
