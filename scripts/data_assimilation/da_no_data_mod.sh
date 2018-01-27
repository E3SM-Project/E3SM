#! /bin/bash

##############################################################################
###
### A stub data assimilation script that prints out information but makes
### no modifications to model data.
### Tests using this script should be BFB with a non-data assimilation run
###
##############################################################################

errcode=0
if [ $# -ne 2 ]; then
  echo "ERROR: Wrong number of arguments, $# (should be 2)"
  errcode=$(( errcode + 1 ))
else
  caseroot=$1
  cycle=$2
  echo "caseroot: ${caseroot}"
  echo "cycle: ${cycle}"
  cd ${caseroot}
  res=$?
  if [ $res -ne 0 ]; then
    echo "ERROR: Unable to cd to caseroot, ${caseroot}"
    errcode=$(( errcode + 1 ))
  else
    dval="`./xmlquery --value COMPSET 2> /dev/null | grep DWAV`"
    ./xmlchange DATA_ASSIMILATION_WAV=TRUE
    res=$?
    if [ $res -ne 0 ]; then
      echo "ERROR: Unable to change DATA_ASSIMILATION_WAV to TRUE"
      errcode=$(( errcode + 1 ))
    fi
    if [ $cycle -gt 0 -a -n "${dval}" ]; then
      rundir="`./xmlquery --value RUNDIR`"
      if [ -n "${rundir}" -a -d "${rundir}" ]; then
        cd ${rundir}
        res=$?
        if [ $res -ne 0 ]; then
          echo "ERROR: Unable to cd to caseroot, ${caseroot}"
          errcode=$(( errcode + 1 ))
        else
          # Check the latest log file for a resume signal
          lfile="`ls -t wav.log.* | head -1`"
          if [ -n "${lfile}" -a -f "${lfile}" ]; then
            if [ -n "`echo ${lfile} | grep 'gz$'`" ]; then
              sig="`zcat ${lfile} | grep 'Resume signal, TRUE' 2> /dev/null`"
            else
              sig="`cat ${lfile} | grep 'Resume signal, TRUE' 2> /dev/null`"
            fi
            if [ -n "${sig}" ]; then
              echo "Post-DA resume signal found for cycle ${cycle}"
            else
              echo "ERROR: No post-DA resume signal found for cycle ${cycle}"
              errcode=$(( errcode + 1 ))
            fi
          else
            echo "ERROR: Unable to find wav log file in `pwd -P`"
            errcode=$(( errcode + 1 ))
          fi
        fi
      else
        echo "ERROR: RUNDIR (${rundir}) is not a valid directory"
        errcode=$(( errcode + 1 ))
      fi
    fi
  fi
fi

exit $errcode
