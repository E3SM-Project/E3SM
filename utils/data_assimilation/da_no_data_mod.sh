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
  cycle=$1
  shift
  total_cycles=$1
  shift
  if [ -n ${CASEROOT} ]; then
    echo "caseroot: ${CASEROOT}"
  else
    echo "ERROR: CASEROOT not defined"
    errcode=$(( errcode + 1 ))
  fi
  echo "cycle: ${cycle}"
  echo "total_cycles: ${total_cycles}"
fi

exit $errcode
