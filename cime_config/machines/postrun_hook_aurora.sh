#!/bin/bash

# This script assumes that DAOS pre-run hook is already
# executed (DFUSE is running)

LOG_FILE="post_run_hook_script.log"

echo "=========== Post-run HOOK for DAOS on Aurora ===============" | tee -a "$LOG_FILE"
if [ -z "${DAOS_POOL}" ]; then
  DAOS_POOL="E3SM_Dec"
else
  DAOS_POOL="${DAOS_POOL}"
fi

if [ -z "${USER}" ]; then
  DAOS_USER="e3sm"
else
  DAOS_USER="${USER}"
fi

if [ -z "${DAOS_CONT}" ]; then
  DAOS_CONT="e3sm"
else
  DAOS_CONT="${DAOS_CONT}"
fi

echo "=========== DAOS OBJECT INFORMATION ===========================" | tee -a "$LOG_FILE"
echo "Querying DAOS Pool : $DAOS_POOL" | tee -a "$LOG_FILE"
daos pool query $DAOS_POOL | tee -a "$LOG_FILE"
echo "Querying DAOS Container, $DAOS_CONT, in $DAOS_POOL" | tee -a "$LOG_FILE"
daos cont query $DAOS_POOL $DAOS_CONT | tee -a "$LOG_FILE"

echo "=========== ls container ===================================" | tee -a "$LOG_FILE"
ls -lRt /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}/${DAOS_USER} | tee -a "$LOG_FILE"
echo "=========== CLEANUP DFUSE =====================================" | tee -a "$LOG_FILE"
clean-dfuse.sh $DAOS_POOL:$DAOS_CONT | tee -a "$LOG_FILE"
echo "===============================================================" | tee -a "$LOG_FILE"
