#!/bin/bash
echo "============== Post-config hook ======================"
echo "============= Loading mods for DAOS on Aurora ============="
module use /soft/modulefiles
module load daos

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

echo "================ DAOS CLEANING UP ===================="
clean-dfuse.sh ${DAOS_POOL}:${DAOS_CONT}
fusermount3 -u /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}

echo "======================================================"
