#!/bin/bash -x

echo "============= ls for DAOS on Aurora ============="
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

echo "=========== DAOS OBJECT INFORMATION ==========================="
echo "Querying DAOS Pool : $DAOS_POOL"
daos pool query $DAOS_POOL
echo "Querying DAOS Container, $DAOS_CONT, in $DAOS_POOL"
daos cont query $DAOS_POOL $DAOS_CONT

echo "=========== LAUNCHING DFUSE ==================================="
mkdir -p /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}
start-dfuse.sh -m /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER} --pool ${DAOS_POOL} --cont ${DAOS_CONT}
mount | grep dfuse
mount | grep e3sm

echo "=============== listing contents of directory ================"
ls -lRt /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}/${DAOS_USER}
#ls -lt /tmp/${DAOS_POOL}/${DAOS_CONT}/run
echo "=============== Rm'ing contents of directory ================"
if [ "$#" -ne 0 ]; then
  cp $@
fi

echo "================ CLEANING UP ======================"
clean-dfuse.sh ${DAOS_POOL}:${DAOS_CONT}
fusermount3 -u /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}
