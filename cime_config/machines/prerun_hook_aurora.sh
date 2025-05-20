#!/bin/bash

LOG_FILE="pre_run_hook_script.log"
echo "Pre-run hook...." | tee -a "$LOG_FILE"
echo "============= Pre-run HOOK for DAOS Support on Aurora =========" | tee -a "$LOG_FILE"
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

echo "=========== LAUNCHING DFUSE ===================================" | tee -a "$LOG_FILE"
#launch-dfuse.sh $DAOS_POOL:$DAOS_CONT | tee -a "$LOG_FILE"
# start dfuse on each compute node
# mount each container on the system

BINDIR=/soft/daos/bin
mountpt="/tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}"
CLUSH_CMD="clush --hostfile=$PBS_NODEFILE -f 208 -B "

${CLUSH_CMD} "mkdir -p ${mountpt} && ${BINDIR}/start-dfuse.sh \
   --pool ${DAOS_POOL} \
   --cont ${DAOS_CONT} \
   -m ${mountpt} \
   --disable-caching \
   --disable-wb-cache " | tee -a "$LOG_FILE"

echo "=========== Setting DAOS INTERCEPTION LIB =====================" | tee -a "$LOG_FILE"
#export LD_PRELOAD="/usr/lib64/libpil4dfs.so"
export LD_PRELOAD="/soft/perftools/darshan/darshan-3.4.7/lib/libdarshan.so:/opt/aurora/24.347.0/spack/unified/0.9.2/install/linux-sles15-x86_64/oneapi-2025.0.5/hdf5-1.14.5-zrlo32i/lib/libhdf5.so:/opt/aurora/24.347.0/spack/unified/0.9.2/install/linux-sles15-x86_64/oneapi-2025.0.5/parallel-netcdf-1.12.3-cszcp66/lib/libpnetcdf.so:/usr/lib64/libpil4dfs.so"
export DARSHAN_LOGFILE="/tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}/${DAOS_USER}/darshan_log.bin"
echo "=========== ls container ======================================" | tee -a "$LOG_FILE"
ls -lRt /tmp/${DAOS_POOL}/${DAOS_CONT}_${DAOS_USER}/${DAOS_USER} | tee -a "$LOG_FILE"
echo "==============================================================" | tee -a "$LOG_FILE"
