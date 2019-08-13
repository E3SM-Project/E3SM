#!/bin/bash

#MSUB -A ${HOMME_PROJID}
#MSUB -l nodes=1:ppn=12,walltime=00:45:00
#MSUB -q pdebug
#MSUB -d ${TEST_DIR}
#MSUB -j
#MSUB -o ${TEST_NAME}.stdout.${SLURM_JOB_ID}
#MSUB -e ${TEST_NAME}.stderr.${SLURM_JOB_ID}

