#!/bin/bash

#MSUB -A ${HOMME_PROJID}
#MSUB -l nodes=1:ppn=12,walltime=00:45:00
#MSUB -q pdebug
#MSUB -d ${testDir}
#MSUB -j
#MSUB -o ${testName}.stdout.${SLURM_JOB_ID}
#MSUB -e ${testName}.stderr.${SLURM_JOB_ID}

