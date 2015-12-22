#!/bin/bash -l

#PBS -A ${HOMME_PROJID}
#PBS -o ${TEST_NAME}.stdout.${PBS_JOBID}
#PBS -e ${TEST_NAME}.stderr.${PBS_JOBID}
#PBS -N ${TEST_NAME}
#PBS -l nodes=1
#PBS -l walltime=2:00:00

