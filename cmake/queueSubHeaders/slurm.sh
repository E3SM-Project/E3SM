#!/bin/bash

#SBATCH -J ${testName} # Job name
#SBATCH -o ${testname}.stdout.%j     # Name of stdout output file(%j expands to jobId)
#SBATCH -e ${testName}.stderr.%j     # Name of stderr output file(%j expands to jobId)
#SBATCH -p normal           # Submit to the 'normal' or 'development' queue
#SBATCH -N 1                # Total number of nodes requested (16 cores/node)
#SBATCH -n 16                # Total number of mpi tasks requested
#SBATCH -t 0:01:00         # Run time (hh:mm:ss) - 1.5 hours

#SBATCH -A ${HOMME_PROJID} # Allocation name to charge job against

