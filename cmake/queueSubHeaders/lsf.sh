#!/bin/bash

#BSUB -a poe
#BSUB -P ${HOMME_PROJID}

#BSUB -q small
#BSUB -W 0:20
#BSUB -x

#BSUB -J ${testName}

#BSUB -o ${testName}.stdout.%J
#BSUB -e ${testName}.stderr.%J

#BSUB -n 16
#BSUB -R "span[ptile=16]" 

