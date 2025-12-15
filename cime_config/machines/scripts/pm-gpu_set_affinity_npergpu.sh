#!/bin/bash
#num_gpus=$(nvidia-smi -L | wc -l)
tasks_per_node=$1
tasks_per_gpu=$(( ${tasks_per_node} / 4 ))
gpu=$(( (${SLURM_LOCALID} / ${tasks_per_gpu}) % 4 ))
export CUDA_VISIBLE_DEVICES=$gpu
#echo “RANK= ${SLURM_PROCID} LOCAL_RANK= ${SLURM_LOCALID} gpu= ${gpu}”
shift
"$@"
