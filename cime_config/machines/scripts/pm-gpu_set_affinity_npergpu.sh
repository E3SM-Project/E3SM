#!/bin/bash

# Runtime launcher wrapper designed to enforce round-robin GPU affinity for high-density MPI jobs (e.g., MPN=64 on a 4-GPU node).
# It ensures optimal resource sharing and prevents device contention by partitioning MPI ranks into subgroups per device.

# example with mpn=tasks_per_node=64 or 64 MPI's per node:
#+------------------+--------------------------+--------------+
#| Local Rank Range | Logic: (Rank / 16) % 4   | Assigned GPU |
#+------------------+--------------------------+--------------+
#| 00 - 15          | 0 / 16 ... 15 / 16 = 0   | dev0         |
#| 16 - 31          | 16 / 16 ... 31 / 16 = 1  | dev1         |
#| 32 - 47          | 32 / 16 ... 47 / 16 = 2  | dev2         |
#| 48 - 63          | 48 / 16 ... 63 / 16 = 3  | dev3         |
#+------------------+--------------------------+--------------+

# Get total MPI tasks per node from first argument
tasks_per_node=$1

# Dynamically detect the number of GPUs on this node
num_gpus=$(nvidia-smi -L | wc -l)
#num_gpus=4

# Calculate how many tasks share each GPU
# If 64 tasks and 4 GPUs, tasks_per_gpu = 16
tasks_per_gpu=$(( ${tasks_per_node} / ${num_gpus} ))

# Use 0 if SLURM_LOCALID is not set
local_id=${SLURM_LOCALID:-0}

# Assign GPU based on Local Rank
# The modulo (%) handles edge cases if tasks_per_node isn't perfectly divisible
gpu=$(( (${local_id} / ${tasks_per_gpu}) % ${num_gpus} ))

export CUDA_VISIBLE_DEVICES=$gpu

#printf '?RANK= %s LOCAL_RANK= %s gpu= %s?\n' ${SLURM_PROCID} ${SLURM_LOCALID} ${gpu}
#echo "num_gpus=${num_gpus} Rank ${SLURM_PROCID} (Local ${SLURM_LOCALID}) assigned to GPU ${gpu}"

# Clean up arguments and launch the application
shift
exec "$@"
