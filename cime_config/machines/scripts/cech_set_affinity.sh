#!/bin/bash
# On Cech (NVIDIA GB200 NVL4):
# Grace Socket 0 (CPUs 0-71, GPUs 0-1) connects to ConnectX-7 HCA mlx5_0
# Grace Socket 1 (CPUs 72-143, GPUs 2-3) connects to ConnectX-7 HCA mlx5_3
#
# CPU LPDDR5X memory is on NUMA nodes 0 and 1.
# GPU HBM3e memory is exposed by the OS as movable system RAM on NUMA nodes 2, 10, 18, 26.
# Bind user-space host memory allocations (malloc/mmap/Kokkos::HostSpace) to CPU NUMA
# nodes 0,1 so the Linux kernel reclaims Lustre page cache (FilePages) on nodes 0,1
# instead of spilling host allocations onto GPU HBM3e NUMA nodes (which triggers
# cudaErrorNvlinkUncorrectable during cudaMemcpyAsync). Note: cudaMalloc allocations
# in the NVIDIA kernel driver use __GFP_THISNODE and still allocate on GPU HBM3e nodes.

export CUDA_DEVICE_MAX_CONNECTIONS=1

# Raise soft locked-memory limit (RLIMIT_MEMLOCK) to match the compute node's
# unlimited hard limit in case Slurm propagated a low limit (e.g. 8192 KB) from
# the login node, which otherwise causes UCX ibv_reg_mr() failures on multi-node runs.
ulimit -Sl unlimited 2>/dev/null || true

NUMA_NODE="0,1"
TASKS_PER_NODE=${SLURM_NTASKS_PER_NODE:-${SLURM_TASKS_PER_NODE%%[(,]*}}
if [ -n "${TASKS_PER_NODE}" ] && [ "${TASKS_PER_NODE}" -gt 4 ] 2>/dev/null; then
  gpu=$(( (${SLURM_LOCALID:-0} * 4) / ${TASKS_PER_NODE} ))
  export KOKKOS_DEVICE_ID=${gpu}
  case ${gpu} in
    0|1)
      export UCX_NET_DEVICES=mlx5_0:1
      NUMA_NODE="0"
      ;;
    *)
      export UCX_NET_DEVICES=mlx5_3:1
      NUMA_NODE="1"
      ;;
  esac
elif [ -n "${SLURM_LOCALID}" ]; then
  case ${SLURM_LOCALID} in
    0|1)
      export UCX_NET_DEVICES=mlx5_0:1
      NUMA_NODE="0"
      ;;
    2|3)
      export UCX_NET_DEVICES=mlx5_3:1
      NUMA_NODE="1"
      ;;
    *)
      # Fallback for non-standard rank layouts: bind based on assigned CPU core ID
      cpu=$(taskset -c -p $$ 2>/dev/null | awk '{print $NF}' | cut -d'-' -f1 | cut -d',' -f1)
      if [ -n "$cpu" ] && [ "$cpu" -ge 72 ] 2>/dev/null; then
        export UCX_NET_DEVICES=mlx5_3:1
        NUMA_NODE="1"
      else
        export UCX_NET_DEVICES=mlx5_0:1
        NUMA_NODE="0"
      fi
      ;;
  esac
fi

if command -v numactl >/dev/null 2>&1; then
  exec numactl --membind=${NUMA_NODE} "$@"
else
  exec "$@"
fi
