#!/bin/bash
#SBATCH --job-name=parallel_python
#SBATCH --nodes=1              # Ensure single node
#SBATCH --ntasks=1             # Run one task (master process)
#SBATCH --exclusive            # Ensure exclusive access to the node (uses all cpus)
##SBATCH --mem=16G              # Request memory (adjust as needed)
#SBATCH --time=00:30:00
#SBATCH --account e3sm
#SBATCH --partition=regular
#SBATCH --constraint=cpu

# Tell Python/OpenMP to use the requested number of threads
export OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE}
export MKL_NUM_THREADS=${SLURM_CPUS_ON_NODE}

# Run the landgen package
# This command will run the __main__.py file inside the landgen package/directory
python -m landgen