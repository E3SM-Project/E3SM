#!/bin/bash 
#
#   Jobscript for launching DCMIP 2016 Test 2 - Tropical Cyclone on NERSC
#
#SBATCH --job-name=dcmip2016-2
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#SBATCH --qos=debug
#SBATCH --nodes=8
#SBATCH --constraint=cpu

#SBATCH --account=e3sm
#SBATCH --mail-type=all

#SBATCH --time=0-00:30:00

# Process distribution
OMP_NUM_THREADS=1
CORES_PER_NODE=128  # Cores per node on Perlmutter
NCORES=$(( $SLURM_NNODES * $CORES_PER_NODE )) # Physical Cores
CORES_PER_TASK=2    # Perlmutter has hyperthreading enabled, so use 2
NVCORES=$(( CORES_PER_TASK * NCORES ))        # Virtual Cores


# Executable location, module for NCL
EXEC=../../../test_execs/theta-l-nlev30/theta-l-nlev30
module load climate-utils/2023q2

# Function to run the test
function run {
local NTASKS=$1
echo "NTASKS = $NTASKS"
namelist=namelist-$prefix.nl
\cp -f $namelist input.nl
srun -K -c $CORES_PER_TASK -n $NTASKS -N $SLURM_NNODES $EXEC < input.nl
date

# Plot initial conditions for: u, T, th, q, p_nh, geo, and ps
ncl plot-tropical-cyclone-init.ncl
# Plot horizontal cross-section an 10 days for: u, v, T, Q, precl, geo, and ps
ncl plot-horiz-crossx.ncl
# Plot surface pressure and lowest level horizontal wind magnitude over time.
ncl plot-intensity-trace.ncl
# Plot surface pressure map at days 0, 8, 9, and 10.
ncl plot-horiz-ps.ncl

# Save output to run-specific files
\mv -f movies/dcmip2016_test21.nc   movies/${prefix}_dcmip2016_test21.nc

\mv -f init.pdf ${prefix}_init.pdf
\mv -f x-sections.pdf ${prefix}_x-sections.pdf
\mv -f wind.pdf ${prefix}_wind.pdf
\mv -f psmap.pdf ${prefix}_psmap.pdf
}

# Max NTASKS is ne*ne*6, with ne specified in the namelist
MAX_NTASKS=$(( 8 * 8 * 6 ))
prefix=r400  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

MAX_NTASKS=$(( 30 * 30 * 6 ))
prefix=r100  ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

MAX_NTASKS=$(( 60 * 60 * 6 ))
prefix=r50   ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))
prefix=r50-h ; run $(($NVCORES>$(( $CORES_PER_TASK * $MAX_NTASKS ))?MAX_NTASKS:NCORES))

