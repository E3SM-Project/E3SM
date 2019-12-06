#!/bin/bash
#
# PBS batch script to run an MPI application
#
#PBS -N pyEnsSum            
#PBS -A STDD0001        
#PBS -l walltime=00:30:00               
#PBS -q regular          
#PBS -j oe     
#PBS -M haiyingx@ucar.edu
#PBS -l select=1:ncpus=12:mpiprocs=12                 

### setenv OMP_NUM_THREADS 1
### run the executable
mkdir -p /glade/scratch/haiyingx/temp
export TMPDIR=/glade/scratch/haiyingx/temp



### setenv MPI_USE_ARRAY=false

#mpiexec_mpt python pyEnsSum.py --verbose --indir /glade/p/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel --esize 151   --tslice 1  --tag cesm1_3_beta11 --sumfile /glade/scratch/haiyingx/cesm1.3.b11.ne30_ne30.FC5_V6.nc --mpi_enable --mach yellowstone --compset FC5 --res ne30_ne30 --jsonfile empty.json --gmonly

mpiexec_mpt python pyEnsSum.py --indir /glade/p/tdd/asap/verification/cesm2_0_beta07/uf.yellow.intel.f19.F2000-DEV --esize 350  --verbose --tslice 1 --tag cesm2_0_beta07 --compset F2000_DEV --res f19_f19 --sumfile /glade/scratch/haiyingx/cesm2_0_beta07-ab-uf-summary.nc --gmonly --mpi_enable  --mach cheyenne --verbose --jsonfile empty.json

