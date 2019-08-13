#!/bin/tcsh

#PBS -A NIOW0001
#PBS -N ensSumPop
#PBS -q regular
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=0:15:00
#PBS -j oe
#PBS -M abaker@ucar.edu

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR


mpiexec_mpt python pyCECT.py --popens --sumfile  /glade/p/cisl/iowa/pop_verification/cesm2_0_beta10/pop.ens.sum.cesm2.0.b10.nc --indir /glade/p/cisl/iowa/pop_verification/cesm2_0_beta10/testcases/C96 --jsonfile pop_ensemble.json --input_glob C96.pop.000.pop.h.0001-12

