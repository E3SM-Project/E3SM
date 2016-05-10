#!/bin/bash
#BSUB -n 12 
#BSUB -q regular
#BSUB -R "span[ptile=1]"
#BSUB -a poe
#BSUB -N
#BSUB -x
#BSUB -o pop.%J.stdout
#BSUB -e pop.%J.stdout
#BSUB -J pop
#BSUB -P STDD0002
#BSUB -W 3:25
#BSUB -u haiyingx@ucar.edu

export MP_LABELIO=yes;
export MP_COREFILE_FORMAT=lite
mpirun.lsf python pyEnsSumPop.py --verbose --indir /glade/scratch/haiyingx/pop_ensemble_data60/ --sumfile /glade/scratch/haiyingx/pop.ens.sum.2.nc --tslice 0 --nyear 1 --nmonth 12 --npert 80 --jsonfile pop_ensemble.json  --mpi_enable --zscoreonly
