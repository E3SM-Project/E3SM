#! /bin/tcsh -f
#BSUB -n 47
#BSUB -R "span[ptile=4]"
#BSUB -q regular
#BSUB -N
#BSUB -o ens.%J.stdout
#BSUB -e ens.%J.stdout
#BSUB -a poe
#BSUB -x
#BSUB -J pyEnsSum
#BSUB -W 00:20
#BSUB -P STDD0002

#mpirun.lsf python pyEnsSum.py --indir /glade/u/tdd/asap/verification/cesm1_3_beta06/intel.yellowstone.151 --esize 151  --verbose --tslice 1  --tag cesm1_3_beta06 --sumfile /glade/scratch/haiyingx/intel.151.beta06.nc --jsonfile beta06_ens_excluded_varlist.json --gmonly --mpi_enable 
mpirun.lsf python pyEnsSum.py --indir /glade/p/work/jshollen/validation/cesm1_4_beta06/ensemble --esize 151  --verbose --tslice 1  --tag cesm1_3_beta06 --sumfile /glade/scratch/haiyingx/intel.151.beta06.nc --jsonfile beta06_ens_excluded_varlist_a.json --gmonly --mpi_enable 
