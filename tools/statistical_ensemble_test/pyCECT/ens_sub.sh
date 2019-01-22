#! /bin/tcsh -f
#BSUB -n 23 
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


#generate a new ens sum of cesm1_3_beta11
#have excluded list
#mpirun.lsf python pyEnsSum.py --indir /glade/p/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel --esize 151   --tslice 1  --tag cesm1_3_beta11 --sumfile /glade/scratch/haiyingx/cesm1.3.b11.ne30_ne30.FC5_V6.nc --mpi_enable --mach yellowstone --compset FC5 --res ne30_ne30 --jsonfile ens_excluded_varlist.json --gmonly

#generate a new ens sum of cesm1_3_beta11
#have included list
mpirun.lsf python pyEnsSum.py --indir /glade/p/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel --esize 151   --tslice 1  --tag cesm1_3_beta11 --sumfile /glade/scratch/haiyingx/cesm1.3.b11.ne30_ne30.FC5_V6.nc --mpi_enable --mach yellowstone --compset FC5 --res ne30_ne30 --jsonfile included_varlist.json --gmonly
