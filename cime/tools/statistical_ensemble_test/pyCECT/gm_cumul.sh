#! /bin/tcsh -f
#BSUB -n 450 
#BSUB -R "span[ptile=4]"
#BSUB -q small
#BSUB -N
#BSUB -o cumul.%J.stdout
#BSUB -e cumul.%J.stdout
#BSUB -a poe
#BSUB -x
#BSUB -J pyEnsSum
#BSUB -W 01:00
#BSUB -P STDD0002

#mpirun.lsf python pyEnsSum.py --indir /glade/scratch/abaker/monthly_allcompilers/ --esize 150  --verbose --tslice 1  --tag cesm1_3_beta06 --sumfile /glade/scratch/haiyingx/gm.cumul.150.nc --jsonfile beta06_ens_excluded_varlist.json --cumul --mpi_enable --regx '(pgi(.)*-(01|02|03|04|05|06|07|08|09|10|11|12))' --compset FC --mach ys --res ne30-ne30
mpirun.lsf python pyEnsSum.py --indir /glade/scratch/abaker/monthly_allcompilers/ --esize 450  --verbose --tslice 1  --tag cesm1_3_beta06 --sumfile /glade/scratch/haiyingx/gm.cumul.450.nc --jsonfile beta06_ens_excluded_varlist.json --cumul --mpi_enable  --compset FC --mach ys --res ne30-ne30 --startMon 1 --endMon 2 --fIndex 151
