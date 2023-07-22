#!/bin/bash -e
#SBATCH --job-name bubble
##SBATCH -p acme-medium
#SBATCH -p compute,debug
##SBATCH -p debug
##SBATCH --account=condo           
#SBATCH -N 8
#SBATCH --time=15:00:00


source ~/load-modules-sky
#make -j16 theta-l-nlev60-native theta-l-nlev128-native theta-l-nlev20-native

#NMAX, TIMESTEP, OUTFREQ

#for 128 nlev, we need only 8 frames, so running till 700 sec would be ok
totlength=800 #in sec  ; for 128 lev it was 1000 sec/100 sec output frequency
freqnumber=100 #output each 100 sec

dtarray=( 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0001)
#dtarray=(  0.01)
#dtarray=( 0.0001 )

#option=cv
option=eamcpdry


for ddt in ${dtarray[@]} ; do

    nlname=new-${option}${ddt}.nl

    nmax=$( echo "(( $totlength / $ddt ))" | bc)
    ofreq=$( echo "(( $freqnumber / $ddt ))" | bc)

    sed -e s/TIMESTEP/${ddt}/ -e s/NMAX/${nmax}/ -e s/OUTFREQ/${ofreq}/ \
        ${option}TEMPLATE > ${nlname}

    echo "dt is ${ddt}"
    echo "nmax is ${nmax}"
    echo "outfreq is ${ofreq}"

srun -n 512 -K --cpu_bind=cores test_execs/theta-l-nlev128-native/theta-l-nlev128-native < ${nlname} > output-${option}-l128-dt${ddt}
#srun -n 512 -K --cpu_bind=cores test_execs/theta-l-nlev20-native/theta-l-nlev20-native < ${nlname} > output-${option}-dt${ddt}
#srun -n 512 -K --cpu_bind=cores test_execs/theta-l-nlev60-native/theta-l-nlev60-native < ${nlname}
#srun -n 512 -K --cpu_bind=cores test_execs/theta-l-nlev128-native/theta-l-nlev128-native < ${nlname}

done







