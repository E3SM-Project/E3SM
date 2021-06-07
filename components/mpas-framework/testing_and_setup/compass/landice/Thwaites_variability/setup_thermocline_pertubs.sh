#!/bin/bash
# bash script to set up a suite of Thwaites variability runs


# ===========================

# Manually set/check these variables before running!

BASE_DIR=`pwd`   # Path in which all the versions will be set up
#TEMPLATE_DIR=/scratch1/scratchdirs/hoffman2/thwaites_variability_ensemble_25pctGLcalc_from_year5_steady_from_bedmap2_D15_kappa11.0_sill663/run_template_base
TEMPLATE_DIR=`pwd`/adjust+000

elevs="-025 -050 -075 -100 -125 -150 -175 -200 -225 -250 -275 -300 +025 +050 +075 +100 +125 +150 +175 +200 +225 +250 +275 +300"

#  ==========================



nlfile="namelist.landice"
# now create all the instances in a flat dir structure

cp ../slurm.edison.run .

for elev in $elevs; do
         cd $BASE_DIR
          
          # build the dir name
          dirname=adjust${elev}

          # === Setting up the run ===
             echo Setting up:  $dirname
             cp -r $TEMPLATE_DIR $dirname

             # update the nl settings
             cd $BASE_DIR/$dirname
             udepth=`grep config_temperature_profile_thermocline_upper_depth $nlfile |cut -d "=" -f 2 |tr -d " " |cut -d "." -f 1`
             echo udepth is $udepth
             udepthnew=`python -c "x=$udepth;y=int('$elev'); print x+y"`
             echo udepthnew is $udepthnew
             sed -i.SEDBACKUP "s/config_temperature_profile_thermocline_upper_depth.*/config_temperature_profile_thermocline_upper_depth = $udepthnew/" $nlfile
             ldepth=`grep config_temperature_profile_thermocline_lower_depth $nlfile |cut -d "=" -f 2 |tr -d " " |cut -d "." -f 1`
             ldepthnew=`python -c "x=$ldepth;y=int('$elev'); print x+y"`
             sed -i.SEDBACKUP "s/config_temperature_profile_thermocline_lower_depth.*/config_temperature_profile_thermocline_lower_depth = $ldepthnew/" $nlfile
             rm $nlfile.SEDBACKUP

             # add this run to the edison bundle
             echo "cd ${BASE_DIR}/${dirname}" >> ../slurm.edison.run
             echo  "srun -n 48 -N 1 ../landice_model " >> ../slurm.edison.run
done


cd $BASE_DIR
