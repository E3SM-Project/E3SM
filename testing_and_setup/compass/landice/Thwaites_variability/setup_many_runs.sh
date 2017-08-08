 #!/bin/bash
'''
bash script to set up a suite of Thwaites variability runs
'''

# ===========================

# Manually set/check these variables before running!

SETUP=1
SUBMIT_RUN=1
SUBMIT_RUN_RESTART=0

BASE_DIR=`pwd`   # Path in which all the versions will be set up 
TEMPLATE_DIR=`pwd`  # Path where the test case infrastructure test was set up

amplitudes='200 500'
periods='02 20 70'
phases="0.0"

#  ==========================



nlfile="namelist.landice"

if [ $SETUP = 1 ]; then
   # First make a 'blank' template copy of the directory
   cp -r $TEMPLATE_DIR/setup_model $BASE_DIR/run_template
   cd $BASE_DIR/run_template
   # remove some unneeded garbage
   rm ./make_graph_file.py ./metis ./setup_model.py
   # remove symlink to albany file
   cp --remove-destination `readlink albany_input.xml` albany_input.xml
   cp --remove-destination `readlink slurm.wolf.run` slurm.wolf.run
fi

# now create all the instances in a flat dir structure
cd $BASE_DIR
for amp in $amplitudes; do
   for per in $periods; do
      for pha in $phases; do
          # create the dir
          dirname=amp${amp}_per${per}_pha${pha}

          if [ $SETUP = 1 ]; then
             echo Setting up:  $dirname
             cp -r $BASE_DIR/run_template $BASE_DIR/$dirname
   
             # update the nl settings
             cd $BASE_DIR/$dirname
             sed -i.SEDBACKUP "s/config_basal_mass_bal_seroussi_amplitude.*/config_basal_mass_bal_seroussi_amplitude = $amp/" $nlfile
             sed -i.SEDBACKUP "s/config_basal_mass_bal_seroussi_period.*/config_basal_mass_bal_seroussi_period = $per/" $nlfile
             sed -i.SEDBACKUP "s/config_basal_mass_bal_seroussi_phase.*/config_basal_mass_bal_seroussi_phase = $pha/" $nlfile
             rm $nlfile.SEDBACKUP
             sed -i.SEDBACKUP "s/^#SBATCH --job-name.*/#SBATCH --job-name=$dirname/" slurm.wolf.run
             rm slurm.wolf.run.SEDBACKUP
          fi

          if [ $SUBMIT_RUN = 1 ]; then
             cd ${BASE_DIR}/${dirname}  # redundant if also setting up, but needed otherwise

             if [ $SUBMIT_RUN_RESTART = 1 ]; then
                sed -i.SEDBACKUP "s/config_do_restart.*/config_do_restart = .true./" $nlfile
                sed -i.SEDBACKUP "s/config_do_restart.*/config_start_time = 'file'/" $nlfile
             fi

             echo submitting job for run $dirname
             sbatch slurm.wolf.run
          fi

          cd $BASE_DIR  # not needed, but seems safer to do 
      done
   done
done

cd $BASE_DIR

