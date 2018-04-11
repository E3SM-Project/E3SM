.. _use-cases-basic-example:

Setting up a multi-year run
----------------------------

This shows all of the steps necessary to do a multi-year simulation.

1. Create a new case named EXAMPLE_CASE in your **$HOME** directory. Use an 1850 control compset at 1-degree resolution (CESM components/resolution).

   ::

   > cd $CIME/scripts
   > ./create_newcase --case ~/EXAMPLE_CASE --compset B1850_CN --res f09_g16

2. Go to the $CASEROOT directory. Edit **env_mach_pes.xml** if you need a different pe-layout. Then set up and build the case as shown.

   ::

   > cd ~/EXAMPLE_CASE
   > ./case.setup
   > ./case.build

3. In your case directory, set the job to run 12 model months, set the wallclock time, and submit the job.

   ::

   > cd ~/EXAMPLE_CASE
   > xmlchange STOP_OPTION=nmonths
   > xmlchange STOP_N=12
   > xmlchange JOB_WALLCLOCK_TIME=06:00
   > ./case.submit 

4. Make sure the run succeeded. Look for the following line at the end of the **cpl.log file** in your run directory.

   ::

   (seq_mct_drv): ===============       SUCCESSFUL TERMINATION OF CPL7-CESM ===============

5. In the same case directory, Set the case to resubmit itself 10 times so it will run a total of 11 years (including the initial year), and resubmit the case. (Note that a resubmit will automatically change the run to be a continuation run).

   ::

   > xmlchange RESUBMIT=10
   > case.submit
