.. _faq-basic-example:

A Basic Example
---------------

This specifies all the steps necessary to create, set up, build, and run a case. The following assumes that ``$CIMEROOT`` is the root directory of CIME.

1. Create a new case named EXAMPLE_CASE in the ~/cesm directory. Use an 1850 control compset at 1-degree resolution on yellowstone.
::

   > cd $CIME/scripts
   > ./create_newcase -case ~/EXAMPLE_CASE -compset B1850_CN -res f09_g16

2. Go to the $CASEROOT directory. Edit env_mach_pes.xml if a different pe-layout is desired first. Then set up and build the case.
::

   > cd ~/EXAMPLE_CASE
   > ./case.setup
   > ./case.build

3. Go back to the case directory, set the job to run 12 model months and increase the job wallclock time
::

   > cd ~/EXAMPLE_CASE
   > xmlchange STOP_OPTION=nmonths
   > xmlchange STOP_N=12
   > xmlchange JOB_WALLCLOCK_TIME=06:00
   > ./case.submit 

4. Make sure the run succeeded. Look for the following line at the end of the cpl.log file in your run directory.:
::

   (seq_mct_drv): ===============       SUCCESSFUL TERMINATION OF CPL7-CESM ===============

5. Set it to resubmit itself 10 times so that it will run a total of 11 years (including the initial year), and resubmit the case. 
(Note that a resubmit will automatically change the run to be a continuation run).:
::

   > xmlchange RESUBMIT=10
   > case.submit
