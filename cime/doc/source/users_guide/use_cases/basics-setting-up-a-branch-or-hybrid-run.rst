.. _setting-up-a-branch-run:

Setting up a branch or hybrid run
---------------------------------

A branch or hybrid run uses initialization data from a previous run. Here is an example in which a valid load-balanced scenario is assumed.

The first step in setting up a branch or hybrid run is to create a new case. 
::

   > cd $CIMEROOT/scripts
   > create_newcase --case ~/EXAMPLE_CASEp --compset B1850 --res f09_g16
   > cd ~/EXAMPLE_CASEp


For a branch run, modify **env_run.xml** to branch from **EXAMPLE_CASE** at year 0001-02-01.
::

   > xmlchange RUN_TYPE=branch
   > xmlchange RUN_REFCASE=EXAMPLE_CASE
   > xmlchange RUN_REFDATE=0001-02-01

For a hybrid run, modify **env_run.xml** to start up from **EXAMPLE_CASE** at year 0001-02-01.
::

   > xmlchange RUN_TYPE=hybrid
   > xmlchange RUN_REFCASE=EXAMPLE_CASE
   > xmlchange RUN_REFDATE=0001-02-01

For a branch run, your **env_run.xml** file for **EXAMPLE_CASEp** should be identical to the file for **EXAMPLE_CASE** except for the ``$RUN_TYPE`` setting. 
Also, modifications introduced into **user_nl_** files in **EXAMPLE_CASE** should be reintroduced in **EXAMPLE_CASEp**.

Next, set up and build your case executable.
::

   > ./case.setup
   > ./case.build

Pre-stage the necessary restart/initial data in ``$RUNDIR``. Assume for this example that it was created in the **/rest/0001-02-01-00000** directory shown here:
::

   > cd $RUNDIR
   > cp /user/archive/EXAMPLE_CASE/rest/0001-02-01-00000/* . 

It is assumed that you already have a valid load-balanced scenario. 
Go back to the case directory, set the job to run 12 model months, and submit the job.
::

   > cd ~/EXAMPLE_CASEp
   > xmlchange STOP_OPTION=nmonths
   > xmlchange STOP_N=12
   > xmlchange JOB_WALLCLOCK_TIME=06:00
   > case.submit

Make sure the run succeeded. Look for the following line at the end of the **cpl.log** file in your run directory.
::

   (seq_mct_drv): ===============       SUCCESSFUL TERMINATION OF CPL7-CCSM ===============

Change the run to a continuation run. Set it to resubmit itself 10 times so it will run a total of 11 years (including the initial year), then resubmit the case.
::

   > xmlchange CONTINUE_RUN=TRUE
   > xmlchange RESUMIT=10
   > case.submit

