.. _running-a-case:

***************
Running a Case
***************

.. _case-submit:

========================
Calling **case.submit**
========================

The script `case.submit <../Tools_user/case.submit.html>`_  will submit your run to the batch queueing system on your machine.
If you do not have a batch queueing system, `case.submit <../Tools_user/case.submit.html>`_ will start the job interactively, given that you have a proper MPI environment defined.
Running `case.submit <../Tools_user/case.submit.html>`_ is the **ONLY** way you should start a job.

To see the options to `case.submit <../Tools_user/case.submit.html>`_, issue the command
::

   > ./case.submit --help

A good way to see what `case.submit <../Tools_user/case.submit.html>`_ will do, is to first call `preview_run <../Tools_user/preview_run.html>`_
::

   > ./preview_run

which will output the environment for your run along with the batch submit and mpirun commands.
As an example, on the NCAR machine, cheyenne, for an A compset at the f19_g17_rx1 resolution, the following is output from `preview_run <../Tools_user/preview_run.html>`_:
::

   CASE INFO:
      nodes: 1
      total tasks: 36
      tasks per node: 36
      thread count: 1

   BATCH INFO:
      FOR JOB: case.run
   ENV:
      module command is /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/libexec/lmod python purge
      module command is /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/libexec/lmod python load ncarenv/1.2 intel/17.0.1 esmf_libs mkl esmf-7.0.0-defio-mpi-O mpt/2.16 netcdf-mpi/4.5.0 pnetcdf/1.9.0 ncarcompilers/0.4.1
      Setting Environment OMP_STACKSIZE=256M
      Setting Environment TMPDIR=/glade/scratch/mvertens
      Setting Environment MPI_TYPE_DEPTH=16
   SUBMIT CMD:
      qsub    -q regular -l walltime=12:00:00 -A P93300606 .case.run

   FOR JOB: case.st_archive
      ENV:
         module command is /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/libexec/lmod python purge
         module command is /glade/u/apps/ch/opt/lmod/7.5.3/lmod/lmod/libexec/lmod python load ncarenv/1.2 intel/17.0.1 esmf_libs mkl esmf-7.0.0-defio-mpi-O mpt/2.16 netcdf-mpi/4.5.0 pnetcdf/1.9.0 ncarcompilers/0.4.1
         Setting Environment OMP_STACKSIZE=256M
         Setting Environment TMPDIR=/glade/scratch/mvertens
         Setting Environment MPI_TYPE_DEPTH=16
         Setting Environment TMPDIR=/glade/scratch/mvertens
         Setting Environment MPI_USE_ARRAY=false
   SUBMIT CMD:
      qsub    -q share -l walltime=0:20:00 -A P93300606  -W depend=afterok:0 case.st_archive

   MPIRUN:
      mpiexec_mpt  -np 36 -p "%g:"  omplace -tm open64  /glade/scratch/mvertens/jim/bld/cesm.exe  >> cesm.log.$LID 2>&1

Each of the above sections is defined in the various **$CASEROOT** xml files and the associated variables can be modified using the
`xmlchange <../Tools_user/xmlchange.html>`_ command (or in the case of tasks and threads, this can also be done with the `pelayout <../Tools_user/pelayout.html>`_ command).

- The PE layout is set by the xml variables **NTASKS**, **NTHRDS** and **ROOTPE**. To see the exact settings for each component, issue the command
  ::

     ./xmlquery NTASKS,NTHRDS,ROOTPE

  To change all of the **NTASKS** settings to say 30 and all of the **NTHRDS** to 4, you can call
  ::

     ./xmlchange NTASKS=30,NTHRDS=4

  To change JUST the ATM NTASKS to 8, you can call
  ::

     ./xmlchange NTASKS_ATM=8

- Submit parameters are set by the xml variables in the file **env_batch.xml**. This file is special in certain xml variables can appear in more than one group.
  NOTE: The groups are the list of jobs that are submittable for a case.
  Normally, the minimum set of groups are  **case.run** and **case.st_archive**.
  We will illustrate how to change an xml variable in **env_batch.xml** using the xml variable ``JOB_WALLCLOCK_TIME``.

  - To change ``JOB_WALLCLOCK_TIME`` for all groups to 2 hours for cheyenne, use
    ::

       ./xmlchange JOB_WALLCLOCK_TIME=02:00:00

  - To change ``JOB_WALLCLOCK_TIME`` to 20 minutes for cheyenne for just **case.run**, use
    ::

       ./xmlchange JOB_WALLCLOCK_TIME=00:20:00 --subgroup case.run


Before you submit the case using `case.submit <../Tools_user/case.submit.html>`_, make sure the batch queue variables are set correctly for your run
In particular, make sure that you have appropriate account numbers (``PROJECT``), time limits (``JOB_WALLCLOCK_TIME``), and queue (``JOB_QUEUE``).

Also modify **$CASEROOT/env_run.xml** for your case using **xmlchange**.

Once you have executed `case.setup <../Tools_user/case.setup.html>`_ and `case.build <../Tools_user/case.build.html>`_ , call `case.submit <../Tools_user/case.submit.html>`_
to submit the run to your machine's batch queue system.
::

   > cd $CASEROOT
   > ./case.submit

---------------------------------
Result of running case.submit
---------------------------------

When called, the `case.submit <../Tools_user/case.submit.html>`_ script will:

- Load the necessary environment.

- Confirm that locked files are consistent with the current xml files.

- Run `preview_namelist <../Tools_user/preview_namelist.html>`_, which in turn will run each component's **cime_config/buildnml** script.

- Run :ref:`check_input_data<input_data>` to verify that the required data are present.

- Submit the job to the batch queue. which in turn will run the `case.run <../Tools_user/case.run.html>`_ script.

Upon successful completion of the run, `case.run <../Tools_user/case.run.html>`_  will:

- Put timing information in **$CASEROOT/timing**.
  See :ref:`model timing data<model-timing-data>` for details.

- Submit the short-term archiver script `case.st_archive <../Tools_user/case.st_archive.html>`_  to the batch queue if ``$DOUT_S`` is TRUE.
  Short-term archiving will copy and move component history, log, diagnostic, and restart files from ``$RUNDIR`` to the short-term archive directory ``$DOUT_S_ROOT``.

- Resubmit `case.run <../Tools_user/case.run.html>`_ if ``$RESUBMIT`` > 0.


---------------------------------
Monitoring case job statuses
---------------------------------

The **$CASEROOT/CaseStatus** file contains a log of all the job states and `xmlchange <../Tools_user/xmlchange.html>`_ commands in chronological order.
Below is an example of status messages:
::

  2017-02-14 15:29:50: case.setup starting
  ---------------------------------------------------
  2017-02-14 15:29:54: case.setup success
  ---------------------------------------------------
  2017-02-14 15:30:58: xmlchange success <command> ./xmlchange STOP_N=2,STOP_OPTION=nmonths  </command>
  ---------------------------------------------------
  2017-02-14 15:31:26: xmlchange success <command> ./xmlchange STOP_N=1  </command>
  ---------------------------------------------------
  2017-02-14 15:33:51: case.build starting
  ---------------------------------------------------
  2017-02-14 15:53:34: case.build success
  ---------------------------------------------------
  2017-02-14 16:02:35: case.run starting
  ---------------------------------------------------
  2017-02-14 16:20:31: case.run success
  ---------------------------------------------------
  2017-02-14 16:20:45: st_archive starting
  ---------------------------------------------------
  2017-02-14 16:20:58: st_archive success
  ---------------------------------------------------

.. note::
  After a successful first run, set the **env_run.xml** variable ``$CONTINUE_RUN`` to ``TRUE`` before resubmitting or the job will not
  progress.

  You may also need to modify the **env_run.xml** variables
  ``$STOP_OPTION``, ``$STOP_N`` and/or ``$STOP_DATE`` as well as
  ``$REST_OPTION``, ``$REST_N`` and/or ``$REST_DATE``, and ``$RESUBMIT``
  before resubmitting.

See the :ref:`basic example<basic_example>` for a complete example of how to run a case.

---------------------------------
Troubleshooting a job that fails
---------------------------------

There are several places to look for information if a job fails.
Start with the **STDOUT** and **STDERR** file(s) in **$CASEROOT**.
If you don't find an obvious error message there, the
**$RUNDIR/$model.log.$datestamp** files will probably give you a
hint.

First, check **cpl.log.$datestamp**, which will often tell you
*when* the model failed. Then check the rest of the component log
files. See :ref:`troubleshooting run-time problems<troubleshooting>` for more information.

.. _input_data:

====================================================
Input data
====================================================

The **check_input_data** script determines if the required data files
for your case exist on local disk in the appropriate subdirectory of
``$DIN_LOC_ROOT``. It automatically downloads missing data required for your simulation.

.. note:: It is recommended that users on a given system share a common ``$DIN_LOC_ROOT`` directory to avoid duplication on
	  disk of large amounts of input data. You may need to talk to your system administrator in order to set this up.

The required input data sets needed for each component are found in the
**$CASEROOT/Buildconf** directory. These files are generated by a call
to **preview_namlists** and are in turn created by each component's
**buildnml** script. For example, for compsets consisting only of data
models (i.e. ``A`` compsets), the following files are created:
::

   cpl.input_data_list
   datm.input_data_list
   dice.input_data_list
   docn.input_data_list
   drof.input_data_list

You can independently verify the presence of the required data by
using the following commands:
::

   > cd $CASEROOT
   > ./check_input_data --help
   > ./check_input_data

If data sets are missing, obtain them from the input data server(s) via the commands:
::

   > cd $CASEROOT
   > ./check_input_data --download

``check_input_data`` is automatically called by the case control
system, when the case is built and submitted.  So manual usage of this
script is optional.

-----------------------------------
Distributed Input Data Repositories
-----------------------------------

CIME has the ability to utilize multiple input data repositories, with
potentially different protocols.  The repositories are defined in the
file **$CIMEROOT/config/$model/config_inputdata.xml**.  The currently
supported server protocols are: ``gridftp``, ``subversion``, ``ftp`` and
``wget``. These protocols may not all be supported on your machine,
depending on software configuration.

.. note:: You now have the ability to create your own input data
          repository and add it to the **config_inputdata.xml**. This
          will permit you to easily collaborate by sharing your
          required inputdata with others.


.. _controlling-start-stop-restart:

====================================================
Starting, Stopping and Restarting a Run
====================================================

The file **env_run.xml** contains variables that may be modified at
initialization or any time during the course of a model run. Among
other features, the variables comprise coupler namelist settings for
the model stop time, restart frequency, coupler history frequency, and
a flag to determine if the run should be flagged as a continuation run.

At a minimum, you will need to set the variables ``$STOP_OPTION`` and
``$STOP_N``. Other driver namelist settings then will have consistent and
reasonable default values. The default settings guarantee that
restart files are produced at the end of the model run.

By default, the stop time settings are:
::

  STOP_OPTION = ndays
  STOP_N = 5
  STOP_DATE = -999

The default settings are appropriate only for initial testing. Before
starting a longer run, update the stop times based on the case
throughput and batch queue limits. For example, if the model runs 5
model years/day, set ``RESUBMIT=30, STOP_OPTION= nyears, and STOP_N=
5``. The model will then run in five-year increments and stop after
30 submissions.

.. _run-type-init:

---------------------------------------------------
Run-type initialization
---------------------------------------------------

The case initialization type is set using the ``$RUN_TYPE`` variable in
**env_run.xml**. A CIME run can be initialized in one of three ways:

``startup``

  In a startup run (the default), all components are initialized using
  baseline states. These states are set independently by each component
  and can include the use of restart files, initial  files, external
  observed data files, or internal initialization (that is, a "cold start").
  In a startup run, the coupler sends the start date to the components
  at initialization. In addition, the coupler does not need an input data file.
  In a startup initialization, the ocean model does not start until the second
  ocean coupling step.

``branch``

  In a branch run, all components are initialized using a consistent
  set of restart files from a previous run (determined by the
  ``$RUN_REFCASE`` and ``$RUN_REFDATE`` variables in **env_run.xml**).
  The case name generally is changed for a branch run, but it
  does not have to be. In a branch run, the ``$RUN_STARTDATE`` setting is
  ignored because the model components obtain the start date from
  their restart data sets. Therefore, the start date cannot be changed
  for a branch run. This is the same mechanism that is used for
  performing a restart run (where ``$CONTINUE_RUN`` is set to TRUE in
  the **env_run.xml** file). Branch runs typically are used when
  sensitivity or parameter studies are required, or when settings for
  history file output streams need to be modified while still
  maintaining bit-for-bit reproducibility. Under this scenario, the
  new case is able to produce an exact bit-for-bit restart in the same
  manner as a continuation run if no source code or component namelist
  inputs are modified. All models use restart files to perform this
  type of run. ``$RUN_REFCASE`` and ``$RUN_REFDATE`` are required for
  branch runs. To set up a branch run, locate the restart tar file or
  restart directory for ``$RUN_REFCASE`` and ``$RUN_REFDATE`` from a
  previous run, then place those files in the ``$RUNDIR``  directory.
  See :ref:`setting up a branch run<setting-up-a-branch-run>`.

``hybrid``

  A hybrid run is initialized like a startup but it uses
  initialization data sets from a previous case. It is similar
  to a branch run with relaxed restart constraints.
  A hybrid run allows users to bring together
  combinations of initial/restart files from a previous case
  (specified by ``$RUN_REFCASE``) at a given model output date
  (specified by ``$RUN_REFDATE``). Unlike a branch run, the starting
  date of a hybrid run (specified by ``$RUN_STARTDATE``) can be
  modified relative to the reference case. In a hybrid run, the model
  does not continue in a bit-for-bit fashion with respect to the
  reference case.  The resulting climate, however, should be
  continuous provided that no model source code or namelists are
  changed in the hybrid run. In a hybrid initialization, the ocean
  model does not start until the second ocean coupling step, and the
  coupler does a "cold start" without a restart file.

The variable ``$RUN_TYPE`` determines the initialization type. This
setting is only important for the initial production run when
the ``$CONTINUE_RUN`` variable is set to FALSE. After the initial
run, the ``$CONTINUE_RUN`` variable is set to TRUE, and the model
restarts exactly using input files in a case, date, and bit-for-bit
continuous fashion.

The variable ``$RUN_STARTDATE`` is the start date (in yyyy-mm-dd format)
for either a startup run or a hybrid run. If the run is targeted to be
a hybrid or branch run, you must specify values for ``$RUN_REFCASE`` and
``$RUN_REFDATE``.

.. _starting_from_a_refcase:

----------------------------------------
Starting from a reference case (REFCASE)
----------------------------------------

There are several xml variables that control how either a branch or a hybrid case can start up from another case.
The initial/restart files needed to start up a run from another case are required to be in $EXEROOT.
The xml variable ``$GET_REFCASE`` is a flag that if set will automatically prestaging the refcase restart data.

- If ``$GET_REFCASE`` is ``TRUE``, then the the values set by ``$RUN_REFDIR``, ``$RUN_REFCASE``, ``$RUN_REFDATE`` and  ``$RUN_TOD`` are
  used to prestage the data by symbolic links to the appropriate path.

  The location of the necessary data to start up from another case is controlled by the xml variable ``$RUN_REFDIR``.

  - If ``$RUN_REFDIR`` is an absolute pathname, then it is expected that initial/restart files needed to start up a model run are in ``$RUN_REFDIR``.

  - If ``$RUN_REFDIR`` is a relative pathname, then it is expected that initial/restart files needed to start up a model run are in a path relative to ``$DIN_LOC_ROOT`` with the absolute pathname  ``$DIN_LOC_ROOT/$RUN_REFDIR/$RUN_REFCASE/$RUN_REFDATE``.

  - If ``$RUN_REFDIR`` is a relative pathname AND is not available in ``$DIN_LOC_ROOT`` then CIME will attempt to download the data from the input data repositories.


- If ``$GET_REFCASE`` is ``FALSE`` then the data is assumed to already exist in ``$EXEROOT``.

.. _controlling-output-data:

=========================
Controlling output data
=========================

During a model run, each model component produces its own output
data sets in ``$RUNDIR`` consisting of history, initial, restart, diagnostics, output
log and rpointer files. Component history files and restart files are
in netCDF format. Restart files are used to either restart the same
model or to serve as initial conditions for other model cases. The
rpointer files are ascii text files that list the component history and
restart files that are required for restart.

Archiving (referred to as short-term archiving here) is the phase of a model run when output data are
moved from ``$RUNDIR`` to a local disk area (short-term archiving).
It has no impact on the production run except to clean up disk space
in the ``$RUNDIR`` which can help manage user disk quotas.

Several variables in **env_run.xml** control the behavior of
short-term archiving. This is an example of how to control the
data output flow with two variable settings:
::

  DOUT_S = TRUE
  DOUT_S_ROOT = /$SCRATCH/$user/$CASE/archive


The first setting above is the default, so short-term archiving is enabled. The second sets where to move files at the end of a successful run.

Also:

- All output data is initially written to ``$RUNDIR``.

- Unless you explicitly turn off short-term archiving, files are
  moved to ``$DOUT_S_ROOT`` at the end of a successful model run.

- Users generally should turn off short-term archiving when developing new code.

Standard output generated from each component is saved in ``$RUNDIR``
in a  *log file*. Each time the model is run, a single coordinated datestamp
is incorporated into the filename of each output log file.
The run script generates the datestamp in the form YYMMDD-hhmmss, indicating
the year, month, day, hour, minute and second that the run began
(ocn.log.040526-082714, for example).

By default, each component also periodically writes history files
(usually monthly) in netCDF format and also writes netCDF or binary
restart files in the ``$RUNDIR`` directory. The history and log files
are controlled independently by each component. History output control
(for example, output fields and frequency) is set in each component's namelists.

The raw history data does not lend itself well to easy time-series
analysis. For example, CAM writes one or more large netCDF history
file(s) at each requested output period. While this behavior is
optimal for model execution, it makes it difficult to analyze time
series of individual variables without having to access the entire
data volume. Thus, the raw data from major model integrations usually
is post-processed into more user-friendly configurations, such as
single files containing long time-series of each output fields, and
made available to the community.

For CESM, refer to the `CESM2 Output Filename Conventions
<http://www.cesm.ucar.edu/models/cesm2.0/cesm/filename_conventions_cesm.html>`_
for a description of output data filenames.

.. _restarting-a-run:

======================
Restarting a run
======================

Active components (and some data components) write restart files
at intervals that are dictated by the driver via the setting of the
``$REST_OPTION`` and ``$REST_N`` variables in **env_run.xml**. Restart
files allow the model to stop and then start again with bit-for-bit
exact capability; the model output is exactly the same as if the model
had not stopped. The driver coordinates the writing of restart
files as well as the time evolution of the model.

Runs that are initialized as branch or hybrid runs require
restart/initial files from previous model runs (as specified by the
variables ``$RUN_REFCASE`` and ``$RUN_REFDATE``). Pre-stage these
iles to the case ``$RUNDIR`` (normally ``$EXEROOT/run``)
before the model run starts. Normally this is done by copying the contents
of the relevant **$RUN_REFCASE/rest/$RUN_REFDATE.00000** directory.

Whenever a component writes a restart file, it also writes a restart
pointer file in the format **rpointer.$component**. Upon a restart, each
component reads the pointer file to determine which file to read in
order to continue the run. These are examples of pointer files created
for a component set using full active model components.
::

  - rpointer.atm
  - rpointer.drv
  - rpointer.ice
  - rpointer.lnd
  - rpointer.rof
  - rpointer.cism
  - rpointer.ocn.ovf
  - rpointer.ocn.restart


If short-term archiving is turned on, the model archives the
component restart data sets and pointer files into
**$DOUT_S_ROOT/rest/yyyy-mm-dd-sssss**, where yyyy-mm-dd-sssss is the
model date at the time of the restart. (See `below for more details
<http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.)

---------------------------------
Backing up to a previous restart
---------------------------------

If a run encounters problems and crashes, you will normally have to
back up to a previous restart. If short-term archiving is enabled,
find the latest **$DOUT_S_ROOT/rest/yyyy-mm-dd-ssss/** directory
and copy its contents into your run directory (``$RUNDIR``).

Make sure that the new restart pointer files overwrite older files in
in ``$RUNDIR`` or the job may not restart in the correct place. You can
then continue the run using the new restarts.

Occasionally, when a run has problems restarting, it is because the
pointer and restart files are out of sync. The pointer files
are text files that can be edited to match the correct dates
of the restart and history files. All of the restart files should
have the same date.

============================
Archiving model output data
============================

The output data flow from a successful run depends on whether or not
short-term archiving is enabled, as it is by default.

-------------
No archiving
-------------

If no short-term archiving is performed, model output data remains
remain in the run directory as specified by ``$RUNDIR``.

---------------------
Short-term archiving
---------------------

If short-term archiving is enabled, component output files are moved
to the short-term archiving area on local disk, as specified by
``$DOUT_S_ROOT``. The directory normally is **$EXEROOT/../archive/$CASE.**
and has the following directory structure: ::

   rest/yyyy-mm-dd-sssss/
   logs/
   atm/hist/
   cpl/hist
   glc/hist
   ice/hist
   lnd/hist
   ocn/hist
   rof/hist
   wav/hist
   ....

The **logs/** subdirectory contains component log files that were
created during the run. Log files are also copied to the short-term
archiving directory and therefore are available for long-term archiving.

The **rest/** subdirectory contains a subset of directories that each contains
a *consistent* set of restart files, initial files and rpointer
files. Each subdirectory has a unique name corresponding to the model
year, month, day and seconds into the day when the files were created.
The contents of any restart directory can be used to create a branch run
or a hybrid run or to back up to a previous restart date.

---------------------
Long-term archiving
---------------------

Users may choose to follow their institution's preferred method for long-term
archiving of model output. Previous releases of CESM provided an external
long-term archiver tool that supported mass tape storage and HPSS systems.
However, with the industry migration away from tape archives, it is no longer
feasible for CIME to support all the possible archival schemes available.

================================================
Data Assimilation and other External Processing
================================================

CIME provides a capability to run a task on the compute nodes either
before or after the model run.  CIME also provides a data assimilation
capability which will cycle the model and then a user defined task for
a user determined number of cycles.


-------------------------
Pre and Post run scripts
-------------------------

Variables ``PRERUN_SCRIPT`` and ``POSTRUN_SCRIPT`` can each be used to name
a script which should be exectuted immediately prior starting or
following completion of the CESM executable within the batch
environment.  The script is expected to be found in the case directory
and will recieve one argument which is the full path to that
directory.  If the script is written in python and contains a
subroutine with the same name as the script, it will be called as a
subroutine rather than as an external shell script.

-------------------------
Data Assimilation scripts
-------------------------

Variables ``DATA_ASSIMILATION``, ``DATA_ASSIMILATION_SCRIPT``, and
``DATA_ASSIMILATION_CYCLES`` may also be used to externally control
model evolution.  If ``DATA_ASSIMILATION`` is true after the model
completes the ``DATA_ASSIMILATION_SCRIPT`` will be run and then the
model will be started again ``DATA_ASSIMILATION_CYCLES`` times.  The
script is expected to be found in the case directory and will recieve
two arguments, the full path to that directory and the cycle number.
If the script is written in python and contains a subroutine with the
same name as the script, it will be called as a subroutine rather than
as an external shell script.

..: A simple example pre run script.

::

   #!/usr/bin/env python
   import sys
   from CIME.case import Case

   def myprerun(caseroot):
       with Case(caseroot) as case:
            print ("rundir is ",case.get_value("RUNDIR"))

    if __name__ == "__main__":
      caseroot = sys.argv[1]
      myprerun(caseroot)
