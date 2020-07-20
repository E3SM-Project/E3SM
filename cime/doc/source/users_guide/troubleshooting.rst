.. _troubleshooting:

Troubleshooting
===============

Troubleshooting case creation
-----------------------------

Generally, `create_newcase  <../Tools_user/create_newcase.html>`_ errors are reported to the terminal and should provide some guidance about what caused them.

If `create_newcase  <../Tools_user/create_newcase.html>`_ fails on a relatively generic error, first check to make sure the command-line arguments match the interface's specification. See the help text to review usage.
::

   > create_newcase --help

Troubleshooting problems in cime scripts
----------------------------------------

If any of the python-based cime scripts are dying in a mysterious way, more information can be obtained by rerunning the script with the ``--debug`` option.

Troubleshooting job submission
-------------------------------

Most problems associated with submission or launch are site-specific.
The batch and run aspects of the `case.submit  <../Tools_user/case.submit.html>`_ script are created by parsing the variables in **$CASEROOT/env_batch.xml** file.

Take these steps to check for problems:

1. Review the batch submission options in **$CASEROOT/env_batch.xml**. Confirm that they are consistent with the site-specific batch environment, and that the queue names, time limits, and hardware processor request make sense and are consistent with the case.

2. Make sure that `case.submit  <../Tools_user/case.submit.html>`_ uses the correct batch job tool for submitting the `case.submit  <../Tools_user/case.submit.html>`_ script. Depending on the batch environment, it might be **bsub**, **qsub** or another command. Also confirm if a redirection "<" character is required. The information for how **case.submit** submits jobs appears at the end of the standard output stream.

Troubleshooting runtime problems
---------------------------------

To see if a run completed successfully, check the last several lines of the **cpl.log** file for a string like ``SUCCESSFUL TERMINATION``. A successful job also usually copies the log files to the **$CASEROOT/logs** directory.

Check these things first when a job fails:

- Did the model time out?

- Was a disk quota limit hit?

- Did a machine go down?

- Did a file system become full?

If any of those things happened, take appropriate corrective action (see suggestions below) and resubmit the job.

If it is not clear that any of the above caused a case to fail, there are several places to look for error messages.

- Check component log files in your run directory (``$RUNDIR``).
  This directory is set in the **env_run.xml** file.
  Each component writes its own log file, and there should be log files for every component in this format: **cpl.log.yymmdd-hhmmss**.
  Check each log file for an error message, especially at or near the end.

- Check for a standard out and/or standard error file in ``$CASEROOT``.
  The standard out/err file often captures a significant amount of extra model output and also often contains significant system output when a job terminates.
  Useful error messages sometimes are found well above the bottom of a large standard out/err file. Backtrack from the bottom in search of an error message.

- Check for core files in your run directory and review them using an appropriate tool.

- Check any automated email from the job about why a job failed. Some sites' batch schedulers send these.

- Check the archive directory: **$DOUT_S_ROOT/$CASE**.   If a case failed, the log files
  or data may still have been archived.

**Common errors**

One common error is for a job to time out, which often produces minimal error messages.
Review the daily model date stamps in the **cpl.log** file and the timestamps of files in your run directory to deduce the start and stop time of a run.
If the model was running fine, but the wallclock limit was reached, either reduce the run length or increase the wallclock setting.

If the model hangs and then times out, that usually indicates an MPI or file system problem or possibly a model problem. If you suspect an intermittent system problem, try resubmitting the job. Also send a help request to local site consultants to provide them with feedback about system problems and to get help.

Another error that can cause a timeout is a slow or intermittently slow node.
The **cpl.log** file normally outputs the time used for every model simulation day. To review that data, grep the **cpl.log** file for the string ``tStamp`` as shown here:
::

     > grep tStamp cpl.log.* | more

The output looks like this:
::

  tStamp_write: model date = 10120 0 wall clock = 2009-09-28 09:10:46 avg dt = 58.58 dt = 58.18
  tStamp_write: model date = 10121 0 wall clock = 2009-09-28 09:12:32 avg dt = 60.10 dt = 105.90


Review the run times at the end of each line for each model day.
The "avg dt =" is  the average time to simulate a model day and "dt = " is the time needed to simulate the latest model day.

The model date is printed in YYYYMMDD format and the wallclock is the local date and time.
In the example, 10120 is Jan 20, 0001, and the model took 58 seconds to run that day.
The next day, Jan 21, took 105.90 seconds.

A wide variation in the simulation time for typical mid-month model days suggests a system problem. However, there are variations in the cost of the model over time.
For instance, on the last day of every simulated month, the model typically writes netcdf files, which can be a significant intermittent cost.
Also, some model configurations read data mid-month or run physics intermittently at a timestep longer than one day.
In those cases, some variability is expected. The time variation typically is quite erratic and unpredictable if the problem is system performance variability.

Sometimes when a job times out or overflows disk space, the restart files will get mangled.
With the exception of the CAM and CLM history files, all the restart files have consistent sizes.

Compare the restart files against the sizes of a previous restart. If they don't match, remove them and move the previous restart into place before resubmitting the job.
See `Restarting a run <http://esmci.github.io/cime/users_guide/running-a-case.html#restarting-a-run>`_.

It is not uncommon for nodes to fail on HPC systems or for access to large file systems to hang. Before you file a bug report, make sure a case fails consistently in the same place.

**Rerunning with additional debugging information**

There are a few changes you can make to your case to get additional information that aids in debugging:

- Increase the value of the run-time xml variable ``INFO_DBUG``: ``./xmlchange INFO_DBUG=2``.
  This adds more information to the cpl.log file that can be useful if you can't tell what component is aborting the run, or where bad coupling fields are originating.
  (This does NOT require rebuilding.)

- Try rebuilding and rerunning with the build-time xml variable ``DEBUG`` set to ``TRUE``: ``./xmlchange DEBUG=TRUE``.

  - This adds various runtime checks that trap conditions such as out-of-bounds array indexing, divide by 0, and other floating point exceptions (the exact conditions checked depend on flags set in ``config_compilers.xml``).

  - The best way to do this is often to create a new case and run ``./xmlchange DEBUG=TRUE`` before running ``./case.build``.
    However, if it is hard for you to recreate your case, then you can run that xmlchange command from your existing case; then you must run ``./case.build --clean-all`` before rerunning ``./case.build``.

  - Note that the model will run significantly slower in this mode, so this may not be feasible if the model has to run a long time before producing the error.
    (Sometimes it works well to run the model until shortly before the error in non-debug mode, have it write restart files, then restart after rebuilding in debug mode.)
    Also note that answers will change slightly, so if the error arises from a rare condition, then it may not show up in this mode.
