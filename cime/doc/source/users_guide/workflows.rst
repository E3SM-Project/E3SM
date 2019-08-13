.. _workflows:

*********
Workflows
*********

Currently, there are three kinds of workflow controls available in CIME.

1. Multiple jobs workflow

   The file, **$CIMEROOT/config/$model/machines/config_batch.xml**, contains a section called ``<batch_jobs>`` which defines the submit script templates and the pre-requisites for running them.
   As an example, in cesm, the default ``<batch_jobs>`` section is the following, with an explanation give by the NOTES additions.

   ::

      <batch_jobs>
         <!-- order matters, with no-batch jobs will be run in the order listed here -->
	 <job name="case.run">
	   <template>template.case.run</template>
	   <prereq>$BUILD_COMPLETE and not $TEST</prereq>
	 </job>
	 <job name="case.test">
	    <template>template.case.test</template>
            <prereq>$BUILD_COMPLETE and $TEST</prereq>
         </job>
         <job name="case.st_archive">
            <template>template.st_archive</template>
            <task_count>1</task_count>
            <walltime>0:20:00</walltime>
            <!-- If DOUT_S is true and case.run (or case.test) exits successfully then run st_archive-->
            <dependency>case.run or case.test</dependency>
            <prereq>$DOUT_S</prereq>
         </job>
      </batch_jobs>

   The elements that can be contained in the ``<batch_jobs>`` section are:

   * <job name="NAME"> : the name of the batch job.

   Batch jobs can contain one or more jobs, and each job has the following elements:

   * <template> : required, the associated template that is used to create the job.
   * <prepreq> : required, the pre-requiste settings of xml variables that determine whether the job will be run.
   * <task_count> : optional, the number of tasks the job should use.
   * <walltime> : optional, the maximum walltime the job should require.
   * <dependency> : optional, the previous job expected to complete before this job can start.

   In the above the scripts for each ``<job>`` entry,
   i.e. ``case.run``, ``case.test`` and ``case.st_archive`` are
   created from their respective templates the corresponding
   ``<template>`` entry. The templates are located in the
   **$CIMEROOT/config/$model/machines** directory.

   Furthermore, in the above ``case.run`` will be run if the xml
   variables ``BUILD_COMPLETE`` is true and ``TEST`` is false.  On the
   other hand, if ``BUILD_COMPLETE`` is true and ``TEST`` is true then
   ``case.test`` will be run.  In addition, if ``DOUT_S`` is true,
   then ``case.st_archive`` will be run **after** either ``case.run``
   or ``case.test`` successfully completes. Also note that
   ``case.st_archive`` is a separate job that only requests a single
   task.

   Users can extend the above workflow with their own custom needs.

2. Pre-run and post-run script

   CIME provides the ability to execute user-defined scripts during
   the execution of ``case.run``. These user-defined scripts can be
   invoked either before and/or after the model is run. The xml variables that controls this capability are:

   * ``PRERUN_SCRIPT``: points to an external script to be run before model execution.

   * ``POSTRUN_SCRIPT``: points to an external script to be run after successful model completion.

   Note, that when these scripts are called, the full processor allocation for the job will be used - even if only 1 processor actually is invoked for the external script.

3. Data assimilation controls

  CIME provides the ability to hook in a data assimilation utility via a set of xml variables:

  * ``DATA_ASSIMILATION_SCRIPT``:  points to an external script to be run **after** model completion

  * ``DATA_ASSIMILATION_CYCLES``: integer that controls the number of data assimilation cycles. The run script
    will loop over these number of data assimilation cycles and for each cycle will run the model and subsequently run the data assimilation script.

  * ``DATA_ASSIMILATION``: if set to TRUE for a given component, then
    a resume signal will be sent to that component at
    initialization. If set, the component will execute special post
    data assimilation logic on initialization.  See the component
    documentation for details. This flag is a bit subtle in that it is a per-component flag, not a model wide flag.

    ::

       To see what the component flags are call
        > ./xmlquery DATA_ASSIMILATION

       The output will look like
     	>   DATA_ASSIMILATION: ['CPL:FALSE', 'ATM:TRUE', 'LND:FALSE', 'ICE:FALSE', 'OCN:FALSE', 'ROF:FALSE', 'GLC:FALSE', 'WAV:FALSE']

	To change the LND value to TRUE issue
	> ./xmlchange DATA_ASSIMILATION_LND=TRUE
