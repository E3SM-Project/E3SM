.. _timers:

===================
Timers and timing
===================

CIME includes a copy of the General Purpose Timing Library (GPTL) and timers are placed throughout the CIME driver.  CIME-driven models typically
also have GPTL timers in their code and very detailed timing information can be obtained.

.. _model-timing-data:

Model timing data
------------------

Every model run produces a summary timing output file. The file is placed in
**$CASEROOT/timing/ccsm_timing.$CASE.$date**, where $date is a datestamp that CIME sets at runtime.
The following describes the most important parts of a timing file:

- CCSM TIMING PROFILE is the first section in the timing output. It summarizes general timing information for the run. The total run time and cost are given in several metrics to facilitate analysis and comparisons with other runs. These metrics includ pe-hrs per simulated year (cost), simulated years per wall day (thoughput), seconds, and seconds per model day. The total run time for each component and the time for initialization of the model also are provided. These times are the aggregate over the total run and do not take into account any temporal or processor load imbalances.

- DRIVER TIMING FLOWCHART is the second section in the timing output. It provides timing information for the driver in sequential order and indicates which processors are involved in the cost. Finally, the timings for the coupler are broken out at the bottom of the timing output file.

- Another file in the timing directory, **ccsm_timing_stats.$date**, summarizes the minimum and maximum of all the model timers.

- Another stream of useful timing information is the **cpl.log.$date** file that every run produces. It contains the run time for each model day during the run and is output during the run. You can search for ``tStamp`` in the cpl.log file to see the information, which is useful for tracking down temporal variability in cost due to inherent model variability or to hardware. The model daily cost generally is pretty constant unless I/O is written intermittently, such as at the end of the month.


Controlling timers
------------------

.. todo:: Add info on how to control timers

Setting the time limits
-----------------------

When you look at the **ccsm_timing.$CASE.$datestamp** file for "Model Throughput", you will find output like this:
 ::

  Overall Metrics:
  Model Cost: 327.14 pe-hrs/simulated_year (scale= 0.50)
  Model Throughput: 4.70 simulated_years/day

The model throughput is the estimated number of model years that you can run in a wallclock day. Based on this, you can maximize your **$CASE.run** queue limit and change ``$STOP_OPTION`` and ``$STOP_N`` in **env_run.xml**. 

For example, say a model's throughput is 4.7 simulated_years/day, and the maximum runtime limit on your machine is 12 hours. 4.7 model years/24 hours * 12 hours = 2.35 years. On the massively parallel computers, there is always some variability in how long it will take a job to run. On some machines, you may need to leave as much as 20% buffer time in your run to guarantee that jobs finish reliably before the time limit. For that reason, set your model to run only one model year/job. In this example, set your wallclock at 12 hours and invoke **xmlchange** in ``CASEROOT`` as shown here:
 ::

  >./xmlchange STOP_OPTION=nyears
  >./xmlchange STOP_N=1 
  >./xmlchange REST_OPTION=nyears
  >./xmlchange REST_N=1 

