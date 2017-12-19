.. _optimizing-processor-layout:

============================
Optimizing Processor Layout
============================

Load balancing is the practice of specifying processor layout for a given model configuration
(compset, grid, and so on) to optimize throughput and efficiency. For a fixed total number of
processors, the goal of optimization to achieve maximum throughput. In contrast, for a given 
configuration across varied processor counts, the purpose is to find several "sweet spots" where 
the model is minimally idle, cost is relatively low, and the throughput is relatively high.
 
As with most models, increasing total processors normally results in both increased throughput 
and increased cost. 
If models scaled linearly, the cost would remain constant across different processor counts, 
but models generally don't scale linearly and the cost increases as processor count increases.

Performing a load-balancing exercise on a proposed model run before undertaking a long production run is recommended practice.

CIME experimental cases have significant flexibility with respect to the layout of components 
across different hardware processors. There are eight unique models (atm, lnd, rof, ocn, ice, 
glc, wav, cpl) that are managed independently by the CIME driver, each with a unique MPI communicator. 
In addition, the driver runs on the union of all processors and controls the sequencing and hardware partitioning.

See :ref:`How to customize the PE layout<changing-the-pe-layout>` for a detailed discussion of how to set processor layouts.

.. _model-timing-data:

Model timing data
------------------

In order to perform a load-balancing exercise, you must first be aware of the different types 
of timing information produced by a model run. 

Every model run produces a summary timing output file. The file is placed in
**$CASEROOT/timing/ccsm_timing.$CASE.$date**, where $date is a datestamp that CIME sets at runtime.
The following describes the most important parts of a timing file:

- CCSM TIMING PROFILE is the first section in the timing output. It summarizes general timing information for the run. The total run time and cost are given in several metrics to facilitate analysis and comparisons with other runs. These metrics includ pe-hrs per simulated year (cost), simulated years per wall day (thoughput), seconds, and seconds per model day. The total run time for each component and the time for initialization of the model also are provided. These times are the aggregate over the total run and do not take into account any temporal or processor load imbalances.

- DRIVER TIMING FLOWCHART is the second section in the timing output. It provides timing information for the driver in sequential order and indicates which processors are involved in the cost. Finally, the timings for the coupler are broken out at the bottom of the timing output file.

- Another file in the timing directory, **ccsm_timing_stats.$date**, summarizes the minimum and maximum of all the model timers.

- Another stream of useful timing information is the **cpl.log.$date** file that every run produces. It contains the run time for each model day during the run and is output during the run. You can search for ``tStamp`` in the cpl.log file to see the information, which is useful for tracking down temporal variability in cost due to inherent model variability or to hardware. The model daily cost generally is pretty constant unless I/O is written intermittently, such as at the end of the month.

Using model timing data
------------------------

Load balancing requires you to consider a number of factors, such as which components are run; their absolute and relative resolution; cost, scaling and processor count sweet spots for each component; and internal load imbalance within a component.
 
It is often best to load balance a system with all significant run-time I/O turned off because it occurs infrequently, typically just one timestep per month. It is best treated as a separate cost as it can otherwise bias interpretation of the overall balance. 
Also, the use of OpenMP threading in some or all of the components is dependent on the hardware/OS support as well as whether the system supports running all MPI and mixed MPI/OpenMP on overlapping processors for different components. 

Finally, decide whether components should run sequentially, concurrently, or in some combination.
 
Typically, a series of short test runs with the desired production configuration can establish a reasonable load balance setup for the production job. The timing output can be used to compare test runs to help determine the optimal load balance.

Changing the pe layout of the model has NO IMPACT on the scientific results. The basic order of operations and calling sequence are hardwired into the driver and do not change with the pe layout. that doesn't change when the pe layout is changed. 
There are some constraints on the ability of either the CESM or E3SM fully active configuration fully concurrent. For example, the atmosphere model always run sequentially with the ice and land models for scientific reasons. As a result, running the atmosphere concurrently with the ice and land will result in idle processors at some point in the timestepping sequence.

For more information about how the driver is implemented, see (Craig, A.P., Vertenstein, M., Jacob, R., 2012: A new flexible coupler for earth system modeling developed for CCSM4 and CESM1.0. International Journal of High Performance Computing Applications, 26, 31-42, 10.1177/1094342011428141). 

**One approach to load balancing**

Carry out 20-day model runs with restarts and history turned off in order to find the layout that has the best load balance for the targeted number of processors. This provides a reasonable performance estimate for the production run for most of the runtime.
 
Treat the end-of-month history and end-of-run restart I/O as a separate cost. 

To set up this test configuration, create your production case, and then edit **env_run.xml** to set ``STOP_OPTION`` to ndays, ``STOP_N`` to 20, and ``RESTART_OPTION`` to never.

Seasonal variation and spin-up costs can change performance over time, so even after a production run has started, review the timing output occasionally to see if any layout changes might improve throughput or decrease cost.

In determining an optimal load balance for a specific configuration, two pieces of information are useful.

- Which components are most expensive.

- How individual components scale. Do they run faster with all MPI or mixed MPI/OpenMP decomposition strategies? What are their optimal decompositions at each processor count? If the cost and scaling of the components are unknown, several short tests with arbitrary component pe counts can help establish component scaling and sweet spots.

**Determining an optimal load balance**

- Start with the most expensive component and a fixed optimal processor count and decomposition for that component.

- Test the systems, varying the sequencing/concurrency of the components and the pe counts of the other components.

- Identify a few potential load balance configurations, then run each a few times to establish run-to-run variability and determine the best layout.

In all cases, review the component run times in the timing output file for both overall throughput and independent component timings. Identify idle processors by considering the component concurrency in conjunction with the component timing.

In general, a few component layout options are most reasonable:

- fully sequential,

- fully sequential except the ocean running concurrently,

- fully concurrent except the atmosphere running sequentially with the ice, rof, and land components.

Finally, run on a subset of the atmosphere processors, either sequentially or concurrently with the land and ice. 

The concurrency is limited in part by hardwired sequencing in the driver. The sequencing is set by scientific constraints, although there may be some addition flexibility with respect to concurrency when running with mixed active and data models.

**Some general rules for finding optimal configurations**

- Make sure you have set a processor layout where each hardware processor is assigned to at least one component. There is rarely a reason to have completely idle processors.

- Make sure your cheapest components keep up with your most expensive components. In other words, a component that runs on 1024 processors should not be waiting on a component running on 16 processors.

- Before running the job, make sure the batch queue settings in the **$CASE.run** script are set correctly for your run. Review the account numbers, queue names and time limits. The ideal time limit, queue and run length are dependent on each other and on the current model throughput.

- Take full advantage of the hardware resources. If you are charged by the 32-way node, you might as well target a total processor count that is a multiple of 32.

- Keep a single component on a single node, if possible, to minimize internal component communication cost.

- Assume that hardware performance can vary due to contention on the interconnect, file systems, or other areas. If you are unsure of a timing result, run cases multiple times.


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

