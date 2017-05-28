.. _optimizing-processor-layout:

============================
Optimizing Processor Layout
============================

Load balancing refers to the optimization of the processor layout for a given model configuration (compset, grid, etc) such that the cost and throughput will be optimal. 
Optimal is a somewhat subjective thing. 
For a fixed total number of processors, it means achieving the maximum throughput. 
For a given configuration across varied processor counts, it means finding several "sweet spots" where the model is minimally idle, the cost is relatively low, and the throughput is relatively high. 
As with most models, increasing total processors normally results in both increased throughput and increased cost. 
If models scaled linearly, the cost would remain constant across different processor counts, but generally, models don't scale linearly and cost increases with increasing processor count. 
It is strongly recommended that a user perform a load-balancing exercise on their proposed model run before undertaking a long production run.

CIME experimental cases have significant flexibility with respect to the layout of components across different hardware processors. 
In general, there are eight unique models (atm, lnd, rof, ocn, ice, glc, wav, cpl) that are managed independently by the CIME driver, each with a unique MPI communicator. 
In addition, the driver runs on the union of all processors and controls the sequencing and hardware partitioning.

See :ref:`now to cusomize the PE layout<changing-the-pe-layout>` for a detailed discussion of how to set processor layouts.

.. _model-timing-data:

Model timing data
------------------

In order to perform a load balancing exercise, you must first be aware of the different types of timing information produced by every model run. 

A summary timing output file is produced after every model run. This file is placed in ``$CASEROOT/timing/ccsm_timing.$CASE.$date``, where $date is a datestamp set by CIME at runtime, and contains a summary of various information. 
The following provides a description of the most important parts of a timing file.

The first section in the timing output, CCSM TIMING PROFILE, summarizes general timing information for the run. 
The total run time and cost is given in several metrics including pe-hrs per simulated year (cost), simulated years per wall day (thoughput), seconds, and seconds per model day. 
This provides general summary information quickly in several units for analysis and comparison with other runs. 
The total run time for each component is also provided, as is the time for initialization of the model. 
These times are the aggregate over the total run and do not take into account any temporal or processor load imbalances.

The second section in the timing output, "DRIVER TIMING FLOWCHART", provides timing information for the driver in sequential order and indicates which processors are involved in the cost. Finally, the timings for the coupler are broken out at the bottom of the timing output file.

Separately, there is another file in the timing directory, ccsm_timing_stats.$date that accompanies the above timing summary. 
This second file provides a summary of the minimum and maximum of all the model timers.

There is one other stream of useful timing information in the cpl.log.$date file that is produced for every run. 
The cpl.log file contains the run time for each model day during the model run. 
This diagnostic is output as the model runs. 
You can search for tStamp in the cpl.log file to see this information. 
This timing information is useful for tracking down temporal variability in model cost either due to inherent model variability cost (I/O, spin-up, seasonal, etc) or possibly due to variability due to hardware. 
The model daily cost is generally pretty constant unless I/O is written intermittently such as at the end of the month.

Using model timing data
------------------------

In practice, load-balancing requires a number of considerations such as which components are run, their absolute and relative resolution; cost, scaling and processor count sweet-spots for each component; and internal load imbalance within a component. 
It is often best to load balance the system with all significant run-time I/O turned off because this occurs very infrequently, typically one timestep per month, and is best treated as a separate cost as it can bias interpretation of the overall model load balance. 
Also, the use of OpenMP threading in some or all of the components is dependent on the hardware/OS support as well as whether the system supports running all MPI and mixed MPI/OpenMP on overlapping processors for different components. 
A final point is deciding whether components should run sequentially, concurrently, or some combination of the two with each other. 
Typically, a series of short test runs is done with the desired production configuration to establish a reasonable load balance setup for the production job. 
The timing output can be used to compare test runs to help determine the optimal load balance.

Changing the pe layout of the model has NO IMPACT on the scientific results. The basic order of operations and calling sequence is hardwired into the driver and that doesn't change when the pe layout is changed. 
There are some constraints on the ability of either the CESM or ACME fully active configuraiton fully concurrent. 
In particular, the atmosphere model always run sequentially with the ice and land for scientific reasons. 
As a result, running the atmosphere concurrently with the ice and land will result in idle processors in these components at some point in the timestepping sequence. 

For more information about how the driver is implemented, see (Craig, A.P., Vertenstein, M., Jacob, R., 2012: A new flexible coupler for earth system modeling developed for CCSM4 and CESM1.0. International Journal of High Performance Computing Applications, 26, 31-42, 10.1177/1094342011428141). 

In general, we normally carry out 20-day model runs with restarts and history turned off in order to find the layout that has the best load balance for the targeted number of processors. 
This provides a reasonable performance estimate for the production run for most of the runtime. 
The end of month history and end of run restart I/O is treated as a separate cost from the load balance perspective. 
To set up this test configuration, create your production case, and then edit env_run.xml and set STOP_OPTION to ndays, STOP_N to 20, and RESTART_OPTION to never. 
Seasonal variation and spin-up costs can change performance over time, so even after a production run has started, it's worthwhile to occasionally review the timing output to see whether any changes might be made to the layout to improve throughput or decrease cost.

In determining an optimal load balance for a specific configuration, two pieces of information are useful.

- Determine which component or components are most expensive.

- Understand the scaling of the individual components, whether they run faster with all MPI or mixed MPI/OpenMP decomposition strategies, and their optimal decompositions at each processor count. If the cost and scaling of the components are unknown, several short tests can be carried out with arbitrary component pe counts just to establish component scaling and sweet spots.

One method for determining an optimal load balance is as follows

- start with the most expensive component and a fixed optimal processor count and decomposition for that component

- test the systems, varying the sequencing/concurrency of the components and the pe counts of the other components

- identify a few best potential load balance configurations and then run each a few times to establish run-to-run variability and to try to statistically establish the faster layout

In all cases, the component run times in the timing output file can be reviewed for both overall throughput and independent component timings. Using the timing output, idle processors can be identified by considering the component concurrency in conjunction with the component timing.

In general, there are only a few reasonable component layout options.

- fully sequential

- fully sequential except the ocean running concurrently

- fully concurrent except the atmosphere run sequentially with the ice, rof, and land components

- finally, it makes best sense for the coupler to run on a subset of the atmosphere processors and that can be sequentially or concurrently run with the land and ice

The concurrency is limited in part by the hardwired sequencing in the driver. This sequencing is set by scientific constraints, although there may be some addition flexibility with respect to concurrency when running with mixed active and data models.

There are some general rules for finding optimal configurations:

- Make sure you have set a processor layout where each hardware processor is assigned to at least one component. There is rarely a reason to have completely idle processors in your layout.

- Make sure your cheapest components keep up with your most expensive components. In other words, a component that runs on 1024 processors should not be waiting on a component running on 16 processors.

- Before running the job, make sure the batch queue settings in the $CASE.run script are set correctly for the specific run being targetted. The account numbers, queue names, time limits should be reviewed. The ideal time limit, queues, and run length are all dependent on each other and on the current model throughput.

- Make sure you are taking full advantage of the hardware resources. If you are charged by the 32-way node, you might as well target a total processor count that is a multiple of 32.

- If possible, keep a single component on a single node. That usually minimizes internal component communication cost. That's obviously not possible if running on more processors than the size of a node.

- And always assume the hardware performance could have variations due to contention on the interconnect, file systems, or other areas. If unsure of a timing result, run cases multiple times.


Setting the time limits
-----------------------
In looking at the ccsm_timing.$CASE.$datestamp files for "Model Throughput", output like the following will be found:

```
Overall Metrics:
Model Cost: 327.14 pe-hrs/simulated_year (scale= 0.50)
Model Throughput: 4.70 simulated_years/day
```

The model throughput is the estimated number of model years that you can run in a wallclock day. Based on this, you can maximize $CASE.run queue limit and change $STOP_OPTION and $STOP_N in ``env_run.xml``. For example, say a model's throughput is 4.7 simulated_years/day. On yellowstone(??), the maximum runtime limit is 6 hours. 4.7 model years/24 hours * 6 hours = 1.17 years. On the massively parallel computers, there is always some variability in how long it will take a job to run. On some machines, you may need to leave as much as 20% buffer time in your run to guarantee that jobs finish reliably before the time limit. For that reason we will set our model to run only one model year/job. Continuing to assume that the run is on yellowstone, in ``$CASE.yellowstone.run set``:

```
#BSUB -W 6:00
```

and ``xmlchange`` should be invoked as follows in ``CASEROOT``:

```
./xmlchange STOP_OPTION=nyears
./xmlchange STOP_N=1 
./xmlchange REST_OPTION=nyears
./xmlchange REST_N=1 
```
