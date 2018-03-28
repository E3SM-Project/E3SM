.. _pesthreads:

==================================
Controlling processors and threads
==================================

Once a compset and resolution for a case has been defined, CIME provides ways to define how many processors and
threads the case will use.


.. _defining-pes:

pe-settings for a case
-------------------------

CIME looks at the xml element ``PES_SPEC_FILE`` in the **config_files.xml** file to determine where
to find the supported out-of-the-box model grids for the target component.

Each component that sets compsets has an associated **config_pes.xml** file that specifies an out-of-the-box pe-layout for those compsets.
The pe-layout might also have dependencies on the model grid and the target machine.
Finally, there might be more than one out-of-the-box pe-layout that could be used for a compset/grid/machine combination: one for a low processor setting and one for a high processor setting.

A typical entry in **config_pes.xml** looks like this:

::

  <grid name="a%T62">
    <mach name="cheyenne">
      <pes pesize="any" compset="DATM%IAF">
      .......
      </pes>
    </mach>
  </grid>

Given the various dependencies, CIME uses an order of precedence to determine the optimal match. This order is as follows:

1. grid match

   CIME first searches the grid nodes for a grid match in **config_grids.xml**.
   The search is based on a regular expression match for the grid longname.
   All nodes that have a grid match are used in the subsequent search. If there is no grid match, all nodes that have ``<grid name="any">`` are used in the subsequent search.


2. machine match

   CIME next uses the list of nodes obtained in the grid match to search for the machine name using the ``<mach>`` nodes. If there is no machine match, then all nodes with ``<machine name="any">`` are used in the subsequent search.


3. pesize and compset match

   CIME next uses the list of nodes obtained in the machine match to search for pesize and compset using the ``<pes>`` nodes. If there is no match, the node with ``<pes pesize="any" compset="any">`` is used.

The **create_newcase** script outputs the matches that are found in determining the best out-of-the-box pe-layout.

Threading control
-------------------------

.. todo:: Add threading control info

SMT
-------------------------

.. todo:: Add SMT info


.. _optimizing-processor-layout:
Optimizing processor layout
----------------------------

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
