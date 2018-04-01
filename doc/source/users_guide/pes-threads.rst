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

Changing default PE layout
--------------------------

Settings in the **env_mach_pes.xml** file determine:

- the number of MPI tasks and OpenMP threads for each component.
- the number of instances of each component.
- the layout of the components across the hardware processors.

Optimizing the throughput and efficiency of a CIME experiment often involves customizing the processor (PE) layout. (See :ref:`load balancing <optimizing-processor-layout>`.)
CIME provides significant flexibility with respect to the layout of components across different hardware processors. In general, the CIME components -- atm, lnd, ocn, and so on -- can run on overlapping or mutually unique processors. While each component is associated with a unique MPI communicator, the CIME driver runs on the union of all processors and controls the sequencing and hardware partitioning.

The component processor layout is determined by the following settings:

- the number of MPI tasks.
- the number of OpenMP threads per task.
- the root MPI task number from the global communicator.
- the maximum number of MPI tasks per node.

The entries in **env_mach_pes.xml** have the following meanings:

.. csv-table:: "Entries in env_mach_pes.xml"
   :header: "xml variable", "description"
   :widths: 25, 75

   "MAX_TASKS_PER_MODE",  "The total number of (MPI tasks) * (OpenMP threads) allowed on a node. This is defined in **config_machines.xml** and therefore given a default setting, but can be user modified."
   "MAX_MPITASKS_PER_NODE", "The maximum number of MPI tasks per node. This is defined in **config_machines.xml** and therefore given a default setting, but can be user modified."
   "NTASKS", "Total number of MPI tasks. A negative value indicates nodes rather than tasks, where MAX_MPITASKS_PER_NODE * -NTASKS equals the number of MPI tasks."
   "NTHRDS", "Number of OpenMP threads per MPI task."
   "ROOTPE", "The global MPI task of the component root task; if negative, indicates nodes rather than tasks."
   "PSTRID", "The stride of MPI tasks across the global set of pes (for now set to 1)."
   "NINST", "The number of component instances, which are spread evenly across NTASKS."

**Example 1**
~~~~~~~~~~~~~~~

If a component has **NTASKS=16**, **NTHRDS=4** and **ROOTPE=32**, it will run on 64 hardware processors using 16 MPI tasks and 4 threads per task starting at global MPI task 32.

Each CIME component has corresponding entries for ``NTASKS``, ``NTHRDS``, ``ROOTPE`` and ``NINST`` in the **env_mach_pes.xml** file.

**Note:**

- ``NTASKS`` must be greater or equal to 1 even for inactive (stub) components.
- ``NTHRDS`` must be greater or equal to 1.
- If ``NTHRDS`` = 1, this generally means threading parallelization will be off for that component.
- ``NTHRDS`` should never be set to zero.
- The total number of hardware processors allocated to a component is ``NTASKS`` * ``NTHRDS``.
- The coupler processor inputs specify the pes used by coupler computation such as mapping, merging, diagnostics, and flux calculation. This is distinct from the driver, which automatically runs on the union of all processors to manage model concurrency and sequencing.
- The root processor is set relative to the MPI global communicator, not the hardware processors counts. An example of this is below.
- The layout of components on processors has no impact on the science.
- If all components have identical ``NTASKS``, ``NTHRDS``, and ``ROOTPE`` settings, all components will run sequentially on the same hardware processors.

The scientific sequencing is hardwired into the driver. Changing processor layouts does not change intrinsic coupling lags or coupling sequencing.

For a **fully active configuration**, the atmosphere component is hardwired in the driver to never run concurrently with the land or ice component. Performance improvements associated with processor layout concurrency therefore are constrained in this case such that there is never a performance reason not to overlap the atmosphere component with the land and ice components. Beyond that constraint, the land, ice, coupler and ocean models can run concurrently, and the ocean model can also run concurrently with the atmosphere model.

An important but often misunderstood point: The root processor for any given component is set relative to the MPI global communicator, not the hardware processor counts. For instance, in the following example, the atmosphere and ocean will run concurrently, each on 64 processors with the atmosphere running on MPI tasks 0-15 and the ocean running on MPI tasks 16-79.
::

   NTASKS(ATM)=6  NTHRRDS(ATM)=4  ROOTPE(ATM)=0
   NTASKS(OCN)=64 NTHRDS(OCN)=1   ROOTPE(OCN)=16

The first 16 tasks are each threaded 4 ways for the atmosphere. CIME ensures that the batch submission script (**$CASE.run**) automatically requests 128 hardware processors, and the first 16 MPI tasks will be laid out on the first 64 hardware processors with a stride of 4. The next 64 MPI tasks are laid out on the second set of 64 hardware processors.

If you had set ``ROOTPE_OCN`` to 64 in this example, a total of 176 processors would be requested, the atmosphere would be laid out on the first 64 hardware processors in 16x4 fashion, and the ocean model would be laid out on hardware processors 113-176. Hardware processors 65-112 would be allocated but completely idle.

**Example 2**
~~~~~~~~~~~~~~~

If a component has **NTASKS=-2**, **NTHRDS=4** and **ROOTPE=0**, **MAX_MPITASKS_PER_NODE=4**, **MAX_TASKS_PER_NODE=4**, it will run on (8 MPI tasks * 4 threads) = 32 hardware processors on 8 nodes.

If you intended 2 nodes INSTEAD of 8 nodes, then you would change **MAX_MPITASKS_PER_NODE=1** (using **xmlchange**).


**Note**: **env_mach_pes.xml** *cannot* be modified after **case.setup** has been invoked without first running the following:
::

   case.setup --clean





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
