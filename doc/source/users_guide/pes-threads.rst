.. _pesthreads:

==================================
Controlling processors and threads
==================================

Once a compset and resolution for a case has been defined, CIME
provides ways to define the processor layout the case will use.

CIME cases have significant flexibility with respect to the layout of
components across different hardware processors. There are up to eight
unique models (atm, lnd, rof, ocn, ice, glc, wav, cpl) that are
managed independently by the CIME driver, each with a unique MPI
communicator.  In addition, the driver runs on the union of all
processors and controls the sequencing and hardware partitioning.

.. _defining-pes:

pe-settings for a case
-------------------------

CIME looks at the xml element ``PES_SPEC_FILE`` in the **$CIMEROOT/config/$model/config_files.xml** file to determine where
to find the supported out-of-the-box model pe-settings for the primary component (See :ref:`Compsets<compsets>` for definition of primary component.)

When your run `create_newcase  <../Tools_user/create_newcase.html>`_, CIME identifies the primary component and the setting of the ``PES_SPEC_FILE`` in the standard output.

By default, each primary component has a **config_pes.xml** file in
its **cime_config** directory.  That file specifies out-of-the-box
pe-layout for compsets that the primary component defines.  Currently,
the pe-layout can have dependencies on the compset, the model grid and
the target machine.  Finally, there might be more than one
out-of-the-box pe-layout that could be used for a compset/grid/machine
combination: one for a low processor setting and one for a high
processor setting.

A typical entry in a **config_pes.xml** looks like this:

::

  <grid name="a%T62">
    <mach name="cheyenne">
      <pes pesize="any" compset="DATM%IAF">
      .......
      </pes>
    </mach>
  </grid>

Currently, the pesize can have values of ``[any,S,M,L,X1,X2]``.

Given the various dependencies, CIME uses an order of precedence to determine the optimal match. This order is as follows:

1. grid match

   | CIME first searches the grid nodes for a grid match in **config_grids.xml**.
   | The search is based on a regular expression match for the grid longname.
   | All grid matches are then used in the subsequent search.
   | If there is no grid match, all nodes that have ``<grid name="any">`` are used in the subsequent search.

2. machine match

   | CIME next uses the list of nodes obtained in the grid match to search for the machine name using the ``<mach>`` nodes.
   | If there is no machine match, then all nodes with ``<machine name="any">`` are used in the subsequent search.

3. pesize and compset match

   | CIME next uses the list of nodes obtained in the machine match to search for pesize and compset using the ``<pes>`` nodes.
   | If there is no match, the node with ``<pes pesize="any" compset="any">`` is used.

When `create_newcase  <../Tools_user/create_newcase.html>`_  is called, it outputs the matches that are found in determining the best out-of-the-box pe-layout.

Setting the PE layout
---------------------

Optimizing the throughput and efficiency of a CIME experiment often
involves customizing the processor (PE) layout. (See :ref:`load
balancing <optimizing-processor-layout>`.)  CIME provides significant
flexibility with respect to the layout of components across different
hardware processors.  In general, the CIME components -- atm, lnd,
ocn, and so on -- can run on overlapping or mutually unique
processors.  While each component is associated with a unique MPI
communicator, the CIME driver runs on the union of all processors and
controls the sequencing and hardware partitioning.

The pe-layout settings are controlled by the ``$CASEROOT`` file
**env_mach_pes.xml** file. Variables in this file determine the number
of MPI tasks and OpenMP threads for each component, the number of
instances of each component and the layout of the components across
the hardware processors. The entries in **env_mach_pes.xml** have the
following meanings:

.. list-table:: Entries in **env_mach_pes.xml**
   :widths: 10 40
   :header-rows: 1

   * - XML variable
     - Description
   * - MAX_MPITASKS_PER_NODE
     - The maximum number of MPI tasks per node. This is defined in **config_machines.xml** and therefore given a default setting, but can be user modified.
   * - MAX_TASKS_PER_NODE
     - The total number of (MPI tasks) * (OpenMP threads) allowed on a node. This is defined in **config_machines.xml** and therefore given a default setting, but can be user modified. Some computational platforms use a special software customized for the target hardware called symmetric multi-threading (SMT). This allows for over-subscription of the hardware cores. In cases where this is beneficial to model performance, the variable ``MAX_TASKS_PER_NODE`` will be greater than the hardware cores per node as specified by ``MAX_MPITASKS_PER_NODE``.
   * - NTASKS
     - Total number of MPI tasks. A negative value indicates nodes rather than tasks, where *MAX_MPITASKS_PER_NODE \* -NTASKS* equals the number of MPI tasks.
   * - NTHRDS
     - Number of OpenMP threads per MPI task. ``NTHRDS`` must be greater than or equal to 1. If ``NTHRDS`` = 1, this generally means threading parallelization will be off for the given component.
   * - ROOTPE
     -  The global MPI task of the component root task; if negative, indicates nodes rather than tasks. The root processor for each component is set relative to the MPI global communicator.
   * - PSTRID
     - The stride of MPI tasks across the global set of pes (for now set to 1). This variable is currently not used and is a placeholder for future development.
   * - NINST
     -  The number of component instances, which are spread evenly across NTASKS.
   * - COST_PER_NODE
     -  The numbers of cores/node used for accounting purposes. The user should not normally need to set this - but it is useful for understanding how you will be charged.

Each CIME component has corresponding entries for ``NTASKS``, ``NTHRDS``, ``ROOTPE`` and ``NINST`` in the **env_mach_pes.xml** file. The layout of components on processors has no impact on the science.
If all components have identical ``NTASKS``, ``NTHRDS``, and ``ROOTPE`` settings, all components will exectute sequentially on the same hardware processors.

.. hint:: To view the current settings, use the `pelayout <../Tools_user/pelayout.html>`_ tool

The time sequencing is hardwired into the driver. Changing
processor layouts does not change intrinsic coupling lags or coupling
sequencing.

The coupler component has its own processor set for doing
computations such as mapping, merging, diagnostics, and flux
calculation.  This is distinct from the driver, which always
runs on the union of all processors to manage model concurrency and
sequencing.

For a **fully active configuration**, the atmosphere component is
hardwired in the driver to never run concurrently with the land or ice
component.  Performance improvements associated with processor layout
concurrency therefore are constrained in this case such that there is
never a performance reason not to overlap the atmosphere component
with the land and ice components.  Beyond that constraint, the land,
ice, coupler and ocean models can run concurrently, and the ocean
model can also run concurrently with the atmosphere model.

.. note:: if **env_mach_pes.xml** is modified after `case.setup <../Tools_user/case.setup.html>`_  has been called, then you must run `case.setup --reset <../Tools_user/case.setup.html>`_ and the call `case.build <../Tools_user/case.build.html>`_.  **case.build** will only recompile any source code that depends on values in **env_mach_pes.xml**

Case Resource Allocation
------------------------

Resources for your case will be allocated according to the following logic.

* ``NTASKS`` * ``NTHRDS`` is the total number of hardware processors allocated to a component.

* The total number of cores that are allocated will be based on the product of (1) and (2) below where

  1. ``MAX(ROOTPE(comp) + NTASKS(comp))`` across all components
  2. ``MAX(NTHRDS)`` across all components

In the following example, the atmosphere and ocean will run concurrently. The atmosphere will use 16 MPI tasks each with 4 threads per task for a total of 64 cores. The ocean will use 16 MPI tasks with 1 thread per task. BUT since the atmosphere has 4 threads, the ocean will use 64 total cores. The total number of cores will be 128. The atmosphere will run on MPI tasks 0-15 and the ocean will run on MPI tasks 16-31 in the global MPI communicators.

  ::

     NTASKS_ATM=16 NTHRDS_ATM=4  ROOTPE_ATM=0
     NTASKS_OCN=16 NTHRDS_OCN=1  ROOTPE_OCN=16

CIME ensures that the batch submission script (`case.submit
<../Tools_user/case.submit.html>`_ ) will automatically requests 128
hardware processors, and the first 16 MPI tasks will be laid out on
the first 64 hardware processors with a stride of 4. The next 16 MPI
tasks are laid out on the second set of 64 hardware processors in the
same manner, even though the ocean is not threaded.  If you had set
``ROOTPE_OCN`` to 64 in this example, a total of 312 processors would
be requested, the atmosphere would be laid out on the first 64
hardware processors in 16x4 fashion, and the ocean model would be laid
out on hardware processors 255-311. Hardware processors 64-254 would
be allocated but completely idle.

We strongly encourage you to use the `preview_run
<../Tools_user/preview_run.html>`_ script to review the environment
and job submit commands for your case.

.. _optimizing-processor-layout:

Optimizing processor layout
----------------------------

Load balancing is the practice of specifying a processor layout for a given model configuration
(compset, grid, and so on) to maximize simulation speed while minimizing processor idle time.
For a fixed total number of processors, the goal of this optimization is to achieve maximum throughput.
For a set of processor counts, the purpose is to find several "sweet spots" where
the model is minimally idle, cost is relatively low, and the throughput is relatively high.

As with most models, increasing total processors normally results in both increased throughput
and increased cost.
If models scaled linearly, the cost would remain constant across different processor counts,
but models generally don't scale linearly and the cost increases as processor count increases.

Performing a load-balancing exercise on a proposed case before
undertaking a long production run is recommended practice.  Load
balancing requires you to consider a number of factors, such as which
components are run; their absolute and relative resolution; cost,
scaling and processor count sweet spots for each component; and
internal load imbalance within a component.

It is often best to load balance a system with all significant
run-time I/O turned off because it occurs infrequently, typically just
one timestep per simulated  month. It is best treated as a separate cost as it
can otherwise bias interpretation of the overall balance.  Also, the
use of OpenMP threading in some or all of the components is dependent
on the hardware/OS support as well as whether the system supports
running all MPI and mixed MPI/OpenMP on overlapping processors for
different components.

Finally, decide whether components should run sequentially, concurrently, or in some combination.

Typically, a series of short test runs with the desired production
configuration can establish a reasonable load balance setup for the
production job. The timing output can be used to compare test runs to
help determine the optimal load balance.

Changing the pe layout of the model has NO IMPACT on the scientific
results. The basic order of operations and calling sequence are
hardwired into the driver and do not change with the pe
layout. However, both CESM and E3SM do impose some contraints in the
tempororal evolution of the components.  For example, the prognostic
atmosphere model always run sequentially with the ice and land models
for scientific reasons. As a result, running the atmosphere
concurrently with the ice and land will result in idle processors at
some point in the timestepping sequence.

.. hint:: If you need to load balance a fully coupled case, use the :ref:`Load Balancing Tool<load_balancing_tool>`

**One approach to load balancing**

Carry out a :ref:`PFS test <testing>`. This test is by default a
20-day model run with restarts and history output turned off. This
should help you find the layout that has the best load balance for the
targeted number of processors. This provides a reasonable performance
estimate for the production run for most of the runtime.

Seasonal variation and spin-up costs can change performance over time,
so even after a production run has started, review the timing output
occasionally to see if any layout changes might improve throughput or
decrease cost.

In determining an optimal load balance for a specific configuration,
two pieces of information are useful.

* Which components are most expensive.

* How individual components scale. Do they run faster with all MPI or
  mixed MPI/OpenMP decomposition strategies? What are their optimal
  decompositions at each processor count? If the cost and scaling of
  the components are unknown, several short tests with arbitrary
  component pe counts can help establish component scaling and sweet
  spots.

**Determining an optimal load balance**

* Start with the most expensive component and a fixed optimal processor count and decomposition for that component.

* Vary the concurrency and pe counts of the other components.

* Identify a few potential load balance configurations, then run each a few times to establish run-to-run variability and determine the best layout.

In all cases, review the component run times in the timing output file for both overall throughput and independent component timings. Identify idle processors by considering the component concurrency in conjunction with the component timing.

In general, a few component layout options are most reasonable:

* fully sequential,
* fully sequential except the ocean running concurrently,
* fully concurrent except the atmosphere running sequentially with the ice, rof, and land components.

The concurrency is limited in part by hardwired sequencing in the
driver. The sequencing is set by scientific constraints, although
there may be some addition flexibility with respect to concurrency
when running with mixed active and data models.

**Some general rules for finding optimal configurations**

- Make sure you have set a processor layout where each hardware processor is assigned to at least one component. There is rarely a reason to have completely idle processors.

- Make sure your cheapest components keep up with your most expensive components. In other words, a component that runs on 1024 processors should not be waiting on a component running on 16 processors.

- Before running the job, make sure the batch queue settings are set correctly for your run. Review the account numbers, queue names and time limits. The ideal time limit, queue and run length are dependent on each other and on the current model throughput.

- Take full advantage of the hardware resources. If you are charged by the 32-way node, you might as well target a total processor count that is a multiple of 32.

- Keep a single component on a single node, if possible, to minimize internal component communication cost.

- Assume that hardware performance can vary due to contention on the interconnect, file systems, or other areas. If you are unsure of a timing result, run cases multiple times.

The pe-layout and the associated timings are found in the  :ref:`timing files <model-timing-data>` generated for your run.
