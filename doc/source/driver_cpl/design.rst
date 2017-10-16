Design
======

--------
Overview
--------
cpl7 is built as a single executable with a single high-level driver.
The driver runs on all processors and handles coupler sequencing, model concurrency, and communication of data between components. 
The driver calls all model components via common and standard interfaces.
The driver also directly calls coupler methods for mapping (interpolation), rearranging, merging, an atmosphere/ocean flux calculation, and diagnostics. 
The model components and the coupler methods can run on subsets of all the processors.
In other words, cpl7 consists of a driver that controls the top level sequencing, the processor decomposition, and communication between components and the coupler while coupler operations such as mapping and merging are running under the driver on a subset of processors as if there were a unique coupler model component.

In general, an active component both needs data from and provides data to the coupler while data models generally read data from I/O and then just provide data to the coupler.
Currently, the atmosphere, land, river, and sea ice models are always tightly coupled to better resolve the diurnal cycle. 
This coupling is typically half-hourly, although at higher resolutions, can be more frequent.
The ocean model coupling is typically once or a few times per day. 
The diurnal cycle of ocean surface albedo is computed in the coupler for use by the atmosphere model.
The looser ocean coupling frequency means the ocean forcing and response is lagged in the system. 
There is an option in cpl7 to run the ocean tightly coupled without any lags, but this is more often used only when running with data ocean components.

--------------------------
Sequencing and Concurrency
--------------------------
The component processor layouts and MPI communicators are derived from namelist input. 
At the present time, there are eight (10) basic processor groups in cpl7.
These are associated with the atmosphere, land, river, ocean, sea ice, land ice, wave, external-system-process, coupler, and global groups, although others could be easily added later. 
Each of the processor groups can be distinct, but that is not a requirement of the system.
A user can overlap processor groups relatively arbitrarily. 
If all processors sets overlap each other in at least one processor, then the model runs sequentially.
If all processor sets are distinct, the model runs as concurrently as science allows. 
The processor sets for each component group are described via 3 basic scalar parameters at the present time; the number of mpi tasks, the number of openmp threads per mpi task, and the global mpi task rank of the root mpi task for that group.
For example, a layout where the number of mpi tasks is 8, the number of threads per mpi task is 4, and the root mpi task is 16 would create a processor group that consisted of 32 hardware processors, starting on global mpi task number 16 and it would contain 8 mpi tasks. 
The global group would have at least 24 tasks and at least 48 hardware processors.
The driver derives all MPI communicators at initialization and passes them to the component models for use. 
More information on the coupler concurrency can be found in the Craig et al IJHPCA 2012 reference mentioned in the top section of this document.

As mentioned above, there are two issues related to whether the component models run concurrently. 
The first is whether unique chunks of work are running on distinct processor sets.
The second is the sequencing of this work in the driver. 
As much as possible, the driver sequencing has been implemented to maximize the potential amount of concurrency of work between different components.
Ideally, in a single coupling step, the forcing for all models would be computed first, the models could then all run concurrently, and then the driver would advance. 
However, scientific requirements such as the coordination of surface albedo and atmosphere radiation computations as well as general computational stability issues prevents this ideal implementation in cpl7.
`Figure 1 <cplug-02.1-figx1.jpg?raw=true>`_ shows the maximum amount of concurrency supported for a fully active system. 
In practice, the scientific constraints mean the active atmosphere model cannot run concurrently with the land, runoff, and sea-ice models.
Again, `figure 1 <cplug-02.1-figx1.jpg?raw=true>`_ does not necessarily represent the optimum processor layout for performance for any configuration, but it provides a practical limit to the amount of concurrency in the system due to scientific constraints. 
Results are bit-for-bit identical regardless of the component sequencing because the scientific lags are fixed by the implementation, not the processor layout.

image:: cplug-02.1-figx1.jpg

Figure 1: Maximum potential processor concurrency designed to support scientific requirements and stability.

--------------------
Component Interfaces
--------------------
The standard cpl7 component model interfaces are based upon the ESMF design.
Each component provides an init, run, and finalize method with consistent arguments. 
The component interface arguments currently consist of Fortran and MCT datatypes.
The physical coupling fields are passed through the interfaces in the init, run, and finalize phases. 
As part of initialization, an MPI communicator is passed from the driver to the component, and grid and decomposition information is passed from the component back to the driver.
The driver/coupler acquires all information about resolution, configurations, and processor layout at run-time from either namelist or from communication with components.


Initialization of the system is relatively straight-forward.
First, the MPI communicators are computed in the driver. 
Then the component model initialization methods are called on the appropriate processor sets, and an mpi communicator is sent, and the grid and decomposition information are passed back to the driver.
Once the driver has all the grid and decomposition information from the components, various rearrangers and mappers are initialized that will move data between processors, decompositions, and grids as needed at the driver level. 
No distinction is made in the coupler implementation for sequential versus concurrent execution.
In general, even for cases where two components have identical grids and processor layouts, often their decomposition is different for performance reasons. 
In cases where the grid, decomposition, and processor layout are identical between components, the mapping or rearranging operation will degenerate to a local data copy.

The interface to the components' run method consists of two distinct bundles of fields. 
One is the data sent to force the model.
The second is data received from the model for coupling to other components. 
The run interface also contains a clock that specifies the current time and the run length for the model and a data type that encapsulates grid, decomposition, and scalar coupling information.
These interfaces generally follow the ESMF design principles.

-------------------------------
MCT, The Model Coupling Toolkit
-------------------------------
In cpl7, the MCT attribute_vector, global_segmap, and general_grid datatypes have been adopted at the highest levels of the driver, and they are used directly in the component init, run, and finalize interfaces. 
In addition, MCT is used for all data rearranging and mapping (interpolation).
The clock used by cpl7 at the driver level is based on the ESMF specification. 
Mapping weights are still generated off-line using the SCRIP or ESMF packages as a preprocessing step.
They are read using a subroutine that reads and distributes the mapping weights in reasonably small chunks to minimize the memory footprint. 
Development of the cpl7 coupler not only relies on MCT, but MCT developers contributed significantly to the design and implementation of the cpl7 driver.
Development of cpl7 coupler resulted from a particularly strong and close collaboration between NCAR and the Department of Energy Argonne National Lab.

------------------------------------
Memory, Parallel IO, and Performance
------------------------------------
Scaling to tens-of-thousands of processors requires reasonable performance scaling of the models, and all components have worked at improving scaling via changes to algorithms, infrastructure, or decompositions.
In particular, decompositions using shared memory blocking, space filling curves, and all three spatial dimensions have been implemented to varying degrees in all components to increase parallelization and improve scalability. 
The Craig et al IJHPCA 2012 reference mentioned in the first section of this document provides a summary of scaling performance of cpl7 for several coupler kernals.

In practice, performance, load balance, and scalability are limited as a result of the size, complexity, and multiple model character of the system. 
Within the system, each component has its own scaling characteristics.
In particular, each may have processor count "sweet-spots" where the individual component model performs particularly well. 
This might occur within a component because of internal load balance, decomposition capabilities, communication patterns, or cache usage.
Second, component performance can vary over the length of the model run. 
This occurs because of seasonal variability of the cost of physics in models, changes in performance during an adjustment (spin-up) phase, and temporal variability in calling certain model operations like radiation, dynamics, or I/O.
Third, the hardware or batch queueing system might have some constraints on the total number of processors that are available. 
For instance, on 16 or 32 way shared memory node, a user is typically charged based on node usage, not processor usage.
So there is no cost savings running on 40 processors versus 64 processors on a 32-way node system. 
As a result of all of these issues, perfect load-balancing is generally not possible.
But to a large degree, if one accepts the limitations, a load balance configuration with acceptable idle-time and reasonably good throughput is nearly always possible to configure.


Load-balancing requires a number of considerations such as which components are run, their absolute resolution, and their relative resolution; cost, scaling and processor count sweet-spots for each component; and internal load imbalance within a component.
It is often best to load balance the system with all significant run-time I/O turned off because this generally occurs very infrequently (typically one timestep per month), is best treated as a separate cost, and can bias interpretation of the overall model load balance. 
Also, the use of OpenMP threading in some or all of the system is dependent on the hardware/OS support as well as whether the system supports running all MPI and mixed MPI/OpenMP on overlapping processors for different components.
Finally, should the components run sequentially, concurrently, or some combination of the two. 
Typically, a series of short test runs is done with the desired production configuration to establish a reasonable load balance setup for the production job.



