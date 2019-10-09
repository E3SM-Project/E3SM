====================================
Initialization and Restart
====================================

The initialization has been developed over the last two decades to meet the scientific goals, minimize the communication required, and ensure a consistent and well defined climate system.
The order of operations is critical. The initialization is basically as follows:

- The ``ccsm_pes`` namelist is read and mpi communicators are initialized.
- The ``seq_infodata`` namelist is read and configuration settings are established.
- The ``prof_inparm`` namelist is read and the timing tool is initialized.
- The ``pio_inparm`` namelist is read and the driver IO is initialized.
- The ``seq_timemgr`` namelist is read and the driver time manager and clocks are initialized.
- The atmosphere init routine is called, the mpi communicator and clock are sent,  and the atmosphere grid is returned.
- The land init routine is called, the mpi communicator and clock are sent, and the land grid is returned.
- The runoff init routine is called, the mpi communicator and clock are sent, and the runoff grid is returned.
- The ocean init routine is called, the mpi communicator and clock are sent, and the ocean grid is returned.
- The ice init routine is called, the mpi communicator and clock are sent,  and the ice grid is returned.
- The land ice init routine is called, the mpi communicator and clock are sent, and the land ice grid is returned.
- The infodata buffer is synchronized across all processors.  This buffer contains many model configuration settings set by the driver but also  sent from the components.
- The atmosphere, land, runoff, ice, land ice, and ocean rearrangers are initialized.  - These rearrangers move component data between the component pes and the coupler pes.
- The Remaining attribute datatypes associated are initialized
- The mapping weights and areas are read.  
- The Component grids are checked using the domain checking method.
- The flux area corrections are initialized on the component pes and applied to the initial fields sent by each component on the component pes. Those initial fields are then rearranged to the coupler pes.
- The fractions are initialized on the coupler pes. 
- The atmosphere/ocean flux computation is initialized and initial ocean albedos are computed on the coupler pes.
- The land, ocean, and ice initial albedos are mapped to the atmosphere grid and merged to generate initial surface albedos.
- The initial atmosphere forcing data (albedos) is rearranged from the coupler pes to the atmosphere pes, and the area corrections are applied.
- The second phase of the atmosphere init method is to initialize the atmosphere radiation from the surface albedos.
- The new atmosphere initial data is area corrected and rearranged to the coupler pes.
- The budget diagnostics are zeroed out.
- The coupler restart file is read.
- Initialization is complete.

