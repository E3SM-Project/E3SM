# The E3SM Atmosphere Model in C++ (EAMxx)

EAMxx is an entirely new atmosphere model for E3SM. It is written in C++ using the Kokkos performance portability library to enable it to run efficiently on CPUs, GPUs, and (hopefully) whatever comes next. Currently only the km-scale explicit-convection version called SCREAM (the Simple Cloud-Resolving E3SM Atmosphere Model) is available, but a low-resolution version is in the works.

Like the documentation for other component models, EAMxx documentation is divided into:

* The [User Guide](user/index.md) which explains how to run EAMxx and describes all the compile-time and run-time options available for modifying a simulation
* The [Developer Guide](developer/index.md) which contains all the information needed to contribute to the development of EAMxx
* The [Technical Guide](technical/index.md) which describes the equations and numerical methods used by each of the specific process implementations and diagnostics available in the current version of EAMxx

Put another way, all information about how to customize runs without changing code is included in the User Guide, general information about software design which is needed for intelligently modifying code goes in the Developer Guide, and details about the specific process implementations in the current model version are included in the Technical Guide. 
