# The E3SM Atmosphere Model in C++ (EAMxx)

##EAMxx
EAMxx is almost completely different in all ways from the atmosphere model used for E3SM versions 1-3. 


EAMxx was built from the ground up using C++ in order to embrace modern software practices and to allow "performance portability" across various supercomputers. The latter goal is achieved by using the Kokkos library. EAMxx is a "clean-start" model with almost no similarity to the E3SM atmosphere model used in versions 1-3.  Currently only the km-scale explicit-convection version called SCREAM (the Simple Cloud-Resolving E3SM Atmosphere Model) is available, but a low-resolution version is in the works.

Like the documentation for other component models, EAMxx documentation is divided into:

* The [User Guide](user/index.md) - info about running EAMxx and all options for modifying a simulation
* The [Developer Guide](developer/index.md) - information needed to contribute to EAMxx development
* The [Technical Guide](technical/index.md) - equations and numerical methods used in EAMxx

Put another way, all information about how to customize runs without changing code is included in the User Guide, general information about software design which is needed for intelligently modifying code goes in the Developer Guide, and details about the specific process implementations in the current model version are included in the Technical Guide. 
