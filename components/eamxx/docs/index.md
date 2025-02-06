# EAMxx

## A High-performance E3SM Atmosphere Model in C++

<!-- EAMxx is almost completely different in all ways from the atmosphere model used for E3SM versions 1-3. -->

EAMxx was designed and built from the ground up employing modern C++ in order to embrace leading-edge software practices and enable "[performance portability](common/glossary.md#performance-portability)"[^perf-port_def] across various leadership-class supercomputers.
The latter goal is achieved by using the Kokkos library.
EAMxx is a "clean-start" model with almost no similarity to the E3SM atmosphere model used in versions 1-3.
Currently only the km-scale explicit-convection version called SCREAM[^eamxx_v_scream] (the Simple Cloud-Resolving E3SM Atmosphere Model) is available, but a low-resolution version is in development.

## Trail Map

Like the documentation for other component models, EAMxx documentation is divided into:

* The [User Guide](user/index.md) - info about running EAMxx and all options for modifying a simulation
* The [Developer Guide](developer/index.md) - information needed to contribute to EAMxx development
* The [Technical Guide](technical/index.md) - equations and numerical methods used in EAMxx
* [Testing](testing/index.md) - a detailed description of the testing philosophy and practices tied to EAMxx development and maintenance
    * A focused treatment of testing, specifically tailored to **Users** or **Developers** appears in the respective guides.
    * ==Move this back to Dev Guide==

### Details

Put another way, all information about how to customize runs without changing code is included in the User Guide, general information about software design required to intelligently modify code goes in the Developer Guide, and details about the specific process implementations in the current model version are included in the Technical Guide.
Finally, we devote a section to detailing our rigorous verification and validation testing because this was integral to our approach of moving quickly and embracing bleeding-edge software practices to design a novel, high-resolution model of Earth's atmosphere.

!!! info "Super cool eamxx figure goes here"


[^perf-port_def]: ***Performance Portability:*** A software-design practice focused on providing the ability to run the same model on differing architectures (cpu, gpu, ...) without modifying the source code and, ideally, without sacrificing performance.
[^eamxx_v_scream]: ***Note:*** Before **EAMxx** version 1, the project, software library, and atmosphere driver were generally referred to as **SCREAM**. However, at this point, we wish to make the distinction between **EAMxx**, the next-gen atmosphere driver and C++ software library, and **SCREAM**, the particular high-resolution configuration of **EAMxx**.
