# EAMxx

## A High-performance E3SM Atmosphere Model in C++

EAMxx was designed and built from the ground up employing modern C++ in order to
embrace leading-edge software practices and enable
"performance portability"
[^perf-port_def] across various
leadership-class supercomputers.
The latter goal is achieved by using the Kokkos library.
EAMxx is a "clean-start" model with almost no similarity to the E3SM atmosphere
model used in versions 1-3.
Currently only the km-scale explicit-convection version called
SCREAM[^eamxx_v_scream] (the Simple Cloud-Resolving E3SM Atmosphere Model) is
available, but a low-resolution version is in development.

## Trail Map

Like the documentation for other component models, EAMxx documentation is
divided into the following sections:

- [User Guide](user/index.md) - info about running EAMxx and all options
for modifying a simulation
- [Developer Guide](developer/dev_quickstart.md) - information needed to contribute
to EAMxx development
- [Technical Guide](technical/index.md) - equations and numerical methods used
in EAMxx

Put another way, the ***User Guide*** provides information about how to
successfully run a model that includes EAMxx and how to customize runs without
changing source code.
The ***Developer Guide*** includes information about EAMxx software design
that will assist those who wish to intelligently modify source code.
Scientific and modeling details about the specific process implementations
in the current model version are provided in the ***Technical Guide***.

<!-- !!! info "Super cool eamxx figure goes here" -->

<!-- ======================================================================= -->

[^perf-port_def]: ***Performance Portability:*** A software-design practice focused on providing the ability to run the same model on differing architectures (cpu, gpu, ...) without modifying the source code and, ideally, without sacrificing performance.
[^eamxx_v_scream]: ***Note:*** Prior to **EAMxx** version 1, this project, software library, and atmosphere driver were generally referred to as **SCREAM**. However, at this point, we wish to make the distinction between **EAMxx**, the next-gen atmosphere driver and C++ software library, and **SCREAM**, the particular high-resolution configuration of **EAMxx**.
