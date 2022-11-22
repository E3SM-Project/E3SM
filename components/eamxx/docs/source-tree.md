# SCREAM's Source Tree

All SCREAM-specific code can be found in `components/eamxx` within the
[SCREAM repo](https://github.com/E3SM-Project/scream). Here's how things are
organized:

+ `cime_config`: Tools and XML files for integrating SCREAM with E3SM via the
  CIME framework.
+ `cmake`: CMake functions and macros used by the configuration/build system.
+ `data`: Data files used by our tests.
+ `docs`: Documentation for the SCREAM project, including design documents,
  instructions for building and testing SCREAM, and this document.
+ `extern`: Source for certain lightweight third-party libraries, embedded
  directly into the repo (and not imported as submodules).
+ `scripts`: Miscellaneous scripts that implement workflows for running tests
  and analyzing performance.
+ `src`: All C++ source code (and any bridges to Fortran) for SCREAM are stored
  here. We describe the contents of this directory in greater detail below.
+ `tests`: Implements standalone, end-to-end tests for various SCREAM
  components (RRTMG, HOMME, P3, SHOC, etc).

In addition, you'll notice the following files in `components/eamxx`:

+ `CMakeLists.txt`: The CMake file that defines SCREAM's configuration/build
  system.
+ `CTestConfig.cmake`: This CTest file contains parameters that determine how
  our test results are reported to the [E3SM CDash Site](http://my.cdash.org/submit.php?project=E3SM).
+ `README.md`: SCREAM's top-level README file, which describes the project and
  its purpose.

## The `src` Directory

Herein lіes the source code for SCREAM. Broadly, here's where things are:

+ `control`: Contains the atmosphere driver and basic tests for it.
+ `dynamics`: Here's where HOMME lives within SCREAM, along with code for
  interfacing with it using SCREAM's data structures.
+ `interface`: Glue code for embedding SCREAM within E3SM as an atmosphere
  component.
+ `physics`: Source code for physics-related atmospheric processes, including
  + `p3`: The C++/Kokkos implementation of P3 microphysics within SCREAM.
  + `shoc`: The C++/Kokkos implementation of SHOC macrophysics within SCREAM.
  + `rrtmgp`: A stub for the radiation processes as represented in SCREAM.
  + `common`: Utilities and data structures common to these processes.
+ `share`: Utilities used by various components within SCREAM. A lot of things
  here will likely end up in `ekat`.

Each of these directories contains a `CMakeLists.txt` file for defining how
things are build, and a `tests/` subdirectory that houses relevant
unit and verification tests.

You'll also see some other files in the `src/` directory itself:

+ `scream_config.f.in`: A template for generating a Fortran include file with
  SCREAM configuration information.
+ `scream_config.h.in`: A template for generating a C++ header file with
  SCREAM configuration information.

