# EAMxx's Source Tree

All EAMxx-specific code can be found in `components/eamxx` within the
[EAMxx repo](https://github.com/E3SM-Project/scream). Here's how things are
organized:

+ `cime_config`: Tools and XML files for integrating EAMxx with E3SM via the
  CIME framework.
+ `cmake`: CMake functions and macros used by the configuration/build system.
+ `data`: Data files used by our tests.
+ `docs`: Documentation for the EAMxx project, including design documents,
  instructions for building and testing EAMxx, and this document.
+ `scripts`: Miscellaneous scripts that implement workflows for running tests
  and analyzing performance.
+ `src`: All C++ source code (and any bridges to Fortran) for EAMxx are stored
  here. We describe the contents of this directory in greater detail below.
+ `tests`: Implements standalone, end-to-end tests for various EAMxx
  components (RRTMG, HOMME, P3, SHOC, etc).

In addition, you'll notice the following files in `components/eamxx`:

+ `CMakeLists.txt`: The CMake file that defines EAMxx's configuration/build
  system.
+ `CTestConfig.cmake`: This CTest file contains parameters that determine how
  our test results are reported to the [E3SM CDash Site](http://my.cdash.org/submit.php?project=E3SM).
+ `README.md`: EAMxx's top-level README file, which describes the project and
  its purpose.
+ `mkdocs.yml`: The configuration file for [mkdocs](https://www.mkdocs.org/),
  the tool we currently use to build and publish our documentation.

## The `src` Directory

Herein lіes the source code for EAMxx. Broadly, here's where things are:

+ `control`: Contains the atmosphere driver and basic tests for it.
+ `dynamics`: Here's where HOMME lives within EAMxx, along with code for
  interfacing with it using EAMxx's data structures.
+ `mct_coupling`: Glue code for embedding EAMxx within E3SM as an atmosphere
  component using the MCT coupler.
+ `physics`: Source code for physics-related atmospheric processes, including
  + `p3`: The C++/Kokkos implementation of P3 microphysics within EAMxx.
  + `shoc`: The C++/Kokkos implementation of SHOC macrophysics within EAMxx.
  + `rrtmgp`: A stub for the radiation processes as represented in EAMxx.
  + `share`: Utilities and data structures common to these processes.
+ `share`: Utilities used by various components within EAMxx. Of note:
  + `io`: EAMxx's interface to the [SCORPIO](https://e3sm.org/scorpio-parallel-io-library/)
    library.
+ `diagnostics`: A collection of simple classes used to compute diagnostic
  quantities.

Each of these directories contains a `CMakeLists.txt` file for defining how
things are build, and a `tests/` subdirectory that houses relevant
unit and verification tests.

You'll also see some other files in the `src/` directory itself, such as

+ `scream_config.h.in`: A template for generating a C++ header file with
  EAMxx configuration information.
