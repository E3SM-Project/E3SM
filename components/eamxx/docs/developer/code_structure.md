# EAMxx Source Code Structure

!!! Info inline "E3SM/components/eamxx/"

    ```
    .
    ├── cime_config/
    ├── cmake/
    │   └── machine-files/
    │── CMakeLists.txt
    │── CTestConfig.cmake
    ├── docs/
    ├── mkdocs.yml
    ├── README.md
    ├── scripts/
    │   ├── atmchange
    │   ├── atmquery
    │   └── test-all-eamxx
    ├── src/
    │   ├── control/
    │   ├── diagnostics/
    │   ├── dynamics/
    │   │   └── homme/
    │   ├── mct_coupling/
    │   ├── physics/
    │   │   ├── cld_fraction/
    │   │   ├── cosp/
    │   │   ├── iop_forcing/
    │   │   ├── mam/
    │   │   ├── ml_correction/
    │   │   ├── nudging/
    │   │   ├── p3/
    │   │   ├── rrtmgp/
    │   │   ├── share/
    │   │   ├── shoc/
    │   │   ├── spa/
    │   │   └── tms/
    │   ├── python/
    │   ├── eamxx_config.h.in
    │   └── share/
    ├── tests/
    │   ├── generic/
    │   ├── meta-tests/
    │   ├── multi-process/
    │   │   ├── dynamics_physics/
    │   │   └── physics_only/
    │   ├── python/
    │   └── single-process/
    └── tpls/
    ```

<!-- <div class="grid" markdown> -->
=== "EAMxx Root Directory"
    All EAMxx-specific code can be found in `components/eamxx` within the
    [E3SM repo](https://github.com/E3SM-Project/E3SM).
    Here's how things are organized:

    - `cime_config`: Tools and XML files for integrating EAMxx with E3SM via the
      CIME framework.
    - `cmake`: CMake functions and macros used by the configuration/build system.
    - `docs`: Documentation for the EAMxx project, including design documents,
      instructions for building and testing EAMxx, and this document.
    - `scripts`: Miscellaneous scripts that implement workflows for running tests
      and analyzing performance.
    - `src`: All C++ source code (and any bridges to Fortran) for EAMxx are stored
      here. We describe the contents of this directory in greater detail below.
    - `tests`: Implements standalone, end-to-end tests for various EAMxx
      components (RRTMG, HOMME, P3, SHOC, etc).
    - `tpls`: Utilities for building EAMxx-specific third-party libraries.[^arcane]
    [^arcane]:
        *Here there be dragons...* These are some of the more arcane EAMxx
        features that even the bravest of developers may never lay eyes on.
    
    In addition, you'll notice the following files in `components/eamxx`:
    
    - `CMakeLists.txt`: The CMake file that defines EAMxx's configuration/build
      system.
    - `CTestConfig.cmake`: This CTest file contains parameters that determine how
      our test results are reported to the [E3SM CDash Site](http://my.cdash.org/submit.php?project=E3SM).
    - `README.md`: EAMxx's top-level README file, which describes the project and
      its purpose.
    - `mkdocs.yml`: The configuration file for [mkdocs](https://www.mkdocs.org/),
      the tool we currently use to build and publish our documentation.

=== "`src/` Directory"

    Herein lies the source code for EAMxx.
    Broadly, here's where things are located:
    
    - `control`: Contains the atmosphere driver and basic tests for it.
    - `diagnostics`: A collection of simple classes used to compute diagnostic
      quantities.
    - `dynamics`: Home the files implementing the necessary interfaces between
      EAMxx and the HOMME non-hydrostatic dycore.
        - ***Note:*** HOMME itself is in components/homme
    - `mct_coupling`: Glue code for embedding EAMxx within E3SM as an atmosphere
      component using the MCT coupler.
    - `physics`: Source code for physics-related atmospheric processes, including
        - `cosp`: Diagnostic-only[^diag] package used for computing
        radiation-related quantities.
        [^diag]:
            ***Diagnostic:*** Refers to quantities or variables that are not
            inputs to any part of the model but capture useful descriptive data.
        - `iop_forcing`: "Intensive Observation Period" package that can be used
        in conjunction with the Doubly-Periodic (DP) configuration of SCREAM.
        - `mam`: Contains the high-performance **M**odal **A**erosol **M**odel,
        parameterized with **4** size modes known as **MAM4xx**.
        - `ml_correction`: This is a product developed as a part of a Lawrence
        Livermore National Laboratory LDRD project in which the goal was to
        improve the results of low-resolution SCREAM runs by employing
        machine-learning algorithms to analyze the output of nudged
        high-resolution runs.
            - This code is not actively supported but may be useful/or
            interesting to some.
        - `nudging`: Contains machinery that can be applied to "nudge"
        quantities in EAMxx toward desired values, typically coming from
        previous model runs or reanalysis data.
            - See the [Nudging in EAMxx](../user/nudging.md) page in the
            User Guide for further details.
        - `p3`: The C++/Kokkos implementation of P3 microphysics within EAMxx.
        - `rrtmgp`: A stub for the radiation processes as represented in EAMxx.
        - `share`: Utilities and data structures common to these processes.
        - `shoc`: The C++/Kokkos implementation of SHOC macrophysics within EAMxx.
    - `python`: Source code for the experimental ***PyEAMxx*** package that
    enables building and running an individual EAMxx atmosphere process using python/[conda](https://docs.conda.io/en/latest/).
        - As of time of writing, this feature is still in development and should
        e considered a prototype.
            - See the [pyEAMxx](../user/pyeamxx.md) page in the User Guide for
            a more detailed description.
    - `share`: Utilities used by various components within EAMxx. Of note:
        - `io`: EAMxx's interface to the [SCORPIO](https://e3sm.org/scorpio-parallel-io-library/)
        library.

    Each of these directories contains a `CMakeLists.txt` file for defining how
    EAMxx is configured and built, and a `tests/` subdirectory that houses relevant
    unit and verification tests.
    
    You'll also see some other files in the `src/` directory itself, such as
    
    - `eamxx_config.h.in`: A template for generating a C++ header file with
      EAMxx configuration information.

=== "`tests/` Directory"

    - `generic/`:
    - `meta-tests/`: 
    - `multi-process/`: 
        - `dynamics_physics/`: 
        - `physics_only/`: 
    - `python/`: 
    - `single-process/`: 
