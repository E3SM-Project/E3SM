# Quick-start Guide for Developers

- Mention full-model, but don't go into detail
- Focus on standalone build (test all scream)
  - maybe high-level description, then link to it
  - why you should use it
  - presets for important configurations
    - dbg, dbg-single, etc.
  - also possible to roll your own cmake build and focus narrowly on what you're doing
  - reach out if you believe you're running into a bug within `test-all-scream`
  - mention baselines/-comparisons
  

==Move Source Tree Here==

## Motivation

The intention of this guide is to provide a new developer with the most straightforward path toward productively building, running, and testing EAMxx.
Specifically, we will discuss building and testing in the EAMxx "standalone" configuration using the `test-all-scream` workflow.
Note that tests related to running the full E3SM model are discussed in [Full Model Testing](dev_testing/full_model_testing.md).

## Stratosphere-level Overview

At the most basic level, EAMxx is built and run using [CMake](https://cmake.org), a widely-used tool to "automatically"[^automagic-cmake] configure, build, and test a project that will run on the system being employed.
Thus, one can certainly build and test EAMxx with their own customized CMake setup.
However, we would strongly advise against doing this unless you know what you're doing and have good reasons for choosing this option.
Instead, we recommend using the EAMxx-designed `test-all-scream` workflow that employs a variety of shell and python scripts to automate the series of tasks common to development.
More details on this are provided below in the [Testing Section](#code-build-test-repeat), and full details may be found in [EAMxx Automated Standalone Testing](dev_testing/test_all_scream.md)

### Note on Computing Platforms

EAMxx is tested and ==(unofficially?) supported== on the [E3SM-supported machines](https://e3sm.org/model/running-e3sm/supported-machines/).
This means that, if you have an account on one of these machines, EAMxx tests can be run with no extra configuration required.
However, many of the standalone EAMxx tests are able to be run on a reasonably-equipped laptop or personal workstation.[^but-mac]
This second choice requires some extra configurations but the effort required may be worth the potential speedup in development.

## Let's Get Started (Quickly)

For the purposes of this quick-start guide, we will employ the working example of adding a new [processes](../common/glossary.md#atmospheric-process) to the MAM4xx aerosol library, specifically the not-real, example process "**Aerosolve**."
The source files for this process will be found at `src/physics/mam/`.
Note that here, and for the remainder of this guide, we will denote all paths relative to the EAMxx root directory `E3SM/components/eamxx/` (see [Code Structure](code_structure.md) for additional details on how EAMxx code is organized).

### Modifying Source Code

The minimum files required for a new process will be those defining the interface to EAMxx.
For our case, and following the current informal naming standard, these will be the `eamxx_aerosolve_process_interface.<x>pp` files, and trivial skeleton examples are provided below.

<div class="grid" markdown>
```c++ title="src/physics/mam/eamxx_aerosolve_process_interface.hpp"
#include <src/physics/mam/all_the_headers.hpp>

namespace scream {
class MAM_Aerosolve final : public scream::AtmosphereProcess {
  using KT = ekat::KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Real>;
  using view_2d = typename KT::template view_2d<Real>;

  view_1d answer("answer", nlevels, nthings);
  view_2d input("input", nthings);

private:
  void aerosolver(const view_2d input, const int indicator, view_1d &sol);

public:
  // Initialize variables
  void initialize_impl(const view_2d inmat) override;
  // Run the process by one time step
  void run_impl(const double dt) override;
  // Finalize--cleanup anything that needs to be manually released
  void finalize_impl(){};
}
}
```

```c++ title="src/physics/mam/eamxx_aerosolve_process_interface.cpp"
#include <src/physics/mam/eamxx_aerosolve_process_interface.hpp>

namespace scream {
void MAM_Aerosolve::initialize_impl(const view_2d inmat) {
  [...]
}
void MAM_Aerosolve::run_impl(const double dt) {
  [...]
  aerosolver(input, i, answer);
}
void MAM_Aerosolve::aerosolver(const view_2d input, const int indicator,
                               view_1d &ans) {
  if (indicator == 42) {
    ans = Kokkos::subview(input, indicator, Kokkos::ALL());
  } else {
    Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, nthings), [&](const int i) {
          ans[i] = input(0, i);
    });
  }
}
```
</div>

#### Add New Process to CMake Build

In order to ensure the code you've added is configured and built by CMake, we need to add to the relevant `CMakeLists.txt` file.
At minimum, we need to add `eamxx_aerosolve_process_interface.cpp` to the EAMxx library named `mam`, as shown below.
Of course, there may be additional configuration commands that are required for a proper build, but this is beyond the scope of this guide.

```cmake title="src/physics/mam/CMakeLists.txt"
[...]
# EAMxx mam4xx-based atmospheric processes
add_library(mam
  eamxx_mam_microphysics_process_interface.cpp
  eamxx_mam_optics_process_interface.cpp
  [...]
  eamxx_aerosolve_process_interface.cpp)
[...]
```

### Adding a Test to Verify/Validate New Process

All new code that is added to EAMxx is required to be tested.
At minimum, a process must be **Validation Tested**.
That is, tested against known, trusted data, and this data can be samples taken by real-world instruments or comparable results generated by other code.
Ideally, any new and sufficiently complex functions should also have **Verification Tests** that test for mathematical or scientific "correctness".[^correct-huh]

#### Single-process Validation Test

First, we add a validation test for (only) our new process to `tests/single-process/mam/aerosolve/`.
Defining the test behavior and pass/fail criteria requires 3 files:
- `CMakeLists.txt`
    - Largely reusable boilerplate code but should be tailored to the proper test parameters.
    - Defines build and run behavior.
    - Enumerates the required NetCDF input files.
- `input.yaml`
    - Defines key parameters for the model and test run.
    - Enumerates the keys (`std::string`) used to access the input files from within the process source code.
    - Provides information about the test's initial condition.
        - Input file names/paths.
        - Other required input values.
- `output.yaml`
    - Defines the output fields used as pass/fail criteria.
    - Defines the tests runtime output behavior.

```cmake title="tests/single-process/mam/aerosolve/CMakeLists.txt"
include (ScreamUtils)

set (TEST_BASE_NAME mam4_aerosolve_standalone)
set (FIXTURES_BASE_NAME ${TEST_BASE_NAME}_generate_output_nc_files)
# Create the test
CreateADUnitTest(${TEST_BASE_NAME}
  LABELS mam4_aerosolve physics
  LIBS mam
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  FIXTURES_SETUP_INDIVIDUAL ${FIXTURES_BASE_NAME}
)

# Set AD configurable options
set (ATM_TIME_STEP 1800)
SetVarDependingOnTestSize(NUM_STEPS 2 5 48)  # 1h 2.5h 24h
set (RUN_T0 2021-10-12-45000)

# Copy (and configure) yaml files needed by tests
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/output.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/output.yaml)

# Ensure test input files are present in the data dir
set (TEST_INPUT_FILES
  scream/init/${aerosolve_IC_file}
  scream/mam4xx/.../aerosolve_input_data_file_ABC.nc
  scream/mam4xx/.../aerosolve_input_data_file_XYZ.nc
  [...]
  <path-to-required-input-file>
)
foreach (file IN ITEMS ${TEST_INPUT_FILES})
  GetInputFile(${file})
endforeach()

# Compare results between runs employing differing numbers of MPI ranks
# to ensure they are bfb
include (CompareNCFiles)
CompareNCFilesFamilyMpi (
  TEST_BASE_NAME ${TEST_BASE_NAME}
  FILE_META_NAME ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.npMPIRANKS.${RUN_T0}.nc
  MPI_RANKS ${TEST_RANK_START} ${TEST_RANK_END}
  LABELS mam4_aerosolve physics
  META_FIXTURES_REQUIRED ${FIXTURES_BASE_NAME}_npMPIRANKS_omp1
)
# Compare one of the output files with the baselines.
# Note: one is enough, since we already check that np1 is BFB with npX
if (SCREAM_ENABLE_BASELINE_TESTS)
  set (OUT_FILE ${TEST_BASE_NAME}_output.INSTANT.nsteps_x1.np${TEST_RANK_END}.${RUN_T0}.nc)
  CreateBaselineTest(${TEST_BASE_NAME} ${TEST_RANK_END} ${OUT_FILE} ${FIXTURES_BASE_NAME}) 
endif()
```

```yaml title="tests/single-process/mam/aerosolve/input.yaml"
%YAML 1.1
---
driver_options:
  atmosphere_dag_verbosity_level: 5
time_stepping:
  time_step: ${ATM_TIME_STEP}
  run_t0: ${RUN_T0}  # YYYY-MM-DD-XXXXX
  number_of_steps: ${NUM_STEPS}
atmosphere_processes:
  atm_procs_list: [mam4_aerosolve]
  mam4_aerosolve:
    aerosolve_ABC_data: ${SCREAM_DATA_DIR}/mam4xx/.../aerosolve_input_data_file_ABC.nc
    aerosolve_XYZ_data: ${SCREAM_DATA_DIR}/mam4xx/.../aerosolve_input_data_file_XYZ.nc
    <identifying-string-used-by-aerosolve_xyz.xpp>: /path/to/.../data_file.nc
grids_manager:
  Type: Mesh Free
  geo_data_source: IC_FILE
  grids_names: [Physics GLL]
  Physics GLL:
    type: point_grid
    aliases: [Physics]
    number_of_global_columns:   218
    number_of_vertical_levels:  72
initial_conditions:
  # The name of the file containing the initial conditions for this test.
  Filename: ${SCREAM_DATA_DIR}/init/${aerosolve_IC_file}
  topography_filename: ${TOPO_DATA_DIR}/${EAMxx_tests_TOPO_FILE}
  
  # other variables to pass as input
  num_ABC : 42.0
  num_XYZ : 1.618
  <input_variable> : <value>
# The parameters for I/O control
Scorpio:
  output_yaml_files: ["output.yaml"]
...
```

```yaml title="tests/single-process/mam/aerosolve/output.yaml"
%YAML 1.1
---
filename_prefix: mam4_aerosolve_standalone_output
Averaging Type: Instant
Fields:
  Physics:
    Field Names:
      - answer_field
output_control:
  Frequency: 1
  frequency_units: nsteps
...
```

#### Function-level Verification Tests

Next, because of the highly-technical nature of our new `MAM_Aerosolve::aerosolver()` function, we add a verification test to the 'tests/' subdirectory that co-located with the source code (`src/physics/mam/`).
This subdirectory contains a `CMakeLists.txt` that we modify to add our new test to the build system.
Any required data files should be placed here, and the preferred format for data is [yaml](https://yaml.org), though [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) can also be used.

```c++ title="src/physics/mam/aerosolver_test.cpp"
#include "catch2/catch.hpp"
#include "physics/mam/eamxx_aerosolve_process_interface.hpp"

namespace {
using namespace scream;
TEST_CASE("verify aerosolver","mam")
{
  // initialize the data and test variables
  tol = 1.0e-12;
  input = [...]
  i = [...]
  [...]
  // call the function we are testing
  MAM_Aerosolveaerosolver(input, i, test_answer);
  // apply the REQUIRE() function from catch2/catch.hpp to test
  // whether the resulting test answer is within a set error tolerance
  // from the reference answer
  REQUIRE(std::abs(test_answer - reference_answer) < tol);
}
```

#### Multi-process Validation Tests

Finally, because the Aerosolve process is not called in isolation, we must validate the results generated when Aerosolve is called in combination with other processes.
At time of writing, these tests are divided into the self-describing directories `physics_only` and `dynamics_physics`.
The subdirectories and test directories are named according to the convention of listing the names of included processes and potentially relevant grid information.
Generally speaking, these tests are configured in the same manner as the single-process tests, via the 3 files `CMakeLists.txt`, `input.yaml`, `output.yaml`, and for this reason we omit further details

==More details? Modifications?==

#### Testing Modifications to Existing Code

Note that if you are merely modifying or augmenting existing code, you may not need to create additional tests.
However, one case for which you ***should*** create a test is when adding a new function to existing code.
In this case, a verification test should be added that checks the new functionality for correctness.

### Code, Build, Test, Repeat...

Now that we've written all of the necessary code and there are no errors, we build EAMxx and run our tests.
This build and test workflow can be conducted most easily by using the `test-all-scream` script, located in the `eamxx/scripts/` directory.
We give an overview of common basic use cases here, and we provide full details on the various capabilities and options of `test-all-scream` in [EAMxx Automated Standalone Testing](dev_testing/test_all_scream.md).

##### Only Configure, No Build or Test

**Useful for:** Mid-development testing when it may be useful to confirm or debug changes to CMakeLists.txt, file dependencies, TPL availability, etc.

```shell
./scripts/test-all-scream <required-args> [optional-args] --config-only
```
- regex for specific test
- quick rerun (failed)
- baseline
- separate dirs per arch
- parallel build/test
- verbose

### Testing on Supported Machines

Easy. Everything should work with zero configuration

### Testing Locally

Requires `-l, --local` flag and correctly-configured files in `~/.cime`.
See [Local Configuration Files](dev_testing/test_all_scream.md#local-configuration-files).

### Submit Code to E3SM

- Link to RJ confluence pages
- screenshots?

<!-- ======================================================================= -->

[^automagic-cmake]: In reality, significant effort is involved in designing a project to build "automagically," but the idea is to take that process out of the user's hands.
[^but-mac]: Mac computers with an Apple silicon CPU (e.g., M1, M2...) are not officially supported and present challenges related to handling floating-point exceptions. However, some have had success building on this architecture, and if you achieve a robust build configuration for Mac please reach out to the development team!
[^correct-huh]: Admittedly, "correctness" is a slippery concept, but at the very least we can always test for properties like mass conservation, non-negativity, etc.
