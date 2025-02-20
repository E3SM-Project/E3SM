# Quick-start Guide for Developers

## TODO

- Mention full-model, but don't go into detail

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
For our case, and following the current informal naming standard, these will be the `eamxx_aerosolve_process_interface.<x>pp` files, located at `src/physics/mam/`, and trivial skeleton examples are provided below.

=== "eamxx_aerosolve_process_interface.hpp"
    ```c++
    #ifndef EAMXX_MAM_AEROSOLVE_HPP
    #define EAMXX_MAM_AEROSOLVE_HPP

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

    #endif  // EAMXX_MAM_AEROSOLVE_HPP
    ```
=== "eamxx_aerosolve_process_interface.cpp"
    ```c++
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

#### Add the New Process to CMake Build

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

### Adding Verification and Validation Tests

All new code that is added to EAMxx is required to be tested.
At minimum, a process must be **Validation Tested**.
That is, tested against known, trusted data, and this data can be samples taken by real-world instruments or comparable results generated by other code.
Ideally, any new and sufficiently complex functions should also have **Verification Tests** that test for mathematical or scientific "correctness".[^correct-huh]

#### Single-process Validation Test

First, we add a validation test for (only) our new process to 
```
tests/single-process/mam/aerosolve/
```
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

=== "CMakeLists.txt"
    ```cmake
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
    
===  "input.yaml"
    ```yaml
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

=== "output.yaml"
    ```yaml
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

Next, because of the highly-technical nature of our new `MAM_Aerosolve::aerosolver()` function, we add a verification test to the 'tests/' subdirectory, co-located with the source code (`src/physics/mam/`).
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

Finally, because the **Aerosolve** process is not called in isolation, we must validate the results generated when **Aerosolve** is called in combination with other processes.
At time of writing, these tests are divided into the self-describing directories `physics_only` and `dynamics_physics`.
The subdirectories and test directories are named according to the convention of listing the names of included processes and potentially relevant grid information.
Generally speaking, these tests are configured in the same manner as the single-process tests, via the 3 files `CMakeLists.txt`, `input.yaml`, `output.yaml`, and for this reason we omit further details

#### Testing Modifications to Existing Code

Note that if you are merely modifying or augmenting existing code, you may not need to create additional tests, provided there is already a test that targets the piece of code that was edited.
However, one case for which you ***should*** create a test is when adding a new function to existing code.
In this case, a verification test should be added that checks the new functionality for correctness.

### Code, Build, Test, Repeat...

Now that we've written all of the necessary code and there are no errors, we build EAMxx, run our tests, and watch them pass with flying colors.
:white_check_mark:

However, for those of us that are mere mortals, we may need to iterate on this process a bit until everything works perfectly.
In this case, the automated build and test workflow provided by `scripts/test-all-scream` is designed to make this cycle as painless as possible.

#### Selected Useful Testing Options

We give an overview of common basic use cases here, and we provide full details on the various capabilities and options of `test-all-scream` in [EAMxx Automated Standalone Testing](dev_testing/test_all_scream.md).

??? Example "Configure Only Without Build or Test"

    **Useful for:** Mid-development testing when debugging compilation, and it may be useful to confirm/debug changes to CMakeLists.txt, file dependencies, TPL availability, etc.
    
    ``` {.shell .copy}
    $ ./scripts/test-all-scream <required-args> [optional-args] --config-only
    ```

??? Tip "Parallel Build and Test Execution"

    **Note:**

    - This is unnecessary when running on [supported machines](https://e3sm.org/model/running-e3sm/supported-machines/) because efficient defaults are already set within `test-all-scream`.
    - If you are testing locally or on an otherwise unsupported machine, it is recommended to set these values manually.
        - The reason for this is that the default value for each case is `0`, which should theoretically choose the optimal number of threads on modern systems, but there is no guarantee this will behave as expected or that it will not bog down your system.
    
    ``` {.shell .copy}
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        # if running multiple test configurations, this will launch those build stacks in parallel \
        --parallel \
        # will compile using M threads (equivalent to a manual 'make -jM') \
        --make-parallel-level <M> \
        # will run tests of the same configuration using N threads \
        # (FIXME: interaction with '--parallel' ?) \
        --ctest-parallel-level <N> \
    ```
??? Quote "Only Run Requested Test Cases"

    As a first step, you are likely to only be concerned with testing your specific modifications, so there is no need to run the full suite of EAMxx tests.
    So, we may restrict the tests that will be run by using an appropriate [regular expression](https://en.wikipedia.org/wiki/Regular_expression) ("regex").[^regex-tool]

    ``` {.shell .copy}
    # This example regex will match the exact substring "aerosolve".
    # That is, it will match "aerosolve_validate_test",
    # "test_suite_aerosolve", and "test_aerosolve_quick" but not
    # "aeroplane_solver_test" or "aerosol_coagulation_test".

    $ test_choice="[^\s]*aerosolve[^\s]*"
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --limit-test-regex "${test_choice}"
    ```

??? Success "Baseline Tests"

    [Baseline tests](dev_testing/test_all_scream.md#baseline-tests) compare the current state of the code to some "baseline" value to ensure that the modifications do not negatively impact the current solutions provided by the model.
    An example of this, provided below, is focused on testing the code with our new **Aerosolve**-related modifications as compared to the code before we started working.

    ***Note***, however, that any new **Aerosolve**-specific tests or tests involving processes that interact with **Aerosolve** may give different answers for perfectly good reasons.
    Regardless, we would like to establish that our modifications do not cause any unexpected interactions, so all remaining baseline tests should pass.

    ``` {.shell .copy}
    # For simplicity, we'll assume we begin at the most up-to-date
    # commit on origin/master
    $ cd "${eamxx_root}"
    $ git checkout master
    $ export bl_loc="${HOME}/eamxx-baselines"
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        # set the location where we will save the baseline files we are about to generate \
        --baseline-dir "${bl_loc}" \ 
        # this flag indicates we will generate new baselines \
        --generate
    ["one eternity later..."]
    # create a new branch to develop from
    $ git checkout -b jane_dev/eamxx/aerosolve_proc
    ["write all the code..."]
    # Note: the presence of the '--baseline-dir' flag indicates we will
    # compare test results to the latest baselines found in that directory
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        # set the location where we saved the previous baselines \
        --baseline-dir "${bl_loc}"
    ["mam4_aerosolve_standalone test fails... no baseline! (so do the others)"]
    ["homme_standalone test passes! (so do other non-associated tests)"]
    ```

??? Tip ""Quick-rerunning" Tests"
    
    Rebuilding EAMxx can take quite some time, so `test-all-scream` has a couple of options to rerun tests quickly and avoiding configure or build costs.

    `--quick-rerun`

      - **Useful for:** rerunning when there's no need to reconfigure.
          - No edits made to any `CMakeLists.txt` files.
          - No relevant files renamed or moved to a different directory.
          - You plan to re-run tests that did not necessarily fail on the previous run.

    ``` {.shell .copy}
    $ test_choice="[^\s]*aerosolve[^\s]*"
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --limit-test-regex "${test_choice}"
    ["tests run..."]
    ["modify aerosolve.cpp..."]
    # no need to reconfigure--let's be quick!
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --quick-rerun
    ```

    `--quick-rerun-failed`

      - **Useful for:** rerunning only failed tests when there's no need to reconfigure or, potentially, no need to rebuild.
          - No edits made to any `CMakeLists.txt` files.
          - No relevant files renamed or moved to a different directory.
          - In contrast to above, you only wish to rerun tests that failed the previous run.

    ``` {.shell .copy}
    # run all of the eamxx tests
    $ ./scripts/test-all-scream <required-args> [optional-args]
    ["tests run..."]
    ["aerosolve-related tests fail... so much output, though :("]
    # no need to reconfigure or rebuild--rerun to make the output easier to parse.
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --quick-rerun
    ```

??? Example "Efficiently Building for Multiple Architectures"

    The default build directory employed by `test-all-scream` is `${eamxx_root}/ctest-build/<test-configuration>`, where `test-configuration` can be `full_debug`, `release`, etc.
    The implication of this is that if we are testing for both CPU and GPU on the same machine, the default configuration will result in a GPU build overwriting the previous CPU build, or vice-versa.
    This will result in long build times and make testing on both architectures inefficient.

    For this reason, it is recommended to use the `--work-dir` flag to set separate build directories for each architecture.

    ``` {.shell .copy}
    # we are on a machine with lots of cpu cores, so we will be extra clever
    # and run simultaneously by using only a portion of cores for each build

    $ echo ${total_cores}
      100

    $ build_loc_cpu="ctest-cpu-build"
    $ build_threads_cpu=20
    $ test_threads_cpu=20
    $ build_loc_gpu="ctest-gpu-build"
    $ build_threads_gpu=40

    # Launch the cpu build/test
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --work-dir "${build_loc_cpu}" \
        --make-parallel-level "${build_threads_cpu}" \
        --ctest-parallel-level "${test_threads_cpu}"

    # launch the gpu build/test
    $ ./scripts/test-all-scream <required-args> [optional-args] \
        --work-dir "${build_loc_gpu}" \
        --make-parallel-level "${build_threads_gpu}"
    ```

#### Testing on Supported Machines

When testing on a supported E3SM machine, the only configuration required is to set the proper "machine ID" flag indicating the machine name, and potentially desired device type (CPU, GPU).
This is because `test-all-scream` configures the run according to a "machine file" that sets an environment known to build and test successfully and efficiently[^mach-file].

??? Example "Example: Running on Perlmutter"

    ``` {.shell .copy}
    # run tests on cpu
    ["start interactive session on cpu node (see below)..."]
    $ ./scripts/test-all-scream --machine pm-cpu
    ["test runs..."]
    ["logout from cpu node..."]
    ["start interactive session on gpu node..."]
    $ ./scripts/test-all-scream --machine pm-gpu
    ```

==FIXME: is the below correct?==

??? Warning "Must Be Run Interactively on a Machine with Batch Scheduling"

    The `test-all-scream` script assumes it is being run on a single, (interactive?) node.
    This can be achieved by launching the job manually from an interactive session.
    To illustrate, if the platform's batch scheduler is `slurm`, we would do

    ```shell
    $ salloc --nodes 1 --qos interactive --time 02:00:00 [...]
    <response indicating interactive node name>
    $ ssh <node-id>
    $ ./scripts/test-all-scream --machine <machine-id>
    ```

[Link to WIP](../WIP.md)

#### Testing Locally

Requires `-l, --local` flag and correctly-configured files in `~/.cime`.
See [Local Configuration Files](dev_testing/test_all_scream.md#local-configuration-files).

#### Suggested (Very Thorough and Conservative) Testing Progression

To fully cover all the potential interactions of the **Aerosolve** process, a reasonable series of tests is provided in the example script below.

**Note:** There is no reason all of this could not be done at once, but we write the script in favor of being clear and explicit.

??? Abstract "Thorough and Careful Testing Script"

    For the purposes of this example, we assume we are on a non-supported local workstation named "bingo-pajama" that has 64 CPU cores and an NVIDIA GPU.
    We also assume that we have working configuration files in `~/.cime/` for both CPU and GPU build configurations.

    ``` {.shell .copy}
    #!/bin/bash

    # machine-related variables
    machine_name="bingo-pajama"
    total_ncores=64
    id_cpu="${machine_name}-cpu"
    id_gpu="${machine_name}-gpu"

    # name of log file for conveniently capturing output
    outfile_prefix="aerosolve_buildNtest"
    ofile_cpu="${outfile_prefix}-${id_cpu}.log"
    ofile_gpu="${outfile_prefix}-${id_gpu}.log"

    # parallelism variables
    make_ncores=60
    test_ncores=60
    pbuild="-p --make-parallel-level ${make_ncores}"
    ptest="--ctest-parallel-level ${test_ncores}"

    # flags to assist with testing/debugging
    vflag="--extra-verbose"
    # limit test configuration to double-precision debug
    test_config="-t dbg"
    bl_loc="${HOME}/eamxx-baselines"

    # navigate to the root directory of eamxx
    cd "${eamxx_root}"
    # switch to the master branch
    git checkout master
    # pull so we are at the latest commit (origin/HEAD)
    git pull origin master
    # update the submodules (the below is typically unnecessary, but thorough)
    git submodule sync --recursive && \
      git submodule update --init --recursive --progress

    # flags common to all tests
    base_flags="${pbuild} ${vflag} ${test_config}"
    # flags common to this stage of tests
    test_base_flags="${base_flags} --baseline-dir ${bl_loc} --generate"
    # general test flags per architecture
    cpu_flags="--work-dir ${PWD}/${id_cpu} --local ${id_cpu} ${ptest}"
    gpu_flags="--work-dir ${PWD}/${id_gpu} --local ${id_gpu}"
    # specific flags for this test per architecture
    cpu_test_flags="${test_base_flags} ${cpu_flags}"
    gpu_test_flags="${test_base_flags} ${gpu_flags}"

    # generate baselines for cpu and gpu
    echo "TEST FLAGS = ${cpu_test_flags}"
    ./scripts/test-all-scream ${cpu_test_flags}
    echo "TEST FLAGS = ${gpu_test_flags}"
    ./scripts/test-all-scream ${gpu_test_flags}

    # Verify the 'aerosolver()' function works as expected by running only that test
    test_choice="[^\s]*aerosolver_test[^\s]*"
    test_base_flags="${base_flags} --limit-test-regex ${test_choice}"
    cpu_test_flags="${test_base_flags} ${cpu_flags}"
    gpu_test_flags="${test_base_flags} ${gpu_flags}"

    # run verification tests
    echo "TEST FLAGS = ${cpu_test_flags}"
    # now that we're running tests that may fail, use 'tee' to
    # duplicate stdout and stderr to file so debugging is slightly easier
    ./scripts/test-all-scream ${cpu_test_flags} |& tee "${ofile_cpu}"
    echo "TEST FLAGS = ${gpu_test_flags}"
    ./scripts/test-all-scream ${gpu_test_flags} |& tee "${ofile_gpu}"

    # Validate the Aerosolve process generates the solutions we expect
    test_choice="[^\s]*mam4_aerosolve_standalone[^\s]*"
    test_base_flags="${base_flags} --limit-test-regex ${test_choice}"
    cpu_test_flags="${test_base_flags} ${cpu_flags}"
    gpu_test_flags="${test_base_flags} ${gpu_flags}"
    # run single-process verification tests
    echo "TEST FLAGS = ${cpu_test_flags}"
    ./scripts/test-all-scream ${cpu_test_flags} |& tee "${ofile_cpu}"
    echo "TEST FLAGS = ${gpu_test_flags}"
    ./scripts/test-all-scream ${gpu_test_flags} |& tee "${ofile_gpu}"

    # validate all the processes that interact with aerosolve generate
    # expected solutions
    aerosolve_multiproc_tests=(
                               # physics_only
                               shoc_mam4_aerosolve
                               p3_mam4_aerosolve
                               # dynamics_physics
                               homme_shoc_cld_p3_mam4_aerosolve_rrtmgp
                               homme_shoc_cld_spa_p3_rrtmgp_mam4_aerosolve
                               )
    for tchoice in ${aerosolve_multiproc_tests[@]}; do
      test_choice="[^\s]*${tchoice}[^\s]*"
      test_base_flags="${base_flags} --limit-test-regex ${test_choice}"
      cpu_test_flags="${test_base_flags} ${cpu_flags}"
      gpu_test_flags="${test_base_flags} ${gpu_flags}"
      # run single-process verification tests
      echo "TEST FLAGS = ${cpu_test_flags}"
      ./scripts/test-all-scream ${cpu_test_flags} |& tee "${ofile_cpu}"
      echo "TEST FLAGS = ${gpu_test_flags}"
      ./scripts/test-all-scream ${gpu_test_flags} |& tee "${ofile_gpu}"
    done
    # NOTE: just setting 'test_choice="[^\s]*_mam4_aerosolve^\s]*"' would also do
    #       the trick in this case

    # run full eamxx test suite against pre-Aerosolve baselines for all
    # test configurations (not just 'dbg')
    base_flags="${pbuild} ${vflag} --baseline-dir ${bl_loc}"
    cpu_test_flags="${base_flags} ${cpu_flags}"
    gpu_test_flags="${base_flags} ${gpu_flags}"
    # run the full suite of eamxx baseline tests
    echo "TEST FLAGS = ${cpu_test_flags}"
    ./scripts/test-all-scream ${cpu_test_flags} |& tee "${ofile_cpu}"
    echo "TEST FLAGS = ${gpu_test_flags}"
    ./scripts/test-all-scream ${gpu_test_flags} |& tee "${ofile_gpu}"
    ```

### Submit Code to E3SM

- Link to RJ confluence pages
- screenshots?

<!-- ======================================================================= -->

[^automagic-cmake]: In reality, significant effort is involved in designing a project to build "automagically," but the idea is to take that process out of the user's hands.
[^but-mac]: Mac computers with an Apple silicon CPU (e.g., M1, M2...) are not officially supported and present challenges related to handling floating-point exceptions. However, some have had success building on this architecture, and if you achieve a robust build configuration for Mac please reach out to the development team!
[^correct-huh]: Admittedly, "correctness" is a slippery concept, but at the very least we can always test for properties like mass conservation, non-negativity, etc.
[^regex-tool]: For those that are not regex ninjas (yet), there are useful tools only to construct or decode regex strings--e.g., [RegExr](https://regexr.com).
[^mach-file]: If you discover that a supported machine does not successfully build and test, please file an Issue ==with the tag ***xyz***?==.
Also, no promises the standard configuration is perfectly efficient.
