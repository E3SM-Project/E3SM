<!-- Line break seems necessary for aesthetic's sake -->
<!-- markdownlint-disable-next-line MD033 -->
# EAMxx Automated Standalone Testing <br> (`test-all-eamxx`)

## Local Configuration Files

<!-- - automatically takes care of spreading out the parallel workload
      - For speed and efficient use of resources
- can build dbg/release at the same time
- single-node testing script
      - Can handle multi-process/thread
      - Assumes everything is on the same node
      - must be launched from where the tests will be run -->

In this section we describe our testing methodology for standalone EAMxx
configurations. We use several types of tests

- **Unit tests** are individual test programs that demonstrate that a small set
  of code performs a single function or a set of related functions. We use
  a C++ unit testing framework called [Catch2](https://catch2-temp.readthedocs.io/en/latest/index.html)
  to implement unit tests.
- **Property (verification) tests** are test programs that configure code that
  demonstrates that a part of EAMxx (for example, an atmospheric physics
  parameterization or the dynamical core) is able to produce an answer that
  satisfies some physical constraint or matches a known solution under specific
  circumstances.
- **Fortran-C++ "bit-for-bit" (BFB) tests** are test programs, often implemented
  as unit tests, that demonstrate that a set of C++ code ported from Fortran
  produces bit-for-bit identical results to its Fortran counterpart, provided
  certain compiler options are enabled (such as "strict" floating-point
  arithmetic).
- **Test Suites** are named collections of tests that can be run on demand using
  the [ctest](https://cmake.org/cmake/help/latest/manual/ctest.1.html) command.

We also support a `test-all-eamxx` configuration that runs all of the
standalone tests for an EAMxx configuration. Note, your current machine
must be known to EAMxx before `test-all-eamxx` will work. A machine can
be made known to EAMxx by editing the eamxx/scripts/machines_specs.py files.
There are some instructions on what to do at the top of this file.

`test-all-eamxx` has an informative help dump, displaying some examples,
important information, and options flags/descriptions
(provided below in the collapsed text box).

``` {.shell .copy}
cd ${e3sm_root}/components/eamxx
./scripts/test-all-eamxx --help
```
<!-- disable long-line warnings since this is the actual output and -->
<!-- doesn't look bad on the page -->
<!-- markdownlint-disable MD013 -->
??? note "Output from `test-all-eamxx --help`"

    <!-- markdownlint-disable-line MD046 -->
    ```shell
    usage: 
    test-all-eamxx <ARGS> [--verbose]
    OR
    test-all-eamxx --help

    EXAMPLES (assumes user is on machine mappy):
        # Run all tests on current machine using the EAMxx-approved
        env for this machine 
        > cd $eamxx_repo/components/eamxx
        > ./scripts/test-all-eamxx -m mappy

        # Run all tests on current machine with default behavior except using your current shell env 
        > cd $eamxx_repo/components/eamxx
        > ./scripts/test-all-eamxx --preserve-env -m mappy

        # Run all tests on current machine with default behavior except using non-default baselines 
        > cd $eamxx_repo/components/eamxx
        > ./scripts/test-all-eamxx -m mappy --baseline-dir=PATH_TO_BASELINES

        # Run all tests on current machine with default behavior except using local baselines 
        > cd $eamxx_repo/components/eamxx
        > ./scripts/test-all-eamxx -m mappy --baseline-dir=LOCAL

        # Run only the dbg test on current machine with default behavior otherwise
        > cd $eamxx_repo/components/eamxx
        > ./scripts/test-all-eamxx -m mappy -t dbg

    Drive `ctest` testing of EAMxx for a complete set of tests. This will be our
    gold standard to determine if the code is working or not on the current platform.
    For batch machines, this script expects to already be on a compute node;
    it does not do batch submissions.

    IMPORTANT: the default behavior of this script *changes your environment*,
               by loading machine-specific modules and setting machine-specific
               env vars. To prevent this behavior, use --preserve-env flag.

    Baselines: By default, test-all-eamxx will not run baseline tests. If you set
    -b AUTO, baseline tests will be done with the pre-existing public baselines in
    the location specified by the machine spec. You can regenerate
    baselines any time by using the -g flag, but be aware this will impact everyone
    if you regenerate the public baselines. You can change the target baseline area
    using -b $path. You can also use the magic word "LOCAL" to have test-all-eamxx
    pick a local directory for you if you want to manage your own baselines within
    the current repo. If -g is provided, no tests will be run; -g means generate only.

    The general workflow for baseline-changing PRs is:
    1) Issue PR
    2) AT will fail with differences in the baseline tests
    3) Review and merge.
    4) Ask JimF or LucaB for a bless to baselines on mappy and weaver.

    If you are developing a branch and baselines tests are failing unexpectedly, it is
    likely that your branch has fallen out of date. You should upstream merge or rebase
    your branch.

    To sum up, the baseline handling for test-all-eamxx should basically match what we
    do for create_test tests, the only difference is that test-all-eamxx does baseline
    comparison tests by default.

    options:
      -h, --help            show this help message and exit
      -cxx, --cxx-compiler CXX_COMPILER
                            C++ compiler (default: None)
      -f90, --f90-compiler F90_COMPILER
                            F90 compiler (default: None)
      --c-compiler C_COMPILER
                            C compiler (default: None)
      -s, --submit          Submit results to dashboard (default: False)
      -p, --parallel        Launch the different build types stacks in parallel
                            (default: False)
      -g, --generate        Instruct test-all-eamxx to generate baselines from
                            current commit. Skips tests (default: False)
      -b, --baseline-dir BASELINE_DIR
                            Directory where baselines should be read from (or
                            written to, if -g is used). Default is None which
                            skips all baseline tests. AUTO means use public
                            baselines. You can also use LOCAL to manage baselines
                            in your local work dir. Lastly, you can provide a path
                            here as well. (default: None)
      -m, --machine MACHINE
                            Provide machine name. This is *always* required. It
                            can, but does not have to, match EAMXX_MACHINE. You
                            can decorate this with compiler info if a machine
                            supports multiple compiler types. This value will
                            be used as the CTEST_SITE for cdash if the tests are
                            submitted. It is expected that an EAMxx machine file
                            exists for this value. (default: None)
      --config-only         In the testing phase, only run config step, skip build
                            and tests (default: False)
      -c, --custom-cmake-opts CUSTOM_CMAKE_OPTS
                            Extra custom options to pass to cmake. Can use
                            multiple times for multiple cmake options. The -D is
                            added for you (default: [])
      -e, --custom-env-vars CUSTOM_ENV_VARS
                            Extra custom environment variables to be used. These
                            will override(if applicable) whatever was found in
                            machine_specs. Each -e flag supports a single env var,
                            so to pass multiple env var, do -e 'KEY1=VALUE1' -e
                            'KEY2=VALUE2' (default: [])
      --preserve-env        Whether to skip machine env setup, and preserve the
                            current user env (useful to manually test new modules)
                            (default: False)
      -t, --test TESTS      Only run specific test configurations, choices='dbg'
                            (debug), 'sp' (debug single precision), 'fpe' (debug
                            pksize=1 floating point exceptions on), 'opt'
                            (release), 'cov' (debug coverage), 'valg' (Release
                            build where tests run through valgrind), 'csm' (debug
                            with compute sanitizer memcheck), 'csr' (debug with
                            compute sanitizer racecheck), 'csi' (debug with
                            compute sanitizer initcheck), 'css' (debug with
                            compute sanitizer synccheck) (default: [])
      -l, --local           Allow to not specify a machine name, and have test-
                            all-eamxx to look for '~/.cime/scream_mach_specs.py'
                            for machine specifications. (default: False)
      -r, --root-dir ROOT_DIR
                            The root directory of the EAMxx src you want to test.
                            Default will be the EAMxx src containing this script.
                            (default: None)
      -w, --work-dir WORK_DIR
                            The work directory where all the building/testing will
                            happen. Defaults to ${root_dir}/ctest-build (default:
                            None)
      --quick-rerun         Do not clean the build dir, and do not reconfigure.
                            Just (incremental) build and test. (default: False)
      --quick-rerun-failed  Do not clean the build dir, and do not reconfigure.
                            Just (incremental) build and retest failed tests only.
                            (default: False)
      --make-parallel-level MAKE_PARALLEL_LEVEL
                            Max number of jobs to be created during compilation.
                            If not provided, use default for given machine.
                            (default: 0)
      --ctest-parallel-level CTEST_PARALLEL_LEVEL
                            Force to use this value for CTEST_PARALLEL_LEVEL. If
                            not provided, use default for given machine. (default:
                            0)
      -x, --extra-verbose   Have ctest run in extra-verbose mode, which should
                            print full output from all tests (default: False)
      --limit-test-regex LIMIT_TEST_REGEX
                            Limit ctest to running only tests that match this
                            regex (default: None)
      --test-level {at,nightly,experimental}
                            Set the testing level. (default: at)
      --test-size {short,medium,long}
                            Set the testing level. Defaults to medium unless the
                            test is cov or mem_check(short is default in those
                            cases). (default: None)

    ```
<!-- markdownlint-enable MD013 -->

If you are unsure of the CMake configuration for you development cycle, one
trick you can use is to run `test-all-eamxx --config-only`
for the `dbg` test and then copy the `cmake` commands from the resulting output.

```shell
cd $eamxx_repo/components/eamxx
./scripts/test-all-eamxx -t dbg -m $machine
# wait for a few seconds*
# Ctrl-C *
# Copy the contents of DCMAKE_COMMAND that was passed to ctest *
# Add "cmake" to beginning of contents and path to eamxx at the end. *
```

Considerations for using `test-all-eamxx`:

- Your machine must be known to our scripts, see above.
- If you try to run commands by-hand (i.e., outside of `test-all-eamxx`,
`cmake`, `make`, `ctest`, etc.), you'll need to
load the eamxx-env into your shell, which can be achieved by running
`cd eamxx/scripts; eval $(./eamxx-env-cmd $machine)`
- `test-all-eamxx` expects to be run from a compute node if you
are on a batch machine.
- You'll need to think about your baseline situation, as many of our
  tests rely on pre-existing baselines. The `-b` option controls the baseline
  location and can have the following values (**Note:** If no `-b` flag is
  passed at all, no baseline testing will be done):
      - `AUTO`
          - A common public baseline area shared by all developers
      - `$path`
          - A specific arbitrary path

## Baseline Tests

***(Work in Progress)***

<!-- ==details go here== -->
<!-- 
## Running EAMxx's Tests with CTest

Before running the tests, generate a baseline file:

```shell
cd $RUN_ROOT_DIR
make baseline
```

The tests will run, automatically using the baseline file, which is located in
the CMake-configurable path `${EAMXX_BASELINES_DIR}`.
By default, this path is set to an invalid string.
If baselines tests are enabled, we check that a valid
path has been provided.

To run all of EAMxx's tests, make sure you're in `$RUN_ROOT_DIR` and type

```shell
ctest -VV
```

This runs everything and reports results in an extra-verbose (`-VV`) manner.

You can also run subsets of the EAMxx tests. For example, to run only the
P3 regression tests (again, from the `$RUN_ROOT_DIR` directory), use

```shell
ctest -R p3_regression
``` -->

<!-- ### Grouping Tests with Labels

We can create groupings of tests by using **labels**. For example, we have a
`driver` label that runs tests for EAMxx's standalone driver. You can see a
list of available labels by typing

```shell
ctest --print-labels
```

To see which tests are associated with a given label (e.g. `driver`), use

```shell
ctest -L driver -N
``` -->

## EAMxx Test Suites

### The `p3_regression` Suite

`p3_regression` uses a baseline file to compare any new or altered
implementations with our P3 Fortran reference implementation. If you're working
on the C++/Kokkos implementation, you can invoke any new tests to the function
`Baseline::run_and_cmp` in
`${EAMxx_SRC_DIR}/components/eamxx/p3/tests/p3_run_and_cmp.cpp`.

If the reference Fortran implementation changes enough that a new baseline file
is required, make sure to let other EAMxx team members know, in order to
minimize disruptions.
