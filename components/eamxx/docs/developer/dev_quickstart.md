# Quick-start Guide for Developers

## Motivation

The intention of this guide is to provide a new developer with the most
straightforward path toward productively building, running, and testing EAMxx.
Specifically, we will discuss building and testing in the EAMxx "standalone"
configuration using the `test-all-eamxx` workflow.
Note that tests related to running the full E3SM model are discussed in
[Full Model Testing](dev_testing/full_model_testing.md).

## Stratosphere-level Overview

At the most basic level, EAMxx is built and run using
[CMake](https://cmake.org), a widely-used tool to "automatically"[^automagic-cmake]
configure, build, and test a project on any arbitrary computing platform.
Thus, one can certainly build and test EAMxx with their own customized
CMake/[CTest](https://cmake.org/cmake/help/latest/manual/ctest.1.html) setup.
However, we would strongly advise against doing this unless you know what
you're doing and have good reasons for choosing this option.

Instead, we recommend using the EAMxx-designed `test-all-eamxx` workflow that
employs a variety of shell and python scripts to automate the series of tasks
common to development.
More details on this are provided below in the section on
[Running EAMxx](#running-eamxx-via-test-all-eamxx).

A detailed, example-based walk-through of modifying and testing new code
may be found in [Testing for Development](dev_testing/testing_for_development.md)
and comprehensive details about the standalone testing tools are located in
[EAMxx Automated Standalone Testing](dev_testing/test_all_eamxx.md)

### Note on Computing Platforms

EAMxx is tested and supported on the
[E3SM-supported machines](https://e3sm.org/model/running-e3sm/supported-machines/),
as well as a handful of other machines.[^supported-machines]
This means that, if you have an account on one of these machines, EAMxx tests
can be run with no extra configuration required.
However, many of the standalone EAMxx tests are able to be run on a
reasonably-equipped laptop or personal workstation.[^but-mac]
This second choice requires some extra configurations[^TPLs-etc] but the effort required
may be worth the potential speedup in development.

## Let's Get Started (Quickly)

!!! Warning "Disclaimer"

    We make use of many code snippets below for ***demonstrative purposes***.
    As such, we provide no guarantee that any of these are working code, and
    the author would be quite surprised if more than a couple of them run
    without edits.
    If you take the time to debug these examples, feel free to submit a
    PR with a working version, and we are sure some future developer will
    appreciate it!

## Running EAMxx via `test-all-eamxx`

The quickest method for a new developer to get EAMxx up and running is to use
the automated configure/build/test workflow provided by `scripts/test-all-eamxx`.
This is a good "smoke test"[^smoke-test-def] to confirm that your development
environment is configured properly and that there are no glaring run-time
errors preventing successful running and testing.

Running via `test-all-eamxx` requires ***at minimum*** 1 to 2 command-line
arguments to work correctly and likely more if your goals are more than modest.
In the following sections, we present the configurations a new developer is
most likely to require.

### Selected Useful Testing Options

We give an overview of common basic use cases here, and we provide full details
on the various capabilities and options of `test-all-eamxx` in
[EAMxx Automated Standalone Testing](dev_testing/test_all_eamxx.md).

!!! Warning "Note"

    For all discussion of running standalone testing via `test-all-eamxx`,
    including the following, we assume the commands are run from the root
    directory of EAMxx. That is, if one begins in the root E3SM directory,
    assume that each code snippet begins with an implicit

    ``` {.shell .copy}
    $ eamxx_root="<some/file/path>/E3SM/components/eamxx"
    $ cd "${eamxx_root}"
    ```

??? Example "Configure Only Without Build or Test"

    **Useful for:** Mid-development testing when debugging compilation, and it
    may be useful to confirm/debug changes to CMakeLists.txt, file dependencies,
    TPL availability, etc.

    ``` {.shell .copy}
    $ ./scripts/test-all-eamxx <required-args> [optional-args] --config-only
    ```

??? Tip "Parallel Build and Test Execution"

    **Note:**

    - This is unnecessary when running on
    [supported
    machines](https://e3sm.org/model/running-e3sm/supported-machines/) because
    efficient defaults are already set within `test-all-eamxx`.
    - If you are testing locally or on an otherwise unsupported machine, it is
    recommended to set these values manually.
        - The reason for this is that the default value for each case is `0`,
        which should theoretically choose the optimal number of threads on
        modern systems, but there is no guarantee this will behave as expected
        or that it will not bog down your system.

    ``` {.shell .copy}
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        # if running multiple test configurations, this will launch those build
        # stacks in parallel \
        --parallel \
        # will compile using M threads (equivalent to a manual 'make -jM') \
        --make-parallel-level <M> \
        # will run tests of the same configuration using N threads \
        --ctest-parallel-level <N> \
    ```

??? Quote "Only Run Requested Test Cases"

    As a first step, you are likely to only be concerned with testing your
    specific modifications, so there is no need to run the full suite of EAMxx tests.
    So, we may restrict the tests that will be run by using an appropriate
    [regular expression](https://en.wikipedia.org/wiki/Regular_expression) ("regex").[^regex-tool]

    ``` {.shell .copy}
    # This example regex will match the exact substring "ash_injection".
    # That is, it will match "ash_injection_validate_test",
    # "test_suite_ash_injection", and "test_ash_injection_quick" but not
    # "ash_in_jection_test" or "ash_injector_test".

    $ test_choice="[^\s]*ash_injection[^\s]*"
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --limit-test-regex "${test_choice}"
    ```

??? Success "Baseline Tests"

    [Baseline tests](dev_testing/test_all_eamxx.md#baseline-tests) compare the
    current state of the code to some "baseline" value to ensure that the
    modifications do not negatively impact the current solutions provided by
    the model.
    An example of this, provided below, is focused on testing the code with our
    new **Ash-Injection**-related modifications as compared to the code before
    we started working.

    ***Note***, however, that any new **Ash-Injection**-specific tests or tests
    involving processes that interact with **Ash-Injection** may give different
    answers for perfectly good reasons.
    Regardless, we would like to establish that our modifications do not cause
    any unexpected interactions, so all remaining baseline tests should pass.

    ``` {.shell .copy}
    # For simplicity, we'll assume we begin at the most up-to-date
    # commit on origin/master
    $ cd "${eamxx_root}"
    $ git checkout master
    $ export bl_loc="${HOME}/eamxx-baselines"
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        # set the location where we will save the baseline files we are about
        # to generate \
        --baseline-dir "${bl_loc}" \
        # this flag indicates we will generate new baselines \
        --generate
    ["one eternity later..."]
    # create a new branch to develop from
    $ git checkout -b jane_dev/eamxx/ash_injection_proc
    ["write all the code..."]
    # Note: the presence of the '--baseline-dir' flag indicates we will
    # compare test results to the latest baselines found in that directory
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        # set the location where we saved the previous baselines \
        --baseline-dir "${bl_loc}"
    ["mam4_ash_injection_standalone test fails... no baseline! (so do the others)"]
    ["homme_standalone test passes! (so do other non-associated tests)"]
    ```

??? Tip ""Quick-rerunning" Tests"

    Rebuilding EAMxx can take quite some time, so `test-all-eamxx` has a couple
    of options to rerun tests quickly and avoid configure or build costs.

    `--quick-rerun`

      - **Useful for:** rerunning when there's no need to reconfigure.
          - No edits made to any `CMakeLists.txt` files.
          - No relevant files renamed or moved to a different directory.
          - You plan to re-run tests that did not necessarily fail on the
          previous run.

    ``` {.shell .copy}
    $ test_choice="[^\s]*ash_injection[^\s]*"
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --limit-test-regex "${test_choice}"
    ["tests run..."]
    ["modify ash_injection_test.cpp..."]
    # no need to reconfigure--let's be quick!
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --quick-rerun
    ```

    `--quick-rerun-failed`

      - **Useful for:** rerunning only failed tests when there's no need to
      reconfigure or, potentially, no need to rebuild.
          - No edits made to any `CMakeLists.txt` files.
          - No relevant files renamed or moved to a different directory.
          - In contrast to above, you only wish to rerun tests that failed the
          previous run.

    ``` {.shell .copy}
    # run all of the eamxx tests
    $ ./scripts/test-all-eamxx <required-args> [optional-args]
    ["tests run..."]
    ["ash_injection-related tests fail... so much output, though :("]
    # no need to reconfigure or rebuild--rerun to make the output easier to parse.
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --quick-rerun
    ```

??? Example "Efficiently Building for Multiple Architectures"

    The default build directory employed by `test-all-eamxx` is
    `${eamxx_root}/ctest-build/<test-configuration>`, where `test-configuration`
    can be `full_debug`, `release`, etc.
    The implication of this is that if we are testing for both CPU and GPU on
    the same machine, the default configuration will result in a GPU build
    overwriting the previous CPU build, or vice-versa.
    This will result in long build times and make testing on both architectures inefficient.

    For this reason, it is recommended to use the `--work-dir` flag to set
    separate build directories for each architecture.

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
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --work-dir "${build_loc_cpu}" \
        --make-parallel-level "${build_threads_cpu}" \
        --ctest-parallel-level "${test_threads_cpu}"

    # launch the gpu build/test
    $ ./scripts/test-all-eamxx <required-args> [optional-args] \
        --work-dir "${build_loc_gpu}" \
        --make-parallel-level "${build_threads_gpu}"
    ```

### Testing on Supported Machines

When testing on a supported E3SM machine, the only configuration required is to
set the proper "machine ID" flag that indicates the machine name and potentially
the desired device to run with (CPU, GPU).
This is because `test-all-eamxx` configures the run according to a
"machine file" that sets an environment known to build and test successfully
and efficiently[^mach-file].

??? Example "Example: Running on Perlmutter"

    ``` {.shell .copy}
    # run tests on cpu
    ["start interactive session on cpu node (see below)..."]
    $ ./scripts/test-all-eamxx --machine pm-cpu
    ["test runs..."]
    ["logout from cpu node..."]
    ["start interactive session on gpu node..."]
    $ ./scripts/test-all-eamxx --machine pm-gpu
    ```

??? Warning "`test-all-eamxx` must be run on a `compute` node"

    The `test-all-eamxx` script assumes it is being run on a compute node
    (i.e., not the `login` or `head` node).

    This can be handled most easily by launching the job manually from an
    interactive session.[^tas-compute-node]

    To illustrate, if the platform's batch scheduler is `slurm`, we could do

    ```shell
    $ salloc --nodes 1 --qos interactive --time 02:00:00 [...]
    <response indicating interactive node name>
    $ ssh <node-id>
    $ ./scripts/test-all-eamxx --machine <machine-id>
    ```

### Testing Locally

Running EAMxx on your local workstation, laptop, or less-trafficked cluster
machine can potentially make your development quicker and less complicated.
However, following the ***No Free Lunches*** principle, there is some fiddly
setup involved in making your life easier.

Testing on a local, or otherwise unknown, machine requires the `--local`
(or `-l`) flag and correctly-configured files in the `~/.cime/` directory.
For a full explanation of how to configure your own machine,
see the pages on [Local Configuration Files](dev_testing/test_all_eamxx.md#local-
configuration-files) and [Third-party Libraries](TPLs.md), but we will briefly summarize here.

??? Example "Example `scream_mach_specs.py` Configuration File"

    The `test-all-eamxx` workflow expects to find a file called
    `scream_mach_specs.py` in the directory `~/.cime/`.
    This file is imported using the utility function `get_all_machines()`
    located in `scripts/machine_specs.py`, and it is used to setup the
    required build environment for your machine. 

    Your local `scream_mach_specs.py` configuration file must be a valid python
    module that defines a `Local` specialization of the `Machine` class,
    that is defined in `machine_specs.py`.
    The parent `Machine` class definition is included below, and we also provide
    a bare-bones example of a local `scream_mach_specs.py` defining the
    `Local` version of the `Machine` class.

    Your local `scream_mach_specs.py` configuration file should specify such
    things as:

    - Compiler paths and/or binary name
    - The number of resources (CPU threads) to be used for building or running
        - `num_bld_res` and `num_run_res`
    - Baselines directory location
    - Commands required to setup the build environment
        - E.g., the paths to the requisite third-party libraries (TPLs).
    - GPU architecture
    - Etc.

    We recommend setting as many of the class variables as you can, because
    (mostly) more information is better than less.

    ```{.shell .copy title="${eamxx_root}/scripts/machine_specs.py:[L58-L88]"}
    ###############################################################################
    class Machine(object):
    ###############################################################################
        """
        Parent class for objects describing a machine to use for EAMxx
        standalone testing.
        """
        concrete = False
        name = ""
        mach_file = ""
        env_setup = []
        gpu_arch = "none"
        num_bld_res = -1
        num_run_res = -1
        batch = ""
        cxx_compiler = "mpicxx"
        c_compiler   = "mpicc"
        ftn_compiler = "mpifort"
        baselines_dir = ""

        @classmethod
        def setup_base(cls,name,num_bld_res=-1,num_run_res=-1):
            cls.name = name
            cls.mach_file = pathlib.Path(EAMXX_DIR) / "cmake" /
                            "machine-files" / (cls.name + ".cmake")
            cls.num_bld_res = num_bld_res if num_bld_res > 0 else get_available_cpu_count()
            cls.num_run_res = num_run_res if num_run_res > 0 else get_available_cpu_count()

        @classmethod
        def uses_gpu (cls):
            return cls.gpu_arch!="none"

    ###############################################################################
    ```

    ```{.shell .copy title="~/.cime/scream_mach_specs.py"}
    from machines_specs import Machine

    ###############################################################################
    class Local(Machine):
    ###############################################################################
        concrete = True
        @classmethod
        def setup(cls):
            # NOTE: do not change this line because `test-all-eamxx -l`
            # specifically looks for a machine named "local".
            super().setup_base("local")

            # NOTE: these are run by the shell (likely bash), so must have
            #       shell syntax, and can use shell environment variables
            #       (e.g, $HOME or ~)
            cls.env_setup = [
                             # source my personal .bashrc 
                             ". ~/.bashrc",
                             # source a custom script to setup my eamxx
                             # environment
                             "source ${HOME}/.cime/eamxx-setup.sh"
                             ]

            # define the python-variable home directory for use below
            home_dir = "/home/arthur_dent"

            # NOTE: the syntax used in the python strings below is called
            #       "formatted string literals" or more commonly "f-strings"

            # provide the path to my baselines directory
            cls.baselines_dir = f"{home_dir}/eamxx-baselines"

            # provide the path to a custom cmake file, used to pass build
            # variables directly to cmake
            cls.mach_file = f"{home_dir}/.cime/eamxx_mach_file.cmake"

            # set a couple of environment variables that are used in the
            # setup script and/or eamxx_mach_file.cmake
            # path to directory containing personally-built libraries
            tpl_path = f"{home_dir}/TPLs"
            # path to personally-built openmpi compilers
            mpi_bindir = f"{tpl_path}/openmpi-4.1.2/install/bin"

            # specific paths to 
            cls.cxx_compiler = f"{mpi_bindir}/mpicxx"
            cls.c_compiler = f"{mpi_bindir}/mpicc"
            cls.ftn_compiler = f"{mpi_bindir}/mpifort"

    ###############################################################################
    ```

<!-- ======================================================================= -->

[^automagic-cmake]: In reality, significant effort is involved in designing a project to build "automagically," but the idea is to take that process out of the user's hands.
[^supported-machines]: If you're feeling brave and want to read some python, you can find details about supported machines and their configurations in [machine_specs.py](https://github.com/E3SM-Project/E3SM/blob/master/components/eamxx/scripts/machines_specs.py) in the `scripts/` directory. Or, for general machine files for E3SM, you can take a look at `E3SM/cime_config/machines/`
[^but-mac]: Mac computers with an Apple silicon CPU (e.g., M1, M2...) are not officially supported and present challenges related to handling floating-point exceptions. However, some have had success building on this architecture, and if you achieve a robust build configuration for Mac please reach out to the development team!
[^TPLs-etc]: In particular, it is likely you will need to manually build/install the required [Third-party libraries](TPLs.md).
[^smoke-test-def]: ***Smoke Test*** is a term used by software developers to describe a type of test that initially indicates whether something is working properly--as in, "flip the switch and see if anything starts smoking." :fingers_crossed:
<!-- doesn't like the separation between footnote and text reference appears
to be triggered by the admonition environment -->
<!-- markdownlint-disable-next-line MD053 -->
[^regex-tool]: For those that are not regex ninjas (yet), there are useful tools only to construct or decode regex strings--e.g., [RegExr](https://regexr.com).
[^mach-file]: If you discover that a supported machine does not successfully build and test, please file an [Issue](https://github.com/E3SM-Project/E3SM/issues) with the tag ***Machine Files***.
Also, no promises the standard configuration is perfectly efficient.
<!-- doesn't like the separation between footnote and text reference appears
to be triggered by the admonition environment -->
<!-- markdownlint-disable-next-line MD053 -->
[^tas-compute-node]: While probably not time-efficient, there is no reason `test-all-eamxx` cannot also be run by submitting to the batch scheduler.
If this is your goal, the exercise is left to the reader.
