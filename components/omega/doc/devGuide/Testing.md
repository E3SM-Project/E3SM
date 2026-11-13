(omega-dev-testing)=

# Testing Code

Omega development requires both **standalone testing** and **coupled testing**:

- **Standalone testing:** CTest unit tests and the Polaris `omega_pr` regression suite.
- **Coupled testing:** CIME system tests and coupled cases.

```{important}
**Coupled testing is just as important as standalone testing.** Omega is the ocean
component of E3SM, and changes that pass standalone testing can still break
coupled simulations (e.g., causing restart `ERS` or `PEM` tests to fail).
Developers are expected to validate both standalone and coupled behavior on at
least one CPU and one GPU machine before submitting pull requests.
```

For instructions on how to build Omega and run CTest unit tests, see
{ref}`omega-dev-quick-start-build-test`.


## Running the Polaris `omega_pr` Test Suite

The CTest unit tests on their own are not sufficient to ensure that changes
to Omega behave as expected.  Developers also need to run the regression tests
in the form of the
[omega_pr](https://docs.e3sm.org/polaris/main/users_guide/ocean/suites.html#omega-pr-suite)
test suite from Polaris.  The tests in the suite are defined in
[omega_pr.txt](https://github.com/E3SM-Project/polaris/blob/main/polaris/suites/ocean/omega_pr.txt)
in the Polaris repository.

You will run the suite twice: once using the Omega version vendored in Polaris
under `e3sm_submodules/Omega` (the official baseline) and once using your
development branch. This means you will need to
{ref}`build Omega <omega-dev-quick-start-build>` from both sources before
running the suite. You may wish to use the {ref}`omega-dev-quick-start-ctest-util`
to automate the build and run the CTests for each.

```{note}
While the Omega `develop` branch often matches the Polaris submodule
bit-for-bit (b4b), the submodule at `e3sm_submodules/Omega` is the
authoritative baseline for regression comparison.
```

### Prerequisites: build Omega for baseline and PR

Before setting up and running the `omega_pr` suite, you must build Omega twice:

- Build the baseline from the Omega sources vendored by Polaris at
  `<path-to-polaris-repo>/e3sm_submodules/Omega`. Record the CMake build
  directory path as `BASELINE_BUILD` (used below with `polaris suite -p`).
- Build from your development (PR) branch. Record the CMake build directory
  path as `PR_BUILD` (also used with `-p`).

The `-p` option passed to `polaris suite` should point to the CMake build
directory for the corresponding branch (the directory where you ran CMake and
invoked `./omega_build.sh`).

```{note}
Run both baseline and PR suite on the same machine, CPU/GPU partition,
compiler, and Omega build type (Debug/Release) to avoid false positive/negative
diffs.
```

### Setting up the suite

If you have not yet set up Polaris, you will need to follow its
[quick start](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html)
including [cloning the repo](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#set-up-a-polaris-repository-for-beginners),
[configuring the software environment](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#supported-machines),
and [activating the environment](https://docs.e3sm.org/polaris/main/developers_guide/quick_start.html#activating-the-environment).

Then, with the Polaris environment active, you can set up the baseline
`omega_pr` suite as follows:
``` bash
BASELINE_BUILD=/path/to/build/from/polaris/e3sm_submodules/Omega
POLARIS_SCRATCH=/path/to/scratch/polaris_testing/
polaris suite -c ocean -t omega_pr --model omega \
  -w $POLARIS_SCRATCH/baseline_omega_pr \
  -p $BASELINE_BUILD
```
where `$POLARIS_SCRATCH` is a directory within your scratch space that is
useful for setting up and running the Polaris test suites and
`$BASELINE_BUILD` is the directory where you built Omega from the Omega
sources located at `<polaris-repo>/e3sm_submodules/Omega`, or where the
Polaris CTest utility built it for you (e.g. `build_omega/build_chrysalis_intel`).

Next, you can set up the `omega_pr` for your development branch, using the
baseline set up above:

``` bash
PR_BUILD=/path/to/omega/pr/branch/build
POLARIS_SCRATCH=/path/to/scratch/polaris_testing/
polaris suite -c ocean -t omega_pr --model omega \
    -b $POLARIS_SCRATCH/baseline_omega_pr \
    -w $POLARIS_SCRATCH/my_branch_omega_pr \
  -p $PR_BUILD
```
Here, `$POLARIS_SCRATCH` is the same scratch directory as above,
`$PR_BUILD` is the directory where you build Omega from your
development branch (or where the CTest utility built for you), and
`my_branch_omega_pr` can be any name that helps you keep track of what you
were testing.

### Running the suite

First, run the baseline, e.g. for Slurm systems:
``` bash
cd $POLARIS_SCRATCH/baseline_omega_pr
sbatch job_script.omega_pr.sh
```
Note the job ID for the next step.  Next, run the suite for the test branch,
waiting until the baseline is done:
``` bash
cd $POLARIS_SCRATCH/my_branch_omega_pr
sbatch --kill-on-invalid-dep=yes --dependency=afterok:$JOBID job_script.omega_pr.sh
```
where `$JOBID` is the JOB ID of the baseline slurm job.

Monitor the jobs (e.g. `squeue -j $JOBID`) and take a look at combined
stdout/stderr log in `polaris_omega_pr.o$JOBID` in each work directory to see
whether the tests passed.  Report the results as part of any (non-trivial)
Omega PR.

```{note}
For other schedulers (e.g. PBS on Aurora), take a look at the machine's
documentation for equivalent commands.
```

### Handling test failures

What should you do if you see test failures?

#### Test execution errors for baseline (Polaris `e3sm_submodules/Omega`)

If you see test execution failures in the baseline run from Polaris
`e3sm_submodules/Omega`,
these need to be reported as [issues](https://github.com/E3SM-Project/Omega/issues)
on the Omega repo, preferably mentioning a relevant Omega developer who may
be able to investigate the cause of the failure.  Please provide the contents
of the relevant log file for the failing test(s) and optionally the location of
the tests (and log files) if other Omega developers will have access.

#### Test execution errors for your branch

If you see execution failures on the tests for your development branch, please
take a look at the log files and attempt to debug the issues yourself first,
then reach out to the Omega team if you need help.  In many cases, the test
failures may related to changes that are required in Polaris, see
{ref}`omega-dev-testing-update-polaris` below.

#### Diffs between your branch and baseline

If diffs are expected:
- Document the diffs (in Omega PR comment)
- Update `e3sm_submodules/Omega` in Polaris (after merge of Omega PR)

If *not* expected:
- Debug and eliminate

If you see errors related to comparisons with the baseline, you will need to
determine if these are expected based on the code changes or not.  If they
are expected, you will need to document (in a PR comment) how large the
differences are in the tests that show differences.  In addition, you will need
to make a pull request to the Polaris repository to update the Omega submodule
in `e3sm_submodules/Omega` once your branch has been merged so that the changes
become the new baseline in Polaris.

```{note}
We do not update `e3sm_submodules/Omega` after every branch is merged into
`develop` on Omega to reduce maintenance burden.
```

If you see unexpected differences between the results for your development
branch and those from the baseline, you should attempt to debug the cause yourself
and then reach out to the Omega team if you need help.

You may want to look at Polaris'
[bisect](https://github.com/E3SM-Project/polaris/blob/main/utils/bisect/README.md)
utility to help debug which commit introduced the unexpected differences.


(omega-dev-testing-update-polaris)=

### Updating Polaris

Many Omega changes are likely to also require Polaris changes.  In such cases,
it makes sense to co-develop your Omega branch and a Polaris branch.  The
Omega branch will need to be merged first, then the `e3sm_submodules/Omega`
submodule will need to be updated to include the Omega changes as part of
the Polaris pull request.

Here are some examples of Omega changes that will require corresponding
Polaris changes.

#### Adding new Omega dimensions or variables

If you add new Omega variables that have a corresponding MPAS-Ocean name,
you need to add the mapping between the names to
[mpaso_to_omega.yaml](https://github.com/E3SM-Project/polaris/blob/main/polaris/ocean/model/mpaso_to_omega.yaml).
For more details, see the Polaris
[Ocean Framework documentation](https://docs.e3sm.org/polaris/main/developers_guide/ocean/framework.html).

For example:
```yaml
dimensions:
  Time: time
  nVertLevels: NVertLayers

variables:
  temperature: Temperature
  salinity: Salinity
  ssh: SshCell
```
In each case, the keys are the MPAS-Ocean names and the values are the
corresponding Omega names.

#### Updates to Default.yml

If you add or remove sections, config options, streams, etc. in
[Default.yml](https://github.com/E3SM-Project/Omega/blob/develop/components/omega/configs/Default.yml),
you may need to make corresponding changes to the YAML files for Polaris
test cases for them to continue to work as expected.

Whenever possible, we try to use a mapping between MPAS-Ocean namelist options
and Omega config options within
[mpaso_to_omega.yaml](https://github.com/E3SM-Project/polaris/blob/main/polaris/ocean/model/mpaso_to_omega.yaml).
For example:
```yaml
config:
- section:
    time_management: TimeIntegration
  options:
    config_start_time: StartTime
    config_stop_time: StopTime
    config_run_duration: RunDuration
    config_calendar_type: CalendarType
```

If you have added new Omega config sections or options that have corresponding
MPAS-Ocean namelist sections or options, please add them to the map. In each
case, the keys are the MPAS-Ocean names and the values are the corresponding
Omega names.

Frequently, there is no MPAS-Ocean equivalent to an Omega config option, and
we currently do not try to map streams from MPAS-Ocean to Omega (their
syntax is too different).  In such cases, we have do define Omega YAML config
options and streams in each test case, and you may need to update each one
with your changes.  Here is an example exerpt from the
[manufactured_solution](https://github.com/E3SM-Project/polaris/blob/main/polaris/tasks/ocean/manufactured_solution/forward.yaml)
test:
```yaml
Omega:
  Tendencies:
    VelDiffTendencyEnable: false
    VelHyperDiffTendencyEnable: false
    UseCustomTendency: true
    ManufacturedSolutionTendency: true
  IOStreams:
    InitialVertCoord:
      Filename: init.nc
    InitialState:
      UsePointerFile: false
      Filename: init.nc
      Contents:
      - NormalVelocity
      - PseudoThickness
    RestartRead: {}
    History:
      Filename: output.nc
      Freq: {{ output_freq }}
      FreqUnits: Seconds
      IfExists: append
      # effectively never
      FileFreq: 9999
      FileFreqUnits: years
      Contents:
      - NormalVelocity
      - PseudoThickness
      - SshCell
```
For more details on updating either the map or the YAML files for individual
tests, see the Polaris
[Ocean Framework documentation](https://docs.e3sm.org/polaris/main/developers_guide/ocean/framework.html).

(omega-dev-cime-testing)=

## Coupled Testing with CIME

Because Omega is the ocean component of E3SM, changes must be validated in
coupled configurations using E3SM's CIME test harness.

```{note}
All CIME creation commands shown below (`create_test`, `create_newcase`, etc.)
are executed from the `cime/scripts/` directory relative to the **top-level
E3SM repository root**, *not* from within the `components/omega/` subdirectory.
```

For general CIME testing concepts, options, and advanced workflows, refer to the
comprehensive [CIME System Testing documentation](https://esmci.github.io/cime/versions/master/html/system_testing.html).
This section focuses on the specific configurations, conventions, and knowledge an
Omega developer needs.

### Why Coupled Testing is Essential

Coupled testing exercises critical behaviors that cannot be validated in
standalone runs:

- **Flux exchanges:** Validates ocean-atmosphere and ocean-ice boundary flux
  exchanges (`SfcCoupling` and the MCT/MOAB bridges).
- **Coupler clock stepping:** Ensures ocean time stepping advances correctly
  under coupler clock alarms and coupling intervals.
- **Exact restart reproducibility (`ERS`):** Verifies bit-for-bit restart
  capability under CIME driver control.
- **PE layout independence (`PEM`/`PET`):** Verifies bit-for-bit results across
  different processor element (PE, i.e., MPI task) layouts and thread counts.

(omega-cime-env)=

### Sourcing the CIME Environment

```{important}
**Required environment for running Omega under CIME:**
Prior to Omega being added to E3SM, developers did not need to source a special
environment before running CIME commands; legacy components relied exclusively on
the host system's base Python interpreter.

However, Omega uses YAML configuration files instead of traditional Fortran namelists.
To generate `omega.yml`, CIME invokes Omega's `components/omega/cime_config/buildnml`
in-process using Python. This requires the `pyyaml` package, which is not part of the
Python standard library.

This generation step occurs during case setup and inside batch jobs during `case.run`
(including between segments of `ERS` and `PEM` tests). Sourcing the machine's environment
script beforehand ensures `pyyaml` remains available at all stages of execution.

A proposed standard developer environment ([E3SM Discussion #8906](https://github.com/E3SM-Project/E3SM/discussions/8906))
aims to make sourcing an environment a standard, required step for running
all E3SM v4 components. Until that is fully adopted in E3SM, developers
**must explicitly source the CIME environment script for their machine before invoking
any CIME command**:
```

| Machine | Environment Script |
| :--- | :--- |
| aurora | `source /lus/flare/projects/E3SMinput/soft/e3sm-unified/load_latest_cime_env.sh` |
| chrysalis | `source /lcrc/soft/climate/e3sm-unified/load_latest_cime_env.sh` |
| frontier | `source /ccs/proj/cli115/software/e3sm-unified/load_latest_cime_env.sh` |
| perlmutter (`pm-cpu`, `pm-gpu`) | `source /global/common/software/e3sm/anaconda_envs/load_latest_cime_env.sh` |

### Omega Compsets and Grid Compatibility

Omega compsets are grouped into three main tiers (defined in
`components/omega/cime_config/config_compsets.xml` and
`cime_config/allactive/config_compsets.xml`):

- **C Cases (Data Atmosphere & Data Sea Ice):** `COMEGA-IAF`, `COMEGA-JRA1p5`
  - Lightweight tests for atmospheric and runoff forcing ingestion and ocean
    state evolution without active ice or atmosphere components.
- **G Cases (Active Sea-Ice MPASSI):** `GOMEGA-IAF`, `GOMEGA-JRA1p5`
  - Tests bidirectional ice-ocean coupling, ice-ocean fluxes, and frazil ice.
- **B Cases (Fully Coupled):** `WCYCL1850NS-OMEGA`
  - Fully active simulation across all Earth system components (EAM, ELM,
    MOSART, MPASSI, OMEGA).

```{important}
**Data atmosphere (DATM) grid and forcing compatibility:**
In C and G compsets, the atmosphere component is a data atmosphere (DATM) reading
prescribed atmospheric forcing files, *not* a prognostic atmosphere component
(such as EAM, which runs on spectral element grids like `ne4pg2` or `ne30pg2`).
For DATM, the atmosphere grid identifier specifies the resolution of the
prescribed input data files:
- `IAF` (Interannual Forcing) compsets use the `T62` atmospheric forcing grid (e.g., `T62_oQU240`).
- `JRA1p5` compsets use the `TL319` atmospheric forcing grid (e.g., `TL319_EC30to60E2r2`).

Attempting to pair `T62` with a `JRA1p5` compset (or `TL319` with `IAF`) is
invalid and will fail during case setup or input data checking. Refer to the
`e3sm_omega_developer` test suite in
[`cime_config/tests.py`](https://github.com/E3SM-Project/Omega/blob/develop/cime_config/tests.py)
for verified compset, resolution, and test modifier combinations (such as
`.omega-jra_1958`).
```

### Running Coupled Tests with `create_test`

Automated coupled tests are run using CIME's
[`create_test`](https://esmci.github.io/cime/versions/master/html/system_testing.html#run-a-individual-tests)
script located in `cime/scripts`.

#### Test Name Syntax

CIME encodes the test configuration into a standardized name format (see
[CIME Test Name Syntax](https://esmci.github.io/cime/versions/master/html/system_testing.html#test-name-syntax)):

```text
TESTTYPE[_MODIFIERS].GRID.COMPSET[.MACHINE_COMPILER][.GROUP-TESTMODS]
```

- **`TESTTYPE`:** The core test type without underscores (e.g., `ERS`, `PEM`, `SMS`).
  See [CIME Test Types](https://esmci.github.io/cime/versions/master/html/system_testing.html#test-type).
- **`_MODIFIERS`:** Optional modifiers prefixed with underscores to alter settings
  (e.g., `_D` enables `DEBUG=TRUE`, `_Ld3` sets run length to 3 days).
  See [CIME Modifiers](https://esmci.github.io/cime/versions/master/html/system_testing.html#modifiers).
- **`GRID`:** Atmosphere and ocean grid alias or longname (e.g., `T62_oQU240`, `ne4pg2_oQU240`).
- **`COMPSET`:** The compset alias (e.g., `COMEGA-IAF`, `GOMEGA-IAF`, `WCYCL1850NS-OMEGA`).
- **`MACHINE_COMPILER`:** Optional machine and compiler joined by an underscore
  (e.g., `chrysalis_intel`, `pm-gpu_gnugpu`). Can alternatively be passed via CLI flags
  (`-m <machine> -c <compiler>`).
- **`GROUP-TESTMODS`:** Optional directory with namelist or XML modifications
  (see [CIME Test Mods](https://esmci.github.io/cime/versions/master/html/system_testing.html#test-mods)).
  Omega's test mods are located under
  `components/omega/cime_config/testdefs/testmods_dirs/omega/`
  (for example, `omega-jra_1958` sets the atmospheric and runoff forcing start/end
  years to 1958 for JRA-55 simulations).

```{note}
**Machine and compiler detection:**
Specifying the machine and compiler is **optional**. If omitted, `create_test`
automatically detects the host machine and uses that machine's default compiler.

The primary exception is platforms like **NERSC Perlmutter**, where both `pm-cpu`
and `pm-gpu` machine profiles are accessible from the same login environment.
In that case, you must deliberately specify the machine (e.g., append `.pm-gpu`
or pass `-m pm-gpu`) to target GPUs rather than CPU nodes.
```

#### Recommended Tests: Exact Restart (`ERS`) and PE-Count (`PEM`)

Rather than running an extensive test suite across multiple machines, the best
bet for developers is to run targeted **`ERS` (Exact Restart)** and **`PEM`
(PE-count reproducibility)** tests:

- **`ERS` (Exact Restart):** Runs a continuous simulation and an equivalent
  segmented run that stops and restarts mid-run, comparing results bit-for-bit.
  This catches issues with state persistence, boundary conditions, and checkpointing.
- **`PEM` (PE-Count Reproducibility):** Runs the same configuration using two
  different processor element counts (PEs, i.e., MPI tasks) and domain
  decompositions, verifying bit-for-bit identity. This catches bugs in domain
  decomposition, halo exchanges, and reductions.
- **Debug modifier (`_D`):** Using debug versions (such as `ERS_D` or `PEM_D`)
  is strongly recommended. The `_D` modifier enables compiler assertions, array
  bounds checking, and floating-point exception (FPE) trapping, catching latent
  numerics or memory errors early.

```{important}
**Critical on GPU machines and compilers:**
Running `ERS` and `PEM` tests (especially with `_D`) is particularly vital on
GPU machines and compilers (e.g. `pm-gpu`, `frontier`, `aurora`). GPU builds
are sensitive to asynchronous kernel synchronization, race conditions, memory
space transfers, and decomposition-dependent reductions. These issues often
pass short smoke runs without symptom, but reliably fail bit-for-bit restart
(`ERS`) or PE-count (`PEM`) checks.
```

To run targeted tests from `cime/scripts` (relative to the E3SM repo root):

```bash
# 1. Source CIME environment on Perlmutter
source /global/common/software/e3sm/anaconda_envs/load_latest_cime_env.sh

# 2. Navigate to cime/scripts from the E3SM repository root
cd cime/scripts

# Exact restart test with debug flags (machine & compiler auto-detected):
./create_test ERS_D.T62_oQU240.COMEGA-IAF

# Explicitly specifying pm-gpu at NERSC (where login node supports both pm-cpu and pm-gpu):
./create_test ERS_D.T62_oQU240.COMEGA-IAF.pm-gpu_gnugpu
# (or equivalently using CLI flags: ./create_test ERS_D.T62_oQU240.COMEGA-IAF -m pm-gpu -c gnugpu)

# PE-count reproducibility test with debug flags (crucial on GPUs):
./create_test PEM_D.T62_oQU240.COMEGA-IAF

# Active sea-ice exact restart test:
./create_test ERS.T62_oQU240.GOMEGA-IAF
```

````{tip}
**Developer workflow: `--no-run` and interactive execution**
Passing `--no-run` to `create_test` sets up and builds the test case without
submitting the run phase to the batch system:

```bash
./create_test ERS_D.T62_oQU240.COMEGA-IAF --no-run
```

This is especially helpful when you want to:
- Make case modifications (e.g., XML settings via `./xmlchange` or Omega
  runtime configuration edits in `user_nl_omega`) before running.
- Run the test from an interactive compute node allocation (such as `salloc`),
  which requires running `./case.submit --no-batch` from within the generated test
  case directory.
````

#### The `e3sm_omega_developer` Test Suite

The `e3sm_omega_developer` test suite defined in
[`cime_config/tests.py`](https://github.com/E3SM-Project/Omega/blob/develop/cime_config/tests.py)
is maintained primarily for automated nightly testing in E3SM and serves as a
living catalog of valid compset, grid, and modifier combinations. Developers
are not expected to run the entire suite for every PR given the variety of
supported platforms. However, it can be run when comprehensive multi-case
coverage is desired on a supported machine:

```bash
cd cime/scripts
./create_test e3sm_omega_developer
# Or explicitly specifying machine and compiler if desired:
./create_test e3sm_omega_developer -m <machine> -c <compiler>
```

#### Monitoring, Diagnosing, and Iterating on Failed Tests

Each test generated by `create_test` creates a standalone case directory under
your configured test scratch root (named after the test, e.g.,
`ERS_D.T62_oQU240.COMEGA-IAF.<test-id>`). Test status across all phases is
tracked in `TestStatus` and `TestStatus.log`.

If a test fails, you do **not** need to set up a new manual case to diagnose or
fix the issue:
1. Navigate directly into the generated test case directory.
2. Inspect `TestStatus` to identify the failing phase (e.g., `MODEL_BUILD`, `RUN`, `COMPARE_base_rest`).
3. Check the logs:
   - For build errors: inspect `bld/ocn.bldlog.*`.
   - For runtime or comparison errors: inspect logs in the test's `run/` directory
     (`e3sm.log.*`, `ocn.log.*` , `cpl.log.*`).
4. Fix the bug directly in the repository source code checkout used to set up the test.
5. Rebuild and resubmit right from within the existing test case directory:
   ```bash
   # Recompile only modified components
   ./case.build -m

   # Rerun the test
   ./case.submit
   ```

### Interactive Coupled Cases

For active feature development or exploring new configurations, creating an
interactive case from scratch on Chrysalis is straightforward:

```bash
# 1. Source the CIME environment on Chrysalis
source /lcrc/soft/climate/e3sm-unified/load_latest_cime_env.sh

# 2. From the E3SM repository root, navigate to cime/scripts and create a case on Chrysalis
cd cime/scripts
./create_newcase --case /lcrc/group/e3sm/ac.username/cases/omega_test --compset COMEGA-IAF --res T62_oQU240 --mach chrysalis

# 3. Setup the case (invokes Omega's buildnml to generate omega.yml)
cd /lcrc/group/e3sm/ac.username/cases/omega_test
./case.setup

# 4. (Optional) Enable debug compiler flags or customize configuration
./xmlchange DEBUG=TRUE

# You can customize Omega configuration by editing user_nl_omega:
# (see user guide: {ref}`omega-user-config`)
./preview_namelists

# 5. Build the case
./case.build

# 6. Submit the case
./case.submit
```

- **Fast rebuilds:** Run `./case.build -m` to recompile only modified components
  without rebuilding external libraries.
- **Resetting setup:** Run `./case.setup --reset` after modifying grid resolutions,
  component layouts, or processor element (PE) counts.

## What to Include in Your PR Testing Comment

Please add a PR comment titled "Testing" that summarizes what you ran and the
outcomes. Both standalone and coupled testing must be reported for non-trivial
PRs:

```{important}
At a minimum, testing must be performed on at least **one CPU and one GPU machine/compiler combination**, covering:
1. CTest unit tests
2. Polaris `omega_pr` regression suite (standalone)
3. Targeted coupled CIME tests (such as `ERS` or `PEM` tests)
```

- **CTest unit tests:**
  - Machine(s) and partition/queue (if applicable)
  - Compiler(s), and build type(s) (Debug/Release)
  - Result: either "All tests passed" or a short list of failing tests with
    a one-line note each
- **Polaris `omega_pr` regression suite (standalone):**
  - Baseline run (Polaris `e3sm_submodules/Omega`) build path (`-p`) and workdir (`-w`)
  - PR branch run build path (`-p`) and workdir (`-w`)
  - Machine/partition, compiler, and build type (should match between baseline and PR)
  - Result: "All tests passed" or a concise list of diffs/failures with a brief note
  - Log paths (e.g., `polaris_omega_pr.o$JOBID`) if issues occur
- **Coupled CIME tests:**
  - Machine(s) and compiler(s) (highlighting GPU platforms if applicable)
  - Test name(s) (e.g., `ERS_D.T62_oQU240.COMEGA-IAF`, `PEM_D.T62_oQU240.COMEGA-IAF`)
  - Result: "All tests passed (PASS)" or note on failures/phases
  - Log paths (e.g., `TestStatus.log`) if issues occur
- **Tests added/modified/impacted:**
  - Briefly list any new or updated CTest, Polaris, or CIME tests and what they cover
- **Performance PRs (if applicable):**
  - Link to the relevant PACE experiment and summarize before/after metrics

Use the following template for your PR comment:

```text
## Testing

CTest unit tests
- Machine: <machine>[, <partition>]
- Compiler: <compiler>
- Build type: <Debug|Release>
- Result: All tests passed | Failures (X of Y): <name> (<one-line note>), ...

Polaris omega_pr regression suite (standalone)
- Baseline build (-p): <path/to/build of Polaris e3sm_submodules/Omega>
- Baseline workdir (-w): <path/to/baseline_omega_pr>
- PR build (-p): <path/to/pr/build>
- PR workdir (-w): <path/to/my_branch_omega_pr>
- Machine/partition: <machine>[, <partition>]
- Compiler/build type: <compiler>, <Debug|Release>
- Result: All tests passed | Diffs/Failures: <brief summary>
- Logs (if applicable): <path/to/polaris_omega_pr.o$JOBID>

Coupled CIME tests
- Machine: <machine>
- Compiler: <compiler>
- Test(s): <e.g. ERS_D.T62_oQU240.COMEGA-IAF, PEM_D.T62_oQU240.COMEGA-IAF>
- Result: All tests passed | Failures (X of Y): <test name> (<phase>), ...
- Logs (if applicable): <path/to/TestStatus.log>

Tests added/modified/impacted
- CTest: <name> — <one-line purpose>
- Polaris: <path or name> — <one-line purpose>
- CIME: <test or compset> — <one-line purpose>

Performance (if applicable)
- PACE: <link>
- Before: <metric(s)>
- After: <metric(s)>
- Delta: <summary>
```
