<!-- markdownlint-disable -->
# AGENTS.md

E3SM (Energy Exascale Earth System Model) is a state-of-the-art fully coupled Earth System Model. This guide helps you work effectively with this codebase.

## Project Overview

E3SM (Energy Exascale Earth System Model) is a state-of-the-art fully coupled Earth System Model developed for DOE Leadership Computing Facilities. It includes biogeochemical and cryospheric processes and is designed for high-performance computing clusters.

- **Website:** https://e3sm.org
- **Documentation:** https://docs.e3sm.org/E3SM/
- **DOI:** 10.11578/E3SM/dc.20240930.1

## Build Commands

E3SM uses the CIME (Common Infrastructure for Modeling the Earth) case control system.
Build is done out-of-source and is not under the source directory.
First a CASEDIR is created using create_newcase.  In the CASEDIR, you can run ./case.build
after running ./case.setup.
Supported machines are in cime_config/machines/config_machines.xml and grouped by hostname.
Do not try and create a case or test unless on a supported machine.

```bash
#list supported machines
cd cime/scripts
./query_config --machines
The string "(current)" will appear next to the machine you are on if it is supported.
```

### Creating and Building a Case
Ask the user to specify a CASENAME. Create the case in the top level project directory unless
the case name includes a absolute or relative path.

```bash
# Navigate to CIME scripts directory
cd cime/scripts

# Create a new case
./create_newcase --case ../../CASENAME --compset COMPSET --res RESOLUTION --mach MACHINE

# Setup the case
cd /path/to/case
./case.setup

# Build the code for the case
./case.build

#Find the directory of the build output
./xmlquery EXEROOT

# Submit the case
./case.submit

# Find the directory with executable output
./xmlquery RUNDIR
```

### Rebuild a case
```bash
#Do not clean just build again
./case.build

#if only components or driver were modified its faster to do
./case.build -m

#if the pelayout is modified, do
./case.setup --reset
./case.build
```

### Rerun a case
```bash
# Submit the case
./case.submit
```

### Common Compsets and Resolutions
Use these when testing changes to code made only in the specified component
- elm: COMPSET I1850CNPRDCTCBCTOP, RESOLUTION ne4pg2_ne4pg2
- eam: COMPSET F1850, RESOLUTION ne4pg2_oQU480
- mosart: COMPSET RMOSGPCC, RESOLUTION r05_r05
- eamxx: COMPSET F2010-SCREAMv1, RESOLUTION ne4pg2_ne4pg2
- mpas-ocean: COMPSET CMPASO-NYF, RESOLUTION T62_oQU120
- mpas-seaice: COMPSET DTESTM, RESOLUTION T62_oQU240

If code has been changed in elm and eam, use the F1850 compeset.
If code has been changed in 2 or more other components, use the
coupled case: COMPSET WCYCL1850NS, RESOLUTION ne4pg2_r05_oQU480


## Compsets

Compsets define the component configuration for a simulation. The longname format is:
`TIME_ATM[%phys]_LND[%phys]_ICE[%phys]_OCN[%phys]_ROF[%phys]_GLC[%phys]_WAV[%phys]`

Compset definitions are in `cime_config/allactive/config_compsets.xml` and component-specific files in `components/*/cime_config/config_compsets.xml`.

```bash
#list all compsets
cd cime/scripts
./query_config --compsets all
```

For a created case, the longname is in the README.case file in the case directory.

### Component Abbreviations in Compset Names
- **ATM:** EAM (atmosphere), EAMXX/SCREAM (next-gen), DATM (data), SATM (stub)
- **LND:** ELM (land), SLND (stub)
- **ICE:** MPASSI (MPAS sea ice), CICE (legacy), DICE (data), SICE (stub)
- **OCN:** MPASO (MPAS ocean), DOCN (data), SOCN (stub)
- **ROF:** MOSART (river routing), SROF (stub)
- **GLC:** MALI (land ice), SGLC (stub)
- **WAV:** WW3 (waves), SWAV (stub)

## Grid Resolutions

Grid definitions are in `cime_config/config_grids.xml`. The grid alias format is:
`a%name_l%name_oi%name_r%name_m%mask_g%name_w%name` (atm_lnd_ocn/ice_rof_mask_glc_wav)

```bash
#list all grids
cd cime/scripts
./query_config --grids
```

## Testing

While doing development, create a single case that uses the component
being modified and run ./case.build and ./case.submit as needed to test
code modifications.

Most of E3SM does not have unit tests.

### EAMxx (SCREAM) Standalone Testing
Code modifications in components/eamxx can be tested
in a standalone mode without CIME or a driver.

```bash
cd components/eamxx

# Run all tests on a machine (requires compute node on batch systems)
./scripts/test-all-eamxx -m MACHINE

# Run specific test type (sp, dbg, fpe, opt)
./scripts/test-all-eamxx -m MACHINE -t dbg

# Preserve current environment instead of loading machine-specific modules
./scripts/test-all-eamxx --preserve-env -m MACHINE

# Run with local baselines
./scripts/test-all-eamxx -m MACHINE --baseline-dir=LOCAL
```
### CIME-Based System Tests

```bash
cd cime/scripts

# Run a single test
./create_test TEST_NAME

# Run a test suite
./create_test SUITE_NAME

# Monitor test progress by going to the test directory and tailing the TestStatus.log file

# Test naming format: TESTTYPE[_modifiers].GRID.COMPSET.MACHINE_COMPILER[.TESTMODS]
# Example: SMS_D_Ln5_P4.ne4pg2_oQU480.F2010.pm-cpu_gnu
# Example with a testmod: ERS_Ld3.ne30pg2_r05_IcoswISC30E3r5_gis4to40.IGELM_MLI.chrysalis_gnu.mali-landiceIG
#
# MACHINE_COMPILER is required when a TESTMODS field follows. Test names listed in
# tests.py omit it, so they cannot be pasted directly into create_test -- insert the
# machine_compiler field or the testmod is parsed as the machine.
```

Test types for various properties of the model:
- **SMS** - Short smoke test
- **ERS** - Exact restart test
- **SMS_D** - Debug smoke test
- **PEM** - Test bit-for-bit when changing mpi task count
- **PET** - Test bit-for-bit when changing thread count
- **ERS_Ld** - Exact restart test run for lenght of d days

Before claiming a change is tested, confirm the code path actually executed.
Several driver features are off by default and are not exercised by the
standard tests -- coupler budgets (`BUDGETS`, default FALSE) and auxiliary
history files are examples. Grep the run log for output the new code produces,
or enable the feature explicitly.

To run a test with modified namelist settings, build it without running,
edit the namelist, then submit:

```bash
./create_test TESTNAME -t TESTID --no-run
# append settings to user_nl_cpl (or user_nl_<component>) in the case dir
cd $CASEDIR && ./case.submit
```

Do not use `create_clone` on a test case to do this -- the clone inherits
`env_test.xml` and `TestStatus` from the original and fails on the test-name
mismatch.

### Comparing model output between runs

Use `cprnc` to compare two NetCDF history files field by field. Its path is a
per-machine CIME variable, so get it from the case rather than hardcoding:

```bash
CPRNC=$(./xmlquery --value CCSM_CPRNC)
$CPRNC -m baseline.nc new.nc
```

The summary line reports how many of N fields had non-zero differences, and
`-m` adds per-field RMS and normalized RMS.

Two rules that are easy to get wrong:

- **Only compare runs of the same test.** Different test types stop at
  different times and write history at different dates, so field counts are
  not comparable between, say, an `SMS` and an `ERS_Ld3`. A field differing at
  the last bit may simply not have diverged yet in the shorter run, which reads
  as an improvement that is not real.
- **Save a baseline before re-running.** Re-running a case overwrites its
  `RUNDIR`. If the current output is a reference for a later comparison, copy
  it somewhere first.

### Test Suites (defined in cime_config/tests.py)
- e3sm_land_developer
- e3sm_mosart_developer
- e3sm_land_exeshare

```bash
#get list of test suites
cd cime/CIME/Tools
./list_e3sm_tests

#list compsets in a test suite
cd cime/CIME/Tools
./list_e3sm_tests -t compsets SUITE_NAME

#list compset longnames in a test suite
./list_e3sm_tests -t compsets -l SUITE_NAME
```

Nightly results from test suites appear on https://my.cdash.org/index.php?project=E3SM

### Testing EAMxx Modifications

When modifying files under `components/eamxx/**`, see `components/eamxx/AGENTS.md`
for standalone testing instructions.

## Architecture

### Primary Languages
- **Fortran** (~4300 files): Core model physics and dynamics
- **Python** (~2100 files): CIME infrastructure, testing, configuration
- **C++**: Modern EAMxx/SCREAM atmosphere component, external libraries

### Major Components (in /components/)

| Component | Description |
|-----------|-------------|
| eam | E3SM Atmosphere Model (Fortran) |
| eamxx | Next-gen atmosphere model (C++/Fortran, SCREAM) |
| elm | E3SM Land Model |
| mosart | River routing model |
| mpas-ocean | Ocean model |
| mpas-seaice | Sea ice model |
| mpas-albany-landice | MALI land ice model |
| homme | Spectral element dynamical core |

### Key Infrastructure Directories

- `/cime/` - CIME submodule (case control system)
- `/cime_config/` - E3SM-specific CIME configuration
  - `machines/` - Machine, compiler, batch specifications
  - `machines/config_machines.xml` - supported machines
  - `machines/cmake_macros/` - Machine-specific CMake settings
  - `tests.py` - Test suite definitions
  - `config_grids.xml` - Grid definitions
  - `config_compilers.xml` - Compiler configurations
- `/driver-mct/` - MCT coupler driver
- `/driver-moab/` - MOAB coupler driver
- `/share/` - Shared utilities (timing, RNG, streams)
- `/externals/` - External library submodules

### Coupler drivers

`driver_default = "moab"` in `cime_config/customize/config.py`, so cases use
driver-moab unless opted out. The `_Vmct` test modifier (e.g. `ERS_Vmct_Ld3...`)
selects driver-mct. When changing code shared by both drivers, or code in one
driver that is meant to reproduce the other, test both -- a passing `_Vmct` test
shows mct still runs, not that its answers are unchanged.

### External Dependencies (Git Submodules)

Key submodules in `/externals/`:
- **ekat** - E3SM Kokkos Application Toolkit
- **scorpio** - Parallel I/O library
- **mam4xx** - Modal Aerosol Model (C++)
- **haero** - Aerosol package
- **YAKL** - Yet Another Kernel Language
- **mct** - Model Coupling Toolkit

After cloning, initialize submodules:
```bash
git submodule update --init --recursive --depth=1
```

### Build System

- **Primary:** CMake 3.18+ (top-level CMakeLists.txt in /components)
- **Secondary:** Make (MPAS components, MCT)
- **GPU Support:** CUDA and HIP via cmake_macros configuration

## CI/CD

GitHub Actions workflows in `/.github/workflows/`:
- `eamxx-sa-testing.yml` - EAMxx standalone (gcc-openmp, gcc-cuda)
- `eamxx-v1-testing.yml` - CIME tests with EAMxx/MAM4xx
- `e3sm-gh-testing.yml`  - CIME tests with E3SM cases

Build types tested: sp (single precision), dbg (debug), fpe (floating point exceptions), opt (optimized)

## Development Workflow

1. Testing changes to E3SM requires running on supported machines
2. E3SM uses gitworkflows but without a 'seen' branch https://www.kernel.org/pub/software/scm/git/docs/gitworkflows.html
3. Development guide: https://e3sm.org/model/running-e3sm/developing-e3sm/
4. PR-based workflow with automated testing
5. Feature branch names should use this pattern: ```<github username>/<source code area or component>/<feature-description>```

### Writing about code

When writing about code in github issues or pull requests, use these rules.

Assume that the audience or readership is an E3SM developer with
familiarity at a high level with the overall structure and intention of the
code base.

Never repeat explanations of existing, unchanged features.
Explanations are appropriate where there is a change, particularly one of a
conceptual nature. Explanations are also appropriate to highlight nuance
that informs the current discussion.

Write less; do not pack the same content into denser sentences. Keep
headings, tables and links. Colleagues mostly write unstructured prose, and
structure is an improvement on it. The problem is length.

- **Lead with the answer.** The first two sentences say what you found,
  changed, or propose. Setup and reproduction go last.
- **One point per paragraph, and few paragraphs.** Colleagues write one to
  three per comment; recent AI-written ones ran to eighteen. That gap is
  the complaint. Say each thing once.
- **Do not narrate the mechanism.** The chain of calls, and why the fix is
  right, go in the commit message. Here, say what broke and where to look.
- **Cut clauses that qualify rather than inform**, and any sentence whose
  only job is to justify the one before it. One clause per sentence where
  one will do.
- **Use backticks about half as often as feels natural.** They are for what
  a reader would type or grep. Code blocks hold artifacts you did not
  write, never authored prose.
- **One document, one decision.** Anything still relevant after this merges
  is an issue, not a comment.
- **For a PR description** Should have two parts.  The fist should be a brief description in plain
  text with no markdown or other formatting. Then add "Fixes" for any issues that are closed.
  Then add the PR signature described below.  Then add a line.
  After the line should be the full description with markdown formatting, tables, etc. 
  The description should use imperative tense and start with a verb like 'Fix' or 'Add'.

### PR signature
When the PR is created by the agent, use this signature.
```
*Created by <agent> on @<user>'s behalf.
```
If the changes grew out of the user's own changes or substantial
direction, sign it instead with:

```
*Created by @<user> and <agent> together*
```
## Writing rules for using collaboration or communication software (Github, Slack, etc.)

When posting a comment (not a PR) on someone's behalf, preface it with this on its own line.

```
*Posted by <agent>*
```

If the post grew out of the user's own draft or substantial
direction, preface it instead with:

```
*Drafted by @<user> and <agent> together*
```

Name the agent harness, not the vendor: `Claude Code`, `Codex`, and so on.
