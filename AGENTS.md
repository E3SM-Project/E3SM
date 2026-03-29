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

# Test naming format: TEST_TYPE.RESOLUTION.COMPSET[.MACHINE][.TESTMOD]
# Example: SMS_D_Ln5_P4.ne4pg2_oQU480.F2010.ghci-oci_gnu
```

Test types for various properties of the model:
- **SMS** - Short smoke test
- **ERS** - Exact restart test
- **SMS_D** - Debug smoke test
- **PEM** - Test bit-for-bit when changing mpi task count
- **PET** - Test bit-for-bit when changing thread count
- **ERS_Ld** - Exact restart test run for lenght of d days

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

Build types tested: sp (single precision), dbg (debug), fpe (floating point exceptions), opt (optimized)

## Development Workflow

1. Developing E3SM requires running on supported machines.
2. E3SM uses gitworkflows but without a 'seen' branch https://www.kernel.org/pub/software/scm/git/docs/gitworkflows.html
3. Development guide: https://e3sm.org/model/running-e3sm/developing-e3sm/
4. PR-based workflow with automated testing
5. Feature branch names should use this pattern: ```<github username>/<source code area or component>/<feature-description>```
6. New features should be coordinated via E3SM management and science plan

## Testing EAMxx Modifications (for AI Agents)

When modifying files under `components/eamxx/**`, use the EAMxx standalone testing
system to verify your changes. This does NOT require a supported HPC machine.

### One-Time Environment Setup

Run the setup script from the repository root:
```bash
source components/eamxx/scripts/setup-copilot-env.sh
```

This script handles everything: installing system packages (apt-get on Debian/Ubuntu,
with spack as fallback), Python packages, and git submodules. It uses `sudo -E` to
preserve proxy environment variables. If SSH is unavailable (common in CI containers),
it automatically switches to HTTPS for GitHub.

If the setup script does not work for your environment, install these manually:
- **System packages:** gfortran, gcc, g++, cmake, make, MPI (openmpi), netcdf-c,
  netcdf-fortran, pnetcdf, boost, yaml-cpp, blas/lapack, hdf5, perl, libcurl, zlib
- **Python packages:** `pip install psutil pyyaml netCDF4 packaging`
- **Git submodules:**
  ```bash
  git submodule update --init --recursive \
      externals/ekat externals/scorpio externals/mam4xx externals/haero \
      components/eam/src/physics/cosp2/external \
      components/eam/src/physics/rrtmgp/external
  git submodule update --init cime
  cd cime && git submodule update --init CIME/non_py/cprnc && cd ..
  ```
- **Environment variable:** `export SCREAM_INPUT_ROOT=$HOME/e3sm-inputdata && mkdir -p $SCREAM_INPUT_ROOT`

### One-Time Configure

Run configure ONCE from the `components/eamxx` directory:
```bash
cd components/eamxx
./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only
```

The default build type is `dbg` (debug with bounds checking). Choose a different
type based on the issue:
- `-t dbg` -- debug with bounds checking (default, catches most bugs)
- `-t sp` -- single-precision debug (for single-precision bugs)
- `-t opt` -- optimized/release build (for release-mode or performance bugs)
- `-t fpe` -- floating-point exception checking (for NaN/Inf issues)

If input data download fails (e.g., `web.lcrc.anl.gov` is unreachable), create
empty placeholder files. Configure only checks that the files exist; the actual data
is only needed when running tests that use it.

### Incremental Build and Test

After each code change, go directly to the build directory and run make/ctest:
```bash
cd components/eamxx/ctest-build/copilot-testing/full_debug
make -j$(nproc)
ctest -j$(nproc)                   # run all tests
ctest -R <test_name_regex>         # run specific test(s)
ctest --rerun-failed               # rerun only failed tests
```

Build type to directory name mapping:
| Flag | Directory name |
|------|---------------|
| `-t dbg` | `full_debug` |
| `-t sp` | `full_sp_debug` |
| `-t opt` | `release` |
| `-t fpe` | `debug_nopack_fpe` |

### Reconfiguring

If you need a different build type, configure again:
```bash
cd components/eamxx
./scripts/test-all-eamxx -m copilot-testing -t sp --config-only
cd ctest-build/copilot-testing/full_sp_debug
make -j$(nproc) && ctest -j$(nproc)
```

If cmake files or build system files changed, reconfigure the existing build type:
```bash
cd components/eamxx
rm -rf ctest-build/copilot-testing/full_debug/CMake*
./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only
```

## Rules for pull request (PR) description
1. PR description should have two parts.  The fist should be a brief description in plain
 text with no markdown or other formatting. Add a line.  The second part should be the full
description with markdown formatting.
2. PR description should use imperative tense and start with a verb like 'Fix' or 'Add'

