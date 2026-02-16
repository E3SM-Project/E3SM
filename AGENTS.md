# AGENTS.md

E3SM (Energy Exascale Earth System Model) is a state-of-the-art fully coupled Earth System Model. This guide helps you work effectiv
ely with this codebase.

## Project Overview

E3SM (Energy Exascale Earth System Model) is a state-of-the-art fully coupled Earth System Model developed for DOE Leadership Computing Facilities. It includes biogeochemical and cryospheric processes and is designed for high-performance computing clusters.

- **Website:** https://e3sm.org
- **Documentation:** https://docs.e3sm.org/E3SM/
- **DOI:** 10.11578/E3SM/dc.20240930.1

## Build Commands

E3SM uses the CIME (Common Infrastructure for Modeling the Earth) case control system. Building requires a supported machine.
Build is done out-of-source and is not under the source directory.

### Creating and Building a Case

```bash
# Navigate to CIME scripts directory
cd cime/scripts

# Create a new case
./create_newcase --case /path/to/case --compset COMPSET --res RESOLUTION --machine MACHINE

# Setup the case
cd /path/to/case
./case.setup

# Build the case
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

#if only components or driver were modified
./case.build -m
```

### Rerun a case
```bash
# Submit the case
./case.submit
```

### Common Compsets and Resolutions
- **Compsets:** WCYCL1850, F2010, I1850ELM
- **Resolutions:** ne30pg2_r05_IcoswISC30E3r5, ne4pg2_oQU480

## Compsets

Compsets define the component configuration for a simulation. The longname format is:
`TIME_ATM[%phys]_LND[%phys]_ICE[%phys]_OCN[%phys]_ROF[%phys]_GLC[%phys]_WAV[%phys]`

Compset definitions are in `cime_config/allactive/config_compsets.xml` and component-specific files in `components/*/cime_config/config_compsets.xml`.

### Fully Coupled (WCYCL - Water Cycle)
| Alias | Description |
|-------|-------------|
| WCYCL1850 | Pre-industrial (1850) fully coupled |
| WCYCL1950 | Mid-century (1950) fully coupled |
| WCYCL20TR | 20th century transient fully coupled |
| WCYCL1850NS | 1850 without spun-up ocean/ice ICs (for testing) |
| WCYCLSSP585/370/245 | Future scenario projections |

### Atmosphere-Only (F compsets)
| Alias | Description |
|-------|-------------|
| F1850 | 1850 atmosphere with prescribed ocean/ice |
| F1950 | 1950 atmosphere with prescribed ocean/ice |
| F2010 | 2010 atmosphere with prescribed ocean/ice |
| F20TR | 20th century transient atmosphere |
| FSCM | Single column model |

### Land-Only (I compsets)
| Alias | Description |
|-------|-------------|
| IELM | Present-day land with data atmosphere |
| I1850ELM | 1850 land-only |
| I1850ELMCN | 1850 land with carbon-nitrogen cycle |
| I20TRELM | 20th century transient land |
| IELMBC | Land with black carbon |
| IERA5ELM | Land forced by ERA5 reanalysis |

### Component Abbreviations in Compset Names
- **ATM:** EAM (atmosphere), EAMXX/SCREAM (next-gen), DATM (data), SATM (stub)
- **LND:** ELM (land), SLND (stub)
- **ICE:** MPASSI (MPAS sea ice), CICE (legacy), DICE (data), SICE (stub)
- **OCN:** MPASO (MPAS ocean), DOCN (data), SOCN (stub)
- **ROF:** MOSART (river routing), SROF (stub)
- **GLC:** MALI (land ice), SGLC (stub)
- **WAV:** WW3 (waves), SWAV (stub)

### Physics Options (% modifiers)
- `%CMIP6` - CMIP6 configuration
- `%CNPRDCTCBCTOP` - Carbon-nitrogen-phosphorus with RDCTC biogeochemistry
- `%SP` / `%SPBC` - Satellite phenology (with/without black carbon)
- `%PRES` - Prescribed mode
- `%DOM` - Data ocean model mode

## Grid Resolutions

Grid definitions are in `cime_config/config_grids.xml`. The grid alias format is:
`a%name_l%name_oi%name_r%name_m%mask_g%name_w%name` (atm_lnd_ocn/ice_rof_mask_glc_wav)

### Atmosphere Grids
| Grid | Description | Approximate Resolution |
|------|-------------|------------------------|
| ne4np4.pg2 | Cubed-sphere, very coarse | ~7.5° (testing only) |
| ne30np4.pg2 | Cubed-sphere, low-res | ~1° (100 km) |
| ne120np4 | Cubed-sphere, high-res | ~0.25° (25 km) |
| 0.9x1.25 (f09) | Finite-volume | ~1° |
| 1.9x2.5 (f19) | Finite-volume | ~2° |
| 4x5 (f45) | Finite-volume | ~4° |

### Ocean/Ice Grids (MPAS)
| Grid | Description |
|------|-------------|
| oQU480 | Quasi-uniform 480 km (testing) |
| oQU240 | Quasi-uniform 240 km |
| oEC60to30v3 | Eddy-closure 60-30 km variable |
| IcoswISC30E3r5 | Icosahedral with ice-shelf cavities 30 km |
| SOwISC12to60E2r4 | Southern Ocean 12-60 km variable |

### River Routing Grids
| Grid | Description |
|------|-------------|
| r05 | 0.5° river routing (standard) |
| r01 | 0.1° river routing (high-res) |
| rx1 | 1° river routing |

### Common Combined Grid Aliases
| Alias | Components | Use Case |
|-------|------------|----------|
| ne30pg2_r05_IcoswISC30E3r5 | 1° atm + 0.5° river + 30km ocean | Production coupled |
| ne4pg2_oQU480 | Very coarse atm + 480km ocean | Testing/debugging |
| f09_g16 | 1° FV atm + 1° POP ocean | Legacy configurations |
| f19_g16 | 2° FV atm + 1° POP ocean | Coarse testing |
| f09_f09 | 1° FV atm-only | Atmosphere-only runs |
| r05_r05 | 0.5° land/river only | Land-only runs |
| hcru_hcru | 0.5° CRU grid | Land with CRU forcing |
| ELM_USRDAT | User-defined | Custom regional grids |

### EAMxx (SCREAM) Standalone Testing

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

## Testing

### CIME-Based Tests

```bash
cd cime/scripts

# Run a single test
./create_test TEST_NAME --wait

# Test naming format: TEST_TYPE.GRID.COMPSET[.MACHINE][.TESTMOD]
# Example: SMS_D_Ln5_P4.ne4pg2_oQU480.F2010.ghci-oci_gnu
```

Test types:
- **SMS** - Short smoke test
- **ERS** - Exact restart test
- **SMS_D** - Debug smoke test
- **ERS_Ld** - Exact restart test run for lenght of d days

### Test Suites (defined in cime_config/tests.py)
- e3sm_land_developer
- e3sm_mosart_developer
- e3sm_land_exeshare

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
git submodule update --init --recursive
```

### Build System

- **Primary:** CMake 3.18+ (top-level CMakeLists.txt in /components)
- **Secondary:** Make (MPAS components, MCT)
- **GPU Support:** CUDA and HIP via cmake_macros configuration

## CI/CD

GitHub Actions workflows in `/.github/workflows/`:
- `eamxx-sa-testing.yml` - EAMxx standalone (gcc-openmp, gcc-cuda)
- `e3sm-gh-ci-cime-tests.yml` - CIME integration tests
- `eamxx-sa-coverage.yml` - Code coverage
- `eamxx-gh-clang-format.yml` - Code formatting

Build types tested: sp (single precision), dbg (debug), fpe (floating point exceptions), opt (optimized)

## Development Workflow

1. Development guide: https://e3sm.org/model/running-e3sm/developing-e3sm/
2. E3SM requires running on supported machines (see https://e3sm.org/model/running-e3sm/supported-machines/)
3. PR-based workflow with automated testing
4. Cosmetic-only changes from non-staff are generally not accepted
5. New features should be coordinated via E3SM management and science plan

## Rules for pull requests
1. PR description should be brief and in plain text with no markdown and with lines wrapped at 80 characters.
2. PR description should use imperative tense and start with a verb like Fix or Add
3. After the plan text, you can add a line and under the line provide more description with markdown formatting and tables, etc.

