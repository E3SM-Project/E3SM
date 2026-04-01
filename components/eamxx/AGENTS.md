# EAMxx (SCREAM) Agent Guide

This guide covers standalone testing of EAMxx modifications. It does NOT
require a supported HPC machine.

## One-Time Environment Setup

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

## One-Time Configure

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

## Incremental Build and Test

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

## Reconfiguring

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
