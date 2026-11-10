<!-- markdownlint-disable -->
# EAMxx (SCREAM) standalone testing

EAMxx builds and tests without CIME and without a supported HPC machine, using
the `copilot-testing` machine entry.

## Setup (once per environment)

From the repository root:

```bash
source components/eamxx/scripts/setup-copilot-env.sh
```

It installs system packages (apt-get, or spack as a fallback) and Python
packages, initializes the needed submodules, and exports `SCREAM_INPUT_ROOT`.
If a step fails, read the script and do that step by hand. The script is the
source of truth for the dependency list.

## Configure (once per build type)

```bash
cd components/eamxx
./scripts/test-all-eamxx -m copilot-testing -t dbg --config-only
```

| `-t` | Use for | Build directory under `ctest-build/copilot-testing/` |
|------|---------|------|
| `dbg` | default, debug with bounds checking | `full_debug` |
| `sp` | single-precision bugs | `full_sp_debug` |
| `opt` | release-mode or performance bugs | `release` |
| `fpe` | NaN/Inf, floating-point exceptions | `debug_nopack_fpe` |

If the input data download fails (for example `web.lcrc.anl.gov` is
unreachable), create empty placeholder files. Configure only checks that the
files exist. The data is needed only by the tests that read it.

## Build and test

Do not rerun `test-all-eamxx` after each edit. Work in the build directory:

```bash
cd components/eamxx/ctest-build/copilot-testing/full_debug
make -j$(nproc)
ctest -R <test_name_regex>      # or: ctest -j$(nproc), ctest --rerun-failed
```

If CMake or build system files changed, remove `CMake*` in the build directory
and run the configure command again.
