# EAMxx–MAM4xx all-aerosol surface-emissions regression

Date: 2026-10-07

Status: **PASS — 7 CTests, zero failures after restoration; required mutation
failed in exactly 8 `num_a1`/`num_a2` assertions.**

Validation mode: **STRICT**. This is test-only evidence for raw surface-flux
composition and the immediate lowest-layer constituent-flux application at a
fixed model state. It is not a concentration-level additivity claim.

## Revisions and isolation

| Item | Pinned value |
|---|---|
| Worktree | `/Users/burr114/CODEX_assistant/worktrees/e3sm-mam4xx-emissions-all-aerosol-test` |
| Branch | `susburrows/eamxx/mam4xx-emissions-all-aerosol-test` |
| E3SM / PR #8799 head | `8641fbdde8587b79c1979c38b39cc73ec90c4673` |
| MAM4xx gitlink / checkout | `b8a1c3cef3060de4817e1f07538e89494b930eef` |
| PR state checked | open, not draft, base `master`; metadata updated `2026-10-06T17:52:45Z` |
| Container image | `sha256:41ff04ef90bfbdd91889d9a5333fb3b7e3ad0e6b2184dd6c1f45aca9a0306af0` |
| Build/log volume | `mam-emissions-all-aerosol-20261007`, directory `/work/run` |

The older `e3sm-mam4xx-emissions-overwrite` worktree was read only. No GitHub
write, push, comment, PR update, or maintainer contact was performed.

External research: quick lookup — Perplexity was unavailable, so the official
GitHub PR metadata was queried directly. It confirmed that PR #8799 remains
open and that `8641fbd...` is still its current head. All implementation and
scientific checks were determined from the pinned repository source and tests.

## Expected-source matrix

The table below is explicit test data. Source class is not derived from the
algorithm under test. It was checked against:

- valid mode/species pairs in MAM4xx `aero_modes.hpp`;
- established slots in MAM4xx `aero_config.hpp`;
- online dust, sea-salt, and marine-organic targets in
  `aero_model_emissions.hpp`; and
- prescribed target/sector configuration in the EAMxx surface-emissions
  process.

| Mode | Tracer | Expected class |
|---:|---|---|
| 1 | `so4_a1` | PrescribedOnly |
| 1 | `pom_a1` | Inactive |
| 1 | `soa_a1` | Inactive |
| 1 | `bc_a1` | Inactive |
| 1 | `dst_a1` | OnlineOnly |
| 1 | `ncl_a1` | OnlineOnly |
| 1 | `mom_a1` | OnlineOnly |
| 1 | `num_a1` | Both |
| 2 | `so4_a2` | PrescribedOnly |
| 2 | `soa_a2` | Inactive |
| 2 | `ncl_a2` | OnlineOnly |
| 2 | `mom_a2` | OnlineOnly |
| 2 | `num_a2` | Both |
| 3 | `dst_a3` | OnlineOnly |
| 3 | `ncl_a3` | OnlineOnly |
| 3 | `so4_a3` | Inactive |
| 3 | `bc_a3` | Inactive |
| 3 | `pom_a3` | Inactive |
| 3 | `soa_a3` | Inactive |
| 3 | `mom_a3` | Inactive |
| 3 | `num_a3` | OnlineOnly |
| 4 | `pom_a4` | PrescribedOnly |
| 4 | `bc_a4` | PrescribedOnly |
| 4 | `mom_a4` | Inactive |
| 4 | `num_a4` | PrescribedOnly |

Counts are 2 `Both`, 8 `OnlineOnly`, 5 `PrescribedOnly`, and 10 `Inactive`.
The pinned source matches all source classes. One repository naming boundary is
intentional and now checked explicitly: surface-flux symbols use `ncl_a*`,
while established EAMxx interstitial/cloudborne field names use `nacl_a*`.

## Logic, units, and fixture

At fixed time and state, the test checks each aerosol tracer `c` at two
consecutive 1800 s steps:

```text
combined[c] ~= online_only[c] + prescribed_only[c]
```

The tolerance is
`10 * epsilon * max(abs(actual), abs(expected), 1)`. Inactive entries and all
gas/aerosol entries in the no-source configuration must be exactly zero.
Every expected online or prescribed component is required to be positive in at
least one represented `Real`; both components must overlap in at least one
column for `num_a1` and `num_a2` at each step.

The 218 controlled columns cycle through land, mixed ocean/land, and ocean;
use nonzero and varying wind, 300 K SST, nonzero marine-organic input data,
four nonzero dust bins with differing strengths, and dust/sea-salt scale
factors 1.5/0.6. Prescribed inputs are the current `c20260730` ne2np4 files;
the marine-organic file is `c20260807`. Zeroed and oracle NetCDF copies exist
only in per-test `/tmp/mam-emissions-*` directories and are removed by the
fixture.

For each of `so4_a1`, `so4_a2`, `bc_a4`, `pom_a4`, `num_a1`, `num_a2`, and
`num_a4`, an independent single-target fixture sets every stored sector to a
controlled constant. The oracle reads the stored float back and computes:

```text
expected = stored_sector_sum * source_scale * 1.65979e-23
           * mam4::gas_chemistry::adv_mass[target_slot - offset]
```

The pinned aerosol source scale is 1. Each oracle requires the intended target
slot to match and every other aerosol slot to remain exactly zero. The existing
number normalization, including the `1.0074` conversion entry, is preserved;
this test neither changes nor scientifically endorses it.

The real `mam4_constituent_fluxes` consumer is then checked for all 25 aerosol
tracers. Interstitial mass and number fields start from distinct positive
baselines. For the bottom layer:

```text
expected = baseline + surface_flux * dt * g / dp
```

All upper interstitial layers and every cloudborne layer must remain exactly at
baseline. Zero-flux/inactive bottom layers must also remain exactly unchanged.
Slots come from `mam4::AeroConfig`; field names come from `mam_coupling`.

## Build configuration and commands

The source and pinned dependency/input mounts were read only. The container was
limited to 2 CPUs, 6 GiB memory, 6 GiB swap, 512 PIDs, and one OpenMP thread.
The new test itself is registered at one rank because it keeps all 218 columns
local. Existing np1/np2 tests retain MPI comparison coverage.

Configuration recorded in `CMakeCache.txt`:

```text
CMAKE_BUILD_TYPE=Debug
CMAKE_CXX_COMPILER=/usr/bin/mpicxx (GNU 11.4.0)
CMAKE_C_COMPILER=/usr/bin/mpicc
CMAKE_Fortran_COMPILER=/usr/bin/mpifort
SCREAM_DOUBLE_PRECISION=TRUE
SCREAM_FPE=ON
SCREAM_FPMODEL=strict
SCREAM_NUM_VERTICAL_LEV=72
SCREAM_PACK_SIZE=16
SCREAM_TEST_MAX_TOTAL_THREADS=2
```

Exact commands executed inside the constrained container were:

```sh
cmake -S /src/components/eamxx -B /work/run/build \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_Fortran_COMPILER=mpifort \
  -DSCREAM_ENABLE_MAM=ON \
  -DSCREAM_NUM_VERTICAL_LEV=72 \
  -DSCREAM_FPE=ON \
  -DSCREAM_DYNAMICS_DYCORE=NONE \
  -DSCREAM_INPUT_ROOT=/inputdata \
  -DSCREAM_ENABLE_BASELINE_TESTS=OFF \
  -DSCREAM_TEST_MAX_TOTAL_THREADS=2 \
  -DEKAT_ENABLE_MPI=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_MODULE_PATH=/src/externals/scorpio/cmake \
  -DNetCDF_C_PATH=/system-prefix \
  -DNetCDF_Fortran_PATH=/system-prefix \
  -DPnetCDF_PATH=/system-prefix \
  -DFETCHCONTENT_SOURCE_DIR_SPDLOG=/src/externals/ekat/extern/spdlog \
  -DFETCHCONTENT_SOURCE_DIR_YAML_CPP=/src/externals/ekat/extern/yaml-cpp \
  -DFETCHCONTENT_SOURCE_DIR_CATCH2=/src/externals/ekat/extern/Catch2

cmake --build /work/run/build --parallel 2 --target \
  mam_surface_emissions_test \
  mam4_srf_online_emiss_standalone \
  mam4_srf_online_emiss_mam4_constituent_fluxes \
  cprnc

ctest --test-dir /work/run/build --output-on-failure -V \
  -R 'mam_surface_emissions_test|mam4_srf_online_emiss'
```

The host-side `docker run` supplied the limits above, mounted this worktree at
`/src` read only, the new volume at `/work`, and these older pinned-cache
subpaths read only: `deps/0..4`, `inputdata`, and `system-prefix` from
`mam-emissions-20260924`. Complete logs remain in
`mam-emissions-all-aerosol-20261007:/work/run/*.log`; they are not committed.

## Executed results

| Test | Result |
|---|---|
| `mam_surface_emissions_test` (np1) | PASS; 1 case, 1,762,948 assertions |
| `mam4_srf_online_emiss_standalone_np1` | PASS |
| `mam4_srf_online_emiss_standalone_np2` | PASS |
| standalone np2-vs-np1 CPRNC | PASS |
| `mam4_srf_online_emiss_mam4_constituent_fluxes_np1` | PASS |
| `mam4_srf_online_emiss_mam4_constituent_fluxes_np2` | PASS |
| coupled np2-vs-np1 CPRNC | PASS |

Final result: **7/7 CTests passed; 0 failed; 2.20 s wall time.** The OpenMPI
container emitted known nonfatal memory-binding warnings.

The first focused run failed one assertion because the initial fixture assumed
the `ncl_a*` flux spelling was also the prognostic-field spelling. The fix did
not change expectations: it made the established `ncl_a* -> nacl_a*` boundary
explicit and verified it. The next run passed. No production change was needed.

## Required overwrite mutation

Only a scratch overlay changed:

```cpp
constituent_fluxes_ispe_srf.update(sector_field_sum_, 1, 1);
```

to the old overwrite behavior:

```cpp
constituent_fluxes_ispe_srf.update(sector_field_sum_, 1, 0);
```

The completed generalized test then failed as required: 1,762,854 assertions
passed and **8 failed**, exclusively the `num_a1` and `num_a2` additivity and
overlap checks at both timesteps. CTest returned nonzero. A subsequently
created fixed overlay was byte-identical to the pinned source, all targets were
rebuilt, and the final seven-test suite passed without an overlay mount.

| Mutation artifact | SHA-256 |
|---|---|
| pinned/fixed production source | `504938f4c7810bfad75eca654344181e2c498e49c0803fe02847414e6e299410` |
| overwrite scratch overlay | `569d1a9fcad1cd34bd8be03b2f538dc572e437c889e2c3895c60eeb19dabcb88` |

No production mutation exists in the Git worktree. A separate mapping mutation
was not added; the explicit 25-row table, unique-slot check, field-name check,
seven single-target oracles, and consumer checks already fail on missing or
misplaced mappings without expanding the production mutation mechanism.

## Source and test hashes

| File | SHA-256 |
|---|---|
| `components/eamxx/src/physics/mam/tests/mam_surface_emissions_test.cpp` | `d25d41435c8c7faf49e60f5989f4d1a7889a04ce9d8e47e7fd7ae367faabd42c` |
| `components/eamxx/src/physics/mam/tests/CMakeLists.txt` | `de9218542550ac2a47d2dda286b97be66b8f920bab2a39c0da7f7ad1d61072f3` |
| surface-emissions production source | `504938f4c7810bfad75eca654344181e2c498e49c0803fe02847414e6e299410` |
| constituent-flux consumer functions | `3e53418942b9a1a6a124f274ffcf3d28077d2c6be2d74830e4b300f675a8b3df` |
| MAM4xx online-emissions source | `460189be88d2a83f2a1bfcd3c150a2e24f5f500e02e6ee394b0d7f34268761f1` |
| MAM4xx aerosol configuration | `f049a73cd0e3134b05c73579969619e0657ddaf640a3f36dc8c4246b52ffd55b` |

## Final review and limitations

The second review re-counted all 25 table entries and all four classes; checked
mass versus number units; traced every slot through `AeroConfig` and every
consumer field through `mam_coupling`; rechecked exact two-step reset/no-source
behavior; verified all lowest-layer, upper-layer, and cloudborne assertions;
confirmed the mutation failures; and inspected the final diff. `git diff
--check` passed, no changed line exceeds 100 columns, and only the two expected
test-source files plus this evidence README differ from the pinned PR head.
The required clang-format version 14 was not installed locally or in the pinned
container, so no claim of an executed clang-format check is made.

Remaining limits:

- The generalized 218-column fixture itself ran only at one rank. The existing
  standalone and coupled tests supplied one-/two-rank and CPRNC comparisons.
- Validation is local Linux ARM64, GNU 11.4, double-precision Debug/FPE only;
  no GPU, single-precision, full coupled model, historical baseline,
  performance, or long-run stability claim is made.
- The test preserves the current number normalization and process order.
- It establishes raw surface-flux composition and immediate lowest-layer
  application only, not downstream concentration additivity.

With those limits, the test-only result is ready to propose to Oscar for
PR #8799 review. It is not pushed or attached to the PR.
