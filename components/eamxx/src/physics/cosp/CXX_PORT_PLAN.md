# EAMxx COSP: C++/Kokkos port plan

Status: planning, no port code yet.
Branch: `whannah/eamxx/cosp-cxx-port`, based on E3SM master `ac6eaa2`.
COSP submodule (`components/eam/src/physics/cosp2/external`): E3SM-Project/COSPv2.0 `afd43be` (branch `e3sm-master`).

## 1. Goal

Enabling COSP in EAMxx added about 25% to the runtime of an ne30pg2 case on 4 PM-GPU nodes. The
coupled production runs planned at ne256pg2 can't eat that. The fix is to make COSP run on the GPU
with the rest of EAMxx, so the ISCCP, MODIS and MISR simulators, plus the subcolumn generator and
optics that feed them, get ported from Fortran to C++/Kokkos. The bigger "who owns the COSP code"
conversation with the community is out of scope here and can happen later.

## 2. What EAMxx actually runs today

All of this was checked against the code on this branch.

### 2.1 The wrapper (`components/eamxx/src/physics/cosp/`)

| File | Lines | Role |
|---|---|---|
| `eamxx_cosp.cpp/.hpp` | 366 / 82 | `Cosp` AtmosphereProcess: field requests, call frequency, z_mid on device, host syncs, night masking |
| `cosp_functions.hpp` | 104 | `CospFunc::main`: copies LayoutRight host views into LayoutLeft host views, calls Fortran, copies back |
| `cosp_c2f.F90` | 1390 | iso_c_binding bridge: fills COSP derived types, runs `subsample_and_optics`, calls `cosp_simulator` |
| `CMakeLists.txt` | 60 | Builds the Fortran `cosp` library and `eamxx_cosp` |

Data flow on every COSP step (`eamxx_cosp.cpp:219-356`):

1. `sync_to_host` on 12 input fields, plus a device deep copy and host sync of the sunlit mask.
2. z_int and z_mid computed **on device** (`:257-282`), then z_mid is synced to host.
3. Host-side transpose of 12 (ncol, nlev) fields into LayoutLeft copies (`cosp_functions.hpp:59-78`).
4. Serial Fortran on one CPU core per rank: subcolumns, optics, ISCCP, MISR, MODIS.
5. Host-side transpose of the three 3D histograms back to LayoutRight, host loop for night fill values.
6. `sync_to_dev` on the 4 outputs, then device kernels to broadcast the sunlit mask.

### 2.2 Which Fortran gets compiled

This one matters for the "our COSP has already deviated" point. **EAMxx compiles the
submodule sources (`cosp2/external/...`), not EAM's modified copies in `cosp2/local/`.** EAM's
`Cosp.cmake` uses `local/cosp.F90`, `local/cosp_config.F90` and `local/modis_simulator.F90`, while
EAMxx's `CMakeLists.txt` uses the `external/src/` versions of all three. The local copies differ
from the submodule (132, 63 and 74 diff lines), so EAM and EAMxx are not running identical COSP.

The submodule itself (E3SM fork, `e3sm-master`) has also drifted from CFMIP upstream master
(`5eb05e5`, Apr 2026). Upstream has since cleaned up single precision literals in `icarus.F90`
and restructured `cosp.F90` and `cosp_stats.F90` (about 650 insertions and 2300 deletions over the
files in scope). Neither version has any OpenACC/OpenMP offload, so there is no ready-made GPU
path to pull in.

### 2.3 Hard-coded configuration (`cosp_c2f.F90:44-185`)

- Simulators on: ISCCP, MODIS, MISR. Everything else (CALIPSO, ATLID, ground lidar, PARASOL,
  CloudSat, RTTOV) is off.
- Outputs actually consumed by EAMxx: `isccp_totalcldarea`, `isccp_fq`,
  `modis_Optical_Thickness_vs_Cloud_Top_Pressure`, `misr_fq`. Nothing else leaves the wrapper.
- `overlap = 3` (max/random), `isccp_topheight = 1`, `isccp_topheight_direction = 2`,
  `emsfc_lw = 0.99` (set in `eamxx_cosp.cpp:305`).
- Convective cloud, convective condensate, and all precipitation inputs are zeroed.

### 2.4 Subcolumns and frequency (`cime_config/namelist_defaults_eamxx.xml:712-719`)

- `cosp_subcolumns` defaults to **10**, but is **1 for ne256np4, ne512np4 and ne1024np4** (the
  selector is a regex prefix match, so `ne256np4.pg2` gets 1). ne30pg2 gets 10.
- `cosp_frequency = 1`, `cosp_frequency_units = hours`.
- With 1 subcolumn, `subsample_and_optics` skips SCOPS and the RNG entirely and just sets
  `frac_out = 1` (`cosp_c2f.F90:1088-1233`). It also converts in-cloud optical depth and
  emissivity to grid-box values first (`cosp_c2f.F90:303-311`).

**So the ne256pg2 production configuration never touches the subcolumn generator or RNG.**

## 3. Corrections to the earlier summary

The earlier summary was a solid starting point, but several things in it don't match the code.

| Earlier claim | What the code shows |
|---|---|
| Overhead is dominated by host/device data movement | Not established, and the numbers make it doubtful (see 3.1). It needs measuring first. |
| Colleague's test used 1 subcolumn | Default at ne30pg2 is 10. Worth confirming what the colleague's case actually set. |
| `prec_scops.F90` is on the active path | It's called, but with zero precipitation rates, and its output only feeds radar/lidar quantities. It has no effect on the ISCCP/MODIS/MISR outputs, so it can be skipped. |
| `cosp_stats.F90` (509 lines) is needed, including vertical regrid | Only `hist2D` (about 30 lines, used by MODIS) is on the path. `cosp_change_vertical_grid` is only used by lidar/radar. |
| `cosp_optics.F90` (~500 lines) | Only `cosp_simulator_optics`, `modis_optics_partition` and `modis_optics` are used, about 150 lines. `lidar_optics` is the bulk and is off. |
| ~4,000 to 5,000 lines in scope | Closer to 2,000 to 2,500 (table in section 4). |
| RNG BFB on GPU is the biggest risk, may require abandoning BFB | Per-column streams already exist, and BFB looks very doable (3.2). Also irrelevant at ne256 (1 subcolumn). |
| Thread-safe reductions needed for binning | Histograms are per-column counts over at most `nsubcol` entries. With parallelism over columns, no atomics are needed. |

### 3.1 Is it really the transfers?

Rough size of the per-call traffic at ne30pg2 on 16 GPUs (21,600 columns, so about 1,350 per rank,
128 levels): 13 input fields of 1,350 x 128 doubles is about 18 MB down, and the outputs
(1 + 49 + 49 + 112 values per column) are about 2 MB back. At PCIe Gen4 rates that's on the order
of a millisecond. The host transposes are a similar size of memory traffic.

The Fortran itself, on the other hand, runs serially on one CPU core per rank, for 1,350 columns
x 10 subcolumns x 128 levels, including the MODIS size retrieval (a two-stream solve against 15
trial radii for every cloudy subcolumn) and the ICARUS brightness temperature iteration. Meanwhile
the other ~15 cores per rank sit idle. My guess is that the Fortran compute is the larger share,
not the copies. The port fixes both either way, but it changes which stopgaps are worth trying
(section 7, M0) and how much to expect from the 1-subcolumn ne256 config.

### 3.2 The RNG (`subsample_and_optics_example/subcol/rng/mo_rng.F90`)

- It's a Marsaglia congruential generator with one 32-bit integer of state per column, already
  seeded per column (`cosp_c2f.F90:1093-1099`). That maps straight onto one GPU thread per column.
- To match bit-for-bit, the C++ has to reproduce: signed 32-bit wraparound in
  `69069*seed+1234567` (do it in `uint32_t` and reinterpret, since signed overflow is UB in C++),
  Fortran `mod` with negative operands (same sign rule as C++ `%`), the single precision product
  `seed*0.931322574615479E-09`, and the `r = r+1; r = r-int(r)` fold for negative values.
- The sign-preserving `mod` means there's no simple jump-ahead, so draws for one column stay
  sequential. Cheap enough: `nsubcol*(nlev+1)` draws, about 1,300 for 10 x L128. Plan is a
  per-column kernel that fills `rand(ncol, ndraw)`, then the rest of SCOPS runs parallel over
  subcolumns.
- SCOPS quirk: at `ilev == 1` it uses `ran(npoints)`, the *last* column's draw, for every column
  (`scops.F90:136`). With `overlap = 3` and no convective cloud this value gets overwritten before
  it affects anything, but the draws still advance the streams, so the draw order has to be kept.
- The seed comes from min/max of rank-local surface pressure (`cosp_c2f.F90:1095-1098`), so COSP
  already isn't BFB across MPI decompositions (the standalone test's CMake comment says as much).
  A port can match this exactly with a `parallel_reduce`. Switching to a decomposition-independent
  seed (global column id plus step) would be nicer, but it's a non-BFB change for later.

## 4. Real scope of the port

| Piece | Source (under `cosp2/external/`) | Approx. lines in scope | Needed for 1 subcol (ne256) |
|---|---|---|---|
| Input prep, grid-box conversion | `cosp_c2f.F90` `cosp_c2f_run` | 90 | yes |
| Subcolumn generator | `subsample_and_optics_example/subcol/scops.F90` | 240 | no |
| RNG | `subsample_and_optics_example/subcol/rng/mo_rng.F90` | 60 | no |
| Subcolumn hydrometeor assignment | `cosp_c2f.F90` `subsample_and_optics` (cloud species only) | 80 | trivial |
| Optics | `cosp_optics.F90`: `cosp_simulator_optics`, `modis_optics_partition`, `modis_optics` | 150 | yes |
| ISCCP | `src/simulator/icarus/icarus.F90` + `cosp_isccp_interface.F90` | 650 | yes |
| MISR | `src/simulator/MISR_simulator/MISR_simulator.F90` + interface | 300 | yes |
| MODIS | `src/simulator/MODIS_simulator/modis_simulator.F90` (minus `*_old` fits) + interface | 800 | yes |
| Histogram | `src/cosp_stats.F90` `hist2D` | 30 | yes |
| Driver | `src/cosp.F90` `cosp_simulator`, ISCCP/MISR/MODIS branches only | 300 | yes |
| Bin tables | `src/cosp_config.F90` | 100 | yes |

Notes on dependencies: MODIS uses the ISCCP cloud-top pressure (`isccp_boxptop`), so ICARUS's
subcolumn part has to run before MODIS. The MODIS column routine computes about 20 outputs, of
which EAMxx uses one (the tau/CTP joint histogram), so the port can skip the rest for now.

## 5. BFB hazards worth knowing up front

1. **Single precision literals and default-real conversions.** The E3SM fork of COSP is full of
   them: `transmax(j) .le. 0.9999999` (a float, so 0.99999988...) and `tb .gt. 260.` in
   `icarus.F90:352-365`, `1.e-20` and `3.41817e-7` in `icarus.F90:268-269`,
   `real(count(...))/ncol` in `icarus.F90:574` and `MISR_simulator.F90:278`, and the `-999.`
   single precision sentinel in MISR. Each one needs an explicit `float` in the C++ to stay BFB. Fixing them is a
   separate, deliberately non-BFB step after the port lands.
2. **FMA on GPU.** nvcc contracts `a*b+c` by default. BFB against Fortran is a CPU-build test.
   GPU builds get compared against CPU C++ with a tolerance.
3. **Histogram outputs are discrete.** A roundoff difference in tau or cloud-top pressure can move
   a subcolumn to the neighbouring bin, so "matches to roundoff" isn't a meaningful check on the
   histograms. CPU comparisons should be exact, and GPU comparisons should count bin flips rather
   than look at max abs diff.
4. **Night handling.** MODIS packs sunlit columns (`pack` in `cosp.F90:722-727`) and writes
   `R_UNDEF` elsewhere, and EAMxx then overwrites night columns with its own fill value. On GPU
   this becomes a predicated loop over all columns with no compaction.
5. **Bin order flips.** ISCCP reverses the pressure axis (`cosp.F90:962`) and MODIS does the same
   (`cosp.F90:1399-1400`). Easy to get wrong when the outputs move to LayoutRight.
6. **Vertical order.** COSP gets levels top to bottom here, same as EAMxx (`cosp_c2f.F90:334-344`
   passes them unflipped), so no flips are needed on input.

## 6. Target design

- **Location and structure.** `components/eamxx/src/physics/cosp/` gets a
  `cosp_functions.hpp` with a templated `CospFunctions<Real, Device>` struct (same pattern as
  P3/SHOC), with kernels in `impl/`: `cosp_rng_impl.hpp`, `cosp_scops_impl.hpp`,
  `cosp_optics_impl.hpp`, `cosp_isccp_impl.hpp`, `cosp_misr_impl.hpp`, `cosp_modis_impl.hpp`,
  `cosp_hist_impl.hpp`. Bin edges and centres from `cosp_config.F90` become `constexpr` arrays,
  and `cosp_c2f_get_bins` goes away.
- **Data.** Use the field views directly (LayoutRight, (ncol, nlev)). Subcolumn work arrays
  (ncol, nsubcol, nlev) are allocated once in `initialize_impl`, later from the ATM buffer.
  `frac_out` becomes an integer view. No derived types: plain views passed to each kernel.
- **Parallelism.** `TeamPolicy` with one team per column and threads over subcolumns, with vector
  ranges over levels where the work is a reduction (column optical depth, brightness temperature
  sums). Histograms are filled by the column's team with no atomics.
- **Only compute what's consumed.** Four outputs today. The structure should make it easy to add
  more MODIS/ISCCP diagnostics later without redesigning.
- **Transition switch.** A runtime parameter, say `cosp_impl: fortran | cxx`, so both paths can run
  in the same build for A/B comparisons in real cases. It goes away once the C++ path is
  validated.
- **Incremental, but device-resident from day one.** Instead of swapping kernels one at a time
  behind the current host interface (which means writing throwaway host/device plumbing for every
  intermediate state), each kernel is unit tested against the Fortran on its own, and the driver
  switches over in one step once the whole 1-subcolumn path is ported. The pipeline is linear and
  only ~2,000 lines, so this is less risky than it sounds and avoids the glue code.

## 7. Milestones

**M0. Measure first (small, no port code).** Add timers in `Cosp::run_impl` around the host syncs,
the transposes in `CospFunc::main`, `cosp_c2f_run`, and the output syncs. Run ne30pg2 on PM-GPU
with `cosp_subcolumns` = 10 and = 1. That tells us how the 25% splits, and what to expect at
ne256. Stopgaps to consider meanwhile, depending on the result: call less often (3-hourly), use 1
subcolumn at ne30 too, or thread the Fortran over column chunks with OpenMP to use the idle cores.

**M1. Test harness and leaf kernels.** RNG, `cosp_simulator_optics`, `hist2D`, the MODIS optics
partition, and `modis_optics`. Each one gets a test that calls Fortran and C++ on the same random
inputs (CPU build) and requires exact equality. The Fortran side is reached via small
`bind(C)` test shims for the leaf routines.

**M2. The 1-subcolumn path (everything ne256 needs).** ICARUS subcolumn and column, MISR subcolumn
and column, MODIS subcolumn and column, then the driver. Flip `run_impl` to the C++ path under the
`cosp_impl` switch, with no host syncs. Compare against Fortran in `cosp_standalone` with
`cosp_subcolumns: 1`.

**M3. SCOPS and the 10-subcolumn path.** RNG pre-draw kernel, SCOPS, subcolumn assignment, and
`cosp_standalone` comparisons with 10 subcolumns.

**M4. Validation and performance.** GPU vs CPU C++ comparison (bin-flip counts), short ne30pg2
run on PM-GPU, timing with COSP on vs off, and an ne256 timing if allocation allows. Make `cxx`
the default.

**M5. Cleanup.** Drop the Fortran path from the EAMxx build (EAM keeps using it), convert the BFB
tests to the baseline-file style P3/SHOC now use (`UnitBase` and `m_baseline_action` in
`src/share/physics/physics_test_data.hpp`), and refresh the `cosp_standalone` baselines.

Optional follow-ups, all non-BFB: fix the single precision literals, decomposition-independent
seeding, and pulling in upstream CFMIP changes.

## 8. Testing notes

- Standalone EAMxx builds and runs on a generic Linux box per `components/eamxx/AGENTS.md`
  (`test-all-eamxx -m copilot-testing`), CPU only. That's enough for all the BFB-vs-Fortran work
  in M1 to M3. GPU checks need Perlmutter.
- `tests/single-process/cosp/` already has a `cosp_standalone` AD test (ne4, 72 levels, 2/5/48
  steps) with a baseline compare. It's the natural place for Fortran vs C++ whole-pipeline
  comparisons.

## 9. Open questions

1. What did the colleague's ne30pg2 case set for `cosp_subcolumns` (default is 10)?
2. Is strict BFB against the current Fortran (on CPU) the bar for the port? It looks achievable
   and keeps validation simple, so that's the recommendation.
3. Port the E3SM fork as-is (BFB, single precision quirks included), then fix the quirks as a
   separate step? Recommended, versus porting from upstream CFMIP directly.
4. OK with a temporary `cosp_impl` runtime switch during the transition?
