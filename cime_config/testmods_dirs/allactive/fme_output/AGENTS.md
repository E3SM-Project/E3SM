# FME Online Output Processing

Testmod that adds online horizontal remapping, vertical coarsening, derived
field computation, and vertical reduction to E3SM for FME (Full Model
Emulation) training data generation (ACE/Samudra). Eliminates offline
postprocessing. Overhead: +1.6% wallclock.

## Source Files

Shared infrastructure (pure computation, no I/O or MPI):
- `share/util/shr_horiz_remap_mod.F90` -- CRS SpMV remap + apply_masked
- `share/util/shr_vcoarsen_mod.F90` -- vertical coarsening (pressure/depth)
- `share/util/shr_derived_mod.F90` -- expression parser (e.g. `STW=Q+CLDICE+CLDLIQ`)

EAM wrappers (in `components/eam/src/control/`):
- `horiz_remap_mod.F90` -- per-tape remap via cam_history (transparent)
- `eam_vcoarsen.F90` -- 8-layer pressure-bounded averaging + column integration
- `eam_derived.F90` -- chained expressions, tendencies, phys/dyn split

MPAS-Ocean (in `components/mpas-ocean/src/`):
- `analysis_members/mpas_ocn_fme_depth_coarsening.F` -- 25-level depth coarsening
- `analysis_members/mpas_ocn_fme_derived_fields.F` -- SST, SSS, heat flux, etc.
- `analysis_members/mpas_ocn_fme_vertical_reduce.F` -- OHC, FWC, KE
- `shared/mpas_ocn_fme_horiz_remap.F` -- direct PIO writes to lat-lon files

MPAS-Sea-Ice (in `components/mpas-seaice/src/`):
- `analysis_members/mpas_seaice_fme_derived_fields.F` -- 9 derived fields
  (iceAreaTotal, iceVolumeTotal, snowVolumeTotal, iceThicknessMean,
  surfaceTemperatureMean, airStressZonal, airStressMeridional,
  uVelocityGeo, vVelocityGeo)
- `shared/mpas_seaice_fme_horiz_remap.F` -- singleton horiz remap module

Testmod (this directory):
- `shell_commands` -- case setup, map file detection, namelist generation
- `verify_eam.py` -- EAM verification dashboard (HTML + figures)
- `verify_mpas.py` -- MPAS-O/SI verification dashboard
- `../fme_legacy_output/shell_commands` -- raw output for offline pipeline comparison

## Build and Test

```bash
cd cime/scripts
./create_test SMS_Ld2_P512.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850NS.pm-cpu_gnu \
    --testmod allactive/fme_output

# EAM verification (requires xarray, matplotlib, cartopy, netCDF4):
micromamba run -n xgns python \
    cime_config/testmods_dirs/allactive/fme_output/verify_eam.py \
    --rundir $RUNDIR --outdir /path/to/output_figs

# MPAS verification:
micromamba run -n xgns python \
    cime_config/testmods_dirs/allactive/fme_output/verify_mpas.py \
    --rundir $RUNDIR --outdir /path/to/output_figs
```

Map files default to `/pscratch/sd/m/mahf708/fme_maps/` (NERSC); override
with the `FME_MAPS_DIR` environment variable on other systems.
- `map_ne30pg2_to_gaussian_180by360_shifted.nc` (EAM, n_a=21600)
- `map_IcoswISC30E3r5_to_gaussian_180by360_shifted.nc` (MPAS, n_a=465044)

## Architecture Gotchas

These are hard-won lessons. Read before modifying FME code.

1. **PIO type mismatch garbles output silently.** `pio_initdecomp(PIO_REAL)` +
   `pio_write_darray(real(r8))` produces garbage with no error. Always convert
   to `real(4)` before writing when decomposition uses `PIO_REAL`.

2. **Fortran sequence association stride bug.** Passing `array(pcols,nlev)` to a
   dummy declared as `(ncol,nlev)` with `ncol < pcols` reads wrong memory at
   levels k>1. Always pass explicit sections: `array(1:ncol,:)`.

3. **apply_masked packs validity mask as level numlev+1.** One `MPI_Alltoallv` +
   one SpMV handles both field and mask. Normalizes by `valid_frac`
   (sum of weight * mask), not `frac_b`. Outputs `SHR_FILL_VALUE` (1e20)
   where `valid_frac < epsilon`.

4. **MPAS framework calls AM init twice.** `register_var` must skip duplicates
   or PIO throws "name in use" errors.

5. **Fill value flow.** All FME code uses `SHR_FILL_VALUE` (1e20) consistently.
   AMs set fill at source -> remap caller zeros field + sets mask=0 ->
   `apply_masked` normalizes by valid_frac -> output gets `SHR_FILL_VALUE`.
   Detection threshold is `SHR_FILL_VALUE * 0.1` (1e19).

6. **MPAS remapped output bypasses the stream manager.** The stream manager only
   supports `nCells`/`nEdges`/`nVertices` decomposed dims (hardcoded in
   `mpas_stream_manager.F:is_decomposed_dim`). Remapped `(lon,lat)` output uses
   direct PIO writes with date-based file rotation via `check_rotate` +
   `register_var` replay.

7. **Map file uses global cell IDs** (n_a = N_total), not ocean-only. The
   `gcol_to_ocean_col` identity branch is taken.

8. **shr_vcoarsen_mod requires monotonically INCREASING coordinates.** The overlap
   formula `max(0, min(c_hi,b_hi) - max(c_lo,b_lo))` assumes `c_lo < c_hi` and
   `b_lo < b_hi`. This covers pressure (Pa) and depth (m).

9. **Depth bound parsing validates at init.** `iostat` error handling on each
   token, then non-negativity and strict monotonicity checks.
   `nCoarsenLevels > nFmeDepthLevels` is a fatal error.

10. **send_local_cell==0 is fatal.** An incompatible map file (wrong mesh) causes
    unmapped source cells. This is caught at init with `MPAS_LOG_CRIT`, not
    silently masked.

11. **MPAS time averaging accumulates in remap space.** Accumulators are
    module-scope arrays of size `n_b_local` (target grid points on this rank).
    Each compute call remaps fields and adds to accumulators. When the output
    interval elapses, values are normalized by `nAccum` and written. Native
    stream output remains instantaneous (last snapshot). The trigger advances
    `avg_last_output_time` by exactly one `avg_output_interval` (not by
    `currTime`) so the schedule does not drift across many windows. On
    restart the baseline is `MPAS_NOW`, not `MPAS_START_TIME`, so the
    trigger doesn't fire every step until catching up; the partial pre-
    restart sample is dropped (the accumulators are remap-space and not
    in the Registry, so MPAS restart streams cannot persist them).
    `write_avg_field` guards against `nAccum == 0` and writes
    `SHR_FILL_VALUE` instead of dividing by zero.

12. **Remapped output files include xtime.** The `xtime(StrLen, Time)` character
    variable records the MPAS date-time string for each output record. CF
    attributes include `axis='X'`/`'Y'` on lon/lat and `calendar='noleap'`
    as a global attribute. PIO `pio_put_var` for xtime requires a 1-element
    character array (`character(len=64) :: buf(1)`), not a scalar — PIO's
    `put_vara_1d_text` interface demands `character(len=*) :: val(:)`.

13. **MPAS coupler fluxes are passed through raw (no ice-fraction recovery).**
    The forcing-pool fields (`shortWaveHeatFlux`, `longWaveHeatFluxDown`,
    `latentHeatFlux`, `sensibleHeatFlux`, `windStressZonal`,
    `windStressMeridional`) are delivered by the coupler as
    `flux_atm * (1 - iceFraction)`. `mask_and_remap` does *land* masking
    only -- it does NOT divide by `(1 - iceFraction)`. This matches
    `timeSeriesStatsCustom` semantics (raw ocean-side fluxes). If
    SamudrACE training needs the atmospheric-side flux, take it from
    EAM (TAUX, TAUY, FLDS, FSDS, ...); the MPAS-O versions are the
    ocean-side ice-weighted values. An earlier revision did the
    recovery via `iceFraction` and cell-fill for `iceFraction>=0.99`;
    that path was reverted to keep parity with the legacy pipeline.
    SST, SSS, and SSH are state variables and were never corrected.

14. **EAM averaging convention.** State variables (T, U, V, STW, PS, TS)
    are instantaneous snapshots (`avgflag_pertape='I'`). Fluxes and
    tendencies (FSDS, FLDS, FLUT, FSUTOA, FLUS, FSUS, LHFLX, SHFLX,
    TAUX, TAUY, PRECT, SOLIN) use per-field `:A` suffix for time
    averaging over the output interval. EAM shortwave radiation fields
    (FSUTOA, FSUS) show sharp gradients at ice edges — this is real
    surface albedo contrast (ice ~55% vs ocean ~3%), not an artifact.

15. **MPAS instantaneous output remaps and writes per variable, not in two
    passes.** The earlier pattern was `loop { remap_to_fld_remap }` then
    `loop { write_inst(name) }` which silently wrote the *last* remapped
    field under every variable name. Each AM now interleaves remap+write
    in instantaneous mode (`inst_remap_and_write` / `remap_and_write_inst`).
    Time-averaged mode is unaffected because it indexes accumulators
    correctly. The production testmod sets `time_averaging=.true.`.

16. **eam_derived tendency flags are per-chunk arrays.** `tend_initialized`
    and `phys_snap_valid` are dimensioned `(begchunk:endchunk)`, not
    scalar; otherwise the OMP-PARALLEL chunk loop in `physpkg.F90:1429`
    races and produces chunk-dependent tendencies on the first timestep
    after init or restart. Each thread now reads/writes only its own
    chunk's slot.

17. **shr_horiz_remap coverage threshold is unified.** Both `apply` (frac_b
    normalization) and `apply_masked` (valid_frac normalization) use the
    same `SHR_COVERAGE_EPS=1e-6` cutoff so the same coastline cells fill
    regardless of which path the caller uses.

18. **shr_derived divide-by-zero produces SHR_FILL_VALUE, not 0.** A `RATIO=A/B`
    expression with `B==0` writes `1e20` (matching the remap fill
    convention) so silent zeros don't mask buggy ratio definitions. The
    only configured expression is `STW=Q+CLDICE+CLDLIQ+RAINQM`, which
    has no division.

## Runtime Configuration

Both `fme_output` and `fme_legacy_output` testmods accept environment variables:

```bash
FME_EAM_OUTPUT_HOURS=6          # EAM output frequency in hours (default: 6)
FME_EAM_AVGFLAG=I               # 'I' = instantaneous (default); fluxes use per-field :A
FME_EAM_MFILT=4                 # samples per output file (default: 4)
FME_MPAS_INTERVAL=00-00-01_00:00:00  # MPAS output/averaging interval (default: daily)
FME_MAPS_DIR=/path/to/fme_maps  # SCRIP map directory (default: NERSC mahf708 scratch)
```

Example: `FME_EAM_OUTPUT_HOURS=24 FME_MPAS_INTERVAL=00-00-05_00:00:00 ./create_test ...`

MPAS FME AMs support time averaging via namelist:
- `config_AM_fme*_time_averaging = .true.` enables accumulation
- `config_AM_fme*_compute_interval = 'dt'` samples every timestep
- `config_AM_fme*_stream_output_interval` sets the averaging window
- The `fme_output` testmod enables averaging by default
- Three ocean AMs enabled: fmeDepthCoarsening, fmeDerivedFields, fmeVerticalReduce
- One sea-ice AM enabled: fmeSeaiceDerivedFields

Verification scripts support cross-comparison:
```bash
# EAM verification
python verify_eam.py --rundir $RUNDIR --outdir /path/to/figs

# MPAS verification with legacy cross-comparison
python verify_mpas.py --rundir $RUNDIR --outdir /path/to/figs \
    --legacy-rundir $LEGACY_RUNDIR
```

## Remaining Work

### Near-term
- Add CI test variant (SMS_Ld2, ne4pg2_oQU480 for fast builds)
- Add ERS restart test and PEM MPI reproducibility test
- Parameterize map file paths via env variables in shell_commands

### Code quality
- 3D PIO decomposition for depth-coarsened fields (single 3D var instead of
  per-level 2D vars like `temperatureCoarsened_0..24`)
- Depth coordinate variable in remapped 3D output
- SSH masking under ice shelf cavities (cavities have maxLevelCell > 0 and
  can have SSH < -1000m; not caught by iceFraction correction)
- Add lonCell/latCell to MPAS-O FME Registry stream definitions for native
  grid diagnostic plots

### Extensions
- Numeric `time` coordinate + `time_bnds` for CF-compliant time axis
- Online spherical harmonics roundtrip (replaces ACE offline SH smoothing)
- 2D/3D wetmask computation
- ELM / MOSART FME integration (deferred -- needs reimplementation)
- Performance: persistent work arrays, O(n) send_local_cell lookup, OpenMP
