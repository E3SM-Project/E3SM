# FME Online Output Processing

Testmod that adds online horizontal remapping, vertical coarsening, derived
field computation, and vertical reduction to E3SM for FME (Full Model
Emulation) training data generation (ACE/Samudra). Eliminates offline
postprocessing. Overhead: +1.6% wallclock.

> **Next major workstream:** see `## PLAN: Coupler-Native FME Output —
> As-Exchanged Merged Forcing (driver-mct)` at the END of this file — a
> self-contained, design-approved brief (2026-06-20, not yet implemented)
> for emitting the merged post-coupler forcing (`x2o`/`x2i`/`xao`) that no
> component tape can provide.

## Source Files

Shared infrastructure (pure computation, no I/O or MPI):
- `share/util/shr_horiz_remap_mod.F90` -- CRS SpMV remap + apply_masked
- `share/util/shr_vcoarsen_mod.F90` -- vertical coarsening (pressure/depth)
- `share/util/shr_derived_mod.F90` -- expression parser (e.g. `STW=Q+CLDICE+CLDLIQ`)

EAM wrappers (in `components/eam/src/control/`):
- `horiz_remap_mod.F90` -- per-tape remap via cam_history (transparent)
- `eam_vcoarsen.F90` -- 8-layer pressure-bounded averaging + column integration
- `eam_derived.F90` -- chained expressions, tendencies, phys/dyn split
- `cam_history.F90` (modified) -- adds `hist_file_storage_type` namelist
  for calendar-based file rotation (`one_month` / `one_year`),
  ported from EAMxx's StorageType pattern

MPAS-Ocean (in `components/mpas-ocean/src/`):
- `analysis_members/mpas_ocn_fme_depth_coarsening.F` -- 25-level depth coarsening
- `analysis_members/mpas_ocn_fme_derived_fields.F` -- SST, SSS, heat flux, etc.
- `analysis_members/mpas_ocn_fme_vertical_reduce.F` -- OHC, FWC, KE
- `shared/mpas_ocn_fme_horiz_remap.F` -- direct PIO writes to lat-lon files

MPAS-Sea-Ice (in `components/mpas-seaice/src/`):
- `analysis_members/mpas_seaice_fme_derived_fields.F` -- 9 derived fields
  (iceAreaTotal, iceVolumeTotal, snowVolumeTotal, iceThicknessMean,
  surfaceTemperatureMean, airStressZonal, airStressMeridional,
  uVelocityGeoCell, vVelocityGeoCell — renamed from uVelocityGeo /
  vVelocityGeo to avoid a name collision with the vertex-located
  geographicalVectors AM in maint-3.0)
- `shared/mpas_seaice_fme_horiz_remap.F` -- singleton horiz remap module

Testmod (this directory):
- `shell_commands` -- case setup, map file detection, namelist generation
- `verify_eam.py` -- EAM verification dashboard (HTML + figures)
- `verify_mpas.py` -- MPAS-O/SI verification dashboard
- `../fme_legacy_output/shell_commands` -- raw output for offline pipeline comparison

Production run script (at repo root):
- `run_e3sm.fme.sh` -- create_newcase wrapper for SamudrACE production.
  Branches from a spun-up v3 piControl, applies the FME testmod via
  `--user-mods-dirs`, sets 5-yr restart cadence aligned with the
  monthly file rotation. PLACEHOLDER refcase fields must be filled
  in before submit; an early-fail guard refuses to proceed otherwise.

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

Map files default to `$DIN_LOC_ROOT/fme/` (the standard E3SM inputdata tree --
`shell_commands` sets `FME_MAPS_DIR=$(xmlquery --value DIN_LOC_ROOT)/fme`);
override with the `FME_MAPS_DIR` environment variable to point elsewhere.
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
   Detection threshold is `SHR_FILL_VALUE * 0.1` (1e19). Note: this only
   handles the *spatial* dimension. Time-averaging accumulators need
   per-cell valid-sample counts to avoid summing fill into the mean —
   see gotcha #36.

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

11. **MPAS time averaging accumulates in remap space; restart-persisted
    via sidecar files.** Accumulators are module-scope arrays of size
    `n_b_local` (target grid points on this rank). Each compute call
    remaps fields and adds to accumulators. When the output interval
    elapses, values are normalized by `nAccum` and written. The trigger
    advances `avg_last_output_time` by exactly one `avg_output_interval`
    (not by `currTime`) so the schedule does not drift across many
    windows.

    **Restart persistence** (added 2026-05-01 to make ERS pass for FME
    streams): each AM writes a `<amName>.fme_accum_restart.nc` sidecar
    in the run dir whenever the coupler-driven restart alarm fires
    (`seq_timemgr_RestartAlarmIsOn(EClock)` branch in `ocn_comp_mct.F`
    and `ice_comp_mct.F`). The sidecar holds `remap_accum`,
    `remap_accum_count` (added 2026-05-04 — gotcha #36),
    `remap_latest_inst` (depth-coarsening only), `nAccum`, and
    `avg_last_output_time` as a global string attribute. PIO writes use
    a dedicated `iodesc_2d_dbl` (PIO_DOUBLE) with `n_avg_fields` as the
    record (UNLIMITED) dim — the variable is 3D `(lon, lat, n_avg_fields)`
    so `pio_setframe` can walk slots; defining it as 2D segfaulted at
    setframe call. On AM init when `config_do_restart=.true.`, the AM
    reads the sidecar and overrides the default `MPAS_NOW` baseline
    with the persisted timestamp. `write_avg_field` guards against
    `nAccum == 0` and writes `SHR_FILL_VALUE` instead of dividing by
    zero.

    **`compute_on_startup=.false.` is required for warm-restart BFB.**
    With the default `=true`, the AM compute runs once at MPAS_NOW
    before timestepping — a sample at start_time on cold runs, at
    restart_time on warm runs. The averaging window straddling the
    restart point sees one extra sample on warm legs that the
    continuous-run baseline doesn't, producing a ~1e-4 RMS diff in
    that record. Disabling `compute_on_startup` for all four FME AMs
    (in `shell_commands` and the `user_nl_*` for testmod cases) makes
    cold and warm runs accumulate identically. Cost: the very first
    averaged record of a cold run loses one sample (1 out of ~N per
    window — invisible in production at daily/30-min cadence).
    `write_on_startup` was already `.false.`.

12. **OBSOLETE 2026-05-04: xtime is no longer written.** Remapped output
    files used to include `xtime(StrLen, Time)` (legacy MPAS character
    timestamp), which required a `(1)`-element character buffer trick to
    work around the PIO `put_vara_1d_text` interface. With the CF-time
    rewrite (gotcha #39), xtime, the StrLen dim, and the PIO_CHAR
    dependency are all gone. Kept the entry as a breadcrumb for anyone
    reading older branches.

13. **MPAS coupler fluxes are passed through raw (no ice-fraction recovery).**
    The forcing-pool fields (`shortWaveHeatFlux`, `longWaveHeatFluxDown`,
    `latentHeatFlux`, `sensibleHeatFlux`, `windStressZonal`,
    `windStressMeridional`) are delivered by the coupler as
    `flux_atm * (1 - iceFraction)`. `mask_and_remap` does *land* masking
    only -- it does NOT divide by `(1 - iceFraction)`. This matches
    `timeSeriesStatsCustom` semantics (raw ocean-side fluxes). If
    SamudrACE training needs the atmospheric-side flux, take it from
    EAM (TAUX, TAUY, FLDS, FSDS, ...); the MPAS-O versions are the
    ocean-side ice-weighted values. SST, SSS, and SSH are state
    variables and were never corrected.

    *History (revisit if needed):* an earlier revision of this branch
    did recover the atmospheric-side flux by dividing by
    `(1 - iceFraction)` in both `mask_and_remap` and the
    `surfaceHeatFluxTotal` cell loop (MPAS-O and MPAS-SI), with cells
    at `iceFraction >= 0.99` falling back to `fillValue` to avoid
    amplifying noise. That path was deliberately reverted on
    2026-04-29 to keep parity with `timeSeriesStatsCustom` and the
    legacy offline pipeline (decision applied while adopting fixes
    from `claude/review-fme-reconciled-E56bM`). To revisit: see the
    `mask_and_remap` body and the `surfaceHeatFluxTotal` cell loop in
    `mpas_ocn_fme_derived_fields.F` and `mpas_seaice_fme_derived_fields.F`
    in branch `claude/fme-reconciled.bak-pre-rebase-maint30` (commit
    `707a185 "Recover atmospheric fluxes from coupler ice-weighting"`).

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

19. **Verify dashboard uses per-layer bathymetry masks for depth-coarsened
    fields.** A naive surface-SST land mask flags every deep
    layer-below-seafloor cell as a false-positive "ocean NaN". The dashboard
    now bootstraps both a surface mask (from SST) and a per-layer mask (from
    `layerThicknessCoarsened_{k}`) in a dedicated first pass before the main
    check loop, then chooses the appropriate mask per variable. The first
    pass also fixes a quiet ordering bug where flux variables sorted before
    `sst` lexicographically were never mask-checked.

20. **Verify cross-compare uses abs_diff for sign-changing / near-zero
    fields.** Global-mean TAUX/TAUY oscillate around zero, so any small
    discrepancy blows up `rel_diff` (e.g. 0.003 N/m^2 of remap noise becomes
    37%). Same for ICEFRAC, FSUS, FLUS over partial-ice cells when the FME
    averaged tape is compared against the legacy daily tape. EAM
    cross-verify defines `GMEAN_ABS_TOL = {TAUX, TAUY, ICEFRAC, FSUS, FLUS}`
    and passes if `abs_diff < tol` even when `rel_diff` is large. MPAS-O
    `check_remap_conservation` uses `ssh_abs_tol = 3 m` because SSH global
    mean depends on the volume-conservation reference choice
    (pressureAdjustedSSH vs raw) and can offset by O(1 m) even when the
    spatial pattern matches.

21. **`surfaceTemperatureMean` is Celsius, not Kelvin.** The MPAS-SI
    `mpas_seaice_fme_derived_fields.F` averages Icepack `surfaceTemperature`
    which is in Celsius. The Registry XML now declares `units="degC"` to
    match. All verify plot vmin/vmax use `(-45, 5)` degC; RANGE_CHECKS
    uses `(-50, 5)`. Do NOT change the Fortran to convert to Kelvin
    without auditing all four sites in `verify_mpas.py`.

22. **Verify fill threshold is `1e18`, not `1e10`.** `oceanHeatContent`
    reaches `~6e11 J/m^2` in deep equatorial Pacific
    (`T_C * h * rho_sw * cp_sw = 30*5000*1025*4000`), so the older
    `fill_thresh = 1e10` clipped genuine physical values and reported
    them as "ocean NaN/fill" (~31000 cells per timestep). All 22 sites
    in `verify_mpas.py` now use `1e18` -- well above any geophysical
    magnitude, well below `SHR_FILL_VALUE = 1e20`. Safe because
    apply_masked produces no intermediate fill artifacts (see #3, #5,
    #17). EAM `verify_eam.py` keeps `1e10` since no EAM field
    approaches that magnitude.

23. **`uVelocityGeoCell`/`vVelocityGeoCell` are FME-AM only.** The
    cell-centered geographic-frame ice velocities live in
    `Registry_seaice_fme_derived_fields.xml` and are computed inside
    `mpas_seaice_fme_derived_fields.F`. They appear in the remapped FME
    output but NOT in the native `timeSeriesStatsCustom` tape. Verify
    has an `ICE_FME_ONLY` set that exempts them from the native
    presence check while still requiring them in the remapped check.

24. **`DEPTH_BOUNDS` matches the production namelist.** The 25-level
    Registry max for `nFmeDepthLevels` is the *allocated* dimension;
    the *active* layer count comes from
    `config_AM_fmeDepthCoarsening_depth_bounds` (production: 19
    layers, 20 boundaries -- matching `verify_mpas.py:DEPTH_BOUNDS`).
    The offline depth-coarsening cross-check now detects the active
    layer count by counting non-fill slabs in `temperatureCoarsened`
    rather than using the dimension. If the namelist changes, update
    both `shell_commands` and `DEPTH_BOUNDS` in lockstep.

25. **Verify ice-presence mask is bootstrapped from `iceAreaTotal`.**
    Sea-ice derived fields (`iceVolumeTotal`, `iceThicknessMean`,
    `surfaceTemperatureMean`, `airStress*`, `u/vVelocityGeoCell`) are
    validly fill where there is no ice -- *both* on land *and* on
    open-ocean cells without ice. The surface SST land mask says
    "valid in all ocean", which false-positives ~31000 open-ocean
    cells per file. The ice mask is `iceAreaTotal is fill OR exactly
    zero`. `iceAreaTotal` itself is fill on land and 0 on open ocean
    (so the ice mask correctly flags it as "valid only with ice
    presence" while CHECK 4 exempts it from "valid-on-land" anyway).

26. **EAM STW vcoarsen linearity uses cos-lat weighting for FME mean.**
    FME output is on a regular lat-lon grid; legacy native (ne30pg2)
    is quasi-equal-area. For fields with a strong meridional gradient
    (STW = Q+CLDLIQ+CLDICE+RAINQM, which drops by orders of magnitude
    from equator to pole), an unweighted FME mean over-represents
    polar cells and biases low by 10-25% vs the area-weighted legacy
    mean -- producing a spurious "STW linearity FAIL" even though
    vcoarsen IS exactly linear. The Test 3 fallback now mirrors Test
    2's cos-lat weighting and uses the same 2% cross-grid threshold.

27. **EAM `hist_file_storage_type` enables calendar-based file rotation.**
    Ported from EAMxx's StorageType pattern (`scream_io_file_specs.hpp:32-71`).
    Per-tape namelist option in `cam_history_nl` with three values:
    - `'num_snapshots'` (default): legacy `mfilt`-based rotation
    - `'one_month'`: rotate at calendar-month boundary, file `%y-%m.nc`
    - `'one_year'`: rotate at calendar-year boundary, file `%y.nc`

    Implementation: `wshist` checks the new record's month/year against
    `hist_storage_curr_idx(t)` *before* the "Starting a new volume"
    branch; on mismatch, closes the file and resets `nfils=0` so the
    existing new-file path opens a fresh one with the new date.
    Lazy-init: a `-1` curr_idx (e.g. fresh run) inherits from the
    current record. After restart, `read_restart_history` parses
    `.YYYY-MM.nc` / `.YYYY.nc` from the open file's name to seed
    `hist_storage_curr_idx`, so calendar-boundary restarts rotate
    correctly. `mfilt` becomes a safety upper bound when calendar
    rotation is in charge -- production sets 1500 (any month has
    at most 124 records at 6-h cadence). FME testmod uses
    `'one_month'` by default via `FME_EAM_STORAGE=one_month`.

28. **MPAS FME streams use `.YYYY-MM.nc` for monthly aggregation.**
    The runtime `streams.ocean` and `streams.seaice` files are
    generated by `cime_config/buildnml` (Python) at case_setup time,
    NOT from the Registry XMLs. Edits to the Registry stream
    declarations are documentation-only -- they don't change
    runtime behavior. The four FME stream blocks live at:
    - `components/mpas-ocean/cime_config/buildnml:1251-1319`
      (fmeDepthCoarsening, fmeDerivedFields, fmeVerticalReduce)
    - `components/mpas-seaice/cime_config/buildnml:827-845`
      (fmeSeaiceDerivedFields)

    Each has `filename_template="...$Y-$M.nc"` and
    `filename_interval="00-01-00_00:00:00"`. The Registry XMLs were
    also updated for consistency, but the buildnml strings are what
    actually drive runtime. The remap path
    (`mpas_ocn_fme_horiz_remap.F:722`, `ocn_fme_remap_build_filename`)
    inherits the template via `MPAS_stream_mgr_get_property` -- so
    no Fortran change in the remap layer either. Daily records
    still accumulate (`output_interval="00-00-01"`), just bundled
    into monthly files. For 100-yr production: 1200 EAM files +
    4800 MPAS files = ~6000 total, vs 182,500 daily files.
    Per-file size: ~1.5-1.8 GB.

    **First-time gotcha (2026-04-30)**: I initially edited only the
    Registry XMLs and the runtime files still emerged as `$Y-$M-$D`.
    The empirical signal was a `Ld4` ERS test producing
    `fmeDepthCoarsening.0001-01-02.remapped.nc` through `...05.nc`
    (still daily). When the EAM monthly file works
    (`eam.h0.0001-01.nc`) but the MPAS files don't, the buildnml
    edits are missing.

    **Live-run gotcha (2026-05-10)**: PIO buffers writes to the
    UNLIMITED time dim. If you `xarray.open_dataset` the *current*
    month's file mid-leg, you'll see `time_dim=0` (and a header-only
    file size, ~3-60 KB depending on stream) even though the AM has
    been writing records for days. The counter doesn't materialize
    on disk until file close — which happens at month rotation,
    restart-write, or run end. Don't panic; the previous-month file
    is already closed and authoritative for what's been written.

29. **MPAS FME restart now works correctly via append-mode reopen +
    accumulator sidecar (FIXED 2026-05-01; ERS_Ld4 PASSes
    COMPARE_base_rest BFB-clean for all FME files).** The original
    behavior was: `check_rotate` always created in `PIO_CLOBBER` mode,
    so leg 2 truncated leg 1's file. AM `*_restart_*` hooks were empty
    stubs, so accumulators reset on warm restart. Combined, ERS would
    produce a 1-record `.rest` (Time=1) vs phase-1's 4-record `.base`,
    failing cprnc.

    **The fix has four pieces, all required:**

    1. **Append-mode reopen** in `mpas_ocn_fme_horiz_remap.F` and
       `mpas_seaice_fme_horiz_remap.F`. `remap_file_open` now
       inquires whether the target filename exists. When it does AND
       `is_restart_run==.true.` (cached at module init from
       `config_do_restart`), it opens via `pio_openfile(PIO_WRITE)`
       instead of `pio_createfile(PIO_CLOBBER)`. Dim/var IDs are
       looked up via `pio_inq_dimid` / `pio_inq_varid` instead of
       being defined; `write_time` skips the "first record only"
       lon/lat seed block.

    2. **Frame tracking via numeric `time`**, gated on STRICT-less-than.
       After reopen, scan the existing CF `time(time)` variable and seed
       `time_record = count(existing_time < MPAS_NOW_in_days)`. The next
       `write_time` then increments to the FIRST frame whose existing
       time equals MPAS_NOW (= the leg-1 leftover record past the
       restart point) and OVERWRITES it. With `<=` instead of `<`, that
       frame would be excluded from the count and leg-2 would append
       past the end, producing a duplicate time — a very subtle false
       PASS where cprnc compares `min(time)` records on each side and
       finds them all matching because `.rest`'s frames 1..N are
       literally leg-1's bytes. (Pre-2026-05-04 this used the legacy
       MPAS `xtime` character var with a string lex-compare; switched
       to numeric `time` when xtime was removed in gotcha #39.)

    3. **Accumulator sidecar files** (see gotcha #11) — required so
       leg-2's first averaging window after restart starts from the
       same `(remap_accum, nAccum, avg_last_output_time)` state that
       leg 1 had at the restart-write moment. Without this, leg-2's
       record at the restart-spanning xtime would average over a
       smaller subset of timesteps than leg 1 did.

    4. **`compute_on_startup=.false.`** for all four FME AMs (see
       gotcha #11 for the asymmetric-extra-sample story). This is
       the warm-restart-BFB closer.

    **Production cadence still benefits from calendar alignment**:
    `REST_OPTION='nyears'` with monthly file rotation ensures leg 2
    opens a fresh `*.YYYY-01.nc` after a year-end restart — the
    append-mode branch is exercised on the rare mid-month-restart
    case but not on the common year-boundary case. Documented in
    `run_e3sm.fme.sh`.

    **Hybrid + RESUBMIT chains work correctly:**
    - Leg 1 (hybrid initial): `config_do_restart=.false.` →
      clobber-create everything fresh.
    - Subsequent legs: warm restart, sidecars loaded, monthly files
      either don't exist yet (year-end alignment → fresh create) or
      get appended to via the `<` xtime overwrite logic
      (mid-window restart, rare).

30. **`eam_derived` only sees state/pbuf/constituents -- not `addfld`
    diagnostics.** `validate_field_name` (`eam_derived.F90:801-866`)
    walks 5 sources: known 3D/2D state vars, prior derived fields,
    constituents (`cnst_get_ind`), and physics buffer (`pbuf_get_index`).
    Fields registered via `addfld` and populated only via `outfld` --
    PRECSC, PRECSL, PRECC, PRECL, FLNT, FLDS, most cam_diagnostics
    outputs -- are NOT in any of those. So
    `derived_fld_defs = 'X=PRECSC+PRECSL'` aborts at init with
    "unknown field PRECSC". Discovered 2026-04-30 trying to define
    `PRECST = PRECSC+PRECSL`. Workaround: add the sum as a real
    `addfld`+`outfld` next to the existing diagnostic. PRECST is
    now natively output via `cam_diagnostics.F90:627,2050` as
    `snowc + snowl`, and `'PRECST:A'` in fincl1 references that
    native field directly (no derived-expression hop). Before
    adding any new derived expression: grep the operand -- if it
    only appears in `addfld`/`outfld`, you cannot reach it via
    eam_derived; write the sum in physics dispatch instead.

31. **`*at2m` / `*at10m` height-interp fields restart bug — FIXED
    and verified 2026-04-30 (ERS_Ld4 PASSes COMPARE_base_rest, all 1478
    EAM fields IDENTICAL incl. Tat2m, Qat2m, STWat2m, Uat10m,
    Vat10m).** Root cause was gotcha #2 (Fortran sequence association
    stride): `eam_vcoarsen.F90` passed `src_field` and `coord_mid`
    (declared `(pcols, pver)`) to `shr_vcoarsen_select_nearest` /
    `shr_vcoarsen_select_index` whose dummies are `(ncol, nlev)`. With
    `ncol < pcols`, levels k>1 read shifted memory; for `i ≤ pcols-ncol`
    at `k=2` the shifted slot lands in the zeroed tail
    `field_out(ncol+1:pcols, 1)` set by `get_state_field` (line 577),
    producing the spurious zeros cprnc reported. The avg path was BFB
    because it already passed explicit slices (`src_field(1:ncol, :)`,
    `coord_iface(1:ncol, :)`). Fix: all three select paths (sel_lev at
    line 452, sel_pres at line 475, sel_height at line 501) now pass
    `src_field(1:ncol, :)`, `coord_mid(1:ncol, :)`, `nlev_max(1:ncol)`,
    and `selected(1:ncol)` -- mirrors the avg pattern exactly.

    **Workaround if fix regresses**: drop Tat2m, Qat2m, STWat2m,
    Uat10m, Vat10m from `fincl1` in `shell_commands`. SamudrACE
    doesn't strictly need them -- TREFHT/QREFHT cover 2 m temperature
    and humidity, and ACE training pipelines often derive 10 m wind
    from level-7 coarsened U/V (boundary layer) anyway.

32. **MPAS native FME files are disabled in production.** As of 2026-04-30,
    the four FME stream blocks in `cime_config/buildnml` are emitted with
    `output_interval="none"` -- this prevents MPAS from writing the native
    ne30 mesh `*.fmeDepthCoarsening.YYYY-MM.nc` etc. files (saving ~3-5×
    the storage of the remapped files), while leaving the stream metadata
    (filename_template) registered so `mpas_ocn_fme_horiz_remap.F`'s
    `MPAS_stream_mgr_get_property` query still succeeds and constructs
    the `.remapped.nc` filename correctly. Native files were a duplicate
    of the remapped output on a different grid -- not needed for the
    SamudrACE training tape. To re-enable for debugging, change
    `output_interval="none"` to a real interval like `"00-00-01_00:00:00"`.

33. **`fmeVerticalReduce` AM is disabled in production.** Its three
    outputs (`oceanHeatContent`, `freshwaterContent`, `kineticEnergy`) are
    column-integrated diagnostics not in the SamudrACE Confluence spec
    (p3ai/6154289880). The testmod sets
    `config_AM_fmeVerticalReduce_enable = .false.`; the Fortran AM
    (`mpas_ocn_fme_vertical_reduce.F`) and the Registry stay in place
    so it can be re-enabled if column-integrated diagnostics become
    useful for climate-drift monitoring or independent ocean-state
    sanity checks. OHC can also be derived post-hoc from
    `temperatureCoarsened_*` + `layerThicknessCoarsened_*` in the
    `fmeDepthCoarsening` tape.

34. **MPAS FME `.remapped.nc` files now compared in CIME tests (ERS etc).**
    `cime_config/config_archive.xml` previously had `exclude_testing="true"`
    on `mpaso` and `mpassi`, which caused `hist_utils.py:79` to skip both
    components entirely whenever `TEST=TRUE`. As of 2026-04-30 that
    attribute is removed. `<hist_file_extension>` is now
    `hist\..*\.nc$` -- *not* the bare `hist` -- because of a subtle CIME
    regex-construction bug: when `_component_compare_copy` invokes
    `copy_histfiles(case, "base", match_suffix="nc")`,
    `archive_base.py:142-148` builds a literal `\.nc$` suffix onto the
    regex *unless* `nc` already appears in it OR the ext ends in `$`.
    Bare `hist` triggered the suffix append, producing
    `mpaso\d?_?(\d{4})?\.hist\.nc$` -- which matches NO real MPAS file
    (FME files end `.remapped.nc`, native streams end
    `.YYYY-MM-DD.nc`). The narrow regex (`hist\.am\.fme\w+\.\S+\.remapped\.nc$`)
    worked by luck because it ended in `$`; the broad regex must too,
    or contain `nc`. Pattern `hist\..*\.nc$` matches everything we want
    (FME remapped + native streams), excludes `.rst`/`.rest`/`.base`/`.cprnc.out`,
    and produces the same regex regardless of `match_suffix`.
    Per-stream `<hist_file_ext_regex>` entries are added for the FME streams
    (`hist\.am\.fmeDepthCoarsening`, `hist\.am\.fmeDerivedFields`,
    `hist\.am\.fmeVerticalReduce`, `hist\.am\.fmeSeaiceDerivedFields`)
    so each FME stream gets its own ext-bucket and all of them survive
    the `archive_base.py:111-118` `latest_files[ext]=hist` collapse.
    Native non-FME streams still collapse into the default `\w+` →
    `"hist"` bucket, so only one (alphabetically last) gets compared --
    deemed acceptable because (a) those streams are upstream MPAS, not
    FME-specific, and (b) ALL of them get archived correctly.
    A case must be re-created after this edit for the env_archive.xml
    in its CASEROOT to pick up the new spec; pre-existing case dirs
    will still skip MPAS comparison.

35. **`_inst` companion writes for fmeDepthCoarsening are SPEC-MANDATED,
    not optional.** The Samudra training spec
    (Confluence p3ai/6154289880, "case run options for v3 LR piControl
    SamudrACE") explicitly lists "1D avg & 1D instant" for all four 3D
    ocean state vars (temperature, salinity, velocityZonal,
    velocityMeridional). `config_AM_fmeDepthCoarsening_write_instantaneous_companion=.true.`
    in `shell_commands:217` is the production setting; the Registry
    default is `.false.` (off for non-FME consumers, on for this
    testmod). Storage cost: ~900 GB extra over a 100-yr run for the
    fmeDepthCoarsening stream (~50% of that stream's volume). Do not
    turn this off without re-checking with the SamudrACE team --
    investigated 2026-04-30 and confirmed the spec demands both modes.
    The `_inst` channel is otherwise unused by `verify_mpas.py` (open
    hardening item: add presence/range check for `*_inst` companions).

36. **MPAS time-averaging needed per-cell valid-sample counts, not a
    global `nAccum`. FIXED 2026-05-04 (sea-ice AM had observable
    pollution; ocean AMs fixed defensively).** The original pattern was

    ```
    remap_accum(:, idx) = remap_accum(:, idx) + fld_remap(:)   ! per timestep
    fld_out             = remap_accum(:, idx) / nAccum         ! at output trigger
    ```

    `apply_masked` (gotcha #5) sets a target cell to `SHR_FILL_VALUE`
    when *all* its source cells are fill. For ocean cells whose source
    mask is static (land), this means the cell is fill on *every*
    timestep, and `nAccum * 1e20 / nAccum` cleanly preserves fill.
    But for **sea ice**, the source-side ice-presence mask flips
    day-to-day at the ice edge: a target cell can be valid on day 1
    (some ice present), fill on day 2 (no ice), valid on day 3, etc.
    The accumulator then ends up with `k * 1e20 + (nAccum - k) * v`,
    and dividing by `nAccum` produces *fractional-of-fill* values
    like `1e20/24 = 4.17e18` — outside the fill-detection threshold
    used by `verify_mpas.py` (`1e18`) but very much not physical.

    **Empirical signature**: in
    `*.fmeSeaiceDerivedFields.YYYY-MM.remapped.nc`, the six
    intensive fields (`iceThicknessMean`, `surfaceTemperatureMean`,
    `airStress{Zonal,Meridional}`, `{u,v}VelocityGeoCell` — all set
    to `fillValue` in the AM when `iceAreaTotal == 0`) had ~3000 cells
    per file with `1e18 < |x| < 1e20` whose values were exact
    multiples of `1e20/nAccum`. The three extensive sums
    (`iceAreaTotal`, `iceVolumeTotal`, `snowVolumeTotal`) were
    pollution-free because the AM keeps them at `0.0` (not
    `fillValue`) when no ice is present.

    **Fix**: parallel `remap_accum_count(:, :)` real(r8) array.
    Each `remap_*` subroutine detects fill in `fld_remap`
    (`abs(val) < SHR_FILL_VALUE * 0.1_r8`) and only accumulates
    valid samples, incrementing the per-cell count. The write path
    divides by the per-cell count and emits fill where count == 0.
    Sidecar restart-state files now also persist the count array
    so warm-restart BFB still holds (gotcha #29) — both
    `ocn_fme_horiz_remap_{read,write}_accum` and the seaice
    counterparts gained an optional `count_arr=` argument.

    **Sidecar format**: the new variable is `remap_accum_count`,
    same shape and PIO decomposition as `remap_accum`. It is
    optional: callers that omit `count_arr=` get the legacy
    behavior (no count read/write). The sea-ice and ocean AMs all
    pass `count_arr=remap_accum_count`.

    **Consequence for cold runs**: cells that flip valid↔fill
    within an averaging window now produce a true average over
    only the valid samples rather than a fill-polluted ratio.
    Spatial coverage of the averaged field can change very
    slightly at the ice edge (an ice-edge cell that previously
    showed `2e19` will now show its actual mean over the valid
    days); this is a correctness improvement, not a regression.

    **Defensive fix on ocean side**: `mpas_ocn_fme_derived_fields.F`,
    `mpas_ocn_fme_depth_coarsening.F`, and
    `mpas_ocn_fme_vertical_reduce.F` were also converted to the
    per-cell-count pattern even though their masks are static
    (land / per-layer bathymetry). This keeps the four AMs
    coherent and protects against future fields whose validity
    flips in time (e.g., ice-shelf cavity exposure if SSH masking
    is added under cavities).

    **What it does NOT fix**: sub-window pollution of the *source*
    field before remap (the AM still sets ice-free cells to
    `fillValue` on the native mesh, and `apply_masked` still
    correctly handles that at the spatial level). The fix is purely
    in the time-axis accumulator step.

37. **Static `mask_2d` / `mask_<k>` written into ocean FME files for
    SamudrACE compatibility (ADDED 2026-05-04).** Each ocean output file
    now ships canonical wet masks alongside the data:

    - `fmeDerivedFields*.remapped.nc` -> `mask_2d(lon, lat)`
    - `fmeDepthCoarsening*.remapped.nc` -> `mask_0(lon, lat)`,
      `mask_1(lon, lat)`, ..., `mask_{nCoarsenLevels-1}(lon, lat)`
      (production: `mask_0` through `mask_18`)

    Values: `1.0` = valid ocean cell, `0.0` = land/below seafloor.
    Stored as `PIO_REAL` (float32). No `Time` dim — the masks are
    static for the AM's lifetime.

    **Why these specific names**: the SamudrACE preprocessing script
    `compute_ocean_dataset_e3sm.py` (lines 46-65) constructs masks with
    exactly this naming -- `mask_2d` for any 2D variable, `mask_<k>` for
    any 3D variable ending in `_<k>`. By emitting the masks in the FME
    output directly, we let the downstream training pipeline read them
    as-is and skip the `~ds[var].isnull().any(dim="time")` derivation
    step. The dataloader's lookup is `re.compile(r"_(\d+)$")` then
    `f"mask_{level}"` for 3D / `"mask_2d"` for 2D
    (`compute_ocean_dataset_e3sm.py:get_mask_for`).

    **How they're computed**: in each AM, on the first compute call
    after `do_horiz_remap` is active, a constant-1 native field is
    remapped through `apply_masked` (gotcha #5) with the correct
    source-side mask:

    - `fmeDerivedFields`: `mask=maxLevelCell` (static land mask).
    - `fmeDepthCoarsening`: per level, `mask=levelMask` derived from
      `layerThicknessCoarsened_k != fillValue` (static bathymetry mask).

    The remapped value is `1.0` on target cells with any valid source
    coverage and `SHR_FILL_VALUE` elsewhere. Threshold to a hard
    `1.0`/`0.0` and cache in module-scope arrays
    (`mask_2d_remap`, `mask_levels_remap`).

    **When they're written**: at the first time record of every fresh
    file, gated on `time_record == 1 .and. .not. existing_file_reopened`
    (warm-restart reopens skip — masks are already in the file from
    leg 1). Helper `write_static_masks_if_needed()` lives inside the
    compute routine of each AM.

    **PIO/file-format pieces**: `ocn_fme_remap_file_t` gained
    `stored_is_static(:)` and `pio_var_is_static(:)` boolean arrays
    plus an optional `is_static=` keyword on `register_var` and
    `def_var`. Static vars are defined with shape `(lon, lat)` only
    (no Time dim), have no `_FillValue` attribute, and `write_var`
    skips `pio_setframe` for them. Append-mode reopen (warm restart)
    looks them up via `pio_inq_varid` like any other var.

    **Storage cost** (production, 100 yr):
    - 1 `mask_2d` × 6000 files × 256 KB = ~1.5 GB total.
    - 19 `mask_k` × 6000 files × 256 KB = ~29 GB total.
    Combined ~30 GB; negligible vs. the ~6 TB total tape volume.

    **Sea-ice and fmeVerticalReduce are not covered** — Elynn's request
    was explicit ocean-only (Slack 2026-05-04). The same mechanism
    transfers cleanly if needed later: in seaice, `mask_2d` would be
    derived from `iceAreaTotal`'s ocean-mask pattern; for vertical
    reduce, from `maxLevelCell`. The horiz_remap helper's
    `is_static` plumbing is ocean-side only — the seaice helper would
    need the parallel addition.

    **What it does NOT include**: `mask_2d` is the *land* mask only.
    Time-varying validity (e.g., sea-ice presence at the ice edge,
    SSH below ice shelves) is encoded in the data variables' fill
    pattern itself, not in the static mask. The per-cell-count
    accumulator (gotcha #36) handles that correctly during time
    averaging.

38. **`totalFreshWaterTemperatureFlux` removed from `fmeDerivedFields`.**
    The SamudrACE Confluence spec strikes through this field. Removed
    2026-04-30 from `mpas_ocn_fme_derived_fields.F` (register_var,
    pointer declaration, mpas_pool_get_array, mask_and_remap call,
    write_avg_field call, inst_remap_and_write call), the Registry
    XML var list, and the `cime_config/buildnml` stream block. The
    accumulator slot at idx=12 in `remap_accum` is left allocated but
    unused for now (renumbering the 24 indices was deemed not worth
    the churn risk). Reusable for any future field at index 12.

39. **MPAS FME files use CF time variables matching EAM (ADDED 2026-05-04;
    epoch plumbing simplified to runtime auto-discovery later that
    same day -- see "Epoch plumbing" subsection below for the why).**
    The legacy `xtime(StrLen, Time)` character timestamp was replaced
    with the EAM-style trio:

    ```
    dimensions:
        time = UNLIMITED ;   nbnd = 2 ;   lon = 360 ;   lat = 180 ;
    variables:
        double time(time) ;
            time:units    = "days since {RUN_STARTDATE} 00:00:00" ;
            time:calendar = "noleap" ;
            time:bounds   = "time_bnds" ;
            time:axis     = "T" ;
        double time_bnds(time, nbnd) ;
        int    date(time) ;       // YYYYMMDD
        int    datesec(time) ;    // seconds of day
    ```

    **Why**: bytewise time alignment with EAM for online ACE coupling.
    `xtime` was redundant given `time + time_bnds` (CF-canonical) and
    a footgun: legacy MPAS `timeSeriesStatsCustom` stamps xtime at the
    START of an averaging window while our FME code stamps at the END,
    so anyone applying MPAS muscle memory to our files got it wrong
    by one window. Dropping xtime eliminates that trap.

    **Time stamp convention**: numeric `time` is the END of the
    averaging window (matches EAM); `time_bnds = [t_lo, t_hi]` carries
    both endpoints. For instantaneous records (`*_inst` companion
    fields and inst-mode AMs), `time_bnds = [t, t]`.

    **Epoch plumbing -- auto-discovered at runtime, no namelist needed.**
    Each FME AM calls `ocn_fme_horiz_remap_init_time_epoch(domain,
    amName)` (or the seaice equivalent) once at init. The first caller
    populates the module-scope epoch state; subsequent AMs see
    `time_epoch_set=.true.` and no-op. Discovery order:

    1. **Sidecar lookup** (`try_read_epoch_from_sidecar`): if the
       per-AM accumulator sidecar `<amName>.fme_accum_restart.nc`
       exists and has a `time_reference_date` global attribute,
       use that. This is the warm-restart path -- the epoch baked
       into the sidecar by the prior cold leg is the case start
       date, and reading it back ensures legs after a restart use
       the same epoch they started with.
    2. **Clock fallback**: query `MPAS_START_TIME` from
       `domain%clock` and convert to a date-time string. On a cold
       run this equals CIME `RUN_STARTDATE` -- verified at
       `ocn_comp_mct.F:611-614` where the driver calls
       `seq_timemgr_EClockGetData(EClock, ECurrTime=...)` for
       `runtype == 'initial'` and stores the result as
       MPAS_START_TIME. On warm restart the same code path runs
       but with `runtype == 'continue'` and the current restart
       time -- which is why the sidecar lookup wins on warm legs.

    The sidecar persistence is automatic: each
    `*_write_accum` call (which fires on every coupler-driven
    restart-write, see gotcha #11) writes the current
    `time_epoch_str` as a global attribute alongside `nAccum` and
    `avg_last_output_time`. Cost: one PIO put_att per restart write,
    negligible.

    **Failure mode** (warm restart with no sidecar): if the
    sidecar is absent on a warm restart (e.g., user nuked the run
    dir between legs), the AM falls back to MPAS_START_TIME =
    leg-start time, NOT case start. The result would be MPAS files
    with `time:units` shifted relative to EAM. Logged as a loud
    `MPAS_LOG_WARN` so the issue surfaces in mpaso/mpassi.log.

    **Why we replaced the namelist plumbing (2026-05-04 evening,
    same day as initial CF time rewrite):** an earlier revision
    used per-AM `config_AM_*_time_reference_date` namelist options
    set by `shell_commands` from `xmlquery RUN_STARTDATE`. That
    plumbing had two real problems: (1) four near-duplicate flags
    that had to agree, validated at runtime; (2) the production
    script `run_e3sm.fme.sh` runs `xmlchange RUN_STARTDATE` AFTER
    `create_newcase` -- which is when shell_commands fires and
    bakes the (then-default) epoch into `user_nl_mpa*`. Production
    runs would have silently shipped MPAS files with a bogus
    epoch (e.g., `0001-01-01` from the WCYCL1850NS compset
    default) while EAM used the correct case start date. The
    auto-discovery design eliminates both problems by not needing
    case-construction-time knowledge of the epoch at all.

    **Numeric conversion**: `MPAS_Time_Type` to days-since-epoch goes
    through `mpas_get_timeInterval(t - epoch, StartTimeIn=epoch,
    dt=seconds)` then `seconds / 86400`. ESMF respects the noleap
    calendar exactly. Helpers are
    `time_to_days_since_epoch` (ocean) and
    `si_time_to_days_since_epoch` (seaice). Both modules also store
    the epoch as a CF-form string (`YYYY-MM-DD HH:MM:SS` with a
    SPACE, not the MPAS underscore) for the `time:units` attribute.

    **CF units underscore-vs-space**: MPAS xtime form is
    `YYYY-MM-DD_HH:MM:SS`; CF "days since" requires
    `YYYY-MM-DD HH:MM:SS`. The internal `set_time_epoch` /
    `set_time_epoch_si` helper substitutes the underscore on the fly
    into the `time_epoch_cf` string used in attributes.

    **Append-mode reopen**: see gotcha #29 point 2 for the parallel
    update -- frame tracking now uses numeric `time` instead of
    string `xtime`.

    **`write_time` signature**: changed from
    `write_time(time_str)` (string-based, end-of-window only) to
    `write_time(time_lo, time_hi)` taking two `MPAS_Time_Type`
    arguments. AM call sites at all 8 locations (avg path + inst
    path × 4 AMs) updated.

40. **Per-variable units / cell_methods + global provenance attrs filled
    in (ADDED 2026-05-05).** A full audit found three CF-attribute
    classes of bug across all four FME streams; all fixed in one pass:

    1. **EAM `eam_vcoarsen.F90` hardcoded `'varies'`** at all 6
       `addfld` register sites — placeholder that was never filled in.
       The vcoarsen avg / sel_lev / sel_pres / sel_height paths are
       mathematically unit-preserving, so the output unit equals the
       source. Replaced with a `lookup_field_units(name)` helper that
       maps each known source field (T→K, U/V→m/s, OMEGA→Pa/s, Z3→m,
       Q/CLDLIQ/CLDICE/RAINQM/SNOWQM/STW→kg/kg, NUM*→1/kg) and warns
       at masterproc with `units="unknown"` for unrecognised inputs.
       The column-integrated `*_INT` and `d*_INT_dt` paths use a
       `compose_int_units` helper that multiplies through by `kg/m^2`
       (and `/s` for the tendency) since vcoarsen integrates dp/g.
       Production: Q_k, T_k, U_k, V_k, STW_k, *at2m/at10m all carry
       the right unit string in eam.h0.YYYY-MM.nc.

    2. **MPAS `mpas_ocn_fme_depth_coarsening.F` empty units** (`''`)
       on the four 3D state vars and their `_inst` companions.
       Replaced with a `coarsen_field_units(name)` helper at the
       module bottom mapping `temperature→degC`, `salinity→PSU`,
       `velocityZonal/velocityMeridional→m s-1`. Falls back to
       `'1'` (CF dimensionless) with a WARN if a future field is
       added without a units mapping.

    3. **No `cell_methods` on any data var.** Added optional
       `cell_methods=` keyword to `register_var` and `def_var` in
       both `mpas_ocn_fme_horiz_remap.F` and
       `mpas_seaice_fme_horiz_remap.F`. The four AMs read their
       `config_AM_*_time_averaging` flag once at register time and
       pass `'time: mean'` (averaged-mode) or `'time: point'`
       (instantaneous-mode); depth-coarsen `_inst` companions
       always get `'time: point'`. Static masks (`mask_2d`,
       `mask_*`) intentionally get no `cell_methods` (no time dim).

    Also fixed in this pass:
    - **`surfaceTemperatureMean` was `'K'` in the seaice AM** while
      the Registry XML and the actual field are degC (gotcha #21).
      Fixed to `'degC'`. Verify dashboards already used the correct
      Celsius range; this was only a metadata mismatch.
    - **`iceAreaTotal` was `'unitless'`** — replaced with `'1'` for
      CF compliance.

    **MPAS global attrs enriched** to match EAM's set. Previously
    only `source / Conventions / calendar / time_reference_date`;
    now also `title / realm / product / institution_id / contact /
    case / username / hostname / git_version / history`. Plumbing:

    - Added `prov_*` module-scope strings + `*_set_provenance(case,
      username, hostname, model_version, history)` setter in both
      `mpas_ocn_fme_horiz_remap.F` and `mpas_seaice_fme_horiz_remap.F`.
    - `ocn_comp_mct.F` and `ice_comp_mct.F` call the setter once
      during init, right after `seq_infodata_GetData(...)` populates
      caseid/model_version/username/hostname. Idempotent — first
      caller wins, subsequent calls no-op (guards against MPAS
      framework double-init).
    - `remap_file_open` writes the prov_* strings via
      `pio_put_att(PIO_GLOBAL)` on every fresh-create. Append-mode
      reopen skips them (already on disk from leg 1).

    **Why this matters**: downstream tools (xarray, ncview, ESMValTool)
    use these attrs for axis labels, plot titles, and provenance
    tracking. CF-aware ML pipelines like SamudrACE will use
    `cell_methods` to distinguish averaged-vs-instantaneous records
    when concatenating tapes.

    **Files touched (10)**: `eam_vcoarsen.F90`, `mpas_ocn_fme_horiz_remap.F`,
    `mpas_seaice_fme_horiz_remap.F`, `ocn_comp_mct.F`, `ice_comp_mct.F`,
    `mpas_ocn_fme_depth_coarsening.F`, `mpas_ocn_fme_derived_fields.F`,
    `mpas_ocn_fme_vertical_reduce.F`, `mpas_seaice_fme_derived_fields.F`,
    plus this AGENTS.md gotcha. Verified end-to-end with ERS_Ld4
    test 52484623 (RUN + COMPARE_base_rest BFB-clean).

41. **Per-layer-interface vertical coordinate scalars (ADDED 2026-05-05).**
    The flat `Q_k`/`T_k`/`temperatureCoarsened_k` data fields have no
    per-layer hint about *which pressure or depth they correspond to*.
    Added scalar metadata fields so downstream consumers can reconstruct
    pressure (atm) or depth (ocn) directly from the file.

    **EAM h0**: `ak_0..ak_N` and `bk_0..bk_N` (N = `n_avg_levs`).
    Production: 9 of each (interfaces of the 8 vcoarsen layers).
    True 0-dim netCDF scalars in `ds.data_vars` (no dim, like `P0`).
    Dimensionless `units='1'`. Pressure at vcoarsen interface k:
    `p = ak_k * P0 + bk_k * PS`. In level-index mode (the production
    setting), values come straight from `hyai`/`hybi` at the native
    interface index `vcoarsen_level_bounds(k)`. In pressure-bounded
    mode, an `interp_native_interface` helper does linear interpolation
    against `p_iface(k) = hyai(k)*P0 + hybi(k)*PS_ref` so the
    coefficients are still PS-independent.

    **EAM plumbing — `add_hist_scalar` API** (NEW in
    `cam_history_support.F90`). Registers a true 0-dim netCDF scalar
    variable that lands in `ds.data_vars` (xarray) — modeled on the
    P0 path (lines 1627-1638 / 1771-1779) but available to any
    caller. Public signature:
    `call add_hist_scalar(name, value, units, long_name)`. A parallel
    `hist_scalars(maxhistscalars=64)` array stores the registered set;
    `define_hist_scalars(File)` and `write_hist_scalars(File)` are
    called from inside the existing `write_hist_coord_attrs` and
    `write_hist_coord_vars` so `cam_history.F90` needs no changes.

    **Why not `add_hist_coord(..., vlen=1)`?** That works (creates a
    1D var of length 1) but consumes one slot per scalar in the 25-cap
    `hist_coords` table. With 18 new scalars (9 ak + 9 bk) atop the
    baseline ~7-9 entries, Crux/ALCF blew the cap with
    `'Too many dimensions in add_hist_coord.'`; Perlmutter was just
    under by luck. The 0-dim scalar path eliminates the cap concern
    entirely AND classifies as `data_vars` (matching `Q_k`, `T_k`)
    rather than `coords`. The previous `add_hist_coord` approach was
    abandoned 2026-05-05 in favor of `add_hist_scalar`.

    **MPAS-O fmeDepthCoarsening**: `idepth_0..idepth_N`
    (N = `nCoarsenLevels`). Production: 20 (interfaces of the 19
    coarsened layers). True 0-dim PIO_DOUBLE scalars in `ds.data_vars`.
    `units='m'`, `positive='down'`, `standard_name='depth'`. Layer k
    spans `idepth_k → idepth_{k+1}`; thickness equals
    `layerThicknessCoarsened_k` (sanity invariant).

    **MPAS plumbing — `is_scalar` shape on `register_var` / `def_var`**.
    `ocn_fme_remap_file_t` gained an `is_scalar=.true.` flag with a
    `scalar_value=` kwarg, parallel to the existing `is_static`.
    Scalars are defined with an empty dim list and written via
    `pio_put_var(file, desc, value)` at define-time; PIO transitions
    from def-mode to data-mode automatically. The `stored_scalar_value`
    array carries the constant across rotation so `check_rotate`
    re-emits each scalar in every fresh file. `write_var` errors
    gracefully if accidentally invoked on a scalar. Sea-ice horiz_remap
    module has no `is_scalar` plumbing (no scalar metadata needed).

    **Why interfaces, not midpoints**: interfaces are what vcoarsen /
    depth-coarsen actually integrate over. With N+1 interface values
    a downstream tool can derive layer means, midpoints, and
    thicknesses; with N midpoint values it cannot recover interfaces.
    Interface form has zero aggregation logic — values come straight
    from the source's interface arrays.

    **What it does NOT cover**: `*at2m` / `*at10m` height-interp
    fields keep the height in their `long_name` only (per user
    direction) — no scalar `height_2m`/`height_10m` companion. Add
    later if a CF-strict reader needs it.

    **Verification example**: `ds.ak_0.values.item() ≈ 0.0001` (TOA),
    `ds.bk_8.values.item() == 1.0` (surface), `ds.idepth_0.item() ==
    0.0`, `ds.idepth_19.item() == 6380.0`. Reconstruct any layer's
    top pressure: `p_top = ds.ak_k.item()*ds.P0 + ds.bk_k.item()*ds.PS`.

42. **Height-interp field rename + 5-day MPAS companion stream
    (ADDED 2026-05-09).** Two related changes landed together; ERS_Ld10
    PASSes COMPARE_base_rest BFB-clean (job 52745629).

    **Variable rename — `_at_z`/`_at_P`/`_at_L` → no underscores.**
    The three `make_sel_*_name` helpers in `eam_vcoarsen.F90` now
    produce names without internal underscores so the only `_<token>`
    suffix on output names is the level-index `_<k>` form
    (`T_0..T_7`). Mapping:

    - `T_at_z2` → `Tat2m`,    `U_at_z10` → `Uat10m` (height-interp,
       `m` suffix indicates meters above surface)
    - `U_at_P850` → `Uat850hPa` (pressure-level select)
    - `U_at_L5` → `UatL5` (level-index select; no unit suffix)

    Production fincl1 in `shell_commands` accordingly: `Tat2m, STWat2m,
    Qat2m, Uat10m, Vat10m`. The `namelist_definition.xml` doc strings,
    AGENTS.md citations (gotchas #31, #40, this entry), and the
    workaround note in #31 were all updated. Only EAM is affected; no
    MPAS field renamed. **No back-compat shim** — older outputs with
    the old names need a downstream `xr.rename({...})` if loaded
    alongside post-2026-05-09 files.

    **5-day-mean MPAS companion stream — `*Output5D`.** SamudrACE wants
    BOTH 1-day and 5-day means as outputs. Each FME AM (depth-coarsening,
    derived-fields ocean, derived-fields seaice) now optionally
    maintains a SECOND accumulator and output file at a longer cadence,
    sharing the same per-step samples as the primary 1D path but writing
    averaged-only (no `_inst` companion). Driven by a new namelist
    option:

    ```
    config_AM_fmeDepthCoarsening_stream_output_interval_5D = '00-00-05_00:00:00'
    config_AM_fmeDerivedFields_stream_output_interval_5D = '00-00-05_00:00:00'
    config_AM_fmeSeaiceDerivedFields_stream_output_interval_5D = '00-00-05_00:00:00'
    ```

    Set to `'none'` (or omit) to disable. The testmod sets all three
    via a single `FME_MPAS_INTERVAL_5D` env var (default `00-00-05`).
    Output filenames: `*.fmeDepthCoarsening5D.YYYY-MM.remapped.nc`
    etc. Sidecar accumulator restart files use distinct amName keys
    (`fmeDepthCoarsening5D`, `fmeDerivedFields5D`,
    `fmeSeaiceDerivedFields5D`) so 1D and 5D state don't collide.

    **Plumbing changes** (touched files):
    - `Registry_fme_depth_coarsening.xml`,
      `Registry_fme_derived_fields.xml`,
      `Registry_seaice_fme_derived_fields.xml`: new `<nml_option
      name="config_AM_*_stream_output_interval_5D">` and a parallel
      `<stream name="*Output5D">` declaration.
    - `namelist_definition_mpaso.xml` /
      `namelist_definition_mpassi.xml` AND
      `namelist_defaults_mpaso.xml` /
      `namelist_defaults_mpassi.xml`: new entry per AM (default
      `'none'`). **Both files matter** — the build-time MPAS namelist
      validator uses these XMLs (separate from Registry XML); a
      missing entry yields `"either config_am_X is not a valid MPASSI
      namelist variable or X = 'val' is not a valid value"` at the
      buildnml SETUP phase. We hit this once and added them.
    - `cime_config/buildnml` (mpas-ocean and mpas-seaice): emit the
      `*Output5D` stream blocks with `output_interval="none"` (same
      pattern as the primary streams — buildnml is authoritative for
      filename_template; native ne30 file is suppressed; the
      `.remapped.nc` is what the AM produces via the file's
      filename_template).
    - **MPAS-O AMs** (`mpas_ocn_fme_depth_coarsening.F`,
      `mpas_ocn_fme_derived_fields.F`): clean — the existing
      `ocn_fme_remap_file_t` already supports multiple instances, so
      we just add a second `remap_outfile_5D` plus parallel
      `remap_accum_5D` / `remap_accum_count_5D` arrays + a duplicate
      output-trigger block in `compute`.
    - **MPAS-SI horiz_remap** (`mpas_seaice_fme_horiz_remap.F`)
      pre-dates the type-encapsulated refactor, so the per-file PIO
      state is still at module level. Rather than refactor 1200 lines
      under time pressure, we added a parallel set of `_5D` state
      vars + `_5D` procedures (`open_file_5D`, `def_var_5D`,
      `write_var_5D`, `write_time_5D`, `close_file_5D`,
      `register_var_5D`, `check_rotate_5D`). The accum-sidecar
      routines are amName-keyed and reused via thin `_5D` aliases.
      **TODO/cleanup**: a future refactor to a
      `seaice_fme_remap_file_t` type would erase the duplication and
      bring the seaice path in line with the ocean side.

    **Cross-check semantics — sample-weighted vs day-weighted means.**
    `verify_production.py` grew a 5D-vs-mean(1D) cross-check, but the
    comparison is subtler than it looks at first. By construction:

    - `5D[j] = sum(per-step samples in window j) /
       count(per-step samples in window j)`  → *sample-weighted*.
    - `mean(1D[i for i in window j]) = (1D[j*5] + ... + 1D[j*5+4]) / 5`
       → *day-weighted*.

    These coincide ONLY when each 1D record in the window has the same
    per-cell sample count. At fill-edge cells (sea-ice edge, ice
    runoff, river runoff, ocean–land transitions) the per-cell count
    flips day-to-day and the two means legitimately differ by O(spike
    × count_imbalance). The output files don't carry the per-cell
    counts, so we cannot reconstruct the sample-weighted 1D mean
    exactly without the AM also writing `count_<var>` companions
    (deferred — see "Extensions" below).

    Practical mitigation in the cross-check: filter cells where any 1D
    record in the window is fill-magnitude, then report `p99(diff)`
    (not `max(diff)`) as a percentage of the variable's `p99(|value|)`.
    p99 across cells skips the few count-flicker cells that legitimately
    diverge. Default tolerance is 1% of p99-scale (WARN), 5% (FAIL).
    On ERS_Ld10 all three streams pass cleanly under this metric.

    **Workaround if the 5D feature regresses**: drop the
    `_stream_output_interval_5D` line from `shell_commands` and the
    AMs silently skip the 5D path (`time_avg_5D_enabled` stays
    `.false.`). The primary 1D path is unchanged in this branch and
    is bit-identical to pre-2026-05-09 behavior — verified by
    inspection of the modified files (no edits to the existing 1D
    code paths beyond the dual-accumulate `if (time_avg_5D_enabled)`
    block in `remap_1level` / `remap_1field`, gated on a flag that
    defaults to `.false.`).

43. **Atmosphere aerosol/cloud/CCN additions + day+night AOD companions +
    CO2 scalar metadata (ADDED 2026-05-24).** Four loosely coupled
    additions to the EAM side of the FME tape:

    **A. `AODVISall` / `AODUVall` / `AODNIRall` — day+night unmasked AODs.**
    `modal_aer_opt.F90:modal_aero_sw` already accumulates AODVIS/AODUV/AODNIR
    for every column (`do i = 1, ncol` at lines 998-1017). The trailing
    `do i = 1, nnite ; aodvis(idxnite(i)) = fillvalue` loops at lines
    1239-1247 (visible) and 1283-1287 (UV/NIR) are *purely cosmetic*
    masking — the optics computation itself is solar-zenith-independent.
    The new `*all` variants snapshot the unmasked arrays via `outfld`
    BEFORE those fill loops run. AODVIS itself is unchanged; downstream
    consumers that depend on the nighttime-masked semantics keep
    working. All three `*all` variants are addfld'd, but only
    `AODVISall` is selected in fincl1 for the production tape
    (`AODUVall` / `AODNIRall` are registered and available — add them
    to fincl1 later if the SamudrACE team wants UV/NIR). The fincl1
    entry carries NO explicit averaging flag, so it inherits the
    per-tape default `avgflag_pertape` (production: `'I'`,
    instantaneous). All three are gated on `list_idx == 0` (climate
    diagnostic list only — no `_1`/`_2` diag list variants).

    Why this matters for ACE/Samudra: every-column AOD lets the training
    pipeline treat aerosol forcing as a continuous field instead of a
    half-sky-fill mask. Production currently emits it at the tape
    default (instantaneous on the radiation-snapshot value); switch the
    fincl1 entry to `'AODVISall:A'` if a time-averaged AOD is wanted
    instead (optics is only recomputed on radiation timesteps, so an
    average folds in the diurnal cycle of aerosol water uptake on
    radiation cadence).

    **B. `do_aerocom_ind3=.true.` enables AEROCOM cloud diagnostics.**
    Off by default in EAM (`phys_control.F90:137`). Turning it on in
    the testmod (`user_nl_eam`) lights up the `output_aerocom_aie`
    register/init path and unlocks: `aerindex` (Angstrom × AOD),
    `angstrm` (Angstrom coefficient), `aod400`/`aod700`, `cdnc` (cloud
    droplet number at cloud top, #/m³), `lwp`/`iwp` (water paths,
    kg/m²), `cod` (cloud optical depth), `clt`/`lcc`/`icc` (cloud
    cover fractions), `cdnumbf`/`icnumbf` (column droplet/ice number
    before microphysics), and `colccn.1`/`colccn.3`/`ccn.1bl`/`ccn.3bl`
    (CCN at S=0.1%/0.3%, column and 1 km AGL).

    Selected for the FME tape (fincl1): `aerindexall`, `cdnc`, `lwp`,
    `ccn.3bl`, `colccn.3`, `lcc`. Note `aerindexall` (the day+night
    companion, section A pattern) is used rather than the
    nighttime-masked `aerindex`; `angstrm`/`angstrmall` are addfld'd
    but NOT selected for the production tape. None of these carry an
    explicit `:A` — they inherit the per-tape `avgflag_pertape`
    (production: `'I'`, instantaneous). The other 15+ aerocom fields
    are addfld'd but not requested in fincl1; expand later if the
    SamudrACE team wants more.

    **RRTMG vs RRTMGP caveat for `aerindex` / `angstrm` / `aod400` /
    `aod700`.** These four are outfld'd only from
    `components/eam/src/physics/rrtmg/radiation.F90:1322-1352`,
    NOT from the RRTMGP path. WCYCL1850 / WCYCL1850NS resolve to
    `eamv3_phys_defaults` which does NOT set `-rad`, so the
    `configure` default `$rad_pkg = 'rrtmg'` (line 1081) wins —
    confirmed via `components/eam/cime_config/config_component.xml`
    (no `-rad` override on `_EAM%CMIP6.*` compsets). If a future user
    switches to `-rad rrtmgp` (e.g., SCREAM or MMF compsets), these
    four fields — AND the day+night `aerindexall`/`angstrmall`
    companions, which also outfld from `radiation.F90` — will be
    silently absent. Since `aerindexall` IS in the production fincl1
    (section B), losing the RRTMG path would drop it from the tape.
    `cdnc`, `lwp`, `lcc`, `colccn.3`, and the `ccn.*` family come from
    microphysics / `output_aerocom_aie` routines that run regardless of
    radiation package, so those are fine either way.

    Also: the legacy `aerindex` and `angstrm` are nighttime-masked
    at `radiation.F90:1341-1346` (same cosmetic-fill pattern as
    AODVIS). The underlying `aer_tau` IS computed for all `1:ncol`
    inside `modal_aero_sw` (it's pure aerosol optics — no solar-zenith
    dependence), so the values at idxnite columns are physically
    valid before the mask. New day+night companions `aerindexall` and
    `angstrmall` snapshot those values BEFORE the fill loop, parallel
    to the `AOD*all` pattern (added 2026-05-24; addfld in
    `output_aerocom_aie.F90:137-141`, outfld in
    `radiation.F90:1341-1348`). Both are gated on the same
    `do_aerocom_ind3=.true.` and `if(icall.eq.0)` (climatology run)
    conditions as the masked siblings, so all four fields appear
    together in eam.h0.

    **C. `cdnc` semantics.** From `output_aerocom_aie.F90`, `cdnc` is
    "Cloud droplet number concentration at cloud top" (#/m³), NOT a
    column or 3D field. It is populated inside `cloud_top_aerocom`
    (the AEROCOM cloud-top sampling routine) and is the canonical
    aerosol-cloud-interaction diagnostic for AEROCOM intercomparisons.
    For 3D droplet number on every level, use `AWNC` (P3 microphysics,
    in-cloud average) instead; for column-integrated, use `CDNUMC` (P3)
    or `cdnum` / `cdnumbf` (aerocom). The fincl1 entry intentionally
    uses the cloud-top `cdnc`.

    **D. CO2 (and friends) are already per-record in `eam.h0` —
    no FME-side code needed.** `cam_history.F90:4948-4953` writes
    `co2vmr(time)`, `ch4vmr(time)`, `n2ovmr(time)`, `f11vmr(time)`,
    `f12vmr(time)`, and `sol_tsi(time)` per record on every output
    interval, sourced from `chem_surfvals_co2_rad(vmr_in=.true.)` /
    `chem_surfvals_get('CH4VMR')` etc. These are 1D time-only doubles
    sitting in the file's coord vars (alongside `time`, `date`,
    `datesec`, etc.); xarray exposes them as `ds.co2vmr.values` of
    shape `(time,)`. Verified empirically in
    `v3.LR.piControl.aigo.eam.h0.0401-01.nc`: 123 records, all
    `2.84317e-4 mol/mol` for piControl 1850.

    **Why this matters**: for transient compsets (CMIP6 historical,
    scenario runs, 1pctCO2, 4xCO2), `flbc_inti` updates `co2vmr` as
    simulation time advances. The cam_history per-record write
    captures that drift automatically — every output record has the
    current CO2 forcing. No special FME plumbing needed.

    **Don't add `add_hist_scalar('co2vmr', ...)`** — an earlier
    revision of this branch tried that in `eam_vcoarsen_register`
    and would have produced a name collision with the existing
    `tape(t)%co2vmrid` (or silently emitted a constant frozen at
    `t=0`). Removed 2026-05-24 after auditing
    `cam_history.F90:4942-4953` and confirming `co2vmr(time)` is
    already in every h0 file we've ever written. The
    `add_hist_scalar` machinery added in gotcha #41 stays the right
    tool for *genuinely constant* per-file metadata (ak/bk, idepth)
    where one value per file is desirable; it's the wrong tool for
    time-varying scalars.

    **Files touched (5)**: `components/eam/src/physics/cam/modal_aer_opt.F90`
    (addfld + outfld for the 3 AOD*all variants),
    `components/eam/src/physics/cam/output_aerocom_aie.F90` (addfld
    for `aerindexall` / `angstrmall`),
    `components/eam/src/physics/rrtmg/radiation.F90` (outfld for
    `aerindexall` / `angstrmall` before the fill loop),
    `cime_config/testmods_dirs/allactive/fme_output/shell_commands`
    (`do_aerocom_ind3=.true.`, fincl1 extension, three more non-FME
    MPAS-O AMs disabled), this AGENTS.md entry.

    **Smoke-test checklist for the next live run**:
    - `AODVISall` shows non-fill values at nighttime grid cells
      (sanity-check global min/max excludes 1e36 outside polar night).
    - `aerindexall` (the fincl1-selected one) shows up only on RRTMG
      compsets (will be absent on `-rad rrtmgp`).
    - `cdnc.values.mean()` is in the O(1e7-1e8) #/m³ range over ocean
      stratus regions, ~0 over deserts.
    - `lwp.values.mean()` is ~O(0.1) kg/m² over ITCZ, ~0 over deserts.
    - `lcc` (liquid cloud cover fraction) is in [0,1]; `colccn.3` is
      a column CCN burden at S=0.3%.
    - `ccn.3bl` is O(1e7-1e9) #/m³ in polluted continental BL,
      O(1e6-1e7) over remote ocean.
    - `ds.co2vmr.values` is per-record (shape `(time,)`) and matches
      ~2.847e-4 for WCYCL1850 piControl 1850 (already in h0 files
      from cam_history.F90:4948 — not FME work).
    - `AODVIS` (legacy) is still nighttime-masked (fill at idxnite) —
      we did NOT break backward compatibility.

44. **MPAS FME remap streams: PIO "invalid IO type" abort at month
    rollover (RESOLVED 2026-05-24).** A continued (warm-restart) leg died
    exactly at the 1940-01 -> 1940-02 boundary with
    `PIO: FATAL ERROR ... invalid IO type (scorpio/.../pio_lists.c:120)`
    on a handful of OCN ranks, then thousands of `srun: ... Killed` lines.
    The kills are just MPI_ABORT teardown — not the cause. The tell is the
    new-month file: `fmeDerivedFields.1940-02.remapped.nc` left at **0
    bytes** (created, died in define mode) while its January file and the
    sibling `fmeDepthCoarsening`/`fmeSeaiceDerivedFields` Feb files were
    fine. `pio_lists.c:120` is `pio_get_file()` asserting
    `iotype_is_valid(cfile->iotype)`; since `pioc_support.c` always stores
    a valid iotype on a successful create, the assert means the looked-up
    ncid resolved to a **stale/orphaned** entry in `pio_file_list`.

    Both hand-rolled remap modules (`mpas_ocn_fme_horiz_remap.F` and
    `mpas_seaice_fme_horiz_remap.F` — these open/close/rotate PIO files
    directly instead of via the stream manager) had three defects:

    **A. Define-mode leak (load-bearing).** `open_file`'s define-error
    path `return`ed WITHOUT `pio_closefile` and without setting
    `file_is_open`, orphaning the just-created ncid; the next rotation
    overwrote `self%pio_file` with a fresh create, permanently leaking the
    old id. Once PIO recycled that id, the linear lookup returned the
    stale struct -> invalid iotype. Fix: close + `pio_file%fh = -1` +
    `file_is_open = .false.` on that path.

    **B. Close didn't reset the handle (defensive).** `close_file` now
    also zeroes `pio_file%fh` and the `dim_*_id` ids so a leftover handle
    can't alias a recycled ncid. Only matters once (A) has happened, but
    cheap insurance.

    **C. `is_restart_run` latched run-global (load-bearing).** It is set
    once at init from `config_do_restart` and never cleared, so EVERY
    mid-run month rotation took the append-reopen branch whenever the
    target file already existed — including a 0-byte orphan left by a
    prior crash, which `pio_openfile` then mishandles. Fix: a per-instance
    fire-once `restart_reopen_pending` latch (one for the primary stream,
    one for each `_5D` companion on the sea-ice side). Only the leg-start
    month appends; every genuinely-new month clobber-creates; and a failed
    reopen now falls through to clobber instead of aborting. This PRESERVES
    the leg-2-reopens-leg-1-leftover intent of gotcha #29 (the first open
    after restart still appends) — it just stops that privilege from
    leaking to later months.

    Edits land at both open/close sites per file (sea-ice has parallel
    primary + `_5D` state, gotcha #42), so six logical fixes total. Files:
    `components/mpas-ocean/src/shared/mpas_ocn_fme_horiz_remap.F`,
    `components/mpas-seaice/src/shared/mpas_seaice_fme_horiz_remap.F`.

    **Verification.** Syntax-only compile of each file by replaying the
    case's `flags.make` (`mpif90 -syntax-only -free -module <tmp>` with the
    `-D`/`-I` flags from
    `build/cmake-bld/mpas-framework/src/CMakeFiles/{ocn,ice}.dir/`) — both
    exit 0. Then a live restart run cleared the 1940-01 -> 02 rollover and
    reached 1940-02-05 with zero `invalid IO type`/`MPI_ABORT`; all six Feb
    `.remapped.nc` files (3 monthly + 3 `_5D`) populate with headers, incl.
    the `fmeDerivedFields.1940-02` that was 0 bytes before.

    **Smoke-test for the next rollover-crossing run**: confirm every
    `*.1940-MM.remapped.nc` is non-zero on the FIRST day of a new month,
    and that a restart leg whose start month already has a leftover file
    logs `clobber-creating instead` (the (C) fallback) rather than
    aborting.

45. **FME output schedule not reset on a branch/rewind restart + first
    post-gap record was a whole-gap average (FIXED 2026-05-28).** Two
    coupled defects surfaced when a run was rewound to an earlier date
    (e.g. `CONTINUE_RUN=FALSE` re-running from 0401 with output through
    0411 already on disk).

    **Bug 1 — schedule never re-fires.** `*_read_restart` restores
    `avg_last_output_time` (and the `_5D` companion) verbatim from the
    accumulator sidecar's `avg_last_output_time` global attribute
    (`ocn_fme_horiz_remap_read_accum`). The restore assumed a restart can
    only move time forward. On a rewind the restored value (0411-01-01)
    is in the FUTURE relative to `MPAS_NOW` (0401-...), so the output
    gate `currTime - avg_last_output_time >= avg_output_interval`
    (`mpas_ocn_fme_derived_fields.F` compute) is negative for the entire
    re-run and never fires — the re-crossed months get **no** remapped
    output; existing files are left byte-for-byte untouched. Output
    silently resumes only once the clock passes the stale time.
    Secondary: when output did resume, the pre-existing monthly file was
    append-reopened (gotcha #29/#44), and `write_static_masks_if_needed()`
    short-circuits on reopened files, so `mask_2d`/`history` kept the OLD
    map while the data records used the new map (the 0411-01 mixed-map
    file).

    **Bug 2 — contaminated first record.** Accumulation runs every
    timestep; the accumulator is reset only inside the flush block. While
    Bug 1 suppressed flushes, the restored sidecar sums kept growing
    across legs, so the first flush after the suppression lifted divided a
    ~10-model-year sum by its per-cell count and stamped it as a single
    1-day (or 5-day) mean. Independent latent defect: ANY mechanism that
    delays a flush past one interval (long stall, misconfig, rewind)
    produces a multi-interval first record.

    **Fix (both bugs, one mechanism).** In every AM `*_read_restart`,
    after restoring the schedule time, detect rewind with the strict test
    `avg_last_output_time > MPAS_NOW` (strict `>` leaves an exact
    same-point warm restart untouched, preserving the gotcha #29
    warm-restart BFB). On a hit: reset `avg_last_output_time = MPAS_NOW`,
    zero `remap_accum` / `remap_accum_count` / `nAccum` (and the `_5D`
    trio), and call the new `ocn_fme_horiz_remap_set_rewind()` /
    `seaice_fme_horiz_remap_set_rewind()`. That setter latches a
    module-scope `rewind_restart` flag that makes the file-open paths
    (`remap_file_open`, seaice `open_file` + `open_file_5D`) **clobber-
    create** instead of append-reopen (extra `.and. .not. rewind_restart`
    on the `may_append` predicate), so the re-crossed files get fresh
    records, a fresh `mask_2d`/`mask_*`, and fresh `history` for the
    current map. The flag is run-global and harmless once the clock
    catches up to the original timeline (the leg-start month is then in
    the future and doesn't yet exist, so it is created fresh anyway).

    **Applied to all FME AM schedules**: ocean `fmeDerivedFields`
    (1D + 5D), `fmeDepthCoarsening` (1D + 5D), `fmeVerticalReduce` (1D),
    and seaice `fmeSeaiceDerivedFields` (1D + 5D). The rewind detection
    lives in a small contained `rewind_reset` helper in the AMs that have
    a `_5D` companion (called once per schedule), inline in
    `fmeVerticalReduce`. Files touched (6):
    `mpas_ocn_fme_horiz_remap.F`, `mpas_seaice_fme_horiz_remap.F`,
    `mpas_ocn_fme_derived_fields.F`, `mpas_ocn_fme_depth_coarsening.F`,
    `mpas_ocn_fme_vertical_reduce.F`, `mpas_seaice_fme_derived_fields.F`.

    **Acceptance / smoke checks.** (1) ERS_Ld4 FME BFB still PASSes
    (strict `>` keeps same-point warm restart byte-identical). (2) A
    branch/rewind to an earlier epoch re-emits the re-crossed months with
    rerun mtimes, whose `mask_2d` zero-count and `history` match the
    current map. (3) The first resumed record is a plausible single-
    interval field — `mean|t1 - t2|` comparable to `mean|t2 - t3|`, not an
    anomalously smooth decadal mean. Log tell on a rewind leg:
    `rewind restart detected ... resetting schedule to MPAS_NOW` followed
    by `rewind restart flagged; re-crossed monthly files will be
    clobber-created`.

    **Optional further hardening (NOT done — flagged for the BFB-cautious).**
    The bug report also suggests guarding the steady-state flush path
    (discard/skip a record whose `currTime - avg_last_output_time`
    exceeds one interval by more than a tolerance). Skipped to avoid a
    tolerance-based false-discard in the hot compute path; the rewind
    reset already covers the only in-practice trigger (a stalled flush
    schedule). Revisit if a non-rewind long-gap scenario ever appears.

    **Forward-restart variant of the secondary defect — found in
    production, ALSO fixed 2026-05-28.** Scanning the live
    `v3.LR.piControl.aigo` run (90 yr, 1081 files/stream, forward
    20-yr-leg restarts only — no rewind) found NO Bug-1 suppression and NO
    Bug-2 contamination (output continuous, 0 duplicate dates; first-record
    spatial-std ratio ~1.0 at every leg boundary across all 6 streams; SST
    glides smoothly across each restart). BUT exactly one file —
    `fmeDerivedFields.0451-01` and its `fmeDepthCoarsening`/`*5D` siblings —
    shipped an **all-zero static mask** (`mask_2d`/`mask_0..18` = 64800
    zeros vs the correct 19170) while its data records were perfectly fine.
    Mechanism (confirmed via the per-leg `history` global attr): every
    20-yr leg ran slightly PAST the year boundary and created the `YY-01`
    file, then the next leg append-reopened it. At 0421/0441/0471 the
    creating leg reached the first flush (wrote the Jan-1 record AND the
    mask → 31 records, correct mask). At 0451 the creating leg defined the
    file but STOPPED BEFORE its first flush, leaving `mask_2d`
    defined-but-unwritten (= all-zero); the reopening leg then skipped the
    mask write (the `existing_file_reopened` early-return in
    `write_static_masks_if_needed`) → 30 records (starts Jan 2), all-zero
    mask. This is the same secondary defect as above but on a FORWARD
    restart, which the rewind-gated clobber does NOT cover.

    Fix: drop the `existing_file_reopened` early-return in all four mask
    writers (`write_static_masks_if_needed` + `_5D` in
    `mpas_ocn_fme_derived_fields.F` and `mpas_ocn_fme_depth_coarsening.F`),
    keeping the `time_record == 1` gate. A normal warm-restart reopen of a
    file that already holds records has `time_record >= 2` there, so its
    correct mask is left untouched (ERS warm-restart BFB preserved); only a
    reopen of a created-but-empty file (time_record reaches 1) rewrites the
    mask, with the correct value (`mask_*_remap` is recomputed each leg
    from the static land mask). Self-heals both the all-zero case and the
    bug report's stale-old-map case. **Existing 0451-01 files on disk are
    not retroactively fixed** — patch post-hoc by copying `mask_2d` /
    `mask_0..18` from any neighbor month (the masks are static), or the
    SamudrACE preprocessor will treat that month as all-land and drop it.

46. **Changing the remap MAP across a warm restart contaminates the
    straddling averaged record (found in production 2026-05-28; NOT yet
    fixed in code).** The `v3.LR.piControl.aigo` run changed the ocean
    remap map at the 0421 restart (old map → `map_IcoswISC30E3r5_to_
    gaussian_180by360_shifted_trintbilin.nc`; the target grid and land
    mask are identical, only the interpolation weights differ). The
    time-averaging accumulators (`remap_accum`, `remap_accum_5D`) live in
    REMAP (lat-lon) space (gotcha #11), and the sidecar (gotcha #29)
    carries them across the restart so leg-2 *continues* leg-1's partial
    averaging window. That continuation is BFB-correct when the map is
    unchanged (why ERS passes and the same-map restarts 0441/0451/0471 are
    bit-clean) — but when the map changes, leg-2 adds NEW-map remapped
    samples onto leg-1's OLD-map partial sums, so the first window that
    straddles the restart is an old+new-map blend.

    **Empirical signature (0421, all averaged streams; same-map restarts
    bit-exact by contrast):**
    - 1D: only the first full window computed across the restart is wrong
      (`fmeDerivedFields.0421-01` Jan-2 record: ocean-mean SST 13.259
      spikes off the smooth 13.190→13.197 trajectory, spatial step-RMS
      ~1.8 vs the ~0.10 baseline). The Jan-1 straddle record is the
      old-map record the prior leg already wrote (preserved, not
      overwritten — see gotcha #45 mechanism), so it is pure OLD map.
    - 5D: the first window `[Jan1,Jan6]` (dated 0421-01-06) differs from
      `mean(1D[Jan2..Jan6])` by 8.5% (SST) / 37% (SW flux); depth-coarsen
      temperature 8.3%; sea-ice area 27%. Every later window is bit-exact
      (rel ~2e-8). Diagnostic: a correct 5-day mean equals the mean of its
      five overlapping daily means to float32 roundoff — see the
      `verify_5d` cross-check in this incident's scripts.

    So at a map-change restart the boundary month holds: 0401-0420 = pure
    OLD map; `YY-01` rec0 (the Dec31→Jan1 straddle) = OLD map; the first
    full new-leg window (1D Jan-2; 5D Jan1-Jan6) = OLD+NEW blend; and
    everything after = pure NEW map. The bulk of the data IS overwritten
    with the new map — only the one straddling window per stream is
    contaminated, plus the single old-map straddle record.

    **Immediate workaround (operational, no code change):** when changing
    the map mid-run, DELETE the `*.fme_accum_restart.nc` sidecars before
    resuming. The averagers then start fresh on the new map; the first
    window after the restart is a clean (if slightly short) pure-new-map
    partial mean instead of a blend.

    **Proper fix (analogous to the gotcha #45 rewind reset, not yet
    implemented):** stamp the map-file identity (filename or n_s weight
    count) into each sidecar as a global attr; in `*_read_restart`, if the
    restored map identity differs from the active map, DISCARD the
    accumulator (zero `remap_accum`/`_count`/`nAccum`, set
    `avg_last_output_time = MPAS_NOW`) exactly like the rewind path —
    rather than blindly continuing an accumulator that is meaningless in
    the new remap space. This makes a deliberate mid-run map change safe.

    **Note on the 5D averager itself:** it is CORRECT in general — mid-leg
    and at all same-map restarts, every 5D record equals the mean of its
    five overlapping daily records to float32 roundoff. The only 5D errors
    found are the map-change straddle windows above, and they are not
    5D-specific (the 1D straddle record is wrong the same way, same cause).

47. **Negative-weight remap maps (TempestRemap 2nd-order FV / "trfv2")
    silently corrupt the masked remap path (found in production
    2026-06-02; FIXED with a self-guarding `apply_masked` + init refusal).**
    The `v3.LR.piControl.aigo` production tape had non-physical extremes in
    the ocean depth-coarsened fields: `temperatureCoarsened_10` up to
    250 degC, `velocityMeridionalCoarsened` to -13 m/s, and the smoking gun,
    `layerThicknessCoarsened` down to **-4.8e5 m** (negative thickness is
    impossible — `accumulate_thickness` only sums `max(0,...)` overlaps).
    ~2400 columns/file had a negative coarsened thickness, 96% of them at the
    column's DEEPEST valid level (the rest within 1-2 levels); ~40 of those
    blew temperature out of `[-3,40]`. Static, identical in `_inst` and
    `avg`, so NOT a time-averaging artifact; SST / level-0 were clean.

    **Root cause.** The FME remap always uses
    `shr_horiz_remap_apply_masked`, which normalizes each target cell by
    `valid_frac = sum_j w_j * mask_j` over only the VALID (in-coverage)
    sources. With a NON-NEGATIVE map this is a convex combination → output
    stays within the source range. The run's map was
    `..._shifted_trfv2.nc`, a TempestRemap 2nd-order FV map that is
    **52.8% negative weights** (min weight -3.28, max +4.60). At the surface
    every source is valid so `valid_frac` = full row sum = 1.0 and the
    negatives balance as designed (surface clean). At a column's deepest
    valid level only a few sources remain valid; their signed weights no
    longer sum to ~1 and can cancel to ~0 or go NEGATIVE, so
    `numerator / valid_frac` explodes / inverts sign → negative thickness,
    250 degC. The `SHR_COVERAGE_EPS=1e-6` floor only catches
    `valid_frac ~ 0`, not small-positive or negative cancellation.

    **Timing — it is the MAP, not the restart.** Corruption spanned
    0401 → 0411-01 (CFS run) and 0401 → 0421-01 (a second run on
    `$PSCRATCH/E3SMv3/...`). The 5-yr restarts WITHIN those spans
    (0406, 0411, 0416, 0421-01) did NOT clear it — proof a restart alone
    doesn't fix it. Each run cleared only when the operator swapped the map
    at a restart: CFS → first-order `conservative` (44892 target cells),
    the pscratch run → `trintbilin` (45630, same footprint), both
    all-positive. Don't mistake the coincidence ("cleared after the
    restart") for causation; the *map swap* (which requires a restart to
    take effect) is what fixed it. Verified exhaustively month-by-month:
    zero negative thickness anywhere after each run's swap. Map weight signs:
    `conservative` and `trintbilin` have 0 negatives; `trfv2` has 53%.
    (The `..._patch.nc` / higher-order maps also have negatives — never use
    them on the masked path.)

    **Fix (`share/util/shr_horiz_remap_mod.F90`):**
    1. `read_mapfile` scans the global weight list and sets
       `has_negative_weights` / `min_weight` on `shr_horiz_remap_t`.
    2. `apply_masked` accumulates `valid_frac_abs = sum |w_j|*mask_j`
       alongside the signed `valid_frac`, and only normalizes when
       `valid_frac >= WELL_COND_FRAC * valid_frac_abs` (0.5), else emits
       `SHR_FILL_VALUE`. For a non-negative map `valid_frac ==
       valid_frac_abs` so this is **BFB** (ERS stays clean); for a
       negative-weight map it fills the cancelling cells and bounds any
       residual output to `<= max|source| / 0.5`.
    3. `mpas_ocn_fme_horiz_remap.F` init checks `has_negative_weights` and
       aborts with `MPAS_LOG_CRIT` naming the map and min weight — so a
       trfv2-type map can never SILENTLY ship a corrupt tape again.
    The shared-module guard protects EAM/sea-ice masked callers too; only
    the MPAS-ocean wrapper got the loud init abort (add to the seaice/EAM
    wrappers if they ever take a user-supplied map).

    **Operational:** the corrupt legs (CFS 0401-0411-01; pscratch
    0401-0421-01) must be masked or re-run; everything after each swap is
    clean. The only post-swap out-of-generous-range values are real model
    features, NOT this bug: SSS ~50-54 PSU at hypersaline basins
    (Persian Gulf 24-25N/50.5E, Gulf of Suez 28-29N/32-33E, Shark Bay
    26S/114.5E) and sea-ice `surfaceTemperatureMean` ~ -50 degC in polar
    winter — both identical before and after the swap.

48. **Code-review follow-up — latent-path correctness + verify-script /
    submission fixes (APPLIED 2026-06-19; full findings in `REVIEW.md` in this
    dir).** A branch-wide re-review found no CRITICAL defects in
    production-exercised paths, but several latent bugs one namelist edit from
    triggering, verify scripts drifted from the production config, and two
    submission-script issues. Fixes applied:

    **EAM latent Fortran (production-unused today, but namelist-reachable):**
    - `physics_buffer.F90.in` `pbuf_get_field_ndims` counted padded trailing
      dims (always returned 6, so eam_derived misrouted 2D pbuf operands into
      OOB reads). Now counts dims with size > 1 (column dim always present).
    - `eam_derived.F90` passed the whole `tmp_fields(pcols,pver,*)` to
      `shr_derived_eval`; for a 2D expression (nlev=1) sequence association
      mapped operands 2..N onto operand 1's upper levels (gotcha #2 class),
      silently evaluating them as zero. Now passes `(:,1:nlev,1:n_operands)`.
    - `eam_vcoarsen.F90` `int_prev_valid` was a module-scope SCALAR used inside
      the OMP chunk loop (race + wrong first-step `d*_INT_dt`). Now a per-chunk
      `(begchunk:endchunk)` array. **This corrects gotcha #16:** the
      scalar-flag race was eliminated in eam_derived but NOT in eam_vcoarsen
      until now.
    - `*at2m`/`*at10m` long_names now disclose that the height select collapses
      to the lowest model level (~10-25 m AGL) when the target is below the
      lowest midpoint.

    **Shared remap (`shr_horiz_remap_mod.F90`):** `build_comm` now counts source
    gcols with no valid owner rank and returns nonzero `ierr` instead of
    silently dropping them (which left `col_recvidx=0` → OOB read at apply with
    `ierr=0`). All three wrappers already abort on a nonzero build_comm ierr.

    **MPAS remap def_var replay:** ocean + seaice (primary and `_5D`)
    `check_rotate` replay loops now fold each `def_var` ierr into `err` and log
    loudly, so a failed (re)definition can't leave a stale descriptor that
    `write_var` then writes through silently. (No slot-name clearing: the
    registered-var and per-file var lists share arrays — the
    `seaice_fme_remap_file_t` refactor in "Extensions" would let `write_var`
    skip dead slots cleanly; deferred.)

    **CIME test integrity (`config_archive.xml`) — corrects gotcha #34.** The
    1D and 5D FME streams collapsed into one `_get_extension` ext-bucket (the
    unanchored 1D regex `fmeDepthCoarsening` also matches `…Coarsening5D…`
    files), so ERS silently restart-BFB-compared only one of {1D,5D} per AM and
    base/rest could pick different ones (flaky NO_COMPARE). Added `…5D`
    `<hist_file_ext_regex>` entries BEFORE the 1D ones and anchored the 1D
    regexes with a trailing `\.`. **Re-create existing test cases** so their
    CASEROOT `env_archive.xml` picks up the new spec (per gotcha #34's note).

    **Verify scripts (now usable on the production config):**
    - `verify_mpas.py` certification made `fmeVerticalReduce` OPTIONAL (disabled
      in production, gotcha #33) so the cert can PASS; it previously hard-required
      that stream and always FAILed + `exit 1`.
    - `verify_mpas.py` ice-presence fill mask is now rebuilt PER RECORD from each
      record's `iceAreaTotal` (was a single first-file snapshot → false FAILs on
      seasonal ice).
    - `verify_eam.py` `cross_verify` SKIP paths now return the 3-tuple shape
      (were bare lists → `ValueError` unpack crash when h0 missing/header-only).
    - `verify_eam.py` Test 1b TOPO_FILL now expects 0% fill in level-index
      coarsening mode (production) via a `VCOARSEN_MODE` constant, instead of a
      pressure-bound expectation that always DIFFed.
    - `verify_production.py` `main()` now aggregates FAIL sanity entries (and
      missing PRIMARY streams; a missing 5D companion is tolerated) into a
      nonzero exit code — was implicit `exit 0` on the normal path (the
      exception handler already exited 1).

    **Submission script (`run_e3sm.fme.sh`):** dropped `-cosp` (no COSP fields
    on the FME tape → 100 yr of wasted simulator runs) and set the default
    `WALLTIME` to the 24 h pm-cpu regular-QOS cap (was 34 h, which Slurm
    rejects).

    **Other doc corrections from the review (no code change needed):** gotcha
    #47's loud negative-weight init abort is NOT ocean-only — the seaice
    (`mpas_seaice_fme_horiz_remap.F:224`) and EAM (`horiz_remap_mod.F90:110`)
    wrappers already have it. Gotcha #9's `nCoarsenLevels > nFmeDepthLevels`
    guard fires at first compute, not init. Gotcha #22 now covers 23
    fill-threshold sites, not 22. `FME_MAPS_DIR` defaults to `$DIN_LOC_ROOT/fme`,
    not NERSC scratch (corrected in the two sections above).

    **NOT changed (deferred, see REVIEW.md):** O2 (no `fmeVerticalReduceOutput`
    buildnml block — only matters if that AM is re-enabled; the stream query
    falls back); gotcha #46 proper map-identity sidecar discard; the seaice
    type-encapsulation refactor; the remaining LOW/robustness items.

## Runtime Configuration

Both `fme_output` and `fme_legacy_output` testmods accept environment variables:

```bash
FME_EAM_OUTPUT_HOURS=6          # EAM output frequency in hours (default: 6)
FME_EAM_AVGFLAG=I               # 'I' = instantaneous (default); fluxes use per-field :A
FME_EAM_STORAGE=one_month       # File rotation policy (default: one_month).
                                # Values: one_month, one_year, num_snapshots.
                                # See gotcha #27.
FME_EAM_MFILT=1500              # Safety bound when calendar rotation is in
                                # charge (default: 1500). Force daily
                                # behavior with FME_EAM_STORAGE=num_snapshots
                                # FME_EAM_MFILT=4.
FME_MPAS_INTERVAL=00-00-01_00:00:00      # MPAS output/averaging interval (default: daily)
FME_MPAS_INTERVAL_5D=00-00-05_00:00:00   # 5-day-mean companion stream cadence
                                         # (set to 'none' to disable; see #42)
FME_MAPS_DIR=/path/to/fme_maps  # SCRIP map directory (default: $DIN_LOC_ROOT/fme)
```

Example: `FME_EAM_OUTPUT_HOURS=24 FME_MPAS_INTERVAL=00-00-05_00:00:00 ./create_test ...`

MPAS FME AMs support time averaging via namelist:
- `config_AM_fme*_time_averaging = .true.` enables accumulation
- `config_AM_fme*_compute_interval = 'dt'` samples every timestep
- `config_AM_fme*_stream_output_interval` sets the averaging window
- The `fme_output` testmod enables averaging by default
- Two ocean AMs enabled in production: fmeDepthCoarsening, fmeDerivedFields
  (fmeVerticalReduce was disabled 2026-04-30 -- see gotcha #33)
- One sea-ice AM enabled: fmeSeaiceDerivedFields
- All four AMs' native streams are gated to `output_interval="none"` --
  only the `.remapped.nc` files are written for the SamudrACE tape (#32)

Verification scripts support cross-comparison:
```bash
# EAM verification
python verify_eam.py --rundir $RUNDIR --outdir /path/to/figs

# MPAS verification with legacy cross-comparison
python verify_mpas.py --rundir $RUNDIR --outdir /path/to/figs \
    --legacy-rundir $LEGACY_RUNDIR

# Production-tape validator: 6 streams (3 1D + 3 5D), CF sanity,
# fill-leak scan, 5D-vs-mean(1D) cross-check (#42), HTML viewer.
# Push to a shared web-viewable area with --share-dir; the curated
# top-level $PSCRATCH/share/index.html is updated additively (a new
# <li> is inserted before the closing </ul> if no link to the case
# already exists; existing entries are NEVER touched).
python verify_production.py --rundir $RUNDIR \
    --outdir /tmp/fme-verify \
    --share-dir $PSCRATCH/share \
    --label v3.LR.piControl.aigo.test5
```

## Verification Dashboard

`verify_eam.py` and `verify_mpas.py` produce HTML+figure dashboards. Each
runs three classes of check:

1. **Per-stream sanity:** file inventory, variable presence, range bounds,
   global stats. Generates summary tables and last-timestep maps.
2. **Cross-stream physics:** radiation budget closure (EAM), heat-flux
   sum closure (MPAS-O), per-layer bathymetry consistency
   (depth-coarsening), temporal monotonicity.
3. **Cross-testmod comparison:** FME (online remap) vs legacy
   (timeSeriesStatsCustom + offline pipeline) global means, side-by-side
   maps, vcoarsen linearity. Requires both `--rundir` and `--legacy-rundir`.

Each script ends with an "ACE Training Readiness Certification" section
that gates production. The certification passes when:
- All expected ACE variables are present
- No NaN/fill in cells expected to be ocean (per-layer bathymetry-aware
  for depth-coarsened fields; surface SST mask for surface fields;
  ice-only fields exempt)
- Physical ranges respected
- Radiation/heat-flux budgets close to within ~5 W/m^2

### Production-readiness assessment (as of 2026-04-30)

After the 2026-04-30 verify-script edits (gotchas #19-#26) plus the
SamudrACE-spec edits later that day (Q vcoarsened to Q_0..Q_7,
Qat2m, TREFHT, QREFHT, PRECST=PRECSC+PRECSL), both dashboards
converge cleanly. EAM PASSes 7/7 with all 22 cross-compare checks
green. STW linearity (Test 3) was the last remaining EAM "failure"
and resolved when the FME global-mean was switched to cos-lat
weighting (#26). MPAS PASSes after the ice-presence mask (#25),
bumped fill threshold (#22), and ICE_FME_ONLY exemption (#23).

**Production file aggregation** (gotchas #27, #28): EAM and MPAS
both produce monthly files (~6000 total over 100 yr instead of
~182,500 daily). EAM via the new `hist_file_storage_type='one_month'`
namelist; MPAS via Registry `filename_template="$Y-$M.nc"`. The
production run script is `run_e3sm.fme.sh` at the repo root with
5-yr restart cadence aligned to the monthly file boundary --
mitigates the file-clobber-on-restart concern (gotcha #29).

**EAM topographic-fill question -- RESOLVED.** The earlier concern
about lower-troposphere vcoarsen layers (5-7) showing zero fill
where pressure-range coarsening would put cells below terrain no
longer applies: production switched to **level-index-based**
coarsening (`vcoarsen_level_bounds = 0,25,38,46,52,56,61,69,80` in
`shell_commands`, defining 8 layers by EAM hybrid level indices),
not pressure ranges. Every column has all 80 model levels by
construction, so every output layer always has valid samples on
every column. No below-terrain cells, no fill expected, no
inconsistency between Fortran and verify. The pressure-bounded
code path in `eam_vcoarsen.F90` is still available for users who
prefer it; production just doesn't use it.

**One open MPAS question (likely physical):**
- Latent/sensible heat flux extremes at Gulf Stream / Kuroshio winter
  outbreaks reach -1000 W/m^2. Range bounds were loosened from
  ±500 W/m^2 to (-1200, 500) W/m^2 for LHFLX and (-800, 500) W/m^2 for
  SHFLX based on reanalysis observations. Spot-check FME vs legacy
  values cell-by-cell at the extreme locations to confirm physical
  (not remap artifact) before declaring done.

## Remaining Work

### Open spec questions to resolve before the 100-yr submit

Cross-reference of `shell_commands` (testmod) vs the SamudrACE
Confluence spec ("case run options for v3 LR piControl SamudrACE",
p3ai/6154289880, last edited 2026-04-25) was done 2026-05-01.
Items needing explicit confirmation from the SamudrACE team before
submit (one resolved 2026-05-05; two remain):

1. **STW formula -- RESOLVED 2026-05-05.** Spec wording "RIME" was
   loose; the actual desired field is `Q + CLDLIQ + CLDICE + RAINQM`
   (i.e., total water including rain, NOT rime mass on snow).
   Confirmed with the SamudrACE team. Our `derived_fld_defs =
   'STW=Q+CLDICE+CLDLIQ+RAINQM'` in `shell_commands` is correct as
   written; no change needed.

2. **`uVelocityGeoCell` (cell-centered) vs spec `uVelocityGeo`
   (vertex-located).** The existing `geographicalVectors` AM
   defines `uVelocityGeo` on `nVertices`, derived from the B-grid
   primary solve at vertices. Our `uVelocityGeoCell` is on
   `nCells`, derived by rotating the pre-aggregated
   `uVelocityCell/vVelocityCell` (one interpolation step removed
   from the solver). For SamudrACE training on a regular lat-lon
   grid, cell-centered is almost certainly what's wanted (it
   matches the cell-based remap framework and all other Samudra
   fields are cell-centered), but the spec literal is ambiguous.
   Confirm with team. If they want vertex-located: substantial
   implementation lift (new vertex-mesh map file + vertex-aware
   remap path; the FME framework only does cell-based remap
   currently).

3. **Sea-ice aggregate naming: spec `iceAreaCell/Volume/snowVolumeCell`
   vs ours `iceAreaTotal/Volume/snowVolumeTotal`.** The data is
   bit-identical (same `sum_iCategory` formula, verified in
   `mpas_seaice_icepack.F:4450` vs `mpas_seaice_fme_derived_fields.F:541`)
   — only the output name differs. Trivial to rename our outputs
   to match spec. Verify this is just a naming-convention
   preference, then rename in `mpas_seaice_fme_derived_fields.F`
   (register_var calls + local target arrays + verify_mpas.py
   expected names).

### Near-term
- **Monthly-boundary smoke test for the 5D companion stream
  (RESOLVED 2026-05-10).** Verified empirically on
  `v3.LR.piControl.aigo.ppl_test6` `PM_2x5_nyears` running from
  `0401-01-01`. All seven streams (EAM h0 + 3 MPAS 1D + 3 MPAS 5D)
  rotated cleanly across both the Jan→Feb and Feb→Mar boundaries.
  Observed counts and the rule that drives them — month-N's file
  holds records whose *closing timestamp* falls in month N, so the
  window straddling the month boundary lands in month N+1:
    - **EAM h0 (6h inst): Jan=123, Feb=112.** Jan loses the t=0
      initial state (first sample is at Jan 1 06:00). Feb is the
      full 28×4=112 because the Feb 1 00:00 sample has
      `date=4010201` → Feb file.
    - **MPAS 1D: Jan=30, Feb=28.** The `[day 30, day 31]` window
      (Jan 31 → Feb 1) end-stamps Feb 1 → lands in Feb. The
      `[day 58, day 59]` window (Feb 28 → Mar 1) end-stamps Mar 1
      → lands in Mar. So month-N's daily file holds windows
      *closing* in month N, which is one window short of "all
      windows whose data is in month N" (the last day's window
      is in month N+1).
    - **MPAS 5D: Jan=6, Feb=5.** Boundary-straddling first-Feb
      record: `time=35`, `time_bnds=[30, 35]`
      (Jan 31 00:00 → Feb 5 00:00), `date=4010205`.
  Continuity Δ between Jan-last and Feb-first time stamps: 0.25
  (EAM 6h), 1.0 (1D), 5.0 (5D) — exact, no gaps, no overlaps.
  `verify_production.py` cross-check deferred until production
  submit.
- Add CI test variant (SMS_Ld2, ne4pg2_oQU480 for fast builds)
- Add SMS_Ld40 monthly-rotation smoke test (verify Jan->Feb file
  rotation, record counts per month match calendar 124/112/124/...,
  no duplicate timestamps across consecutive files)
- PEM MPI reproducibility test
- ERS_Ld4 PASSes COMPARE_base_rest BFB-clean as of 2026-05-01 with
  the append-mode + accumulator-sidecar + frame-tracking +
  compute_on_startup=.false. fix bundle (see gotchas #11 and #29).
  Verified independently across two case directories: 766 fields in
  fmeDepthCoarsening / 98 in fmeDerivedFields / 42 in
  fmeSeaiceDerivedFields all IDENTICAL. EAM h0 and cpl.hi remain
  IDENTICAL. The non-FME AM blockers (timeSeriesStatsDaily restart
  bug; vertex-velocity diffs) are sidestepped by disabling
  globalStats / regionalStatistics / timeSeriesStats* in the
  testmod — those AMs aren't part of the SamudrACE training tape.
- A follow-up ERS at production cadence (REST_OPTION='nyears',
  monthly file rotation) should run cleanly given the above fixes
  exercise the same code paths -- still worth scheduling for
  paranoia.

### Verify-script hardening (open after 2026-04-30 punch list pass)
These are robustness improvements that don't block production but
would catch latent issues during the 100-yr run:
- Calendar-aligned file inventory: assert all EAM/MPAS files match
  the `*.YYYY-MM.nc` pattern; per-month record count matches
  noleap calendar (Jan=124, Feb=112, ...); no duplicate xtime
  values across adjacent files.
- Time-axis-on-grid check (assert `time mod nominal_interval ≈ 0` now
  that `avg_last_output_time + avg_output_interval` advances exactly,
  per gotcha #11).
- `_inst` companion variable presence/range check (production sets
  `config_AM_fmeDepthCoarsening_write_instantaneous_companion=.true.`
  but verify never reads `*_inst` fields). Suggested: at minimum, a
  one-file presence/range check in the depth-coarsening readiness path.
- Coverage-eps coastline parity: compare fill mask of
  `temperatureCoarsened_0` (apply_masked path) vs `sst` (apply path);
  both should now agree to within ~one cell after the unified
  `SHR_COVERAGE_EPS=1e-6` (gotcha #17). Currently no test asserts this.

### Code quality
- 3D PIO decomposition for depth-coarsened fields (single 3D var instead of
  per-level 2D vars like `temperatureCoarsened_0..24`)
- Depth coordinate variable in remapped 3D output
- SSH masking under ice shelf cavities (cavities have maxLevelCell > 0 and
  can have SSH < -1000m; not caught by iceFraction correction)
- Add lonCell/latCell to MPAS-O FME Registry stream definitions for native
  grid diagnostic plots

### Extensions
- **Coupler-native FME output (merged forcing).** Design approved, not yet
  implemented — full brief in `## PLAN: Coupler-Native FME Output` at the
  end of this file. Emits `x2o`/`x2i`/`xao` (the post-merge forcing each
  component integrates), which no component tape exposes.
- Online spherical harmonics roundtrip (replaces ACE offline SH smoothing)
- ELM / MOSART FME integration (deferred -- needs reimplementation)
- Performance: persistent work arrays, O(n) send_local_cell lookup, OpenMP
- **Write per-cell sample counts to MPAS FME outputs** so the
  5D-vs-mean(1D) cross-check (#42) can be exact rather than
  p99-tolerant. The AM already maintains `remap_accum_count` /
  `remap_accum_count_5D`; just emit them as a `count_<var>` companion
  to each averaged var. Then verify_production can compute
  `weighted_mean_1D = sum(count_i * v_i) / sum(count_i)` and compare
  bit-by-bit. Cost: doubles the var count in the file (~2x storage
  overhead — not free for 100-yr SamudrACE tape) so probably only
  worth it as a sanity-only stream variant, not production.
- **Refactor `mpas_seaice_fme_horiz_remap.F` to a type-encapsulated
  `seaice_fme_remap_file_t`** (mirroring `ocn_fme_remap_file_t`).
  Erases the duplicated `_5D` state vars and procedures added in
  #42 — those landed via parallel state because a clean refactor
  was deemed too risky under time pressure. ~1200 lines to touch
  but net code reduction once done.
- **`verify_production.py --quick` mode** that re-renders the HTML +
  re-runs sanity / cross-check WITHOUT regenerating the ~580 PNGs
  (timeseries + climatology per var). Plotting takes the bulk of
  the runtime; iterating on the cross-check / sanity logic doesn't
  need fresh figures.

## PLAN: Coupler-Native FME Output — As-Exchanged Merged Forcing (driver-mct)

**Status:** design approved 2026-06-20. **IMPLEMENTED 2026-06-20/06-21 — SEVEN
bundles. Ocean-mesh bundles each 1D + 5D (10 streams): imports x2o/x2i, exports
o2x/i2x, plus coupler xao. Atm-grid bundles x2a(in)/a2x(out) each a SINGLE
configurable sub-daily stream (`cpl_fme_atm_interval`, default 6 h) = 2 streams.
12 streams total.** Purely additive: the EAM h0 tape and the MPAS-O/SI FME AMs
stay exactly as-is. All file:line refs were verified against *this* checkout's
`driver-mct` before coding.

> ### ⇨ HANDOFF — read this first (state as of 2026-06-21 PM, 4th update)
> **Where it stands: ✅ ACCEPTANCE MET — coupler-native FME is BFB-validated
> end-to-end on ERS_Ld4 with all 12 streams on.** Fixed: `cime_final` segfault
> (gotcha #49); archive bug that compared only a sidecar (gotcha #50); the
> merged-import warm-restart BFB bug (gotcha #51). **Progression of ERS_Ld4
> COMPARE_base_rest runs:**
> - run `...185103_coawyr` (pre-#50): "PASS" but only 1 sidecar compared — false
>   confidence.
> - run `...143731_xn3ikq` (post-#50, pre-#51): all 7 outputs compared → 5
>   IDENTICAL, **2 FAILED (`x2a`,`x2i`)** — diagnosed gotcha #51.
> - run `...153319_od1ku8` (post-#51): **COMPARE_base_rest PASS, ALL 7 coupler
>   outputs IDENTICAL** (x2o/x2i/o2x/i2x/x2a/a2x/xao). x2a came back clean too,
>   so no residual upstream-driver issue. Full test green (RUN, COMPARE, MEMLEAK,
>   SHORT_TERM_ARCHIVER).
>
> **Remaining work is NOT correctness — it's housekeeping:**
> 1. **Commit** the validated tree (currently the `tmp` commit `e15fe601a3` +
>    uncommitted rename/archive/#51 changes). Suggested split: (a) coupler module
>    `cpl_fme_mod.F90` + driver hooks + namelist + verify_cpl.py; (b) archive fix
>    + ocn/ice→x2o/x2i rename; (c) the #51 phase-hook fix. (User is handling the
>    commit.)
> 2. **5D outputs are still unexercised** — a 4-day ERS never flushes a 5-day
>    window, so only `.accum.nc` sidecars exist for the 5D streams. Run a longer
>    test (≥ ERS_Ld11, or production cadence) to BFB-validate the 5D output
>    flush + its append-mode reopen.
> 3. Optional: add `cpl_fme_mod.F90` + `verify_cpl.py` to the "Source Files"
>    section of this doc.
>
> **Streams (renamed 2026-06-21, gotcha #50 "Stream rename"):** ocean/ice
> IMPORTS are `x2o`/`x2i` (were `ocn`/`ice`). The 6 component bundles are
> symmetric by direction: `x2o`/`o2x`, `x2i`/`i2x`, `x2a`/`a2x` (+ coupler
> `xao`). Files `<case>.cpl.fme.<stream>{1D,5D}.YYYY-MM.nc`, stream ∈
> {x2o,x2i,o2x,i2x,xao} (1D+5D) and {x2a,a2x} (single 6 h). Per-bundle enable
> flags `cpl_fme_<stream>_enable`/`_5d_enable`; SHARED ocean-mesh map/window vars
> stay `cpl_fme_ocn_{map,interval,5d_interval}` on purpose.
>
> **To re-validate / re-run from scratch:** `cd cime/scripts && ./create_test
> ERS_Ld4_P1024.ne30pg2_r05_IcoswISC30E3r5.WCYCL1850.pm-cpu_gnu --testmod
> allactive/fme_output` — MUST be a fresh case-create (the `config_archive.xml`
> spec only reaches a case via its CASEROOT `env_archive.xml` at create time, not
> a rebuild). Acceptance = COMPARE_base_rest PASS with all 7 `cpl.fme.*1D`/`x2a`/
> `a2x` outputs IDENTICAL. If `x2a`/`x2i` ever regress, check only the
> post-restart records (day 3.25 for x2a, day 4 for x2i) — that points back at
> the merged-import sampling phase (gotcha #51).
>
> **Config (testmod defaults, all `.true.`):** all 12 streams,
> x2a/a2x at 6 h, the ocean/ice bundles at 1D + 5D. Atm maps confirmed present:
> `$DIN_LOC_ROOT/fme/map_{ne30pg2,IcoswISC30E3r5}_to_gaussian_180by360_shifted_trintbilin.nc`.
>
> **Known-unfixed:** gotcha #46 (map-change-across-restart accumulator
> contamination) applies to the coupler sidecars too in principle, but the
> coupler accumulates in SOURCE space (map-independent) so it is **immune** —
> see the "KEY DESIGN DECISION" note below.

**Coupler data model (why these names):** the MCT mediator stores, per
component, an IMPORT (`x2c`: x2a/x2o/x2i — what the component integrates) and an
EXPORT (`c2x`: a2x/o2x/i2x — what it produces). There is NO `a2o`/`a2i` bundle;
the atm→ocean piece is the transient `a2x_ox` consumed by the `x2o` merge. So
the FME coverage is per-component in/out: atm = x2a(in)+a2x(out), ocn =
x2o(in)+o2x(out), ice = x2i(in)+i2x(out), plus xao (coupler bulk fluxes). x2a
AND a2x ride the ATM grid (ne30pg2 → `cpl_fme_atm_map` = the EAM map); all
others ride the ocean mesh (`cpl_fme_ocn_map`; sea-ice shares it). Note a2x is
ALSO available on the EAM h0 tape, but the coupler-native a2x is emitted on the
same lat-lon grid / cadence as x2a so SamudrACE gets atm forcing both IN and OUT
in one self-consistent place.

**Average vs. flux (matches the coupler) — IMPLEMENTED.** Every field in every
bundle is accumulated (running source-space SUM over the output window) and
divided by `nAccum` at flush → a true time-mean, tagged `cell_methods='time:
mean'`. This mirrors the coupler's own treatment: x2o is sampled from the
`x2oacc` accumulator the coupler itself sums-then-averages; the other bundles
have no native accumulator so we average coupling-cadence snapshots. Crucially,
**no field is summed instead of averaged** because the coupler exchanges flux
DENSITIES (rates: W m⁻², kg m⁻² s⁻¹), not per-step accumulated totals — so the
mean rate is the correct aggregation for both state and flux variables, exactly
as the coupler/EAM-h0 do. There is no quantity in these bundles that is a true
per-step running total, so uniform time-mean is right across the board.

**Atm cadence (x2a/a2x) — configurable.** `cpl_fme_atm_interval` (days, default
0.25 = 6 h) sets the atm window; sampling happens once per ocean coupling step
(the `cpl_fme_accum` hook), so the interval must be ≥ and a multiple of the
coupling step. For default WCYCL (ATM_NCPL=OCN_NCPL=ICE_NCPL=48 → 30-min step)
a 6-h window averages 12 exact snapshots. If a config ever couples the ocean
coarser than the atm, x2a/a2x would undersample — only then is a second hook in
`cime_run_atmocn_fluxes` warranted (the still-live "second hook" open question).

### Implementation notes (2026-06-20)
New module `driver-mct/main/cpl_fme_mod.F90` + 5 hooks in `cime_comp_mod.F90`
(use:133; `cpl_fme_init` at the end of `init_readrestart` ~2438;
`cpl_fme_accum` right after `prep_ocn_accum_avg` ~4014; `cpl_fme_restart_write`
after `seq_rest_write` ~5155; `cpl_fme_final` in `cime_final` ~3605).
**12 config-driven streams across 7 bundles:** ocean-mesh `x2o{1D,5D}` (ocean
merged import), `xao{1D,5D}` (Faox_*), `x2i{1D,5D}` (sea-ice import, ice gsmap —
sea-ice shares the ocean mesh so the SAME map file), `o2x{1D,5D}` (ocean
export), `i2x{1D,5D}` (ice export); atm-grid `x2a` (merged import) and `a2x`
(atm export), each a single
stream at `cpl_fme_atm_interval` (default 6 h). x2o is sampled from the
cpl-restart-safe `x2oacc` window-mean (PLAN enabler #1); all other bundles have
no native accumulator so we average coupling-cadence snapshots (sampled at the
ocean hook). `cpl_fme_init` takes `(infodata, EClock, ocn, ocn_present, ice,
ice_present, atm, atm_present, read_restart)` — atm args are REQUIRED (they feed
x2a/a2x). Config: `cpl_fme_inparm` group in `namelist_definition_drv.xml`
(auto-written to `drv_in` by the existing nmlgen — no buildnml change; ocn map +
atm map + ocn 1D/5D windows + atm window + per-bundle enable flags), `FME_CPL_*`
env + `user_nl_cpl` block in the testmod `shell_commands`, per-stream
`<hist_file_extension>fme\.<stream>\.\d.*\.nc$</hist_file_extension>` (dated
outputs only, NOT the `.accum.nc` sidecars — gotcha #50) plus ONE
`<hist_file_ext_regex>fme\.\w+</hist_file_ext_regex>` for the per-stream COMPARE
bucket, for `{x2o,x2i,xao,o2x,i2x}{1D,5D}` plus `x2a` and `a2x`, in
`driver-mct/cime_config/config_archive.xml` so ERS actually COMPAREs the new
files (gotchas #34/#48/#50 — **re-create test cases** to pick up the spec).
Readiness check: `verify_cpl.py` (CF time/grid/field-presence/fill-leak +
5D-vs-mean(1D) cross-check for the ocean-mesh bundles; xgns→share flow).

**KEY DESIGN DECISION — accumulate in SOURCE space, not remap space.** The
super-accumulator is an `mct_aVect` of the merged forcing on the ocean-cell
decomposition (built/averaged with the same MCT ops as `x2oacc_ox` itself);
fields are remapped via `shr_horiz_remap` *only at flush*. For a linear masked
remap with a static land mask this is bit-equivalent to remap-then-accumulate,
but it (a) reuses the proven `x2oacc` persistence path for restart-BFB trust,
(b) is **immune to gotcha #46** (map-change contamination — the accumulator is
map-independent), and (c) is cheaper (1 remap/field/window vs per coupling
step). Consequence: the per-cell-count machinery of gotcha #36 is unnecessary
here (coupler bundles carry no time-varying source fill; the only fill is the
static land mask introduced by `apply_masked` at the spatial level). The
restart sidecar is `<case>.cpl.fme.<stream>.accum.nc` (per-field source-space
sums + `nAccum` + drift-free schedule + epoch), written on the coupler restart
alarm and read at init; rewind handled per gotcha #45 (strict `>` keeps warm-
restart BFB; future schedule ⇒ reset + clobber-create). Honors #11 (advance
`next_output_time` by exactly one interval — no drift), #29/#44 (first-open
append + strict-`<` frame seeding to overwrite the straddling leftover record),
#1 (real(4)/PIO_REAL), #3/#5 (apply_masked + SHR_FILL_VALUE), #47 (init refuses
negative-weight maps). Output files: `<case>.cpl.fme.x2o1D.YYYY-MM.nc` and
`...x2o5D...` (renamed from `ocn1D`/`ocn5D` 2026-06-21 — see gotcha #50), CF
time = END of window (matches EAM/#39), reusing the existing
`map_IcoswISC30E3r5_to_gaussian_180by360_shifted*.nc` (no new map files).

**Acceptance bar: MET (2026-06-21).** Fresh ERS_Ld4 `...153319_od1ku8` with all
fixes + 12 streams PASSED `COMPARE_base_rest` with every per-stream 1D output
(`cpl.fme.{x2o,x2i,o2x,i2x,xao}1D.*` + `x2a`/`a2x`) IDENTICAL — see the HANDOFF
box at the top of this PLAN. **Remaining housekeeping:** (1) commit the validated
tree (user is doing this); (2) a longer run (≥ ERS_Ld11 / production cadence) to
BFB-exercise the 5D output flush, which a 4-day ERS never reaches. The default
field lists are the authoritative sets in `cpl_fme_mod.F90`
(`x2o_default`/`xao_default`/`x2i_default`/`o2x_default`/`i2x_default`/
`x2a_default`/`a2x_default`); absent fields are skipped+warned at setup for
compset safety. **Source files for the FME testmod doc list:** add
`cpl_fme_mod.F90` (coupler) and `verify_cpl.py` (this dir) when next updating
the "Source Files" section.

### gotcha #49 — coupler FME segfault in `cime_final`: free PIO decomp BEFORE close flushes the buffered darray (FIXED 2026-06-21)
First ERS_Ld4 (built from the 3-bundle ocn/xao/ice tree) ran the entire leg —
wrote the restart at day 4, all `fme.*` output + `.accum.nc` sidecars on disk —
then **SIGSEGV at the very end** on ~all cpl ranks. The GNU backtrace was
misleading (line numbers skewed by `-O2` inlining; the SIGSEGV itself was a
*secondary* fault during Kokkos `atexit` teardown after an `MPI_Abort`).
`addr2line` on the executable gave the true chain:
`pio_closefile → PIOc_closefile → sync_file → flush_buffer →
PIOc_write_darray_multi → pio_err → piodie`.
**Root cause:** `pio_write_darray` (in `cpl_fme_file_write_var`) BUFFERS field
writes and PIO flushes them inside `pio_closefile`. `cpl_fme_final` was freeing
`iodesc_dst` (`pio_freedecomp`) *before* `cpl_fme_file_close`, so the close-time
flush dereferenced a freed decomposition → PIO aborts. (The sidecar path didn't
crash because it frees via the *file* overload `pio_freedecomp(pf, iod)` AFTER
its own write, and closes after — order was fine there.) **Fix:** in
`cpl_fme_final`, **close the file FIRST** (flush buffered darrays with the decomp
still valid), **then** `pio_freedecomp`. One-spot reorder; canonical PIO
lifecycle. NB: the run also exposed a compile-blocking typo — `cpl_fme_init`
declared `atm`/`atm_present` but omitted them from the SUBROUTINE arg list (the
3-bundle tree that built predates the atm additions); fixed alongside.
**Lesson:** never free a PIO decomposition while any buffered (unflushed)
`write_darray` still references it; close/sync first. Diagnose PIO segfaults
with `addr2line -f -e bld/e3sm.exe <addr>`, not the inlined Fortran line nums.

### gotcha #50 — coupler FME files weren't actually BFB-compared; sidecars must NOT be archived (FIXED 2026-06-21)
The first passing ERS_Ld4 (WCYCL1850, COMPARE_base_rest PASS) hid a gap: of all
`cpl.fme.*` files, cprnc compared exactly ONE — `cpl.fme.xao5D.accum.nc`, a
restart *sidecar*, not even a data output. Two distinct CIME-archive bugs, both
in the `drv` `comp_archive_spec` (`driver-mct/cime_config/config_archive.xml`):

1. **Bucket collapse (the gotcha #34 mechanism, again).** `config_archive.xml`
   had 12 per-stream `<hist_file_extension>` entries but NO
   `<hist_file_ext_regex>`. `archive_base._get_extension` uses the latter to
   key the "latest file per ext" COMPARE dedup; with none, it falls back to
   `\w+` and collapses EVERY `cpl.fme.*` file into one `"fme"` bucket → only the
   alphabetically-last survives (`xao5D.accum.nc`). **Fix:** add one
   `<hist_file_ext_regex>fme\.\w+</hist_file_ext_regex>` → each stream gets its
   own bucket (`fme.x2o1D`, `fme.x2a`, ...). Verified by simulating
   `get_all_hist_files`+`_get_extension` on the real rundir: 3/3 enabled 1D
   outputs now bucketed separately.

   **Stream rename (2026-06-21, same pass):** the ocean/ice IMPORT streams were
   renamed `ocn`→`x2o` and `ice`→`x2i` so all six component bundles read
   symmetrically by direction: `x2o`/`o2x`, `x2i`/`i2x`, `x2a`/`a2x` (+ `xao`).
   Files are now `<case>.cpl.fme.{x2o,x2i}{1D,5D}.YYYY-MM.nc`. The rename touched,
   in lockstep: the `add_stream` name strings + per-bundle enable flags
   `cpl_fme_{x2o,x2i}_enable`/`_5d_enable` in `cpl_fme_mod.F90`; the matching
   `<entry id=...>` in `namelist_definition_drv.xml`; the `user_nl_cpl` heredoc
   in `shell_commands`; the `config_archive.xml` globs + tfiles; and
   `verify_cpl.py` (`EXPECTED` keys + `STREAMS`). The SHARED ocean-mesh
   map/window namelist vars stay named `cpl_fme_ocn_{map,interval,5d_interval}`
   on purpose — they are mesh-level (used by x2o/o2x/x2i/i2x/xao), not
   bundle-specific.

2. **Sidecars were being swept in as hist (latent PRODUCTION restart bug).**
   The cpl accumulator sidecars are named `<case>.cpl.fme.<stream>.accum.nc` —
   they carry the casename prefix, so unlike the MPAS `*.fme_accum_restart.nc`
   sidecars (which dodge archiving by NOT starting with the casename) they
   matched the bare `fme\.<stream>\..*\.nc$` glob. With `disposition="move"`,
   st_archive would MOVE them out of the run dir at end-of-leg → the next
   RESUBMIT leg's `cpl_fme_read_sidecar` (reads from the run dir at init) finds
   nothing → cold-starts the accumulator (WARN) → breaks warm-restart BFB across
   the archive boundary. Didn't surface in ERS (no post-archive leg). **Fix:**
   tighten each glob to dated OUTPUTS only — `fme\.<stream>\.\d.*\.nc$` (the
   `\.\d` requires the YYYY-MM date right after the stream name, excluding
   `.accum.nc`). Sidecars now stay in the run dir for restart and are simply not
   compared — which is fine: their BFB correctness is already exercised by the
   dated 1D output comparison (same accumulator + append-mode-reopen path).
   **Consequence:** a short ERS (<5 d) compares only the 1D outputs; the 5D
   outputs need a test long enough to flush a 5-day window.

**Acceptance bar still requires a case RE-CREATE** (not just rebuild): the
`env_archive.xml` baked into an existing CASEROOT is a stale copy of this spec.
The 2026-06-20/21 case was created from the 3-bundle tree, so its
`env_archive.xml` predates even the o2x/i2x/x2a/a2x streams.

### gotcha #51 — merged-import streams (x2i, x2a) not warm-restart BFB: sample each import AFTER its own merge, not at the ocean hook (FIXED + CONFIRMED 2026-06-21)
With gotcha #50's archive fix the coupler OUTPUTS finally got cprnc-compared on
a fresh case (`...20260621_143731_xn3ikq`), and ERS_Ld4 COMPARE_base_rest
**FAILED on exactly two streams — `x2a` and `x2i`** — while `x2o`, `o2x`,
`i2x`, `xao`, `a2x` were all IDENTICAL. The two failures are precisely the
merged IMPORTS (`component_get_x2c_cx`).

**Diagnosis (the arithmetic nails it).** ERS restarts at end of day 3, so only
post-restart records can differ:
- `x2i` (1-day window, 48 coupling steps/window): its single leg-2 record
  (day 4) is off by mean ~5.87 K on `Sa_tbot` (Sa_tbot≈282 K). **282/48 =
  5.875 K** = exactly ONE zero-valued sample dragging down a 48-sample mean.
- `x2a` (6-h window, 12 steps/window): differs ONLY on the first post-restart
  record (day 3.25); later records identical. `Sx_tref` normalized RMS ≈ 8.7%
  ≈ **1/12** = one zero sample in 12. Only the dynamic ocean/ice-derived fields
  (`Sx_*`,`So_*`,`Faxx_*`) differ; near-static `Sf_*` fractions don't (a single
  zero barely moves a constant). For `x2i` the WHOLE bundle is stale so even
  `Sa_*`/`Faxa_*` differ.

**Root cause.** The single accum hook was at the OCEAN phase (right after
`prep_ocn_accum_avg`). In the driver step order: `prep_atm_mrg` (merges x2a) →
ocean hook → `prep_ice_mrg` (merges x2i). So at that hook **x2i is sampled one
merge too early** — it still holds the PREVIOUS step's value, which is a
re-initialized ZERO on the first step of a warm-restart leg. x2a is merged
before the hook but from an ocean→atm contribution that is stale/zero on the
first restart step. Either way the first post-restart sample is ~0 and pollutes
that window's mean by `field/nAccum`. The other five streams are immune: `x2o`
reads the restart-persisted `x2oacc` window-mean; `xao` is recomputed each step;
`o2x`/`i2x`/`a2x` are component EXPORTS straight from the (BFB) models.

**Why skipping the bad sample does NOT work.** Dropping the first post-restart
sample would leave the restart window with `nAccum`-1 samples vs the base
window's `nAccum` → different mean → still not BFB. The only correct fix is to
sample the import at a point where it is validly merged in BOTH runs.

**Fix.** Sample each merged import right after its own merge — the same point
`component_diag(flow='x2c')` reads it for the (BFB) cpl history. `cpl_fme_accum`
gained a `phase` argument (`CPL_FME_PHASE_{OCN,ICE,ATM}`); each stream carries a
`phase` (set in `add_stream`: x2i→ICE, x2a→ATM, all else→OCN) and is processed
only at its phase. `cime_comp_mod.F90` now calls the hook three times per step:
OCN after `prep_ocn_accum_avg`, ATM right after `prep_atm_mrg`, ICE right after
`prep_ice_mrg`. The driver clock advances once per step (`seq_timemgr_clockAdvance`,
well before any prep), so all three phases see the same `curr_time` and the
drift-free flush schedule (gotcha #11) stays consistent. Files: `cpl_fme_mod.F90`
(+phase param/field, `add_stream` select, `cpl_fme_accum(EClock, phase)` filter)
and `cime_comp_mod.F90` (use-list + 3 hook calls). Minimal/surgical — only x2i
and x2a moved; the 5 already-BFB streams stay at the ocean hook.

**STATUS: CONFIRMED.** Fresh ERS_Ld4 `...153319_od1ku8` (all fixes #49/#50/#51,
12 streams on) PASSED COMPARE_base_rest with **all 7 coupler outputs IDENTICAL**
— `x2a` and `x2i` now bit-for-bit alongside the 5 that already passed. `x2a`
came back clean, so the feared upstream-driver residue (ocean→atm export not
restored before the first restart-step atm merge) did NOT materialize: sampling
right after `prep_atm_mrg` is sufficient. Full test green (RUN, COMPARE_base_rest,
MEMLEAK, SHORT_TERM_ARCHIVER).

### Why (the one quantity component tapes cannot provide)
The EAM h0 tape and the MPAS FME AMs emit each component's *internal*
fields. They do NOT emit the **merged, post-coupler forcing each component
actually integrates** — that exists only on the mediator:
- `x2o` — merged forcing the ocean imports (atm precip/radiation +
  coupler-computed turbulent fluxes `Foxx_*` + ice melt `Fioi_*`).
- `x2i` — merged forcing sea-ice imports.
- `xao` — the coupler's own atm–ocean bulk fluxes (`Faox_taux/tauy/lat/
  sen/lwup/evap/swdn/swup`), recomputed on the exchange grid every atm
  step — distinct from EAM's surface fluxes.
For a *coupled* ACE/Samudra emulator the training target is precisely this
exchanged forcing. That is the SOLE justification for a coupler path.
**`a2x` is NOT a reason** — EAM already emits it (the `zbot..swnet` block
added to `cam_diagnostics::diag_export` 2026-06-20 mirrors `atm_export`
bit-for-bit). So skip the 6h atm coupler stream unless you want
config consolidation; it duplicates EAM h0.

### Key enabler #1 — the merged ocean forcing is already accumulated AND restart-safe
`x2oacc_ox` (`driver-mct/main/prep_ocn_mod.F90:104-105,247-252`) is a running
SUM of merged `x2o` over the ocean coupling window:
- accumulated every atm step in `prep_ocn_accum` (prep_ocn_mod.F90:394-423),
- averaged in place (÷count) in `prep_ocn_accum_avg` (prep_ocn_mod.F90:427-455),
  which copies the window-mean into `x2o_ox = component_get_x2c_cx(ocn)`
  and resets the counter (`x2oacc_ox_cnt = 0`),
- **persisted to the cpl restart** (`x2oacc_ox` + `x2oacc_ox_cnt`:
  `seq_rest_mod.F90:271-272` read / `607-609` write) — so the window-mean
  base is BFB across a mid-window stop with ZERO new code.
- exposed through **existing public getters** `prep_ocn_get_x2oacc_ox()` /
  `prep_ocn_get_x2oacc_ox_cnt()` (prep_ocn_mod.F90:50-51, 1481-1489) — the
  same ones `seq_rest_mod` uses (:54-55).
⇒ The "as-exchanged ocean forcing" is `x2o_ox` immediately after
`prep_ocn_accum_avg`. Cleanest tap: `use prep_ocn_mod, only:
prep_ocn_get_x2oacc_ox` and read the accumulator pointer directly (it holds
the in-place-averaged window-mean after `prep_ocn_accum_avg`) — this
decouples `cpl_fme` from call-ordering instead of relying on `x2o_ox`
still holding the mean. There is **no** native
accumulator for ice (`x2i`) or atm (`x2a`) — those are sampled at
coupling cadence. `xao` (`Faox_*`) is recomputed every atm step in
`seq_flux_atmocn_mct` (`driver-mct/main/seq_flux_mct.F90:1316`; exchange-
grid variant `seq_flux_atmocnexch_mct` at :960 — confirm which is active).

### Key enabler #2 — no new map files
Coupler bundles live on decompositions we already have maps for (see
`$DIN_LOC_ROOT/fme/`):
- atm-grid (a2x/x2a) → `map_ne30pg2_to_gaussian_180by360_shifted.nc` (n_a=21600)
- ocn/ice-grid (x2o/o2x/x2i/i2x/xao; **sea-ice shares the ocean mesh**) →
  `map_IcoswISC30E3r5_to_gaussian_180by360_shifted.nc` (n_a=465044)
The shared remap `share/util/shr_horiz_remap_mod.F90` is component-agnostic
— pack `x2o`/`x2i` as a flat array sized to the gsmap exactly the way
`components/eam/src/control/horiz_remap_mod.F90` packs CAM chunks.

### Design
New module `driver-mct/main/cpl_fme_mod.F90`, mirroring the MPAS FME AM
pattern (init → accumulate → remap → direct-PIO monthly write + sidecar
restart). N independent streams, each with own alarm / field list / map /
accumulator. `seq_hist_mod.F90:seq_hist_writeaux` (1060-1365) is the
structural template for multi-stream/multi-cadence/optional-averaging —
but route the averaged AV through `shr_horiz_remap`, not raw `seq_io`.

Streams (priority order):
1. **ocn-merged-1D** — `x2o` (the prize) + `o2x` ocean state, daily mean, IcoswISC30E3r5 map.
2. **ocn/ice-merged-5D** — `x2o` + `x2i` + `xao`, 5-day mean, IcoswISC30E3r5 map.
3. *(optional, redundant — default OFF)* **atm-6h** — `a2x`, ne30pg2 map.

Accumulation:
- ocean streams: accumulate the window-mean `x2o_ox` (BFB base) into a
  1D / 5D super-accumulator. ONLY this super-accumulator needs a sidecar.
- ice/xao: accumulate `x2i_ix` / `xao_ox` at coupling cadence (fluxes ⇒
  average, not snapshot).

### Integration points (`driver-mct/main/cime_comp_mod.F90`)
- `cpl_fme_init` — near existing FME/AM + histaux setup, after gsmaps/doms.
- `cpl_fme_accum(...)` — call right after `prep_ocn_accum_avg` (so
  `x2o_ox` holds the window-mean) and after the ice merge each ocn/ice
  coupling step. Co-locate with the histaux call block (~4912-5078).
- `cpl_fme_write(EClock)` — on each stream's alarm; remap + monthly PIO write.
- `cpl_fme_restart_write/read` — co-locate with `cime_run_write_restart`
  (~5101-5148) and the `seq_rest` read.

### Config
Mirror the `FME_EAM_*`/`FME_OCN_*` pattern: add `FME_CPL_*` env vars in the
testmod `shell_commands` → translate to a `cpl_fme` namelist group in
`drv_in` via the driver `cime_config/buildnml` +
`driver-mct/cime_config/namelist_definition_drv.xml` (the `histaux_*` defs
at 1211-1346 are the template: per-stream enable, cadence option/N, field
list, map file, storage type).

### Packaging
Monthly files with sub-window snapshots, via the MPAS `$Y-$M` filename
rotation + `check_rotate` / `register_var` replay (gotcha #6) — direct
PIO, NOT the CAM history path.

### Restart / BFB — reuse the proven patterns, do not reinvent
- Port the MPAS accumulator sidecar
  (`mpas_ocn_fme_horiz_remap.F` write_accum 1527-1643 / read_accum
  1660-1787) → `*.cpl.fme_accum_restart.nc` for the 1D/5D super-accumulators.
- Honor gotcha **#11** (advance `avg_last_output_time + interval` exactly —
  no drift), **#29** (append-mode reopen + frame tracking), **#36**
  (per-cell valid-sample counts so fill never enters the mean), **#45**
  (rewind-restart: clobber-create + strict `>` to stay BFB).
- The base `x2oacc` is already cpl-restart-safe — do NOT re-accumulate it;
  only the 1D/5D layer is ours.
- Acceptance bar: ERS_Ld4 PASSes COMPARE_base_rest with the cpl FME streams
  ON (same bar as the existing EAM/MPAS path).

### Remap details (reuse)
- Output `real(4)` with `PIO_REAL` decomp (gotcha #1).
- `apply_masked` with level-numlev+1 mask packing (gotcha #3) and
  `SHR_FILL_VALUE` flow (gotcha #5).

### Verify
Extend `verify_production.py` with a coupler-stream readiness check (CF
time/bnds, fill-leak scan, range sanity, 5D-vs-mean(1D) cross-check #42 vs
the daily stream), or a parallel `verify_cpl.py`; add to the xgns run block.

### Field selection (authoritative lists in `driver-mct/shr/seq_flds_mod.F90`)
- **x2o** merged forcing: `Foxx_taux/tauy/sen/lat/lwup/evap`, `Foxx_swnet`,
  `Faxa_lwdn`, `Faxa_rain/snow`, `Fioi_melth/meltw/salt`, `Sa_pslv`,
  `Si_ifrac` (x2o assembled at :3934; fluxes 1307-1682).
- **xao**: `Faox_taux/tauy/lat/sen/lwup/evap/swdn/swup`,
  `So_tref/qref/u10/ustar/duu10n` (:1176-1491, 2642-2649).
- **x2i**: atm state `Sa_*` + `Faxa_*` + ocean state `So_*` + `Fioo_q`
  (:627-1744).
Start with `x2o` (the exact ocean-integrated forcing); add `xao`/`x2i`
once the path is proven.

### Task checklist
1. `cpl_fme_mod.F90` skeleton: types, init (read maps via shr_horiz_remap), per-stream state.
2. Accum hook after `prep_ocn_accum_avg` + ice merge; super-accumulator + per-cell counts.
3. Remap + direct-PIO monthly write (port from `mpas_ocn_fme_horiz_remap.F`).
4. Sidecar restart read/write; wire into `seq_rest` / `cime_run_write_restart`.
5. `namelist_definition_drv.xml` + driver buildnml `cpl_fme` group; `shell_commands` `FME_CPL_*` block.
6. ERS_Ld4 COMPARE_base_rest BFB validation.
7. `verify_production.py` (or `verify_cpl.py`) coupler-stream check.

### Open questions
- Ocean/ice coupling cadence vs the 1D/5D windows: confirm the
  window-mean→day super-accumulation count math (reuse #11 drift-free advance).
- `x2i` has no native accumulator: snapshot-at-cadence vs our-own-average
  (fluxes ⇒ average).
- Whether to ship the redundant 6h `a2x` stream (default: NO).
- Which `seq_flux_atmocn*` variant is active in WCYCL (decides `xao` cadence).
