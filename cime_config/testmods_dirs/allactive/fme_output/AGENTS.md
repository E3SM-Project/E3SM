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

31. **`*_at_z2` / `*_at_z10` height-interp fields restart bug — FIXED
    and verified 2026-04-30 (ERS_Ld4 PASSes COMPARE_base_rest, all 1478
    EAM fields IDENTICAL incl. T_at_z2, Q_at_z2, STW_at_z2, U_at_z10,
    V_at_z10).** Root cause was gotcha #2 (Fortran sequence association
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

    **Workaround if fix regresses**: drop T_at_z2, Q_at_z2, STW_at_z2,
    U_at_z10, V_at_z10 from `fincl1` in `shell_commands`. SamudrACE
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
FME_MPAS_INTERVAL=00-00-01_00:00:00  # MPAS output/averaging interval (default: daily)
FME_MAPS_DIR=/path/to/fme_maps  # SCRIP map directory (default: NERSC mahf708 scratch)
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
Q_at_z2, TREFHT, QREFHT, PRECST=PRECSC+PRECSL), both dashboards
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
- Online spherical harmonics roundtrip (replaces ACE offline SH smoothing)
- ELM / MOSART FME integration (deferred -- needs reimplementation)
- Performance: persistent work arrays, O(n) send_local_cell lookup, OpenMP
