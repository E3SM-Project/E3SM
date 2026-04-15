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
- `analysis_members/mpas_seaice_fme_derived_fields.F` -- 7 derived fields
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

Map files (NERSC): `/pscratch/sd/m/mahf708/fme_maps/`
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
    stream output remains instantaneous (last snapshot).

12. **Remapped output files include xtime.** The `xtime(StrLen, Time)` character
    variable records the MPAS date-time string for each output record. CF
    attributes include `axis='X'`/`'Y'` on lon/lat and `calendar='noleap'`
    as a global attribute. PIO `pio_put_var` for xtime requires a 1-element
    character array (`character(len=64) :: buf(1)`), not a scalar — PIO's
    `put_vara_1d_text` interface demands `character(len=*) :: val(:)`.

13. **MPAS coupler fluxes are ice-fraction-weighted.** The forcing pool fields
    (`shortWaveHeatFlux`, `longWaveHeatFluxDown`, `latentHeatFlux`,
    `sensibleHeatFlux`, `windStressZonal`, `windStressMeridional`) represent
    the flux received by the ocean surface, NOT the atmospheric flux.
    Under heavy ice coverage (iceFraction >= 0.99), these values are near
    zero and not physically meaningful. The FME derived fields AM masks them
    to `fillValue` before remapping. `surfaceHeatFluxTotal` (the sum of
    coupler fluxes) is likewise masked. SST/SSS remain valid under ice.

14. **EAM averaging convention.** State variables (T, U, V, STW, PS, TS)
    are instantaneous snapshots (`avgflag_pertape='I'`). Fluxes and
    tendencies (FSDS, FLDS, FLUT, FSUTOA, FLUS, FSUS, LHFLX, SHFLX,
    TAUX, TAUY, PRECT, SOLIN) use per-field `:A` suffix for time
    averaging over the output interval. EAM shortwave radiation fields
    (FSUTOA, FSUS) show sharp gradients at ice edges — this is real
    surface albedo contrast (ice ~55% vs ocean ~3%), not an artifact.

## Runtime Configuration

Both `fme_output` and `fme_legacy_output` testmods accept environment variables:

```bash
FME_EAM_OUTPUT_HOURS=6          # EAM output frequency in hours (default: 6)
FME_EAM_AVGFLAG=I               # 'I' = instantaneous (default); fluxes use per-field :A
FME_EAM_MFILT=4                 # samples per output file (default: 4)
FME_MPAS_INTERVAL=00-00-01_00:00:00  # MPAS output/averaging interval (default: daily)
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
- Add FSUTOA, FSNTOA to legacy testmod fincl1 for native-vs-interpolated comparison

### Code quality
- 3D PIO decomposition for depth-coarsened fields (single 3D var instead of
  per-level 2D vars like `temperatureCoarsened_0..24`)
- Depth coordinate variable in remapped 3D output
- SSH masking under ice shelf cavities (currently only masked by iceFraction;
  cavities have maxLevelCell > 0 and can have SSH < -1000m)

### Extensions
- Numeric `time` coordinate + `time_bnds` for CF-compliant time axis
- Online spherical harmonics roundtrip (replaces ACE offline SH smoothing)
- 2D/3D wetmask computation
- ELM / MOSART FME integration (deferred -- needs reimplementation)
- Performance: persistent work arrays, O(n) send_local_cell lookup, OpenMP
