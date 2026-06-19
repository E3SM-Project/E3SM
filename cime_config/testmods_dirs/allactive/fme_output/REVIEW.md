# FME Branch Review — `claude/fme-reconciled` vs `upstream/maint-3.0`

- **Date:** 2026-06-12
- **Branch tip reviewed:** `0a38867bb29307b00ba015182a2e7b610c47361c`
- **Base (upstream/maint-3.0):** `6dfcad8d7a6984ddd42717784a2de4df4948c415`
- **Diff size:** 65 files, ~25,300 insertions / 53 deletions, 64 commits
- **Method:** six parallel subsystem reviews (shared utils, EAM, MPAS-O, MPAS-SI,
  Python verifiers, config/scripts/MOSART), all AGENTS.md gotchas treated as
  hypotheses and re-verified against code. Every CRITICAL/HIGH finding below was
  independently re-confirmed by reading the cited code a second time before
  inclusion. Line numbers refer to the branch tip above.

---

## Executive summary

**No CRITICAL defects were found in the production-exercised paths.** The
restart/rotation/rewind/accumulator state machines in both MPAS remap modules
check out against the documented invariants (gotchas #29/#36/#39/#44/#45/#47
all verified implemented as described, with the exceptions noted below). The
feared sea-ice `_5D` copy-paste divergence does **not** exist — verified by
normalized diff *and* targeted identifier sweeps in both directions.

The real exposure falls in four buckets:

1. **Latent Fortran bugs in shipped, namelist-reachable but
   production-unused paths** — three HIGH findings in EAM (`pbuf` ndims
   misclassification, a surviving sequence-association bug in `eam_derived`
   2D eval, a module-scope-scalar OMP race in `eam_vcoarsen` column
   integration) and one HIGH in the shared remap module (silently dropped
   unowned source columns → out-of-bounds read).
2. **Verification scripts have drifted from the production config** — the
   MPAS readiness certification can **never PASS** the current testmod
   (hard-requires the disabled `fmeVerticalReduce` stream), `verify_eam`'s
   cross-verify **crashes** on its own SKIP path, and `verify_production`
   **always exits 0** even on FAIL.
3. **One test-infrastructure gap that silently weakens ERS** — the
   `config_archive.xml` ext-regexes predate the 5D streams; 1D and 5D files
   collapse into one comparison bucket, so one of them is never
   restart-BFB-checked (and base/rest can pick different ones → flaky
   NO_COMPARE).
4. **Production-submission script issues** — `-cosp` wasted compute,
   34 h walltime vs the 24 h pm-cpu QOS cap, dead PLACEHOLDER guard.

**Must-fix before the 100-yr submit:** V1 (cert dead-on-arrival), C1
(config_archive 5D buckets + case re-creation), V3/V4 (verify crash / exit-0),
S1 (`-cosp`), S2 (walltime). **Should-fix soon:** E1–E3, U1 (latent
correctness bugs one namelist edit away from triggering).

Severity scale: CRITICAL (wrong data in production now) · HIGH (wrong
data/crash on a reachable path, or a test that cannot do its job) · MEDIUM
(robustness gap with data-affecting consequence) · LOW (hardening) · NOTE
(record-keeping).

---

## 1. Shared utilities (`share/util/`)

### U1. HIGH — `build_comm` silently drops needed source columns with no valid owner → out-of-bounds read with `ierr=0`
`share/util/shr_horiz_remap_mod.F90:336-342, 355-362, 399`

A needed gcol whose `gcol_to_rank` entry is outside `[0,nprocs)` (the MPAS
wrappers initialize `gcol_to_rank(:) = -1` and fill only owned cells —
`mpas_ocn_fme_horiz_remap.F:352`, `mpas_seaice_fme_horiz_remap.F:300`) is
silently excluded from the request lists. Its `gcol_to_recvpos` stays 0, so
`col_recvidx = 0`, and at apply time
`ws_recv_buf((src_idx-1)*numlev + k)` with `src_idx = 0` reads index
`k - numlev ≤ 0` — silent garbage or segfault, `ierr = 0` throughout. The
wrappers' `send_local_cell == 0` fatal check (gotcha #10) only covers the
*send* side and cannot catch this: an unowned gcol is never requested from
anyone, so it never appears in any rank's send list. Trigger: any map whose
`col` entries reference cells not owned by any rank (map/mesh mismatch — the
exact scenario the init checks claim to be fatal on). **Fix:** count dropped
gcols in build_comm step 1 and return nonzero `ierr`.

### U2. MEDIUM — map-file `col` values unvalidated → out-of-bounds write at init
`shr_horiz_remap_mod.F90:197-199`. `gcol_marker(rd%tmp_src_gid(i)) = 1` with
`gcol_marker` sized `n_a` and `tmp_src_gid` read straight from the file's
`col` variable. A corrupt or 0-based map writes out of bounds (heap
corruption, no error). `row` is implicitly safe (filtered to
`[row_start,row_end]`). Also: `n_s == 0` yields an all-fill output with
`ierr=0` instead of a loud refusal.

### U3. MEDIUM — 32-bit overflow in target partition arithmetic at scale
`shr_horiz_remap_mod.F90:159-160`. `row_start = myrank * n_b / nprocs + 1`
in default integer. For `n_b = 64800` this overflows at ~33k ranks; for a
0.25° target (~1.04e6) at ~2k ranks. Silent garbage partitioning. Fix:
compute the product in `integer(i8)`.

### U4. MEDIUM-LOW — double-init early return leaves `send_gcol_list` unallocated with `ierr=0`
`shr_horiz_remap_mod.F90:330-331`. The `intent(out), allocatable` dummy is
auto-deallocated on entry before `if (rd%initialized) return`. Currently
unreachable (all three wrappers guard with their own `is_initialized` flag),
but it is a trap for any future caller.

### U5. LOW — assorted hardening gaps
- Failed init is not re-entrant (allocate-on-allocated on retry); no destroy
  routine (`:171, :289`).
- `pio_get_var` / `pio_inq_dimlen` return codes unchecked (`:239-281`) —
  truncated map file yields uninitialized data consumed as weights.
- lon-fastest destination-grid ordering assumed, unvalidated (`:289-293`) —
  lat-major map produces silently wrong coordinate metadata.
- `shr_vcoarsen_select_nearest` reads OOB if `nlev_max(i) > nlev`
  (`shr_vcoarsen_mod.F90:318-324`); `select_index` validates, this path
  doesn't.
- `shr_vcoarsen_cat_sum` has no fill handling (`:383-388`) and `cat_wgtavg`
  doesn't skip fill field values with positive weight (`:420-424`) — safe
  today only by caller convention (sea-ice categories are fill-in-all or
  0-in-empty).
- `shr_derived_mod` exponent heuristic swallows `+`/`-` operators after field
  names ending digit+[eEdD] (`X=Q1D+X` parses as one bogus field name —
  fails loud downstream, but with a mystifying message) (`:135-147`).
  Silent truncation at the 256-char rhs / 34-char token buffers, no ierr
  (`:114, :117, :164`). Stale header comment says divide-by-zero "produces
  0.0" — contradicts the implemented fill behavior (`:223`).

### Verified correct (traced, not assumed)
- **Gotcha #47 BFB claim is exact**: for non-negative maps `valid_frac` and
  `valid_frac_abs` accumulate identical values in identical order, so the
  `WELL_COND_FRAC=0.5` guard degenerates bit-for-bit to the pre-guard
  `> 1e-6` test. Negative-weight detection is rank-consistent (every rank
  scans the full global `S`), and all three wrappers check
  `has_negative_weights` at init.
- Gotcha #3 (mask packed as level numlev+1, can never leak into output),
  #17 (unified `SHR_COVERAGE_EPS`), #18 (div-by-zero → 1e20), #8
  (increasing-coordinate precondition degrades to all-fill, not garbage).
- MPI alltoall/alltoallv count-displacement math, ordering invariants, CRS
  bucket sort, zero-cell-rank buffer sizing, and file-close-on-all-paths in
  `read_pio_contents` all traced correct.
- The unmasked `apply` (frac_b path) has **zero callers** on this branch —
  dead at runtime; keep or delete deliberately.

---

## 2. EAM (atmosphere)

### E1. HIGH — `pbuf_get_field_ndims` always returns 6; 2D-pbuf classification in `eam_derived` is dead code that misroutes 2D pbuf fields into out-of-bounds reads
`components/eam/src/physics/cam/physics_buffer.F90.in:769-795` (re-verified)

The function's contract claims `(pcols, 0, 0, ...)` → ndims=1, but
`pbuf_register_field_int` (`physics_buffer.F90.in:617-618`) does
`dimsizes(:,col_type_grid) = 1` then overwrites the leading `dimcnt` entries —
trailing dims are **1, never 0**. The loop `if (dimsizes(i) > 0)` counts all
six → **ndims = 6 for every registered pbuf field**. Consequences in
`eam_derived.F90`: `get_field` (`:654`) never takes the 2D branch, so a
horizontal-only pbuf operand (e.g. `PBLH`) is read as `(pcols, pver)` —
79 levels past the end; `get_field_ndims` (`:737`) and `is_field_2d` (`:794`)
are equally broken. Trigger: `derived_fld_defs='X=PBLH*2.0'` or
`tend_flds='PBLH'` — accepted by the validator, then garbage or DEBUG abort.
Not triggered by production (constituents only). **Fix:** compare `> 1` for
dims ≥ 2, or store the registered `dimcnt` in the header.

### E2. HIGH — sequence-association bug (gotcha #2 class) survives in `eam_derived` for 2D expressions with ≥2 field operands
`components/eam/src/control/eam_derived.F90:421, 454` vs
`share/util/shr_derived_mod.F90:229` (re-verified)

`tmp_fields(pcols, pver, max_operands)` is passed **whole** to a dummy
declared `fields(max_col, nlev, n_operands)`. For a 2D expression
(`nlev=1`), storage association maps dummy `fields(i,1,2)` onto actual
`tmp_fields(i,2,1)` — operand 1's (zeroed) level 2 — instead of
`tmp_fields(i,1,2)` where operand 2's data lives. Any 2D expression with two
or more field operands silently evaluates operands 2..N as **zero**. No
crash, no warning. Production's `STW` is 3D (`nlev=pver`, exactly
conforming), which is why this was never observed. This is precisely the
stride bug that gotchas #2/#31 document as fixed in `eam_vcoarsen` — it was
not fixed here. **Fix:** pass
`tmp_fields(:, 1:nlev, 1:expressions(i)%n_operands)`.

### E3. HIGH — `eam_vcoarsen` column-integration "previous valid" flag is a module-scope scalar: OMP race plus wrong-even-serial first-timestep tendencies
`components/eam/src/control/eam_vcoarsen.F90:104, 592-606` (re-verified)

`logical :: int_prev_valid = .false.` is shared across chunks while
`eam_vcoarsen_write` runs inside the `!$OMP PARALLEL DO` chunk loop
(`physpkg.F90:1419-1452`). Even serially, the first chunk processed on step 1
sets the flag and every later chunk computes `(integrated - 0)/dt` — garbage
`d*_INT_dt` on all chunks but one; with OMP, *which* chunks is
timing-dependent (non-reproducible). `int_prev` is also not
restart-persisted, so the same garbage recurs on the first step of every
restart leg. This is **exactly the scalar-flag pattern gotcha #16 claims was
eliminated** — it was fixed in `eam_derived` (per-chunk arrays) but not here.
Latent (production sets no `vcoarsen_int_flds`). **Fix:** make it
`logical, allocatable :: int_prev_valid(:)` over `begchunk:endchunk`,
mirroring `tend_initialized`.

### E4. MEDIUM — tendency state not restart-persisted; derived-operand phys tendencies emit garbage on the first post-restart step
`eam_derived.F90:122-123, 328, 396-403, 481-493`. `tend_prev` /
`tend_initialized` live only in memory → `d{NAME}_dt` outputs 0 for the first
post-restart step (breaks ERS for any `d*_dt` field on tape; pollutes one
sample of an averaged record). Worse: `eam_derived_stage` snapshots derived
operands from `derived_cache`, which is **all zeros** on the first step of a
leg, yet still sets `phys_snap_valid=.true.` → `d{DERIVED}_dt_phys =
(curr - 0)/dt`, a spuriously huge value with no gating. (`d*_dt_dyn` is
correctly suppressed.) Latent in production (no `tend_flds`).

### E5. MEDIUM — height/pressure select silently collapses out-of-range targets to the boundary level, with misleading metadata
`shr_vcoarsen_mod.F90:326-338`; long_names at `eam_vcoarsen.F90:397-399`.
The lowest EAM midpoint is ~10–25 m AGL, so **`Tat2m`/`Qat2m`/`STWat2m` are
identically the lowest-model-level values for essentially every column**, not
2 m values, while the long_name claims "linearly interpolated to 2.0 m above
surface". `Uat10m`/`Vat10m` genuinely interpolate only where the lowest
midpoint sits below 10 m. The namelist doc discloses the collapse; the
per-field metadata does not. ACE consumers comparing `Tat2m` to `TREFHT`
(both on the tape) will see a systematic, unexplained discrepancy. Pressure
select (unused) returns bottom-level values below ground instead of fill.
**Fix:** at minimum, correct the long_names; consider documenting in the
SamudrACE handoff that `*at2m` ≈ lowest model level.

### E6. MEDIUM — calendar rotation not excluded for satellite tapes
`cam_history.F90:4829-4833` gates rotation on `.not. is_initfile` and
`.not. restart` but not `is_satfile(t)`. Sat tapes force `nfils=1` in
`h_define` and are never date-stamped, so storage_type on a sat tape would
close/recreate files with untested semantics. One-line fix.

### E7. LOW — assorted
- `one_year` storage + `nhtfrq=0` picks the `%y-%m` filename spec
  (`cam_history.F90:777-796`): file named `*.0001-01.nc` containing Jan–Dec.
- After a restart-filename parse failure, the first post-restart write
  lazy-inits **without rotating** (`:4839-4840`) — only reachable with a
  custom `hfilename_spec`.
- `interp_native_interface` (`eam_vcoarsen.F90:915-924`): out-of-range
  `p_target` interpolates across the whole column (defaults `1, pverp`) —
  silently wrong ak/bk metadata. No check that
  `vcoarsen_level_bounds(n+1) <= pver`; a bound of 90 on an 80-level model
  indexes `hyai(91)` OOB at register time.
- `horiz_remap_mod.F90:248-249` checks fill on **level 1 only** (relies on
  EAM's z-invariant-fill convention); no guard that remapped-tape fields live
  on the physics grid (a dyn-grid fincl field would index `hbuf` with wrong
  coordinates); `pio_write_darray` ierr discarded (`:314, :337`).
- `add_hist_scalar`: `maxhistscalars=64` caps at n_avg_levs ≥ 32 while the
  namelist allows 50 (clean endrun); restart-history files carry
  defined-but-never-written `ak_*/bk_*` vars (cosmetic).
- Adiabatic physics + `tend_stages` set → `eam_derived_stage` dereferences
  unallocated `phys_snap` (segfault rather than message) — register is gated
  on `moist_physics`, the stage call (`cam_comp.F90:258-261`) is not.
- On `CONTINUE_RUN`, the restart copy of `horiz_remap_file` silently
  overrides any user change (`cam_history.F90:1848-1860`) — defensible
  (prevents an EAM analog of gotcha #46) but undocumented.

### Notes
- `Z3` in eam_derived/eam_vcoarsen maps to `state%zm` (height above
  *surface*), not the h-file `Z3` (geopotential above sea level) — a naming
  trap for anyone writing a derived expression.
- On an RRTMGP compset the production fincl1 would hit the FLDLST `endrun`
  on `FSUS`/`FLUS` (addfld'd only in `rrtmg/radiation.F90:644-646,700-702`)
  — a hard abort, not the "silently absent" that gotcha #43 predicts for the
  aerocom fields (those *are* registered unconditionally and would be
  fill/zero).
- Namelist XML closure verified: all 15 new Fortran namelist vars are in
  `namelist_definition.xml` with matching types/extents/valid_values;
  `vcoarsen_int_flds` lacks a `namelist_defaults_eam.xml` entry (harmless
  inconsistency). `restartvarcnt` 38→39 matches the one added restart var.
- `DTENDTTW`: `total_water_ac` is a `'global'` pbuf field → restart-persisted
  → BFB-safe. Correct.

### Verified correct
Gotcha #31 fix (all three select paths pass explicit `(1:ncol,:)` slices);
#16 per-chunk arrays in eam_derived (the race claim is false only for
eam_vcoarsen — see E3); #27 rotation ordering (check before "new volume"
branch, stamp on every non-restart write, restart parse seeds curr_idx,
hybrid lazy-init clean, mfilt is a pure safety net, rotation decisions
rank-deterministic); #30 validator coverage and native `PRECST`;
#41 ak/bk scalar plumbing incl. init-order (`hycoef_init` precedes
`phys_register`); #42 renames complete with no stale references; #43 A–D
(AOD*all/aerindexall outfld before the fill loops, no `add_hist_scalar('co2vmr')`
anywhere, per-record co2vmr confirmed at `cam_history.F90:4948-4953`); #47
EAM wrapper hard-aborts on negative weights. Remap write path always
PIO_DOUBLE + r8 (gotcha #1 not triggered); 3D `idof` ordering matches the
buffer layout.

---

## 3. MPAS-Ocean

### O1. MEDIUM — append-reopen schema-mismatch path: uninitialized `time_len`, writes proceed through a garbage frame index
`mpas_ocn_fme_horiz_remap.F:779, 832-897, 1389`

In the append branch, `pio_inq_*` failures accumulate into `err` but the code
continues: a failed `inq_dimid('time')` leaves `time_len` (no initializer)
uninitialized, `allocate(existing_time(int(time_len)))` can crash or OOM, and
`time_record` is seeded from garbage. Lines 884-886 set `file_is_open=.true.`
and `existing_file_reopened=.true.` **even when `err /= 0`**; `check_rotate`
returns early on the error without replaying def_vars, but the AMs gate
writes only on `file_is_open` — so `write_time` runs `pio_setframe` at an
arbitrary frame. The reopen-*open*-failure path correctly falls through to
clobber (`go to 100`); the inq-failure path should do the same. Trigger is
low-probability (an openable file with the wrong schema at the leg-start
filename — e.g. a legacy xtime-era file), but the consequence is undefined
behavior rather than a clean error.

### O2. MEDIUM — `buildnml` emits no `fmeVerticalReduceOutput` stream block
`components/mpas-ocean/cime_config/buildnml:1256-1345` (re-verified: zero
matches for `fmeVerticalReduce` in the file)

buildnml is the runtime truth (gotcha #28); the Registry stream is inert. If
`config_AM_fmeVerticalReduce_enable=.true.` is ever set, the
`MPAS_stream_mgr_get_property` query in `ocn_fme_remap_build_filename`
(`mpas_ocn_fme_horiz_remap.F:1450-1457`) fails and silently falls back to
`fmeVerticalReduceOutput.$Y-$M-$D.remapped.nc` — daily files, no casename,
wrong rotation, and they don't match the CIME archive regexes. **AGENTS.md
gotcha #28 is wrong** in claiming the buildnml blocks cover fmeVerticalReduce.

### O3. MEDIUM — `check_rotate` ignores all errors from the def_var replay; a failed def leaves an invalid `var_desc` slot live
`mpas_ocn_fme_horiz_remap.F:1395-1405`; `remap_file_def_var:1078-1083,
1101-1106`. The per-var `ierr` is overwritten each iteration and never
folded into `err`. On a `pio_def_var`/`pio_inq_varid` failure the slot's
`n_pio_vars` is already incremented with an invalid descriptor; a later
`write_var` matches by name and writes through it. Concrete scenario: add a
new output field to an AM, warm-restart mid-month — the reopened leg-1 file
lacks the var, inq fails, every write of that var that month is garbage.
**Fix:** accumulate errors; mark failed slots dead and skip them in
`write_var`.

### O4. MEDIUM — epoch sidecar lookup not gated on `is_restart_run`
`mpas_ocn_fme_horiz_remap.F:607, 644-675`. On a **cold start** in a run dir
containing a stale `<amName>.fme_accum_restart.nc` (CONTINUE_RUN=FALSE rerun,
reused dir), the epoch silently comes from the stale sidecar instead of
`MPAS_START_TIME`. Benign for same-case reruns (epoch identical); silently
wrong `time:units` relative to EAM if a dir is reused across cases with
different start dates. The loud WARN exists only for the opposite case.
**Fix:** gate the sidecar lookup on `is_restart_run` (or WARN when a cold
start adopts a sidecar epoch).

### O5. MEDIUM — non-short-circuit `.and.` referencing an absent optional dummy
`mpas_ocn_fme_horiz_remap.F:551`:
`if (present(mask) .and. mask(icell) < k)`. Fortran does not short-circuit;
referencing `mask` when absent is UB. Currently latent — every call site in
all three AMs passes `mask=` — but the argument is optional and the *other*
validity check in the same routine handles absence carefully. Nest the
conditions.

### O6. LOW — assorted
- A failed rotate/open silently drops the whole averaging window: writes are
  skipped but accumulators are zeroed and the schedule advanced (all AMs,
  e.g. `mpas_ocn_fme_depth_coarsening.F:890-938`). Deliberate trade-off
  (avoids a gotcha-#45-style multi-interval blend); deserves a loud ERROR
  log per lost window.
- `build-namelist` missing `add_default` for six FME options (they fall back
  to Registry defaults, so `mpaso_in` just doesn't show them);
  `_horiz_remap_file` add_default'd twice in two blocks.
- Whole compute (remap+accumulate+flush) runs **per block** inside the block
  loop; >1 block/rank would double-accumulate and double-write, and
  `send_local_cell` stores block-local indices with no offset. E3SM runs
  1 block/rank — add a fatal guard on `block%next` being associated.
- No finalize for the horiz_remap module (iodescs never freed — run-lifetime,
  not per-rotation); AM finalize deallocates the 5D arrays but not
  `remap_accum`/`remap_accum_count`.
- `fmeVerticalReduce` writes `kineticEnergy = 0` (not fill) when velocity
  diagnostics are unassociated (`:484-489`).
- Standalone MPAS `src/shared/Makefile` build broken (`shr_*`/`pio` deps not
  available outside the E3SM cmake build); irrelevant for E3SM.
- `nCoarsenLevels > nFmeDepthLevels` guard lives in compute, not init
  (`mpas_ocn_fme_depth_coarsening.F:744-751`) — still fatal-before-output,
  later than gotcha #9 documents.
- Gotcha #46 (map-change accumulator blend) confirmed **still unfixed**; the
  sidecar carries no map identity. Same hazard class: changing `depth_bounds`
  values while keeping the layer count across a restart is silently accepted.
- buildnml flips `timeSeriesStatsCustomOutput` to a 5-day interval and
  rewrites its var list — intentional for cross-verification, but it alters
  legacy tapes for any case built from this branch with that AM enabled.

### Verified correct
Gotchas #44 A/B/C (define-error close+reset on the fresh-create path; close
zeroes handles; fire-once `restart_reopen_pending` consumed unconditionally,
failed reopen falls through to clobber; the ocean 5D companion gets its own
latch via a separate `ocn_fme_remap_file_t` instance); #29 pt 2 (strict `<`
frame seeding, empty-file reseed of lon/lat); #45 (strict `>` rewind test in
all four schedules, full accumulator zeroing incl. the `_5D` trio,
`rewind_restart` in the `may_append` predicate; `remap_latest_inst` not
zeroed but overwritten before any flush — safe); the forward-restart mask fix
(the `existing_file_reopened` early-return is gone from all four mask
writers, `time_record == 1` gate kept); #36/#42 (per-cell counts at every
accumulate site, fill at count==0, distinct 5D sidecar keys, full
registry/namelist closure for the 5D options); full sidecar read/write
symmetry; #39 epoch ordering with `ocn_comp_mct.F:611-619` confirmed; #1/#3/
#5/#4/#7/#9/#10/#13/#15/#32/#37/#40/#41/#47 all verified as documented. The
ERS timeline (leg-1 flush → sidecar → leg-2 reopen → BFB overwrite of the
straddle frame) was traced end-to-end and is self-consistent.

---

## 4. MPAS-Sea-Ice

**Headline: the `_5D` duplication is clean.** Every primary/`_5D` procedure
pair was diffed after normalization, and — because normalization cannot catch
missed renames — the entire 5D region (lines 1329-1724) was swept for bare
primary-state identifiers and vice versa. Zero cross-contamination. The only
intentionally shared state is the PIO decompositions, remap data, epoch, and
provenance, exactly as the header comment declares.

### S1. MEDIUM — `cellMask` does not exist in the sea-ice Registry; the non-identity-map branch dereferences a null pointer
`mpas_seaice_fme_horiz_remap.F:259-272`. The ocean-only-indexed-map branch
(`n_a /= N_total`) does `mpas_pool_get_array(meshPool, 'cellMask', ...)` —
no such field exists anywhere in `components/mpas-seaice/src/` → unassociated
pointer → segfault. Latent (production maps are full-mesh, identity branch),
but a user supplying a masked SCRIP map gets an unexplained crash instead of
an error. **Fix:** remove the branch and `MPAS_LOG_CRIT` on
`n_a /= N_total`, or use a real field.

### S2. MEDIUM — `check_rotate` (primary and `_5D`) ignores per-variable `def_var` errors → stale-descriptor writes, gotcha-#44-class
`mpas_seaice_fme_horiz_remap.F:1015-1019` (primary), `1679-1683` (5D). Same
defect as O3, with a sharper consequence on the append path: a failed
`pio_inq_varid` leaves `pio_var_descs(iv)` holding the **previous file's**
descriptor while the name stays valid — the next `write_var` writes through
the stale varid into the new file. The gotcha #44 open/close fixes were
applied; this replay loop was not hardened. Identical in both copies.

### S3. LOW/MEDIUM — partially-failed append-reopen leads to clobber of leg-1 records on the next rotation check
`mpas_seaice_fme_horiz_remap.F:696-705` (and 5D `1435-1438`). On an inq
failure during reopen, `file_is_open=.true.` is still set with nonzero err;
`check_rotate` returns without setting `current_filename`; the next compute
sees a filename mismatch → close + reopen, but `restart_reopen_pending` was
already consumed → **clobber-create**, destroying leg-1's records in the
leg-start month. Acceptable as self-heal for a corrupt file; also fires on a
transient PIO error.

### S4. LOW — `uVelocityGeoCell`/`vVelocityGeoCell` silently zero if `config_use_high_frequency_coupling=.false.`
`mpas_seaice_fme_derived_fields.F:709-710, 769-785`. `uVelocityCell/
vVelocityCell` are populated only under HF coupling
(`mpas_seaice_icepack.F:1664-1668`). Default is `.true.`, so production is
fine, but there is no init-time guard; disabling it ships plausible-looking
rotated zeros in the training tape. Add a WARN/CRIT at AM init.

### S5. NOTES
- Flush discards the window on a failed file open (same design as ocean O6).
- Multi-block decomposition structurally unsupported (same as ocean O6).
- `pio_var_names` dual-indexed by `n_registered_vars` and `n_pio_vars`;
  correct only because check_rotate replays registration order — fragile
  invariant in both copies (the type-encapsulation refactor in AGENTS.md's
  Extensions list would fix it).
- `read_accum`/`try_read_epoch_from_sidecar_si` open sidecars with
  `PIO_WRITE` intent (`:514, :1226`) — fails on read-only files.
- A pre-#36 sidecar lacking `remap_accum_count` restores sums but leaves
  counts 0 → first post-restart flush emits fill (cross-version restart
  only).
- buildnml's native seaice stream blocks omit `u/vVelocityGeoCell` that the
  Registry lists — harmless (streams inert), inconsistent.

### Verified correct
Gotchas #44 A/B/C in both copies (incl. separate `restart_reopen_pending` /
`_5D` latches); #45 (`rewind_restart` gates both `may_append` predicates;
strict `>`; both accumulator trios zeroed; bonus: `write_time` reseeds
lon/lat on a reopened zero-record file — the seaice analog of the ocean mask
self-heal); #29 (strict `<` frame seed; sidecar hook fires from the coupler
restart-alarm branch in `ice_comp_mct.F` before stream writes; distinct
amName keys); #36 (per-cell counts both paths; extensive sums stay 0.0,
intensive fields get fill — the exact #25 contract); #39 CF time incl.
underscore/space pairing; #21 (`surfaceTemperatureMean` is the iceArea-
weighted category mean, degC everywhere); the geographic rotation call was
checked against both the `geographical_vectors` AM and — decisively — the
coupler export path (`ice_comp_mct.F:2720-2731`), same routine, same
argument order, no sign flip, so FME's airStress matches `Faii_taux`
semantics; #47 (the seaice wrapper **has** the loud negative-weight abort
that AGENTS.md #47 says is ocean-only — AGENTS.md understates the coverage);
#4, #1, #15; full namelist closure for all 9 options across Registry /
definition / defaults XMLs and both group lists.

---

## 5. Python verification scripts

### V1. HIGH — MPAS readiness certification can never PASS the production testmod
`verify_mpas.py:4055-4070, 4078-4079` (re-verified). The certification
streams dict hard-requires `fmeVerticalReduce`, and
`REMAPPED_VERTREDUCE_VARS` are folded into `expected_vars` — but the testmod
disables that AM (`shell_commands:303`, gotcha #33, since 2026-04-30).
**Every run of the current testmod produces TRAINING READINESS: FAIL and
`sys.exit(1)`.** AGENTS.md's "both dashboards converge cleanly" claim predates
the disable. **Fix:** drop the stream from certification or make it
conditional.

### V2. HIGH — ice-presence mask is a single-time snapshot applied to all timesteps
`verify_mpas.py:4256-4268, 4312-4348`. The ice mask is built from
`iceAreaTotal` at one record of the first file, then applied to every
timestep of every file. Ice cover is strongly seasonal: any cell icy at mask
time and melted later is fill at the later record → flagged "NaN/fill in
ocean cells". Fine for Ld10 smoke tests; guaranteed false FAILs on a
multi-season production tape — the stated purpose of certification. The mask
must be rebuilt per-record from the same record's `iceAreaTotal`.

### V3. HIGH — `verify_eam.py` cross_verify returns the wrong shape on SKIP paths → crash
`verify_eam.py:1364, 1372` (`return issues` — a bare list) vs the normal
3-tuple return (`:1961`) and the 3-way unpack at `:3407` (re-verified).
Triggers whenever h0 files are missing in either rundir **or the first h0
file is header-only (`time` dim = 0)** — exactly the mid-leg PIO-buffering
state of gotcha #28. `ValueError: not enough values to unpack` kills the
dashboard before HTML is written.

### V4. HIGH — `verify_production.py` always exits 0
`verify_production.py:1163-1169` catches every exception, prints a traceback,
and falls off the end (implicit exit 0); `main()` never aggregates `FAIL`
sanity entries into an exit code; a stream with zero files renders only as
"*no files found*" with no FAIL entry (`:1027-1030, :1142-1146`). A missing
stream, a fill-leak FAIL, a 5D cross-check FAIL, and a crash are all
indistinguishable from success to any caller or CI.

### V5. HIGH — `verify_eam.py` Test 1b TOPO_FILL gate is stale w.r.t. level-index coarsening
`verify_eam.py:1574-1601`. Expected fill is computed from pressure bounds
(`PS < VCOARSEN_PBOUNDS[k]`), but production uses level-index coarsening
(`shell_commands:176`), for which fill is 0 by construction (AGENTS.md
"production-readiness" section). Expected ~5-10% vs actual ~0% → `TOPO_FILL
DIFF` issue → exit 1 on every cross-verify against production output. Should
assert fill == 0 in level-index mode.

### V6. MEDIUM — assorted
- **5D-vs-1D cross-check thresholds**: FAIL is hard-coded at 1% while the CLI
  help and gotcha #42 say "1% WARN / 5% FAIL"; with the default
  `rel_tol=1e-2` the WARN branch is unreachable; three inconsistent default
  values for one knob (`verify_production.py:136-144, 265, 390-393`).
- **EAM certification variable list is stale vs production fincl1**: none of
  `Q_0..Q_7`, `Tat2m`/`Qat2m`/`STWat2m`/`Uat10m`/`Vat10m`, `TREFHT`,
  `QREFHT`, `PRECST`, `AODVISall`, `aerindexall`, `cdnc`, `lwp`, `lcc`,
  `ccn.3bl`, `colccn.3` appear in any check list in `verify_eam.py` — the
  entire 2m/10m + aerosol addition is unverified; those fields could be
  missing or garbage and certification still PASSes (`verify_eam.py:63-122`).
- **Radiation-budget closure uses an unweighted lat-lon mean** (the exact
  gotcha-#26 error mode) and gates at 50 W/m² instead of the advertised
  ~5 W/m² — both wrong mean and near-vacuous gate
  (`verify_eam.py:961-981, 3246-3266`).
- **Test 3 `same_grid` formula is wrong** (`stw_2d.size / N_VCOARSEN_LAYERS`
  is meaningless; true same-grid case gets the 2000×-looser cross-grid
  threshold) and `stw_2d` can be unbound → NameError
  (`verify_eam.py:1859-1861`).
- **Depth-layer auto-detect can collapse to 0 expected layers** → vacuous
  presence pass if `temperatureCoarsened_0` regresses
  (`verify_mpas.py:4082-4095`); fall back to 19 or FAIL.
- **Corrupt files silently reclassified as "legacy format"** (WARN only) in
  verify_production (`:165-174`).
- **Scaling**: eager per-file `xr.open_dataset` + `xr.concat` holds handles
  and materializes arrays — ~1200 EAM files likely trips the 1024-FD ulimit
  on the full 100-yr tape; enforce/document `--year` (`:968-985`).
- `push_to_share` docstring claims rebuild-from-filesystem; the code is
  additive insert (the insert itself is correct and idempotent per AGENTS.md,
  but it's an unguarded read-modify-write — concurrent pushes can lose an
  entry).

### V7. LOW — assorted
Dead `last_time_index()` (returns −1 on every branch, both scripts);
`safe_open` silently swallows xarray decode errors; potential `NameError` on
leftover `lons`/`lats` bindings (`verify_eam.py:841-851`); `RANGE_CHECKS`
covers only T/Q/PS/TS/PHIS so most EAM "physical range" certification is
vacuous; cross-grid legacy-vs-FME comparison has no functioning pass/fail
gate (NaN diff_rms never trips `>1e3`); 13 of the 23 ocean derived-fields
variables (the freshwater/heat companion fluxes) are never presence/range
checked by verify_mpas; truncated 5D windows at leg boundaries compared as if
complete; `user_nl_eam` embedded in HTML unescaped; exact float equality
against `SHR_FILL_VALUE` (verified safe for the current dtypes).

### Verified correct
Gotchas #19 (mask bootstrap first-pass, ordering-bug fix real and complete),
#20 (`GMEAN_ABS_TOL` set, `ssh_abs_tol=3 m`), #21 (degC at all four sites),
#22 (all fill-threshold sites at 1e18 — now 23 sites, not 22), #23
(`ICE_FME_ONLY`), #24 (`DEPTH_BOUNDS` byte-identical to `shell_commands`,
active-layer autodetect), #25 (ice mask = fill OR exactly 0), #26 (Test 3
cos-lat fallback — modulo the `same_grid` bug above), #42 (sample-vs-day
weighting documented, fill-window filtering, p99/p99 normalization, 1D globs
correctly exclude 5D companions, zero stale `totalFreshWaterTemperatureFlux`
or `T_at_z2`-style references anywhere). All three scripts use the Agg
backend, close every figure, and use `decode_times=False` consistently;
header-only files are handled everywhere except the V3 path.

---

## 6. Configuration, scripts, MOSART

### C1. HIGH — 1D and 5D FME streams collapse into one ERS comparison bucket
`cime_config/config_archive.xml:88, 107-109` (re-verified: no 5D entries
exist). The per-stream `<hist_file_ext_regex>` entries are unanchored, so
`*.fmeDepthCoarsening5D.0001-01.remapped.nc` matches the 1D regex and lands
in the same bucket (`archive_base.py:_get_extension` keys on the matched
group). `get_latest_hist_files` keeps **one file per bucket**, winner
determined by directory-listing order. Consequences for ERS: (a) one of
{1D, 5D} per AM is silently never restart-BFB-compared; (b) `.base` and
`.rest` are bucketed in separate scans, so base can pick 1D while rest picks
5D → spurious NO_COMPARE flake. Affects all three AMs. Gotcha #34's "each
FME stream gets its own ext-bucket" was true on 2026-04-30 and went stale
when the 5D streams landed 2026-05-09. **Fix:** add
`hist\.am\.fmeDepthCoarsening5D` (etc.) entries *before* the 1D ones, or
anchor the 1D regexes with a trailing `\.`; re-create test cases afterward
(per gotcha #34's note). Archiving itself is unaffected.

### C2. MEDIUM — production script issues (`run_e3sm.fme.sh`)
- **`-cosp` appended to CAM_CONFIG_OPTS** (`:353`) with zero COSP fields on
  the FME tape: `docosp=.true.` runs the simulator every few radiation steps
  for 100 years producing nothing. Drop it.
- **`WALLTIME=34:00:00`** (`:135`) exceeds the 24 h pm-cpu regular-QOS cap;
  CIME only warns, Slurm will likely reject at submission. The 5-leg
  RESUBMIT design absorbs a shorter walltime anyway.
- **`PROJECT` fallback fragile** (`:39`): on non-Slurm machines the `sacctmgr`
  substitution fails silently under `readonly`, the unquoted `--project
  ${PROJECT}` collapses, and create_newcase misparses the next flag. Breaks
  the advertised `MACHINE=chrysalis` path.
- **PLACEHOLDER guard is dead code as shipped** (`:161-171`): all defaults
  are real values so it can never fire; it checks neither `RUN_REFDATE` nor
  `START_DATE` despite the error message naming RUN_REFDATE;
  case-sensitive match.
- `if [ $? != 0 ]` after create_newcase unreachable under `bash -fe`
  (`:323-328`); `BRANCH="mahf708/fme/aigo"` cloned from E3SM-Project — no
  guarantee it equals this reviewed branch (provenance drift); stale
  header comments.

### C3. MEDIUM — embedded map-regeneration doc produces the wrong files
`fme_output/shell_commands:50-91` documents plain `ncremap` (default
algorithm) producing `..._shifted.nc`, but detection (`:107-108`) expects
`..._shifted_trintbilin.nc` — and with the gotcha-#47 negative-weight
refusal, regenerating with the wrong algorithm and renaming aborts at init.
The doc needs `--alg_typ` and the `_trintbilin` names.

### C4. LOW
- `FME_MPAS_INTERVAL_5D=` (empty) silently re-enables the 5D stream — `:-`
  substitutes on null; only `'none'` disables (`shell_commands:26-29`).
  AGENTS.md says "'none' (or omit)", the inline comment wrongly adds
  "or empty".
- mpassi leaves `write_on_startup=.true.` (mpaso sets `.false.`) — harmless
  (driver skips on `output_stream=='none'`) but logs a WARNING every leg.
- `FME_EAM_OUTPUT_HOURS=0` → `nhtfrq=0` → monthly averages, silently very
  different from intent; no env-var validation.
- `fme_legacy_output` sets `cosp_lite=.true.` without `-cosp` (inert);
  empty `user_nl_mosart` heredoc is dead code.
- Stale docs: AGENTS.md says `FME_MAPS_DIR` defaults to NERSC scratch — the
  actual default is `$DIN_LOC_ROOT/fme` (`shell_commands:40`); gotcha #35's
  `shell_commands:217` is now line 285.

### MOSART — verified correct
The `ntapes==0` restart-read guard (`RtmHistFile.F90:1404-1416`) covers both
`locfnh` and `locfnhr` reads plus the strip-null loop in one guarded block;
downstream `do t=1,ntapes` loops no-op naturally. `rtmhist_empty_htapes`
declared/read/bcast correctly with a matching XML entry; semantics mirror ELM
exactly — fincl-requested fields always create their tape, so **no wanted
history file can be discarded**. Restart *write* with ntapes==0 still defines
a zero-length dim + vars (empirically fine, ERS-validated; worth one check
per PIO iotype). MOSART remains `exclude_testing` in config_archive, so the
guard is exercised by crash/no-crash, not cprnc — consistent with intent.

### Config — verified correct
Gotcha #34's regex math holds end-to-end against `archive_base.py:139-148`
(the `\.nc$`-suffix-append trap, `.base`/`.rest` suffix forms, no match on
`.rst.*`, sidecars excluded by the casename-prefix filter, rpointer specs
untouched) — modulo the C1 5D gap. All spot-checked fincl1 fields exist with
correct registration paths (incl. `PRECST` native per gotcha #30, `FLUS`/
`FSUS` on RRTMG per #43). Full namelist plumbing closure for EAM and both
MPAS components. Depth bounds byte-identical between `shell_commands` and
`verify_mpas.py`. All six FME stream blocks (3 AMs × 1D+5D) emitted by
buildnml with `$Y-$M` templates and `output_interval="none"` — except the
missing fmeVerticalReduce block (O2). Production cadence math checks out
(5 × 20 yr legs, REST_N=5 on year boundaries → append path unexercised in
production, as designed). `RUN_STARTDATE` after build is safe (case.submit
re-runs buildnml; branch runs are governed by RUN_REFDATE). `_inst`
companions live inside the primary files — no extra archive entries needed.

---

## 7. AGENTS.md corrections needed

1. **Gotcha #28** claims the buildnml stream blocks cover
   "fmeDepthCoarsening, fmeDerivedFields, fmeVerticalReduce" — there is no
   fmeVerticalReduce block (O2).
2. **Gotcha #34** ("each FME stream gets its own ext-bucket and all of them
   survive") is stale since the 5D streams landed (C1).
3. **Gotcha #47** says the loud init abort is ocean-only ("add to the
   seaice/EAM wrappers if..."), but both the seaice wrapper
   (`mpas_seaice_fme_horiz_remap.F:224-233`) and the EAM wrapper
   (`horiz_remap_mod.F90:110-119`) already have it. Coverage is better than
   documented.
4. **Gotcha #16** ("the OMP race was eliminated") is true for `eam_derived`
   but false for `eam_vcoarsen`'s `int_prev_valid` (E3).
5. **Gotcha #9** — the `nCoarsenLevels > nFmeDepthLevels` guard fires at
   first compute, not init.
6. **Gotcha #22** — 23 fill-threshold sites now, not 22.
7. **Runtime Configuration** — `FME_MAPS_DIR` default is `$DIN_LOC_ROOT/fme`,
   not NERSC scratch; gotcha #35's line citation is stale.
8. The "Production-readiness assessment / both dashboards converge cleanly"
   section predates the fmeVerticalReduce disable and the fincl1 expansion;
   per V1/V2/V5/M2 the dashboards currently cannot converge on the
   production config.

---

## 8. Recommended fix order

**Blockers for the 100-yr submit (test/verification integrity + submission):**
1. **C1** — config_archive.xml 5D ext-regex entries (two-line XML edit;
   re-create cases). Until fixed, ERS does not actually verify 5D restart
   BFB.
2. **V1** — drop/conditionalize fmeVerticalReduce in verify_mpas
   certification (cert is dead-on-arrival).
3. **V3** — fix cross_verify SKIP-path return shape (crash).
4. **V4** — verify_production exit-code aggregation (CI blindness).
5. **S1/S2 (script)** — drop `-cosp`; fix walltime to ≤24 h.

**High-value correctness (one namelist edit from triggering):**
6. **E3** — `int_prev_valid` per-chunk (5-line fix, removes
   nondeterminism).
7. **E1** — `pbuf_get_field_ndims` trailing-1s fix.
8. **E2** — explicit section in the `shr_derived_eval` call.
9. **U1** — error out on dropped unowned gcols in build_comm.
10. **O3/S2 (seaice)** — harden the check_rotate def_var replay in both
    remap modules (stale-descriptor write path).

**Before any config variation in production:**
11. **O2** — add the fmeVerticalReduceOutput buildnml block (or a loud
    init abort when the stream query fails).
12. **V2** — per-timestep ice mask in certification (false FAILs on
    multi-season tapes).
13. **V5/M2/M3 (python)** — level-index TOPO_FILL gate; extend EAM cert
    lists to the production fincl1; cos-lat weighting in the radiation
    budget.
14. **E5** — fix `*at2m` long_names (metadata honesty for the training
    pipeline).
15. **Gotcha #46 proper fix** — stamp map identity into the sidecars and
    reset accumulators on mismatch (already designed in AGENTS.md, not yet
    implemented).

Everything else is hardening that can land opportunistically.
