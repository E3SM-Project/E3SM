# Campaign: Water Isotope Infrastructure for EAMxx

## Metadata

```yaml
campaign_id: water-isotope-infrastructure
title: Water Isotope Infrastructure for EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-10
status: draft
base_branch: wiso-dev
repository_path: /Users/rfiorella/code/E3SM/EAMXX-wiso
branch_prefix: rfiorella/eamxx/water-isotope-infra
execution_mode: continuous
estimated_total_hours: 26
```

## Campaign Overview

This campaign implements the foundational infrastructure required for water isotope and tracer capability in EAMxx. It establishes the process registration framework, multi-component tracer definitions, and surface flux extensions needed to track water isotopes (HDO, H2-18O, HTO) through the atmospheric component of E3SM.

The campaign consists of three specs executed in sequence with controlled dependencies:
1. Register water tracer processes in the EAMxx physics factory, with namelist-defaults plumbing so the processes can be placed in the atm DAG
2. Define multi-component tracer fields with runtime-configurable component count, registered as advected tracers
3. Extend surface evaporation flux to support isotope-specific boundary conditions

This is the first campaign in a multi-phase effort to add full water isotope capability to E3SM following conventions used in other ESMs (CESM, GISS ModelE, ECHAM).

## Execution Plan

```yaml
execution_mode: continuous
pause_between_specs: false
create_prs: draft
pr_labels: []
pr_size_limits:
  max_files_per_pr: 25
  max_lines_per_pr: 1000  # raised from 600: spec 002 includes standalone test scaffolding
```

## Specs

### Spec 001: Register Water Tracer Processes

```yaml
spec_id: 001-register-water-tracer-processes
spec_path: specs/001-register-water-tracer-processes.md
order: 1
estimated_hours: 8
dependencies: []
branch_name: rfiorella/eamxx/water-isotope-infra-001-register-processes
pr_target: wiso-dev
pr_title: "Register water_tracers and water_isotopes processes in EAMxx"
pr_draft: true
validation_tier: 0
critical_path: true
```

**Summary:** Creates process stubs for `water_tracers` and `water_isotopes` (`WaterIsotopes` subclasses `WaterTracers` per campaign instructions: isotopes are a special case of tracers) and registers them in the EAMxx AtmosphereProcessFactory. Establishes CMake build configuration following the P3 pattern (PUBLIC `EAMXX_HAS_*` compile definitions, link into `eamxx_physics`) and adds the `<water_tracers>` namelist-defaults block (with `tracer_count`) so the processes can be added to `atm_procs_list` via atmchange. The default process list is NOT modified.

**Key Deliverables:**
- `components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.{hpp,cpp}` (new directory — greenfield)
- `components/eamxx/src/physics/water_tracers/eamxx_water_isotopes_process_interface.{hpp,cpp}`
- `components/eamxx/src/physics/water_tracers/CMakeLists.txt`
- Updated `components/eamxx/src/physics/CMakeLists.txt` (add_subdirectory)
- Updated `components/eamxx/src/physics/register_physics.hpp`
- Updated `components/eamxx/cime_config/namelist_defaults_eamxx.xml` (defaults blocks only)

**Success Gates:**
- CMake build completes for the water_tracers library target
- Both processes registered in `register_physics.hpp` with `EAMXX_HAS_*` guards
- Namelist defaults block present so atmchange accepts the processes
- Full EAMxx case build completes with the library linked

### Spec 002: Multi-Component Water Tracer Field Definitions

```yaml
spec_id: 002-multi-component-tracer-fields
spec_path: specs/002-multi-component-tracer-fields.md
order: 2
estimated_hours: 12
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
branch_name: rfiorella/eamxx/water-isotope-infra-002-tracer-fields
pr_target: wiso-dev
pr_title: "Add multi-component water tracer fields (qv_iso, qc_iso, qi_iso, qr_iso families)"
pr_draft: true
validation_tier: 0
critical_path: true
```

**Summary:** Defines four tracer field families (`qv_iso`, `qc_iso`, `qi_iso`, `qr_iso` — mirroring P3 mass tracers; EAMxx names cloud liquid `qc`, and `qm` is ice rime mass, folded into qi for now) with component count runtime-configurable via the `tracer_count` parameter. Because EAMxx's advected "tracers" group is a monolithic (col, Q, lev) array and `add_tracer` only accepts scalar 3D fields, each component is registered as an individual scalar tracer (gaining HOMME/turbulence advection automatically) and grouped per family to provide the (col, cmp, lev) view. Includes a standalone single-process test under `tests/single-process/water_tracers/`.

**Key Deliverables:**
- `components/eamxx/src/physics/water_tracers/field_registry.cpp`
- Updated `components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.{hpp,cpp}`
- Updated `components/eamxx/src/physics/water_tracers/CMakeLists.txt`
- `components/eamxx/tests/single-process/water_tracers/{CMakeLists.txt,input.yaml,output.yaml}`

**Success Gates:**
- All four field families registered with kg/kg units and long_name metadata
- Component count runtime-configurable (no hardcoded sizes)
- Standalone single-process test builds, runs, and writes fields to history output
- BFB regression: 5-day F2010-SCREAMv1 run with tracer_count=1 is bit-for-bit identical to baseline b58a0fbea6 (cprnc)
- Case build completes

**Dependency Rationale:** Requires process registration and the tracer_count namelist plumbing from spec 001.

### Spec 003: Surface Flux Tracer Extension

```yaml
spec_id: 003-surface-flux-tracer-extension
spec_path: specs/003-surface-flux-tracer-extension.md
order: 3
estimated_hours: 6
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
branch_name: rfiorella/eamxx/water-isotope-infra-003-surface-flux
pr_target: wiso-dev
pr_title: "Extend surface evaporation flux to multi-component tracer field"
pr_draft: true
validation_tier: 0
critical_path: true
```

**Summary:** Extends the surface coupling importer to define `surf_evap_iso` as a (col, cmp) field following the importer's existing vector2d/vector4d pattern (`surf_mom_flux`, `dstflx`). Registered only when `tracer_count > 0`, zero-initialized with timestamp set in `initialize_impl` (the import loop only writes coupler-named fields, so an unguarded Computed field would risk timestamp postcondition failures in every EAMxx run).

**Key Deliverables:**
- Updated `components/eamxx/src/control/atmosphere_surface_coupling_importer.{hpp,cpp}`

**Success Gates:**
- `surf_evap_iso` registered with (col, cmp) layout, kg/m2/s units, guarded on tracer_count
- Build completes in both configurations (tracer_count absent/0 and =1)
- BFB regression vs baseline b58a0fbea6 in both configurations — the importer runs in every EAMxx case, making this the campaign's largest BFB risk

**Dependency Rationale:** Requires the tracer_count parameter plumbing from spec 001 but is independent of 3D tracer definitions (spec 002). Surface fluxes are 2D boundary conditions.

**Parallel Execution Note:** Specs 002 and 003 can be executed in parallel after spec 001 completes, as they touch non-overlapping code regions and have no mutual dependencies.

## Dependency Graph

```
001-register-water-tracer-processes
├── 002-multi-component-tracer-fields
└── 003-surface-flux-tracer-extension
```

**Critical Path:** 001 → 002 (longest path at 20 hours)

**Parallel Opportunities:** Specs 002 and 003 are independent after 001 completes. If parallel execution is enabled, they can run simultaneously to reduce total wall time from 26 hours to 20 hours.

## Integration Points

### Component Scope
- **EAMxx only:** All changes confined to `components/eamxx/`
- No modifications to dynamics, coupler, land model, or other E3SM components
- No changes to CIME infrastructure or case control (namelist_defaults_eamxx.xml is EAMxx-owned config)

### File Ownership
Each spec touches distinct files with no overlap:
- **Spec 001:** `water_tracers/` directory creation, `register_physics.hpp`, `physics/CMakeLists.txt`, `namelist_defaults_eamxx.xml`
- **Spec 002:** `water_tracers/` interior files (field_registry, process implementation), `tests/single-process/water_tracers/`
- **Spec 003:** `control/atmosphere_surface_coupling_importer.*`

### Build System Impact
- New CMake library target: `water_tracers`, linked into the `eamxx_physics` INTERFACE target
- New preprocessor guards: `EAMXX_HAS_WATER_TRACERS`, `EAMXX_HAS_WATER_ISOTOPES` — set automatically as PUBLIC compile definitions by the library's CMakeLists (P3 pattern); they are NOT CIME XML variables or user-facing CMake options
- Runtime enablement is per-case via `atmchange atm_procs_list+=water_tracers`; default configurations unaffected
- No changes to existing targets or dependencies

### Testing Strategy
- **Validation Tier 0:** compile + standalone single-process test + BFB regression vs baseline
- **Target Compset:** F2010-SCREAMv1
- **Target Resolution:** ne4pg2_ne4pg2
- **Test Platform:** pm-cpu or other supported E3SM machine. The development laptop (macOS) is not a supported machine; all case.build/case.submit/test-all-eamxx gates must run remotely.
- **Baseline Tag:** b58a0fbea6 — baseline generated per the regression-baseline procedure with a BASELINE.txt manifest (parent SHA, date, compset, resolution, run length); shared by specs 002 and 003
- **BFB Contract:** Infrastructure-only changes must leave existing answers bit-for-bit unchanged (cprnc). Follow-on campaigns that extend bulk water arrays inherit the slice-1-unchanged contract (bulk-H2O slice BFB vs pre-extension baseline). Conservation closure tests (per tracer-conservation tolerance defaults) begin only when physics writes nonzero tracer values.

## Out of Scope

```yaml
out_of_scope:
  - Full implementation of water isotope fractionation physics
  - Integration with P3/SHOC microphysics kernels
  - Separate q*_trace (non-fractionating) array family (deferred; per-component runtime naming permits adding it later without structural change)
  - Rime-specific isotope tracking (qm); rime isotope content folded into qi until P3 integration
  - Adding water_tracers/water_isotopes to any default atm_procs_list
  - Coupling to land model or coupler for isotope-specific boundary conditions
  - Input data preparation or initial condition files
  - Performance optimization or GPU-specific tuning
  - Conservation closure validation (meaningless on zero-valued fields; starts with physics integration)
  - Multi-tracer configurations beyond single-component testing
  - Dynamics modifications (components join the standard advected tracer group; no dycore changes)
  - Documentation or user guide updates
```

## Success Criteria for Campaign Completion

```yaml
campaign_success_criteria:
  - id: CSC1
    description: All three specs pass their individual success criteria
    verification: Review each spec's validation output

  - id: CSC2
    description: Full EAMxx build succeeds and both processes can be enabled in a case
    verification: |
      cd cime/scripts
      ./create_newcase --case ../../test_water_iso --compset F2010-SCREAMv1 --res ne4pg2_ne4pg2 --mach <machine>
      cd ../../test_water_iso
      ./case.setup
      ./atmchange atm_procs_list+=water_tracers
      ./atmchange water_tracers::tracer_count=1
      ./case.build
      # Note: EAMXX_HAS_* are compile definitions set automatically when the
      # library links (P3 pattern); they are not xmlchange-able CIME variables.

  - id: CSC3
    description: No regression in existing EAMxx tests
    verification: |
      cd components/eamxx
      ./scripts/test-all-eamxx -m <machine> -t dbg

  - id: CSC4
    description: All PRs successfully merge to wiso-dev without conflicts
    verification: Manual review of PR merge history

  - id: CSC5
    description: BFB regression holds — default-configuration 5-day run identical to baseline b58a0fbea6 after all three specs merge
    verification: cprnc of final history file vs baseline (see regression-baseline procedure); spec 002/003 gates re-run on merged wiso-dev

  - id: CSC6
    description: Standalone water_tracers single-process test passes in CI-style invocation
    verification: ctest -R water_tracers in EAMxx standalone build
```

## Post-Campaign Integration

After campaign completion, the following integration work is required (tracked in separate specs/campaigns):

1. **Physics Integration:** Implement actual fractionation physics in water_isotopes process (unit tests vs published alpha values per tracer-conservation reference-data conventions)
2. **Microphysics Coupling:** Connect to P3/SHOC for phase-change fractionation; revisit rime (qm) isotope treatment
3. **Boundary Conditions:** Coupler field-list extension plus land/ocean isotope sources for realistic surface fluxes
4. **q*_trace family:** Non-fractionating tagged-water tracers, if separate arrays are preferred over alpha=1 components
5. **Validation:** Tier 1-3 testing with observational datasets; conservation closure per tracer-conservation tolerances
6. **Documentation:** User guide, namelist documentation, example configurations

## Risk Assessment

```yaml
risks:
  - risk: Monolithic per-family subgroup view (col, cmp, lev) proves infeasible in the field manager
    mitigation: Spec 002 Phase 1 investigation resolves go/no-go early; fallback is iteration over named per-component subfields, which preserves all functionality
    severity: medium

  - risk: Adding components to the advected tracers group changes dycore qsize handling and breaks BFB even with zero-valued tracers
    mitigation: BFB gate in spec 002 catches it immediately; ask-before item requires consultation before any tolerance-based fallback
    severity: medium

  - risk: Importer change perturbs existing runs (importer executes in every EAMxx case)
    mitigation: tracer_count>0 guard plus dual-configuration BFB gate in spec 003
    severity: medium

  - risk: No supported machine access during execution (development host is macOS, unsupported)
    mitigation: All build/run gates scheduled on pm-cpu or equivalent; checkpoint execution lets implementation proceed locally with remote validation batched
    severity: medium

  - risk: Field manager API changes between spec creation and execution
    mitigation: Investigation phase in spec 002 detects API mismatches early
    severity: low

  - risk: Parallel execution of specs 002 and 003 causes merge conflicts
    mitigation: File ownership is disjoint; no shared files modified
    severity: low

  - risk: Baseline tag b58a0fbea6 becomes stale during campaign execution
    mitigation: Follow regression-baseline refresh policy (diff parent SHAs; regenerate to a -rN suffixed directory only if relevant code paths changed; never overwrite)
    severity: low
```

## Author Notes

This campaign establishes infrastructure only. No physics is implemented. All fields are initialized to zero. The goal is to enable subsequent campaigns to add fractionation physics and microphysics integration without revisiting basic infrastructure.

Design decisions of record (2026-06-10 revision):
- `WaterIsotopes` subclasses `WaterTracers`, honoring the campaign instruction that isotopes are a special case of water tracers.
- Components are individual advected scalar tracers grouped per family — not bespoke (col, cmp, lev) fields outside the tracer group, which would never be advected (add_tracer hardcodes the scalar 3D layout; the tracers group is monolithic (col, Q, lev) consumed by HOMME).
- Field families mirror P3 mass tracers: qv_iso, qc_iso, qi_iso, qr_iso. EAMxx uses qc for cloud liquid (not ql); qm is ice rime mass (not mixed-phase) and is deferred.
- The `_iso` suffix (rather than `_HDO`) supports future multi-isotope configurations without field renames.
- Validation replaces the original trivially-passing "mass conservation on zero fields" gates with the BFB-vs-baseline contract plus a standalone single-process test; conservation closure testing begins when physics produces nonzero values.

All specs use checkpoint execution mode to enable human review at natural breakpoints (investigation reports, build verification, BFB results).

## Campaign Execution Timeline

**Estimated Duration:** 26 hours total
- Spec 001: 8 hours
- Spec 002: 12 hours (can run parallel with 003 after 001)
- Spec 003: 6 hours (can run parallel with 002 after 001)

**Critical Path:** 001 → 002 = 20 hours

**Parallel Execution:** If enabled, specs 002 and 003 run simultaneously after 001, reducing total time to 20 hours.

**Continuous Mode:** Campaign runs without pause between specs. Each spec's success gates must pass before proceeding to dependents.

## Validation Status

Campaign manifest revised 2026-06-10 against codebase audit:

- **C1 (Metadata):** ✓ All required fields present
- **C2 (Specs List):** ✓ All three specs listed with paths, order, hours, dependencies
- **C3 (Dependencies):** ✓ Dependency graph clear; 001 → {002, 003} with rationale
- **C4 (Integration):** ✓ File ownership, build impact, testing strategy documented; build-guard mechanics corrected to match P3 pattern
- **C5 (PR Strategy):** ✓ Branch naming, PR titles, draft status, size limits specified
- **C6 (Success Criteria):** ✓ Six campaign-level criteria; CSC2 corrected (atmchange, not xmlchange); CSC5 now a BFB gate; CSC6 covers the standalone test
- **C7 (Out of Scope):** ✓ Campaign-level exclusions documented, including deferred q*_trace family and qm
- **C8 (Timeline):** ✓ Estimated hours per spec, critical path, parallel opportunities identified
- **C9 (Context):** ✓ Campaign overview explains purpose and relationship to broader effort

**Status:** Ready for execution with `/run-campaign campaigns/2026-06-09-water-isotope-infrastructure.md` (validation gates require a supported E3SM machine)
