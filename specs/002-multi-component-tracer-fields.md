# Spec: Multi-Component Water Tracer Field Definitions in EAMxx

## Metadata

```yaml
spec_id: 002-multi-component-tracer-fields
title: Multi-Component Water Tracer Field Definitions in EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-10
status: draft
spec_type: model-e3sm-eamxx
campaign_id: water-isotope-infrastructure
priority: high
estimated_effort_hours: 12
execution_mode: checkpoint
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
    description: Requires water_tracers and water_isotopes processes to be registered, plus the namelist defaults block for tracer_count
```

## Model-Specific Configuration

```yaml
model_specific:
  target_compset: F2010-SCREAMv1
  target_resolution: ne4pg2_ne4pg2
  validation_tier: 0  # compile + standalone test + BFB regression vs baseline (5-day run)
  requires_input_data: true  # BFB regression run needs F2010-SCREAMv1 inputdata
  platform: pm-cpu or other supported E3SM machine (local macOS cannot run CIME builds)
  baseline_tag: b58a0fbea6
  stop_n: 5
  stop_option: ndays
  affects_components:
    - eamxx
  affects_compsets:
    - F2010-SCREAMv1
```

## Inputs

```yaml
inputs:
  - path: components/eamxx/src/share/field/field_manager.hpp
    description: EAMxx field manager for field registration and allocation
    format: C++ header
    required: true

  - path: components/eamxx/src/share/field/field_request.hpp
    description: FieldRequest/GroupRequest structures, including monolithic group allocation
    format: C++ header
    required: true

  - path: components/eamxx/src/share/atm_process/atmosphere_process.hpp
    description: add_field/add_tracer/add_group helpers; note add_tracer hardcodes the scalar 3D layout (lines ~375-393)
    format: C++ header
    required: true

  - path: components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    description: Example of tracer registration via add_tracer (qv, qc, qi, qr, qm at lines ~78-86)
    format: C++ source
    required: true

  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.hpp
    description: Water tracers process interface (from spec 001)
    format: C++ header
    required: true

  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    description: Water tracers process implementation (from spec 001)
    format: C++ source
    required: true

  - path: components/eamxx/tests/single-process/p3
    description: Example standalone single-process test layout to mirror for water_tracers
    format: directory
    required: true
```

## Deliverables

```yaml
deliverables:
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.hpp
    description: Updated process interface with tracer registration for multi-component water tracers
    format: C++ header

  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    description: Updated create_requests() registering per-component tracers for qv_iso, qc_iso, qi_iso, qr_iso families
    format: C++ source

  - path: components/eamxx/src/physics/water_tracers/field_registry.cpp
    description: Helper that generates per-component tracer names and registers them with metadata (units, long_name) and group membership
    format: C++ source

  - path: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    description: Updated to include field_registry.cpp in build
    format: CMake

  - path: components/eamxx/tests/single-process/water_tracers/CMakeLists.txt
    description: Standalone single-process test for water_tracers (mirrors p3 test layout)
    format: CMake

  - path: components/eamxx/tests/single-process/water_tracers/input.yaml
    description: Test input with tracer_count=1
    format: YAML

  - path: components/eamxx/tests/single-process/water_tracers/output.yaml
    description: Output stream config exercising history output of the new fields
    format: YAML
```

## Success Criteria

```yaml
success_criteria:
  - id: SC1
    phase: implementation
    description: All four isotope tracer field families are registered in the process
    type: shell
    command: |
      grep -c -E "(qv_iso|qc_iso|qi_iso|qr_iso)" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/field_registry.cpp
    assertion: count >= 4
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Missing field registrations mean tracers cannot be allocated or tracked
    note: Cheap textual pre-check only; SC3/SC5 are the real gates

  - id: SC2
    phase: implementation
    description: Component count is runtime-configurable (tracer_count read from ParameterList; no hardcoded component sizes)
    type: shell
    command: |
      grep -q "tracer_count" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    failure_mode: Component count baked in at compile time defeats the runtime-configuration requirement

  - id: SC3
    phase: validation
    description: Standalone single-process test builds, runs, and writes all four field families to history output
    type: shell
    command: |
      cd ${EAMXX_BUILD_DIR} && ctest -R water_tracers --output-on-failure
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/tests/single-process/water_tracers/CMakeLists.txt
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Process cannot be instantiated by factory, fields fail allocation, or output registration broken

  - id: SC4
    phase: implementation
    description: Field metadata includes units (kg/kg) and long_name for all isotope field families
    type: shell
    command: |
      grep -q "kg/kg" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/field_registry.cpp && \
      grep -q "long_name" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/field_registry.cpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Missing metadata causes I/O errors or unintelligible output files

  - id: SC5
    phase: validation
    description: BFB regression — enabling water_tracers with tracer_count=1 (zero-initialized, no physics) leaves existing model answers bit-for-bit unchanged vs baseline b58a0fbea6
    type: comparison
    command: |
      cd ${CASE_DIR} && \
      ./xmlchange STOP_OPTION=ndays,STOP_N=5 && \
      ./case.submit && \
      cprnc ${BASELINE_DIR}/h0.final.nc ${RUNDIR}/h0.final.nc
    expect: BFB
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    failure_mode: New tracers perturb existing physics (e.g., via tracer-group resizing affecting dycore state); refactor is functional, not structural

  - id: SC6
    phase: validation
    description: Build completes with multi-component tracer fields enabled
    type: shell
    command: |
      cd ${CASE_DIR} && ./case.build
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/CMakeLists.txt
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Compilation errors from incorrect field dimensions or template instantiation failures
```

## Out of Scope

```yaml
out_of_scope:
  - Implementation of fractionation physics (deferred to integration spec)
  - Custom tracer advection schemes (components join the standard advected tracers group; HOMME/SHOC handle transport)
  - Coupling to microphysics parameterizations (P3/SHOC integration - separate spec)
  - Boundary condition specification for isotope tracers
  - Input data preparation or initial condition files
  - Performance optimization or GPU-specific memory layouts
  - Separate q*_trace array family (see resolved decision; deferred)
  - Mass-conservation closure tests (meaningless on zero-valued fields; deferred to physics integration specs per tracer-conservation conventions)
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: Register each component as an individual scalar tracer via add_tracer, with names generated at runtime (e.g., qv_iso_0 .. qv_iso_{N-1}); the (col, cmp, lev) view is recovered through a field group per family
    rationale: add_tracer hardcodes the scalar 3D layout (atmosphere_process.hpp:375-393) and the advected "tracers" group is a monolithic (col, Q, lev) array consumed by HOMME. A standalone (col, cmp, lev) field would never be advected. Per-component scalar tracers get HOMME + turbulence advection for free; grouping preserves multi-component indexing. Investigation phase confirms whether a monolithic subgroup allocation is feasible or whether the process iterates over named subfields.
    date: 2026-06-10

  - decision: Four field families — qv_iso, qc_iso, qi_iso, qr_iso — mirroring P3 mass tracers (qv, qc, qi, qr)
    rationale: EAMxx names cloud liquid qc (not ql). qm is ice rime mass (not "mixed-phase"); rime isotope content is treated as part of qi for now and revisited at P3 integration if rime-specific isotope tracking is needed.
    date: 2026-06-10

  - decision: Defer separate q*_trace (non-fractionating tracer) array family
    rationale: Per-component runtime naming supports adding a parallel _trace family later without structural change; isotopes-first keeps this spec bounded. Original campaign instructions list both; revisit at physics-integration campaign.
    date: 2026-06-10

  - decision: Use generic _iso names (NOT qv_HDO)
    rationale: Component index maps to species (HDO, H2-18O, ...) via configuration; allows future expansion without field name changes
    date: 2026-06-09

  - decision: Runtime-configurable component count via tracer_count parameter (defaults block added in spec 001)
    rationale: Number of tracers controlled via YAML/atmchange; allows 1-N tracers without recompilation
    date: 2026-06-09

  - decision: Set tracer_count to 1 for initial testing
    rationale: Single tracer simplifies validation; multi-tracer testing deferred to integration phase
    date: 2026-06-09

  - decision: Primary validation is BFB vs pre-extension baseline plus a standalone single-process test (replaces earlier mass-conservation criterion)
    rationale: Fields are zero-initialized with no physics, so any conservation test passes trivially and proves nothing. The real risk of this change is perturbing existing physics through tracer-group changes; cprnc BFB vs baseline b58a0fbea6 is the proof (regression-baseline convention). History/restart I/O in EAMxx is configured via YAML output streams, not C++ calls — the earlier add_to_history/add_to_restart criterion referenced an API that does not exist.
    date: 2026-06-10
```

## Context

This spec builds on spec 001 (process registration) to define the multi-component tracer storage required for water isotope tracking in EAMxx. Water isotopes (HDO, H2-18O, HTO, etc.) must track through the water mass phases P3 carries as tracers: vapor (qv), cloud liquid (qc), cloud ice (qi), and rain (qr). (P3 also carries qm = ice rime mass and number concentrations; isotopic rime content is folded into qi for now.)

Key architectural constraint (verified on this branch): the advected "tracers" group in EAMxx is monolithically allocated as (col, Q, lev) and consumed by the dycore; `add_tracer` only accepts scalar 3D fields. Therefore each isotope component is registered as its own scalar tracer — gaining HOMME and turbulence advection automatically — and the per-family multi-component view is assembled via field groups. This satisfies the campaign requirement of (col, cmp, lev)-indexed access with runtime-configurable cmp while remaining advection-correct.

EAMxx I/O: history output is configured through YAML output streams (OutputManager); restart handling is automatic for fields required by registered processes. No per-field C++ "add to history/restart" call exists.

The primary tier-0 validation is the regression-baseline contract: with tracer_count=1, zero-valued tracers, and no-op physics, a 5-day F2010-SCREAMv1 run must be bit-for-bit identical to the pre-extension baseline (b58a0fbea6). Any non-BFB result means the structural change leaked into existing physics. Later specs that extend bulk arrays inherit the slice-1-unchanged contract (bulk-H2O slice BFB vs pre-extension baseline).

## Approach

### Phase 1: Investigation - Tracer Group and Field Group Mechanics

**Subagent goals:**
- Confirm how per-component scalar tracers join the "tracers" group and how the dycore picks them up (GroupRequest, MonolithicAlloc)
- Determine whether a per-family subgroup (e.g., group "qv_iso") can obtain a contiguous (col, cmp, lev) view, or whether iteration over named subfields is the supported pattern
- Identify metadata attachment (units, long_name) for tracer fields
- Confirm tracer-count effects on dycore state sizing (qsize) and any namelist implications

**Deliverable:** Investigation report with the exact registration code pattern and a go/no-go on monolithic subgroup views.

### Phase 2: Investigation - I/O Streams and Restart

**Subagent goals:**
- Document YAML output-stream configuration needed for the new fields (history)
- Confirm restart behavior for runtime-named tracer fields
- Identify constraints on field names in output (runtime-generated names vs fixed names)

**Deliverable:** Investigation report; informs output.yaml in the standalone test.

### Phase 3: Field Registry Implementation

1. Create `field_registry.cpp`:
   - Helper generating component tracer names (`qv_iso_<k>`, `qc_iso_<k>`, `qi_iso_<k>`, `qr_iso_<k>` for k in [0, tracer_count))
   - Register each via `add_tracer<Updated>` with units kg/kg and per-family group membership
   - Attach long_name metadata per field

2. Update `CMakeLists.txt` to compile field_registry.cpp

### Phase 4: Update Process Interface

1. `eamxx_water_tracers_process_interface.hpp`: members for tracer_count and field/group handles
2. `eamxx_water_tracers_process_interface.cpp`:
   - Read tracer_count from ParameterList (registered in defaults by spec 001)
   - `create_requests()`: call field registry helper; skip entirely when tracer_count == 0
   - `initialize_impl()`: retrieve fields, deep_copy(0), stamp timestamps
   - `run_impl()`: no-op with comment (advection by dynamics; fractionation in later campaign)

### Phase 5: Standalone Single-Process Test

1. Create `tests/single-process/water_tracers/` mirroring the p3 test layout:
   - input.yaml with `water_tracers: tracer_count: 1` in atm_procs
   - output.yaml streaming the new fields
   - CMakeLists.txt registering the ctest
2. Test asserts: process instantiates via factory, fields allocate with expected extents, output file contains the fields

### Phase 6: Build + BFB Verification (supported machine)

1. Create F2010-SCREAMv1 / ne4pg2_ne4pg2 case; `atmchange atm_procs_list+=water_tracers` and `atmchange water_tracers::tracer_count=1`
2. `./case.build` (SC6)
3. Generate/locate baseline from b58a0fbea6 per regression-baseline procedure (record BASELINE.txt manifest: parent SHA, date, compset, resolution, stop_n/stop_option)
4. Run 5 days (`./xmlchange STOP_OPTION=ndays,STOP_N=5 && ./case.submit`)
5. `cprnc` final history file vs baseline; require BFB (SC5)
6. If non-BFB: investigate tracer-group sizing effects on dycore before any tolerance-based fallback; a non-BFB result here means the change is functional, not structural

## Ask-Before Items

```yaml
ask_before:
  - Adding new dependencies to CMakeLists.txt beyond field manager/share
  - Modifying existing field dimensions or registration for standard water fields (qv, qc, qi, qr, qm)
  - Any change that alters dycore qsize handling for runs without water_tracers enabled
  - Running simulations longer than 5 days for validation
  - Falling back to a tolerance criterion if BFB fails (requires proxy-BFB justification per regression-baseline convention)
```

## Checkpoint Strategy

```yaml
checkpoints:
  - phase: Phase 1-2 complete (Investigation)
    validation: |
      - Tracer-group registration pattern confirmed with code references
      - Monolithic subgroup view feasibility resolved (go/no-go)
      - I/O stream + restart approach documented

  - phase: Phase 3-4 complete
    validation: |
      - field_registry.cpp exists; CMakeLists updated
      - tracer_count plumbed from ParameterList; tracer_count==0 is a clean no-op
      - Code compiles

  - phase: Phase 5 complete
    validation: |
      - Standalone test exists and passes locally (SC3)

  - phase: Phase 6 complete
    validation: |
      - Case builds (SC6); 5-day run BFB vs baseline (SC5)
      - All success criteria pass
```

## Notes

- Each component is a real advected tracer; "cmp" indexing is provided by per-family groups, not a bespoke (col, cmp, lev) field outside the tracer group
- tracer_count=1 for initial testing; multi-tracer validation deferred to integration specs
- BFB vs baseline b58a0fbea6 is the primary tier-0 gate; conservation closure tests start only when physics writes nonzero values (see tracer-conservation tolerance defaults for that phase)
- Fractionation physics and microphysics coupling are out of scope for this spec
- CIME gates require a supported machine (pm-cpu or equivalent)
