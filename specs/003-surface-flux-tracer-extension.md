# Spec: Extend Surface Evaporation Flux to Multi-Component Tracer Field

## Metadata

```yaml
spec_id: 003-surface-flux-tracer-extension
title: Extend Surface Evaporation Flux to Multi-Component Tracer Field
author: rfiorella
created: 2026-06-09
updated: 2026-06-10
status: draft
spec_type: model-e3sm-eamxx
campaign_id: water-isotope-infrastructure
priority: high
estimated_effort_hours: 6
execution_mode: checkpoint
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
    description: Requires the tracer_count parameter plumbing (namelist defaults block) from spec 001
```

## Model-Specific Configuration

```yaml
model_specific:
  target_compset: F2010-SCREAMv1
  target_resolution: ne4pg2_ne4pg2
  validation_tier: 0  # compile + BFB regression vs baseline (5-day run)
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
  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: Surface coupling importer; surf_evap is scalar2d (line ~44). Multi-component 2D precedents already exist here — surf_mom_flux is vector2d (line ~45) and dstflx is vector4d (line ~61). Import loop only writes fields named by the coupler via SCDataManager (setup_surface_coupling_data, lines ~67-107).
    format: C++ source
    required: true

  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.hpp
    description: Header for surface coupling importer
    format: C++ header
    required: true

  - path: components/eamxx/src/share/field/field_identifier.hpp
    description: Field identifier and metadata structures for understanding field layouts
    format: C++ header
    required: true

  - path: components/eamxx/cime_config/namelist_defaults_eamxx.xml
    description: tracer_count parameter source (water_tracers block from spec 001); sc_import block where importer params live
    format: XML
    required: true
```

## Deliverables

```yaml
deliverables:
  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: Updated to define surf_evap_iso with (col, cmp) dimensions, registered only when tracer_count > 0; zero-initialized with timestamp set in initialize_impl
    format: C++ source

  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.hpp
    description: Updated to declare surf_evap_iso field handle and tracer_count member
    format: C++ header
```

## Success Criteria

```yaml
success_criteria:
  - id: SC1
    phase: implementation
    description: surf_evap_iso field is requested with (col, cmp) dimensions, guarded on tracer_count > 0
    type: shell
    command: |
      grep -q "surf_evap_iso" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp && \
      grep -q "tracer_count" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    failure_mode: Field not created with correct dimensions or registered unconditionally (perturbing every EAMxx run)
    note: Cheap textual pre-check; SC3/SC4 are the real gates

  - id: SC2
    phase: implementation
    description: surf_evap_iso has correct units (kg/m2/s) consistent with surf_evap
    type: shell
    command: |
      grep -B 2 -A 2 "surf_evap_iso" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp | grep -q "kg/m2/s"
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    failure_mode: Missing or incorrect units cause I/O errors or physical inconsistencies

  - id: SC3
    phase: validation
    description: EAMxx builds successfully with surf_evap_iso field
    type: shell
    command: |
      cd ${CASE_DIR} && ./case.build
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.hpp
    failure_mode: Compilation errors from field dimension mismatch or missing declarations

  - id: SC4
    phase: validation
    description: BFB regression — importer change leaves a standard run (water tracers disabled, tracer_count=0) bit-for-bit identical to baseline b58a0fbea6; a second run with tracer_count=1 is also BFB since surf_evap_iso stays zero and is never consumed
    type: comparison
    command: |
      cd ${CASE_DIR} && \
      ./xmlchange STOP_OPTION=ndays,STOP_N=5 && \
      ./case.submit && \
      cprnc ${BASELINE_DIR}/h0.final.nc ${RUNDIR}/h0.final.nc
    expect: BFB
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    failure_mode: Importer change perturbs import indexing or field timestamps for existing fields; the importer runs in EVERY EAMxx case, so this is the campaign's largest BFB risk
```

## Out of Scope

```yaml
out_of_scope:
  - Modification of coupler interface or data transfer between components (cpl field lists, scream_cpl_indices Fortran side)
  - Implementation of isotope-specific surface flux physics
  - Changes to land model (ELM) boundary conditions
  - Extension of other surface fluxes (sensible heat, momentum, etc.)
  - Coupling to ocean or sea ice isotope tracers
  - Conservation closure tests (field is zero-valued; deferred to coupling integration)
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: Use field name surf_evap_iso (not surf_evap_HDO)
    rationale: Generic naming consistent with qv_iso pattern from spec 002; allows future multi-isotope extension
    date: 2026-06-09

  - decision: Dimensions are (col, cmp) following the existing vector2d pattern in the importer (surf_mom_flux); dstflx (vector4d) shows the multi-component precedent
    rationale: Surface fluxes are 2D; the importer already supports vector-valued 2D fields — reuse that pattern rather than invent one
    date: 2026-06-10

  - decision: Register surf_evap_iso only when tracer_count > 0; zero-init and set timestamp in initialize_impl
    rationale: The importer is instantiated in every EAMxx run. The import loop only writes fields the coupler names via SCDataManager, so surf_evap_iso would otherwise be a Computed field that is never written — risking field-timestamp postcondition failures and breaking BFB for all users. Guarding on the parameter keeps existing runs untouched.
    date: 2026-06-10

  - decision: tracer_count single source is the water_tracers namelist block from spec 001; the importer reads the same parameter (plumbing path resolved in Phase 1 investigation)
    rationale: Two independent copies of the count (process vs importer) would drift; one authoritative parameter
    date: 2026-06-10

  - decision: Hard dependency only on spec 001, not spec 002
    rationale: Surface flux extension is independent of 3D tracer definitions; only requires the tracer_count parameter plumbing
    date: 2026-06-09

  - decision: EAMxx-only changes, no coupler modifications
    rationale: Initial implementation holds zeros; real isotope fluxes require coupler + land/ocean work deferred to integration phase. Until then downstream processes can reference the field but values are placeholders.
    date: 2026-06-09

  - decision: Primary validation is BFB vs baseline (replaces earlier "flux conservation" criterion)
    rationale: A zero-valued, never-consumed flux trivially "conserves"; the earlier criterion also referenced a verify_flux_conservation.py script that does not exist. The meaningful tier-0 proof is that the importer change is invisible to existing runs (regression-baseline convention).
    date: 2026-06-10
```

## Context

This spec extends the surface evaporation boundary condition to support multi-component water isotope tracers. The current `surf_evap` field is a scalar 2D field (`add_field<Computed>("surf_evap", scalar2d, kg/m2/s, ...)`, importer line ~44). For water isotope tracking, `surf_evap_iso` adds a component dimension to hold isotope-specific evaporation fluxes.

In the full isotope-enabled model, surface evaporation will include equilibrium fractionation, kinetic fractionation, and source-specific isotopic signatures (ocean vs. land vs. vegetation). This spec only defines the field structure; fractionation physics and coupler transport are deferred.

Importer architecture (verified on this branch):
- The importer registers candidate fields in `create_requests()`; the actual import set is dictated by the coupler through `SCDataManager` (field names, cpl indices, vector components — `setup_surface_coupling_data`).
- Fields not named by the coupler are never written by the import loop. A `Computed` field that is never written can fail field-initialization/timestamp checks — hence the explicit zero-init + timestamp in `initialize_impl` and the `tracer_count > 0` guard.
- Multi-component 2D fields already exist here: `surf_mom_flux` (vector2d) and `dstflx` (vector4d). `surf_evap_iso` follows the same layout pattern with extent = tracer_count.

The real coupling integration point — adding isotope fluxes to the coupler field list and the Fortran-side surface-coupling indices — is explicitly out of scope and lands in a later campaign together with land/ocean isotope sources.

## Approach

### Phase 1: Examine Existing Pattern and Parameter Plumbing

1. Document how `surf_evap` and the vector-valued fields (`surf_mom_flux`, `dstflx`) are registered in `create_requests()`
2. Trace how the importer's ParameterList is populated (sc_import block in namelist defaults) and determine the cleanest path for it to read `tracer_count` from the water_tracers block (options: duplicate param under sc_import set by buildnml, or query via shared driver params) — record the chosen mechanism
3. Confirm behavior of Computed-but-never-written fields (timestamp/property checks) to validate the zero-init requirement

### Phase 2: Update Header File

1. `atmosphere_surface_coupling_importer.hpp`:
   - Add `m_tracer_count` member (default 0)
   - Add field handle for surf_evap_iso
   - Follow naming/style of existing flux fields

### Phase 3: Implement Multi-Component Surface Flux Field

1. `atmosphere_surface_coupling_importer.cpp`:
   - Read tracer_count in constructor (default 0 — absent parameter must behave exactly as today)
   - In `create_requests()`, when tracer_count > 0: `add_field<Computed>("surf_evap_iso", <(COL,CMP) layout, extent=tracer_count>, kg/m2/s, grid_name)` following the vector2d pattern
2. In `initialize_impl()`:
   - When registered: deep_copy(0) and update the field timestamp so postcondition checks pass
   - Comment: "Isotope-specific surface fluxes will be provided by coupler in future work"

### Phase 4: Build Verification (supported machine)

1. F2010-SCREAMv1 / ne4pg2_ne4pg2 case; build with no water-tracer settings (tracer_count absent ⇒ 0 ⇒ field not registered)
2. Second configuration with `atmchange water_tracers::tracer_count=1`; verify surf_evap_iso appears in field manager with extents (ncols, 1)
3. Both builds compile and link clean (SC3)

### Phase 5: BFB Regression

1. Use the baseline from b58a0fbea6 (shared with spec 002; BASELINE.txt manifest per regression-baseline procedure)
2. Run 5 days with tracer_count=0 (default path); cprnc vs baseline — must be BFB (this is the gate protecting every existing EAMxx user)
3. Run 5 days with tracer_count=1; cprnc vs baseline — expected BFB since the field is zero and unconsumed; investigate any diff before proceeding
4. Document results (SC4)

## Ask-Before Items

```yaml
ask_before:
  - Modifying coupler interface or data structures (SCDataManager, Fortran cpl indices)
  - Adding new dependencies beyond field manager
  - Creating new YAML/namelist parameters beyond reusing tracer_count
  - Registering surf_evap_iso unconditionally (without the tracer_count guard)
  - Running simulations longer than 5 days for validation
```

## Checkpoint Strategy

```yaml
checkpoints:
  - phase: Phase 1 complete
    validation: |
      - Vector-field registration pattern documented
      - tracer_count plumbing path into importer decided and recorded
      - Never-written-Computed-field behavior confirmed

  - phase: Phase 2-3 complete
    validation: |
      - Header + implementation updated; guard on tracer_count present
      - Code compiles

  - phase: Phase 4 complete
    validation: |
      - Both configurations build (SC3)
      - surf_evap_iso present with correct extents when enabled; absent when disabled

  - phase: Phase 5 complete
    validation: |
      - BFB vs baseline in both configurations (SC4)
      - All success criteria pass
```

## Notes

- Independent of spec 002 (3D tracers); only shares the tracer_count parameter from spec 001
- The importer runs in every EAMxx case — the tracer_count guard plus BFB gate is what makes this change safe to merge
- Real flux values require coupler + surface-component work (later campaign); until then the field is a zero-valued placeholder with correct shape and metadata
- CIME gates require a supported machine (pm-cpu or equivalent)
