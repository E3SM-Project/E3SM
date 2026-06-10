# Spec: Multi-Component Water Tracer Field Definitions in EAMxx

## Metadata

```yaml
spec_id: 002-multi-component-tracer-fields
title: Multi-Component Water Tracer Field Definitions in EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-09
status: draft
spec_type: model-e3sm-eamxx
campaign_id: water-isotope-infrastructure
priority: high
estimated_effort_hours: 12
execution_mode: checkpoint
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
    description: Requires water_tracers and water_isotopes processes to be registered
```

## Model-Specific Configuration

```yaml
model_specific:
  target_compset: F2010-SCREAMv1
  target_resolution: ne4pg2_ne4pg2
  validation_tier: 0  # compile/lint only
  requires_input_data: false
  platform: any supported EAMxx platform
  baseline_tag: b58a0fbea6
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
  
  - path: components/eamxx/src/share/field/field_manager.cpp
    description: Field manager implementation
    format: C++ source
    required: true
  
  - path: components/eamxx/src/share/field/field_identifier.hpp
    description: Field identifier and metadata structures
    format: C++ header
    required: true
  
  - path: components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    description: Example of field requests and multi-dimensional field handling
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
```

## Deliverables

```yaml
deliverables:
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.hpp
    description: Updated process interface with field requests for multi-component tracers
    format: C++ header
    
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    description: Updated with field registration in create_requests() for qv_iso, qi_iso, ql_iso, qr_iso, qm_iso
    format: C++ source
    
  - path: components/eamxx/src/physics/water_tracers/field_registry.cpp
    description: Field registration logic for multi-component tracers with metadata
    format: C++ source
    
  - path: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    description: Updated to include field_registry.cpp in build
    format: CMake
```

## Success Criteria

```yaml
success_criteria:
  - id: SC1
    phase: implementation
    description: All five isotope tracer fields are requested with correct dimensions
    type: shell
    command: |
      grep -E "(qv_iso|qi_iso|ql_iso|qr_iso|qm_iso)" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp | wc -l
    assertion: count >= 5
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    failure_mode: Missing field requests mean tracers cannot be allocated or tracked
    
  - id: SC2
    phase: implementation
    description: Fields have (col, cmp, lev) dimensions with cmp runtime-configurable
    type: shell
    command: |
      grep -q "FieldLayout.*col.*cmp.*lev" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp || \
      grep -q "FieldLayout.*COL.*CMP.*LEV" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    failure_mode: Incorrect dimensions prevent proper indexing and cause runtime array bounds errors
    
  - id: SC3
    phase: validation
    description: Mass conservation test for each tracer component passes within tolerance
    type: tolerance
    command: |
      cd ${CASE_DIR} && \
      ./xmlchange --id RUN_LENGTH --val "5 days" && \
      ./case.submit && \
      python3 ${SOURCE_ROOT}/components/eamxx/scripts/verify_mass_conservation.py \
        --field qv_iso --rtol 1e-12 --atol 1e-14 ${RUNDIR}/output.nc
    rtol: 1.0e-12
    atol: 1.0e-14
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Mass leaks indicate tracer not properly advected or numerics unstable
    
  - id: SC4
    phase: implementation
    description: Field metadata includes units and long_name for all isotope fields
    type: shell
    command: |
      for field in qv_iso qi_iso ql_iso qr_iso qm_iso; do
        grep -q "units.*kg/kg" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/field_registry.cpp || exit 1
        grep -q "long_name" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/field_registry.cpp || exit 1
      done
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/field_registry.cpp
    failure_mode: Missing metadata causes I/O errors or unintelligible output files
    
  - id: SC5
    phase: implementation
    description: Fields are added to history and restart I/O streams
    type: shell
    command: |
      grep -E "(add_to_history|add_to_restart)" ${SOURCE_ROOT}/components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp | wc -l
    assertion: count >= 5
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    failure_mode: Fields not written to output; cannot verify results or restart simulations
    
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
  - Multi-tracer advection schemes (handled by dynamics)
  - Coupling to microphysics parameterizations (P3/SHOC integration - separate spec)
  - Boundary condition specification for isotope tracers
  - Input data preparation or initial condition files
  - Performance optimization or GPU-specific memory layouts
  - Validation beyond mass conservation (tier 1+ testing deferred)
  - YAML configuration beyond tracer_count parameter
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: Use field names qv_iso, qi_iso, ql_iso, qr_iso, qm_iso (NOT qv_HDO)
    rationale: Generic naming allows future expansion to other isotope species without field name changes
    date: 2026-06-09
    
  - decision: Dimensions are (col, cmp, lev) with runtime-configurable cmp
    rationale: Component dimension size controlled via YAML parameter; allows 1-N tracers per field
    date: 2026-06-09
    
  - decision: Set tracer_count to 1 for initial testing
    rationale: Single tracer simplifies validation; multi-tracer testing deferred to integration phase
    date: 2026-06-09
    
  - decision: Use investigation phase with subagents
    rationale: Field manager architecture and multi-component handling patterns need exploration before implementation
    date: 2026-06-09
    
  - decision: Include metadata (units, long_name) for all fields
    rationale: Required for CF-compliant I/O and output file usability
    date: 2026-06-09
    
  - decision: Add fields to history and restart I/O streams
    rationale: Necessary for verification and simulation continuity
    date: 2026-06-09
    
  - decision: No boundary conditions in this spec
    rationale: Boundary handling deferred to physics integration; initial work uses zero initialization
    date: 2026-06-09
    
  - decision: Mass conservation tolerance rtol=1e-12, atol=1e-14
    rationale: Standard EAMxx defaults for tracer conservation tests
    date: 2026-06-09
```

## Context

This spec builds on spec 001 (process registration) to define the multi-component field arrays required for water isotope tracking in EAMxx. Water isotopes (HDO, H2-18O, HTO, etc.) must track through all water phases: vapor (qv), ice (qi), liquid (ql), rain (qr), and mixed-phase (qm).

Unlike standard scalar tracers, each water phase field requires a "component" dimension to hold multiple isotope species simultaneously. This architecture allows runtime configuration of the number of isotope tracers without recompilation, supporting both single-isotope (e.g., HDO only) and multi-isotope (HDO + H2-18O) configurations.

The field names use generic `_iso` suffixes rather than species-specific names to enable future extension. The component index will map to specific isotope species (HDO=1, H2-18O=2, etc.) via a separate configuration mechanism.

EAMxx field registration requires:
1. Requesting fields in `create_requests()` with `FieldRequest` objects specifying dimensions and metadata
2. Allocation happens automatically via the field manager after `create_requests()` returns
3. Fields are accessed in `initialize_impl()` and `run_impl()` via the grids manager

This spec includes an investigation phase to explore:
- Field manager API for multi-dimensional field requests
- YAML parameter handling for runtime configuration (tracer_count)
- Pattern for grouping related fields with shared dimensions
- I/O stream registration for history and restart output

Mass conservation is the primary validation criterion. Each tracer component must conserve mass within machine precision (rtol=1e-12, atol=1e-14) over a 5-day simulation.

## Approach

### Phase 1: Investigation - Field Manager Architecture

**Subagent goals:**
- Understand `FieldRequest` API for multi-dimensional fields
- Identify how to specify runtime-configurable dimension sizes
- Find examples of (col, cmp, lev) or similar layouts in existing code
- Determine metadata attachment mechanism (units, long_name)

**Deliverable:** Investigation report documenting:
- Code patterns for field request with 3D layout
- YAML parameter plumbing for tracer_count
- Metadata registration approach

### Phase 2: Investigation - I/O Stream Integration

**Subagent goals:**
- Find how fields are added to history output streams
- Find restart I/O registration mechanism
- Identify any special handling for multi-component fields in I/O

**Deliverable:** Investigation report documenting:
- I/O stream registration API
- Any dimension ordering constraints for output
- Restart file field requirements

### Phase 3: Field Registry Implementation

1. Create `field_registry.cpp` with:
   - Helper function to register all five isotope fields
   - Metadata specification: units="kg/kg", appropriate long_name for each
   - Dimension specifications with runtime cmp size
   
2. Add to `CMakeLists.txt`:
   - Include field_registry.cpp in library source list
   - Verify no additional dependencies needed

### Phase 4: Update Process Interface - Field Requests

1. Edit `eamxx_water_tracers_process_interface.hpp`:
   - Add member variables for tracer_count parameter
   - Add member variables to hold field references
   
2. Edit `eamxx_water_tracers_process_interface.cpp`:
   - In constructor, read tracer_count from ParameterList
   - In `create_requests()`, call field registry helper
   - Request all five fields: qv_iso, qi_iso, ql_iso, qr_iso, qm_iso
   - Specify FieldLayout as (COL, CMP, LEV) with CMP size = tracer_count
   - Add fields to history and restart I/O streams

### Phase 5: Initialization Logic

1. In `initialize_impl()`:
   - Retrieve field references from grids manager
   - Initialize all fields to zero (simple starting condition)
   - Add comment: "Initial conditions will be set by CAM infrastructure"
   
2. In `run_impl()`:
   - Add comment: "Tracer advection handled by dynamics; fractionation physics to be added in integration spec"
   - No-op implementation for now

### Phase 6: Build Verification

1. Create test case with F2010-SCREAMv1 / ne4pg2_ne4pg2
2. Add to YAML configuration:
   ```yaml
   water_tracers:
     tracer_count: 1
   ```
3. Build with water_tracers enabled
4. Verify compilation succeeds
5. Check field dimensions in build artifacts or debug output

### Phase 7: Mass Conservation Validation

1. Run 5-day simulation
2. Write or adapt Python script to verify mass conservation:
   - Load qv_iso from output netCDF
   - Compute global integral at each timestep
   - Check integral change < rtol * mean_mass + atol
3. Verify all five fields pass conservation test
4. Document results

## Ask-Before Items

```yaml
ask_before:
  - Adding new dependencies to CMakeLists.txt beyond field manager
  - Modifying existing field dimensions for standard water fields
  - Creating YAML configuration template files
  - Running simulations longer than 5 days for validation
```

## Checkpoint Strategy

```yaml
checkpoints:
  - phase: Phase 1-2 complete (Investigation)
    validation: |
      - Investigation reports document field manager API
      - YAML parameter handling approach identified
      - I/O stream registration approach identified
    
  - phase: Phase 3 complete
    validation: |
      - field_registry.cpp exists
      - CMakeLists.txt updated
      - Code compiles in isolation
    
  - phase: Phase 4 complete
    validation: |
      - Process interface updated with field requests
      - Grep finds all five field names in create_requests()
      - Code compiles
    
  - phase: Phase 5 complete
    validation: |
      - initialize_impl() retrieves and initializes fields
      - run_impl() contains appropriate comments
      - Code compiles
    
  - phase: Phase 6 complete
    validation: |
      - Test case builds successfully
      - Success criteria SC1, SC2, SC4, SC5, SC6 pass
    
  - phase: Phase 7 complete
    validation: |
      - Mass conservation test (SC3) passes
      - All success criteria pass
```

## Notes

- Field dimensions follow (col, cmp, lev) convention where cmp is the component/tracer index
- tracer_count=1 for initial testing; multi-tracer validation deferred to integration specs
- Mass conservation test is the primary validation for tier 0
- Fractionation physics and microphysics coupling are out of scope for this spec
- Field initialization uses zero values; proper initial conditions will be set by CAM infrastructure in later work
