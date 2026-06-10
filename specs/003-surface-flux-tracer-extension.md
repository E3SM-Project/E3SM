# Spec: Extend Surface Evaporation Flux to Multi-Component Tracer Field

## Metadata

```yaml
spec_id: 003-surface-flux-tracer-extension
title: Extend Surface Evaporation Flux to Multi-Component Tracer Field
author: rfiorella
created: 2026-06-09
updated: 2026-06-09
status: draft
spec_type: model-e3sm-eamxx
campaign_id: water-isotope-infrastructure
priority: high
estimated_effort_hours: 6
execution_mode: checkpoint
dependencies:
  - spec_id: 001-register-water-tracer-processes
    relationship: hard
    description: Requires water_tracers process registration infrastructure
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
  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: Surface coupling importer where surf_evap is currently defined as scalar2d
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
  
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    description: Reference for multi-component field pattern (from spec 001)
    format: C++ source
    required: true
```

## Deliverables

```yaml
deliverables:
  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    description: Updated to define surf_evap_iso with (col, cmp) dimensions
    format: C++ source
    
  - path: components/eamxx/src/control/atmosphere_surface_coupling_importer.hpp
    description: Updated to declare surf_evap_iso field and related member variables
    format: C++ header
```

## Success Criteria

```yaml
success_criteria:
  - id: SC1
    phase: implementation
    description: surf_evap_iso field is requested with (col, cmp) dimensions
    type: shell
    command: |
      grep -q "surf_evap_iso" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp && \
      grep -E "(col.*cmp|COL.*CMP)" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    failure_mode: Field not created with correct dimensions; cannot store per-component surface fluxes
    
  - id: SC2
    phase: implementation
    description: surf_evap_iso has correct units and metadata
    type: shell
    command: |
      grep -A 5 "surf_evap_iso" ${SOURCE_ROOT}/components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp | grep -q "kg/m2/s"
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
    description: Mass conservation test for surf_evap_iso within tolerance
    type: tolerance
    command: |
      cd ${CASE_DIR} && \
      ./xmlchange --id RUN_LENGTH --val "5 days" && \
      ./case.submit && \
      python3 ${SOURCE_ROOT}/components/eamxx/scripts/verify_flux_conservation.py \
        --field surf_evap_iso --rtol 1e-12 --atol 1e-14 ${RUNDIR}/output.nc
    rtol: 1.0e-12
    atol: 1.0e-14
    verifies:
      deliverable: components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp
    failure_mode: Flux accumulation errors or incorrect dimension indexing cause mass imbalance
```

## Out of Scope

```yaml
out_of_scope:
  - Modification of coupler interface or data transfer between components
  - Implementation of isotope-specific surface flux physics
  - Changes to land model (ELM) boundary conditions
  - Extension of other surface fluxes (sensible heat, momentum, etc.)
  - Validation beyond compilation and mass conservation (tier 1+ deferred)
  - Coupling to ocean or sea ice isotope tracers
  - YAML configuration for surface flux parameters
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: Use field name surf_evap_iso (not surf_evap_HDO)
    rationale: Generic naming consistent with qv_iso pattern from spec 002; allows future multi-isotope extension
    date: 2026-06-09
    
  - decision: Dimensions are (col, cmp) not (col, cmp, lev)
    rationale: Surface fluxes are 2D quantities at the surface boundary; no vertical structure
    date: 2026-06-09
    
  - decision: Runtime-configurable cmp dimension via tracer_count parameter
    rationale: Consistent with field definitions in spec 002; enables flexible isotope configurations
    date: 2026-06-09
    
  - decision: Hard dependency only on spec 001, not spec 002
    rationale: Surface flux extension is independent of 3D field definitions; only requires process registration infrastructure
    date: 2026-06-09
    
  - decision: EAMxx-only changes, no coupler modifications
    rationale: Initial implementation uses zero or dummy values; full coupling deferred to integration phase
    date: 2026-06-09
    
  - decision: Validation tier 0 with mass conservation test
    rationale: Same validation approach as specs 001-002; tier 1+ requires full physics implementation
    date: 2026-06-09
    
  - decision: Mass conservation tolerance rtol=1e-12, atol=1e-14
    rationale: Standard EAMxx defaults for flux conservation tests
    date: 2026-06-09
```

## Context

This spec extends the surface evaporation flux field to support multi-component water isotope tracers. The current `surf_evap` field is a scalar 2D field (columns only) that represents total surface evaporation flux. For water isotope tracking, we need `surf_evap_iso` with an additional component dimension to hold isotope-specific evaporation fluxes.

Surface evaporation is a critical boundary condition for atmospheric water tracers. In the full isotope-enabled model, surface evaporation will include:
1. Equilibrium fractionation (vapor pressure effects)
2. Kinetic fractionation (diffusion-limited transport)
3. Source-specific isotopic signatures (ocean vs. land vs. vegetation)

However, this spec focuses solely on defining the multi-component field structure in the atmosphere surface coupling importer. The physics of isotopic fractionation during evaporation is deferred to later integration specs.

The `atmosphere_surface_coupling_importer` is responsible for importing surface boundary conditions from the coupler into EAMxx's field manager. Currently it handles standard fluxes like `surf_evap`, sensible heat, momentum, etc. This spec adds `surf_evap_iso` alongside the existing `surf_evap` field.

Key architectural points:
- Surface fluxes are 2D (column-only), not 3D like atmospheric fields
- The component dimension (cmp) must be runtime-configurable via the same tracer_count parameter used in spec 002
- Initial implementation will use zero or placeholder values since full coupling requires land model changes (out of scope)
- The field must be registered in the field manager so downstream processes can access it

## Approach

### Phase 1: Examine Existing Surface Flux Pattern

1. Read `components/eamxx/src/control/atmosphere_surface_coupling_importer.cpp` to understand:
   - How `surf_evap` is currently defined and allocated
   - Field request pattern for 2D surface fluxes
   - Metadata specification (units, long_name)
   - Integration with field manager
   
2. Read `components/eamxx/src/control/atmosphere_surface_coupling_importer.hpp` to identify:
   - Member variables holding field references
   - Constructor and initialization patterns
   - Any YAML parameter handling

3. Document the current architecture for creating and registering surface fluxes

### Phase 2: Update Header File

1. Edit `atmosphere_surface_coupling_importer.hpp`:
   - Add member variable for tracer_count parameter
   - Add member variable to hold surf_evap_iso field reference
   - Add declaration for surf_evap_iso in the class interface
   - Follow naming and style conventions from existing flux fields

### Phase 3: Implement Multi-Component Surface Flux Field

1. Edit `atmosphere_surface_coupling_importer.cpp`:
   - In constructor, read tracer_count from ParameterList (if not already available)
   - Add field request for surf_evap_iso with:
     - FieldLayout: (COL, CMP) where CMP size = tracer_count
     - Units: "kg/m2/s" (consistent with surf_evap)
     - Long name: "Surface evaporation flux for water isotope tracers"
   - Follow the pattern used for existing surface flux fields
   
2. Add initialization logic:
   - Retrieve field reference from field manager
   - Initialize to zero (placeholder until coupling implemented)
   - Add comment: "Isotope-specific surface fluxes will be provided by coupler in future work"

3. Ensure field is properly registered so downstream processes can access it

### Phase 4: Build Verification

1. Create test case with F2010-SCREAMv1 / ne4pg2_ne4pg2
2. Add to YAML configuration:
   ```yaml
   water_tracers:
     tracer_count: 1
   ```
3. Build with water_tracers enabled
4. Verify:
   - Compilation succeeds
   - surf_evap_iso field appears in field manager
   - Dimensions are correct (columns x 1)
   - No linker errors

### Phase 5: Mass Conservation Validation

1. Run 5-day simulation with surf_evap_iso initialized to zero
2. Write or adapt Python script to verify conservation:
   - Load surf_evap_iso from history output
   - Verify values remain zero or stable
   - Check no spurious flux accumulation
   - Verify dimensions match expected (ncols, tracer_count)
3. Document results and confirm SC4 passes

## Ask-Before Items

```yaml
ask_before:
  - Modifying coupler interface or data structures
  - Adding new dependencies beyond field manager
  - Creating new YAML configuration parameters beyond tracer_count
  - Running simulations longer than 5 days for validation
```

## Checkpoint Strategy

```yaml
checkpoints:
  - phase: Phase 1 complete
    validation: |
      - Current surface flux architecture documented
      - Field request pattern identified
      - Metadata requirements understood
    
  - phase: Phase 2 complete
    validation: |
      - Header file updated with surf_evap_iso declarations
      - Code compiles (header only)
    
  - phase: Phase 3 complete
    validation: |
      - Implementation complete in .cpp file
      - grep finds surf_evap_iso in field requests
      - Code compiles
    
  - phase: Phase 4 complete
    validation: |
      - Test case builds successfully
      - Success criteria SC1, SC2, SC3 pass
    
  - phase: Phase 5 complete
    validation: |
      - Mass conservation test (SC4) passes
      - All success criteria pass
```

## Notes

- This spec is independent of spec 002 (3D field definitions) - surface fluxes are 2D only
- Hard dependency on spec 001 for process registration infrastructure
- No coupler modifications in this spec; full coupling deferred to integration phase
- Initial values are zero/placeholder; proper coupling requires land model work
- Field naming follows generic "_iso" pattern for future multi-isotope extension
- Mass conservation test with rtol=1e-12, atol=1e-14 matches specs 001-002 validation approach
