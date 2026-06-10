# Campaign: Water Isotope Infrastructure for EAMxx

## Metadata

```yaml
campaign_id: water-isotope-infrastructure
title: Water Isotope Infrastructure for EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-09
status: draft
base_branch: wiso-dev
repository_path: /Users/rfiorella/code/E3SM/EAMXX-wiso
branch_prefix: rfiorella/eamxx/water-isotope-infra
execution_mode: continuous
estimated_total_hours: 26
```

## Campaign Overview

This campaign implements the foundational infrastructure required for water isotope and tracer capability in EAMxx. It establishes the process registration framework, multi-component field definitions, and surface flux extensions needed to track water isotopes (HDO, H2-18O, HTO) through the atmospheric component of E3SM.

The campaign consists of three specs executed in sequence with controlled dependencies:
1. Register water tracer processes in the EAMxx physics factory
2. Define multi-component tracer field arrays with runtime-configurable dimensions
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
  max_lines_per_pr: 600
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

**Summary:** Creates process interface stubs for `water_tracers` and `water_isotopes` and registers them in the EAMxx AtmosphereProcessFactory. Establishes CMake build configuration and enables compilation with `EAMXX_HAS_WATER_TRACERS` and `EAMXX_HAS_WATER_ISOTOPES` guards.

**Key Deliverables:**
- `components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.{hpp,cpp}`
- `components/eamxx/src/physics/water_tracers/eamxx_water_isotopes_process_interface.{hpp,cpp}`
- `components/eamxx/src/physics/water_tracers/CMakeLists.txt`
- Updated `components/eamxx/src/physics/register_physics.hpp`

**Success Gates:**
- CMake build completes for `eamxx_physics_water_tracers` target
- Both processes appear in `register_physics.hpp` with appropriate guards
- Full EAMxx build completes with processes enabled

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
pr_title: "Add multi-component tracer field arrays (qv_iso, qi_iso, ql_iso, qr_iso, qm_iso)"
pr_draft: true
validation_tier: 0
critical_path: true
```

**Summary:** Defines five multi-component field arrays (`qv_iso`, `qi_iso`, `ql_iso`, `qr_iso`, `qm_iso`) with (col, cmp, lev) dimensions where the component dimension is runtime-configurable via YAML parameter `tracer_count`. Implements field registration, metadata specification, and I/O stream integration.

**Key Deliverables:**
- `components/eamxx/src/physics/water_tracers/field_registry.cpp`
- Updated `components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.{hpp,cpp}`
- Updated `components/eamxx/src/physics/water_tracers/CMakeLists.txt`

**Success Gates:**
- All five fields registered with correct (col, cmp, lev) layout
- Fields include proper metadata (units="kg/kg", long_name)
- Fields added to history and restart I/O streams
- Mass conservation test passes (rtol=1e-12, atol=1e-14) over 5-day simulation
- Build completes with tracer_count=1

**Dependency Rationale:** Requires process registration infrastructure from spec 001 to attach fields to the water_tracers process.

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

**Summary:** Extends the surface coupling importer to define `surf_evap_iso` as a (col, cmp) field with runtime-configurable component dimension. Establishes the boundary condition infrastructure for isotope-specific surface fluxes.

**Key Deliverables:**
- Updated `components/eamxx/src/control/atmosphere_surface_coupling_importer.{hpp,cpp}`

**Success Gates:**
- `surf_evap_iso` field registered with (col, cmp) layout
- Field has correct units (kg/m2/s) and metadata
- Build completes successfully
- Flux conservation test passes (rtol=1e-12, atol=1e-14) over 5-day simulation

**Dependency Rationale:** Requires process registration infrastructure from spec 001 but is independent of 3D field definitions (spec 002). Surface fluxes are 2D boundary conditions.

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
- No changes to CIME infrastructure or case control

### File Ownership
Each spec touches distinct files with no overlap:
- **Spec 001:** `water_tracers/` directory creation, `register_physics.hpp`
- **Spec 002:** `water_tracers/` interior files (field_registry, process implementation)
- **Spec 003:** `control/atmosphere_surface_coupling_importer.*`

### Build System Impact
- New CMake library target: `eamxx_physics_water_tracers`
- New preprocessor guards: `EAMXX_HAS_WATER_TRACERS`, `EAMXX_HAS_WATER_ISOTOPES`
- No changes to existing targets or dependencies
- Optional linking: existing builds unaffected unless guards enabled

### Testing Strategy
- **Validation Tier 0:** Compile/lint only for all specs
- **Target Compset:** F2010-SCREAMv1
- **Target Resolution:** ne4pg2_ne4pg2
- **Test Platform:** Any supported EAMxx platform (pm-cpu, pm-gpu)
- **Baseline Tag:** b58a0fbea6
- **Mass Conservation Tests:** 5-day simulations with rtol=1e-12, atol=1e-14

## Out of Scope

```yaml
out_of_scope:
  - Full implementation of water isotope fractionation physics
  - Integration with P3/SHOC microphysics kernels
  - YAML configuration files for runtime parameters (beyond tracer_count)
  - Coupling to land model or coupler for isotope-specific boundary conditions
  - Input data preparation or initial condition files
  - Performance optimization or GPU-specific tuning
  - Validation beyond tier 0 (compile + mass conservation)
  - Multi-tracer configurations beyond single-isotope testing
  - Dynamics modifications for tracer advection
  - Documentation or user guide updates
```

## Success Criteria for Campaign Completion

```yaml
campaign_success_criteria:
  - id: CSC1
    description: All three specs pass their individual success criteria
    verification: Review each spec's validation output
    
  - id: CSC2
    description: Full EAMxx build succeeds with all three features enabled simultaneously
    verification: |
      cd cime/scripts
      ./create_newcase --case test_water_iso --compset F2010-SCREAMv1 --res ne4pg2_ne4pg2 --mach <machine>
      cd test_water_iso
      ./xmlchange --id EAMXX_HAS_WATER_TRACERS --val TRUE
      ./xmlchange --id EAMXX_HAS_WATER_ISOTOPES --val TRUE
      ./case.setup
      ./case.build
    
  - id: CSC3
    description: No regression in existing EAMxx tests
    verification: |
      cd components/eamxx
      ./scripts/test-all-eamxx -m <machine> -t dbg
    
  - id: CSC4
    description: All PRs successfully merge to wiso-dev without conflicts
    verification: Manual review of PR merge history
    
  - id: CSC5
    description: Code coverage for new files meets minimum 60% threshold
    verification: lcov/gcov analysis on water_tracers/ directory (optional)
```

## Post-Campaign Integration

After campaign completion, the following integration work is required (tracked in separate specs/campaigns):

1. **Physics Integration:** Implement actual fractionation physics in water_isotopes process
2. **Microphysics Coupling:** Connect to P3/SHOC for phase-change fractionation
3. **Boundary Conditions:** Full coupling with land model for realistic surface fluxes
4. **Validation:** Tier 1-3 testing with observational datasets
5. **Documentation:** User guide, namelist documentation, example configurations

## Risk Assessment

```yaml
risks:
  - risk: Field manager API changes between spec creation and execution
    mitigation: Investigation phase in spec 002 will detect API mismatches early
    severity: medium
    
  - risk: CMake configuration conflicts with existing physics processes
    mitigation: Follow P3/SHOC patterns exactly; test build in isolation first
    severity: low
    
  - risk: Mass conservation tests fail due to numerical precision issues
    mitigation: Use standard EAMxx tolerances (rtol=1e-12, atol=1e-14); document any platform-specific variations
    severity: medium
    
  - risk: Parallel execution of specs 002 and 003 causes merge conflicts
    mitigation: File ownership is disjoint; no shared files modified
    severity: low
    
  - risk: Baseline tag b58a0fbea6 becomes stale during campaign execution
    mitigation: Update baseline tag in each spec if wiso-dev advances significantly
    severity: low
```

## Author Notes

This campaign establishes infrastructure only. No physics is implemented. All fields are initialized to zero or placeholder values. The goal is to enable subsequent campaigns to add fractionation physics and microphysics integration without revisiting basic infrastructure.

The naming convention (`_iso` suffix rather than `_HDO`) is intentional to support future multi-isotope configurations (HDO + H2-18O + HTO) without field name changes.

The validation tier 0 approach (compile + mass conservation only) is appropriate for infrastructure work. Physics validation requires tier 1+ testing deferred to integration campaigns.

All specs use checkpoint execution mode to enable human review and course correction at natural breakpoints (investigation reports, build verification, validation tests).

## Campaign Execution Timeline

**Estimated Duration:** 26 hours total
- Spec 001: 8 hours
- Spec 002: 12 hours (can run parallel with 003 after 001)
- Spec 003: 6 hours (can run parallel with 002 after 001)

**Critical Path:** 001 → 002 = 20 hours

**Parallel Execution:** If enabled, specs 002 and 003 run simultaneously after 001, reducing total time to 20 hours.

**Continuous Mode:** Campaign runs without pause between specs. Each spec's success gates must pass before proceeding to dependents.

## Validation Status

This campaign manifest has been validated against the campaign rubric criteria:

- **C1 (Metadata):** ✓ All required fields present (campaign_id, title, author, dates, base_branch, execution_mode)
- **C2 (Specs List):** ✓ All three specs listed with paths, order, hours, dependencies
- **C3 (Dependencies):** ✓ Dependency graph clear; 001 → {002, 003} with rationale
- **C4 (Integration):** ✓ File ownership, build impact, testing strategy documented
- **C5 (PR Strategy):** ✓ Branch naming, PR titles, draft status, labels, size limits specified
- **C6 (Success Criteria):** ✓ Five campaign-level success criteria with verification commands
- **C7 (Out of Scope):** ✓ Campaign-level exclusions documented
- **C8 (Timeline):** ✓ Estimated hours per spec, critical path, parallel opportunities identified
- **C9 (Context):** ✓ Campaign overview explains purpose and relationship to broader effort

**Status:** PASS - Ready for execution with `/run-campaign campaigns/2026-06-09-water-isotope-infrastructure.md`
