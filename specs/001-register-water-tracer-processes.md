# Spec: Register Water Tracer Processes in EAMxx

## Metadata

```yaml
spec_id: 001-register-water-tracer-processes
title: Register Water Tracer Processes in EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-09
status: draft
spec_type: model-e3sm-eamxx
campaign_id: water-isotope-infrastructure
priority: high
estimated_effort_hours: 8
execution_mode: checkpoint
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
  - path: components/eamxx/src/physics/register_physics.hpp
    description: Registry for physics processes in EAMxx
    format: C++ header
    required: true
  
  - path: components/eamxx/src/physics/p3/eamxx_p3_process_interface.hpp
    description: Example process interface for pattern matching
    format: C++ header
    required: true
  
  - path: components/eamxx/src/physics/p3/eamxx_p3_process_interface.cpp
    description: Example process implementation for pattern matching
    format: C++ source
    required: true
  
  - path: components/eamxx/src/physics/p3/CMakeLists.txt
    description: Example CMake configuration for physics process
    format: CMake
    required: true
  
  - path: components/eamxx/src/physics/water_tracers/
    description: Directory that already contains some implementation work
    format: directory
    required: true
```

## Deliverables

```yaml
deliverables:
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.hpp
    description: Process interface header for water_tracers process
    format: C++ header
    
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_tracers_process_interface.cpp
    description: Process implementation for water_tracers process
    format: C++ source
    
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_isotopes_process_interface.hpp
    description: Process interface header for water_isotopes process (subset of water_tracers)
    format: C++ header
    
  - path: components/eamxx/src/physics/water_tracers/eamxx_water_isotopes_process_interface.cpp
    description: Process implementation for water_isotopes process
    format: C++ source
    
  - path: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    description: CMake configuration to build water_tracers physics process
    format: CMake
    
  - path: components/eamxx/src/physics/register_physics.hpp
    description: Updated to register water_tracers and water_isotopes processes
    format: C++ header
```

## Success Criteria

```yaml
success_criteria:
  - id: SC1
    phase: implementation
    description: CMake build completes without errors for water_tracers physics library
    type: shell
    command: |
      cd ${BUILD_DIR}/components/eamxx && \
      cmake --build . --target eamxx_physics_water_tracers -- -j4
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    failure_mode: Build system cannot find headers, malformed CMake syntax, or missing dependencies
    
  - id: SC2
    phase: implementation
    description: water_tracers process is registered in AtmosphereProcessFactory
    type: shell
    command: |
      grep -q "water_tracers" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp && \
      grep -q "EAMXX_HAS_WATER_TRACERS" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/register_physics.hpp
    failure_mode: Process not added to registry, cannot be used in yaml configuration files
    
  - id: SC3
    phase: implementation
    description: water_isotopes process is registered in AtmosphereProcessFactory
    type: shell
    command: |
      grep -q "water_isotopes" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp && \
      grep -q "EAMXX_HAS_WATER_ISOTOPES" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/register_physics.hpp
    failure_mode: Process not added to registry, isotope-specific fractionation physics cannot be invoked
    
  - id: SC4
    phase: validation
    description: Full EAMxx build completes with water_tracers enabled
    type: shell
    command: |
      cd ${CASE_DIR} && ./case.build
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/register_physics.hpp
      deliverable: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    failure_mode: Linker errors, missing symbols, or header inclusion issues
```

## Out of Scope

```yaml
out_of_scope:
  - Full implementation of water tracer physics (separate spec)
  - Definition of water tracer field arrays (qv_trace, qi_trace, etc.) - separate spec
  - Integration with P3/SHOC microphysics kernels (already partially implemented)
  - YAML configuration file setup for runtime parameters
  - Modification of any process other than water_tracers and water_isotopes
  - Changes to dynamics, coupler, or other components outside eamxx/src/physics
  - Performance optimization or GPU-specific tuning
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: water_isotopes is a separate process from water_tracers
    rationale: Isotopes require special fractionation physics; keeping separate allows independent enabling/disabling
    date: 2026-06-09
    
  - decision: Follow existing EAMxx process registration pattern
    rationale: Consistency with p3, shoc, mam processes; use AtmosphereProcess base class
    date: 2026-06-09
    
  - decision: Both processes live in components/eamxx/src/physics/water_tracers/ directory
    rationale: Shared implementation code for tracer handling; directory already exists with partial implementation
    date: 2026-06-09
    
  - decision: Validation tier 0 (compile/lint only)
    rationale: Initial infrastructure work; full physics testing deferred to integration specs
    date: 2026-06-09
```

## Context

The EAMxx component of E3SM requires infrastructure to track water isotopes and general water tracers through the atmosphere. This work is the first step in a campaign to add full water isotope capability to EAMxx following conventions used in other ESMs.

Water tracers are additional water species (e.g., tagged water from different source regions) that track through the model but do not undergo fractionation. Water isotopes (HDO, H2-18O, HTO) are a special case of water tracers that undergo equilibrium and kinetic fractionation during phase changes.

The `components/eamxx/src/physics/water_tracers/` directory already contains implementation work for fractionation physics and kinetic fractionation hooks, but the processes are not yet registered in the EAMxx process factory. This spec focuses solely on registering the processes so they can be invoked by the atmosphere driver.

EAMxx uses a factory pattern to register physics processes. Each process:
1. Inherits from `AtmosphereProcess` base class
2. Implements `create_requests()`, `initialize_impl()`, and `run_impl()` methods
3. Is registered in `register_physics.hpp` with a preprocessor guard (e.g., `EAMXX_HAS_WATER_TRACERS`)
4. Has a corresponding CMake target that can be optionally linked

This spec does NOT include defining the field arrays (qv_trace, qi_trace, etc.) or YAML configuration - those are separate tasks in the campaign.

## Approach

### Phase 1: Examine Existing Process Pattern

1. Read `components/eamxx/src/physics/p3/eamxx_p3_process_interface.hpp` and `.cpp` to understand the process interface pattern
2. Read `components/eamxx/src/physics/p3/CMakeLists.txt` to understand build configuration
3. Identify required methods: constructor, `create_requests()`, `initialize_impl()`, `run_impl()`
4. Note the pattern for requesting fields from grids_manager

### Phase 2: Create Water Tracers Process

1. Create `eamxx_water_tracers_process_interface.hpp`:
   - Define `WaterTracers` class inheriting from `AtmosphereProcess`
   - Declare required methods
   - Add minimal member variables (grids, number of columns/levels)
   
2. Create `eamxx_water_tracers_process_interface.cpp`:
   - Implement constructor taking `ekat::Comm` and `ekat::ParameterList`
   - Implement `create_requests()` - initially empty or minimal
   - Implement `initialize_impl()` - initially no-op
   - Implement `run_impl()` - initially no-op with comment "TODO: tracer transport"
   
3. Add namespace: use `scream` namespace consistent with other processes

### Phase 3: Create Water Isotopes Process

1. Create `eamxx_water_isotopes_process_interface.hpp`:
   - Define `WaterIsotopes` class inheriting from `AtmosphereProcess`
   - Similar structure to `WaterTracers` but separate class
   
2. Create `eamxx_water_isotopes_process_interface.cpp`:
   - Similar stub implementation to `WaterTracers`
   - Add comment indicating this will apply fractionation physics

### Phase 4: Create CMake Configuration

1. Create `components/eamxx/src/physics/water_tracers/CMakeLists.txt`:
   - Define library target `eamxx_physics_water_tracers`
   - Add source files: both `.cpp` files
   - Link dependencies: `ekat`, `share`, any existing water_tracers utilities
   - Use pattern from `p3/CMakeLists.txt`
   
2. Verify no preprocessor guards preventing build

### Phase 5: Register Processes

1. Edit `components/eamxx/src/physics/register_physics.hpp`:
   - Add `#ifdef EAMXX_HAS_WATER_TRACERS` block with include for water_tracers header
   - Add `#ifdef EAMXX_HAS_WATER_ISOTOPES` block with include for water_isotopes header
   - In `register_physics()` function, add registration calls:
     - `proc_factory.register_product("water_tracers", &create_atmosphere_process<WaterTracers>);`
     - `proc_factory.register_product("water_isotopes", &create_atmosphere_process<WaterIsotopes>);`
   - Follow alphabetical ordering convention if present

### Phase 6: Build Verification

1. Configure a test case with EAMxx enabled
2. Build with `-DEAMXX_HAS_WATER_TRACERS=ON -DEAMXX_HAS_WATER_ISOTOPES=ON`
3. Verify no compilation or linking errors
4. Check that symbols for both processes exist in build artifacts

## Ask-Before Items

```yaml
ask_before:
  - Modifying any process registration file outside water_tracers directory
  - Adding new dependencies beyond what P3/SHOC already use
  - Changing the process naming convention (currently "water_tracers" and "water_isotopes")
```

## Checkpoint Strategy

```yaml
checkpoints:
  - phase: Phase 2 complete
    validation: |
      - WaterTracers .hpp and .cpp exist
      - Files compile in isolation
    
  - phase: Phase 3 complete
    validation: |
      - WaterIsotopes .hpp and .cpp exist
      - Files compile in isolation
    
  - phase: Phase 4 complete
    validation: |
      - CMakeLists.txt exists
      - cmake configure succeeds
    
  - phase: Phase 5 complete
    validation: |
      - register_physics.hpp updated
      - grep finds both process names
    
  - phase: Phase 6 complete
    validation: |
      - Full build succeeds
      - Success criteria SC1-SC4 pass
```

## Notes

- The water_tracers directory already contains fractionation physics implementation; this spec only adds process registration
- Field definitions (qv_trace, qi_trace dimensions) deferred to spec 002
- YAML configuration deferred to spec 003
- Integration testing with actual tracer transport deferred to post-campaign validation
