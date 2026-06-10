# Spec: Register Water Tracer Processes in EAMxx

## Metadata

```yaml
spec_id: 001-register-water-tracer-processes
title: Register Water Tracer Processes in EAMxx
author: rfiorella
created: 2026-06-09
updated: 2026-06-10
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
  validation_tier: 0  # compile/lint only; no model runs required by this spec
  requires_input_data: false
  platform: pm-cpu or other supported E3SM machine (local macOS cannot run CIME builds)
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

  - path: components/eamxx/src/physics/CMakeLists.txt
    description: Parent CMake file that adds physics subdirectories and owns the eamxx_physics INTERFACE target
    format: CMake
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
    description: Example CMake configuration for physics process (see target_compile_definitions PUBLIC EAMXX_HAS_P3 at line ~81 and conditional link into eamxx_physics at lines ~115-117)
    format: CMake
    required: true

  - path: components/eamxx/cime_config/namelist_defaults_eamxx.xml
    description: Namelist defaults; every runnable atm process needs a defaults block here so atmchange/buildnml accept its parameters
    format: XML
    required: true
```

Note: `components/eamxx/src/physics/water_tracers/` does **not** exist yet.
This spec creates it from scratch. (An earlier draft of this spec claimed the
directory already contained fractionation physics; that was incorrect.)

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
    description: Process interface header for water_isotopes process (subclass of WaterTracers)
    format: C++ header

  - path: components/eamxx/src/physics/water_tracers/eamxx_water_isotopes_process_interface.cpp
    description: Process implementation for water_isotopes process
    format: C++ source

  - path: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    description: CMake configuration to build water_tracers physics library, define EAMXX_HAS_WATER_TRACERS/EAMXX_HAS_WATER_ISOTOPES, and link into eamxx_physics
    format: CMake

  - path: components/eamxx/src/physics/CMakeLists.txt
    description: Updated with add_subdirectory(water_tracers)
    format: CMake

  - path: components/eamxx/src/physics/register_physics.hpp
    description: Updated to register water_tracers and water_isotopes processes
    format: C++ header

  - path: components/eamxx/cime_config/namelist_defaults_eamxx.xml
    description: New <water_tracers> and <water_isotopes> defaults blocks (with tracer_count parameter) so atmchange can add the processes to atm_procs_list; default atm_procs_list is NOT changed
    format: XML
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
      cmake --build . --target water_tracers -- -j4
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/water_tracers/CMakeLists.txt
    failure_mode: Build system cannot find headers, malformed CMake syntax, or missing dependencies

  - id: SC2
    phase: implementation
    description: water_tracers process is registered in AtmosphereProcessFactory behind its guard
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
    description: water_isotopes process is registered in AtmosphereProcessFactory behind its guard
    type: shell
    command: |
      grep -q "water_isotopes" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp && \
      grep -q "EAMXX_HAS_WATER_ISOTOPES" ${SOURCE_ROOT}/components/eamxx/src/physics/register_physics.hpp
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/register_physics.hpp
    failure_mode: Process not added to registry, isotope-specific fractionation physics cannot be invoked

  - id: SC4
    phase: implementation
    description: Namelist defaults blocks exist so buildnml/atmchange accept the new processes
    type: shell
    command: |
      grep -q "<water_tracers" ${SOURCE_ROOT}/components/eamxx/cime_config/namelist_defaults_eamxx.xml && \
      grep -q "tracer_count" ${SOURCE_ROOT}/components/eamxx/cime_config/namelist_defaults_eamxx.xml
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/cime_config/namelist_defaults_eamxx.xml
    failure_mode: atmchange rejects water_tracers parameters; process cannot be added to the atm DAG at runtime

  - id: SC5
    phase: validation
    description: Full EAMxx build completes with water_tracers library linked
    type: shell
    command: |
      cd ${CASE_DIR} && ./case.build
    assertion: exit_code == 0
    verifies:
      deliverable: components/eamxx/src/physics/register_physics.hpp
      deliverable: components/eamxx/src/physics/water_tracers/CMakeLists.txt
      deliverable: components/eamxx/src/physics/CMakeLists.txt
    failure_mode: Linker errors, missing symbols, or header inclusion issues
```

## Out of Scope

```yaml
out_of_scope:
  - Full implementation of water tracer physics (separate spec)
  - Definition of water tracer field arrays - separate spec (002)
  - Integration with P3/SHOC microphysics kernels
  - Adding water_tracers/water_isotopes to the default atm_procs_list (user opt-in via atmchange only)
  - Modification of any process other than water_tracers and water_isotopes
  - Changes to dynamics, coupler, or other components outside eamxx
  - Performance optimization or GPU-specific tuning
```

## Resolved Decisions

```yaml
resolved_decisions:
  - decision: WaterIsotopes subclasses WaterTracers (not a sibling process)
    rationale: Campaign instructions require isotopes to be a special case of water tracers; subclassing shares all field handling and lets isotopes override only fractionation hooks. Both names still registered separately in the factory so either can be enabled alone.
    date: 2026-06-10

  - decision: EAMXX_HAS_WATER_TRACERS / EAMXX_HAS_WATER_ISOTOPES are PUBLIC compile definitions set by the library's own CMakeLists (P3 pattern, p3/CMakeLists.txt:81), not user-facing CMake options or CIME XML variables
    rationale: Matches every existing EAMxx physics process; the guard becomes active automatically when the library is linked into eamxx_physics
    date: 2026-06-10

  - decision: Runtime enablement is via atm_procs_list (atmchange), requiring a defaults block in namelist_defaults_eamxx.xml; this spec adds the defaults block but does not change the default process list
    rationale: Factory registration alone never places a process in the execution DAG; the namelist defaults block is the missing plumbing. Keeping the default list unchanged guarantees no impact on existing runs.
    date: 2026-06-10

  - decision: Follow existing EAMxx process registration pattern
    rationale: Consistency with p3, shoc, mam processes; use AtmosphereProcess base class
    date: 2026-06-09

  - decision: Both processes live in components/eamxx/src/physics/water_tracers/ directory (created by this spec)
    rationale: Shared implementation code for tracer handling
    date: 2026-06-09

  - decision: Validation tier 0 (compile/lint only)
    rationale: This spec is build-infrastructure only; nothing runs because the processes are not in atm_procs_list by default. BFB regression proof against baseline b58a0fbea6 is performed in specs 002/003 where runtime behavior could change.
    date: 2026-06-10
```

## Context

The EAMxx component of E3SM requires infrastructure to track water isotopes and general water tracers through the atmosphere. This work is the first step in a campaign to add full water isotope capability to EAMxx following conventions used in other ESMs.

Water tracers are additional water species (e.g., tagged water from different source regions) that track through the model but do not undergo fractionation. Water isotopes (HDO, H2-18O, HTO) are a special case of water tracers that undergo equilibrium and kinetic fractionation during phase changes. The class design mirrors this: `WaterIsotopes` derives from `WaterTracers`.

EAMxx process plumbing (verified against this branch):
1. Each process inherits from `AtmosphereProcess` and implements `create_requests()` (pure virtual, called by `set_grids()` — see `share/atm_process/atmosphere_process.hpp:108-112`), `initialize_impl()`, and `run_impl()`, plus `name()` and `type()` (returning `AtmosphereProcessType::Physics`).
2. The process library's CMakeLists sets `target_compile_definitions(<lib> PUBLIC EAMXX_HAS_<NAME>)` and conditionally links itself into the `eamxx_physics` INTERFACE target (P3 pattern).
3. `register_physics.hpp` includes the header and calls `proc_factory.register_product(...)` inside the matching `#ifdef` guard.
4. To actually run, the process name must appear in `atm_procs_list`, which requires a defaults block in `cime_config/namelist_defaults_eamxx.xml`. Without that block, `atmchange` rejects the process and its parameters.

This spec covers all four layers but leaves the default `atm_procs_list` untouched; enabling the processes is a per-case user action.

This spec does NOT include defining the field arrays - those are a separate task in the campaign (spec 002).

## Approach

### Phase 1: Examine Existing Process Pattern

1. Read `components/eamxx/src/physics/p3/eamxx_p3_process_interface.hpp` and `.cpp` to understand the process interface pattern
2. Read `components/eamxx/src/physics/p3/CMakeLists.txt`, noting `target_compile_definitions(p3 PUBLIC EAMXX_HAS_P3)` and the `if (TARGET eamxx_physics) ... target_link_libraries(eamxx_physics INTERFACE p3)` block
3. Identify required overrides: constructor, `name()`, `type()`, `create_requests()`, `initialize_impl()`, `run_impl()`
4. Read the `<p3>` and `<sc_import>` blocks in `namelist_defaults_eamxx.xml` to mirror the defaults-block structure

### Phase 2: Create Water Tracers Process

1. Create `eamxx_water_tracers_process_interface.hpp`:
   - Define `WaterTracers` class inheriting from `AtmosphereProcess`
   - Declare `name()`, `type()` (Physics), `create_requests()`, `initialize_impl()`, `run_impl()`
   - Member hooks that subclasses may override (e.g., a virtual no-op fractionation hook) — keep minimal
   - Add minimal member variables (grid pointer, ncols, nlevs, tracer_count)

2. Create `eamxx_water_tracers_process_interface.cpp`:
   - Implement constructor taking `ekat::Comm` and `ekat::ParameterList`; read `tracer_count` (default 0)
   - Implement `create_requests()` - initially empty (fields arrive in spec 002)
   - Implement `initialize_impl()` - initially no-op
   - Implement `run_impl()` - initially no-op with comment "TODO: tracer physics (spec 002+)"

3. Use `scream` namespace consistent with other processes

### Phase 3: Create Water Isotopes Process

1. Create `eamxx_water_isotopes_process_interface.hpp`:
   - Define `WaterIsotopes` class inheriting from `WaterTracers`
   - Override `name()`; inherit all field handling
   - Declare override points for fractionation hooks (no-op for now)

2. Create `eamxx_water_isotopes_process_interface.cpp`:
   - Stub implementation; comment indicating fractionation physics lands in a later campaign

### Phase 4: Create CMake Configuration

1. Create `components/eamxx/src/physics/water_tracers/CMakeLists.txt`:
   - Define library target `water_tracers` with both `.cpp` files
   - `target_compile_definitions(water_tracers PUBLIC EAMXX_HAS_WATER_TRACERS EAMXX_HAS_WATER_ISOTOPES)`
   - Link `scream_share` (and `eamxx_physics_share` if needed), following `p3/CMakeLists.txt`
   - `if (TARGET eamxx_physics) target_link_libraries(eamxx_physics INTERFACE water_tracers) endif()`

2. Edit `components/eamxx/src/physics/CMakeLists.txt`: add `add_subdirectory(water_tracers)`

### Phase 5: Register Processes

1. Edit `components/eamxx/src/physics/register_physics.hpp`:
   - Add `#ifdef EAMXX_HAS_WATER_TRACERS` block with include + `proc_factory.register_product("water_tracers", &create_atmosphere_process<WaterTracers>);`
   - Add `#ifdef EAMXX_HAS_WATER_ISOTOPES` block with include + `proc_factory.register_product("water_isotopes", &create_atmosphere_process<WaterIsotopes>);`

### Phase 6: Namelist Defaults Block

1. Edit `cime_config/namelist_defaults_eamxx.xml`:
   - Add `<water_tracers inherit="atm_proc_base">` block with `<tracer_count type="integer" doc="...">0</tracer_count>`
   - Add `<water_isotopes inherit="atm_proc_base">` block (inherits tracer_count semantics or references the water_tracers block; investigation during implementation decides the cleanest inheritance)
   - Do NOT modify any `<atm_procs_list>` entries

2. Verify with buildnml dry-run or `atmchange atm_procs_list+=water_tracers` on a test case (on a supported machine)

### Phase 7: Build Verification

1. On a supported machine, create an F2010-SCREAMv1 / ne4pg2_ne4pg2 case
2. `./case.setup && ./case.build` — library builds and links by default (guards become defined; processes registered but inert since not in atm_procs_list)
3. Check that symbols for both processes exist in the built library

## Ask-Before Items

```yaml
ask_before:
  - Modifying any process registration file outside water_tracers directory (register_physics.hpp and physics/CMakeLists.txt edits above are pre-approved)
  - Adding new dependencies beyond what P3/SHOC already use
  - Changing the process naming convention (currently "water_tracers" and "water_isotopes")
  - Adding either process to any default atm_procs_list
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
      - WaterIsotopes .hpp and .cpp exist; class derives from WaterTracers
      - Files compile in isolation

  - phase: Phase 4 complete
    validation: |
      - CMakeLists.txt exists; add_subdirectory added
      - cmake configure succeeds

  - phase: Phase 5-6 complete
    validation: |
      - register_physics.hpp updated; grep finds both process names
      - namelist_defaults_eamxx.xml has water_tracers block with tracer_count

  - phase: Phase 7 complete
    validation: |
      - Full build succeeds on supported machine
      - Success criteria SC1-SC5 pass
```

## Notes

- The water_tracers directory is created from scratch by this spec; no prior implementation exists on wiso-dev
- Field definitions deferred to spec 002; surface flux extension to spec 003
- Because the default atm_procs_list is unchanged and run_impl is a no-op, this spec cannot change model answers; build gates are sufficient (BFB proof happens in specs 002/003)
- CIME build gates (SC5) require a supported machine; plan execution on pm-cpu or equivalent
