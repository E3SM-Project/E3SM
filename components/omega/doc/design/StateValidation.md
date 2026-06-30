<!--- OMEGA Ocean State Validation requirements and design ------------------->

(omega-design-state-validation)=

# Ocean State Validation

## 1 Overview

The Ocean State Validation module provides a mechanism for checking the
physical plausibility of the ocean prognostic state and selected auxiliary
variables at runtime. After each timestep (or at user-defined intervals) the
model can call `validateOceanState` to scan every owned cell and vertical layer
for NaN (Not-a-Number) values and values that lie outside a pre-defined
physically meaningful range. Any detected anomaly is reported through the
OMEGA logging infrastructure as a critical error message.

Whether the run is aborted on detection is controlled by two runtime
configuration options under `Omega::State::Validation`:

- **`AbortOnNan`** (default: `true`) — when `true`, the presence of any NaN
  values in a checked field causes the run to be terminated via `MPI_Abort`
  (with a full stack backtrace printed first).
- **`AbortOnOOB`** (default: `false`) — when `true`, out-of-bounds values
  trigger the same hard abort.

When both flags are `false` and errors are detected, a critical-level log
message is emitted and the run continues rather than being terminated.

## 2 Requirements

### 2.1 Requirement: Check for NaN values

The validation function must detect NaN values in the ocean state fields.
NaNs indicate a numerical instability or a programming error, and their
presence in the prognostic state makes continued time-stepping meaningless.
Every element of each validated array must be tested.

### 2.2 Requirement: Check for out-of-bounds values

In addition to NaN detection, the validation function must check that every
element of each validated field lies within a physically plausible range.
Values outside these ranges indicate catastrophic model failure and should
halt the simulation before corrupted output is written to disk.

### 2.3 Requirement: Validate PseudoThickness

`PseudoThickness` must be validated for each owned cell over all vertical
layers. The valid range is $[10^{-10},\, 1000]$ m.
Negative or near-zero layer thicknesses indicate numerical collapse of the
column and must be caught immediately.

### 2.4 Requirement: Validate KineticEnergyCell

`KineticEnergyCell` from the kinematic auxiliary state must be validated for
each owned cell over all vertical layers. The valid range is $[0,\, 10]$
m$^2$ s$^{-2}$.  Negative kinetic energies are unphysical, and values
exceeding 10 m$^2$ s$^{-2}$ correspond to current speeds above
$\sim 4.5$ m s$^{-1}$, which are unrealistic for the open ocean.

### 2.5 Requirement: Validate Temperature tracer

Ocean Conservative Temperature must be validated for each owned cell over
all vertical layers. The valid range is $[-10,\, 50]$ °C.
This broad range accommodates all realistic oceanographic regimes including
polar and hydrothermal vent environments.

### 2.6 Requirement: Validate Salinity tracer

Ocean Absolute Salinity must be validated for each owned cell over all
vertical layers. The valid range is $[-2,\, 60]$ g kg$^{-1}$.
Values below $-2$ g kg$^{-1}$ are unphysical and values above 60 g kg$^{-1}$
are outside the valid domain of the TEOS-10 equation of state.

### 2.7 Requirement: GPU/CPU portability

All validation kernels must execute on both CPU and GPU hardware using the
Kokkos parallel programming model and therefore must be expressed as Kokkos
parallel reductions.

### 2.8 Requirement: Informative error reporting

On detection of any failure the module must log a critical-level message that
identifies the field name, the nature of the problem (NaN or out-of-bounds),
and the number of offending elements.

### 2.9 Requirement: Graceful handling of absent tracers

If the Temperature or Salinity tracer is not present in the tracer registry
(e.g. in configurations that do not activate active tracers) the
corresponding check must be skipped silently rather than causing an error.

### 2.10 Requirement: Configurable abort-on-error behavior

The user must be able to independently control whether NaN errors and
out-of-bounds errors cause the run to abort.  Two boolean YAML options
(`AbortOnNan` and `AbortOnOOB`, both nested under
`Omega::State::Validation`) govern this:

- When `AbortOnNan: true`, any detected NaN value triggers a hard abort via
  `MPI_Abort` (with a stack backtrace printed beforehand).
- When `AbortOnOOB: true`, any detected out-of-bounds value triggers the same
  hard abort.
- When both are `false`, errors are logged at critical severity and the run
  continues.

The defaults are `AbortOnNan: true` and `AbortOnOOB: false`.

### 2.11 Desired: Configurable valid ranges

In the future it may be desirable to allow the user to override the default
valid ranges through the OMEGA configuration system (e.g. for idealised
process studies that intentionally use non-oceanic parameter values).

### 2.12 Desired: Configurable validation frequency

In the future it may be desirable to allow the user to control whether
validation is performed every timestep, every N timesteps, or only at
specific points in the run (e.g. after restart reads).

## 3 Algorithmic Formulation

No complex numerical algorithms are required. Each field is checked with two
independent `parallelReduce` passes over the domain, restricted to active
cells where `CellMask(i, k) > 0` (inactive cells, such as land cells, are
skipped):

1. **NaN pass** – counts elements for which `Kokkos::isnan(val)` is `true`.
2. **Bounds pass** – counts elements that are finite yet lie outside
   $[\text{MinVal},\, \text{MaxVal}]$.

Separating the two passes avoids potentially undefined behaviour when
comparing NaN values with `<` or `>`.

For a 2-D field $f_{i,k}$ with $i \in [0,\, N_\text{cells})$ and
$k \in [0,\, N_\text{vert})$ the two counts are:

$$
N_\text{NaN}    = \sum_{i,k} M_{i,k}\,\mathbf{1}[\,\text{isnan}(f_{i,k})\,]
$$

$$
N_\text{OOB}    = \sum_{i,k} M_{i,k}\,\mathbf{1}[\,\lnot\,\text{isnan}(f_{i,k})
                   \land (f_{i,k} < f_\text{min} \lor f_{i,k} > f_\text{max})\,]
$$

where $M_{i,k}$ is `CellMask(i, k)` (1 for active cell-layer, 0 for
inactive).

For a 3-D tracer array $T_{n,i,k}$ the same expressions are applied at a
fixed tracer index $n$.

The validation traverses only the `NCellsOwned` cells (excluding halo cells)
to avoid double-counting in the parallel decomposition.

## 4 Design

The module is implemented as a free function (`validateOceanState`) plus two
file-local helper functions. It does not introduce a class or persistent state.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

The valid ranges for each field are compile-time constants embedded in the
implementation:

| Field                | MinVal    | MaxVal |
|----------------------|-----------|--------|
| `PseudoThickness`    | 1×10⁻¹⁰   | 1000   |
| `KineticEnergyCell`  | 0         | 10     |
| `Temperature`        | −10       | 50     |
| `Salinity`           | −2        | 60     |

The abort-on-error behavior is governed by two runtime YAML options nested
under `Omega::State::Validation`:

| YAML key      | Type | Default | Meaning                                      |
|---------------|------|---------|----------------------------------------------|
| `AbortOnNan`  | bool | `true`  | Abort via `MPI_Abort` when NaNs are detected |
| `AbortOnOOB`  | bool | `false` | Abort via `MPI_Abort` when OOBs are detected |

When both flags are `false` and errors are found, a CRITICAL log message is
emitted and the run continues.

#### 4.1.2 Class/structs/data types

The module introduces one file-local plain struct:

```c++
struct ValidationAbortConfig {
   bool AbortOnOOB = false;
   bool AbortOnNan = true;
};
```

This struct is populated by `getValidationAbortConfig` (see §4.2.5) and
consumed by `validateOceanState`. All other types are the existing
`OceanState`, `AuxiliaryState`, `VertCoord`, and `Tracers` types from the
OMEGA ocean component.

### 4.2 Methods

#### 4.2.1 `checkOceanState` (public)

Performs all field checks and returns separate NaN and out-of-bounds counts:

```c++
std::pair<I4, I4> checkOceanState(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  const VertCoord *VCoord,
                                  I4 TimeLevel);
```

Checks all fields described in Section 2, skipping inactive cells
(`CellMask == 0`). Logs critical messages for each type of error. Returns
`{TotalNaNs, TotalOOBs}` where both values are 0 if all checks pass. Does
**not** abort. Suitable for calling from tests.

#### 4.2.2 `validateOceanState` (public)

Production entry-point that applies the configured abort policy on failure:

```c++
void validateOceanState(const OceanState *State,
                        const AuxiliaryState *AuxState,
                        const VertCoord *VCoord,
                        I4 TimeLevel);
```

Calls `checkOceanState` to obtain `{NaNs, OOBs}`, then reads
`ValidationAbortConfig` via `getValidationAbortConfig`. Aborts via
`abortWithMessage` (which logs at CRITICAL severity, prints a stack
backtrace, and calls `MPI_Abort`) if:
- `NaNs > 0` **and** `AbortOnNan` is `true`, **or**
- `OOBs > 0` **and** `AbortOnOOB` is `true`.

If neither abort condition is met but errors were detected, a CRITICAL log
message is emitted and the function returns normally.

#### 4.2.3 `checkArray2D` (file-local helper)

```c++
static std::pair<I4, I4> checkArray2D(const Array2DReal &Arr,
                                      I4 NRows, I4 NCols,
                                      Real MinVal, Real MaxVal,
                                      bool CheckMin,
                                      const Array2DReal &CellMask);
```

Performs the NaN and bounds counts for a 2-D device array over the first
`NRows` rows and `NCols` columns, restricted to active cells
(`CellMask(Row, Col) > 0`). When `CheckMin` is `false` only the upper
bound is enforced (not needed for any current field but useful for future
extension). Returns `{NaNCount, OutOfRangeCount}`.

#### 4.2.4 `checkTracerArray` (file-local helper)

```c++
static std::pair<I4, I4> checkTracerArray(const Array3DReal &Tracers3D,
                                          I4 TracerIdx,
                                          I4 NCells, I4 NVert,
                                          Real MinVal, Real MaxVal,
                                          const Array2DReal &CellMask);
```

Performs the NaN and bounds counts for a single tracer slice (identified by
`TracerIdx`) of the 3-D tracer array, restricted to active cells
(`CellMask(Cell, K) > 0`). Returns `{NaNCount, OutOfRangeCount}`.

#### 4.2.5 `getValidationAbortConfig` (file-local helper)

```c++
static ValidationAbortConfig getValidationAbortConfig();
```

Reads the two boolean options `AbortOnNan` and `AbortOnOOB` from the OMEGA
config under `Omega::State::Validation` and returns a populated
`ValidationAbortConfig`. If the config node is absent the struct's
default-initialised values (`AbortOnNan = true`, `AbortOnOOB = false`) are
returned unchanged.

#### 4.2.6 `abortWithMessage` (file-local helper)

```c++
static void abortWithMessage(const std::string &Msg);
```

Logs `Msg` at critical severity, prints a stack backtrace using `cpptrace`,
then calls `MPI_Abort` on the default OMEGA communicator with error code
`ErrorCode::Critical`.

## 5 Verification and Testing

### 5.1 Test: Valid state passes without abort

A unit test constructs a minimal OMEGA environment (MachEnv, Decomp,
HorzMesh, VertCoord, Tracers, OceanState, AuxiliaryState) using the standard
test mesh `OmegaMesh.nc`.  All state arrays are filled with physically
plausible values:

- `PseudoThickness` = 100 m (valid range [1×10⁻¹⁰, 1000])
- `NormalVelocity` = 0 m s⁻¹ (not directly validated but required for
  `KineticEnergyCell` to be zero)
- `Temperature` = 10 °C (valid range [−10, 50])
- `Salinity` = 35 g kg⁻¹ (valid range [−2, 60])

`AuxiliaryState::computeAll` is called to populate `KineticEnergyCell`
before `validateOceanState` is invoked.

The test passes if `validateOceanState` returns without calling `MPI_Abort`.

Tests requirements: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.9.

### 5.2 Negative tests: invalid values are detected

The public `checkOceanState` function is used for negative tests so that
errors can be detected without triggering `MPI_Abort`. Each sub-test:
1. Resets the state to valid values via `restoreValidState`.
2. Injects a single type of invalid value (NaN or OOB) into one field using
   a `parallelFor` kernel that overwrites all owned cell-layer entries.
3. Calls `checkOceanState`, which returns `{NaNs, OOBs}`, and verifies that
   the sum `NaNs + OOBs` is greater than zero.

The following sub-tests are implemented:

| Sub-test                      | Injected value         | Field              |
|-------------------------------|------------------------|--------------------|
| `testNaNPseudoThickness`      | NaN                    | PseudoThickness    |
| `testOOBHighPseudoThickness`  | 2000 m (> max 1000 m)  | PseudoThickness    |
| `testOOBLowPseudoThickness`   | −1 m (< min 1×10⁻¹⁰)   | PseudoThickness    |
| `testNaNKineticEnergy`        | NaN                    | KineticEnergyCell  |
| `testOOBKineticEnergy`        | 9999 J kg⁻¹ (> max 10) | KineticEnergyCell  |
| `testNaNTemperature`          | NaN                    | Temperature        |
| `testOOBTemperature`          | 9999 °C (> max 50)     | Temperature        |
| `testNaNSalinity`             | NaN                    | Salinity           |
| `testOOBSalinity`             | 9999 g kg⁻¹ (> max 60) | Salinity           |

Tests requirements: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.10.
