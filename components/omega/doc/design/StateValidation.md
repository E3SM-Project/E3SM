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
OMEGA logging infrastructure as a critical error together with a full stack
backtrace; the run is then terminated via `MPI_Abort` so that corrupted output
is not silently written to disk.

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

### 2.3 Requirement: Validate LayerThickness

`LayerThickness` must be validated for each owned cell over all vertical
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
and the number of offending elements.  After all fields are checked the
module must additionally print a stack backtrace to assist with debugging,
then abort the run via `MPI_Abort`.

### 2.9 Requirement: Graceful handling of absent tracers

If the Temperature or Salinity tracer is not present in the tracer registry
(e.g. in configurations that do not activate active tracers) the
corresponding check must be skipped silently rather than causing an error.

### 2.10 Desired: Configurable valid ranges

In the future it may be desirable to allow the user to override the default
valid ranges through the OMEGA configuration system (e.g. for idealised
process studies that intentionally use non-oceanic parameter values).

### 2.11 Desired: Configurable validation frequency

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
| `LayerThickness`     | 1×10⁻¹⁰  | 1000   |
| `KineticEnergyCell`  | 0         | 10     |
| `Temperature`        | −10       | 50     |
| `Salinity`           | −2        | 60     |

#### 4.1.2 Class/structs/data types

No new classes or data types are introduced. The module uses the existing
`OceanState`, `AuxiliaryState`, `VertCoord`, and `Tracers` types from the
OMEGA ocean component.

### 4.2 Methods

#### 4.2.1 `validateOceanState` (public)

The sole public interface of the module:

```c++
void validateOceanState(const OceanState *State,
                        const AuxiliaryState *AuxState,
                        const VertCoord *VCoord,
                        I4 TimeLevel);
```

Validates all fields described in Section 2, skipping inactive cells
(`CellMask == 0`). Logs critical errors for each failed check and aborts the
run if any check fails. `VCoord` supplies the `CellMask` array. `TimeLevel`
specifies the time-level index within the state arrays to validate (0 =
current timestep).

#### 4.2.2 `checkArray2D` (file-local helper)

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

#### 4.2.3 `checkTracerArray` (file-local helper)

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

#### 4.2.4 `abortWithMessage` (file-local helper)

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

- `LayerThickness` = 100 m (valid range [1×10⁻¹⁰, 1000])
- `NormalVelocity` = 0 m s⁻¹ (not directly validated but required for
  `KineticEnergyCell` to be zero)
- `Temperature` = 10 °C (valid range [−10, 50])
- `Salinity` = 35 g kg⁻¹ (valid range [−2, 60])

`AuxiliaryState::computeAll` is called to populate `KineticEnergyCell`
before `validateOceanState` is invoked.

The test passes if `validateOceanState` returns without calling `MPI_Abort`.

Tests requirements: 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.9.

### 5.2 Test: NaN in LayerThickness triggers abort

A future test should set a subset of `LayerThickness` entries to `NaN` and
verify that `validateOceanState` detects and logs the error and calls
`MPI_Abort`. Because `MPI_Abort` terminates the process this test would
need to be run in a separate executable or with a death-test framework.

Tests requirement: 2.1, 2.3, 2.8.

### 5.3 Test: Out-of-bounds value in Temperature triggers abort

A future test should set a subset of `Temperature` entries to a value
outside [−10, 50] (e.g. 999 °C) and verify that `validateOceanState`
detects and logs the error.

Tests requirement: 2.2, 2.5, 2.8.

### 5.4 Test: Missing tracer is skipped gracefully

A future test should invoke `validateOceanState` in a configuration where
neither Temperature nor Salinity tracers are registered and verify that the
function completes without error.

Tests requirement: 2.9.
