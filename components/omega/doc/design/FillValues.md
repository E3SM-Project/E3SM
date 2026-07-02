(omega-design-fill-values)=
# Fill Values

**Table of Contents**
1. [Overview](#1-overview)
2. [Requirements](#2-requirements)
3. [Algorithmic Formulation](#3-algorithmic-formulation)
4. [Design](#4-design)
5. [Verification and Testing](#5-verification-and-testing)

## 1 Overview

Omega supports output of ocean fields to NetCDF files via I/O streams. Each field
carries a declared fill value stored as metadata under the `_FillValue` attribute,
which marks entries that are undefined or outside the active domain. Currently, fill
value constants are declared locally in each module with inconsistent values
(typically `-9.99e30` for real-valued fields, `-999` for integers). These values do
not match the NetCDF-C standard fill values that analysis tools and visualization
software commonly expect. Additionally, the Kokkos arrays backing these fields are
never explicitly initialized with fill values: inactive ocean layers (model layers
below `MaxLayerCell`) hold uninitialized memory rather than a well-defined sentinel,
and the distinction between "no data" and "zero" is lost.

This design standardizes fill values across all Omega fields to match the NetCDF-C
standard, ensures that field arrays are automatically initialized with the correct
fill value when data is attached to a `Field`, and verifies through CTests that the
resulting behavior is correct for inactive layers and for edge/vertex fields at
domain boundaries.

## 2 Requirements

### 2.1 Requirement: Centralized fill value constants

Fill value constants must be defined in a single location in Omega that is accessible
to all modules. Scatter of locally-declared constants has produced inconsistencies
(e.g., `-9.99e30` vs `-9.99E+30`) and makes future maintenance error-prone.

### 2.2 Requirement: Standard fill values matching NetCDF-C

Fill value constants must exactly match the NetCDF-C standard values (`NC_FILL_*`)
as exposed by the SCORPIO library (`PIO_FILL_*` in `pio.h`). Using standard values
ensures compatibility with NetCDF-aware tools (ncview, Xarray, NCO, etc.) that
recognize and handle fill values automatically. The required values are:

| Type | Omega name       | Value                         |
|------|-----------------|-------------------------------|
| I4   | `FillValueI4`   | -2147483647 (`PIO_FILL_INT`)  |
| I8   | `FillValueI8`   | `PIO_FILL_INT64`              |
| R4   | `FillValueR4`   | 9.9692099683868690e+36f (`PIO_FILL_FLOAT`) |
| R8   | `FillValueR8`   | 9.9692099683868690e+36 (`PIO_FILL_DOUBLE`) |
| Real | `FillValueReal` | Alias for `FillValueR4` or `FillValueR8` depending on build |

### 2.3 Requirement: Automatic array initialization at attach time

When `Field::attachData()` is called, the attached Kokkos array must be
automatically initialized with the standard fill value for the array's element
type before any compute routine runs. The fill value is deduced from
`typename T::value_type` at `attachData()` time using the constants in
`FillValues.h` (`FillValueI4`, `FillValueI8`, `FillValueR4`, or `FillValueR8`),
so the caller does not need to specify it. This guarantees that inactive ocean
layers (below `MaxLayerCell`) and any other undefined entries contain a
well-defined sentinel in output, without requiring each module to manually
perform the fill.

### 2.4 Requirement: No fill value argument in Field::create()

The `FillValue` argument is removed from `Field::create()`. Because the fill
value is always one of the standard constants determined by the array element
type, specifying it at field creation time is redundant. The fill value and the
`_FillValue` metadata entry are set automatically when `attachData()` is called.

### 2.5 Requirement: Inactive layers contain fill values in output

After compute routines execute, model layers with index `k > MaxLayerCell(ICell)`
must contain the declared fill value in output. With automatic initialization at
attach time (Requirement 2.3) and compute routines that only write to active layers,
this is satisfied without changes to individual compute functions.

### 2.6 Desired: Valid values at domain boundary edges and vertices

Edge and vertex fields at valid boundaries adjacent to land or to bathymetry must
contain computed valid values, not fill values. Switching from zero-initialization
to fill-value initialization must not silently break any boundary computation that
previously assumed zero defaults. This mirrors a problem encountered in MPAS-Ocean
when initializing `NormalVelocity` with a fill value rather than zero: code that
implicitly assumed zero at land boundaries produced incorrect results. During
implementation, all compute routines for edge/vertex fields must be audited to
confirm that every valid boundary entry receives an explicit computed value. As an
example of what such an audit can reveal: `TracerHorzAdvOnCell` accumulated its
internal working array `HighOrderFlxHorz` over a cell's full active layer range
without clamping to each contributing edge's active range, reading fill values from
edges shallower than the cell — see Section 4.5 for the resolution.

### 2.7 Requirement: CTests

CTests must verify Requirements 2.1–2.3, 2.5, and Desired Requirement 2.6. The
three-zone boundary behavior (2.6) is exercised for `NormalVelocity` via the edge
mask; the two-zone cell and vertex masks are verified with synthetic cell and
vertex fields.

## 3 Algorithmic Formulation

No new numerical algorithm is introduced. The key behavioral change is in
`Field::attachData<T>()`, which gains a fill step equivalent to:

```c++
Kokkos::deep_copy(InDataArray, fill_scalar);
```

where `fill_scalar` is the typed fill value extracted from the field's `FieldMeta`
map (stored as `std::any` under key `"_FillValue"`), and the element type
`typename T::value_type` of the Kokkos view determines the `std::any_cast` type.
`Kokkos::deep_copy` with a scalar source broadcasts that value to every element of
the view and dispatches correctly for both host and device memory spaces.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Constants: FillValue definitions in FillValues.h

The fill value constants are defined in a new lightweight header
`components/omega/src/base/FillValues.h`. The values match the NetCDF-C `NC_FILL_*`
constants (which are identical to the SCORPIO `PIO_FILL_*` constants), but are
written as numeric literals so that the header can be included from any context
without pulling in `pio.h` or `mpi.h`, which carry strict include-order requirements.

A type-indexed variable template `FillValue<T>` is the primary definition. Its
primary template is intentionally left undefined so that instantiation with an
unsupported type produces a link error. Explicit specializations are provided for
each supported scalar type:

```c++
template <typename T>
constexpr T FillValue;  // primary template intentionally undefined

template <> inline constexpr I4 FillValue<I4> = -2147483647;             // NC_FILL_INT
template <> inline constexpr I8 FillValue<I8> = -9223372036854775806LL;  // NC_FILL_INT64
template <> inline constexpr R4 FillValue<R4> = 9.9692099683868690e+36f; // NC_FILL_FLOAT
template <> inline constexpr R8 FillValue<R8> = 9.9692099683868690e+36;  // NC_FILL_DOUBLE
```

Named aliases are provided for readability in comparison and test code:

```c++
constexpr I4   FillValueI4   = FillValue<I4>;
constexpr I8   FillValueI8   = FillValue<I8>;
constexpr R4   FillValueR4   = FillValue<R4>;
constexpr R8   FillValueR8   = FillValue<R8>;
#if defined(SINGLE_PRECISION)
constexpr Real FillValueReal = FillValue<R4>;
#else
constexpr Real FillValueReal = FillValue<R8>;
#endif
```

`FillValues.h` is included from `Field.h`. Because every module that creates or
uses fields already includes `Field.h` (directly or transitively), no additional
include directives are needed in module source files.

#### 4.1.2 Class/struct changes

No new classes or structs are introduced. The changes are confined to:
- new header `FillValues.h` with constants (Section 4.1.1)
- `FillValues.h` included from `Field.h`
- the `attachData` template methods in `Field.h` / `Field.cpp` (Section 4.2)
- removal of local fill value declarations in module files (Section 4.3)

### 4.2 Methods

#### 4.2.1 Field::attachData() auto-fill

The `attachData<T>()` instance method and the `attachFieldData<T>()` static helper
in `components/omega/src/infra/Field.h` gain a private fill step controlled by an
optional `FillOnAttach` parameter (default `true`):

```c++
template <typename T>
void attachData(const T &InDataArray, bool FillOnAttach = true) {
    OMEGA_ASSERT(isKokkosArray<T>,
                 "Field::attachData requires Kokkos array as input");
    DataArray = std::make_shared<T>(InDataArray);
    DataType  = checkArrayType<T>();
    MemLoc    = findArrayMemLoc<T>();
    if (FillOnAttach)
        fillWithValue<T>(InDataArray); // auto-fill with type-deduced FillValue
}

template <typename T>
static void attachFieldData(const std::string &FieldName,
                            const T &InDataArray, bool FillOnAttach = true);
```

The default `FillOnAttach = true` preserves the fill behavior for the normal
initialization workflow (fresh arrays attached before any data is written). Pass
`FillOnAttach = false` when re-attaching an existing data-filled array where the
fill must not overwrite computed values — see Section 4.2.3 for the canonical use case.

The private helper `fillWithValue<T>()` performs the following steps:

1. Derives the element type as `typename T::non_const_value_type` (using the
   non-const form so that const view types resolve to the same underlying type).
2. Looks up the fill value via `FillValue<ValType>`, a constexpr variable template
   specialization from `FillValues.h`. Supporting a new type requires only adding
   one specialization there; no changes to `fillWithValue` are needed.
3. Stores the typed fill value in `FieldMeta["_FillValue"]` (the NetCDF/CF
   standard attribute name), making it available to IO and metadata queries.
4. Calls `Kokkos::deep_copy(InDataArray, scalar)` to broadcast the scalar to every
   element of the view. This works for host views (`HostArray*`), device views
   (`Array*`), and all dimensionalities (1D–5D) without additional dispatch.
5. Only the view being attached is filled; no assumption is made about whether the
   other memory space (host or device) has been attached.

Ordering invariants preserved by this design:
- `Field::create()` is always called before `attachData()`, so the field exists
  in `AllFields` before fill time.
- The `"_FillValue"` metadata entry is set at `attachData()` time, which precedes
  any IO write operation in the initialization workflow.
- Restart reads happen after `attachData()` and overwrite fill entries with data
  from the restart file.
- Compute routines run after initialization and overwrite active entries, leaving
  inactive entries (below `MaxLayerCell`) at the fill value.

The fill value also serves as the "not read" sentinel for optional reads. A field
marked with `Field::setOptionalRead(true)` that is absent from an input file is
skipped by the IOStream read and keeps the fill value from `attachData()`; the
owning module detects the fill value afterward and substitutes a default. For
example, `VertCoord::SurfacePressure` defaults to zero when it is not present in
the initial-condition or restart file.

#### 4.2.2 No per-module manual initialization required

Because `attachData()` handles initialization automatically, no module needs to
call an explicit fill routine. Modules may still call `Kokkos::deep_copy` or
similar to set specific values, but this is independent of fill-value management.

#### 4.2.3 Ordering invariant and anti-patterns

**Invariant**: `attachData()` with `FillOnAttach = true` (the default) must always
be called *before* any real values are written to the array. Any data written before
such a call is silently overwritten by the auto-fill step and is dead code.

Three anti-patterns discovered and corrected during implementation:

**Class A — Wasted pre-fill**: A `deepCopy` (or equivalent write) is applied to an
array before `defineFields()` calls `attachData()` on it. Because `attachData()`
overwrites the entire array, the earlier write has no effect.

```c++
// WRONG: initial value is overwritten by attachData()
deepCopy(SpecVol, 1.0_Real / RhoSw);
defineFields(); // calls attachData(), fills SpecVol with FillValueReal
```

**Class B — Compute-time re-initialization**: A `deepCopy(arr, 0)` (or any scalar)
at the top of a recurring compute method resets the entire array on every timestep,
overwriting the correct fill values in inactive layers with zero.

```c++
// WRONG: destroys fill values in inactive layers every timestep
void Foo::compute(...) {
    deepCopy(LocArr, 0);
    parallelFor({NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
        // writes only active layers KMin..KMax
    });
}
```

The correct pattern for Classes A and B is to do nothing: `attachData()` sets the
fill value once at initialization, and compute kernels that iterate only over active
layers naturally leave inactive layers at the fill value.

**Class C — Time-level pointer update**: `attachData()` is called each timestep to
re-point the IO field to the newly-current time-level array (e.g., in
`OceanState::updateTimeLevels()` and `Tracers::updateTimeLevels()`). The array
already contains the just-computed state. Calling `attachData()` with the default
`FillOnAttach = true` here would silently destroy the computed values before IO
writes them. This is the most dangerous anti-pattern because the model runs to
completion without aborting, but all output fields contain fill values.

```c++
// WRONG: destroys computed state before IO can write it
void OceanState::updateTimeLevels() {
    CurTimeIndex = (CurTimeIndex + 1) % NTimeLevels;
    Field::attachFieldData<Array2DReal>(NormalVelocityFldName,
                                        NormalVelocity[CurTimeIndex]); // fills!
}

// CORRECT: pointer update only, computed values preserved
void OceanState::updateTimeLevels() {
    CurTimeIndex = (CurTimeIndex + 1) % NTimeLevels;
    Field::attachFieldData<Array2DReal>(NormalVelocityFldName,
                                        NormalVelocity[CurTimeIndex], false);
}
```

A related instance occurs in the IO write path when the time-coordinate array is
populated before `attachFieldData("time", OutTime)` is called. The correct pattern
is to attach first (which fills the single-element array with the fill value), then
assign the elapsed time to `OutTime(0)`. Because Kokkos host views share backing
memory, the assignment is visible through the stored pointer before the field is
written to the file.

### 4.3 Module-level changes: removing FillValue from Field::create() calls

The `FillValue` parameter is removed from `Field::create()`. All ~50 call sites
across Omega drop the fill value argument. The pattern is:

```c++
// BEFORE — explicit fill value argument
Field::create(FieldName, Description, Units, StdName,
              ValidMin, ValidMax, FillValueReal, NDims, DimNames);

// AFTER — fill value is deduced automatically from the array type at attachData()
Field::create(FieldName, Description, Units, StdName,
              ValidMin, ValidMax, NDims, DimNames);
```

The `Tracers::define()` helper in `src/ocn/Tracers.h` / `Tracers.cpp` and
`TracerDefs.inc` also have their `FillValue` parameter and arguments removed,
since `Tracers::define()` calls `Field::create()` internally.

### 4.4 Boundary edge and vertex handling

During implementation, all compute routines that populate fields defined on mesh
edges or vertices must be audited to confirm that every valid boundary entry —
including edges and vertices adjacent to land cells and edges at the base of the
water column — receives an explicit computed value rather than relying on
zero-initialization. Where such reliance is found, an explicit assignment must be
added. The CTest in Section 5.5 detects regressions after the change.

### 4.5 Edge flux field masking

Edge fields that represent fluxes normal to the edge face — such as
`NormalVelocityTend`, `VelocityDel2Aux.Del2Edge`, and `NormalVelocity` itself —
must carry zeros (not fill values) in the boundary layer sub-ranges
`[MinLayerEdgeTop, MinLayerEdgeBot)` and `(MaxLayerEdgeTop, MaxLayerEdgeBot]`,
where one neighboring cell is active but the other is not. Two helper methods
are added to `VertCoord` (`components/omega/src/ocn/VertCoord.h/.cpp`):

**`zeroEdgeField(Array2DReal &Arr, I4 NEdgesAll)`** — zeros all layers in
`[MinLayerEdgeTop, MaxLayerEdgeBot]`. Used before recomputing a flux-type edge
field each time step so that boundary edges show 0, not `FillValueReal`. Layers
outside the valid range retain their fill value from `attachData()`. Called for:
- `NormalVelocityTend` at the start of `Tendencies::computeVelocityTendenciesOnly()`
- `VelocityDel2Aux.Del2Edge` before the `edgeAuxState2` kernel in
  `AuxiliaryState::computeMomAux()`, because the Del4 hyperdiffusion kernel
  accumulates `Del2Edge` over all edges sharing a cell or vertex; fill values in
  boundary layers of a neighbouring edge would corrupt `Del2DivCell` and
  `Del2RelVortVertex`.

**`applyEdgeLayerMask(Array2DReal &Arr, I4 NEdgesAll)`** — enforces the full
three-zone mask on an edge field after IC or restart read: layers outside
`[MinLayerEdgeTop, MaxLayerEdgeBot]` are set to `FillValueReal`; layers inside
`[MinLayerEdgeTop, MaxLayerEdgeBot]` but outside the active range
`[MinLayerEdgeBot, MaxLayerEdgeTop]` are set to 0; active layers are left
unchanged.

**`applyCellLayerMask(Array2DReal &Arr, I4 NCellsAll)`** — enforces the
two-zone mask on a cell field after IC or restart read: layers outside the
active range `[MinLayerCell, MaxLayerCell]` are set to `FillValueReal`; active
layers are left unchanged. Cell fields have no boundary (zero) zone because a
cell column is either active or inactive at a given layer — there is no
neighbor-dependent partial-activity range as there is for edges and vertices.

**`applyVertexLayerMask(Array2DReal &Arr, I4 NVerticesAll)`** — enforces the
two-zone mask on a vertex field: layers outside `[MinLayerVertexTop,
MaxLayerVertexBot]` are set to `FillValueReal`; active layers are left
unchanged. Unlike an edge, a vertex does **not** have a zeroed boundary zone: a
boundary vertex with one or more active surrounding cells holds valid, generally
non-zero data (for example, relative vorticity at such a vertex is computed from
its active surrounding edges — see `VorticityAuxVars::computeVarsOnVertex`, which
loops over the full `[MinLayerVertexTop, MaxLayerVertexBot]` range). Zeroing that
boundary band would discard real signal, so the whole valid range is kept. There
is currently no vertex-based state field read from IC/restart, so this method has
no production caller yet; it is provided for completeness and is exercised by the
CTest (Section 5.6).

All three `apply*LayerMask` methods use the inclusive `[Min, Max]` layer
convention consistent with `computeGeomZHeight`/`computePressure`. They are
driven through `OceanState::applyLayerMasks(TimeLevel)`
(`src/ocn/OceanState.cpp`), which is called from `ocnInit` (`OceanInit.cpp`)
after `exchangeHalo` and before `copyToHost`. `applyLayerMasks` applies the
edge mask to `NormalVelocity` and the cell mask to `PseudoThickness`,
`Temperature`, and `Salinity` — the state and tracer fields read from
IC/restart.

Fields that are *not* fluxes — `MeanPseudoThickEdge`, `FluxPseudoThickEdge` —
are genuinely undefined at boundary edges (interpolating thickness when one
neighboring cell is land has no physical meaning) and therefore correctly retain
`FillValueReal` in those layers; no zeroing is applied to them.
`VorticityAux` on edges is already computed over the full valid range
`[MinLayerEdgeTop, MaxLayerEdgeBot]`, so it requires no additional masking.

**`HighOrderFlxHorz` in `TracerHorzAdvOnCell`** — the 3-D working array
`(NTracers, NEdgesSize, NVertLayers)` that accumulates horizontal tracer flux
contributions per edge requires separate treatment. The edge-pass overload of
`TracerHorzAdvOnCell::operator()` writes `HighOrderFlxHorz(L, IEdge, K)` only
for layers `K ∈ [0, MaxLayerEdgeTop(IEdge)]`; layers beyond that index retain
their fill values from `attachData()`. The cell-pass overload then reads
`HighOrderFlxHorz` over the cell's own active range `[MinLayerCell, MaxLayerCell]`,
which for cells adjacent to shallower edges extends past `MaxLayerEdgeTop` of
those edges — reading fill values and corrupting active-layer tracer tendencies.
The fix is to add `MinLayerEdgeBot` and `MaxLayerEdgeTop` as members of
`TracerHorzAdvOnCell` and clamp the cell-pass inner loop to
`[max(KStart, MinLayerEdgeBot(IEdge)), min(KEnd, MaxLayerEdgeTop(IEdge)+1))`.
This follows the identical pattern already used by `PseudoThicknessFluxDivOnCell`,
`TracerDiffOnCell`, `TracerAuxVars::computeVarsOnCells`, and every other
cell-level flux-accumulation kernel, and avoids the need to zero
`HighOrderFlxHorz` before the edge pass each timestep.

## 5 Verification and Testing

A test `test/ocn/FillValueTest.cpp` is added and registered as `FILL_VALUE_TEST`
with 8 MPI tasks (building `testFillValue.exe`) using `add_omega_test()` in the
test `CMakeLists.txt`. The test performs a full ocean initialization (matching
the `StateTest` pattern: `MachEnv`, logging, Config, `TimeStepper`, IO, `Field`,
`IOStream`, `Decomp`, `Halo`, `HorzMesh`, `VertCoord`, `Tracers`,
`AuxiliaryState`, `PressureGrad`, `Tendencies`, `OceanState`), reads the initial
state, exchanges halos, and applies `applyEdgeLayerMask` to `NormalVelocity`
before running its assertions. Each test counts mismatches and, on any failure,
accumulates an `Error(ErrorCode::Fail, ...)` into a top-level `Error` whose
final state sets the process return code.

### 5.1 Test: fill constant values

Verify that each Omega fill value constant exactly equals its NetCDF-C
counterpart (`FillValueI4 == NC_FILL_INT`, `FillValueI8 == NC_FILL_INT64`,
`FillValueR4 == NC_FILL_FLOAT`, `FillValueR8 == NC_FILL_DOUBLE`), comparing
directly against the `<netcdf.h>` `NC_FILL_*` macros.

Tests requirements: 2.1, 2.2, 2.7

### 5.2 Test: attachData auto-fill

Create a `Field`, allocate a Kokkos host array whose elements are all set to a
distinct sentinel (`0`), then call `attachData()`. Verify that after the call
every element of the array equals `FillValueR8`.

Tests requirements: 2.3, 2.7

### 5.3 Test: inactive layers contain fill values

After `VertCoord` initialization, verify that the cell field `GeomZMid` carries
`FillValueReal` in every inactive layer (`k >= MaxLayerCell(ICell)`) for all
owned cells. This confirms that the auto-fill applied at `attachData()` time
survives initialization for inactive layers.

Tests requirements: 2.5, 2.7

### 5.4 Test: NormalVelocity three-zone fill/zero/real pattern

After full state initialization, IC read, halo exchange, and
`VertCoord::applyEdgeLayerMask()`, verify that `NormalVelocity` satisfies the
three-zone invariant for every owned edge:

- **Fully inactive layers** (outside `[MinLayerEdgeTop, MaxLayerEdgeBot]`):
  equal `FillValueReal`.
- **Boundary layers** (inside `[MinLayerEdgeTop, MaxLayerEdgeBot]` but outside
  `[MinLayerEdgeBot, MaxLayerEdgeTop]`): equal `0`.
- **Active layers** (`[MinLayerEdgeBot, MaxLayerEdgeTop]`): not equal to
  `FillValueReal` (for a quiescent IC the value is `0`, which is valid).

This guards against two regressions: (a) the original zero-initialization
assumption (boundary layers silently filled with zero before `attachData()`
was changed) and (b) IC files that store zeros in inactive layers, which
`applyEdgeLayerMask` must replace with `FillValueReal`.

Tests requirements: 2.5, 2.6, 2.7

### 5.5 Test: applyCellLayerMask two-zone pattern

Allocate a synthetic cell field, fill it with a sentinel value (≠ 0 and
≠ `FillValueReal`), apply `VertCoord::applyCellLayerMask()`, and verify for every
owned cell that:

- **Inactive layers** (outside `[MinLayerCell, MaxLayerCell]`): equal
  `FillValueReal`.
- **Active layers** (`[MinLayerCell, MaxLayerCell]`): retain the sentinel value.

Using a sentinel rather than IC data gives an exact check that the mask
overwrites only the inactive layers.

Tests requirements: 2.5, 2.7

### 5.6 Test: applyVertexLayerMask two-zone pattern

Allocate a synthetic vertex field, fill it with the same sentinel value, apply
`VertCoord::applyVertexLayerMask()`, and verify for every owned vertex that:

- **Inactive layers** (outside `[MinLayerVertexTop, MaxLayerVertexBot]`): equal
  `FillValueReal`.
- **Active layers** (`[MinLayerVertexTop, MaxLayerVertexBot]`): retain the
  sentinel value.

Vertices have no zeroed boundary zone (a boundary vertex holds valid data), so
the whole valid range is checked for the sentinel. Because no vertex-based state
field is read from IC/restart, this synthetic test is the primary coverage for
`applyVertexLayerMask`.

Tests requirements: 2.5, 2.7
