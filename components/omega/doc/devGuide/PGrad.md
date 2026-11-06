(omega-dev-pgrad)=

# Pressure Gradient (PGrad)

Omega includes a `PressureGrad` class that computes horizontal pressure gradient
tendencies for the non-Boussinesq momentum equation. The implementation supports a
layer-integrated finite-volume scheme as the default and a centered difference
scheme as an option. The class follows the same factory pattern used by other
Omega modules.

## PressureGradType enum

An enumeration of the available pressure gradient schemes is defined in `PGrad.h`:

```c++
enum class PressureGradType {
   Centered,    ///< existing 2nd-order Montgomery scheme
   FiniteVolume ///< layer-integrated finite-volume scheme
   // , <FutureVariant>  ///< e.g. a 6th-order option, added when implemented
};
```

This is used to select which pressure gradient method is applied at runtime.

## Initialization

An instance of `PressureGrad` requires both a [`HorzMesh`](#omega-dev-horz-mesh) and
a [`VertCoord`](#omega-dev-vert-coord), so these classes and all of their dependencies
must be initialized before `PressureGrad` can be initialized. The static method:

```c++
OMEGA::PressureGrad::init();
```

initializes the default `PressureGrad` instance using the default `HorzMesh` and
`VertCoord` instances and the global Omega configuration. A pointer to the default
instance can be retrieved at any time using:

```c++
OMEGA::PressureGrad* DefPGrad = OMEGA::PressureGrad::getDefault();
```

## Creating additional instances

Additional named instances can be created using:

```c++
OMEGA::PressureGrad* MyPGrad =
    OMEGA::PressureGrad::create("MyPGrad", Mesh, VCoord, Options);
```

where `Mesh` is a pointer to a `HorzMesh`, `VCoord` is a pointer to a `VertCoord`,
and `Options` is a pointer to the `Config`. A named instance can be retrieved later
using:

```c++
OMEGA::PressureGrad* MyPGrad = OMEGA::PressureGrad::get("MyPGrad");
```

## Constructor behavior

The constructor reads the `PressureGrad` section from the configuration, stores
references to mesh and vertical coordinate data needed for computation, and enables
the appropriate functor based on the configured `PressureGradType`. It also allocates
placeholder arrays for tidal potential and self-attraction/loading, which are intended
to be populated by a future tidal forcing module. Currently these arrays are initialized
to zero.

## Computing the pressure gradient

To compute pressure gradient tendencies and accumulate them into a tendency array:

```c++
PGrad->computePressureGrad(Tend, PressureMid, PressureInterface, SpecVol,
                           GeomZInterface, PseudoThick, ConservTemp,
                           AbsSalinity, EqState);
```

where:
- `Tend` is a 2D array `(NEdgesAll × NVertLayers)` that the pressure gradient
  tendency is accumulated into
- `PressureMid`, `PressureInterface` and `GeomZInterface` come from `VertCoord`
- `SpecVol` is the specific volume field from `Eos`
- `PseudoThick` is the pseudo-thickness at the desired time level
- `ConservTemp` and `AbsSalinity` are the layer-mean tracer fields, used by the
  finite-volume scheme's vertical reconstruction and ignored by the centered one
- `EqState` supplies the specific volume derivatives, likewise used only by the
  finite-volume scheme

The method uses hierarchical Kokkos parallelism: an outer `parallelForOuter` loop
iterates over edges, and an inner `parallelForInner` loop iterates over vertical
chunks. The appropriate functor is dispatched based on `PressureGradChoice`. The
finite-volume branch additionally runs two kernels ahead of that dispatch: one
forming the per-cell reconstruction slopes, and the per-edge column scan.

Note that `AuxiliaryState::computeMomVertAux` calls `Eos::computeSpecVolAndDerivs`
in place of `Eos::computeSpecVol` when the finite-volume scheme is selected. That
call fills `SpecVol` as well, so it replaces rather than accompanies the simpler
one, and the branch keeps a centered run from paying for the derivative
arithmetic.

## Functors

### PressureGradCentered

This functor implements a centered difference approximation of the pressure gradient
tendency. For each edge, it first computes the layer-invariant tidal and
self-attraction/loading contribution:

```
GradGeoPot = grad(TidalPotential) + grad(SelfAttractionLoading)
```

Then, for each vertical layer `K`, it computes three terms:

1. **Montgomery potential gradient**: The average of the horizontal gradients of the
   Montgomery potential ($\alpha p + g z$) at the top (interface `K`) and bottom
   (interface `K+1`) of the layer. This compactly represents the combined effect
   of the pressure gradient and the geopotential contribution from tilted coordinate
   surfaces.

2. **Specific volume correction**: A correction term equal to the edge-averaged
   pressure at mid-layer multiplied by the horizontal gradient of specific volume.
   This accounts for horizontal density variations that are not captured by the
   Montgomery potential form.

3. **Tidal and geopotential forcing** (`GradGeoPot`): The external geopotential
   contribution from tidal forcing and self-attraction/loading, applied uniformly
   across all layers at an edge.

The tendency update for each layer is:

```
Tend(IEdge, K) += EdgeMask(IEdge, K) * (-GradMontPot + PGradAlpha - GradGeoPot)
```

where `EdgeMask` is applied to enforce land boundary conditions. The functor operator
signature is:

```c++
KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                const Array2DReal &PressureMid,
                                const Array2DReal &PressureInterface,
                                const Array2DReal &GeomZInterface,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading,
                                const Array2DReal &SpecVol) const;
```

### PressureGradFiniteVolume

This functor assembles the finite-volume pressure gradient tendency. The whole
horizontal pressure gradient is the geopotential compared at **fixed pressure**,
and almost all of the work of forming that comparison happens before the functor
runs, in the column scan described below. What remains per layer is:

```
DeltaPress = edge average of the two columns' layer pressure thicknesses
LayerMean  = DeltaZFixedP(IEdge, K+1) + DeltaZMoment(IEdge, K) / DeltaPress
Tend(IEdge, K) += EdgeMask(IEdge, K) * (-Gravity / DcEdge * LayerMean - GradGeoPot)
```

`LayerMean` is the layer average of the fixed-pressure height difference,
recovered from its value at the layer's bottom interface and the first moment of
the integrand over the layer. Its signature therefore takes the column scan's two
output arrays rather than the state the scan consumed:

```c++
KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                const Array2DReal &PressureInterface,
                                const Array2DReal &DeltaZFixedP,
                                const Array2DReal &DeltaZMoment,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading) const;
```

## The per-edge column scan

`PressureGrad::computeColumnScan` fills `DeltaZFixedP`, the fixed-pressure height
difference at each edge-layer interface, and `DeltaZMoment`, the first moment of
the integrand over each layer.

**This cannot live in the functor.** It is a prefix sum down each column with
edge-dependent coefficients, so it is not expressible as an independent
per-vertical-chunk operation. It runs as a `parallelForOuter` over edges: a
`parallelForInner` over the edge's active layers forms both integrals, and after a
team barrier a `parallelScanInner` accumulates them up the column, in the same
shape as `VertCoord::computeGeomZHeight`. It is the one structural addition the
finite-volume scheme makes beyond the per-edge, per-chunk pattern the centered
scheme uses.

The scan is anchored at the **sea floor**. `VertCoord` builds geometric height
upward from a prescribed bathymetry, so at the bottom interface the cross-edge
height difference is exact input and vanishes identically for a flat floor, where
at the surface it would be the small residual of two column-length accumulations.
There is a second reason: `VertCoord` accumulates a midpoint rule over each
column's own layers, and on a curved profile two columns with different layer
partitions give sums that differ at second order in layer thickness. Anchored at
the surface that discrepancy would enter the height difference directly; anchored
at the sea floor it never enters, because the scheme integrates its own
reconstruction rather than accumulating `VertCoord`'s height.

The anchor is computed, not assumed: each column's height at the deepest shared
interface is shifted to the common pressure by integrating its own reconstruction
over half the cross-edge pressure difference, and both short integrals vanish
where the two columns' interface pressures agree.

Those short integrals use one first-order expansion of the equation of state,
shared by the two bottom layers, so they are only accurate across a pressure
difference of about a layer. That covers a flat floor, and most of a sloping
one cut by partial bottom cells beneath level interior layers. The floors are
just over a layer apart wherever a thin bottom cell lies entirely above the
deeper column's layer of the same index, and terrain-following layers over a
steep slope put one floor several layers deeper than the other. An expansion
stretched across that whole offset leaves an error that grows with the cube of
the offset and is the same at every depth. So where the floors are more than a
layer apart, the shallower floor is compared instead with the deeper column's own
interface just below its pressure, taking the deeper column's height there from
`VertCoord` and sharing the expansion between the two layers that meet there.
The anchor is then carried down to the edge's deepest interface with the
matched-pressure integrand, which vanishes on the exact set like every other
increment of the scan. Where the floors are within one layer of each other the
deeper column's interface is the shared one, and the anchor reduces exactly to
the comparison above.

The scheme has a known instability on terrain-following (sigma) layers over a
steep seamount that the centered scheme does not share. It is described, with
the tests that rule out the anchor, extrapolation below the floor and the tracer
limiter as causes, in the
[user guide](omega-user-pgrad-sigma-instability).

## Each column's state is looked up by pressure, not by layer index

At each quadrature point of edge layer `K`, the integrand needs the
reconstruction of whichever of *that column's* layers contains that pressure,
which under coordinate tilt is generally **not** layer `K`. At a tilt of 50 m/km
with 64 m layers the two columns' layer `K` are offset by nearly three layer
thicknesses and do not overlap in pressure at all.

`findLayerForPress` in `PGradRecon.h` performs that lookup. Within the column
scan the answer advances monotonically with `K`, so passing the previous answer
as a hint makes it a pair of incremented cursors rather than a search; the result
does not depend on the hint. A pressure outside the column clamps to the
outermost valid layer and extrapolates its reconstruction, which is the rule
where the edge control volume extends past a column's own floor.

**This is the single most important thing not to get wrong, and no answer-level
test can catch it.** Replacing the lookup with a layer-index lookup passes every
exactness and accuracy check in the test suite: on a profile linear in pressure
every layer's mean-preserving reconstruction is that same line, so looking up the
wrong layer costs nothing. It is pinned instead by direct property tests in
`PGradTest.cpp`, which assert that the returned layer brackets the pressure, that
it differs from the edge layer index under tilt, that the answer does not depend
on the starting hint, and that a pressure outside the column clamps. Those tests
are not optional.

## Supporting headers

| Header | Contents |
| ------ | -------- |
| `PGradRecon.h` | the mean-preserving linear reconstruction of temperature and salinity in pressure, and the per-column pressure lookup |
| `PGradFiniteVolume.h` | the edge-shared equation-of-state expansion, the matched-pressure integrand, and the Gauss-Legendre rule |

The equation-of-state expansion is **shared across each edge**: the coefficients
and the expansion state are averaged from the edge's two cells and used for
*both* columns. That one set multiplies both columns is what the exactness rests
on -- it is what makes the constant and pressure-derivative terms cancel in the
matched-pressure difference. Which set is used is an ordinary accuracy question
and cannot break exactness, because the coefficients end up multiplying a
quantity that is identically zero on profiles the reconstruction resolves.

Specific volume is never integrated directly and the equation of state is never
evaluated inside any integral. The four expansion coefficients come from the
`SpecVol`, `SpecVolDCt`, `SpecVolDSa` and `SpecVolDP` fields `Eos` already
computes, one evaluation per cell per layer, and the pressure gradient adds none
of its own.

## Configuration

The pressure gradient type is selected in the input YAML file:

```yaml
PressureGrad:
   PressureGradType: 'FiniteVolume'
```

Valid options for `PressureGradType` are:
- `'centered'` or `'Centered'`: centered difference approximation
- `'finiteVolume'` or `'FiniteVolume'`: the layer-integrated finite-volume method
  (default)

The `FiniteVolume` scheme reads three further keys from the same group, all
optional: `HorzOrder` (default 2), `VerticalReconstruction` (default `'linear'`)
and `QuadraturePoints` (default 2). Only `HorzOrder: 2` and
`VerticalReconstruction: 'linear'` are implemented; other values are rejected
with an error rather than falling back, so that a configuration requesting them
cannot quietly run as the implemented scheme. `QuadraturePoints` is an accuracy
setting only: the integrand is zero at every point on a resolved profile, so no
quadrature rule can affect the exactness.

An unrecognized value is a fatal error rather than a silent fallback to the
centered scheme, so that a typo -- or a configuration naming a scheme that no
longer exists -- cannot produce a run that looks like a passing centered run.

## Data members

The `PressureGrad` class stores the following key data:

| Member | Type | Description |
| ------ | ---- | ----------- |
| `NEdgesAll` | `I4` | Total number of edges including halo |
| `NChunks` | `I4` | Number of vertical chunks for vectorization |
| `NVertLayers` | `I4` | Number of vertical layers |
| `NVertLayersP1` | `I4` | Number of vertical layers plus one |
| `MinLayerEdgeBot` | `Array1DI4` | Minimum active layer index for each edge |
| `MaxLayerEdgeTop` | `Array1DI4` | Maximum active layer index for each edge |
| `TidalPotential` | `Array1DReal` | Tidal potential (placeholder, currently zero) |
| `SelfAttractionLoading` | `Array1DReal` | Self-attraction and loading term (placeholder, currently zero) |
| `CenteredPGrad` | `PressureGradCentered` | Centered pressure gradient functor |
| `FiniteVolumePGrad` | `PressureGradFiniteVolume` | Finite-volume pressure gradient functor |
| `PressureGradChoice` | `PressureGradType` | Selected pressure gradient method |
| `ReconSlopeCt` | `Array2DReal` | Reconstruction slope of temperature in pressure, per cell and layer |
| `ReconSlopeSa` | `Array2DReal` | Reconstruction slope of salinity in pressure |
| `DeltaZIncr` | `Array2DReal` | Per-layer integral of the matched-pressure integrand |
| `DeltaZMoment` | `Array2DReal` | Its first moment about the layer's top interface |
| `DeltaZFixedP` | `Array2DReal` | Fixed-pressure height difference at edge-layer interfaces |

The last five are allocated only when the `FiniteVolume` scheme is selected, so a
centered run pays no memory for them. The reconstruction slopes are computed once
per cell and layer and reused across each cell's edges rather than recomputed per
edge.

## Removal

To remove all `PressureGrad` instances:

```c++
OMEGA::PressureGrad::clear();
```

To remove a specific named instance:

```c++
OMEGA::PressureGrad::erase("MyPGrad");
```
