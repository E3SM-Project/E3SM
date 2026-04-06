(omega-dev-pgrad)=

# Pressure Gradient (PGrad)

Omega includes a `PressureGrad` class that computes horizontal pressure gradient
tendencies for the non-Boussinesq momentum equation. The implementation supports a
centered difference scheme as the default, with a placeholder for future high-order
methods. The class follows the same factory pattern used by other Omega modules.

## PressureGradType enum

An enumeration of the available pressure gradient schemes is defined in `PGrad.h`:

```c++
enum class PressureGradType { Centered, HighOrder1, HighOrder2 };
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
PGrad->computePressureGrad(Tend, State, VCoord, EqState, TimeLevel);
```

where:
- `Tend` is a 2D array `(NEdgesAll × NVertLayers)` that the pressure gradient
  tendency is accumulated into
- `State` is the current `OceanState`, from which layer thickness is extracted
  at the given `TimeLevel`
- `VCoord` provides pressure, interface height, and geopotential fields
- `EqState` provides the specific volume field
- `TimeLevel` selects which time level of the state to use

The method uses hierarchical Kokkos parallelism: an outer `parallelForOuter` loop
iterates over edges, and an inner `parallelForInner` loop iterates over vertical
chunks. The appropriate functor is dispatched based on `PressureGradChoice`.

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
                                const Array2DReal &ZInterface,
                                const Array1DReal &TidalPotential,
                                const Array1DReal &SelfAttractionLoading,
                                const Array2DReal &SpecVol) const;
```

### PressureGradHighOrder

This functor is a placeholder for a future high-order pressure gradient implementation
suitable for ice shelf cavities and complex bathymetry. Currently it performs no
computation (a no-op).

## Configuration

The pressure gradient type is selected in the input YAML file:

```yaml
PressureGrad:
   PressureGradType: 'centered'
```

Valid options for `PressureGradType` are:
- `'centered'` or `'Centered'`: centered difference approximation (default)
- `'HighOrder1'`: first high-order method (placeholder, future implementation)

If an unrecognized value is provided, the implementation falls back to the centered
scheme and logs an informational message.

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
| `HighOrderPGrad` | `PressureGradHighOrder` | High-order pressure gradient functor |
| `PressureGradChoice` | `PressureGradType` | Selected pressure gradient method |

## Removal

To remove all `PressureGrad` instances:

```c++
OMEGA::PressureGrad::clear();
```

To remove a specific named instance:

```c++
OMEGA::PressureGrad::erase("MyPGrad");
```
