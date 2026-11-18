(omega-dev-submeso-eddies)=

# Submesoscale Eddy Parameterization

Omega includes a `SubmesoEddies` class (`src/ocn/SubmesoEddies.{h,cpp}`) that implements the
submesoscale mixed-layer instability closure used to produce an eddy-induced
normal transport velocity. It provides methods to compute buoyancy gradients, and eddy velocity.
The current implementation follows the Fox-Kemper et al. (2011) (FK11) closure,
as described in the
{ref}`omega-design-submesoscale-eddies` design document.

## Initialization

The `SubmesoEddies` class is implemented as a singleton. Before creating it,
[`HorzMesh`](#omega-dev-horz-mesh) and [`VertCoord`](#omega-dev-vert-coord)
must be initialized. Create the instance with the static method
```c++
SubmesoEddies::init();
```
Retrieve the pointer at any time with:
```c++
SubmesoEddies* DefSubEddies = SubmesoEddies::getInstance();
```

## Data and algorithms

Core public fields include:

- `GradBuoyEdgeInterface`
- `EddyVelocity`

Key implementation details:

- `computeBuoyGrad` computes horizontal buoyancy gradients at edges and adds
  the tilted-coordinate correction using `BruntVaisalaFreqSq`.
- `computeEddyVelocity` forms mixed layer averaged buoyancy and stratification
  terms, evaluates frontal-width limits (`LfMin`, `DsMax`), computes a
  streamfunction with `shapeFunction()`, and takes its vertical divergence to
  produce edge-normal eddy velocity.

## Computation of buoyancy gradient
To compute buoyancy gradient `GradBuoyEdgeInterface` from specific volume
`SpecVol`, mean pseudo-thickness on edges `MeanPseudoThickEdge`, mid-layer
`z` coordinate `GeomZMid`, and squared Brunt-Vaisala frequency
`BruntVaisalaFreqSq`, use
```c++
SubEddies.computeBuoyGrad(SpecVol, MeanPseudoThickEdge, GeomZMid, BruntVaisalaFreqSq);
```

## Computation of eddy velocity
To compute eddy velocity array `EddyVelocity` from mixed layer depth `DenMixLayerDepth`,
mixed layer index `DenMixLayerIndex`, squared Brunt-Vaisala
frequency field `BruntVaisalaFreqSq` and mean pseudo-thickness on edges
`MeanPseudoThickEdge`, use
```c++
SubEddies.computeEddyVelocity(DenMixLayerDepth, DenMixLayerIndex, BruntVaisalaFreqSq, MeanPseudoThickEdge);
```

## Finalization
To clear the singleton instance, use the static method
```c++
SubmesoEddies::destroyInstance();
```
