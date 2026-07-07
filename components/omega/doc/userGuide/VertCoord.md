(omega-user-vert-coord)=

## Vertical Coordinate

### Overview

Omega uses pseudo height, $\tilde{z} = \frac{p}{\rho_0 g}$,  as the vertical coordinate ([V1 governing equation document](omega-design-governing-eqns-omega1)).
In practice, $\tilde{z}$ is not computed directly.
Instead, the model tracks and evolves the pseudo-thickness, $\tilde{h} = \Delta \tilde{z}$, defined as the difference between adjacent layer interfaces.
The pseudo height is essentially a normalized pressure coordinate, with the advantage that it has units of meters.
The `VertCoord` class contains variables and functions relevant to keeping track of:
 - the maximum possible number of layers in a cell, i.e. the extent of the vertical dimension (read in from the mesh file)
 - the bottom depth of each cell (read in from the mesh file)
 - the location of active vertical layers (used to set extents of vertical loop bounds):
   - the max and min indices of active layers in each cell (read in from the mesh file)
   - the max and min indices of active layers for edges and vertices at the bottom and top of the water column (computed from min/max cell layers)
 - pressure (computed from the pseudo-thickness and surface pressure)
 - geometric height (computed from the bottom depth, specific volume, and pseudo-thickness)
 - geopotential (computed from the geometric height and tidal forcing)
 - desired vertical interface locations (computed from pressure, reference layer pseudo-thickness, and user-specified weights)

Multiple instances of the vertical coordinate class can be created and accessed by a unique name.

### Variables

| Variable Name | Description | Units |
| ------------- | ----------- | ----- |
| NVertLayers   | maximum number of vertical layers | - |
| NVertLayersP1 | maximum number of vertical layers plus 1 | - |
| PressureInterface | pressure at layer interfaces | force per unit area at layer interfaces | kg m$^{-1}$ s$^{-2}$ |
| PressureMid | pressure at layer mid points | force per unit area at layer mid point | kg m$^{-1}$ s$^{-2}$ |
| GeomZInterface | geometric height of layer interfaces | m |
| GeomZMid | geometric height of layer midpoint | m |
| GeopotentialMid | geopotential at layer mid points | m$^2$/s$^2$|
| PseudoThicknessTarget | desired pseudo-thickness based on total perturbation from the reference pseudo-thickness | - |
| MinLayerCell | first active layer for cell | - |
| MaxLayerCell | last active layer for cell | - |
| MinLayerEdgeTop | min of the first active layers for cells on edge | - |
| MaxLayerEdgeTop | min of the last active layer for cells on edge | - |
| MinLayerEdgeBot | max of the first active layer for cells on edge | - |
| MaxLayerEdgeBot | max of the last active layer for cells on edge | - |
| MinLayerVertexTop | min of the first active layer for cells on vertex | - |
| MaxLayerVertexTop | min of the last active layer for cells on vertex | - |
| MinLayerVertexBot | max of the first active layer for cells on vertex | - |
| MaxLayerVertexBot | max of the last active layer for cells on vertex | - |
| VertCoordMovementWeights | weights to specify how total column thickness changes are distributed across layers | - |
| RefPseudoThickness | reference pseudo-thickness used to distribute total column thickness changes | m |
| BottomGeomDepth | positive down distance from the reference geoid to the bottom | m |

### Inactive and boundary layers in output

Each output field declares a fill value (matching the NetCDF-C standard `_FillValue`) that marks entries outside the active domain.
In a cell column only the active layers (`[MinLayerCell, MaxLayerCell]`) carry valid data; deeper, inactive layers contain the fill value.
Edge and vertex fields have an intermediate boundary range where some, but not all, of the surrounding cells are active.
For flux-type edge fields (such as `NormalVelocity`), the flux through a face bordering land vanishes, so those boundary layers carry `0`; layers fully outside the valid range carry the fill value.
Vertex fields, by contrast, keep their computed (generally non-zero) values throughout the valid range, since a boundary vertex still has active surrounding cells contributing to quantities like vorticity; only layers outside the valid range carry the fill value.
After an initial-condition or restart read, `VertCoord` enforces this pattern on the state and tracer fields so that analysis tools (ncview, Xarray, NCO, etc.) recognize and mask the inactive entries automatically.
