(omega-user-horz-mesh)=

## Horizontal Mesh

Omega uses the MPAS mesh specification found
[here](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf). The names
of the mesh variables have been retained, with the caveat that they now begin
with a capital letter.

The Mesh class is meant to be a container for all mesh variables local to a
decomposed sub-domain that can be easily passed among the dycore routines. It
depends on a given [Decomp](#omega-user-decomp) and reproduces the
cell/edge/vertex totals and connectivity information from that class. The Mesh
class also creates parallel I/O decompositions that are used to read in the
additional mesh variables, which are not required for Decomp.

Currently, the Mesh class reads in all variables from the MPAS mesh
specification except those read by [Decomp](#omega-dev-decomp).
This includes the following variables:

| Variable Name | Description | Units |
| ------------- | ----------- | ----- |
| XCell, YCell, ZCell | Cartesian coordinates of cell centers | m |
| XEdge, YEdge, ZEdge | Cartesian coordinates of edge centers | m |
| XVertex, YVertex, ZVertex | Cartesian coordinates of vertices | m |
| BottomDepth | Depth of the bottom of the ocean at cell centers | m |
| FCell, FEdge, FVertex | Coriolis parameter at cell centers/edges/vertices | radians/s |
| LonCell, LatCell | Longitude/latitude coordinates of cell centers | radians |
| LonEdge, LatEdge | Longitude/latitude coordinates of edge centers | radians |
| LonVertex, LatVertex | Longitude/latitude coordinates of vertices | radians |
| AreaCell | Area of each cell | m^2 |
| AreaTriangle | Area of each triangle in the dual grid | m^2 |
| KiteAreasOnVertex | Area of the portions of each dual cell that are part of each cellsOnVertex | m^2 |
| DvEdge | Length of each edge, computed as the distance between verticesOnEdge | m |
| DcEdge | Length of each edge, computed as the distance between CellsOnEdge | m |
| AngleEdge | Angle the edge normal makes with local eastward direction | radians |
| MeshDensity | Value of density function used to generate a particular mesh at cell centers | - |
| WeightsOnEdge | Reconstruction weights associated with each of the edgesOnEdge | - |

In the future, the Mesh class will optionally compute the mesh variables that
are dependent on the Cartesian mesh coordinates internally.
This includes the various areas, lengths, angles, and weights needed for the
TRiSK discretization (e.g. rows 5-11 in the table above).
