(omega-dev-decomp)=

## Domain Decomposition (Decomp)

In order to run across nodes in a parallel computer, OMEGA subdivides the
horizontal domain into subdomains that are distributed across the machine
and communicate using message passing via the Message Passing Interface (MPI).
To decompose the domain, we utilize the
[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) library.
An OMEGA mesh is fully described by the
[MPAS Mesh Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf)
which we will reproduce here eventually. The Decomp class decomposes
the domain based on the index space and number of MPI tasks (currently
one subdomain per task). Once decomposed, the Decomp class holds all
of the index space information as described below.

Within OMEGA, a default decomposition of the mesh is first created
with the call:
```c++
 Decomp::init();
```
This must be called very early in the init process, just after initializing
the MachEnv, Config and IO.  Mesh information is first read using parallel
IO into an equally-spaced linear decomposition, then partitioned by METIS
into a more optimal decomposition. The parallel METIS implementation
(ParMETIS) will eventually be used to reduce the memory footprint for
high-resolution configurations, but currently we use the serial METIS library
for the partitioning.

METIS requires information about the connectivity in the mesh. In particular,
it needs the total number of cells and edges in the mesh and the connectivity
given by the CellsOnCell array that stores the indices of neighboring cells
for each cell. Once the cells have been partitioned, we determine multiple
layers of halo cells that will be needed to avoid excessive communcation with
remote neighboring subdomains. The CellsOnCell, EdgesOnCell
and VerticesOnCell arrays are then distributed according to this cell
partitioning. Edges and Vertices are then partitioned. Ownership of each
edge/vertex is determined by the first (valid) cell index in the CellsOnEdge
or CellsOnVertex for that edge/vertex. The remaining connectivity arrays
(CellsOnEdge, EdgesOnEdge, VerticesOnEdge, CellsOnVertex, EdgesOnVertex)
are the redistributed to match the edge and vertex distributions. Halos
are filled to ensure all necessary edge and vertex information for the
cell decomposition (and cell halos) are present in the subdomain.

After the call to the Decomp initialization routine, a Decomp named
Default has been created and can be retrieved with
```c++
OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
```
Once retrieved all Decomp members are public and can be accessed using
```c++
OMEGA::I4 NCells = DefDecomp->NCells;
OMEGA::HostArray1DI4 CellIDH = DefDecomp->CellIDH;
```
Decomp is a container for all mesh index and connectivity arrays as
described in the mesh specification above. In particular, it contains
  - NCellsGlobal: the total number of cells
  - NCellsAll: the total number of cells in the local partition
  - NCellsSize: the length of the Cell array axis (typicall NCellsAll+1)
  - NCellsOwned: the number of cells owned by this task
  - NCellsHalo(i): the number of owned+halo cells for each halo layer
  - Analogous size variables for Edges and Vertices
  - MaxEdges: the max number of edges on a cell (and array size)
  - VertexDegree: the number of cells/edges at each vertex
  - CellID(NCellsSize): the global index for each local cell
  - CellLoc(NCellsSize,2): the task and local index for every local cell
    (redundant for owned cells but gives the remote location for halo cells)
  - Analogous arrays for Edges, Vertices
  - CellsOnCell(NCellsSize,MaxEdges): the local index for neighbor cells
    across each cell edge
  - EdgesOnCell(NCellsSize,MaxEdges): the local edge index for each cell edge
  - VerticesOnCell(NCellsSize,MaxEdges): the local index for each cell vertex
  - CellsOnEdge(NEdgesSize,2): the index of the 2 cells sharing an edge
  - VerticesOnEdge(NEdgesSize,2): the vertex index at each edge endpoint
  - EdgesOnEdge(NEdgesSize,MaxEdges x 2): the index for all edges contributing
    to an edge (all edges of the 2 cells sharing an edge)
  - CellsOnVertex(NVerticesSize,VertexDegree): indices for all cells meeting
    at a vertex
  - EdgesOnVertex(NVerticesSize,VertexDegree): indices for all edges meeting
    at a vertex
  - NEdgesOnCell(NCellsSize): the number of actual edges on each cell
  - NEdgesOnEdge(NEdgesSize): the number of actual edges on each edge

For each of the arrays above, there is a copy of the array on the host and
device (GPU) with the host array named with an extra H on the end
(eg CellsOnCellH). All are Kokkos arrays so are accessed with (index) rather
than [index] and for some of the arrays noted above are multi-dimensional.
A typical host loop might then look something like:
```c++
for (int Cell = 0; Cell < NCellsOwned; ++Cell) {
   for (int Edge = 0; Edge < NEdgesOnCell(Cell); ++Edge) {
      VarH(Cell) = VarH(Cell) + FluxH(Cell,Edge);
   }
}
```
And on the device, we use `OMEGA::parallelFor` in place of `Kokkos::parallel_for`;
```c++
OMEGA::parallelFor( {NCellsOwned,MaxEdges},
                       KOKKOS_LAMBDA (int Cell, int Edge) {
  if (Edge < NEdgesOnCell(Cell)) {
      Var(Cell) = Var(Cell) + Flux(Cell,Edge);
  }
});
```

Any defined decomposition can be removed by name using
```c++
Decomp::erase(Name);
```
and all decompositions *must* be removed before the Kokkos finalize call using
```c++
Decomp::clear();
```
which destroys all host and device arrays before Kokkos finalizes and removes
the memory pool in which all the arrays are allocated. Failure to call clear
before `Kokkos::finalize()` will result in an error.
