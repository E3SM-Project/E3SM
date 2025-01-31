(omega-user-decomp)=

## Domain Decomposition (Decomp)

Omega is designed to be run in parallel across multiple nodes and cores
in a clustered architecture. As in many Earth System Models, we decompose
the horizontal domain into subdomains that are distributed across the
system. Communication between the domains is accomplished using the
Message Passing Interface standard. To reduce the communication between
subdomains, a halo of points is filled with information from remote neighbor
cells, edges and vertices so that a time step can largely be completed
without the need to communicate.

The partitioning itself is performed using the
[METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) library that
can partition unstructured domains in a way that also optimizes communication.
More details on the mesh, connectivity and partitioning can be found in
the [Developer's Guide](#omega-dev-decomp).

There are three parameters that are set by the user in the input configuration
file. These are:
```yaml
Decomp:
   HaloWidth: 3
   MeshFileName: OmegaMesh.nc
   DecompMethod: MetisKWay
```
(until the config module is complete, these are currently hardwired to
the defaults above). The HaloWidth is set to be able to compute all of the
baroclinic terms in a timestep without communication and for higher-order
tracer advection terms, this currently must be at least 3. The MeshFileName
should include the complete path and filename to a standard Omega mesh file
that contains at a minimum
  - the total number of cells, edges and vertices (NCells, NEdges, NVertices)
  - the mesh connectivity contained in the arrays CellsOnCell, EdgesOnCell
    VerticesOnCell, CellsOnEdge, EdgesOnEdge, CellsOnVertex, EdgesOnVertex.
Again, a full description of the mesh is given in the
[Developer's Guide](#omega-dev-decomp).
The mesh information is read via parallel IO into an initial linear domain
decomposition and then is partitioned by METIS and rearranged into the
final METIS parallel decomposition.

METIS and ParMETIS support a number of partitioning schemes, but MetisKWay
is currently the only supported decomposition method for Omega and is
generally the better option.

Once the mesh is decomposed, all of the mesh index arrays are stored in
a Decomp named Default which can be retrieved as described in the
Developer guide. In the future, additional decompositions associated
with processor subsets (as described in MachEnv) but this capability is
not yet supported.
