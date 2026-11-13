(omega-dev-horz-mesh)=

## Horizontal Mesh

The Omega horizontal mesh uses the [MPAS Mesh
Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf).
A mesh object is created by the `init` method, which assumes that `Decomp` has
already been initialized.
```c++
OMEGA::Decomp::init();
OMEGA::HorzMesh::init(ModelClock);
```
The model clock should be the ocean model clock and is required by the IOStream
routines that read the HorzMesh variables from the input stream
[HorzMeshIn](#omega-user-horz-mesh).
The constructor replicates the subdomain mesh cell/edge/vertex counts and
connectivity information from Decomp so this information can be passed among the
computational routines, alongside the other local mesh information.  It then
reads most of the remaining mesh variables (coordinates, lengths, areas)
from the input HorzMeshIn IOStream. A few variables, like masks and operator
scaling are computed internally; the mesh scaling is relative to
the global median cell width, which is computed during initialization.
After reading, the init routine fills all
halos and synchronizes the host and device copies of the mesh variables.
These tasks are organized into several private methods. Eventually,
dependent mesh variables will be computed from the minimum set of required mesh
information.

After initialization, the default mesh object can be retrieved via:
```
OMEGA::HorzMesh *HMesh = OMEGA::HorzMesh::getDefault();
```
Once this retrieval has been performed, public member variables can be accessed
using:
```
OMEGA::I4 NCellsOwned = HMesh->NCellsOwned;
OMEGA::Array1DR8 AreaCell = HMesh->AreaCell;
```

The HorzMesh is meant to be a container that allows the mesh information to be
passed to the PDE solver routines:
```
void computeFluxTendency(OMEGA::HorzMesh *HMesh, ...) {
OMEGA::parallelFor({HMesh->NCellsOwned,HMesh->MaxEdges},
                      KOKKOS_LAMBDA (int Cell, int Edge) {
  if (Edge < HMesh->NEdgesOnCell(Cell)) {
      Var(Cell) = Var(Cell) + Flux(Cell,Edge);
  }
}
```

`MaxEdges` is read from the mesh file and is therefore only known at run time.
Kernels that need a fixed-size per-thread array indexed by edge should size it
with the compile-time bound instead:
```
OMEGA::Real Weights[OMEGA::HorzMesh::MaxEdgesBound];
```
`MaxEdgesBound` is a single shared upper bound on `MaxEdges`; meshes exceeding
it are rejected. Do not introduce a local copy of this limit.

For member variables that are host arrays, variable names are appended with an
`H`.  Array variable names not ending in `H` are device arrays.

The device arrays are deallocated by the `HorzMesh::clear()` method, which is
necessary before calling `Kokkos::finalize`. Because Field and IOStream also
have references to HorzMesh variables, the `IOStream::clear()` and
`Field::clear()` should also be called to ensure all device arrays are removed.
