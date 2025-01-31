(omega-dev-horz-mesh)=

## Horizontal Mesh

The Omega horizontal mesh uses the [MPAS Mesh
Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf).
A mesh object is created by the `init` method, which assumes that `Decomp` has
already been initialized.
```c++
Err = OMEGA::Decomp::init();
Err = OMEGA::HorzMesh::init();
```
The constructor replicates the subdomain mesh cell/edge/vertex counts and
connectivity information from Decomp so this information can be passed among the
computational routines, alongside the other local mesh information.  It then
creates several parallel I/O decompositions and reads in the remaining subdomain
mesh information.  Finally, any mesh information needed on the device is copied
from the host to a device Kokkos array. Arrays such as the coordinate variables,
which are not involved in tendency calculations, are not transferred to the
device. These tasks are organized into several private methods. Eventually,
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

For member variables that are host arrays, variable names are appended with an
`H`.  Array variable names not ending in `H` are device arrays.  The copy from
host to device array is performed in the constructor via:
```c++
AreaCell = OMEGA::createDeviceMirrorCopy(AreaCellH);
```

The device arrays are deallocated by the `HorzMesh::clear()` method, which is
necessary before calling `Kokkos::finalize`.
