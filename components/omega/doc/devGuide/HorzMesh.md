(omega-dev-horz-mesh)=

## Horizontal Mesh

The OMEGA horizontal mesh uses the [MPAS Mesh
Specification](https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf).
A mesh object is created by using the constructor, which requires a Decomp object:
```c++
Err = OMEGA::Decomp::init();
OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
OMEGA::HorzMesh Mesh(DefDecomp);
```
The constructor replicates the subdomain mesh cell/edge/vertex counts and
connectivity information from Decomp so this information can be passed among the
computational rountines.  It then creates several parallel I/O decompositions
and reads in the remaining subdomain mesh information.  Finally, any mesh
information needed on the device to a device YAKL array from the host.  These
tasks are orgainzed into several private methods.  The variable names for host
arrays hare appended with in `H`, array variable names not ending in `H` are
device arrays.  The copy from host to device array is performed via:
```c++
AreaCell = AreaCellH.createDeviceCopy();
```
The device arrays are deallocated by the `HorzMesh::clear()` method, which is
necessary before calling `yakl::finalize`.
