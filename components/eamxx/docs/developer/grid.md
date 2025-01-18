# Grids and Remappers

In EAMxx, the `AbstractGrid` is an interface used to access information regarding
the horizontal and vertical discretization. The most important information that
the grid stores is:

* the number of local/global DOFs: these are the degrees of freedom of the
  horizontal grid only. Here, local/global refers to the MPI partitioning.
* the DOFs global IDs (GIDs): a list of GIDs of the DOFs on the current MPI rank,
  stored as a Field
* the local IDs (LIDs) to index list: this list maps the LID of a DOF (that is,
  the position of the DOF in the GID list) to a "native" indexing system for that
  DOF. For instance, a `PointGrid` (a class derived from `AbstractGrid`) is a
  simple collection of points, so the "native" indexing system coincides with the
  LIDs. However, for a `SEGrid` (a derived class, for spectral element grids),
  the "native" indexing is a triplet `(ielem,igp,jgp)`, specifying the element
  index, and the two indices of the Gauss point within the element.
* geometry data: stored as a `std::map<std::string,Field>`, this represent any
  data that is intrinsically linked to the grid (either along the horizontal or
  vertical direction), such as lat/lon coordinates, vertical coordinates, area
  associated with the DOF.

Grids can also be used to retrieve the layout of a 2d/3d scalar/vector field,
which allows certain downstream classes to perform certain operations without
assuming anything on the horizontal grid.

In general, grid objects are passed around the different parts of EAMxx as const
objects (read-only). The internal data can only be modified during construction,
which usually is handled by a `GridsManager` object.
