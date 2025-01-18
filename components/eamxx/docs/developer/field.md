# Field

In EAMxx, a `Field` is a data structure holding two things: pointers to the
data and pointers to metadata. Both the data and metadata are stored in
`std::shared_ptr` instances, to ensure consistency across all copies of
the field. This allows for fast shallow copy semantic for this class.

The data is stored on both CPU and device memory (these may be the same,
depending on the Kokkos backend). In EAMxx, we always assume and guarantee
that the device data is up to date. That implies that the data must be
explicitly synced to host before using it on host, and explicitly synced
to device after host manipulation, in order to ensure correctness.
In order to access the data, users must use the `get_view`/
`get_strided_view` methods, which takes two template arguments:
the data type, and an enum specifying whether CPU or device data is needed.
The data type is used to reinterpret the generic pointer stored inside
to a view of the correct scalar type and layout. It is a possibly
const-qualified type, and if the field was marked as "read-only",
the method ensures that the provided data type is const. A read-only field
can be created via the `getConst` method, which returns a shallow copy of
the field, but marked as read-only. The enum specifying host or device data
is optional, with device being the default.

The metadata is a collection of information on the field, such as name, layout, units,
allocation size, and more. Part of the metadata is immutable after creation (e.g.,
name, units, or layout), while some metadata can be partially or completely modified.
The metadata is contained in the `FieldHeader` data structure, which contains four
parts:

* `FieldIdentifier`: stores the field's name, layout, units, data type,
  and name of the grid where it's defined. These information are condensed
  in a single string, that can be used to uniquely identify a field, allowing
  to distinguish between different version of the same field.
  The layout is stored in the `FieldLayout` data structure, which includes:
  * the field tags: stored as a `std::vector<FieldTag>`, they give context to the
    field's extents.
  * the field dims: stored both as a `std::vector<int>`, as well as a 1d `Kokkos::View`.
* `FieldTracking`: stores information on the usage of the field, as well as its
  possible connections to other fields. In particular, the tracked items are:
  * the field time stamp: the time stamp when the field was last updated.
  * the field accumulation start time: used for fields that are accumulated over
    several time steps (or time step subcycles). For instance, it allows to
    reconstruct fluxes from raw accumulations.
  * the providers/customers: lists of atmosphere processes (see below) that
    respectively require/compute the field in their calculations.
  * the field groups: a list of field groups that this field belongs too. Field groups
    are used to access a group of fields without explicit prior knowledge about the
    number and/or names of the fields.
* `FieldAllocProp`: stores information about the allocation. While the field is not
  yet allocated, users can request special allocations for the field, for instance
  to accommodate packing (for SIMD), which may require padding. Upon allocation,
  this information is then used by the Field structure to extract the actual data,
  wrapped in a properly shaped `Kokkos::View`. The alloc props are also
  responsible of tracking additional information in case the field is a "slice" of
  a higher-dimensional one, a fact that can affect how the data is accessed.
* Extra data: stored as a `std::map<std::string,ekat::any>`, allows to catch any
  metadata that does not fit in the above structures. This is a last resort structure,
  intended to accommodate the most peculiar corner cases, and should be used sparingly.
