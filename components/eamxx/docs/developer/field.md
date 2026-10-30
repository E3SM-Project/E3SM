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

## Looping over the entries of fields

Writing a kernel that touches the data of a field is less straightforward than it
looks. The `Kokkos::View` that `get_view` returns has a type that depends on the
field's rank, on its scalar type, and on its layout (a `LayoutRight` view if the field
is contiguous in memory, a `LayoutStride` one if it is, for instance, a slice of a
larger field). These are runtime properties of a field, but a kernel needs the
view type at compile time. To bridge the two, EAMxx has `field_pfor`, which turns
a field into a view of the right type and runs a `Kokkos::parallel_for` over all
its entries, for any rank in `[0,6]`:

```cpp
#include "share/field/field_pfor.hpp"

// A kernel is a functor with a templated operator(). It receives the views first,
// then one index per dimension of the fields
struct AddKernel {
  template<typename YView, typename XView, typename... Idx>
  KOKKOS_INLINE_FUNCTION
  void operator() (const YView& y, const XView& x, const Idx&... i) const {
    y(i...) += x(i...);
  }
};

// y += x, whatever the rank and layout of the two fields
using details::view_as;
details::field_pfor(AddKernel(), view_as<Real>(y_field), view_as<const Real>(x_field));
```

`view_as<T>(field)` says which scalar type to use for that field. It must match the
field's data type, and it can be `const`-qualified for fields that the kernel only
reads. The kernel is called as `func(v_1,...,v_K,i_0,...,i_{N-1})`, with `v_k` the view of
the `k`-th field, and `N` the rank of the fields. Because the same kernel
runs for any rank, write the indexing as `y(i...)`: for a rank-0 field, the pack `i...`
is empty. The members of the kernel struct are how you pass parameters (e.g., scalars
or other views), and they are copied to device. A kernel can still branch on the
number of indices with `if constexpr (sizeof...(i)==2)`.

Some things to keep in mind:

* All fields must have the same rank and extents (the tags are not compared), and
  at most three fields can be passed.
* Like `Kokkos::parallel_for`, the call is asynchronous: it does not fence. Kernels launched on
  the same execution space run in order, so consecutive calls can be chained without fencing.
  Fence (e.g., with `Kokkos::fence()`, or by calling `sync_to_host()`) before using the results on host.
* The kernel must be a **functor struct, not a `KOKKOS_LAMBDA`**. nvcc does not allow
  generic extended `__host__ __device__` lambdas (needed here, since the view types are
  not known when the lambda is written), nor extended lambdas whose captures are first used inside an
  `if constexpr`. A struct with a templated `operator()` has neither restriction. Two more
  nvcc restrictions on extended lambdas are easy to run into: they cannot be defined inside
  a generic lambda, nor inside a private or protected member function. (In unit tests, we also
  define lambdas before, and not inside, the argument of a `REQUIRE`, since Catch may evaluate that
  argument in an unevaluated context, where nvcc does not accept them.)
* The kernel is instantiated for every rank in `[0,6]` and for every combination of
  contiguous/strided views of the fields. With `K` fields, that is `7 * 2^K` instantiations
  for each call site, so keep the kernel small (if needed, have it forward to a templated
  helper). The `Field` update methods that use `field_pfor` are explicitly instantiated in the
  ETI files (`share/field/eti/`), so that this cost is paid once, rather than in every file
  that includes the header.

For examples, see the update kernels in `share/field/field_update_impl.hpp` (`UpdateKernel`),
`field_update_masked.hpp`, and `field_update_fill_aware.hpp`, and the unit test
`share/field/tests/field_pfor_tests.cpp`.

### Under the hood

Kokkos calls the functor of a `parallel_for` with the indices only, so the views cannot be
passed directly to the kernel. `field_pfor` uses a few generic utilities (in `share/util`,
none of which depends on fields) to bridge the gap:

* `BoundFunctor` (`eamxx_bound_functor.hpp`): a device-friendly version of `std::bind_front`.
  `bind_functor(func,a,b)` stores `func` and `(a,b)`, and when called with `(i,j)`
  it calls `func(a,b,i,j)`. The call arguments are perfectly forwarded, so it also works
  with `parallel_reduce`, where Kokkos passes the reduction result after the indices.
  `field_pfor` binds the visited views to the kernel.
* `Tuple` (`eamxx_tuple.hpp`, in `scream::details`): a minimal tuple which can be used
  in device code (`std::tuple` and `std::apply` are not guaranteed to work in kernels),
  where `BoundFunctor` keeps the bound arguments.
* `parallel_for_nd` (`eamxx_parallel_for_nd.hpp`): launches the right Kokkos policy
  (a `RangePolicy` for rank 0 and 1, an `MDRangePolicy` otherwise) for the extents of the fields.
  For rank 0, where Kokkos has no 0-dimensional policy, it drops the dummy index.

`field_pfor` is for elementwise kernels. There is no equivalent yet for reductions or for
kernels that use hierarchical parallelism.
