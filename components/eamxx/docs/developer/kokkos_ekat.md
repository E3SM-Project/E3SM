# Kokkos and EKAT

## Kokkos

EAMxx uses Kokkos for performance portable abstractions for parallel execution
of code and data management to various HPC platforms, including OpenMP, Cuda,
HIP, and SYCL. Here we give a brief overview of the important concepts for
understanding Kokkos in EAMxx. For a more in depth description, see the
[Kokkos wiki](https://kokkos.org/kokkos-core-wiki).

### Kokkos::Device

`Kokkos::Device` is a struct which contain the type definitions for two main
Kokkos concepts: execution space (`Kokkos::Device::execution_space`), the place
on-node where parallel operations (like for-loops, reductions, etc.) are
executed, and the memory space (`Kokkos::Device::memory_space`), the memory
location on-node where data is stored. Given your machine architecture, Kokkos
defines a default "device" space, given by

```cpp
Kokkos::Device<Kokkos::DefaultExecutionSpace,
         Kokkos::DefaultExecutionSpace::memory_space>
```

where all performance critical code should be executed (e.g., on an NVIDIA
machine, this device would be the GPU accelerators) and a default "host" space,
given by

```c++
Kokkos::Device<Kokkos::DefaultHostExecutionSpace,
         Kokkos::DefaultHostExecutionSpace::memory_space>
```

where data can be accessed by the CPU cores and is necessary for I/O
interfacing, for example. Currently, these default spaces are the ones used by
EAMxx. On CPU-only machines, host and device represent the same space.

### Kokkos Views

The main data struct provided by Kokkos used in EAMxx in the `Kokkos::View`.
This is a multi-dimensional data array that can live on either device or host
memory space. These Views are necessary when running on GPU architectures as
data structures like `std::vector` and `std::array` will be unavailable on
device.

Views are constructed in EAMxx most commonly with the following template and
input arguments

```cpp
Kokkos::View<DataType, LayoutType, DeviceType>(const std::string& label,
                         int dim0, int dim1, ...)
```

where

- `DataType`: scalar type of the view, given as `ScalarType`+`*`(x's number of
  run-time dimensions). E.g., a 2D view of doubles will have `DataType =
  double**`. There is also an ability to define compile-time dimensions by
  using `[]`, see
  [Kokkos wiki section on views]
  (<https://kokkos.org/kokkos-core-wiki/API/core/view/view.html).
- `LayoutType`: mapping of indices into the underlying 1D memory storage. Types
  are:
  - `LayoutRight` (used in EAMxx): strides increase from the right most to the
  left most dimension, right-most dimension is contiguous
  - `LayoutLeft`: strides increase from the left most to the right most
  dimension, left-most dimension is contiguous
  - `LayoutStride`: strides can be arbitrary for each dimension
- `DeviceType`: provides space where data live, defaults to the default device

The following example defines a view "temperature" which has dimensions columns
and levels:

```cpp
Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultDevice> temperature(
  "temperature", ncols, nlevs);
```

### Deep Copy

Kokkos provides `Kokkos::deep_copy(dst, src)` which copies data between views
of the same dimensions, or a scalar values into a view. Common uses

```cpp
Kokkos::deep_copy(view0, view1); // Copy all data from view1 into view0
Kokkos::deep_copy(view0, 5); // Set all values of view0 to 5
```

As seen in the next section, we can use `deep_copy()` to copy data between host
and device.

### Mirror Views

We will often need to have memory allocation the resides on device (for
computation), and then need that identical data on host (say, for output).
Kokkos has a concept of mirror views, where data can be copied from host to
device and vice versa.

Here is an example using the device view `temperature` from above

```cpp
// Allocate view on host that exactly mirrors the size of layout of the device
view
auto host_temperature = Kokkos::create_mirror_view(temperature);

// Copy all data from device to host
Kokkos::deep_copy(host_temperature, temperature);
```

Kokkos also offers an all-in-one option

```cpp
// Note: must hand the host device instance as first argument
auto host_temperature = Kokkos::create_mirror_view_and_copy(
  Kokkos::DefaultHostDevice(), temperature);
```

### Parallel Execution

The most basic parallel execution pattern used by EAMxx is the
`Kokkos::parallel_for` which defines a for-loop with completely independent
iterations. The `parallel_for` takes in an optional label for debugging, an
execution policy, which defines a range and location (host or device) for the
code to be run, and a lambda describing the body of code to be executed. The
following are execution policies used in EAMxx

- `int count`: 1D iteration range `[0, count)`
- `RangePolicy<ExecSpace>(int beg, int end)`: 1D iteration range for indices
  `[beg, end)`
- `MDRangePolicy<ExecSpace, Kokkos::Rank<N>>(int[N] beg, int[N] end)`: multi-
  dimensional iteration range `[beg, end)`
- `TeamPolicy<ExecSpace>(int league_size, int team_size, int vector_size)`: 1D
  iteration over `league_size`, assigned to thread teams of size `team_size`,
  each with `vector_size` vector lanes. Both `team_size` and `vector_size` can
  be given `Kokkos::AUTO` as input for Kokkos to automatically compute.

If no `ExecSpace` template is given, the default execution space is used.

For lambda capture, use `KOKKOS_LAMBDA` macro which sets capture automatically
based on architecture.

Example using `RangePolicy` to initialize a view

```cpp
Kokkos::View<double**, Kokkos::LayoutRight> temperature("temperature", ncols,
                            nlevs);
Kokkos::parallel_for("Init_temp",
           Kokkos::RangePolicy(0, ncols*nlevs),
           KOKKOS_LAMBDA (const int idx) {
  int icol = idx/nlevs;
  int ilev = idx%nlevs;

  temperature(icol, ilev) = 0;
});
```

Same example with `TeamPolicy`

```cpp
Kokkos::parallel_for("Init_temp",
           Kokkos::TeamPolicy(ncols*nlevs, Kokkos::AUTO, Kokkos::AUTO),
           KOKKOS_LAMBDA (const TeamPolicy::member_type& team) {
  // league_rank() gives the index for this team
  int icol = team.league_rank()/nlevs;
  int ilev = team.league_rank()%nlevs;

  temperature(icol, ilev) = 0;
});
```

### Hierarchical Parallelism

Using `TeamPolicy`, we can have up to three nested levels of parallelism: team
parallelism, thread parallelism, vector parallelism. These nested policies can
be called within the lambda body using the following execution policies

- `TeamThreadRange(team, begin, end)`: execute over threads of a team
- `TeamVectorRange(team, begin, end)`: execute over threads and vector lanes of
  a team
- `ThreadVectorRange(team, begin, end)`: execute over vector lanes of a thread

An example of using these policies

```cpp
Kokkos::View<double***> Q("tracers", ncols, ntracers, nlevs);
Kokkos::parallel_for(Kokkos::TeamPolicy(ncols, Kokkos::AUTO),
           KOKKOS_LAMBDA (TeamPolicy::member_type& team) {
  int icol = team.league_rank();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&](int ilev) {
  temperature(icol, ilev) = 0;
  });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlevs), [&](int ilev) {
  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, ntracers), [&](int iq) {
    Q(icol, iq, ilev) = 0;
  });
  });
});
```

IMPORTANT! Nested policies cannot be used in arbitrary order. `ThreadVectorRange`
must be used inside a `TeamThreadRange`, and `TeamVectorRange` must be the only
level of nested parallelism.

```cpp
Kokkos::parallel_for(TeamPolicy(...), ... {
   // OK
   Kokkos::parallel_for(TeamThreadRange, ... {

   });

   // OK
   Kokkos::parallel_for(TeamThreadRange, ... {
   Kokkos::parallel_for(ThreadVectorRange, ... {

   });
   });

   // OK
   Kokkos::parallel_for(TeamVectorRange, ...{

   });

   // WRONG,ThreadVectorRange must be nested in TeamThreadRange
   Kokkos::parallel_for(ThreadVectorRange, ...{

   });

   // WRONG, a TeamVectorRange must be the only nested level
   Kokkos::parallel_for(TeamVectorRange, ...{
   Kokkos::parallel_for(ThreadVectorRange, ... {

   });
   });
});
```

Using these incorrectly can be very tricky to debug as the code almost certainly
will _not_ error out, but race conditions will exist among threads.

## EKAT

### KokkosTypes

### ExeSpaceUtils

### Vectorization: Packs

### Scratch Memory: WorspaceManager

### Algorithms
