# RDycore Code Structure and Organization

In this section, we describe the organization of RDycore's source files. We have
attempted to group functions and data types together into subsystems that can
be easily understood separately from one another.

Because RDycore's behavior is determined entirely by input parameters specified
in the [YAML input format](../common/input.md), it's a good idea to familiarize
yourself with this format as you read this section.

## A Brief Tour of RDycore's Source Tree

The root level of the [RDycore](https://github.com/RDycore/RDycore) source tree
contains a `README.md` file, a `CMakeLists.txt` file that defines the build
system, and some configuration files for documentation and tooling.

There are several folders at the root of the source tree:

```
+ RDyore +- .github/workflows
         |
         +- cmake
         |
         +- config
         |
         +- docs
         |
         +- driver
         |
         +- external
         |
         +- include
         |
         +- share
         |
         +- src
         |
         +- tools
```

* [`.github/workflows`](https://github.com/RDycore/RDycore/tree/main/.github/workflows): workflows
  that support our [GitHub Continuous Integration
  Environment](development.md#GitHub-Continuous-Integration-Environment)
* [`cmake`](https://github.com/RDycore/RDycore/tree/main/cmake): CMake scripts
  that support the build system, testing, etc.
* [`config`](https://github.com/RDycore/RDycore/tree/main/config): shell scripts
  to help with running RDycore on our target platforms
* [`docs`](https://github.com/RDycore/RDycore/tree/main/docs): Markdown source
  files for our [mkdocs](https://squidfunk.github.io/mkdocs-material/)-based
  documentation
* [`driver`](https://github.com/RDycore/RDycore/tree/main/driver): C and Fortran
  driver programs, including the standalone RDycore drivers and a few other
  development tools
* [`external`](https://github.com/RDycore/RDycore/tree/main/external): third-
  party libraries used by RDycore that aren't provided by PETSc
* [`include`](https://github.com/RDycore/RDycore/tree/main/include): the
  `rdycore.h.in` template that generates the `rdycore.h` API header file, and
  all private header files (located within the `private` subfolder)
* [`share`](https://github.com/RDycore/RDycore/tree/main/share): data files for
  initial conditions, material properties, and unstructured grids (used mainly
  for testing)
* [`src`](https://github.com/RDycore/RDycore/tree/main/src): the source code for
  RDycore
* [`tools`](https://github.com/RDycore/RDycore/tree/main/src): miscellaneous
  tools and scripts for building and deploying Docker images

Take a look at each of these folders to familiarize yourself with their
contents. In particular, the `src` folder has a few important subfolders:

* `src/f90-mod`: the RDycore Fortran module that mirrors the C library's
  capabilities
* `src/swe`: functions and data structures specific to the solution of the
  2D shallow water equations
* `src/tests`: unit tests for subsystems within RDycore

The private headers in the `include` folder contain definitions of opaque types
and functions that are helpful throughout RDycore but are not part of the API.

The `driver` folder also has its own `tests` subfolder that defines tests
that run the standalone drivers and perform convergence tests for selected
problems.

## The RDycore Object and its Lifecycle

RDycore (the "river dynamical core") is represented in code by a data structure
called `RDy`, declared as an opaque type in [include/rdycore.h](https://github.com/RDycore/RDycore/tree/main/include/rdycore.h).
It is defined in [include/private/rdycoreimpl.h](https://github.com/RDycore/RDycore/tree/main/include/private/rdycoreimpl.h).
In this section, we describe how to manipulate `RDy` to set up and run
simulations. It may be helpful to refer to the [standalone C driver](https://github.com/RDycore/RDycore/tree/main/driver/main.c)
(and or the corresponding [Fortran driver](https://github.com/RDycore/RDycore/tree/main/driver/main.F90))
as you read this section.

1. **Initialization**

Before RDycore can be used in any program, the supporting systems must be
initialized with a call to `RDyInit`:

```
  // initialize RDycore subsystems
  PetscCall(RDyInit(argc, argv, helpString));
```

This function accepts the `argc` and `argv` command line argument parameters
accepted by any C program, plus a string printed to display help/usage
information.

Next, create an `RDy` object, passing it an MPI communicator, the path to a
[YAML input file](../common/input.md), and the pointer to the desired `RDy`
object:

```
  // create an RDy on the given communicator with the given input
  RDy rdy;
  PetscCall(RDyCreate(MPI_COMM_WORLD, "my-input.yaml", &rdy));
```

This only creates the object--it does not read input or allocate any resources
for the problem described within the input.

2. **Setup**

To read the input file and ready your `RDy` object to run the simulation
defined therein, call `RDySetup`:

```
  // read input and set up the dycore
  PetscCall(RDySetup(rdy));
```

After this call, you're ready to run a simulation. You can also all any query
or utility functions on your `RDy` object to do whatever you need for your own
simulation. For example, you might need information about the unstructured grid,
or perhaps you are specifying specific boundary values for a Dirichlet boundary
condition. See the API header file for all the possibilities.

3. **Timestepping**

The simplest way to start running your simulation after setup is to run
`RDyAdvance` (which takes a single step) within a `while` loop that uses
`RDyFinished` as termination condition:

```
  // run the simulation to completion
  while (!RDyFinished(rdy)) {
    // ... pre-step logic here ...

    // advance the simulation by a single step
    PetscCall(RDyAdvance(rdy));

    // ... post-step logic here ...
  }
```

If you look in the standalone C driver, you can see various pre- and post-step
logic to accommodate data transfer, restarts, etc.

4. **Finalization**

At the end of an RDycore-enabled program, you must destroy the `RDy` object with
a call to `RDyDestroy`:

```
  // destroy the dycore
  PetscCall(RDyDestroy(&rdy));
```

Finally, call `RDyFinalize` to reclaim all resources used by RDycore and its
subsystems:

```
  // clean up
  PetscCall(RDyFinalize());
```

## Computational Domain

RDycore solves equations on a topologically-two-dimensional domain with
topography represented by an elevation coordinate $z(x, y)$.

RDycore uses two objects to represent the computional domain:

* [DMPlex](https://petsc.org/release/manual/dmplex/): a PETSc data structure
  that stores the minimum set of information needed to define an unstructured
  grid
* [RDyMesh](https://github.com/RDycore/RDycore/blob/main/include/private/rdymeshimpl.h):
  An RDycore-specific data structure (implementation [here](https://github.com/RDycore/RDycore/blob/main/src/rdymesh.c))
  that stores data for cells, edges, and vertices, constructed using the
  `DMPlex` object

Both of these objects are created by a call to [RDySetup](https://github.com/RDycore/RDycore/blob/main/src/rdysetup.c).
The `RDyMesh` object is a simple container with

* metadata like numbers of cells, edges, vertices
* a `cells` field that stores cell data in an `RDyCells` data structure
* an `edges` field that stores edge data in an `RDyEdges` data structure
* a `vertices` field that stores vertex data in an `RDyVertices` data structure
* metadata useful for output/visualization

The domain is partitioned into disjoint _regions_, each of which consists of a
set of contiguous _cells_, usually triangles or quads. Additionally, the domain
is bounded by one or more _boundaries_, each of which consists of a set of
edges (which can be thought of as two-dimensional "faces"). We describe how
regions and boundaries work below.

### Regions

When a mesh file is read to create a `DMPlex` object during [RDySetup](https://github.com/RDycore/RDycore/blob/main/src/rdysetup.c),
[DMLabel](https://petsc.org/release/manualpages/DMLabel/) objects are created
that represent disjoint sets of cells. Each cell set represents a region. From
each label/cell set, RDycore constructs an [RDyRegion](https://github.com/RDycore/RDycore/blob/main/include/private/rdycoreimpl.h#L20),
which is basically a named array of local cell IDs.

RDycore then uses information in the [YAML input file](../common/input.md) to
associate initial condition and source data with each region by name. This
process is implemented in the `InitRegions` function (https://github.com/RDycore/RDycore/blob/main/src/rdysetup.c#L182).


### Boundaries

Similarly to regions, the `DMPlex` object created during [RDySetup](https://github.com/RDycore/RDycore/blob/main/src/rdysetup.c)
constructs `DMLabel` objects represents disjoint sets of edges, each of which
represents a boundary. From each label/edge set, RDycore constructs an
[RDyBoundary](https://github.com/RDycore/RDycore/blob/main/include/private/rdycoreimpl.h#L29),
a named array of local edge IDs.

RDycore then uses information in the [YAML input file](../common/input.md) to
associate boundary condition data with each boundary by name. This process is
implemented in the `InitBoundaries` function, which is called in `RDySetup`.

## Operators

For purposes of our discussion, an "operator" is a mathematical construct for
computing interior and boundary fluxes and source terms. Each operator is
represented by a specific data structure associated with functions that
implement its operations.

_**NOTE**: at the time of writing, RDycore solves the shallow water equations
exclusively, so everything here refers to operators pertaining only to these
equations._

### Shallow Water Equations (SWE) operators
All shallow-water-equations-specific code lives within the [swe source subfolder](https://github.com/RDycore/RDycore/blob/main/src/swe/)
and is initialized by the `InitSWE` function, which is called in `RDySetup`.
The interface for the "SWE operator" is visible in the private header
[rdysweimpl.h](https://github.com/RDycore/RDycore/blob/main/include/private/rsdysweimpl.h).

All operator code has two parts:

* "host" code that runs on the CPU and dispatches calls to the "device" (a CPU
  or GPU, depending on your runtime configuration) to create, update,
  manipulate, and destroy the operator
* "device" code that performs the mathematical operators associated with the
  operator using the [libceed](https://ceed.exascaleproject.org/libceed/)
  exascale library.

The SWE operators are created and configured in the `InitSWE` function, which
lives in [`swe/swe.c`](https://github.com/RDycore/RDycore/blob/main/src/swe/swe.c).

The `libceed`-based "device" code for these operators is defined within [`swe/swe_ceed.c`](https://github.com/RDycore/RDycore/blob/main/src/swe/swe_ceed.c)
and the header file [`swe/swe_ceed_impl.h`](https://github.com/RDycore/RDycore/blob/main/src/swe/swe_ceed_impl.h).
The latter file defines inline [`CeedQFunction`](https://libceed.org/en/latest/api/CeedQFunction/#ceedqfunction)s
that run on the device.

We also implement a host-only version of the operators in [`swe/swe_petsc.c`](https://github.com/RDycore/RDycore/blob/main/src/swe/swe_petsc.c),
which we use mostly as a simple mechanism for troubleshoots our CEED
implementation.

We implement the following SWE operators:

* **Riemann Flux operator**: computes finite-volume fluxes between interior
  cells and on the boundaries, according to specified boundary conditions

* **Source operator**: computes the source term for the shallow water equations

In the CEED implementation, these two operators are combined into a
[composite operator](https://libceed.org/en/latest/api/CeedOperator/#c.CeedCompositeOperatorCreate)
that is called using in-place device data. In the PETSc version, these operators
are implemented by CPU code that is applied sequentially to relevant data.

## Input

The parsing of input is dispatched from [RDySetup](https://github.com/RDycore/RDycore/blob/main/src/rdysetup.c#L244)
and implemented in [ReadConfigFile](https://github.com/RDycore/RDycore/blob/main/src/yaml_input.c#L1173).
At the top of `yaml_input.c`, several structs are defined that represent a YAML
"schema". In this "schema," each struct represents a YAML mapping object, with
fields in the mapping corresponding to fields in the struct. Accordinglhy, YAML
arrays within mappings are represented by array fields.

For a more detailed explanation of how this YAML schema is parsed, see the
[libcyaml documentation](https://github.com/tlsa/libcyaml/blob/main/docs/guide.md).

## Output

RDycore can produce output for visualization, as well as diagnostic output
helpful for troubleshooting and debugging.

### Visualization

RDycore supports two visualizable output formats:

* [XDMF](https://www.xdmf.org/index.php/Main_Page): an ancient,
  marginally-supported format, popular in the earth science community, that
  stores data in HDF5 files and metadata in XML files
* [CFD General Notational System (CGNS)](https://cgns.github.io/): a standard
  format with lots of traction in the CFD community but less visibility within
  the earth science community
* PETSc's binary output format, which is better than nothing but probably not
  very full-featured.

The writing of XDMF files is implemented in [`xdmf_output.c`](https://github.com/RDycore/RDycore/blob/main/src/xdmf_output.c)
in the `WriteXDMFOutput` function. CGNS and binary output are handled mostly by
PETSc's built-in output machinery.

All output functions are registered using PETSc's
[TSMonitor](https://petsc.org/release/manualpages/TS/TSMonitor) feature and
called in `RDyAdvance` (within [rdyadvance.c](https://github.com/RDycore/RDycore/blob/main/src/rdyadvance.c))
at the frequency specified in the YAML input file.

### Diagnostics: time series

RDycore can write out time series data useful for debugging simulations.
Currently, the only quantity that can be recorded as a time series is the
"boundary flux"--the total flux through the domain boundary, which is needed to
account for mass non-conservation in simulations with open boundaries.

Time series data are accumulated and written by logic in [`time_series.c`](https://github.com/RDycore/RDycore/blob/main/src/time_series.c).
The function `InitTimeSeries` initializes this subsystem, and is called in
`RDyAdvance`.

These time series diagnostics are invasive, and can negatively affect simulation
performance, particularly when using GPUs. Therefore they are intended only for
debugging and troubleshooting.

## Checkpoints and Restarts

The writing of _checkpoint files_ (which store all state data necessary to
restart a simulation) and the process of restarting a previously-running
simulation are largely handled by PETSc:

* `InitCheckpoints` (defined in [checkpoint.c](https://github.com/RDycore/RDycore/blob/main/src/checkpoint.c)
  sets up a [TSMonitor](https://petsc.org/release/manualpages/TS/TSMonitor) that
  calls the `WriteCheckpoint` function at an interval specified in the YAML
  input file. This function is called in `RDySetup`.
* `ReadCheckpointFile` (also defined in `checkpoint.c`) reads the checkpoint
  file (in a format specified within the YAML input. This function is called
  in `RDySetup` if the input file specifieѕ that a simulation should be
  restarted.

## Running Ensembles

RDycore can run _ensembles_, which are sets of concurrently-running simulations
("ensemble members") with variations in selected input parameters. RDycore uses
an extremely simple mechanism to run these simulations within a single MPI job:

* on initialization, RDycore splits its "world" communicator into a number of
  smaller, communicators of equal size
* each of these smaller communicators is assigned to a single ensemble member
* each ensemble member is configured according to the ensemble specification
  in the [YAML input file](../common/input.md)
* all ensemble members are run concurrently to completion

Currently, ensemble calculations only run ensemble members, possibly generating
member-specific output. No post-processing or analysis has been implemented,
though it would likely be simple to do so.

### Ensemble member configuration

The data for each ensemble member is listed explicitly in the [ensemble section](../common/input.md#ensemble)
of the YAML input specification. All parameters for ensemble members override
the corresponding parameters in other sections of the file. This logic is
implemented in the `ConfigEnsembleMember` function within [ensemble.c](https://github.com/RDycore/RDycore/blob/main/src/ensemble.c).

## Support for Fortran

A Fortran version of the public C interface is available in an [`rdycore` Fortran module](https://github.com/RDycore/RDycore/blob/main/src/f90-mod/rdycore.F90)
in the `src/f90-mod` directory. This module is hand-written and uses the
[Fortran 2003 `iso_c_binding`](https://fortranwiki.org/fortran/show/iso_c_binding)
standard intrinsic module to map C functions to equivalent Fortran subroutines
and functions, and defines appropriate data types.

The mapping from a C function to a Fortran subroutine (or function) is
accomplished in two parts:

1. A C function is made available to Fortran by defining a Fortran function
   within an `interface` block that uses a `bind(C)` annotation specifying the
   case-sensitive name of the C function. All data types in this Fortran
   function must correspond to supported C data types via the `iso_c_binding`
   module. The function returns an integer corresponding to the `PetscErrorCode`
   return type of the C function. We refer to such a function as a "C-bound
   Fortran function".

2. A Fortran "wrapper" subroutine is defined that calls the C-bound Fortran
   function defined in item 1. This subroutine follows the PETSc convention in
   which the last argument (`ierr`) gets the return value of the C-bound Fortran
   function.

For example, let's take a look at the Fortran function for `RDyCreate`, which
calls the C function `RDyCreateF90`. This separate C function is required
because `RDyCreate` accepts an `MPI_Comm` input parameter, and Fortran uses
integers for MPI communicators.

First, we create a C-bound Fortran function `rdycreate_` and associate it with
the C function `RDyCreateF90` within the `interface` block near the top of
[`f90-mod/rdycore.F90`](https://github.com/RDycore/RDycore/blob/main/src/f90-mod/rdycore.F90):

```
  interface
    ...
    integer(c_int) function rdycreate_(comm, filename, rdy) bind(c, name="RDyCreateF90")
      use iso_c_binding, only: c_int, c_ptr
      integer,            intent(in)  :: comm
      type(c_ptr), value, intent(in)  :: filename
      type(c_ptr),        intent(out) :: rdy
    end function
    ...
  end interface
```

Then we define a subroutine `RDyCreate` that calls `rdycreate_` and sets `ierr`
to its return value:

```
  subroutine RDyCreate(comm, filename, rdy_, ierr)
    use iso_c_binding, only: c_null_char
    character(len=1024), intent(in) :: filename
    integer,   intent(in)  :: comm
    type(RDy), intent(out) :: rdy_
    integer,   intent(out) :: ierr

    integer                      :: n
    character(len=1024), pointer :: config_file

    n = len_trim(filename)
    allocate(config_file)
    config_file(1:n) = filename(1:n)
    config_file(n+1:n+1) = c_null_char
    ierr = rdycreate_(comm, c_loc(config_file), rdy_%c_rdy)
    deallocate(config_file)
  end subroutine
```

Notice that we have to do some things to construct a NULL-terminated C string
for the filename with a character array and `c_null_char`, and then pass a
pointer to this array with `c_loc`. This is how things are done with Fortran's
`iso_c_binding` module, which you can learn about in the link at the top of this
section.

Whenever a new C function is added to the public interface in `rdycore.h`, these
two corresponding Fortran items must be added to the [template from which it is generated](https://github.com/RDycore/RDycore/blob/main/include/rdycore.h.in)
to support its use in Fortran.

### Special considerations

* **C pointers**: Perhaps counterintuitively, C pointers must be passed by value
  to Fortran-bound C functions. In other words, any argument of type `type(c_ptr)`
  must have the `value` attribute and `intent(in)`. You can think of a pointer
  as a memory location, which is morally equivalent to an integer of the
  appropriate size. The `intent(in)` expresses that the pointer itself remains
  unchanged even if the data it points to is modified by the function.
  **NOTE**: an `RDy` C object is expressed as a C pointer in the Fortran module.
* **C primitives**: Because C passes arguments to functions by value and Fortran
  by reference, it is necessary to add the `value` attribute to any parameter
  in a C-bound Fortran function that has a non-pointer C type. This includes
  most `intent(in)` primitive parameters in C-bound Fortran functions. Note
  that `intent(out)` parameters in these functions must necessarily be pointers
  in C functions, so they must not have the `value` attribute.
* **Enumerated types**: An enumerated type in C can be mapped to Fortran as a
  set of related integer parameters. See, for example, the way [time units](https://github.com/RDycore/RDycore/blob/main/src/f90-mod/rdycore.F90#L35)
  are expressed in the `rdycore` Fortran module.
* **PETSc types**: PETSc types like `Vec` are passed to C-bound Fortran
  functions as `PetscFortranAddr` with the `value` attribute and `intent(in)`.
  This is very similar to the way we treat C pointers, with some magic PETSc
  pixie dust sprinkled on it to satisfy conditions required by PETSc.
