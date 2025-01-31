(omega-design-decomp)=

# Decomp


## 1 Overview

For parallel execution, OMEGA will utilize a hierarchical parallelism with the
coarse-grained parallelism achieved via message passing. An OMEGA mesh will be
decomposed into subdomains so that computations on each sub domain can proceed
in parallel and messages will be passed to communicate data between domains.
This document describes the requirements and design of this domain
decomposition.


## 2 Requirements

### 2.1 Requirement: Mesh partition

The primary requirement is to be able to partition a mesh across a
distributed memory machine. Many standard tools exist to perform
this partitioning (eg metis, parmetis, zoltan) and their use is
encouraged. This normally requires an input file with the primary
mesh indices and adjacency (see Requirement 2.4).

### 2.2 Requirement: Metis option

For backward compatibility and to facilitate bit-for-bit comparisons
with MPAS, at least one partition option must use the Metis library.

### 2.3 Required: Online partitioning

We require the partitioning to occur online (ie as part of model
initialization rather than a separate pre-processing step).
This is meant to avoid the past difficulties in maintaining a large
number of partition input files for each MPI configuration. Experience
has shown that even high-resolution partitioning is very inexpensive
and can be performed as part of the model run and does not need this
off-line step.

### 2.4 Required: Input mesh information

As input into the partitioning, we require a file with information
on cell adjacency. There are at least two potential design approaches
for this.  The first option is to use the global `graph.info` file
which is already in Metis input format with a line for each cell
listing the global indices of that cell's neighbors. A second option
is to use the input mesh netCDF file which has the same information
in the global cellsOnCell array. This second option would enable us to
eliminate all graph.info files after the mesh file has been created.

### 2.5 Required: All mesh locations

While the primary partition is based on the cells and cell
adjacency, the decomposition must also partition the edge and vertex
index spaces.

### 2.6 Required: Halo

For parallel efficiency, a halo containing copies of nearest neighbor
points will be defined. A separate design document describes the
capability needed for updating those halo points for field data.
Within the decomp type, we simply keep the halo index information.

### 2.7 Required: Local to global index mapping

Once the mesh is decomposed, each parallel process has a set of
local indices renumbered for the local space. The decomposition
must provide the ability to return the global index given the
processor's local index. It must also be able to
provide the location (MPI rank and local index) for all the
neighbor or halo points to facilitate later setup of halos and
other infrastructure. Some of this information may only be needed
for initialization (see Requirement 2.8 on release of memory).

### 2.8 Required: Local sorting

At a minimum, the local addresses must be sorted with owned
indices first and each level of halo following. This is needed
for both optimal communication and for selectively computing
only what is needed. Additional sorting within those blocks
may be desireable for optimizing memory locality.

### 2.9 Required: Release memory

Once the partition information has been used to set up halos,
mesh variables, I/O and other quantities that depend on the
partition, all or parts of the partitioning data may be
released to free memory. Note that some of the information
may be transferred to mesh variables as part of later
mesh initialization.

### 2.10 Desired: Weighted partitioning option

For better load-balancing, support for supplying or computing weights
for each cell index based on estimated workload is desireable.

### 2.11 Desired: Multiple partitions of same mesh

In the future, it may be desireable to perform parts of the calculation
(eg the communication-dominated barotropic mode) on a smaller
partition. We may need to support multiple partitions of the
same mesh on different numbers of MPI ranks.

### 2.12 Desired: Multiple domains, sub-blocking or coloring

For some future capabilities (eg local timestepping in regionally
focused mesh locations), it may be desireable to decompose the
global domain into smaller subdomains and over-subscribe subdomains
to a rank, combining low- and high-resolution subdomains for load
balancing. It is unclear whether this is best achieved with a multiple
sub-block structure (as in older MPAS) or accomodating this within a
weighted partitioning. Identifying the mesh regions using a coloring
of the mesh may be helpful.

## 3 Algorithmic Formulation

Most of the algorithms required are embedded and documented in the
relevant packages (eg Metis multi-level k-way partitioning).

## 4 Design

The design presented here is initially for decomposing the global
domain into a single block per MPI rank and does not yet include
weighted decompositions. However, it should not prevent these
options from being added later. The parmetis package will be
used initially to provide metis back compatibility while reducing
memory use over a serial metis approach. The cell mesh indices and
adjacency will be read from the non-decomposed graph.info files in
metis format, allowing the adjacency info to be read serially
before the parallel IO is set up.

Once the primary (cell) mesh is partitioned, the edge and vertex
index spaces will be assigned based on the adjacency to the local
cells, as in the current MPAS model. Because much of this information
will be later replicated in a mesh class with Kokkos arrays, this
decomposition structure will be destroyed after initialization and
any related setup of halos and other infrastructure.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

The decomposition includeds a public enum to define the
supported partitioning method. Initially, this will support
the default Metis method, but others can be added:

```c++
    enum partMethod {
       partMethodUnknown,        ///< undefined for error checking
       partMethodMetisKWay,      ///< default KWay method in metis
       }
```

Other methods, like the Metis GeomKWay (combination of space-filling
curve with KWay) or Zoltan options may be added later.

Eventually, an enum defining options for computing weights for
weighted partitions may also be required, but will not appear in
this initial implementation.

Another parameter is the `haloWidth` that defines the width
of the halo needed to provide neighbor information for all
local operators. The minimum is typically 3 for current
algorithms.

These parameters and the mesh input file will be read as part of
the input configuration file in a decomp configuration group:

```yaml
    decomp:
       meshInputFilename: 'omegaMeshFile.nc'
       partitionMethod: 'MetisKWay'
       haloWidth: 3
```

The Metis KWay partition method will be default and the haloWidth
will initially default to 3.

#### 4.1.2 Class/structs/data types

There will be a Decomp class that stores and defines a partitioning
of the address space.

```c++

    class Decomp {

       private:
          partMethod method;   ///< method to use for partitioning


       public:

          I4 haloDepth;    ///< depth of the halo

          I4 nCellsGlobal; ///< num cells in the global domain
          I4 nCellsAll;    ///< num cells on local rank (own+halos)
          I4 nCellsOwned;  ///< num cells exclusively owned by local rank

          std::vector<I4>[] nCellsHalo[]; ///< num of owned+halo points
                                          ///< for each level of halo

          std::vector<I4>[] cellGlobalID[]; ///< global index for each
                                            ///< local cell (own+halo)

          /// The adjacency (neighbor) info will be stored in the
          /// compact form that Metis uses. In this form the neighbor
          /// cell location for each owned and ghost cell in local cell
          /// order is stored as a linear vector of integer pairs
          /// (partition/processor number and local address). A second
          /// vector stores the starting address in this list for the
          /// neighbors associated with the local cell address.

          std::vector<I4>[] nbrStartCell ///< start addr in nbr array for cell
          std::vector<I4>[] nbrLocCell   ///< list of nbr addresses

          [ repeat identical variables for edge and vertex index spaces ]


          // methods described below

    }
```

### 4.2 Methods

Because most class members are public, only constructors and
destructors are really needed.

#### 4.2.1 Constructor

The main method will be a constructor that decomposes the mesh and
creates the public variables. The default decomposition will require
no arguments and will read options from the input config file:

```c++
    Decomp mainDecomp;
```

For multiple decomposition of the same mesh, a second constructor
will be needed that starts from the default decomposition and
uses arguments for determining decomposition options. The exact form
will be determined later, but might look something like:

```c++
    Decomp newDecomp(oldDecomp, nParts, method, haloWidth);
```

#### 4.2.2 Destructor

After most of the initialization is complete and the decomposition
information has been used to initialize meshes, halo, I/O, etc.,
we will supply a destructor to release the space in memory.

```c++
    delete myDecomp;
```

## 5 Verification and Testing

### 5.1 Test Metis

We will create a small sample graph/mesh for partitioning,
similar to examples in Metis documention and create a partition.
The resulting partition will be compared to that generated
by the standalone gpmetis partition.
  * tests requirement 2.1, 2.2, 2.3, 2.4

### 5.2 Back compatibility

We will create a decomposition based on an MPAS input mesh
and verify that the decomposition is the same as the older
MPAS Fortran code. This test should only be needed after initial
development and will not need to be repeated as part of
routine testing

### 5.3 Check global IDs

To make sure all cells, etc. are accounted for, we can do a global
sum of the number of cells as well as a sum of global IDs across
the partition and ensure they the expected sum.
