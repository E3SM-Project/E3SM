(omega-design-halo)=
# Halo

## 1 Overview

In a distributed memory parallel environment, it is necessary to regularly
exchange data on the interfaces between adjacent partitions of the mesh
(i.e. halo elements). These halo exchanges will be accomplished utilizing
the MPI library.

## 2 Requirements

### 2.1 Requirement: Compatible with all types of arrays

Handle halo exchanges for arrays of each supported data type and dimensionality.

### 2.2 Requirement: Exchange at all types of mesh elements

Meshes are composed of three types of elements (cells, edges, and vertices)
at which values can be defined, we must be able to handle exchanges at each
type of element.

### 2.3 Requirement: Variable halo depth

The default minimum required halo depth for most of the algorithms we will
use is three cells. We should also be able to handle different halo depths.

### 2.4 Requirement: Exchange lists

Each MPI rank needs a send list and receive list for every neighboring rank and
each type of mesh element (cells, edges, and vertices) containing the local
indices for each owned element to be packed and sent to neighboring ranks and
each element owned by neighboring ranks to be received and unpacked into local
field arrays. Exchange lists should be organized by halo layer.

### 2.5 Requirement: Send and receive buffers

Each MPI rank needs two arrays for every neighboring rank to serve as a send
buffer and a receive buffer to be passed to the MPI send and receive routines.
On each rank, one buffer per neighbor will be sent concurrently, so there needs
to be one send and one receive buffer array for each neighbor to facilitate
communication.

### 2.6 Requirement: Non-blocking communication

To efficiently handle exchanges between each rank and all of its neighbors, we
will use the non-blocking MPI routines `MPI_Isend` and `MPI_Irecv`. The
`MPI_Test` routine will be used to determine when each communication is
completed. Boolean variables will be needed to track when messages have been
received and when buffers have been unpacked.

### 2.7 Requirement: Compatible with device-resident arrays

In environments with GPU devices, many arrays will reside solely on the device
and we will need to be able to manage halo exchanges for these arrays, either
via GPU-aware MPI or handling array transfers to host.

### 2.8 Desired: Communicate subset of halo layers

We may want the ability to communicate individual halo layers as opposed to
always exchanging the entire halo for a given array. Exchange lists should be
organized by halo layer to facilitate this.

### 2.9 Desired: Exchange multiple arrays simultaneously

Communications in MPAS are generally latency limited, so it may be advantageous
to exchange multiple arrays at the same time, packing the halo elements of
multiple arrays into the same buffer to reduce the number of messages.

### 2.10 Desired: Multiple environments/decompositions

OMEGA runs may contain multiple communication environments or mesh
decompositions, would need to be able to perform halo exchanges in these
circumstances.

### 2.11 Desired: OpenMP threading

If we allow OpenMP threading, halo exchanges would need to be carried out in a
thread-safe manner, and a process to carry out local exchanges may be needed.

### 2.12 Desired: Minimize buffer memory allocation/deallocation

It might be advantageous to have persistent buffer arrays that are large enough
to handle the largest halo exchange in a particular run, as opposed to
calculating the needed buffer size and reallocating the buffer arrays for each
halo exchange. Only a subset of the buffer array would need to be communicated
for smaller exchanges. A method for determining the largest necessary buffer
size for each neighbor would be needed, and would depend on which arrays are
active, the number of halo elements in each index space, the number of vertical
layers and tracers.

## 3 Algorithmic Formulation

No specific algorithms needed outside of those provided by standard MPI library

## 4 Design

The mesh decomposition will be performed during initialization via the Decomp
class based on the input mesh file and options defined in the configuration
file. A set of class objects described below will take the Decomp object and
store information about the decomposition as well as allocate buffer memory
necessary for performing halo exchanges.

### 4.1 Data types and parameters

#### 4.1.1 Parameters

An enum defining mesh element types would be helpful to control which exchange
lists to use for a particular field array:

```c++
enum meshElement{onCell, onEdge, onVertex};
```

An MPI datatype ```MPI_RealKind``` will be defined based on whether the default
real data type is single or double precision via the SINGLE_PRECISION compile
time switch:

```c++
#ifdef SINGLE_PRECISION
MPI_Datatype MPI_RealKind = MPI_FLOAT;
#else
MPI_Datatype MPI_RealKind = MPI_DOUBLE;
#endif
```

Halo exchanges will depend on the ```haloWidth``` parameter defined in the
decomp configuration group.

#### 4.1.2 Class/structs/data types

The ExchList class will hold both send and receive lists. Lists are separated by
halo layer to allow for the possibility of exchanging individual layers. Buffer
offsets are needed to pack/unpack all the halo layers into/from the same buffer.

```c++
class ExchList {

   private:

      I4 nList[haloWidth];        ///< number of mesh elements in each
                                  ///<   layer of list
      I4 nTot;                    ///< number of elements summed over layers
      I4 offsets[haloWidth];      ///< offsets for each halo layer

      std::vector<I4> indices[haloWidth];   ///< list of local indices

   friend class Neighbor;
   friend class Halo;
};
```

The Neighbor class contains all the information and buffer memory needed for
carrying out each type of halo exchange with one neighboring rank. It contains
the ID of the neighboring MPI rank, an object of the ExchList class for sends
and receives for each mesh index space, a send and a receive buffer, MPI request
handles that are returned by `MPI_Irecv` and `MPI_Isend` which are
needed to test for completion of an individual message transfer, and boolean
switches to control the progress of a Halo exchange. To avoid the need for
separate integer and floating-point buffers, during the packing process integer
values can be recast as reals in a bit-preserving manner using
`reinterpret_cast` and then recast as integers on the receiving end.

```c++
class Neighbor {

   private:

      I4 rankID;                           ///< ID of neighboring MPI rank
      ExchList sendLists[3], recvLists[3]; ///< 0 = onCell, 1 = onEdge,
                                           ///< 2 = onVertex
      std::vector<Real> sendBuffer, recvBuffer;
      MPI_Request rReq, sReq;              ///< MPI request handles
      bool received = false;
      bool unpacked = false;

   friend class Halo;
};
```

The Halo class collects all Neighbor objects needed by an MPI rank to perform
a full halo exchange with each of its neighbors. The total number of neighbors,
the local rank, and the MPI communicator are saved here. The mesh element type
for the current array being transferred is also stored here. This class will be
the user interface for halo exchanges, it is a friend class to the subordinate
classes so that it has access to all private data.

```c++
class Halo {

   private:

      I4 nNghbr;                           ///< number of neighboring ranks
      I4 myRank                            ///< local MPI rank ID
      MPI_Comm myComm;                     ///< MPI communicator handle
      meshElement elemType;                ///< index space of current array
      std::vector<Neighbor> neighbors;

   public:

      // methods

};
```

### 4.2 Methods

#### 4.2.1 Constructors

The constructor for the Halo class will be the interface for declaring an
instance of the Halo class and instances of all associated member classes,
and requires info from the Decomp and MachEnv objects.

```c++
Halo(Decomp inDecomp, MachEnv inEnv);
```

The constructors for Neighbor and ExchList will be called from within the Halo
constructor and will create instances of these classes based on info from the
Decomp object.

#### 4.2.2 Array halo exchange

The primary use for the Halo class will be a member function called
`exchangeFullArrayHalo` which will exchange halo elements for the input
array with each of its neighbors across all layers of the halo. The index space
the array is located in (cells, edges, or vertices) needs to be fetched from the
metadata associated with the array stored in the Halo object to determine which
exchange lists to use. It will return an integer error code to catch errors
(likewise for each subprocess below).

```c++
int exchangeFullArrayHalo(ArrayLocDDTT &array, meshElement elemType);
```

This will be an interface for different exchange funcitons for each type of
ArrayLocDDTT (where Loc is the array location, i.e. device or host, DD is the
dimensionality and TT is the data type) as defined in the DataTypes doc. If the
array is located on a GPU device, this function may tranfer arrays from device
to host. The ordering of steps for completing a halo exchange is:
1. start receives
2. pack buffers
3. start sends
4. unpack buffers

#### 4.2.3 Start receive/send

A `startReceives` function will loop over all member Neighbor objects and
call `MPI_Irecv` for each neighbor. It takes no arguments because all the
info needed is already contained in the Halo object and its member objects.

```c++
int StartReceives();
```

Likewise, `startSends` will loop over all Neighbor objects and call
`MPI_Isend` to send the buffers to each neighbor.

```c++
int StartSends();
```

#### 4.2.4 Buffer pack/unpack

For each type of ArrayDDTT, there will be a buffer pack function aliased to a
`packBuffer` interface:

```c++
int packBuffer(ArrayDDTT array, int iNeighbor);
```

where halo elements of the potentially multidimensional array ArrayDDTT are
packed into the 1D send buffer for a particular Neighbor in the member
`std::vector` neighbors using the associated ExchList based on the
meshElement the array is defined on. The `exchangeFullArrayHalo` will loop
over the member neighbors and call `packBuffer` for each Neighbor.

Similarly, buffer unpack functions for each ArrayDDTT type will be aliased to
an `unpackBuffer` interface:

```c++
int unpackBuffer(ArrayDDTT &array, int iNeighbor);
```

## 5 Verification and Testing

### 5.1 Halo exchange tests

A unit test where a set of halo exchanges are performed for each type of
ArrayDDTT and each index space could verify all requirements, except for 2.3.
The test should use a mesh decomposition across multiple MPI ranks, and define
an array for each type of ArrayDDTT and index space. The array elements could be
filled based on the global IDs for the cells, edges, or vertices. After an
exchange, all the elements in an array (owned+halo elements) would be checked
to ensure they have the expected value to verify a successful test. A similar
test utilizing a mesh decomposition with HaloWidth other than the default value
would satisfy requirement 2.3 as well.
