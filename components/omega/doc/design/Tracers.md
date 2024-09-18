(omega-design-tracers)=

# Tracer Infrastructure

## 1 Overview

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents). This document
describes a tracer module that manages the definition, description, storage,
and access of tracers throughout OMEGA. For reasons outlined in the following
requirements, specific methods or algorithms for transporting or mixing
tracers are not addressed here. These algorithms will be covered in other
design documents.

## 2 Requirements

### 2.1 Requirement: Tracer definition and metadata

The tracer module must provide a centralized location where all supported
tracers are defined. This definition includes all metadata associated with
the tracers—see related documents for metadata and I/O. Where possible,
CF-compliant metadata must be used.

### 2.2 Requirement: Tracer identification

The tracer module must provide a simple way to identify a tracer, either
by name or by index. The latter may be preferred for performance reasons.

### 2.3 Requirement: Tracer groups

It is often convenient to select or define properties for tracers as
a group. It must be possible to assign tracers to particular groups.
Properties of the group should apply to all tracers assigned to that group.
Tracers can belong to multiple groups, as long as properties of each group
do not conflict. A group called "All" will include all selected tracers.
For performance reasons, it may be desirable for some non-intersecting
groups to be assigned contiguous memory or contiguous index space.

### 2.4 Requirement: Tracer selection

At the start of a simulation, users must be able to select which of the
supported tracers are enabled for the experiment via the configuration
input file. To simplify the configuration, tracers should be selectable
by group as well as individually.

### 2.5 Requirement: Tracer restart and I/O

With few exceptions, tracers will generally be included in all restarts and
often in other outputs for post-processing and analysis. Thus, it is
essential for tracers to be registered as available I/O fields (see related
I/O designs) so that they can be read from or written to a file.

### 2.6 Requirement: Time levels

Time integration schemes often require the storage of multiple time levels.
Tracer storage must support multiple time levels, and a mechanism for
specifying which time level is needed for a given calculation must be
provided. Ideally, changing time levels at the end of a time step should
avoid a memory copy and instead involve a simple change of time index or
pointer.

### 2.7 Requirement: Acceleration or supercycling

Some tracers or tracer groups can be evolved using a longer time step. For
performance reasons, it should be possible to advance these tracers with a
longer time step and avoid computing tendencies at every dynamical time step.

### 2.8 Desired: Per-tracer/group algorithmic requirements

Different tracers may require different methods for their evolution to
preserve desired properties. For example, temperature and salinity should
maintain consistency with the continuity equation. Many tracers will require
monotonic or positive-definite advection schemes to prevent negative
concentrations. We may even use a different time-stepping scheme for
supercycled tracers.

## 3 Algorithmic Formulation

The algorithms used to advect, mix, and evolve tracers in time are described
elsewhere in the respective modules responsible for computing these tendencies.

## 4 Design

We will store all tracers in a vector of 3-dimensional device arrays.

```c++
  std::vector<Array3DReal> TracerArrays; // multiple time levels
```

Each device array in the vector has dimensions corresponding to the number
of tracers, the number of all cells in the rank, and the vertical length.
This allows parallelization over tracer and cell dimensions, with the vertical
dimension either parallel or vectorized, depending on the architecture and
computation involved. Many of the requirements for grouping, selecting, and
swapping time levels can be achieved through the collection and manipulation
of array indices. Individual tracers can still be extracted as arrays with
contiguous memory.

All tracers will be defined in a set of code blocks (described below) in
a file named `TracerDefs.inc`, which will be directly included in the
initialization routine to define all supported tracers and their metadata.
Tracer groups are defined as sets of tracer indices in the OMEGA YAML
configuration file.

### 4.1 Data Types and Parameters

#### 4.1.1 Definition

All supported tracers are defined in a file called `TracerDefs.inc`, which
includes the code blocks needed to define each tracer and its metadata. This
file must include all tracers potentially used in the model. The file will
look something like this:

```c++
// At the top of the file, list all available tracers in a comment
// to provide a sort of table of contents.
// Tracers:
//   Temp
//   Salt
//   [Others...]

define(
       "Temp",                  // Name of variable
       "Potential Temperature", // Long name or description
       "degree_C",              // Units
       "sea_water_potential_temperature", // CF standard name
       -273.15,                 // Min valid field value
       100.0,                   // Max valid field value
       1.e33                    // Fill value for undefined entries
);
```

#### 4.1.2 Configuration

Users will select desired tracers either individually or by group. If a group
is listed, it implies the entire group, but users can choose a subset of
a group by specifying individual tracer names and omitting the group name.

```yaml
omega:
  Tracers:
    Base: [Temp, Salt]
    Debug: [Debug1, Debug2, Debug3]
    [other individual tracers or groups as needed]
```

For requirement 2.7 (per-tracer algorithm choice), the configuration might
look something like this:

```yaml
omega:
  TracerAdv:
    # Start with default values for all tracers
    - All:
         Method: standard
         SomeParameter: 0.5
    # Specialized options - override default by group or tracer name
    - Ecosys:
         Method: MyMonotoneMethod
    - [Other tracers or tracer groups, only if needed]:
         Method: MyAdditionalOption
         SomeParameter: 1.0
```

#### 4.1.3 Classes/Structs/Data Types

There will be a single class that contains the static array for storing
all tracers, along with several containers for managing tracer groups and
indices. Methods are described in Section 4.2 below.

```c++
/// The Tracers class provides a container for tracer arrays and methods
class Tracers {

 private:

   // Total number of tracers defined at initialization
   static I4 NumTracers;

   // Static storage of the tracer arrays for device and host memory spaces
   // The arrays are 3-dimensional, with dimensions for Tracer, Cell, and
   // Vertical layers
   static std::vector<Array3DReal>
       TracerArrays; ///< Time levels -> [Tracer, Cell, Vertical]
   static std::vector<HostArray3DReal>
       TracerArraysH; ///< Time levels -> [Tracer, Cell, Vertical]

   // Maps for managing tracer groups
   // The key of this map is a group name and the value is a pair of
   // GroupStartIndex and GroupLength
   static std::map<std::string, std::pair<I4, I4>> TracerGroups;

   // Maps for matching tracer names with indices (in both directions)
   static std::map<std::string, I4> TracerIndexes;
   static std::map<I4, std::string> TracerNames;

   // Halo exchange
   static Halo *MeshHalo;

   // Tracer dimension names
   static std::vector<std::string> TracerDimNames;

   // Current time index
   // This index is circular, so it returns to index 0 when it exceeds the max
   // index
   static I4 CurTimeIndex; ///< Time dimension array index for the current level

   // Pack tracer field name
   static std::string packTracerFieldName(const std::string &TracerName);

   // Defines all tracers locally without allocating memory
   static I4
   define(const std::string &Name,        ///< [in] Name of the tracer
          const std::string &Description, ///< [in] Long name or description
          const std::string &Units,       ///< [in] Units
          const std::string &StdName,     ///< [in] CF standard name
          const Real ValidMin,            ///< [in] Min valid field value
          const Real ValidMax,            ///< [in] Max valid field value
          const Real FillValue            ///< [in] Value for undefined entries
   );

  public:
    [methods defined below]
};
```
Here is the corrected version with improved grammar and clarity:

---

### 4.2 Methods

#### 4.2.1 `init`

The initialization method defines all the tracers by including the
`TracerDefs.inc` file. Once all fields are defined, it allocates the full
tracer array based on the input configuration of selected tracers. The
interface is as follows:

```c++
static int Err = init();
```

Note that this function does not initialize the tracer values. The values
will be filled either by input I/O streams or by individual modules
associated with specific tracers and tracer groups.

#### 4.2.2 `define`

All tracers must be defined. This routine takes the name of the tracer and
some required metadata as input. It uses this information to create
a `MetaData` entry and an `IOField` entry, but only if the tracer is
selected in the input configuration file. It also increments the number of
tracers and assigns a tracer index if the tracer is selected. The function
returns a non-zero error code if any of several possible error conditions
are encountered.

```c++
static int define(
   const std::string &Name,        ///< [in] Name of the tracer
   const std::string &Description, ///< [in] Long name or description
   const std::string &Units,       ///< [in] Units
   const std::string &StdName,     ///< [in] CF standard name
   OMEGA::Real ValidMin,           ///< [in] Minimum valid field value
   OMEGA::Real ValidMax,           ///< [in] Maximum valid field value
   OMEGA::Real FillValue           ///< [in] Fill value for undefined entries
);
```

#### 4.2.3 `getNumTracers`

This function retrieves the total number of tracers activated for the current
simulation.

```c++
static I4 getNumTracers();
```

#### 4.2.4 `getIndex`

This function retrieves the tracer index based on the tracer name.

```c++
// Retrieves tracer index from tracer name
static I4 getIndex(I4 &TracerIndex,              ///< [out] Tracer index
                   const std::string &TracerName ///< [in] Tracer name

);
```

#### 4.2.5 `getAll`

This function retrieves the full array of all tracers at a specified time
level. The current time level is specified by an integer value of 0, the
previous time level is -1, and so on. If the requested time level does not
exist, the function returns a negative integer.

```c++
// get a device array for all tracers
static I4 getAll(Array3DReal &TracerArray, ///< [out] tracer device array
                 const I4 TimeLevel        ///< [in] time level index
);
```

#### 4.2.6 `getByIndex` (Single Tracer)

This function retrieves a single tracer array by index and time level.
Similar to `getAll`, the current time level is specified by an integer value
of 0, and the previous time level is -1. If the requested time level does
not exist, the function returns a negative integer.

```c++
// get a device array by tracer index
static I4 getByIndex(Array2DReal &TracerArray, ///< [out] tracer device array
                     const I4 TimeLevel,       ///< [in] time level index
                     const I4 TracerIndex      ///< [in] global tracer index
);
```

#### 4.2.7 `getGroupRange`

This function retrieves a pair of the group's start index and group length.

```c++
// Retrieves a pair of (group start index, group length)
static I4
getGroupRange(std::pair<I4, I4> &GroupRange, ///< [out] Group range
              const std::string &GroupName   ///< [in] Group name

);
```

A typical use of this function might look like:

```c++
// Get group range
std::pair<I4, I4> GroupRange;
I4 Err = Tracers::getGroupRange(GroupRange, GroupName);

// Unpack group range
auto [StartIndex, GroupLength] = GroupRange;

// Get all tracers at the current time level (0)
Array3DReal TracerArray;
Err = OMEGA::Tracers::getAll(TracerArray, 0);
if (Err != 0)
   LOG_ERROR("getAll returns an error code: {}.", Err);

OMEGA::parallelFor(
   "ComputeGroupTendency",
   {GroupLength, NCells, NVertLayers},
   KOKKOS_LAMBDA(Int iIndex, Int iCell, Int iVert) {
      int iTracer = TracerArray[iIndex + StartIndex];
      // Perform operations on TracerArray(iTracer, iCell, iVert);
   });
```

#### 4.2.8 `isGroupMemberByIndex`

When looping over tracers, this function returns `true` if the current
tracer index belongs to a group with a given name.

```c++
// Checks if a tracer is a member of a group by tracer index
static bool
isGroupMemberByIndex(const I4 TracerIndex,      ///< [in] Tracer index
                     const std::string GroupName ///< [in] Group name
);
```

#### 4.2.9 Managing Time Levels

The exact process for managing time integration is not yet clear. We
assume a time integration module will assign and swap time indices (e.g.,
CurTime, PrevTime) to point to the appropriate array index in all
prognostic arrays like tracers. The current time is represented by
an integer value of "0", while the previous time is represented as "-1",
and so on.

#### 4.2.10 `clear`

This function clears all defined tracers and cleans up memory. It returns
an error code if it is unable to deallocate or destroy all memory.

```c++
// Clears all defined tracers and deallocates memory.
static int clear();
```

Here is the corrected version with proper grammar:

---

## 5 Verification and Testing

### 5.1 Test All

A sample `TracerDefs.inc` file will be located in the source directory,
containing the base tracers and three tracers in a Debug group. We will
assume two time levels (Cur: 0, Prev: -1). A unit test driver will read an
input configuration that requests these two tracer groups. The tracers will
be initialized with artificial data, and basic arithmetic will be performed
on the device. In one loop, we will perform arithmetic on one group using
the vector of group tracer indices. In a second loop, we will perform
operations on the second group using the `isMembers` function to select tracers.

This test covers requirements 2.1–2.4 and 2.6.

All other requirements will need to be tested in later system tests once
I/O streams and time integration schemes have been implemented.
