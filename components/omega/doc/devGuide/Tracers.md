(omega-dev-tracers)=

# Tracers

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents).

To manage tracer definitions, users will update the `TracerDefs.inc` file
located in the `omega/src/ocn` directory of E3SM. The file contains
`define` C++ function calls, each of which defines a tracer as shown below:

```c++
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

To add a new tracer, simply call the `define` function with the appropriate
arguments.

The following sections explain the concept and implementation of OMEGA tracers.

## Concepts

### Arrays in Host and Device Memory

All tracers are stored in a vector of 3-dimensional device and host arrays.

```c++
  std::vector<Array3DReal> TracerArrays;   // Device: multiple time levels
  std::vector<HostArray3DReal> TracerArraysH; // Host: multiple time levels
```

Each device and host array in the vector has dimensions corresponding to the
number of tracers, the number of cells in the rank, and the vertical levels.

The host tracer array (`TracerArraysH`) is internally managed, so users
should not directly modify it in most cases.

### Tracer Groups

Each tracer is assigned to a tracer group defined in the OMEGA YAML
configuration file, such as `Default.yml`. A tracer should be assigned to only one group.

To access the member tracers of a group in the code, users can use the
`getGroupRange` function as shown below:

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

### Time Levels

During the initialization of Tracers, the number of time levels is determined,
and the length of the tracer vectors (`TracerArrays` and `TracerArraysH`)
is set accordingly. The Tracers class internally manages the index of the
current time level in the arrays. To access a specific time level in the
arrays, users will use the integer "0" or negative integers. For example,
the following code retrieves the device tracer arrays for the current and
previous time levels:

```c++
// Get all tracers at the current time level (0)
Array3DReal CurrentTracerArray;
I4 Err1 = OMEGA::Tracers::getAll(CurrentTracerArray, 0);

// Get all tracers at the previous time level (-1)
Array3DReal PreviousTracerArray;
I4 Err2 = OMEGA::Tracers::getAll(PreviousTracerArray, -1);
```

## Initialization and Finalization

To use Tracers, the class should first be initialized by calling the
`init()` function and finalized by calling `finalize()`. These function
calls are typically handled within the initialization and finalization of
OMEGA itself, so users may not need to call them separately.

## Key APIs

### `getAll` and `getAllHost`

These functions return all device and host tracer arrays, respectively. If
the specified `TimeLevel` does not exist, they return a negative integer.

```c++
static HostArray3DReal getAllHost(
   HostArray3DReal &TracerArrayH, ///< [out] tracer host array
   const I4 TimeLevel ///< [in] Time level index
);

static Array3DReal getAll(
   Array3DReal &TracerArray, ///< [out] tracer device array
   const I4 TimeLevel ///< [in] Time level index
);
```

### `getGroupRange`

`getGroupRange` returns a pair of two `I4` values representing the group
start index and the group length. If the group is not found, it returns
negative integers.

```c++
static I4 getGroupRange(
   std::pair<I4, I4> &GroupRange, ///< [out] Group range
   const std::string &GroupName   ///< [in] Group name
);
```

### `getByIndex` and `getHostByIndex`

These functions return device and host tracer arrays, respectively, based
on the `TracerIndex`. If the specified `TimeLevel` and/or `TracerIndex`
does not exist, they return a negative integer.

```c++
static Array2DReal getByIndex(
   Array2DReal &TracerArray, ///< [out] tracer device array
   const I4 TimeLevel,  ///< [in] Time level index
   const I4 TracerIndex ///< [in] Global tracer index
);

static HostArray2DReal getHostByIndex(
   HostArray2DReal &TracerArrayH, ///< [out] tracer host array
   const I4 TimeLevel,  ///< [in] Time level index
   const I4 TracerIndex ///< [in] Global tracer index
);
```

### `getIndex`

`getIndex` returns the index of the tracer specified by the `TracerName`
argument. If the tracer is not found, it returns a negative integer.

```c++
static I4 getIndex(
   I4 &TracerIndex,              ///< [out] Tracer index
   const std::string &TracerName ///< [in] Tracer name
);
```

### `getFieldByIndex`

`getFieldByIndex` returns the `Field` object associated with the tracer
specified by the `TracerIndex`. If not found, it returns `nullptr`.

```c++
// Returns a field by tracer index. If it does not exist, returns nullptr
static std::shared_ptr<Field>
getFieldByIndex(const I4 TracerIndex ///< [in] Global tracer index
);
```

### `updateTimeLevels`

`updateTimeLevels` increments the current time level by one. If the current
time level exceeds the number of time levels, it returns to `0`. It also
exchanges the halo cells of all tracers at the current time level and updates
the associated field with the new tracer arrays.

```c++
/// Increment time levels
static I4 updateTimeLevels();
```

### `getNumTracers`

`getNumTracers` returns the total number of tracers used in the simulation.

```c++
// Get total number of tracers
static I4 getNumTracers();
```
