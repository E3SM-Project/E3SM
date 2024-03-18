(omega-design-tracers)=
# Tracer infrastructure

## 1 Overview

Tracers are either heat (temperature) or material carried by a
fluid parcel (eg salt, chemical or biological constituents). Here
we describe a tracer module that manages how tracers are defined,
described, stored and accessed throughout OMEGA. For reasons
outlined in the requirements below, we do not address
the specific methods or algorithms for transporting or mixing
tracers here. Those algorithms will be covered in other
design documents.

## 2 Requirements

### 2.1 Requirement: Tracer definition and metadata

The tracer module must provide a single location in which all
supported tracers are defined. This definition includes all
metadata associated with the tracer - see related documents
for metadata and IO. Where possible, CF compliant metadata
must be used.

### 2.2 Requirement: Tracer identification

The tracer module must provide a simple way to identify a tracer,
either by name or by index. The latter may be preferred for
performance reasons.

### 2.3 Requirement: Tracer groups

It is often convenient to select or define properties of tracers
as a group. It must be possible to assign tracers to a particular
group or groups. Properties of the group should apply to all
tracers assigned to that group. Tracers can be members of multiple
groups as long as properties of each group do not conflict. There
will be a tracer group for All that will include all selected
tracers. For performance reasons, it may be desirable for some
non-intersecting groups to be assigned contiguous memory or with
a contiguous index space.

### 2.4 Requirement: Tracer selection

At the start of the simulation, users must be able to select
which of the supported tracers are enabled for the experiment
through the configuration input file. For a more concise
configuration, it should be possible to select tracers by group
as well as individually.

### 2.5 Requirement: Tracer restart and IO

With few exceptions, tracers will generally be part of all
restarts and will often be included in other output for
postprocessing and analysis. Practically, this is a
requirement for tracers to be registered as available IO
fields (see related IO designs) so that they can be read
from or written to a file.

### 2.6 Requirement: Time levels

Time integration schemes often require the storage of multiple
time levels. Tracer storage must support multiple time levels and
a mechanism for specifying which time level is needed for a
given calculation. Ideally the change of time levels at the end
of a time step should avoid a memory copy and be a simple change
of time index or change of pointer.

### 2.7 Requirement: Acceleration or supercycling

Some tracers or tracer groups can be evolved using a longer
time step. For performance, it should be possible to advance
these tracers with a longer time step and avoid computing
tendencies at every dynamical time step.

### 2.7 Desired: Per-tracer/group algorithmic requirements

Different tracers may require different methods for their evolution
to preserve desired properties. For example, temperature and
salinity should preserve tracer consistency with the continuity
equation. Many tracers will require monotone or positive-definite
advection schemes to prevent negative concentrations. We might
even use a different time stepping scheme for supercycled tracers.


## 3 Algorithmic Formulation

Algorithms used to advect, mix and evolve tracers in time are
described elsewhere in the respective modules that compute these
tendencies.

## 4 Design

We will store all tracers in a single static array on the device:
```c++
  Array4DReal AllTracers(MaxTimeLvls, NumTracers, NCellsSize,
                         VertSize);
```
This will allow us to parallelize over tracer and cell dimensions
with the vertical dimension either parallel or vectorized depending
on architecture and computation involved. Many of the requirements
for grouping, selecting and swapping time levels can be achieved
through collection and manipulation of array indices. Individual
tracers can still be extracted as arrays with contiguous memory.

All tracers will be defined in a set of code blocks (described
below) in a file `TracerDefs.inc` that will be directly included
in the init routine to define all supported tracers and their
metadata. Tracer groups are defined as a set of tracer indices
that belong to the group. A set of primary groups will be defined
in `TracerDefs.inc` that define membership for the purposes of
selecting tracers. Additional groups may be added by other
modules as needed (eg if we decide to do per-tracer algorithm
choices, the tracer advection module may define groups that
assign tracers to a particular method).

### 4.1 Data types and parameters

#### 4.1.1 Definition

All supported tracers are defined in a file called `TracerDefs.inc`
with the code blocks needed to define the tracer and its metadata.
This must include all tracers potentially used in the model.
The file will look something like:
```c++
// At the top of the file, list all available tracers in a
// comment to provide a sort of table of contents.
// Tracers:
//   Temp
//   Salt
//   [Others...]
// Tracer Groups:
//   Base (or Active?)
//   Debug
//   Ecosys
//   Analysis

// Individual tracer entries will look like:
// Tracer: Temp

   Err = OMEGA::Tracer::define(
          "Temp",                  // Name of var
          "Potential Temperature", // long Name or description
          "degree_C",              // units
          "sea_water_potential_temperature", /// CF standard Name
          -273.15,                 // min valid field value
          100.0,                   // max valid field value
          1.e33                    // fill value for undef entries
   );
   Err = OMEGA::Tracer::addToGroup("Base", "Temp");
```

#### 4.1.2 Configuration

Users will select desired tracers either individually or by
group. If the group is in the list, it implies the entire
group, but the user can choose a subset of a group by
supplying the individual tracer names and omitting the group name.
```yaml
omega:
  Tracers:
     - Base (or Active or ? - basically the required T,S)
     - Debug
     - Ecosys
     - IdealAge
     - [other individuals or groups as needed]

```

For requirement 2.7 (per-tracer algorithm choice), the configuration
could look something like:
```yaml
omega:
  TracerAdv:
    # start with default values for all tracers
    - All
         Method: standard
         SomeParameter: 0.5
    # specialized options - override default by group or tracer name
    - Ecosys:
         Method: MyMonotoneMethod
    - [Other tracers or tracer groups, only if needed]
         Method: My additional option
         SomeParameter: 1.0
```


#### 4.1.3 Class/structs/data types

There is a single class that contains the static array for storing
all tracers and a number of containers to manage tracer groups and
indices. Methods are described in section 4.2 below.

```c++
class Tracers {

  private:
    // Total number of tracers activated for this simulation
    static int NumTracers

    // static storage of the tracer array for all tracers
    static Array4DReal AllTracers;

    // maps for matching tracer names with indices (both directions)
    static std::map<std::string, int> Index
    static std::map<int, std::string> Name

    // maps for managing tracer groups
    /// This map associates a group name with a list of tracer
    /// indices that are assigned to the group
    static std::map<std::string, std::vector<int>> Groups

    /// This map matches the group name with a boolean vector
    /// of NumTracers that flags a tracer index as a member
    static std::map<std::string, std::vector<int>> MemberFlag;

  public:
    [methods defined below]
}
```

### 4.2 Methods

#### 4.2.1 init

The initialization method is where all the tracers are defined
by including the `TracerDefs.inc` file. Once all the fields have
been defined, it allocates the full tracer array based on the
input configuration of selected tracers. The interface is simply:
```c++
static int Err = init();
```
Note that this function does not actually initialize the
tracer values. The values will be filled either by input IO streams
or by individual modules associated with some tracers and tracer
groups.

#### 4.2.2 define

All tracers must be defined. This routine takes as input the
name of the tracer and some required metadata. It uses these
inputs to create a MetaData entry and an IOField entry, but only
if the tracer is selected in the input configuration file. It also
increments the number of tracers and assigns a tracer index if the
tracer is selected. It returns an error code that is non-zero if
any of a number of error conditions are encountered.
```c++
static int define(
   const std::string &Name,        ///< [in] Name of tracer
   const std::string &Description, ///< [in] Long name or description
   const std::string &Units,       ///< [in] Units
   const std::string &StdName,     ///< [in] CF standard Name
   OMEGA::Real ValidMin,           ///< [in] min valid field value
   OMEGA::Real ValidMax,           ///< [in] max valid field value
   OMEGA::Real FillValue           ///< [in] value for undef entries
);
```

#### 4.2.3 addToGroup

Tracers can be added to a group with this method:
```c++
static int addToGroup(
   const std::string &GroupName, ///< [in] name of group
   const std::string &TracerName ///< [in] name of tracer to add
);
```
and this routine is typically called right after the tracer is
defined and included as part of the `TracerDefs.inc` file. Tracers
can be added to more than one group so multiple calls are possible.
If the group does not yet exist, one will be created with the name
provided. If the tracer already belongs to the group, no action is
taken. At this point, not all information on tracers and groups are
available since other fields in a group may not have been defined
or added yet. Once all fields have been defined and added to groups
(i.e. at the end of the tracers init routine), the `pruneGroups`
function is called to ensure the groups and group members are
consistent with the tracers selected in the input config file.

#### 4.2.4 organize

During the definition phase, groups are defined that may not
have been fully selected and some information cannot be inferred.
A final pass of the tracers and groups will be necessary to sort
tracer indices for more optimal access (eg contiguous indices),
prune the groups, and finalize group members. This should only be
called at the end of the `Tracers::init` so is defined as a
private function:
```c++
private:
  static int organize();
```
An error code is returned if any error conditions are encountered.

#### 4.2.5 getNumTracers

Retrieves the total number of selected tracers activated for the
current simulation.
```c++
static int getNumTracers();
```

#### 4.2.6 getIndex

This function retrieves the tracer index given a tracer name.
```c++
static int getIndex(const std::string &TracerName);
```

#### 4.2.7 getAll

This function retrieves the full array of all tracers at a
given time level.
```c++
static Array3DReal getAll(int TimeIndx ///< [in] time level index
);
```

#### 4.2.8 get (single tracer and time level)

To retrieve a single tracer array by index and time level, we
provide the function:
```c++
static Array2dReal get(int TrcrIndx, ///< [in] global tracer index
                       int TimeIndx  ///< [in] time level index
);
```

#### 4.2.9 getGroup

This function retrieves a vector of tracer indices that belong
to a group. The group name must be provided.
```c++
static std::vector<int> getGroup(const std::string &GroupName);
```
A typical use of this function might look like:
```c++
   Array4DReal TrcrArray = OMEGA::Tracers::getAll();
   std::vector<int> EcoTracerIndices =
                    OMEGA::Tracers::getGroup("Ecosys")

   int NGrpMembers = EcoTracerIndices.size();
   OMEGA::parallelFor(
        "ComputeGroupTendency",
        {NGrpMembers, NCells, NVertLayers},
        KOKKOS_LAMBDA(Int iTracer, Int iCell, Int k) {
           int TrcrIndx = EcoTracerIndices[iTracer];
           do stuff with TrcrArray(TrcrIndx, iCell, k);
    });

```

#### 4.2.10 isMember

When looping over tracers, this function returns true if the
current tracer index belongs to a group of a given name.
```c++
static bool isMember(
         int Index,   ///< [in] global tracer index
         const std::string &GroupName);
```
This provides a complement to the getGroup function so the
loop in the previous section could also look like:
```c++
   int NTracers = OMEGA::Tracers::getNumTracers;
   Array4DReal TrcrArray = OMEGA::Tracers::getAll();

   OMEGA::parallelFor(
        "ComputeGroupTendency",
        {NTracers, NCells, NVertLayers},
        KOKKOS_LAMBDA(Int iTracer, Int iCell, Int k) {
           if (isMember(iTracer,"GrpName"){
              do stuff with TrcrArray(TrcrIndx, iCell, k);
           }
    });
```

#### 4.2.11 Managing time levels

It is not clear yet how time integration will be managed. We assume
here that a time integration module will be assigning and swapping
time indices (eg CurTime, NewTime, OldTime?, etc.) to point to
appropriate array index in all prognostic arrays like tracers.
If these need to be managed within the tracer module, we may need
a function like swapTimeLevels to identify and modify time level
indices.

#### 4.2.12 clear

A function is needed to clear all defined tracers and clean up
memory. An error code is returned if it is unable to deallocate
or destroy all memory.
```c++
static int clear();
```

## 5 Verification and Testing

### 5.1 Test All

A sample TracerDefs.inc file will be created in the test directory
with the base tracers (T,S) and two tracers in a Debug group.
We will assume two time levels (Cur, New). A unit test driver
will read an input configuration that requests these two
tracer groups. The tracers will be initialized with artificial
data and basic arithmetic will be performed on the device. In one
loop nest, we will perform arithmetic on one group using the
vector of group tracer indices. In a second loop nest, we will
perform operations on the second group using the isMembers function
to select.
  - tests requirements 2.1-2.4, 2.6

All other requirements will need to be tested in later system
tests once IOStreams and time integration schemes have been
implemented.
