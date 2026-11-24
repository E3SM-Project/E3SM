(omega-dev-analysis)=

# Analysis

The Analysis module provides in-situ computation of analysis (diagnostic)
fields from the ocean model state during a simulation. Analysis fields are
computed on the fly at user-specified intervals and written to output
streams, avoiding the need to store large volumes of raw model output for
offline post-processing. User-level configuration is described in the
[User Guide](#omega-user-analysis) and the algorithmic formulation is
described in the [Design Document](#omega-design-analysis). This document
describes the implementation: the classes, their interfaces, and how to
extend the module.

The module is built around small, composable *operators*. Each operator
performs a single, well-defined transformation of one or more input
[Fields](#omega-dev-field) into one or more output Fields. Operators are
combined into *chains* (e.g. a spatial mean followed by a time mean), and the
set of all chains forms a directed acyclic graph in which nodes are operators
and edges are data dependencies. Bundled *analysis groups* (such as
`GlobalStats`) read the configuration file and construct the appropriate
chains and output [IOStreams](#omega-dev-iostreams) automatically.

In the initial implementation, analysis is produced by pre-defined bundled
groups. Support for fully user-defined, composable operator chains specified
directly in the configuration file is planned for a future version.

```{note}
This module is under active development, and a number of the implementation
details described below are temporary and expected to be refactored in future
updates. In particular, the name-based dependency resolution, the
insertion-order operator storage, the hard-coded group dispatch in the
`Analysis` constructor, and the currently limited set of supported array
ranks are all first-version choices that may change. Such areas are flagged
throughout this document. The public interfaces (`init`, `create`, `get`,
`computeAll`, `finalize`) and the operator/group extension points are
expected to remain stable.
```

## Architecture

The module lives in `src/analysis/` and is organized as follows:

- `Analysis.{h,cpp}` — the top-level `Analysis` orchestrator class.
- `AnalysisOperator.{h,cpp}` — the `AnalysisOperator` abstract base class and
  the `opParam`/`makeOpConfig` configuration helpers.
- `AnalysisOpFactory.{h,cpp}` — the `AnalysisOpFactory` used to create the
  correct templated operator specialization at run time.
- `AnalysisGroup.{h,cpp}` — the `AnalysisGroup` abstract base class for
  bundled groups.
- `operators/` — the built-in operator classes, plus `Ops.h` (a convenience
  header that includes every operator) and `Ops.cpp` (registers all
  operators with the factory).
- `analysisGroups/` — the built-in group classes, plus `Groups.h` (a
  convenience header that includes every group).

The core classes are:

- **`AnalysisOperator`**: abstract base class for all operators. Concrete
  operators are class templates parameterized on the Kokkos array type of
  their primary input Field.
- **`AnalysisOpFactory`**: a static factory that maps operator type names to
  constructor functions and dispatches to the correct templated
  specialization based on the runtime metadata of the input Field.
- **`OperatorNode`**: an internal structure that wraps an operator instance
  together with its upstream dependencies, output stream names, and compute
  alarms. It represents a single node in the dependency graph.
- **`AnalysisGroup`**: abstract base class for bundled groups that parse
  configuration, build operator chains, and create output streams. Includes
  helper structures (`OpChainInfo`, `StreamInfo`, `StreamParams`) for managing
  chain metadata and stream configuration parameters.
- **`Analysis`**: the top-level orchestrator that owns the operator graph,
  resolves dependencies, sets up alarm-based scheduling, and drives
  computation each timestep.

## Initialization and lifecycle

The Analysis module depends on the [`HorzMesh`](#omega-dev-horz-mesh),
[`VertCoord`](#omega-dev-vert-coord), and [`TimeStepper`](#omega-dev-time-stepping)
classes, so these must be initialized before Analysis. The default instance
is created with the static `init` method:
```c++
OMEGA::Analysis::init();
```
This registers all built-in operators with the factory, retrieves the
default mesh, vertical coordinate, machine environment, model clock, and
configuration, and then constructs the default `Analysis` instance. A pointer
to the default instance can be retrieved at any time with:
```c++
OMEGA::Analysis *DefAnalysis = OMEGA::Analysis::getDefault();
```

A non-default, named instance can be created from a name, machine
environment, mesh, vertical coordinate, model clock, and configuration:
```c++
OMEGA::Analysis *NewAnalysis = OMEGA::Analysis::create(Name, Env, Mesh,
                                                       VCoord, ModelClock,
                                                       Options);
```
All instances are tracked in a static map `AllAnalysisObjects`, and `create`
returns `nullptr` if an instance with the given name already exists.

Once per timestep, the main driver loop calls:
```c++
DefAnalysis->computeAll();
```
which computes every analysis field whose compute alarm is ringing at the
current time (see [Alarm-based scheduling](#alarm-based-scheduling) below).
The associated output streams are written separately through the normal
[IOStream](#omega-dev-iostreams) `writeAll` mechanism.

At the end of a simulation all instances are removed with:
```c++
OMEGA::Analysis::finalize();
```

## The AnalysisOperator base class

Every operator derives from `AnalysisOperator`. Because operators are
templated on the array type of their primary input, the base class provides a
common, non-templated interface used by the orchestrator. The key public
methods are:
```c++
const std::string getOperatorType();              // e.g. "SpatialMean"
const std::string getName();                       // unique instance name
const std::vector<std::string> getInputFieldNames();
const std::vector<std::string> getOutputFieldNames();
bool isCacheValid(const TimeInstant &TimeStamp);
virtual void initialize(const MachEnv *Env, const HorzMesh *Mesh,
                        const VertCoord *VCoord, Config Options);
virtual void setPeriodAlarm(Alarm *Alarm);
virtual void compute(const TimeInstant &TimeStamp) = 0;
```

The base class stores protected members shared by all operators, including
pointers to the mesh (`Mesh`) and vertical coordinate (`VCoord`), the MPI
communicator (`Comm`), the operator type name, the unique instance name, the
input and output Field name lists, and the caching state (`LastComputed`
timestamp and `FieldComputed` flag).

A concrete derived operator is responsible for three things:

1. **Construction**: pass the operator type name to the base constructor,
   record the input Field names, construct the output Field name(s), allocate
   the output data array, and register the output [Field](#omega-dev-field)
   with the appropriate metadata. Output Field names follow the convention
   `<input>_<OperatorSuffix>` (for example, `Temperature_SpatialMax`). The
   output's CF metadata is derived from the input's, as described under
   [Output metadata](#omega-dev-analysis-output-metadata).
2. **Initialization**: `initialize()` is called by the orchestrator after all
   Fields exist. The base implementation stores the mesh, vertical
   coordinate, and communicator; operators may override to perform additional
   setup.
3. **Computation**: `compute()` retrieves the typed input array(s) from the
   Field registry by name, performs the transformation, writes the result to
   the operator-owned output array, and updates `LastComputed`/`FieldComputed`
   for caching. Temporal operators additionally override `setPeriodAlarm()`
   so they know which alarm triggers finalization.

All derived operator classes must implement a constructor with the signature:
```c++
OperatorName(const std::vector<std::string> &UpstreamNames, Config Options)
```
This standardized interface enables the factory to instantiate operators
uniformly, passing upstream field names and configuration parameters
regardless of operator type.

The `isCacheValid()` method returns `true` when the operator has already been
computed for a given timestamp. This prevents redundant work when several
downstream operators share the same upstream intermediate result.

### Configuration helpers

Operators receive parameters through a `Config` object. The `Config` class
provides a generic, type-safe key-value container that can hold parameters of
different types (strings, integers, floats, etc.) without requiring explicit
knowledge of the parameter set at compile time, making operator constructors
flexible and extensible. To make it easy to build a `Config` inline, whether
an operator is created from user configuration or programmatically by a
bundled group, `AnalysisOperator.h` provides the `opParam` and `makeOpConfig`
helpers:
```c++
Config Opts = makeOpConfig(opParam("Period", std::string("1Day")));
```
`opParam(Key, Value)` builds a key/value pair and `makeOpConfig(...)`
assembles those pairs into a `Config` via variadic template expansion.

## Built-in operators

The following operators are currently implemented. All operators are class
templates parameterized on the input array type. The output Field name is the
input Field name with the operator suffix appended. The class name for each
operator is the operator type name with `Op` appended (e.g.
`SpatialMaxOp<ArrayType>`, `TimeMeanOp<ArrayType>`). Operators compute over all
active layers of owned mesh entities (cells, edges, or vertices) of the input
field(s).

| Operator | Inputs | Outputs | Output type | Output suffix | Config parameters | Description |
| -------- | ------ | ------- | ----------- | ------------- | ----------------- | ----------- |
| `SpatialMin` | 1 | 1 | scalar with same type as input (`Array1D<InType>`, dimension `Scalar`) | `_SpatialMin` | — | Global minimum of the input field. |
| `SpatialMax` | 1 | 1 | scalar with same type as input (`Array1D<InType>`, dimension `Scalar`) | `_SpatialMax` | — | Global maximum of the input field. |
| `SpatialMean` | 1 | 1 | scalar (`Array1DReal`, dimension `Scalar`) | `_SpatialMean` | — | Global weighted mean of the input field: by area for a horizontal field, by mass for a layered one (see [Spatial weights](#omega-dev-analysis-spatial-weights)). |
| `SpatialStdDev` | 2 (the field and its `_SpatialMean`) | 1 | scalar (`Array1DReal`, dimension `Scalar`) | `_SpatialStdDev` | — | Global weighted standard deviation of the input field, with the same weights as `SpatialMean`. Requires the field's `SpatialMean` as an upstream input, which is added to its input list automatically. |
| `TimeMean` | 1 | 1 | same rank and dimensions as the input (`Real`) | `_TimeMean<period>` | `Period` (string, e.g. `"1Day"`); `SamplingInterval` (string, optional) | Time average of the input field over a configurable period (e.g. `1Day`). By default accumulates every time step; the optional `SamplingInterval` instead gates accumulation to a coarser cadence (see below). Finalizes the mean when the period alarm rings. Output name embeds the period, e.g. `_TimeMean1Day`. |
| `BinaryMultiply` | 2 | 1 | same as first input | `_BinaryMultiply(<field2>)` | — | Element-wise multiplication of two fields. Supports vertical expansion (1D field replicated across vertical layers when multiplied with 2D/3D field). |
| `HorzMean` | 1 | 1 | 1D Real (`Array1DReal`) over the input's vertical dimension | `_HorzMean` | `WeightField` (string, optional, default `AreaCell`) | Area-weighted horizontal mean per vertical layer/interface. |
| `BinnedAccumulator` | 2 (data field, bin index field) | 1 | replaces horizontal dimension with `NumBins` | `_BinnedAccumulator(<binfield>)` | `NumBins` (I4) | Accumulates field values into spatial bins. Uses MPI reduction for global totals. Automatically applies regional mask if attached to input Field. |
| `CoordinateBinning` | 1 | 1 | 1D integer (`Array1DI4`) | `_BinIndex` | `NumBins` (I4), `MinBin` (Real, optional), `MaxBin` (Real, optional) | Assigns mesh entities to bins based on coordinate values (e.g., latitude). Computed once during initialization and cached. |
| `ExtractRegion` | 1 | 1 | same as input | `_<RegionName>` | `MaskName` (string) | Applies regional mask to field by multiplying values with mask. Attaches mask to output Field for downstream operators. Supports mask intersection. |
| `PrefixSum` | 1 | 1 | same as input | `_PrefixSum` | `Dimension` (I4), `Reverse` (bool, default false) | Cumulative summation along specified dimension. Forward (start→end) or reverse (end→start) scan. |
| `PseudoToGeometric` | 1 | 1 | same as input | `_Geometric` | — | Converts pseudo-height coordinates to geometric coordinates using specific volume. Handles vertical grid staggering automatically. |
| `ScalarMultiply` | 1 | 1 | same as input | `_ScalarMultiply(<value>)` | `Scalar` (string, parsed to Real) | Multiplies field by scalar constant. Can be used for unit conversions (e.g., ×1e-6 for Sverdrups). |
| `TransectAccumulator` | 3 (data field, mask field, sign field) | 1 | 1D vertical profile | `_TransectAccumulator(<transect>)` | `TransectName` (string) | Accumulates transport across transect edges. Only supports 2D Real arrays for the data field. Uses MPI reduction for global totals. |

(omega-dev-analysis-sampling-gate)=
### Sampling gate for temporal reductions

Temporal reduction operators (currently `TimeMeanOp`) inherit a reusable
sampling gate from the `AnalysisOperator` base class. The base class provides:

- `TimeInterval SamplingInterval`, `bool HasSamplingInterval`, and
  `Alarm *SampleAlarm` members (sampling disabled by default).
- `Alarm *getSampleAlarm()`, which returns `SampleAlarm` (null when sampling
  every timestep).
- `bool shouldSample(const TimeInstant &)`: returns `true` on every call when
  no sampling interval is configured (`HasSamplingInterval == false`), so the
  default path is bit-for-bit with the prior every-step behavior; otherwise it
  returns `true` only on steps where `SampleAlarm` is ringing and resets the
  alarm to schedule the next sample.

`TimeMeanOp` wraps its accumulate block in `if (shouldSample(TimeStamp))`; the
period-finalization logic runs **outside** this gate (every call) so a period
boundary that lands on a non-sample step is still honored. The operator reads
the `SamplingInterval` string from its own `Config` and, in `setPeriodAlarm()`
(where the reduction period and model timestep are both available), normalizes
the three accepted forms to a single `TimeInterval`:

- `"<n>perPeriod"` → `PeriodInterval / n`, divided in seconds for a period
  in days and in whole months for a period in months or years, since
  `TimeInterval` division truncates calendar intervals
- `"<n>steps"` → `n * ModelTimeStep`
- otherwise parsed directly as a `TimeInterval` string (e.g. `"6hour"`)

It then enforces the single constraint
`PeriodInterval.isDivisibleBy(SamplingInterval)` (`ABORT_ERROR` on failure) and
builds an operator-owned recurring `SampleAlarm` **anchored to the absolute
clock start** (`Clock::getStartTime()`), attaching it to the model clock; the
operator's destructor detaches it with `Clock::detachAlarm()`. The
alarm is then seeded to the current time via `SampleAlarm->reset(getCurrentTime())`,
which walks its boundaries forward to the first interval after "now". On a
restart, `Clock::setCurrentTime()` resets it again, with the clock's other
periodic alarms, once the restart time has been read. Because
the boundary sequence `ClockStart + n*SamplingInterval` is a deterministic walk
from a fixed origin, it is reconstructed identically on fresh starts and
restarts, making arbitrary `SamplingInterval` values ERS-safe.

No new restart state is needed: the sample schedule is reconstructed purely from
the absolute clock start, `SamplingInterval`, and the current time, so a restart
lands on the same boundary a straight-through run would hold.

`Analysis::setComputeAlarms()` uses the sample alarm, in place of the
every-step accumulation alarm, as the operator's accumulation alarm. Since
compute alarms propagate upstream, the operators that feed a sampled reduction
(for example, the MOC binning or a spatial mean) are also computed only on
sample steps, as with MPAS-Ocean's `compute_interval`.

(omega-dev-analysis-spatial-weights)=
### Spatial weights

`SpatialMean` and `SpatialStdDev` weight each mesh entity by how much of the
ocean it represents, so that the statistic does not depend on the mesh
resolution. The weights are built by the `SpatialWeights` class in
`analysis/SpatialWeights.h` from the input field's dimension names:

- A horizontal field (rank 1) is weighted by area: `AreaCell` for cells,
  `DcEdge * DvEdge` for edges and `AreaTriangle` for vertices.
- A layered field (last dimension `NVertLayers`) is weighted by mass, the
  area times the current `PseudoThickness` of the state, which is the layer's
  renormalized mass per unit area. The pseudo-thickness at an edge is the mean
  of its two cells and at a vertex the kite-area-weighted mean of its cells,
  as in the auxiliary variables.
- An interface field (last dimension `NVertLayersP1`) is weighted by half the
  mass of each adjacent layer, so the weights sum to the column mass and the
  result is the mass-weighted mean of the field taken as piecewise linear
  between interfaces. The surface and bottom interfaces have one adjacent
  layer each.

Only active entries take part. The weights are built from the active masks
of the vertical coordinate by selection rather than multiplication, and the
reductions skip inactive entries, so the fill values in inactive layers are
never read. Any leading extra dimensions of the field (e.g. tracers) share
the same weights.

An operator uses the weights by calling `init()` with the input field name
in `initialize()`, `update()` in `compute()` (mass weights follow the state,
so they are recomputed each time), and then `weightedSum()` or
`weightedSquaredDeviation()` with the input array and `weightSum()` for the
normalization. The sums accumulate each owned entity's terms in a fixed
order and reduce them with `globalSum`, so they are as reproducible as
Omega's other global sums.

(omega-dev-analysis-output-metadata)=
### Output metadata

An output Field inherits its CF metadata from the input it is derived from.
The base class provides the helpers that do this:

- `inheritMetadata(InputName, CellMethod)` returns an `InheritedMetadata`
  struct holding the input's `units`, `standard_name` and `cell_methods`,
  each empty if the input has none, with `CellMethod` appended to the cell
  methods. CF requires `cell_methods` to list the reductions in the order
  they were applied, so a time mean of a spatial mean carries
  `area: depth: mean time: mean`.
- `spatialCellMethod(InputName, Method)` returns the cell method for a
  reduction over all owned mesh entities and layers of the input:
  `area: <Method>` for a horizontal (1D) field and `area: depth: <Method>`
  for a field with a vertical dimension. `Method` is a CF cell method name
  (`mean`, `minimum`, `maximum`, `standard_deviation`, `sum`).
- `createOutputField(...)` calls `Field::create()` with the inherited units
  and standard name and adds the `cell_methods` attribute when it is not
  empty.
- `getFieldAttribute(FieldName, AttName)` returns any string attribute of a
  Field, or an empty string if the Field has none.

The built-in reductions (mean, minimum, maximum, standard deviation, time
mean) all have the units and standard name of the field they reduce, so
they pass `inheritMetadata()` straight to `createOutputField()`:
```c++
auto CellMethod  = spatialCellMethod(InputNames[0], "mean");
auto Meta        = inheritMetadata(InputNames[0], CellMethod);
auto OutputField = createOutputField(OutputNames[0],
                                     "Spatial mean of " + InputNames[0],
                                     Meta, ValidMin, ValidMax, NDims, DimNames);
```
An operator whose result has different units derives them with
[CFUnits](#omega-dev-cfunits) before creating the output, and clears the
standard name if the result is no longer the same physical quantity. For
example, an area-weighted sum multiplies the input's units by those of the
cell area:
```c++
auto Meta  = inheritMetadata(InputNames[0],
                             spatialCellMethod(InputNames[0], "sum"));
Meta.Units = CFUnits::multiply(Meta.Units, "m2");
Meta.StdName.clear();
```
The other operators derive their units the same way:

- `BinaryMultiply` multiplies the units of its two inputs (a velocity in
  `m s-1` times an area in `m2` is `m3 s-1`) and has no standard name or cell
  methods, since the product is a new quantity.
- `ScalarMultiply` scales the units by the reciprocal of its scalar with
  `CFUnits::scale()`, so a flux in `m3 s-1` multiplied by `1e-6` is in
  `1e+06 m3 s-1` (Sverdrups). A negative scalar also clears the standard
  name, since it reverses the quantity.
- `HorzMean` inherits the units and standard name and appends `area: mean`.
- `BinnedAccumulator`, `TransectAccumulator` and `PrefixSum` sum their input,
  and `ExtractRegion` and `PseudoToGeometric` multiply it by a mask or a
  dimensionless factor, so all of them inherit its units. `PrefixSum` aborts
  if its boundary condition field has different units from its input.
- `CoordinateBinning` outputs a bin index, which is dimensionless (`1`).

A field without units yields outputs without units, never with invented ones,
and a malformed units string aborts the run when the operator is constructed.

## Operator factory and type dispatch

Fields carry runtime metadata describing their scalar type, rank, and memory
location, but operators are compile-time templates. `AnalysisOpFactory`
bridges this gap. It maintains a static registry (a Meyer's singleton) that
maps a fully-qualified key to a constructor function. Keys have the form
`<BaseName>_<ArrayType>_<MemLoc>`, for example `SpatialMax_Array2DR8_Device`.

All variants of an operator template are registered with a single call:
```c++
AnalysisOpFactory::registerAllArrayVariants<SpatialMaxOp>("SpatialMax");
```
This expands the `OMEGA_ANALYSIS_ARRAY_TYPES` macro over the supported
combinations of scalar type (`I4`, `I8`, `R4`, `R8`), rank, and memory
location, registering a creator lambda for each. The currently enabled ranks
are 1 through 3 (the 4D and 5D entries are present but commented out because
Omega currently utilizes no arrays with rank $>$ 3, and each registered
operator variant increases compile time and binary size).

For operators that only support specific array ranks, the factory provides
specialized registration methods using SFINAE (Substitution Failure Is Not An
Error) to prevent invalid template instantiations at compile time:

- `register1DVariants` — registers only 1D array variants (horizontal fields without vertical structure)
- `register2DVariants` — registers only 2D array variants (horizontal × vertical)
- `register2DRealVariants` — registers only 2D Real array variants (excludes integer types)

These methods use compile-time type checking to filter out incompatible array
types, reducing compile time and binary size while catching type errors early.
The SFINAE helpers (`register1DVariantHelper`, etc.) expand over all array types
but only register the matching variants, with no-op branches for mismatched types.

At creation time, the orchestrator calls:
```c++
auto Op = AnalysisOpFactory::createOp(OpType, UpstreamNames, Options);
```
`createOp` inspects the metadata of the primary upstream Field
(`UpstreamNames[0]`), builds the fully-qualified key, looks up the matching
creator, and invokes it. If no matching variant is registered, it aborts.

All built-in operators are registered in one place,
`Analysis::registerAllBaseAnalysisOperators()` (defined in
`operators/Ops.cpp`), which is called during `Analysis::init()`.

## Operator chains, dependencies, and caching

A chain is expressed as an underscore-delimited string beginning with a
source Field name followed by operator tokens, for example
`Temperature_SpatialMean_TimeMean1Day`. The orchestrator turns a chain string
into operator instances with:
```c++
DefAnalysis->parseChainAndBuildOps(ChainStr);
```
This walks the chain left to right, building the progressively longer output
name at each stage. If a Field with that name already exists, the stage is
skipped and reused, which is how common intermediate results (such as a
shared `SpatialMean`) are naturally deduplicated across chains. Otherwise the
operator kind is identified from the token (a `Spatial` prefix denotes a
spatial reduction; a `Time` prefix denotes a temporal reduction, with the
period parsed from the first digit in the token) and the operator is created
via the factory and wrapped in an `OperatorNode`.

**Parenthesized arguments:** The parser supports operators with parenthesized
arguments for passing field names, scalar values, or other parameters. The
parser respects parenthesis depth when splitting on underscores, so
`Field_BinaryMultiply(OtherField)_ScalarMultiply(1.0e-6)` correctly parses as
three tokens: `Field`, `BinaryMultiply(OtherField)`, and `ScalarMultiply(1.0e-6)`.
Different operators interpret parenthesized arguments differently:
- `BinaryMultiply(FieldName)` — second input field name
- `ScalarMultiply(Value)` — scalar value (parsed from string to Real)
- `ExtractRegion(MaskName)` — mask field name passed via Config
- `BinnedAccumulator(BinIndexField)` — bin index field name
- `TransectAccumulator(TransectName)` — transect name for mask lookup

An optional `Config` parameter can be passed to `parseChainAndBuildOps()` to
provide operator-specific configuration that applies to all operators in the
chain (e.g., binning parameters, integration direction).

After all groups have registered their chains, the constructor resolves the
graph edges with `buildOperatorDependencies()`, which matches each node's
input Field names against the output Field names of every other node and
records the matches as upstream pointers.

Each timestep, `computeAll()` checks every node's compute alarms and, for any
node whose alarm is ringing, calls `computeRecursive()`. This method first
checks the operator's cache (`isCacheValid()`) and returns immediately on a
hit; otherwise it recursively computes all upstream dependencies before
computing the node itself. Combined with the timestamp cache, this guarantees
that each operator runs at most once per timestep even when shared by
multiple chains.

```{note}
In this version, dependency resolution uses post-hoc name matching and the
operators are stored in insertion order rather than a formally sorted
topological order. The recursive, cache-guarded evaluation ensures upstream
operators are always computed before their downstream consumers. A
signature-based deduplication and explicit ordering are planned for a future
version.
```

## Alarm-based scheduling

Each `OperatorNode` carries a list of compute alarms in `ComputeAlarms`.
`setComputeAlarms()` assigns these as follows:

- **Terminal operators** (those that write to one or more output streams)
  borrow the write alarm from each associated stream via
  `IOStream::getAlarm()`.
- **Temporal reduction operators** additionally receive an accumulation alarm.
  With a `SamplingInterval`, it is the operator's sample alarm (see
  [Sampling gate](#omega-dev-analysis-sampling-gate)); otherwise it rings
  every timestep and is created and owned by the `Analysis` instance and
  attached to the model clock. The stream's write alarm serves as the period
  (output) alarm;
  it is passed to the operator with `setPeriodAlarm()` so the operator knows
  when to finalize and reset its accumulator. Thus a temporal operator has
  two alarms: one that triggers accumulation and one that triggers
  finalization at the end of the period.
- **Instantaneous output operators** receive only the stream's alarm, so they
  are evaluated only when a snapshot is written.

Finally, `propagateAlarmsUpstream()` performs a fixed-point iteration that
copies each node's alarms to all of its upstream dependencies. This ensures
that an intermediate operator is computed whenever any downstream consumer
needs fresh data.

## Analysis groups

An `AnalysisGroup` encapsulates the configuration parsing, chain
construction, and output stream creation for a named, bundled set of analysis
outputs. The base class provides helper structures and methods:

- `OpChainInfo` records, for a single chain, the chain string, its output
  frequency, and whether it performs temporal reduction or instantaneous
  output.
- `StreamInfo` aggregates multiple operator chains that share the same output
  frequency and type (temporal reduction vs. instantaneous) for grouping into
  a single IOStream. Contains operator names, frequency string, and reduction
  type flag.
- `StreamParams` holds default values for all IOStream options, allows them
  to be overridden from the group configuration via `apply()`, and converts
  them to a `Config` with `toConfig()`.
- `createAnalysisGroupStreams()` groups the chains by their output
  characteristics (period and reduction type), creates the associated
  IOStreams, and associates operator output Fields with those streams.

Concrete groups are dispatched by name in the `Analysis` constructor, which
reads the `Analysis` configuration node, iterates over its child group nodes,
and constructs each enabled group.

```{note}
The group dispatch is currently a hard-coded name check in the `Analysis`
constructor, and enabling a group that is not a recognized bundled group
aborts. This is a first-version choice; the dispatch mechanism and support
for user-defined composable groups are expected to be generalized in a future
update.
```

### Restart persistence for temporal reductions

A temporal reduction may span a model restart (for example, a monthly mean
when restarts are written every ten days). To keep the running average
correct across such a boundary, each `TimeMeanOp` checkpoints its mid-period
accumulation state in the restart file.

Two pieces of state are persisted per operator, both as data Fields added to
the existing `Restart` FieldGroup (the same mechanism other modules such as
`OceanState`, `Tracers`, and `VertCoord` already use):

1. The partial accumulated sum — the operator's existing output array is
   attached directly to a restart Field (no copy), since that array already
   holds the running sum mid-period.
2. `NumAccum` — a one-element `I4` data Field holding the sample count. A data
   Field (not metadata) is required because Field-level metadata is written
   but not read back on restart.

The restart Fields are created at the end of the `TimeMeanOp` constructor via
the operator-private, non-virtual `createRestartFields()`; there is no
base-class hook for it. Both Fields use `OptionalRead = true` for backward
compatibility, so restart files written before this feature read cleanly with
`NumAccum` left at `0` (a fresh period, i.e. the pre-existing behavior).

`IsNewPeriod` is not stored; it is reconstructed as `(NumAccum == 0)`. The
backing count array is kept in sync inline inside `TimeMeanOp::compute()`
(the only place `NumAccum` is mutated), so the checkpoint is consistent
regardless of when the restart write occurs. After a period finalizes the
count is normalized to `0` to mark that no partial period is in progress.

On restart, `Analysis::applyRestartState()` iterates the operator nodes and,
for each operator reporting `needsRestart()`, derives `IsNewPeriod` from the
restored `NumAccum`. It is a restart-only step, invoked in `OceanInit.cpp`
after the restart read and gated on `ReadRestart` for standalone runs and
`CoupledReadRestart` in coupled runs.

### GlobalStats

`GlobalStats` is the first concrete group. It computes global spatial
statistics of a configurable set of fields, optionally followed by temporal
reduction. Its configuration lists the fields, the spatial statistics, the
temporal reduction periods, and/or the instantaneous output intervals:
```yaml
GlobalStats:
  Enable: true
  Fields: [Temperature, Salinity]
  SpatialStats: [Mean, Max, Min]
  ReductionPeriod: [1day, 1month]
  SnapshotPeriod: [6hour]
```
For each (field, statistic) pair, the constructor builds a spatial reduction
chain (for example `Temperature_SpatialMean`). For each reduction period it
appends a `TimeMean` operator to produce a time-averaged chain (for example
`Temperature_SpatialMean_TimeMean1Day`); for each snapshot period it records
the discrete-sampling chain for instantaneous output. It then calls
`createAnalysisGroupStreams()` to create the corresponding output streams. At
least one of `ReductionPeriod` or `SnapshotPeriod` must be present.

### MOC

`MOC` computes the Meridional Overturning Circulation (MOC) streamfunction
using two methods: latitude-binned regional MOC and transect-based MOC.
Configuration specifies binning parameters, region/transect names, and output
frequencies:

```yaml
MOC:
  Enable: true
  NumBins: 180                  # Number of latitude bins
  MinLat: -90.0                 # Minimum latitude in degrees
  MaxLat: 90.0                  # Maximum latitude in degrees
  Regions: [Global, Atlantic]   # Regional MOC (latitude × depth)
  Transects: [Drake, Atlantic26N]  # Transect MOC (depth only)
  ReductionPeriod: [1Month]
  SnapshotPeriod: [1Day]
```

**Latitude-binned MOC** computes the streamfunction as a function of latitude
and depth for each specified region. The constructor builds a complex operator
chain for each region:
1. `LatCell_CoordinateBinning` — assigns cells to latitude bins (initialization only, shared across regions)
2. `VerticalPseudoVelocity_PseudoToGeometric` — converts to geometric coordinates
3. `VerticalVelocity_BinaryMultiply(AreaCell)` — computes vertical flux
4. `[Optional] VerticalFlux_ExtractRegion(RegionMask)` — applies regional mask (for non-global regions)
5. `VerticalFlux_BinnedAccumulator(LatCell_BinIndex)` — accumulates flux into latitude bins
6. `BinnedFlux_PrefixSum` — horizontal integration (south to north)
7. `MOC_ScalarMultiply(1.0e-6)` — converts to Sverdrups (Sv)

**Transect-based MOC** computes the streamfunction across specific transects
as a function of depth only. The constructor builds a chain for each transect:
1. `PseudoThickness_PseudoToGeometric` — converts to geometric layer thickness
2. `LayerThickness_BinaryMultiply(NormalVelocity)` — thickness × velocity
3. `EdgeTransport_BinaryMultiply(DvEdge)` — transport × edge width
4. `TransportField_TransectAccumulator(TransectName)` — accumulates across transect edges
5. `TransectTransport_PrefixSum` — vertical integration (bottom to top, reverse)
6. `TransectMOC_ScalarMultiply(1.0e-6)` — converts to Sverdrups

For each chain, the constructor stores an operator configuration (`Config`)
containing parameters like `NumBins`, `Dimension`, `Reverse`, and `TransectName`
in a `ChainConfigs` vector. These configs are passed to `buildTemporalChains()`,
which appends temporal operators and applies custom `IOName` metadata to
produce user-friendly output variable names (e.g., `MOC_streamfunction_Global`
instead of the full operator chain string). At least one of `Regions` or
`Transects` must be specified.

```{note}
Regional and transect masks are not yet implemented in the Omega infrastructure.
The `Regions` configuration currently only supports `[Global]`, and `Transects`
is a placeholder. The operator infrastructure (ExtractRegionOp,
TransectAccumulatorOp) is fully implemented and ready for use once mask fields
become available.
```

## Extensibility

The module is designed so that new operators and new groups plug in through
the registration and dispatch mechanisms without changes to the core compute
loop.

### Adding a new operator

To add a new operator:

1. **Create the operator class** in `operators/` as a template on the input
   array type, deriving from `AnalysisOperator`. In the constructor, pass the
   operator type name to the base constructor, set `InputNames`, build the
   output name as `<input>_<Suffix>`, set `OutputNames` and `InstanceName`,
   allocate the output data array, and register the output Field with
   `createOutputField()` plus `attachData()`, deriving its units, standard
   name and cell methods from the input as described under
   [Output metadata](#omega-dev-analysis-output-metadata). If the operator
   depends on another operator's output (as `SpatialStdDev` depends on
   `SpatialMean`), append that dependency to `InputNames`.
2. **Implement `compute()`** to retrieve the input array(s) from the Field
   registry, perform the transformation, write the output array(s), and update
   `LastComputed`/`FieldComputed`.
3. **Implement `setPeriodAlarm()`** only if the operator is a temporal
   reduction that finalizes on a period boundary.
4. **Include the header** in `operators/Ops.h`.
5. **Register the operator** by adding a
   `AnalysisOpFactory::registerAllArrayVariants<NewOp>("NewName")` call to
   `Analysis::registerAllBaseAnalysisOperators()` in `operators/Ops.cpp`.
6. **Update the chain parser** if your operator introduces a new naming
   convention. The current parser in `Analysis::parseChainAndBuildOps()`
   identifies operators by prefix (`Spatial` for spatial reductions, `Time`
   for temporal reductions). If your operator uses a different naming
   pattern, add logic to recognize and dispatch the new token. The parser is
   preliminary and expected to evolve as new operator categories are added.
7. **Add unit tests** in `test/analysis/AnalysisOperatorTest.cpp` to verify the
   new operator implementation. Tests should validate correct computation on
   representative input data.

If your operator follows an existing naming convention, no parser changes are
needed. No changes to the dependency resolution, alarm propagation, or compute
logic are required.

### Adding a new analysis group

To add a new bundled group:

1. **Create the group class** in `analysisGroups/`, deriving from
   `AnalysisGroup`. In the constructor, read the group's configuration
   options, build chain stem strings (field names with spatial operators but
   without temporal operators), and call `buildTemporalChains()` to handle
   temporal operator construction, operator instantiation via
   `parseChainAndBuildOps()`, and automatic `OpChainInfo` population.
   Alternatively, for groups with custom chain building logic, manually call
   `AnalysisManager->parseChainAndBuildOps()` for each chain and record an
   `OpChainInfo` for each chain describing its output frequency and whether it
   performs temporal reduction or instantaneous output.
2. **Create the output streams** by calling `createAnalysisGroupStreams()`
   at the end of the constructor.
3. **Include the header** in `analysisGroups/Groups.h`.
4. **Add a dispatch branch** so the `Analysis` constructor constructs the new
   group when it is enabled in the configuration.

The guiding principle, as with the rest of the module, is that extensions
should not require restructuring the core: new capabilities are added by
registering operators and dispatching groups, leaving the graph construction,
scheduling, and compute loop unchanged.

## Verification and testing

The Analysis module is covered by two standalone C++ test executables, both
registered in `test/CMakeLists.txt` and run with eight MPI ranks.
`AnalysisOperatorTest.cpp` verifies the accuracy of the operators, while
`AnalysisSystemTest.cpp` tests the individual components of the Analysis
framework.
