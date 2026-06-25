<!--- OMEGA Analysis Module design document ---------------------------------->

(omega-design-analysis)=

# Analysis

## 1 Overview

The Omega Analysis module provides in-situ computation of desired analysis
fields from the ocean model state. Analysis fields are computed on-the-fly
during simulation runtime and written to output streams at user-specified
intervals, providing an alternative to extensive offline post-processing.

The framework is built on a composable operator architecture where operators
can be chained together to produce analysis outputs. This approach enables
user flexibility, avoids the proliferation of hard-coded analysis routines,
and supports future extensibility without architecture changes. The initial
delivery (v1) provides a set of bundled `AnalysisGroup` types; full
user-configurable operator composition is planned for subsequent updates.

## 2 Requirements

### 2.1 Requirement: Composable operator framework

The Analysis system depends on simple, composable operators where each
operator performs a single, well-defined transformation. This enables:
- New analysis outputs via configuration rather than new code
- Testing of individual operations in isolation
- Reuse of common operations (spatial and temporal reductions, binary
  operations) across analysis computations

### 2.2 Requirement: Availability of all model variables

All simulation variables produced by the model and available for I/O in
Omega must be available to the Analysis module. Variables produced by the
Analysis system should also be available for further Analysis computation.

### 2.3 Requirement: Field access via dependency declaration

Operators must declare their input field dependencies at construction time.
During initialization, the orchestrator resolves dependencies and provides
operators with persistent pointers/references to input fields (from simulation
model fields or upstream operators). Operators retain these references and
access fields directly during `compute()`.

### 2.4 Requirement: Operator registration and factory

New operators must be registerable via a factory pattern, without changes to
the core analysis architecture. Operators self-register during initialization,
and the orchestrator queries the factory for operators by name. This
facilitates future extensibility; new operators integrate into the analysis
framework without modifying orchestration code.

### 2.5 Requirement: Multi-input and multi-output operators

Operators must be able to accept multiple input fields and produce multiple
output fields. Multi-input capability enables operators that combine fields
(e.g., binary operations, vector operations requiring multiple components).
Multi-output capability allows operators to simultaneously return separable
results (e.g., components of a vector field, or the components of a spatial
gradient).

### 2.6 Requirement: Computation caching

When multiple output streams or analysis fields depend on the same
intermediate result, that result must be computed once per timestep and
cached. Timestamp-based cache validation prevents stale results.

### 2.7 Requirement: Time operators

Time-based operations (mean, min, max over a period) must be regular
operators within the analysis framework, enabling composition with spatial
operations. Time period specification should be flexible (not limited to
hard-coded groups).

### 2.8 Requirement: Stream integration

Analysis fields must be integrated into the Omega output stream framework.
Configurable output stream parameters (filename, precision, period, etc.)
must be provided for fields produced by the analysis system. Fields will be
written to output with associated metadata.

### 2.9 Requirement: Polaris compatibility

Output from the Analysis module must be compatible with Polaris for
post-processing.

### 2.10 Requirement: Requested initial analysis capability

Initial delivery of the Analysis system will supply operators necessary for
computing a specified set of Analysis outputs:
- Global stats: global reduction to mean, min, max, and standard deviation
  of configurable fields
- AMOC: stream function for Atlantic meridional overturning circulation
- Eddy stats

## 3 Algorithmic Formulation

### 3.1 Operator Composition and Dependency Resolution

The Analysis system represents computations as a directed acyclic graph (DAG)
where nodes are operators and edges represent data dependencies. A single
Analysis field computation is defined by a string name that may expand into
multiple operators forming a chain.

#### 3.1.1 Operator dependencies

Each operator $\mathcal{O}_i$ produces one or more output fields and requires
zero or more input fields:

$$ \{\mathcal{O}_i^{\text{out},1}, \mathcal{O}_i^{\text{out},2}, \ldots\} =
f_i(\mathcal{I}_{i,1}, \mathcal{I}_{i,2}, \ldots, \mathcal{I}_{i,k}) $$

where each input $\mathcal{I}_{i,j}$ is either:
- A simulation field from the model (terminal node, no incoming operator edge)
- An output of another operator $\mathcal{O}_j$ (creating dependency edge
  $\mathcal{O}_j \to \mathcal{O}_i$)

**Operator chains:** A single Analysis field name $a$ may parse into an
ordered sequence of operators:

$$ a \xmapsto{\text{parse}} \{\mathcal{O}_1, \mathcal{O}_2, \ldots,
\mathcal{O}_m\} $$

where intermediate operators produce fields consumed by subsequent operators
in the chain, and only the terminal operator $\mathcal{O}_m$ writes to the
output stream.

**Shared intermediates:** When multiple Analysis fields require the same
intermediate result, the dependency resolver identifies structurally
equivalent operators via signature matching:

$$ \text{sig}(\mathcal{O}) = (\text{type}(\mathcal{O}),
\{\mathcal{I}_1, \mathcal{I}_2, \ldots\}) $$

Two operators with identical signatures are merged into a single node in the
DAG, preventing redundant computation.

> **v1 implementation note:** The full DAG construction algorithm below is
> the target design. The v1 implementation uses a simpler approximation:
> operator chains are parsed left-to-right and nodes are appended in natural
> dependency order; dependency edges are resolved post-hoc by matching
> operator input names against other operators' output names. Signature-based
> deduplication, cycle detection, and formal topological sort are planned for
> subsequent updates.

#### 3.1.2 Dependency graph construction

**Algorithm**: $\texttt{Analysis::buildDependencyGraph}$

Input: Set of requested Analysis field names $\mathcal{A} = \{a_1, a_2,
\ldots, a_n\}$ from all output streams

Output: Directed acyclic graph $\mathcal{G} = (\mathcal{V}, \mathcal{E})$
where $\mathcal{V}$ are operator nodes and $\mathcal{E}$ are data dependency
edges, with topological ordering $\pi : \mathcal{V} \to \mathbb{N}$

**Phase 1**: Parse and expand operator chains
1. Initialize: $\mathcal{V} \leftarrow \emptyset$, $\mathcal{E} \leftarrow
   \emptyset$, $\Sigma \leftarrow \emptyset$ (signature cache)
2. **For** each analysis field $a \in \mathcal{A}$:
    - Parse string into chain of operators: $\{\mathcal{O}_1, \ldots,
      \mathcal{O}_m\} \leftarrow \texttt{parseOperatorChain}(a)$
    - **For** $i = 1$ to $m$:
        - Compute signature: $s \leftarrow \text{sig}(\mathcal{O}_i)$
        - **If** $s \in \Sigma$ (operator already exists):
            - Retrieve existing node: $v \leftarrow \Sigma[s]$
            - **If** $i = m$ (final operator): Add $a$ to $v$'s output list
        - **Else** (create new node):
            - Create node: $v \leftarrow \text{OperatorNode}(\mathcal{O}_i)$
            - **If** $i = m$ (final operator): add output for node $v$ to
              stream for $a$, set alarm period
            - **Else** (intermediate operator): no stream output, computed
              on-demand when downstream alarm rings
            - Add to graph: $\mathcal{V} \leftarrow \mathcal{V} \cup \{v\}$
            - Cache signature: $\Sigma \leftarrow \Sigma \cup \{(s, v)\}$

**Phase 2**: Resolve dependencies
1. **For** each operator node $v \in \mathcal{V}$:
    - Let $\mathcal{I}(v) = \{\mathcal{I}_1, \ldots, \mathcal{I}_n\}$ be
      input fields for $v$
    - **For** each required input $\mathcal{I}_j \in \mathcal{I}(v)$:
        - **If** $\mathcal{I}_j$ is a simulation field from the model:
          terminal dependency, no edge needed
        - **Else if** $\exists\ u \in \mathcal{V}$ such that $\mathcal{I}_j
          \in \text{outputs}(u)$:
            - Add dependency edge: $\mathcal{E} \leftarrow \mathcal{E} \cup
              \{(u, v)\}$
            - Propagate alarms: **For** each $\text{Alarm} \in
              v.\text{ComputeAlarms}$:
                - **If** $\text{alarm} \notin u.\text{ComputeAlarms}$:
                    - $u.\text{ComputeAlarms} \leftarrow
                      u.\text{ComputeAlarms} \cup \{\text{Alarm}\}$
                     (upstream nodes observe all downstream alarms)
        - **Else** (field not found): ERROR

**Phase 3**: Validate acyclicity
1. Detect cycles using depth-first search with recursion stack:
$$ \text{hasCycle}(\mathcal{G}) = \begin{cases} \texttt{true} & \text{if }
\exists \text{ path } v_1 \to v_2 \to \cdots \to v_n \to v_1 \\
\texttt{false} & \text{otherwise} \end{cases} $$
    - **If** cycle detected: ERROR

**Phase 4**: Topological sort
1. Compute topological ordering $\pi : \mathcal{V} \to \{0, 1, \ldots,
   |\mathcal{V}|-1\}$ using Kahn's algorithm:
    - $\text{inDegree}(v) \leftarrow |\{u \in \mathcal{V} : (u,v) \in
      \mathcal{E}\}|$ for all $v$
    - $Q \leftarrow \{v \in \mathcal{V} : \text{inDegree}(v) = 0\}$
    - **While** $Q \neq \emptyset$:
        - Remove $v$ from $Q$; assign $\pi(v) \leftarrow \text{order}$;
          increment order; append $v$ to sorted list
        - For each $(v, w) \in \mathcal{E}$: decrement
          $\text{inDegree}(w)$; if zero, add $w$ to $Q$
    - **If** $|\text{sorted}| \neq |\mathcal{V}|$: ERROR (cycle)
2. Return $\mathcal{G}$ with ordering $\pi$

### 3.2 Operator Factory and Registration

The operator factory provides a runtime registry that maps operator type
names to constructor functions. This enables:
- **Decentralized registration**: Operators register themselves via a
  template helper before `main()` executes
- **Dynamic instantiation**: The orchestrator creates operators by name
  without hard-coded switch statements
- **Type-safe dispatch**: The factory selects the correct templated
  specialization based on the input field's runtime metadata (scalar type,
  rank, memory location)
- **Extensibility**: New operators can be added without modifying
  orchestration code

#### 3.2.1 Templated operator specializations

Analysis operators are class templates parameterized on the concrete Kokkos
array type `ArrayT` of their primary input field:

```c++
template<typename ArrayT>
class SpatialMaxOp : public AnalysisOperator { ... };
```

The factory registers all combinations of scalar type (I4/I8/R4/R8), rank
(1–5), and memory location (Device/Host/Both) for each operator template:

```c++
AnalysisOpFactory::registerAllArrayVariants<SpatialMaxOp>("SpatialMax");
```

At operator creation time, the factory inspects the primary upstream Field's
metadata to select the matching specialization.

#### 3.2.2 Registration

**Algorithm**: `AnalysisOpFactory::registerAllArrayVariants<OpT>(BaseName)`

1. Expand `OMEGA_ANALYSIS_ARRAY_TYPES` macro over all (DType, Rank, MemLoc,
   ArrayT) combinations
2. For each combination, call `registerOperator` with key
   `baseName + "_" + ArrayT + "_" + memloc` and a lambda that constructs
   `OpT<ArrayT>(UpstreamNames, Options)`
3. Validate key is unique (abort on duplicate)

All base analysis operators are registered at program startup by
`registerAllBaseAnalysisOperators()`, called from `Analysis::init()`.

#### 3.2.3 Factory operator creation

**Algorithm**: `AnalysisOpFactory::createOp`

**Input**: Operator type name, upstream field names, configuration options

**Output**: `unique_ptr<AnalysisOperator>` with the correct typed specialization

1. Retrieve the primary upstream Field from the Field registry
2. Extract `ArrayDataType`, rank, and `ArrayMemLoc` from Field metadata
3. Build fully-qualified type key: `OpType + "_" + ArrayTypeName + "_" + MemLoc`
4. Look up constructor in registry; abort if not found
5. Invoke constructor with `(UpstreamNames, Options)` and return result

### 3.3 Runtime Dispatch

The main Analysis computational loop is executed every timestep.
Dependencies are traversed recursively so that upstream operators are always
fresh when a downstream operator needs them. Caching prevents redundant work
when multiple downstream operators share an upstream.

> **v1 implementation note:** The loop over `SortedOperators` below assumes
> a topological ordering computed by `buildDependencyGraph`. In v1, nodes
> are iterated in insertion order (which is naturally dependency-correct for
> linearly chained operators). The full topological sort and
> `computeRecursive` are the target design.

**Algorithm**: `Analysis::computeAll`

**Input**: Topologically sorted operator list, current timestamp

**Output**: Updated Analysis fields written to registered output streams

1. **For** each $\texttt{Op} \in \texttt{SortedOperators}$:
    - **If** any alarm in $\texttt{Op.ComputeAlarms}$ is ringing:
        - $\texttt{computeRecursive(Op, TimeStamp)}$

2. $\texttt{computeRecursive(Op, TimeStamp)}$:
    - **If** $\texttt{Op.FieldComputed}$ **AND**
      $\texttt{Op.LastComputed == TimeStamp}$: return (cache hit)
    - **For** each $\texttt{UpstreamOp} \in \texttt{Op.Upstreams}$:
        - $\texttt{computeRecursive(UpstreamOp, TimeStamp)}$
    - $\texttt{Op.compute(TimeStamp)}$
    - $\texttt{Op.LastComputed} \leftarrow \texttt{TimeStamp}$;
      $\texttt{Op.FieldComputed} \leftarrow \texttt{true}$

### 3.4 Alarm Model

Each `OperatorNode` holds a vector of non-owning alarm pointers
(`vector<Alarm*> ComputeAlarms`). An operator is triggered when any of its
alarms rings.

- **Discrete-sampling (non-temporal-reduction) terminal operators**: borrow a
  raw pointer to the write alarm of the associated output stream. The stream
  owns this alarm.
- **Temporal reduction terminal operators**: require two alarms. An
  **accumulation alarm** controls how frequently a sample is added to the
  running sum; its interval is a user-configurable `AccumulationInterval`
  parameter (defaulting to every timestep in v1). An **output alarm**
  (borrowed from the associated stream, as for discrete-sampling operators)
  controls when the accumulated sum is divided by the sample count and
  written to output. Accumulation alarms are owned by the `Analysis` object
  as `vector<unique_ptr<Alarm>> AccumulationAlarms`. Each temporal reduction
  operator's `ComputeAlarms` vector contains two pointers: a raw pointer to
  its accumulation alarm and a raw pointer to its output alarm.
- **Intermediate (non-terminal) operators**: receive alarm pointers
  propagated from their downstream operators. Propagation is performed by
  `Analysis::propagateAlarmsUpstream()`, which iterates until no further
  changes occur.

This design ensures alarms have a clear single owner (a stream or the Analysis
class) while allowing any number of operator nodes to observe them.

> **v1 constraint:** Temporal reduction periods must be evenly divisible into the
> restart interval. This is validated during `createAnalysisGroupStreams()` to
> ensure proper checkpoint/restart behavior.

### 3.5 AnalysisGroup Configuration

Each child node of the `Analysis:` group in the configuration YAML file
represents an `AnalysisGroup`. The orchestrator iterates over these nodes during
initialization and dispatches to the appropriate handler:

- **Named pre-defined groups** (e.g. `GlobalStats`): Dispatched by name to
  a derived `AnalysisGroup` subclass. The subclass reads its own config
  parameters and constructs the appropriate operator chains and output
  streams internally.
- **Custom user-defined groups** (future): Config nodes not matching a
  pre-defined name will be parsed as user-defined groups of composable
  operator chains using the full chain-parsing and DAG machinery described
  in Section 3.1.

Example config structure:

```yaml
Omega:
  Analysis:
    GlobalStats:              # pre-defined bundled group
      Enable: true
      Fields: ["NormalVelocity", "Temperature", "Salinity"]
      SpatialStats: ["Max", "Min", "Mean", "StdDev"]
      ReductionPeriod: ["1Day", "1Month"]
      SampleFreq: ["1Hour"]
      Filename: global.stats.$Y
      Stream:                 # define optional stream parameters
        FileFreq: 1
        FileFreqUnits: years
    MyCustomGroup:            # future: user-defined composable group
      Enable: false
      OperatorChains:         # final DSL syntax to be determined
        - "FieldA_Op1_Op2(FieldB)"
        - "Op3(FieldC,FieldD)_Op4"
      Filename: custom.analysis.$Y.$M
```

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Configuration

The `Analysis` config node is a map of group names to group-specific config
sub-nodes. Each group sub-node must contain at minimum an `Enable` boolean
key. Additional keys are group-specific. The `AnalysisGroup` base class
provides a `StreamParams` helper for translating group config options into
`IOStream::create` arguments. Each group may generate multiple output
streams depending on its configuration; for example, a `GlobalStats` group
with multiple reduction periods (e.g., `["1Day", "1Month"]`) will create
separate streams for each period, grouping operator chains by their output
frequency and whether they perform temporal reduction (e.g., TimeMeanOp) or
discrete sampling (i.e., instantaneous snapshots).

#### 4.1.2 Classes

##### AnalysisOperator

The `AnalysisOperator` class is the abstract base class from which all
concrete operators are derived. It is parameterized on the Kokkos array
type `ArrayT` in derived classes. Output field data arrays are allocated as
members of the derived class and created in the constructor; the Field
registry entry is also created at construction time. The `initialize()`
method is called after all fields exist, primarily to store mesh/env
pointers needed by `compute()`.

```c++
// Temporal operators have an accumulation phase and an operation/output phase

class AnalysisOperator {
 public:
   AnalysisOperator();

   ~AnalysisOperator();

   /// Return name for this operator type
   const std::string getOperatorType();

   /// Return unique name for this operator instance.
   /// Derived from the concatenated upstream field names and operator type,
   /// e.g. "Temperature_SpatialMean_TimeMean1Day"
   const std::string getName();

   /// Return names of fields required by this operator
   const std::vector<std::string> getInputFieldNames();

   /// Return names of output fields produced by this operator
   const std::vector<std::string> getOutputFieldNames();

   /// Returns true if the output field has already been computed for TimeStamp
   bool isCacheValid(const TimeInstant &TimeStamp);

   /// Initialize operator: store mesh/env pointers needed by compute().
   virtual void initialize(const MachEnv *InEnv,
                           const HorzMesh *Mesh,
                           const VertCoord *VCoord,
                           Config Options);

   /// Set period alarm for temporal reduction operators
   /// Default implementation does nothing (non-temporal operators ignore this)
   virtual void setPeriodAlarm(Alarm *Alarm);

   /// Perform computation of Analysis fields. Retrieves input data from the
   /// Field registry using input field names. Writes to operator-owned output
   /// arrays which are attached to the Field registry.
   virtual void compute(const TimeInstant &TimeStamp) = 0;

 protected:
   std::string OperatorTypeName;
   std::string InstanceName;
   std::vector<std::string> InputNames;
   std::vector<std::string> OutputNames;

   TimeInstant LastComputed;
   bool FieldComputed;
};
```

Helper utilities for building operator `Config` objects inline:

```c++
// Create a Config from key-value pairs
// Usage: makeOpConfig(opParam("Period", "1day"), opParam("Layer", 10))
template<typename T>
OpParam<T> opParam(std::string Key, T&& Value);

template<typename T, typename... Args>
Config makeOpConfig(const std::pair<std::string, T>& Param, Args... OtherArgs);
```

These helpers enable in-code construction of `Config` objects for passing
parameters to operator constructors, using the same YAML-based `Config`
interface that reads from configuration files. This provides a uniform
parameter-passing mechanism: operators receive a `Config` object whether
instantiated from user config or programmatically by a bundled
`AnalysisGroup`. The pattern avoids constructor signature proliferation as
operators gain parameters, maintains type safety via `Config::get<T>()`, and
allows operator-specific validation and defaults to be centralized in the
constructor.

##### Example derived operator — SpatialMaxOp

```c++
template<typename ArrayT>
class SpatialMaxOp : public AnalysisOperator {
 public:
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructor: sets InputNames, creates output Field and data array.
   /// InstanceName = UpstreamNames[0] + "_SpatialMax"
   SpatialMaxOp(const std::vector<std::string> &UpstreamNames,
                Config Options);

   /// Retrieves typed input array from the Field registry and calls
   /// globalMaxVal() to compute the MPI-global maximum.
   void compute(const TimeInstant &TimeStamp) override;

 private:
   const HorzMesh *Mesh;
   const VertCoord *VCoord;
   MPI_Comm Comm;

   /// Output data — one scalar value stored as a 1D Array of length 1
   typename Array1D<ScalarT>::type OutputData;
   ScalarT SpatialMax;
};
```

##### AnalysisOpFactory

Factory class for creating `AnalysisOperator` instances. The class
itself is a singleton with all static methods; internally it maintains a
Meyer's singleton registry map. The factory dispatches to the correct
templated specialization at runtime by inspecting the primary upstream
Field's metadata.

```c++
class AnalysisOpFactory {
 public:
   using CreatorFunc = std::function<std::unique_ptr<AnalysisOperator>(
       const std::vector<std::string> &UpstreamNames, Config Options)>;

   /// Register a single operator variant by string label
   static void registerOperator(const std::string &Label,
                                CreatorFunc Creator);

   /// Create an operator instance. Inspects Field metadata of UpstreamNames[0]
   /// to select the correct templated specialization.
   static std::unique_ptr<AnalysisOperator> createOp(
       const std::string &OpType,
       const std::vector<std::string> &UpstreamNames,
       Config Options
   );

   /// Register all scalar type × rank × memory location variants of a
   /// templated operator class.
   /// Usage: registerAllArrayVariants<SpatialMaxOp>("SpatialMax");
   template<template<typename> class OperatorTemplate>
   static void registerAllArrayVariants(const std::string &BaseName);

   /// Check if operator type is registered
   static bool hasOperator(const std::string &Type);

 private:
   static std::map<std::string, CreatorFunc>& registry(); // Meyer's singleton
   static std::string getArrayTypeName(ArrayDataType DType,
                                       I4 Rank,
                                       ArrayMemLoc MemLoc);
};
```

All base analysis operators are registered by calling:

```c++
void Analysis::registerAllBaseAnalysisOperators();
```

from `Analysis::init()` before any operators are instantiated.

##### OperatorNode

Internal representation of a node in the Analysis operator graph.

```c++
struct OperatorNode {
   std::unique_ptr<AnalysisOperator> Op; ///< Operator instance (owned)
   std::vector<OperatorNode*> Upstreams; ///< Upstream dependencies (non-owning)
   std::vector<std::string> StreamNames; ///< Associated output stream names
   std::vector<Alarm*> ComputeAlarms;    ///< Alarms triggering compute (non-owning)
};
```

Operators with a non-empty `StreamNames` vector are terminal nodes whose
output is written to one or more output streams. Operators with an empty
`StreamNames` vector are intermediate nodes computed on demand when a
downstream alarm rings.

##### AnalysisGroup

`AnalysisGroup` is the abstract base class for bundled analysis groups. In
v1, concrete derived classes (e.g. `GlobalStats`) encapsulate the config
parsing, operator construction, and stream creation for a named analysis
group. In the future, the same base class will support user-defined custom
groups specified entirely in config, where the user supplies composable
operator chains within the group's config node.

The base class provides a `StreamParams` helper for translating group config
into `IOStream::create` arguments, and `createAnalysisGroupStreams()` which
groups operator chains by their output period and type, validates temporal
reduction periods against the restart interval, and creates the associated
`IOStream` objects.

```c++
class AnalysisGroup {
 public:
   virtual ~AnalysisGroup() = default;

   std::string getName();

   /// Groups operator chains by stream characteristics, creates IOStream
   /// objects, associates operator output fields with streams, and stores
   /// AnalysisStream metadata on the Analysis orchestrator.
   void createAnalysisGroupStreams(
       const std::string &GroupName,
       Config &AnalysisGroupOptions,
       Analysis *AnalysisMgr
   );

 protected:
   /// Metadata about a single operator chain within this group
   struct OpChainInfo {
      std::string ChainStr; ///< Operator instance name (output field name)
      std::string FreqStr;  ///< Period/frequency string, e.g. "1day", "6hour"
      bool IsTimeReduction; ///< true = temporal reduction; false = discrete sample
   };

   /// Template for constructing an IOStream config for this group's output
   struct StreamParams {
      StreamParams();  // default values for all IOStream options
      void apply(const std::map<std::string, std::string> &Overrides);
      Config toConfig() const;
      std::map<std::string, std::string> Params;
   };

   std::string GroupName;
   std::vector<OpChainInfo> OpChainInfos; ///< All operator chains in this group
};
```

##### GlobalStats (derived AnalysisGroup)

`GlobalStats` is the first concrete `AnalysisGroup` subclass. It reads
`Fields`, `SpatialStats`, `ReductionPeriod`, and `SampleFreq` from the group
config and constructs a matrix of spatial-reduction operator chains, each
optionally followed by a temporal reduction operator. The `ReductionPeriod`
parameter specifies temporal reduction intervals (e.g., "1Day", "1Month")
for outputs computed by temporal reduction operators such as `TimeMeanOp`,
while the `SampleFreq` parameter specifies discrete sampling intervals for
instantaneous snapshots of the analysis fields.

```c++
class GlobalStats : public AnalysisGroup {
 public:
   GlobalStats(const std::string &GroupName,
               Config &AnalysisGroupOptions,
               Analysis *AnalysisMgr);
   ~GlobalStats() = default;
};
```

For each `(field, stat, period)` combination, the constructor builds a chain
string of the form `FieldName_SpatialStat_TimeMeanPeriod` and calls
`AnalysisMgr->parseChainAndBuildOps()`. For each
`(field, stat, samplefreq)` combination, it builds `FieldName_SpatialStat` chains.
After all chains are registered, it calls `createAnalysisGroupStreams()`.

##### Analysis

`Analysis` is the top-level orchestrator class. It owns the `OperatorNode`
list, the accumulation alarms for temporal reduction operators. It is
responsible for reading the config, constructing `AnalysisGroup` instances,
resolving operator dependencies, and scheduling compute calls via the
alarm system.

```c++
class Analysis {
 public:
   /// Initialize the Analysis module: register all base operators,
   /// retrieve mesh/vertcoord/clock, create the Default Analysis instance.
   /// Must be called after HorzMesh, VertCoord, and TimeStepper are initialized.
   static void init();

   /// Create a named Analysis instance
   static Analysis *create(const std::string &Name,
                           const MachEnv *Env,
                           const HorzMesh *Mesh,
                           const VertCoord *VCoord,
                           Clock *ModelClock,
                           Config *Options);

   /// Called each timestep to trigger all operators whose alarms are ringing
   void computeAll();

   /// Parse an underscore-delimited operator chain string and register all
   /// operators in the chain that do not yet exist as Fields
   void parseChainAndBuildOps(const std::string &OpChainStr);

   /// Instantiate a single operator and append it as an OperatorNode
   void registerAnalysisOp(const std::string &OpName,
                            const std::vector<std::string> &UpstreamNames,
                            Config Options);

   /// Get a pointer to the model clock (used by AnalysisGroup for stream creation)
   Clock *&getModelClock();

   /// Check whether a node with FullOpName is already registered
   bool OpNodeExists(const std::string &FullOpName);

   static Analysis *getDefault();
   static void finalize();
   ~Analysis();

 private:
   /// Accumulation alarms owned by Analysis for temporal reduction operators
   std::vector<std::unique_ptr<Alarm>> AccumulationAlarms;

   static Analysis *DefAnalysis;
   static std::map<std::string, std::unique_ptr<Analysis>> AllAnalysisObjects;

   Analysis(const std::string &Name,
            const MachEnv *Env,
            const HorzMesh *Mesh,
            const VertCoord *VCoord,
            Clock *ModelClock,
            Config *Options);

   std::string Name;
   Clock *ModelClock;
   const HorzMesh *Mesh;
   const VertCoord *VCoord;

   /// All registered operator nodes
   std::vector<std::unique_ptr<OperatorNode>> OpNodes;

   // Private Methods

   /// Register all built-in operator types with the AnalysisOpFactory
   static void registerAllBaseAnalysisOperators();

   /// Post-hoc dependency resolution: match input field names against
   /// other nodes' output field names to populate Upstreams vectors.
   void buildOperatorDependencies();

   /// Set ComputeAlarms on terminal nodes and propagate alarms upstream.
   void setComputeAlarms();

   /// Iteratively propagate downstream alarms to upstream nodes
   void propagateAlarmsUpstream();

   Analysis(const Analysis &) = delete;
   Analysis(Analysis &&)      = delete;
};
```

### 4.2 Operator chain string convention

Operator instance names (and the names of the Fields they produce) follow
the convention that each component is separated by an underscore character:

```
FieldName_Op1[Params]_Op2[Params]...
```

Examples:
- `Temperature_SpatialMax` — spatial maximum of Temperature
- `NormalVelocity_SpatialMean_TimeMean1day` — 1-day time average of the
  spatial mean of NormalVelocity
- `PseudoThickness_SpatialStdDev` — spatial standard deviation of
  PseudoThickness (implicitly requires `PseudoThickness_SpatialMean` as a
  shared intermediate)

The `parseChainAndBuildOps()` method splits on `_`, reconstructs the running
prefix at each node, and creates an operator only if the corresponding output
Field does not already exist — enabling natural sharing of intermediate
results without an explicit signature cache.

> **Note on operator chain syntax**: The exact form of operator chain strings
> shown in examples throughout this document represents a preliminary syntax
> for the v1 implementation. The final syntax for fully composable
> user-defined operator chains will be refined in future versions. The current
> v1 implementation focuses on pre-defined bundled groups (e.g.,
> `GlobalStats`) with group-specific configuration parameters.

## 5 Verification and Testing

### 5.1 Test: Individual operator correctness

For each operator type (SpatialMax, SpatialMin, SpatialMean, SpatialStdDev,
TimeMean in the first batch), construct a small test mesh with analytic field
values. Call `compute()` directly and verify output against a known-answer
solution. For TimeMean specifically, verify accumulation over multiple
timesteps, verify correct mean calculation at period end, and test with
different `AccumulationInterval` settings. This unit test validates each
operator in isolation before integration testing.

### 5.2 Test: Dependency resolution and execution order

Create configurations with shared intermediate operators (e.g.,
`Field_SpatialMean_TimeMean1day` and `Field_SpatialStdDev` both requiring
`Field_SpatialMean`). Verify that `buildOperatorDependencies()` correctly
populates the `Upstreams` vectors, that intermediate results are computed
exactly once per timestep (cache validation), and that upstream operators
complete before downstream operators execute (correct execution order). This
test verifies DAG construction and cache-based deduplication.

### 5.3 Test: Alarm system

Create operators with multiple downstream consumers at different frequencies.
Verify that `propagateAlarmsUpstream()` correctly propagates alarms from
terminal nodes to all upstream dependencies. Verify that `setPeriodAlarm()`
correctly injects period alarms for temporal reduction operators. Verify that
TimeMeanOp correctly accumulates samples during the accumulation phase and
finalizes when the period alarm rings. Verify that operators with multiple
alarms in `ComputeAlarms` trigger when ANY alarm rings. Verify that
intermediate (non-terminal) operators with empty `StreamNames` are computed
on-demand when downstream alarms ring and do not create output files. This
test verifies the alarm-driven scheduling mechanism.

### 5.4 Test: Factory registration and type dispatch

Verify that all base analysis operators register correctly via
`registerAllBaseAnalysisOperators()`. Verify that the factory can instantiate
operators for all supported array types (I4/I8/R4/R8, ranks 1-5,
Device/Host/Both). Verify that `AnalysisOpFactory::createOp()` correctly
inspects upstream Field metadata (scalar type, rank, memory location) and
selects the matching template specialization. Verify that appropriate errors
are produced when requesting unregistered operator types or array type
combinations. This test verifies the extensibility mechanism and type-safe
dispatch.

### 5.5 Test: Configuration parsing and validation

Verify that `parseChainAndBuildOps()` correctly handles valid operator chain
strings and reuses existing intermediate Fields rather than creating
duplicates. Verify that `parseChainAndBuildOps()` produces informative error
messages for unrecognized operator names or missing input fields. Verify that
`makeOpConfig()` and `opParam()` helper functions correctly construct Config
objects for inline parameter passing. Verify that operator constructors
correctly extract and validate parameters from Config objects, with appropriate
error handling for missing required parameters or invalid types. Verify that
`createAnalysisGroupStreams()` correctly groups operator chains by period and
type, validates temporal reduction periods against the restart interval via
`TimeInterval::isDivisibleBy()`, and creates the expected set of IOStream
objects. Verify that `StreamParams::apply()` correctly overrides default stream
parameters with group-specific configuration. This test verifies the user
interface and configuration system.

### 5.6 Test: End-to-end integration

Complete system test exercising all components from configuration parsing
through NetCDF output for global statistics. Advance the clock through one
or more output periods, and verify that output files contain the expected
fields with correct values. This test validates the complete workflow with
real mesh and I/O.

### 5.7 Test: Advanced DAG features (future)

Once the full DAG construction algorithm is implemented, create configurations
with circular dependencies and verify that cycle detection produces appropriate
errors. Test signature-based deduplication to ensure structurally equivalent
operators are merged into single nodes. Verify formal topological sort produces
correct execution ordering for complex DAGs. This test validates future
enhancements to dependency resolution.
