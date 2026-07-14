(omega-user-analysis)=

# Analysis

The Analysis module provides in-situ computation of diagnostic fields from the
ocean model state during a simulation. Analysis fields are computed on the fly
at user-specified frequencies and written to output streams, eliminating the
need to store large volumes of raw model output for offline post-processing.

In the initial delivery, analysis is configured through pre-defined bundled
**AnalysisGroups** such as `GlobalStats`. These bundles provide convenient
shorthands for commonly-used analysis patterns and are configured entirely
through the YAML input file. Future updates will add user-defined composable
operator chains for fully customizable analysis, though the capability of
using pre-defined bundled groups will persist as convenient defaults.

## Background

The Analysis module supports two output modes:

- **Temporal reduction:** Time-averaged output over a specified period (e.g.,
 daily or monthly means). Fields are accumulated each timestep; the reduction
 is finalized and written at the end of each period. Temporal reduction is
 configured with the `ReductionPeriod` parameter. When using temporal
 reduction, the reduction period must be evenly divisible into the restart
 interval. This ensures that temporal reduction states can be properly
 checkpointed and restored. For calendars with variable month lengths, restart
 intervals based on months can only be divided by reduction periods based on
 months or periods of 1 day or whole-number fractions of a day, e.g., 6 hours.
 The Analysis system validates that the reduction period divides evenly into
 the restart interval during initialization and will abort if this requirement
 is violated.

- **Instantaneous output:** Snapshots of analysis fields written at specified
 frequencies (e.g., every 6 hours). Instantaneous output is configured with the
 `SnapshotPeriod` parameter.

## Configuration Overview

Analysis groups are configured in the YAML input file under `Omega: Analysis:`.
Each group can be independently enabled or disabled. The initial delivery
provides the `GlobalStats` group.

### Common Configuration Parameters

All AnalysisGroups support the following parameters:

- **Enable:** Required boolean (`true` or `false`). Enables or disables the
  group.

- **Filename:** Required string. Template for output filenames with optional
  time-based variables (`$Y`, `$M`, `$D`, `$h`, `$m`, `$s`). The group
  automatically appends frequency and type information to the filename.

- **Stream:** Optional map of IOStream parameters to customize output behavior.
  See [Stream Parameters](#stream-parameters) below for details.

## Analysis Groups

Each AnalysisGroup encapsulates a set of related analysis operations. Groups
inherit from the `AnalysisGroup` base class, which provides common
functionality for operator chain construction and stream management.

When initialized, each group:
1. Parses its configuration to determine which operator chains to build
2. Registers operator chains with the Analysis orchestrator
3. Groups chains by output frequency and type (time reduction vs. instantaneous
output)
4. Creates IOStream objects for each unique temporal configuration

All groups automatically handle intermediate operator sharing—if multiple
chains require the same intermediate result, it is computed once and reused.
Groups also validate that temporal reduction periods divide evenly into the
restart interval to ensure proper checkpoint/restart behavior.

### Stream Parameters

All AnalysisGroups support the optional `Stream` sub-node for customizing
IOStream output parameters. If omitted, defaults are used.

**Default values set by the Analysis module:**
- `UsePointerFile: false`
- `IfExists: append`
- `Precision: double`
- `UseStartEnd: false`

All other Stream parameters (`FileFreq`, `FileFreqUnits`, `PointerFilename`,
`StartTime`, `EndTime`) are empty or disabled by default.

**Note:** The `Filename`, `Freq`, and `FreqUnits` stream parameters are set
automatically by `AnalysisGroup` configuration. Any user-specified values for
these parameters will be overwritten.

For complete IOStream documentation, see the
[IOStreams](#omega-user-iostreams) section of the User's Guide.

### GlobalStats

Computes global spatial statistics (mean, min, max, standard deviation) for
user-specified fields.

**Example:**

```yaml
Omega:
  Analysis:
    GlobalStats:
      Enable: true
      Fields: [Temperature, GeomZInterface, Salinity, NormalVelocity]
      SpatialStats: [Mean, Max, Min, StdDev]
      ReductionPeriod: [1Day, 1Month]
      SnapshotPeriod: [6Hour]
      Filename: global.stats.$Y
      Stream:
        FileFreq: 1
        FileFreqUnits: years
        Precision: single
```

**Group-Specific Parameters:**

- **Fields:** Required list of field names to analyze. The simulation will
  abort if a requested field does not exist.

- **SpatialStats:** Required list of spatial operators.

- **ReductionPeriod:** Optional list of time periods for temporal reduction
  (e.g., `1Hour`, `1Day`, `1Month`, `1Year`). Each period must divide evenly
  into the restart interval. At least one of `ReductionPeriod` or
  `SnapshotPeriod` must be specified.

- **SnapshotPeriod:** Optional list of intervals for instantaneous output
  (e.g., `1Hour`, `6Hour`, `1Day`). At least one of `ReductionPeriod` or
  `SnapshotPeriod` must be specified.

**Output fields:**
For each `(Field, SpatialStat)` combination:
 - Instantaneous output: `FieldName_SpatialStat` (e.g.,
 `Temperature_SpatialMean`)
 - Temporal reduction: `FieldName_SpatialStat_TimeMeanPeriod` (e.g.,
 `Temperature_SpatialMean_TimeMean1Day`)

**Output streams:**
Automatically created and named:
 - Time reduction: `GlobalStats_FreqTimeStats` (e.g.,
 `GlobalStats_1DayTimeStats`)
 - Instantaneous output: `GlobalStats_FreqInstants` (e.g.,
 `GlobalStats_6HourInstants`)

## Usage Notes

### Temporal Reduction Period Constraint

When using temporal reduction (`ReductionPeriod`), the period must divide
evenly into the restart interval to ensure proper checkpoint/restart behavior.
The Analysis module validates this during initialization. For example, if the
restart interval is 6 months, valid periods include `6hours`, `1Day`, `1Month`,
`2Months`, `3Months`, and `6Months`, but `5Days` or `7Months` are invalid.

### Field Availability

The Analysis module initializes after the model state and validates that all
requested fields exist. If a requested field is not found, the simulation
aborts during initialization with an error message.

### Intermediate Results

Some analysis operators require intermediate computations. For example,
computing the standard deviation requires the mean as input. Currently, when
configuring the `SpatialStats` list, users must include both `Mean` and
`StdDev` if standard deviation is desired. The mean computation will be
automatically shared if multiple statistics require it.

## Future Extensions

The Analysis framework is designed for extensibility. Future versions will
include:

- **User-defined composable operator chains:** Full control over operator
  composition specified directly in the configuration file, enabling custom
  analysis workflows beyond the pre-defined bundles.

- **Additional pre-defined analysis groups:** Bundled groups for common
  analysis patterns such as AMOC stream function, eddy statistics, regional
  diagnostics, and more.

- **Additional operators:** Vertical operations, horizontal gradients, spatial
  filtering, and custom transformations.

The bundled `AnalysisGroup` types like `GlobalStats` will persist as convenient
shorthands for commonly-used analysis outputs, even as the composable operator
capability becomes available.

For detailed information about the Analysis architecture, operator composition,
dependency resolution, and extending the module with new operators or groups,
see the [Analysis](#omega-dev-analysis) section of the Developer's Guide.
