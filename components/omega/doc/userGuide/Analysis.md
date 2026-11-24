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
 configured with the `ReductionPeriod` parameter. A reduction period may span a
 restart: the mid-period accumulation state (the partial sum and sample count)
 is checkpointed in the restart file and restored on continuation, so the final
 average is identical to an uninterrupted run. There is no requirement that the
 reduction period divide evenly into the restart interval.

- **Instantaneous output:** Snapshots of analysis fields written at specified
 frequencies (e.g., every 6 hours). Instantaneous output is configured with the
 `SnapshotPeriod` parameter.

(omega-user-analysis-sampling-interval)=
### Sampling interval for temporal reductions

The `SamplingInterval` parameter set the temporal reduction accumulation
cadence. The values `dt` and `timestep` are accepted as keys to accumulate
a sample on **every** model timestep. Setting `SamplingInterval` to one of the
accepted forms below instead collects samples at a coarser cadence, so the time
average is formed from fewer, evenly-spaced samples. This is the in-model
equivalent of MPAS-Ocean's driver-level `compute_interval` for the
`time_series_stats` analysis member.

`SamplingInterval` accepts three forms, distinguished by suffix:

| Form               | Example        | Meaning                              |
|--------------------|----------------|--------------------------------------|
| TimeInterval units | `6hour`        | a fixed time interval between samples |
| Timestep count     | `4steps`       | every 4 model timesteps              |
| Samples per period | `8perPeriod`   | 8 evenly-spaced samples per reduction period |

Behavior and constraints:

- **Omitted/empty:** sample every timestep.
- **Phasing:** the first sample falls one `SamplingInterval` after the start of
  each reduction period (samples at `PeriodStart + n*SamplingInterval`).
- **Snap-to-step:** like all model alarms, sampling fires on the first timestep
  at or after each target time. A `SamplingInterval` that is not a whole
  multiple of the timestep therefore snaps to the next step, giving a bounded
  sub-timestep phase error (zero when the interval is a timestep multiple).
- **Samples per period:** `<n>perPeriod` divides a period given in days or
  shorter units into `n` equal intervals. Months vary in length, so for a
  period given in months or years `n` must divide the number of months, as
  in `3perPeriod` for `1Year`; otherwise use a fixed interval such as `1day`.
- **Period-divisibility (required):** the reduction period must be evenly
  divisible by the sampling interval, so the simple count-mean remains an
  unbiased estimate of the time average. A violation aborts the run. **This
  check is unrelated to the restart interval** — temporal reductions may still
  span restarts with no restart-interval divisibility requirement.
- **Accuracy caveat:** the averaged value stays unbiased, but when the interval
  is not a timestep multiple the per-period sample count may differ by ±1
  between periods (a boundary-adjacent sample can snap to either side). This
  per-period count variation is **reproducible across restarts**: the sample
  schedule is anchored to the absolute clock start, so an exact-restart (ERS)
  run reconstructs the identical boundaries and reproduces the straight-through
  result bit-for-bit. When the interval is a timestep multiple, the count is
  constant and this does not arise at all.


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
Temporal reduction operators automatically checkpoint their mid-period
accumulation state to the restart file, so reduction periods may span restarts
without any divisibility requirement.

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
user-specified fields. The mean and standard deviation are weighted so that
they do not depend on the mesh resolution: by the area of each cell, edge or
vertex for a horizontal field, and by mass (area times pseudo-thickness) for
a field with a vertical dimension, with an interface field weighted by half
the mass of each adjacent layer. Only active layers take part. Each statistic
is written with the `units` and `standard_name` of the field it reduces and a
CF `cell_methods` attribute describing the reduction, for example
`area: depth: mean` for the mean of a layered field and
`area: depth: mean time: mean` for its time mean.

**Example:**

```yaml
Omega:
  Analysis:
    GlobalStats:
      Enable: true
      Fields: [Temperature, GeomZInterface, Salinity, NormalVelocity]
      SpatialStats: [Mean, Max, Min, StdDev]
      ReductionPeriod: [1Day, 1Month]
      SamplingInterval: "dt"
      SnapshotPeriod: [6Hour]
      Filename: global.stats.$Y.nc
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
  (e.g., `1Hour`, `1Day`, `1Month`, `1Year`). At least one of `ReductionPeriod`
  or `SnapshotPeriod` must be specified.

- **SamplingInterval:** Cadence at which the temporal reductions accumulate
  samples; see
  [Sampling interval for temporal reductions](#omega-user-analysis-sampling-interval).
  Empty, `dt`, or `timestep` samples every timestep.

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

### MOC

Computes the Meridional Overturning Circulation (MOC) streamfunction using two
methods: latitude-binned regional MOC and transect-based MOC. The MOC
represents zonally integrated meridional mass transport as a function of
latitude and depth (regional MOC) or depth alone (transect MOC). Output is in
Sverdrups (Sv), where 1 Sv = 10⁶ m³/s, so its `units` attribute is
`1e+06 m3 s-1`, the udunits spelling of a Sverdrup (udunits reads `Sv` as the
sievert).

**Example:**

```yaml
Omega:
  Analysis:
    MOC:
      Enable: true
      NumBins: 180
      MinLat: -90.0
      MaxLat: 90.0
      Regions: [Global]
      Transects: []
      ReductionPeriod: [1Month]
      SamplingInterval: "dt"
      SnapshotPeriod: [1Day]
      Filename: moc.$Y-$M.nc
      Stream:
        FileFreq: 1
        FileFreqUnits: months
```

**Group-Specific Parameters:**

- **NumBins:** Optional integer specifying the number of latitude bins for
  regional MOC (default: 180, approximately 1-degree resolution).

- **MinLat:** Optional minimum latitude in degrees for binning (default: -90.0).

- **MaxLat:** Optional maximum latitude in degrees for binning (default: 90.0).

- **Regions:** Optional list of region names for regional MOC computation
  (default: `[Global]`). Regional MOC computes the streamfunction as a function
  of latitude and depth. Currently only `Global` is supported; region masks are
  not yet implemented.

- **Transects:** Optional list of transect names for transect-based MOC
  computation (default: `[]`). Transect MOC computes the streamfunction as a
  function of depth only by accumulating transport across specified transect
  edges. Transect masks are not yet implemented.

- **ReductionPeriod:** Optional list of time periods for temporal reduction.
  At least one of `ReductionPeriod` or `SnapshotPeriod` must be specified.

- **SamplingInterval:** Cadence at which the temporal reductions accumulate
  samples; see
  [Sampling interval for temporal reductions](#omega-user-analysis-sampling-interval).
  Empty, `dt`, or `timestep` samples every timestep.

- **SnapshotPeriod:** Optional list of intervals for instantaneous output. At
  least one of `ReductionPeriod` or `SnapshotPeriod` must be specified.

**Output fields:**

- Regional MOC: `MOC_streamfunction_RegionName` (2D: latitude bins × depth)
- Transect MOC: `MOC_streamfunction_transect_TransectName` (1D: depth only)

For temporal reduction, `_TimeMeanPeriod` is appended (e.g.,
`MOC_streamfunction_Global_TimeMean1Month`).

**Output streams:**
Automatically created and named:
 - Time reduction: `MOC_FreqTimeStats` (e.g., `MOC_1MonthTimeStats`)
 - Instantaneous output: `MOC_FreqInstants` (e.g., `MOC_1DayInstants`)

```{note}
Regional and transect masks are not yet implemented. The `Regions` parameter
currently only supports `[Global]`, and `Transects` is a placeholder for future
functionality. The underlying operator infrastructure is complete and ready for
use once the mask fields are implemented in Omega.
```

### MonthlyAverages

Computes monthly time averages of a set of requested global fields. Unlike
`GlobalStats`, no spatial reduction is applied: each requested field is passed
directly to a temporal-reduction (`TimeMean`) operator. The averaging period
defaults to one month, but may be overridden with the optional
`ReductionPeriod` parameter. The `MonghtlyAverages` group is specifically for
temporal reductions, any defined `SnapshotPeriod` parameter is ignored. This
group mirrors the `timeSeriesStatsMonthly` analysis member in MPAS-Ocean, which
delivers monthly means of prognostic and derived fields.

**Example:**

```yaml
Omega:
  Analysis:
    MonthlyAverages:
      Enable: true
      Fields: [Temperature, Salinity, NormalVelocity, PseudoThickness]
      ReductionPeriod: [1Month]   # optional; defaults to [1Month]
      SamplingInterval: "dt"
      Filename: monthly.averages.$Y.nc
      Stream:
        FileFreq: 1
        FileFreqUnits: years
```

**Group-Specific Parameters:**

- **Fields:** Required list of field names to average. The simulation will abort
  if a requested field does not exist.

- **ReductionPeriod:** Optional list of time periods for temporal reduction
  (e.g., `1Day`, `1Month`, `1Year`). Defaults to `[1Month]` when not specified.
  Any `SnapshotPeriod` supplied for this group is ignored.

- **SamplingInterval:** Cadence at which the temporal reductions accumulate
  samples; see
  [Sampling interval for temporal reductions](#omega-user-analysis-sampling-interval).
  Empty, `dt`, or `timestep` samples every timestep.

**Output fields:**
For each requested field, a time-averaged field is produced:
 - `FieldName_TimeMeanPeriod` (e.g., `Temperature_TimeMean1Month`)

**Output streams:**
Automatically created and named:
 - Time reduction: `MonthlyAverages_FreqTimeStats` (e.g.,
 `MonthlyAverages_1MonthTimeStats`)

```{note}
The default field list in `Default.yml` includes the core prognostic fields
(`Temperature`, `Salinity`, `NormalVelocity`, `PseudoThickness`). Additional
registered global fields may be requested by adding them to the `Fields` list.
A reduction period may span a restart; see
[Temporal Reduction and Restarts](#temporal-reduction-and-restarts) below.
```

## Usage Notes

### Temporal Reduction and Restarts

Temporal reductions (`ReductionPeriod`) may span a model restart. There is no
requirement that the reduction period divide evenly into the restart interval,
so any combination of periods and restart intervals is valid (for example, a
`1Month` average with a 10-day restart interval). When a restart is written in
the middle of a reduction period, the mid-period accumulation state (the
partial sum and the sample count) is checkpointed in the restart file and
restored on continuation, so the final written average is identical to that of
an uninterrupted run. Restart files written before this capability was added
remain readable; a reduction that lacks checkpointed state simply begins a
fresh accumulation period on restart.

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
