# User Guide

E3SM is not just one earth system model but a modeling system that allows
many different configurations of atmosphere, ocean, land and other
components with both full model and data model options. Also, the configurations of model
components can run at different resolutions.  Some configurations
can run easily on a laptop.  Other's require the most powerful
supercomputers.

Using the model requires first deciding what configuration to use and creating
a case with that configuration.
The configuration options are managed by the Case Control System (CCS)
within the
[Community Infrastructure for Modeling the Earth](https://esmci.github.io/cime/versions/master/html/what_cime/index.html) (CIME).

Before reading the rest of this guide, you should become familiar with
cases, compsets and grids by reading the
[Case Control System Basic Usage](https://esmci.github.io/cime/versions/master/html/users_guide/index.html#case-control-system-part-1-basic-usage)

A [step-by-step guide](https://docs.e3sm.org/running-e3sm-guide/) for running E3SM with a run script
is available.

## Calendar Selection

`NO_LEAP` remains the default calendar. To opt a new case into the proleptic
Gregorian calendar, run `./xmlchange CALENDAR=GREGORIAN` in its case directory
before setup. Model days have 86400 seconds; UTC leap seconds are not modeled.
Select `./xmlchange CALENDAR=NO_LEAP` to use the default instead. For supported
CIME-managed cases, `CALENDAR` is the single calendar control: MPAS sea-ice,
ocean, and land-ice namelist defaults and the coupled WW3 runtime calendar
follow it. Coupled EAMxx
also follows it without a separate `SCREAM_HAS_LEAP_YEAR` setting. Native
standalone executables retain their own calendar configuration.
Do not reinterpret an existing NO_LEAP restart by changing its calendar.

Surface albedos may deliberately predict the next atmosphere shortwave event.
An albedo's valid time can therefore differ from a flux-history timestamp.
Conservation checks must pair current fluxes with the retained matching
properties, using actual averaging intervals and timestep area weights.
The opt-in [solar conservation suite](../dev-guide/index.md#solar-conservation-suite)
checks TOA/heating closure, actual imported-albedo consumption, surface
exchange and slower-ocean delivery inventories under both calendars,
including restart and scheduled-time consistency. These radiative partition
and transfer checks are not a full prognostic thermodynamic storage budget;
supported-machine coupled runtime qualification remains pending.
For the EAM/MCT scheduling timeline and diagnostic interpretation, see
[timing and diagnostic records](../dev-guide/index.md#timing-and-diagnostic-records);
the developer guide also provides
[failure triage](../dev-guide/index.md#local-validation-and-failure-triage).

NO_LEAP orbital phase retains the day-80.5 vernal-equinox reference.
Independently, `orb_mode = 'fixed_year'` anchors the reference to calendar day
80.5 each year under either calendar. Gregorian fixed-year phase advances
continuously between successive annual day-80.5 anchors, using the actual
365- or 366-day interval between them. Before the anchor, this is the previous
civil year's length; after it, this is the current year's length. Day 80.5 is
March 21 at noon in a common year and March 20 at noon in a Gregorian leap year.
Phase does not jump at January 1. Its rate can change at an anchor.

Gregorian variable-year and explicit-parameter modes retain continuous elapsed
model time relative to 2000-01-01 00:00, with an idealized 365.2422-day period
and the legacy January-1 phase calibration. Their phase does not reset at
civil-year boundaries. Fixed-year and explicit orbital elements remain fixed;
variable-year coefficients use the configured year alignment and fractional
civil year. Perpetual-season geometry uses the configured reference date.
These are Berger-compatible model conventions, not precision ephemerides;
the truncated orbital series need not give exactly zero declination at the
reference for nonzero eccentricity.

This fixed-year convention changes Gregorian results relative to both the
earlier elapsed-time convention and discontinuous civil-year scaling. Use fresh
initialization and new baselines. Gregorian fixed-year exact restarts require
the continuous day-80.5 phase marker; an absent or older convention is rejected.

### CMIP Compatibility

These changes address independent aspects of solar forcing and calendar
support:

- **Equation of time:** Correct solar hour angles to improve the seasonally
	varying timing of radiation, regardless of calendar or CMIP participation.
	This correction is intended to be consistent with CMIP7; it does not shift
	orbital seasonal phase.
- **Orbital year:** Correct the blanket 1990 default by deriving the orbital
	year from the compset, such as 1850. CMIP configurations must use the
	applicable protocol-prescribed year. This bug fix is independent of calendar
	choice and equation-of-time support.
- **Gregorian calendar:** Provide an optional proleptic Gregorian calendar for
	observational consistency in seasonal-to-decadal (S2D) simulations, while
	retaining CMIP orbital conventions and respecting experiment-specific
	calendar requirements.
- **Variable orbits:** Offer evolving orbital parameters as an explicit S2D
	option, not for CMIP-protocol forced simulations requiring fixed orbital
	settings. Solar System and Earth-Moon barycenter corrections are outside
	the present scope; they may be considered for additional precision in
	future non-CMIP applications.
- **Vernal equinox:** Preserve the calendar-day-80.5 reference for fixed-year
	orbits under either calendar as the intended CMIP convention. Do not change
	this reference when using fixed-year orbital settings.

These are intended compatibility choices, not certification of CMIP
compliance. Confirmation requires the applicable CMIP/CMIP7 experiment and
forcing guidance, together with model validation; local implementation checks
alone do not establish compliance.

### Orbital Year Selection

Orbital mode is independent of `CALENDAR`. Both calendars default to
`orb_mode = 'fixed_year'` in the coupled driver. Compsets whose time token
starts with a numeric year use that fixed orbital year, including transient
ranges and suffixed tokens such as `1850-2014` and `1950SOI`. Historical
`20TR` compsets use a fixed 1850 orbit. An ambiguous time token, such as
`SSP370`, requires an explicit `orb_iyear` or complete `fixed_parameters`
settings in `user_nl_cpl`; namelist generation otherwise fails. There is no
1990 fallback, and `RUN_STARTDATE` is not used to infer the forcing epoch.
The orbital year does not change merely because forcing varies or the
model clock advances.

Set explicit orbital options in `user_nl_cpl`; these override the defaults.
To retain the previous orbital choice, set `orb_mode = 'fixed_year'` and
`orb_iyear = 1990`. Changing the default orbital year changes scientific
results for both calendars; an unchanged NO_LEAP calendar alone does not
preserve the old 1990-orbit baseline.

To vary the orbit with the model year, set
`orb_mode = 'variable_orbit_year'`. The old spelling `variable_year` remains
a supported alias, not a separate logical switch. The orbital year is
`orb_iyear + model_year - orb_iyear_align`; both year settings default to the
compset orbital epoch above. For an ambiguous epoch, both must be set
explicitly. Set `orb_iyear_align` explicitly when the model clock and orbital
epoch differ. For example, `orb_iyear = 1850` and
`orb_iyear_align = 1` associate model year 0001 with the 1850 orbit.
`orb_mode = 'fixed_parameters'` continues to use explicit `orb_eccen`,
`orb_mvelp`, and `orb_obliq` values instead of an orbital year.

### Equation of Time

The equation-of-time correction defaults on for supported MCT and MOAB cases,
under both NO_LEAP and Gregorian calendars and every orbital mode. To retain
the mean-solar-time convention, explicitly disable it in `user_nl_cpl`:

```fortran
orb_use_equation_of_time = .false.
```

The correction is apparent solar time minus mean solar time, in seconds.
A positive value advances apparent solar time, so apparent noon occurs earlier
in mean time. It is derived from the existing orbital solution's mean longitude
and solar right ascension, not a modern-year empirical fit. It changes hour
angle, not the declination or Earth-Sun distance formulas. It does not shift
the model clock, forcing timestamps, forcing-year mapping, or coupling schedule.

This option is independent of `CALENDAR` and `orb_mode`. Fixed orbital years and
fixed elements still give a seasonally varying correction. Perpetual-season
calculations use the same reference date for declination and equation of time.
Fixed or uniform zenith-angle overrides retain precedence.

| Path | Enabled Behavior |
| --- | --- |
| EAM radiation and supplied photochemistry angles | Correct shared zenith geometry. |
| EAMxx RRTMGP and MAM photochemistry | Compute each process's correction at its own evaluation time. |
| ELM albedo, terrain zenith/azimuth, and FATES | Use the full scheduled next-shortwave timestamp. |
| Coupled CICE and MPAS-Seaice/Icepack | Use scheduled declination and equation of time from the driver. |
| DATM cosine-zenith interpolation | Recompute the correction for instantaneous and normalization samples. |
| Coupler ocean albedo | Use scheduled declination, EOT and the atmosphere's shortwave averaging interval. |
| MOSART stream heating | Consume ELM's corrected exported zenith angle. |

The coupled driver's equation-of-time setting is authoritative in EAMxx,
including when it is disabled. Leave
process-local `use_equation_of_time` unset to inherit it; an explicitly
contradictory process setting is rejected. Standalone EAMxx defaults on;
set `use_equation_of_time: false` in each participating RRTMGP/MAM process's
parameters to opt out. Keep those standalone process settings consistent. EAMxx's exact
coupled solar scheduling requires unsubcycled radiation and enclosing groups.

Computed ocean albedo uses the existing albedo formula evaluated at the same
forward interval-averaged cosine as atmospheric shortwave. EAM supplies its
actual averaging width, including zero when averaging is disabled; EAMxx supplies
radiation frequency times its atmosphere timestep. DATM supplies zero for its
instantaneous solar interpolation. Unknown intervals are rejected for computed
ocean albedos rather than silently treated as instantaneous.

A sunrise-crossing interval with a positive mean cosine therefore receives
daylight albedos and positive absorption even when its starting instant is dark.
This is the existing formula applied to the mean cosine, not a new
irradiance-weighted mean-albedo parameterization. Received downward shortwave is
not discarded because a saved albedo equals one; full reflection preserves
equal and opposite gross fluxes. Current fluxes retain their saved albedos, and
the later albedo update prepares only the next radiation interval.

Instantaneous solar geometry uses a common geometric horizon at cosine zero.
Shared cosine-interpolated forcing no longer discards small positive daylight
cosines below `0.001`; a genuinely all-dark normalization interval still gives
zero flux without division by zero. Coupled Icepack also uses scheduled time,
rather than its local clock, for NO_LEAP/EOT-off hour angles.

ELM coupler-bypass solar interpolation now uses driver-authoritative orbital,
calendar and EOT geometry. Its midpoint samples retain their full dates and
seconds across midnight, leap days and year changes, replacing the independent
365-day solar approximation. This does not adapt the bypass meteorological or
aerosol readers' fixed-calendar indexing. Prescribed prepartitioned fluxes and
the existing interpolation ratio cap retain their separate behavior.

These corrections do not force identical radiation masks: instantaneous geometry
and an interval's daylight average can legitimately have different signs near
sunrise or sunset. Radiation/chemistry averaging windows and scheduled land/ice
albedo sampling remain unchanged. The fixes deliberately change affected
sea-ice and forcing/bypass baselines even with EOT disabled; cross-version saved
forcing/albedo states are not guaranteed bit-for-bit continuation.

Computed-ocean exact restarts require interval-aware ocean-solar provenance
version 1 and a finite nonnegative saved averaging width. Missing or unsupported
provenance and changed averaging widths are rejected before saved albedos are
used. Older affected restarts require startup/hybrid initialization. No-ocean
and constant-albedo configurations do not require this ocean-specific marker.

Enabled legacy ESMF coupling and uncoupled Icepack solar geometry are rejected.
MPAS-Ocean's `config_sw_absorption_type = 'ohlmann00'` is also rejected when
enabled, because its prescribed `zenithAngle` data have no compatible convention
marker. The `jerlov` and `none` paths remain available. GOTM turbulence does not
call the separate standalone GOTM air-sea solar calculation. Prescribed solar
data are not re-timed or corrected by guessing their original convention.
Components without independent solar hour-angle calculations need no clock shift.
The enabled default does not remove these restrictions. Unsupported enabled
configurations require an explicit opt-out or a supported configuration.

Driver restarts record the selected convention. Exact restarts reject a changed
setting; a missing marker means legacy/off and is rejected when enabling the
correction. No numerical equation-of-time offset is restored or reused across
timestamps. Use startup/hybrid initialization for an intentional convention change.
To continue an old off-convention exact restart, explicitly set
`orb_use_equation_of_time = .false.` rather than adopting the new default.

The enabled default changes scientific results and requires new baselines,
including for NO_LEAP. Disabling EOT retains the mean-solar-time convention,
but does not disable the ocean interval correction or gross-flux preservation.
Those corrections deliberately change affected dawn/dusk absorption on both
calendars. Coupled EAMxx RRTMGP and MAM use shared driver-authoritative geometry
even with EOT disabled; RRTMGP also uses frequency-aware scheduling. These
corrections can change legacy NO_LEAP results, including MAM's former local-orbit
fallback with EOT off.
Full-model
bit-for-bit, component agreement, forcing-energy, and continuous-versus-restart
checks still require a supported machine; local compiled procedure checks are
not certification. Use identical driver orbital settings in comparisons.
EAM refreshes its saved orbital context at each coupling entry. In particular,
NO_LEAP variable-year/EOT-on runs no longer retain initialization-year orbital
parameters or acquire different parameters merely by restarting in a later year.

With the default input-calendar policy, DATM retains its existing forcing-date
selection, including repeated February 28 data on February 29 when reading
NO_LEAP forcing. The optional
[no-leap input compatibility policy](#using-no-leap-forcing-with-leap-years)
changes this selection for supported readers. Gregorian solar geometry uses
actual model dates instead. Cosine-zenith interpolation normalizes over model-step
samples across each contiguous model-time application of a forcing record, with
matching daylight cutoffs in the instantaneous and averaged geometry. Dark days
within an illuminated multi-day interval share that interval's normalization;
they do not independently discard part of the record's energy. A complete
interval without resolved daylight stays dark. Repeated forcing days do not
imply unchanged annual energy totals.

New driver restarts store the full next-shortwave date, seconds, and validity.
An old Gregorian restart with an ambiguous solar timestamp is rejected unless
the component scheduler can establish a unique replacement. There is no silent
fallback to the legacy Gregorian orbit.

Gregorian driver coupling supports `NCPL_BASE_PERIOD=hour` or `day`, but not
`year` or `decade`. Changing the base period also requires preserving the
component coupling intervals; it is not a calendar-only adjustment.
Coupled Gregorian EAMxx radiation requires unsubcycled RRTMGP and enclosing
atmosphere groups; unsupported schedules fail explicitly. Standalone Gregorian CICE and
uncoupled Gregorian Icepack solar calculations are not supported by this
interface. External CICE ESMF hosts must supply `nextsw_ymd`, `nextsw_tod`,
`nextsw_valid`, `nextsw_declination`, and `nextsw_eccf` attributes alongside
`nextsw_cday`, using the same scheduled instant and orbital context.

These Gregorian changes are experimental pending scientific calibration review
and supported-machine validation. Required checks include NO_LEAP bit-for-bit
comparisons, leap and century boundaries, restart equivalence, component solar
agreement, forcing energy normalization, and internal/external ESMF consistency.

## Using No-Leap Forcing With Leap Years

!!! warning "Experimental Input Compatibility, Not Calendar Conversion"
		`fudge_noleap_for_gregorian` is a special, default-off driver option for
		using no-leap forcing in a Gregorian simulation. It changes input time
		interpretation and scientific results. It does not conserve annual forcing
		totals, convert restart state, or make every input reader calendar-aware.
		Read the coverage table before enabling it. Full supported-machine builds,
		scientific comparisons and split-run restart validation remain required.

### Enable the Policy

For a new case, select the Gregorian calendar before setup:

```bash
./xmlchange CALENDAR=GREGORIAN
```

Set the following in the case's `user_nl_cpl`:

```fortran
fudge_noleap_for_gregorian = .true.
```

Generate the case namelists through the normal case setup/build workflow.
Check that `drv_in` contains the option in `seq_infodata_inparm`, then inspect
the coupler and affected component logs at initialization. Both MCT and MOAB
drivers recognize the option. It is not a new `xmlchange` variable, and coupled
components must not set independent, contradictory versions of the policy.
Enabling it with `CALENDAR=NO_LEAP` is an error. Standalone component input
configuration is not enabled by this driver option.

The default is `.false.`. Leaving the option unset preserves the existing input
handling, including the shared-stream February 28 fallback. A policy change is
not permitted on an exact continuation. To make an intentional change, use the
component-supported fresh startup or hybrid initialization workflow rather
than relabeling an existing restart.

### How Input Dates Are Mapped

The reader first decodes each numeric offset using the **source calendar**.
It then places that civil month, day and time on the Gregorian model timeline.
Existing year alignment and cycling still select the corresponding model year
where supported. The source file is not modified and its calendar attribute
is not replaced with the model calendar.

For example, `59 days since 2000-01-01 00:00:00` in a `noleap` file denotes
March 1, not February 29. The mapped record remains March 1. Adding that numeric
offset directly to a Gregorian clock would assign the wrong date.

For linear point samples, interpolation uses elapsed Gregorian seconds between
the mapped records. Consider daily source records with these values:

| Gregorian Model Time | Value | Origin |
| --- | --- | --- |
| February 28, 00:00 | 10 | Source record. |
| February 29, 00:00 | 12 | Equal-weight interpolation. |
| February 29, 12:00 | 13 | Interpolation, three quarters of the interval. |
| March 1, 00:00 | 14 | Source record. |

This inserts a gap in the input timeline; it does not duplicate February 28 or
stretch every date in the year. Gregorian century rules apply: 2000 and 2400
are leap years, while 1900 and 2100 are not. Outside leap-spanning intervals,
record dates and ordinary interpolation remain unchanged.

For subdaily data, **all original record times are retained**. A final record
at February 28, 18:00 and the next at March 1, 00:00 are 30 Gregorian hours
apart in a leap year. February 29, 00:00 is one fifth of the way through that
interval. This is not a same-hour blend of the two neighboring days and does
not reconstruct a missing diurnal cycle. Monthly inputs likewise interpolate
between their actual mapped record dates, with the longer Gregorian interval
where it contains leap day.

### Input Metadata and Quantities

- Supply explicit calendar metadata on each input time coordinate. Recognized
	no-leap names include `noleap`, `no_leap`, `365_day` and `365day`. Gregorian
	input is not forcibly relabeled or adapted. Missing or empty calendar
	metadata requires an explicit override with the option enabled, as described
	below. The default-off legacy fallback is not used to guess an enabled
	input's calendar. `standard` and `gregorian` share the mixed-calendar flavor;
	`proleptic_gregorian` retains a distinct source flavor. Internal whitespace
	does not turn malformed tokens into valid aliases.
- Keep the calendar consistent across a file sequence. Shared streams and the
	EAMxx file interpolator reject mixed-calendar sequences. Chemistry validates
	individual time-coordinate objects; this does not certify mixed-calendar
	sequences managed by their callers. Unsupported calendars, such as
	`360_day`, are not covered by this option.
- Keep the reader's existing epoch and time-unit requirements. In particular,
	EAM chemistry's numeric coordinate path supports its existing `days since`
	format. The EAMxx file interpolator retains its existing `case_t0` epoch and
	supported units, and requires finite, nonnegative source offsets when
	adapting no-leap files. This is not a general CF time-reader replacement.
- State quantities and instantaneous rates can be linearly interpolated as
	point samples. Missing `cell_methods` means the point-sample assumption is
	reported for the EAMxx file interpolator; declared interval data are not
	silently reclassified as points.
- Declared `time: sum`, unsupported time reductions, and repeated time methods
	require a quantity-aware adapter and are rejected by the adapted readers.
	Do not remove metadata to bypass these checks. Accumulations cannot be
	made safe simply by interpolating their numbers.
- EAM chemistry maps both endpoints of time bounds in the source calendar
	before using them in Gregorian time. Declared means require usable bounds.
	The existing interval-mean selection behavior is retained rather than
	silently replaced with point interpolation.
- Shared streams reject declared time means and ambiguous bounded data until
	a bounds-aware adapter is available. Cosine-zenith forcing retains its
	existing solar normalization over mapped record intervals, but that does
	not establish support for arbitrary CF interval means. Do not remove
	`time: mean` or time bounds to bypass this restriction.

### Identifying an Input Calendar

Identify the calendar of each input coordinate and the reader that consumes it,
not just the receiving component. The model's `CALENDAR` setting does not prove
the source calendar. These checks apply to both compatibility directions.

For CF-style NetCDF input, inspect the relevant time coordinate's `calendar`
and `units`, including its reference epoch. For example:

```text
time:units = "days since 1850-01-01"
time:calendar = "noleap"
```

| Calendar Declaration | Meaning |
| --- | --- |
| `noleap`, `365_day` | Every year has 365 days. The resolver also accepts `no_leap` and `365day`. |
| `proleptic_gregorian` | Gregorian leap-year rules throughout the timeline. |
| `standard`, `gregorian` | Mixed Julian/Gregorian calendar under CF conventions. |
| `360_day`, `all_leap`, `366_day`, `julian` | Other calendars, unsupported by these compatibility options. |

The resolver retains both the `GREGORIAN` model-calendar category and the
scientifically distinct source flavor. Shared streams, EAM chemistry
`time_coordinate`, EAMxx `TimeInterpolation` and preflight reject mixed
`standard`/`gregorian` inputs whose epoch, records or time bounds extend before
1582-10-15. This includes negative offsets crossing the cutoff and early epochs
whose records describe modern dates. The runtime gate also applies with both
compatibility options off; the legacy missing-metadata fallback remains unchanged.
Explicit `proleptic_gregorian` inputs retain historical support within the
reader's existing date/unit limits. A proleptic override cannot relabel contrary
mixed-calendar metadata, and the override restart signature preserves this
distinction. Full mixed Julian/Gregorian decoding is not implemented. See the
[CF calendar conventions](https://cfconventions.org/cf-conventions/cf-conventions.html#calendar).

Inspect the coordinate's `bounds` or `climatology` attribute and the forcing
variable's `cell_methods` as well. Instantaneous samples, interval means,
accumulations and climatologies have different temporal semantics. Calendar
identification alone does not establish that a reader supports those semantics.

**CMIP and CMIP6 are not calendar declarations.** Check model-produced CMIP
datasets and prescribed CMIP forcing products independently, including every
file in a sequence. Twelve monthly records, a filename, an experiment name,
or an absent February 29 do not establish no-leap: a Gregorian common year
can have the same records. Check provider documentation and preprocessing
history when metadata is missing, ambiguous or potentially mislabeled.

| Input Path | Calendar Evidence and Current Limits |
| --- | --- |
| DATM meteorology, DOCN ocean forcing and DICE ice forcing through shared streams | Read `calendar` on the configured time coordinate. Stream year alignment and cycling are separate from the source declaration. |
| MPAS prescribed sea ice through shared streams | Resolve source metadata separately from the model's `config_calendar_type`; this is not the native MPAS forcing-group path. |
| Native MPAS ocean and sea-ice forcing groups | Use native clocks, configured start times, intervals and cycling. Active groups reject either enabled compatibility policy; they are guarded, not adapted. |
| EAM aircraft emissions using `time_coordinate` | Read `time:calendar` through the adapted chemistry reader, subject to its units, date and quantity restrictions. |
| EAM prescribed aerosols using `tracer_data` | Reject time-varying input under either enabled policy. Fixed inputs are not evidence of time-varying reader support. |
| EAM `VOLC_CMIP6` volcanic aerosols | A separate reader converts integer `date` values through model-time routines without the strict source-calendar resolver. The CMIP6 label does not establish compatibility. |
| EAMxx SPA aerosols | Use `DataInterpolation`, not the adapted `TimeInterpolation`. The inspected path constructs model timestamps from reference dates and numeric day offsets without the strict resolver. |
| ELM nitrogen/phosphorus deposition and other shared streams | Use shared-stream source-calendar resolution; this does not establish coverage for independent land input paths. |
| ELM `CPL_BYPASS` meteorology and aerosol deposition | Meteorology uses fixed 365-day indexing; aerosol deposition uses positional monthly records and a fixed month table. Those reader assumptions do not establish the files' actual calendars. |

The CMIP6 volcanic, EAMxx SPA and ELM bypass paths above require separate
qualification. They must not be presented as covered by the adapted-reader
tests or by the enabled-policy startup banner. This source audit is not an
inspection or certification of actual forcing datasets.

Coupler-delivered fields have already been evaluated upstream; the receiving
component does not re-detect the original source-file calendar. Static meshes,
maps and parameter data generally have no forcing-time interpolation to adapt.

For ordinary CF coordinates without an explicit custom calendar, an omitted
`calendar` implies `standard`. That convention differs from the shared-stream
legacy fallback to NO_LEAP when both compatibility options are off. With either
option enabled, adapted readers deliberately require supported explicit metadata
or an exact asserted override. An override records external knowledge; it does
not discover a calendar, convert a dataset or adapt an unsupported reader.

### Explicit Metadata Overrides

For legacy files with absent or empty calendar metadata, supply paired arrays
in `user_nl_cpl`:

```fortran
fudge_noleap_for_gregorian = .true.
input_calendar_override_selectors = 'file:/input/forcing.nc', 'stream:forcing.txt'
input_calendar_override_calendars = 'noleap', '365_day'
```

Both arrays default to empty; at most 128 rules are supported. Nonempty rules
require the experimental option. These are explicit assertions about the
source data, not inferred calendars, and do not rewrite the files. A rule may
fill missing metadata or confirm an equivalent declaration. It cannot replace
an unsupported declaration or contradict an existing calendar. Duplicate
selectors and disagreeing file/stream rules are errors.

Matching is literal and case-sensitive for identifiers, without wildcards,
basename matching, path resolution or symlink normalization. Copy the file
identifier printed by the reader's `NOLEAP_INPUT` diagnostics. Shared streams
use the path/name opened by the reader; a `stream:` rule uses the literal
stream-description identifier in the generated component namelist, commonly
a relative filename. Chemistry uses its original requested input filepath,
before local staging. EAMxx uses its requested filename and supports `file:`
rules, not `stream:` rules. Programmatically constructed shared streams without
a description identifier require file rules. Calendars are case-insensitive
and accept the aliases above.

Enabled readers report the declaration, canonical resolved calendar and origin:
`metadata`, `override` or `metadata+override`. An explicit override produces an
additional warning that no calendar was inferred. The rule registry is immutable
after initialization. Exact continuation requires the same canonical rule set;
reordering entries or using equivalent calendar aliases is allowed. Restart
files retain the canonical sorted rules. A legacy restart without this record
can continue only with an empty rule set, subject to the existing calendar and
enabled-policy restart checks.

Neither metadata nor an override proves that a mislabeled dataset was actually
produced with that calendar. Coordinate counts and an absent February 29 do not
establish source intent. Enabled readers validate finite, strictly increasing
source coordinates and source-calendar dates before mapping. Reader-specific
units, epoch, quantity and native-reader restrictions still apply.

### Input Calendar Preflight

Run the project-owned command from the run directory before submitting a case:

```bash
python /path/to/E3SM/cime_config/scripts/input_calendar_preflight.py \
	--policy-namelist drv_in --namelist datm_in \
	--input chemistry /input/chemistry.nc
```

Repeat `--namelist` for generated component namelists that contain shared
`streams` arrays, `--stream` for an exact stream-description identifier, and
`--input READER FILE` for explicitly identified native inputs. Adapted reader
labels are `stream`, `chemistry` and `eamxx`. Active unsupported paths can be
inventoried as `mpas`, `ww3`, `flbc`, `tracer` or `eamxx-iop`; they report an
error rather than silently claiming support. To read `user_nl_cpl` instead of
generated `drv_in`, add `--groupless`. Use `--rundir` when running elsewhere.
The policy must be enabled. CIME's actual parser reads namelists; the optional
`netCDF4` and NumPy dependencies are needed only when auditing inputs.

The JSON report lists file/stream identity, reader, declared and resolved
calendar, provenance, matched rules, coordinate span and errors. It checks
finite and strictly increasing numeric time, CF units/epoch decoding, explicit
date/datesec consistency where present, shared-stream calendar consistency,
and unsupported reader/quantity declarations. Missing units, unresolved remote
staging or unexpanded paths require explicit local inventory; no dataset is
downloaded or guessed. Unmatched rules are errors so misspelled overrides are
visible. Errors produce a nonzero exit status.

**The inventory is always explicitly marked partial.** Only supplied inputs and
configured streams are audited. Exit zero means those audited checks passed,
not that every native input was discovered or that a model run is qualified.
The preflight does not certify interpolation algorithms, cycling, spatial
transforms, conservation or restart behavior. It never substitutes for the
runtime guards or dependency-backed tests below.

### Reader Coverage

Coverage is by **input path**, not by component name. An eligible reader must
still satisfy its metadata and algorithm requirements.

| Input Path | Enabled Behavior |
| --- | --- |
| Shared streams with `linear` interpolation | Keep source date labels and compute Gregorian weights for point samples. Includes eligible data-component streams, ELM streams and MPAS-Seaice prescribed shared streams. Declared means and ambiguous bounded data are rejected. |
| Shared streams with `coszen` interpolation | Keep mapped record intervals and Gregorian solar normalization; do not repeat February 28. Declared means and ambiguous bounded data still require an adapter. |
| Shared streams with `lower`, `upper`, `nearest` or other algorithms | Reject no-leap input rather than silently substitute a different algorithm. |
| EAMxx file-backed `TimeInterpolation` | Decode no-leap offsets before forming file timestamps; interpolate point samples across mapped gaps. Declared time means, sums or ambiguous bounds require another adapter. |
| EAM chemistry `time_coordinate` | Map numeric source timestamps and usable bounds before model-time lookup; reject unsupported quantity semantics. |
| EAM legacy FLBC and `tracer_data` | Reject active time-varying input under this policy. Fixed input is not temporally interpolated. |
| EAMxx IOP input | Reject when the policy is enabled; this is a separate, unadapted reader. |
| Native MPAS forcing groups | Reject when a group is initialized under the enabled policy. Coupled ocean, sea ice and land ice pass the policy to the framework; standalone MPAS is unchanged. |
| Coupled WW3 | Reject the enabled policy until its native input path has a calendar adapter. |
| Coupler-exported fields and static meshes, maps and parameters | No input-calendar interpolation is applied by this option. |

Additional specialized readers are not certified by this table. In particular,
sharing a component with an adapted reader does not establish coverage for its
nudging, initialization or other independent input paths. Do not treat the
startup warning banner as proof of universal reader support.

If an algorithm or reader is rejected, use a supported input path or prepare
compatible Gregorian data and disable the option. Do not change a stream's
algorithm merely to suppress an error: rates, interval means and accumulated
totals have different semantics.

### Initial Conditions and Restarts

The model clock, timestep, coupling schedule, output calendar, orbital policy
and solar-geometry convention remain Gregorian. No leap seconds are added.
Static inputs are not re-timed, and a startup state snapshot is not interpolated
like a forcing time series.

Use only a component-supported startup import that initializes a new Gregorian
clock and rebuilds its scheduling state. ELM's startup route initializes its
time manager separately from continuation restart loading. This does not certify
arbitrary state imports for other components. Choose a physically appropriate
snapshot explicitly; the option does not invent a February 29 initial state
from no-leap restart files.

New driver restarts record the policy and model calendar. Exact continuation
requires the same policy. A missing policy marker means legacy/off. An enabled
restart also requires a Gregorian calendar marker. ELM checks its component
restart calendar before allowing it to replace its time-manager calendar.
The same rules apply to driver branch runs that restore restart state;
branching alone is not a cross-calendar conversion.

A Gregorian run using supported no-leap forcing may write and read normal
Gregorian restarts with the same option. This is different from continuing a
NO_LEAP simulation as Gregorian. Restart equivalence across February 29 still
requires supported-machine validation, including reconstructed forcing bounds,
record selection, solar caches and component scheduling.

### Warning Output and Troubleshooting

Search the coupler and component logs for `NOLEAP_INPUT`. Warnings are emitted
on the relevant reader's logging rank, not independently by every MPI rank.
They are deliberately visible even when ordinary informational output is
reduced. They are event-based rather than repeated every timestep.

For example, the driver emits:

```text
WARNING NOLEAP_INPUT: fudge_noleap_for_gregorian is ENABLED
WARNING NOLEAP_INPUT: preserve source dates; interpolate leap-day gaps
WARNING NOLEAP_INPUT: annual forcing totals are NOT conserved
WARNING NOLEAP_INPUT: this is not a cross-calendar restart converter
```

| Event | Expected Diagnostic |
| --- | --- |
| Driver initialization | Enabled-policy banner, date-preservation policy, non-conservation warning and restart-conversion restriction. |
| Source inspection | File or stream identity, source/model calendars and the interpretation applied. |
| Leap-spanning interpolation interval | Shared-stream bounding records; EAMxx point-interpolation bounds and weights at the reported sample. |
| February 29 entry and exit | A warning identifying the affected stream or chemistry file; EAMxx includes the model timestamp. |
| Reader cleanup | A summary of adapted inputs or leap-day events for that execution segment, where the reader's cleanup/finalize method is invoked. |
| Restart initialization | Validation of the policy/calendar markers and, for ELM, its component restart calendar. |

For a missing-calendar error, identify the true source calendar from dataset
provenance and obtain correctly described data. For an unsupported quantity,
reader or algorithm error, use a supported adapter or compatible forcing, not
a metadata override. For a restart-policy error, retain the original setting
for exact continuation or use a documented fresh initialization workflow.

Source warnings may recur when files are reopened or a restart segment begins.
Segment summary counts are not accumulated across jobs. A run that never
reaches a leap-spanning interval can still report adapted input files, and an
aborted run may not reach reader cleanup.

### Scientific and Validation Limits

Adding a day of interpolated flux does not preserve the source annual integral.
Precipitation, radiation, emissions, deposition and freshwater or heat fluxes
can therefore change annual budgets. The option provides calendar compatibility,
not a conservative remapper or a physically reconstructed missing day.

Before production use, check daily, subdaily and monthly inputs across February
28 through March 1, year alignment/cycling, source-file boundaries and leap versus
non-leap centuries. Compare uninterrupted and restarted runs with identical
policy settings. Verify forcing budgets and scientific suitability for the
experiment. Source parsing, namelist checks and ordinary calendar tests alone
do not establish model bit-for-bit behavior or scientific certification.

The [no-leap compatibility test matrix](../dev-guide/index.md#no-leap-input-compatibility-tests)
documents the opt-in smoke, short, negative and cycle suites, numerical reference
values, warning checks and outstanding reader/guard gates. Local fixture and
harness tests are not full reader qualification: coupled runs and the new
dependency-backed tests still need execution, and some native-reader/guard
harnesses are not yet implemented. Do not interpret these test definitions as
support for every input dataset.

## Using Gregorian Forcing With No-Leap Models

!!! warning "Experimental Input Compatibility"

	`fudge_gregorian_to_noleap` is default-off and requires `CALENDAR=NO_LEAP`.
	It retains Gregorian source interpolation, including February 29 records.
	It does not convert model clocks or restart state, conserve annual forcing
	totals, or qualify every input reader. Dependency-backed reader tests and
	supported-machine integration runs remain required.

### Enable the Reverse Policy

For a new case, select the model calendar before setup:

```bash
./xmlchange CALENDAR=NO_LEAP
```

Set these driver options in `user_nl_cpl`:

```fortran
fudge_gregorian_to_noleap = .true.
fudge_noleap_for_gregorian = .false.
```

Both MCT and MOAB recognize the option in `seq_infodata_inparm`. Generate
namelists normally and verify the resulting `drv_in`. The two options are
mutually exclusive; enabling both is an error. Neither option is an
`xmlchange` variable. The correctly spelled reverse option has no
`gregorial` alias. Unset or explicit `.false.` preserves legacy behavior.

| Policy | Model Calendar | Mismatched Source | Interpolation |
| --- | --- | --- | --- |
| `fudge_noleap_for_gregorian` | Gregorian | NO_LEAP | Preserve source civil labels and interpolate Gregorian leap-day gaps. |
| `fudge_gregorian_to_noleap` | NO_LEAP | Gregorian | Preserve Gregorian source intervals and February 29 brackets; skip only the model visit. |

The reverse policy is not the mathematical inverse of the forward policy.
Compatible NO_LEAP forcing stays on its ordinary NO_LEAP input path.

### February 29 Remains a Source Record

The model advances from February 28 to March 1 without a February 29 model
timestamp. Input decoding and weights nevertheless use the Gregorian source
calendar. A leap-day source record can therefore influence a February 28
query. It is not dropped, compressed out, copied onto another date, or
redistributed among neighboring days.

For Gregorian source year 2000, consider:

| Source Record | Value |
| --- | --- |
| February 28, 00:00 | 10 |
| February 29, 00:00 | 100 |
| March 1, 00:00 | 14 |

The NO_LEAP model query February 28, 12:00 returns **55**, and March 1,
00:00 returns **14**. There is no model/history sample on February 29.
If the leap-day record is absent and only the February 28 and March 1
points remain, February 28, 12:00 returns **11**, not the compressed-calendar
value 12: twelve hours have elapsed in a 48-hour Gregorian bracket.

For subdaily points February 28, 18:00 = 10 and February 29, 00:00 = 14,
February 28, 21:00 returns 12. If the second point is instead March 1,
00:00 = 14, the same query returns **10.4**, using three hours of a 30-hour
bracket. Missing February 29 records are not automatically missing coverage.
Monthly points likewise retain Gregorian interval lengths. Gregorian source
years 2000 and 2400 are leap years; 1900 and 2100 are not.

### Metadata, Readers and Cycling

First identify each source using
[Identifying an Input Calendar](#identifying-an-input-calendar), including its
reader-specific gaps. Those limits apply in both compatibility directions.

Supply truthful, explicit source calendars. The aliases and exact
[metadata override rules](#explicit-metadata-overrides) apply in both
directions. For a legacy Gregorian file with missing metadata, for example:

```fortran
fudge_gregorian_to_noleap = .true.
input_calendar_override_selectors = 'file:/input/forcing.nc'
input_calendar_override_calendars = 'gregorian'
```

An override cannot contradict declared metadata or rescue an unsupported
calendar. In particular, do not relabel Gregorian input as NO_LEAP.
Finite, strictly increasing source coordinates and valid source dates are
required. Numeric Gregorian offsets 59 and 60 from January 1, 2000 denote
distinct February 29 and March 1 records. Gregorian leap-day source epochs
are valid even though they are not valid NO_LEAP model timestamps.

| Input Path | Reverse Behavior |
| --- | --- |
| Shared `linear` streams | Gregorian source-time point interpolation. Year alignment/cycling retains valid source-year queries and brackets instead of constructing February 29 in a common model year. |
| Shared `coszen` streams | Existing normalization with NO_LEAP model solar time; actual-reader/coszen qualification remains pending. |
| EAM chemistry `time_coordinate` | Explicit Gregorian source-time context for numeric/date coordinates, bounds and current/fixed model queries; existing input-time shift and bounded-mean selection semantics remain. |
| EAMxx file-backed `TimeInterpolation` | Internal Gregorian-seconds input axis for ordering, coverage, bracket lookup and weights; the model `TimeStamp` remains NO_LEAP. Existing `case_t0`, supported units and nonnegative offset requirements apply. |
| Shared `nearest`, `lower`, `upper` and other unsupported algorithms | Reject adapted Gregorian streams rather than substitute an algorithm. |
| Time-varying legacy EAM FLBC/tracer, EAMxx IOP, native MPAS forcing and coupled WW3 | Reject either enabled policy; fixed/inactive paths retain their existing restrictions. |

The source-calendar view is separate from the model clock. Shared-stream
cycling selects the configured source year and preserves its Gregorian
intervals; cycle-edge positions use the configured source cycle's elapsed
duration. It does not manufacture February 29 in a common source year.
This does not add cycling to readers that lack it.

Shared streams and EAMxx retain their point-sample quantity restrictions:
declared means, sums, repeated/unsupported time methods and ambiguous bounds
are rejected where no quantity-aware adapter exists. Chemistry retains its
existing bounded-mean path with valid bounds, not annual-total redistribution.
Static meshes, parameters and coupling fields receive no time interpolation.
Other specialized native readers are not certified by this table.

### Preflight and Diagnostics

From the run directory, audit explicitly configured inputs:

```bash
python /path/to/E3SM/cime_config/scripts/input_calendar_preflight.py \
	--policy-namelist drv_in --model-calendar NO_LEAP --namelist datm_in \
	--input chemistry /input/chemistry.nc
```

For `user_nl_cpl`, add `--groupless` and supply `--model-calendar NO_LEAP`.
A generated calendar declaration is also accepted and must agree with the
policy and explicit context. The JSON reports direction, model calendar,
resolved source/interpolation calendars and leap-day record counts. Its
inventory remains **partial**; exit zero does not certify undiscovered inputs,
interpolation, cycling, restarts or scientific suitability.

Look for `GREGORIAN_INPUT` activation, source/model identity, retained-bracket
and non-conservation warnings. Shared metadata/rejection diagnostics can
still contain the historical `NOLEAP_INPUT` prefix. Reverse mode does not
announce entry into February 29 model time, because that model time does
not exist. Files may be reported again when reopened or after restart.
MPI-root ownership and warning-frequency qualification remain pending.

### Restarts and Scientific Limits

Model clocks, outputs, coupling schedules and solar/orbital conventions
remain NO_LEAP. Four model years contain 1460 days. Neither policy converts
initial-condition or restart state. Use only supported startup imports.

Exact continuation requires identical policy direction, NO_LEAP model
calendar and canonical override rules. New driver restarts store a separate
reverse-policy marker. A missing marker means reverse disabled, so an
enabled reverse continuation rejects an old/unmarked restart. Invalid
markers, changed flags, changed rules or incompatible component calendars
are errors. ELM checks stored calendar state before restoring it. Ordinary
same-policy NO_LEAP restart use still needs full split-run validation.

Skipping a model visit does not conserve annual flux integrals. Precipitation,
radiation, emissions and other budgets can change, and a discontinuity can
occur when the model jumps to March 1. Check the experiment's scientific
suitability rather than treating this option as a conservative remapper.
The [reverse test matrix](../dev-guide/index.md#gregorian-input-for-no-leap-tests)
lists independent reference values, opt-in suites and outstanding gates.

## Supported Coupled Compsets

A *fully coupled* compset is one which has active components for at least the atmosphere, ocean, land surface, ocean and
sea-ice all interacting.  Each compset is associated with a specific forcing condition.
Coupled compsets in E3SM are developed  for three science-driven simulation campaigns:  **water cycle change and impacts**, **human-earth system feedbacks**, and **polar processes, sea-level rise and coastal impacts**. The standard coupled configurations -- which consist of prognostic atmosphere, land, river, ocean and sea-ice components -- form the base physical coupled system and are mainly designed for `water cycle change and impacts` simulation campaign.
Below list the standard configuration compsets supported in the current version of E3SM:

|Compset alias | Description |
|:-----------  |:----------- |
|`WCYCL1850` | Standard configuration with pre-industrial climatological forcings |
|`WCYCL1850-4xCO2` | Same as `WCYCL1850` except with abrupt (then persistent) 4xCO2 forcing. |
|`WCYCL1850-1pctCO2` | Same as `WCYCL1850` except with 1 percent per year increase of CO2 concentration |
|`WCYCL1950` | Standard configuration with perpetual 1950 forcings |
|`WCYCL20TR` | Standard configuration with prescribed transient forcings over the historical period (1850-2014) |
|`WCYCLSSP245` | Standard configuration with prescribed SSP-245 forcings |
|`WCYCLSSP370` | Standard configuration with prescribed SSP-370 forcings |
|`WCYCLSSP585` | Standard configuration with prescribed SSP-585 forcings |
| | **Compsets for single-forcings simulations of the historical period (1850-2014)** |
|`WCYCL20TR-GHG` | Greenhouse gases only configuration (`GHGs`)|
|`WCYCL20TR-aer` | Anthropogenic aerosols and precursors only configuration (`aer`)|
|`WCYCL20TR-xGHG-xaer` | All forcings except GHGs and aer (`solar irradiance, stratospheric ozone and volcanic emissions, land use land cover`) |
|`WCYCL20TR-xaer` | All forcings except aer (`GHGs, solar irradiance, stratospheric ozone and volcanic emissions, land use land cover`) |
|`WCYCL20TR-nat` | Natural-only configuration (`solar irradiance, stratospheric volcanic emissions, land use land cover`) |
|`WCYCL20TR-ozone` | Stratospheric ozone only configuration |
|`WCYCL20TR-lulc` | Land use land cover only configuration |
|`WCYCL20TR-volc` | Stratospheric volcanic emissions only configuration |

The table below specifies which forcing category is fixed at 1850 conditions and which are allowed to vary over the historical period
for each of the historical ("20TR") compsets including the single-forcing compsets.

|Compset alias  |  GHGs    | Aerosol & precursors | Oxidants | Ozone (CI & Linoz) | Volcano | Solar   |  Land Use & ndep/popdensa |
|:------------  |:-----:   | :---:                | :---:    | :---:              | :---:   | :---:   | :---:   |
|`WCYCL20TR`    | varying       | varying         | varying  | varying            | varying | varying | varying |
|`WCYCL20TR-GHG`| varying       | 1850            | 1850     | 1850               | 1850    | 1850    | 1850    |
|`WCYCL20TR-aer`| 1850          | varying         | varying  | 1850               | 1850    | 1850    | 1850    |
|`WCYCL20TR-xGHG-xaer`| 1850    | 1850            | 1850     | varying            | varying | varying | varying |
|`WCYCL20TR-xaer`| varying      | 1850            | 1850     | varying            | varying | varying | varying |
|`WCYCL20TR-nat`| 1850          | 1850            | 1850     | 1850               | varying | varying | 1850    |
|`WCYCL20TR-ozone`| 1850        | 1850            | 1850     | varying            | 1850    | 1850    | 1850    |
|`WCYCL20TR-lulc`| 1850         | 1850            | 1850     | 1850               | 1850    | 1850    | varying |
|`WCYCL20TR-volc`| 1850         | 1850            | 1850     | 1850               | varying | 1850    | 1850    |

<!-- markdownlint-disable-next-line MD033 -->
- <sub> *Volcano* refers to stratospheric volcanic SO2 emissions; *1850* for *Volcano* refers to background (average) stratospheric volcanic emissions used in pre-industrial control experiments</sub>
<!-- markdownlint-disable-next-line MD033 -->
- <sub> *Oxidants* always follow *Aerosol & precursors* for fixed or varying.</sub>

The compsets for the other two science simulation campaigns are being finalized, with additional components and/or features.
The compset naming follows the same convention, e.g., `CRYO1850` and `CRYO1850-4xCO2` are with prognostic ice-shelf melt fluxes for the `polar processes` simulation campaign.

Compsets are also available for standalone component model configurations, See the User Guides for the components for more information.

## Supported resolution

Currently two grid sets are supported for the above compsets, including a nominal low-resosluton confiuguration and one regionally refined mesh. Additional regionally refined meshes and a high-resolution grid will become available in the near future.

| Grid Alias  |  Description  |
| ----------- |  ------------ |
|ne30pg2_r05_IcoswISC30E3r5 | For this grid set, the atmosphere is on the ne30pg2 cubed-sphere mesh with approximately 100km resolution, the land and river are on a 0.5deg x 0.5deg structured grid, and the ocean and sea ice are on a hexagonal mesh dervied from the dual of a 30km resolution icosahedral mesh with ice shelf cavities (wISC) around Antarctica.|
|northamericax4v1pg2_r025_IcoswISC30E3r5 | The atmosphere for this grid set uses North America regionally refined mesh from about 110 km down to 25 km over the refined region. The land and river are on 0.25deg x 0.25deg structured grid. The ocean and sea ice are on the same icosahedral mesh as for `ne30pg2_r05_IcoswISC30E3r5`.|

## Input data

Inputdata for coupled compsets at component model levels are the same as for the standalone component configurations
for a given forcing scenario (e.g., `1850` for the pre-industrial period,  `20TR` for the historical period, `2010`
for present-day condition, and `SSPs` for Shared Socioeconomic Pathways).
Between the coupled compsets, the differences are in the prescribed solar forcing, volcanic emissions,
atmospheric forcing data, and land use and land cover. The required inputdata for the pre-industrial and the historical periods
as well as the present-day condition are described in [the EAM User's Guide](https://e3sm-project.github.io/E3SM/EAM) and
[the ELM's User's Guide](https://e3sm-project.github.io/E3SM/ELM). Below are the prescribed forcing data for the SSP scenarios.

### [SSP245 forcing data](ssp245-forcings.md)

### [SSP370 forcing data](ssp370-forcings.md)

### [SSP585 forcing data](ssp585-forcings.md)
