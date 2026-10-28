# Development Guide

Development of E3SM is carefully planned and tied to goals in the proposals that fund
work on E3SM.

## Code development

Most code development occurs within one component and developers should consult the
relevant component's development guide.

More information on code development.

+ [Development Big Picture](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/7997024/Development+Big+Picture)
+ [Getting Started](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/1868455/Development+Getting+Started+Guide)
+ [Development Reference](https://acme-climate.atlassian.net/wiki/spaces/DOC/pages/2523163/Development+Reference)

## Creating new input

E3SM's public input data directory contains initial and boundary conditions for many cases.
Creating new initial or boundary conditions typically starts with creating a new grid.

+ [Adding Support for New Grids](adding-grid-support/index.md)

## Calendar and solar integration tests

The following opt-in CIME suites exercise calendar, orbital-mode, and
equation-of-time (EOT) combinations on economical grids. They are not inherited
by the standard developer suites.

| Suite | Tests | Coverage |
| --- | ---: | --- |
| `e3sm_calendar_smoke` | 16 | Gregorian fixed/EOT-on leap day and variable/EOT-on year rollover |
| `e3sm_calendar_short` | 128 | Both short boundaries for all eight profiles |
| `e3sm_calendar_cycle` | 64 | Four-year restart runs for all eight profiles |

On a supported machine, run from `cime/scripts`:

```bash
./create_test e3sm_calendar_smoke
./create_test e3sm_calendar_short
./create_test e3sm_calendar_cycle
```

Run the smoke suite first. The full short suite is intended for periodic
integration testing; the four-year suite is an expensive qualification test.
Use machine-appropriate processor layouts and walltimes. These suites do not
force build sharing across different component configurations.

### Component configurations

Each suite includes these eight configurations:

| Family | Compset | Grid | Main integration coverage |
| --- | --- | --- | --- |
| A | `A` | `ne4pg2_oQU480_rx1` | Driver and prescribed atmosphere, ocean, and ice; stub land |
| B | `WCYCL1850NS` | `ne4pg2_oQU480` | Fully coupled atmosphere, land, ocean, and ice |
| C | `CMPASO-NYF` | `T62_oQU120` | Ocean with prescribed forcing |
| G | `GMPAS-NYF` | `T62_oQU240` | Ocean and ice with prescribed forcing |
| D | `DTESTM` | `T62_oQU240` | Sea ice with prescribed atmosphere and slab ocean |
| I | `I1850ELM` | `f19_g16` | Land with prescribed atmospheric forcing |
| EAM F | `F1850` | `ne4pg2_oQU480` | EAM radiation and surface coupling |
| EAMxx F | `F2010-SCREAMv1` | `ne4pg2_ne4pg2` | EAMxx radiation and surface coupling |

These coarse grids reduce integration-test cost; they do not replace
production-resolution coupling and decomposition tests. Resource allocations
should match the small grids rather than reuse production node counts.

### Profiles and boundaries

Profile names encode calendar (`n`: NO_LEAP, `g`: Gregorian), orbital mode
(`f`: fixed, `v`: variable), and EOT (`0`: off, `1`: on).

| Profile | Calendar | Orbital mode | EOT |
| --- | --- | --- | --- |
| `nf0` | NO_LEAP | Fixed year | Off |
| `nf1` | NO_LEAP | Fixed year | On |
| `nv0` | NO_LEAP | Variable orbit year | Off |
| `nv1` | NO_LEAP | Variable orbit year | On |
| `gf0` | Gregorian | Fixed year | Off |
| `gf1` | Gregorian | Fixed year | On |
| `gv0` | Gregorian | Variable orbit year | Off |
| `gv1` | Gregorian | Variable orbit year | On |

Fixed profiles leave
the orbital year unset to exercise the compset default. Variable profiles
explicitly align orbital year 1999 with model year 1999; they are utility
experiments, not CMIP-protocol configurations. All profiles explicitly set EOT,
so these runs do not independently check its default-on policy.

| Boundary modification | Start date | Length | Restart interval |
| --- | --- | --- | --- |
| `calendar_feb` | `2000-02-27` | 6 days | 1 day |
| `calendar_rollover` | `1999-12-29` | 6 days | 1 day |
| `calendar_cycle` | `1999-01-01` | 1461 Gregorian or 1460 NO_LEAP days | 730 days |

Short tests run six days from February 27, 2000 or December 29, 1999 and restart
after one day. Cycle tests start January 1, 1999, restart after 730 days, and
run 1461 Gregorian days or 1460 NO_LEAP days, ending January 1, 2003. Day-based
lengths avoid CIME's fixed-365-day arithmetic for year-based restart pointers.
Including the restarted suffix, each long test simulates approximately six
model-years. The 730-day restart falls on December 31, 2000 in Gregorian and
January 1, 2001 in NO_LEAP; both runs end on January 1, 2003.

To run individual tests rather than a whole suite, use the composed profile and
boundary modifications. For example, from `cime/scripts`:

```bash
./create_test ERS_Ld6.ne4pg2_oQU480.F1850.allactive-calendar_gf1--allactive-calendar_feb
./create_test ERS_Ld1461.ne4pg2_ne4pg2.F2010-SCREAMv1.allactive-calendar_gv1--allactive-calendar_cycle
```

The first checks EAM Gregorian fixed-year restart behavior across leap day.
The second runs the four-year Gregorian variable-orbit EAMxx restart test.
The smoke suite selects `gf1` with `calendar_feb` and `gv1` with
`calendar_rollover` for each configuration. The short suite combines every
profile with both short boundaries; the cycle suite combines every profile
with `calendar_cycle`.

### Lightweight and numerical checks

Lightweight regression checks require Python and the checked-out CIME submodule,
but no model case or input datasets. Run from the repository root:

```bash
python3 share/test/unit/shr_cal_test/test_calendar_config.py
```

Append `-v` to display test method names and their results.

| Test method | Checks |
| --- | --- |
| `test_epoch_defaults_and_eot_default` | 1850, 2010, and 20TR epoch defaults, fixed-year default, and default-on EOT under both calendars and drivers |
| `test_profiles` | All eight profile namelists and their orbital-year/EOT settings for both drivers |
| `test_explicit_epoch_and_legacy_alias` | Explicit orbital-year/alignment overrides and the `variable_year` alias |
| `test_suite_modifications` | Suite counts, unique entries, and resolution of composed modifications |
| `test_restart_pointers` | Calendar-specific restart timestamps for both short boundaries and the four-year runs |
| `test_noleap_input_defaults` | Compatibility option defaults off under both calendars and drivers |
| `test_noleap_input_user_overrides` | Actual parsing of true/false user-namelist overrides for both drivers |
| `test_noleap_input_suites` | Compatibility-suite counts, uniqueness and modifier resolution |
| `test_noleap_system_test_loader` | Actual CIME discovery of the project-owned NLF/NLN classes |

These checks use both drivers' actual namelist generators and orbital-default
helpers. They check compset epochs, explicit overrides, the legacy variable-year
alias, default-on EOT, all eight profile namelists, suite modification resolution,
and CIME restart-pointer calculations. The existing shared calendar pFUnit tests
also cover aligned variable-year Gregorian dispatch, century leap rules, and
pre-1582 continuity. Run those tests in the shared-library unit-test environment.

The shared numerical test source is
`share/test/unit/shr_cal_test/test_shr_cal.pf`. In addition to its existing solar
geometry and EOT tests, it includes:

| Test | Checks |
| --- | --- |
| `variable_orbital_year_alignment` | Gregorian model years 1999-2003 aligned to orbital year 1850 at model year 1999, EOT on/off, and detection of frozen orbital parameters |
| `proleptic_gregorian_boundaries` | February boundaries in 1500, 1582, 1600, 1900, 2000, and 2100; continuity across October 1582; and the 1461-day cycle |
| `fixed_year_rollover_continuity` | Decreasing-step fixed-year continuity across common/leap/century rollovers, EOT on/off, and perpetual-date stability |

When `gfortran` is available, the existing Python configuration test also builds
the complete shared orbital module and existing orbital tester using minimal
kind/constants/log/abort stubs. It checks finite physical ranges, annual anchors,
decreasing-step rollover continuity, and unchanged 365-day arithmetic. Other
compiled probes execute both drivers' actual ocean-albedo geometry blocks at
dawn/dusk with both EOT signs, test the full EAMxx Fortran bridge against shared
dispatch, and read ordinary and maximum-length override strings through the
actual restart character routines with I/O stubs. These are not dependency-backed
model, Scorpio or PIO tests.

Positive forcing fixtures select `coszen` for solar-only streams and retain
linear interpolation for other streams; mixed-quantity streams do not establish
solar-only coverage. Instantaneous history already exports
`xaoo_Faox_swdn`, `xaoo_Faox_swup`, `xaoo_So_anidr` and `domo_mask`.
Solar fixture runs check these for missing/non-finite active values, albedo
range, the negative-upward flux convention and excessive magnitudes. The default
Earth-fixture downward-SW qualification limit is 2000 W m-2, configurable through
the fixture manifest's `solar_flux_limit`; it does not clip production fields.
Valid inactive masks and genuine night zeros remain allowed. Albedo and flux
history use different scheduled evaluation times, so albedo alone is not used
as a daylight oracle for a recorded flux. Coupled daylight/zero-flux, forcing
energy and cyclic/restart comparisons still need supported-machine execution.

### Solar conservation suite

`e3sm_solar_conservation` is an opt-in coupled shortwave conservation suite. It
does not inherit into the ordinary developer suites. Its `SCV` system test
extends ERS and checks both the continuous run and restarted suffix, rather
than equating successful execution or bit-for-bit restart with conservation.

**Qualification limits:** SCV implements TOA/heating, actual imported-albedo
consumption, global surface exchange and unequal-ocean-cadence inventory
checks for EAM/MCT. It does not prove full prognostic thermodynamic storage
closure or local land/ice-column closure in a running model. Unequal land/ice
cadence, EAMxx, MOAB and variable-orbit cases remain outside this suite.
Global surface residuals can hide compensating local errors even though
atmospheric columns and ocean transitions also have pointwise assertions.
The compiled probes use infrastructure stubs, not a dependency-backed EAM/MCT
build. Supported-machine builds, MPI/log collection and coupled continuous/
restart runs are still required for runtime qualification. No such coupled
run has been performed on the development host.

For scientific configuration and restart migration, see the user guide's
[calendar selection](../user-guide/index.md#calendar-selection) discussion.
The [failure triage](#local-validation-and-failure-triage) below maps diagnostic
failures to the first quantities to inspect.

#### Coverage and execution

| Suite suffix | Cases | Coverage |
| --- | ---: | --- |
| `smoke` | 4 | Two days near the March equinox |
| `boundaries` | 28 | Six days across February 2000, February 2100, year rollover, both equinoxes and both solstices |
| `cadence` | 4 | Six days with shortwave radiation every second atmosphere timestep |
| `unequal` | 12 | Six days with ocean delivery every two or four atmosphere steps, including four-step ocean delivery with two-step radiation |

Each group crosses `NO_LEAP`/`GREGORIAN` with EOT off/on and uses fixed-year
orbits. All cases use prognostic atmosphere, land, ocean and sea ice through
`WCYCL1850NS` at `ne4pg2_oQU480`. Prescribed-ocean or prescribed-ice F cases
are intentionally not used as whole-system conservation examples. Each case
restarts after one day. The base modification sets atmosphere, land, ice,
ocean and runoff coupling to 48 exchanges/day, atmosphere `dtime=1800`,
`use_rad_dt_cosz=.true.`, and `iradsw=1`; the cadence modification sets
`iradsw=2`. Solar averaging widths are therefore 1800 or 3600 seconds. The
unequal-cadence modifications reduce ocean coupling to 24 or 12 exchanges/day;
land and ice remain at 48/day. There are 48 cases in the complete family.

Run only on a supported machine, for example:

```bash
cd cime/scripts
./create_test e3sm_solar_conservation_smoke --machine MACHINE --compiler COMPILER
./create_test e3sm_solar_conservation
```

#### Timing and diagnostic records

The base modification enables `BUDGETS` and `budget_inst=4`. This new MCT
diagnostic level adds full-precision `SOLAR_EXCHANGE_V1`,
`SOLAR_SCHEDULE_V1`, `SOLAR_COLUMN_V1`, `SOLAR_CONSUME_V1` and
`SOLAR_BUFFER_V1` records. Levels
0 through 3 are unchanged. The exchange record is five signed global means
in W/m2: atmospheric net surface SW (negative), land absorption, ocean
receipt, northern ice retention and southern ice retention. Existing budget
reductions supply instantaneous area-weighted products, not products of
averaged areas and fluxes. Ocean receipt already includes ice transmission;
ice retention subtracts it. Never add transmission again to that global sum.

The checker integrates each record over its actual 1800-second interval and
reports component energies and residuals in J/m2 of Earth's surface. Driver
budget timestamps are interval ends. Missing, duplicated, reordered,
non-finite, wrongly signed, or wholly zero-throughput records fail. Every
instantaneous residual must satisfy `abs(residual) <= 1e-7 + 1e-10*scale`
W/m2, where `scale` is the sum of absolute component terms. The integrated
test uses the corresponding duration-scaled absolute bound. These are
explicit initial qualification limits, not calibrated solver/remapping
error estimates: supported-machine runs must establish whether they are
appropriate. Investigate failures rather than silently relaxing limits.

The independent civil-date checker validates calendar days, averaging widths,
restart provenance, and upcoming scheduled events. For this equal-timestep
EAM configuration the schedule is two driver timesteps ahead: the atmosphere
interface checks its internal clock minus one timestep against the driver
clock, then asks radiation for the following internal timestep. In this
equal-step EAM interface, the timeline is:

```text
driver exchange label:       t_driver
actual EAM internal clock:   t_atm   = t_driver + dt_atm
next valid SW event:         t_valid = t_atm + dt_atm  (when valid)
forward geometry interval:   [t_valid, t_valid + dt_avg]
```

The scheduling lead is distinct from the averaging width `dt_avg`; neither
is a universal offset for other components or configurations. Every-step radiation requires a
valid event at every exchange; two-step radiation alternates valid events
after startup. A future albedo timestamp is not itself an error. The restart
suffix must reproduce the continuous run's overlapping exchange, schedule,
column, consumption and buffer records exactly, excluding the deliberate
fresh-run/restored provenance-marker transition. Only fresh content from all
`*.log.*` files is inspected, including compressed logs. Atmosphere ranks can
write to the main executable log rather than the root atmosphere log; all
rank output must be retained for complete spatial coverage. Results are written
as `solar-conservation.base.json` and `solar-conservation.rest.json` in RUNDIR.
Log snapshots retain decoded lengths and SHA-256 hashes, not full raw text.
The checker verifies existing prefixes, streams new lines through one record
dispatch pass, and rejects conflicting compressed/uncompressed copies rather
than silently choosing one. Parsed records are retained for restart comparison.

#### Atmosphere and surface budgets

EAM captures the four actual imported albedo channels after `atm_import`.
Their order is visible direct, visible diffuse, near-infrared direct and
near-infrared diffuse; values are dimensionless and columns are physics columns.
On each shortwave radiation call, it checks that the consumed channels are
unchanged and that their scheduled civil date and time equal EAM's actual
internal clock. Held-radiation steps do not consume new albedos. These checks
are independent of the published-schedule oracle; a future history label is
not treated as stale merely because it differs from driver time.

Both RRTMG and RRTMGP check each column before cached heating is converted
for reuse. The modeled-column absorption is `sum(qrs*pdel)/gravity`, with
`qrs` in W/kg and pressure thickness in Pa. `FSNT` is at the model top, not
the top of the extra radiation layer. The diagnostics therefore separately
integrate that layer's heating and check it against the TOA-to-model-top
flux difference. Full column closure uses true TOA net SW minus surface net
SW minus modeled-column and extra-layer absorption. Extra-layer absorption
is radiative absorption, not a claim of heating applied to prognostic model
levels. Net TOA SW is not expected to vanish.

`SOLAR_COLUMN_V1` records date, time, chunk, column count, maximum pointwise
residual, area-weighted TOA/surface/total-heating terms, maximum column scale
and summed area fraction. Physics-grid areas are in radians squared and are
normalized by `4*pi`. The checker requires complete time and chunk coverage,
unit global area and agreement with the actual atmospheric surface exchange.
`SOLAR_CONSUME_V1` records each checked date, time, chunk and column count.
The opt-in `SOLAR_TOP` global physics-buffer field saves TOA flux and
extra-layer absorption in native EAM restarts, including held-radiation steps.
Registration occurs before component initialization; runtime consumption
tracing begins only in `cime_run`, after initialization radiation.

#### Ocean inventories and full-chain closure

Ocean instrumentation surrounds the actual MCT accumulation, averaging and
copy to the component input. `SOLAR_BUFFER_V1` records date, time, instance,
`ADD`/`SEND`, actual sample count and before/incoming/delivered/after energy
inventories. Inventories use the actual accumulated field and restored
counter, driver interval, ocean-domain area and ocean-plus-ice fraction.
Each cell and global transition must satisfy
`before + incoming - delivered - after = 0`. Counts, continuity, interval
coverage and incoming energy are checked independently against the coupled
exchange records. Initial and final pending inventories are retained, not
assumed zero. A nonzero delivery with zero accumulated samples fails as an
unfunded initialization transfer; it is not silently discarded.

The integrated full-chain test requires TOA energy minus atmospheric
absorption, land absorption, retained ice absorption, ocean-delivered energy
and the change in pending ocean energy to close. Both calendars use an
independent civil-date oracle. Slower ocean coupling is supported only when
its frequency divides 48/day; a component cadence faster than the atmosphere
is rejected because it changes the driver sampling interval.

#### Local validation and failure triage

Local validation includes actual CIME class/modifier resolution, compiled
real-source atmosphere and ocean-accumulator probes, independent Python
calendar fixtures, negative controls and a mocked ERS/log lifecycle at equal,
two-step and four-step ocean cadences. Controls reject lost heating, partial
spatial coverage, shifted or changed albedo consumption, missing/duplicated
buffer samples and an internally closing but incorrectly weighted ledger.
The mocked runner carries nonempty initial/final inventories across restart
and rejects a reset ledger even when its individual transitions close.
An additional compiled probe executes Icepack's actual `icepack_prep_radiation`
routine with changed spectra, fractional ice, surface/interior/transmitted
partitions, layer absorption, darkness and a zero-reference fully reflecting
state. It tests broadband conservation, not unsupported per-band interior
absorption invariants. These probes are not coupled model qualification.

| Failure | First checks |
| --- | --- |
| Stale or changed consumed albedo | Compare the reported chunk, actual EAM clock and scheduled valid time; inspect import order and the reported channel before considering history labels. |
| TOA/heating or upper-layer closure | Use the reported column, residual, scale and tolerance; check W/kg and Pa conversions, `FSNT` versus true TOA, and independent extra-layer absorption. |
| Missing columns or consumption events | Retain every rank's executable/atmosphere log; check global area and actual radiation cadence before changing coverage requirements. |
| Ocean buffer or full-chain closure | Inspect ADD/SEND counts, area-times-fraction weights, incoming energy and pending inventories; do not zero pending state or add ice transmission twice. |
| Restart overlap mismatch | Compare native `SOLAR_TOP`, ocean accumulator/counter and scheduled-event state with the continuous run; a locally closing reset ledger is still wrong. |
| Replaced or conflicting logs | Preserve append-only run output and matching gzip/plain copies; do not mix old and new segments or discard rank logs. |

Preserve failing records and reports before changing code or tolerances. Local
probes isolate arithmetic and lifecycle defects; only supported-machine runs
can establish behavior with the real coupler, dependencies and MPI layout.

#### Forward albedo references

The archived CCSM3 sea-ice scientific description in
`components/cice/docs/SciGuide/ice_scidoc.tex` explicitly explains computing
albedos last so the next atmosphere radiation and ice thermodynamics use
the same albedos. That historical contract supports preserving forward
scheduling; it is not proof of every current runtime cadence.

The [CESM sea-ice history guide](https://cesmcice.readthedocs.io/en/latest/users_guide/ice_history.html)
warns against constructing budgets from unweighted per-ice-area history
averages. Consult the
[E3SM MPAS sea-ice technical guide](https://docs.e3sm.org/E3SM/MPAS-seaice/tech-guide/)
and [E3SM Icepack documentation](https://e3sm-icepack.readthedocs.io/en/latest/)
for the host/column boundary and radiation scheme. Prescribed ice does not
provide a physical storage-conservation oracle.
[Dang, Zender and Flanner (2019)](https://doi.org/10.5194/tc-13-2325-2019)
discuss reflectance, transmission, layer absorption and compensating spectral
errors; their paper does not define the exact coupled albedo event schedule.
Briegleb and Light (2007), NCAR/TN-472+STR,
[doi:10.5065/D6B27S71](https://doi.org/10.5065/D6B27S71), and Holland et al.
(2012), [doi:10.1175/JCLI-D-11-00078.1](https://doi.org/10.1175/JCLI-D-11-00078.1),
are relevant references, but their full texts were not retrieved during this
review; no additional timing claim is inferred from their citations alone.

### Ocean shortwave intervals

EAMxx RRTMGP and MAM use `physics::solar_geometry` from shared physics rather
than selecting calendar/orbital/EOT branches inside each process. The returned
`SolarGeometry` provides the eccentricity factor and host cosine-zenith
evaluation for a caller-supplied averaging interval. RRTMGP retains its radiation
interval and fixed-angle/irradiance overrides; MAM retains its timestep interval.
The shared library owns the Fortran orbital bridge, and coupled driver geometry
is authoritative for both processes even when EOT is disabled. This also corrects
MAM's former NO_LEAP/EOT-off process-local orbital fallback.

The atmospheric `nextsw_dt_avg` metadata specifies seconds of forward cosine
averaging from the scheduled shortwave timestamp. Zero means instantaneous;
the internal unknown value is not a fallback. EAM exports its initialized
`dt_avg`, EAMxx exports radiation frequency times the actual atmosphere timestep,
and DATM explicitly exports zero. Both driver metadata APIs broadcast the width
with the schedule. A timestamp update without a width invalidates the old width.
Other atmosphere producers must declare a supported interval before using
computed ocean albedos.

Both ocean routines retain the current saved-albedo spectral reflection and
absorption pairing before preparing the next albedos. Merge/accumulation must
remain before that update. The existing albedo formula uses the atmospheric
mean cosine, allowing absorption in sunrise-crossing intervals without discarding
received gross fluxes. This intentionally changes affected Gregorian and NO_LEAP
dawn/dusk baselines, independently of EOT. It does not align every component's
separate daylight mask or averaging window.

All three restart writers persist `seq_infodata_nextsw_dt_avg` and
`seq_infodata_ocean_solar_version=1`. Computed-ocean continuation validates saved
provenance during initialization and checks the newly published width against
the saved width before reuse. No-ocean and constant-albedo paths are exempt from
this ocean-specific compatibility requirement.

Existing Python regressions compile both real ocean albedo/flux blocks and
metadata validators, check each spectral band and full-reflection closure, and
retain the merge-before-update ordering guard. They also compile actual EAM
interval getters and link the actual EAMxx scheduler to its complete Fortran
interface, exercising radiation cadence, disabled radiation, missing RRTMGP,
and coupled geometry with EOT off. The same test links the actual shared C++
solar API to the real Fortran orbital routines and checks both calendars,
fixed/current-year and explicit-element orbits, EOT settings, averaging widths,
and callback failures/conflicts using clock/parameter infrastructure stubs.
Infrastructure-stub probes do not qualify
real MPI/PIO/Scorpio exchange, coupled energy budgets, or uninterrupted versus
exact-restart model runs; those remain supported-machine requirements.

### Instantaneous terminator geometry

Coupled Icepack uses the scheduled calendar day for both declination and hour
angle in NO_LEAP/EOT-off runs; the local sea-ice clock no longer supplies a
different hour angle. Gregorian/EOT timestamp validation and standalone behavior
are unchanged. ELM refreshes its driver orbital context before forcing import.
Its coupler-bypass solar interpolation uses `elm_solar_cosine`, not the legacy
365-day sinusoidal declination approximation. Target and normalization midpoints
use calendar-aware signed offsets, retaining seconds and fractional midpoint
seconds across midnight, leap-day and year boundaries. This does not change the
bypass readers' separate forcing-calendar/indexing limitations.

Shared cosine interpolation and bypass solar scaling use the geometric horizon
at cosine zero, not the former `0.001` daylight cutoff. Normalization averages
include nonnegative cosines; all-dark means produce zero scaling without
division by zero. Shared averaging rejects nonpositive intervals and timesteps
on both calendars and weights a partial final sample only over the remaining
interval. Bypass retains its midpoint quadrature, existing ratio cap and
prepartitioned prescribed-flux handling.

Common instantaneous geometry does not require common flux masks. Preserve EAM's
optional radiation averaging and instantaneous chemistry, EAMxx's radiation and
MAM timestep intervals, and land/ice instantaneous scheduled albedos. Ocean
albedos retain the atmosphere-matched interval described above. Terrain and
fixed-zenith overrides and prescribed fluxes are not instantaneous horizontal
daylight masks.

Existing numerical subtests compile the actual Icepack routine, ELM helper,
bypass interpolation block, forcing cosine/scaling blocks and shared averaging
routine against the real orbital implementation. Calendar/context/time-service
fixtures isolate these procedures; golden date/offset cases check sample
forwarding, not full ESMF integration. These repairs intentionally change
affected sea-ice and forcing/bypass baselines, including NO_LEAP with EOT off.
Cross-version saved forcing/albedo states are not claimed to continue bit for
bit. Full component builds, actual-reader normalization/energy budgets and
uninterrupted/restart qualification remain supported-machine requirements.

### Qualification limits

These are integration test definitions, not certified input-data combinations.
Before qualification, verify initial-condition compatibility and forcing-year
alignment, cycling, and leap-day interpolation for each configuration. Do not
relabel forcing calendars to make a test pass. Archive profile-specific
baselines; physically different profiles are not expected to agree bit for bit.

ERS checks continuous/restarted reproducibility within each profile. It does
not prove the orbital default is correct, variable orbital parameters advance,
or EOT has the correct sign and magnitude. The lightweight checks validate the
configuration defaults, and the shared numerical tests validate Gregorian
variable-year dispatch against explicit orbital calculations. Independent EOT
reference values, NO_LEAP driver orbital updates, and checks of component output
calendars and forcing interpolation remain scientific qualification requirements.
Passing configuration checks alone does not certify the numerical tests or full
model runs. Production builds and model execution must occur on supported machines.

## No-leap input compatibility tests

The [source-calendar identification audit](../user-guide/index.md#identifying-an-input-calendar)
distinguishes file metadata from model-calendar settings and documents reader
coverage by input path. Neither a CMIP label nor an aerosol/component name
establishes coverage. EAM `VOLC_CMIP6`, EAMxx SPA `DataInterpolation`, and ELM
`CPL_BYPASS` meteorology/aerosol deposition require separate qualification;
passing shared-stream or `TimeInterpolation` tests does not qualify those paths.
Mixed CF `standard`/`gregorian` epochs, records and bounds before 1582-10-15 are
rejected in the audited readers and preflight, even for modern record dates with
an ancient epoch. Explicit proleptic Gregorian remains distinct and supported
within each reader's existing limits. Inventory actual forcing files and their
provenance; the preflight inventory remains explicitly partial.

These opt-in tests exercise `fudge_noleap_for_gregorian`, described in
[Using No-Leap Forcing With Leap Years](../user-guide/index.md#using-no-leap-forcing-with-leap-years).
The original calendar suites and their counts remain unchanged. No compatibility
suite is inherited by a standard developer suite.

!!! warning "Coverage is not yet complete"
	The Python configuration, harness and generated-NetCDF checks run locally.
	Coupled NLF/NLN, EAMxx and shared pFUnit tests still require dependency-backed
	execution. The native-reader and guard gaps below remain qualification
	blockers. Passing a smoke test does not certify every E3SM input reader,
	existing input dataset, conservative flux handling or cross-calendar state
	conversion.

### Source-calendar Resolution Checks

`share/util/shr_input_calendar_mod.F90` owns canonical aliases, explicit
metadata/override conflict handling, exact rule lookup, immutable registration,
finite/ordered coordinate checks, source date checks, restart provenance and
the C binding used by coupled EAMxx. Shared streams and chemistry use the shared
Fortran wrapper. Standalone EAMxx retains equivalent metadata-only resolution;
it does not acquire a dependency on the coupler's Fortran library.

The dependency-free contract test can run on a local Fortran installation:

```bash
build_dir=$(mktemp -d)
gfortran -ffree-line-length-none -fcheck=all -Wall -Werror \
	-J "$build_dir" share/util/shr_input_calendar_mod.F90 \
	share/test/unit/shr_cal_test/test_input_calendar.F90 \
	-o "$build_dir/test_input_calendar"
"$build_dir/test_input_calendar"
python -m unittest discover -s share/test/unit/shr_cal_test -p 'test_*.py'
```

The standalone CMake target is `shr_input_calendar_exe`, registered as the
`shr_input_calendar` test alongside the existing shared calendar test. It uses
a separate module output directory to avoid competing `.mod` writes.

The Fortran test executes the real resolver, validation, registry and C bridge,
including missing metadata, conflicts, exact matching and immutability. Python
tests use both actual driver namelist generators and on-disk NetCDF fixtures
for preflight. Cross-language contract checks can load the compiled C binding
to compare canonical normalization with preflight. These do not replace the
actual Scorpio/MCT/PIO reader tests or a restart-file round trip. In particular,
both drivers' character restart I/O now accommodates the full canonical rule
record; dependency-backed qualification must verify that long record survives
writing and reading without truncation.

The contract cases distinguish mixed and proleptic metadata/override flavors,
reject the year-0001 mixed epoch even for modern records, retain historical
proleptic input, and reject negative offsets and bounds crossing the reform
cutoff. The C callback returns source flavor as well as category and provenance;
its Fortran and C++ declarations must change together. Gregorian fixed-year
restarts also carry `seq_infodata_solar_phase_version=1`; both MOAB writer paths
emit it. Exact continuation rejects an absent/old marker for that convention
without imposing it on unaffected NO_LEAP runs.

The user guide describes explicit override identifiers and the partial-inventory
preflight command. Keep its aliases/conflict rules synchronized with the shared
resolver. Do not infer calendars from filenames, record counts or missing dates.

### Suite Matrix

| Suite | Tests | Definition |
| --- | ---: | --- |
| `e3sm_calendar_noleap_smoke` | 2 | A configuration, daily/subdaily point forcing across February 29 |
| `e3sm_calendar_noleap_short` | 56 | Four configurations x two short boundaries x seven fixture/control choices |
| `e3sm_calendar_noleap_negative` | 11 | Specific metadata, interpolation, model-calendar and restart-policy rejections in A |
| `e3sm_calendar_noleap_cycle` | 2 | A and I1850ELM, one common forcing year cycled over a 1461-day Gregorian run |

The short configurations are `A` on `ne4pg2_oQU480_rx1`, `CMPASO-NYF` on
`T62_oQU120`, `DTESTM` on `T62_oQU240`, and `I1850ELM` on `f19_g16`.
The seven choices are daily, subdaily, monthly, Gregorian-source control,
explicitly disabled policy, omitted-option default, and common-year cycling.
Both `calendar_feb` and `calendar_rollover` are included. Gregorian `gf1`
supplies fixed-year orbital settings and EOT on; these suites are not another
orbital-profile cross product. The cycle tests use `calendar_cycle`.

`NLF` inherits ERS: run continuously, restart from the configured checkpoint,
and perform the normal history comparison. Before **each** run leg, generate
namelists, replace the generated prescribed streams with test-owned fixtures,
then skip subsequent namelist regeneration for that run. Preserve spatial
schemas and domain files. Temperature varies in time but is spatially constant;
other forcing fields retain their first source slice. Point metadata describe
these synthetic records, not relabeled production interval means.

Every run leg must produce fresh instantaneous `*.cpl*.hi.*.nc` files containing
`a2x_Sa_tbot`. The MCT writer constructs the name from prefix `a2x` and attribute
`Sa_tbot`; time is Gregorian elapsed days since the run's reference start.
The driver advances clocks, runs/receives DATM, then writes history. DATM queries
the shared stream using its current component date and time. History is scheduled
at `24 / gcd(24, ATM_NCPL)` hours, so snapshots coincide with atmosphere coupling
boundaries without changing the component's coupling frequency. Only positive
day-based coupling frequencies are supported by this harness.

The export comparison uses a `1e-4 K` absolute tolerance, reads actual CF time
units/calendar rather than inferring time from filenames, and rejects stale
first-leg files. Supported DATM transformations are `COPYALL`, `CLMNCEP` and
`CORE2_NYF`. For CORE2, compare only latitudes within 60 degrees of the equator:
polar temperature corrections occur **after** interpolation. Anomaly temperature
forcing and unverified DATM modes fail closed. There is no independent numeric
assertion on every DOCN, DICE or prescribed MPAS-Seaice field.

Enabled runs require the global policy banner and an activation-description
warning for **each** staged no-leap stream. Disabled runs reject an enabled
banner. This checks presence and stream identity, not MPI-root exclusivity or
the absence of repeated warning messages. No annual-integral conservation or
log-byte equivalence is asserted.

`NLN` succeeds only if a run fails with its specified diagnostic in new/changed
logs. A successful run, missing input, build failure or unrelated abort cannot
stand in for the expected rejection. Restart-policy cases first complete a
successful leg, change `user_nl_cpl`, regenerate/stage the second leg, and require
the restart-policy mismatch. Staging occurs inside `run_indv`, after ERS resets
its namelist flag. The negative cases are:

| Fixture | Required rejection |
| --- | --- |
| `missing_calendar`, `empty_calendar` | Missing/empty explicit source calendar |
| `bad_method`, `mean` | `time: sum` / `time: mean` without a supported interval adapter |
| `ambiguous_bounds` | Real time bounds with no field cell-method declaration |
| `bad_algorithm`, `bad_lower`, `bad_upper` | Nearest/lower/upper instead of the allowed adapted algorithms |
| `bad_model_calendar` | Enabled policy with NO_LEAP model calendar |
| `restart_on`, `restart_off` | Policy changes from off to on or on to off |

### Local checks

From the repository root:

```bash
python3 share/test/unit/shr_cal_test/test_calendar_config.py -v
python3 -m unittest discover -s share/test/unit/shr_cal_test -p 'test_*.py' -v
```

The first command requires only Python and CIME. File-backed checks in the second
also require NumPy and `netCDF4`; otherwise they are reported as skipped, not
verified. They generate temporary fixtures, check missing/bounded/mean metadata,
test stale logs and exports, corrupt an export to prove the numeric validator
fails, and check restart staging and coupling-aligned history setup. These are
harness tests, not executions of the Fortran or C++ readers.

Independent numeric references include 10 and 14 at February 28 00:00 and
March 1 00:00, giving 12 on February 29 00:00 and 13 at noon. Moving the first
record to February 28 18:00 gives 10.8 at February 29 00:00, because the bracket
is 30 hours. Monthly references use actual elapsed Gregorian time. Disabled
shared-stream references repeat February 28 and use a 365-day source clock.
Do not use the production mapping routine to calculate expected results.

### Dependency-backed runs

In a configured EAMxx standalone build, build target `time_interpolation`, then
run `ctest --test-dir <eamxx-build> -L calendar_noleap --output-on-failure`.
The label selects real Scorpio-backed `TimeInterpolation` tests at the configured
MPI rank counts. The ordinary interpolation registration excludes this tag to
avoid simultaneous tests writing the same fixtures. Cases cover daily,
subdaily, monthly and multiyear offsets; common/leap century years; aliases;
Gregorian and policy-off controls; missing/unsupported/mixed calendars;
nonfinite/negative offsets; duplicates; unsupported cell methods; ambiguous
bounds; and lower/upper coverage errors. These tests are added but not yet
compiled or executed in the current environment.

The existing shared pFUnit target `shr_cal_exe` also includes input cell-method
parsing and calendar-alias checks. It does not initialize/mutate the immutable
policy or exercise stream I/O. Run it in the configured shared-library unit-test
environment. The historical stream tester's obsolete Makefile is not a
replacement for a current dependency-backed reader target.

On a supported machine, run from `cime/scripts`:

```bash
./create_test e3sm_calendar_noleap_smoke
./create_test e3sm_calendar_noleap_short
./create_test e3sm_calendar_noleap_negative
./create_test e3sm_calendar_noleap_cycle
```

Run smoke before larger suites; use machine-appropriate layouts and walltimes.
The cycle suite is expensive. Python fixture tooling must be available on the
execution host. Existing input datasets are still required for domains, spatial
schemas and component initial conditions. Do not create model cases on an
unsupported workstation. No coupled runs or scientific baselines have been
generated for these suites here.

### Outstanding reader and guard gates

| Reader or contract | Current coverage | Remaining gate |
| --- | --- | --- |
| Shared prescribed point streams | NLF/NLN definitions and locally checked fixtures/oracles | Execute each configuration; independently test all advertised output fields |
| Shared cycling/cache/file transitions | Cycle-suite definitions; EAMxx file sorting tests | Real shared-reader multi-file/cache and year-alignment values |
| Shared coszen path | Not exercised: synthetic streams select linear interpolation | Gregorian solar normalization and mapped interval duration |
| EAMxx file interpolation | Real-reader Catch2 cases added | Compile/run with Scorpio, one/multiple MPI ranks |
| EAM chemistry `time_coordinate` | No component-linked reader test yet | Numeric offsets, date/datesec and try_dates, both mean-bound endpoints, finite-time checks |
| Time-varying EAM FLBC / tracer_data | Production rejection guards, no automated guard harness | Invoke real initializers; fixed-input controls; require their specific diagnostics |
| EAMxx IOP | Production rejection guard, no automated guard harness | Activate IOP through the real driver and require its adapter diagnostic |
| Native MPAS forcing | Production rejection guard, no automated guard harness | Each coupled owner, active forcing-group rejection and inactive controls |
| Coupled WW3 | Production rejection guard, no automated guard harness | Reach real wave initialization and require its native-input diagnostic |
| Restart provenance | NLF same-policy and NLN bidirectional policy-change definitions | Execute; missing/corrupt markers, immutable-policy mutation, ELM calendar-restore mismatch |
| Default-off / Gregorian-source controls | Separate NLF definitions and local default/override checks | Paired omitted-versus-explicit-off BFB and unchanged Gregorian-input comparison |
| Restart placement and warnings | Existing boundary checkpoints and warning-presence checks | Restart on/after February 29; root-only and bounded warning-frequency assertions |

Qualification remains blocked until these gaps are closed and the relevant
tests execute successfully. On 2026-10-27 the local host was not recognized by
CIME; CMake, MPI tools and NetCDF-Fortran were unavailable. Source inspection,
Python test success and clean editor diagnostics do not remove those gates.

## Gregorian Input for No-Leap Tests

The [source-calendar identification audit](../user-guide/index.md#identifying-an-input-calendar)
and its specialized-reader qualification gaps also apply to this reverse policy.

This experimental reverse policy uses `fudge_gregorian_to_noleap`, default
false, with `CALENDAR=NO_LEAP`. It is mutually exclusive with the forward
flag. See [the user contract](../user-guide/index.md#using-gregorian-forcing-with-no-leap-models).
The policy retains Gregorian interpolation, including February 29 source
brackets. No leap-day record compression or annual conservation is intended.

### Source-Time Architecture

`shr_cal_mod` owns immutable forward/reverse policy flags and separate
forward, reverse and any-policy queries. Both drivers initialize/broadcast
the reverse flag and write/read its own restart marker alongside the model
calendar and canonical override signature. Missing reverse provenance
permits only reverse-disabled continuation. Malformed markers and policy,
calendar or rule changes are rejected before normal component reads.

Shared streams decode Gregorian coordinates normally. Enabled reverse
queries, cached bounds and interpolation factors use the aligned source
year. `shr_stream_inputDate` supplies the same mapping to cache checks and
weights. Cycle-edge brackets shift by the Gregorian duration of the entire
source cycle, avoiding invalid February 29 model-year labels. Raw source
dates and file indices are unchanged. Coszen normalization uses NO_LEAP
model-time samples and remains a dependency-backed qualification gate.

Chemistry's `time_coordinate` records an explicit Gregorian input view for
numeric/date decoding, bounds, fixed/current queries, searches, copies and
teardown. Source coordinates and model civil queries use the shared explicit
Gregorian elapsed-time helper. Input `dtime` is applied in the input-time
view, not by changing the model clock.

EAMxx resolves source metadata before parsing `case_t0`. Its reverse file
path uses scalar Gregorian seconds for ordering, duplicate detection,
coverage, cached bracket updates and weights. Gregorian source leap days
never need a NO_LEAP model `TimeStamp`. Coupled EAMxx receives an explicit
Gregorian-seconds callback backed by `shr_cal_mod`; standalone input uses
an explicit Gregorian fallback. Both test files contain matching known-date
fixtures, including century years. Model timestamps and field-to-field
interpolation retain their existing semantics. No global calendar toggles
are used for file decoding.

Existing unsupported-reader guards now recognize either policy. ELM checks
NO_LEAP component restart state before allowing reverse continuation.
These source changes are not proof that native guards or restart I/O have
executed successfully.

### Independent References

| Gregorian Source Points | NO_LEAP Query | Expected Value |
| --- | --- | --- |
| Feb 28 00:00=10, Feb 29 00:00=100, Mar 1 00:00=14 in 2000 | Feb 28 12:00 | 55 |
| Same three records | Mar 1 00:00 | 14 |
| Feb 28 00:00=10, Mar 1 00:00=14 in 2000 | Feb 28 12:00 | 11, not compressed-calendar 12 |
| Feb 28 18:00=10, Feb 29 00:00=14 in 2000 | Feb 28 21:00 | 12 |
| Feb 28 18:00=10, Mar 1 00:00=14 in 2000 | Feb 28 21:00 | 10.4, not 12 |

Century controls include 1900, 1999, 2000, 2004, 2100 and 2400. Gregorian
offsets 59/60 from January 1, 2000 remain distinct February 29/March 1
records. A leap-day source epoch is valid; a February 29 model timestamp
is not. Four NO_LEAP model years contain exactly 1460 days.

### Opt-In Suites

| Suite | Tests | Matrix |
| --- | --- | --- |
| `e3sm_calendar_gregorian_to_noleap_smoke` | 2 | A, daily/subdaily across February. |
| `e3sm_calendar_gregorian_to_noleap_short` | 56 | Four configurations x two boundaries x seven fixture/control choices. |
| `e3sm_calendar_gregorian_to_noleap_negative` | 18 | Specific metadata, quantity, algorithm, policy and restart rejections in A. |
| `e3sm_calendar_gregorian_to_noleap_cycle` | 2 | A and I1850ELM, leap-year source 2000 cycled over 1460 NO_LEAP model days. |

The configurations are `ne4pg2_oQU480_rx1.A`, `T62_oQU120.CMPASO-NYF`,
`T62_oQU240.DTESTM` and `f19_g16.I1850ELM`. Short boundaries are February
and year rollover; variants are daily, subdaily, monthly, matching-NO_LEAP,
explicit off, unset default and cycling. All use the existing `nf1` profile.
None inherits into normal developer suites. Existing forward and original
calendar suite memberships/counts remain unchanged.

The 18 negatives are missing/empty/unsupported/mixed metadata, sum, mean,
ambiguous bounds, nearest/lower/upper algorithms, wrong model calendar,
both flags, contradictory override, restart enable/disable, changed restart
override rules, missing reverse marker and corrupt reverse marker.
NLN must observe the specified diagnostic in fresh logs; an unrelated
failure does not pass. Marker mutations are limited to test-owned run files
identified by the generated coupler restart pointer.

NLF/NLN reuse the existing harness and actual generated stream schemas.
Reverse fixtures include a distinct February 29 temperature sentinel.
The reverse oracle evaluates source Gregorian intervals, projecting cycling
queries into the configured source year. Export validation decodes the
actual CF calendar, requires NO_LEAP output, rejects February 29 model
samples and retains the established DATM transformation/CORE2 latitude
restrictions. Synthetic date enumeration uses NO_LEAP decoding, not Python
Gregorian timedeltas across leap day. NLF retains ERS comparisons and fresh
atmosphere-cadence forcing-export checks; those model runs are not executed
by the local Python tests.

On a supported machine:

```bash
cd cime/scripts
./create_test e3sm_calendar_gregorian_to_noleap_smoke
./create_test e3sm_calendar_gregorian_to_noleap_short
./create_test e3sm_calendar_gregorian_to_noleap_negative
./create_test e3sm_calendar_gregorian_to_noleap_cycle
```

For an existing configured EAMxx standalone build, use
`ctest --test-dir <eamxx-build> -L calendar_gregorian_to_noleap --output-on-failure`.
The opt-in actual-Scorpio Catch2 cases cover retained/sparse leap-day
brackets, daily/subdaily/monthly inputs, source epochs, reversed file order,
matching-NO_LEAP input and coverage failures at one/multiple MPI ranks.
The registered C++ target and added shared pFUnit callback fixtures still
need compilation/execution with their real dependencies.

### Verification Status

Local checks exercise actual CIME namelist generators, modifier resolution,
preflight and NetCDF fixture/export/restart-marker logic. Full-module
preprocessing and Fortran syntax parsing are separate checks, not semantic
builds. EAMxx editor diagnostics are not C++ compilation. The standalone
metadata resolver is independently compiled/run, but is not a reader test.

Still required: supported-machine MCT/MOAB and component builds; all four
reverse suites; actual shared/chemistry/Scorpio reader values and cache/file
transitions; coszen with NO_LEAP model solar context; standalone/coupled
Gregorian-seconds parity; long restart override-record I/O; split-run BFB;
paired off/unset BFB; real native-reader/ELM guard executions; root-only,
bounded warning behavior; and a full documentation build. Local success does
not establish production readiness, universal input coverage or scientific
budget conservation. Supported-machine work remains explicitly deferred.
