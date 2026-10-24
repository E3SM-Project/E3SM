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

NO_LEAP orbital phase retains the day-80.5 vernal-equinox reference.
Independently, `orb_mode = 'fixed_year'` anchors the reference to calendar day
80.5 each year under either calendar. For Gregorian fixed-year runs, orbital
speed uses the actual 365- or 366-day year. Day 80.5 is March 21 at noon in a
common year and March 20 at noon in a Gregorian leap year. This annual convention
does not preserve continuous orbital phase across changes in civil-year length.

Gregorian variable-year and explicit-parameter modes retain continuous elapsed
model time relative to 2000-01-01 00:00, with an idealized 365.2422-day period
and the legacy January-1 phase calibration. Their phase does not reset at
civil-year boundaries. Fixed-year and explicit orbital elements remain fixed;
variable-year coefficients use the configured year alignment and fractional
civil year. Perpetual-season geometry uses the configured reference date.
These are Berger-compatible model conventions, not precision ephemerides;
the truncated orbital series need not give exactly zero declination at the
reference for nonzero eccentricity.

The annual fixed-year convention changes Gregorian results relative to the
earlier continuous-phase implementation. Use fresh initialization and new
baselines rather than assuming compatibility with those Gregorian restarts.

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
| MOSART stream heating | Consume ELM's corrected exported zenith angle. |

The coupled driver's equation-of-time setting is authoritative in EAMxx,
including when it is disabled. Leave
process-local `use_equation_of_time` unset to inherit it; an explicitly
contradictory process setting is rejected. Standalone EAMxx defaults on;
set `use_equation_of_time: false` in each participating RRTMGP/MAM process's
parameters to opt out. Keep those standalone process settings consistent. EAMxx's exact
coupled solar scheduling requires unsubcycled radiation and enclosing groups.

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
including for NO_LEAP. Disabled numerical branches retain their original
arithmetic; legacy bit-for-bit comparisons require an explicit opt-out. Full-model
bit-for-bit, component agreement, forcing-energy, and continuous-versus-restart
checks still require a supported machine; local compiled procedure checks are
not certification. Use identical orbital settings in comparisons, including
EAMxx's legacy NO_LEAP process-local orbit versus the enabled driver context.

DATM retains its existing forcing-date selection, including repeated February 28
data on February 29 when reading NO_LEAP forcing. Gregorian solar geometry uses
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
