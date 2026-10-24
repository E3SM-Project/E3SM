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

Append `-v` to display the five test method names and their results.

| Test method | Checks |
| --- | --- |
| `test_epoch_defaults_and_eot_default` | 1850, 2010, and 20TR epoch defaults, fixed-year default, and default-on EOT under both calendars and drivers |
| `test_profiles` | All eight profile namelists and their orbital-year/EOT settings for both drivers |
| `test_explicit_epoch_and_legacy_alias` | Explicit orbital-year/alignment overrides and the `variable_year` alias |
| `test_suite_modifications` | Suite counts, unique entries, and resolution of composed modifications |
| `test_restart_pointers` | Calendar-specific restart timestamps for both short boundaries and the four-year runs |

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
