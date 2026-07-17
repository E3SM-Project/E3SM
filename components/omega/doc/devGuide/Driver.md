(omega-dev-driver)=

## Driver

The Omega Driver runs Omega either as a standalone model or as the ocean
component of E3SM (currently via the MCT coupler). In both modes execution
is divided into init, run, and finalize phases, but the coupled driver
splits init into two phases and calls run once per coupling interval
instead of once for the whole simulation. All entry points below are
declared in `OceanDriver.h`.

### Standalone driver

`ocnInit` reads `omega.yml`, calls `initOmegaModules` to initialize every
Omega module (including the default `TimeStepper`, which owns the model
`Clock` and `EndAlarm`), reads the `InitialState` or `RestartRead` IOStream,
and updates halo/host arrays for the resulting state:
```c++
int ocnInit(MPI_Comm Comm);
```
`ocnRun` advances the model, calling the default `TimeStepper`'s `doStep`
once per ocean time step, until the `EndAlarm` attached to the model
`Clock` rings:
```c++
int ocnRun(TimeInstant &CurrTime);
```
`ocnFinalize` cleans up all Omega objects; `CurrTime` is passed in case a
final restart needs to be written:
```c++
int ocnFinalize(const TimeInstant &CurrTime);
```
The standalone `main` (`src/drivers/standalone/OceanDriver.cpp`) wraps
these three calls with `MPI_Init`/`Kokkos::initialize` and their shutdown
counterparts.

### Coupled driver

The coupled driver (`src/drivers/coupled/`) runs Omega as the `ocn`
component of an MCT-based E3SM system. `ocn_comp_mct.F90` is the MCT cap,
exposing `ocn_init_mct`/`ocn_run_mct`/`ocn_final_mct` to the coupler. It
calls `bind(c)` interfaces declared in `omega_f2cxx_interface.F90` and
implemented in `omega_cxx2f_interface.cpp`, which call the C++ entry
points below.

#### Init

Init is split into two phases because the coupler needs Omega's
decomposition (`NCellsOwned`) before it can size and allocate the MCT
attribute vectors (`x2o`/`o2x`) that `SfcCoupling` attaches to:
```c++
int ocnInit1(MPI_Comm Comm, const int OcnId, const std::string &ConfigFile,
             const std::string &LogFile, const StartType StartType,
             const TimeInitParams &TimeParams,
             const CouplingInitParams &CouplingParams);
int ocnInit2(const Real *CplToOcnData, Real *OcnToCplData);
```
`ocnInit1` initializes every Omega module, using the coupled overload of
`initOmegaModules(Comm, TimeParams, CouplingParams)`, which also calls
`SfcCoupling::init`. It then reads the `InitialState` or `RestartRead`
IOStream based on `StartType`:
```c++
enum class StartType { StartUp, Continue, Branch };
```
converted from the coupler's integer start-type code by
`safeIntToStartType`. On restart/branch, the simulation time is taken from
the restart file rather than the coupler (which only knows the case start
time). On a cold start, `ocnInit1` advances the model `Clock` to the first
coupling-interval boundary so it is in sync with the coupler's clock.

`ocnInit2` runs once the coupler has allocated `x2o`/`o2x`: it attaches
them to `SfcCoupling` (`SfcCoupling::attachData`), does the initial
export/import exchange with the coupler, and updates halo/host arrays.

#### Run

The coupled overload of `ocnRun` advances the model one coupling interval,
rather than to the end of the simulation: it imports coupler fields into
`Forcing` at the start of the interval, steps until `SfcCoupling`'s
`CouplingAlarm` rings (updating `SfcCoupling`'s export accumulators every
ocean time step), and exports the accumulated fields back to the coupler.
`WriteRestart` is set by the coupler (from its own restart alarm) to force
a restart write at the end of the interval:
```c++
int ocnRun(TimeInstant &CurrTime, bool WriteRestart);
```

#### Finalize

`ocnFinalize` is shared with the standalone driver; it additionally clears
all `SfcCoupling` instances.

#### Field name/index mapping

`omega_cpl_indices.F90` builds the import (`x2o`) and export (`o2x`) field
name and coupler-column-index arrays from a dummy MCT attribute vector,
before Omega (and its `HorzMesh`) exists. `omega_ocn_init1` passes these
arrays through to `omega_cxx2f_interface.cpp`'s `buildFieldIndexMap`, which
builds the `std::map<std::string, int>` (`ImportIdx`/`ExportIdx`) required
by `CouplingInitParams`; see [Surface Coupling](#omega-dev-sfc-coupling)
for how those maps are used.

#### Mesh/decomposition queries

Before the coupler-owned buffers exist, the Fortran bridge needs Omega's
decomposition to build its MCT `gsMap` and domain. These bridge-only
queries expose that information and have no dependency on `SfcCoupling`:
```c++
int omega_get_ncells_local();
int omega_get_ncells_global();
void omega_get_index_to_cell_id(int *CellID);
void omega_get_lonlat_cell(double *LonCell, double *LatCell);
void omega_get_area_cell(double *AreaCell);
```
