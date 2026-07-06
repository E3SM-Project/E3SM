(omega-design-ocn-coupler)=
# Surface Coupling

## 1 Overview

The `SfcCoupling` class manages all fields exchanged between Omega and
the coupled system (via MCT or MOAB drivers). It owns import and export
coupling field arrays, unpacks driver attribute-vector data into typed
Omega arrays, applies imported fields to internal forcing state,
accumulates export fields over each coupling interval, and repacks them
for the driver. The design cleanly separates driver-specific operations
(handled in a Fortran/C++ bridge layer) from the C++ class, is agnostic
to the underlying driver architecture, and is extensible for future
fields (BGC, waves, land-ice).

## 2 Requirements

### 2.1 Requirement: Own coupling field storage

The class must own all coupling fields as `Array1DReal`/`HostArray1DReal`
arrays indexed over ocean cells, consistent with other Omega field
storage.

### 2.2 Requirement: Clean driver boundary

The driver interface layer (Fortran/C++ bridge) must handle all
driver-specific operations (index lookups, attribute vector layouts,
MCT/MOAB protocol). The C++ class receives and returns plain arrays of
doubles with no direct dependencies on MCT, MOAB, or other coupling
infrastructure.

### 2.3 Requirement: Field system integration

All coupling fields must be registerable with Omega's Field system for
diagnostic output. Not yet implemented in the initial pass.

### 2.4 Requirement: Per-coupling-interval operations

The class must provide methods to import fields, pass the received
import fields to the appropriate forcing class member variables,
accumulate export quantities each ocean timestep, and export
accumulated fields. Accumulators are reset as part of the export step,
so no separate reset call is required in the run loop.

`applyImportFields` writes directly into `Forcing` member arrays, the
same arrays populated by file-based input in standalone runs. This
keeps `Forcing` agnostic to whether its data originates from the
coupler or from a file, and `SfcCoupling` itself has no user-facing
runtime configuration options; all configuration comes from the
coupler at initialization (see 4.1.1).

### 2.5 Requirement: Driver architecture agnosticism

The class must support both MCT and MOAB layouts without conditional
compilation. The data layout (MCT column-major vs MOAB field-major) is
configured at runtime through the bridge layer.

### 2.6 Requirement: Coupling interval tracking

The class must track the coupling interval using Omega's `Alarm` system,
creating an internal periodic alarm at initialization to manage when
per-coupling-interval operations should occur.

### 2.7 Requirement: Coupler conversion layer

Unit conversions and state variable transformations between Omega's
internal representation and coupler expectations (e.g., conservative
temperature → in situ temperature, absolute → practical salinity) must
be supported, and must be done once per coupling interval rather than
every ocean timestep, since evaluating the EOS polynomial is expensive.
Converted and unconverted representations must not both be reachable
from outside the class, to prevent unit confusion. Temperature
conversion is implemented (4.2.5); practical salinity conversion is a
TODO (see 2.8).

### 2.8 Requirement: Extensible design

Adding new coupling fields must require only localized changes: add
array members and update import/export logic. No restructuring of class
interfaces or method signatures should be needed.

## 3 Algorithmic Formulation

Export fields are accumulated on-device each ocean timestep using Welford's
online running-average algorithm, rather than a running sum divided at
export time:
```cpp
KOKKOS_INLINE_FUNCTION Real updateAverage(const Real OldAvg,
                                          const Real NewValue,
                                          const I4 NAccumSteps) {
   return OldAvg + (NewValue - OldAvg) / (NAccumSteps + 1);
}
```
This produces the same interval-mean result as a naive sum-then-divide, but
avoids growing partial sums. Instantaneous fields (e.g. `SSH`) bypass
accumulation and are read directly from their source at export time.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

An enum class specifies the coupled driver layout:
```cpp
enum class CouplingLayout { MCT, MOAB };
```

All initialization inputs are supplied by the coupler at runtime and
bundled into a single struct, rather than passed as separate arguments:
```cpp
struct CouplingInitParams {
   int NImportFields;
   int NExportFields;
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;
   TimeInterval CouplingTimeStep;
   CouplingLayout Layout;
};
```
`SfcCoupling` has no user-facing runtime configuration options.

#### 4.1.2 Class/structs/data types

`SfcCoupling` lives in `src/ocn/` alongside `AuxiliaryState` and
`OceanState`. It has no direct driver dependencies. The bridge layer in
`src/drivers/coupled/` is the only code that calls it.

Import (x2o) and export (o2x) fields are grouped into two small container
classes, rather than flat members directly on `SfcCoupling`. Each is
constructed with a name `Suffix` and `Mesh` so multiple named `SfcCoupling`
instances don't collide on Kokkos view labels:

```cpp
// x2o: Coupler to Ocean. Host-only; applyImportFields() copies to the
// device arrays owned by Forcing.
class CplToOcnFields {
 public:
   HostArray1DReal SfcStressZonal; ///< Foxx_taux  [N m⁻²]
   HostArray1DReal SfcStressMerid; ///< Foxx_tauy  [N m⁻²]

   CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh);
};

// o2x: Ocean to Coupler. Averaged fields keep a private device array
// (updated each ocean timestep, native Omega units) plus a public host
// mirror (converted units, packed into the driver buffer). Device
// arrays are private so external code cannot read them in native units
// and mistake them for the (unit-converted) host mirrors.
class OcnToCplFields {
 public:
   HostArray1DReal AvgSfcTemperatureH;     ///< So_t   [K], in situ approx
   HostArray1DReal AvgSfcSalinityH;        ///< So_s   [g kg⁻¹], TODO: practical
   HostArray1DReal AvgSfcVelocityZonalH;   ///< So_u   [m s⁻¹]
   HostArray1DReal AvgSfcVelocityMeridH;   ///< So_v   [m s⁻¹]
   HostArray1DReal InstSshCellH;           ///< So_ssh [m], instantaneous

   void updateFields(const OceanState *State, const Array3DReal &TracerArray,
                     I4 NAccumSteps, I4 NCellsOwned);
   void copyToHost();
   void resetFields();

   OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh);

 private:
   Array1DReal AvgSfcTemperature;      ///< [°C], conservative temperature
   Array1DReal AvgSfcSalinity;         ///< [g kg⁻¹], absolute salinity
   Array1DReal AvgSfcVelocityZonal;
   Array1DReal AvgSfcVelocityMerid;

   // Scratch buffer for the in situ/Kelvin conversion in copyToHost()
   Array1DReal InSituTempScratch;
};
```

`updateFields` (accumulation) and the conversion in `copyToHost` both
need direct access to the device arrays, so that logic lives on
`OcnToCplFields` itself rather than on `SfcCoupling`;
`SfcCoupling::updateExportFields` is now a thin wrapper that just calls
`OcnToCpl.updateFields(...)` and increments `NAccumSteps`.

Only the fields above are wired up in the initial implementation. The
remaining fields from the original field list are deferred to follow-on
work and added the same way (4.5): SW/LW heat fluxes, latent/sensible
heat, snow/rain/evaporation, sea-ice and iceberg heat/freshwater fluxes,
`SeaIceSaltFlux` (renamed from `SeaIceSalinityFlux`), river/ice runoff,
ice fraction, `WindSpeed10m` (needed for KPP), sea-ice/atmospheric
pressure, SSH gradients, sea-ice formation heat (to be renamed once
sign-definite, e.g. melt potential), frazil ice mass, and freshwater heat
flux.

```cpp
class SfcCoupling {
 public:
   std::string Name;
   I4 NCellsOwned;      ///< Number of cells owned by this task
   I4 NImportFields;    ///< Num of fields in the x2o pointer array
   I4 NExportFields;    ///< Num of fields in the o2x pointer array

   CplToOcnFields CplToOcn; ///< Coupler to Ocean (x2o)
   OcnToCplFields OcnToCpl; ///< Ocean to Coupler (o2x)
   Alarm CouplingAlarm;     ///< Alarm for the coupling interval

   // Unmanaged, strided views over the driver's raw x2o/o2x buffers
   Kokkos::View<const Real **, Kokkos::LayoutStride, Kokkos::HostSpace,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
       CplToOcnView;
   Kokkos::View<Real **, Kokkos::LayoutStride, Kokkos::HostSpace,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
       OcnToCplView;

   static SfcCoupling *create(const std::string &Name, const HorzMesh *Mesh,
                              int NImportFields, int NExportFields,
                              const std::map<std::string, int> &ImportIdx,
                              const std::map<std::string, int> &ExportIdx,
                              TimeStepper *Stepper,
                              const TimeInterval &CouplingTimeStep,
                              const CouplingLayout &Layout);
   static int init(const CouplingInitParams &Params);
   ~SfcCoupling();
   static void clear();
   static void erase(const std::string Name);
   static SfcCoupling *getDefault();
   static SfcCoupling *get(const std::string Name);
   I4 getNAccumSteps() const;

   void attachData(const Real *CplToOcnData, Real *OcnToCplData);
   void importFromCoupler();
   void exportToCoupler();
   void applyImportFields(Forcing *Forcing);
   void updateExportFields(const OceanState *State,
                           const Array3DReal &TracerArray);

 private:
   I4 NAccumSteps = 0;
   CouplingLayout Layout;
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;

   template <class View> auto ownedSubView(const View &V) const;
};
```

Temperature conversion (2.7) is implemented in `copyToHost` (4.2.5);
practical-salinity conversion is a TODO, so `AvgSfcSalinityH` is still
absolute salinity.

### 4.2 Methods

#### 4.2.1 Creation and Initialization

`init` retrieves the default `HorzMesh` and `TimeStepper`, validates that
the coupling interval is not shorter than, and is evenly divisible by, the
ocean timestep (aborting otherwise), then calls `create` to build the
default instance:
```cpp
static int SfcCoupling::init(const CouplingInitParams &Params);
```
`create` allocates the `CplToOcn`/`OcnToCpl` field containers sized to
`NCellsOwned`, stores the index maps and `Layout`, and constructs
`CouplingAlarm` on `Stepper`'s `Clock`:
```cpp
static SfcCoupling *SfcCoupling::create(
    const std::string &Name, const HorzMesh *Mesh, int NImportFields,
    int NExportFields, const std::map<std::string, int> &ImportIdx,
    const std::map<std::string, int> &ExportIdx, TimeStepper *Stepper,
    const TimeInterval &CouplingTimeStep, const CouplingLayout &Layout);
```

#### 4.2.2 Retrieval

```cpp
SfcCoupling *SfcCoupling::getDefault();
SfcCoupling *SfcCoupling::get(const std::string &Name);
```

Other portions of the code can inquire whether it is running in coupled
mode by doing:

```cpp
if (!SfcCoupling::getDefault()) {
    // no instance of SfcCoupling, therefore running in standalone mode
}
```

#### 4.2.3 Attaching coupler data

Before import/export, `attachData` wraps the driver's raw `x2o`/`o2x`
pointers in unmanaged, layout-strided Kokkos views, computing the
layout-dependent stride once rather than on every element access: MCT
lays fields out as `(NCellsOwned, NImportFields)` (field index strides
fastest), MOAB as `(NImportFields, NCellsOwned)` (cell index strides
fastest).
```cpp
void SfcCoupling::attachData(const Real *CplToOcnData, Real *OcnToCplData);
```

#### 4.2.4 Import, Apply, Update, Export

`importFromCoupler` looks up each field's column index from `ImportIdx`
and copies it out of `CplToOcnView` into the corresponding `CplToOcn`
array.

`applyImportFields` deep-copies `CplToOcn` arrays into the matching
`Forcing` arrays (e.g. `Forcing->SfcStressForcing.ZonalStressCell`),
restricted to the owned-cell subview since `SfcCoupling` has no halo
information; `Forcing` is responsible for the halo exchange. Called once
per coupling interval, directly after `importFromCoupler`, outside the
ocean timestep loop.

`updateExportFields` is called *within* the ocean timestep loop. It
updates each `OcnToCpl` running average in place using Welford's
algorithm (3), then increments `NAccumSteps`. Velocity currently uses a
placeholder constant (`1e-4`) pending vector reconstruction, to avoid
producing zero (and therefore divide-by-zero/infinite flux) velocities
downstream in the coupler; this is a known limitation to revisit before
production use.

`exportToCoupler` calls `OcnToCpl.copyToHost()`, which converts and
copies device arrays to their host mirrors (SSH is read directly from
`VertCoord`'s host SSH array, since it is instantaneous and not
accumulated), then packs the host mirrors into `OcnToCplView` at their
export indices, then resets `OcnToCpl` and `NAccumSteps` for the next
interval. Resetting is folded into `exportToCoupler`; there is no
separate `resetAccumulators` call in the run loop.

A rough sketch of how/when these functions are called within `OcnRun`:

```cpp
// fetch default OceanState, TimeStepper, and SfcCoupling
OceanState *DefOceanState   = OceanState::getDefault();
TimeStepper *DefTimeStepper = TimeStepper::getDefault();
SfcCoupling *DefSfcCoupling = SfcCoupling::getDefault();

DefSfcCoupling->attachData(CplToOcnData, OcnToCplData);

// these two could be wrapped into a single call
DefSfcCoupling->importFromCoupler();
DefSfcCoupling->applyImportFields(DefForcing);

while (Err == 0 && !(DefSfcCoupling->CouplingAlarm.isRinging())) {

   DefTimeStepper->doStep(DefOceanState, SimTime);

   DefSfcCoupling->updateExportFields(DefOceanState, TracerArray);

}

// if coupler tells us, force write restart stream

DefSfcCoupling->exportToCoupler(); // also resets accumulators
```

#### 4.2.5 Conversion Methods

Conversion happens in `copyToHost`, once per coupling interval rather
than every ocean timestep, since the EOS conversion polynomial is
expensive. Conservative temperature is converted to potential
temperature via a local `Teos10Eos` instance's `calcPtFromCt` (only
under `EosType::Teos10Eos`; otherwise passed through unconverted), then
shifted to Kelvin (`+ TkFrz`) as an in situ approximation. The result is
written to a private scratch buffer (`InSituTempScratch`) rather than
in place, so the device `AvgSfcTemperature` array driving the running
average (3) always stays in native (deg C, conservative) units. A
local `Teos10Eos` is constructed rather than reusing `Eos::getInstance()`
since `calcPtFromCt` needs no config-derived state, avoiding exposure of
`Eos` internals. Salinity conversion (absolute → practical) is a TODO;
`AvgSfcSalinityH` is currently absolute salinity, copied through
unconverted.

The device arrays backing the conversion (`AvgSfcTemperature`,
`AvgSfcSalinity`, ...) are private members of `OcnToCplFields`; only the
post-conversion host mirrors are public. This prevents external code
from reading the native-unit device arrays and confusing them with the
converted host mirrors -- a hazard motivated by `copyToHost` itself
needing to reach into `VertCoord::SshCell` (a raw device array) to
populate `InstSshCellH`, which showed how easy it is to mix up
device/host and native/converted units without this guard.

#### 4.2.6 Destruction and removal

```cpp
static void SfcCoupling::erase(const std::string Name);
static void SfcCoupling::clear();
```
The destructor does not detach `CouplingAlarm` from the shared `Clock`;
callers must ensure `SfcCoupling` instances outlive the `Clock`, consistent
with other `Alarm` owners in Omega.

### 4.3 Driver Interface Bridge

The bridge layer in `src/drivers/coupled/` is not yet implemented; it will
provide `extern "C"` entry points, called from the Fortran driver
(`ocn_comp_mct.F90`), that build a `CouplingInitParams` and call
`SfcCoupling::init`. The coupling interval and ocean timestep are
independent: the driver calls `ocn_run_mct` once per coupling interval;
the ocean can take multiple timesteps within it.

#### 4.3.1 Name-based Field Index Mapping

At startup, the Fortran bridge will populate parallel arrays of Omega
field names and driver column indices and pass them into
`CouplingInitParams::ImportIdx`/`ExportIdx`. This approach is
driver-agnostic: MCT and MOAB differ only in how column indices are
obtained.

### 4.4 Interval Accumulation

Export fields accumulate a running average each ocean timestep via
Welford's algorithm (3), in native Omega units; `SSH` is read directly
from its source at export time instead of being accumulated. Unit
conversion (4.2.5) is deferred to `copyToHost`/`exportToCoupler`, so it
runs once per interval rather than once per accumulation step.
`exportToCoupler` resets all accumulators after packing, so no separate
reset call is needed in the run loop. Interval summation (as opposed to
averaging), needed for flux-like fields, is not yet implemented; it can
reuse `NAccumSteps` (average x steps) or a dedicated accumulation path
once flux fields are added.

### 4.5 Extensibility

To add a new coupling field:
1. Add an `Array1DReal`/`HostArray1DReal` member to `CplToOcnFields` or
   `OcnToCplFields` in `SfcCoupling.h`, and initialize it in the
   corresponding constructor
2. Add an unpack/pack call in `importFromCoupler`/`exportToCoupler` (and,
   for export fields, an update in `updateExportFields`)
3. Add the name/index entry in the Fortran bridge file
4. Fill the field in `applyImportFields()` or `updateExportFields()`

No changes to method signatures or ordering constraints are needed.

## 5 Verification and Testing

### 5.1 Test: Import round-trip

Attach synthetic `x2o` buffers with known, cell-varying values, call
`importFromCoupler`, and verify the `CplToOcn` arrays match expected
values. Test both MCT and MOAB layouts.
  - tests requirement 2.1, 2.2, 2.5

### 5.2 Test: Apply imported fields to Forcing

Populate `CplToOcn` arrays directly, call `applyImportFields`, and verify
the owned-cell subview of the matching `Forcing` arrays matches.
  - tests requirement 2.1, 2.4

### 5.3 Test: Running-average accumulation and conversion

Advance the model clock through a known number of ocean timesteps within
one coupling interval, calling `updateExportFields` with cell- and
step-varying tracer values each timestep (using the coupling `Alarm` to
end the loop, as in the sketch in 4.2.4). Verify `NAccumSteps` matches
the step count and, after `copyToHost` (the only sanctioned way to read
the averages), the resulting host-mirror averages match the analytic
mean within round-off tolerance. `EosType` is forced to `constant` for
the test suite so `calcPtFromCt` is an identity, isolating the `+TkFrz`
Kelvin shift; salinity is checked unconverted (practical conversion is a
TODO).
  - tests requirement 2.4, 2.6, 2.7

### 5.4 Test: Export round-trip and reset

Populate `OcnToCpl` averages via `updateExportFields` (seeded to the
target value using `NAccumSteps == 0`), call `exportToCoupler`, and
verify the packed `o2x` buffer matches expected values (temperature
offset by `TkFrz`) for both MCT and MOAB layouts. Also verify `OcnToCpl`
host mirrors and `NAccumSteps` are zeroed afterward.
  - tests requirement 2.1, 2.4, 2.5, 2.7

### 5.5 Test: Object lifecycle

Create a non-default, named `SfcCoupling`, verify it can be retrieved
with `get`, erase it, and verify retrieval afterward fails.
  - tests requirement 2.8
