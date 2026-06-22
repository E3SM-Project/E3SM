(omega-design-ocn-coupler)=
# Surface Coupling

## 1 Overview

The `SurfaceCoupling` class manages all fields exchanged between Omega and
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

The class must own all coupling fields as `Array1DReal` arrays indexed
over ocean cells, consistent with other Omega field storage.

### 2.2 Requirement: Clean driver boundary

The driver interface layer (Fortran/C++ bridge) must handle all
driver-specific operations (index lookups, attribute vector layouts,
MCT/MOAB protocol). The C++ class receives and returns plain arrays of
doubles with no direct dependencies on MCT, MOAB, or other coupling
infrastructure.

### 2.3 Requirement: Field system integration

All coupling fields must be registerable with Omega's Field system for
diagnostic output.

### 2.4 Requirement: Per-coupling-interval operations

The class must provide methods to import fields, pass the received
import fields to the appropriate forcing class member variables,
accumulate export quantities each ocean timestep, export
accumulated fields, and reset accumulators.

### 2.5 Requirement: Driver architecture agnosticism

The class must support both MCT and MOAB layouts without conditional
compilation. The data layout (MCT column-major vs MOAB field-major) is
configured at runtime through the bridge layer.

### 2.6 Requirement: Coupling interval tracking

The class must track the coupling interval using Omega's `Alarm` system,
creating an internal periodic alarm at initialization to manage when
per-coupling-interval operations should occur.

### 2.7 Requirement: Coupler conversion layer

Private conversion methods must handle unit conversions and state
variable transformations between Omega's internal representation and
coupler expectations (e.g., conservative temperature ↔ in situ
temperature).

### 2.8 Requirement: Extensible design

Adding new coupling fields must require only localized changes: add
array members and update import/export logic. No restructuring of class
interfaces or method signatures should be needed.

## 3 Algorithmic Formulation

No algorithms are required.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

An enum class specifies the coupled driver layout:
```cpp
enum class CouplingLayout { MCT, MOAB };
```

The initial implementation will not have any runtime configuration options.

#### 4.1.2 Class/structs/data types

`SurfaceCoupling` lives in `src/ocn/` alongside `AuxiliaryState` and
`OceanState`. It has no direct driver dependencies. The bridge layer in
`src/drivers/coupled/` is the only code that calls it.

```cpp
class SurfaceCoupling {
 public:
   // Import fields (x2o): coupler → ocean
   Array1DReal SurfaceStressZonal;      ///< Foxx_taux  [N m⁻²]
   Array1DReal SurfaceStressMeridional; ///< Foxx_tauy  [N m⁻²]
   Array1DReal SWHeatFlux;              ///< Foxx_swnet [W m⁻²]
   Array1DReal LatentHeatFlux;          ///< Foxx_lat   [W m⁻²]
   Array1DReal SensibleHeatFlux;        ///< Foxx_sen   [W m⁻²]
   Array1DReal LWHeatFluxUp;            ///< Foxx_lwup  [W m⁻²]
   Array1DReal LWHeatFluxDown;          ///< Faxa_lwdn  [W m⁻²]
   Array1DReal SeaIceHeatFlux;          ///< Fioi_melth [W m⁻²]
   Array1DReal IcebergHeatFlux;         ///< Fioi_bergh [W m⁻²]
   Array1DReal SnowFlux;                ///< Faxa_snow  [kg m⁻² s⁻¹]
   Array1DReal RainFlux;                ///< Faxa_rain  [kg m⁻² s⁻¹]
   Array1DReal EvaporationFlux;         ///< Foxx_evap  [kg m⁻² s⁻¹]
   Array1DReal SeaIceFreshWaterFlux;    ///< Fioi_meltw [kg m⁻² s⁻¹]
   Array1DReal IcebergFreshWaterFlux;   ///< Fioi_bergw [kg m⁻² s⁻¹]
   Array1DReal SeaIceSalinityFlux;      ///< Fioi_salt  [kg m⁻² s⁻¹]
   Array1DReal RiverRunoffFlux;         ///< Foxx_rofl  [kg m⁻² s⁻¹]
   Array1DReal IceRunoffFlux;           ///< Foxx_rofi  [kg m⁻² s⁻¹]
   Array1DReal IceFraction;             ///< Si_ifrac   [-]
   Array1DReal SeaIcePressure;          ///< Si_bpress  [Pa]
   Array1DReal AtmosphericPressure;     ///< Sa_pslv    [Pa]

   // Export fields (o2x): ocean → coupler
   Array1DReal SurfaceTemperature;        ///< So_t    [°C]
   Array1DReal SurfaceSalinity;           ///< So_s    [g kg⁻¹]
   Array1DReal SurfaceVelocityZonal;      ///< So_u    [m s⁻¹]
   Array1DReal SurfaceVelocityMeridional; ///< So_v    [m s⁻¹]
   Array1DReal SSH;                       ///< So_ssh  [m]
   Array1DReal SSHGradientZonal;          ///< So_dhdx [m m⁻¹]
   Array1DReal SSHGradientMeridional;     ///< So_dhdy [m m⁻¹]
   Array1DReal SeaIceFormationHeat;       ///< Fioo_q      [W m⁻²]
   Array1DReal FrazilIceMass;             ///< Fioo_frazil [kg m⁻² s⁻¹]
   Array1DReal FreshWaterHeatFlux;        ///< Faoo_h2otemp[W m⁻²]

   // Public methods
   int importFromCoupler(const Real *x2oData, int NFields, int NCells);
   int exportToCoupler(Real *o2xData, int NFields, int NCells);
   int applyImportFields(AuxiliaryState *AuxState);
   int accumulateExportFields(const OceanState *State,
                             const Array3DReal &TracerArray,
                             int TimeLevel);
   void resetAccumulators();

 private:
   int NAccumSteps = 0;
   Alarm CouplingAlarm;
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;
   CouplingLayout DataLayout;

   static Real potTempToConservTemp(Real PotTemp, Real Salinity, Real Pressure);
   static Real conservTempToPotTemp(Real ConsTemp, Real Salinity, Real Pressure);
};
```

### 4.2 Methods

#### 4.2.1 Creation

The `create` method allocates import/export arrays sized to `NCellsOwned`, stores
the field index maps and `DataLayout`, creates the coupling `Alarm`, and calls
`registerFields`. After creation, `create` internally calls
`accumulateExportState` once to populate the initial export state.

```cpp
SurfaceCoupling *SurfaceCoupling::create(
   const std::string &Name, int OcnId, const HorzMesh *Mesh,
   const std::map<std::string, int> &ImportIdx,
   const std::map<std::string, int> &ExportIdx,
   const TimeInterval &CplInterval, CouplingLayout Layout);
```

#### 4.2.2 Retrieval

```cpp
SurfaceCoupling *SurfaceCoupling::getDefault();
SurfaceCoupling *SurfaceCoupling::get(const std::string &Name);
```

Other portions of the code can inquire whether running in coupled mode by
doing:

```cpp
if (!SurfaceCoupling::getDefault) {
    // no instance of SurfaceCoupling, therefore running in standalone mode
}
```

#### 4.2.3 Import, Apply, Accumulate, Export

`importFromCoupler` unpacks the driver array into typed member arrays using
the name→column-index map. The stride calculation is layout-dependent: MCT
uses `col*NCells + i`, MOAB uses `i*NFields + col`. Unit conversions (e.g.,
in situ → conservative temperature) and corrections (e.g. forcing shortwave
radiation to be positive) are applied inline. The corrections are applied to
address numerical issues, not and physical one, in the off chance monotonicity
is not preserved during remapping.

`applyImportFields` copies imported fields from the `SurfaceCoupling` class
into the appropriate arrays in the `Forcing` class. This method will be called
directly after `importFromCoupler` is called, but outside of the coupled loop.

`accumulateExportFields` is called *within* the coupled loop, and adds the
current-step values to running sums each ocean timestep. The accumulation
step in agnostic of whether the fields are passed back to the coupler as
sums or averages.

`exportToCoupler` packs export arrays into the driver array. This function is
called directly after the coupled loop. State fields are divided by
`NAccumSteps` (interval mean); flux fields are packed directly (interval total);
and `SSH` is passed as an instantaneous field.

`resetAccumulators` zeroes all export arrays and sets `NAccumSteps = 0`.

A rough sketch of how/when these function will be called within `OcnRun` looks
like:

```cpp
// fetch default OceanState, TimeStepper, and SurfaceCoupler
OceanState *DefOceanState   = OceanState::getDefault();
TimeStepper *DefTimeStepper = TimeStepper::getDefault();
SurfaceCoupler *DefSurfaceCoupler = SurfaceCoupler::getDefault();

// get coupling alarm
Alarm *CouplingAlarm = DefTimeStepper->getCouplingAlarm();

// these two could be wrapped into a single call
DefSurfaceCoupler->importFromCoupler(...);
DefSurfaceCoupler->applyImportFields(...);

while (Err == 0 && !(CouplingAlarm->isRinging())) {

   DefTimeStepper->doStep(DefOceanState, SimTime);

   DefSurfaceCoupler->accumulateExportFields(...);

}

// if coupler tells us, force write restate stream

DefSurfaceCoupler->exportToCoupler(...);

// Reset accumulators
DefSurfaceCoupler->resetAccumulators(...);
```

#### 4.2.4 Conversion Methods

Temperature conversion between conservative (Omega internal) and in situ
(coupler expectation) is handled as private static methods. Similar handling
will occur for salinity with conversion between absolute (Omega internal) and
practical (MPAS-Seaice expectation).

#### 4.2.5 Destruction and removal

```cpp
void SurfaceCoupling::erase(const std::string &Name);
void SurfaceCoupling::clear();
```

### 4.3 Driver Interface Bridge

The bridge layer in `src/drivers/coupled/` provides `extern "C"` entry points
called from the Fortran driver (`ocn_comp_mct.F90`). The coupling interval
and ocean timestep are independent: the driver calls `ocn_run_mct` once per
coupling interval; the ocean *can* take multiple timesteps within the coupling
interval.

#### 4.3.1 Name-based Field Index Mapping

At startup, the Fortran bridge populates parallel arrays of Omega field names
and driver column indices, then passes them to `omega_ocn_init`. C++ builds
`std::map<std::string, int>` objects for import/export lookups. This approach
is driver-agnostic: MCT and MOAB differ only in how column indices are
obtained. Field names are passed as fixed-length null-terminated character
arrays (standard Fortran-C interop), parsed by stepping through the flat
array at stride `FieldNameLen`.

### 4.4 Interval Accumulation

Export state fields are accumulated each ocean timestep and divided by
`NAccumSteps` at export time (interval mean). Export flux fields are
accumulated and packed directly without dividing (interval total).
`resetAccumulators` is called by `omega_ocn_export` after packing.

### 4.5 Extensibility

To add a new coupling field:
1. Add an `Array1DReal` member to `SurfaceCoupling.h`
2. Add the field registration in `registerFields()`
3. Add an unpack/pack call in `importFromCoupler`/`exportToCoupler`
4. Add the name/index entry in the Fortran bridge file
5. Fill the field in `applyImportFields()` or `accumulateExportFields()`

No changes to method signatures or ordering constraints are needed.

## 5 Verification and Testing

### 5.1 Test: Import/export round-trip

Create synthetic coupling arrays with known values, call
`importFromCoupler`, verify member arrays match expected values. Pack
export arrays and verify the output matches. Test both MCT and MOAB
layouts.
  - tests requirement 2.1, 2.2

### 5.2 Test: Interval averaging for state fields

Accumulate the same state values over multiple simulated timesteps,
export, and verify the output equals the mean.
  - tests requirement 2.4

### 5.3 Test: Interval summation for flux fields

Accumulate flux values over multiple timesteps, export, and verify the
output equals the total sum (not divided).
  - tests requirement 2.4

### 5.4 Test: Name-keyed field mapping

Build the name→index map with different fill orders and verify
import/export produces correct results regardless of order. Request an
unknown field name and verify error logging.
  - tests requirement 2.2, 2.8

### 5.5 Test: Coupling alarm

Create a `SurfaceCoupling` with a known coupling interval. Advance the
simulation clock and verify the coupling alarm rings at the expected
interval boundaries and does not ring between them.
  - tests requirement 2.6

### 5.6 Test: Coupling loop integration

Simulate a complete coupling interval: import, run multiple ocean
timesteps with accumulation, export, and reset. Verify the full
import → accumulate → export → reset cycle produces correct results and
that the accumulators are properly zeroed for the next interval.
  - tests requirement 2.4, 2.6
