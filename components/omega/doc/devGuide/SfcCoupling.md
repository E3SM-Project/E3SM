(omega-dev-sfc-coupling)=

# Surface Coupling

The `SfcCoupling` class manages the variables exchanged to (`o2x`) and from
(`x2o`) the coupler. It handles import/export of the raw coupler data,
application of imported fields to the [`Forcing`](#omega-dev-forcing) object,
and accumulation of export fields over the coupling interval. It is possible
to have multiple `SfcCoupling` instances; every instance has a name and is
tracked in a static C++ map called `AllSfcCoupling`.

## Initialization

The static method:
```c++
OMEGA::SfcCoupling::init(CouplingInitParams);
```
initializes the default `SfcCoupling` given a `CouplingInitParams` struct
(number of import/export fields, name-to-index maps, coupling time step, and
coupling `Layout`, i.e. `MCT` or `MOAB`) supplied by the coupler. A pointer to
the default instance can be retrieved at any time using:
```c++
OMEGA::SfcCoupling* DefSfcCoupling = OMEGA::SfcCoupling::getDefault();
```

## Creation of non-default surface coupling objects

A non-default instance can be created with the static `create` method,
which returns a pointer to the newly created object:
```c++
OMEGA::SfcCoupling* NewSfcCoupling = OMEGA::SfcCoupling::create(
    Name, Mesh, NImportFields, NExportFields, ImportIdxMap, ExportIdxMap,
    Stepper, CouplingTimeStep, Layout);
```
Given its name, a pointer to a named instance can be obtained at any time:
```c++
OMEGA::SfcCoupling* NewSfcCoupling = OMEGA::SfcCoupling::get(Name);
```

## Data exchange

Before import/export, unmanaged views must be attached to the raw coupler
data pointers:
```c++
SfcCoupling.attachData(CplToOcnData, OcnToCplData);
```
To copy data from the coupler into the `CplToOcn` fields, and from the
`OcnToCpl` fields back out to the coupler:
```c++
SfcCoupling.importFromCoupler();
SfcCoupling.exportToCoupler();
```
Imported fields are applied to a `Forcing` object with:
```c++
SfcCoupling.applyImportFields(ForcingPtr);
```
Export fields are accumulated (running average) each ocean time step with:
```c++
SfcCoupling.updateExportFields(State, TracerArray);
```

## Units

`OcnToCplFields`'s averaged temperature is Conservative Temperature (deg C)
on device (`AvgSfcTemperature`), invariant at all times. `copyToHost()`
converts it to in-situ temperature in Kelvin (via TEOS-10, when in use)
into the host mirror `AvgSfcTemperatureH`, once per coupling interval.
All other averaged fields keep the same units and value on host and device.

To keep this invariant, `OcnToCplFields`'s averaged device arrays are
private: they can only be written via `updateAverages()` (called each
ocean time step from `updateExportFields()`) and read via `copyToHost()`.

## Removal of surface coupling objects

To erase a specific named instance use `erase`:
```c++
OMEGA::SfcCoupling::erase(Name);
```
To clear all instances do:
```c++
OMEGA::SfcCoupling::clear();
```
