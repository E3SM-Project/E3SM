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

## Coupled energy/heat correction (aka energy fixer)

The export field `Faoo_h2otemp` carries the extra energy flux associated with
freshwater mass exchange, in W m$^{-2}$. It communicates the energy already
included in the ocean's surface heat budget so that the coupled system can
account for the corresponding energy as an atmospheric adjustment. It is not an
additional ocean temperature tendency or a total heat-flux export.

The instantaneous field is owned by
`Forcing::TracerForcing.ExtraEnergyCell`. Its calculation, physical assumptions,
and reset behavior are described in the
[forcing guide](#omega-dev-forcing-extra-energy).

## Units

`OcnToCplFields`'s averaged temperature is Conservative Temperature (deg C)
on device (`AvgSfcTemperature`), invariant at all times. `copyToHost()`
converts it to in-situ temperature in Kelvin (via TEOS-10, when in use)
into the host mirror `AvgSfcTemperatureH`, once per coupling interval.
All other averaged fields keep the same units and value on host and device.

To keep this invariant, `OcnToCplFields`'s averaged device arrays are
private: they can only be written via `updateAverages()` (called each
ocean time step from `updateExportFields()`) and read via `copyToHost()`.

## Frazil Sign Conventions

Frazil exports follow the existing MPAS-O / E3SM convention:

| Coupler field | Units | Positive direction | Description |
| --- | --- | --- | --- |
| `Fioo_frazilm` | $\mathrm{kg\ m^{-2}\ s^{-1}}$ | Out of ocean | Frazil mass flux |
| `Fioo_frazils` | $\mathrm{kg\ m^{-2}\ s^{-1}}$ | Out of ocean | Frazil salt flux |
| `Fioo_frazilh` | $\mathrm{W\ m^{-2}}$ | Into ocean | Heat flux from frazil formation |

These mass and salt flux signs are an exception to the standard coupler
convention used by other components. Omega adopts this inverse convention
to interface with existing E3SM coupling plumbing. Standardizing this across
E3SM is tracked in [E3SM-Project/E3SM#8928](https://github.com/E3SM-Project/E3SM/issues/8928).

Because newly formed frazil has negative enthalpy, extracting frazil increases
ocean heat energy. `Fioo_frazilh` is exported as this positive-into-ocean heat
flux, which sea ice uses to reconstruct negative ice enthalpy.

`Fioo_q` ($\mathrm{W\ m^{-2}}$) is separate, nonpositive melt potential. It is
set to zero when frazil is exported; otherwise it represents energy available
at the surface to melt ice.

## Removal of surface coupling objects

To erase a specific named instance use `erase`:
```c++
OMEGA::SfcCoupling::erase(Name);
```
To clear all instances do:
```c++
OMEGA::SfcCoupling::clear();
```
