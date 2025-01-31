(omega-user-tracers)=

# Tracers

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents).

## Updating tracer defintions in `TracerDefs.inc` file

To manage tracer definitions, users will update the `TracerDefs.inc` file
located in the `omega/src/ocn` directory of E3SM. The file contains
Tracer index variables and `defineAllTracers` C++ function defintion that
contains the calls to `define` function per each tracer as shown below:

```c++
inline static I4 IndxTemp             = Tracers::IndxInvalid;
inline static I4 IndxSalt             = Tracers::IndxInvalid;
inline static I4 IndxMyBGCTracer      = Tracers::IndxInvalid;

// Tracer definitions packaged in a defineAllTracers function
static void defineAllTracers() {

   define("Temp",                            ///< [in] Name of tracer
          "Potential Temperature",           ///< [in] Long name or description
          "degree_C",                        ///< [in] Units
          "sea_water_potential_temperature", ///< [in] CF standard Name
          -273.15,                           ///< [in] min valid field value
          100.0,                             ///< [in] max valid field value
          1.e33,                             ///< [in] value for undef entries
          IndxTemp);                         ///< [out] (optional) static index

   define("Salt", "Salinity", "psu", "sea_water_salinity", 0.0, 50.0, 1.e33,
          IndxSalt);
   define("Debug1", "Debug Tracer 1", "none", "none", 0.0, 100.0, 1.e33);
   define("Debug2", "Debug Tracer 2", "none", "none", 0.0, 100.0, 1.e33);
   define("Debug3", "Debug Tracer 3", "none", "none", 0.0, 100.0, 1.e33);
}
```

To add a new tracer, simply call the `define` function with the appropriate
arguments. Index argument is optional one that allows to access the tracer
data using the given tracer index variable.

## Selecting tracers using YAML configuration file

Note that not all tracers defined in `TracerDefs.inc` will be used during
a simulation. To select tracers and groups of tracers, users will configure
YAML files in OMEGA, such as `Default.yml` in the `omega/configs` directory,
as shown below:

```yaml
omega:
  Tracers:
    Base: [Temp, Salt]
    Debug: [Debug1, Debug2, Debug3]
    [other individual tracers or groups as needed]
```

In the above example, two tracer groups (Base and Debug) are selected. The
Base group includes the `Temp` and `Salt` tracers, while the Debug group
includes `Debug1`, `Debug2`, and `Debug3`.
