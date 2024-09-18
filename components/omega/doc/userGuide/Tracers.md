(omega-user-tracers)=

# Tracers

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents).

To manage tracer definitions, users will update the `TracerDefs.inc` file
located in the `omega/src/ocn` directory of E3SM. The file contains
`define` C++ function calls, each of which defines a tracer as shown below:

```c++
define(
       "Temp",                  // Name of variable
       "Potential Temperature", // Long name or description
       "degree_C",              // Units
       "sea_water_potential_temperature", // CF standard name
       -273.15,                 // Min valid field value
       100.0,                   // Max valid field value
       1.e33                    // Fill value for undefined entries
);
```

To add a new tracer, simply call the `define` function with the appropriate
arguments.

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
