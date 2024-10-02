(omega-user-tracers)=

# Tracers

Tracers refer to either heat (temperature) or material carried by a fluid
parcel (e.g., salt, chemical, or biological constituents).

To manage tracer definitions, users will update the `TracerDefs.inc` file
located in the `omega/src/ocn` directory of E3SM. The file contains
`define` C++ function calls, each of which defines a tracer as shown below:

```c++
static I4
define(const std::string &Name,        ///< [in] Name of tracer
       const std::string &Description, ///< [in] Long name or description
       const std::string &Units,       ///< [in] Units
       const std::string &StdName,     ///< [in] CF standard Name
       const Real ValidMin,            ///< [in] min valid field value
       const Real ValidMax,            ///< [in] max valid field value
       const Real FillValue,           ///< [in] value for undef entries
       I4 &Index = IndxInvalid         ///< [out] (optional) index value
);
```

To add a new tracer, simply call the `define` function with the appropriate
arguments. Index argument is optional one that allows to access the tracer
data using the given tracer index variable.

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
