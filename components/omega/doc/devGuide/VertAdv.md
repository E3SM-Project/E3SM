(omega-dev-vert-adv)=

## Vertical Advection

The `VertAdv` class contains methods and arrays for calculating vertical
velocities and the tendencies of thickness, horizontal velocity, and tracers due
to vertical advection.

### Initialization

The default `VertAdv` instance is created by calling `VertAdv::init()` method.
The default `HorzMesh`, `VertCoord`, and `Tracers` objects must be initialized
first. The default instance is retrieved by:
```
VertAdv *DefVertAdv = VertAdv::getDefault();
```

Additional instances can be created with the `create` method, which calls the
constructor and places the new instance in a map identified with the label
`Name`:
```
VertAdv *VertAdv::create(const std::string &Name, ///< [in] name for new VertAdv
                         const HorzMesh *Mesh,    ///< [in] associated HorzMesh
                         const VertCoord *VCoord, ///< [in] associated VertCoord
                         Config *Options          ///< [in] configuration options
)
```
This instance is retrieved with the get method:
```
VertAdv *NewVertAdv = VertAdv::get("Name");
```
The constructor reads configuration options, stores some info from the
`HorzMesh`, `VertCoord`, and `Tracers` objects, allocates needed arrays,
and defines `Field` metadata.

### Variables

The following arrays are stored as members of the `VertAdv` class.

| Variable Name | Type | Dimensions |
| ------------- | ---- | ---------- |
| VerticalVelocity | Real | NCellsSize, NVertLayersP1 |
| TotalVerticalVelocity | Real | NCellsSize, NVertLayersP1 |
| VertFlux | Real | NTracers, NCellsSize, NVertLayersP1 |
| LowOrderVertFlux | Real | NTracers, NCellsSize, NVertLayersP1 |

The `VerticalVelocity` represents the raw vertical pseudovelocity through the
top interfaces of cell layers, computed from the divergence of the horizontal
velocity. The `TotalVerticalVelocity` is the transport velocity and includes
corrections applied to the `VerticalVelocity`. The `VertFlux` and
`LowOrderVertFlux` store tracer fluxes at layer interfaces. The specific
numerical algorithms to compute these fluxes are chosen by the user through
configuration options.

### Methods

The `VerticalVelocity` and `TotalVerticalVelocity` arrays are computed by
calling the `computeVerticalVelocity` method:
```
VertAdv::computeVerticalVelocity(NormalVelocity, FluxPseudoThickEdge);
```
This method takes as input the `NormalVelocity` field from the `OceanState`
object, and the `FluxPseudoThickEdge` field from the `AuxiliaryState`. At
present, `TotalVerticalVelocity` is equivalent to `VerticalVelocity`;
additional corrections will be added in subsequent updates. The
`TotalVerticalVelocity` array is used by each of the tendency methods,
therefore `computeVerticalVelocity` must be called before any tendency
computations.

There are three methods for computing vertical advection tendencies of
thickness, horizontal velocity and tracers that are called from the time
stepper: `computeThicknessVAdvTend`, `computeVelocityVAdvTend`, and
`computeTracerVAdvTend`. Each method updates the corresponding tendency array
in place (inout).

Only the thickness tendency array is passed to the `computeThicknessVAdvTend`
method:
```
VertAdv::computeThicknessVAdvTend(ThickTend);
```

The `computeVelocityVAdvTend` method requires the `NormalVelocity` and
`FluxPseudoThickEdge` fields, in addition to the velocity tendency array:
```
VertAdv::computeVelocityVAdvTend(VelTend, NormalVelocity, FluxPseudoThickEdge);
```

The tracer vertical advection tendency depends on the configured settings.
The `computeTracerVAdvTend` method takes as arguments the full 3D array of
tracers, a thickness array (selected based on the configured advection
algorithm) and a `TimeInterval` representing the time step, along with the
tracer tendency array:
```
VertAdv::computeTracerVAdvTend(TracerTend, TracerArray, Thickness, TimeStep);
```
For the standard advection algorithm, `PseudoThickness` from the `OceanState` is
passed as the thickness argument; for the FCT algorithm, the `ProvPseudoThickness`
from the `AuxiliaryState` is used instead.

This method first calls `computeVerticalFluxes` to compute the tracer fluxes at
layer interfaces, using the order of accuracy specified in the configuration
file. It then dispatches to the selected advection algorithm (standard or FCT)
to compute the tracer tendencies.
