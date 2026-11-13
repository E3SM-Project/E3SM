(omega-user-tendencies)=

## Tendencies

The `Tendencies` class provides a container for the [tendency terms](#omega-user-tend-terms) in OMEGA.
Upon creation of an `Tendencies` instance, these functors are initialized and arrays for the
accumulated tendencies are allocated.
There are no user-configurable options beyond those for the tendency term functors.

Frazil is one of the configurable tendency contributions. It is controlled by
`Omega.Tendencies.FrazilTendencyEnable` together with the `Omega.Frazil`
configuration block (`FrazilType`, `LayerMassFracMax`, `Phi`, `DepthLimit`, and
`ConservationCheck`). For operational guidance and options, see
[Frazil](Frazil.md).
