(omega-user-vert-adv)=

## Vertical Advection

### Overview

In Omega, vertical advection represents the transport of mass, momentum, and
tracer quantities along the vertical coordinate due to vertical motion within a
water column. In the context of Omega's
[governing equations](omega-design-governing-eqns-omega1), vertical advection
contributes to the tendencies of pseudo-thickness, horizontal velocity, and
tracer fields. Vertical motion is inferred from the divergence of horizontal
flow through the continuity equation and is necessary for accurately
representing vertical transport in three-dimensional ocean simulations.

The `VertAdv` class contains variables and functions for computing the vertical
velocity and the tendencies of `PseudoThickness`, `NormalVelocity`, and
`Tracers`. The algorithms used are determined by options specified in the
configuration file.

### Configuration Options

Within the `Tendencies` group of the configuration file, flags are used to
enable or disable the contribution of vertical advection to the tendencies of
thickness, velocity, and tracers:
| Configuration name | Options |
| ------------------ | ------- |
| ThicknessVertAdvTendencyEnable | true / false |
| VelocityVertAdvTendencyEnable | true / false |
| TracerVertAdvTendencyEnable | true / false |

Within the `Advection` group, options are provided to configure the algorithms
used to compute the tracer tendencies:
| Configuration name | Description | Options |
| ------------------ | ----------- | ------- |
| VerticalTracerFluxOrder | order of accuracy used for the tracer flux calculation | 2, 3, or 4 |
| VerticalTracerFluxLimiterEnable | selects the tracer tendency algorithm | true / false |
| Coef3rdOrder | coefficient used for blending 3rd- and 4th-order fluxes when `VerticalTracerFluxOrder == 3` | 0.0 - 1.0 |

For `VerticalTracerFluxLimiterEnable`, `true` enables flux-corrected transport,
and `false` selects the standard algorithm. For `Coef3rdOrder`, `1` corresponds
to purely 3rd-order, `0` purely 4th-order.

The `VertAdv` class provides four fields that can be written to output by
adding them to the contents of an output stream in the configuration file:
`VerticalPseudoVelocity`, `TotalVerticalPseudoVelocity`, `VertFlux`, and
`LowOrderVertFlux`. These fields collectively make up the `VertAdv` output group.
