(omega-dev-rayleigh-damping)=

# Rayleigh Damping

Rayleigh damping (see the [user guide](#omega-user-rayleigh-damping)) is not a separate class or tendency term. It is a diagonal contribution to the tridiagonal system solved by the implicit velocity vertical mixing in `VertMix::applyVelVertMixImplicit`.

## Implementation

The `VelVertMixSetupOnEdge` functor in `VertMix.h`, which sets up the tridiagonal system for each edge column, has two public members for Rayleigh damping, alongside the members for implicit bottom drag:

```c++
bool RayleighDampingEnabled; ///< Enable Rayleigh damping flag
Real RayleighDampingCoeff;   ///< Rayleigh damping coefficient (s^-1)
```

The constructor sets them to `false` and `0`. They are read from the `RayleighDamping` config group in `VertMix::init()`, which aborts if `DampingCoeff` is negative.

The tridiagonal system solved for the new normal velocity $Y_k$ in each active layer is

$$
-G_{k-1} Y_{k-1} + \left( G_{k-1} + G_k + H_k \right) Y_k - G_k Y_{k+1} = X_k .
$$

The backward-Euler Rayleigh term adds $\Delta t \, c_R \, \tilde{h}_k$ to the diagonal-only coefficient $H_k$, where $\tilde{h}_k$ is the pseudo-thickness of layer $k$ interpolated to the edge. With bottom drag, $H_k$ becomes

$$
H_k = \tilde{h}_k \left( 1 + \Delta t \, c_R \right)
+ \delta_{k,k_{\mathrm{bot}}} \, \Delta t \, C_D \, \frac{\rho}{\rho_0} \left| \mathbf{u} \right| .
$$

The off-diagonal coefficients $G_k$ and the right-hand side $X_k = \tilde{h}_k u^n_k$ are unchanged. The functor is only called for active layers, so damping is applied to exactly those layers. The contribution is guarded by `RayleighDampingEnabled`, so the disabled path is unchanged.

Because the damping is part of the implicit solve, `Tendencies` aborts at initialization if Rayleigh damping is enabled but `VelVertMixTendencyEnable` is `false`, as it already does for implicit bottom drag. `VertMix::init()` must therefore be called before `Tendencies::init()`.

All time steppers apply the damping through `VertMix::VertMixImplicit`. No fields, halo exchanges or I/O are added.

## Testing

`VertMixTest` calls `VertMix::applyVelVertMixImplicit` with a test `OceanState` and checks:

- With no vertical viscosity or bottom drag, every active layer decays by exactly $1/(1 + \Delta t \, c_R)$, both for a typical coefficient and for $c_R \Delta t = 10^3$. Velocity must keep its sign and not grow. Inactive layers must be unchanged.
- With viscosity and bottom drag, disabled damping and enabled damping with $c_R = 0$ give bit-for-bit identical results.
- With nonuniform viscosity, implicit bottom drag and damping all enabled, the result matches a host Thomas-algorithm solve of the system above.
