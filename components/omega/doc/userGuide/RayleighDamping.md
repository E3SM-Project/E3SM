(omega-user-rayleigh-damping)=

# Rayleigh Damping

Rayleigh damping adds a linear drag on the normal velocity to the momentum equation:

$$
\left. \frac{\partial u}{\partial t} \right|_{\mathrm{Rayleigh}} = -c_R u ,
$$

where $c_R$ (s$^{-1}$) is the damping coefficient and $1/c_R$ is the e-folding time scale of the damping. The damping acts on every active layer of every edge, but not on pseudo-thickness or tracers.

Rayleigh damping is a spin-up tool, not part of the physics of a production simulation. It is used for *dynamic adjustment*: a short run, starting from an initial condition that is not in dynamical balance (e.g. one built from observational climatologies), in which strong damping removes the large barotropic and inertia-gravity waves generated as the model adjusts. A dynamic adjustment is typically run as a sequence of jobs (restarting from the previous one) in which the damping coefficient is reduced and the time step increased from one job to the next until damping can be turned off and the time step has reached its production value. Coefficients of order $10^{-4}$ s$^{-1}$ (a time scale of about 2.8 hours) are typical for the first, most aggressive job.

Rayleigh damping is configured with the `RayleighDamping` group in the yaml configuration file:

```yaml
RayleighDamping:
  Enable: false        # Enables Rayleigh damping
  DampingCoeff: 1.0e-4 # Damping coefficient (s^-1)
```

Rayleigh damping is disabled by default, in which case `DampingCoeff` has no effect and results are bit-for-bit identical to results without the feature. `DampingCoeff` must be non-negative and is read only at initialization, so changing it requires a new job.

The damping term is treated implicitly in time as part of the implicit vertical mixing of velocity, so it is stable for any value of $c_R \Delta t$. Rayleigh damping therefore requires implicit velocity vertical mixing to be enabled:

```yaml
Tendencies:
  VelVertMixTendencyEnable: true
```

which is the default. Omega aborts during initialization if Rayleigh damping is enabled and `VelVertMixTendencyEnable` is `false`.

See the [Rayleigh damping design document](#omega-design-rayleigh-damping) for the details of the formulation.
