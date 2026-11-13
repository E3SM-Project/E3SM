(omega-user-split-time-stepping)=

# Split Time Stepping

Ocean motion spans two very different speeds. The barotropic (depth-averaged)
mode carries fast external gravity waves, while the baroclinic
(depth-varying) mode evolves more slowly. An unsplit explicit time stepper
must use a time step short enough for the fast barotropic waves, even though
most of the ocean changes much more slowly. A split-explicit time stepper
separates the two modes: the baroclinic velocity, thickness and tracers are
advanced with the long model time step, while the barotropic mode is subcycled
with a much shorter barotropic time step.

Omega provides two time steppers built on the same RK2 scheme:

| Config option name | Scheme |
| ------------------- | ------- |
| SplitExplicitRK2 | mode-split RK2, with the barotropic mode subcycled by a forward-backward predictor-corrector scheme |
| UnsplitRK2 | the same RK2 scheme without the mode split |

`UnsplitRK2` advances the barotropic mode with the full time step, so it needs
a time step short enough for the barotropic mode. It is useful as a reference
for assessing the effect of mode splitting in `SplitExplicitRK2`.

`UnsplitRK2` and `RungeKutta2` both use an unsplit RK2 approach.
`UnsplitRK2` uses the split-explicit framework with mode splitting disabled
and computes the Coriolis term separately, while `RungeKutta2` includes it in
the potential-vorticity tendency. Both require a time step short enough for
external gravity waves.

## Configuration

The time stepper is selected in the `TimeIntegration` section of the
configuration file, together with a `ModeSplitShare` subgroup holding the
options for `SplitExplicitRK2` and `UnsplitRK2`:

```yaml
  TimeIntegration:
    TimeStepper: SplitExplicitRK2
    TimeStep: 0000_00:10:00
    ModeSplitShare:
      BtrTimeStepper: Predictor-Corrector
      BtrTimeStep: 0000_00:00:20
      NTimeStepIteration: 2
      NBclCoriolisIteration: 2
      ReinitSplitVelocity: false
```

The remaining `TimeIntegration` options are described in the [Time
stepping](#omega-user-time-stepping) section, and the `ModeSplitShare` values
shown are those of the default configuration. The `ModeSplitShare` subgroup is
required whenever one of the two steppers is selected. Within it, only
`BtrTimeStep` is required, and only for `SplitExplicitRK2`; the other options
fall back to their defaults:

| Configuration name | Description | Default |
| ------------------ | ----------- | ------- |
| BtrTimeStepper | scheme for the barotropic subcycle; `Predictor-Corrector` is the only option | Predictor-Corrector |
| BtrTimeStep | barotropic time step, in the same format as `TimeStep`; must be positive and shorter than it | required |
| NTimeStepIteration | positive integer number of predictor-corrector iterations per time step | 2 |
| NBclCoriolisIteration | positive integer number of Coriolis iterations in the baroclinic velocity update; forced to 1 for `UnsplitRK2` | 2 |
| ReinitSplitVelocity | recompute the velocity split at every time step; has no effect for `UnsplitRK2` | false |

`BtrTimeStepper` and `BtrTimeStep` are ignored by `UnsplitRK2`, which also
always uses a single Coriolis iteration. If `NBclCoriolisIteration` is
supplied, it must still be a positive integer before it is overridden. Keep
`NTimeStepIteration: 2` for the RK2 predictor-corrector scheme. A single
iteration retains only the predictor and reduces temporal accuracy to first
order; more iterations do not make this a higher-order Runge Kutta method. The
name `SplitExplicitRK2` describes its RK2 structure, rather than a guarantee of
second-order convergence of the complete split scheme: the current time-stepper
test checks first-order convergence for its split configuration.


For `UnsplitRK2`, a configuration block with all applicable default values
explicitly specified is as follows:

```yaml
  TimeIntegration:
    TimeStepper: UnsplitRK2
    TimeStep: 0000_00:00:20
    ModeSplitShare:
      NTimeStepIteration: 2
```

The time steps in these examples are illustrative; stability depends on the
mesh and the simulated state.

### Barotropic time step

The subcycle count parameter `NBtrSubcycles` is `ceil(TimeStep/BtrTimeStep)`,
and the barotropic time step actually used is `TimeStep` divided by that
number, so it is never longer than `BtrTimeStep`. The resolved values are
written to the log at initialization. For the `SplitExplicitRK2` example above,
the `omega.log` shows:

```
SplitExplicitRK2: TimeStep=600 s, BtrTimeStep=20 s, NBtrSubcycles=30, BtrDt=20 s, NTimeStepIteration=2, NBclCoriolisIteration=2
```

`BtrTimeStep` must be short enough (approximately 20~30x shorter than `TimeStep`)
to resolve the external gravity waves on the
mesh, while `TimeStep` must satisfy the stability limits of the remaining
explicit processes. In each time-step iteration, the complete barotropic
subcycling sequence spans twice the model time step. Each individual barotropic
step has duration `TimeStep/NBtrSubcycles`, so a model time step takes `2 *
NTimeStepIteration * NBtrSubcycles` barotropic steps, 120 in the example above.

### Other requirements

`SplitExplicitRK2` places two requirements on other sections of the
configuration file:

- `HaloWidth` in the `Decomp` section must be at least 3, which is the default.
  The barotropic subcycle uses the extra halo layers to communicate only once
  per subcycle. This is checked when the barotropic update is first called; see [Domain Decomposition](#omega-user-decomp).
- `SSHTendencyEnable` in the `Tendencies` section must be false, which is also
  the default. The split time stepper computes the barotropic pressure gradient
  itself. This is checked at initialization; see [Tendency Terms](#omega-user-tend-terms).

## State variables and restart

The split time steppers carry three state variables in addition to
`NormalVelocity` and `PseudoThickness` (for both `SplitExplicitRK2` and
`UnsplitRK2`):

| Field name | Description | Units | Dimensions |
| ---------- | ----------- | ----- | ---------- |
| NormalBaroclinicVelocity | baroclinic velocity component normal to edge | m s-1 | NEdges, NVertLayers |
| NormalBarotropicVelocity | barotropic velocity component normal to edge | m s-1 | NEdges |
| BarotropicPressureAnomaly | barotropic pressure anomaly | Pa | NCells |

These fields exist only when `SplitExplicitRK2` or `UnsplitRK2` is selected.
They are part of the `Restart` field group but not the `State` group, so they
are written to restart files but are not needed in the initial-condition file.
For `SplitExplicitRK2`, on a cold start the velocity split is computed from the
initial `NormalVelocity` and `PseudoThickness`, and the barotropic pressure
anomaly from the initial pressure. On a restart, `SplitExplicitRK2` reads all
three from the restart file instead. `UnsplitRK2` initializes its baroclinic
velocity from the full `NormalVelocity` and sets its barotropic velocity to
zero on both cold starts and restarts; its barotropic pressure anomaly is
unused.
