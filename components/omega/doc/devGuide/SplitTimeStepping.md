(omega-dev-split-time-stepping)=

# Split Time Stepping

The `SplitExplicitRK2Stepper` class implements the mode-split RK2 time stepping
scheme. It is used by two time stepper types that share all of their code and
select different splitting and Coriolis-iteration settings:

| Enum value | Split factor | Barotropic subcycle |
| ---------- | ------------ | ------------------- |
| TimeStepperType::SplitExplicitRK2 | 1 | yes, the barotropic time stepping is subcycled |
| TimeStepperType::UnsplitRK2 | 0 | no, the barotropic time stepping is skipped |

The type is used only by the constructor, to decide how the options are read:

```c++
SEConfig(SplitExplicitInit::readConfigOptions(
    InTimeStep, InType == TimeStepperType::UnsplitRK2))
```

The algorithm then uses the resolved `SEConfig` options. For `UnsplitRK2`,
`SplitFactor` is zero, so the barotropic velocity and the barotropic forcing
are zero, the full velocity is carried in the baroclinic velocity array.
Everything below applies to both types unless it is marked as
split-only.

Both `UnsplitRK2` and `RungeKutta2` steppers use an unsplit RK2 approach.
`UnsplitRK2` uses `SplitExplicitRK2Stepper` with `SplitFactor = 0` and computes
the linear Coriolis term separately with `CoriolisTendMode::Separate` and one
Coriolis iteration. `RungeKutta2` uses `RungeKutta2Stepper`, where the Coriolis
term is included in the potential-vorticity tendency (`CoriolisTendMode::PVFlux`).

## Source layout

| File | Contents |
| ---- | -------- |
| `SplitExplicitTypes.h` | `SplitExplicitConfig` and `SplitExplicitScratch` |
| `SplitExplicitInit.{h,cpp}` | config reading, scratch allocation, velocity split |
| `SplitExplicitRK2Stepper.{h,cpp}` | Stages 1 and 3, and the time-step iteration |
| `SplitExplicitBarotropicPCStepper.{h,cpp}` | Stage 2, the barotropic subcycle |

## Configuration

The options are read from the `ModeSplitShare` subgroup of `TimeIntegration`,
which is required whenever one of the two steppers is selected, and are stored
in a `SplitExplicitConfig` member of the stepper.

| Option | Default | Notes |
| ------ | ------- | ----- |
| BtrTimeStepper | Predictor-Corrector | only choice; ignored by `UnsplitRK2` |
| BtrTimeStep | none | required for `SplitExplicitRK2`, must satisfy `0 < BtrTimeStep < TimeStep` |
| NTimeStepIteration | 2 | positive integer iterations over Stages 1-3; 1 retains only the first-order predictor |
| NBclCoriolisIteration | 2 | must be positive if supplied; then forced to 1 for `UnsplitRK2` |
| ReinitSplitVelocity | false | recompute the velocity split every step; no effect for `UnsplitRK2` |

`NBtrSubcycles` is not read from the config; it is derived as
`ceil(TimeStep/BtrTimeStep)`, and the resolved options are reported with
`LOG_INFO` at initialization.

Two further requirements are enforced elsewhere. `Tendencies::setModeSplit`
aborts if `SSHTendencyEnable` is true together with a nonzero split factor, and
the barotropic subcycle aborts unless `Decomp::HaloWidth` is at least 3.

## Initialization

`finalizeInit` is the split-explicit part of `TimeStepper::init2`. After the
usual pointer checks it configures the shared velocity tendency, allocates the
scratch arrays and initializes the barotropic substepper:

```c++
Tend->setModeSplit(CoriolisTendMode::Separate, SEConfig.SplitFactor);
SplitExplicitInit::allocateScratch(SEScratch, Mesh, VCoord, Name);
BarotropicPCStepper.init(AuxState, &SEScratch, &SEConfig, Mesh, MeshHalo, VCoord);
```

The first call tells the velocity tendency to leave the linear Coriolis
acceleration out of the vorticity flux, so that Stage 1 can iterate it, and
which barotropic weight to use for the surface pressure gradient.

`isSplit()` returns true for both types implemented by this class, including
`UnsplitRK2`; it indicates that the additional state arrays are required.
`OceanState::init` queries it on the default time stepper and only allocates
`NormalBaroclinicVelocity`, `NormalBarotropicVelocity` and
`BarotropicPressureAnomaly` when it is true.

Once the initial or restart file has been read, the driver calls
`initializeStateFromInput`. For a cold start it initializes the barotropic
pressure anomaly and splits `NormalVelocity` into its two parts:

```c++
SplitExplicitInit::initializeBarotropicPressure(SEScratch, State, Mesh, VCoord, CurLevel);
SplitExplicitInit::computeVelocitySplit(State, Mesh, VCoord, CurLevel);
```

while a zero split factor instead copies the full velocity into the baroclinic
array at both time levels:

```c++
SplitExplicitInit::computeUnsplitVelocitySplit(State, Mesh, VCoord, CurLevel, NextLevel);
```

On a split restart, the fields read from the restart file are retained instead
of recomputing the initial split and pressure. Either way it finishes by
copying the current level into the next level with `initializeNextState`, then
exchanging the halo and mirroring the state to the host, since the split is
computed over all edges.

## One time step

`doStep` copies the current level into the next level and then repeats three
stages `NTimeStepIteration` times:

```c++
doBaroclinicVelocityUpdate(State, NextTracerArray, CurLevel, NextLevel, VelStageTime, TimeStep);
if (SEConfig.SplitFactor != 0._Real) {
   BarotropicPCStepper.doBarotropicVelocityUpdate(State, CurLevel, NextLevel, TimeStep);
}
computeTransportVelocity(State, NextLevel);
doThicknessTracerUpdate(State, CurTracerArray, NextTracerArray, CurLevel, NextLevel, StageTime, TimeStep, FinalIteration);
```

The first outer iteration is the predictor and evaluates the momentum right-hand side
at time `n`; later iterations see the midpoint state left by their predecessor,
so `VelStageTime` is `StageTime + 0.5 * TimeStep` for them. The last iteration
is flagged with `FinalIteration`. The baroclinic velocity halo is exchanged
after Stage 1, and the whole state and tracer halos after Stage 3 whenever
another iteration follows.

### Stage 1: baroclinic velocity

`doBaroclinicVelocityUpdate` computes the velocity tendencies with
`Tend->computeVelocityTendencies` and then calls
`doBaroclinicCoriolisIteration`, which repeats the following
`NBclCoriolisIteration` times:

```c++
Tend->computeCoriolisAccelerationOnEdge(IterVelocityTend, BaseVelocityTend, NormalBclVelEdge, FEdge);
updateBaroclinicVelocityWithBarotropicForcing(State, CurLevel, NextLevel, StageTimeStep);
```

The base tendency in `Tend->NormalVelocityTend` is left untouched so that every
iteration can re-read it; the Coriolis acceleration of the current baroclinic
velocity is added to it in the `IterVelocityTend` scratch array instead.
Intermediate iterations exchange the baroclinic velocity halo; after the loop
the barotropic forcing halo is exchanged for the split case.

`updateBaroclinicVelocityWithBarotropicForcing` makes one pass over each edge
column. It first forms the barotropic forcing, the thickness-weighted column
mean of the provisional full-step baroclinic velocity divided by the time step,
and stores it in `BarotropicForcing` for Stage 2. For each wet edge column,
with $h_k$ from `MeanPseudoThickEdge`, $u^n_{bcl,k}$ from `NormalBclVelCur`,
$R_k$ from `IterVelocityTend`, and split factor $s$,

$$
F = \frac{s}{\Delta t}
\frac{\sum_k h_k (u^n_{bcl,k} + \Delta t R_k)}{\sum_k h_k}.
$$

This includes any column mean remaining in the old baroclinic velocity under
the current thickness weights. In the code, `FluxSum` is the numerator:

```c++
Forcing = LocSplitFactor * (FluxSum / ThicknessSum) * InvDtSeconds;
```

and then advances the baroclinic velocity to the stage midpoint with that
forcing removed:

```c++
NormalBclVelNext(IEdge, K) = NormalBclVelCur(IEdge, K) + HalfDtSeconds * (IterVelocityTend(IEdge, K) - Forcing);
```

The forcing block is skipped when the split factor is zero, so the `UnsplitRK2`
update is the plain half-step advance of the full velocity.

### Stage 2: barotropic velocity (split only)

`SplitExplicitBarotropicPCStepper::doBarotropicVelocityUpdate` integrates the
barotropic mode with a forward-backward predictor-corrector scheme, subcycled
over `2*NBtrSubcycles` steps of `TimeStep/NBtrSubcycles`, with the Stage 1
barotropic forcing held fixed. Each subcycle applies a velocity predictor, a
pressure predictor, a velocity corrector and a pressure corrector, blended with
the feedback weights `Gamma1 = 0.5333`, `Gamma2 = 0.5333` and `Gamma3 = 1.0`,
and then rotates the `Cur`/`Cor` scratch buffers. The corrector velocity and
flux are accumulated, so `NormalBarotropicVelocity` at the next level ends up
holding the subcycle mean rather than an instantaneous value, and
`BarotropicFlux` holds the mean barotropic transport used by Stage 3. The
velocity mean includes the initial velocity and all `2*NBtrSubcycles` corrector
velocities; the flux mean averages the `2*NBtrSubcycles` corrector fluxes. Each
Stage 2 call restarts from the current model time level, integrating over a
window of length `2*TimeStep`.

To avoid halo exchanges between predictor and corrector kernels, the two `Cur`
buffers are exchanged, once at the top of each subcycle. Each kernel then
computes over one halo layer fewer than its inputs cover, which is why the loop
ranges use `NEdgesHaloH`/`NCellsHaloH` and why a halo width of at least 3 is
required.

### Transport velocity

`computeTransportVelocity` sums the barotropic and baroclinic velocities into
`NormalVelocity` at the next level and stores, in the `TransportVelocityAdd`
scratch array, the per-column correction that makes the column transport match
the time-averaged barotropic flux from Stage 2:

```c++
VelCorrection = (BtrFlux(IEdge) / RhoGravity - TransportSum) / ThickSum;
```

The correction is left at zero when the split factor is zero.

### Stage 3: pseudo-thickness and tracers

`doThicknessTracerUpdate` computes the thickness and tracer auxiliary variables
with `TransportVelocityAdd` applied, computes their tendencies at the new time
level, and then updates depending on the iteration. The final iteration retains
the conservative update over the full time step:

```c++
updateThicknessByTend(State, NextLevel, State, CurLevel, StageTimeStep);
updateTracersByTend(NextTracerArray, CurTracerArray, State, NextLevel, State, CurLevel, StageTimeStep);
```

while earlier iterations construct the RK2 midpoint state, in which the tracer
concentration is the average of its old and provisional full-step values rather
than a conservative update over half a time step:

```c++
updateThicknessByTend(State, NextLevel, State, CurLevel, 0.5 * StageTimeStep);
updateTracersToMidpoint(NextTracerArray, CurTracerArray, State, CurLevel, StageTimeStep);
```

`finalizeTimeStepIterationState` closes the iteration. On the final iteration
`reconstructNormalVelocity` extrapolates the baroclinic velocity from `n+1/2`
to `n+1` and adds the barotropic velocity back in:

```c++
NormalVelNext(IEdge, K) = 2._Real * NormalBclVelNext(IEdge, K) - NormalBclVelCur(IEdge, K) + NormalBtrVelNext(IEdge);
```

For the split stepper, on every iteration it also recomputes the total
pseudo-thickness and resets the barotropic pressure anomaly from it, so that
pressure and thickness stay consistent:

```c++
BtrPressAnomalyNext(ICell) = RhoGravity * (TotalPseudoThickness(ICell) - BottomGeomDepth(ICell));
```

### After the iterations

`doStep` rotates the state and tracer time levels, refreshes the kinetic
auxiliary variables from the completed `n+1` velocity, applies implicit
vertical mixing with a further halo exchange when it is enabled, validates the
state with `validateOceanState`, and advances the clock.

## Scratch arrays

The scratch arrays are allocated once by `SplitExplicitInit::allocateScratch`
and are held in a `mutable SplitExplicitScratch` member, since `doStep` is
const.

| Array | Location | Purpose |
| ----- | -------- | ------- |
| NormalBarotropicVelocitySubcycle{Cur,Pre,Cor} | edge | subcycle state, predictor and corrector |
| BarotropicPressureAnomalySubcycle{Cur,Pre,Cor} | cell | subcycle state, predictor and corrector |
| BarotropicForcing | edge | column-mean forcing passed from Stage 1 to Stage 2 |
| BarotropicFlux | edge | subcycle-mean barotropic transport |
| BaroclinicPseudoThicknessEdge | edge | column sum of the flux pseudo-thickness |
| IterVelocityTend | edge, layer | base plus Coriolis tendency of the current iteration |
| TransportVelocityAdd | edge, layer | transport velocity correction for Stage 3 |

## Testing

Both types are exercised by `TimeStepperTest`, which measures the convergence
rate against an exact exponentially decaying velocity solution with fixed
thickness. `UnsplitRK2` is checked for second order; `SplitExplicitRK2` is
currently checked for first order in that test configuration. This test
expectation does not establish second-order accuracy of the full split scheme
or isolate the source of its error. The test also exercises both types without
a stop time/end alarm and checks that the step count increments across calls to
`doStep`.
