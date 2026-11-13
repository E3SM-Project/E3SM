(omega-user-time-stepping)=

# Time stepping
Time stepper refers to the means by which a simulation is advanced in time.
The configuration of a simulation is set by the ``TimeIntegration`` section
of the Omega configuration file:
```yaml
  TimeIntegration:
    CalendarType: No Leap
    TimeStepper: Forward-Backward
    TimeStep: 0000_00:10:00
    StartType: StartUp
    StartTime: 0001-01-01_00:00:00
    StopType: AtTime
    StopCriterion: 0001-01-01_02:00:00
    ModeSplitShare:
      BtrTimeStepper: Predictor-Corrector
      BtrTimeStep: 0000_00:00:20
      NTimeStepIteration: 2
      NBclCoriolisIteration: 2
      ReinitSplitVelocity: false
```
This configuration refers to the default time stepping used for the model
dynamics (momentum and continuity equations). Additional time steppers can
be used for other portions of the model (eg the barotropic mode or tracer
transport) that may have different time steps and use different algorithms.

The ``CalendarType`` choice is described in the
[Time Management](#omega-user-time-manager) section.

The ``TimeStepper`` option refers to the numerical scheme used to advance the
model in time. Omega implements a number of time-stepping schemes. The user
can select the scheme they want in the configuration file.
The following time steppers are currently available:
| Config option name | Scheme |
| ------------------- | ------- |
| Forward-Backward | forward-backward |
| RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
| SplitExplicitRK2 | mode-split RK2 with a subcycled barotropic mode |
| UnsplitRK2 | the same RK2 scheme with mode splitting disabled |

`SplitExplicitRK2` and `UnsplitRK2` use the `ModeSplitShare` configuration
subgroup shown above. See [Split time stepping](#omega-user-split-time-stepping)
for descriptions of these schemes and their configuration options.

The ``TimeStep`` refers to the main model time step used to advance the solution
forward. The time step is specified as a formatted string and can be provided
in any of the following forms:

- ``DDDD_HH:MM:SS(.sss...)``
- ``HH:MM:SS(.sss...)``
- ``MM:SS(.sss...)``
- ``SS(.sss...)``

Days, hours and minutes are optional but must be in order if included.
Fractional seconds are optional.

The ``StartOption`` can be one of five choices. The ``StartUp`` option is for
starting a solution from scratch from an initial state file. The ``Continue``
option is for continuing a simulation from a restart file. The ``Branch``
option will branch from an existing simulation by reading from the restart
file and is identical to ``Continue`` except that some other simulation
metadata (eg case name or simulation name) can be changed. The ``Hybrid``
option will continue a branch from a restart file but will reset the start
time and should only be used for the first step of a hybrid run, similar to the
``StartUp`` option. It is equivalent to using the ``StartUp`` option and
specifying the restart file as an init file. It is only provided in the event
a future mult-level time stepping scheme is used in which starting from a
restart with multiple time levels is an advantage over a fresh start. Finally,
a ``Coupled`` option specifies that the ``StartType`` is provided by the
coupler in a coupled configuration. The actual ``StartType`` used internally
is created by the E3SM CIME tools from a combination of the ``RUN_TYPE`` and
``CONTINUE`` options.

The ``StartTime`` refers to the starting time for the full simulation (not the
current leg of an ongoing simulation). It is in the
format ``yyyy-mm-day_hh:mm:ss`` for year, month, day, hour, minute, second.
If this is a continuation of an existing simulation, the current time for this
leg of the full simulation will be set by the restart file.

A ``StopType`` determines (with the ``StopCriterion`` below) how the simulation
will be stopped. There are three options. The ``AtTime`` option will stop the
simulation at a specific time and the ``StopCriterion`` holds that specific time
as described below. The ``AfterDuration`` option runs the simulation for a fixed
time interval and the ``StopCriterion`` is used to define that interval. A final
option called ``Coupled`` is primarily for coupled simulations where the
simulation will stop after receiving a signal from the coupler.

For the ``AtTime`` stop type, the ``StopCriterion`` must be a time instant in
the format ``yyyy-mm-dd_hh:mm:ss``. If the ``StopType`` is ``AfterDuration``,
the ``StopCriterion`` is a time interval in the format described above for the
time step (but typically ``dddd_hh:mm:ss``). For the ``Coupled`` option the
``StopCriterion`` is ignored.

```{toctree}
:hidden:

SplitTimeStepping
```
