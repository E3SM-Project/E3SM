(omega-dev-time-stepping)=

# Time stepping

Omega supports choice of multiple time stepping schemes. The time stepping
module provides implementations of different time stepping schemes and general
infrastructure for managing time steppers. To achieve modularity, each time
stepper is implemented as a separate class. To ensure common interface and
enable sharing of common code and data, every time stepper inherits from the
`TimeStepper` base class.

## Time stepper type
An enumeration listing all implemented schemes is provided. It needs to be extended every time a new time stepper
is added. It is used to identify a time stepper at run time.
```c++
enum class TimeStepperType { ForwardBackward, RungeKutta4, RungeKutta2 };
```

## TimeStepper base class
The `TimeStepper` class is a base class for all Omega time steppers. It stores
as members common data needed by every time stepper. These include a name and
time stepping method as well as some time tracking variables, including
StartTime, StopTime, TimeStep, a Clock and an EndAlarm. Its `doStep` method
defines the interface that every time stepper needs to implement. In addition to this
method, it provides a number of other methods that can divided into two groups.
The first group is for general time stepper management. The second group is
utility functions that provide common functionality that can be re-used in
implementation of different time steppers.

### Do step method
The `doStep` method is the main method of every time stepper class. It exists
in the base class as a pure virtual method, meaning it defines an interface but
has no implementation. Its signature is
```c++
    void TimeStepper::doStep(OceanState *State, TimeInstant &SimTime) const;
```
It is meant to advance `State` by one time step, from `Time` to `Time +
TimeStep`. Classes implementing concrete time stepping schemes need to provide
implementations of this method. As part of the time step, this routine advances
the Clock associated with the time stepper and returns the current simulation
time. For the default time stepper, this is advances the primary model clock.

### Time stepper management methods

#### Initialization

Two static methods:
```c++
TimeStepper::init1();
TimeStepper::init2();
```
initialize the default time stepper in two phases. The first phase must be
called early in the initialization process. It creates and time stepper
instance, fills it with various time quantities (StartTime, EndTime, TimeStep)
and creates a Clock (the Model Clock for the simulation) and an EndAlarm to
ring when this portion of the simulation should stop. The second phase is
called after the default tendency, auxiliary state, horizontal mesh, and halo
objects have been initialized. It attaches pointers to these objects for later
use during the time stop phase. A pointer to the default time stepper can be
retrieved at any time using:
```c++
TimeStepper* DefTimeStepper = TimeStepper::getDefault();
```

#### Creation of non-default time steppers

A non-default time stepper can be created from a string `Name`, time stepper
type `Type`, `TimeStep`, `StartTime`, `EndTime`, tendencies `Tend`,
auxiliary state `AuxState`, horizontal mesh `Mesh`, and halo layer `MyHalo`
```c++
TimeStepper*  NewTimeStepper = TimeStepper::create(Name, Type, StartTime,
         EndTime, TimeStep, Tend, AuxState, Mesh, MyHalo);
```
For convenience, this returns a pointer to the newly created time stepper.
Given its name, a pointer to a named time stepper can be obtained at any time
by calling the static `get` method:
```c++
TimeStepper* NewTimeStepper = TimeStepper::get(Name);
```

#### Changing the time step
The time step is set on creation, but if the time step must change, the
function
```c++
Stepper->changeTimeStep(TimeStep);
```
can be used where `TimeStep` is an instance of `TimeInterval` class.

#### Getters
Given a pointer to a `TimeStepper` you can obtain various class members.
Pointers to the associated Clock and EndAlarm can also be retrieved.
```c++
TimeStepperType Type = Stepper->getType();
int NTimeLevels = Stepper->getNTimeLevels();
std::string Name = Stepper->getName();
TimeInterval TimeStep = Stepper->getTimeStep();
TimeInstant StartTime = Stepper->getStartTime();
TimeInstant StopTime = Stepper->getStopTime();
Clock *ModelClock = Stepper->getClock();
Alarm *EndAlarm = Stepper->getEndAlarm();
```

#### Removal of time steppers
To erase a specific named time stepper use `erase`
```c++
TimeStepper::erase(Name);
```
To clear all instances do:
```c++
TimeStepper::clear();
```

### Utility methods

In addition to the above management methods, the base class provides common
functions that are useful in implementing different time steppers. Given two
pointers to `OceanState` objects `State1` and `State2`, two integers
representing two time levels `TimeLevel1` and `TimeLevel2`, and coefficient
`Coeff` of type `TimeInterval` calling
```c++
updateStateByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
conceptually performs the following update
`State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend`
element-wise for every state variable. To perform this update only for
layer thickness do
```c++
updateThicknessByTend(State1, TimeLevel1, *State2, TimeLevel2, Coeff);
```
and only for normal velocity do
```c++
updateVelocityByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
Similarly for tracers, the method for updating tracers is
```c++
updateTracersByTend(NextTracers, CurTracers, State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
Additionally, for the RK4 method, there are methods for initializing, accumulating, and
finalizing the tracer integration for a given time step
```c++
weightTracers(NextTracers, CurTracers, CurState, TimeLevel1);
accumulateTracersUpdate(AccumTracer, Coeff);
finalizeTracersUpdate(NextTracers, State, TimeLevel);
```

## Implemented time steppers
The following time steppers are currently implemented
| Class name | Enum value | Scheme |
| ---------- | -------------------------- | ------ |
| ForwardBackwardStepper | TimeStepperType::ForwardBackward | forward-backward |
| RungeKutta2Stepper | TimeStepperType::RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4Stepper | TimeStepperType::RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
