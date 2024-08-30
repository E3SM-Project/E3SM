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
as members common data needed by every times stepper. Its `doStep` method
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
    void TimeStepper::doStep(OceanState *State, TimeInstant Time) const;
```
It is meant to advance `State` by one time step, from `Time` to `Time +
TimeStep`. Classes implementing concrete time stepping schemes need to provide
implementations of this method.

### Time stepper management methods

#### Initialization

The static method:
```c++
TimeStepper::init();
```
initializes the default time stepper. A pointer to it can be retrieved at any time using:
```c++
TimeStepper* DefTimeStepper = TimeStepper::getDefault();
```

#### Creation of non-default time steppers

A non-default time stepper can be created from a string `Name`, time stepper
type `Type`, tendencies `Tend`, auxiliary state `AuxState`, horizontal mesh
`Mesh`, and halo layer `MyHalo`
```c++
TimeStepper*  NewTimeStepper = TimeStepper::create(Name, Type, Tend, AuxState, Mesh, MyHalo);
```
For convenience, this returns a pointer to the newly created time stepper.
Given its name, a pointer to a named time stepper can be obtained at any time
by calling the static `get` method:
```c++
TimeStepper* NewTimeStepper = TimeStepper::get(Name);
```

#### Setting the time step
After creating a time stepper, the time step can be set by
```c++
Stepper->setTimeStep(TimeStep);
```
where `TimeStep` is an instance of `TimeInterval` class.

#### Getters
Given a pointer to a `TimeStepper` you can obtain its type, name, number of time levels,
or time step by calling
```c++
TimeStepperType Type = Stepper->getType();
int NTimeLevels = Stepper->getNTimeLevels();
std::string Name = Stepper->getName();
TimeInterval TimeStep = Stepper->getTimeStep();
```
respectively.


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

In addition to the above management methods, the base class provides common functions that
are useful in implementing different time steppers. Given two pointers to
`OceanState` objects `State1` and `State2`, two integers representing two time
levels `TimeLevel1` and `TimeLevel2`, and coefficient `Coeff` of type
`TimeInterval` calling
```c++
updateStateByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```
conceptually performs the following update `State1(TimeLevel1) = State2(TimeLevel2) + Coeff * Tend`
element-wise for every state variable. To perform this update only for
layer thickness do
```c++
updateThicknessByTend(State1, TimeLevel1, *State2, TimeLevel2, Coeff);
```
and only for normal velocity do
```c++
updateVelocityByTend(State1, TimeLevel1, State2, TimeLevel2, Coeff);
```

## Implemented time steppers
The following time steppers are currently implemented
| Class name | Enum value | Scheme |
| ---------- | -------------------------- | ------ |
| ForwardBackwardStepper | TimeStepperType::ForwardBackward | forward-backward |
| RungeKutta2Stepper | TimeStepperType::RungeKutta2 | second-order two-stage midpoint Runge Kutta method |
| RungeKutta4Stepper | TimeStepperType::RungeKutta4 | classic fourth-order four-stage Runge Kutta method |
