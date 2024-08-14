(omega-dev-aux-state)=

# Auxiliary State

The `AuxiliaryState` class groups together all of Omega [auxiliary
variables](#omega-dev-aux-vars). It is responsible for managing their lifetime,
registering them with IOStreams, and providing functions that compute them,
to be used in the tendency evaluation. It is possible
to have multiple `AuxiliaryState` instances in Omega. Every instance has a name
and is tracked in a static C++ map called `AllAuxStates`.

## Initialization

An instance of the `AuxiliaryState` class requires a [`HorzMesh`](#omega-dev-horz-mesh), so
the mesh class and all of its dependencies need to be initialized before the `AuxiliaryState` class
can be. The static method:
```c++
OMEGA::AuxiliaryState::init();
```
initializes the default `AuxiliaryState`. A pointer to it can be retrieved at any time using:
```c++
OMEGA::AuxiliaryState* DefAuxState = OMEGA::AuxiliaryState::getDefault();
```

## Creation of non-default auxiliary states

A non-default auxiliary state can be created from a string `Name`, horizontal mesh `Mesh`, and number of
vertical levels `NVertLevels`:
```c++
OMEGA::AuxiliaryState*  NewAuxState = OMEGA::AuxiliaryState::create(Name, Mesh, NVertLevels);
```
For conveniece, this returns a pointer to the newly created state. Given its name, a pointer to a named auxiliary state
can be obtained at any time by calling the static `get` method:
```c++
OMEGA::AuxiliaryState* NewAuxState = OMEGA::AuxiliaryState::get(Name);
```

## Computation of auxiliary variables
To compute all auxiliary variables stored in an auxiliary state `AuxState`,
given ocean state `State` and time level `TimeLevel` do:
```c++
AuxState.computeAll(State, TimeLevel);
```

## Removal of auxiliary states
To erase a specific named auxiliary state use `erase`
```c++
OMEGA::AuxiliaryState::erase(Name);
```
To clear all instances do:
```c++
OMEGA::AuxiliaryState::clear();
```
