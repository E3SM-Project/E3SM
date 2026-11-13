(omega-design-state)=
# State

## 1 Overview

The state class will provide a container for the prognostic (normal velocity and layer thickness) variables on a subdomain within OMEGA.
It will have functionality for registering fields with IOStreams and performing time level updates.
This will also be a data type that can be used for tendency terms and provis variables within the timestepping scheme.

## 2 Requirements

### 2.1 Requirement: State should register prognostic variables with IOStreams
The state class will have a method to register the progostic variables with IOStreams to provide initial condition input and solution output.

### 2.2 Requirement: Update timestep
A method will be responsible for advancing the prognostic variables prior to the beginning of a new timestep.
The left-most index of the prognostic variables will be a time index, which will facilitate flexibility with timestepping schemes.
The time level update will include a halo exchange and an exchange of time level.

### 2.3 Requirement: Can be used as a common data type for tendency and provis variables
The timestepping scheme will require storage for the tendency terms and the intermediate stage. The state class can also be used to fulfill this need.

## 3 Algorithmic Formulation

No algorithms are required

## 4 Design

### 4.1 Data types and parameters
#### 4.1.1 Parameters
No parameters are required.


#### 4.1.2 Class/structs/data types
The State class will contain public Kokkos host and device views for the normal velocity and layer thickness variables.

```c++
class State{

public:

  Array3DR8 LayerThickness;
  Array3DR8 NormalVelocity;

}
```

### 4.2 Methods

There will be a constructor and destructor for the class and several public methods. Some of these methods are only applicable to the prognostic state and not to the tendency and provis contexts.

#### 4.2.1 Constructor
The constructor will also be responsible for:
  * allocating prognostic variables
  * registering fields and metadata with the I/O infrastructure

```c++
State(Decomp *DefDecomp, int TimeLevels);
```

#### 4.2.2 Destructor
A destructor will be available to release memory.

#### 4.2.3 Time level update
The time level update method will perform a halo exchange for the prognostic variables and advance the state variables for the next timestep.

## 5 Verification and Testing

### 5.1 Test read

The sample domain used for the Decomp test will be used to test obtaining the correct local values from the IO infrastructure

### 5.2 Test time level update

There will be a test to ensure the time level update functions correctly.
