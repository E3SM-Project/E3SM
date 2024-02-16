(omega-design-state)=
# State

## 1 Overview

The state class will provide a container for the prognostic (normal velocity and layer thickness) variables on a subdomain within OMEGA.
It will have functionality for reading initial conditions, managing host/device data transfer, and performing time level updates.
This will also be a data type that can be used for tendency terms and provis variables within the timestepping scheme.

## 2 Requirements

### 2.1 Requirement: State should read in initial conditions for prognostic variables
The state class will have a method to read the subdomain initial conditions in parallel using information from Decomp.

### 2.2 Requirement: Manage host/device data transfer
Initial condition data read on host will be transferred to device and the prognostic variables computed on the device will be transferred for output.

### 2.3 Requirement: Update timestep
A method will be responsible for advancing the prognostic variables prior to the beginning of a new timestep. This will include a halo exchange a deep copy from the current state variable to an "old" variable.

### 2.4 Requirement: Can be used as a common data type for tendency and provis variables
The timestepping scheme will require storage for the tendency terms and the intermediate stage. The state class can also be used to fulfill this need.

## 3 Algorithmic Formulation

No algorithms are required

## 4 Design

### 4.1 Data types and parameters
#### 4.1.1 Parameters
No parameters are required.


#### 4.1.2 Class/structs/data types
The State class will contain public Kokkos host and device views for the normal velocity and layer thickness variables. There will also be an additional device view for the "old" prognostic variables on the device. 

```c++
class State{

public:

  Array2DR8 LayerThickness;
  ArrayHost2DR8 LayerThicknessH;
  Array2DR8 LayerThicknessOld;

  Array2DR8 NormalVelocity;
  ArrayHost2DR8 NormalVelocityH;
  Array2DR8 NormalVelocityOld;

}
```

### 4.2 Methods

There will be a constructor and destructor for the class and several public methods. Some of these methods are only applicable to the prognostic state and not to the tendency and provis contexts (e.g. read).

#### 4.2.1 Constructor
The constructor will also be responsible for:
  * allocating prognostic variables
  * getting global IDs for potential read
  * registering metadata with the I/O infrastructure

```c++
State(Decomp *DefDecomp);
```

#### 4.2.2 Destructor
A destructor will be available to release memory.

#### 4.2.3 Read
The read method will read the state variable information from a initial condition or restart file in parallel.

#### 4.2.4 Copy to device
A method will provide functionality to deep copy host data to the device for the prognostic variables.

#### 4.2.5 Copy to host
There will also be a method to perform a deep copy of the device state to the host.

#### 4.2.6 Metadata registration
The metadata associated with each state variable will be registered within the I/O infrastructure in a private method called by the constructor.

#### 4.2.7 Time level update
The time level update method will perform a halo exchange for the prognostic variables and advance the state variables for the next timestep.

## 5 Verification and Testing

### 5.1 Test read

The sample domain used for the Decomp test will be used to test obtaining the correct local values from the read routine.

### 5.2 Test time level update

There will be a test to ensure the time level update functions correctly.
