(omega-design-aux-state)=
# AuxiliaryState

## 1 Overview

The auxiliary state class will provide a container for the auxiliary variables
on a subdomain within OMEGA. It will provide methods for computing the
variables and it will register fields with IOStreams.

## 2 Requirements

### 2.1 Requirement: Auxiliary state should register auxiliary variables with IOStreams
The auxiliary state class will register the auxiliary variables with IOStreams
in its constructor to facilitate their output.

### 2.2 Requirement: Computation of auxiliary variables
Methods to compute the auxiliary variables will be provided. At minimum, a
method to compute all of the variables is needed, but convenience functions
that compute just a subset of them might also be useful.

### 2.3 Requirement: Multiple auxiliary states can be present
There will be methods to create and remove non-default auxiliary states, that
can be associated with different meshes.

## 3 Algorithmic Formulation

No algorithms are required

## 4 Design

### 4.1 Data types and parameters
#### 4.1.1 Parameters
No parameters are required.


#### 4.1.2 Class/structs/data types
The auxiliary state class will contain objects representing auxiliary variable
groups, strings for names of the auxiliary state and its metadata group, and a
pointer to the associated mesh. Additionally, it will have a static map of all
auxiliary states and a pointer to the default auxiliary state for quick
retrieval.

```c++
class AuxiliaryState {
 public:
   std::string Name;
   std::string GroupName;
   KineticAuxVars KineticAux;
   LayerThicknessAuxVars LayerThicknessAux;
   VorticityAuxVars VorticityAux;
   VelocityDel2AuxVars VelocityDel2Aux;
   // more variables in the future
 private:
   const HorzMesh *Mesh;
   static AuxiliaryState *DefaultAuxState;
   static std::map<std::string, std::unique_ptr<AuxiliaryState>> AllAuxStates;
};
```

### 4.2 Methods

There will be a constructor and destructor for the class and several public
methods. The constructor will be be private to make sure that every auxiliary
state is properly stored in the map. However, a static method to create new
auxiliary states, enforcing this requirement, will be provided.

#### 4.2.1 Creation
The constructor will be responsible for:
  * constructing auxiliary variables
  * registering fields and metadata with the I/O infrastructure

```c++
AuxiliaryState(const std::string &Name, const HorzMesh *Mesh, int NVertLevels);
```

The create method will take the same arguments as the constructor, use it to
create a new auxiliary state, and put it in the static map of all auxiliary
states. It will return a pointer to the newly created state.
```c++
AuxiliaryState* AuxiliaryState::create(const std::string &Name, const HorzMesh *Mesh, int NVertLevels);
```

#### 4.2.2 Initialization
The init method will create the default auxiliary state and return an error code:
```c++
int AuxiliaryState::init();
```

#### 4.2.3 Retrieval
There will be methods for getting the default and non-default auxiliary states:
```c++
AuxiliaryState *AuxiliaryState::getDefault();
AuxiliaryState *AuxiliaryState::get(const std::string &Name);
```

#### 4.2.4 Computation
The 'computeAll' method will compute all auxiliary variables using the ocean
state from a given time level:
```c++
void AuxiliaryState::computeAll(const OceanState *State, int TimeLevel);
```

#### 4.2.5 Destruction and removal
The destructor will unregister fields and remove metadata. The erase method
will remove a named auxiliary state, whereas the clear method will remove all of
them. Both will call the destructor in the process.
```c++
void AuxiliaryState::erase(const std::string &Name);
void AuxiliaryState::clear();
```

## 5 Verification and Testing

### 5.1 Test multiple auxiliary states

There will be a test checking that multiple auxiliary state object can be
created and that removing them works properly.

### 5.2 Test computation

There will be a test to ensure all of the variables are correctly computed by
the compute method.
