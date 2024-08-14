(omega-design-tendencies)=
# Tendencies

## 1 Overview

The tendencies class will provide a container for the tendency terms
on a subdomain within OMEGA. It will provide methods for computing the
terms for the layer thickness and normal velocity tendencies.

## 2 Requirements

### 2.1 Requirement: Computation of tendency terms
Methods to compute the tendency terms will be provided. This will include
functions to compute the layer thickness tendencies, normal velocity tendencies,
and all tendencies. This will also involve calling the compute function for the
auxiliary state.

### 2.2 Requirement: Multiple tendencies can be present
There will be methods to create and remove non-default tendencies instances, that
can be associated with different meshes.

## 3 Algorithmic Formulation

No algorithms are required

## 4 Design

### 4.1 Data types and parameters
#### 4.1.1 Parameters
No parameters are required.


#### 4.1.2 Class/structs/data types
The tendencies class is an object representing a group of tendency terms and
arrays for the accumulated tendency calculations. Additionally, it will have
a static map of all tendencies and a pointer to the default tendencies instance for quick
retrieval.

```c++
class Tendencies{
 public:
   Array2DReal LayerThicknessTend;
   Array2DReal NormalVelocityTend;
   ThicknessFluxDivOnCell ThicknessFluxDiv;
   PotentialVortHAdvOnEdge PotientialVortHAdv;
   KEGradOnEdge KEGrad;
   SSHGradOnEdge SSHGrad;
   VelocityDiffusionOnEdge VelocityDiffusion;
   VelocityHyperDiffOnEdge VelocityHyperDiff;
 private:
   static Tendencies *DefaultTendencies;
   static std::map<std::string, std::unique_ptr<Tendencies>> AllTendencies;
};
```

### 4.2 Methods

There will be a constructor and destructor for the class and several public
methods. The constructor will be be private to make sure that every tendencies instance
is properly stored in the map. However, a static method to create new
tendencies instance, enforcing this requirement, will be provided.

#### 4.2.1 Creation
The constructor will be responsible for:
  * initializing tendency term functors
  * allocating tendency arrays

```c++
Tendencies(const std::string &Name, const HorzMesh *Mesh, int NVertLevels, Config *Options);
```

The create method will take the same arguments as the constructor, use it to
create a new tendencies instance, and put it in the static map of all tendencies.
It will return a pointer to the newly created object.
```c++
Tendencies* Tendencies::create(const std::string &Name, const HorzMesh *Mesh, int NVertLevels, Config *Options);
```

#### 4.2.2 Initialization
The init method will create the default tendencies and return an error code:
```c++
int Tendencies::init();
```

#### 4.2.3 Retrieval
There will be methods for getting the default and non-default tendencies instances:
```c++
Tendencies *Tendencies::getDefault();
Tendencies *Tendencies::get(const std::string &Name);
```

#### 4.2.4 Computation
The 'computeAll' method will compute and accumulate all layer thickness and normal velocity tendencies using the ocean
state and auxiliary state from a given time level:
```c++
void Tendencies::computeAllTendencies(const OceanState *State, const AuxilaryState *AuxState, int TimeLevel);
```
The layer thickness tendencies will be computed with a method:
```c++
void Tendencies::computeThicknessTendencies(const OceanState *State, const AuxilaryState *AuxState, int TimeLevel);
```
The normal velocity tendencies will be computed with a method:
```c++
void Tendencies::computeVelocityTendencies(const OceanState *State, const AuxilaryState *AuxState, int TimeLevel);
```


#### 4.2.5 Destruction and removal
No operations are needed in the destructor as the Kokkos arrays are removed
when they are no longer in scope. The erase method
will remove a named tendencies instance, whereas the clear method will remove all of
them. Both will call the destructor in the process.
```c++
void Tendencies::erase(const std::string &Name);
void Tendencies::clear();
```

## 5 Verification and Testing

### 5.1 Test multiple tendencies instances

There will be a test checking that multiple tendencies objects can be
created and that removing them works properly.

### 5.2 Test computation

There will be a test to ensure all of the tendencies are correctly computed by
the compute method.
