(omega-design-auxvars)=
# Auxiliary variables

## 1 Overview

Auxiliary variables are variables that are needed to compute the tendency terms
and advance the model state. While they are functions of the model state, they
differ from variables needed only for diagnostic purposes. Auxiliary variables
can often be bundled together into groups of variables that are logically
related, or can be efficiently computed together. Each auxiliary variable group
is implemented as a class containing array(s) storing the variable(s) and member
functions that compute them. The member functions compute the variable (or
variable group) over a chunk of vertical levels for a particular mesh location.
There may be more than one compute function, since some groups may contain
variables defined on different mesh elements. This approach allows flexible
groupings of computational work within larger cell/edge/vertex loops. This type
of flexibility will be important for optimizing the code across different device
architectures.

## 2 Requirements

### 2.1 Requirement: Functions that compute auxiliary variables should have simple, intuitive calling arguments
Auxiliary variables take in constant data as constructor arguments, which
simplifies the arguments used to call the compute methods.

### 2.2 Requirement: Auxiliary variables computations should be able to be grouped under different outer loops over mesh locations
The ability to fuse/separate auxiliary variables computations will allow for
optimization across different device architectures. Member functions will
implement a given variable (group) computation on a specific mesh location
(cell, edge, vertex).

### 2.3 Requirement: Vectorization on CPU architectures
Auxiliary variables computations have inner loops over a chunk of vertical
levels. The chunk size will be set to the vector length on CPU machines and 1
for GPUs. This will allow for the possibility of vectorization on CPUs.

### 2.4 Requirement: Configuration options
Some auxiliary variables need to be computed differently based on configuration
options. The class constructor will retrieve possible options from the Config
and store them as enums. The compute functions will do different things based on
the stored options.

### 2.5 Requirement: Output
Auxiliary variables can also be useful as diagnostics. Each class will have a
method to register the variables with IOStreams to provide diagnostic output.

### 2.6 Desired: Small memory footprint
Some auxiliary variables are used only to compute another auxiliary variable,
and are not needed afterwards. It may be desirable to find a memory management
strategy that would allow to minimize the amount of persistent memory needed.

## 3 Algorithmic Formulation

The algorithmic formulation for each auxiliary variable can be found in the
algorithms design document.

## 4 Design

### 4.1 Data types and parameters
#### 4.1.1 Parameters

Each configurable auxiliary variable will define its parameters as an enum. Example:
```c++
enum FluxThickType { Center, Upwind };
```
These parameters will be read from the input configuration file
```yaml
    advection:
       fluxThicknessType: 'Center'
```
and stored inside the class.

#### 4.1.2 Class/structs/data types
The auxiliary variables group is a class that contains arrays for storing the
variables and methods to compute them. Configuration choices will also be
stored. The class also contains private member variables for any constant data
that are defined by the constructor.

```c++
class LayerThicknessAuxVars {
 public:
   Array2DReal FluxLayerThickEdge;
   Array2DReal MeanLayerThickEdge;
   FluxThickType FluxThickChoice;

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge, int KChunk,
                                          const Array2DReal &LayerThickness,
                                          const Array2DReal &NormalVelocity) const;

 private:
   Array2DI4 CellsOnEdge;
};
```

### 4.2 Methods

#### 4.2.1 Constructor
The constructor will be responsible for:
  * allocating auxiliary variables
  * retrieving and storing the configuration options if applicable
  * registering fields and metadata with the I/O infrastructure.

```c++
   LayerThicknessAuxVars(const HorzMesh *mesh, int NVertLevels, const Config *Options)
       : FluxLayerThickEdge("FluxLayerThickEdge", mesh->NEdgesSize,
                            NVertLevels),
         MeanLayerThickEdge("MeanLayerThickEdge", mesh->NEdgesSize,
                            NVertLevels),
         CellsOnEdge(mesh->CellsOnEdge) {

             std::string FluxThickTypeStr;
             Config->get("fluxThicknessType", FluxThickTypeStr);
             switch (FluxThickTypeStr) {
                 case "Center":
                 FluxThickType = Center;
                 break;
                 case "Upwind":
                 FluxThickType = Upwind;
                 break;
                 default:
                 FluxThickType = Upwind;
                 break;
             }

             defineIOFields();
    }
```

#### 4.2.2 Compute methods
Compute methods implement the auxiliary variables computations for a chunk of
vertical levels at a given horizontal mesh location. The mesh location is
indicated in the method name. There may be more than one compute method to
compute different groups of variables over different mesh locations. Any
configurable computation options are handled inside the compute method. The
inner loop over a chunk of vertical levels enables CPU vectorization.

```c++
   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &LayerThickCell,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = KChunk * VecLength;
      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         MeanLayerThickEdge(IEdge, K) =
             0.5_Real * (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
      }

      switch (FluxThickEdgeChoice) {
      case Center:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            FluxLayerThickEdge(IEdge, K) =
                0.5_Real *
                (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
         }
         break;
      case Upwind:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            if (NormalVelEdge(IEdge, K) > 0) {
               FluxLayerThickEdge(IEdge, K) = LayerThickCell(JCell0, K);
            } else if (NormalVelEdge(IEdge, K) < 0) {
               FluxLayerThickEdge(IEdge, K) = LayerThickCell(JCell1, K);
            } else {
               FluxLayerThickEdge(IEdge, K) = Kokkos::max(
                   LayerThickCell(JCell0, K), LayerThickCell(JCell1, K));
            }
         }
         break;
      }
   }
```

## 5 Verification and Testing

### 5.1 Unit testing
Each auxiliary variable will be used in a CTest which will compare its result
against a reference value.

### 5.2 Convergence testing
The auxiliary variables for the shallow water equations will be tested on
existing test cases implemented in Polaris to verify they produce the correct
convergence rates.
