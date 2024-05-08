(omega-design-tendency)=
# Tendency

## 1 Overview

The tendency terms in OMEGA are implemented as functors, which define an operation over a number of vertical levels for particular a cell, edge, or vertex.
Tendency functors take information which remains constant during the forward simulation (such as `HorzMesh` and `Config`)  as constructor arguments.
The `operator()` method is overloaded with the relevant discrete /parameterization.
This approach allows for a modularization of the tendency terms that enables flexible groupings of work within larger cell/edge/vertex loops. 
This type of flexibility will be important for optimizing the code across different device architectures.

## 2 Requirements

### 2.1 Requirement: Tendencies should be able to be grouped under different outer loops over mesh locations.
The ability to fuse/separate tendency computations will allow for optimization across different device architectures.
Each functor will implement a given tendency operation on a specific mesh location (cell, edge, vertex).

### 2.2 Requirement: Tendency operator calls should have simple, intuitive calling arguments.
Tendency functors take in constant data as constructor arguments, which simplifies the arguments used to call the operator method.

### 2.3 Requirement: Tendencies must allow for vectorization on CPU architectures.
Tendency operations have inner loops over a chunk of vertical levels.
The chunk size will be set to the vector length on CPU machines and 1 for GPUs.
This will allow for the possibility of vectorization on CPUs.

## 3 Algorithmic Formulation

The algorithmic formulation for each tendency can be found in the algorithms design document.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters

#### 4.1.2 Class/structs/data types
The tendency functor is a class the overrides the `operator()` method.
The class contains private member variables for any constant data that are defined by the constructor. 
Each tendency term will have a `bool` variable which can be set to enable/disable the computation of the tendency.

```c++
class ThicknessFluxDivergenceOnCell {
  public:
 
    bool enabled = false;
 
    ThicknessFluxDivergenceOnCell(const HorzMesh *Mesh, Config *Options);
  
    KOKKOS_INLINE_FUNCTION Real operator()(int ICell,
                                           int KChunk,
                                           const OceanState *State,
                                           const OceanAuxState *AuxState,
                                           Array2DReal &Tend);
  
  private:
    Array1DI4 NEdgesOnCell;
    Array2DI4 EdgesOnCell;
    Array1DR8 DvEdge;
    Array1DR8 AreaCell;
    Array2DR8 EdgeSignOnCell;
};
```

### 4.2 Methods

#### 4.2.1 Constructor
Tendency functor constructors are responsible for initializing the private member variables:

```c++
ThicknessFluxDivergenceOnCell::ThicknessFluxDivergenceOnCell(HorzMesh const *Mesh, Config *Options)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) {

    enabled = Options->get("ThicknessFluxTendencyEnable")
}
```

#### 4.2.2 operator
The operator method implements the tendency computation for a chunk of vertical levels at a given horizontal mesh location.
The inner loop over a chunk of vertical levels enables CPU vectorization.

```c++
KOKKOS_INLINE_FUNCTION Real ThicknessFluxDivergenceOnCell::operator()(int ICell,
                                                                      int KChunk,
                                                                      const OceanState *State,
                                                                      const OceanAuxState *AuxState,
                                                                      Array2DReal &Tend)  const {
   const Real InvAreaCell = 1. / AreaCell(iCell);
   for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
      const int JEdge = EdgesOnCell(ICell, J);
      for (int K = KChunk * LevelsPerChunk; K < (KChunk + 1) * LevelsPerChunk; ++K) {
         Tend(JEdge,K) -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) * AuxState->ThicknessFlux(JEdge, K) * InvAreaCell;
      }
   }
}
```

## 5 Verification and Testing

### 5.1 Unit testing
Each tendency operator will be used in a CTest which will compare its result against a reference value.

### 5.2 Convergence testing
The tendency terms for the shallow water equations will be tested on existing test cases implemented in Polaris to verify they produce the correct convergence rates.

