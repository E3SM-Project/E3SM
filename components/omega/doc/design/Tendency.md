(omega-design-tendency)=
# Tendency

## 1 Overview

The tendency terms in OMEGA are implemented as functors, which define an operation on a particular vertical level of a cell, edge, or vertex.
Tendency functors take information which remains constant during the forward simulation (such as the `HorzMesh` and config options)  as constructor arguments.
The `operator()` method is overloaded with the relevant discrete /parameterization.
This approach allows for a modularization of the tendency terms that enables flexible groupings of work within larger cell/edge/vertex loops. 
This type of flexibility will be important for optimizing the code across different device architectures.
An "orchestrator" will be responsible for grouping tendency calculations inside loops over difference mesh locations.

## 2 Requirements

### 2.1 Requirement: Tendencies should be able to be grouped under different outer loops over mesh locations.
The ability to fuse/separate tendency computations will allow for optimization across different device architectures.
Each functor will implement a given tendency operation on a specific mesh location (cell, edge, vertex).

### 2.2 Requirement: Tendency operator calls should have simple, intuitive calling arguments.
Tendency functors take in constant data as constructor arguments, which simplifies the arguments used to call the operator method.

## 3 Algorithmic Formulation

The algorithmic formulation for each tendency can be found in the algorithms design document.

## 4 Design

### 4.1 Data types and parameters

#### 4.1.1 Parameters
Each tendency functor will have a `name` and `locaion` parameter associated with it which will be used to identify it within the orchestrator.

#### 4.1.2 Class/structs/data types
The tendency functor is a class the overrides the `operator()` method.
The class contains private member variables for any constant data that are defined by the constructor. 

```c++
class ThicknessFluxDivergenceOnCell {
  public:
    std::string Name;
    std::string MeshLocation = "Cell";
  
    ThicknessFluxDivergenceOnCell(HorzMesh const *mesh, std::string TendName="ThicknessFluxDivergenceOnCell");
  
    KOKKOS_INLINE_FUNCTION Real operator()(int ICell,
                                           int KLevel,
                                           const Array2DReal &ThicknessFlux);
  
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
ThicknessFluxDivergenceOnCell::ThicnessFluxDivergenceOnCell(HorzMesh const *mesh, std::string TendName="ThicknessFluxDivergenceOnCell")
    : NEdgesOnCell(mesh->NEdgesOnCell), EdgesOnCell(mesh->EdgesOnCell),
      DvEdge(mesh->DvEdge), AreaCell(mesh->AreaCell),
      EdgeSignOnCell(mesh->EdgeSignOnCell), Name(TendName) {}
```

#### 4.2.2 operator
The operator method implements the tendency computation for a given vertical level/horizontal mesh location

```c++
KOKKOS_INLINE_FUNCTION Real ThicknessFluxDivergenceOnCell::operator()(int ICell,
                                                                      int KLevel,
                                                                      const Array2DReal &ThicknessFlux) const {
   Real DivCell = 0;
   for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
      const int JEdge = EdgesOnCell(ICell, J);
      DivCell -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) * ThicknessFlux(JEdge, KLevel);
   }
   const Real InvAreaCell = 1. / AreaCell(ICell);
   DivCell *= InvAreaCell;
   return DivCell;
}
```

## 5 Verification and Testing

### 5.1 Unit testing
Each tendency operator will be used in a CTest which will compare its result against a reference value.

### 5.2 Convergence testing
The tendency terms for the shallow water equations will be tested on existing test cases implemented in Polaris to verify they produce the correct convergence rates.

