(omega-design-tendency)=
# Tendency

## 1 Overview

The tendency terms in OMEGA are implemented as functors, which define an operation over a number of vertical layers for particular a cell, edge, or vertex.
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
Tendency operations have inner loops over a chunk of vertical layers.
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
class ThicknessFluxDivOnCell {
  public:

    bool enabled = false;

    ThicknessFluxDivOnCell(const HorzMesh *Mesh, Config *Options);

    KOKKOS_FUNCTION void operator()(Array2DReal &Tend
                                    int ICell,
                                    int KChunk,
                                    const Array2DReal &ThicknessFlux,
                                    const Array2DReal &NormalVelEdge);

  private:
    Array1DI4 NEdgesOnCell;
    Array2DI4 EdgesOnCell;
    Array1DReal DvEdge;
    Array1DReal AreaCell;
    Array2DReal EdgeSignOnCell;
};
```

### 4.2 Methods

#### 4.2.1 Constructor
Tendency functor constructors are responsible for initializing the private member variables:

```c++
ThicknessFluxDivOnCell::ThicknessFluxDivOnCell(HorzMesh const *Mesh)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) { }
```

#### 4.2.2 operator
The operator method implements the tendency computation for a chunk of vertical layers at a given horizontal mesh location.
The inner loop over a chunk of vertical layers enables CPU vectorization.

```c++
KOKKOS_FUNCTION void ThicknessFluxDivOnCell::operator()(Array2DReal &Tend
                                                        int ICell,
                                                        int KChunk,
                                                        const Array2DReal &ThicknessFlux,
                                                        const Array2DReal &NormalVelEdge)  const {
    const I4 KStart        = KChunk * VecLength;
    const Real InvAreaCell = 1._Real / AreaCell(ICell);
    Real DivTmp[VecLength] = {0};
    for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
       const I4 JEdge = EdgesOnCell(ICell, J);
       for (int KVec = 0; KVec < VecLength; ++KVec) {
          const I4 K = KStart + KVec;
          DivTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                          ThicknessFlux(JEdge, K) * NormalVelEdge(JEdge, K) *
                          InvAreaCell;
       }
    }
    for (int KVec = 0; KVec < VecLength; ++KVec) {
       const I4 K = KStart + KVec;
       Tend(ICell, K) -= DivTmp[KVec];
    }
}
```

## 5 Verification and Testing

### 5.1 Unit testing
Each tendency operator will be used in a CTest which will compare its result against a reference value.

### 5.2 Convergence testing
The tendency terms for the shallow water equations will be tested on existing test cases implemented in Polaris to verify they produce the correct convergence rates.
