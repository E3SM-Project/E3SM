(omega-dev-horz-operators)=

# Horizontal Operators

The TRiSK scheme horizontal operators are implemented in C++ as functors, which
are classes that can be called with input arguments in the same way as functions
can. However, in contrast to plain functions, functors can contain internal
state. In the case of Omega operators they contain shallow copies of various
`HorzMesh` arrays that are needed to compute the operator, but are not strictly
inputs of the operator.

Operators are constructed from an instance of the `HorzMesh` class, for example
```c++
    auto mesh = OMEGA::HorzMesh::getDefault();
    DivergenceOnCell DivOnCell(mesh);
```
which sets up the previously mentioned internal state.

Each Omega operator provides a C++ call operator, which computes its value on a
vertical chunk of mesh elements given the element index, the chunk index,
and operator-specific input arrays. The first argument to an operator is
an array that gets updated with the computed values. Typically, operators are
created outside of a parallel region and are used inside a parallel loop
over mesh elements, for example
```c++
    auto mesh = OMEGA::HorzMesh::getDefault();
    DivergenceOnCell DivOnCell(mesh);
    parallelFor({mesh->NCellsOwned, NVertLevels / VecLength}, KOKKOS_LAMBDA(int ICell, int KChunk) {
        // computes divergence of Vec for cells with indices (ICell, KChunk:KChunk+VecLength-1)
        // stores the result in DivVec
        DivOnCell(DivVec, ICell, KChunk, Vec);
    });
```

Currently, the following operators are implemented:
- `DivergenceOnCell`
- `GradientOnEdge`
- `CurlOnVertex`
- `TangentialReconOnEdge`

Some tendency terms in the Omega PDE solver could in principle be constructed
using these operators as building blocks. However, very often tendency terms
require evaluation of slightly modified operators. Moreover, there is a
potential performance cost of nesting operators within classes. Hence, the
primary motivation for introducing these classes is to provide reference
implementation of the basic TRiSK operators and for diagnostic and debugging
purposes.
