(omega-dev-tend-terms)=

# Tendency Terms

The tendency terms are each implemented as a separate functor. They contain
private member variables to store any constant data needed for computation.
This data is initialized during construction and comes from either an instance
of the `HorzMesh` class or a `Config` object, which are passed as arguments to
the constructor:
```c++
   auto Mesh = OMEGA::HorzMesh::getDefault();
   auto Options = OMEGA::Config::getOmegaConfig();
   OMEGA::ThicknessFluxDivOnCell ThickFluxDivOnC(Mesh, Options));
```

Each functor is called like a function to compute tendency values for a range
of vertical layers at a given mesh element. So, the functors all take as
input the mesh element index and the index for the vertical chunk, along with
arrays computed by [Auxillary Variables](#omega-dev-aux-vars) classes. In each
case, the first argument is a tendency array that gets updated by the
contributon to the tendency for the given term. The functors are designed to be
used inside a Kokkos parallel loop:
```c++
   OMEGA::parallelFor(
       {mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int KChunk) {
          ThickFluxDixOnC(ThicknessTend, ICell, KChunk, ThickFluxEdge);
       });
```

The following tendency terms for the shallow water equations are currently
implemented:
- `ThicknessFluxDivOnCell`
- `PotentialVortHAdvOnEdge`
- `KEGradOnEdge`
- `SSHGradOnEdge`
- `VelocityDiffusionOnEdge`
- `VelocityHyperDiffOnEdge`
