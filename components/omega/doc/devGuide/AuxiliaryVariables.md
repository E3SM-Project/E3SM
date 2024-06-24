(omega-dev-aux-vars)=

# Auxiliary Variables

Auxiliary variables are variables that need to be computed at every model time
step in order to advance the model state. Similar auxiliary variables are
grouped together in auxiliary variable groups. An auxiliary variable group is
implemented as a C++ class.

## Basic class structure

Each auxiliary group class has the same basic structure. It is a light-weight
class, since it needs to be copied to the GPU to perform computations. Each
class contains:
 - a constructor
 - public arrays for the auxiliary variables
 - public member function(s) that compute the variables
 - member functions for adding metadata and registering with IO
 - static strings that contain the names used for metadata and IO

Additionally, there may be enum members that store user-configurable options.

As an example, the class `KineticAuxVars` that groups kinetic energy and
velocity divergence has the following structure:
```c++
class KineticAuxVars {
 public:
   Array2DReal KineticEnergyCell;
   Array2DReal VelocityDivCell;

   KineticAuxVars(const HorzMesh *mesh, int NVertLevels);
   void addMetaData() const;
   void defineIOFields() const;

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int KChunk,
                     const Array2DReal &NormalVelEdge) const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;

   // names used in defining MetaData and IOFields
   static const std::string KineticEnergyCellName;
   static const std::string VelocityDivCellName;
};
```

## Configuration enums

Some auxiliary variables have user-configurable options. In that case, for each
option an enum is defined in the auxiliary group header file. For example, the
flux thickness choice is represented using
```c++
enum FluxThickEdgeOption { Center, Upwind };
```

## Constructor

The constructor takes in a mesh, and does the following:
- allocates the auxiliary variable arrays
- grabs the neccessary data from the mesh
- retrieves configuration options
- calls the `addMetaData()` and `defineIOFields()` helper functions

## Compute member functions

The compute functions operate on a vertical chunk of mesh elements and take in
the element index, the chunk index, and input arrays, which may contain state
variables or other auxiliary variables. There is no output argument, they write
their results directly to the member arrays. There might be more than one
compute function, if variables located on different mesh elements are grouped
together. For example, the `VorticityAuxVars` class provides both
```c++
void computeVarsOnVertex(int IVertex, int KChunk, const Array2DReal &LayerThickCell, const Array2DReal &NormalVelEdge);
void computeVarsOnEdge(int IEdge, int KChunk);
```
because the vorticity variables are first computed on vertices and then
averaged to edges. Note that the edge compute function doesn't take any array
arguments, because it operates only on the class member arrays.

## Implemented groups
The following auxiliary variable groups are currently implemented:
| Group | Auxiliary Variable | Available options |
| ----- | ------------------- | ------- |
| KineticAuxVars | KineticEnergyCell ||
|| VelocityDiv ||
| LayerThicknessAuxVars | FluxLayerThickEdge | Center or Upwind|
|| MeanLayerThickEdge ||
| VorticityAuxVars |  RelVortVertex ||
||  NormRelVortVertex ||
||  NormPlanetVortVertex ||
||  NormRelVortEdge ||
||  NormPlanetVortEdge ||
