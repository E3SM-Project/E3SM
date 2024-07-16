#include "VorticityAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

VorticityAuxVars::VorticityAuxVars(const std::string &AuxStateSuffix,
                                   const HorzMesh *Mesh, int NVertLevels)
    : RelVortVertex("RelVortVertex" + AuxStateSuffix, Mesh->NVerticesSize,
                    NVertLevels),
      NormRelVortVertex("NormRelVortVertex" + AuxStateSuffix,
                        Mesh->NVerticesSize, NVertLevels),
      NormPlanetVortVertex("NormPlanetVortVertex" + AuxStateSuffix,
                           Mesh->NVerticesSize, NVertLevels),
      NormRelVortEdge("NormRelVortEdge" + AuxStateSuffix, Mesh->NEdgesSize,
                      NVertLevels),
      NormPlanetVortEdge("NormPlanetVortEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, NVertLevels),
      VertexDegree(Mesh->VertexDegree), CellsOnVertex(Mesh->CellsOnVertex),
      EdgesOnVertex(Mesh->EdgesOnVertex),
      EdgeSignOnVertex(Mesh->EdgeSignOnVertex), DcEdge(Mesh->DcEdge),
      KiteAreasOnVertex(Mesh->KiteAreasOnVertex),
      AreaTriangle(Mesh->AreaTriangle), FVertex(Mesh->FVertex),
      VerticesOnEdge(Mesh->VerticesOnEdge) {}

void VorticityAuxVars::registerFields(const std::string &AuxGroupName) const {
   addMetaData(AuxGroupName);
   defineIOFields();
}

void VorticityAuxVars::unregisterFields() const {
   IOField::erase(RelVortVertex.label());
   IOField::erase(NormRelVortVertex.label());
   IOField::erase(NormPlanetVortVertex.label());
   IOField::erase(NormRelVortEdge.label());
   IOField::erase(NormPlanetVortEdge.label());
   MetaData::destroy(RelVortVertex.label());
   MetaData::destroy(NormRelVortVertex.label());
   MetaData::destroy(NormPlanetVortVertex.label());
   MetaData::destroy(NormRelVortEdge.label());
   MetaData::destroy(NormPlanetVortEdge.label());
}

void VorticityAuxVars::addMetaData(const std::string &AuxGroupName) const {
   auto VertexDim    = MetaDim::get("NVertices");
   auto EdgeDim      = MetaDim::get("NEdges");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get(AuxGroupName);

   const Real FillValue = -9.99e30;

   // Relative vorticity on vertices
   auto RelVortVertexMeta = ArrayMetaData::create(
       RelVortVertex.label(),
       "curl of horizontal velocity, defined at vertices", /// long Name or
                                                           /// description
       "s^-1",                                             /// units
       "ocean_relative_vorticity",                         /// CF standard Name
       std::numeric_limits<Real>::min(),                   /// min valid value
       std::numeric_limits<Real>::max(),                   /// max valid value
       FillValue,           /// scalar used for undefined entries
       2,                   /// number of dimensions
       {VertexDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(RelVortVertex.label());

   // Normalized relative vorticity on vertices
   auto NormRelVortVertexMeta = ArrayMetaData::create(
       NormRelVortVertex.label(),
       "curl of horizontal velocity divided by layer thickness", /// long Name
                                                                 /// or
                                                                 /// description
       "m^-1 s^-1",                                              /// units
       "",                               /// CF standard Name
       std::numeric_limits<Real>::min(), /// min valid value
       std::numeric_limits<Real>::max(), /// max valid value
       FillValue,                        /// scalar used for undefined entries
       2,                                /// number of dimensions
       {VertexDim, VertDim}              /// dim pointers
   );
   AuxMetaGroup->addField(NormRelVortVertex.label());

   // Normalized planetary vorticity on vertices
   auto NormPlanetVortVertexMeta = ArrayMetaData::create(
       NormPlanetVortVertex.label(),
       "earth's rotational rate (Coriolis parameter, f) divided by layer "
       "thickness",                      /// long Name or description
       "m^-1 s^-1",                      /// units
       "",                               /// CF standard Name
       std::numeric_limits<Real>::min(), /// min valid value
       std::numeric_limits<Real>::max(), /// max valid value
       FillValue,                        /// scalar used for undefined entries
       2,                                /// number of dimensions
       {VertexDim, VertDim}              /// dim pointers
   );
   AuxMetaGroup->addField(NormPlanetVortVertex.label());

   // Normalized relative vorticity on edges
   auto NormRelVortEdgeMeta = ArrayMetaData::create(
       NormRelVortEdge.label(),
       "curl of horizontal velocity divided by layer thickness, averaged from "
       "vertices to edges",              /// long Name or description
       "m^-1 s^-1",                      /// units
       "",                               /// CF standard Name
       std::numeric_limits<Real>::min(), /// min valid value
       std::numeric_limits<Real>::max(), /// max valid value
       FillValue,                        /// scalar used for undefined entries
       2,                                /// number of dimensions
       {EdgeDim, VertDim}                /// dim pointers
   );
   AuxMetaGroup->addField(NormRelVortEdge.label());

   // Normalized planetary vorticity on edges
   auto NormPlanetVortEdgeMeta = ArrayMetaData::create(
       NormPlanetVortEdge.label(),
       "earth's rotational rate (Coriolis parameter, f) divided by layer "
       "thickness, averaged from vertices to edges", /// long Name or
                                                     /// description
       "m^-1 s^-1",                                  /// units
       "",                                           /// CF standard Name
       std::numeric_limits<Real>::min(),             /// min valid value
       std::numeric_limits<Real>::max(),             /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {EdgeDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(NormPlanetVortEdge.label());
}

void VorticityAuxVars::defineIOFields() const {
   int Err;

   // Relative vorticity on vertices
   Err = IOField::define(RelVortVertex.label());
   Err = IOField::attachData(RelVortVertex.label(), RelVortVertex);

   // Normalized relative vorticity on vertices
   Err = IOField::define(NormRelVortVertex.label());
   Err = IOField::attachData(NormRelVortVertex.label(), NormRelVortVertex);

   // Normalized planetary vorticity on vertices
   Err = IOField::define(NormPlanetVortVertex.label());
   Err =
       IOField::attachData(NormPlanetVortVertex.label(), NormPlanetVortVertex);

   // Normalized relative vorticity on edges
   Err = IOField::define(NormRelVortEdge.label());
   Err = IOField::attachData(NormRelVortEdge.label(), NormRelVortEdge);

   // Normalized planetary vorticity on edges
   Err = IOField::define(NormPlanetVortEdge.label());
   Err = IOField::attachData(NormPlanetVortEdge.label(), NormPlanetVortEdge);
}

} // namespace OMEGA
