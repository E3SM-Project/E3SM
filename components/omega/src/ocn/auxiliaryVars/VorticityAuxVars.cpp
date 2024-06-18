#include "VorticityAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

const std::string VorticityAuxVars::RelVortVertexName     = "RelVortVertex";
const std::string VorticityAuxVars::NormRelVortVertexName = "NormRelVortVertex";
const std::string VorticityAuxVars::NormPlanetVortVertexName =
    "NormPlanetVortVertex";
const std::string VorticityAuxVars::NormRelVortEdgeName = "NormRelVortEdge";
const std::string VorticityAuxVars::NormPlanetVortEdgeName =
    "NormPlanetVortEdge";

VorticityAuxVars::VorticityAuxVars(const HorzMesh *mesh, int NVertLevels)
    : RelVortVertex("RelVortVertex", mesh->NVerticesSize, NVertLevels),
      NormRelVortVertex("NormRelVortVertex", mesh->NVerticesSize, NVertLevels),
      NormPlanetVortVertex("NormPlanetVortVertex", mesh->NVerticesSize,
                           NVertLevels),
      NormRelVortEdge("NormRelVortEdge", mesh->NEdgesSize, NVertLevels),
      NormPlanetVortEdge("NormPlanetVortEdge", mesh->NEdgesSize, NVertLevels),
      VertexDegree(mesh->VertexDegree), CellsOnVertex(mesh->CellsOnVertex),
      EdgesOnVertex(mesh->EdgesOnVertex),
      EdgeSignOnVertex(mesh->EdgeSignOnVertex), DcEdge(mesh->DcEdge),
      KiteAreasOnVertex(mesh->KiteAreasOnVertex),
      AreaTriangle(mesh->AreaTriangle), FVertex(mesh->FVertex),
      VerticesOnEdge(mesh->VerticesOnEdge) {

   addMetaData();
   defineIOFields();
}

void VorticityAuxVars::addMetaData() const {
   auto VertexDim    = MetaDim::get("NVertices");
   auto EdgeDim      = MetaDim::get("NEdges");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get("auxiliaryVars");

   const Real FillValue = -9.99e30;

   // Relative vorticity on vertices
   auto RelVortVertexMeta = ArrayMetaData::create(
       RelVortVertexName,
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
   AuxMetaGroup->addField(RelVortVertexName);

   // Normalized relative vorticity on vertices
   auto NormRelVortVertexMeta = ArrayMetaData::create(
       NormRelVortVertexName,
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
   AuxMetaGroup->addField(NormRelVortVertexName);

   // Normalized planetary vorticity on vertices
   auto NormPlanetVortVertexMeta = ArrayMetaData::create(
       NormPlanetVortVertexName,
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
   AuxMetaGroup->addField(NormPlanetVortVertexName);

   // Normalized relative vorticity on edges
   auto NormRelVortEdgeMeta = ArrayMetaData::create(
       NormRelVortEdgeName,
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
   AuxMetaGroup->addField(NormRelVortEdgeName);

   // Normalized planetary vorticity on edges
   auto NormPlanetVortEdgeMeta = ArrayMetaData::create(
       NormPlanetVortEdgeName,
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
   AuxMetaGroup->addField(NormPlanetVortEdgeName);
}

void VorticityAuxVars::defineIOFields() const {
   int Err;

   // Relative vorticity on vertices
   Err = IOField::define(RelVortVertexName);
   Err = IOField::attachData(RelVortVertexName, RelVortVertex);

   // Normalized relative vorticity on vertices
   Err = IOField::define(NormRelVortVertexName);
   Err = IOField::attachData(NormRelVortVertexName, NormRelVortVertex);

   // Normalized planetary vorticity on vertices
   Err = IOField::define(NormPlanetVortVertexName);
   Err = IOField::attachData(NormPlanetVortVertexName, NormPlanetVortVertex);

   // Normalized relative vorticity on edges
   Err = IOField::define(NormRelVortEdgeName);
   Err = IOField::attachData(NormRelVortEdgeName, NormRelVortEdge);

   // Normalized planetary vorticity on edges
   Err = IOField::define(NormPlanetVortEdgeName);
   Err = IOField::attachData(NormPlanetVortEdgeName, NormPlanetVortEdge);
}

} // namespace OMEGA
