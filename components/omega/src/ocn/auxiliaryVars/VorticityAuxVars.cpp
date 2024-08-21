#include "VorticityAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

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

void VorticityAuxVars::registerFields(const std::string &AuxGroupName,
                                      const std::string &MeshName) const {

   int Err = 0; // error code for some calls

   // Create fields with metadata
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NVertices" + DimSuffix; // for first three fields
   DimNames[1] = "NVertLevels";           // same for all fields below

   // Relative vorticity on vertices
   auto RelVortVertexField = Field::create(
       RelVortVertex.label(),                              // field name
       "curl of horizontal velocity, defined at vertices", // long name/describe
       "s^-1",                                             // units
       "ocean_relative_vorticity",                         // CF standard Name
       std::numeric_limits<Real>::min(),                   // min valid value
       std::numeric_limits<Real>::max(),                   // max valid value
       FillValue, // scalar for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   // Normalized relative vorticity on vertices
   auto NormRelVortVertexField = Field::create(
       NormRelVortVertex.label(),                                // field name
       "curl of horizontal velocity divided by layer thickness", // long Name
       "m^-1 s^-1",                                              // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar used for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Normalized planetary vorticity on vertices
   auto NormPlanetVortVertexField = Field::create(
       NormPlanetVortVertex.label(), // field name
       "earth's rotational rate (Coriolis parameter, f) divided by layer "
       "thickness",                      // long Name or description
       "m^-1 s^-1",                      // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar used for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Normalized relative vorticity on edges
   DimNames[0]               = "NEdges" + DimSuffix; // for last two fields
   auto NormRelVortEdgeField = Field::create(
       NormRelVortEdge.label(), // field name
       "curl of horizontal velocity divided by layer thickness, averaged from "
       "vertices to edges",              // long Name or description
       "m^-1 s^-1",                      // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar used for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Normalized planetary vorticity on edges
   auto NormPlanetVortEdgeField = Field::create(
       NormPlanetVortEdge.label(), // field name
       "earth's rotational rate (Coriolis parameter, f) divided by layer "
       "thickness, averaged from vertices to edges", // long Name or description
       "m^-1 s^-1",                                  // units
       "",                                           // CF standard Name
       std::numeric_limits<Real>::min(),             // min valid value
       std::numeric_limits<Real>::max(),             // max valid value
       FillValue, // scalar used for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   // Add fields to Aux field group
   Err = FieldGroup::addFieldToGroup(RelVortVertex.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", RelVortVertex.label(),
                AuxGroupName);

   Err = FieldGroup::addFieldToGroup(NormRelVortVertex.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", NormRelVortVertex.label(),
                AuxGroupName);

   Err =
       FieldGroup::addFieldToGroup(NormPlanetVortVertex.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}",
                NormPlanetVortVertex.label(), AuxGroupName);

   Err = FieldGroup::addFieldToGroup(NormRelVortEdge.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", NormRelVortEdge.label(),
                AuxGroupName);

   Err = FieldGroup::addFieldToGroup(NormPlanetVortEdge.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", NormPlanetVortEdge.label(),
                AuxGroupName);

   // Attach data to fields
   Err = RelVortVertexField->attachData<Array2DReal>(RelVortVertex);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", RelVortVertex.label());

   Err = NormRelVortVertexField->attachData<Array2DReal>(NormRelVortVertex);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", NormRelVortVertex.label());

   Err =
       NormPlanetVortVertexField->attachData<Array2DReal>(NormPlanetVortVertex);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}",
                NormPlanetVortVertex.label());

   Err = NormRelVortEdgeField->attachData<Array2DReal>(NormRelVortEdge);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", NormRelVortEdge.label());

   Err = NormPlanetVortEdgeField->attachData<Array2DReal>(NormPlanetVortEdge);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", NormPlanetVortEdge.label());
}

void VorticityAuxVars::unregisterFields() const {
   int Err = 0;

   Err = Field::destroy(RelVortVertex.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", RelVortVertex.label());

   Err = Field::destroy(NormRelVortVertex.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", NormRelVortVertex.label());

   Err = Field::destroy(NormPlanetVortVertex.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", NormPlanetVortVertex.label());

   Err = Field::destroy(NormRelVortEdge.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", NormRelVortEdge.label());

   Err = Field::destroy(NormPlanetVortEdge.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", NormPlanetVortEdge.label());
}

} // namespace OMEGA
