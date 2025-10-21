#include "VorticityAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

VorticityAuxVars::VorticityAuxVars(const std::string &AuxStateSuffix,
                                   const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : RelVortVertex("RelVortVertex" + AuxStateSuffix, Mesh->NVerticesSize,
                    VCoord->NVertLayers),
      NormRelVortVertex("NormRelVortVertex" + AuxStateSuffix,
                        Mesh->NVerticesSize, VCoord->NVertLayers),
      NormPlanetVortVertex("NormPlanetVortVertex" + AuxStateSuffix,
                           Mesh->NVerticesSize, VCoord->NVertLayers),
      NormRelVortEdge("NormRelVortEdge" + AuxStateSuffix, Mesh->NEdgesSize,
                      VCoord->NVertLayers),
      NormPlanetVortEdge("NormPlanetVortEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, VCoord->NVertLayers),
      VertexDegree(Mesh->VertexDegree), CellsOnVertex(Mesh->CellsOnVertex),
      EdgesOnVertex(Mesh->EdgesOnVertex),
      EdgeSignOnVertex(Mesh->EdgeSignOnVertex), DcEdge(Mesh->DcEdge),
      KiteAreasOnVertex(Mesh->KiteAreasOnVertex),
      AreaTriangle(Mesh->AreaTriangle), FVertex(Mesh->FVertex),
      VerticesOnEdge(Mesh->VerticesOnEdge) {}

void VorticityAuxVars::registerFields(const std::string &AuxGroupName,
                                      const std::string &MeshName) const {

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
   DimNames[1] = "NVertLayers";           // same for all fields below

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
   FieldGroup::addFieldToGroup(RelVortVertex.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(NormRelVortVertex.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(NormPlanetVortVertex.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(NormRelVortEdge.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(NormPlanetVortEdge.label(), AuxGroupName);

   // Attach data to fields
   RelVortVertexField->attachData<Array2DReal>(RelVortVertex);
   NormRelVortVertexField->attachData<Array2DReal>(NormRelVortVertex);
   NormPlanetVortVertexField->attachData<Array2DReal>(NormPlanetVortVertex);
   NormRelVortEdgeField->attachData<Array2DReal>(NormRelVortEdge);
   NormPlanetVortEdgeField->attachData<Array2DReal>(NormPlanetVortEdge);
}

void VorticityAuxVars::unregisterFields() const {
   Field::destroy(RelVortVertex.label());
   Field::destroy(NormRelVortVertex.label());
   Field::destroy(NormPlanetVortVertex.label());
   Field::destroy(NormRelVortEdge.label());
   Field::destroy(NormPlanetVortEdge.label());
}

} // namespace OMEGA
