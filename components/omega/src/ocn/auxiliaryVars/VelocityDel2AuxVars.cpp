#include "VelocityDel2AuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

VelocityDel2AuxVars::VelocityDel2AuxVars(const std::string &AuxStateSuffix,
                                         const HorzMesh *Mesh, int NVertLevels)
    : Del2Edge("VelDel2Edge" + AuxStateSuffix, Mesh->NEdgesSize, NVertLevels),
      Del2DivCell("VelDel2DivCell" + AuxStateSuffix, Mesh->NCellsSize,
                  NVertLevels),
      Del2RelVortVertex("VelDel2RelVortVertex" + AuxStateSuffix,
                        Mesh->NVerticesSize, NVertLevels),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgesOnVertex(Mesh->EdgesOnVertex), CellsOnEdge(Mesh->CellsOnEdge),
      VerticesOnEdge(Mesh->VerticesOnEdge),
      EdgeSignOnVertex(Mesh->EdgeSignOnVertex),
      AreaTriangle(Mesh->AreaTriangle), VertexDegree(Mesh->VertexDegree) {}

void VelocityDel2AuxVars::registerFields(const std::string &AuxGroupName,
                                         const std::string &MeshName) const {

   int Err = 0; // Error flag for some calls

   // Create/define fields
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[1] = "NVertLevels";
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Del2 vel on edges
   DimNames[0]        = "NEdges" + DimSuffix;
   auto Del2EdgeField = Field::create(
       Del2Edge.label(),                            // field name
       "laplacian of horizontal velocity on edges", // long Name or description
       "m^-1 s^-1",                                 // units
       "",                                          // CF standard Name
       std::numeric_limits<Real>::min(),            // min valid value
       std::numeric_limits<Real>::max(),            // max valid value
       FillValue,                                   // scalar for undef entries
       NDims,                                       // number of dimensions
       DimNames                                     // dimension names
   );

   // Del2 vel on cells
   DimNames[0] = "NCells" + DimSuffix;
   auto Del2DivCellField =
       Field::create(Del2DivCell.label(), // field name
                     "divergence of laplacian of horizontal velocity "
                     "on cells",  // long Name or description
                     "m^-2 s^-1", // units
                     "",          // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar used for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dimension names
       );

   // Del2 of relative vorticity on vertices
   DimNames[0]                 = "NVertices" + DimSuffix;
   auto Del2RelVortVertexField = Field::create(
       Del2RelVortVertex.label(),                     // field name
       "laplacian of relative vorticity at vertices", // long name, description
       "m^-2 s^-1",                                   // units
       "",                                            // CF standard Name
       std::numeric_limits<Real>::min(),              // min valid value
       std::numeric_limits<Real>::max(),              // max valid value
       FillValue, // scalar used for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   // Add fields to Aux Field group
   Err = FieldGroup::addFieldToGroup(Del2Edge.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", Del2Edge.label(),
                AuxGroupName);

   Err = FieldGroup::addFieldToGroup(Del2DivCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", Del2DivCell.label(),
                AuxGroupName);

   Err = FieldGroup::addFieldToGroup(Del2RelVortVertex.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", Del2RelVortVertex.label(),
                AuxGroupName);

   // Attach data to fields
   Err = Del2EdgeField->attachData<Array2DReal>(Del2Edge);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", Del2Edge.label());

   Err = Del2DivCellField->attachData<Array2DReal>(Del2DivCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", Del2DivCell.label());

   Err = Del2RelVortVertexField->attachData<Array2DReal>(Del2RelVortVertex);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", Del2RelVortVertex.label());
}

void VelocityDel2AuxVars::unregisterFields() const {
   int Err = 0;

   Err = Field::destroy(Del2Edge.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", Del2Edge.label());

   Err = Field::destroy(Del2DivCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", Del2DivCell.label());

   Err = Field::destroy(Del2RelVortVertex.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", Del2RelVortVertex.label());
}

} // namespace OMEGA
