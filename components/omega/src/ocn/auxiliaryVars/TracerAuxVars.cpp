#include "TracerAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

TracerAuxVars::TracerAuxVars(const std::string &AuxStateSuffix,
                             const HorzMesh *Mesh, const I4 NVertLevels,
                             const I4 NTracers)
    : HTracersOnEdge("ThickTracersOnEdge" + AuxStateSuffix, NTracers,
                     Mesh->NEdgesSize, NVertLevels),
      Del2TracersOnCell("Del2TracerOnCell" + AuxStateSuffix, NTracers,
                        Mesh->NCellsSize, NVertLevels),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell) {}

void TracerAuxVars::registerFields(const std::string &AuxGroupName,
                                   const std::string &MeshName) const {

   int Err = 0; // error code

   // Create fields
   const Real FillValue = -9.99e30;
   int NDims            = 3;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Thickness-weighted tracers on Edge
   DimNames[0]            = "NTracers";
   DimNames[1]            = "NEdges" + DimSuffix;
   DimNames[2]            = "NVertLevels";
   auto HTracersEdgeField = Field::create(
       HTracersOnEdge.label(), // field name
       "thickness-weighted tracers at edges. May be centered, upwinded, or a "
       "combination of the two.",        // long name or description
       "",                               // units
       "",                               // CF standard name
       0,                                // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Del2 tracers on Cell
   DimNames[1]               = "NCells" + DimSuffix;
   auto Del2TracersCellField = Field::create(
       Del2TracersOnCell.label(),                                // field name
       "laplacian of thickness-weighted tracers at cell center", // long name or
                                                                 // description
       "",                                                       // units
       "",                               // CF standard name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Add fields to Aux Field group
   Err = FieldGroup::addFieldToGroup(HTracersOnEdge.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", HTracersOnEdge.label(),
                AuxGroupName);

   Err = FieldGroup::addFieldToGroup(Del2TracersOnCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", Del2TracersOnCell.label(),
                AuxGroupName);

   // Attach data to fields
   Err = HTracersEdgeField->attachData<Array3DReal>(HTracersOnEdge);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", HTracersOnEdge.label());

   Err = Del2TracersCellField->attachData<Array3DReal>(Del2TracersOnCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", Del2TracersOnCell.label());
}

void TracerAuxVars::unregisterFields() const {
   int Err = 0;

   Err = Field::destroy(HTracersOnEdge.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", HTracersOnEdge.label());

   Err = Field::destroy(Del2TracersOnCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", Del2TracersOnCell.label());
}

} // namespace OMEGA
