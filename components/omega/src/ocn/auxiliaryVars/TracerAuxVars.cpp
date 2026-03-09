#include "TracerAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

TracerAuxVars::TracerAuxVars(const std::string &AuxStateSuffix,
                             const HorzMesh *Mesh, const VertCoord *VCoord,
                             const I4 NTracers)
    : Del2TracersCell("Del2TracerCell" + AuxStateSuffix, NTracers,
                      Mesh->NCellsSize, VCoord->NVertLayers),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

void TracerAuxVars::registerFields(const std::string &AuxGroupName,
                                   const std::string &MeshName) const {

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

   // Del2 tracers on Cell
   DimNames[0]               = "NTracers";
   DimNames[1]               = "NCells" + DimSuffix;
   DimNames[2]               = "NVertLayers";
   auto Del2TracersCellField = Field::create(
       Del2TracersCell.label(),                                  // field name
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
   FieldGroup::addFieldToGroup(Del2TracersCell.label(), AuxGroupName);

   // Attach data to fields
   Del2TracersCellField->attachData<Array3DReal>(Del2TracersCell);
}

void TracerAuxVars::unregisterFields() const {
   Field::destroy(Del2TracersCell.label());
}

} // namespace OMEGA
