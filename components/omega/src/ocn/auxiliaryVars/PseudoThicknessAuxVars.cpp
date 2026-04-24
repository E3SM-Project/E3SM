#include "PseudoThicknessAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

PseudoThicknessAuxVars::PseudoThicknessAuxVars(
    const std::string &AuxStateSuffix, const HorzMesh *Mesh,
    const VertCoord *VCoord)
    : FluxPseudoThickEdge("FluxPseudoThickEdge" + AuxStateSuffix,
                          Mesh->NEdgesSize, VCoord->NVertLayers),
      MeanPseudoThickEdge("MeanPseudoThickEdge" + AuxStateSuffix,
                          Mesh->NEdgesSize, VCoord->NVertLayers),
      ProvPseudoThickness("ProvPseudoThickness" + AuxStateSuffix,
                          Mesh->NCellsSize, VCoord->NVertLayers),
      AreaCell(Mesh->AreaCell), DvEdge(Mesh->DvEdge),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), CellsOnEdge(Mesh->CellsOnEdge),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

void PseudoThicknessAuxVars::registerFields(const std::string &AuxGroupName,
                                            const std::string &MeshName) const {

   // Create/define fields
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NEdges" + DimSuffix;
   DimNames[1] = "NVertLayers";

   // Flux pseudo-thickness on edges
   auto FluxPseudoThickEdgeField = Field::create(
       FluxPseudoThickEdge.label(), // field name
       "pseudo-thickness used for fluxes through edges. May be centered, "
       "upwinded, or a combination of the two.", // long Name or description
       "m",                                      // units
       "",                                       // CF standard Name
       0,                                        // min valid value
       std::numeric_limits<Real>::max(),         // max valid value
       FillValue,                                // scalar for undefined entries
       NDims,                                    // number of dimensions
       DimNames                                  // dimension names
   );

   // Mean pseudo-thickness on edges
   auto MeanPseudoThickEdgeField = Field::create(
       MeanPseudoThickEdge.label(),                           // field name
       "pseudo-thickness averaged from cell center to edges", // long Name or
                                                              // description
       "m",                                                   // units
       "",                               // CF standard Name
       0,                                // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValue,                        // scalar used for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   // Provisional Thickness
   auto ProvPseudoThicknessField = Field::create(
       ProvPseudoThickness.label(),              // field name
       "pseudo-thickness after horizontal flux", // long Name or description
       "m",                                      // units
       "",                                       // CF standard Name
       0,                                        // min valid value
       std::numeric_limits<Real>::max(),         // max valid value
       FillValue,                                // scalar for undefined entries
       NDims,                                    // number of dimensions
       DimNames                                  // dimension names
   );

   // Add fields to Aux field group
   FieldGroup::addFieldToGroup(FluxPseudoThickEdge.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(MeanPseudoThickEdge.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(ProvPseudoThickness.label(), AuxGroupName);

   // Attach field data
   FluxPseudoThickEdgeField->attachData<Array2DReal>(FluxPseudoThickEdge);
   MeanPseudoThickEdgeField->attachData<Array2DReal>(MeanPseudoThickEdge);
   ProvPseudoThicknessField->attachData<Array2DReal>(ProvPseudoThickness);
}

void PseudoThicknessAuxVars::unregisterFields() const {
   Field::destroy(FluxPseudoThickEdge.label());
   Field::destroy(MeanPseudoThickEdge.label());
   Field::destroy(ProvPseudoThickness.label());
}

} // namespace OMEGA
