#include "LayerThicknessAuxVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

LayerThicknessAuxVars::LayerThicknessAuxVars(const std::string &AuxStateSuffix,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : FluxLayerThickEdge("FluxLayerThickEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, VCoord->NVertLayers),
      MeanLayerThickEdge("MeanLayerThickEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, VCoord->NVertLayers),
      SshCell("SshCell" + AuxStateSuffix, Mesh->NCellsSize,
              VCoord->NVertLayers),
      ProvThickness("ProvThickness" + AuxStateSuffix, Mesh->NCellsSize,
                    VCoord->NVertLayers),
      AreaCell(Mesh->AreaCell), DvEdge(Mesh->DvEdge),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), CellsOnEdge(Mesh->CellsOnEdge),
      BottomDepth(VCoord->BottomDepth),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

void LayerThicknessAuxVars::registerFields(const std::string &AuxGroupName,
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

   // Flux layer thickness on edges
   auto FluxLayerThickEdgeField = Field::create(
       FluxLayerThickEdge.label(), // field name
       "layer thickness used for fluxes through edges. May be centered, "
       "upwinded, or a combination of the two.", // long Name or description
       "m",                                      // units
       "",                                       // CF standard Name
       0,                                        // min valid value
       std::numeric_limits<Real>::max(),         // max valid value
       FillValue,                                // scalar for undefined entries
       NDims,                                    // number of dimensions
       DimNames                                  // dimension names
   );

   // Mean layer thickness on edges
   auto MeanLayerThickEdgeField = Field::create(
       MeanLayerThickEdge.label(),                           // field name
       "layer thickness averaged from cell center to edges", // long Name or
                                                             // description
       "m",                                                  // units
       "",                                                   // CF standard Name
       0,                                                    // min valid value
       std::numeric_limits<Real>::max(),                     // max valid value
       FillValue, // scalar used for undefined entries
       NDims,     // number of dimensions
       DimNames   // dimension names
   );

   // Sea surface height
   DimNames[0]       = "NCells" + DimSuffix;
   auto SshCellField = Field::create(
       SshCell.label(),                     // field name
       "sea surface height at cell center", // long Name or description
       "m",                                 // units
       "sea_surface_height",                // CF standard Name
       0,                                   // min valid value
       std::numeric_limits<Real>::max(),    // max valid value
       FillValue,                           // scalar for undefined entries
       NDims,                               // number of dimensions
       DimNames                             // dimension names
   );

   // Provisional Thickness
   auto ProvThicknessField = Field::create(
       ProvThickness.label(),                   // field name
       "layer thickness after horizontal flux", // long Name or description
       "m",                                     // units
       "",                                      // CF standard Name
       0,                                       // min valid value
       std::numeric_limits<Real>::max(),        // max valid value
       FillValue,                               // scalar for undefined entries
       NDims,                                   // number of dimensions
       DimNames                                 // dimension names
   );

   // Add fields to Aux field group
   FieldGroup::addFieldToGroup(FluxLayerThickEdge.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(MeanLayerThickEdge.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(SshCell.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(ProvThickness.label(), AuxGroupName);

   // Attach field data
   FluxLayerThickEdgeField->attachData<Array2DReal>(FluxLayerThickEdge);
   MeanLayerThickEdgeField->attachData<Array2DReal>(MeanLayerThickEdge);
   SshCellField->attachData<Array2DReal>(SshCell);
   ProvThicknessField->attachData<Array2DReal>(ProvThickness);
}

void LayerThicknessAuxVars::unregisterFields() const {
   Field::destroy(FluxLayerThickEdge.label());
   Field::destroy(MeanLayerThickEdge.label());
   Field::destroy(SshCell.label());
   Field::destroy(ProvThickness.label());
}

} // namespace OMEGA
