#include "LayerThicknessAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

LayerThicknessAuxVars::LayerThicknessAuxVars(const std::string &AuxStateSuffix,
                                             const HorzMesh *Mesh,
                                             int NVertLevels)
    : FluxLayerThickEdge("FluxLayerThickEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, NVertLevels),
      MeanLayerThickEdge("MeanLayerThickEdge" + AuxStateSuffix,
                         Mesh->NEdgesSize, NVertLevels),
      SshCell("SshCell" + AuxStateSuffix, Mesh->NCellsSize, NVertLevels),
      CellsOnEdge(Mesh->CellsOnEdge), BottomDepth(Mesh->BottomDepth) {}

void LayerThicknessAuxVars::registerFields(
    const std::string &AuxGroupName) const {
   addMetaData(AuxGroupName);
   defineIOFields();
}

void LayerThicknessAuxVars::unregisterFields() const {
   IOField::erase(FluxLayerThickEdge.label());
   IOField::erase(MeanLayerThickEdge.label());
   IOField::erase(SshCell.label());
   MetaData::destroy(FluxLayerThickEdge.label());
   MetaData::destroy(MeanLayerThickEdge.label());
   MetaData::destroy(SshCell.label());
}

void LayerThicknessAuxVars::addMetaData(const std::string &AuxGroupName) const {
   auto EdgeDim      = MetaDim::get("NEdges");
   auto CellDim      = MetaDim::get("NCells");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get(AuxGroupName);

   const Real FillValue = -9.99e30;

   // Flux layer thickness on edges
   auto FluxLayerThickEdgeMeta = ArrayMetaData::create(
       FluxLayerThickEdge.label(),
       "layer thickness used for fluxes through edges. May be centered, "
       "upwinded, or a combination of the two.", /// long Name or description
       "m",                                      /// units
       "",                                       /// CF standard Name
       0,                                        /// min valid value
       std::numeric_limits<Real>::max(),         /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {EdgeDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(FluxLayerThickEdge.label());

   // Mean layer thickness on edges
   auto MeanLayerThickEdgeMeta = ArrayMetaData::create(
       MeanLayerThickEdge.label(),
       "layer thickness averaged from cell center to edges", /// long Name or
                                                             /// description
       "m",                                                  /// units
       "",                               /// CF standard Name
       0,                                /// min valid value
       std::numeric_limits<Real>::max(), /// max valid value
       FillValue,                        /// scalar used for undefined entries
       2,                                /// number of dimensions
       {EdgeDim, VertDim}                /// dim pointers
   );
   AuxMetaGroup->addField(MeanLayerThickEdge.label());

   // Sea surface height
   auto SshCellMeta = ArrayMetaData::create(
       SshCell.label(),
       "sea surface height at cell center", /// long Name or description
       "m",                                 /// units
       "",                                  /// CF standard Name
       0,                                   /// min valid value
       std::numeric_limits<Real>::max(),    /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {CellDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(SshCell.label());
}

void LayerThicknessAuxVars::defineIOFields() const {
   int Err;

   // Flux layer thickness on edges
   Err = IOField::define(FluxLayerThickEdge.label());
   Err = IOField::attachData(FluxLayerThickEdge.label(), FluxLayerThickEdge);

   // Mean layer thickness on edges
   Err = IOField::define(MeanLayerThickEdge.label());
   Err = IOField::attachData(MeanLayerThickEdge.label(), MeanLayerThickEdge);

   // Sea surface height
   Err = IOField::define(SshCell.label());
   Err = IOField::attachData(SshCell.label(), SshCell);
}

} // namespace OMEGA
