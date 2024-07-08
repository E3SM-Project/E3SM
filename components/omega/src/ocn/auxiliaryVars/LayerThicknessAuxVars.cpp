#include "LayerThicknessAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

const std::string LayerThicknessAuxVars::FluxLayerThickEdgeName =
    "FluxLayerThickEdge";
const std::string LayerThicknessAuxVars::MeanLayerThickEdgeName =
    "MeanLayerThickEdge";

LayerThicknessAuxVars::LayerThicknessAuxVars(const HorzMesh *Mesh,
                                             int NVertLevels)
    : FluxLayerThickEdge("FluxLayerThickEdge", Mesh->NEdgesSize, NVertLevels),
      MeanLayerThickEdge("MeanLayerThickEdge", Mesh->NEdgesSize, NVertLevels),
      CellsOnEdge(Mesh->CellsOnEdge) {
   addMetaData();
   defineIOFields();
}

void LayerThicknessAuxVars::addMetaData() const {
   auto EdgeDim      = MetaDim::get("NEdges");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get("auxiliaryVars");

   const Real FillValue = -9.99e30;

   // Flux layer thickness on edges
   auto FluxLayerThickEdgeMeta = ArrayMetaData::create(
       FluxLayerThickEdgeName,
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
   AuxMetaGroup->addField(FluxLayerThickEdgeName);

   // Mean layer thickness on edges
   auto MeanLayerThickEdgeMeta = ArrayMetaData::create(
       MeanLayerThickEdgeName,
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
   AuxMetaGroup->addField(MeanLayerThickEdgeName);
}

void LayerThicknessAuxVars::defineIOFields() const {
   int Err;

   // Flux layer thickness on edges
   Err = IOField::define(FluxLayerThickEdgeName);
   Err = IOField::attachData(FluxLayerThickEdgeName, FluxLayerThickEdge);

   // Mean layer thickness on edges
   Err = IOField::define(MeanLayerThickEdgeName);
   Err = IOField::attachData(MeanLayerThickEdgeName, MeanLayerThickEdge);
}

} // namespace OMEGA
