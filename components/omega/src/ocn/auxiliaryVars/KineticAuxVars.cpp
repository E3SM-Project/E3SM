#include "KineticAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

const std::string KineticAuxVars::KineticEnergyCellName = "KineticEnergyCell";
const std::string KineticAuxVars::VelocityDivCellName   = "VelocityDivCell";

KineticAuxVars::KineticAuxVars(const HorzMesh *mesh, int NVertLevels)
    : KineticEnergyCell("KineticEnergyCell", mesh->NCellsSize, NVertLevels),
      VelocityDivCell("VelocityDivCell", mesh->NCellsSize, NVertLevels),
      NEdgesOnCell(mesh->NEdgesOnCell), EdgesOnCell(mesh->EdgesOnCell),
      EdgeSignOnCell(mesh->EdgeSignOnCell), DcEdge(mesh->DcEdge),
      DvEdge(mesh->DvEdge), AreaCell(mesh->AreaCell) {
   addMetaData();
   defineIOFields();
}

void KineticAuxVars::addMetaData() const {
   auto CellDim      = MetaDim::get("NCells");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get("auxiliaryVars");

   const Real FillValue = -9.99e30;

   // Kinetic energy on cells
   auto KineticEnergyCellMeta = ArrayMetaData::create(
       KineticEnergyCellName,
       "kinetic energy of horizontal velocity on cells", /// long Name or
                                                         /// description
       "m^2 s^-2",                                       /// units
       "",                                               /// CF standard Name
       0,                                                /// min valid value
       std::numeric_limits<Real>::max(),                 /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {CellDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(KineticEnergyCellName);

   // Velocity divergence on cells
   auto VelocityDivCellMeta = ArrayMetaData::create(
       VelocityDivCellName,
       "divergence of horizontal velocity", /// long Name or description
       "s^-1",                              /// units
       "",                                  /// CF standard Name
       std::numeric_limits<Real>::min(),    /// min valid value
       std::numeric_limits<Real>::max(),    /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {CellDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(VelocityDivCellName);
}

void KineticAuxVars::defineIOFields() const {
   int Err;

   // Kinetic energy on cells
   Err = IOField::define(KineticEnergyCellName);
   Err = IOField::attachData(KineticEnergyCellName, KineticEnergyCell);

   // Velocity divergence on cells
   Err = IOField::define(VelocityDivCellName);
   Err = IOField::attachData(VelocityDivCellName, VelocityDivCell);
}

} // namespace OMEGA
