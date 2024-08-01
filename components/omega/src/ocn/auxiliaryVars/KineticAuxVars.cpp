#include "KineticAuxVars.h"
#include "IOField.h"
#include "MetaData.h"

#include <limits>

namespace OMEGA {

KineticAuxVars::KineticAuxVars(const std::string &AuxStateSuffix,
                               const HorzMesh *Mesh, int NVertLevels)
    : KineticEnergyCell("KineticEnergyCell" + AuxStateSuffix, Mesh->NCellsSize,
                        NVertLevels),
      VelocityDivCell("VelocityDivCell" + AuxStateSuffix, Mesh->NCellsSize,
                      NVertLevels),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell) {}

void KineticAuxVars::registerFields(const std::string &AuxGroupName) const {
   addMetaData(AuxGroupName);
   defineIOFields();
}

void KineticAuxVars::unregisterFields() const {
   IOField::erase(KineticEnergyCell.label());
   IOField::erase(VelocityDivCell.label());
   MetaData::destroy(KineticEnergyCell.label());
   MetaData::destroy(VelocityDivCell.label());
}

void KineticAuxVars::addMetaData(const std::string &AuxGroupName) const {
   auto CellDim      = MetaDim::get("NCells");
   auto VertDim      = MetaDim::get("NVertLevels");
   auto AuxMetaGroup = MetaGroup::get(AuxGroupName);

   const Real FillValue = -9.99e30;

   // Kinetic energy on cells
   auto KineticEnergyCellMeta = ArrayMetaData::create(
       KineticEnergyCell.label(),
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
   AuxMetaGroup->addField(KineticEnergyCell.label());

   // Velocity divergence on cells
   auto VelocityDivCellMeta = ArrayMetaData::create(
       VelocityDivCell.label(),
       "divergence of horizontal velocity", /// long Name or description
       "s^-1",                              /// units
       "",                                  /// CF standard Name
       std::numeric_limits<Real>::min(),    /// min valid value
       std::numeric_limits<Real>::max(),    /// max valid value
       FillValue,         /// scalar used for undefined entries
       2,                 /// number of dimensions
       {CellDim, VertDim} /// dim pointers
   );
   AuxMetaGroup->addField(VelocityDivCell.label());
}

void KineticAuxVars::defineIOFields() const {
   int Err;

   // Kinetic energy on cells
   Err = IOField::define(KineticEnergyCell.label());
   Err = IOField::attachData(KineticEnergyCell.label(), KineticEnergyCell);

   // Velocity divergence on cells
   Err = IOField::define(VelocityDivCell.label());
   Err = IOField::attachData(VelocityDivCell.label(), VelocityDivCell);
}

} // namespace OMEGA
