#include "KineticAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

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

void KineticAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   int Err = 0; // error flag for some calls

   // Create fields
   const Real FillValue = -9.99e30;
   int NDims            = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Kinetic energy on cells
   DimNames[0]                 = "NCells" + DimSuffix;
   DimNames[1]                 = "NVertLevels";
   auto KineticEnergyCellField = Field::create(
       KineticEnergyCell.label(),                        // field name
       "kinetic energy of horizontal velocity on cells", // long name/describe
       "m^2 s^-2",                                       // units
       "specific_kinetic_energy_of_sea_water",           // CF standard Name
       0,                                                // min valid value
       std::numeric_limits<Real>::max(),                 // max valid value
       FillValue, // scalar for undefined entries
       2,         // number of dimensions
       DimNames   // dim names
   );

   // Velocity divergence on cells
   auto VelocityDivCellField = Field::create(
       VelocityDivCell.label(),             // field name
       "divergence of horizontal velocity", // long Name or description
       "s^-1",                              // units
       "",                                  // CF standard Name
       std::numeric_limits<Real>::min(),    // min valid value
       std::numeric_limits<Real>::max(),    // max valid value
       FillValue,                           // scalar used for undefined entries
       NDims,                               // number of dimensions
       DimNames                             // dimension names
   );

   // Add fields to FieldGroup
   Err = FieldGroup::addFieldToGroup(KineticEnergyCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", KineticEnergyCell.label(),
                AuxGroupName);
   Err = FieldGroup::addFieldToGroup(VelocityDivCell.label(), AuxGroupName);
   if (Err != 0)
      LOG_ERROR("Error adding field {} to group {}", VelocityDivCell.label(),
                AuxGroupName);

   // Attach data
   Err = KineticEnergyCellField->attachData<Array2DReal>(KineticEnergyCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", KineticEnergyCell.label());

   Err = VelocityDivCellField->attachData<Array2DReal>(VelocityDivCell);
   if (Err != 0)
      LOG_ERROR("Error attaching data to field {}", VelocityDivCell.label());
}

void KineticAuxVars::unregisterFields() const {
   int Err = 0;
   Err     = Field::destroy(KineticEnergyCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", KineticEnergyCell.label());
   Err = Field::destroy(VelocityDivCell.label());
   if (Err != 0)
      LOG_ERROR("Error destroying field {}", VelocityDivCell.label());
}

} // namespace OMEGA
