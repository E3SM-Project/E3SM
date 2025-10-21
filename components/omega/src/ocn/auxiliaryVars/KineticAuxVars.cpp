#include "KineticAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

KineticAuxVars::KineticAuxVars(const std::string &AuxStateSuffix,
                               const HorzMesh *Mesh, const VertCoord *VCoord)
    : KineticEnergyCell("KineticEnergyCell" + AuxStateSuffix, Mesh->NCellsSize,
                        VCoord->NVertLayers),
      VelocityDivCell("VelocityDivCell" + AuxStateSuffix, Mesh->NCellsSize,
                      VCoord->NVertLayers),
      NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), DcEdge(Mesh->DcEdge),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell) {}

void KineticAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

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
   DimNames[1]                 = "NVertLayers";
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
   FieldGroup::addFieldToGroup(KineticEnergyCell.label(), AuxGroupName);
   FieldGroup::addFieldToGroup(VelocityDivCell.label(), AuxGroupName);

   // Attach data
   KineticEnergyCellField->attachData<Array2DReal>(KineticEnergyCell);
   VelocityDivCellField->attachData<Array2DReal>(VelocityDivCell);
}

void KineticAuxVars::unregisterFields() const {
   Field::destroy(KineticEnergyCell.label());
   Field::destroy(VelocityDivCell.label());
}

} // namespace OMEGA
