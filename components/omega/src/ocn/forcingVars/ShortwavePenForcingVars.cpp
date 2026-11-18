#include "ShortwavePenForcingVars.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

ShortwavePenForcingVars::ShortwavePenForcingVars(const std::string &Suffix,
                                                 const HorzMesh *Mesh)
    : ExtinctionCoeffRedCell("ExtinctionCoeffRedCell" + Suffix,
                             Mesh->NCellsSize),
      ExtinctionCoeffBlueCell("ExtinctionCoeffBlueCell" + Suffix,
                              Mesh->NCellsSize) {}

void ShortwavePenForcingVars::registerFields(
    const std::string &GroupName, const std::string &MeshName) const {
   const int NDims = 1;
   std::string DimSuffix;
   if (MeshName != "Default") {
      DimSuffix = MeshName;
   }
   const std::vector<std::string> DimNames = {"NCells" + DimSuffix};

   auto ExtinctionCoeffRedCellField = Field::create(
       ExtinctionCoeffRedCell.label(), "red-band extinction coefficient",
       "m-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);
   auto ExtinctionCoeffBlueCellField = Field::create(
       ExtinctionCoeffBlueCell.label(), "blue-band extinction coefficient",
       "m-1", "", std::numeric_limits<Real>::lowest(),
       std::numeric_limits<Real>::max(), NDims, DimNames);

   FieldGroup::addFieldToGroup(ExtinctionCoeffRedCell.label(), GroupName);
   FieldGroup::addFieldToGroup(ExtinctionCoeffBlueCell.label(), GroupName);

   ExtinctionCoeffRedCellField->attachData<Array1DReal>(ExtinctionCoeffRedCell);
   ExtinctionCoeffBlueCellField->attachData<Array1DReal>(
       ExtinctionCoeffBlueCell);
}

void ShortwavePenForcingVars::unregisterFields() const {
   Field::destroy(ExtinctionCoeffRedCell.label());
   Field::destroy(ExtinctionCoeffBlueCell.label());
}

} // namespace OMEGA
