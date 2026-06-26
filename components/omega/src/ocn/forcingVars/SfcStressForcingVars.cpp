#include "SfcStressForcingVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

SfcStressForcingVars::SfcStressForcingVars(const std::string &Suffix,
                                           const HorzMesh *Mesh)
    : NormalStressEdge("NormalStressEdge" + Suffix, Mesh->NEdgesSize),
      ZonalStressCell("SfcStressZonal" + Suffix, Mesh->NCellsSize),
      MeridStressCell("SfcStressMeridional" + Suffix, Mesh->NCellsSize),
      CellsOnEdge(Mesh->CellsOnEdge), AngleEdge(Mesh->AngleEdge), Interp(Mesh) {
}

void SfcStressForcingVars::registerFields(
    const std::string &MeshName // name of horizontal mesh
) const {

   const Real FillValue = -9.99e30;
   int NDims            = 1;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NCells" + DimSuffix;
   auto ZonalStressCellField =
       Field::create(ZonalStressCell.label(),          // field name
                     "zonal surface stress",           // long name/describe
                     "N m^{-2}",                       // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     NDims,                            // number of dimensions
                     DimNames                          // dim names
       );

   auto MeridStressCellField =
       Field::create(MeridStressCell.label(),     // field name
                     "meridional surface stress", // long Name or description
                     "N m^{-2}",                  // units
                     "",                          // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     NDims,                            // number of dimensions
                     DimNames                          // dimension names
       );

   FieldGroup::addFieldToGroup(ZonalStressCell.label(), "Forcing");
   FieldGroup::addFieldToGroup(MeridStressCell.label(), "Forcing");

   ZonalStressCellField->attachData<Array1DReal>(ZonalStressCell);
   MeridStressCellField->attachData<Array1DReal>(MeridStressCell);
}

void SfcStressForcingVars::unregisterFields() const {
   if (Field::exists(ZonalStressCell.label())) {
      Field::destroy(ZonalStressCell.label());
   }
   if (Field::exists(MeridStressCell.label())) {
      Field::destroy(MeridStressCell.label());
   }
}

} // namespace OMEGA
