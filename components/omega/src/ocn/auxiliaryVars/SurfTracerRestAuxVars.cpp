#include "SurfTracerRestAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

SurfTracerRestAuxVars::SurfTracerRestAuxVars(const std::string &AuxStateSuffix,
                                             const HorzMesh *Mesh,
                                             const I4 NTracers)
    : TracersMonthlySurfClimoCell("TracersMonthlySurfClimoCell" +
                                      AuxStateSuffix,
                                  NTracers, Mesh->NCellsSize) {}

void SurfTracerRestAuxVars::registerFields(
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

   auto TracersMonthlySurfClimoCellField =
       Field::create(TracersMonthlySurfClimoCell.label(),
                     "monthly surface tracer climatology", // long name/describe
                     "tracer units",                       // units
                     "",                                   // CF standard Name
                     std::numeric_limits<Real>::min(),     // min valid value
                     std::numeric_limits<Real>::max(),     // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames); // dim names

   // Add fields to FieldGroup
   FieldGroup::addFieldToGroup(TracersMonthlySurfClimoCell.label(),
                               AuxGroupName);

   // Attach data
   TracersMonthlySurfClimoCellField->attachData<Array2DReal>(
       TracersMonthlySurfClimoCell);
}

void SurfTracerRestAuxVars::unregisterFields() const {
   Field::destroy(TracersMonthlySurfClimoCell.label());
}

} // namespace OMEGA
