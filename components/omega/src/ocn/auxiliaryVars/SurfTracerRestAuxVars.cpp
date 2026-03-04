#include "SurfTracerRestAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

SurfTracerRestAuxVars::SurfTracerRestAuxVars(const std::string &AuxStateSuffix,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             const I4 NTracers)
    : SurfTracerRestValuesCell("SurfTracerRestValuesCell" + AuxStateSuffix,
                               NTracers, Mesh->NCellsSize),
      TracersMonthlySurfClimoCell("TracersMonthlySurfClimoCell" +
                                      AuxStateSuffix,
                                  NTracers, Mesh->NCellsSize),
      MinLayerCell(VCoord->MinLayerCell) {}

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

   // Tracers surface restoring values
   DimNames[0] = "NTracers";
   DimNames[1] = "NCells" + DimSuffix;
   auto SurfTracerRestValuesCellField =
       Field::create(SurfTracerRestValuesCell.label(), // field name
                     "tracer surface restoring value", // long name/describe
                     "tracer units",                   // units
                     "",                               // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValue, // scalar for undefined entries
                     NDims,     // number of dimensions
                     DimNames   // dim names
       );

   // Add fields to FieldGroup
   FieldGroup::addFieldToGroup(SurfTracerRestValuesCell.label(), AuxGroupName);

   // Attach data
   SurfTracerRestValuesCellField->attachData<Array2DReal>(
       SurfTracerRestValuesCell);
}

void SurfTracerRestAuxVars::unregisterFields() const {
   Field::destroy(SurfTracerRestValuesCell.label());
}

} // namespace OMEGA
