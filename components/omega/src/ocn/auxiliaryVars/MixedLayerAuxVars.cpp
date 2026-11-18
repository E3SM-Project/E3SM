#include "MixedLayerAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

MixedLayerAuxVars::MixedLayerAuxVars(const std::string &AuxStateSuffix,
                                     const HorzMesh *Mesh,
                                     const VertCoord *VCoord)
    : DenMixLayerDepth("DenMixLayerDepth" + AuxStateSuffix, Mesh->NCellsSize),
      DenMixLayerIndex("DenMixLayerIndex" + AuxStateSuffix, Mesh->NCellsSize),
      ReferencePressure("ReferencePressure" + AuxStateSuffix, Mesh->NCellsSize,
                        VCoord->NVertLayers),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      GeomZInterface(VCoord->GeomZInterface), GeomZMid(VCoord->GeomZMid) {
   deepCopy(ReferencePressure, ReferencePressureVal);
}

void MixedLayerAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   // Create and add mixed layer depth field
   {
      int NDims = 1;
      std::vector<std::string> DimNames(NDims);
      DimNames[0] = "NCells" + DimSuffix;

      auto DenMixLayerDepthField =
          Field::create(DenMixLayerDepth.label(),         // Field name
                        "Mixed Layer Depth",              // Long Name
                        "m",                              // Units
                        "",                               // CF-ish Name
                        0.0,                              // Min valid value
                        std::numeric_limits<Real>::max(), // Max valid value
                        NDims,   // Number of dimensions
                        DimNames // Dimension names
          );

      FieldGroup::addFieldToGroup(DenMixLayerDepth.label(), AuxGroupName);
      DenMixLayerDepthField->attachData<Array1DReal>(DenMixLayerDepth);
   }
}

void MixedLayerAuxVars::unregisterFields() const {
   Field::destroy(DenMixLayerDepth.label());
}

} // namespace OMEGA
