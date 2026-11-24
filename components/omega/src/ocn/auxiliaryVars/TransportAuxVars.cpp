#include "TransportAuxVars.h"
#include "DataTypes.h"
#include "Field.h"

#include <limits>

namespace OMEGA {

TransportAuxVars::TransportAuxVars(const std::string &AuxStateSuffix,
                                   const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : NormalTransportVelocity("NormalTransportVelocity" + AuxStateSuffix,
                              Mesh->NEdgesSize, VCoord->NVertLayers) {}

void TransportAuxVars::registerFields(
    const std::string &AuxGroupName, // name of Auxiliary field group
    const std::string &MeshName      // name of horizontal mesh
) const {

   // Create/define fields
   int NDims = 2;
   std::vector<std::string> DimNames(NDims);
   std::string DimSuffix;
   if (MeshName == "Default") {
      DimSuffix = "";
   } else {
      DimSuffix = MeshName;
   }

   DimNames[0] = "NEdges" + DimSuffix;
   DimNames[1] = "NVertLayers";

   auto NormalTransportVelocityField =
       Field::create(NormalTransportVelocity.label(), // field name
                     "horizontal velocity used to transport pseudo-thickness "
                     "and tracers", // long Name or description
                     "m s-1",       // units
                     "",            // CF standard Name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     NDims,   // number of dimensions
                     DimNames // dimension names
       );

   // Add fields to FieldGroup
   FieldGroup::addFieldToGroup(NormalTransportVelocity.label(), AuxGroupName);

   // Attach data
   NormalTransportVelocityField->attachData<Array2DReal>(
       NormalTransportVelocity);
}

void TransportAuxVars::unregisterFields() const {
   Field::destroy(NormalTransportVelocity.label());
}

} // namespace OMEGA
