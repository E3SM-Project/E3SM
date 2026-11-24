#ifndef OMEGA_AUX_TRANSPORT_H
#define OMEGA_AUX_TRANSPORT_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class TransportAuxVars {
 public:
   Array2DReal NormalTransportVelocity;

   TransportAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                    const VertCoord *VCoord);

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;
};

} // namespace OMEGA
#endif
