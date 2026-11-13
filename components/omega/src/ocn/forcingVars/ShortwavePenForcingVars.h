#ifndef OMEGA_FORCING_SHORTWAVE_PEN_H
#define OMEGA_FORCING_SHORTWAVE_PEN_H

#include "DataTypes.h"
#include "HorzMesh.h"

#include <string>

namespace OMEGA {

class ShortwavePenForcingVars {
 public:
   Array1DReal ExtinctionCoeffRedCell;
   Array1DReal ExtinctionCoeffBlueCell;

   ShortwavePenForcingVars(const std::string &Suffix, const HorzMesh *Mesh);

   void registerFields(const std::string &GroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;
};

} // namespace OMEGA

#endif
