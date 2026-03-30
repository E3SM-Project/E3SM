#ifndef OMEGA_AUX_SURF_REST_H
#define OMEGA_AUX_SURF_REST_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class SurfTracerRestAuxVars {
 public:
   Array2DReal TracersMonthlySurfClimoCell;

   SurfTracerRestAuxVars(const std::string &AuxStateSuffix,
                         const HorzMesh *Mesh, const I4 NTracers);

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;
};

} // namespace OMEGA
#endif
