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
   Array2DReal SurfTracerRestValuesCell;
   Array2DReal TracersMonthlySurfClimoCell;
   Real MaxDiff = 100.0; // maximum allowed difference for restoring
   /// Need to add under sea ice restoring option when that is available

   SurfTracerRestAuxVars(const std::string &AuxStateSuffix,
                         const HorzMesh *Mesh, const VertCoord *VCoord,
                         const I4 NTracers);

   KOKKOS_FUNCTION void
   computeVarsOnCells(int L, int ICell, const Array3DReal &TracerCell) const {

      const I4 K = MinLayerCell(ICell);

      if ((TracersMonthlySurfClimoCell(L, ICell) - TracerCell(L, ICell, K)) >
          MaxDiff) {
         SurfTracerRestValuesCell(L, ICell) = MaxDiff;
      } else if ((TracersMonthlySurfClimoCell(L, ICell) -
                  TracerCell(L, ICell, K)) < -MaxDiff) {
         SurfTracerRestValuesCell(L, ICell) = -MaxDiff;
      } else {
         SurfTracerRestValuesCell(L, ICell) =
             TracersMonthlySurfClimoCell(L, ICell) - TracerCell(L, ICell, K);
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 MinLayerCell;
};

} // namespace OMEGA
#endif
