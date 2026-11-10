#ifndef OMEGA_STRESS_H
#define OMEGA_STRESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class SfcStressForcingVars {
 public:
   Array1DReal NormalStressEdge;
   Array1DReal ZonalStressCell;
   Array1DReal MeridStressCell;
   InterpCellToEdgeOption InterpChoice;

   SfcStressForcingVars(const std::string &Suffix, const HorzMesh *Mesh);

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge) const {
      const Real ZonalStressEdge = Interp(IEdge, ZonalStressCell, InterpChoice);
      const Real MeridStressEdge = Interp(IEdge, MeridStressCell, InterpChoice);

      NormalStressEdge(IEdge) =
          Kokkos::cos(AngleEdge(IEdge)) * ZonalStressEdge +
          Kokkos::sin(AngleEdge(IEdge)) * MeridStressEdge;
   }

   void registerFields(const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   InterpCellToEdge Interp;
   Array2DI4 CellsOnEdge;
   Array1DReal AngleEdge;
};

} // namespace OMEGA
#endif
