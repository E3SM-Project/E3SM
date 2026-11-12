#ifndef OMEGA_AUX_MIX_LAYER_H
#define OMEGA_AUX_MIX_LAYER_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

// Generic linear interpolation routine. This should be in OmegaMath.h or
// something like that once it exists.
KOKKOS_INLINE_FUNCTION Real linearInterp(Real x, Real y1, Real x1, Real y2,
                                         Real x2) {
   if (x1 == x2)
      return y1;
   const Real A = (y1 - y2) / (x1 - x2);
   const Real B = y1 - A * x1;
   return A * x + B;
}

class MixedLayerAuxVars {
 public:
   // Reference depth
   Real ReferenceDepth = 10;

   // Reference pressure value
   Real ReferencePressureVal = 1e5;

   // Density threshold for determining the mixed layer depth
   Real DenThreshold = 0.03;

   // Mixed layer index
   Array1DI4 DenMixLayerIndex;

   // Mixed layer depth
   Array1DReal DenMixLayerDepth;

   // Reference pressure array for computation of displaced spec volume
   Array2DReal ReferencePressure;

   MixedLayerAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                     const VertCoord *VCoord);

   // Compute mixed layer index and  depth based on the density difference
   // criterion
   KOKKOS_FUNCTION void
   computeVarsOnCell(const TeamMember &Team, int ICell,
                     const Array2DReal &SpecVolDisp) const {

      const Real SSH = GeomZInterface(ICell, MinLayerCell(ICell));

      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);

      // Find first interface where depth >= reference depth
      int KRef;
      parallelSearchInner(
          Team, Range{KMin + 1, KMax},
          INNER_LAMBDA(int K) {
             const Real Depth = SSH - GeomZInterface(ICell, K);
             return Depth >= ReferenceDepth;
          },
          KRef);

      // Not found, set mixed layer depth to to the depth of the deepest layer
      // and return early
      if (KRef == -1) {
         DenMixLayerIndex(ICell) = KMax;
         DenMixLayerDepth(ICell) = SSH - GeomZMid(ICell, KMax);
         return;
      }

      const int KRefM1 = Kokkos::max(KRef - 1, MinLayerCell(ICell));

      const Real DepthKRef   = SSH - GeomZMid(ICell, KRef);
      const Real DepthKRefM1 = SSH - GeomZMid(ICell, KRefM1);

      const Real ReferenceSpecVol =
          linearInterp(ReferenceDepth, SpecVolDisp(ICell, KRef), DepthKRef,
                       SpecVolDisp(ICell, KRefM1), DepthKRefM1);

      // Start searching from reference level - 1
      int KDen;
      parallelSearchInner(
          Team, Range{KRefM1, KMax},
          INNER_LAMBDA(int K) {
             return (ReferenceSpecVol / SpecVolDisp(ICell, K) - 1) >=
                    DenThreshold * ReferenceSpecVol;
          },
          KDen);

      // Not found. Setting to the depth of the deepest layer
      if (KDen == -1) {
         DenMixLayerIndex(ICell) = KMax;
         DenMixLayerDepth(ICell) = SSH - GeomZMid(ICell, KMax);
      } else { // Found
         const int KDenM1 = Kokkos::max(KDen - 1, MinLayerCell(ICell));

         const Real DepthKDen   = SSH - GeomZMid(ICell, KDen);
         const Real DepthKDenM1 = SSH - GeomZMid(ICell, KDenM1);

         const Real FactorKDen =
             ReferenceSpecVol / SpecVolDisp(ICell, KDen) - 1;
         const Real FactorKDenM1 =
             ReferenceSpecVol / SpecVolDisp(ICell, KDenM1) - 1;

         Real MixedLayerDepth =
             linearInterp(DenThreshold * ReferenceSpecVol, DepthKDen,
                          FactorKDen, DepthKDenM1, FactorKDenM1);

         // guarantee MLD between DepthKDenM1 and DepthKDen
         // this can happen because density difference in the first layer
         // can already be above the threshold
         MixedLayerDepth =
             Kokkos::clamp(MixedLayerDepth, DepthKDenM1, DepthKDen);

         DenMixLayerIndex(ICell) = KDen;
         DenMixLayerDepth(ICell) = MixedLayerDepth;
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal GeomZInterface;
   Array2DReal GeomZMid;
};

} // namespace OMEGA
#endif
