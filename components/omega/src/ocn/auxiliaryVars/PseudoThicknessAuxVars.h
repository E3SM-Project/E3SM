#ifndef OMEGA_AUX_THICKNESS_H
#define OMEGA_AUX_THICKNESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

enum class FluxThickEdgeOption { Center, Upwind };

class PseudoThicknessAuxVars {
 public:
   Array2DReal FluxPseudoThickEdge;
   Array2DReal MeanPseudoThickEdge;
   Array2DReal ProvPseudoThickness;

   FluxThickEdgeOption FluxThickEdgeChoice;

   PseudoThicknessAuxVars(const std::string &AuxStateSuffix,
                          const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &PseudoThickCell,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K = KStart + KVec;
         MeanPseudoThickEdge(IEdge, K) =
             0.5_Real *
             (PseudoThickCell(JCell0, K) + PseudoThickCell(JCell1, K));
      }

      switch (FluxThickEdgeChoice) {
      case FluxThickEdgeOption::Center:
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            FluxPseudoThickEdge(IEdge, K) =
                0.5_Real *
                (PseudoThickCell(JCell0, K) + PseudoThickCell(JCell1, K));
         }
         break;
      case FluxThickEdgeOption::Upwind:
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            if (NormalVelEdge(IEdge, K) > 0) {
               FluxPseudoThickEdge(IEdge, K) = PseudoThickCell(JCell0, K);
            } else if (NormalVelEdge(IEdge, K) < 0) {
               FluxPseudoThickEdge(IEdge, K) = PseudoThickCell(JCell1, K);
            } else {
               FluxPseudoThickEdge(IEdge, K) = Kokkos::max(
                   PseudoThickCell(JCell0, K), PseudoThickCell(JCell1, K));
            }
         }
         break;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnCells(int ICell, int KChunk,
                                           const Array2DReal &PseudoThickCell,
                                           const Array2DReal &NormalVelEdge,
                                           const Real Dt) const {

      // Temporary for stacked shallow water
      const int KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      Real TmpProv[VecLength] = {0.};

      Real DtInvAreaCell = Dt / AreaCell(ICell);
      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);
         const Real Factor =
             DtInvAreaCell * DvEdge(JEdge) * EdgeSignOnCell(ICell, J);
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            TmpProv[KVec] += Factor * FluxPseudoThickEdge(JEdge, K) *
                             NormalVelEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K = KStart + KVec;
         ProvPseudoThickness(ICell, K) =
             PseudoThickCell(ICell, K) + TmpProv[KVec];
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DReal AreaCell;
   Array1DReal DvEdge;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DReal EdgeSignOnCell;
   Array2DI4 CellsOnEdge;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
