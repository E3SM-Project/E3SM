#ifndef OMEGA_AUX_TRACER_H
#define OMEGA_AUX_TRACER_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "auxiliaryVars/LayerThicknessAuxVars.h"

#include <string>

namespace OMEGA {

class TracerAuxVars {
 public:
   Array3DReal HTracersOnEdge;
   Array3DReal Del2TracersOnCell;

   FluxThickEdgeOption TracersOnEdgeChoice = Center;

   TracerAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                 const I4 NVertLevels, const I4 NTracers);

   KOKKOS_FUNCTION void computeVarsOnEdge(int L, int IEdge, int KChunk,
                                          const Array2DReal &NormalVelEdge,
                                          const Array2DReal &HCell,
                                          const Array3DReal &TrCell) const {
      const int KStart = KChunk * VecLength;
      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      switch (TracersOnEdgeChoice) {
      case Center:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            HTracersOnEdge(L, IEdge, K) =
                0.5_Real * (HCell(JCell0, K) * TrCell(L, JCell0, K) +
                            HCell(JCell1, K) * TrCell(L, JCell1, K));
         }
         break;
      case Upwind:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            if (NormalVelEdge(IEdge, K) > 0) {
               HTracersOnEdge(L, IEdge, K) =
                   HCell(JCell0, K) * TrCell(L, JCell0, K);
            } else if (NormalVelEdge(IEdge, K) < 0) {
               HTracersOnEdge(L, IEdge, K) =
                   HCell(JCell1, K) * TrCell(L, JCell1, K);
            } else {
               HTracersOnEdge(L, IEdge, K) =
                   Kokkos::max(HCell(JCell0, K) * TrCell(L, JCell0, K),
                               HCell(JCell1, K) * TrCell(L, JCell1, K));
            }
         }
         break;
      }
   }

   KOKKOS_FUNCTION void
   computeVarsOnCells(int L, int ICell, int KChunk,
                      const Array2DReal &LayerThickEdgeMean,
                      const Array3DReal &TrCell) const {

      const int KStart       = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real Del2TrCellTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);

         const int JCell0 = CellsOnEdge(JEdge, 0);
         const int JCell1 = CellsOnEdge(JEdge, 1);

         const Real DvDcEdge = DvEdge(JEdge) / DcEdge(JEdge);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K           = KStart + KVec;
            const Real TracerGrad = TrCell(L, JCell1, K) - TrCell(L, JCell0, K);
            Del2TrCellTmp[KVec] -= EdgeSignOnCell(ICell, J) * DvDcEdge *
                                   LayerThickEdgeMean(JEdge, K) * TracerGrad;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K                    = KStart + KVec;
         Del2TracersOnCell(L, ICell, K) = Del2TrCellTmp[KVec] * InvAreaCell;
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;
};

} // namespace OMEGA
#endif
