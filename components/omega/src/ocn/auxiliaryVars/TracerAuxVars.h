#ifndef OMEGA_AUX_TRACER_H
#define OMEGA_AUX_TRACER_H

#include "DataTypes.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

enum class FluxTracerEdgeOption { Center, Upwind };

class TracerAuxVars {
 public:
   Array3DReal HTracersEdge;
   Array3DReal Del2TracersCell;

   FluxTracerEdgeOption TracersOnEdgeChoice;

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
      case FluxTracerEdgeOption::Center:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            HTracersEdge(L, IEdge, K) =
                0.5_Real * (HCell(JCell0, K) * TrCell(L, JCell0, K) +
                            HCell(JCell1, K) * TrCell(L, JCell1, K));
         }
         break;
      case FluxTracerEdgeOption::Upwind:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            if (NormalVelEdge(IEdge, K) > 0) {
               HTracersEdge(L, IEdge, K) =
                   HCell(JCell0, K) * TrCell(L, JCell0, K);
            } else if (NormalVelEdge(IEdge, K) < 0) {
               HTracersEdge(L, IEdge, K) =
                   HCell(JCell1, K) * TrCell(L, JCell1, K);
            } else {
               HTracersEdge(L, IEdge, K) =
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
         const int K                  = KStart + KVec;
         Del2TracersCell(L, ICell, K) = Del2TrCellTmp[KVec] * InvAreaCell;
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
};

} // namespace OMEGA
#endif
