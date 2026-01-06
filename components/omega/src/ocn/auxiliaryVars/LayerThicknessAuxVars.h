#ifndef OMEGA_AUX_THICKNESS_H
#define OMEGA_AUX_THICKNESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

enum class FluxThickEdgeOption { Center, Upwind };

class LayerThicknessAuxVars {
 public:
   Array2DReal FluxLayerThickEdge;
   Array2DReal MeanLayerThickEdge;
   Array2DReal SshCell;
   Array2DReal ProvThickness;

   FluxThickEdgeOption FluxThickEdgeChoice;

   LayerThicknessAuxVars(const std::string &AuxStateSuffix,
                         const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &LayerThickCell,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = chunkStart(KChunk, MinLayerEdgeBot(IEdge));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerEdgeTop(IEdge));

      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K = KStart + KVec;
         MeanLayerThickEdge(IEdge, K) =
             0.5_Real * (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
      }

      switch (FluxThickEdgeChoice) {
      case FluxThickEdgeOption::Center:
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            FluxLayerThickEdge(IEdge, K) =
                0.5_Real *
                (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
         }
         break;
      case FluxThickEdgeOption::Upwind:
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            if (NormalVelEdge(IEdge, K) > 0) {
               FluxLayerThickEdge(IEdge, K) = LayerThickCell(JCell0, K);
            } else if (NormalVelEdge(IEdge, K) < 0) {
               FluxLayerThickEdge(IEdge, K) = LayerThickCell(JCell1, K);
            } else {
               FluxLayerThickEdge(IEdge, K) = Kokkos::max(
                   LayerThickCell(JCell0, K), LayerThickCell(JCell1, K));
            }
         }
         break;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnCells(int ICell, int KChunk,
                                           const Array2DReal &LayerThickCell,
                                           const Array2DReal &NormalVelEdge,
                                           const Real Dt) const {

      // Temporary for stacked shallow water
      const int KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const int KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K       = KStart + KVec;
         SshCell(ICell, K) = LayerThickCell(ICell, K) - BottomDepth(ICell);
      }

      Real TmpProv[VecLength] = {0.};

      Real DtInvAreaCell = Dt / AreaCell(ICell);
      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);
         const Real Factor =
             DtInvAreaCell * DvEdge(JEdge) * EdgeSignOnCell(ICell, J);
         for (int KVec = 0; KVec < KLen; ++KVec) {
            const int K = KStart + KVec;
            TmpProv[KVec] +=
                Factor * FluxLayerThickEdge(JEdge, K) * NormalVelEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K             = KStart + KVec;
         ProvThickness(ICell, K) = LayerThickCell(ICell, K) + TmpProv[KVec];
      }

      /*
      Real TotalThickness = 0.0;
      for (int K = 0; K < NVertLayers; K++) {
         TotalThickness += LayerThickCell(ICell, K);
      }

      SshCell(ICell) = TotalThickness - BottomDepth(ICell);
      */
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
   Array1DReal BottomDepth;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
