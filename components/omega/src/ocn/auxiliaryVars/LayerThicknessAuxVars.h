#ifndef OMEGA_AUX_THICKNESS_H
#define OMEGA_AUX_THICKNESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

enum FluxThickEdgeOption { Center, Upwind };

class LayerThicknessAuxVars {
 public:
   Array2DReal FluxLayerThickEdge;
   Array2DReal MeanLayerThickEdge;
   Array2DReal SshCell;

   FluxThickEdgeOption FluxThickEdgeChoice;

   LayerThicknessAuxVars(const std::string &AuxStateSuffix,
                         const HorzMesh *Mesh, int NVertLevels);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &LayerThickCell,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = KChunk * VecLength;
      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         MeanLayerThickEdge(IEdge, K) =
             0.5_Real * (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
      }

      switch (FluxThickEdgeChoice) {
      case Center:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            FluxLayerThickEdge(IEdge, K) =
                0.5_Real *
                (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
         }
         break;
      case Upwind:
         for (int KVec = 0; KVec < VecLength; ++KVec) {
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

   KOKKOS_FUNCTION void
   computeVarsOnCells(int ICell, int KChunk,
                      const Array2DReal &LayerThickCell) const {

      // Temporary for stacked shallow water
      const int KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K       = KStart + KVec;
         SshCell(ICell, K) = LayerThickCell(ICell, K) - BottomDepth(ICell);
      }

      /*
      Real TotalThickness = 0.0;
      for (int K = 0; K < NVertLevels; K++) {
         TotalThickness += LayerThickCell(ICell, K);
      }

      SshCell(ICell) = TotalThickness - BottomDepth(ICell);
      */
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal BottomDepth;
};

} // namespace OMEGA
#endif
