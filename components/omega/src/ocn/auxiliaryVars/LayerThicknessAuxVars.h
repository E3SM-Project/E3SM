#ifndef OMEGA_AUX_THICKNESS_H
#define OMEGA_AUX_THICKNESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

namespace OMEGA {

class LayerThicknessAuxVars {
 public:
   Array2DReal FluxLayerThickEdge;
   Array2DReal MeanLayerThickEdge;

   LayerThicknessAuxVars(const HorzMesh *mesh, int NVertLevels)
       : FluxLayerThickEdge("FluxLayerThickEdge", mesh->NEdgesSize,
                            NVertLevels),
         MeanLayerThickEdge("MeanLayerThickEdge", mesh->NEdgesSize,
                            NVertLevels),
         CellsOnEdge(mesh->CellsOnEdge) {}

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &LayerThickCell,
                     const Array2DReal &NormalVelEdge) const {
      const int KStart = KChunk * VecLength;
      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         // TODO: add upwind option
         FluxLayerThickEdge(IEdge, K) =
             0.5_Real * (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));

         MeanLayerThickEdge(IEdge, K) =
             0.5_Real * (LayerThickCell(JCell0, K) + LayerThickCell(JCell1, K));
      }
   }

 private:
   Array2DI4 CellsOnEdge;
};

} // namespace OMEGA
#endif
