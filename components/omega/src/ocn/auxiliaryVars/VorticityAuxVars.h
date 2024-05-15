#ifndef OMEGA_AUX_VORTICITY_H
#define OMEGA_AUX_VORTICITY_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

namespace OMEGA {

class VorticityAuxVars {
 public:
   Array2DReal RelVortVertex;
   Array2DReal NormRelVortVertex;
   Array2DReal NormPlanetVortVertex;

   Array2DReal NormRelVortEdge;
   Array2DReal NormPlanetVortEdge;

   VorticityAuxVars(const HorzMesh *mesh, int NVertLevels)
       : RelVortVertex("RelVortVertex", mesh->NVerticesSize, NVertLevels),
         NormRelVortVertex("NormRelVortVertex", mesh->NVerticesSize,
                           NVertLevels),
         NormPlanetVortVertex("NormPlanetVortVertex", mesh->NVerticesSize,
                              NVertLevels),
         NormRelVortEdge("NormRelVortEdge", mesh->NEdgesSize, NVertLevels),
         NormPlanetVortEdge("NormPlanetVortEdge", mesh->NEdgesSize,
                            NVertLevels),
         VertexDegree(mesh->VertexDegree), CellsOnVertex(mesh->CellsOnVertex),
         EdgesOnVertex(mesh->EdgesOnVertex),
         EdgeSignOnVertex(mesh->EdgeSignOnVertex), DcEdge(mesh->DcEdge),
         KiteAreasOnVertex(mesh->KiteAreasOnVertex),
         AreaTriangle(mesh->AreaTriangle), FVertex(mesh->FVertex),
         VerticesOnEdge(mesh->VerticesOnEdge) {}

   KOKKOS_FUNCTION void
   computeVarsOnVertex(int IVertex, int KChunk,
                       const Array2DReal &LayerThickCell,
                       const Array2DReal &NormalVelEdge) const {
      const int KStart           = KChunk * VecLength;
      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      Real LayerThickVertex[VecLength] = {0};
      Real RelVortVertexTmp[VecLength] = {0};

      for (int J = 0; J < VertexDegree; ++J) {
         const int JCell = CellsOnVertex(IVertex, J);
         const int JEdge = EdgesOnVertex(IVertex, J);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            LayerThickVertex[KVec] += InvAreaTriangle *
                                      KiteAreasOnVertex(IVertex, J) *
                                      LayerThickCell(JCell, K);
            RelVortVertexTmp[KVec] += InvAreaTriangle * DcEdge(JEdge) *
                                      EdgeSignOnVertex(IVertex, J) *
                                      NormalVelEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K                    = KStart + KVec;
         const Real InvLayerThickVertex = 1._Real / LayerThickVertex[KVec];

         RelVortVertex(IVertex, K) = RelVortVertexTmp[KVec];
         NormRelVortVertex(IVertex, K) =
             RelVortVertexTmp[KVec] * InvLayerThickVertex;
         NormPlanetVortVertex(IVertex, K) =
             FVertex(IVertex) * InvLayerThickVertex;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge, int KChunk) const {
      const int KStart   = KChunk * VecLength;
      const int JVertex0 = VerticesOnEdge(IEdge, 0);
      const int JVertex1 = VerticesOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         NormRelVortEdge(IEdge, K) =
             0.5_Real *
             (NormRelVortVertex(JVertex0, K) + NormRelVortVertex(JVertex1, K));

         NormPlanetVortEdge(IEdge, K) =
             0.5_Real * (NormPlanetVortVertex(JVertex0, K) +
                         NormPlanetVortVertex(JVertex1, K));
      }
   }

 private:
   I4 VertexDegree;
   Array2DI4 CellsOnVertex;
   Array2DI4 EdgesOnVertex;
   Array2DR8 EdgeSignOnVertex;
   Array1DR8 DcEdge;
   Array2DR8 KiteAreasOnVertex;
   Array1DR8 AreaTriangle;
   Array2DI4 VerticesOnEdge;
   Array1DR8 FVertex;
};

} // namespace OMEGA
#endif
