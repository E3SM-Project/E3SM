#ifndef OMEGA_AUX_VORTICITY_H
#define OMEGA_AUX_VORTICITY_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class VorticityAuxVars {
 public:
   Array2DReal RelVortVertex;
   Array2DReal NormRelVortVertex;
   Array2DReal NormPlanetVortVertex;

   Array2DReal NormRelVortEdge;
   Array2DReal NormPlanetVortEdge;

   VorticityAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                    const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnVertex(int IVertex, int KChunk,
                       const Array2DReal &PseudoThickCell,
                       const Array2DReal &NormalVelEdge) const {

      const int KStartVertex = chunkStart(KChunk, MinLayerVertexTop(IVertex));
      const int KLenVertex =
          chunkLength(KChunk, KStartVertex, MaxLayerVertexBot(IVertex));
      const int KEndVertex = KStartVertex + KLenVertex - 1;

      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      Real PseudoThickVertex[VecLength] = {0};
      Real RelVortVertexTmp[VecLength]  = {0};

      for (int J = 0; J < VertexDegree; ++J) {
         const int JCell = CellsOnVertex(IVertex, J);

         const int KStartCell = Kokkos::max(KStartVertex, MinLayerCell(JCell));
         const int KEndCell   = Kokkos::min(KEndVertex, MaxLayerCell(JCell));

         for (int K = KStartCell; K <= KEndCell; ++K) {
            const int KVec = K - KStartVertex;
            PseudoThickVertex[KVec] += InvAreaTriangle *
                                       KiteAreasOnVertex(IVertex, J) *
                                       PseudoThickCell(JCell, K);
         }

         const int JEdge = EdgesOnVertex(IVertex, J);
         const int KStartEdge =
             Kokkos::max(KStartVertex, MinLayerEdgeTop(JEdge));
         const int KEndEdge = Kokkos::min(KEndVertex, MaxLayerEdgeBot(JEdge));

         for (int K = KStartEdge; K <= KEndEdge; ++K) {
            const int KVec = K - KStartVertex;
            RelVortVertexTmp[KVec] += InvAreaTriangle * DcEdge(JEdge) *
                                      EdgeSignOnVertex(IVertex, J) *
                                      NormalVelEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < KLenVertex; ++KVec) {
         const int K                     = KStartVertex + KVec;
         const Real InvPseudoThickVertex = 1._Real / PseudoThickVertex[KVec];

         RelVortVertex(IVertex, K) = RelVortVertexTmp[KVec];
         NormRelVortVertex(IVertex, K) =
             RelVortVertexTmp[KVec] * InvPseudoThickVertex;
         NormPlanetVortVertex(IVertex, K) =
             FVertex(IVertex) * InvPseudoThickVertex;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnEdge(int IEdge, int KChunk) const {
      const int KStart   = chunkStart(KChunk, MinLayerEdgeTop(IEdge));
      const int KLen     = chunkLength(KChunk, KStart, MaxLayerEdgeBot(IEdge));
      const int JVertex0 = VerticesOnEdge(IEdge, 0);
      const int JVertex1 = VerticesOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const int K = KStart + KVec;
         NormRelVortEdge(IEdge, K) =
             0.5_Real *
             (NormRelVortVertex(JVertex0, K) + NormRelVortVertex(JVertex1, K));

         NormPlanetVortEdge(IEdge, K) =
             0.5_Real * (NormPlanetVortVertex(JVertex0, K) +
                         NormPlanetVortVertex(JVertex1, K));
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   I4 VertexDegree;
   Array2DI4 CellsOnVertex;
   Array2DI4 EdgesOnVertex;
   Array2DReal EdgeSignOnVertex;
   Array1DReal DcEdge;
   Array2DReal KiteAreasOnVertex;
   Array1DReal AreaTriangle;
   Array2DI4 VerticesOnEdge;
   Array1DReal FVertex;

   Array1DI4 MinLayerVertexTop;
   Array1DI4 MaxLayerVertexBot;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeBot;
};

} // namespace OMEGA
#endif
