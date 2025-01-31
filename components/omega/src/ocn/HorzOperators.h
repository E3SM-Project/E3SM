#ifndef OMEGA_HORZOPERATORS_H
#define OMEGA_HORZOPERATORS_H

#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

class DivergenceOnCell {
 public:
   DivergenceOnCell(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &DivCell, int ICell,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart       = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivCellTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            DivCellTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                                VecEdge(JEdge, K) * InvAreaCell;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K       = KStart + KVec;
         DivCell(ICell, K) = DivCellTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
};

class GradientOnEdge {
 public:
   GradientOnEdge(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &GradEdge, int IEdge,
                                   int KChunk,
                                   const Array2DReal &ScalarCell) const {
      const int KStart     = KChunk * VecLength;
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);
      const auto JCell0    = CellsOnEdge(IEdge, 0);
      const auto JCell1    = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         GradEdge(IEdge, K) =
             InvDcEdge * (ScalarCell(JCell1, K) - ScalarCell(JCell0, K));
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
};

class CurlOnVertex {
 public:
   CurlOnVertex(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &CurlVertex, int IVertex,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart           = KChunk * VecLength;
      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      Real CurlVertexTmp[VecLength] = {0};

      for (int J = 0; J < VertexDegree; ++J) {
         const int JEdge = EdgesOnVertex(IVertex, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            CurlVertexTmp[KVec] += DcEdge(JEdge) *
                                   EdgeSignOnVertex(IVertex, J) *
                                   VecEdge(JEdge, K) * InvAreaTriangle;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K            = KStart + KVec;
         CurlVertex(IVertex, K) = CurlVertexTmp[KVec];
      }
   }

 private:
   I4 VertexDegree;
   Array2DI4 EdgesOnVertex;
   Array1DReal DcEdge;
   Array1DReal AreaTriangle;
   Array2DReal EdgeSignOnVertex;
};

class TangentialReconOnEdge {
 public:
   TangentialReconOnEdge(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &ReconEdge, int IEdge,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart = KChunk * VecLength;

      Real ReconEdgeTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const int JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            ReconEdgeTmp[KVec] += WeightsOnEdge(IEdge, J) * VecEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K         = KStart + KVec;
         ReconEdge(IEdge, K) = ReconEdgeTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
};

} // namespace OMEGA
#endif
