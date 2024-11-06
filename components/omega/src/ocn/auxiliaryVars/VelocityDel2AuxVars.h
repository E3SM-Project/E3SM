#ifndef OMEGA_AUX_VELDEL2_H
#define OMEGA_AUX_VELDEL2_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class VelocityDel2AuxVars {
 public:
   Array2DReal Del2Edge;
   Array2DReal Del2DivCell;
   Array2DReal Del2RelVortVertex;

   VelocityDel2AuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                       int NVertLevels);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(int IEdge, int KChunk, const Array2DReal &VelocityDivCell,
                     const Array2DReal &RelVortVertex) const {
      const int KStart = KChunk * VecLength;

      const int JCell0   = CellsOnEdge(IEdge, 0);
      const int JCell1   = CellsOnEdge(IEdge, 1);
      const int JVertex0 = VerticesOnEdge(IEdge, 0);
      const int JVertex1 = VerticesOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);
      const Real InvDvEdge =
          1._Real / Kokkos::max(DvEdge(IEdge), 0.25_Real * DcEdge(IEdge));

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         const Real GradDiv =
             (VelocityDivCell(JCell1, K) - VelocityDivCell(JCell0, K)) *
             InvDcEdge;
         const Real CurlVort =
             -(RelVortVertex(JVertex1, K) - RelVortVertex(JVertex0, K)) *
             InvDvEdge;
         Del2Edge(IEdge, K) = GradDiv + CurlVort;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnCell(int ICell, int KChunk) const {
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      const int KStart       = KChunk * VecLength;

      Real Del2DivCellTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge     = EdgesOnCell(ICell, J);
         const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            Del2DivCellTmp[KVec] -= DvEdge(JEdge) * InvAreaCell *
                                    EdgeSignOnCell(ICell, J) *
                                    Del2Edge(JEdge, K);
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K           = KStart + KVec;
         Del2DivCell(ICell, K) = Del2DivCellTmp[KVec];
      }
   }

   KOKKOS_FUNCTION void computeVarsOnVertex(int IVertex, int KChunk) const {
      const int KStart           = KChunk * VecLength;
      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      Real Del2RelVortVertexTmp[VecLength] = {0};

      for (int J = 0; J < VertexDegree; ++J) {
         const int JEdge = EdgesOnVertex(IVertex, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            Del2RelVortVertexTmp[KVec] += InvAreaTriangle * DcEdge(JEdge) *
                                          EdgeSignOnVertex(IVertex, J) *
                                          Del2Edge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K                   = KStart + KVec;
         Del2RelVortVertex(IVertex, K) = Del2RelVortVertexTmp[KVec];
      }
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DReal EdgeSignOnCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DI4 EdgesOnVertex;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array2DReal EdgeSignOnVertex;
   Array1DReal AreaTriangle;
   I4 VertexDegree;
};

} // namespace OMEGA
#endif
