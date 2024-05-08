#ifndef OMEGA_HORZOPERATORS_H
#define OMEGA_HORZOPERATORS_H

#include "DataTypes.h"
#include "HorzMesh.h"

namespace OMEGA {

class DivergenceOnCell {
 public:
   DivergenceOnCell(HorzMesh const *mesh);

   KOKKOS_FUNCTION Real operator()(int ICell,
                                   const Array1DReal &VecEdge) const {
      Real DivCell = 0;
      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);
         DivCell -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) * VecEdge(JEdge);
      }
      const Real InvAreaCell = 1. / AreaCell(ICell);
      DivCell *= InvAreaCell;
      return DivCell;
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;
   Array2DR8 EdgeSignOnCell;
};

class GradientOnEdge {
 public:
   GradientOnEdge(HorzMesh const *mesh);

   KOKKOS_FUNCTION Real operator()(int IEdge,
                                   const Array1DReal &ScalarCell) const {
      const auto JCell0    = CellsOnEdge(IEdge, 0);
      const auto JCell1    = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1. / DcEdge(IEdge);
      const Real GradEdge =
          InvDcEdge * (ScalarCell(JCell1) - ScalarCell(JCell0));
      return GradEdge;
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DR8 DcEdge;
};

class CurlOnVertex {
 public:
   CurlOnVertex(HorzMesh const *mesh);

   KOKKOS_FUNCTION Real operator()(int IVertex,
                                   const Array1DReal &VecEdge) const {
      Real CurlVertex = 0;
      for (int J = 0; J < VertexDegree; ++J) {
         const int JEdge = EdgesOnVertex(IVertex, J);
         CurlVertex +=
             DcEdge(JEdge) * EdgeSignOnVertex(IVertex, J) * VecEdge(JEdge);
      }
      const Real InvAreaTriangle = 1. / AreaTriangle(IVertex);
      CurlVertex *= InvAreaTriangle;
      return CurlVertex;
   }

 private:
   I4 VertexDegree;
   Array2DI4 EdgesOnVertex;
   Array1DR8 DcEdge;
   Array1DR8 AreaTriangle;
   Array2DR8 EdgeSignOnVertex;
};

class TangentialReconOnEdge {
 public:
   TangentialReconOnEdge(HorzMesh const *mesh);

   KOKKOS_FUNCTION Real operator()(int IEdge,
                                   const Array1DReal &VecEdge) const {
      Real ReconEdge = 0;
      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const int JEdge = EdgesOnEdge(IEdge, J);
         ReconEdge += WeightsOnEdge(IEdge, J) * VecEdge(JEdge);
      }
      return ReconEdge;
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DR8 WeightsOnEdge;
};

} // namespace OMEGA
#endif
