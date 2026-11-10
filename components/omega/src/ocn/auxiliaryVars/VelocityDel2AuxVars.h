#ifndef OMEGA_AUX_VELDEL2_H
#define OMEGA_AUX_VELDEL2_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class VelocityDel2AuxVars {
 public:
   Array2DReal Del2Edge;
   Array2DReal Del2DivCell;
   Array2DReal Del2RelVortVertex;

   VelocityDel2AuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                       const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(const TeamMember &Team, int IEdge,
                     const Array2DReal &VelocityDivCell,
                     const Array2DReal &RelVortVertex) const {

      const int JCell0   = CellsOnEdge(IEdge, 0);
      const int JCell1   = CellsOnEdge(IEdge, 1);
      const int JVertex0 = VerticesOnEdge(IEdge, 0);
      const int JVertex1 = VerticesOnEdge(IEdge, 1);

      const Real InvDcEdge = 1._Real / DcEdge(IEdge);
      const Real InvDvEdge =
          1._Real / Kokkos::max(DvEdge(IEdge), 0.25_Real * DcEdge(IEdge));

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             const Real GradDiv =
                 (VelocityDivCell(JCell1, K) - VelocityDivCell(JCell0, K)) *
                 InvDcEdge;
             const Real CurlVort =
                 -(RelVortVertex(JVertex1, K) - RelVortVertex(JVertex0, K)) *
                 InvDvEdge;
             Del2Edge(IEdge, K) = EdgeMask(IEdge, K) * GradDiv + CurlVort;
          });
   }

   KOKKOS_FUNCTION void computeVarsOnCell(const TeamMember &Team,
                                          int ICell) const {
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal Del2DivCellTmp(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { Del2DivCellTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge     = EdgesOnCell(ICell, J);
         const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                Del2DivCellTmp(K) -= DvEdge(JEdge) * InvAreaCell *
                                     EdgeSignOnCell(ICell, J) *
                                     Del2Edge(JEdge, K);
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell},
          INNER_LAMBDA(int K) { Del2DivCell(ICell, K) = Del2DivCellTmp(K); });
   }

   KOKKOS_FUNCTION void computeVarsOnVertex(const TeamMember &Team,
                                            int IVertex) const {
      // Compute over the full vertex valid range [MinLayerVertexTop,
      // MaxLayerVertexBot] so that boundary-vertex layers (where only some
      // surrounding cells are active) receive a valid value. Each edge's
      // contribution is clamped to that edge's valid range [MinLayerEdgeTop,
      // MaxLayerEdgeBot], where Del2Edge has been computed or zeroed; this
      // matches VorticityAuxVars::computeVarsOnVertex and avoids reading
      // uninitialized (fill-value) layers of Del2Edge for deeper edges.

      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      ScratchArray1DReal Del2RelVortVertexTmp(teamScratch(Team), NVertLayers);

      parallelForInner(
          Team, NVertLayers,
          INNER_LAMBDA(int K) { Del2RelVortVertexTmp(K) = 0; });

      const int MinLyrVertexTop = MinLayerVertexTop(IVertex);
      const int MaxLyrVertexBot = MaxLayerVertexBot(IVertex);

      for (int J = 0; J < VertexDegree; ++J) {
         const int JEdge = EdgesOnVertex(IVertex, J);

         const int MinLyrEdgeTop = MinLayerEdgeTop(JEdge);
         const int MaxLyrEdgeBot = MaxLayerEdgeBot(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeTop, MaxLyrEdgeBot}, INNER_LAMBDA(int K) {
                Del2RelVortVertexTmp(K) += InvAreaTriangle * DcEdge(JEdge) *
                                           EdgeSignOnVertex(IVertex, J) *
                                           Del2Edge(JEdge, K);
             });
      }

      parallelForInner(
          Team, Range{MinLyrVertexTop, MaxLyrVertexBot}, INNER_LAMBDA(int K) {
             Del2RelVortVertex(IVertex, K) = Del2RelVortVertexTmp(K);
          });
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
   Array2DReal EdgeMask;
   I4 VertexDegree;
   I4 NVertLayers;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeBot;
   Array1DI4 MinLayerVertexTop;
   Array1DI4 MaxLayerVertexBot;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
