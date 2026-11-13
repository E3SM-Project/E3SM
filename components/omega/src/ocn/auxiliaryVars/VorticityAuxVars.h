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
   computeVarsOnVertex(const TeamMember &Team, int IVertex,
                       const Array2DReal &PseudoThickCell,
                       const Array2DReal &NormalVelEdge) const {

      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      ScratchArray1DReal PseudoThickVertex(teamScratch(Team), NVertLayers);
      ScratchArray1DReal RelVortVertexTmp(teamScratch(Team), NVertLayers);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) {
             PseudoThickVertex(K) = 0;
             RelVortVertexTmp(K)  = 0;
          });

      for (int J = 0; J < VertexDegree; ++J) {
         const int JCell      = CellsOnVertex(IVertex, J);
         const int MinLyrCell = MinLayerCell(JCell);
         const int MaxLyrCell = MaxLayerCell(JCell);

         parallelForInner(
             Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
                PseudoThickVertex(K) += InvAreaTriangle *
                                        KiteAreasOnVertex(IVertex, J) *
                                        PseudoThickCell(JCell, K);
             });

         const int JEdge         = EdgesOnVertex(IVertex, J);
         const int MinLyrEdgeTop = MinLayerEdgeTop(JEdge);
         const int MaxLyrEdgeBot = MaxLayerEdgeBot(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeTop, MaxLyrEdgeBot}, INNER_LAMBDA(int K) {
                RelVortVertexTmp(K) += InvAreaTriangle * DcEdge(JEdge) *
                                       EdgeSignOnVertex(IVertex, J) *
                                       NormalVelEdge(JEdge, K);
             });
      }

      const int MinLyrVertexTop = MinLayerVertexTop(IVertex);
      const int MaxLyrVertexBot = MaxLayerVertexBot(IVertex);

      parallelForInner(
          Team, Range{MinLyrVertexTop, MaxLyrVertexBot}, INNER_LAMBDA(int K) {
             const Real InvPseudoThickVertex = 1._Real / PseudoThickVertex(K);

             RelVortVertex(IVertex, K) = RelVortVertexTmp(K);
             NormRelVortVertex(IVertex, K) =
                 RelVortVertexTmp(K) * InvPseudoThickVertex;
             NormPlanetVortVertex(IVertex, K) =
                 FVertex(IVertex) * InvPseudoThickVertex;
          });
   }

   KOKKOS_FUNCTION void computeVarsOnEdge(const TeamMember &Team,
                                          int IEdge) const {
      const int JVertex0 = VerticesOnEdge(IEdge, 0);
      const int JVertex1 = VerticesOnEdge(IEdge, 1);

      const int KMin = MinLayerEdgeTop(IEdge);
      const int KMax = MaxLayerEdgeBot(IEdge);

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             NormRelVortEdge(IEdge, K) =
                 0.5_Real * (NormRelVortVertex(JVertex0, K) +
                             NormRelVortVertex(JVertex1, K));

             NormPlanetVortEdge(IEdge, K) =
                 0.5_Real * (NormPlanetVortVertex(JVertex0, K) +
                             NormPlanetVortVertex(JVertex1, K));
          });
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

   I4 NVertLayers;
   Array1DI4 MinLayerVertexTop;
   Array1DI4 MaxLayerVertexBot;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeTop;
   Array1DI4 MaxLayerEdgeBot;
};

} // namespace OMEGA
#endif
