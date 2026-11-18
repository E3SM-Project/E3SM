#ifndef OMEGA_AUX_THICKNESS_H
#define OMEGA_AUX_THICKNESS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

enum class FluxThickEdgeOption { Center, Upwind };

class PseudoThicknessAuxVars {
 public:
   Array2DReal FluxPseudoThickEdge;
   Array2DReal MeanPseudoThickEdge;
   Array2DReal ProvPseudoThickness;

   FluxThickEdgeOption FluxThickEdgeChoice;

   PseudoThicknessAuxVars(const std::string &AuxStateSuffix,
                          const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnEdge(const TeamMember &Team, int IEdge,
                     const Array2DReal &PseudoThickCell,
                     const Array2DReal &NormalVelEdge) const {

      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      const int KMin = MinLayerEdgeBot(IEdge);
      const int KMax = MaxLayerEdgeTop(IEdge);

      // Real(0.5) rather than 0.5_Real: nvcc 12.9 on pm-gpu fails to find the
      // _Real literal operator here
      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             MeanPseudoThickEdge(IEdge, K) =
                 Real(0.5) *
                 (PseudoThickCell(JCell0, K) + PseudoThickCell(JCell1, K));
          });

      switch (FluxThickEdgeChoice) {
      case FluxThickEdgeOption::Center:
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                FluxPseudoThickEdge(IEdge, K) =
                    0.5_Real *
                    (PseudoThickCell(JCell0, K) + PseudoThickCell(JCell1, K));
             });
         break;
      case FluxThickEdgeOption::Upwind:
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                if (NormalVelEdge(IEdge, K) > 0) {
                   FluxPseudoThickEdge(IEdge, K) = PseudoThickCell(JCell0, K);
                } else if (NormalVelEdge(IEdge, K) < 0) {
                   FluxPseudoThickEdge(IEdge, K) = PseudoThickCell(JCell1, K);
                } else {
                   FluxPseudoThickEdge(IEdge, K) = Kokkos::max(
                       PseudoThickCell(JCell0, K), PseudoThickCell(JCell1, K));
                }
             });
         break;
      }
   }

   KOKKOS_FUNCTION void computeVarsOnCells(const TeamMember &Team, int ICell,
                                           const Array2DReal &PseudoThickCell,
                                           const Array2DReal &NormalVelEdge,
                                           const Real Dt) const {

      // Temporary for stacked shallow water
      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      ScratchArray1DReal TmpProv(teamScratch(Team), NVertLayers);
      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { TmpProv(K) = 0; });

      Real DtInvAreaCell = Dt / AreaCell(ICell);
      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         const Real Factor =
             DtInvAreaCell * DvEdge(JEdge) * EdgeSignOnCell(ICell, J);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                TmpProv(K) += Factor * FluxPseudoThickEdge(JEdge, K) *
                              NormalVelEdge(JEdge, K);
             });
      }

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             ProvPseudoThickness(ICell, K) =
                 PseudoThickCell(ICell, K) + TmpProv(K);
          });
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
   I4 NVertLayers;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
