#ifndef OMEGA_AUX_TRACER_H
#define OMEGA_AUX_TRACER_H

#include "DataTypes.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class TracerAuxVars {
 public:
   Array3DReal Del2TracersCell;

   TracerAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                 const VertCoord *VCoord, const I4 NTracers);

   KOKKOS_FUNCTION void
   computeVarsOnCells(const TeamMember &Team, int L, int ICell,
                      const Array2DReal &MeanPseudoThickEdge,
                      const Array3DReal &TrCell) const {

      const auto LTrCell =
          subviewUnmanaged(TrCell, L, Kokkos::ALL, Kokkos::ALL);
      const auto LDel2TracersCell =
          subviewUnmanaged(Del2TracersCell, L, Kokkos::ALL, Kokkos::ALL);

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal Del2TrCellTmp(teamScratch(Team), NVertLayers);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) { Del2TrCellTmp(K) = 0; });

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);

         const int JCell0 = CellsOnEdge(JEdge, 0);
         const int JCell1 = CellsOnEdge(JEdge, 1);

         const Real DvDcEdge = DvEdge(JEdge) / DcEdge(JEdge);

         const int MinLyrEdgeBot = MinLayerEdgeBot(JEdge);
         const int MaxLyrEdgeTop = MaxLayerEdgeTop(JEdge);

         parallelForInner(
             Team, Range{MinLyrEdgeBot, MaxLyrEdgeTop}, INNER_LAMBDA(int K) {
                const Real TracerGrad = LTrCell(JCell1, K) - LTrCell(JCell0, K);
                Del2TrCellTmp(K) -= EdgeMask(JEdge, K) *
                                    EdgeSignOnCell(ICell, J) * DvDcEdge *
                                    MeanPseudoThickEdge(JEdge, K) * TracerGrad;
             });
      }

      const int MinLyrCell = MinLayerCell(ICell);
      const int MaxLyrCell = MaxLayerCell(ICell);

      parallelForInner(
          Team, Range{MinLyrCell, MaxLyrCell}, INNER_LAMBDA(int K) {
             LDel2TracersCell(ICell, K) = Del2TrCellTmp(K) * InvAreaCell;
          });
   }

   void registerFields(const std::string &AuxGroupName,
                       const std::string &MeshName) const;
   void unregisterFields() const;

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeMask;
   I4 NVertLayers;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
