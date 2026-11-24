#ifndef OMEGA_AUX_KINETIC_H
#define OMEGA_AUX_KINETIC_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

#include <string>

namespace OMEGA {

class KineticAuxVars {
 public:
   Array2DReal KineticEnergyCell;
   Array2DReal VelocityDivCell;

   KineticAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                  const VertCoord *VCoord);

   KOKKOS_FUNCTION void
   computeVarsOnCell(const TeamMember &Team, int ICell,
                     const Array2DReal &NormalVelEdge) const {

      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      ScratchArray1DReal KineticEnergyCellTmp(teamScratch(Team), NVertLayers);
      ScratchArray1DReal VelocityDivCellTmp(teamScratch(Team), NVertLayers);

      parallelForInner(
          Team, NVertLayers, INNER_LAMBDA(int K) {
             KineticEnergyCellTmp(K) = 0;
             VelocityDivCellTmp(K)   = 0;
          });

      const int KMin = MinLayerCell(ICell);
      const int KMax = MaxLayerCell(ICell);

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge     = EdgesOnCell(ICell, J);
         const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);
         parallelForInner(
             Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                KineticEnergyCellTmp(K) += AreaEdge * 0.5_Real * InvAreaCell *
                                           NormalVelEdge(JEdge, K) *
                                           NormalVelEdge(JEdge, K);
                VelocityDivCellTmp(K) -= DvEdge(JEdge) * InvAreaCell *
                                         EdgeSignOnCell(ICell, J) *
                                         NormalVelEdge(JEdge, K);
             });
      }

      parallelForInner(
          Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
             KineticEnergyCell(ICell, K) = KineticEnergyCellTmp(K);
             VelocityDivCell(ICell, K)   = VelocityDivCellTmp(K);
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
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

} // namespace OMEGA
#endif
