#ifndef OMEGA_AUX_KINETIC_H
#define OMEGA_AUX_KINETIC_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"

#include <string>

namespace OMEGA {

class KineticAuxVars {
 public:
   Array2DReal KineticEnergyCell;
   Array2DReal VelocityDivCell;

   KineticAuxVars(const std::string &AuxStateSuffix, const HorzMesh *Mesh,
                  int NVertLevels);

   KOKKOS_FUNCTION void
   computeVarsOnCell(int ICell, int KChunk,
                     const Array2DReal &NormalVelEdge) const {
      const Real InvAreaCell = 1._Real / AreaCell(ICell);
      const int KStart       = KChunk * VecLength;

      Real KineticEnergyCellTmp[VecLength] = {0};
      Real VelocityDivCellTmp[VecLength]   = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge     = EdgesOnCell(ICell, J);
         const Real AreaEdge = 0.5_Real * DvEdge(JEdge) * DcEdge(JEdge);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            KineticEnergyCellTmp[KVec] += AreaEdge * 0.5_Real * InvAreaCell *
                                          NormalVelEdge(JEdge, K) *
                                          NormalVelEdge(JEdge, K);
            VelocityDivCellTmp[KVec] -= DvEdge(JEdge) * InvAreaCell *
                                        EdgeSignOnCell(ICell, J) *
                                        NormalVelEdge(JEdge, K);
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K                 = KStart + KVec;
         KineticEnergyCell(ICell, K) = KineticEnergyCellTmp[KVec];
         VelocityDivCell(ICell, K)   = VelocityDivCellTmp[KVec];
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
};

} // namespace OMEGA
#endif
