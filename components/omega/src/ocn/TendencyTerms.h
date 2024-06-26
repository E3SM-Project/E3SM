#ifndef OMEGA_TENDENCYTERMS_H
#define OMEGA_TENDENCYTERMS_H

#include "Config.h"
#include "HorzMesh.h"
#include "MachEnv.h"

namespace OMEGA {

class ThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   ThicknessFluxDivOnCell(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell, I4 KChunk,
                                   const Array2DR8 &ThicknessFlux) const {

      const I4 KStart        = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            DivTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                            ThicknessFlux(JEdge, K) * InvAreaCell;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(ICell, K) -= DivTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;
   Array2DR8 EdgeSignOnCell;
};

class PotentialVortHAdvOnEdge {
 public:
   bool Enabled = false;

   PotentialVortHAdvOnEdge(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &NormRVortEdge,
                                   const Array2DR8 &NormFEdge,
                                   const Array2DR8 &LayerThickEdge,
                                   const Array2DR8 &NormVelEdge) const {

      const I4 KStart         = KChunk * VecLength;
      Real VortTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         I4 JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K    = KStart + KVec;
            Real NormVort = (NormRVortEdge(IEdge, K) + NormFEdge(IEdge, K) +
                             NormRVortEdge(JEdge, K) + NormFEdge(JEdge, K)) *
                            0.5_Real;

            VortTmp[KVec] += WeightsOnEdge(IEdge, J) *
                             LayerThickEdge(JEdge, K) * NormVelEdge(JEdge, K) *
                             NormVort;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += VortTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DR8 WeightsOnEdge;
};

class KEGradOnEdge {
 public:
   bool Enabled = false;

   KEGradOnEdge(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &KECell) const {

      const I4 KStart      = KChunk * VecLength;
      const I4 JCell0      = CellsOnEdge(IEdge, 0);
      const I4 JCell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -= (KECell(JCell1, K) - KECell(JCell0, K)) * InvDcEdge;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DR8 DcEdge;
};

class SSHGradOnEdge {
 public:
   bool Enabled = false;

   SSHGradOnEdge(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &HCell) const {

      const I4 KStart      = KChunk * VecLength;
      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -=
             Grav * (HCell(ICell1, K) - HCell(ICell0, K)) * InvDcEdge;
      }
   }

 private:
   R8 Grav;
   Array2DI4 CellsOnEdge;
   Array1DR8 DcEdge;
};

class VelocityDiffusionOnEdge {
 public:
   bool Enabled = false;

   VelocityDiffusionOnEdge(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &DivCell,
                                   const Array2DR8 &RVortVertex) const {

      const I4 KStart = KChunk * VecLength;
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             ((DivCell(ICell1, K) - DivCell(ICell0, K)) * DcEdgeInv -
              (RVortVertex(IVertex1, K) - RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * ViscDel2 * MeshScalingDel2(IEdge) * Del2U;
      }
   }

 private:
   R8 ViscDel2;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 MeshScalingDel2;
   Array2DR8 EdgeMask;
};

class VelocityHyperDiffOnEdge {
 public:
   bool Enabled = false;

   VelocityHyperDiffOnEdge(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &Del2DivCell,
                                   const Array2DR8 &Del2RVortVertex) const {

      const I4 KStart = KChunk * VecLength;
      const I4 ICell0 = CellsOnEdge(IEdge, 0);
      const I4 ICell1 = CellsOnEdge(IEdge, 1);

      const I4 IVertex0 = VerticesOnEdge(IEdge, 0);
      const I4 IVertex1 = VerticesOnEdge(IEdge, 1);

      const Real DcEdgeInv = 1._Real / DcEdge(IEdge);
      const Real DvEdgeInv = 1._Real / DvEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         const Real Del2U =
             ((Del2DivCell(ICell1, K) - Del2DivCell(ICell0, K)) * DcEdgeInv -
              (Del2RVortVertex(IVertex1, K) - Del2RVortVertex(IVertex0, K)) *
                  DvEdgeInv);

         Tend(IEdge, K) -=
             EdgeMask(IEdge, K) * ViscDel4 * MeshScalingDel4(IEdge) * Del2U;
      }
   }

 private:
   Real ViscDel4;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 MeshScalingDel4;
   Array2DR8 EdgeMask;
};

class TracerHorzAdvOnCell {
 public:
   bool Enabled = false;

   TracerHorzAdvOnCell(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk, const Array2DR8 &VEdge,
                                   const Array3DR8 &NormTrCell,
                                   const Array2DR8 &HFluxEdge) const {

      const I4 KStart  = KChunk * VecLength;
      Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HAdvTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            const Real NormTrEdge =
                (NormTrCell(L, JCell0, K) + NormTrCell(L, JCell1, K)) *
                0.5_Real;
            HAdvTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                             HFluxEdge(JEdge, K) * NormTrEdge *
                             VEdge(JEdge, K) * InvAreaCell;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) += HAdvTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DvEdge;
   Array1DR8 AreaCell;
};

class TracerDiffOnCell {
 public:
   bool Enabled = false;

   TracerDiffOnCell(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk, const Array3DR8 &NormTrCell,
                                   const Array2DR8 &HMeanEdge) const {

      const I4 KStart  = KChunk * VecLength;
      Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DiffTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real InvDcEdge = 1._Real / DcEdge(JEdge);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            const Real GradTrEdge =
                (NormTrCell(L, JCell1, K) - NormTrCell(L, JCell0, K)) *
                InvDcEdge;

            DiffTmp[KVec] += EddyDiff2 * DvEdge(JEdge) *
                             EdgeSignOnCell(ICell, J) * HMeanEdge(JEdge, K) *
                             MeshScalingDel2(JEdge) * GradTrEdge * InvAreaCell;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) += DiffTmp[KVec];
      }
   }

 private:
   Real EddyDiff2;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DvEdge;
   Array1DR8 DcEdge;
   Array1DR8 AreaCell;
   Array1DR8 MeshScalingDel2;
};

class TracerHyperDiffOnCell {
 public:
   bool Enabled = false;

   TracerHyperDiffOnCell(const HorzMesh *Mesh, Config *Options);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk,
                                   const Array3DR8 &TrDel2Cell) const {

      const I4 KStart  = KChunk * VecLength;
      Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HypTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real InvDcEdge = 1._Real / DcEdge(JEdge);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            const Real GradTrDel2Edge =
                (TrDel2Cell(L, JCell1, K) - TrDel2Cell(L, JCell0, K)) *
                InvDcEdge;

            HypTmp[KVec] -= EddyDiff4 * DvEdge(JEdge) *
                            EdgeSignOnCell(ICell, J) * MeshScalingDel4(JEdge) *
                            GradTrDel2Edge * InvAreaCell;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) += HypTmp[KVec];
      }
   }

 private:
   Real EddyDiff4;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DR8 EdgeSignOnCell;
   Array1DR8 DvEdge;
   Array1DR8 DcEdge;
   Array1DR8 AreaCell;
   Array1DR8 MeshScalingDel4;
};

} // namespace OMEGA
#endif
