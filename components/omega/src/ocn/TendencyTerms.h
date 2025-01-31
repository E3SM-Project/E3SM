#ifndef OMEGA_TENDENCYTERMS_H
#define OMEGA_TENDENCYTERMS_H
//===-- ocn/TendencyTerms.h - Tendency Terms --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating tendency terms
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate tendencies used to update state variables.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OceanState.h"

#include <functional>
#include <memory>

namespace OMEGA {

/// Divergence of thickness flux at cell centers, for updating layer thickness
/// arrays
class ThicknessFluxDivOnCell {
 public:
   bool Enabled;

   /// constructor declaration
   ThicknessFluxDivOnCell(const HorzMesh *Mesh);

   /// The functor takes cell index, vertical chunk index, and thickness flux
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 ICell, I4 KChunk,
                                   const Array2DReal &ThicknessFlux,
                                   const Array2DReal &NormalVelEdge) const {

      const I4 KStart        = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            DivTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                            ThicknessFlux(JEdge, K) * NormalVelEdge(JEdge, K) *
                            InvAreaCell;
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
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
};

/// Horizontal advection of potential vorticity defined on edges, for
/// momentum equation
class PotentialVortHAdvOnEdge {
 public:
   bool Enabled;

   /// constructor declaration
   PotentialVortHAdvOnEdge(const HorzMesh *Mesh);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// normalized relative vorticity, normalized planetary vorticity, layer
   /// thickness on edges, and normal velocity on edges as inputs,
   /// outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &NormRVortEdge,
                                   const Array2DReal &NormFEdge,
                                   const Array2DReal &FluxLayerThickEdge,
                                   const Array2DReal &NormVelEdge) const {

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
                             FluxLayerThickEdge(JEdge, K) *
                             NormVelEdge(JEdge, K) * NormVort;
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
   Array2DReal WeightsOnEdge;
};

/// Gradient of kinetic energy defined on edges, for momentum equation
class KEGradOnEdge {
 public:
   bool Enabled;

   /// constructor declaration
   KEGradOnEdge(const HorzMesh *Mesh);

   /// The functor takes edge index, vertical chunk index, and kinetic energy
   /// array as inputs, outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &KECell) const {

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
   Array1DReal DcEdge;
};

/// Gradient of sea surface height defined on edges multipled by gravitational
/// acceleration, for momentum equation
class SSHGradOnEdge {
 public:
   bool Enabled;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh);

   /// The functor takes edge index, vertical chunk index, and array of
   /// layer thickness/SSH, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &SshCell) const {

      const I4 KStart      = KChunk * VecLength;
      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) -=
             Grav * (SshCell(ICell1, K) - SshCell(ICell0, K)) * InvDcEdge;
      }
   }

 private:
   Real Grav = 9.80665_Real;
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
};

/// Laplacian horizontal mixing, for momentum equation
class VelocityDiffusionOnEdge {
 public:
   bool Enabled;

   Real ViscDel2;

   /// constructor declaration
   VelocityDiffusionOnEdge(const HorzMesh *Mesh);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// divergence of horizontal velocity (defined at cell centers) and relative
   /// vorticity (defined at vertices), outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &DivCell,
                                   const Array2DReal &RVortVertex) const {

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
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel2;
   Array2DReal EdgeMask;
};

/// Biharmonic horizontal mixing, for momentum equation
class VelocityHyperDiffOnEdge {
 public:
   bool Enabled;

   Real ViscDel4;

   /// Constructor declaration
   VelocityHyperDiffOnEdge(const HorzMesh *Mesh);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// the laplacian of divergence of horizontal velocity and the laplacian of
   /// the relative vorticity, outputs tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &Del2DivCell,
                                   const Array2DReal &Del2RVortVertex) const {

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
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal MeshScalingDel4;
   Array2DReal EdgeMask;
};

// Tracer horizontal advection term
class TracerHorzAdvOnCell {
 public:
   bool Enabled;

   TracerHorzAdvOnCell(const HorzMesh *Mesh);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk, const Array2DReal &NormVelEdge,
                                   const Array3DReal &HTracersOnEdge) const {

      const I4 KStart        = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HAdvTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            HAdvTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                             HTracersOnEdge(L, JEdge, K) *
                             NormVelEdge(JEdge, K) * InvAreaCell;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) -= HAdvTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
};

// Tracer horizontal diffusion term
class TracerDiffOnCell {
 public:
   bool Enabled;

   Real EddyDiff2;

   TracerDiffOnCell(const HorzMesh *Mesh);

   KOKKOS_FUNCTION void
   operator()(const Array3DReal &Tend, I4 L, I4 ICell, I4 KChunk,
              const Array3DReal &TracerCell,
              const Array2DReal &MeanLayerThickEdge) const {

      const I4 KStart        = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DiffTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel2(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            const Real TracerGrad =
                (TracerCell(L, JCell1, K) - TracerCell(L, JCell0, K));

            DiffTmp[KVec] -= EdgeSignOnCell(ICell, J) * RTemp *
                             MeanLayerThickEdge(JEdge, K) * TracerGrad;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) += EddyDiff2 * DiffTmp[KVec] * InvAreaCell;
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel2;
};

// Tracer biharmonic horizontal mixing term
class TracerHyperDiffOnCell {
 public:
   bool Enabled;

   Real EddyDiff4;

   TracerHyperDiffOnCell(const HorzMesh *Mesh);

   KOKKOS_FUNCTION void operator()(const Array3DReal &Tend, I4 L, I4 ICell,
                                   I4 KChunk,
                                   const Array3DReal &TrDel2Cell) const {

      const I4 KStart        = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real HypTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const I4 JEdge = EdgesOnCell(ICell, J);

         const I4 JCell0 = CellsOnEdge(JEdge, 0);
         const I4 JCell1 = CellsOnEdge(JEdge, 1);

         const Real RTemp =
             MeshScalingDel4(JEdge) * DvEdge(JEdge) / DcEdge(JEdge);

         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const I4 K = KStart + KVec;
            const Real Del2TrGrad =
                (TrDel2Cell(L, JCell1, K) - TrDel2Cell(L, JCell0, K));

            HypTmp[KVec] -= EdgeSignOnCell(ICell, J) * RTemp * Del2TrGrad;
         }
      }
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(L, ICell, K) -= EddyDiff4 * HypTmp[KVec] * InvAreaCell;
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnEdge;
   Array2DReal EdgeSignOnCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AreaCell;
   Array1DReal MeshScalingDel4;
};

} // namespace OMEGA
#endif
