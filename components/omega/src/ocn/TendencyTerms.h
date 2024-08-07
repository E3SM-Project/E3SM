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

namespace OMEGA {

/// Divergence of thickness flux at cell centers, for updating layer thickness
/// arrays
class ThicknessFluxDivOnCell {
 public:
   bool Enabled = false;

   /// constructor declaration
   ThicknessFluxDivOnCell(const HorzMesh *Mesh, Config *Options);

   /// The functor takes cell index, vertical chunk index, and thickness flux
   /// array as inputs, outputs the tendency array
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

/// Horizontal advection of potential vorticity defined on edges, for
/// momentum equation
class PotentialVortHAdvOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   PotentialVortHAdvOnEdge(const HorzMesh *Mesh, Config *Options);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// normalized relative vorticity, normalized planetary vorticity, layer
   /// thickness on edges, and normal velocity on edges as inputs,
   /// outputs the tendency array
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DR8 &NormRVortEdge,
                                   const Array2DR8 &NormFEdge,
                                   const Array2DR8 &FluxLayerThickEdge,
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
   Array2DR8 WeightsOnEdge;
};

/// Gradient of kinetic energy defined on edges, for momentum equation
class KEGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   KEGradOnEdge(const HorzMesh *Mesh, Config *Options);

   /// The functor takes edge index, vertical chunk index, and kinetic energy
   /// array as inputs, outputs the tendency array
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

/// Gradient of sea surface height defined on edges multipled by gravitational
/// acceleration, for momentum equation
class SSHGradOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   SSHGradOnEdge(const HorzMesh *Mesh, Config *Options);

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
   R8 Grav = 9.80665_Real;
   Array2DI4 CellsOnEdge;
   Array1DR8 DcEdge;
};

/// Laplacian horizontal mixing, for momentum equation
class VelocityDiffusionOnEdge {
 public:
   bool Enabled = false;

   /// constructor declaration
   VelocityDiffusionOnEdge(const HorzMesh *Mesh, Config *Options);

   /// The functor takes edge index, vertical chunk index, and arrays for
   /// divergence of horizontal velocity (defined at cell centers) and relative
   /// vorticity (defined at vertices), outputs tendency array
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
   R8 ViscDel2 = 1._Real;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 MeshScalingDel2;
   Array2DR8 EdgeMask;
};

/// Biharmonic horizontal mixing, for momentum equation
class VelocityHyperDiffOnEdge {
 public:
   bool Enabled = false;

   /// Constructor declaration
   VelocityHyperDiffOnEdge(const HorzMesh *Mesh, Config *Options);

   /// The functor takes the edge index, vertical chunk index, and arrays for
   /// the laplacian of divergence of horizontal velocity and the laplacian of
   /// the relative vorticity, outputs tendency array
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
   Real ViscDel4 = 1._Real;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array1DR8 DcEdge;
   Array1DR8 DvEdge;
   Array1DR8 MeshScalingDel4;
   Array2DR8 EdgeMask;
};

/// A class that can be used to calculate the thickness and
/// velocity tendencies within the timestepping algorithm.
class Tendencies {
 public:
   // Arrays for accumulating tendencies
   Array2DReal LayerThicknessTend;
   Array2DReal NormalVelocityTend;

   // Instances of tendency terms
   ThicknessFluxDivOnCell ThicknessFluxDiv;
   PotentialVortHAdvOnEdge PotientialVortHAdv;
   KEGradOnEdge KEGrad;
   SSHGradOnEdge SSHGrad;
   VelocityDiffusionOnEdge VelocityDiffusion;
   VelocityHyperDiffOnEdge VelocityHyperDiff;

   // Methods to compute tendency groups
   // TODO Add AuxilaryState as calling argument
   void computeThicknessTendencies(const OceanState *State, const AuxiliaryState *AuxState,
                                   int TimeLevel, Real Time);
   void computeVelocityTendencies(const OceanState *State, const AuxiliaryState *AuxState,
                                  int TimeLevel, Real Time);
   void computeAllTendencies(const OceanState *State, const AuxiliaryState *AuxState,
                             int TimeLevel, Real Time);

   // Create a non-default tendencies
   static Tendencies *
   create(const std::string &Name, ///< [in] Name for tendencies
          const HorzMesh *Mesh,    ///< [in] Horizontal mesh
          int NVertLevels,         ///< [in] Number of vertical levels
          Config *Options          ///< [in] Configuration options
   );

   // Destructor
   ~Tendencies();

   // Initialize Omega tendencies
   static int init();

   // Deallocates arrays
   static void clear();

   // Remove tendencies object by name
   static void erase(const std::string &Name ///< [in]
   );

   // get default tendencies
   static Tendencies *getDefault();

   // get tendencies by name
   static Tendencies *get(const std::string &Name ///< [in]
   );

 private:
   // Mesh sizes
   I4 NCellsOwned; ///< Number of cells owned by this task
   I4 NEdgesOwned; ///< Number of edges owned by this task
   I4 NChunks;     ///< Number of vertical level chunks

   // Construct a new tendency object
   Tendencies(const std::string &Name, ///< [in] Name for tendencies
              const HorzMesh *Mesh,    ///< [in] Horizontal mesh
              int NVertLevels,         ///< [in] Number of vertical levels
              Config *Options          ///< [in] Configuration options
   );

   // Pointer to default tendencies
   static Tendencies *DefaultTendencies;

   // Map of all tendency objects
   static std::map<std::string, std::unique_ptr<Tendencies>> AllTendencies;

}; // end class Tendencies

} // namespace OMEGA
#endif
