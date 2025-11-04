#ifndef OMEGA_PGRAD_H
#define OMEGA_PGRAD_H
//===-- ocn/PGrad.h - Pressure Gradient -----------------*- C++ -*-===//
///
/// Implements the PressureGrad class which provides a centered and
/// high-order pressure gradient option and dispatches computations to
/// functor objects. This follows the patterns used in Eos.h/Eos.cpp.
//
//===----------------------------------------------------------------------===//

#include "Config.h"
#include "Eos.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"
#include <memory>

namespace OMEGA {

enum class PressureGradType { Centered, HighOrder1, HighOrder2 };

/// Centered pressure gradient functor
class PressureGradCentered {
 public:
   bool Enabled;

   /// constructor declaration
   PressureGradCentered(const HorzMesh *Mesh);

   // Compute centered pressure gradient contribution for given edge and
   // vertical chunk. This appends results into the Tend array (in-place).
   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &PressureMid,
                                   const Array2DReal &PressureInterface,
                                   const Array2DReal &Geopotential,
                                   const Array2DReal &LayerThick,
                                   const Array2DReal &ZInterface,
                                   const Array2DReal &SpecVol,
                                   const Array2DReal &SpecVolInterface) const {
      const I4 KStart = KChunk * VecLength;

      const I4 ICell0      = CellsOnEdge(IEdge, 0);
      const I4 ICell1      = CellsOnEdge(IEdge, 1);
      const Real InvDcEdge = 1.0_Real / DcEdge(IEdge);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;

         // Average quantities to edge
         const Real PAlphaEdgeK =
             0.5_Real *
             (PressureInterface(ICell1, K) * SpecVolInterface(ICell1, K) +
              PressureInterface(ICell0, K) * SpecVolInterface(ICell0, K));
         const Real PAlphaEdgeKP1 =
             0.5_Real * (PressureInterface(ICell1, K + 1) *
                             SpecVolInterface(ICell1, K + 1) +
                         PressureInterface(ICell0, K + 1) *
                             SpecVolInterface(ICell0, K + 1));
         const Real InvLayerThickEdge =
             2.0_Real / (LayerThick(ICell1, K) + LayerThick(ICell0, K));

         Real GeoTerm =
             (Geopotential(ICell1, K) - Geopotential(ICell0, K)) * InvDcEdge;
         Real PresTerm = (LayerThick(ICell1, K) * SpecVol(ICell1, K) *
                              PressureMid(ICell1, K) -
                          LayerThick(ICell0, K) * SpecVol(ICell0, K) *
                              PressureMid(ICell0, K)) *
                         InvDcEdge * InvLayerThickEdge;
         Real InterfaceTerm =
             (PAlphaEdgeK * (ZInterface(ICell1, K) - ZInterface(ICell0, K)) -
              PAlphaEdgeKP1 *
                  (ZInterface(ICell1, K + 1) - ZInterface(ICell0, K + 1))) *
             InvDcEdge * InvLayerThickEdge;

         Tend(IEdge, K) +=
             EdgeMask(IEdge, K) * (-GeoTerm - PresTerm + InterfaceTerm);
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
};

/// High-order pressure gradient functor (placeholder)
class PressureGradHighOrder {
 public:
   bool Enabled;

   /// constructor declaration
   PressureGradHighOrder(const HorzMesh *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &Tend, I4 IEdge, I4 KChunk,
                                   const Array2DReal &Pressure,
                                   const Array2DReal &Geopotential,
                                   const Array2DReal &SpecVol) const {
      // Placeholder: for now, no-op (future high-order implementation)
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         Tend(IEdge, K) += 0.0_Real;
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array2DReal EdgeMask;
};

/// Pressure gradient manager class
class PressureGrad {
 public:
   /// Initialize the default instance
   static void init();

   /// Deallocates arrays and deletes instance
   static void clear();

   /// Remove pressure gradient object by name
   static void erase(const std::string &Name ///< [in]
   );

   /// Get the default instance
   static PressureGrad *getDefault();

   /// Get instance by name
   static PressureGrad *get(const std::string &Name ///< [in]
   );

   // Destructor
   ~PressureGrad();

   /// Compute pressure gradient tendencies and add into Tend array
   void computePressureGrad(Array2DReal Tend, const OceanState *State,
                            const VertCoord *VCoord, const Eos *EqState,
                            const int TimeLevel);

   static PressureGrad *create(const std::string &Name, const HorzMesh *Mesh,
                               const VertCoord *VCoord, Config *Options);

 private:
   // Construct a new pressure gradient object
   PressureGrad(const HorzMesh *Mesh, const VertCoord *VCoord, Config *Options);

   // forbid copy and move construction
   PressureGrad(const PressureGrad &) = delete;
   PressureGrad(PressureGrad &&)      = delete;

   // Pointer to default pressure gradient object
   static PressureGrad *DefaultPGrad;

   // Mesh-related sizes
   I4 NEdgesAll   = 0;
   I4 NChunks     = 0;
   I4 NVertLayers = 0;

   // Data required for computation (stored copies of mesh/VCoord arrays)
   Array2DI4 CellsOnEdge;      ///< cells surrounding each edge
   Array1DReal DvEdge;         ///< distance between cell centers across edge
   Array2DReal EdgeSignOnCell; ///< orientation of edge relative to cell

   // Instances of functors
   PressureGradCentered CenteredPGrad;
   PressureGradHighOrder HighOrderPGrad;

   // Choice from config
   PressureGradType PressureGradChoice = PressureGradType::Centered;

   // Map of all pressure gradient objects by name
   static std::map<std::string, std::unique_ptr<PressureGrad>> AllPGrads;

}; // end class PressureGrad

} // namespace OMEGA
#endif
