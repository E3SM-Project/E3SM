#ifndef OMEGA_VERTMIX_H
#define OMEGA_VERTMIX_H
//===-- ocn/VertMix.cpp - Vertical Mixing Coefficients -----------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating vertical diffusivity and viscosity
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the vertical diffusivity and viscosity
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include <string>

namespace OMEGA {

enum class VertMixType {
   PP, /// Pacanowski and Philander (1981)
   KPP /// K-Profile Parameterization (Large et al., 1994)
};

class ConvectiveMix {
   public:
      bool Enabled = true; ///< Enable convective mixing flag

      // Convective mixing parameters
      Real ConvDiff       = 1.0; ///< Convective vertical viscosity and diffusivity (m^2 s^-1)
      Real ConvTriggerBVF = 0.0; ///< Reference density (kg m^-3) at (T,S)=(0,0)

      /// Constructor for ConvectiveMix
      ConvectiveMix();

      KOKKOS_FUNCTION void operator()(Array2DReal VertDiff, Array2DReal VertVisc,
                                      I4 ICell, I4 K,
                                      const Array2DReal &BruntVaisalaFreq) const {

         if (BruntVaisalaFreq(ICell, K) < ConvTriggerBVF) {
            VertDiff(ICell, K) = VertDiff(ICell, K) + ConvDiff;
            VertVisc(ICell, K) = VertVisc(ICell, K) + ConvDiff;
         }
      }
};

class PPShearMix {
   public:
      bool Enabled = true; ///< Enable shear mixing flag

      // Shear mixing parameters
      Real ShearNuZero   = 0.005; ///< Numerator of Pacanowski and Philander (1981) Eq (1).
      Real ShearAlpha    = 5.0;   ///< Alpha value used in Pacanowski and Philander (1981) Eqs (1) and (2).
      Real ShearExponent = 2.0;   /// Exponent value used in Pacanowski and Philander (1981) Eqs (1).

      /// Constructor for PPShearMix
      PPShearMix(const HorzMesh *Mesh);

      KOKKOS_FUNCTION void operator()(Array2DReal VertDiff, Array2DReal VertVisc,
                                      I4 ICell, I4 K,
                                      const Array2DReal &NormalVelocity,
                                      const Array2DReal &TangentialVelocity,
                                      const Array2DReal &BruntVaisalaFreq) const {
         Real ShearViscVal = 0.0;

         Real InvAreaCell = 1.0 / AreaCell(ICell);
         Real ShearSquared = 0.0;
         for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
            I4 JEdge = EdgesOnCell(ICell, J);
            Real Factor = 0.5 * DcEdge(JEdge) * DvEdge(JEdge) * InvAreaCell;
            Real DelU2 = 
               Kokkos::pow(NormalVelocity(JEdge,K-1) - NormalVelocity(JEdge,K), 2.0) + 
               Kokkos::pow(TangentialVelocity(JEdge,K-1) - TangentialVelocity(JEdge,K), 2.0);
            ShearSquared = ShearSquared + Factor * DelU2;
         }

         Real RichardsonNum = BruntVaisalaFreq(ICell, K) / 
                              Kokkos::max(1.0e-12_Real, ShearSquared);

         ShearViscVal = ShearNuZero / 
            Kokkos::pow(1.0_Real + ShearAlpha * RichardsonNum, ShearExponent);
         VertVisc(ICell, K) = VertVisc(ICell, K) + ShearViscVal;
         VertDiff(ICell, K) = VertDiff(ICell, K) + ShearViscVal / (1.0_Real + ShearAlpha * RichardsonNum);
      }

   private:
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   const int NVertLevels;
};

//class KPP {
//   public:

      /// Constructor for KPP
//      KPP();

//      KOKKOS_FUNCTION void operator()(Array2DReal VertDiff, Array2DReal VertVisc,
//                                      I4 ICell, I4 K,
//                                      const Array2DReal &NormalVelocity,
//                                      const Array2DReal &TangentialVelocity,
//                                      const Array2DReal &BruntVaisalaFreq) const {

//         VertDiff(ICell, K) = 0.0;
//         VertVisc(ICell, K) = 0.0;
//      }
//};

/// Class for Vertical Mixing Coefficient (VertMix) calculations
class VertMix {
 public:
   /// Get instance of VertMix
   static VertMix *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   VertMixType VertMixChoice;        ///< Current VertMix type in use
   Array2DReal VertDiff;             ///< Vertical diffusivity field (m^2 s^-1)
   Array2DReal VertVisc;             ///< Vertical viscosity field (m^2 s^-1)

   std::string VertDiffFldName;          ///< Field name for vertical diffusivity
   std::string VertViscFldName;          ///< Field name for vertical viscosity
   std::string VertMixGroupName;         ///< VertMix group name (for config)
   std::string Name;                     ///< Name of this VertMix instance

   // Background mixing parameters
   Real BackDiff = 1.0e-5; ///< Background vertical diffusivity (m^2 s^-1)
   Real BackVisc = 1.0e-4; ///< Background vertical viscosity (m^2 s^-1)

   ConvectiveMix ComputeVertMixConv;   ///< Functor for PP VertMix calculation
   PPShearMix ComputeVertMixShear;   ///< Functor for PP VertMix calculation

   /// Compute vertical diffusivity and viscosity for all cells/levels
   void computeVertMix(const Array2DReal &NormalVelocity,
                       const Array2DReal &TangentialVelocity,
                       const Array2DReal &BruntVaisalaFreq);

   /// Initialize VertMix from config and mesh
   static void init();

 private:
   /// Private constructor
    VertMix(const std::string &Name, const HorzMesh *Mesh, int NVertLevels);
    
   /// Private destructor
    ~VertMix();

    static VertMix* Instance; ///< Instance pointer

    // Delete copy and move constructors and assignment operators
    VertMix(const VertMix&) = delete;
    VertMix& operator=(const VertMix&) = delete;
    VertMix(VertMix&&) = delete;
    VertMix& operator=(VertMix&&) = delete;

   I4 NCellsAll;         ///< Number of horizontal cells
   I4 NVertLevels;       ///< Number of vertical levels

   // Define fields and metadata
   void defineFields();

}; // End class VertMix

} // namespace OMEGA
#endif