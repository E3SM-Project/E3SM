#ifndef OMEGA_EOS_H
#define OMEGA_EOS_H
//===-- ocn/Eos.h - Equation of State --------------------*- C++ -*-===//
//
/// \file
/// \brief Contains functors for calculating specific volume
///
/// This header defines functors to be called by the time-stepping scheme
/// to calculate the specific volume based on the choice of EOS
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

enum class EosType {
   LinearEos, /// Linear equation of state
   Teos10Eos  /// Roquet et al. 2015 75 term expansion
};

/// TEOS10 75-term Polynomial Equation of State
class Teos10Eos {
 public:
   Array2DReal SpecVolPCoeffs;

   /// constructor declaration
   Teos10Eos(int NVertLayers);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   specific volume according to the Roquet et al. 2015 75 term expansion.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   I4 KDisp) const {

      OMEGA_SCOPE(LocSpecVolPCoeffs, SpecVolPCoeffs);
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;

         /// Calculate the local specific volume polynomial pressure
         /// coefficients
         calcPCoeffs(LocSpecVolPCoeffs, KVec, ConservTemp(ICell, K),
                     AbsSalinity(ICell, K));

         /// Calculate the specific volume at the given pressure
         /// If KDisp is 0, we use the current pressure, otherwise we use the
         /// displaced pressure (K + KDisp)
         /// Note: KDisp is only used for TEOS-10, for Linear EOS it
         /// is always 0.
         if (KDisp == 0) {
            // No displacement
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, K)) +
                calcDelta(LocSpecVolPCoeffs, KVec, Pressure(ICell, K));
         } else {
            // Displacement, use the displaced pressure
            I4 KTmp = Kokkos::min(K + KDisp, NVertLayers - 1);
            KTmp    = Kokkos::max(0, KTmp);
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, KTmp)) +
                calcDelta(LocSpecVolPCoeffs, KVec, Pressure(ICell, KTmp));
         }
      }
   }

   /// TEOS-10 helpers
   /// Calculate pressure polynomial coefficients for TEOS-10
   KOKKOS_FUNCTION void calcPCoeffs(Array2DReal SpecVolPCoeffs, const I4 K,
                                    const Real Ct, const Real Sa) const {
      constexpr Real SAu    = 40.0 * 35.16504 / 35.0;
      constexpr Real CTu    = 40.0;
      constexpr Real DeltaS = 24.0;
      Real Ss               = Kokkos::sqrt((Sa + DeltaS) / SAu);
      Real Tt               = Ct / CTu;

      /// Coefficients for the polynomial expansion
      constexpr Real V000 = 1.0769995862e-03;
      constexpr Real V100 = -3.1038981976e-04;
      constexpr Real V200 = 6.6928067038e-04;
      constexpr Real V300 = -8.5047933937e-04;
      constexpr Real V400 = 5.8086069943e-04;
      constexpr Real V500 = -2.1092370507e-04;
      constexpr Real V600 = 3.1932457305e-05;
      constexpr Real V010 = -1.5649734675e-05;
      constexpr Real V110 = 3.5009599764e-05;
      constexpr Real V210 = -4.3592678561e-05;
      constexpr Real V310 = 3.4532461828e-05;
      constexpr Real V410 = -1.1959409788e-05;
      constexpr Real V510 = 1.3864594581e-06;
      constexpr Real V020 = 2.7762106484e-05;
      constexpr Real V120 = -3.7435842344e-05;
      constexpr Real V220 = 3.5907822760e-05;
      constexpr Real V320 = -1.8698584187e-05;
      constexpr Real V420 = 3.8595339244e-06;
      constexpr Real V030 = -1.6521159259e-05;
      constexpr Real V130 = 2.4141479483e-05;
      constexpr Real V230 = -1.4353633048e-05;
      constexpr Real V330 = 2.2863324556e-06;
      constexpr Real V040 = 6.9111322702e-06;
      constexpr Real V140 = -8.7595873154e-06;
      constexpr Real V240 = 4.3703680598e-06;
      constexpr Real V050 = -8.0539615540e-07;
      constexpr Real V150 = -3.3052758900e-07;
      constexpr Real V060 = 2.0543094268e-07;
      constexpr Real V001 = -1.6784136540e-05;
      constexpr Real V101 = 2.4262468747e-05;
      constexpr Real V201 = -3.4792460974e-05;
      constexpr Real V301 = 3.7470777305e-05;
      constexpr Real V401 = -1.7322218612e-05;
      constexpr Real V501 = 3.0927427253e-06;
      constexpr Real V011 = 1.8505765429e-05;
      constexpr Real V111 = -9.5677088156e-06;
      constexpr Real V211 = 1.1100834765e-05;
      constexpr Real V311 = -9.8447117844e-06;
      constexpr Real V411 = 2.5909225260e-06;
      constexpr Real V021 = -1.1716606853e-05;
      constexpr Real V121 = -2.3678308361e-07;
      constexpr Real V221 = 2.9283346295e-06;
      constexpr Real V321 = -4.8826139200e-07;
      constexpr Real V031 = 7.9279656173e-06;
      constexpr Real V131 = -3.4558773655e-06;
      constexpr Real V231 = 3.1655306078e-07;
      constexpr Real V041 = -3.4102187482e-06;
      constexpr Real V141 = 1.2956717783e-06;
      constexpr Real V051 = 5.0736766814e-07;
      constexpr Real V002 = 3.0623833435e-06;
      constexpr Real V102 = -5.8484432984e-07;
      constexpr Real V202 = -4.8122251597e-06;
      constexpr Real V302 = 4.9263106998e-06;
      constexpr Real V402 = -1.7811974727e-06;
      constexpr Real V012 = -1.1736386731e-06;
      constexpr Real V112 = -5.5699154557e-06;
      constexpr Real V212 = 5.4620748834e-06;
      constexpr Real V312 = -1.3544185627e-06;
      constexpr Real V022 = 2.1305028740e-06;
      constexpr Real V122 = 3.9137387080e-07;
      constexpr Real V222 = -6.5731104067e-07;
      constexpr Real V032 = -4.6132540037e-07;
      constexpr Real V132 = 7.7618888092e-09;
      constexpr Real V042 = -6.3352916514e-08;
      constexpr Real V003 = -3.8088938393e-07;
      constexpr Real V103 = 3.6310188515e-07;
      constexpr Real V203 = 1.6746303780e-08;
      constexpr Real V013 = -3.6527006553e-07;
      constexpr Real V113 = -2.7295696237e-07;
      constexpr Real V023 = 2.8695905159e-07;
      constexpr Real V004 = 8.8302421514e-08;
      constexpr Real V104 = -1.1147125423e-07;
      constexpr Real V014 = 3.1454099902e-07;
      constexpr Real V005 = 4.2369007180e-09;

      SpecVolPCoeffs(5, K) = V005;
      SpecVolPCoeffs(4, K) = V014 * Tt + V104 * Ss + V004;
      SpecVolPCoeffs(3, K) =
          (V023 * Tt + V113 * Ss + V013) * Tt + (V203 * Ss + V103) * Ss + V003;
      SpecVolPCoeffs(2, K) =
          (((V042 * Tt + V132 * Ss + V032) * Tt + (V222 * Ss + V122) * Ss +
            V022) *
               Tt +
           ((V312 * Ss + V212) * Ss + V112) * Ss + V012) *
              Tt +
          (((V402 * Ss + V302) * Ss + V202) * Ss + V102) * Ss + V002;
      SpecVolPCoeffs(1, K) =
          ((((V051 * Tt + V141 * Ss + V041) * Tt + (V231 * Ss + V131) * Ss +
             V031) *
                Tt +
            ((V321 * Ss + V221) * Ss + V121) * Ss + V021) *
               Tt +
           (((V411 * Ss + V311) * Ss + V211) * Ss + V111) * Ss + V011) *
              Tt +
          ((((V501 * Ss + V401) * Ss + V301) * Ss + V201) * Ss + V101) * Ss +
          V001;
      SpecVolPCoeffs(0, K) =
          (((((V060 * Tt + V150 * Ss + V050) * Tt + (V240 * Ss + V140) * Ss +
              V040) *
                 Tt +
             ((V330 * Ss + V230) * Ss + V130) * Ss + V030) *
                Tt +
            (((V420 * Ss + V320) * Ss + V220) * Ss + V120) * Ss + V020) *
               Tt +
           ((((V510 * Ss + V410) * Ss + V310) * Ss + V210) * Ss + V110) * Ss +
           V010) *
              Tt +
          (((((V600 * Ss + V500) * Ss + V400) * Ss + V300) * Ss + V200) * Ss +
           V100) *
              Ss +
          V000;

      // could insert a check here (abs(value)> 0 value or <e+33)
   }

   /// Evaluate pressure polynomial delta for TEOS-10
   KOKKOS_FUNCTION Real calcDelta(const Array2DReal &SpecVolPCoeffs, const I4 K,
                                  const Real P) const {

      Real Pp = P * Pa2Db;

      Real Delta = ((((SpecVolPCoeffs(5, K) * Pp + SpecVolPCoeffs(4, K)) * Pp +
                      SpecVolPCoeffs(3, K)) *
                         Pp +
                     SpecVolPCoeffs(2, K)) *
                        Pp +
                    SpecVolPCoeffs(1, K)) *
                       Pp +
                   SpecVolPCoeffs(0, K);
      return Delta;
   }

   /// Calculate reference profile for TEOS-10
   KOKKOS_FUNCTION Real calcRefProfile(Real P) const {
      constexpr Real V00 = -4.4015007269e-05;
      constexpr Real V01 = 6.9232335784e-06;
      constexpr Real V02 = -7.5004675975e-07;
      constexpr Real V03 = 1.7009109288e-08;
      constexpr Real V04 = -1.6884162004e-08;
      constexpr Real V05 = 1.9613503930e-09;
      Real Pp            = P * Pa2Db;

      Real V0 =
          (((((V05 * Pp + V04) * Pp + V03) * Pp + V02) * Pp + V01) * Pp + V00) *
          Pp;
      return V0;
   }

 private:
   const int NVertLayers;
};

/// Linear Equation of State
class LinearEos {
 public:
   /// Coefficients for LinearEos (overwritten by config file if set there)
   Real DRhodT  = -0.2;  ///< Thermal expansion coefficient (kg m^-3 degC^-1)
   Real DRhodS  = 0.8;   ///< Haline contraction coefficient (kg m^-3)
   Real RhoT0S0 = RhoFw; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// constructor declaration
   LinearEos();

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   linear specific volume.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity) const {
      const I4 KStart = KChunk * VecLength;
      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const I4 K = KStart + KVec;
         SpecVol(ICell, K) =
             1.0_Real / (RhoT0S0 + (DRhodT * ConservTemp(ICell, K) +
                                    DRhodS * AbsSalinity(ICell, K)));
      }
   }
};

/// Class for Equation of State (EOS) calculations
class Eos {
 public:
   /// Get instance of Eos
   static Eos *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   EosType EosChoice;            ///< Current EOS type in use
   Array2DReal SpecVol;          ///< Specific volume field
   Array2DReal SpecVolDisplaced; ///< Displaced specific volume field

   std::string SpecVolFldName; ///< Field name for specific volume
   std::string
       SpecVolDisplacedFldName; ///< Field name for displaced specific volume
   std::string EosGroupName;    ///< EOS group name (for config)
   std::string Name;            ///< Name of this EOS instance

   /// Compute specific volume for all cells/layers
   void computeSpecVol(const Array2DReal &ConservTemp,
                       const Array2DReal &AbsSalinity,
                       const Array2DReal &Pressure);

   /// Compute displaced specific volume (for vertical displacement)
   void computeSpecVolDisp(const Array2DReal &ConservTemp,
                           const Array2DReal &AbsSalinity,
                           const Array2DReal &Pressure, I4 KDisp);

   /// Initialize EOS from config and mesh
   static void init();

 private:
   /// Private constructor
   Eos(const std::string &Name, const HorzMesh *Mesh, int NVertLayers);

   /// Private destructor
   ~Eos();

   /// Instance pointer
   static Eos *Instance;

   /// Delete copy and move constructors and assignment operators
   Eos(const Eos &)            = delete;
   Eos &operator=(const Eos &) = delete;
   Eos(Eos &&)                 = delete;
   Eos &operator=(Eos &&)      = delete;

   I4 NCellsAll; ///< Number of horizontal cells
   I4 NChunks;   ///< Number of vertical chunks (for vectorization)

   Teos10Eos ComputeSpecVolTeos10; ///< TEOS-10 specific volume calculator
   LinearEos ComputeSpecVolLinear; ///< Linear specific volume calculator

   // Define fields and metadata
   void defineFields();

}; // End class Eos

} // namespace OMEGA
#endif
