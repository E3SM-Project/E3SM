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
   LinearEos,  /// Linear equation of state
   Teos10Eos,  /// Roquet et al. 2015 75 term expansion
   ConstantEos /// Constant specific volume equation of state
};

/// TEOS10 75-term Polynomial Equation of State
class Teos10Eos {
 public:
   /// constructor declaration
   Teos10Eos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   specific volume according to the Roquet et al. 2015 75 term expansion.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   I4 KDisp) const {

      Real SpecVolPCoeffs[6 * VecLength];
      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         /// Calculate the local specific volume polynomial pressure
         /// coefficients with cell center values
         calcPCoeffs(SpecVolPCoeffs, KVec, ConservTemp(ICell, K),
                     AbsSalinity(ICell, K));

         /// Calculate the specific volume at the given pressure
         /// If KDisp is 0, we use the current pressure, otherwise
         /// we use the displaced pressure (K + KDisp)
         /// Note: KDisp is only used for TEOS-10, for Linear EOS it
         /// is always 0.
         if (KDisp == 0) {
            // No displacement
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, K) * Pa2Db) +
                calcDelta(SpecVolPCoeffs, KVec, Pressure(ICell, K) * Pa2Db);
         } else {
            // Displacement, use the displaced pressure
            I4 KTmp = Kokkos::min(K + KDisp, MaxLayerCell(ICell));
            KTmp    = Kokkos::max(MinLayerCell(ICell), KTmp);
            SpecVol(ICell, K) =
                calcRefProfile(Pressure(ICell, KTmp) * Pa2Db) +
                calcDelta(SpecVolPCoeffs, KVec, Pressure(ICell, KTmp) * Pa2Db);
         }
      }
   }

   /// TEOS-10 helpers
   /// Calculate pressure polynomial coefficients for TEOS-10
   KOKKOS_FUNCTION void calcPCoeffs(Real (&SpecVolPCoeffs)[6 * VecLength],
                                    const I4 KVec, const Real Ct,
                                    const Real Sa) const {
      constexpr Real SaNorm = 40.0 * 35.16504 / 35.0;
      constexpr Real CtNorm = 40.0;
      constexpr Real DeltaS = 24.0;
      Real Ss               = Kokkos::sqrt((Sa + DeltaS) / SaNorm);
      Real Tt               = Ct / CtNorm;

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

      SpecVolPCoeffs[5 + 6 * KVec] = V005;
      SpecVolPCoeffs[4 + 6 * KVec] = V014 * Tt + V104 * Ss + V004;
      SpecVolPCoeffs[3 + 6 * KVec] =
          (V023 * Tt + V113 * Ss + V013) * Tt + (V203 * Ss + V103) * Ss + V003;
      SpecVolPCoeffs[2 + 6 * KVec] =
          (((V042 * Tt + V132 * Ss + V032) * Tt + (V222 * Ss + V122) * Ss +
            V022) *
               Tt +
           ((V312 * Ss + V212) * Ss + V112) * Ss + V012) *
              Tt +
          (((V402 * Ss + V302) * Ss + V202) * Ss + V102) * Ss + V002;
      SpecVolPCoeffs[1 + 6 * KVec] =
          ((((V051 * Tt + V141 * Ss + V041) * Tt + (V231 * Ss + V131) * Ss +
             V031) *
                Tt +
            ((V321 * Ss + V221) * Ss + V121) * Ss + V021) *
               Tt +
           (((V411 * Ss + V311) * Ss + V211) * Ss + V111) * Ss + V011) *
              Tt +
          ((((V501 * Ss + V401) * Ss + V301) * Ss + V201) * Ss + V101) * Ss +
          V001;
      SpecVolPCoeffs[0 + 6 * KVec] =
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
   }

   /// Evaluate pressure polynomial delta for TEOS-10
   KOKKOS_FUNCTION Real calcDelta(const Real (&SpecVolPCoeffs)[6 * VecLength],
                                  const I4 KVec, const Real P) const {

      constexpr Real PNorm = 1e-4;
      Real Pp              = P * PNorm;

      Real Delta = ((((SpecVolPCoeffs[5 + 6 * KVec] * Pp +
                       SpecVolPCoeffs[4 + 6 * KVec]) *
                          Pp +
                      SpecVolPCoeffs[3 + 6 * KVec]) *
                         Pp +
                     SpecVolPCoeffs[2 + 6 * KVec]) *
                        Pp +
                    SpecVolPCoeffs[1 + 6 * KVec]) *
                       Pp +
                   SpecVolPCoeffs[0 + 6 * KVec];

      return Delta;
   }

   /// Calculate reference profile for TEOS-10
   KOKKOS_FUNCTION Real calcRefProfile(Real P) const {
      constexpr Real PNorm = 1e-4;
      constexpr Real V00   = -4.4015007269e-05;
      constexpr Real V01   = 6.9232335784e-06;
      constexpr Real V02   = -7.5004675975e-07;
      constexpr Real V03   = 1.7009109288e-08;
      constexpr Real V04   = -1.6884162004e-08;
      constexpr Real V05   = 1.9613503930e-09;
      Real Pp              = P * PNorm;

      Real V0 =
          (((((V05 * Pp + V04) * Pp + V03) * Pp + V02) * Pp + V01) * Pp + V00) *
          Pp;
      return V0;
   }

   /// Calculate 2nd derivative of Gibbs wrt pot temp at ref P for TEOS-10
   KOKKOS_FUNCTION Real calcGibbsDerivPt0Pt0(Real Sa, Real P) const {
      Real x2 = Sfac * Sa;
      Real x  = Kokkos::sqrt(x2);
      Real y  = P * 0.025;

      Real g03 =
          -24715.571866078 +
          y * (4420.4472249096725 +
               y * (-1778.231237203896 +
                    y * (1160.5182516851419 +
                         y * (-569.531539542516 + y * 128.13429152494615))));

      Real g08 =
          x2 *
          (1760.062705994408 +
           x * (-86.1329351956084 +
                x * (-137.1145018408982 +
                     y * (296.20061691375236 +
                          y * (-205.67709290374563 + 49.9394019139016 * y))) +
                y * (-60.136422517125 + y * 10.50720794170734)) +
           y * (-1351.605895580406 +
                y * (1097.1125373015109 +
                     y * (-433.20648175062206 + 63.905091254154904 * y))));

      return ((g03 + g08) * 0.000625);
   }

   /// Calculate Pot Temmperature from Conservative Temp
   KOKKOS_FUNCTION Real calcPtFromCt(Real Sa, Real Ct) const {
      constexpr Real a0 = -1.446013646344788e-2;
      constexpr Real a1 = -3.305308995852924e-3;
      constexpr Real a2 = 1.062415929128982e-4;
      constexpr Real a3 = 9.477566673794488e-1;
      constexpr Real a4 = 2.166591947736613e-3;
      constexpr Real a5 = 3.828842955039902e-3;
      constexpr Real b0 = 1.000000000000000e0;
      constexpr Real b1 = 6.506097115635800e-4;
      constexpr Real b2 = 3.830289486850898e-3;
      constexpr Real b3 = 1.247811760368034e-6;
      Real a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
      Real pt, pt_old, ptm, dpt_dct, s1;

      s1 = Sa / Psu2Gpkg;

      a5ct = a5 * Ct;
      b3ct = b3 * Ct;

      ct_factor = (a3 + a4 * s1 + a5ct);
      pt_num    = a0 + s1 * (a1 + a2 * s1) + Ct * ct_factor;
      pt_recden = 1.0 / (b0 + b1 * s1 + Ct * (b2 + b3ct));
      pt        = pt_num * pt_recden;

      dpt_dct = pt_recden * (ct_factor + a5ct - (b2 + b3ct + b3ct) * pt);

      ct_diff = calcCtFromPt(Sa, pt) - Ct;
      pt_old  = pt;
      pt      = pt_old - ct_diff * dpt_dct;
      ptm     = 0.5 * (pt + pt_old);

      dpt_dct = -Cp0Sw / ((ptm + TkFrz) * calcGibbsDerivPt0Pt0(Sa, ptm));

      pt      = pt_old - ct_diff * dpt_dct;
      ct_diff = calcCtFromPt(Sa, pt) - Ct;
      pt_old  = pt;
      return (pt_old - ct_diff * dpt_dct);
   }

   /// Calculate Conservative Temmperature from Potential Temp
   KOKKOS_FUNCTION Real calcCtFromPt(Real Sa, Real Pt) const {
      Real x2, x, y, pot_enthalpy;

      x2 = Sfac * Sa;
      x  = Kokkos::sqrt(x2);
      y  = Pt * 0.025e0; /*! normalize for F03 and F08 */
      pot_enthalpy =
          61.01362420681071e0 +
          y * (168776.46138048015e0 +
               y * (-2735.2785605119625e0 +
                    y * (2574.2164453821433e0 +
                         y * (-1536.6644434977543e0 +
                              y * (545.7340497931629e0 +
                                   (-50.91091728474331e0 -
                                    18.30489878927802e0 * y) *
                                       y))))) +
          x2 *
              (268.5520265845071e0 +
               y * (-12019.028203559312e0 +
                    y * (3734.858026725145e0 +
                         y * (-2046.7671145057618e0 +
                              y * (465.28655623826234e0 +
                                   (-0.6370820302376359e0 -
                                    10.650848542359153e0 * y) *
                                       y)))) +
               x * (937.2099110620707e0 +
                    y * (588.1802812170108e0 + y * (248.39476522971285e0 +
                                                    (-3.871557904936333e0 -
                                                     2.6268019854268356e0 * y) *
                                                        y)) +
                    x * (-1687.914374187449e0 +
                         x * (246.9598888781377e0 +
                              x * (123.59576582457964e0 -
                                   48.5891069025409e0 * x)) +
                         y * (936.3206544460336e0 +
                              y * (-942.7827304544439e0 +
                                   y * (369.4389437509002e0 +
                                        (-33.83664947895248e0 -
                                         9.987880382780322e0 * y) *
                                            y))))));

      return (pot_enthalpy / Cp0Sw);
   }

   /// Calculates freezing Conservative Temperature using TEOS-10 polynomial
   /// (polynomial error in [-5e-4, 6e-4] K, from GSW package)
   KOKKOS_FUNCTION Real calcCtFreezing(const Real Sa, const Real P,
                                       const Real SaturationFract) const {
      constexpr Real Sso = 35.16504;
      constexpr Real C0  = 0.017947064327968736;
      constexpr Real C1  = -6.076099099929818;
      constexpr Real C2  = 4.883198653547851;
      constexpr Real C3  = -11.88081601230542;
      constexpr Real C4  = 13.34658511480257;
      constexpr Real C5  = -8.722761043208607;
      constexpr Real C6  = 2.082038908808201;
      constexpr Real C7  = -7.389420998107497;
      constexpr Real C8  = -2.110913185058476;
      constexpr Real C9  = 0.2295491578006229;
      constexpr Real C10 = -0.9891538123307282;
      constexpr Real C11 = -0.08987150128406496;
      constexpr Real C12 = 0.3831132432071728;
      constexpr Real C13 = 1.054318231187074;
      constexpr Real C14 = 1.065556599652796;
      constexpr Real C15 = -0.7997496801694032;
      constexpr Real C16 = 0.3850133554097069;
      constexpr Real C17 = -2.078616693017569;
      constexpr Real C18 = 0.8756340772729538;
      constexpr Real C19 = -2.079022768390933;
      constexpr Real C20 = 1.596435439942262;
      constexpr Real C21 = 0.1338002171109174;
      constexpr Real C22 = 1.242891021876471;

      // Note: a = 0.502500117621 / Sso
      constexpr Real A = 0.014289763856964;
      constexpr Real B = 0.057000649899720;

      Real Sar = Sa * 1.0e-2;
      Real X   = Kokkos::sqrt(Sar);
      Real Pr  = P * 1.0e-4;

      Real CtFreez =
          C0 + Sar * (C1 + X * (C2 + X * (C3 + X * (C4 + X * (C5 + C6 * X))))) +
          Pr * (C7 + Pr * (C8 + C9 * Pr)) +
          Sar * Pr *
              (C10 + Pr * (C12 + Pr * (C15 + C21 * Sar)) +
               Sar * (C13 + C17 * Pr + C19 * Sar) +
               X * (C11 + Pr * (C14 + C18 * Pr) +
                    Sar * (C16 + C20 * Pr + C22 * Sar)));

      /* Adjust for the effects of dissolved air */
      CtFreez = CtFreez - SaturationFract * (1e-3) * (2.4 - A * Sa) *
                              (1.0 + B * (1.0 - Sa / Sso));

      return CtFreez;
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Linear Equation of State
class LinearEos {
 public:
   /// Coefficients for LinearEos (overwritten by config file if set there)
   Real DRhodT  = -0.2;  ///< Thermal expansion coefficient (kg m^-3 degC^-1)
   Real DRhodS  = 0.8;   ///< Haline contraction coefficient (kg m^-3)
   Real RhoT0S0 = RhoFw; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// constructor declaration
   LinearEos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and the ocean tracers (conservative)
   //   temperature, and (absolute) salinity as inputs, and outputs the
   //   linear specific volume.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         SpecVol(ICell, K) =
             1.0_Real / (RhoT0S0 + (DRhodT * ConservTemp(ICell, K) +
                                    DRhodS * AbsSalinity(ICell, K)));
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Constant Equation of State
class ConstantEos {
 public:
   /// constructor declaration
   ConstantEos(const VertCoord *VCoord);

   //   The functor takes the full arrays of specific volume (inout),
   //   the indices ICell and KChunk, and returns a constant specific volume
   //   value for all active layers.
   KOKKOS_FUNCTION void operator()(Array2DReal SpecVol, I4 ICell, I4 KChunk,
                                   const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));
      (void)ConservTemp;
      (void)AbsSalinity;

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K        = KStart + KVec;
         SpecVol(ICell, K) = 1.0_Real / RhoSw;
      }
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Functor for calculating the squared Brunt-Vaisala frequency using TEOS-10
class Teos10BruntVaisalaFreqSq {
 public:
   /// Constructor for BruntVaisalaFreqSq
   Teos10BruntVaisalaFreqSq(const VertCoord *VCoord);

   //   The functor takes the full arrays of squared Brunt-Vaisala frequency
   //   (inout) the index ICell, and the ocean tracers (conservative)
   //   temperature, (absolute) salinity, pressure, specific volume as inputs,
   //   and outputs the squared Brunt-Vaisala frequency.
   KOKKOS_FUNCTION void operator()(Array2DReal BruntVaisalaFreqSq, I4 ICell,
                                   I4 KChunk, const Array2DReal &ConservTemp,
                                   const Array2DReal &AbsSalinity,
                                   const Array2DReal &Pressure,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         if (K == 0) {
            // No Brunt-Vaisala frequency at surface
            BruntVaisalaFreqSq(ICell, K) = 0.0_Real;
         } else {
            // Calculate squared Brunt-Vaisala frequency
            Real CtInt =
                0.5_Real * (ConservTemp(ICell, K) + ConservTemp(ICell, K - 1));
            Real SaInt =
                0.5_Real * (AbsSalinity(ICell, K) + AbsSalinity(ICell, K - 1));
            Real PInt =
                0.5_Real * (Pressure(ICell, K) + Pressure(ICell, K - 1));
            Real SpInt = 0.5_Real * (SpecVol(ICell, K) + SpecVol(ICell, K - 1));
            Real AlphaInt = calcAlpha(SaInt, CtInt, PInt * Pa2Db, SpInt);
            Real BetaInt  = calcBeta(SaInt, CtInt, PInt * Pa2Db, SpInt);
            Real DSa      = AbsSalinity(ICell, K) - AbsSalinity(ICell, K - 1);
            Real DCt      = ConservTemp(ICell, K) - ConservTemp(ICell, K - 1);
            Real DP       = Pressure(ICell, K) - Pressure(ICell, K - 1);

            BruntVaisalaFreqSq(ICell, K) = Gravity * Gravity *
                                           (BetaInt * DSa - AlphaInt * DCt) /
                                           (SpInt * DP);
         }
      }
   }

   /// Calculate alpha values for the squared Brunt-Vaisala frequency
   KOKKOS_FUNCTION Real calcAlpha(Real Sa, Real Ct, Real P, Real Sp) const {
      constexpr Real Factor = 0.0248826675584615;
      constexpr Real Offset = 5.971840214030754e-1;
      constexpr Real PNorm  = 1.0e-4;
      Real Ss               = Kokkos::sqrt(Factor * Sa + Offset);
      Real Tt               = 0.025_Real * Ct;
      Real Pp               = P * PNorm;

      constexpr Real A000 = -1.56497346750e-5;
      constexpr Real A001 = 1.85057654290e-5;
      constexpr Real A002 = -1.17363867310e-6;
      constexpr Real A003 = -3.65270065530e-7;
      constexpr Real A004 = 3.14540999020e-7;
      constexpr Real A010 = 5.55242129680e-5;
      constexpr Real A011 = -2.34332137060e-5;
      constexpr Real A012 = 4.26100574800e-6;
      constexpr Real A013 = 5.73918103180e-7;
      constexpr Real A020 = -4.95634777770e-5;
      constexpr Real A021 = 2.37838968519e-5;
      constexpr Real A022 = -1.38397620111e-6;
      constexpr Real A030 = 2.76445290808e-5;
      constexpr Real A031 = -1.36408749928e-5;
      constexpr Real A032 = -2.53411666056e-7;
      constexpr Real A040 = -4.02698077700e-6;
      constexpr Real A041 = 2.53683834070e-6;
      constexpr Real A050 = 1.23258565608e-6;
      constexpr Real A100 = 3.50095997640e-5;
      constexpr Real A101 = -9.56770881560e-6;
      constexpr Real A102 = -5.56991545570e-6;
      constexpr Real A103 = -2.72956962370e-7;
      constexpr Real A110 = -7.48716846880e-5;
      constexpr Real A111 = -4.73566167220e-7;
      constexpr Real A112 = 7.82747741600e-7;
      constexpr Real A120 = 7.24244384490e-5;
      constexpr Real A121 = -1.03676320965e-5;
      constexpr Real A122 = 2.32856664276e-8;
      constexpr Real A130 = -3.50383492616e-5;
      constexpr Real A131 = 5.18268711320e-6;
      constexpr Real A140 = -1.65263794500e-6;
      constexpr Real A200 = -4.35926785610e-5;
      constexpr Real A201 = 1.11008347650e-5;
      constexpr Real A202 = 5.46207488340e-6;
      constexpr Real A210 = 7.18156455200e-5;
      constexpr Real A211 = 5.85666925900e-6;
      constexpr Real A212 = -1.31462208134e-6;
      constexpr Real A220 = -4.30608991440e-5;
      constexpr Real A221 = 9.49659182340e-7;
      constexpr Real A230 = 1.74814722392e-5;
      constexpr Real A300 = 3.45324618280e-5;
      constexpr Real A301 = -9.84471178440e-6;
      constexpr Real A302 = -1.35441856270e-6;
      constexpr Real A310 = -3.73971683740e-5;
      constexpr Real A311 = -9.76522784000e-7;
      constexpr Real A320 = 6.85899736680e-6;
      constexpr Real A400 = -1.19594097880e-5;
      constexpr Real A401 = 2.59092252600e-6;
      constexpr Real A410 = 7.71906784880e-6;
      constexpr Real A500 = 1.38645945810e-6;

      Real Rval =
          A000 +
          Ss * (A100 + Ss * (A200 + Ss * (A300 + Ss * (A400 + A500 * Ss)))) +
          Tt * (A010 + Ss * (A110 + Ss * (A210 + Ss * (A310 + A410 * Ss))) +
                Tt * (A020 + Ss * (A120 + Ss * (A220 + A320 * Ss)) +
                      Tt * (A030 + Ss * (A130 + A230 * Ss) +
                            Tt * (A040 + A140 * Ss + A050 * Tt)))) +
          Pp * (A001 + Ss * (A101 + Ss * (A201 + Ss * (A301 + A401 * Ss))) +
                Tt * (A011 + Ss * (A111 + Ss * (A211 + A311 * Ss)) +
                      Tt * (A021 + Ss * (A121 + A221 * Ss) +
                            Tt * (A031 + A131 * Ss + A041 * Tt))) +
                Pp * (A002 + Ss * (A102 + Ss * (A202 + A302 * Ss)) +
                      Tt * (A012 + Ss * (A112 + A212 * Ss) +
                            Tt * (A022 + A122 * Ss + A032 * Tt)) +
                      Pp * (A003 + A103 * Ss + A013 * Tt + A004 * Pp)));

      return 0.025_Real * Rval / Sp;
   }

   /// Calculate beta values for the squared Brunt-Vaisala frequency
   KOKKOS_FUNCTION Real calcBeta(Real Sa, Real Ct, Real P, Real Sp) const {
      constexpr Real Factor = 0.0248826675584615;
      constexpr Real Offset = 5.971840214030754e-1;
      constexpr Real PNorm  = 1.0e-4;
      Real Ss               = Kokkos::sqrt(Factor * Sa + Offset);
      Real Tt               = 0.025_Real * Ct;
      Real Pp               = P * PNorm;

      constexpr Real B000 = -3.10389819760e-4;
      constexpr Real B003 = 3.63101885150e-7;
      constexpr Real B004 = -1.11471254230e-7;
      constexpr Real B010 = 3.50095997640e-5;
      constexpr Real B013 = -2.72956962370e-7;
      constexpr Real B020 = -3.74358423440e-5;
      constexpr Real B030 = 2.41414794830e-5;
      constexpr Real B040 = -8.75958731540e-6;
      constexpr Real B050 = -3.30527589000e-7;
      constexpr Real B100 = 1.33856134076e-3;
      constexpr Real B103 = 3.34926075600e-8;
      constexpr Real B110 = -8.71853571220e-5;
      constexpr Real B120 = 7.18156455200e-5;
      constexpr Real B130 = -2.87072660960e-5;
      constexpr Real B140 = 8.74073611960e-6;
      constexpr Real B200 = -2.55143801811e-3;
      constexpr Real B210 = 1.03597385484e-4;
      constexpr Real B220 = -5.60957525610e-5;
      constexpr Real B230 = 6.85899736680e-6;
      constexpr Real B300 = 2.32344279772e-3;
      constexpr Real B310 = -4.78376391520e-5;
      constexpr Real B320 = 1.54381356976e-5;
      constexpr Real B400 = -1.05461852535e-3;
      constexpr Real B410 = 6.93229729050e-6;
      constexpr Real B500 = 1.91594743830e-4;
      constexpr Real B001 = 2.42624687470e-5;
      constexpr Real B011 = -9.56770881560e-6;
      constexpr Real B021 = -2.36783083610e-7;
      constexpr Real B031 = -3.45587736550e-6;
      constexpr Real B041 = 1.29567177830e-6;
      constexpr Real B101 = -6.95849219480e-5;
      constexpr Real B111 = 2.22016695300e-5;
      constexpr Real B121 = 5.85666925900e-6;
      constexpr Real B131 = 6.33106121560e-7;
      constexpr Real B201 = 1.12412331915e-4;
      constexpr Real B211 = -2.95341353532e-5;
      constexpr Real B221 = -1.46478417600e-6;
      constexpr Real B301 = -6.92888744480e-5;
      constexpr Real B311 = 1.03636901040e-5;
      constexpr Real B401 = 1.54637136265e-5;
      constexpr Real B002 = -5.84844329840e-7;
      constexpr Real B012 = -5.56991545570e-6;
      constexpr Real B022 = 3.91373870800e-7;
      constexpr Real B032 = 7.76188880920e-9;
      constexpr Real B102 = -9.62445031940e-6;
      constexpr Real B112 = 1.09241497668e-5;
      constexpr Real B122 = -1.31462208134e-6;
      constexpr Real B202 = 1.47789320994e-5;
      constexpr Real B212 = -4.06325568810e-6;
      constexpr Real B302 = -7.12478989080e-6;

      Real Rval =
          B000 +
          Ss * (B100 + Ss * (B200 + Ss * (B300 + Ss * (B400 + B500 * Ss)))) +
          Tt * (B010 + Ss * (B110 + Ss * (B210 + Ss * (B310 + B410 * Ss))) +
                Tt * (B020 + Ss * (B120 + Ss * (B220 + B320 * Ss)) +
                      Tt * (B030 + Ss * (B130 + B230 * Ss) +
                            Tt * (B040 + B140 * Ss + B050 * Tt)))) +
          Pp * (B001 + Ss * (B101 + Ss * (B201 + Ss * (B301 + B401 * Ss))) +
                Tt * (B011 + Ss * (B111 + Ss * (B211 + B311 * Ss)) +
                      Tt * (B021 + Ss * (B121 + B221 * Ss) +
                            Tt * (B031 + B131 * Ss + B041 * Tt))) +
                Pp * (B002 + Ss * (B102 + Ss * (B202 + B302 * Ss)) +
                      Tt * (B012 + Ss * (B112 + B212 * Ss) +
                            Tt * (B022 + B122 * Ss + B032 * Tt)) +
                      Pp * (B003 + B103 * Ss + B013 * Tt + B004 * Pp)));

      return -0.5_Real * Rval * Factor / (Sp * Ss);
   }

 private:
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Linear squared Brunt-Vaisala frequency calculator
class LinearBruntVaisalaFreqSq {
 public:
   /// Coefficients from LinearEos (overwritten by config file if set there)
   Real RhoT0S0 = 1000.0; ///< Reference density (kg m^-3) at (T,S)=(0,0)

   /// constructor declaration
   LinearBruntVaisalaFreqSq(const VertCoord *VCoord);

   //   The functor takes the full arrays of squared Brunt-Vaisala frequency
   //   (inout), the index ICell, and the specific volume and pseudo-thickness
   //   as inputs, and outputs the squared Brunt-Vaisala frequency.
   KOKKOS_FUNCTION void operator()(Array2DReal BruntVaisalaFreqSq, I4 ICell,
                                   I4 KChunk,
                                   const Array2DReal &SpecVol) const {

      const I4 KStart = chunkStart(KChunk, MinLayerCell(ICell));
      const I4 KLen   = chunkLength(KChunk, KStart, MaxLayerCell(ICell));

      for (int KVec = 0; KVec < KLen; ++KVec) {
         const I4 K = KStart + KVec;
         if (K == 0) {
            /// No Brunt-Vaisala frequency at the top level
            BruntVaisalaFreqSq(ICell, K) = 0.0_Real;
         } else {
            /// Calculate squared Brunt-Vaisala frequency at mid-point between
            /// K-1 and K Do not need to use displaced specific volume here
            /// since only the linear EOS is used with this BVF formulation.
            BruntVaisalaFreqSq(ICell, K) =
                -(Gravity / RhoT0S0) *
                ((1.0_Real / SpecVol(ICell, K - 1)) -
                 (1.0_Real / SpecVol(ICell, K))) /
                (GeomZMid(ICell, K - 1) - GeomZMid(ICell, K));
         }
      }
   }

 private:
   Array2DReal GeomZMid;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
};

/// Class for Equation of State (EOS) calculations
class Eos {
 public:
   /// Get instance of Eos
   static Eos *getInstance();

   /// Destroy instance (frees Kokkos views)
   static void destroyInstance();

   EosType EosChoice;              ///< Current EOS type in use
   Array2DReal SpecVol;            ///< Specific volume field at level centers
   Array2DReal SpecVolDisplaced;   ///< Displaced specific volume field
   Array2DReal BruntVaisalaFreqSq; ///< Squared Brunt-Vaisala frequency field

   std::string SpecVolFldName; ///< Field name for specific volume
   std::string
       SpecVolDisplacedFldName; ///< Field name for displaced specific volume
   std::string BruntVaisalaFreqSqFldName; ///< Field name for squared
                                          ///< Brunt-Vaisala frequency
   std::string EosGroupName;              ///< EOS group name (for config)
   std::string Name;                      ///< Name of this EOS instance

   /// Compute specific volume for all cells/layers
   void computeSpecVol(const Array2DReal &ConservTemp,
                       const Array2DReal &AbsSalinity,
                       const Array2DReal &Pressure);

   /// Compute displaced specific volume (for vertical displacement)
   void computeSpecVolDisp(const Array2DReal &ConservTemp,
                           const Array2DReal &AbsSalinity,
                           const Array2DReal &Pressure, I4 KDisp);

   /// Compute squared Brunt-Vaisala frequency for all cells/layers
   void computeBruntVaisalaFreqSq(const Array2DReal &ConservTemp,
                                  const Array2DReal &AbsSalinity,
                                  const Array2DReal &Pressure,
                                  const Array2DReal &SpecVol);

   /// Convert Conservative Temperature to potential temperature
   Real calcPtFromCt(const Real &Sa, const Real &Ct) const;

   /// Convert potential temperature to Conservative Temperature
   Real calcCtFromPt(const Real &Sa, const Real &Pt) const;

   /// Initialize EOS from config and mesh
   static void init();

 private:
   /// Private constructor
   Eos(const std::string &Name, const HorzMesh *Mesh, const VertCoord *VCoord);

   /// Private destructor
   ~Eos();

   /// Instance pointer
   static Eos *Instance;

   /// Delete copy and move constructors and assignment operators
   Eos(const Eos &)            = delete;
   Eos &operator=(const Eos &) = delete;
   Eos(Eos &&)                 = delete;
   Eos &operator=(Eos &&)      = delete;

   const HorzMesh *Mesh;    ///< Horizontal mesh
   const VertCoord *VCoord; ///< Vertical coordinate

   Teos10Eos ComputeSpecVolTeos10;     ///< TEOS-10 specific volume calculator
   LinearEos ComputeSpecVolLinear;     ///< Linear specific volume calculator
   ConstantEos ComputeSpecVolConstant; ///< Constant specific volume calculator
   Teos10BruntVaisalaFreqSq
       ComputeBruntVaisalaFreqSqTeos10; ///< TEOS-10 squared Brunt-Vaisala
                                        ///< calculator
   LinearBruntVaisalaFreqSq
       ComputeBruntVaisalaFreqSqLinear; ///< Linear squared Brunt-Vaisala
                                        ///< calculator

   // Define fields and metadata
   void defineFields();

}; // End class Eos

} // namespace OMEGA
#endif
