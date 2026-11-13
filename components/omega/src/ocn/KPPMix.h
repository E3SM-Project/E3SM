#ifndef OMEGA_KPP_MIX_H
#define OMEGA_KPP_MIX_H
//===-- ocn/KPPMix.h - K-Profile Parameterization --------*- C++ -*-===//
//
/// \file
/// \brief K-Profile Parameterization (KPP) boundary layer mixing scheme
///
/// This header defines the KPPMix class for computing ocean boundary layer
/// mixing coefficients using the K-Profile Parameterization scheme.
/// Follows Large et al. (1994) formulation with optional Langmuir turbulence
/// enhancement.
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "KPPConstants.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "SurfaceFlux.h"
#include "TimeMgr.h"
#include "VertCoord.h"
#include <string>

namespace OMEGA {

class Forcing;
class OceanState;

/// @brief enum sets which matching criterion is used at the OSBL base.
/// SimpleShapes does no matching, MatchBoth matches KPP diagnosed
/// diffusivity to the interior mixing coefficient.
/// Also selects the shape used for the non-local flux, which follows the
/// scalar diffusivity profile.
enum class KPPMatchType : I4 {
   SimpleShapes = 0, ///< Unmatched Large et al. (1994) cubic shapes
   MatchBoth    = 1  ///< Match the interior coefficient at the OSBL base
};

/// @brief Langmuir turbulence enhancement factor applied to buoyancy forcing
class KPPLangmuirFactor {
 public:
   bool UseLangmuirTurbulence = true; ///< Apply wave enhancement
   /// Disable Langmuir above this ice fraction for theory wave model
   /// This should be disabled for active wave configurations
   Real IceFracThresholdForLangmuir = KPP::IceFracThresh;

   KOKKOS_FUNCTION void operator()(const Array1DReal &LangmuirFactor, I4 ICell,
                                   const Array1DReal &IceFraction,
                                   const Array1DReal &SurfaceFrictionVelocity,
                                   const Array1DReal &WindSpeed10m) const {

      const Real IceFrac = IceFraction(ICell);
      if (UseLangmuirTurbulence && IceFrac < IceFracThresholdForLangmuir) {
         const Real UStar = SurfaceFrictionVelocity(ICell);
         const Real Wind10m =
             (WindSpeed10m.extent(0) > 0) ? WindSpeed10m(ICell) : 0.0_Real;
         LangmuirFactor(ICell) =
             KPP::computeLangmuirEnhancement(Wind10m, UStar, 50.0);
      } else {
         LangmuirFactor(ICell) = 1.0_Real;
      }
   }
};

class PotentialDensityOnCell {
 public:
   KOKKOS_FUNCTION void
   computeRefPressure(const Array2DReal &RefPressure, I4 ICell, I4 K,
                      const Array1DReal &SurfacePressure) const {
      RefPressure(ICell, K) = SurfacePressure(ICell);
   }

   KOKKOS_FUNCTION void operator()(const Array2DReal &PotentialDensity,
                                   I4 ICell, I4 K,
                                   const Array2DReal &SpecVolDisplaced) const {
      PotentialDensity(ICell, K) =
          1.0_Real / Kokkos::max(1.0e-12_Real, SpecVolDisplaced(ICell, K));
   }
};

class KPPSurfaceForcingOnCell {
 public:
   bool UseTracerForcing = false;

   KPPSurfaceForcingOnCell(const HorzMesh *Mesh, const VertCoord *VCoord,
                           const Eos *EosInst);

   KOKKOS_FUNCTION void operator()(
       const Array1DReal &FrictionVelocity, const Array1DReal &BuoyancyFlux,
       I4 ICell, const Array2DReal &ConservTemp, const Array2DReal &AbsSalinity,
       const Array2DReal &PressureMid, const Array2DReal &SpecVol,
       const Array1DReal &ZonalStress, const Array1DReal &MeridStress,
       const Array1DReal &LatentHeatFluxEvap,
       const Array1DReal &SensibleHeatFlux,
       const Array1DReal &LongWaveHeatFluxUp,
       const Array1DReal &LongWaveHeatFluxDown,
       const Array1DReal &SeaIceHeatFlux, const Array1DReal &ShortWaveHeatFlux,
       const Array1DReal &SnowFlux, const Array1DReal &RainFlux,
       const Array1DReal &EvaporationFlux,
       const Array1DReal &SeaIceFreshWaterFlux,
       const Array1DReal &IceRunoffFlux, const Array1DReal &RiverRunoffFlux,
       const Array1DReal &SeaIceSaltFlux) const {
      const I4 KTop = MinLayerCell(ICell);
      if (KTop < 0 || KTop >= NVertLayers)
         return;

      const Real TauMag = Kokkos::sqrt(ZonalStress(ICell) * ZonalStress(ICell) +
                                       MeridStress(ICell) * MeridStress(ICell));
      FrictionVelocity(ICell) =
          Kokkos::sqrt(Kokkos::max(0.0_Real, TauMag / RhoSw));
      BuoyancyFlux(ICell) = 0.0_Real;
      if (!UseTracerForcing)
         return;

      const Real SaTop    = AbsSalinity(ICell, KTop);
      const Real CtTop    = ConservTemp(ICell, KTop);
      const Real PTopDb   = PressureMid(ICell, KTop) * Pa2Db;
      const Real HeatFlux = sfcHeatFluxWithoutMassEnthalpy(
          ICell, SaTop, PTopDb, EosChoice, LongWaveHeatFluxUp,
          LongWaveHeatFluxDown, ShortWaveHeatFlux, SensibleHeatFlux,
          SeaIceHeatFlux, SeaIceFreshWaterFlux, SeaIceSaltFlux,
          LatentHeatFluxEvap, SnowFlux, IceRunoffFlux);
      const Real FreshWaterFlux = sfcFreshWaterFlux(
          ICell, SnowFlux, RainFlux, EvaporationFlux, SeaIceFreshWaterFlux,
          IceRunoffFlux, RiverRunoffFlux);
      const Real TempFlux = HeatFlux * HFluxFac;
      const Real SaltFlux =
          SeaIceSaltFlux(ICell) * SFluxFac - FreshWaterFlux * SaTop / RhoSw;
      const Real SpVol = Kokkos::max(1.0e-12_Real, SpecVol(ICell, KTop));
      // LinearDRhodT/LinearDRhodS are ignored unless EosChoice is LinearEos.
      const Real Alpha = Eos::computeAlpha(EosChoice, SaTop, CtTop, PTopDb,
                                           SpVol, LinearDRhodT);
      const Real Beta = Eos::computeBeta(EosChoice, SaTop, CtTop, PTopDb, SpVol,
                                         LinearDRhodS);
      BuoyancyFlux(ICell) = Gravity * (Alpha * TempFlux - Beta * SaltFlux);
   }

 private:
   Real LinearDRhodT;
   Real LinearDRhodS;
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   EosType EosChoice;
};

/// @brief Stage 1 kernel: Compute OSBL depth from the bulk Richardson number
/// criterion
///
/// The OSBL depth is determined by searching for the depth at which the bulk
/// Richardson number exceeds the critical value.  In KPP the bulk Richardson
/// Number is augmented by a turbulent velocity which is designed to
/// guarantee the entrainment heat flux is -0.2 times the surface buoyancy flux.
///
/// Velocities live on edges (C-grid), so the shear entering Ri is a
/// kite-area weighted average over the edges of each cell.
class KPPOSBLDepthSearch {
 public:
   ///< Ri_crit for determining OSBL base
   Real CriticalRichardson = KPP::CriticalRi;
   Real SurfaceLayerExtent = KPP::SurfaceLayerExtent; ///< Frac of OSBL depth
   /// Apply minimum OSBL depth above this ice fraction
   Real IceFracThresholdForMinimumOSBL = KPP::IceSuppressThresh;
   /// Min OSBL depth under sea ice (m)
   Real MinimumOSBLUnderSeaIce = KPP::MinOSBLUnderIce;
   /// Only search past the base of the OSBL if KPP diagnostics are requested
   bool FullRiProfile = false;

   /// Constructor for KPPOSBLDepthSearch
   KPPOSBLDepthSearch(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(
       const Array1DReal &OSBLDepth, const Array1DI4 &OSBLDepthIndex,
       const Array2DReal &BulkRichardsonNumber,
       const Array2DReal &BulkRichardsonShear,
       const Array2DReal &UnresolvedShear, const Array2DReal &BuoyancyJump,
       I4 ICell, const Array2DReal &PotentialDensity,
       const Array2DReal &NormalVelocity, const Array2DReal &TangentialVelocity,
       const Array1DReal &SurfaceFrictionVelocity,
       const Array1DReal &SurfaceBuoyancyFlux,
       const Array2DReal &BruntVaisalaFreqSq, const Array1DReal &IceFraction,
       const Array1DReal &LangmuirFactor) const {

      const Real UStar = SurfaceFrictionVelocity(ICell);
      // The buoyancy flux defined here is the direct flux only, it does not
      // include any enthalpy from mass fluxes.
      const Real BuoyFlux = SurfaceBuoyancyFlux(ICell);

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);
      if (KMin < 0 || KMin >= NVertLayers || KMax < KMin) {
         return;
      }
      const I4 KIntTop  = KMin + 1;
      const I4 KIntDeep = KMax + 1;

      const Real IceFrac = IceFraction(ICell);

      // KPP depths are measured below the free surface, so geometric
      // heights must be offset by the sea surface height.
      const Real Ssh = SshCell(ICell);

      // Default to the full water column; overwritten if Ri crosses.
      Real OBLDepth         = Ssh - ZInterface(ICell, KIntDeep);
      I4 KCross             = -1;
      const Real RiCritical = CriticalRichardson;
      // Bulk Ri is evaluated at cell centers while the reference average spans
      // the top epsilon*d; this factor corrects for that offset.
      const Real RiScaling   = 1.0_Real - 0.5_Real * SurfaceLayerExtent;
      const Real BuoyFluxEff = BuoyFlux * LangmuirFactor(ICell);

      // Unresolved shear coefficient, Large et al. (1994) Eq. (23):
      // Vt^2 = Cv * sqrt(-beta_T/(c_s*eps)) / (kappa^2 * Ri_crit) * d*N*w_s
      // CSUnres is c_s for the strongly-unstable scalar branch and VtCoef
      // collects the constant prefactor.
      const Real CSUnres = 24.0_Real * Kokkos::sqrt(17.0_Real);
      const Real VtCoef =
          Kokkos::sqrt(0.2_Real / Kokkos::fmax(KPP::NumericalTolerance,
                                               CSUnres * SurfaceLayerExtent)) /
          (VonKar * VonKar);

      // ----------------------------------------------------------------
      // Edge weights are the MPAS triangle areas formed between the
      // cell center, and vertices, with area 0.25*dc*dv, normalized by
      // the cell area.
      // ----------------------------------------------------------------
      const I4 NEdges                 = NEdgesOnCell(ICell);
      bool EdgeValid[MaxEdgesBound]   = {};
      Real EdgeWeights[MaxEdgesBound] = {};
      // NEdgesOnCell is always >= 3 for a valid MPAS mesh cell.
      const Real InvAreaCell = 1.0_Real / AreaCell(ICell);
      for (I4 J = 0; J < NEdges; ++J) {
         const I4 IEdge = EdgesOnCell(ICell, J);
         const I4 KEMin = MinLayerEdgeBot(IEdge);
         const I4 KEMax = MaxLayerEdgeTop(IEdge);
         EdgeValid[J]   = (KEMax >= KEMin && KEMin >= 0 && KEMin < NVertLayers);
         if (EdgeValid[J]) {
            EdgeWeights[J] =
                0.25_Real * DcEdge(IEdge) * DvEdge(IEdge) * InvAreaCell;
         }
      }

      // ----------------------------------------------------------------
      // Bulk Richardson search, Large et al. (1994) Eq. (21):
      //   Ri_b(d) = (B_r - B(d)) d / (|V_r - V(d)|^2 + Vt^2(d))
      // where the reference values B_r, V_r are averaged over the top
      // epsilon*d (d is the candidate OSBL depth) of the column. Since
      // epsilon*d grows monotonically with
      // the trial depth, the averaging window only ever extends downward,
      // so the running sums are carried across trial depths rather than
      // rebuilt from the surface at each one.
      // The OSBL base is the first d at which Ri_b reaches RiCritical, and
      // the search stops there unless the full profile is wanted for the
      // Ri diagnostics.
      // ----------------------------------------------------------------

      // Cell surface-layer density average
      I4 KSurfaceAvg = KMin;
      const Real ThickTop =
          Kokkos::abs(ZInterface(ICell, KMin + 1) - ZInterface(ICell, KMin));
      Real SumThickness = ThickTop;
      Real SumRho       = PotentialDensity(ICell, KMin) * SumThickness;

      // Per-edge surface-layer velocity averages
      // These are updated as the surface layer depth grows with increasing
      // trial depth
      I4 KSurfE[MaxEdgesBound]      = {};
      Real SumThickE[MaxEdgesBound] = {};
      Real SumUnE[MaxEdgesBound]    = {};
      Real SumVtE[MaxEdgesBound]    = {};

      for (I4 J = 0; J < NEdges; ++J) {
         if (!EdgeValid[J]) {
            continue;
         }
         const I4 IEdge  = EdgesOnCell(ICell, J);
         const I4 JCell  = CellsOnCell(ICell, J);
         const I4 KEMin  = MinLayerEdgeBot(IEdge);
         KSurfE[J]       = KEMin;
         const I4 KIntE0 = Kokkos::min(KEMin + 1, NVertLayers);
         const Real ZEdgeTop =
             0.5_Real * (ZInterface(ICell, KEMin) + ZInterface(JCell, KEMin));
         const Real ZEdgeBot =
             0.5_Real * (ZInterface(ICell, KIntE0) + ZInterface(JCell, KIntE0));
         const Real Thick0 = Kokkos::abs(ZEdgeBot - ZEdgeTop);
         SumThickE[J]      = Thick0;
         const I4 KE0      = Kokkos::min(KEMin, NVertLayers - 1);
         SumUnE[J]         = NormalVelocity(IEdge, KE0) * SumThickE[J];
         SumVtE[J]         = TangentialVelocity(IEdge, KE0) * SumThickE[J];
      }

      for (I4 K = KMin; K <= KMax; ++K) {
         const I4 KInt      = K + 1;
         const Real ZDepth  = Ssh - ZInterface(ICell, KInt);
         const Real ZCenter = Ssh - ZMid(ICell, K);
         if (ZDepth < KPP::NumericalTolerance)
            continue;

         const Real SurfLayerDepth = SurfaceLayerExtent * ZDepth;

         // Advance cell surface average for density
         while (KSurfaceAvg < K &&
                (Ssh - ZInterface(ICell, KSurfaceAvg + 1)) < SurfLayerDepth) {
            ++KSurfaceAvg;
            const I4 KSA      = Kokkos::min(KSurfaceAvg, NVertLayers - 1);
            const Real DZ     = Kokkos::abs(ZInterface(ICell, KSurfaceAvg + 1) -
                                            ZInterface(ICell, KSurfaceAvg));
            const Real ThickK = DZ;
            SumThickness += ThickK;
            SumRho += PotentialDensity(ICell, KSA) * ThickK;
         }

         // Advance per-edge surface averages for velocity
         for (I4 J = 0; J < NEdges; ++J) {
            if (!EdgeValid[J]) {
               continue;
            }
            const I4 IEdge     = EdgesOnCell(ICell, J);
            const I4 JCell     = CellsOnCell(ICell, J);
            const I4 KEMax     = MaxLayerEdgeTop(IEdge);
            const Real SshEdge = 0.5_Real * (SshCell(ICell) + SshCell(JCell));
            while (KSurfE[J] < K &&
                   (SshEdge - 0.5_Real * (ZInterface(ICell, KSurfE[J] + 1) +
                                          ZInterface(JCell, KSurfE[J] + 1))) <
                       SurfLayerDepth) {
               ++KSurfE[J];
               const I4 KE = Kokkos::min(
                   Kokkos::max(KSurfE[J], MinLayerEdgeBot(IEdge)), KEMax);
               const Real ZEdgeTop = 0.5_Real * (ZInterface(ICell, KSurfE[J]) +
                                                 ZInterface(JCell, KSurfE[J]));
               const Real ZEdgeBot =
                   0.5_Real * (ZInterface(ICell, KSurfE[J] + 1) +
                               ZInterface(JCell, KSurfE[J] + 1));
               const Real DZ     = Kokkos::abs(ZEdgeBot - ZEdgeTop);
               const Real ThickK = DZ;
               SumThickE[J] += ThickK;
               SumUnE[J] += NormalVelocity(IEdge, KE) * ThickK;
               SumVtE[J] += TangentialVelocity(IEdge, KE) * ThickK;
            }
         }

         // Compute the average density in the surface layer
         const Real InvSumThickness = 1.0_Real / SumThickness;
         const Real RhoAvgSurf      = SumRho * InvSumThickness;

         // Buoyancy jump B_r - B(d), positive for stable stratification
         const Real RhoK     = PotentialDensity(ICell, K);
         const Real DeltaRho = RhoK - RhoAvgSurf;
         const Real DeltaB   = Gravity * DeltaRho / RhoSw;
         // The following field is for diagnostic purposes only
         // and does not affect the computation
         BuoyancyJump(ICell, KInt) = DeltaB;

         // Resolved shear |V_r - V(d)|^2, averaged over the cell edges
         Real DeltaVSq      = 0.0_Real;
         Real EdgeWeightSum = 0.0_Real;
         for (I4 J = 0; J < NEdges; ++J) {
            if (!EdgeValid[J]) {
               continue;
            }
            const I4 IEdge       = EdgesOnCell(ICell, J);
            const I4 KEMin       = MinLayerEdgeBot(IEdge);
            const I4 KEMax       = MaxLayerEdgeTop(IEdge);
            const I4 KE          = Kokkos::min(Kokkos::max(K, KEMin), KEMax);
            const Real InvThickE = 1.0_Real / SumThickE[J];
            const Real UnAvg     = SumUnE[J] * InvThickE;
            const Real VtAvg     = SumVtE[J] * InvThickE;
            const Real UnK       = NormalVelocity(IEdge, KE);
            const Real VtK       = TangentialVelocity(IEdge, KE);
            const Real DUn       = UnK - UnAvg;
            const Real DVt       = VtK - VtAvg;
            const Real EdgeShear = DUn * DUn + DVt * DVt;
            DeltaVSq += EdgeWeights[J] * EdgeShear;
            EdgeWeightSum += EdgeWeights[J];
         }
         if (EdgeWeightSum > 0.0_Real) {
            DeltaVSq /= EdgeWeightSum;
         }
         BulkRichardsonShear(ICell, KInt) = DeltaVSq;

         // Turbulent scalar velocity scale w_s at the surface-layer depth
         Real WTurb = 0.0_Real;
         if (UStar > KPP::NumericalTolerance) {
            const Real U3   = UStar * UStar * UStar;
            const Real Zeta = SurfaceLayerExtent * ZDepth * VonKar *
                              BuoyFluxEff / Kokkos::fmax(U3, KPP::Tiny);
            // kppPhiInvScalar is positive over its whole domain.
            WTurb = VonKar * UStar * KPP::kppPhiInvScalar(Zeta);
         } else if (BuoyFluxEff < 0.0_Real) {
            // Free convection limit: u* drops out and w_s ~ (c_s d B_0)^1/3
            const Real CS = KPP::CMoS;
            const Real WS3 =
                -CS * SurfaceLayerExtent * ZDepth * VonKar * BuoyFluxEff;
            WTurb = VonKar * Kokkos::pow(Kokkos::max(WS3, 0.0_Real),
                                         1.0_Real / 3.0_Real);
         }

         // Unresolved turbulent shear Vt^2 (m^2/s^2), Large et al. Eq.
         // (23). Cv ramps from 2.1 to 1.7 as stratification strengthens.
         const Real NInt = Kokkos::sqrt(
             Kokkos::max(0.0_Real, BruntVaisalaFreqSq(ICell, KInt)));
         const Real Cv =
             (NInt < 0.002_Real) ? (2.1_Real - 200.0_Real * NInt) : 1.7_Real;
         const Real Vt2 = Kokkos::fmax(
             KPP::MinUnresolvedShearSq,
             // CriticalRichardson is clamped positive once in KPPMix::init().
             Cv * VtCoef * ZCenter * NInt * WTurb / RiCritical);
         UnresolvedShear(ICell, KInt) = Vt2;

         // Vt2's own floor (MinUnresolvedShearSq) already exceeds
         // NumericalTolerance, so VelScaleSq can't underflow the divide.
         const Real VelScaleSq = DeltaVSq + Vt2;

         const Real RiBulk = RiScaling * DeltaB * ZCenter / VelScaleSq;
         BulkRichardsonNumber(ICell, KInt) = RiBulk;

         if (KCross < 0 && RiBulk > RiCritical) {
            KCross = K;
            // Levels below the crossing only feed the Ri diagnostics
            if (!FullRiProfile)
               break;
         }
      }

      if (KCross >= KMin) {
         if (KCross > KMin) {
            // Ri values are defined at cell centers, so interpolate on
            // center depths to keep the abscissa consistent.
            // Define quantities to conduct the interpolation of Ri
            // between the two model levels to determine the final
            // boundary layer depth.
            const I4 KAbove    = KCross - 1;
            const I4 KBelow    = KCross;
            const I4 KAboveRi  = KAbove + 1;
            const I4 KBelowRi  = KBelow + 1;
            const Real ZAbove  = Ssh - ZMid(ICell, KAbove);
            const Real ZBelow  = Ssh - ZMid(ICell, KBelow);
            const Real RiAbove = BulkRichardsonNumber(ICell, KAboveRi);
            const Real RiBelow = BulkRichardsonNumber(ICell, KBelowRi);

            const Real H = ZBelow - ZAbove;
            if (H > KPP::NumericalTolerance) {
               // CVMix-style QUAD interpolation for OSBL crossing:
               // - first interior crossing uses zero slope at top point
               // - deeper crossings use upstream slope.
               // Quadratic interpolation is recommended by Danabosglu et al
               // (2006)
               // https://journals.ametsoc.org/view/journals/clim/19/11/jcli3739.1.pdf
               Real SlopeAbove = 0.0_Real;
               if (KCross > KMin + 1) {
                  const I4 KPrev    = KAbove - 1;
                  const I4 KPrevRi  = KPrev + 1;
                  const Real ZPrev  = Ssh - ZMid(ICell, KPrev);
                  const Real RiPrev = BulkRichardsonNumber(ICell, KPrevRi);
                  const Real DZPrev = ZAbove - ZPrev;
                  if (Kokkos::abs(DZPrev) > KPP::NumericalTolerance) {
                     SlopeAbove = (RiAbove - RiPrev) / DZPrev;
                  }
               }

               // In local coordinate T = z - ZAbove:
               // Ri(T) = QuadA T^2 + SlopeAbove T + RiAbove, with QuadA
               // fixed by requiring Ri(H) = RiBelow. The OSBL base is the
               // root of Ri(T) = RiCritical.
               const Real QuadA =
                   (RiBelow - RiAbove - SlopeAbove * H) / (H * H);
               const Real QuadC = RiAbove - RiCritical;

               Real TCross = H;
               if (Kokkos::abs(QuadA) < 1.0e-14_Real) {
                  // Degenerate quadratic -> linear fallback.
                  const Real DRi = RiBelow - RiAbove;
                  if (Kokkos::abs(DRi) > KPP::NumericalTolerance) {
                     const Real Frac = Kokkos::fmax(
                         0.0_Real,
                         Kokkos::fmin(1.0_Real, (RiCritical - RiAbove) / DRi));
                     TCross = Frac * H;
                  }
               } else {
                  const Real Disc =
                      SlopeAbove * SlopeAbove - 4.0_Real * QuadA * QuadC;
                  if (Disc >= 0.0_Real) {
                     const Real SqrtDisc = Kokkos::sqrt(Disc);
                     const Real T1 =
                         (-SlopeAbove + SqrtDisc) / (2.0_Real * QuadA);
                     const Real T2 =
                         (-SlopeAbove - SqrtDisc) / (2.0_Real * QuadA);

                     const bool T1Ok = (T1 >= 0.0_Real && T1 <= H);
                     const bool T2Ok = (T2 >= 0.0_Real && T2 <= H);
                     if (T1Ok && T2Ok) {
                        // Both roots lie in the interval; prefer the one
                        // nearest mid-interval, as CVMix does.
                        const Real Mid = 0.5_Real * H;
                        TCross =
                            (Kokkos::abs(T1 - Mid) <= Kokkos::abs(T2 - Mid))
                                ? T1
                                : T2;
                     } else if (T1Ok) {
                        TCross = T1;
                     } else if (T2Ok) {
                        TCross = T2;
                     } else {
                        TCross = H;
                     }
                  }
               }

               TCross   = Kokkos::fmax(0.0_Real, Kokkos::fmin(H, TCross));
               OBLDepth = ZAbove + TCross;
            } else {
               OBLDepth = ZBelow;
            }
         } else {
            // Match center-based OSBL convention when crossing occurs in
            // the top interval.
            OBLDepth = Ssh - ZMid(ICell, KMin);
         }
      } else {
         OBLDepth = Ssh - ZInterface(ICell, KIntDeep);
      }

      const Real TopLayerThickness =
          Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
      const Real MinOBLDepth = TopLayerThickness;
      const Real MaxOBLDepth = Ssh - ZMid(ICell, KMax);
      // Impose chosen limits on the depth of the OSBL
      OBLDepth = KPP::kppClampOSBLDepth(
          OBLDepth, MinOBLDepth, MaxOBLDepth,
          IceFrac > IceFracThresholdForMinimumOSBL, MinimumOSBLUnderSeaIce);

      const I4 KFinal =
          KPP::kppOSBLIndex(ZInterface, ICell, KMin, KMax, Ssh, OBLDepth);

      OSBLDepth(ICell)      = OBLDepth;
      OSBLDepthIndex(ICell) = KFinal;
   }

 private:
   /// Compile-time bound for the per-thread edge scratch arrays
   static constexpr I4 MaxEdgesBound = HorzMesh::MaxEdgesBound;

   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array1DI4 MinLayerEdgeBot;
   Array1DI4 MaxLayerEdgeTop;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnCell;
   Array1DReal AreaCell;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
};

/// @brief Area-weighted smoothing of the OSBL depth over each cell
/// and its neighbors. This is a laplacian smoothing operation. This suppresses
/// the grid-scale noise that the discrete Ri crossing search and no time
/// history in the OSBL depth calculation can introduce.
class KPPOSBLSmooth {
 public:
   /// Constructor for KPPOSBLSmooth
   KPPOSBLSmooth(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array1DReal &OSBLDepthSmooth, I4 ICell,
                                   const Array1DReal &OSBLDepth) const {
      const I4 NEdges = NEdgesOnCell(ICell);
      Real AreaSum    = 0.0_Real;
      Real BLDSum     = 0.0_Real;
      I4 EdgeCount    = 0;

      for (I4 J = 0; J < NEdges; ++J) {
         const I4 INeighbor = CellsOnCell(ICell, J);
         if (INeighbor == NCellsAll) {
            continue;
         }

         const Real NbrArea = AreaCell(INeighbor);
         BLDSum += 2.0_Real * NbrArea * OSBLDepth(INeighbor);
         AreaSum += 2.0_Real * NbrArea;
         ++EdgeCount;
      }

      if (EdgeCount > 0) {
         const Real SelfArea = AreaCell(ICell);
         BLDSum += OSBLDepth(ICell) * static_cast<Real>(EdgeCount) * SelfArea;
         AreaSum += static_cast<Real>(EdgeCount) * SelfArea;
      }

      if (AreaSum > 0.0_Real) {
         OSBLDepthSmooth(ICell) = BLDSum / AreaSum;
      } else {
         OSBLDepthSmooth(ICell) = OSBLDepth(ICell);
      }
   }

 private:
   I4 NCellsAll;
   Array1DI4 NEdgesOnCell;
   Array2DI4 CellsOnCell;
   Array1DReal AreaCell;
};

/// @brief Re-clamp and re-index the OSBL depth after smoothing
class KPPOSBLCommit {
 public:
   /// Constructor for KPPOSBLCommit
   KPPOSBLCommit(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array1DReal &OSBLDepth,
                                   const Array1DI4 &OSBLDepthIndex, I4 ICell,
                                   const Array1DReal &OSBLDepthSmooth) const {

      const I4 KMin  = MinLayerCell(ICell);
      const I4 KMax  = MaxLayerCell(ICell);
      const Real Ssh = SshCell(ICell);

      const I4 KIntTop = KMin + 1;
      const Real TopLayerThickness =
          Kokkos::abs(ZInterface(ICell, KIntTop) - ZInterface(ICell, KMin));
      const Real MinOBLDepth = TopLayerThickness;
      const Real MaxOBLDepth = Ssh - ZMid(ICell, KMax);

      // The sea-ice minimum is deliberately not reapplied here; it was
      // already enforced before smoothing.
      const Real OBLDepthValue = KPP::kppClampOSBLDepth(
          OSBLDepthSmooth(ICell), MinOBLDepth, MaxOBLDepth, false, 0.0_Real);

      const I4 KFinal =
          KPP::kppOSBLIndex(ZInterface, ICell, KMin, KMax, Ssh, OBLDepthValue);

      OSBLDepth(ICell)      = OBLDepthValue;
      OSBLDepthIndex(ICell) = KFinal;
   }

 private:
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
};

/// @brief Initialize the coefficient arrays with zero KPP contribution, or
/// with precomputed interior mixing for matched-coefficient construction.
class KPPCoeffsInit {
 public:
   bool UseInteriorMix = false; ///< Interior coefficients were supplied

   KOKKOS_FUNCTION void operator()(const Array2DReal &VertDiff,
                                   const Array2DReal &VertVisc,
                                   const Array2DReal &VertNonLocalFlux,
                                   const Array2DReal &TurbulentVelocityScale,
                                   I4 ICell, I4 K,
                                   const Array2DReal &InteriorVertDiff,
                                   const Array2DReal &InteriorVertVisc) const {

      VertDiff(ICell, K) =
          UseInteriorMix ? InteriorVertDiff(ICell, K) : 0.0_Real;
      VertVisc(ICell, K) =
          UseInteriorMix ? InteriorVertVisc(ICell, K) : 0.0_Real;
      VertNonLocalFlux(ICell, K)       = 0.0_Real;
      TurbulentVelocityScale(ICell, K) = 0.0_Real;
   }
};

/// @brief Stage 2 kernel: KPP profile-based diffusivity, viscosity and
/// non-local flux, with optional enhanced mixing at the OSBL base.
class KPPMixingCoeffs {
 public:
   Real SurfaceLayerExtent   = KPP::SurfaceLayerExtent; ///< Frac of OSBL depth
   bool UseEnhancedDiffusion = true;  ///< Apply enhanced mixing at OSBL base
   bool UseInteriorMix       = false; ///< Interior coefficients were supplied
   bool UseMatchedShapes     = false; ///< Match interior value at the OSBL base
   /// Non-local flux normalization; depends only on SurfaceLayerExtent, so the
   /// caller sets it once per call rather than every cell recomputing it.
   Real NonLocalCs = 0.0_Real;

   /// Constructor for KPPMixingCoeffs
   KPPMixingCoeffs(const HorzMesh *Mesh, const VertCoord *VCoord);

   KOKKOS_FUNCTION void operator()(const Array2DReal &VertDiff,
                                   const Array2DReal &VertVisc,
                                   const Array2DReal &VertNonLocalFlux,
                                   const Array2DReal &TurbulentVelocityScale,
                                   I4 ICell, const Array1DReal &OSBLDepth,
                                   const Array1DI4 &OSBLDepthIndex,
                                   const Array1DReal &SurfaceFrictionVelocity,
                                   const Array1DReal &SurfaceBuoyancyFlux,
                                   const Array1DReal &LangmuirFactor,
                                   const Array2DReal &InteriorVertDiff,
                                   const Array2DReal &InteriorVertVisc) const {

      const Real H = OSBLDepth(ICell);

      const I4 KMin = MinLayerCell(ICell);
      const I4 KMax = MaxLayerCell(ICell);
      if (KMin < 0 || KMin >= NVertLayers || KMax < KMin) {
         return;
      }
      const I4 KMatch = Kokkos::min(KMax + 1, OSBLDepthIndex(ICell) + 1);

      // KPP depths are measured below the free surface, so geometric
      // heights must be offset by the sea surface height.
      const Real Ssh = SshCell(ICell);

      const Real UStar    = SurfaceFrictionVelocity(ICell);
      const Real BuoyFlux = SurfaceBuoyancyFlux(ICell);
      // Same Langmuir enhancement applied to the Stage 1 OSBL search, so the
      // mixing magnitude is consistent with the diagnosed boundary layer.
      const Real BuoyFluxEff = BuoyFlux * LangmuirFactor(ICell);

      for (I4 K = KMin; K <= KMax + 1; ++K) {
         const Real ZDepth = Ssh - ZInterface(ICell, K);

         // Check if within OSBL using depth below the free surface.
         if (ZDepth <= H && H > 0.0_Real) {
            // Normalized depth in Omega sign convention: sigma in [-1,0].
            Real Sigma = -ZDepth / H;
            Sigma      = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));

            // CVMix-style turbulent scales: w = kappa*u*/phi in general,
            // with explicit free-convection limits when u*=0. The scales
            // are frozen at the surface-layer depth below the surface
            // layer, so SigmaLoc is capped at SurfaceLayerExtent.
            const Real SigmaCoord = -Sigma; // [0,1]
            const Real SigmaLoc   = Kokkos::fmin(
                SurfaceLayerExtent, Kokkos::fmax(0.0_Real, SigmaCoord));

            Real WMTurb = 0.0_Real;
            Real WSTurb = 0.0_Real;
            KPP::kppTurbScales(UStar, BuoyFluxEff, H, SigmaLoc, VonKar, WMTurb,
                               WSTurb);

            // For MatchBoth, the shape value the KPP profile must reach at
            // the OSBL base so that it joins the interior coefficient there.
            const Real MatchViscShape =
                UseMatchedShapes
                    ? KPP::kppMatchShape(InteriorVertVisc(ICell, KMatch), H,
                                         WMTurb)
                    : 0.0_Real;
            const Real MatchDiffShape =
                UseMatchedShapes
                    ? KPP::kppMatchShape(InteriorVertDiff(ICell, KMatch), H,
                                         WSTurb)
                    : 0.0_Real;

            // Momentum mixing contribution.
            const Real ShapeM =
                UseMatchedShapes ? KPP::kppShapeMatched(Sigma, MatchViscShape)
                                 : KPP::kppShapeMomentum(Sigma);
            VertVisc(ICell, K) = H * WMTurb * ShapeM;

            // Tracer mixing contribution.
            const Real ShapeS =
                UseMatchedShapes ? KPP::kppShapeMatched(Sigma, MatchDiffShape)
                                 : KPP::kppShapeScalar(Sigma);
            VertDiff(ICell, K)               = H * WSTurb * ShapeS;
            TurbulentVelocityScale(ICell, K) = WSTurb;

            // Non-local flux: C_s * G(sigma).
            // C_s = C* * kappa * (c_s * kappa * epsilon)^(1/3)
            // per Large et al. (1994) eq. 20 (~6.33 with default constants).
            // The non-local shape is always the unmatched scalar shape,
            // independent of MatchTechnique, so that gamma vanishes at the
            // OSBL base. The matched shape is non-zero there by construction,
            // which would leave a non-local flux at the base that jumps to
            // zero just below it. CVMix likewise keeps matching (a
            // diffusivity choice) separate from the non-local shape.
            const Real NonLocalShape = KPP::kppShapeScalar(Sigma);

            // Match CVMix behavior: apply non-local term only when
            // surface buoyancy forcing is unstable/neutral.
            if (BuoyFlux <= 0.0_Real) {
               VertNonLocalFlux(ICell, K) = NonLocalCs * NonLocalShape;
            } else {
               VertNonLocalFlux(ICell, K) = 0.0_Real;
            }

         } else {
            // Below OSBL: preserve interior values for MatchBoth, otherwise
            // no KPP contribution.
            VertDiff(ICell, K) =
                UseInteriorMix ? InteriorVertDiff(ICell, K) : 0.0_Real;
            VertVisc(ICell, K) =
                UseInteriorMix ? InteriorVertVisc(ICell, K) : 0.0_Real;
            VertNonLocalFlux(ICell, K)       = 0.0_Real;
            TurbulentVelocityScale(ICell, K) = 0.0_Real;
         }
      }

      // Optional enhanced diffusion/viscosity treatment at OSBL base.
      // Match Large et al (1994) Appendix D weighting at the interface nearest
      // H: the OSBL base rarely lands on an interface, so the coefficient at
      // the neighboring interface KTarget is replaced by a quadratic blend
      // of the KPP value extrapolated from KKtup and the value already
      // there, weighted by where H falls between the two cell centers.
      if (UseEnhancedDiffusion && H > 0.0_Real) {
         const I4 KOBL =
             Kokkos::max(KMin, Kokkos::min(OSBLDepthIndex(ICell), KMax));
         const Real ZMidOBL = Ssh - ZMid(ICell, KOBL);

         const bool TargetOutsideOSBL = H >= ZMidOBL;
         const I4 KKtup =
             TargetOutsideOSBL ? KOBL : Kokkos::max(KMin, KOBL - 1);
         const I4 KTarget = TargetOutsideOSBL ? Kokkos::min(KOBL + 1, KMax + 1)
                                              : Kokkos::max(KMin + 1, KOBL);

         const Real ZKtup = Ssh - ZMid(ICell, KKtup);
         const Real ZNext = (KKtup < KMax)
                                ? (Ssh - ZMid(ICell, KKtup + 1))
                                : (Ssh - ZInterface(ICell, KKtup + 1));
         const Real Delta = Kokkos::fmax(
             0.0_Real, Kokkos::fmin(1.0_Real, (H - ZKtup) / (ZNext - ZKtup)));
         const Real OneMinusDelta = 1.0_Real - Delta;

         Real SigmaKtup = -ZKtup / H;
         SigmaKtup = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, SigmaKtup));
         const Real SigmaCoord = -SigmaKtup;
         const Real SigmaLoc   = Kokkos::fmin(SurfaceLayerExtent,
                                              Kokkos::fmax(0.0_Real, SigmaCoord));

         Real WMKtup = 0.0_Real;
         Real WSKtup = 0.0_Real;
         KPP::kppTurbScales(UStar, BuoyFluxEff, H, SigmaLoc, VonKar, WMKtup,
                            WSKtup);

         const Real MatchViscShape =
             UseMatchedShapes ? KPP::kppMatchShape(
                                    InteriorVertVisc(ICell, KMatch), H, WMKtup)
                              : 0.0_Real;
         const Real MatchDiffShape =
             UseMatchedShapes ? KPP::kppMatchShape(
                                    InteriorVertDiff(ICell, KMatch), H, WSKtup)
                              : 0.0_Real;

         const Real ViscKtup =
             H * WMKtup *
             (UseMatchedShapes ? KPP::kppShapeMatched(SigmaKtup, MatchViscShape)
                               : KPP::kppShapeMomentum(SigmaKtup));
         const Real DiffKtup =
             H * WSKtup *
             (UseMatchedShapes ? KPP::kppShapeMatched(SigmaKtup, MatchDiffShape)
                               : KPP::kppShapeScalar(SigmaKtup));

         const Real ViscProfile = VertVisc(ICell, KTarget);
         const Real DiffProfile = VertDiff(ICell, KTarget);

         const Real EnhVisc = OneMinusDelta * OneMinusDelta * ViscKtup +
                              Delta * Delta * ViscProfile;
         const Real EnhDiff = OneMinusDelta * OneMinusDelta * DiffKtup +
                              Delta * Delta * DiffProfile;

         const Real OldVisc =
             UseInteriorMix ? InteriorVertVisc(ICell, KTarget) : 0.0_Real;
         const Real OldDiff =
             UseInteriorMix ? InteriorVertDiff(ICell, KTarget) : 0.0_Real;
         const Real NewVisc = OneMinusDelta * OldVisc + Delta * EnhVisc;
         const Real NewDiff = OneMinusDelta * OldDiff + Delta * EnhDiff;

         VertVisc(ICell, KTarget) = NewVisc;
         VertDiff(ICell, KTarget) = NewDiff;

         // Keep the non-local term consistent with the rescaled diffusivity
         if (!TargetOutsideOSBL && DiffProfile != 0.0_Real) {
            VertNonLocalFlux(ICell, KTarget) =
                VertNonLocalFlux(ICell, KTarget) * NewDiff / DiffProfile;
         } else if (!TargetOutsideOSBL) {
            VertNonLocalFlux(ICell, KTarget) = 0.0_Real;
         }
      }
   }

 private:
   I4 NVertLayers;
   Array1DI4 MinLayerCell;
   Array1DI4 MaxLayerCell;
   Array2DReal ZInterface;
   Array2DReal ZMid;
   Array1DReal SshCell;
};

/// @brief KPP Boundary Layer Mixing Scheme
///
/// Implements the K-Profile Parameterization following Large et al. (1994)
/// with optional Langmuir turbulence enhancement. Computes vertical
/// diffusivity, viscosity, and non-local flux coefficients for the OSBL.
///
/// Two-stage computation:
/// 1. Stage 1: Compute OSBL depth from bulk Richardson criterion and smooth
/// 2. Stage 2: Compute viscosity and diffusivity coefficients
///
class KPPMix {

 public:
   /// @brief Singleton instance management
   static KPPMix *getInstance();
   static void init();
   static void destroyInstance();

   /// @brief Prepare the current KPP state, surface forcing, and coefficients
   /// from canonical ocean state, EOS, and forcing fields.
   /// InteriorVertDiff/InteriorVertVisc supply a non-KPP interior mixing
   /// estimate for MatchBoth; ignored otherwise.
   void update(const Array3DReal &TracerArray, I4 TempTracerIndex,
               I4 SaltTracerIndex, const Array2DReal &NormalVelocity,
               Eos *EqState, const Forcing *ForcingState, bool UseTracerForcing,
               const Array2DReal &InteriorVertDiff = Array2DReal(),
               const Array2DReal &InteriorVertVisc = Array2DReal());

   /// @brief Prepare KPP fields for the current time step from canonical
   /// ocean state. Looks up Temperature/Salinity tracer indices, gathers the
   /// MatchBoth interior-mixing snapshot from VertMix, and calls update().
   /// Must be called once per step before the first tendency evaluation.
   void computeKPPFields(const OceanState *State,
                         const Array3DReal &TracerArray, int ThickTimeLevel,
                         int VelTimeLevel, bool UseTracerForcing);

   /// @brief Main computation routine
   /// Calls Stage 1 and Stage 2 computation in sequence
   ///
   /// Input arrays should be pre-populated with current state.
   /// NormalVelocity and TangentialVelocity are edge-based quantities
   /// (C-grid convention): dimensions [NEdges × NVertLayers].
   /// Output arrays are computed in-place.
   void computeKPPMix(
       const Array2DReal
           &PotentialDensity, ///< Density (kg/m³) [NCells×NLevels]
       const Array2DReal &NormalVelocity,     ///< Normal vel on edges (m/s)
       const Array2DReal &TangentialVelocity, ///< Tangential vel on edges (m/s)
       const Array1DReal &SurfaceFrictionVelocity, ///< u* (m/s)
       const Array1DReal &SurfaceBuoyancyFlux,     ///< B_0 (m²/s³)
       const Array2DReal &BruntVaisalaFreqSq,      ///< N² (s⁻²)
       const Array1DReal &IceFraction,             ///< Sea ice cover (0-1)
       const Array1DReal &WindSpeed10m =
           Array1DReal(), ///< Wind for Langmuir (m/s)
       const Array2DReal &InteriorVertDiff =
           Array2DReal(), ///< Interior diffusivity for MatchBoth
       const Array2DReal &InteriorVertVisc =
           Array2DReal() ///< Interior viscosity for MatchBoth
   );

   // =======================================================================
   // Output Fields
   // =======================================================================

   /// @brief Vertical diffusivity at layer interfaces (m²/s)
   /// Size: [nCells][nLevels+1]
   Array2DReal VertDiff;

   /// @brief Vertical viscosity at layer interfaces (m²/s)
   /// Size: [nCells][nLevels+1]
   Array2DReal VertVisc;

   /// @brief Boundary layer depth (m)
   /// Size: [nCells]
   Array1DReal OSBLDepth;

   /// @brief OSBL depth as layer index
   /// Size: [nCells]
   Array1DI4 OSBLDepthIndex;

   /// @brief Non-local flux coefficient profile G(σ) (dimensionless)
   /// Size: [nCells][nLevels+1]
   /// Applied to surface tracer fluxes to compute vertical transport
   Array2DReal VertNonLocalFlux;

   /// @brief Bulk Richardson number profile used in OSBL search (dimensionless)
   /// Size: [nCells][nLevels+1]
   Array2DReal BulkRichardsonNumber;

   /// @brief Shear contribution to bulk Richardson denominator (m^2/s^2)
   /// Size: [nCells][nLevels+1]
   Array2DReal BulkRichardsonShear;

   /// @brief Unresolved shear contribution Vt^2 (m^2/s^2)
   /// Size: [nCells][nLevels+1]
   Array2DReal UnresolvedShear;

   /// @brief Buoyancy jump (density anomaly converted to buoyancy) (m/s²)
   /// Size: [nCells][nLevels+1]
   /// Captures delta_b = g * delta_rho / rho_sw at each layer during OSBL
   /// search
   Array2DReal BuoyancyJump;

   /// @brief Turbulent velocity scale profile (m/s), tracer branch
   /// Size: [nCells][nLevels+1]
   Array2DReal TurbulentVelocityScale;

   /// @brief Potential density used by KPP OSBL search (kg/m^3)
   /// Size: [nCells][nLevels]
   Array2DReal PotentialDensity;

   /// @brief Surface friction velocity u* (m/s)
   /// Size: [nCells]
   Array1DReal SurfaceFrictionVelocity;

   /// @brief Surface buoyancy flux B_0 (m²/s³)
   /// Size: [nCells]
   Array1DReal SurfaceBuoyancyFlux;

   // =======================================================================
   // Configuration Parameters
   // =======================================================================

   bool Enabled = true; ///< Enable/disable KPP mixing

   // Defaults below may be overridden from the Config file; where a value also
   // appears in KPPConstants.h, that is the authoritative default.
   Real CriticalRichardson = KPP::CriticalRi;         ///< Ri_crit for OSBL base
   Real SurfaceLayerExtent = KPP::SurfaceLayerExtent; ///< Frac of OSBL depth

   bool UseLangmuirTurbulence = true;  ///< Apply wave enhancement
   bool DebugDiagnostics      = false; ///< Print per-step KPP diagnostics

   // Ice/Langmuir controls (kept configurable to match reference semantics)
   /// Disable Langmuir above this ice fraction
   Real IceFractionThresholdForLangmuir = KPP::IceFracThresh;
   /// Apply minimum OSBL depth above this ice fraction
   Real IceFractionThresholdForMinimumOSBL = KPP::IceSuppressThresh;
   /// Min OSBL depth under sea ice (m)
   Real MinimumOSBLUnderSeaIce = KPP::MinOSBLUnderIce;

   // KPP matching/profile controls (CVMix-style semantics)
   KPPMatchType MatchTechnique = KPPMatchType::SimpleShapes;
   std::string InterpType2Str  = "LMD94"; ///< Linear, Quadratic, Cubic, LMD94
   bool UseEnhancedDiffusion   = true; ///< Apply enhanced mixing at OSBL base
   bool UseOSBLSmoothing = true; ///< Apply MPAS-style OSBL horizontal smoothing

   // Field names for I/O
   std::string BuoyancyJumpFldName;
   std::string VertDiffFldName;
   std::string VertViscFldName;
   std::string OSBLDepthFldName;
   std::string OSBLDepthIndexFldName;
   std::string NonLocalFluxFldName;
   std::string BulkRichardsonFldName;
   std::string BulkRichardsonShearFldName;
   std::string UnresolvedShearFldName;
   std::string TurbulentVelScaleFldName;
   std::string PotentialDensityFldName;
   std::string SurfFricVelFldName;
   std::string SurfBuoyFluxFldName;
   std::string LangmuirFactorFldName;
   std::string Name;

 private:
   /// @brief Private constructor for singleton pattern
   KPPMix(const std::string &InName, const HorzMesh *InMesh,
          const VertCoord *InVCoord);

   /// @brief Private destructor
   ~KPPMix();

   /// @brief Static singleton instance
   static KPPMix *Instance;

   /// @brief Mesh and coordinate references
   const HorzMesh *Mesh;
   const VertCoord *VCoord;

   /// @brief KPP-only workspaces reused across time steps
   Array2DReal RefPressure;
   Array2DReal TangentialVelocity;
   Array1DReal LangmuirFactor;
   Array1DReal OSBLDepthSmooth;

 public:
   /// @brief Stage 1: Compute OSBL depth using edge-based velocity shear
   void computeOSBLDepth(const Array2DReal &PotentialDensity,
                         const Array2DReal &NormalVelocity,
                         const Array2DReal &TangentialVelocity,
                         const Array1DReal &SurfaceFrictionVelocity,
                         const Array1DReal &SurfaceBuoyancyFlux,
                         const Array2DReal &BruntVaisalaFreqSq,
                         const Array1DReal &IceFraction,
                         const Array1DReal &WindSpeed10m);

   /// @brief Stage 2: Compute KPP mixing contribution or matched coefficients
   void computeMixingCoefficients(
       const Array1DReal &SurfaceFrictionVelocity,
       const Array1DReal &SurfaceBuoyancyFlux,
       const Array2DReal &InteriorVertDiff = Array2DReal(),
       const Array2DReal &InteriorVertVisc = Array2DReal());

 private:
   /// @brief Print targeted diagnostics for KPP troubleshooting
   void logDiagnostics(const Array2DReal &PotentialDensity,
                       const Array1DReal &SurfaceFrictionVelocity,
                       const Array1DReal &SurfaceBuoyancyFlux,
                       const Array1DReal &WindSpeed10m);

   /// @brief Register fields with I/O system
   void defineFields();

   // Delete copy and move constructors/assignment
   KPPMix(const KPPMix &)            = delete;
   KPPMix &operator=(const KPPMix &) = delete;
   KPPMix(KPPMix &&)                 = delete;
   KPPMix &operator=(KPPMix &&)      = delete;

}; // class KPPMix

} // namespace OMEGA

#endif // OMEGA_KPP_MIX_H
