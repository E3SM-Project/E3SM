//===-- ocn/KPPMix.cpp - KPP Boundary Layer Mixing Implementation --*- C++
//-*-===//
//
/// \file
/// \brief Implementation of KPP boundary layer mixing scheme
///
/// This file implements the KPPMix class for computing ocean boundary layer
/// mixing following Large et al. (1994) with optional Langmuir enhancement.
//
//===----------------------------------------------------------------------===//

#include "KPPMix.h"
#include "DataTypes.h"
#include "Eos.h"
#include "Error.h"
#include "FillValues.h"
#include "Forcing.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "KPPConstants.h"
#include "Logging.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "VertMix.h"
#include <limits>

namespace OMEGA {

// Singleton instance
KPPMix *KPPMix::Instance = nullptr;

namespace {

bool matchTypeFromString(const std::string &Name, KPPMatchType &Type) {
   if (Name == "SimpleShapes") {
      Type = KPPMatchType::SimpleShapes;
      return true;
   }
   if (Name == "MatchBoth") {
      Type = KPPMatchType::MatchBoth;
      return true;
   }
   return false;
}

const char *matchTypeName(KPPMatchType Type) {
   switch (Type) {
   case KPPMatchType::MatchBoth:
      return "MatchBoth";
   default:
      return "SimpleShapes";
   }
}

} // namespace

KPPSurfaceForcingOnCell::KPPSurfaceForcingOnCell(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord,
                                                 const Eos *EosInst)
    : LinearDRhodT(EosInst->getLinearDRhodT()),
      LinearDRhodS(EosInst->getLinearDRhodS()),
      NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      EosChoice(EosInst->EosChoice) {}

/// Constructor for KPPMix
KPPMix::KPPMix(const std::string &InName, const HorzMesh *InMesh,
               const VertCoord *InVCoord)
    : Name(InName), Mesh(InMesh), VCoord(InVCoord) {

   // Allocate output arrays
   VertDiff  = Array2DReal("VertDiff", Mesh->NCellsSize, VCoord->NVertLayersP1);
   VertVisc  = Array2DReal("VertVisc", Mesh->NCellsSize, VCoord->NVertLayersP1);
   OSBLDepth = Array1DReal("OSBLDepth", Mesh->NCellsSize);
   OSBLDepthIndex = Array1DI4("OSBLDepthIndex", Mesh->NCellsSize);
   VertNonLocalFlux =
       Array2DReal("VertNonLocalFlux", Mesh->NCellsSize, VCoord->NVertLayersP1);
   BulkRichardsonNumber = Array2DReal("BulkRichardsonNumber", Mesh->NCellsSize,
                                      VCoord->NVertLayersP1);
   BulkRichardsonShear  = Array2DReal("BulkRichardsonShear", Mesh->NCellsSize,
                                      VCoord->NVertLayersP1);
   UnresolvedShear =
       Array2DReal("UnresolvedShear", Mesh->NCellsSize, VCoord->NVertLayersP1);
   BuoyancyJump =
       Array2DReal("BuoyancyJump", Mesh->NCellsSize, VCoord->NVertLayersP1);
   TurbulentVelocityScale = Array2DReal(
       "TurbulentVelocityScale", Mesh->NCellsSize, VCoord->NVertLayersP1);
   PotentialDensity =
       Array2DReal("PotentialDensity", Mesh->NCellsSize, VCoord->NVertLayers);
   SurfaceFrictionVelocity =
       Array1DReal("SurfaceFrictionVelocity", Mesh->NCellsSize);
   SurfaceBuoyancyFlux = Array1DReal("SurfaceBuoyancyFlux", Mesh->NCellsSize);
   RefPressure =
       Array2DReal("KPP-RefPressure", Mesh->NCellsSize, VCoord->NVertLayers);
   TangentialVelocity = Array2DReal("KPP-TangentialVelocity", Mesh->NEdgesSize,
                                    VCoord->NVertLayers);
   LangmuirFactor     = Array1DReal("KPP-LangmuirFactor", Mesh->NCellsSize);
   OSBLDepthSmooth    = Array1DReal("KPP-OSBLDepthSmooth", Mesh->NCellsSize);

   // Default to a no-op enhancement so Stage 2 is well-defined even if
   // called before Stage 1 (e.g. standalone unit tests) has computed it.
   deepCopy(LangmuirFactor, 1.0_Real);

   // Set field names
   VertDiffFldName            = "VertDiff";
   VertViscFldName            = "VertVisc";
   OSBLDepthFldName           = "OSBLDepth";
   OSBLDepthIndexFldName      = "OSBLDepthIndex";
   NonLocalFluxFldName        = "VertNonLocalFlux";
   BulkRichardsonFldName      = "BulkRichardsonNumber";
   BulkRichardsonShearFldName = "BulkRichardsonShear";
   UnresolvedShearFldName     = "UnresolvedShear";
   BuoyancyJumpFldName        = "BuoyancyJump";
   TurbulentVelScaleFldName   = "TurbulentVelocityScale";
   PotentialDensityFldName    = "PotentialDensity";
   SurfFricVelFldName         = "SurfaceFrictionVelocity";
   SurfBuoyFluxFldName        = "SurfaceBuoyancyFlux";
   LangmuirFactorFldName      = "KPPLangmuirFactor";

   if (Name != "Default") {
      VertDiffFldName.append(Name);
      VertViscFldName.append(Name);
      OSBLDepthFldName.append(Name);
      OSBLDepthIndexFldName.append(Name);
      NonLocalFluxFldName.append(Name);
      BulkRichardsonFldName.append(Name);
      BulkRichardsonShearFldName.append(Name);
      UnresolvedShearFldName.append(Name);
      BuoyancyJumpFldName.append(Name);
      TurbulentVelScaleFldName.append(Name);
      PotentialDensityFldName.append(Name);
      SurfFricVelFldName.append(Name);
      SurfBuoyFluxFldName.append(Name);
      LangmuirFactorFldName.append(Name);
   }

   defineFields();
}

/// Destructor for KPPMix
KPPMix::~KPPMix() {}

/// Get singleton instance
KPPMix *KPPMix::getInstance() { return Instance; }

/// Destroy singleton instance
void KPPMix::destroyInstance() {
   delete Instance;
   Instance = nullptr;
}

/// Initialize KPPMix from configuration
void KPPMix::init() {
   if (!Instance) {
      Instance = new KPPMix("Default", HorzMesh::getDefault(),
                            VertCoord::getDefault());
   }

   Error Err;
   KPPMix *DefKPPMix   = KPPMix::getInstance();
   Config *OmegaConfig = Config::getOmegaConfig();

   // Get VertMix config group
   Config VertMixConfig("VertMix");
   Err += OmegaConfig->get(VertMixConfig);
   CHECK_ERROR_ABORT(Err, "KPPMix::init: VertMix group not found in Config");

   // Get KPP config subgroup
   Config KPPConfig("KPP");
   Err += VertMixConfig.get(KPPConfig);
   if (Err.isFail()) {
      LOG_WARN("KPPMix::init: KPP subgroup not found, using defaults");
      return; // Continue with defaults
   }

   // Read KPP parameters
   bool Enable = true;
   Err += KPPConfig.get("Enable", Enable);
   DefKPPMix->Enabled = Enable;

   Err += KPPConfig.get("CriticalBulkRichardsonNumber",
                        DefKPPMix->CriticalRichardson);
   // Clamp once here so kernels can trust CriticalRichardson is positive
   DefKPPMix->CriticalRichardson =
       Kokkos::fmax(KPP::NumericalTolerance, DefKPPMix->CriticalRichardson);
   Err += KPPConfig.get("SurfaceLayerExtent", DefKPPMix->SurfaceLayerExtent);
   // Clamp once here so kernels can trust SurfaceLayerExtent is in [0,1]
   DefKPPMix->SurfaceLayerExtent = Kokkos::fmax(
       0.0_Real, Kokkos::fmin(1.0_Real, DefKPPMix->SurfaceLayerExtent));

   // KPP matching/profile semantics.
   std::string MatchStr = "SimpleShapes";
   Error MatchErr       = KPPConfig.get("MatchTechnique", MatchStr);
   if (!MatchErr.isSuccess()) {
      MatchErr.reset();
   }
   Error InterpErr = KPPConfig.get("InterpType2", DefKPPMix->InterpType2Str);
   if (!InterpErr.isSuccess()) {
      InterpErr.reset();
   }
   Error EnhancedErr =
       KPPConfig.get("UseEnhancedDiffusion", DefKPPMix->UseEnhancedDiffusion);
   if (!EnhancedErr.isSuccess()) {
      EnhancedErr.reset();
   }
   Error BLDSmoothErr =
       KPPConfig.get("UseOSBLSmoothing", DefKPPMix->UseOSBLSmoothing);
   if (!BLDSmoothErr.isSuccess()) {
      BLDSmoothErr.reset();
   }

   if (!matchTypeFromString(MatchStr, DefKPPMix->MatchTechnique)) {
      ABORT_ERROR("KPPMix::init: Invalid MatchTechnique='{}', must be "
                  "SimpleShapes or MatchBoth",
                  MatchStr);
   }

   // Wave and flux options
   Err +=
       KPPConfig.get("UseLangmuirTurbulence", DefKPPMix->UseLangmuirTurbulence);
   Err += KPPConfig.get("IceFractionThresholdForLangmuir",
                        DefKPPMix->IceFractionThresholdForLangmuir);
   Err += KPPConfig.get("IceFractionThresholdForMinimumOSBL",
                        DefKPPMix->IceFractionThresholdForMinimumOSBL);
   Err += KPPConfig.get("MinimumOSBLUnderSeaIce",
                        DefKPPMix->MinimumOSBLUnderSeaIce);
   Error DebugErr =
       KPPConfig.get("DebugDiagnostics", DefKPPMix->DebugDiagnostics);
   if (!DebugErr.isSuccess()) {
      DebugErr.reset();
      DefKPPMix->DebugDiagnostics = false;
   }
   if (DefKPPMix->DebugDiagnostics) {
      LOG_WARN("KPP debug diagnostics enabled");
   }

   LOG_WARN("KPPMix::init: KPP initialized enabled={} debugDiagnostics={} "
            "match={}",
            DefKPPMix->Enabled, DefKPPMix->DebugDiagnostics,
            matchTypeName(DefKPPMix->MatchTechnique));
}

void KPPMix::update(const Array3DReal &TracerArray, I4 TempTracerIndex,
                    I4 SaltTracerIndex, const Array2DReal &NormalVelocity,
                    Eos *EqState, const Forcing *ForcingState,
                    bool UseTracerForcing, const Array2DReal &InteriorVertDiff,
                    const Array2DReal &InteriorVertVisc) {

   OMEGA_REQUIRE(EqState, "KPPMix::update: null Eos pointer");
   OMEGA_REQUIRE(ForcingState, "KPPMix::update: null Forcing pointer");

   const I4 NCellsAll   = Mesh->NCellsAll;
   const I4 NVertLayers = VCoord->NVertLayers;

   Array2DReal ConservTemp =
       Kokkos::subview(TracerArray, TempTracerIndex, Kokkos::ALL, Kokkos::ALL);
   Array2DReal AbsSalinity =
       Kokkos::subview(TracerArray, SaltTracerIndex, Kokkos::ALL, Kokkos::ALL);

   OMEGA_SCOPE(LocPressureMid, VCoord->PressureMid);

   PotentialDensityOnCell PotentialDensityCalc;
   OMEGA_SCOPE(LocPotentialDensityCalc, PotentialDensityCalc);
   OMEGA_SCOPE(LocRefPressure, RefPressure);
   OMEGA_SCOPE(LocSurfacePressure, VCoord->SurfacePressure);
   parallelFor(
       "KPP-PotentialDensityPressure", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          LocPotentialDensityCalc.computeRefPressure(LocRefPressure, ICell, K,
                                                     LocSurfacePressure);
       });
   EqState->computeSpecVolDisp(ConservTemp, AbsSalinity, LocRefPressure, 0);

   OMEGA_SCOPE(LocSpecVolPotential, EqState->SpecVolDisplaced);
   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   parallelFor(
       "KPP-PotentialDensity", {NCellsAll, NVertLayers},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          LocPotentialDensityCalc(LocPotentialDensity, ICell, K,
                                  LocSpecVolPotential);
       });

   {
      TangentialReconOnEdge TanReconEdge(Mesh);
      OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
      OMEGA_SCOPE(LocMinLayerEdgeBot, VCoord->MinLayerEdgeBot);
      OMEGA_SCOPE(LocMaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    TanReconEdge(LocTangentialVelocity, IEdge, KChunk,
                                 NormalVelocity);
                 });
          });
   }

   const auto &SfcStress     = ForcingState->SfcStressForcing;
   const auto &TracerForcing = ForcingState->TracerForcing;
   KPPSurfaceForcingOnCell SurfaceForcing(Mesh, VCoord, EqState);
   SurfaceForcing.UseTracerForcing = UseTracerForcing;

   OMEGA_SCOPE(LocSurfaceForcing, SurfaceForcing);
   OMEGA_SCOPE(LocFrictionVelocity, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocBuoyancyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(LocSpecVol, EqState->SpecVol);
   OMEGA_SCOPE(ZonalStress, SfcStress.ZonalStressCell);
   OMEGA_SCOPE(MeridStress, SfcStress.MeridStressCell);
   OMEGA_SCOPE(LatentHeatFluxEvap, TracerForcing.LatentHeatFluxEvapCell);
   OMEGA_SCOPE(SensibleHeatFlux, TracerForcing.SensibleHeatFluxCell);
   OMEGA_SCOPE(LongWaveHeatFluxUp, TracerForcing.LongWaveHeatFluxUpCell);
   OMEGA_SCOPE(LongWaveHeatFluxDown, TracerForcing.LongWaveHeatFluxDownCell);
   OMEGA_SCOPE(SeaIceHeatFlux, TracerForcing.SeaIceHeatFluxCell);
   OMEGA_SCOPE(ShortWaveHeatFlux, TracerForcing.ShortWaveHeatFluxCell);
   OMEGA_SCOPE(SnowFlux, TracerForcing.SnowFluxCell);
   OMEGA_SCOPE(RainFlux, TracerForcing.RainFluxCell);
   OMEGA_SCOPE(EvaporationFlux, TracerForcing.EvaporationFluxCell);
   OMEGA_SCOPE(SeaIceFreshWaterFlux, TracerForcing.SeaIceFreshWaterFluxCell);
   OMEGA_SCOPE(IceRunoffFlux, TracerForcing.IceRunoffFluxCell);
   OMEGA_SCOPE(RiverRunoffFlux, TracerForcing.RiverRunoffFluxCell);
   OMEGA_SCOPE(SeaIceSaltFlux, TracerForcing.SeaIceSaltFluxCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);

   if (UseTracerForcing) {
      I4 FreshPositiveSaltNegativeCount = 0;
      I4 FreshNegativeSaltPositiveCount = 0;
      parallelReduce(
          "KPP-MixedSeaIceFlux-FreshPositiveSaltNegative", {Mesh->NCellsOwned},
          KOKKOS_LAMBDA(I4 ICell, I4 & Count) {
             const I4 KTop = MinLayerCell(ICell);
             if (KTop >= 0 && KTop < NVertLayers &&
                 SeaIceFreshWaterFlux(ICell) > 0.0_Real &&
                 SeaIceSaltFlux(ICell) < 0.0_Real) {
                ++Count;
             }
          },
          Kokkos::Sum<I4>(FreshPositiveSaltNegativeCount));
      parallelReduce(
          "KPP-MixedSeaIceFlux-FreshNegativeSaltPositive", {Mesh->NCellsOwned},
          KOKKOS_LAMBDA(I4 ICell, I4 & Count) {
             const I4 KTop = MinLayerCell(ICell);
             if (KTop >= 0 && KTop < NVertLayers &&
                 SeaIceFreshWaterFlux(ICell) < 0.0_Real &&
                 SeaIceSaltFlux(ICell) > 0.0_Real) {
                ++Count;
             }
          },
          Kokkos::Sum<I4>(FreshNegativeSaltPositiveCount));

      const I4 MixedSignCount =
          FreshPositiveSaltNegativeCount + FreshNegativeSaltPositiveCount;
      if (MixedSignCount > 0) {
         const Real OwnedCellCount = static_cast<Real>(Mesh->NCellsOwned);
         LOG_INFO("KPPMix: mixed-sign sea-ice fluxes in {} owned cells "
                  "({:.3f}%): freshwater>0/salt<0={} ({:.3f}%), "
                  "freshwater<0/salt>0={} ({:.3f}%)",
                  MixedSignCount, 100.0_Real * MixedSignCount / OwnedCellCount,
                  FreshPositiveSaltNegativeCount,
                  100.0_Real * FreshPositiveSaltNegativeCount / OwnedCellCount,
                  FreshNegativeSaltPositiveCount,
                  100.0_Real * FreshNegativeSaltPositiveCount / OwnedCellCount);
      }
   }

   // The individual surface flux components are only valid on owned cells
   // (they are not halo exchanged), so compute u* and B0 on owned cells and
   // exchange them for the halo cells used by the OSBL smoothing and by the
   // edge-based implicit momentum mixing.
   parallelFor(
       "KPP-SurfaceForcing", {Mesh->NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          LocSurfaceForcing(
              LocFrictionVelocity, LocBuoyancyFlux, ICell, ConservTemp,
              AbsSalinity, LocPressureMid, LocSpecVol, ZonalStress, MeridStress,
              LatentHeatFluxEvap, SensibleHeatFlux, LongWaveHeatFluxUp,
              LongWaveHeatFluxDown, SeaIceHeatFlux, ShortWaveHeatFlux, SnowFlux,
              RainFlux, EvaporationFlux, SeaIceFreshWaterFlux, IceRunoffFlux,
              RiverRunoffFlux, SeaIceSaltFlux);
       });

   Halo *MeshHalo = Halo::getDefault();
   OMEGA_REQUIRE(MeshHalo, "KPPMix::update: null default Halo pointer");
   I4 HaloErr = 0;
   HaloErr += MeshHalo->exchangeFullArrayHalo(SurfaceFrictionVelocity, OnCell);
   HaloErr += MeshHalo->exchangeFullArrayHalo(SurfaceBuoyancyFlux, OnCell);
   OMEGA_REQUIRE(HaloErr == 0,
                 "KPPMix::update: halo exchange of surface forcing failed");

   computeKPPMix(PotentialDensity, NormalVelocity, TangentialVelocity,
                 SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                 EqState->BruntVaisalaFreqSq, ForcingState->IceFractionCell,
                 ForcingState->WindSpeed10mCell, InteriorVertDiff,
                 InteriorVertVisc);
}

/// Prepare KPP fields for the current time step from canonical ocean state.
void KPPMix::computeKPPFields(const OceanState *State,
                              const Array3DReal &TracerArray,
                              int ThickTimeLevel, int VelTimeLevel,
                              bool UseTracerForcing) {

   Eos *EqState = Eos::getInstance();
   if (!EqState || !Enabled)
      return;

   Pacer::start("KPP:computeKPPFields", 1);

   I4 TempIdx = -1;
   I4 SaltIdx = -1;
   if (Tracers::getIndex(TempIdx, "Temperature") != 0 ||
       Tracers::getIndex(SaltIdx, "Salinity") != 0) {
      LOG_WARN("KPPMix::computeKPPFields: Temperature/Salinity "
               "tracers not found, skipping KPP update");
      Pacer::stop("KPP:computeKPPFields", 1);
      return;
   }

   const auto *ForcingState = Forcing::getDefault();
   if (!ForcingState) {
      LOG_WARN("KPPMix::computeKPPFields: Forcing has not "
               "been initialized, skipping KPP update");
      Pacer::stop("KPP:computeKPPFields", 1);
      return;
   }

   // MatchBoth needs a non-KPP interior estimate to join its profile to;
   // this is the previous step's snapshot from VertMix (one-step lag).
   Array2DReal InteriorVertDiff;
   Array2DReal InteriorVertVisc;
   if (MatchTechnique == KPPMatchType::MatchBoth) {
      VertMix *VMixInstance = VertMix::getInstance();
      if (VMixInstance) {
         InteriorVertDiff = VMixInstance->InteriorVertDiff;
         InteriorVertVisc = VMixInstance->InteriorVertVisc;
      }
   }

   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);
   update(TracerArray, TempIdx, SaltIdx, NormalVelEdge, EqState, ForcingState,
          UseTracerForcing, InteriorVertDiff, InteriorVertVisc);

   Pacer::stop("KPP:computeKPPFields", 1);
}

/// Main computation routine
void KPPMix::computeKPPMix(const Array2DReal &PotentialDensity,
                           const Array2DReal &NormalVelocity,
                           const Array2DReal &TangentialVelocity,
                           const Array1DReal &SurfaceFrictionVelocity,
                           const Array1DReal &SurfaceBuoyancyFlux,
                           const Array2DReal &BruntVaisalaFreqSq,
                           const Array1DReal &IceFraction,
                           const Array1DReal &WindSpeed10m,
                           const Array2DReal &InteriorVertDiff,
                           const Array2DReal &InteriorVertVisc) {

   if (!Enabled) {
      return; // Skip if disabled
   }

   // Retain PotentialDensity for diagnostics/stream output.
   deepCopy(this->PotentialDensity, PotentialDensity);

   // =======================================================================
   // Stage 1: Compute OSBL depth
   // =======================================================================
   computeOSBLDepth(PotentialDensity, NormalVelocity, TangentialVelocity,
                    SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                    BruntVaisalaFreqSq, IceFraction, WindSpeed10m);

   // =======================================================================
   // Stage 2: Compute Mixing Coefficients
   // =======================================================================
   computeMixingCoefficients(SurfaceFrictionVelocity, SurfaceBuoyancyFlux,
                             InteriorVertDiff, InteriorVertVisc);

   if (DebugDiagnostics) {
      logDiagnostics(PotentialDensity, SurfaceFrictionVelocity,
                     SurfaceBuoyancyFlux, WindSpeed10m);
   }
}

void KPPMix::logDiagnostics(const Array2DReal &PotentialDensity,
                            const Array1DReal &SurfaceFrictionVelocity,
                            const Array1DReal &SurfaceBuoyancyFlux,
                            const Array1DReal &WindSpeed10m) {

   using namespace KPP;

   const auto MinLayerCellH = createHostMirrorCopy(VCoord->MinLayerCell);
   const auto MaxLayerCellH = createHostMirrorCopy(VCoord->MaxLayerCell);
   const auto ZInterfaceH   = createHostMirrorCopy(VCoord->GeomZInterface);
   const auto SshCellH      = createHostMirrorCopy(VCoord->SshCell);
   const auto DensityH      = createHostMirrorCopy(PotentialDensity);
   const auto UStarH        = createHostMirrorCopy(SurfaceFrictionVelocity);
   const auto B0H           = createHostMirrorCopy(SurfaceBuoyancyFlux);
   const auto OSBLDepthH    = createHostMirrorCopy(OSBLDepth);
   const auto OSBLIndexH    = createHostMirrorCopy(OSBLDepthIndex);
   const auto VertDiffH     = createHostMirrorCopy(VertDiff);
   const auto VertViscH     = createHostMirrorCopy(VertVisc);

   const int NCellsAll = Mesh->NCellsAll;
   if (NCellsAll <= 0) {
      return;
   }

   // Domain-wide diagnostic to avoid misleading single-cell checks.
   Real MaxAbsB0       = 0.0_Real;
   Real MaxAbsUStar    = 0.0_Real;
   Real MaxVertDiff    = 0.0_Real;
   Real MaxVertVisc    = 0.0_Real;
   int MaxAbsB0Cell    = -1;
   int MaxAbsUStarCell = -1;
   int MaxVertDiffCell = -1;
   int MaxVertDiffK    = -1;
   int MaxVertViscCell = -1;
   int MaxVertViscK    = -1;
   for (int C = 0; C < NCellsAll; ++C) {
      const Real AbsB0    = Kokkos::abs(B0H(C));
      const Real AbsUStar = Kokkos::abs(UStarH(C));
      if (AbsB0 > MaxAbsB0) {
         MaxAbsB0     = AbsB0;
         MaxAbsB0Cell = C;
      }
      if (AbsUStar > MaxAbsUStar) {
         MaxAbsUStar     = AbsUStar;
         MaxAbsUStarCell = C;
      }
      const int KCMin = MinLayerCellH(C);
      const int KCMax = MaxLayerCellH(C) + 1;
      for (int K = KCMin; K <= KCMax; ++K) {
         const Real Diff = VertDiffH(C, K);
         const Real Visc = VertViscH(C, K);
         if (Diff > MaxVertDiff) {
            MaxVertDiff     = Diff;
            MaxVertDiffCell = C;
            MaxVertDiffK    = K;
         }
         if (Visc > MaxVertVisc) {
            MaxVertVisc     = Visc;
            MaxVertViscCell = C;
            MaxVertViscK    = K;
         }
      }
   }

   const int ICell       = 0;
   const int KMin        = MinLayerCellH(ICell);
   const int KMax        = MaxLayerCellH(ICell);
   const int NVertLayers = VCoord->NVertLayers;

   if (KMin > KMax) {
      return;
   }

   const int KSurf     = Kokkos::min(KMin, NVertLayers - 1);
   const Real RhoSurf  = DensityH(ICell, KSurf);
   const Real UStar    = UStarH(ICell);
   const Real UStarEff = Kokkos::fmax(KPP::MinUStar, UStar);
   const Real BuoyFlux = B0H(ICell);
   Real Wind10m        = 0.0_Real;
   if (WindSpeed10m.extent(0) > 0) {
      const auto Wind10mH = createHostMirrorCopy(WindSpeed10m);
      Wind10m             = Wind10mH(ICell);
   }
   const Real LangmuirFactor =
       UseLangmuirTurbulence
           ? computeLangmuirEnhancement(Wind10m, UStarEff, 50.0)
           : 1.0_Real;
   const Real BuoyFluxEff = BuoyFlux * LangmuirFactor;

   const int KOblIface = Kokkos::min(
       NVertLayers, Kokkos::max(KMin, static_cast<int>(OSBLIndexH(ICell)) + 1));

   const int KTop   = Kokkos::min(KMax, KMin + 3);
   const int KOSBL  = OSBLIndexH(ICell);
   const Real HOSBL = OSBLDepthH(ICell);

   for (int K = KMin; K <= KTop; ++K) {
      const int KCell   = Kokkos::min(K, NVertLayers - 1);
      const int KIface  = Kokkos::min(K + 1, NVertLayers);
      const Real ZDepth = SshCellH(ICell) - ZInterfaceH(ICell, KIface);

      const Real RhoK     = DensityH(ICell, KCell);
      const Real DeltaRho = RhoK - RhoSurf;
      const Real DeltaB   = Gravity * DeltaRho / RhoSw;
      const Real WTurb =
          computeTurbVelocityScale(UStarEff, BuoyFluxEff, ZDepth);
      const Real RiBulk =
          DeltaB * ZDepth / (WTurb * WTurb + KPP::NumericalTolerance);

      Real Sigma = 0.0_Real;
      if (K <= KOSBL) {
         Sigma = -1.0_Real * static_cast<Real>(K - KMin) /
                 static_cast<Real>(KOSBL - KMin + 1);
         Sigma = Kokkos::fmax(-1.0_Real, Kokkos::fmin(0.0_Real, Sigma));
      }

      // Monin-Obukhov coordinate zeta = d/L at this depth
      const Real ZLocal = -Sigma * HOSBL;
      Real Zeta         = 0.0_Real;
      const Real Denom  = VonKar * BuoyFlux;
      if (Kokkos::abs(Denom) > 1.0e-16_Real) {
         const Real LMoninObukhov = (UStarEff * UStarEff * UStarEff) / Denom;
         if (Kokkos::abs(LMoninObukhov) > 1.0e-16_Real) {
            Zeta = ZLocal / LMoninObukhov;
         }
      }

      const Real PhiInvM = KPP::kppPhiInvMomentum(Zeta);
      const Real PhiInvS = KPP::kppPhiInvScalar(Zeta);
   }
}

/// Stage 1: Compute OSBL depth using bulk Richardson search with edge-based
/// velocity shear following the MPAS CVMix reference implementation
/// (mpas_ocn_vmix_cvmix.F)
void KPPMix::computeOSBLDepth(const Array2DReal &PotentialDensity,
                              const Array2DReal &NormalVelocity,
                              const Array2DReal &TangentialVelocity,
                              const Array1DReal &SurfaceFrictionVelocity,
                              const Array1DReal &SurfaceBuoyancyFlux,
                              const Array2DReal &BruntVaisalaFreqSq,
                              const Array1DReal &IceFraction,
                              const Array1DReal &WindSpeed10m) {

   const I4 NVertLayers = VCoord->NVertLayers;

   // =======================================================================
   // Compute Langmuir enhancement factors if wind speed is available
   // =======================================================================
   KPPLangmuirFactor LangmuirCalc;
   LangmuirCalc.UseLangmuirTurbulence       = UseLangmuirTurbulence;
   LangmuirCalc.IceFracThresholdForLangmuir = IceFractionThresholdForLangmuir;

   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocIceFraction, IceFraction);
   OMEGA_SCOPE(LocSurfFricVel, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfBuoyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(LocWindSpeed10m, WindSpeed10m);

   parallelFor(
       "KPP-Langmuir", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          LangmuirCalc(LocLangmuirFactor, ICell, LocIceFraction, LocSurfFricVel,
                       LocWindSpeed10m);
       });

   // =======================================================================
   // Stage 1: Compute OSBL depth using edge-based velocity shear
   // =======================================================================
   KPPOSBLDepthSearch OSBLDepthCalc(Mesh, VCoord);
   OSBLDepthCalc.CriticalRichardson = CriticalRichardson;
   OSBLDepthCalc.SurfaceLayerExtent = SurfaceLayerExtent;
   OSBLDepthCalc.IceFracThresholdForMinimumOSBL =
       IceFractionThresholdForMinimumOSBL;
   OSBLDepthCalc.MinimumOSBLUnderSeaIce = MinimumOSBLUnderSeaIce;
   OSBLDepthCalc.FullRiProfile          = DebugDiagnostics;

   OMEGA_SCOPE(LocPotentialDensity, PotentialDensity);
   OMEGA_SCOPE(LocNormalVelocity, NormalVelocity);
   OMEGA_SCOPE(LocTangentialVelocity, TangentialVelocity);
   OMEGA_SCOPE(LocBruntVaisalaFreqSq, BruntVaisalaFreqSq);
   OMEGA_SCOPE(LocOSBLDepth, OSBLDepth);
   OMEGA_SCOPE(LocOSBLDepthIndex, OSBLDepthIndex);
   OMEGA_SCOPE(LocBulkRichardson, BulkRichardsonNumber);
   OMEGA_SCOPE(LocBulkRichardsonShear, BulkRichardsonShear);
   OMEGA_SCOPE(LocUnresolvedShear, UnresolvedShear);
   OMEGA_SCOPE(LocBuoyancyJump, BuoyancyJump);
   OMEGA_SCOPE(LocMinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, VCoord->MaxLayerCell);

   parallelFor(
       "KPP-OSBLDiagnostics-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const I4 KMin = LocMinLayerCell(ICell);
          const I4 KMax = LocMaxLayerCell(ICell);
          if (K >= KMin && K <= KMax + 1) {
             LocBulkRichardson(ICell, K)      = 0.0_Real;
             LocBulkRichardsonShear(ICell, K) = 0.0_Real;
             LocUnresolvedShear(ICell, K)     = 0.0_Real;
             LocBuoyancyJump(ICell, K)        = 0.0_Real;
          }
       });

   parallelFor(
       "KPP-OSBLDepth", {Mesh->NCellsOwned}, KOKKOS_LAMBDA(I4 ICell) {
          OSBLDepthCalc(LocOSBLDepth, LocOSBLDepthIndex, LocBulkRichardson,
                        LocBulkRichardsonShear, LocUnresolvedShear,
                        LocBuoyancyJump, ICell, LocPotentialDensity,
                        LocNormalVelocity, LocTangentialVelocity,
                        LocSurfFricVel, LocSurfBuoyFlux, LocBruntVaisalaFreqSq,
                        LocIceFraction, LocLangmuirFactor);
       });

   // Exchange OSBLDepth (and OSBLDepthIndex when smoothing is disabled) so
   // halo cells hold exact unsmoothed values prior to Laplacian smoothing.
   Halo *MeshHalo = Halo::getDefault();
   OMEGA_REQUIRE(MeshHalo,
                 "KPPMix::computeOSBLDepth: null default Halo pointer");
   I4 HaloErr = MeshHalo->exchangeFullArrayHalo(OSBLDepth, OnCell);
   if (!UseOSBLSmoothing) {
      HaloErr += MeshHalo->exchangeFullArrayHalo(OSBLDepthIndex, OnCell);
   }
   OMEGA_REQUIRE(HaloErr == 0,
                 "KPPMix::computeOSBLDepth: halo exchange of OSBLDepth failed");

   if (UseOSBLSmoothing) {
      OMEGA_SCOPE(LocOSBLDepthSmooth, OSBLDepthSmooth);

      KPPOSBLSmooth BLDSmoothCalc(Mesh, VCoord);
      KPPOSBLCommit BLDCommitCalc(Mesh, VCoord);

      parallelFor(
          "KPP-OSBLDepth-Smooth", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
             BLDSmoothCalc(LocOSBLDepthSmooth, ICell, LocOSBLDepth);
          });

      parallelFor(
          "KPP-OSBLDepth-CommitSmooth", {Mesh->NCellsAll},
          KOKKOS_LAMBDA(I4 ICell) {
             BLDCommitCalc(LocOSBLDepth, LocOSBLDepthIndex, ICell,
                           LocOSBLDepthSmooth);
          });
   }
}

/// Stage 2: Compute KPP mixing contribution or matched coefficients
void KPPMix::computeMixingCoefficients(
    const Array1DReal &SurfaceFrictionVelocity,
    const Array1DReal &SurfaceBuoyancyFlux, const Array2DReal &InteriorVertDiff,
    const Array2DReal &InteriorVertVisc) {

   I4 NVertLayers = VCoord->NVertLayers;

   const bool LocUseInteriorMix =
       InteriorVertDiff.data() != nullptr && InteriorVertVisc.data() != nullptr;
   const bool LocUseMatchedShapes =
       LocUseInteriorMix && MatchTechnique == KPPMatchType::MatchBoth;

   OMEGA_SCOPE(LocOSBLDepth, OSBLDepth);
   OMEGA_SCOPE(LocOSBLDepthIndex, OSBLDepthIndex);
   OMEGA_SCOPE(LocVertDiff, VertDiff);
   OMEGA_SCOPE(LocVertVisc, VertVisc);
   OMEGA_SCOPE(LocVertNonLocalFlux, VertNonLocalFlux);
   OMEGA_SCOPE(LocTurbulentVelocityScale, TurbulentVelocityScale);
   OMEGA_SCOPE(LocSurfFricVel, SurfaceFrictionVelocity);
   OMEGA_SCOPE(LocSurfBuoyFlux, SurfaceBuoyancyFlux);
   OMEGA_SCOPE(LocLangmuirFactor, LangmuirFactor);
   OMEGA_SCOPE(LocInteriorVertDiff, InteriorVertDiff);
   OMEGA_SCOPE(LocInteriorVertVisc, InteriorVertVisc);
   OMEGA_SCOPE(LocMinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, VCoord->MaxLayerCell);

   // =======================================================================
   // Initialize with zero KPP contribution, or precomputed interior mixing for
   // matched-coefficient construction.
   // =======================================================================
   KPPCoeffsInit CoeffsInit;
   CoeffsInit.UseInteriorMix = LocUseInteriorMix;

   parallelFor(
       "KPP-Coeffs-Init", {Mesh->NCellsAll, NVertLayers + 1},
       KOKKOS_LAMBDA(I4 ICell, I4 K) {
          const I4 KMin = LocMinLayerCell(ICell);
          const I4 KMax = LocMaxLayerCell(ICell);
          const bool IsValidCell =
              KMin >= 0 && KMin < NVertLayers && KMax >= KMin;
          if (IsValidCell && K >= KMin && K <= KMax + 1) {
             CoeffsInit(LocVertDiff, LocVertVisc, LocVertNonLocalFlux,
                        LocTurbulentVelocityScale, ICell, K,
                        LocInteriorVertDiff, LocInteriorVertVisc);
          }
       });

   // =======================================================================
   // Stage 2: Compute KPP profile-based mixing coefficients
   // =======================================================================
   KPPMixingCoeffs MixingCoeffsCalc(Mesh, VCoord);
   MixingCoeffsCalc.SurfaceLayerExtent   = SurfaceLayerExtent;
   MixingCoeffsCalc.UseEnhancedDiffusion = UseEnhancedDiffusion;
   MixingCoeffsCalc.UseInteriorMix       = LocUseInteriorMix;
   MixingCoeffsCalc.UseMatchedShapes     = LocUseMatchedShapes;
   MixingCoeffsCalc.NonLocalCs = KPP::kppNonLocalCs(VonKar, SurfaceLayerExtent);

   parallelFor(
       "KPP-MixingCoeffs", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          MixingCoeffsCalc(LocVertDiff, LocVertVisc, LocVertNonLocalFlux,
                           LocTurbulentVelocityScale, ICell, LocOSBLDepth,
                           LocOSBLDepthIndex, LocSurfFricVel, LocSurfBuoyFlux,
                           LocLangmuirFactor, LocInteriorVertDiff,
                           LocInteriorVertVisc);
       });
}

/// Constructor for KPPOSBLDepthSearch
KPPOSBLDepthSearch::KPPOSBLDepthSearch(const HorzMesh *Mesh,
                                       const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      ZInterface(VCoord->GeomZInterface), ZMid(VCoord->GeomZMid),
      SshCell(VCoord->SshCell), NEdgesOnCell(Mesh->NEdgesOnCell),
      EdgesOnCell(Mesh->EdgesOnCell), CellsOnCell(Mesh->CellsOnCell),
      AreaCell(Mesh->AreaCell), DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge) {}

/// Constructor for KPPOSBLSmooth
KPPOSBLSmooth::KPPOSBLSmooth(const HorzMesh *Mesh, const VertCoord *)
    : NCellsAll(Mesh->NCellsAll), NEdgesOnCell(Mesh->NEdgesOnCell),
      CellsOnCell(Mesh->CellsOnCell), AreaCell(Mesh->AreaCell) {}

/// Constructor for KPPOSBLCommit
KPPOSBLCommit::KPPOSBLCommit(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell), ZInterface(VCoord->GeomZInterface),
      ZMid(VCoord->GeomZMid), SshCell(VCoord->SshCell) {}

/// Constructor for KPPMixingCoeffs
KPPMixingCoeffs::KPPMixingCoeffs(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NVertLayers(VCoord->NVertLayers), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell), ZInterface(VCoord->GeomZInterface),
      ZMid(VCoord->GeomZMid), SshCell(VCoord->SshCell) {}

/// Register fields with I/O system
void KPPMix::defineFields() {
   // OSBLDepth on cells
   std::vector<std::string> CellDims(1);
   CellDims[0] = "NCells";
   auto OSBLDepthField =
       Field::create(OSBLDepthFldName,                     // field name
                     "ocean surface boundary layer depth", // long name
                     "m",                                  // units
                     "",                                   // CF standard name
                     0.0,                                  // min valid value
                     std::numeric_limits<Real>::max(),     // max valid value
                     1,                                    // number of dims
                     CellDims);

   auto OSBLDepthIndexField =
       Field::create(OSBLDepthIndexFldName,                      // field name
                     "ocean surface boundary layer depth index", // long name
                     "",                                         // units
                     "",                             // CF standard name
                     -1,                             // min valid value
                     std::numeric_limits<I4>::max(), // max valid value
                     1,                              // number of dims
                     CellDims);

   // KPP non-local tracer flux profile on cell-layer interfaces
   std::vector<std::string> FluxDims(2);
   FluxDims[0] = "NCells";
   FluxDims[1] = "NVertLayersP1";
   auto NonLocalFluxField =
       Field::create(NonLocalFluxFldName,                 // field name
                     "KPP non-local tracer flux profile", // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     FluxDims);

   auto BulkRichardsonField =
       Field::create(BulkRichardsonFldName,               // field name
                     "bulk Richardson number",            // long name
                     "1",                                 // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     FluxDims);

   auto BulkRichardsonShearField =
       Field::create(BulkRichardsonShearFldName,       // field name
                     "bulk Richardson shear term",     // long name
                     "m2 s-2",                         // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   auto UnresolvedShearField =
       Field::create(UnresolvedShearFldName,           // field name
                     "KPP unresolved shear term Vt2",  // long name
                     "m2 s-2",                         // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   auto BuoyancyJumpField =
       Field::create(BuoyancyJumpFldName,                   // field name
                     "KPP buoyancy jump (density anomaly)", // long name
                     "m s-2",                               // units
                     "",                                    // CF standard name
                     std::numeric_limits<Real>::lowest(),   // min valid value
                     std::numeric_limits<Real>::max(),      // max valid value
                     2,                                     // number of dims
                     FluxDims);

   auto TurbulentVelScaleField =
       Field::create(TurbulentVelScaleFldName,         // field name
                     "KPP turbulent velocity scale",   // long name
                     "m s-1",                          // units
                     "",                               // CF standard name
                     0.0,                              // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     2,                                // number of dims
                     FluxDims);

   std::vector<std::string> LayerDims(2);
   LayerDims[0] = "NCells";
   LayerDims[1] = "NVertLayers";
   auto PotentialDensityField =
       Field::create(PotentialDensityFldName,             // field name
                     "KPP potential density",             // long name
                     "kg m-3",                            // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     2,                                   // number of dims
                     LayerDims);

   // Group KPP-specific outputs for convenient stream selection.
   auto KPPGroup = FieldGroup::create("KPPMix");
   KPPGroup->addField(OSBLDepthFldName);
   KPPGroup->addField(OSBLDepthIndexFldName);
   KPPGroup->addField(NonLocalFluxFldName);
   KPPGroup->addField(BulkRichardsonFldName);
   KPPGroup->addField(BulkRichardsonShearFldName);
   KPPGroup->addField(UnresolvedShearFldName);
   KPPGroup->addField(BuoyancyJumpFldName);
   KPPGroup->addField(TurbulentVelScaleFldName);
   KPPGroup->addField(PotentialDensityFldName);
   KPPGroup->addField(SurfFricVelFldName);
   KPPGroup->addField(SurfBuoyFluxFldName);

   OSBLDepthField->attachData<Array1DReal>(OSBLDepth);
   OSBLDepthIndexField->attachData<Array1DI4>(OSBLDepthIndex);
   NonLocalFluxField->attachData<Array2DReal>(VertNonLocalFlux);
   BulkRichardsonField->attachData<Array2DReal>(BulkRichardsonNumber);
   BulkRichardsonShearField->attachData<Array2DReal>(BulkRichardsonShear);
   UnresolvedShearField->attachData<Array2DReal>(UnresolvedShear);
   BuoyancyJumpField->attachData<Array2DReal>(BuoyancyJump);
   TurbulentVelScaleField->attachData<Array2DReal>(TurbulentVelocityScale);
   PotentialDensityField->attachData<Array2DReal>(PotentialDensity, false);

   // Surface forcing fields on cells
   auto SurfFricVelField =
       Field::create(SurfFricVelFldName,                 // field name
                     "KPP surface friction velocity u*", // long name
                     "m s-1",                            // units
                     "",                                 // CF standard name
                     0.0,                                // min valid value
                     std::numeric_limits<Real>::max(),   // max valid value
                     1,                                  // number of dims
                     CellDims);
   SurfFricVelField->attachData<Array1DReal>(SurfaceFrictionVelocity, false);

   auto SurfBuoyFluxField =
       Field::create(SurfBuoyFluxFldName,                 // field name
                     "KPP surface buoyancy flux",         // long name
                     "m2 s-3",                            // units
                     "",                                  // CF standard name
                     std::numeric_limits<Real>::lowest(), // min valid value
                     std::numeric_limits<Real>::max(),    // max valid value
                     1,                                   // number of dims
                     CellDims);
   SurfBuoyFluxField->attachData<Array1DReal>(SurfaceBuoyancyFlux, false);

   auto LangmuirFactorField =
       Field::create(LangmuirFactorFldName,             // field name
                     "KPP Langmuir enhancement factor", // long name
                     "1",                               // units
                     "",                                // CF standard name
                     1.0,                               // min valid value
                     2.0,                               // max valid value
                     1,                                 // number of dims
                     CellDims);
   LangmuirFactorField->attachData<Array1DReal>(LangmuirFactor, false);

   KPPGroup->addField(LangmuirFactorFldName);

   // OSBLDepth, OSBLDepthIndex, NonLocalFlux, BulkRichardson,
   // BulkRichardsonShear, UnresolvedShear, BuoyancyJump, and
   // TurbulentVelScale attached with FillOnAttach=true above, which already
   // records _FillValue metadata; only fields attached with
   // FillOnAttach=false need it set explicitly here.
   PotentialDensityField->addMetadata("_FillValue", FillValueReal);
   SurfFricVelField->addMetadata("_FillValue", FillValueReal);
   SurfBuoyFluxField->addMetadata("_FillValue", FillValueReal);
   LangmuirFactorField->addMetadata("_FillValue", FillValueReal);

   LOG_INFO("KPPMix::defineFields: registered {}, {}, {}, {}, {}, {}, {}, {}, "
            "{}, {}, {}, {}",
            OSBLDepthFldName, OSBLDepthIndexFldName, NonLocalFluxFldName,
            BulkRichardsonFldName, BulkRichardsonShearFldName,
            UnresolvedShearFldName, BuoyancyJumpFldName,
            TurbulentVelScaleFldName, PotentialDensityFldName,
            SurfFricVelFldName, SurfBuoyFluxFldName, LangmuirFactorFldName);
}

} // namespace OMEGA
