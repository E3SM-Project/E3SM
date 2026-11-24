//===-- SplitExplicitInit.cpp - split-explicit initialization --*- C++ -*-===//
//
// Utilities for initializing split-explicit time stepping state.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitInit.h"
#include "Config.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"

#include <algorithm>
#include <cmath>

namespace OMEGA {

//------------------------------------------------------------------------------
SplitExplicitConfig
SplitExplicitInit::readConfigOptions(const TimeInterval &TimeStep,
                                     bool IsUnsplit) {

   SplitExplicitConfig SEConfig;

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TimeIntConfig("TimeIntegration");
   Error Err = OmegaConfig->get(TimeIntConfig);
   CHECK_ERROR_ABORT(Err, "TimeIntegration group not found in Config");

   Config ModeSplitShareConfig("ModeSplitShare");
   Error ShareErr = TimeIntConfig.get(ModeSplitShareConfig);
   CHECK_ERROR_ABORT(ShareErr,
                     "TimeIntegration::ModeSplitShare group not found in "
                     "Config; it is required by the split time steppers");

   // For the unsplit schemes, SplitFactor is 0.
   if (IsUnsplit) {
      SEConfig.SplitFactor = 0._Real;
   }

   if (!IsUnsplit) {
      std::string BtrTimeStepperStr;
      if (ModeSplitShareConfig.get("BtrTimeStepper", BtrTimeStepperStr)
              .isSuccess()) {
         SEConfig.BtrTimeStepper = getBtrTimeStepperFromStr(BtrTimeStepperStr);
      }

      std::string BtrTimeStepStr;
      Error BtrErr = ModeSplitShareConfig.get("BtrTimeStep", BtrTimeStepStr);
      CHECK_ERROR_ABORT(BtrErr,
                        "TimeIntegration::ModeSplitShare::BtrTimeStep not "
                        "found in Config; it is required by SplitExplicitRK2");
      SEConfig.BtrTimeStep = TimeInterval(BtrTimeStepStr);
   }

   I4 NTimeStepIteration = SEConfig.NTimeStepIteration;
   if (ModeSplitShareConfig.get("NTimeStepIteration", NTimeStepIteration)
           .isSuccess()) {
      if (NTimeStepIteration < 1) {
         ABORT_ERROR("NTimeStepIteration must be greater than zero");
      }
      SEConfig.NTimeStepIteration = NTimeStepIteration;
   }

   I4 NBclCoriolisIteration = SEConfig.NBclCoriolisIteration;
   if (ModeSplitShareConfig.get("NBclCoriolisIteration", NBclCoriolisIteration)
           .isSuccess()) {
      if (NBclCoriolisIteration < 1) {
         ABORT_ERROR("NBclCoriolisIteration must be greater than zero");
      }
      SEConfig.NBclCoriolisIteration = NBclCoriolisIteration;
   }

   ModeSplitShareConfig.get("ReinitSplitVelocity",
                            SEConfig.ReinitSplitVelocity);

   // For the unsplit schemes, the linear Coriolis iteration is 1.
   if (IsUnsplit) {
      SEConfig.NBclCoriolisIteration = 1;
   }

   // Report the resolved options
   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   if (IsUnsplit) {
      LOG_INFO("UnsplitRK2: NTimeStepIteration={}, TimeStep={} s",
               SEConfig.NTimeStepIteration, TimeStepSeconds);
   } else {
      SEConfig.NBtrSubcycles =
          computeSubcycleCount(TimeStep, SEConfig.BtrTimeStep);

      R8 BtrTimeStepSeconds;
      SEConfig.BtrTimeStep.get(BtrTimeStepSeconds, TimeUnits::Seconds);
      LOG_INFO("SplitExplicitRK2: TimeStep={} s, BtrTimeStep={} s, "
               "NBtrSubcycles={}, BtrDt={} s, NTimeStepIteration={}, "
               "NBclCoriolisIteration={}",
               TimeStepSeconds, BtrTimeStepSeconds, SEConfig.NBtrSubcycles,
               TimeStepSeconds / SEConfig.NBtrSubcycles,
               SEConfig.NTimeStepIteration, SEConfig.NBclCoriolisIteration);
   }

   return SEConfig;
}

//------------------------------------------------------------------------------
SplitExplicitBarotropicStepperType
SplitExplicitInit::getBtrTimeStepperFromStr(const std::string &InString) {

   if (InString == "Predictor-Corrector") {
      return SplitExplicitBarotropicStepperType::PredictorCorrector;
   }

   ABORT_ERROR("BtrTimeStepper should be 'Predictor-Corrector' but got {}",
               InString);
   return SplitExplicitBarotropicStepperType::Invalid;
}

//------------------------------------------------------------------------------
I4 SplitExplicitInit::computeSubcycleCount(const TimeInterval &TimeStep,
                                           const TimeInterval &BtrTimeStep) {

   R8 TimeStepSeconds;
   R8 BtrTimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   BtrTimeStep.get(BtrTimeStepSeconds, TimeUnits::Seconds);

   if (BtrTimeStepSeconds <= 0.) {
      ABORT_ERROR("BtrTimeStep must be greater than zero");
   }

   if (BtrTimeStepSeconds >= TimeStepSeconds) {
      ABORT_ERROR("BtrTimeStep ({} s) must be smaller than the baroclinic "
                  "TimeStep ({} s) so that the barotropic mode is subcycled",
                  BtrTimeStepSeconds, TimeStepSeconds);
   }

   return std::max<I4>(
       1, static_cast<I4>(std::ceil(TimeStepSeconds / BtrTimeStepSeconds)));
}

//------------------------------------------------------------------------------
void SplitExplicitInit::allocateScratch(SplitExplicitScratch &Scratch,
                                        const HorzMesh *Mesh,
                                        const VertCoord *VCoord,
                                        const std::string &Name) {

   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");

   Scratch.NormalBarotropicVelocitySubcycleCur = Array1DReal(
       "NormalBarotropicVelocitySubcycleCur" + Name, Mesh->NEdgesSize);
   Scratch.NormalBarotropicVelocitySubcyclePre = Array1DReal(
       "NormalBarotropicVelocitySubcyclePre" + Name, Mesh->NEdgesSize);
   Scratch.NormalBarotropicVelocitySubcycleCor = Array1DReal(
       "NormalBarotropicVelocitySubcycleCor" + Name, Mesh->NEdgesSize);
   Scratch.BarotropicPressureAnomalySubcycleCur = Array1DReal(
       "BarotropicPressureAnomalySubcycleCur" + Name, Mesh->NCellsSize);
   Scratch.BarotropicPressureAnomalySubcyclePre = Array1DReal(
       "BarotropicPressureAnomalySubcyclePre" + Name, Mesh->NCellsSize);
   Scratch.BarotropicPressureAnomalySubcycleCor = Array1DReal(
       "BarotropicPressureAnomalySubcycleCor" + Name, Mesh->NCellsSize);
   Scratch.BarotropicForcing =
       Array1DReal("BarotropicForcing" + Name, Mesh->NEdgesSize);
   Scratch.BarotropicFlux =
       Array1DReal("BarotropicFlux" + Name, Mesh->NEdgesSize);
   Scratch.BaroclinicPseudoThicknessEdge =
       Array1DReal("BaroclinicPseudoThicknessEdge" + Name, Mesh->NEdgesSize);
   Scratch.IterVelocityTend = Array2DReal(
       "IterVelocityTend" + Name, Mesh->NEdgesSize, VCoord->NVertLayers);
   Scratch.TransportVelocityAdd = Array2DReal(
       "TransportVelocityAdd" + Name, Mesh->NEdgesSize, VCoord->NVertLayers);

   parallelFor(
       "initializeCell1D", {Mesh->NCellsAll}, KOKKOS_LAMBDA(I4 ICell) {
          Scratch.BarotropicPressureAnomalySubcycleCur(ICell) = 0._Real;
          Scratch.BarotropicPressureAnomalySubcyclePre(ICell) = 0._Real;
          Scratch.BarotropicPressureAnomalySubcycleCor(ICell) = 0._Real;
       });

   parallelFor(
       "initializeEdge1D", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(I4 IEdge) {
          Scratch.NormalBarotropicVelocitySubcycleCur(IEdge) = 0._Real;
          Scratch.NormalBarotropicVelocitySubcyclePre(IEdge) = 0._Real;
          Scratch.NormalBarotropicVelocitySubcycleCor(IEdge) = 0._Real;
          Scratch.BarotropicForcing(IEdge)                   = 0._Real;
          Scratch.BarotropicFlux(IEdge)                      = 0._Real;
          Scratch.BaroclinicPseudoThicknessEdge(IEdge)       = 0._Real;
       });

   deepCopy(Scratch.IterVelocityTend, 0.);
   deepCopy(Scratch.TransportVelocityAdd, 0.);
}

//------------------------------------------------------------------------------
void SplitExplicitInit::computeVelocitySplit(OceanState *State,
                                             const HorzMesh *Mesh,
                                             const VertCoord *VCoord,
                                             I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocity  = State->getNormalVelocity(TimeLevel);
   Array2DReal PseudoThickness = State->getPseudoThickness(TimeLevel);
   Array2DReal NormalBaroclinicVelocity =
       State->getNormalBaroclinicVelocity(TimeLevel);
   Array1DReal NormalBarotropicVelocity =
       State->getNormalBarotropicVelocity(TimeLevel);

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(MinLayerEdgeTop, VCoord->MinLayerEdgeTop);
   OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);

   deepCopy(NormalBaroclinicVelocity, 0.);
   deepCopy(NormalBarotropicVelocity, 0.);

   parallelForOuter(
       "SplitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);
          const I4 KTop = MinLayerEdgeTop(IEdge);
          const I4 KBot = MaxLayerEdgeBot(IEdge);

          Real BarotropicVelocity = 0._Real;

          if (KMax >= KMin) {
             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             Real ThicknessSum = 0._Real;
             Real FluxSum      = 0._Real;

             parallelReduceInner(
                 Team, Range{KMin, KMax},
                 INNER_LAMBDA(const int K, Real &ThickAccum, Real &FluxAccum) {
                    const Real ThickEdge =
                        0.5_Real *
                        (PseudoThickness(Cell0, K) + PseudoThickness(Cell1, K));

                    ThickAccum += ThickEdge;
                    FluxAccum += ThickEdge * NormalVelocity(IEdge, K);
                 },
                 ThicknessSum, FluxSum);

             BarotropicVelocity = FluxSum / ThicknessSum;

             Kokkos::single(
                 PerTeam(Team), INNER_LAMBDA() {
                    NormalBarotropicVelocity(IEdge) =
                        BarotropicVelocity * EdgeMask(IEdge, KMin);
                 });

          } else {

             Kokkos::single(
                 PerTeam(Team),
                 INNER_LAMBDA() { NormalBarotropicVelocity(IEdge) = 0._Real; });
          }

          parallelForInner(
              Team, Range{KTop, KBot}, INNER_LAMBDA(const int K) {
                 NormalBaroclinicVelocity(IEdge, K) =
                     NormalVelocity(IEdge, K) - BarotropicVelocity;
              });
       });
}

//------------------------------------------------------------------------------
void SplitExplicitInit::computeUnsplitVelocitySplit(OceanState *State,
                                                    const HorzMesh *Mesh,
                                                    const VertCoord *VCoord,
                                                    I4 CurLevel, I4 NextLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");

   Array2DReal NormalVelocityCur = State->getNormalVelocity(CurLevel);
   Array2DReal NormalVelocityNew = State->getNormalVelocity(NextLevel);
   Array2DReal NormalBaroclinicVelocityCur =
       State->getNormalBaroclinicVelocity(CurLevel);
   Array2DReal NormalBaroclinicVelocityNew =
       State->getNormalBaroclinicVelocity(NextLevel);
   Array1DReal NormalBarotropicVelocityCur =
       State->getNormalBarotropicVelocity(CurLevel);
   Array1DReal NormalBarotropicVelocityNew =
       State->getNormalBarotropicVelocity(NextLevel);

   deepCopy(NormalBaroclinicVelocityCur, 0.);
   deepCopy(NormalBaroclinicVelocityNew, 0.);

   deepCopy(NormalBarotropicVelocityCur, 0.);
   deepCopy(NormalBarotropicVelocityNew, 0.);

   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   parallelForOuter(
       "UnsplitVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(I4 K) {
                 NormalBaroclinicVelocityCur(IEdge, K) =
                     NormalVelocityCur(IEdge, K);
                 NormalBaroclinicVelocityNew(IEdge, K) =
                     NormalVelocityNew(IEdge, K);
              });
       });
}

//------------------------------------------------------------------------------
void SplitExplicitInit::initializeBarotropicPressure(
    SplitExplicitScratch &Scratch, OceanState *State, const HorzMesh *Mesh,
    const VertCoord *VCoord, I4 TimeLevel) {

   if (!State)
      LOG_CRITICAL("Invalid State");
   if (!Mesh)
      LOG_CRITICAL("Invalid mesh");
   if (!VCoord)
      LOG_CRITICAL("Invalid vertical coordinate");

   Array1DReal BtrPressAnomaly = State->getBarotropicPressureAnomaly(TimeLevel);
   Array1DReal SurfacePressure = VCoord->SurfacePressure;
   Array2DReal PressureInterface = VCoord->PressureInterface;
   Array1DReal BottomGeomDepth   = VCoord->BottomGeomDepth;
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   parallelFor(
       "initializeBarotropicPressure", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell) {
          const I4 KMax = MaxLayerCell(ICell);

          const Real Pressure =
              PressureInterface(ICell, KMax + 1) - SurfacePressure(ICell);
          BtrPressAnomaly(ICell) =
              Pressure - RhoSw * Gravity * BottomGeomDepth(ICell);
          Scratch.BarotropicPressureAnomalySubcycleCur(ICell) =
              BtrPressAnomaly(ICell);
          Scratch.BarotropicPressureAnomalySubcyclePre(ICell) =
              BtrPressAnomaly(ICell);
          Scratch.BarotropicPressureAnomalySubcycleCor(ICell) =
              BtrPressAnomaly(ICell);
       });
}

} // namespace OMEGA
