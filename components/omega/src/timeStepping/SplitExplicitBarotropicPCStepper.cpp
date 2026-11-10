//===-- SplitExplicitBarotropicPCStepper.cpp - SE stage 2 -----*- C++ -*-===//
//
// Framework for the explicitly subcycled forward-backward
// predictor-corrector barotropic velocity update.
//
//===----------------------------------------------------------------------===//

#include "SplitExplicitBarotropicPCStepper.h"
#include "AuxiliaryState.h"
#include "Eos.h"
#include "GlobalConstants.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"

#include <utility>

namespace OMEGA {

namespace {

/// Blends an old and a new subcycle value with the forward-backward feedback
/// weight Gamma, i.e. (1-Gamma)*Old + Gamma*New.
KOKKOS_INLINE_FUNCTION Real feedbackBlend(Real Gamma, Real Old, Real New) {
   return (1._Real - Gamma) * Old + Gamma * New;
}

/// Edge value of the barotropic pressure deviation from the reference column,
/// centered or upwinded to match the flux-thickness choice used elsewhere.
KOKKOS_INLINE_FUNCTION Real deltaBtrPressureEdge(Real Delta0, Real Delta1,
                                                 Real NormalBtrVelEdge,
                                                 bool Upwind) {
   if (!Upwind)
      return 0.5_Real * (Delta0 + Delta1);
   if (NormalBtrVelEdge > 0._Real)
      return Delta0;
   if (NormalBtrVelEdge < 0._Real)
      return Delta1;
   return Kokkos::max(Delta0, Delta1);
}

} // anonymous namespace

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::init(const AuxiliaryState *InAuxState,
                                            SplitExplicitScratch *InScratch,
                                            const SplitExplicitConfig *InConfig,
                                            const HorzMesh *InMesh,
                                            Halo *InMeshHalo,
                                            const VertCoord *InVCoord) {

   if (!InAuxState)
      LOG_CRITICAL("Invalid auxiliary state");
   if (!InScratch)
      LOG_CRITICAL("Invalid split-explicit scratch data");
   if (!InConfig)
      LOG_CRITICAL("Invalid split-explicit config");
   if (!InMesh)
      LOG_CRITICAL("Invalid mesh");
   if (!InMeshHalo)
      LOG_CRITICAL("Invalid MeshHalo");
   if (!InVCoord)
      LOG_CRITICAL("Invalid vertical coordinate");
   if (InConfig->NBtrSubcycles < 1)
      LOG_CRITICAL("Invalid split-explicit barotropic subcycle count");

   AuxState = InAuxState;
   Scratch  = InScratch;
   SEConfig = InConfig;
   Mesh     = InMesh;
   MeshHalo = InMeshHalo;
   VCoord   = InVCoord;
}

//------------------------------------------------------------------------------
void SplitExplicitBarotropicPCStepper::doBarotropicVelocityUpdate(
    OceanState *State, I4 CurLevel, I4 NextLevel,
    const TimeInterval &StageTimeStep) const {

   if (!Scratch)
      LOG_CRITICAL("Split-explicit barotropic time stepper not initialized");

   Eos *EqState = Eos::getInstance();

   R8 StageTimeStepSeconds;
   StageTimeStep.get(StageTimeStepSeconds, TimeUnits::Seconds);

   const I4 NBtrSubcycles = 2 * SEConfig->NBtrSubcycles;
   const Real BtrDt       = StageTimeStepSeconds / SEConfig->NBtrSubcycles;
   const Real InvBtrVelAvgCount  = 1._Real / (NBtrSubcycles + 1);
   const Real InvBtrFluxAvgCount = 1._Real / NBtrSubcycles;
   // Forward-backward feedback weights for the predictor-corrector subcycle.
   constexpr Real Gamma1     = 0.5333_Real;
   constexpr Real Gamma2     = 0.5333_Real;
   constexpr Real Gamma3     = 1.0_Real;
   constexpr Real RhoGravity = RhoSw * Gravity;

   Pacer::start("SE-RK2:stage2BtrPC", 2);

   // Get array required
   Array1DReal NormalBtrVelCur = State->getNormalBarotropicVelocity(CurLevel);
   Array1DReal BtrPressAnomalyCur =
       State->getBarotropicPressureAnomaly(CurLevel);
   Array1DReal NormalBtrVelNext = State->getNormalBarotropicVelocity(NextLevel);
   Array1DReal BtrPressAnomalyNext =
       State->getBarotropicPressureAnomaly(NextLevel);

   // Cur: the state at the start of the subcycle
   // Pre: the predictor output
   // Cor: the corrector output
   Array1DReal NormalBtrVelSubcycleCur =
       Scratch->NormalBarotropicVelocitySubcycleCur;
   Array1DReal NormalBtrVelSubcyclePre =
       Scratch->NormalBarotropicVelocitySubcyclePre;
   Array1DReal NormalBtrVelSubcycleCor =
       Scratch->NormalBarotropicVelocitySubcycleCor;
   Array1DReal BtrPressAnomalySubcycleCur =
       Scratch->BarotropicPressureAnomalySubcycleCur;
   Array1DReal BtrPressAnomalySubcyclePre =
       Scratch->BarotropicPressureAnomalySubcyclePre;
   Array1DReal BtrPressAnomalySubcycleCor =
       Scratch->BarotropicPressureAnomalySubcycleCor;

   Array1DReal BtrForcing = Scratch->BarotropicForcing;
   Array1DReal BtrFlux    = Scratch->BarotropicFlux;

   Array1DReal BaroclinicPseudoThicknessEdge =
       Scratch->BaroclinicPseudoThicknessEdge;

   const Array2DReal FluxPseudoThickEdge =
       AuxState->PseudoThicknessAux.FluxPseudoThickEdge;
   const bool UpwindFluxThickness =
       AuxState->PseudoThicknessAux.FluxThickEdgeChoice ==
       FluxThickEdgeOption::Upwind;

   OMEGA_SCOPE(CellsOnEdge, Mesh->CellsOnEdge);
   OMEGA_SCOPE(NEdgesOnCell, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(EdgesOnCell, Mesh->EdgesOnCell);
   OMEGA_SCOPE(EdgeSignOnCell, Mesh->EdgeSignOnCell);
   OMEGA_SCOPE(NEdgesOnEdge, Mesh->NEdgesOnEdge);
   OMEGA_SCOPE(EdgesOnEdge, Mesh->EdgesOnEdge);
   OMEGA_SCOPE(WeightsOnEdge, Mesh->WeightsOnEdge);
   OMEGA_SCOPE(DcEdge, Mesh->DcEdge);
   OMEGA_SCOPE(DvEdge, Mesh->DvEdge);
   OMEGA_SCOPE(AreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(FEdge, Mesh->FEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(DepthMeanSpecVol, EqState->DepthMeanSpecificVolume);

   // Initialize barotropic vars
   deepCopy(NormalBtrVelSubcycleCur, NormalBtrVelCur);
   deepCopy(BtrPressAnomalySubcycleCur, BtrPressAnomalyCur);
   deepCopy(NormalBtrVelNext, NormalBtrVelSubcycleCur);
   deepCopy(BtrFlux, 0._Real);

   // Compute the baroclinic pseudo thickness on edges
   parallelForOuter(
       "computeBtrBaroclinicPseudoThicknessEdge", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
          const I4 KMin = MinLayerEdgeBot(IEdge);
          const I4 KMax = MaxLayerEdgeTop(IEdge);

          Real PseudoThicknessSum = 0._Real;
          parallelReduceInner(
              Team, Range{KMin, KMax},
              INNER_LAMBDA(I4 K, Real & Sum) {
                 Sum += FluxPseudoThickEdge(IEdge, K);
              },
              PseudoThicknessSum);

          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 BaroclinicPseudoThicknessEdge(IEdge) = PseudoThicknessSum;
              });
       });

   //-----------------------------------------------------------------------
   // Communication-avoiding barotropic subcycling:
   //   - Each kernel can fill one halo layer less than its inputs cover,
   //     so the ranges below shrink in turn, and the last of them sets the
   //     floor on HaloWidth. Halo exchange is conducted once at the top
   //     of each subcycle.
   const I4 HaloWidth = static_cast<I4>(Mesh->NCellsHaloH.extent(0));

   // Every kernel downstream loses the layers its stencil reaches past
   // what its inputs cover.
   constexpr I4 VelPredShrink   = 0;
   constexpr I4 PressPredShrink = 1;
   constexpr I4 VelCorrShrink   = 1;
   constexpr I4 PressCorrShrink = 2;

   constexpr I4 MinHaloWidth =
       1 + Kokkos::max(Kokkos::max(VelPredShrink, PressPredShrink),
                       Kokkos::max(VelCorrShrink, PressCorrShrink));

   if (HaloWidth < MinHaloWidth)
      ABORT_ERROR("Split-explicit barotropic subcycling needs a Decomp "
                  "HaloWidth of at least {}, got {}",
                  MinHaloWidth, HaloWidth);

   // A halo deeper than the minimum just leaves every range deeper.
   const I4 HaloAll            = HaloWidth - 1;
   const I4 VelPredEdgeRange   = Mesh->NEdgesHaloH(HaloAll - VelPredShrink);
   const I4 PressPredCellRange = Mesh->NCellsHaloH(HaloAll - PressPredShrink);
   const I4 VelCorrEdgeRange   = Mesh->NEdgesHaloH(HaloAll - VelCorrShrink);
   const I4 PressCorrCellRange = Mesh->NCellsHaloH(HaloAll - PressCorrShrink);

   // Start the barotrolic subcycling
   for (I4 Subcycle = 0; Subcycle < NBtrSubcycles; ++Subcycle) {

      // Restore buffers. This is the only communication inside the subcycle.
      MeshHalo->exchangeFullArrayHalo(NormalBtrVelSubcycleCur, OnEdge);
      MeshHalo->exchangeFullArrayHalo(BtrPressAnomalySubcycleCur, OnCell);

      // Barotropic velocity predictor
      parallelFor(
          "btrVelocityPredictor", {VelPredEdgeRange}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcyclePre(IEdge) = 0._Real;
                return;
             }

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real Mask = EdgeMask(IEdge, KMin);

             Real CoriolisTend = 0._Real;
             for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
                const I4 JEdge = EdgesOnEdge(IEdge, J);
                CoriolisTend += WeightsOnEdge(IEdge, J) *
                                NormalBtrVelSubcycleCur(JEdge) * FEdge(JEdge);
             }

             const Real MeanSpecVolEdge =
                 0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));

             const Real BtrPressAnomalyGrad =
                 (BtrPressAnomalySubcycleCur(Cell1) -
                  BtrPressAnomalySubcycleCur(Cell0)) /
                 DcEdge(IEdge);

             NormalBtrVelSubcyclePre(IEdge) =
                 Mask * (NormalBtrVelSubcycleCur(IEdge) +
                         BtrDt * (CoriolisTend -
                                  MeanSpecVolEdge * BtrPressAnomalyGrad +
                                  BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly predictor
      parallelFor(
          "btrPressurePredictor", {PressPredCellRange},
          KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);
                const Real NormalBtrVelEdge =
                    feedbackBlend(Gamma1, NormalBtrVelSubcycleCur(JEdge),
                                  NormalBtrVelSubcyclePre(JEdge));

                const Real DeltaBtrPressure0 =
                    BtrPressAnomalySubcycleCur(Cell0) -
                    BtrPressAnomalyNext(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressAnomalySubcycleCur(Cell1) -
                    BtrPressAnomalyNext(Cell1);
                const Real DeltaBtrPressureEdge =
                    deltaBtrPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                         NormalBtrVelEdge, UpwindFluxThickness);

                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;
                const Real PredictorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;

                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * PredictorBtrFlux;
             }

             constexpr Real SurfaceFreshwaterFlux = 0._Real;
             BtrPressAnomalySubcyclePre(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * (BtrFluxDivTend / AreaCell(ICell) +
                          RhoGravity * SurfaceFreshwaterFlux);
          });

      // Barotropic velocity corrector
      parallelFor(
          "btrVelocityCorrector", {VelCorrEdgeRange}, KOKKOS_LAMBDA(I4 IEdge) {
             const I4 KMin = MinLayerEdgeBot(IEdge);
             if (MaxLayerEdgeTop(IEdge) < KMin) {
                NormalBtrVelSubcycleCor(IEdge) = 0._Real;
                return;
             }

             const Real Mask = EdgeMask(IEdge, KMin);
             const I4 Cell0  = CellsOnEdge(IEdge, 0);
             const I4 Cell1  = CellsOnEdge(IEdge, 1);

             Real CoriolisTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnEdge(IEdge); ++J) {
                const I4 JEdge = EdgesOnEdge(IEdge, J);
                CoriolisTend += WeightsOnEdge(IEdge, J) *
                                NormalBtrVelSubcyclePre(JEdge) * FEdge(JEdge);
             }

             const Real BtrPressure0 =
                 feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                               BtrPressAnomalySubcyclePre(Cell0));

             const Real BtrPressure1 =
                 feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                               BtrPressAnomalySubcyclePre(Cell1));

             const Real MeanSpecVolEdge =
                 0.5_Real * (DepthMeanSpecVol(Cell0) + DepthMeanSpecVol(Cell1));

             const Real BtrPressureGrad =
                 (BtrPressure1 - BtrPressure0) / DcEdge(IEdge);

             NormalBtrVelSubcycleCor(IEdge) =
                 Mask *
                 (NormalBtrVelSubcycleCur(IEdge) +
                  BtrDt * (CoriolisTend - MeanSpecVolEdge * BtrPressureGrad +
                           BtrForcing(IEdge)));
          });

      // Barotropic pressure anomaly corrector
      parallelFor(
          "btrPressureCorrector", {PressCorrCellRange},
          KOKKOS_LAMBDA(I4 ICell) {
             Real BtrFluxDivTend = 0._Real;

             for (I4 J = 0; J < NEdgesOnCell(ICell); ++J) {
                const I4 JEdge = EdgesOnCell(ICell, J);
                if (MaxLayerEdgeTop(JEdge) < MinLayerEdgeBot(JEdge))
                   continue;

                const I4 Cell0 = CellsOnEdge(JEdge, 0);
                const I4 Cell1 = CellsOnEdge(JEdge, 1);

                const Real BtrPressure0 =
                    feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                                  BtrPressAnomalySubcyclePre(Cell0));

                const Real BtrPressure1 =
                    feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                                  BtrPressAnomalySubcyclePre(Cell1));

                const Real NormalBtrVelEdge =
                    feedbackBlend(Gamma3, NormalBtrVelSubcycleCur(JEdge),
                                  NormalBtrVelSubcycleCor(JEdge));

                const Real DeltaBtrPressure0 =
                    BtrPressure0 - BtrPressAnomalyNext(Cell0);
                const Real DeltaBtrPressure1 =
                    BtrPressure1 - BtrPressAnomalyNext(Cell1);
                const Real DeltaBtrPressureEdge =
                    deltaBtrPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                         NormalBtrVelEdge, UpwindFluxThickness);

                const Real BtrPressureEdge =
                    RhoGravity * BaroclinicPseudoThicknessEdge(JEdge) +
                    DeltaBtrPressureEdge;

                const Real CorrectorBtrFlux =
                    BtrPressureEdge * NormalBtrVelEdge;
                BtrFluxDivTend +=
                    EdgeSignOnCell(ICell, J) * DvEdge(JEdge) * CorrectorBtrFlux;
             }

             BtrPressAnomalySubcycleCor(ICell) =
                 BtrPressAnomalySubcycleCur(ICell) +
                 BtrDt * BtrFluxDivTend / AreaCell(ICell);
          });

      // Accumulate the corrector barotropic velocity and flux for the
      // time-averaged transport used by the barotropic-baroclinic velocity
      // correction.
      parallelFor(
          "btrCorrectorAccumulate", {Mesh->NEdgesOwned},
          KOKKOS_LAMBDA(I4 IEdge) {
             NormalBtrVelNext(IEdge) += NormalBtrVelSubcycleCor(IEdge);

             if (MaxLayerEdgeTop(IEdge) < MinLayerEdgeBot(IEdge))
                return;

             const I4 Cell0 = CellsOnEdge(IEdge, 0);
             const I4 Cell1 = CellsOnEdge(IEdge, 1);

             const Real BtrPressure0 =
                 feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell0),
                               BtrPressAnomalySubcycleCor(Cell0));

             const Real BtrPressure1 =
                 feedbackBlend(Gamma2, BtrPressAnomalySubcycleCur(Cell1),
                               BtrPressAnomalySubcycleCor(Cell1));

             const Real NormalBtrVelEdge =
                 feedbackBlend(Gamma3, NormalBtrVelSubcycleCur(IEdge),
                               NormalBtrVelSubcycleCor(IEdge));

             const Real DeltaBtrPressure0 =
                 BtrPressure0 - BtrPressAnomalyNext(Cell0);
             const Real DeltaBtrPressure1 =
                 BtrPressure1 - BtrPressAnomalyNext(Cell1);
             const Real DeltaBtrPressureEdge =
                 deltaBtrPressureEdge(DeltaBtrPressure0, DeltaBtrPressure1,
                                      NormalBtrVelEdge, UpwindFluxThickness);

             const Real BtrPressureEdge =
                 RhoGravity * BaroclinicPseudoThicknessEdge(IEdge) +
                 DeltaBtrPressureEdge;

             BtrFlux(IEdge) += BtrPressureEdge * NormalBtrVelEdge;
          });

      // Rotate the subcycle buffers: the corrector output becomes the state
      // at the start of the next subcycle.
      std::swap(NormalBtrVelSubcycleCur, NormalBtrVelSubcycleCor);
      std::swap(BtrPressAnomalySubcycleCur, BtrPressAnomalySubcycleCor);
   } // Subcycle

   parallelFor(
       "btrSubcycleAverage", {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(I4 IEdge) {
          NormalBtrVelNext(IEdge) *= InvBtrVelAvgCount;
          BtrFlux(IEdge) *= InvBtrFluxAvgCount;
       });

   deepCopy(BtrPressAnomalyNext, BtrPressAnomalySubcycleCur);
   MeshHalo->exchangeFullArrayHalo(NormalBtrVelNext, OnEdge);
   MeshHalo->exchangeFullArrayHalo(BtrFlux, OnEdge);
   MeshHalo->exchangeFullArrayHalo(BtrPressAnomalyNext, OnCell);

   Pacer::stop("SE-RK2:stage2BtrPC", 2);
}

} // namespace OMEGA
