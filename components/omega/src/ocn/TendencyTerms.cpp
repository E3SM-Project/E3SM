//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//
#include <iomanip>
#include <iostream>

#include "DataTypes.h"
#include "Eos.h"
#include "Error.h"
#include "HorzMesh.h"
#include "HorzOperators.h"
#include "TendencyTerms.h"
#include "Tracers.h"

namespace OMEGA {

PseudoThicknessFluxDivOnCell::PseudoThicknessFluxDivOnCell(
    const HorzMesh *Mesh, const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), NVertLayers(VCoord->NVertLayers),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

CoriolisAccelerationOnEdge::CoriolisAccelerationOnEdge(const HorzMesh *Mesh,
                                                       const VertCoord *VCoord)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge), NVertLayers(VCoord->NVertLayers),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 const VertCoord *VCoord)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SfcStressForcingOnEdge::SfcStressForcingOnEdge(const HorzMesh *Mesh,
                                               const VertCoord *VCoord)
    : Enabled(false), EdgeMask(VCoord->EdgeMask),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

BottomDragOnEdge::BottomDragOnEdge(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : Enabled(false), Coeff(0), CellsOnEdge(Mesh->CellsOnEdge),
      NVertLayers(VCoord->NVertLayers), EdgeMask(VCoord->EdgeMask),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SfcThicknessForcingOnCell::SfcThicknessForcingOnCell(const HorzMesh *Mesh,
                                                     const VertCoord *VCoord)
    : MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

SfcTracerForcingOnCell::SfcTracerForcingOnCell(const HorzMesh *Mesh,
                                               const VertCoord *VCoord,
                                               I4 TempTracerIndex,
                                               I4 SaltTracerIndex,
                                               const Eos *EosInst)
    : TempIndex(TempTracerIndex), SaltIndex(SaltTracerIndex),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell),
      EosChoice(EosInst->EosChoice) {}

PenetratingShortwaveOnCell::PenetratingShortwaveOnCell(const HorzMesh *Mesh,
                                                       const VertCoord *VCoord,
                                                       I4 TempTracerIndex)
    : TempIndex(TempTracerIndex), NVertLayers(VCoord->NVertLayers),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {}

TracerHorzAdvOnCell::TracerHorzAdvOnCell(const HorzMesh *Mesh,
                                         const VertCoord *VCoord,
                                         const VertAdv *VAdv)
    : HorzontalMesh(Mesh), VerticalCoord(VCoord),
      NVertLayers(VCoord->NVertLayers),
      NAdvCellsForEdge("NumberOfCellsContribToAdvectionAtEdge",
                       Mesh->NEdgesAll),
      AdvCellsForEdge("IndexOfCellsContributingToAdvection", Mesh->NEdgesAll,
                      Mesh->MaxEdges2 + 2),
      AdvMaskHighOrder("MaskForHighOrderAdvectionTerms", Mesh->NEdgesAll,
                       VCoord->NVertLayers),
      CellsOnCell(Mesh->CellsOnCell),
      AdvCoefs("CommonAdvectionCoefficients", Mesh->MaxEdges2 + 2,
               Mesh->NEdgesAll),
      AdvCoefs3rd("CommonAdvectionCoeffsForHighOrder", Mesh->MaxEdges2 + 2,
                  Mesh->NEdgesAll),
      HighOrderFlxHorz("HigherOrderHorizontalFlux", Tracers::getNumTracers(),
                       Mesh->NEdgesAll, VCoord->NVertLayers),
      TracerCur(), NEdgesOnCell(Mesh->NEdgesOnCell),
      EdgesOnCell(Mesh->EdgesOnCell), CellsOnEdge(Mesh->CellsOnEdge),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop),
      EdgeSignOnCell(Mesh->EdgeSignOnCell), DvEdge(Mesh->DvEdge),
      AreaCell(Mesh->AreaCell), TotalVerticalTransportPseudoVelocity(
                                    VAdv->TotalVerticalTransportPseudoVelocity),
      HProvInv(), HNewInv(), HProv(), TracerMax(), TracerMin(), HighOrderFlx(),
      LowOrderFlx(), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell), WorkTend(), FlxIn(), FlxOut(),
      ActiveTracerHorizontalAdvectionEdgeFlux(),
      ActiveTracerHorizontalAdvectionTendency() {
   deepCopy(HighOrderFlxHorz, 0);
}

TracerDiffOnCell::TracerDiffOnCell(const HorzMesh *Mesh,
                                   const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel2(Mesh->MeshScalingDel2), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

KPPNonLocalTracerFluxOnCell::KPPNonLocalTracerFluxOnCell(
    const HorzMesh *Mesh, const VertCoord *VCoord)
    : Enabled(false), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell) {}

TracerHyperDiffOnCell::TracerHyperDiffOnCell(const HorzMesh *Mesh,
                                             const VertCoord *VCoord)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      CellsOnEdge(Mesh->CellsOnEdge), EdgeSignOnCell(Mesh->EdgeSignOnCell),
      DvEdge(Mesh->DvEdge), DcEdge(Mesh->DcEdge), AreaCell(Mesh->AreaCell),
      MeshScalingDel4(Mesh->MeshScalingDel4), NVertLayers(VCoord->NVertLayers),
      EdgeMask(VCoord->EdgeMask), MinLayerCell(VCoord->MinLayerCell),
      MaxLayerCell(VCoord->MaxLayerCell),
      MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

SurfaceTracerRestoringOnCell::SurfaceTracerRestoringOnCell(
    const HorzMesh *Mesh) {}

FrazilOnCell::FrazilOnCell(const HorzMesh *Mesh, const VertCoord *VCoord)
    : NCellsAll(Mesh->NCellsAll), TempTracerIndex(-1), SaltTracerIndex(-1),
      MinLayerCell(VCoord->MinLayerCell), MaxLayerCell(VCoord->MaxLayerCell) {
   Tracers::getIndex(TempTracerIndex, "Temperature");
   Tracers::getIndex(SaltTracerIndex, "Salinity");

   OMEGA_REQUIRE(TempTracerIndex >= 0,
                 "FrazilOnCell: Temperature tracer index is undefined");
   OMEGA_REQUIRE(SaltTracerIndex >= 0,
                 "FrazilOnCell: Salinity tracer index is undefined");
}

void FrazilOnCell::operator()(const Array2DReal &PseudoThicknessTend,
                              const Array3DReal &TracerTend,
                              const Array3DReal &TracerArray,
                              const Array2DReal &PressureMid,
                              const Array2DReal &PseudoThickness,
                              const TimeInterval &TimeStep,
                              const Real FinalUpdateWeight) const {
   auto *Frazil = Frazil::getDefault();
   if (!Enabled || !Frazil) {
      return;
   }

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   deepCopy(Frazil->FrazilTTend, 0.0_Real);
   deepCopy(Frazil->FrazilSTend, 0.0_Real);
   deepCopy(Frazil->FrazilHTend, 0.0_Real);
   deepCopy(Frazil->AccMIce, 0.0_Real);
   deepCopy(Frazil->AccEIce, 0.0_Real);
   deepCopy(Frazil->AccMLiq, 0.0_Real);
   deepCopy(Frazil->AccELiq, 0.0_Real);
   deepCopy(Frazil->AccMSalt, 0.0_Real);

   const auto ConservTemp =
       Kokkos::subview(TracerArray, TempTracerIndex, Kokkos::ALL, Kokkos::ALL);
   const auto AbsSalinity =
       Kokkos::subview(TracerArray, SaltTracerIndex, Kokkos::ALL, Kokkos::ALL);

   Frazil->computeFrazil(ConservTemp, AbsSalinity, PressureMid,
                         PseudoThickness);
   Frazil->accumulateOcnStepFluxes(FinalUpdateWeight, TimeStepSeconds);

   const auto FrazilHTend = Frazil->FrazilHTend;
   const auto FrazilTTend = Frazil->FrazilTTend;
   const auto FrazilSTend = Frazil->FrazilSTend;
   const I4 TempIndex     = TempTracerIndex;
   const I4 SaltIndex     = SaltTracerIndex;

   OMEGA_SCOPE(LocPseudoThicknessTend, PseudoThicknessTend);
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocFrazilHTend, FrazilHTend);
   OMEGA_SCOPE(LocFrazilTTend, FrazilTTend);
   OMEGA_SCOPE(LocFrazilSTend, FrazilSTend);
   OMEGA_SCOPE(LocMinLayerCell, MinLayerCell);
   OMEGA_SCOPE(LocMaxLayerCell, MaxLayerCell);

   parallelForOuter(
       {NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin = LocMinLayerCell(ICell);
          const int KMax = LocMaxLayerCell(ICell);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 LocPseudoThicknessTend(ICell, K) +=
                     LocFrazilHTend(ICell, K) / TimeStepSeconds;
                 LocTracerTend(TempIndex, ICell, K) +=
                     LocFrazilTTend(ICell, K) / TimeStepSeconds;
                 LocTracerTend(SaltIndex, ICell, K) +=
                     LocFrazilSTend(ICell, K) / TimeStepSeconds;
              });
       });
}

void TracerHorzAdvOnCell::init() {
   const HorzMesh *Mesh    = this->HorzontalMesh;
   const VertCoord *VCoord = this->VerticalCoord;
   const auto MaxEdges2    = Mesh->MaxEdges2;
   const auto NEdgesAll    = Mesh->NEdgesAll;
   const auto NCellsAll    = Mesh->NCellsAll;
   const auto NCellsSize   = Mesh->NCellsSize;
   const auto NEdgesSize   = Mesh->NEdgesSize;
   // Allocate Kokkos arrays in member data

   if (ForceLowOrder) {
      // Return when the 2nd-order tracer horz adv
      deepCopy(NAdvCellsForEdge, 0);
      deepCopy(AdvMaskHighOrder, 0);
      return;
   }

   SecondDerivativeOnCell secondDerivativeOnCell(Mesh);
   Array3DReal DerivTwo("DerivTwo", MaxEdges2 + 2, 2, NEdgesAll);
   parallelFor(
       {NCellsAll},
       KOKKOS_LAMBDA(int ICell) { secondDerivativeOnCell(DerivTwo, ICell); });
   // Compute masks and coefficients
   MasksAndCoefficients masksAndCoefficients(
       Mesh, VCoord, DerivTwo, NAdvCellsForEdge, AdvCellsForEdge,
       AdvMaskHighOrder, AdvCoefs, AdvCoefs3rd);
   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) { masksAndCoefficients(IEdge); });
   if (FCT) {
      const int NVertsFCT = NVertLayers + 1;
      HProvInv =
          Array2DReal("FCTProvesionalLayerThickness", NCellsSize, NVertsFCT);
      HNewInv = Array2DReal("FCTProvesionalNewInverse", NCellsSize, NVertsFCT);
      HProv   = Array2DReal("FCTProvesionalThickness", NCellsSize, NVertsFCT);
      TracerCur    = Array2DReal("TracerCur", NCellsSize, NVertsFCT),
      TracerMax    = Array2DReal("FCTTracerMax", NCellsSize, NVertsFCT);
      TracerMin    = Array2DReal("FCTTracerMin", NCellsSize, NVertsFCT);
      HighOrderFlx = Array2DReal("FCTHighOrderFlx", NEdgesSize, NVertsFCT);
      LowOrderFlx  = Array2DReal("FCTLowOrderFlx", NEdgesSize, NVertsFCT);
      WorkTend     = Array2DReal("WorkTend", NCellsSize, NVertsFCT);
      FlxIn        = Array2DReal("FlxIn", NCellsSize, NVertsFCT);
      FlxOut       = Array2DReal("FlxOut", NCellsSize, NVertsFCT);
      deepCopy(HProvInv, 0.0);
      deepCopy(HNewInv, 0.0);
      deepCopy(HProv, 0.0);
      deepCopy(TracerCur, 0.0);
      deepCopy(TracerMax, 0.0);
      deepCopy(TracerMin, 0.0);
      deepCopy(HighOrderFlx, 0.0);
      deepCopy(LowOrderFlx, 0.0);
      deepCopy(WorkTend, 0.0);
      deepCopy(FlxIn, 0.0);
      deepCopy(FlxOut, 0.0);
      if (ComputeBudgets) {
         // Sized like the NEdges and NCells dimensions they are written with
         const int NTracers = Tracers::getNumTracers();
         ActiveTracerHorizontalAdvectionEdgeFlux =
             Array3DReal("FCTActiveTracerHorizontalAdvectionEdgeFlux", NTracers,
                         NEdgesSize, NVertLayers);
         ActiveTracerHorizontalAdvectionTendency =
             Array3DReal("FCTActiveTracerHorizontalAdvectionTendency", NTracers,
                         NCellsSize, NVertLayers);
         deepCopy(ActiveTracerHorizontalAdvectionEdgeFlux, 0.0);
         deepCopy(ActiveTracerHorizontalAdvectionTendency, 0.0);
         const std::string GroupName = "AuxiliaryState";
         // Both budget arrays hold every tracer on every layer
         std::vector<std::string> FluxDimNames{"NTracers", "NEdges",
                                               "NVertLayers"};
         std::vector<std::string> TendDimNames{"NTracers", "NCells",
                                               "NVertLayers"};

         auto BudgetAdvectionEdgeFlux = Field::create(
             ActiveTracerHorizontalAdvectionEdgeFlux.label(),    // field name
             "Tracer FCT Horizontal Advection Edge Flux Budget", // long name
             "",                                                 // units
             "",                                  // CF standard Name
             std::numeric_limits<Real>::lowest(), // min valid value
             std::numeric_limits<Real>::max(),    // max valid value
             3,                                   // number of dimensions
             FluxDimNames                         // dimension names
         );
         BudgetAdvectionEdgeFlux->attachData<Array3DReal>(
             ActiveTracerHorizontalAdvectionEdgeFlux);
         FieldGroup::addFieldToGroup(
             ActiveTracerHorizontalAdvectionEdgeFlux.label(), GroupName);

         auto BudgetAdvectionCellTend = Field::create(
             ActiveTracerHorizontalAdvectionTendency.label(),      // field name
             "Tracer FCT Horizontal Advection Cell Flux Tendency", // long name
             "",                                                   // units
             "",                                  // CF standard Name
             std::numeric_limits<Real>::lowest(), // min valid value
             std::numeric_limits<Real>::max(),    // max valid value
             3,                                   // number of dimensions
             TendDimNames                         // dimension names
         );
         BudgetAdvectionCellTend->attachData<Array3DReal>(
             ActiveTracerHorizontalAdvectionTendency);
         FieldGroup::addFieldToGroup(
             ActiveTracerHorizontalAdvectionTendency.label(), GroupName);
      }
   }
}
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
