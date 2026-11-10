#include "AuxiliaryState.h"
#include "Config.h"
#include "Error.h"
#include "Field.h"
#include "Logging.h"
#include "PGrad.h"
#include "Pacer.h"
#include "SubmesoEddies.h"
#include "Tendencies.h"
#include "TimeStepper.h"

namespace OMEGA {

// create the static class members
AuxiliaryState *AuxiliaryState::DefaultAuxState = nullptr;
std::map<std::string, std::unique_ptr<AuxiliaryState>>
    AuxiliaryState::AllAuxStates;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

// Constructor. Constructs the member auxiliary variables and registers their
// fields with IOStreams
AuxiliaryState::AuxiliaryState(const std::string &Name, const HorzMesh *Mesh,
                               Halo *MeshHalo, VertCoord *VCoord, VertAdv *VAdv,
                               int NTracers, TimeInterval TimeStep)
    : Mesh(Mesh), MeshHalo(MeshHalo), VCoord(VCoord), VAdv(VAdv),
      Name(stripDefault(Name)), KineticAux(stripDefault(Name), Mesh, VCoord),
      PseudoThicknessAux(stripDefault(Name), Mesh, VCoord),
      VorticityAux(stripDefault(Name), Mesh, VCoord),
      VelocityDel2Aux(stripDefault(Name), Mesh, VCoord),
      SurfTracerRestAux(stripDefault(Name), Mesh, NTracers),
      TracerAux(stripDefault(Name), Mesh, VCoord, NTracers),
      TransportAux(stripDefault(Name), Mesh, VCoord),
      MixedLayerAux(stripDefault(Name), Mesh, VCoord), TimeStep(TimeStep) {

   GroupName = "AuxiliaryState";
   if (Name != "Default") {
      GroupName.append(Name);
   }
   std::string AuxMeshName = Mesh->MeshName;

   auto AuxGroup = FieldGroup::create(GroupName);

   KineticAux.registerFields(GroupName, AuxMeshName);
   PseudoThicknessAux.registerFields(GroupName, AuxMeshName);
   VorticityAux.registerFields(GroupName, AuxMeshName);
   VelocityDel2Aux.registerFields(GroupName, AuxMeshName);
   SurfTracerRestAux.registerFields(GroupName, AuxMeshName);
   TracerAux.registerFields(GroupName, AuxMeshName);
   TransportAux.registerFields(GroupName, AuxMeshName);
   MixedLayerAux.registerFields(GroupName, AuxMeshName);
}

// Destructor. Unregisters the fields with IOStreams and destroys this auxiliary
// state field group.
AuxiliaryState::~AuxiliaryState() {
   KineticAux.unregisterFields();
   PseudoThicknessAux.unregisterFields();
   VorticityAux.unregisterFields();
   VelocityDel2Aux.unregisterFields();
   SurfTracerRestAux.unregisterFields();
   TracerAux.unregisterFields();
   TransportAux.unregisterFields();
   MixedLayerAux.unregisterFields();

   FieldGroup::destroy(GroupName);
}

// Compute auxiliary variables for vertical dynamics
void AuxiliaryState::computeMomVertAux(const OceanState *State,
                                       const Array3DReal &TracerArray,
                                       int ThickTimeLevel) const {

   Pacer::start("AuxState:computeMomVertAux", 2);

   Eos *EosInstance = Eos::getInstance();

   // get pseudo-thickness
   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);

   // get temperature and salinity
   I4 ConservTempIdx;
   I4 AbsSalinityIdx;
   Tracers::getIndex(ConservTempIdx, "Temperature");
   Tracers::getIndex(AbsSalinityIdx, "Salinity");

   const auto ConservTemp =
       Kokkos::subview(TracerArray, ConservTempIdx, Kokkos::ALL, Kokkos::ALL);
   const auto AbsSalinity =
       Kokkos::subview(TracerArray, AbsSalinityIdx, Kokkos::ALL, Kokkos::ALL);

   // compute pressure
   const auto &SurfacePressure = VCoord->SurfacePressure;
   VCoord->computePressure(PseudoThickCell, SurfacePressure);

   // compute specific volume. The FiniteVolume pressure gradient also needs
   // the specific volume derivatives, which computeSpecVolAndDerivs supplies
   // from the same single equation-of-state evaluation. It fills SpecVol as
   // well, so it replaces rather than accompanies computeSpecVol; the branch
   // is here so that a Centered run does not pay for the derivative
   // arithmetic. The default PressureGrad instance is queried because the
   // scheme is a model-wide choice and AuxiliaryState holds no PressureGrad.
   const auto &PressureMid   = VCoord->PressureMid;
   const PressureGrad *PGrad = PressureGrad::getDefault();
   if (PGrad && PGrad->getType() == PressureGradType::FiniteVolume) {
      EosInstance->computeSpecVolAndDerivs(ConservTemp, AbsSalinity,
                                           PressureMid);
   } else {
      EosInstance->computeSpecVol(ConservTemp, AbsSalinity, PressureMid);
   }

   // compute geometric height
   VCoord->computeGeomZHeight(PseudoThickCell, EosInstance->SpecVol);

   EosInstance->computeDepthMeanSpecificVolume(PseudoThickCell);

   // compute Brunt-Vaisala freqency squared
   EosInstance->computeBruntVaisalaFreqSq(ConservTemp, AbsSalinity, PressureMid,
                                          EosInstance->SpecVol);

   // compute target thickness
   VCoord->computeTargetThickness();

   // compute displaced spec volume for mixed layer depth
   EosInstance->computeSpecVolDisp(ConservTemp, AbsSalinity,
                                   MixedLayerAux.ReferencePressure, 0);

   // compute mixed layer depth
   const auto &SpecVolDisplaced = EosInstance->SpecVolDisplaced;
   OMEGA_SCOPE(LocMixedlLayerAux, this->MixedLayerAux);
   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocMixedlLayerAux.computeVarsOnCell(Team, ICell, SpecVolDisplaced);
       });

   Pacer::stop("AuxState:computeMomVertAux", 2);
}

// Compute transport velocity for pseudo-thickness and tracers
void AuxiliaryState::computeTransportVelocity(
    const OceanState *State, int VelTimeLevel,
    const Array2DReal &TransportVelocityAdd) const {
   Pacer::start("AuxState:computeTransportVelocity", 2);

   Array2DReal NormalVelocity          = State->getNormalVelocity(VelTimeLevel);
   const auto &NormalTransportVelocity = TransportAux.NormalTransportVelocity;

   OMEGA_SCOPE(NVertLayers, VCoord->NVertLayers);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   auto *SubEddies = SubmesoEddies::getInstance();
   Array2DReal EddyVelocity;
   if (SubEddies && SubEddies->Enable) {

      Eos *EosInstance = Eos::getInstance();

      const auto &MeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;
      const auto &SpecVol             = EosInstance->SpecVol;
      const auto &BVFreqSq            = EosInstance->BruntVaisalaFreqSq;
      const auto &GeomZMid            = VCoord->GeomZMid;
      const auto &MinLayerEdgeBot     = VCoord->MinLayerEdgeBot;
      const auto &MaxLayerEdgeTop     = VCoord->MaxLayerEdgeTop;
      const auto &DenMixLayerDepth    = MixedLayerAux.DenMixLayerDepth;
      const auto &DenMixLayerIndex    = MixedLayerAux.DenMixLayerIndex;

      SubEddies->computeBuoyGrad(SpecVol, MeanPseudoThickEdge, GeomZMid,
                                 BVFreqSq);
      SubEddies->computeEddyVelocity(DenMixLayerDepth, DenMixLayerIndex,
                                     BVFreqSq, MeanPseudoThickEdge);

      EddyVelocity = SubEddies->EddyVelocity;
   }

   parallelForOuter(
       "computeTransportVelocity", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          // copy NormalVelocity using full vertical extent to preserve
          // the 3-zone masking enforced by applyEdgeLayerMask
          parallelForInner(
              Team, NVertLayers, INNER_LAMBDA(int K) {
                 NormalTransportVelocity(IEdge, K) = NormalVelocity(IEdge, K);
              });

          teamBarrier(Team);

          const int KMin = MinLayerEdgeBot(IEdge);
          const int KMax = MaxLayerEdgeTop(IEdge);

          parallelForInner(
              Team, Range{KMin, KMax}, INNER_LAMBDA(int K) {
                 if (TransportVelocityAdd.data()) {
                    NormalTransportVelocity(IEdge, K) +=
                        TransportVelocityAdd(IEdge, K);
                 }
                 if (EddyVelocity.data()) {
                    NormalTransportVelocity(IEdge, K) += EddyVelocity(IEdge, K);
                 }
              });
       });

   Pacer::stop("AuxState:computeTransportVelocity", 2);
}

// Compute the auxiliary variables needed for pseudo-thickness equation
void AuxiliaryState::computePseudoThicknessAux(
    const OceanState *State, const Array3DReal &TracerArray, int ThickTimeLevel,
    int VelTimeLevel, const TimeInterval ProjDt) const {

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(LocPseudoThicknessAux, PseudoThicknessAux);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   Pacer::start("AuxState:computePseudoThickAux", 2);
   parallelForOuter(
       "computePseudoThickAux", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnEdge(Team, IEdge, PseudoThickCell,
                                                  NormalVelEdge);
       });
   Pacer::stop("AuxState:computePseudoThickAux", 2);

   auto *SubEddies = SubmesoEddies::getInstance();

   computeMomVertAux(State, TracerArray, ThickTimeLevel);

   computeTransportVelocity(State, VelTimeLevel);

   Pacer::start("AuxState:computeVerticalPseudoVelocity", 2);

   const auto &FluxPseudoThickEdge     = PseudoThicknessAux.FluxPseudoThickEdge;
   const auto &NormalTransportVelocity = TransportAux.NormalTransportVelocity;
   VAdv->computeVerticalTransportPseudoVelocity(NormalTransportVelocity,
                                                FluxPseudoThickEdge,
                                                PseudoThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeVerticalPseudoVelocity", 2);
}

// Compute the auxiliary variables needed for momentum equation
void AuxiliaryState::computeMomAux(const OceanState *State,
                                   const Array3DReal &TracerArray,
                                   int ThickTimeLevel, int VelTimeLevel,
                                   const TimeInterval ProjDt) const {

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);

   OMEGA_SCOPE(LocKineticAux, KineticAux);
   OMEGA_SCOPE(LocPseudoThicknessAux, PseudoThicknessAux);
   OMEGA_SCOPE(LocVorticityAux, VorticityAux);
   OMEGA_SCOPE(LocVelocityDel2Aux, VelocityDel2Aux);

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerVertexTop, VCoord->MinLayerVertexTop);
   OMEGA_SCOPE(MaxLayerVertexBot, VCoord->MaxLayerVertexBot);
   OMEGA_SCOPE(MinLayerEdgeTop, VCoord->MinLayerEdgeTop);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   // Sanity checks for the time steps
   OMEGA_REQUIRE(TimeStepSeconds > 0, "TimeStepSeconds has to be positive");
   OMEGA_REQUIRE(ProjDtSeconds > 0, "ProjDtSeconds has to be positive");

   Pacer::start("AuxState:computeMomAux", 1);

   computeMomVertAux(State, TracerArray, ThickTimeLevel);

   // Restrict to vertices fully surrounded by local cells in NCellAll
   // to avoid division by zero in PsuedoThickVertex on outer halo vertices.
   const auto HaloWidth = Mesh->NVerticesHaloH.extent_int(0);
   Pacer::start("AuxState:vertexAuxState1", 2);
   parallelForOuter(
       "vertexAuxState1",
       LaunchConfig({Mesh->NVerticesHaloH(HaloWidth - 1)},
                    TeamScratch<Real>(2 * VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          LocVorticityAux.computeVarsOnVertex(Team, IVertex, PseudoThickCell,
                                              NormalVelEdge);
       });
   Pacer::stop("AuxState:vertexAuxState1", 2);

   Pacer::start("AuxState:cellAuxState1", 2);
   parallelForOuter(
       "cellAuxState1",
       LaunchConfig({Mesh->NCellsAll},
                    TeamScratch<Real>(2 * VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocKineticAux.computeVarsOnCell(Team, ICell, NormalVelEdge);
       });
   Pacer::stop("AuxState:cellAuxState1", 2);

   const auto &VelocityDivCell = KineticAux.VelocityDivCell;
   const auto &RelVortVertex   = VorticityAux.RelVortVertex;

   // Zero Del2Edge at boundary layers before accumulating over neighbouring
   // edges in computeVarsOnCell/Vertex. Del4 hyperdiffusion reads Del2Edge at
   // all edges sharing a cell or vertex; fill values in boundary layers would
   // corrupt Del2DivCell and Del2RelVortVertex.
   VCoord->zeroEdgeField(VelocityDel2Aux.Del2Edge, Mesh->NEdgesAll);

   Pacer::start("AuxState:edgeAuxState2", 2);
   parallelForOuter(
       "edgeAuxState2", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnEdge(Team, IEdge, PseudoThickCell,
                                                  NormalVelEdge);

          LocVelocityDel2Aux.computeVarsOnEdge(Team, IEdge, VelocityDivCell,
                                               RelVortVertex);
       });

   parallelForOuter(
       "edgeAuxState2", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocVorticityAux.computeVarsOnEdge(Team, IEdge);
       });
   Pacer::stop("AuxState:edgeAuxState2", 2);

   Pacer::start("AuxState:vertexAuxState2", 2);
   parallelForOuter(
       "vertexAuxState2",
       LaunchConfig({Mesh->NVerticesAll},
                    TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          LocVelocityDel2Aux.computeVarsOnVertex(Team, IVertex);
       });
   Pacer::stop("AuxState:vertexAuxState2", 2);

   Pacer::start("AuxState:cellAuxState2", 2);
   parallelForOuter(
       "cellAuxState2",
       LaunchConfig({Mesh->NCellsAll}, TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocVelocityDel2Aux.computeVarsOnCell(Team, ICell);
       });
   Pacer::stop("AuxState:cellAuxState2", 2);

   Pacer::start("AuxState:computeVerticalPseudoVelocity", 2);

   const auto &FluxPseudoThickEdge = PseudoThicknessAux.FluxPseudoThickEdge;
   VAdv->computeVerticalPseudoVelocity(NormalVelEdge, FluxPseudoThickEdge,
                                       PseudoThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeVerticalPseudoVelocity", 2);

   Pacer::stop("AuxState:computeMomAux", 1);
}

// Compute the auxiliary variables needed for tracer equation
void AuxiliaryState::computeTracerAux(const OceanState *State,
                                      const Array3DReal &TracerArray,
                                      int ThickTimeLevel, int VelTimeLevel,
                                      const TimeInterval ProjDt) const {

   OMEGA_SCOPE(LocPseudoThicknessAux, PseudoThicknessAux);
   OMEGA_SCOPE(LocTracerAux, TracerAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);

   Pacer::start("AuxState:computePseudoThickAux", 2);
   parallelForOuter(
       "computePseudoThickAux", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnEdge(Team, IEdge, PseudoThickCell,
                                                  NormalVelEdge);
       });
   Pacer::stop("AuxState:computePseudoThickAux", 2);

   const int NTracers              = Tracers::getNumTracers();
   const auto &MeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;

   Pacer::start("Tend:computeTracerAuxCell", 2);
   parallelForOuter(
       "computeTracerAuxCell",
       LaunchConfig({NTracers, Mesh->NCellsAll},
                    TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          LocTracerAux.computeVarsOnCells(Team, LTracer, ICell,
                                          MeanPseudoThickEdge, TracerArray);
       });
   Pacer::stop("Tend:computeTracerAuxCell", 2);

   auto *SubEddies = SubmesoEddies::getInstance();

   computeMomVertAux(State, TracerArray, ThickTimeLevel);

   computeTransportVelocity(State, VelTimeLevel);

   const auto &NormalTransportVelocity = TransportAux.NormalTransportVelocity;

   Pacer::start("AuxState:cellThicknessAux", 2);
   parallelForOuter(
       "cellThicknessAux",
       LaunchConfig({Mesh->NCellsAll}, TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnCells(Team, ICell, PseudoThickCell,
                                                   NormalTransportVelocity,
                                                   TimeStepSeconds);
       });
   Pacer::stop("AuxState:cellThicknessAux", 2);

   Pacer::start("AuxState:computeVerticalPseudoVelocity", 2);

   const auto &FluxPseudoThickEdge = PseudoThicknessAux.FluxPseudoThickEdge;
   VAdv->computeVerticalTransportPseudoVelocity(NormalTransportVelocity,
                                                FluxPseudoThickEdge,
                                                PseudoThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeVerticalPseudoVelocity", 2);
}

// Compute the auxiliary variables
void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray,
                                int ThickTimeLevel, int VelTimeLevel,
                                const TimeInterval ProjDt) const {
   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);

   const int NTracers = TracerArray.extent_int(0);

   OMEGA_SCOPE(LocPseudoThicknessAux, PseudoThicknessAux);
   OMEGA_SCOPE(LocTracerAux, TracerAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   Pacer::start("AuxState:computeAll", 1);

   computeMomAux(State, TracerArray, ThickTimeLevel, VelTimeLevel, ProjDt);

   const auto &MeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;

   Pacer::start("AuxState:cellTracerAux", 2);
   parallelForOuter(
       "tracerCellAux",
       LaunchConfig({NTracers, Mesh->NCellsAll},
                    TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          LocTracerAux.computeVarsOnCells(Team, LTracer, ICell,
                                          MeanPseudoThickEdge, TracerArray);
       });
   Pacer::stop("AuxState:cellTracerAux", 2);

   computeTransportVelocity(State, VelTimeLevel);

   const auto &FluxPseudoThickEdge     = PseudoThicknessAux.FluxPseudoThickEdge;
   const auto &NormalTransportVelocity = TransportAux.NormalTransportVelocity;

   Pacer::start("AuxState:cellThickAux", 2);
   parallelForOuter(
       "thickCellAux",
       LaunchConfig({Mesh->NCellsAll}, TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnCells(Team, ICell, PseudoThickCell,
                                                   NormalTransportVelocity,
                                                   TimeStepSeconds);
       });
   Pacer::stop("AuxState:cellThickAux", 2);

   VAdv->computeVerticalTransportPseudoVelocity(NormalTransportVelocity,
                                                FluxPseudoThickEdge,
                                                PseudoThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeAll", 1);
}

void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray, int TimeLevel,
                                const TimeInterval ProjDt) const {
   computeAll(State, TracerArray, TimeLevel, TimeLevel, ProjDt);
}

void AuxiliaryState::computePseudoThicknessTracerAux(
    const OceanState *State, const Array3DReal &TracerArray, int ThickTimeLevel,
    int VelTimeLevel, const TimeInterval ProjDt,
    const Array2DReal &TransportVelocityAdd) const {

   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(VelTimeLevel);

   const int NTracers = TracerArray.extent_int(0);

   OMEGA_SCOPE(LocPseudoThicknessAux, PseudoThicknessAux);
   OMEGA_SCOPE(LocTracerAux, TracerAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   Pacer::start("AuxState:computePseudoThicknessTracerAux", 1);

   computeMomVertAux(State, TracerArray, ThickTimeLevel);

   Pacer::start("AuxState:edgeThicknessTracerAux", 2);
   parallelForOuter(
       "edgeThicknessTracerAux", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnEdge(Team, IEdge, PseudoThickCell,
                                                  NormalVelEdge);
       });
   Pacer::stop("AuxState:edgeThicknessTracerAux", 2);

   const auto &MeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;

   Pacer::start("AuxState:cellTracerAux", 2);
   parallelForOuter(
       "cellTracerAux",
       LaunchConfig({NTracers, Mesh->NCellsAll},
                    TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          LocTracerAux.computeVarsOnCells(Team, LTracer, ICell,
                                          MeanPseudoThickEdge, TracerArray);
       });
   Pacer::stop("AuxState:cellTracerAux", 2);

   computeTransportVelocity(State, VelTimeLevel, TransportVelocityAdd);

   const auto &NormalTransportVelocity = TransportAux.NormalTransportVelocity;

   Pacer::start("AuxState:cellThicknessAux", 2);
   parallelForOuter(
       "cellThicknessAux",
       LaunchConfig({Mesh->NCellsAll}, TeamScratch<Real>(VCoord->NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          LocPseudoThicknessAux.computeVarsOnCells(Team, ICell, PseudoThickCell,
                                                   NormalTransportVelocity,
                                                   TimeStepSeconds);
       });
   Pacer::stop("AuxState:cellThicknessAux", 2);

   Pacer::start("AuxState:computeVerticalPseudoVelocity", 2);

   const auto &FluxPseudoThickEdge = PseudoThicknessAux.FluxPseudoThickEdge;
   VAdv->computeVerticalTransportPseudoVelocity(NormalTransportVelocity,
                                                FluxPseudoThickEdge,
                                                PseudoThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeVerticalPseudoVelocity", 2);

   Pacer::stop("AuxState:computePseudoThicknessTracerAux", 1);
}

// Create a non-default auxiliary state
AuxiliaryState *AuxiliaryState::create(const std::string &Name,
                                       const HorzMesh *Mesh, Halo *MeshHalo,
                                       VertCoord *VCoord, VertAdv *VAdv,
                                       const int NTracers,
                                       TimeInterval TimeStep) {
   OMEGA_REQUIRE(
       Mesh, "Null HorzMesh pointer in AuxiliaryState::create with Name = {}",
       Name);
   OMEGA_REQUIRE(MeshHalo,
                 "Null Halo pointer in AuxiliaryState::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(
       VCoord,
       "Null VertCoord pointer in AuxiliaryState::create with Name = {}", Name);
   OMEGA_REQUIRE(
       VAdv, "Null VertAdv pointer in AuxiliaryState::create with Name = {}",
       Name);

   if (AllAuxStates.find(Name) != AllAuxStates.end()) {
      LOG_ERROR("Attempted to create a new AuxiliaryState with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewAuxState = new AuxiliaryState(Name, Mesh, MeshHalo, VCoord, VAdv,
                                          NTracers, TimeStep);
   AllAuxStates.emplace(Name, NewAuxState);

   return NewAuxState;
}

// Create the default auxiliary state. Assumes that HorzMesh, VertCoord,
// VertAdv, and Halo have been initialized.
void AuxiliaryState::init() {
   const HorzMesh *DefMesh = HorzMesh::getDefault();
   OMEGA_REQUIRE(DefMesh,
                 "Null default HorzMesh pointer in AuxiliaryState::init");
   Halo *DefHalo = Halo::getDefault();
   OMEGA_REQUIRE(DefHalo, "Null default Halo pointer in AuxiliaryState::init");
   VertCoord *DefVCoord = VertCoord::getDefault();
   OMEGA_REQUIRE(DefVCoord,
                 "Null default VertCoord pointer in AuxiliaryState::init");
   VertAdv *DefVAdv = VertAdv::getDefault();
   OMEGA_REQUIRE(DefVAdv,
                 "Null default VertAdv pointer in AuxiliaryState::init");
   const TimeStepper *DefTimeStepper = TimeStepper::getDefault();
   OMEGA_REQUIRE(DefTimeStepper,
                 "Null default TimeStepper pointer in AuxiliaryState::init");

   int NTracers          = Tracers::getNumTracers();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();

   AuxiliaryState::DefaultAuxState = AuxiliaryState::create(
       "Default", DefMesh, DefHalo, DefVCoord, DefVAdv, NTracers, TimeStep);

   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig,
                 "Null OmegaConfig pointer in AuxiliaryState::init");
   DefaultAuxState->readConfigOptions(OmegaConfig);
}

// Get the default auxiliary state
AuxiliaryState *AuxiliaryState::getDefault() {
   return AuxiliaryState::DefaultAuxState;
}

// Get auxiliary state by name
AuxiliaryState *AuxiliaryState::get(const std::string &Name) {
   // look for an instance of this name
   auto it = AllAuxStates.find(Name);

   // if found, return the pointer
   if (it != AllAuxStates.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("AuxiliaryState::get: Attempt to retrieve non-existent "
                "auxiliary state:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
}

// Remove auxiliary state by name
void AuxiliaryState::erase(const std::string &Name) {
   AllAuxStates.erase(Name);
}

// Remove all auxiliary states
void AuxiliaryState::clear() {
   AllAuxStates.clear();
   DefaultAuxState = nullptr; // prevent dangling pointer
}

// Read and set config options
void AuxiliaryState::readConfigOptions(Config *OmegaConfig) {

   Error Err; // error code

   Config AdvectConfig("Advection");
   Err += OmegaConfig->get(AdvectConfig);
   CHECK_ERROR_ABORT(Err, "AuxiliaryState: Advection group not in Config");

   std::string FluxThickTypeStr;
   Err += AdvectConfig.get("FluxThicknessType", FluxThickTypeStr);
   CHECK_ERROR_ABORT(
       Err, "AuxiliaryState: FluxThicknessType not found in AdvectConfig");

   if (FluxThickTypeStr == "Center") {
      this->PseudoThicknessAux.FluxThickEdgeChoice =
          FluxThickEdgeOption::Center;
   } else if (FluxThickTypeStr == "Upwind") {
      this->PseudoThicknessAux.FluxThickEdgeChoice =
          FluxThickEdgeOption::Upwind;
   } else {
      ABORT_ERROR("AuxiliaryState: Unknown FluxThicknessType requested");
   }
}

//------------------------------------------------------------------------------
// Perform auxiliary state halo exchange
// Note that only non-computed auxiliary variables needs to be exchanged
I4 AuxiliaryState::exchangeHalo() {
   I4 Err = 0;

   // Performing halo exchange on individual tracers because full halo exchange
   // on a 2D array assumes the first dimension is the vertical
   const I4 NTracers =
       SurfTracerRestAux.TracersMonthlySurfClimoCell.extent_int(0);
   for (I4 LTracer = 0; LTracer < NTracers; ++LTracer) {
      auto TracerSurfClimoCell = Kokkos::subview(
          SurfTracerRestAux.TracersMonthlySurfClimoCell, LTracer, Kokkos::ALL);
      Err += MeshHalo->exchangeFullArrayHalo(TracerSurfClimoCell, OnCell);
   }

   return Err;

} // end exchangeHalo

} // namespace OMEGA
