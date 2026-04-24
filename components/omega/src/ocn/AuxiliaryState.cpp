#include "AuxiliaryState.h"
#include "Config.h"
#include "Field.h"
#include "Logging.h"
#include "Pacer.h"
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
      WindForcingAux(stripDefault(Name), Mesh),
      SurfTracerRestAux(stripDefault(Name), Mesh, NTracers),
      TracerAux(stripDefault(Name), Mesh, VCoord, NTracers) {

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
   WindForcingAux.registerFields(GroupName, AuxMeshName);
   SurfTracerRestAux.registerFields(GroupName, AuxMeshName);
   TracerAux.registerFields(GroupName, AuxMeshName);
}

// Destructor. Unregisters the fields with IOStreams and destroys this auxiliary
// state field group.
AuxiliaryState::~AuxiliaryState() {
   KineticAux.unregisterFields();
   PseudoThicknessAux.unregisterFields();
   VorticityAux.unregisterFields();
   VelocityDel2Aux.unregisterFields();
   WindForcingAux.unregisterFields();
   SurfTracerRestAux.unregisterFields();
   TracerAux.unregisterFields();

   FieldGroup::destroy(GroupName);
}

// Compute auxiliary variables for vertical dynamics
void AuxiliaryState::computeMomVertAux(const OceanState *State,
                                       const Array3DReal &TracerArray,
                                       int ThickTimeLevel,
                                       int VelTimeLevel) const {

   Pacer::start("AuxState:computeMomVertAux", 2);

   Eos *EosInstance = Eos::getInstance();

   // get pseudo-thickness
   Array2DReal PseudoThickCell = State->getPseudoThickness(ThickTimeLevel);
   // get normal velocity
   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);

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

   // compute specific volume
   const auto &PressureMid = VCoord->PressureMid;
   EosInstance->computeSpecVol(ConservTemp, AbsSalinity, PressureMid);

   // compute geometric height
   VCoord->computeZHeight(PseudoThickCell, EosInstance->SpecVol);

   // compute target thickness
   VCoord->computeTargetThickness();

   Pacer::stop("AuxState:computeMomVertAux", 2);
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
   OMEGA_SCOPE(LocWindForcingAux, WindForcingAux);

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerVertexBot, VCoord->MinLayerVertexBot);
   OMEGA_SCOPE(MinLayerVertexTop, VCoord->MinLayerVertexTop);
   OMEGA_SCOPE(MaxLayerVertexBot, VCoord->MaxLayerVertexBot);
   OMEGA_SCOPE(MaxLayerVertexTop, VCoord->MaxLayerVertexTop);
   OMEGA_SCOPE(MinLayerEdgeTop, VCoord->MinLayerEdgeTop);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeBot, VCoord->MaxLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);
   R8 ProjDtSeconds;
   ProjDt.get(ProjDtSeconds, TimeUnits::Seconds);

   Pacer::start("AuxState:computeMomAux", 1);

   computeMomVertAux(State, TracerArray, ThickTimeLevel, VelTimeLevel);

   Pacer::start("AuxState:vertexAuxState1", 2);
   parallelForOuter(
       "vertexAuxState1", {Mesh->NVerticesAll},
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          const int KMin   = MinLayerVertexTop(IVertex);
          const int KMax   = MaxLayerVertexBot(IVertex);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocVorticityAux.computeVarsOnVertex(
                     IVertex, KChunk, PseudoThickCell, NormalVelEdge);
              });
       });
   Pacer::stop("AuxState:vertexAuxState1", 2);

   Pacer::start("AuxState:cellAuxState1", 2);
   parallelForOuter(
       "cellAuxState1", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocKineticAux.computeVarsOnCell(ICell, KChunk, NormalVelEdge);
              });
       });
   Pacer::stop("AuxState:cellAuxState1", 2);

   const auto &VelocityDivCell = KineticAux.VelocityDivCell;
   const auto &RelVortVertex   = VorticityAux.RelVortVertex;

   Pacer::start("AuxState:edgeAuxState1", 2);
   parallelFor(
       "edgeAuxState1", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          LocWindForcingAux.computeVarsOnEdge(IEdge);
       });
   Pacer::stop("AuxState:edgeAuxState1", 2);

   Pacer::start("AuxState:edgeAuxState2", 2);
   parallelForOuter(
       "edgeAuxState2", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocPseudoThicknessAux.computeVarsOnEdge(
                     IEdge, KChunk, PseudoThickCell, NormalVelEdge);
                 LocVelocityDel2Aux.computeVarsOnEdge(
                     IEdge, KChunk, VelocityDivCell, RelVortVertex);
              });
       });

   parallelForOuter(
       "edgeAuxState2", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeTop(IEdge);
          const int KMax   = MaxLayerEdgeBot(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocVorticityAux.computeVarsOnEdge(IEdge, KChunk);
              });
       });
   Pacer::stop("AuxState:edgeAuxState2", 2);

   Pacer::start("AuxState:vertexAuxState2", 2);
   parallelForOuter(
       "vertexAuxState2", {Mesh->NVerticesAll},
       KOKKOS_LAMBDA(int IVertex, const TeamMember &Team) {
          const int KMin   = MinLayerVertexBot(IVertex);
          const int KMax   = MaxLayerVertexTop(IVertex);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocVelocityDel2Aux.computeVarsOnVertex(IVertex, KChunk);
              });
       });
   Pacer::stop("AuxState:vertexAuxState2", 2);

   Pacer::start("AuxState:cellAuxState2", 2);
   parallelForOuter(
       "cellAuxState2", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocVelocityDel2Aux.computeVarsOnCell(ICell, KChunk);
              });
       });
   Pacer::stop("AuxState:cellAuxState2", 2);

   Pacer::start("AuxState:computeVerticalVelocity", 2);

   const auto &FluxPseudoThickEdge = PseudoThicknessAux.FluxLayerThickEdge;
   VAdv->computeVerticalVelocity(NormalVelEdge, FluxPseudoThickEdge,
                                 LayerThickCell, ProjDtSeconds);

   Pacer::stop("AuxState:computeVerticalVelocity", 2);

   Pacer::stop("AuxState:computeMomAux", 1);
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

   Pacer::start("AuxState:computeAll", 1);

   computeMomAux(State, TracerArray, ThickTimeLevel, VelTimeLevel, ProjDt);

   Pacer::start("AuxState:cellAuxState3", 2);
   parallelForOuter(
       "cellAuxState3", {Mesh->NCellsAll},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocPseudoThicknessAux.computeVarsOnCells(
                     ICell, KChunk, PseudoThickCell, NormalVelEdge,
                     TimeStepSeconds);
              });
       });
   Pacer::stop("AuxState:cellAuxState3", 2);

   const auto &MeanPseudoThickEdge = PseudoThicknessAux.MeanPseudoThickEdge;

   Pacer::start("AuxState:cellAuxState4", 2);
   parallelForOuter(
       "cellAuxState4", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LocTracerAux.computeVarsOnCells(
                     LTracer, ICell, KChunk, MeanPseudoThickEdge, TracerArray);
              });
       });
   Pacer::stop("AuxState:cellAuxState4", 2);

   Pacer::stop("AuxState:computeAll", 1);
}

void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray, int TimeLevel,
                                const TimeInterval ProjDt) const {
   computeAll(State, TracerArray, TimeLevel, TimeLevel, ProjDt);
}

// Create a non-default auxiliary state
AuxiliaryState *AuxiliaryState::create(const std::string &Name,
                                       const HorzMesh *Mesh, Halo *MeshHalo,
                                       VertCoord *VCoord, VertAdv *VAdv,
                                       const int NTracers,
                                       TimeInterval TimeStep) {
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
   const HorzMesh *DefMesh           = HorzMesh::getDefault();
   Halo *DefHalo                     = Halo::getDefault();
   VertCoord *DefVCoord              = VertCoord::getDefault();
   VertAdv *DefVAdv                  = VertAdv::getDefault();
   const TimeStepper *DefTimeStepper = TimeStepper::getDefault();

   int NTracers          = Tracers::getNumTracers();
   TimeInterval TimeStep = DefTimeStepper->getTimeStep();

   AuxiliaryState::DefaultAuxState = AuxiliaryState::create(
       "Default", DefMesh, DefHalo, DefVCoord, DefVAdv, NTracers, TimeStep);

   Config *OmegaConfig = Config::getOmegaConfig();
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

   Config WindStressConfig("WindStress");
   Err += OmegaConfig->get(WindStressConfig);

   std::string WindStressInterpTypeStr;
   Err += WindStressConfig.get("InterpType", WindStressInterpTypeStr);
   CHECK_ERROR_ABORT(
       Err, "AuxiliaryState: InterpType not found in WindStressConfig");

   if (WindStressInterpTypeStr == "Isotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (WindStressInterpTypeStr == "Anisotropic") {
      this->WindForcingAux.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      ABORT_ERROR("AuxiliaryState: Unknown InterpType requested");
   }
}

//------------------------------------------------------------------------------
// Perform auxiliary state halo exchange
// Note that only non-computed auxiliary variables needs to be exchanged
I4 AuxiliaryState::exchangeHalo() {
   I4 Err = 0;

   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.ZonalStressCell, OnCell);
   Err +=
       MeshHalo->exchangeFullArrayHalo(WindForcingAux.MeridStressCell, OnCell);

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
