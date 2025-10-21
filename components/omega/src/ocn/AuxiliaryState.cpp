#include "AuxiliaryState.h"
#include "Config.h"
#include "Field.h"
#include "Logging.h"
#include "Pacer.h"

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
                               Halo *MeshHalo, const VertCoord *VCoord,
                               int NTracers)
    : Mesh(Mesh), MeshHalo(MeshHalo), VCoord(VCoord), Name(stripDefault(Name)),
      KineticAux(stripDefault(Name), Mesh, VCoord),
      LayerThicknessAux(stripDefault(Name), Mesh, VCoord),
      VorticityAux(stripDefault(Name), Mesh, VCoord),
      VelocityDel2Aux(stripDefault(Name), Mesh, VCoord),
      WindForcingAux(stripDefault(Name), Mesh),
      TracerAux(stripDefault(Name), Mesh, VCoord, NTracers) {

   GroupName = "AuxiliaryState";
   if (Name != "Default") {
      GroupName.append(Name);
   }
   std::string AuxMeshName = Mesh->MeshName;

   auto AuxGroup = FieldGroup::create(GroupName);

   KineticAux.registerFields(GroupName, AuxMeshName);
   LayerThicknessAux.registerFields(GroupName, AuxMeshName);
   VorticityAux.registerFields(GroupName, AuxMeshName);
   VelocityDel2Aux.registerFields(GroupName, AuxMeshName);
   WindForcingAux.registerFields(GroupName, AuxMeshName);
   TracerAux.registerFields(GroupName, AuxMeshName);
}

// Destructor. Unregisters the fields with IOStreams and destroys this auxiliary
// state field group.
AuxiliaryState::~AuxiliaryState() {
   KineticAux.unregisterFields();
   LayerThicknessAux.unregisterFields();
   VorticityAux.unregisterFields();
   VelocityDel2Aux.unregisterFields();
   WindForcingAux.unregisterFields();
   TracerAux.unregisterFields();

   FieldGroup::destroy(GroupName);
}

// Compute the auxiliary variables needed for momentum equation
void AuxiliaryState::computeMomAux(const OceanState *State, int ThickTimeLevel,
                                   int VelTimeLevel) const {
   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   State->getLayerThickness(LayerThickCell, ThickTimeLevel);
   State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

   const int NVertLayers = LayerThickCell.extent_int(1);
   const int NChunks     = NVertLayers / VecLength;

   OMEGA_SCOPE(LocKineticAux, KineticAux);
   OMEGA_SCOPE(LocLayerThicknessAux, LayerThicknessAux);
   OMEGA_SCOPE(LocVorticityAux, VorticityAux);
   OMEGA_SCOPE(LocVelocityDel2Aux, VelocityDel2Aux);
   OMEGA_SCOPE(LocWindForcingAux, WindForcingAux);

   Pacer::start("AuxState:computeMomAux", 1);

   Pacer::start("AuxState:vertexAuxState1", 2);
   parallelFor(
       "vertexAuxState1", {Mesh->NVerticesAll, NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVorticityAux.computeVarsOnVertex(IVertex, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });
   Pacer::stop("AuxState:vertexAuxState1", 2);

   Pacer::start("AuxState:cellAuxState1", 2);
   parallelFor(
       "cellAuxState1", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocKineticAux.computeVarsOnCell(ICell, KChunk, NormalVelEdge);
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
   parallelFor(
       "edgeAuxState2", {Mesh->NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LocVorticityAux.computeVarsOnEdge(IEdge, KChunk);
          LocLayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                                 NormalVelEdge);
          LocVelocityDel2Aux.computeVarsOnEdge(IEdge, KChunk, VelocityDivCell,
                                               RelVortVertex);
       });
   Pacer::stop("AuxState:edgeAuxState2", 2);

   Pacer::start("AuxState:vertexAuxState2", 2);
   parallelFor(
       "vertexAuxState2", {Mesh->NVerticesAll, NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnVertex(IVertex, KChunk);
       });
   Pacer::stop("AuxState:vertexAuxState2", 2);

   Pacer::start("AuxState:cellAuxState2", 2);
   parallelFor(
       "cellAuxState2", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnCell(ICell, KChunk);
       });
   Pacer::stop("AuxState:cellAuxState2", 2);

   Pacer::start("AuxState:cellAuxState3", 2);
   parallelFor(
       "cellAuxState3", {Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocLayerThicknessAux.computeVarsOnCells(ICell, KChunk,
                                                  LayerThickCell);
       });
   Pacer::stop("AuxState:cellAuxState3", 2);

   Pacer::stop("AuxState:computeMomAux", 1);
}

// Compute the auxiliary variables
void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray,
                                int ThickTimeLevel, int VelTimeLevel) const {
   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   State->getLayerThickness(LayerThickCell, ThickTimeLevel);
   State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

   const int NVertLayers = LayerThickCell.extent_int(1);
   const int NChunks     = NVertLayers / VecLength;
   const int NTracers    = TracerArray.extent_int(0);

   OMEGA_SCOPE(LocTracerAux, TracerAux);

   Pacer::start("AuxState:computeAll", 1);

   computeMomAux(State, ThickTimeLevel, VelTimeLevel);

   Pacer::start("AuxState:edgeAuxState4", 2);
   parallelFor(
       "edgeAuxState4", {NTracers, Mesh->NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int IEdge, int KChunk) {
          LocTracerAux.computeVarsOnEdge(LTracer, IEdge, KChunk, NormalVelEdge,
                                         LayerThickCell, TracerArray);
       });
   Pacer::stop("AuxState:edgeAuxState4", 2);

   const auto &MeanLayerThickEdge = LayerThicknessAux.MeanLayerThickEdge;

   Pacer::start("AuxState:cellAuxState4", 2);
   parallelFor(
       "cellAuxState4", {NTracers, Mesh->NCellsAll, NChunks},
       KOKKOS_LAMBDA(int LTracer, int ICell, int KChunk) {
          LocTracerAux.computeVarsOnCells(LTracer, ICell, KChunk,
                                          MeanLayerThickEdge, TracerArray);
       });
   Pacer::stop("AuxState:cellAuxState4", 2);

   Pacer::stop("AuxState:computeAll", 1);
}

void AuxiliaryState::computeAll(const OceanState *State,
                                const Array3DReal &TracerArray,
                                int TimeLevel) const {
   computeAll(State, TracerArray, TimeLevel, TimeLevel);
}

// Create a non-default auxiliary state
AuxiliaryState *AuxiliaryState::create(const std::string &Name,
                                       const HorzMesh *Mesh, Halo *MeshHalo,
                                       const VertCoord *VCoord,
                                       const int NTracers) {
   if (AllAuxStates.find(Name) != AllAuxStates.end()) {
      LOG_ERROR("Attempted to create a new AuxiliaryState with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewAuxState =
       new AuxiliaryState(Name, Mesh, MeshHalo, VCoord, NTracers);
   AllAuxStates.emplace(Name, NewAuxState);

   return NewAuxState;
}

// Create the default auxiliary state. Assumes that HorzMesh, VertCoord and
// Halo have been initialized.
void AuxiliaryState::init() {
   const HorzMesh *DefMesh    = HorzMesh::getDefault();
   Halo *DefHalo              = Halo::getDefault();
   const VertCoord *DefVCoord = VertCoord::getDefault();

   int NVertLayers = VertCoord::getDefault()->NVertLayers;
   int NTracers    = Tracers::getNumTracers();

   AuxiliaryState::DefaultAuxState =
       AuxiliaryState::create("Default", DefMesh, DefHalo, DefVCoord, NTracers);

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
void AuxiliaryState::clear() { AllAuxStates.clear(); }

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
      this->LayerThicknessAux.FluxThickEdgeChoice = FluxThickEdgeOption::Center;
   } else if (FluxThickTypeStr == "Upwind") {
      this->LayerThicknessAux.FluxThickEdgeChoice = FluxThickEdgeOption::Upwind;
   } else {
      ABORT_ERROR("AuxiliaryState: Unknown FluxThicknessType requested");
   }

   std::string FluxTracerTypeStr;
   Err += AdvectConfig.get("FluxTracerType", FluxTracerTypeStr);
   CHECK_ERROR_ABORT(
       Err, "AuxiliaryState: FluxTracerType not found in AdvectConfig");

   if (FluxTracerTypeStr == "Center") {
      this->TracerAux.TracersOnEdgeChoice = FluxTracerEdgeOption::Center;
   } else if (FluxTracerTypeStr == "Upwind") {
      this->TracerAux.TracersOnEdgeChoice = FluxTracerEdgeOption::Upwind;
   } else {
      ABORT_ERROR("AuxiliaryState: Unknown FluxTracerType requested");
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

   return Err;

} // end exchangeHalo

} // namespace OMEGA
