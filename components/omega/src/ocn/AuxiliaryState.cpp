#include "AuxiliaryState.h"
#include "Logging.h"
#include "MetaData.h"

namespace OMEGA {

// create the static class members
AuxiliaryState *AuxiliaryState::DefaultAuxState = nullptr;
std::map<std::string, std::unique_ptr<AuxiliaryState>>
    AuxiliaryState::AllAuxStates;

// Constructor. Constructs the member auxiliary variables and registers their
// fields with IOStreams
AuxiliaryState::AuxiliaryState(const std::string &Name, const HorzMesh *Mesh,
                               int NVertLevels)
    : Mesh(Mesh), Name(Name), KineticAux(Name, Mesh, NVertLevels),
      LayerThicknessAux(Name, Mesh, NVertLevels),
      VorticityAux(Name, Mesh, NVertLevels),
      VelocityDel2Aux(Name, Mesh, NVertLevels) {

   GroupName = "AuxiliaryState";
   if (Name != "Default") {
      GroupName.append(Name);
   }

   MetaGroup::create(GroupName);

   KineticAux.registerFields(GroupName);
   LayerThicknessAux.registerFields(GroupName);
   VorticityAux.registerFields(GroupName);
   VelocityDel2Aux.registerFields(GroupName);
}

// Destructor. Unregisters the fields with IOStreams and destroys this auxiliary
// state meta group.
AuxiliaryState::~AuxiliaryState() {
   KineticAux.unregisterFields();
   LayerThicknessAux.unregisterFields();
   VorticityAux.unregisterFields();
   VelocityDel2Aux.unregisterFields();

   MetaGroup::destroy(GroupName);
}

// Compute the auxiliary variables
void AuxiliaryState::computeAll(const OceanState *State, int TimeLevel) const {
   const Array2DReal &LayerThickCell = State->LayerThickness[TimeLevel];
   const Array2DReal &NormalVelEdge  = State->NormalVelocity[TimeLevel];

   const int NVertLevels = LayerThickCell.extent_int(1);
   const int NChunks     = NVertLevels / VecLength;

   OMEGA_SCOPE(LocKineticAux, KineticAux);
   OMEGA_SCOPE(LocLayerThicknessAux, LayerThicknessAux);
   OMEGA_SCOPE(LocVorticityAux, VorticityAux);
   OMEGA_SCOPE(LocVelocityDel2Aux, VelocityDel2Aux);

   parallelFor(
       "vertexAuxState1", {Mesh->NVerticesHaloH(2), NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVorticityAux.computeVarsOnVertex(IVertex, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });

   parallelFor(
       "cellAuxState1", {Mesh->NCellsHaloH(2), NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocKineticAux.computeVarsOnCell(ICell, KChunk, NormalVelEdge);
       });

   const auto &VelocityDivCell = KineticAux.VelocityDivCell;
   const auto &RelVortVertex   = VorticityAux.RelVortVertex;

   parallelFor(
       "edgeAuxState1", {Mesh->NEdgesHaloH(1), NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LocVorticityAux.computeVarsOnEdge(IEdge, KChunk);
          LocLayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                                 NormalVelEdge);
          LocVelocityDel2Aux.computeVarsOnEdge(IEdge, KChunk, VelocityDivCell,
                                               RelVortVertex);
       });

   parallelFor(
       "vertexAuxState2", {Mesh->NVerticesHaloH(0), NChunks},
       KOKKOS_LAMBDA(int IVertex, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnVertex(IVertex, KChunk);
       });

   parallelFor(
       "cellAuxState2", {Mesh->NCellsHaloH(0), NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocVelocityDel2Aux.computeVarsOnCell(ICell, KChunk);
       });

   parallelFor(
       "cellAuxState3", {Mesh->NCellsHaloH(1), NChunks},
       KOKKOS_LAMBDA(int ICell, int KChunk) {
          LocLayerThicknessAux.computeVarsOnCells(ICell, KChunk,
                                                  LayerThickCell);
       });
}

// Create a non-default auxiliary state
AuxiliaryState *AuxiliaryState::create(const std::string &Name,
                                       const HorzMesh *Mesh, int NVertLevels) {
   if (AllAuxStates.find(Name) != AllAuxStates.end()) {
      LOG_ERROR("Attempted to create a new AuxiliaryState with name {} but it "
                "already exists",
                Name);
      return nullptr;
   }

   auto *NewAuxState = new AuxiliaryState(Name, Mesh, NVertLevels);
   AllAuxStates.emplace(Name, NewAuxState);

   return NewAuxState;
}

// Create the default auxiliary state. Assumes that HorzMesh has been
// initialized.
int AuxiliaryState::init() {
   int Err                 = 0;
   const HorzMesh *DefMesh = HorzMesh::getDefault();

   // These hard-wired variable needs to be updated
   // with retrivals/config options
   const int NVertLevels = 60;

   AuxiliaryState::DefaultAuxState =
       AuxiliaryState::create("Default", DefMesh, NVertLevels);
   return Err;
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

} // namespace OMEGA
