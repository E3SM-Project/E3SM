//===-- ocn/TendencyTerms.cpp - Tendency Terms ------------------*- C++ -*-===//
//
// The tendency terms that update state variables are implemented as functors,
// i.e. as classes that act like functions. This source defines the class
// constructors for these functors, which initialize the functor objects using
// the Mesh objects and info from the Config. The function call operators () are
// defined in the corresponding header file.
//
//===----------------------------------------------------------------------===//

#include "TendencyTerms.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "HorzMesh.h"
#include "OceanState.h"

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh as alread been initialized.
int Tendencies::init() {
   int Err = 0;

   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   Config *TendConfig;

   int NVertLevels = 60;

   Tendencies::DefaultTendencies =
       Tendencies::create("Default", DefHorzMesh, NVertLevels, TendConfig);

   return Err;

} // end init

//------------------------------------------------------------------------------
// Destroys the tendencies
Tendencies::~Tendencies() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Removes all tendencies instances before exit
void Tendencies::clear() { AllTendencies.clear(); } // end clear

//------------------------------------------------------------------------------
// Removes tendencies from list by name
void Tendencies::erase(const std::string &Name) {

   AllTendencies.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get default tendencies
Tendencies *Tendencies::getDefault() {

   return Tendencies::DefaultTendencies;

} // end get default

//------------------------------------------------------------------------------
// Get tendencies by name
Tendencies *Tendencies::get(const std::string &Name ///< [in] Name of tendencies
) {

   auto it = AllTendencies.find(Name);

   if (it != AllTendencies.end()) {
      return it->second.get();
   } else {
      LOG_ERROR(
          "Tendencies::get: Attempt to retrieve non-existent tendencies:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get tendencies

//------------------------------------------------------------------------------
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       Config *Options, ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : ThicknessFluxDiv(Mesh, Options), PotientialVortHAdv(Mesh, Options),
      KEGrad(Mesh, Options), SSHGrad(Mesh, Options),
      VelocityDiffusion(Mesh, Options), VelocityHyperDiff(Mesh, Options),
      CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend) {

   // Tendency arrays
   LayerThicknessTend =
       Array2DReal("LayerThicknessTend", Mesh->NCellsSize, NVertLevels);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, NVertLevels);

   // Array dimension lengths
   NCellsAll = Mesh->NCellsAll;
   NEdgesAll = Mesh->NEdgesAll;
   NChunks   = NVertLevels / VecLength;

} // end constructor

Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       int NVertLevels, ///< [in] Number of vertical levels
                       Config *Options) ///< [in] Configuration options
    : Tendencies(Name, Mesh, NVertLevels, Options, CustomTendencyType{},
                 CustomTendencyType{}) {}

//------------------------------------------------------------------------------
// Compute tendencies for layer thickness equation
void Tendencies::computeThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int TimeLevel,            ///< [in] Time level
    Real Time                 ///< [in] Time
) {

   OMEGA_SCOPE(LocLayerThicknessTend, LayerThicknessTend);
   OMEGA_SCOPE(LocNCellsAll, NCellsAll);
   OMEGA_SCOPE(LocNChunks, NChunks);
   OMEGA_SCOPE(LocThicknessFluxDiv, ThicknessFluxDiv);
   const Array2DReal &NormalVelEdge = State->NormalVelocity[TimeLevel];

   deepCopy(LocLayerThicknessTend, 0);

   // Compute thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;

   if (LocThicknessFluxDiv.Enabled) {
      parallelFor(
          {LocNCellsAll, LocNChunks}, KOKKOS_LAMBDA(int ICell, int KChunk) {
             LocThicknessFluxDiv(LocLayerThicknessTend, ICell, KChunk,
                                 ThickFluxEdge, NormalVelEdge);
          });
   }

   if (CustomThicknessTend) {
      CustomThicknessTend(LocLayerThicknessTend, State, AuxState, TimeLevel,
                          Time);
   }

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int TimeLevel,            ///< [in] Time level
    Real Time                 ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocNEdgesAll, NEdgesAll);
   OMEGA_SCOPE(LocNChunks, NChunks);
   OMEGA_SCOPE(LocPotientialVortHAdv, PotientialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);

   deepCopy(LocNormalVelocityTend, 0);

   // Compute potential vorticity horizontal advection
   const Array2DReal &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   const Array2DReal &NormVelEdge   = State->NormalVelocity[TimeLevel];
   if (LocPotientialVortHAdv.Enabled) {
      parallelFor(
          {LocNEdgesAll, LocNChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocPotientialVortHAdv(LocNormalVelocityTend, IEdge, KChunk,
                                   NormRVortEdge, NormFEdge, FluxLayerThickEdge,
                                   NormVelEdge);
          });
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      parallelFor(
          {LocNEdgesAll, LocNChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocKEGrad(LocNormalVelocityTend, IEdge, KChunk, KECell);
          });
   }

   // Compute sea surface height gradient
   const Array2DReal &SSHCell = AuxState->LayerThicknessAux.SshCell;
   if (LocSSHGrad.Enabled) {
      parallelFor(
          {LocNEdgesAll, LocNChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocSSHGrad(LocNormalVelocityTend, IEdge, KChunk, SSHCell);
          });
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      parallelFor(
          {LocNEdgesAll, LocNChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityDiffusion(LocNormalVelocityTend, IEdge, KChunk, DivCell,
                                  RVortVertex);
          });
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      parallelFor(
          {LocNEdgesAll, LocNChunks}, KOKKOS_LAMBDA(int IEdge, int KChunk) {
             LocVelocityHyperDiff(LocNormalVelocityTend, IEdge, KChunk,
                                  Del2DivCell, Del2RVortVertex);
          });
   }

   if (CustomVelocityTend) {
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, TimeLevel,
                         Time);
   }

} // end velocity tendency compute

void Tendencies::computeThicknessTendencies(
    OceanState *State,        ///< [in] State variables
    AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int TimeLevel,            ///< [in] Time level
    Real Time                 ///< [in] Time
) {
   // only need LayerThicknessAux on edge
   OMEGA_SCOPE(LayerThicknessAux, AuxState->LayerThicknessAux);
   OMEGA_SCOPE(LayerThickCell, State->LayerThickness[TimeLevel]);
   OMEGA_SCOPE(NormalVelEdge, State->NormalVelocity[TimeLevel]);

   parallelFor(
       "computeLayerThickAux", {NEdgesAll, NChunks},
       KOKKOS_LAMBDA(int IEdge, int KChunk) {
          LayerThicknessAux.computeVarsOnEdge(IEdge, KChunk, LayerThickCell,
                                              NormalVelEdge);
       });

   computeThicknessTendenciesOnly(State, AuxState, TimeLevel, Time);
}

void Tendencies::computeVelocityTendencies(
    OceanState *State,        ///< [in] State variables
    AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int TimeLevel,            ///< [in] Time level
    Real Time                 ///< [in] Time
) {
   AuxState->computeAll(State, TimeLevel);
   computeVelocityTendenciesOnly(State, AuxState, TimeLevel, Time);
}

//------------------------------------------------------------------------------
// Compute both layer thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,       ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int TimeLevel,            ///< [in] Time level
    Real Time                 ///< [in] Time
) {

   AuxState->computeAll(State, TimeLevel);
   computeThicknessTendenciesOnly(State, AuxState, TimeLevel, Time);
   computeVelocityTendenciesOnly(State, AuxState, TimeLevel, Time);

} // end all tendency compute

// TODO: Implement Config options for all constructors
ThicknessFluxDivOnCell::ThicknessFluxDivOnCell(const HorzMesh *Mesh,
                                               Config *Options)
    : NEdgesOnCell(Mesh->NEdgesOnCell), EdgesOnCell(Mesh->EdgesOnCell),
      DvEdge(Mesh->DvEdge), AreaCell(Mesh->AreaCell),
      EdgeSignOnCell(Mesh->EdgeSignOnCell) {

   // Options->get("ThicknessFluxTendencyEnable", Enabled);
}

PotentialVortHAdvOnEdge::PotentialVortHAdvOnEdge(const HorzMesh *Mesh,
                                                 Config *Options)
    : NEdgesOnEdge(Mesh->NEdgesOnEdge), EdgesOnEdge(Mesh->EdgesOnEdge),
      WeightsOnEdge(Mesh->WeightsOnEdge) {

   // Options->get("PVTendencyEnable", Enabled);
}

KEGradOnEdge::KEGradOnEdge(const HorzMesh *Mesh, Config *Options)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge) {

   // Options->get("KETendencyEnable", Enabled);
}

SSHGradOnEdge::SSHGradOnEdge(const HorzMesh *Mesh, Config *Options)
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge) {

   // Options->get("SSHTendencyEnable", Enabled);
}

VelocityDiffusionOnEdge::VelocityDiffusionOnEdge(const HorzMesh *Mesh,
                                                 Config *Options)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel2(Mesh->MeshScalingDel2), EdgeMask(Mesh->EdgeMask) {

   // Options->get("VelDiffTendencyEnable", Enabled);
   // Options->get("ViscDel2", ViscDel2);
}

VelocityHyperDiffOnEdge::VelocityHyperDiffOnEdge(const HorzMesh *Mesh,
                                                 Config *Options)
    : CellsOnEdge(Mesh->CellsOnEdge), VerticesOnEdge(Mesh->VerticesOnEdge),
      DcEdge(Mesh->DcEdge), DvEdge(Mesh->DvEdge),
      MeshScalingDel4(Mesh->MeshScalingDel4), EdgeMask(Mesh->EdgeMask) {

   // Options->get("VelHyperDiffTendencyEnable", Enabled);
   // Options->get("ViscDel4", ViscDel4);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
