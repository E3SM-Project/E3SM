//===-- ocn/Tendencies.cpp - Tendencies ------------------*- C++ -*-===//
//
// The Tendencies class is responsible for managing tendencies of state
// variables and tracers. It contains arrays that store the tendency data
// and provides methods for computing different tendency groups for use
// within the timestepping algorithm. At initialization, it determines which
// tendency terms are enabled.
//
//===----------------------------------------------------------------------===//

#include "Tendencies.h"
#include "CustomTendencyTerms.h"
#include "Error.h"
#include "Pacer.h"
#include "Tracers.h"
#include <string>

namespace OMEGA {

Tendencies *Tendencies::DefaultTendencies = nullptr;
std::map<std::string, std::unique_ptr<Tendencies>> Tendencies::AllTendencies;

//------------------------------------------------------------------------------
// Initialize the tendencies. Assumes that HorzMesh and VertCoord has alread
// been initialized.
void Tendencies::init() {
   Error Err; // error code

   HorzMesh *DefHorzMesh   = HorzMesh::getDefault();
   VertCoord *DefVertCoord = VertCoord::getDefault();

   I4 NTracers = Tracers::getNumTracers();

   // Get TendConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Tendencies: Tendencies group not found in Config");

   // Check if use the customized tendencies. If it is not found in the
   // config, we assume it is not used (false)
   bool UseCustomTendency = false;
   Err += TendConfig.get("UseCustomTendency", UseCustomTendency);

   /// Instances of custom tendencies - empty by default
   CustomTendencyType CustomThickTend;
   CustomTendencyType CustomVelTend;

   if (UseCustomTendency) {
      // Check if use manufactured tendency terms if it is not found in
      // the config file, we will assume it is not used (false)
      bool ManufacturedTend = false;
      Error ManufacturedTendErr =
          TendConfig.get("ManufacturedSolutionTendency", ManufacturedTend);

      if (ManufacturedTend) {
         ManufacturedSolution ManufacturedSol;
         ManufacturedSol.init();

         CustomThickTend = ManufacturedSol.ManufacturedThickTend;
         CustomVelTend   = ManufacturedSol.ManufacturedVelTend;

      } // if ManufacturedTend

   } // end if UseCustomTendency

   // Ceate default tendencies
   Tendencies::DefaultTendencies =
       create("Default", DefHorzMesh, DefVertCoord, NTracers, &TendConfig,
              CustomThickTend, CustomVelTend);

   DefaultTendencies->readTendConfig(&TendConfig);
   if (DefaultTendencies->TracerHighOrderHorzAdv.Enabled)
      DefaultTendencies->TracerHighOrderHorzAdv.init();
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
// read and set config options
void Tendencies::readTendConfig(
    Config *TendConfig ///< [in] Tendencies subconfig
) {
   Error Err; // error code

   Err += TendConfig->get("ThicknessFluxTendencyEnable",
                          this->ThicknessFluxDiv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: ThicknessFluxTendencyEnable not found in TendConfig");

   Err += TendConfig->get("PVTendencyEnable", this->PotientialVortHAdv.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: PVTendencyEnable not found in TendConfig");

   Err += TendConfig->get("KETendencyEnable", this->KEGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: KETendencyEnable not found in TendConfig");

   Err += TendConfig->get("SSHTendencyEnable", this->SSHGrad.Enabled);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: SSHTendencyEnable not found in TendConfig");

   Err += TendConfig->get("VelDiffTendencyEnable",
                          this->VelocityDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelDiffTendencyEnable not found in TendConfig");

   Err += TendConfig->get("VelHyperDiffTendencyEnable",
                          this->VelocityHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: VelHyperDiffTendencyEnable not found in TendConfig");

   if (this->VelocityDiffusion.Enabled) {
      Err += TendConfig->get("ViscDel2", this->VelocityDiffusion.ViscDel2);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel2 not found in TendConfig");
   }

   if (this->VelocityHyperDiff.Enabled) {
      Err += TendConfig->get("ViscDel4", this->VelocityHyperDiff.ViscDel4);
      CHECK_ERROR_ABORT(Err, "Tendencies: ViscDel4 not found in TendConfig");
      Err += TendConfig->get("DivFactor", this->VelocityHyperDiff.DivFactor);
      CHECK_ERROR_ABORT(Err, "Tendencies: DivFactor not found in TendConfig");
   }

   if (TendConfig->existsVar("TracerHorzAdvTendencyEnable")) {
      I4 Order = 1;
      if (TendConfig->existsVar("TracerHorzAdvTendencyOrder"))
         TendConfig->get("TracerHorzAdvTendencyOrder", Order);
      if (Order == 1)
         this->TracerHorzAdv.Enabled = true;
      else if (Order == 2)
         this->TracerHighOrderHorzAdv.Enabled = true;
      else {
         const std::string msg =
             "TracerHorzAdvTendencyOrder: Only values are 1 and 2, found " +
             std::to_string(Order);
         ABORT_ERROR(msg);
      }
   }
   Err += TendConfig->get("TracerDiffTendencyEnable",
                          this->TracerDiffusion.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerDiffTendencyEnable not found in TendConfig");

   Err +=
       TendConfig->get("WindForcingTendencyEnable", this->WindForcing.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: WindForcingTendencyEnable not found in TendConfig");

   Err += TendConfig->get("BottomDragTendencyEnable", this->BottomDrag.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: BottomDragTendencyEnable not found in TendConfig");

   Err += TendConfig->get("BottomDragCoeff", this->BottomDrag.Coeff);
   CHECK_ERROR_ABORT(Err,
                     "Tendencies: BottomDragCoeff not found in TendConfig");

   Err += TendConfig->get("TracerHorzAdvTendencyEnable",
                          this->TracerHorzAdv.Enabled);
   CHECK_ERROR_ABORT(
       Err, "Tendencies: TracerHorzAdvTendencyEnable not found in TendConfig");

   if (this->TracerDiffusion.Enabled) {
      Err += TendConfig->get("EddyDiff2", this->TracerDiffusion.EddyDiff2);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff2 not found in TendConfig");
   }

   Err += TendConfig->get("TracerHyperDiffTendencyEnable",
                          this->TracerHyperDiff.Enabled);
   CHECK_ERROR_ABORT(
       Err,
       "Tendencies: TracerHyperDiffTendencyEnable not found in TendConfig");

   if (this->TracerHyperDiff.Enabled) {
      Err += TendConfig->get("EddyDiff4", this->TracerHyperDiff.EddyDiff4);
      CHECK_ERROR_ABORT(Err, "Tendencies: EddyDiff4 not found in TendConfig");
   }
}

//------------------------------------------------------------------------------
// Construct a new group of tendencies
Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       const VertCoord *VCoord, ///< [in] Vertical coordinate
                       int NTracersIn,          ///< [in] Number of tracers
                       Config *Options,         ///< [in] Configuration options
                       CustomTendencyType InCustomThicknessTend,
                       CustomTendencyType InCustomVelocityTend)
    : Mesh(Mesh), VCoord(VCoord), ThicknessFluxDiv(Mesh, VCoord),
      PotientialVortHAdv(Mesh, VCoord), KEGrad(Mesh, VCoord),
      SSHGrad(Mesh, VCoord), VelocityDiffusion(Mesh, VCoord),
      VelocityHyperDiff(Mesh, VCoord), WindForcing(Mesh, VCoord),
      BottomDrag(Mesh, VCoord), TracerHorzAdv(Mesh, VCoord),
      TracerDiffusion(Mesh, VCoord), TracerHyperDiff(Mesh, VCoord),
      TracerHighOrderHorzAdv(Mesh),
      CustomThicknessTend(InCustomThicknessTend),
      CustomVelocityTend(InCustomVelocityTend) {

   // Tendency arrays
   LayerThicknessTend =
       Array2DReal("LayerThicknessTend", Mesh->NCellsSize, VCoord->NVertLayers);
   NormalVelocityTend =
       Array2DReal("NormalVelocityTend", Mesh->NEdgesSize, VCoord->NVertLayers);
   TracerTend = Array3DReal("TracerTend", NTracersIn, Mesh->NCellsSize,
                            VCoord->NVertLayers);

   NTracers = NTracersIn;

} // end constructor

Tendencies::Tendencies(const std::string &Name, ///< [in] Name for tendencies
                       const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                       const VertCoord *VCoord, ///< [in] Vertical coordinate
                       int NTracersIn,          ///< [in] Number of tracers
                       Config *Options)         ///< [in] Configuration options
    : Tendencies(Name, Mesh, VCoord, NTracersIn, Options, CustomTendencyType{},
                 CustomTendencyType{}) {}

//------------------------------------------------------------------------------
// Compute tendencies for layer thickness equation
void Tendencies::computeThicknessTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocLayerThicknessTend, LayerThicknessTend);
   OMEGA_SCOPE(LocThicknessFluxDiv, ThicknessFluxDiv);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   Array2DReal NormalVelEdge = State->getNormalVelocity(VelTimeLevel);

   Pacer::start("Tend:computeThicknessTendenciesOnly", 1);

   parallelForOuter(
       {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const int K                     = KMin + KChunk;
                 LocLayerThicknessTend(ICell, K) = 0;
              });
       });

   // Compute thickness flux divergence
   const Array2DReal &ThickFluxEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;

   if (LocThicknessFluxDiv.Enabled) {
      Pacer::start("Tend:thicknessFluxDiv", 2);
      parallelForOuter(
          {Mesh->NCellsAll}, KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocThicknessFluxDiv(LocLayerThicknessTend, ICell, KChunk,
                                        ThickFluxEdge, NormalVelEdge);
                 });
          });
      Pacer::stop("Tend:thicknessFluxDiv", 2);
   }

   if (CustomThicknessTend) {
      Pacer::start("Tend:customThicknessTend", 2);
      CustomThicknessTend(LocLayerThicknessTend, State, AuxState,
                          ThickTimeLevel, VelTimeLevel, Time);
      Pacer::stop("Tend:customThicknessTend", 2);
   }

   Pacer::stop("Tend:computeThicknessTendenciesOnly", 1);

} // end thickness tendency compute

//------------------------------------------------------------------------------
// Compute tendencies for normal velocity equation
void Tendencies::computeVelocityTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {

   OMEGA_SCOPE(LocNormalVelocityTend, NormalVelocityTend);
   OMEGA_SCOPE(LocPotientialVortHAdv, PotientialVortHAdv);
   OMEGA_SCOPE(LocKEGrad, KEGrad);
   OMEGA_SCOPE(LocSSHGrad, SSHGrad);
   OMEGA_SCOPE(LocVelocityDiffusion, VelocityDiffusion);
   OMEGA_SCOPE(LocVelocityHyperDiff, VelocityHyperDiff);
   OMEGA_SCOPE(LocWindForcing, WindForcing);
   OMEGA_SCOPE(LocBottomDrag, BottomDrag);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computeVelocityTendenciesOnly", 1);

   parallelForOuter(
       {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRange(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const int K                     = KMin + KChunk;
                 LocNormalVelocityTend(IEdge, K) = 0;
              });
       });

   // Compute potential vorticity horizontal advection
   const Array2DReal &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   const Array2DReal &NormRVortEdge = AuxState->VorticityAux.NormRelVortEdge;
   const Array2DReal &NormFEdge     = AuxState->VorticityAux.NormPlanetVortEdge;
   Array2DReal NormVelEdge          = State->getNormalVelocity(VelTimeLevel);
   if (LocPotientialVortHAdv.Enabled) {
      Pacer::start("Tend:potientialVortHAdv", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocPotientialVortHAdv(LocNormalVelocityTend, IEdge, KChunk,
                                          NormRVortEdge, NormFEdge,
                                          FluxLayerThickEdge, NormVelEdge);
                 });
          });
      Pacer::stop("Tend:potientialVortHAdv", 2);
   }

   // Compute kinetic energy gradient
   const Array2DReal &KECell = AuxState->KineticAux.KineticEnergyCell;
   if (LocKEGrad.Enabled) {
      Pacer::start("Tend:KEGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocKEGrad(LocNormalVelocityTend, IEdge, KChunk, KECell);
                 });
          });
      Pacer::stop("Tend:KEGrad", 2);
   }

   // Compute sea surface height gradient
   const Array2DReal &SSHCell = AuxState->LayerThicknessAux.SshCell;
   if (LocSSHGrad.Enabled) {
      Pacer::start("Tend:SSHGrad", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocSSHGrad(LocNormalVelocityTend, IEdge, KChunk, SSHCell);
                 });
          });
      Pacer::stop("Tend:SSHGrad", 2);
   }

   // Compute del2 horizontal diffusion
   const Array2DReal &DivCell     = AuxState->KineticAux.VelocityDivCell;
   const Array2DReal &RVortVertex = AuxState->VorticityAux.RelVortVertex;
   if (LocVelocityDiffusion.Enabled) {
      Pacer::start("Tend:velocityDiffusion", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocVelocityDiffusion(LocNormalVelocityTend, IEdge, KChunk,
                                         DivCell, RVortVertex);
                 });
          });
      Pacer::stop("Tend:velocityDiffusion", 2);
   }

   // Compute del4 horizontal diffusion
   const Array2DReal &Del2DivCell = AuxState->VelocityDel2Aux.Del2DivCell;
   const Array2DReal &Del2RVortVertex =
       AuxState->VelocityDel2Aux.Del2RelVortVertex;
   if (LocVelocityHyperDiff.Enabled) {
      Pacer::start("Tend:velocityHyperDiff", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocVelocityHyperDiff(LocNormalVelocityTend, IEdge, KChunk,
                                         Del2DivCell, Del2RVortVertex);
                 });
          });
      Pacer::stop("Tend:velocityHyperDiff", 2);
   }

   // Compute wind forcing
   const auto &NormalStressEdge = AuxState->WindForcingAux.NormalStressEdge;
   const auto &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   if (LocWindForcing.Enabled) {
      Pacer::start("Tend:windForcing", 2);
      parallelForOuter(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = MinLayerEdgeBot(IEdge);
             const int KMax   = MaxLayerEdgeTop(IEdge);
             const int KRange = vertRangeChunked(KMin, KMax);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocWindForcing(LocNormalVelocityTend, IEdge, KChunk,
                                   NormalStressEdge, MeanLayerThickEdge);
                 });
          });
      Pacer::stop("Tend:windForcing", 2);
   }

   // Compute bottom drag
   if (LocBottomDrag.Enabled) {
      Pacer::start("Tend:bottomDrag", 2);
      parallelFor(
          {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
             LocBottomDrag(LocNormalVelocityTend, IEdge, NormVelEdge, KECell,
                           MeanLayerThickEdge);
          });
      Pacer::stop("Tend:bottomDrag", 2);
   }

   if (CustomVelocityTend) {
      Pacer::start("Tend:customVelocityTend", 2);
      CustomVelocityTend(LocNormalVelocityTend, State, AuxState, ThickTimeLevel,
                         VelTimeLevel, Time);
      Pacer::stop("Tend:customVelocityTend", 2);
   }

   Pacer::stop("Tend:computeVelocityTendenciesOnly", 1);

} // end velocity tendency compute

void Tendencies::computeTracerTendenciesOnly(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   OMEGA_SCOPE(LocTracerTend, TracerTend);
   OMEGA_SCOPE(LocTracerHorzAdv, TracerHorzAdv);
   OMEGA_SCOPE(LocTracerHighOrderHorzAdv, TracerHighOrderHorzAdv);
   OMEGA_SCOPE(LocTracerDiffusion, TracerDiffusion);
   OMEGA_SCOPE(LocTracerHyperDiff, TracerHyperDiff);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);

   Pacer::start("Tend:computeTracerTendenciesOnly", 1);

   parallelForOuter(
       {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRange(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const int K                = KMin + KChunk;
                 LocTracerTend(L, ICell, K) = 0;
              });
       });

   // compute tracer horizotal advection
   const Array2DReal &NormalVelEdge = State->NormalVelocity[VelTimeLevel];
   const Array3DReal &HTracersEdge  = AuxState->TracerAux.HTracersEdge;
   const Array2DReal &FluxLayerThickEdge =
       AuxState->LayerThicknessAux.FluxLayerThickEdge;
   if (LocTracerHorzAdv.Enabled) {
      Pacer::start("Tend:tracerHorzAdv", 2);
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerHorzAdv(LocTracerTend, L, ICell, KChunk,
                                     NormalVelEdge, HTracersEdge);
                 });
          });
      Pacer::stop("Tend:tracerHorzAdv", 2);
   }
   if (LocTracerHighOrderHorzAdv.Enabled) {
      parallelFor(
          {NTracers, NEdgesAll, NChunks},
          KOKKOS_LAMBDA(int L, int IEdge, int KChunk) {
             LocTracerHighOrderHorzAdv(L, IEdge, KChunk, TracerArray,
                                       FluxLayerThickEdge, NormalVelEdge);
          });
      parallelFor(
          {NTracers, NCellsAll, NChunks},
          KOKKOS_LAMBDA(int L, int ICell, int KChunk) {
             LocTracerHighOrderHorzAdv(LocTracerTend, L, ICell, KChunk);
          });
   }

   // compute tracer diffusion
   const Array2DReal &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   if (LocTracerDiffusion.Enabled) {
      Pacer::start("Tend:tracerDiffusion", 2);
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerDiffusion(LocTracerTend, L, ICell, KChunk,
                                       TracerArray, MeanLayerThickEdge);
                 });
          });
      Pacer::stop("Tend:tracerDiffusion", 2);
   }

   // compute tracer hyperdiffusion
   const Array3DReal &Del2TracersCell = AuxState->TracerAux.Del2TracersCell;
   if (LocTracerHyperDiff.Enabled) {
      Pacer::start("Tend:tracerHyperDiff", 2);
      parallelForOuter(
          {NTracers, Mesh->NCellsAll},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const int KMin   = MinLayerCell(ICell);
             const int KMax   = MaxLayerCell(ICell);
             const int KRange = vertRangeChunked(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocTracerHyperDiff(LocTracerTend, L, ICell, KChunk,
                                       Del2TracersCell);
                 });
          });
      Pacer::stop("Tend:tracerHyperDiff", 2);
   }

   Pacer::stop("Tend:computeTracerTendenciesOnly", 1);
} // end tracer tendency compute

void Tendencies::computeThicknessTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   // only need LayerThicknessAux on edge
   Array2DReal LayerThick = State->getLayerThickness(ThickTimeLevel);
   Array2DReal NormVel    = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(LayerThicknessAux, AuxState->LayerThicknessAux);
   OMEGA_SCOPE(LayerThickCell, LayerThick);
   OMEGA_SCOPE(NormalVelEdge, NormVel);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computeThicknessTendencies", 1);

   Pacer::start("Tend:computeLayerThickAux", 2);
   parallelForOuter(
       "computeLayerThickAux", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 LayerThicknessAux.computeVarsOnEdge(
                     IEdge, KChunk, LayerThickCell, NormalVelEdge);
              });
       });
   Pacer::stop("Tend:computeLayerThickAux", 2);

   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);

   Pacer::stop("Tend:computeThicknessTendencies", 1);
}

void Tendencies::computeVelocityTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   Pacer::start("Tend:computeVelocityTendencies", 1);

   AuxState->computeMomAux(State, ThickTimeLevel, VelTimeLevel);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);

   Pacer::stop("Tend:computeVelocityTendencies", 1);
}

void Tendencies::computeTracerTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   Array2DReal LayerThickCell = State->getLayerThickness(ThickTimeLevel);
   Array2DReal NormalVelEdge  = State->getNormalVelocity(VelTimeLevel);
   OMEGA_SCOPE(TracerAux, AuxState->TracerAux);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);

   Pacer::start("Tend:computeTracerTendencies", 1);

   Pacer::start("Tend:computeTracerAuxEdge", 2);
   parallelForOuter(
       "computeTracerAuxEdge", {NTracers, Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int LTracer, int IEdge, const TeamMember &Team) {
          const int KMin   = MinLayerEdgeBot(IEdge);
          const int KMax   = MaxLayerEdgeTop(IEdge);
          const int KRange = vertRangeChunked(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 TracerAux.computeVarsOnEdge(LTracer, IEdge, KChunk,
                                             NormalVelEdge, LayerThickCell,
                                             TracerArray);
              });
       });
   Pacer::stop("Tend:computeTracerAuxEdge", 2);

   const auto &MeanLayerThickEdge =
       AuxState->LayerThicknessAux.MeanLayerThickEdge;
   Pacer::start("Tend:computeTracerAuxCell", 2);
   parallelForOuter(
       "computeTracerAuxCell", {NTracers, Mesh->NCellsAll},
       KOKKOS_LAMBDA(int LTracer, int ICell, const TeamMember &Team) {
          const int KMin   = MinLayerCell(ICell);
          const int KMax   = MaxLayerCell(ICell);
          const int KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 TracerAux.computeVarsOnCells(LTracer, ICell, KChunk,
                                              MeanLayerThickEdge, TracerArray);
              });
       });
   Pacer::stop("Tend:computeTracerAuxCell", 2);

   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);

   Pacer::stop("Tend:computeTracerTendencies", 1);
}

//------------------------------------------------------------------------------
// Compute both layer thickness and normal velocity tendencies
void Tendencies::computeAllTendencies(
    const OceanState *State,        ///< [in] State variables
    const AuxiliaryState *AuxState, ///< [in] Auxilary state variables
    const Array3DReal &TracerArray, ///< [in] Tracer array
    int ThickTimeLevel,             ///< [in] Time level
    int VelTimeLevel,               ///< [in] Time level
    TimeInstant Time                ///< [in] Time
) {
   AuxState->computeAll(State, TracerArray, ThickTimeLevel, VelTimeLevel);
   computeThicknessTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                  Time);
   computeVelocityTendenciesOnly(State, AuxState, ThickTimeLevel, VelTimeLevel,
                                 Time);
   computeTracerTendenciesOnly(State, AuxState, TracerArray, ThickTimeLevel,
                               VelTimeLevel, Time);
} // end all tendency compute

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
