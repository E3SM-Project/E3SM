//===-- ocn/PGrad.cpp - Pressure Gradient Term -----------------*- C++ -*-===//
//
// Implements the PGrad manager and two discretizations: Centered and
// HighOrder.
//
//===----------------------------------------------------------------------===//

#include "PGrad.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "VertCoord.h"

namespace OMEGA {

PressureGrad *PressureGrad::DefaultPGrad = nullptr;
std::map<std::string, std::unique_ptr<PressureGrad>> PressureGrad::AllPGrads;

//------------------------------------------------------------------------------
// Initialize the PressureGrad. Assumes that HorzMesh and VertCoord have already
// been initialized.
void PressureGrad::init() {

   // Retrieve default mesh and vertical coordinate
   HorzMesh *DefMesh    = HorzMesh::getDefault();
   VertCoord *DefVCoord = VertCoord::getDefault();

   // Retrieve omega config
   Config *OmegaConfig = Config::getOmegaConfig();

   // Create the default PressureGrad and set pointer to it
   PressureGrad::DefaultPGrad =
       PressureGrad::create("Default", DefMesh, DefVCoord, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// Create a new PressureGrad object and add it to the map
PressureGrad *
PressureGrad::create(const std::string &Name, /// [in] Name for PressureGrad
                     const HorzMesh *Mesh,    ///< [in] Horizontal mesh
                     const VertCoord *VCoord, ///< [in] Vertical coordinate
                     Config *Options) {       ///< [in] Configuration options

   // Check to see if a PressureGrad of the same name already exists and
   // if so, exit with an error
   if (AllPGrads.find(Name) != AllPGrads.end()) {
      LOG_ERROR("Attempted to create a PressureGrad with name {} but a "
                "PressureGrad of that name already exists",
                Name);
      return nullptr;
   }

   // create a new PressureGrad on the heap and put it in a map of
   // unique_ptrs, which will manage its lifetime
   auto *NewPGrad = new PressureGrad(Mesh, VCoord, Options);
   AllPGrads.emplace(Name, NewPGrad);

   return NewPGrad;

} // end create

//------------------------------------------------------------------------------
// Get the default pressure gradient instance
PressureGrad *PressureGrad::getDefault() {

   return DefaultPGrad;

} // end get default

//------------------------------------------------------------------------------
// Constructor for PressureGrad
PressureGrad::PressureGrad(
    const HorzMesh *Mesh,    ///< [in] Horizontal mesh
    const VertCoord *VCoord, ///< [in] Vertical coordinate
    Config *Options)         ///< [in] Configuration options
    : MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop), CenteredPGrad(Mesh, VCoord),
      HighOrderPGrad(Mesh, VCoord) {

   // store mesh sizes
   NEdgesAll     = Mesh->NEdgesAll;
   NVertLayers   = VCoord->NVertLayers;
   NVertLayersP1 = NVertLayers + 1;
   NChunks       = NVertLayers / VecLength;

   // Read config options for PressureGrad type
   // and enable the appropriate functor
   Config PGradConfig("PressureGrad");
   Error Err;
   Err += Options->get(PGradConfig);
   CHECK_ERROR_ABORT(Err,
                     "PressureGrad: PressureGrad group not found in Config");
   std::string PGradTypeStr;
   Err += PGradConfig.get("PressureGradType", PGradTypeStr);

   if (PGradTypeStr == "centered" || PGradTypeStr == "Centered") {
      PressureGradChoice          = PressureGradType::Centered;
      this->CenteredPGrad.Enabled = true;
   } else if (PGradTypeStr == "HighOrder1") {
      PressureGradChoice           = PressureGradType::HighOrder1;
      this->HighOrderPGrad.Enabled = true;
   } else {
      LOG_INFO(
          "PGrad: Unknown PressureGradType in config, defaulting to centered");
   }

   // Temporary: initialization of tidal potential and SAL
   TidalPotential = Array1DReal("TidalPotential", Mesh->NCellsSize);
   SelfAttractionLoading =
       Array1DReal("SelfAttractionLoading", Mesh->NCellsSize);

} // end constructor

//------------------------------------------------------------------------------
// Destructor for PressureGrad
PressureGrad::~PressureGrad() {

   // No operations needed, Kokkos arrays removed when no longer in scope

} // end destructor

//------------------------------------------------------------------------------
// Remove PressureGrad instances before exit
void PressureGrad::clear() { AllPGrads.clear(); } // end clear

//------------------------------------------------------------------------------
// Remove PressureGrad from list by name
void PressureGrad::erase(const std::string &Name) {

   AllPGrads.erase(Name);

} // end erase

//------------------------------------------------------------------------------
// Get pressure gradient instance by name
PressureGrad *PressureGrad::get(const std::string &Name ///< [in] Name of
) {

   auto it = AllPGrads.find(Name);

   if (it != AllPGrads.end()) {
      return it->second.get();
   } else {
      LOG_ERROR("PressureGrad::get: Attempt to retrieve non-existent "
                "PressureGrad:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get pressure gradient

//------------------------------------------------------------------------------
// Compute pressure gradient tendencies and add into Tend array
void PressureGrad::computePressureGrad(Array2DReal &Tend,
                                       const OceanState *State,
                                       const VertCoord *VCoord,
                                       const Eos *EqState,
                                       const int TimeLevel) const {

   OMEGA_SCOPE(LocCenteredPGrad, CenteredPGrad);
   OMEGA_SCOPE(LocHighOrderPGrad, HighOrderPGrad);
   OMEGA_SCOPE(LocMinLayerEdgeBot, MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTidalPotential, TidalPotential);
   OMEGA_SCOPE(LocSelfAttractionLoading, SelfAttractionLoading);

   const Array2DReal &PressureMid       = VCoord->PressureMid;
   const Array2DReal &PressureInterface = VCoord->PressureInterface;
   const Array2DReal &Geopotential      = VCoord->GeopotentialMid;
   const Array2DReal &SpecVol           = EqState->SpecVol;
   const Array2DReal &ZInterface        = VCoord->ZInterface;
   const Array2DReal &ZMid              = VCoord->ZMid;
   Array2DReal LayerThick;
   State->getLayerThickness(LayerThick, TimeLevel);

   if (PressureGradChoice == PressureGradType::Centered) {

      // computes centered geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-centered", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocCenteredPGrad(Tend, IEdge, KChunk, PressureMid,
                                     PressureInterface, ZInterface,
                                     LocTidalPotential,
                                     LocSelfAttractionLoading, SpecVol);
                 });
          });

   } else {

      // computes high-order geopotential and pressure gradient tendency
      parallelForOuter(
          "pgrad-highorder", {NEdgesAll},
          KOKKOS_LAMBDA(I4 IEdge, const TeamMember &Team) {
             const int KMin   = LocMinLayerEdgeBot(IEdge);
             const int KMax   = LocMaxLayerEdgeTop(IEdge);
             const int KRange = vertRange(KMin, KMax);

             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    LocHighOrderPGrad(Tend, IEdge, KChunk, PressureMid,
                                      PressureInterface, ZInterface,
                                      LocTidalPotential,
                                      LocSelfAttractionLoading, SpecVol);
                 });
          });
   }
} // end compute pressure gradient

//------------------------------------------------------------------------------
// Constructor for centered pressure gradient functor
PressureGradCentered::PressureGradCentered(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

//------------------------------------------------------------------------------
// Constructor for high order pressure gradient functor
PressureGradHighOrder::PressureGradHighOrder(
    const HorzMesh *Mesh,   ///< [in] Horizontal mesh
    const VertCoord *VCoord ///< [in] Vertical coordinate
    )
    : CellsOnEdge(Mesh->CellsOnEdge), DcEdge(Mesh->DcEdge),
      EdgeMask(VCoord->EdgeMask), MinLayerEdgeBot(VCoord->MinLayerEdgeBot),
      MaxLayerEdgeTop(VCoord->MaxLayerEdgeTop) {}

} // namespace OMEGA
