//===-- ocn/VertAdv.cpp - vertical advection --------------------*- C++ -*-===//
//
// The VertAdv class contains methods for calculating the vertical
// pseudo-velocity and the tendencies of thickness, horizontal velocity, and
// tracers due to vertical advective transport, and member variables for storing
// the vertical pseudo-velocity and the fluxes of tracers at vertical interfaces
// in each column.
//
//===----------------------------------------------------------------------===//

#include "VertAdv.h"
#include "Field.h"
#include "Tracers.h"

#include <limits>

namespace OMEGA {

// create static class members
VertAdv *VertAdv::DefaultVertAdv = nullptr;
std::map<std::string, std::unique_ptr<VertAdv>> VertAdv::AllVertAdvs;

//------------------------------------------------------------------------------
// create the default VertAdv, requires prior initialization of default
// HorzMesh, VertCoord, and Tracers
void VertAdv::init() {

   auto Mesh   = HorzMesh::getDefault();
   auto VCoord = VertCoord::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();

   VertAdv::DefaultVertAdv = create("Default", Mesh, VCoord, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// constructor
VertAdv::VertAdv(const std::string &Name_,  //< [in] name for new VertAdv
                 const HorzMesh *InMesh,    //< [in] associated HorzMesh
                 const VertCoord *InVCoord, //< [in] associated VertCoord
                 Config *Options            //< [in] configuration options
) {

   // Read options from Config object
   readConfigOptions(Options);

   // Store name suffix
   Name = Name_;

   // Store pointers to Mesh and VertCoord objects
   Mesh   = InMesh;
   VCoord = InVCoord;

   // Store dimension sizes
   NVertLayers   = VCoord->NVertLayers;
   NVertLayersP1 = VCoord->NVertLayersP1;
   NCellsOwned   = Mesh->NCellsOwned;
   NCellsHalo0   = Mesh->NCellsHaloH(0);
   NCellsAll     = Mesh->NCellsAll;
   NCellsSize    = Mesh->NCellsSize;
   NEdgesOwned   = Mesh->NEdgesOwned;
   NEdgesAll     = Mesh->NEdgesAll;
   NEdgesSize    = Mesh->NEdgesSize;
   NTracers      = Tracers::getNumTracers();

   // Allocate member arrays
   VerticalPseudoVelocity =
       Array2DReal("VerticalPseudoVelocity", NCellsSize, NVertLayersP1);
   TotalVerticalPseudoVelocity =
       Array2DReal("TotalVerticalPseudoVelocity", NCellsSize, NVertLayersP1);
   VertFlux = Array3DReal("VertFlux", NTracers, NCellsSize, NVertLayersP1);

   // Allocate host copies
   VerticalPseudoVelocityH = createHostMirrorCopy(VerticalPseudoVelocity);
   TotalVerticalPseudoVelocityH =
       createHostMirrorCopy(TotalVerticalPseudoVelocity);
   VertFluxH = createHostMirrorCopy(VertFlux);

   // Low-order flux array only needed for flux-corrected transport
   if (VertAdvChoice == VertAdvOption::FCT) {
      LowOrderVertFlux =
          Array3DReal("LowOrderVertFlux", NTracers, NCellsSize, NVertLayersP1);
      LowOrderVertFluxH = createHostMirrorCopy(LowOrderVertFlux);
   }

   defineFields();

} // end constructor

//------------------------------------------------------------------------------
// create a new VertAdv instance
VertAdv *VertAdv::create(const std::string &Name, //< [in] name for new VertAdv
                         const HorzMesh *Mesh,    //< [in] associated HorzMesh
                         const VertCoord *VCoord, //< [in] associated VertCoord
                         Config *Options          //< [in] configuration options
) {
   // Check to see if a VertAdv of the same name already exists and, if so,
   // exit with an error
   if (AllVertAdvs.find(Name) != AllVertAdvs.end()) {
      LOG_ERROR("Attempted to create a VertAdv with name {} but a VertAdv "
                "of that name already exists",
                Name);
      return nullptr;
   }

   // create a new VertAdv on the heap and put it in a map of unique_ptrs,
   // which will manage its lifetime
   auto *NewVertAdv = new VertAdv(Name, Mesh, VCoord, Options);
   AllVertAdvs.emplace(Name, NewVertAdv);

   return NewVertAdv;

} // end create

void VertAdv::defineFields() {

   // set field names (append Name if not default)
   VerticalPseudoVelocityFldName  = "VerticalPseudoVelocity";
   TotalVertPseudoVelocityFldName = "TotalVerticalPseudoVelocity";
   VertFluxFldName                = "VertFlux";
   LowOrderVertFluxFldName        = "LowOrderVertFlux";

   if (Name != "Default") {
      VerticalPseudoVelocityFldName.append(Name);
      TotalVertPseudoVelocityFldName.append(Name);
      VertFluxFldName.append(Name);
      LowOrderVertFluxFldName.append(Name);
   }

   const Real FillValueReal = -9.99e30;
   I4 NDims                 = 2;
   std::vector<std::string> DimNames(NDims);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLayersP1";

   auto VerticalPseudoVelocityField = Field::create(
       VerticalPseudoVelocityFldName,                             // field name
       "Vertical pseudo-velocity across a pseudo-height surface", // long name
                                                                  // or
                                                                  // description
       "m s^-1",                                                  // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValueReal,                    // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   auto TotalVertPseudoVelocityField =
       Field::create(TotalVertPseudoVelocityFldName, // field name
                     "Total vertical pseudo-velocity across a moving, tilted "
                     "pseudo-height surface", // long name or description
                     "m s^-1",                // units
                     "",                      // CF standard Name
                     std::numeric_limits<Real>::min(), // min valid value
                     std::numeric_limits<Real>::max(), // max valid value
                     FillValueReal, // scalar for undefined entries
                     NDims,         // number of dimensions
                     DimNames       // dimension names
       );

   NDims = 3;
   DimNames.resize(NDims);
   DimNames[2] = "NTracers";

   auto VertFluxField = Field::create(
       VertFluxFldName,                                           // field name
       "Vertical flux of tracers across a pseudo-height surface", // long name
                                                                  // or
                                                                  // description
       "",                                                        // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValueReal,                    // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   auto LowOrderVertFluxField = Field::create(
       LowOrderVertFluxFldName, // field name
       "Low-order vertical flux of tracers across a pseudo-height "
       "surface",                        // long name or description
       "",                               // units
       "",                               // CF standard Name
       std::numeric_limits<Real>::min(), // min valid value
       std::numeric_limits<Real>::max(), // max valid value
       FillValueReal,                    // scalar for undefined entries
       NDims,                            // number of dimensions
       DimNames                          // dimension names
   );

   GroupName = "VertAdv";
   if (Name != "Default") {
      GroupName.append(Name);
   }

   // Create a field group for VertAdv fields
   auto VertAdvGroup = FieldGroup::create(GroupName);

   VertAdvGroup->addField(VerticalPseudoVelocityFldName);
   VertAdvGroup->addField(TotalVertPseudoVelocityFldName);
   VertAdvGroup->addField(VertFluxFldName);
   VertAdvGroup->addField(LowOrderVertFluxFldName);

   // Associate Fields with data
   VerticalPseudoVelocityField->attachData<Array2DReal>(VerticalPseudoVelocity);
   TotalVertPseudoVelocityField->attachData<Array2DReal>(
       TotalVerticalPseudoVelocity);
   VertFluxField->attachData<Array3DReal>(VertFlux);
   LowOrderVertFluxField->attachData<Array3DReal>(LowOrderVertFlux);

} // end defineFields

//------------------------------------------------------------------------------
// destructor
VertAdv::~VertAdv() {

   if (FieldGroup::exists(GroupName)) {
      Field::destroy(VerticalPseudoVelocityFldName);
      Field::destroy(TotalVertPseudoVelocityFldName);
      Field::destroy(VertFluxFldName);
      Field::destroy(LowOrderVertFluxFldName);
      FieldGroup::destroy(GroupName);
   }

} // end destructor

//------------------------------------------------------------------------------
// Removes all VertAdvs to clean up before exit
void VertAdv::clear() {

   AllVertAdvs.clear(); // removes all VertAdvs from the list and in the
                        // process, calls the destructors for each

} // end clear

//------------------------------------------------------------------------------
// Removes a VertAdv from map by name
void VertAdv::erase(std::string Name) {
   AllVertAdvs.erase(Name); // removes the VertAdv from the list and in the
                            // process, calls the destructor
} // end erase

//------------------------------------------------------------------------------
// Perform deepCopy for each variable array from device to host
void VertAdv::copyToHost() {

   deepCopy(VerticalPseudoVelocityH, VerticalPseudoVelocity);
   deepCopy(TotalVerticalPseudoVelocityH, TotalVerticalPseudoVelocity);
   deepCopy(VertFluxH, VertFlux);
   if (VertAdvChoice == VertAdvOption::FCT) {
      deepCopy(LowOrderVertFluxH, LowOrderVertFlux);
   }
}

//------------------------------------------------------------------------------
// Perform deepCopy for each variable array from host to device
void VertAdv::copyToDevice() {

   deepCopy(VerticalPseudoVelocity, VerticalPseudoVelocityH);
   deepCopy(TotalVerticalPseudoVelocity, TotalVerticalPseudoVelocityH);
   deepCopy(VertFlux, VertFluxH);
   if (VertAdvChoice == VertAdvOption::FCT) {
      deepCopy(LowOrderVertFlux, LowOrderVertFluxH);
   }
}

//------------------------------------------------------------------------------
// Get default VertAdv
VertAdv *VertAdv::getDefault() { return VertAdv::DefaultVertAdv; }

//------------------------------------------------------------------------------
// Get VertAdv by name
VertAdv *VertAdv::get(const std::string Name //< [in] Name of VertAdv
) {

   // look for an instance of this name
   auto it = AllVertAdvs.find(Name);

   // if found, return the VertAdv pointer
   if (it != AllVertAdvs.end()) {
      return it->second.get();

      // otherwise print error and return null pointer
   } else {
      LOG_ERROR("VertAdv::get: Attempt to retrieve non-existant VertAdv:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }

} // end get VertAdv

//------------------------------------------------------------------------------
// Read and set config options
void VertAdv::readConfigOptions(Config *OmegaConfig) {

   Error Err; // Error code

   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "VertAdv: Tendencies group not in Config");

   Err += TendConfig.get("ThicknessVertAdvTendencyEnable",
                         this->ThickVertAdvEnabled);
   CHECK_ERROR_ABORT(
       Err, "VertAdv: ThicknessVertAdvTendencyEnable not found in TendConfig");

   Err +=
       TendConfig.get("VelocityVertAdvTendencyEnable", this->VelVertAdvEnabled);
   CHECK_ERROR_ABORT(
       Err, "VertAdv: ThicknessVertAdvTendencyEnable not found in TendConfig");

   Err += TendConfig.get("TracerVertAdvTendencyEnable",
                         this->TracerVertAdvEnabled);
   CHECK_ERROR_ABORT(
       Err, "VertAdv: ThicknessVertAdvTendencyEnable not found in TendConfig");

   Config AdvectConfig("Advection");
   Err += OmegaConfig->get(AdvectConfig);
   CHECK_ERROR_ABORT(Err, "VertAdv: Advection group not in Config");

   bool FluxLimiterOn;
   Err += AdvectConfig.get("VerticalTracerFluxLimiterEnable", FluxLimiterOn);
   CHECK_ERROR_ABORT(
       Err,
       "VertAdv: VerticalTracerFluxLimiterEnable not found in AdvectConfig");
   if (FluxLimiterOn) {
      VertAdvChoice = VertAdvOption::FCT;
   } else {
      VertAdvChoice = VertAdvOption::Standard;
   }

   I4 VertFluxOrder;
   Err += AdvectConfig.get("VerticalTracerFluxOrder", VertFluxOrder);
   CHECK_ERROR_ABORT(
       Err, "VertAdv: VerticalTracerFluxOrder not found in AdvectConfig");

   switch (VertFluxOrder) {
   case (2):
      VertFluxChoice = VertFluxOption::Second;
      break;
   case (3):
      VertFluxChoice = VertFluxOption::Third;
      break;
   case (4):
      VertFluxChoice = VertFluxOption::Fourth;
      break;
   default:
      ABORT_ERROR("VertAdv: Invalid option for VerticalTracerFluxOrder found "
                  "in AdvectConfig. Must be 2, 3, or 4");
   }

   Err += AdvectConfig.get("Coef3rdOrder", Coef3rdOrder);
   CHECK_ERROR_ABORT(Err, "VertAdv: Coef3rdOrder not found in AdvectConfig");

} // end readConfigOptions

//------------------------------------------------------------------------------
// Compute VerticalPseudoVelocity and TotalVerticalPseudoVelocity from the
// horizontal velocity (NormalVelocity), the pseudo-thickness used for fluxes
// through edges (FluxPseudoThickEdge), and the cell-based PseudoThickness.
void VertAdv::computeVerticalPseudoVelocity(
    const Array2DReal &NormalVelocity,      //< [in] horizontal velocity
    const Array2DReal &FluxPseudoThickEdge, //< [in] pseudo-thickness at edges
    const Array2DReal &PseudoThickness,     //< [in] pseudo-thickness of layer
    const Real Dt                           //< [in] time interval
) {

   // Return if mesh only has a single vertical layer
   if (NVertLayers == 1)
      return;

   OMEGA_SCOPE(LocVertVel, VerticalPseudoVelocity);
   OMEGA_SCOPE(LocTotVertVel, TotalVerticalPseudoVelocity);
   OMEGA_SCOPE(LocPseudoThickTarget, VCoord->PseudoThicknessTarget);
   OMEGA_SCOPE(LocNVertLayers, NVertLayers);
   OMEGA_SCOPE(LocAreaCell, Mesh->AreaCell);
   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocNEOnC, Mesh->NEdgesOnCell);
   OMEGA_SCOPE(LocEOnC, Mesh->EdgesOnCell);
   OMEGA_SCOPE(LocDvE, Mesh->DvEdge);
   OMEGA_SCOPE(LocESOnC, Mesh->EdgeSignOnCell);

   // Loop over all cells owned by the task
   parallelForOuter(
       "computeVerticalPseudoVelocity",
       LaunchConfig({NCellsHalo0}, TeamScratch<Real>(NVertLayers)),
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          ScratchArray1DReal DivHU(teamScratch(Team), LocNVertLayers);

          const Real InvAreaCell = 1._Real / LocAreaCell(ICell);

          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          I4 KRange     = vertRangeChunked(KMin, KMax);

          // Compute thickness-weighted divergence of horizontal velocity
          // in each layer
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 Real DivHUTmp[VecLength] = {0};
                 const I4 KStart          = chunkStart(KChunk, KMin);
                 const I4 KLen            = chunkLength(KChunk, KStart, KMax);

                 for (int J = 0; J < LocNEOnC(ICell); ++J) {
                    const I4 JEdge = LocEOnC(ICell, J);
                    for (int KVec = 0; KVec < KLen; ++KVec) {
                       const I4 K = KStart + KVec;
                       DivHUTmp[KVec] -= LocDvE(JEdge) * LocESOnC(ICell, J) *
                                         FluxPseudoThickEdge(JEdge, K) *
                                         NormalVelocity(JEdge, K) * InvAreaCell;
                    }
                 }
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    DivHU(K)   = DivHUTmp[KVec];
                 }
              });

          teamBarrier(Team);

          // Set velocity through top and bottom interfaces to zero
          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 LocVertVel(ICell, KMin)        = 0.;
                 LocVertVel(ICell, KMax + 1)    = 0.;
                 LocTotVertVel(ICell, KMin)     = 0.;
                 LocTotVertVel(ICell, KMax + 1) = 0.;
              });
          KRange = vertRange(KMin + 1, KMax);

          // Prefix sum of divergence to determine velocity through each
          // interface
          parallelScanInner(
              Team, KRange, INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KRev = KMax - K;
                 Accum -= DivHU(KRev);

                 if (IsFinal) {
                    LocVertVel(ICell, KRev) = Accum;
                 }
              });

          // Add contribution to transport pseudo-velocity from movement of
          // layer interfaces, store in TotalVerticalPseudoVelocity.
          parallelScanInner(
              Team, KRange, INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KRev      = KMax - K;
                 const Real AleTerm = (LocPseudoThickTarget(ICell, KRev) -
                                       PseudoThickness(ICell, KRev)) /
                                      Dt;

                 Accum -= DivHU(KRev) + AleTerm;

                 if (IsFinal) {
                    LocTotVertVel(ICell, KRev) = Accum;
                 }
              });
       });

} // end computeVerticalPseudoVelocity

//------------------------------------------------------------------------------
// Compute thickness tendency due to vertical advection
void VertAdv::computeThicknessVAdvTend(
    const Array2DReal &ThickTend //< [inout] thickness tendency
) {

   // Return if vertical advection thickness tendency not enabled
   if (!ThickVertAdvEnabled)
      return;

   // Return if mesh only has a single vertical layer
   if (NVertLayers == 1)
      return;

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocTotVertVelocity, TotalVerticalPseudoVelocity);

   // Loop over every owned cell, pseudo-thickness tendency is simply
   // difference in pseudo velocity between bottom and top interface for
   // each layer
   parallelForOuter(
       "computeThicknessVAdvTend", {NCellsOwned},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          const I4 KMin   = MinLayerCell(ICell);
          const I4 KMax   = MaxLayerCell(ICell);
          const I4 KRange = vertRangeChunked(KMin, KMax);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    ThickTend(ICell, K) += LocTotVertVelocity(ICell, K + 1) -
                                           LocTotVertVelocity(ICell, K);
                 }
              });
       });

} // end computeThicknessVAdvTend

//------------------------------------------------------------------------------
// Compute velocity tendency due to vertical advection
void VertAdv::computeVelocityVAdvTend(
    const Array2DReal &VelTend,        //< [inout] horizontal velocity tendency
    const Array2DReal &NormalVelocity, //< [in] horizontal velocity
    const Array2DReal &FluxPseudoThickEdge //< [in] pseudo-thickness at edges
) {

   // Return if vertical advection velocity tendency not enabled
   if (!VelVertAdvEnabled)
      return;

   // Return if mesh only has a single vertical layer
   if (NVertLayers == 1)
      return;

   OMEGA_SCOPE(LocCOnE, Mesh->CellsOnEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTotVertVelocity, TotalVerticalPseudoVelocity);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);

   // Loop over every owned edge
   parallelForOuter(
       "computeVelocityVAdvTend",
       LaunchConfig({NEdgesOwned}, TeamScratch<Real>(NVertLayersP1)),
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const I4 Cell1 = LocCOnE(IEdge, 0);
          const I4 Cell2 = LocCOnE(IEdge, 1);
          const I4 KMin  = MinLayerEdgeBot(IEdge);
          const I4 KMax  = MaxLayerEdgeTop(IEdge);
          I4 KRange      = vertRangeChunked(KMin + 1, KMax);

          // Allocate scratch space for W times Du/Dz at vertical interfaces
          // between edges
          ScratchArray1DReal WDuDzEdge(teamScratch(Team), LocNVertLayersP1);

          // Flux is zero at top and bottom
          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 WDuDzEdge(KMin)     = 0._Real;
                 WDuDzEdge(KMax + 1) = 0._Real;
              });

          // Average vertical velocities from cell centers to edges and multiply
          // by derivative of horizontal velocity to obtain flux of horizontal
          // velocity at the vertical interfaces between edges
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin + 1);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);

                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K      = KStart + KVec;
                    const Real WAvg = 0.5_Real * (LocTotVertVelocity(Cell1, K) +
                                                  LocTotVertVelocity(Cell2, K));
                    WDuDzEdge(K) =
                        WAvg *
                        (NormalVelocity(IEdge, K - 1) -
                         NormalVelocity(IEdge, K)) /
                        (0.5_Real * (FluxPseudoThickEdge(IEdge, K - 1) +
                                     FluxPseudoThickEdge(IEdge, K)));
                 }
              });

          teamBarrier(Team);

          KRange = vertRangeChunked(KMin, KMax);
          // Average W*Du/Dz from interfaces to layer midpoints
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);

                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    VelTend(IEdge, K) -= EdgeMask(IEdge, K) * 0.5_Real *
                                         (WDuDzEdge(K) + WDuDzEdge(K + 1));
                 }
              });
       });

} // end computeVelocityVAdvTend

//------------------------------------------------------------------------------
// Compute tracer tendency due to vertical advection, TimeStep is only needed
// as an argument for flux-corrected transport
void VertAdv::computeTracerVAdvTend(
    const Array3DReal &TracerTend,      //< [inout] tracer tendencies
    const Array3DReal &Tracers,         //< [in] tracer array
    const Array2DReal &PseudoThickness, //< [in] pseudo-thickness
    const TimeInterval TimeStep         //< [in] time step
) {

   // Return if vertical advection tracer tendency not enabled
   if (!TracerVertAdvEnabled)
      return;

   // Return if mesh only has a single vertical layer
   if (NVertLayers == 1)
      return;

   // Compute tracer fluxes at the interfaces
   computeVerticalFluxes(Tracers, PseudoThickness);

   // Dispatch to appropriate algorithm based on configuration settings
   switch (VertAdvChoice) {
   case VertAdvOption::Standard:
      computeStdVAdvTend(TracerTend);
      break;
   case VertAdvOption::FCT:
      R8 Dt;
      TimeStep.get(Dt, TimeUnits::Seconds);
      computeFCTVAdvTend(TracerTend, Tracers, PseudoThickness, Dt);
      break;
   }

} // end computeTracerVAdvTend

//------------------------------------------------------------------------------
// Compute tracer fluxes due to vertical advection, the particular scheme used
// is chosen via configuration settings
void VertAdv::computeVerticalFluxes(
    const Array3DReal &Tracers,        //< [in] tracer array
    const Array2DReal &PseudoThickness //< [in] pseudo-thickness
) {

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocTotVertVel, TotalVerticalPseudoVelocity);
   OMEGA_SCOPE(LocVertFlux, VertFlux);
   OMEGA_SCOPE(LocLowOrderVertFlux, LowOrderVertFlux);
   OMEGA_SCOPE(LocCoef3rdOrder, Coef3rdOrder);

   // Compute the fluxes used for the standard VAdv scheme, or the high-order
   // fluxes used for the FCT VAdv scheme, store in VertFlux member array
   switch (VertFluxChoice) {
   // 2nd-order centered fluxes
   case VertFluxOption::Second:
      parallelForOuter(
          "computeVerticalFluxes-Second", {NTracers, NCellsOwned},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const I4 KMin   = MinLayerCell(ICell);
             const I4 KMax   = MaxLayerCell(ICell);
             const I4 KRange = vertRangeChunked(KMin + 2, KMax - 1);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const I4 KStart = chunkStart(KChunk, KMin + 2);
                    const I4 KLen   = chunkLength(KChunk, KStart, KMax - 1);
                    for (int KVec = 0; KVec < KLen; ++KVec) {
                       const I4 K = KStart + KVec;
                       const Real InvPseudoThickSum =
                           1._Real / (PseudoThickness(ICell, K - 1) +
                                      PseudoThickness(ICell, K));
                       const Real VerticalWeightK =
                           PseudoThickness(ICell, K - 1) * InvPseudoThickSum;
                       const Real VerticalWeightKm1 =
                           PseudoThickness(ICell, K) * InvPseudoThickSum;
                       LocVertFlux(L, ICell, K) =
                           LocTotVertVel(ICell, K) *
                           (VerticalWeightK * Tracers(L, ICell, K) +
                            VerticalWeightKm1 * Tracers(L, ICell, K - 1));
                    }
                 });
          });
      break;
   // 3rd-order upwind fluxes
   case VertFluxOption::Third:
      parallelForOuter(
          "computeVerticalFluxes-Third", {NTracers, NCellsOwned},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const I4 KMin   = MinLayerCell(ICell);
             const I4 KMax   = MaxLayerCell(ICell);
             const I4 KRange = vertRangeChunked(KMin + 2, KMax - 1);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const I4 KStart = chunkStart(KChunk, KMin + 2);
                    const I4 KLen   = chunkLength(KChunk, KStart, KMax - 1);
                    for (int KVec = 0; KVec < KLen; ++KVec) {
                       const I4 K = KStart + KVec;
                       LocVertFlux(L, ICell, K) =
                           (LocTotVertVel(ICell, K) *
                                (7._Real * (Tracers(L, ICell, K) +
                                            Tracers(L, ICell, K - 1)) -
                                 (Tracers(L, ICell, K + 1) +
                                  Tracers(L, ICell, K - 2))) -
                            LocCoef3rdOrder *
                                std::abs(LocTotVertVel(ICell, K)) *
                                ((Tracers(L, ICell, K + 1) -
                                  Tracers(L, ICell, K - 2)) -
                                 3._Real * (Tracers(L, ICell, K) -
                                            Tracers(L, ICell, K - 1)))) /
                           12._Real;
                    }
                 });
          });
      break;
   // 4th-order centered fluxes
   case VertFluxOption::Fourth:
      parallelForOuter(
          "computeVerticalFluxes-Fourth", {NTracers, NCellsOwned},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const I4 KMin   = MinLayerCell(ICell);
             const I4 KMax   = MaxLayerCell(ICell);
             const I4 KRange = vertRangeChunked(KMin + 2, KMax - 1);
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const I4 KStart = chunkStart(KChunk, KMin + 2);
                    const I4 KLen   = chunkLength(KChunk, KStart, KMax - 1);
                    for (int KVec = 0; KVec < KLen; ++KVec) {
                       const I4 K = KStart + KVec;
                       LocVertFlux(L, ICell, K) =
                           LocTotVertVel(ICell, K) *
                           (7._Real * (Tracers(L, ICell, K) +
                                       Tracers(L, ICell, K - 1)) -
                            (Tracers(L, ICell, K + 1) +
                             Tracers(L, ICell, K - 2))) /
                           12._Real;
                    }
                 });
          });
      break;
   }

   // Handle fluxes on interfaces near top and bottom layers. Fluxes are 0 on
   // top-most (KMin) and bottom-most (KMax + 1) interfaces, use second order
   // for fluxes on next-to-top (KMin + 1) and next-to-bottom (KMax) interfaces.
   parallelFor(
       "computeVerticalFluxes-TopBot", {NTracers, NCellsOwned},
       KOKKOS_LAMBDA(int L, int ICell) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          for (int K : {KMin, KMax + 1}) {
             LocVertFlux(L, ICell, K) = 0._Real;
          }
          for (int K : {KMin + 1, KMax}) {
             const Real InvPseudoThickSum =
                 1._Real /
                 (PseudoThickness(ICell, K - 1) + PseudoThickness(ICell, K));
             const Real VerticalWeightK =
                 PseudoThickness(ICell, K - 1) * InvPseudoThickSum;
             const Real VerticalWeightKm1 =
                 PseudoThickness(ICell, K) * InvPseudoThickSum;
             LocVertFlux(L, ICell, K) =
                 LocTotVertVel(ICell, K) *
                 (VerticalWeightK * Tracers(L, ICell, K) +
                  VerticalWeightKm1 * Tracers(L, ICell, K - 1));
          }
       });

   // If using FCT scheme, compute 1st-order upwind fluxes for low-order and
   // remove low-order flux from high-order flux
   if (VertAdvChoice == VertAdvOption::FCT) {
      parallelForOuter(
          "computeVerticalFluxes-LowOrder", {NTracers, NCellsOwned},
          KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
             const I4 KMin   = MinLayerCell(ICell);
             const I4 KMax   = MaxLayerCell(ICell);
             const I4 KRange = vertRangeChunked(KMin + 1, KMax);

             for (int K : {KMin, KMax + 1}) {
                LocLowOrderVertFlux(L, ICell, K) = 0._Real;
             }
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KChunk) {
                    const I4 KStart = chunkStart(KChunk, KMin + 1);
                    const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                    for (int KVec = 0; KVec < KLen; ++KVec) {
                       const I4 K = KStart + KVec;
                       LocLowOrderVertFlux(L, ICell, K) =
                           Kokkos::min(0._Real, LocTotVertVel(ICell, K)) *
                               Tracers(L, ICell, K - 1) +
                           Kokkos::max(0._Real, LocTotVertVel(ICell, K)) *
                               Tracers(L, ICell, K);

                       LocVertFlux(L, ICell, K) -=
                           LocLowOrderVertFlux(L, ICell, K);
                    }
                 });
          });
   }

} // end computeVerticalFluxes

//------------------------------------------------------------------------------
// Compute tracer tendencies due to vertical advection using standard advection
// scheme
void VertAdv::computeStdVAdvTend(
    const Array3DReal &TracerTend //< [inout] tracer tendencies
) {

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocVertFlux, VertFlux);

   // Loop over owned cells, tracer tendency in each layer is computed from
   // difference between fluxes through bottom and top interfaces
   parallelForOuter(
       "computeStdVAdvTend", {NTracers, NCellsOwned},
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const I4 KMin   = MinLayerCell(ICell);
          const I4 KMax   = MaxLayerCell(ICell);
          const I4 KRange = vertRangeChunked(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    TracerTend(L, ICell, K) +=
                        LocVertFlux(L, ICell, K + 1) - LocVertFlux(L, ICell, K);
                 }
              });
       });

} // end computeStdVAdvTend

//------------------------------------------------------------------------------
// Compute tracer tendencies due to vertical advection using flux-corrected
// transport scheme. ProvPseudoThickness input is provisional pseudo-thickness
// after horizontal pseudo-thickness flux
void VertAdv::computeFCTVAdvTend(
    const Array3DReal &TracerTend, //< [inout] tracer tendencies
    const Array3DReal &Tracers,    //< [in] tracer array
    const Array2DReal
        &ProvPseudoThickness, //< [in] provisional pseudo-thickness
    const Real Dt             //< [in] time step
) {

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocTotVertVel, TotalVerticalPseudoVelocity);
   OMEGA_SCOPE(LocVertFlux, VertFlux);
   OMEGA_SCOPE(LocLowOrderVertFlux, LowOrderVertFlux);
   OMEGA_SCOPE(LocNVertLayers, NVertLayers);
   OMEGA_SCOPE(LocEps, Eps);

   parallelForOuter(
       "computeFCTVAdvTend",
       LaunchConfig({NTracers, NCellsOwned},
                    TeamScratch<Real>(5 * NVertLayers + 1)),
       KOKKOS_LAMBDA(int L, int ICell, const TeamMember &Team) {
          const I4 KMin = MinLayerCell(ICell);
          const I4 KMax = MaxLayerCell(ICell);
          I4 KRange     = vertRangeChunked(KMin, KMax);

          ScratchArray1DReal InvNewProvThick(teamScratch(Team), LocNVertLayers);
          ScratchArray1DReal WorkTend(teamScratch(Team), LocNVertLayers);
          ScratchArray1DReal FlxIn(teamScratch(Team), LocNVertLayers);
          ScratchArray1DReal FlxOut(teamScratch(Team), LocNVertLayers);
          ScratchArray1DReal RescaledFlux(teamScratch(Team),
                                          LocNVertLayers + 1);

          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);

                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    InvNewProvThick(K) =
                        1._Real / (ProvPseudoThickness(ICell, K) +
                                   Dt * (LocTotVertVel(ICell, K + 1) -
                                         LocTotVertVel(ICell, K)));
                    Real TracerMax;
                    Real TracerMin;
                    // Determine bounds on tracer from neighbor values for
                    // limiting
                    if (K == KMin) {
                       TracerMax = Kokkos::max(Tracers(L, ICell, K),
                                               Tracers(L, ICell, K + 1));
                       TracerMin = Kokkos::min(Tracers(L, ICell, K),
                                               Tracers(L, ICell, K + 1));
                    } else if (K == KMax) {
                       TracerMax = Kokkos::max(Tracers(L, ICell, K - 1),
                                               Tracers(L, ICell, K));
                       TracerMin = Kokkos::min(Tracers(L, ICell, K - 1),
                                               Tracers(L, ICell, K));
                    } else {
                       TracerMax =
                           Kokkos::max(Tracers(L, ICell, K - 1),
                                       Kokkos::max(Tracers(L, ICell, K),
                                                   Tracers(L, ICell, K + 1)));
                       TracerMin =
                           Kokkos::min(Tracers(L, ICell, K - 1),
                                       Kokkos::min(Tracers(L, ICell, K),
                                                   Tracers(L, ICell, K + 1)));
                    }

                    // Accumulate upwind flux in WorkTend
                    WorkTend(K) = LocLowOrderVertFlux(L, ICell, K + 1) -
                                  LocLowOrderVertFlux(L, ICell, K);
                    // Accumulate remaining high-order flux into layer
                    FlxIn(K) =
                        Kokkos::max(0._Real, LocVertFlux(L, ICell, K + 1)) -
                        Kokkos::min(0._Real, LocVertFlux(L, ICell, K));
                    // Accumulate remaining high-order flux out of layer
                    FlxOut(K) =
                        Kokkos::min(0._Real, LocVertFlux(L, ICell, K + 1)) -
                        Kokkos::max(0._Real, LocVertFlux(L, ICell, K));
                    // Build scale factors to limit flux for FCT using the
                    // bounds determined above and bounds on newly updated
                    // values. Factors are stored in FlxIn and FlxOut
                    // scratch space.
                    Real TracerMinNew =
                        (Tracers(L, ICell, K) * ProvPseudoThickness(ICell, K) +
                         Dt * (WorkTend(K) + FlxOut(K))) *
                        InvNewProvThick(K);
                    Real TracerMaxNew =
                        (Tracers(L, ICell, K) * ProvPseudoThickness(ICell, K) +
                         Dt * (WorkTend(K) + FlxIn(K))) *
                        InvNewProvThick(K);
                    Real TracerUpwindNew =
                        (Tracers(L, ICell, K) * ProvPseudoThickness(ICell, K) +
                         Dt * WorkTend(K)) *
                        InvNewProvThick(K);
                    Real ScaleFactor =
                        (TracerMax - TracerUpwindNew) /
                        (TracerMaxNew - TracerUpwindNew + LocEps);
                    FlxIn(K) =
                        Kokkos::min(1._Real, Kokkos::max(0._Real, ScaleFactor));
                    ScaleFactor = (TracerUpwindNew - TracerMin) /
                                  (TracerUpwindNew - TracerMinNew + LocEps);
                    FlxOut(K) =
                        Kokkos::min(1._Real, Kokkos::max(0._Real, ScaleFactor));
                 }
              });

          teamBarrier(Team);

          KRange = vertRangeChunked(KMin + 1, KMax);

          // Rescale the high-order vertical flux
          RescaledFlux(KMin)     = 0._Real;
          RescaledFlux(KMax + 1) = 0._Real;
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin + 1);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);

                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;

                    RescaledFlux(K) =
                        Kokkos::max(0._Real, LocVertFlux(L, ICell, K)) *
                            Kokkos::min(FlxOut(K), FlxIn(K - 1)) +
                        Kokkos::min(0._Real, LocVertFlux(L, ICell, K)) *
                            Kokkos::min(FlxOut(K - 1), FlxIn(K));
                 }
              });

          teamBarrier(Team);

          // Accumulate total FCT vertical advection tendency
          KRange = vertRangeChunked(KMin, KMax);
          parallelForInner(
              Team, KRange, INNER_LAMBDA(int KChunk) {
                 const I4 KStart = chunkStart(KChunk, KMin);
                 const I4 KLen   = chunkLength(KChunk, KStart, KMax);

                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    WorkTend(K) += RescaledFlux(K + 1) - RescaledFlux(K);
                    TracerTend(L, ICell, K) += WorkTend(K);
                 }
              });
          // TODO: Monotonicity and diagnostic checks
       });

} // end computeFTCVAdvTend

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
