//===-- ocn/VertAdv.cpp - vertical advection --------------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//

#include "VertAdv.h"
#include "Tracers.h"

namespace OMEGA {

// create static class members
VertAdv *VertAdv::DefaultVertAdv = nullptr;
std::map<std::string, std::unique_ptr<VertAdv>> VertAdv::AllVertAdvs;

//------------------------------------------------------------------------------
// create the default VertAdv, requires prior initialization of default
// HorzMesh and VertCoord
void VertAdv::init() {

   auto Mesh   = HorzMesh::getDefault();
   auto VCoord = VertCoord::getDefault();

   Config *OmegaConfig = Config::getOmegaConfig();

   VertAdv::DefaultVertAdv = create("Default", Mesh, VCoord, OmegaConfig);

} // end init

//------------------------------------------------------------------------------
// constructor
VertAdv::VertAdv(const std::string &Name_,  ///< [in] name for new VertAdv
                 const HorzMesh *InMesh,    ///< [in] associated HorzMesh
                 const VertCoord *InVCoord, ///< [in] associated VertCoord
                 Config *Options            ///< [in] configuration options
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
   NCellsAll     = Mesh->NCellsAll;
   NCellsSize    = Mesh->NCellsSize;
   NEdgesOwned   = Mesh->NEdgesOwned;
   NEdgesAll     = Mesh->NEdgesAll;
   NEdgesSize    = Mesh->NEdgesSize;
   NTracers      = Tracers::getNumTracers();

   // Allocate member arrays
   VerticalVelocity =
       Array2DReal("VerticalVelocity", NCellsSize, NVertLayersP1);
   TotalVerticalVelocity =
       Array2DReal("TotalVerticalVelocity", NCellsSize, NVertLayersP1);
   VertFlux = Array3DReal("VertFlux", NTracers, NCellsSize, NVertLayersP1);

   // Low-order flux array only needed for flux-corrected transport
   if (VertAdvChoice == VertAdvOption::FCT) {
      LowOrderVertFlux =
          Array3DReal("LowOrderVertFlux", NTracers, NCellsSize, NVertLayersP1);
   }

} // end constructor

//------------------------------------------------------------------------------
// create a new VertAdv instance
VertAdv *VertAdv::create(const std::string &Name, ///< [in] name for new VertAdv
                         const HorzMesh *Mesh,    ///< [in] associated HorzMesh
                         const VertCoord *VCoord, ///< [in] associated VertCoord
                         Config *Options ///< [in] configuration options
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

//------------------------------------------------------------------------------
// destructor
VertAdv::~VertAdv() {} // end destructor

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
// Get default VertAdv
VertAdv *VertAdv::getDefault() { return VertAdv::DefaultVertAdv; }

//------------------------------------------------------------------------------
// Get VertAdv by name
VertAdv *VertAdv::get(const std::string Name ///< [in] Name of VertAdv
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
   Err += AdvectConfig.get("VerticalTracerFluxLimiterEnabled", FluxLimiterOn);
   CHECK_ERROR_ABORT(
       Err,
       "VertAdv: VerticalTracerFluxLimiterEnabled not found in AdvectConfig");
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
// Compute VerticalVelocity and TotalVerticalVelocity from the horizontal
// velocity (NormalVelocity) and the layer thickness used for fluxes through
// edges (FluxLayerThickEdge)
void VertAdv::computeVerticalVelocity(
    const Array2DReal &NormalVelocity,    ///< [in] horizontal velocity
    const Array2DReal &FluxLayerThickEdge ///< [in] layer thickness at edges
) {

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
       "computeVerticalVelocity", {NCellsOwned},
       KOKKOS_LAMBDA(int ICell, const TeamMember &Team) {
          RealScratchArray DivHU(Team.team_scratch(0), LocNVertLayers);

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
                                         FluxLayerThickEdge(JEdge, K) *
                                         NormalVelocity(JEdge, K) * InvAreaCell;
                    }
                 }
                 for (int KVec = 0; KVec < KLen; ++KVec) {
                    const I4 K = KStart + KVec;
                    DivHU(K)   = DivHUTmp[KVec];
                 }
              });

          Team.team_barrier();

          // Set velocity through top and bottom interfaces to zero
          Kokkos::single(
              PerTeam(Team), INNER_LAMBDA() {
                 VerticalVelocity(ICell, KMin)     = 0.;
                 VerticalVelocity(ICell, KMax + 1) = 0.;
              });
          KRange = vertRange(KMin + 1, KMax);

          // Prefix sum of divergence to determine velocity through each
          // interface
          parallelScanInner(
              Team, KRange, INNER_LAMBDA(int K, Real &Accum, bool IsFinal) {
                 const I4 KRev = KMax - K;
                 Accum -= DivHU(KRev);

                 if (IsFinal) {
                    VerticalVelocity(ICell, KRev) = Accum;
                 }
              });
       },
       NVertLayers);

   // TODO: currently assuming TotalVerticalVelocity = VerticalVelocity, i.e.
   //  purely from divergence of horizontal velocity. Need to add optional
   //  corrections to transport velocity from other contributions, e.g.
   //  movement of vertical interfaces, contribution of horizontal velocity
   //  through tilted interface.
   deepCopy(TotalVerticalVelocity, VerticalVelocity);

} // end computeVerticalVelocity

//------------------------------------------------------------------------------
// Compute thickness tendency due to vertical advection
void VertAdv::computeThicknessVAdvTend(
    const Array2DReal &ThickTend ///< [inout] thickness tendency
) {

   // Return if vertical advection thickness tendency not enabled
   if (!ThickVertAdvEnabled)
      return;

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocTotVertVelocity, TotalVerticalVelocity);

   // Loop over every owned cell, pseudo thickness tendency is simply
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
    const Array2DReal &VelTend,        ///< [inout] horizontal velocity tendency
    const Array2DReal &NormalVelocity, ///< [in] horizontal velocity
    const Array2DReal &FluxLayerThickEdge ///< [in] layer thickness at edges
) {

   // Return if vertical advection velocity tendency not enabled
   if (!VelVertAdvEnabled)
      return;

   OMEGA_SCOPE(LocCOnE, Mesh->CellsOnEdge);
   OMEGA_SCOPE(MinLayerEdgeBot, VCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(MaxLayerEdgeTop, VCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocTotVertVelocity, TotalVerticalVelocity);
   OMEGA_SCOPE(EdgeMask, VCoord->EdgeMask);
   OMEGA_SCOPE(LocNVertLayersP1, NVertLayersP1);

   // Loop over every owned edge
   parallelForOuter(
       "computeVelocityVAdvTend", {NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
          const I4 Cell1 = LocCOnE(IEdge, 0);
          const I4 Cell2 = LocCOnE(IEdge, 1);
          const I4 KMin  = MinLayerEdgeBot(IEdge);
          const I4 KMax  = MaxLayerEdgeTop(IEdge);
          I4 KRange      = vertRangeChunked(KMin + 1, KMax);

          // Allocate scratch space for W times Du/Dz at vertical interfaces
          // between edges
          RealScratchArray WDuDzEdge(Team.team_scratch(0), LocNVertLayersP1);

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
                        (0.5_Real * (FluxLayerThickEdge(IEdge, K - 1) +
                                     FluxLayerThickEdge(IEdge, K)));
                 }
              });

          Team.team_barrier();

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
       },
       NVertLayersP1);

} // end computeVelocityVAdvTend

//------------------------------------------------------------------------------
// Compute tracer tendency due to vertical advection, TimeStep is only needed
// as an arugement for flux-corrected transport
void VertAdv::computeTracerVAdvTend(
    const Array3DReal &TracerTend,     ///< [inout] tracer tendencies
    const Array3DReal &Tracers,        ///< [in] tracer array
    const Array2DReal &LayerThickness, ///< [in] layer thickness
    const TimeInterval TimeStep        ///< [in] (optional) time step
) {

   computeVerticalFluxes(Tracers, LayerThickness);

   // dispatch to appropriate algorithm based on configuration settings
   switch (VertAdvChoice) {
   case VertAdvOption::Standard:
      computeStdVAdvTend(TracerTend);
      break;
   case VertAdvOption::FCT:
      R8 Dt;
      TimeStep.get(Dt, TimeUnits::Seconds);
      computeFCTVAdvTend(TracerTend, Tracers, LayerThickness, Dt);
      break;
   }

} // end computeTracerVAdvTend

//------------------------------------------------------------------------------
// Compute tracer fluxes due to vertical advection, the particular scheme used
// is chosen via configuration settings
void VertAdv::computeVerticalFluxes(
    const Array3DReal &Tracers,       ///< [in] tracer array
    const Array2DReal &LayerThickness ///< [in] layer thickness
) {

   OMEGA_SCOPE(MinLayerCell, VCoord->MinLayerCell);
   OMEGA_SCOPE(MaxLayerCell, VCoord->MaxLayerCell);
   OMEGA_SCOPE(LocTotVertVel, TotalVerticalVelocity);
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
                       const Real InvLayerThickSum =
                           1._Real / (LayerThickness(ICell, K - 1) +
                                      LayerThickness(ICell, K));
                       const Real VerticalWeightK =
                           LayerThickness(ICell, K - 1) * InvLayerThickSum;
                       const Real VerticalWeightKm1 =
                           LayerThickness(ICell, K) * InvLayerThickSum;
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
             const Real InvLayerThickSum =
                 1._Real /
                 (LayerThickness(ICell, K - 1) + LayerThickness(ICell, K));
             const Real VerticalWeightK =
                 LayerThickness(ICell, K - 1) * InvLayerThickSum;
             const Real VerticalWeightKm1 =
                 LayerThickness(ICell, K) * InvLayerThickSum;
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
    const Array3DReal &TracerTend ///< [inout] tracer tendencies
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

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
