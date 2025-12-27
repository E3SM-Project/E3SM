//===-- ocn/VertAdv.cpp - vertical advection --------------------*- C++ -*-===//
//
//===----------------------------------------------------------------------===//

#include "Tracers.h"
#include "VertAdv.h"

namespace OMEGA {

// create static class members
VertAdv *VertAdv::DefaultVertAdv = nullptr;
std::map<std::string, std::unique_ptr<VertAdv>> VertAdv::AllVertAdvs;

//------------------------------------------------------------------------------
// create the default VertAdv, requires prior initialization of default
// HorzMesh and VertCoord
void VertAdv::init() {

   auto Mesh     = HorzMesh::getDefault();
   auto VCoord   = VertCoord::getDefault();

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

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
