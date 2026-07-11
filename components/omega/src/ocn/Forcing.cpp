//===-- ocn/Forcing.cpp - Forcing ------------------*- C++ -*-===//
//
// The Forcing class manages the external forcing (from data or coupled
//  components). For now, it only includes ocean surface stress forcing
// but will include surface restoring and surface thermodynamical forcing.
// For now, it contains
// surface stress data on cells and provides methods to compute
// edge-normal stress components, manage halo exchanges, and handle IO.
//
//===--------------------------------------------------------------===//

#include "Forcing.h"
#include "Field.h"
#include "IOStream.h"
#include "Logging.h"
#include "Pacer.h"

namespace OMEGA {

// Static member initialization
Forcing *Forcing::DefaultForcing = nullptr;
std::map<std::string, std::unique_ptr<Forcing>> Forcing::AllForcing;

static std::string stripDefault(const std::string &Name) {
   return Name != "Default" ? Name : "";
}

// Constructor. Initializes surface stress forcing variables and stores
// mesh/halo.
Forcing::Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo)
    : Name(stripDefault(Name)), SfcStressForcing(stripDefault(Name), Mesh),
      Mesh(Mesh), MeshHalo(MeshHalo) {}

// Destructor. Unregisters fields from IO streams.
Forcing::~Forcing() { unregisterFields(); }

// Register surface stress fields with IO streams for a given mesh.
void Forcing::registerFields(const std::string &MeshName) const {
   SfcStressForcing.registerFields(MeshName);
}

// Unregister surface stress fields from IO streams.
void Forcing::unregisterFields() const { SfcStressForcing.unregisterFields(); }

// Create and register a non-default forcing instance.
Forcing *Forcing::create(const std::string &Name, const HorzMesh *Mesh,
                         Halo *MeshHalo) {
   if (AllForcing.find(Name) != AllForcing.end()) {
      LOG_ERROR("Attempted to create new Forcing with name {} but it already "
                "exists",
                Name);
      return nullptr;
   }

   auto *NewForcing = new Forcing(Name, Mesh, MeshHalo);
   AllForcing.emplace(Name, NewForcing);

   return NewForcing;
}

// Initialize the default forcing instance and read configuration.
// Reads the forcing fields from streams at startup (for now).
void Forcing::init() {
   if (DefaultForcing != nullptr) {
      return;
   }

   FieldGroup::create("Forcing");

   const HorzMesh *DefMesh = HorzMesh::getDefault();
   Halo *DefHalo           = Halo::getDefault();

   DefaultForcing = Forcing::create("Default", DefMesh, DefHalo);

   if (DefaultForcing == nullptr) {
      ABORT_ERROR("Forcing: failed to initialize default forcing state");
   }

   DefaultForcing->registerFields(DefMesh->MeshName);

   Config *OmegaConfig = Config::getOmegaConfig();
   DefaultForcing->readConfigOptions(OmegaConfig);
   // for now, forcing fields are read at start-up only.
   // to be extended to include switch from standalone to coupled.
   // to be moved to a Forcing->prepareForStep(SimTime) method later.
   DefaultForcing->readStreamIntoArrays();
}

// Return the default forcing instance.
Forcing *Forcing::getDefault() { return DefaultForcing; }

// Return a forcing instance by name, or null if not found.
Forcing *Forcing::get(const std::string &Name) {
   auto it = AllForcing.find(Name);
   if (it != AllForcing.end()) {
      return it->second.get();
   }

   LOG_ERROR("Forcing::get: Attempt to retrieve non-existent forcing state:");
   LOG_ERROR("{} has not been defined or has been removed", Name);
   return nullptr;
}

// Return true if a forcing instance with the given name exists.
bool Forcing::exists(const std::string &Name) {
   return AllForcing.find(Name) != AllForcing.end();
}

// Remove a forcing instance by name from the registry.
void Forcing::erase(const std::string &Name) { AllForcing.erase(Name); }

// Clear all forcing instances and destroy the field group.
void Forcing::clear() {
   AllForcing.clear();
   DefaultForcing = nullptr;
   if (FieldGroup::exists("Forcing")) {
      FieldGroup::destroy("Forcing");
   }
}

// Read SfcStress configuration and set interpolation method.
void Forcing::readConfigOptions(Config *OmegaConfig) {
   Error Err;

   Config SfcStressConfig("SfcStress");
   Err += OmegaConfig->get(SfcStressConfig);

   std::string SfcStressInterpTypeStr;
   Err += SfcStressConfig.get("InterpType", SfcStressInterpTypeStr);
   CHECK_ERROR_ABORT(Err, "Forcing: InterpType not found in SfcStressConfig");

   if (SfcStressInterpTypeStr == "Isotropic") {
      this->SfcStressForcing.InterpChoice = InterpCellToEdgeOption::Isotropic;
   } else if (SfcStressInterpTypeStr == "Anisotropic") {
      this->SfcStressForcing.InterpChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      ABORT_ERROR("Forcing: Unknown InterpType requested");
   }
}

// Compute all forcing variables (dispatches to specific computations).
void Forcing::computeAll() const { computeSfcStressForcingOnEdge(); }

// Compute edge-normal stress from cell-center zonal and meridional components.
void Forcing::computeSfcStressForcingOnEdge() const {
   OMEGA_SCOPE(LocSfcStressForcing, SfcStressForcing);

   Pacer::start("Forcing:edge1", 2);
   parallelFor(
       "Forcing:edge1", {Mesh->NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          LocSfcStressForcing.computeVarsOnEdge(IEdge);
       });
   Pacer::stop("Forcing:edge1", 2);
}

// Exchange halo for surface stress cell fields.
I4 Forcing::exchangeHalo() const {
   I4 Err = 0;

   Err += MeshHalo->exchangeFullArrayHalo(SfcStressForcing.ZonalStressCell,
                                          OnCell);
   Err += MeshHalo->exchangeFullArrayHalo(SfcStressForcing.MeridStressCell,
                                          OnCell);

   return Err;
}

// Read forcing fields from input stream.
// To be extended with time indexing later.
void Forcing::readStreamIntoArrays() {
   Error Err;

   std::string StreamName = "Forcing";

   // Attempt to read stream; if unavailable, log and fall back to zero forcing.
   Err = IOStream::read(StreamName);
   if (Err.isFail()) {
      LOG_INFO("Forcing: Error while reading {} stream, using zero forcing",
               StreamName);
      deepCopy(SfcStressForcing.ZonalStressCell, 0._Real);
      deepCopy(SfcStressForcing.MeridStressCell, 0._Real);
   }

   I4 HaloErr = exchangeHalo();
   if (HaloErr != 0) {
      ABORT_ERROR("Forcing: Error exchanging halo for startup forcing fields");
   }

   computeAll();
}

} // namespace OMEGA
