//===-- ocn/Forcing.cpp - Forcing ------------------*- C++ -*-===//
//
// The Forcing class manages external forcing from data or coupled components,
// including ocean surface stress, tracer forcing, and shortwave extinction
// coefficients.
//
//===--------------------------------------------------------------===//

#include "Forcing.h"
#include "Error.h"
#include "Field.h"
#include "FillValues.h"
#include "IOStream.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "Pacer.h"

#include <limits>
#include <vector>

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
      TracerForcing(stripDefault(Name), Mesh),
      ShortwavePenForcing(stripDefault(Name), Mesh),
      WindSpeed10mCell("WindSpeed10m" + stripDefault(Name), Mesh->NCellsSize),
      IceFractionCell("IceFraction" + stripDefault(Name), Mesh->NCellsSize),
      Mesh(Mesh), MeshHalo(MeshHalo) {
   deepCopy(WindSpeed10mCell, 0.0_Real);
   deepCopy(IceFractionCell, 0.0_Real);
}

// Destructor. Unregisters fields from IO streams.
Forcing::~Forcing() { unregisterFields(); }

// Register forcing fields with IO streams for a given mesh.
void Forcing::registerFields(const std::string &MeshName) const {
   ShortwavePenForcing.registerFields("ShortwaveExtinction", MeshName);
   ShortwavePenFieldsRegistered = true;

   if (SfcStressFieldsEnabled) {
      SfcStressForcing.registerFields(MeshName);
   }
   if (TracerForcingFieldsEnabled) {
      TracerForcing.registerFields(MeshName);
   }
   if (WindSpeed10mFieldEnabled) {
      std::vector<std::string> DimNames(1);
      DimNames[0] = MeshName == "Default" ? "NCells" : "NCells" + MeshName;
      auto WindSpeed10mField =
          Field::create(WindSpeed10mCell.label(), "10 m wind speed", "m s-1",
                        "", 0.0, std::numeric_limits<Real>::max(), 1, DimNames);
      WindSpeed10mField->setOptionalRead(true);
      FieldGroup::addFieldToGroup(WindSpeed10mCell.label(), "Forcing");
      WindSpeed10mField->attachData<Array1DReal>(WindSpeed10mCell);
   }
   if (IceFractionFieldEnabled) {
      std::vector<std::string> DimNames(1);
      DimNames[0] = MeshName == "Default" ? "NCells" : "NCells" + MeshName;
      auto IceFractionField =
          Field::create(IceFractionCell.label(), "sea ice fraction", "1", "",
                        0.0, 1.0, 1, DimNames);
      IceFractionField->setOptionalRead(true);
      FieldGroup::addFieldToGroup(IceFractionCell.label(), "Forcing");
      IceFractionField->attachData<Array1DReal>(IceFractionCell);
   }
}

// Unregister forcing fields from IO streams.
void Forcing::unregisterFields() const {
   if (ShortwavePenFieldsRegistered) {
      ShortwavePenForcing.unregisterFields();
      ShortwavePenFieldsRegistered = false;
   }

   if (SfcStressFieldsEnabled) {
      SfcStressForcing.unregisterFields();
   }
   if (TracerForcingFieldsEnabled) {
      TracerForcing.unregisterFields();
   }
   if (WindSpeed10mFieldEnabled && Field::exists(WindSpeed10mCell.label())) {
      Field::destroy(WindSpeed10mCell.label());
   }
   if (IceFractionFieldEnabled && Field::exists(IceFractionCell.label())) {
      Field::destroy(IceFractionCell.label());
   }
}

// Create and register a non-default forcing instance.
Forcing *Forcing::create(const std::string &Name, const HorzMesh *Mesh,
                         Halo *MeshHalo) {
   OMEGA_REQUIRE(
       Mesh, "Null HorzMesh pointer in Forcing::create with Name = {}", Name);
   OMEGA_REQUIRE(MeshHalo,
                 "Null Halo pointer in Forcing::create with Name = {}", Name);

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
   FieldGroup::create("ShortwaveExtinction");

   const HorzMesh *DefMesh = HorzMesh::getDefault();
   OMEGA_REQUIRE(DefMesh, "Null default HorzMesh pointer in Forcing::init");
   Halo *DefHalo = Halo::getDefault();
   OMEGA_REQUIRE(DefHalo, "Null default Halo pointer in Forcing::init");

   DefaultForcing = Forcing::create("Default", DefMesh, DefHalo);

   if (DefaultForcing == nullptr) {
      ABORT_ERROR("Forcing: failed to initialize default forcing state");
   }

   Config *OmegaConfig = Config::getOmegaConfig();
   OMEGA_REQUIRE(OmegaConfig, "Null OmegaConfig pointer in Forcing::init");
   DefaultForcing->readConfigOptions(OmegaConfig);
   DefaultForcing->registerFields(DefMesh->MeshName);
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
   if (FieldGroup::exists("ShortwaveExtinction")) {
      FieldGroup::destroy("ShortwaveExtinction");
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

   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(Err, "Forcing: Tendencies group not found in Config");

   Err +=
       TendConfig.get("SfcStressForcingTendencyEnable", SfcStressFieldsEnabled);
   CHECK_ERROR_ABORT(Err, "Forcing: SfcStressForcingTendencyEnable not found "
                          "in Tendencies config");

   bool SfcThicknessForcingEnabled = false;
   Err += TendConfig.get("SfcThicknessForcingTendencyEnable",
                         SfcThicknessForcingEnabled);
   CHECK_ERROR_ABORT(Err, "Forcing: SfcThicknessForcingTendencyEnable not "
                          "found in Tendencies config");

   bool SfcTracerForcingEnabled = false;
   Err += TendConfig.get("SfcTracerForcingTendencyEnable",
                         SfcTracerForcingEnabled);
   CHECK_ERROR_ABORT(Err, "Forcing: SfcTracerForcingTendencyEnable not found "
                          "in Tendencies config");

   bool KPPNonLocalTracerFluxEnabled = false;
   Error KPPNonLocalTracerFluxErr    = TendConfig.get(
       "KPPNonLocalTracerFluxTendencyEnable", KPPNonLocalTracerFluxEnabled);
   if (!KPPNonLocalTracerFluxErr.isSuccess()) {
      KPPNonLocalTracerFluxErr.reset();
   }

   Config KPPConfig("KPP");
   Config VertMixConfig("VertMix");
   Error KPPConfigErr = OmegaConfig->get(VertMixConfig);
   KPPConfigErr += VertMixConfig.get(KPPConfig);
   bool KPPEnabled            = false;
   bool UseLangmuirTurbulence = false;
   if (KPPConfigErr.isSuccess()) {
      KPPConfig.get("Enable", KPPEnabled);
      KPPConfig.get("UseLangmuirTurbulence", UseLangmuirTurbulence);
   }
   WindSpeed10mFieldEnabled = KPPEnabled && UseLangmuirTurbulence;
   // Ice fraction feeds both Langmuir suppression and the minimum-OSBL-
   // under-ice clamp, so it is needed whenever KPP itself is enabled.
   IceFractionFieldEnabled          = KPPEnabled;
   bool PenetratingShortwaveEnabled = false;
   Config PenSwConfig("PenetratingShortwaveTendency");
   Err += TendConfig.get(PenSwConfig);
   CHECK_ERROR_ABORT(Err, "Forcing: PenetratingShortwaveTendency group not "
                          "found in Tendencies config");
   Err += PenSwConfig.get("Enable", PenetratingShortwaveEnabled);
   CHECK_ERROR_ABORT(Err, "Forcing: PenetratingShortwaveTendency Enable not "
                          "found in PenetratingShortwaveTendency config");

   TracerForcingFieldsEnabled = SfcThicknessForcingEnabled ||
                                SfcTracerForcingEnabled ||
                                PenetratingShortwaveEnabled;
}

// Compute all forcing variables (dispatches to specific computations).
void Forcing::computeAll() const {
   exchangeHalo();
   computeSfcStressForcingOnEdge();
}

// Reset forcing arrays so omitted optional fields remain zero after read.
void Forcing::resetArrays() {
   if (SfcStressFieldsEnabled) {
      deepCopy(SfcStressForcing.NormalStressEdge, 0.0_Real);
      deepCopy(SfcStressForcing.ZonalStressCell, 0.0_Real);
      deepCopy(SfcStressForcing.MeridStressCell, 0.0_Real);
   }

   if (TracerForcingFieldsEnabled) {
      deepCopy(TracerForcing.SnowFluxCell, 0.0_Real);
      deepCopy(TracerForcing.RainFluxCell, 0.0_Real);
      deepCopy(TracerForcing.EvaporationFluxCell, 0.0_Real);
      deepCopy(TracerForcing.SeaIceFreshWaterFluxCell, 0.0_Real);
      deepCopy(TracerForcing.IceRunoffFluxCell, 0.0_Real);
      deepCopy(TracerForcing.RiverRunoffFluxCell, 0.0_Real);
      deepCopy(TracerForcing.LatentHeatFluxEvapCell, 0.0_Real);
      deepCopy(TracerForcing.SensibleHeatFluxCell, 0.0_Real);
      deepCopy(TracerForcing.LongWaveHeatFluxUpCell, 0.0_Real);
      deepCopy(TracerForcing.LongWaveHeatFluxDownCell, 0.0_Real);
      deepCopy(TracerForcing.SeaIceHeatFluxCell, 0.0_Real);
      deepCopy(TracerForcing.ShortWaveHeatFluxCell, 0.0_Real);
      deepCopy(TracerForcing.SeaIceSaltFluxCell, 0.0_Real);
      deepCopy(TracerForcing.SurfaceTracerFluxCell, FillValueReal);
   }

   if (WindSpeed10mFieldEnabled) {
      deepCopy(WindSpeed10mCell, 0.0_Real);
   }

   if (IceFractionFieldEnabled) {
      deepCopy(IceFractionCell, 0.0_Real);
   }
}

void Forcing::setSurfaceTracerFlux(const Array2DReal &Flux) {
   auto &SurfaceTracerFlux = TracerForcing.SurfaceTracerFluxCell;
   OMEGA_REQUIRE(Flux.extent(0) == SurfaceTracerFlux.extent(0),
                 "Forcing::setSurfaceTracerFlux: tracer dimension mismatch");
   OMEGA_REQUIRE(Flux.extent(1) == SurfaceTracerFlux.extent(1),
                 "Forcing::setSurfaceTracerFlux: cell dimension mismatch");
   deepCopy(SurfaceTracerFlux, Flux);

   I4 HaloErr = 0;
   for (I4 L = 0; L < SurfaceTracerFlux.extent(0); ++L) {
      auto TracerFlux = subviewUnmanaged(SurfaceTracerFlux, L, Kokkos::ALL);
      HaloErr += MeshHalo->exchangeFullArrayHalo(TracerFlux, OnCell);
   }
   OMEGA_REQUIRE(HaloErr == 0,
                 "Forcing::setSurfaceTracerFlux: halo exchange failed");
}

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

// Exchange halo for surface stress cell fields. Only needed for variables that
// need information beyond cell-centered values.
I4 Forcing::exchangeHalo() const {
   I4 Err = 0;

   if (SfcStressFieldsEnabled) {
      Err += MeshHalo->exchangeFullArrayHalo(SfcStressForcing.ZonalStressCell,
                                             OnCell);
      Err += MeshHalo->exchangeFullArrayHalo(SfcStressForcing.MeridStressCell,
                                             OnCell);
   }

   if (TracerForcingFieldsEnabled) {
      const I4 NTracers = TracerForcing.SurfaceTracerFluxCell.extent(0);
      for (I4 L = 0; L < NTracers; ++L) {
         auto TracerFlux = subviewUnmanaged(TracerForcing.SurfaceTracerFluxCell,
                                            L, Kokkos::ALL);
         Err += MeshHalo->exchangeFullArrayHalo(TracerFlux, OnCell);
      }
   }

   if (WindSpeed10mFieldEnabled) {
      Err += MeshHalo->exchangeFullArrayHalo(WindSpeed10mCell, OnCell);
   }

   if (IceFractionFieldEnabled) {
      Err += MeshHalo->exchangeFullArrayHalo(IceFractionCell, OnCell);
   }
   return Err;
}

// Read forcing fields from input stream.
// To be extended with time indexing later.
void Forcing::readStreamIntoArrays() {
   Error Err;

   std::string StreamName = "Forcing";

   resetArrays();

   // Nothing to read if no forcing tendencies or KPP wind/ice inputs are
   // enabled.
   if (!SfcStressFieldsEnabled && !TracerForcingFieldsEnabled &&
       !WindSpeed10mFieldEnabled && !IceFractionFieldEnabled) {
      return;
   }

   // Attempt to read stream; if unavailable, log and fall back to zero forcing.
   Err = IOStream::read(StreamName);
   if (Err.isFail()) {
      LOG_INFO("Forcing: Error while reading {} stream, using zero forcing",
               StreamName);
      resetArrays();
   }

   I4 HaloErr = exchangeHalo();
   if (HaloErr != 0) {
      ABORT_ERROR("Forcing: Error exchanging halo for startup forcing fields");
   }

   computeAll();
}

} // namespace OMEGA
