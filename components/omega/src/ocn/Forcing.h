#ifndef OMEGA_FORCING_H
#define OMEGA_FORCING_H
//===-- ocn/Forcing.h - Forcing -------------------------*- C++ -*-===//
//
/// \file
/// \brief Manages forcing for ocean dynamics
///
/// For now, the Forcing class only contains
///  surface stress forcing data and the methods to
/// compute the surface stress contribution on edges necessary
/// to compute momentum tendencies.
//
//===-------------------------------------------------------------===//

#include "Config.h"
#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "forcingVars/SfcStressForcingVars.h"
#include "forcingVars/ShortwavePenForcingVars.h"
#include "forcingVars/TracerForcingVars.h"

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

/// A class that manages forcing data for ocean dynamics.
/// The Forcing class handles field registration, IO, and halo exchanges.
class Forcing {
 public:
   std::string Name; ///< Name identifier for this forcing instance

   SfcStressForcingVars SfcStressForcing; ///< Surface stress forcing variables
   TracerForcingVars TracerForcing; ///< Tracer forcing vars (thickness and T,S)
   ShortwavePenForcingVars ShortwavePenForcing; ///< Shortwave extinction coeffs

   /// 10 m wind speed (m/s), used for the KPP Langmuir enhancement estimate.
   /// Only read/registered when KPP Langmuir turbulence is enabled.
   Array1DReal WindSpeed10mCell;

   /// Sea ice fraction (0-1), used by KPP for Langmuir suppression and the
   /// minimum boundary layer depth under ice. Only read/registered when KPP
   /// is enabled.
   Array1DReal IceFractionCell;

   ~Forcing();

   /// Initialize the default forcing instance
   static void init();

   /// Create a non-default forcing instance
   static Forcing *create(const std::string &Name, const HorzMesh *Mesh,
                          Halo *MeshHalo);

   /// Get the default forcing instance
   static Forcing *getDefault();

   /// Get forcing instance by name
   static Forcing *get(const std::string &Name);

   /// Check if forcing instance exists by name
   static bool exists(const std::string &Name);

   /// Remove forcing instance by name
   static void erase(const std::string &Name);

   /// Remove all forcing instances
   static void clear();

   /// Read and apply configuration options for surface stress
   void readConfigOptions(Config *OmegaConfig);

   /// Register forcing fields with IO streams
   void registerFields(const std::string &MeshName) const;

   /// Unregister surface stress fields from IO streams
   void unregisterFields() const;

   /// Read forcing fields from input stream at startup
   void readStreamIntoArrays();

   /// Reset all forcing arrays to zero before reading optional fields
   void resetArrays();

   /// Replace the complete current-step KPP non-local tracer surface flux.
   /// The caller must provide every non-temperature/salinity tracer and local
   /// cell entry explicitly; KPP overwrites the temperature and salinity slots.
   void setSurfaceTracerFlux(const Array2DReal &Flux);

   /// Compute all forcing variables
   void computeAll() const;

   /// Compute surface stress on edges from cell-center stress components
   void computeSfcStressForcingOnEdge() const;

   /// Exchange halo for surface stress fields
   I4 exchangeHalo() const;

 private:
   Forcing(const std::string &Name, const HorzMesh *Mesh, Halo *MeshHalo);

   Forcing(const Forcing &) = delete;
   Forcing(Forcing &&)      = delete;

   const HorzMesh *Mesh;
   Halo *MeshHalo;
   bool SfcStressFieldsEnabled               = false;
   bool TracerForcingFieldsEnabled           = false;
   bool WindSpeed10mFieldEnabled             = false;
   bool IceFractionFieldEnabled              = false;
   mutable bool ShortwavePenFieldsRegistered = false;

   static Forcing *DefaultForcing;
   static std::map<std::string, std::unique_ptr<Forcing>> AllForcing;
};

} // namespace OMEGA

#endif
