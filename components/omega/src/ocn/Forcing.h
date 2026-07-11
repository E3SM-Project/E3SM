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

#include <map>
#include <memory>
#include <string>

namespace OMEGA {

/// A class that manages surface stress forcing for ocean dynamics.
/// The Forcing class provides surface stress data and handles IO,
/// halo exchanges, and the computation of stress forcing on edges.
class Forcing {
 public:
   std::string Name; ///< Name identifier for this forcing instance

   SfcStressForcingVars SfcStressForcing; ///< Surface stress forcing variables

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

   /// Register surface stress fields with IO streams
   void registerFields(const std::string &MeshName) const;

   /// Unregister surface stress fields from IO streams
   void unregisterFields() const;

   /// Read forcing fields from input stream at startup
   void readStreamIntoArrays();

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

   static Forcing *DefaultForcing;
   static std::map<std::string, std::unique_ptr<Forcing>> AllForcing;
};

} // namespace OMEGA

#endif
