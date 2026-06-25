#ifndef OMEGA_SURFACECOUPLING_H
#define OMEGA_SURFACECOUPLING_H
//===-- ocn/SurfaceCouling.h - surface coupling ----------------*- C++ -*-===//
//
/// \file
/// \brief Contains the coupling variables exchanged with the coupler
///
/// The SurfaceCouling class contains the variables exchanged with the coupler
/// for a sub-domain of the global horizontal mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "HorzMesh.h"
#include "TimeMgr.h"
#include "TimeStepper.h"

#include <string>

namespace OMEGA {

enum class CouplingLayout { MCT, MOAB };

// Parameters needed to initialize a SfcCoupling object. The information
// needed to initialize these parameters is provided by the coupler.
struct CouplingInitParams {
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;
   TimeInterval CouplingTimeStep;
   CouplingLayout Layout;
};

// x2o: Coupler to Ocean
class CplToOcnFields {
 public:
   // x2o fields only need to be stored on the host.
   // The SfcCoupling::applyImportFields() method will handle copying the
   // data to the device.
   HostArray1DReal SurfaceStressZonal;      ///< Foxx_taux  [N m^-2]
   HostArray1DReal SurfaceStressMeridional; ///< Foxx_tauy  [N m^-2]

   CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh);
};

// o2x: Ocean to Coupler
class OcnToCplFields {
 public:
   ///< So_t    [deg C]
   Array1DReal SurfaceTemperature;
   HostArray1DReal SurfaceTemperature_H;

   ///< So_u    [m s^-1]
   Array1DReal SurfaceVelocityZonal;
   HostArray1DReal SurfaceVelocityZonal_H;

   ///< So_v    [m s^-1]
   Array1DReal SurfaceVelocityMeridional;
   HostArray1DReal SurfaceVelocityMeridional_H;

   OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh);
};

/// A class for interfacing with the coupler

/// The SfcCoupling class provides a container for the variables exchanged
/// to (o2x) and from (x2o) the coupler. It containes methods to handle the
/// import and export of raw data from the coupler, unit conversion of said
/// data, application of that data to the model state, and accumulation of the
/// data for avergaing or sumation over multiple ocean time steps.
class SfcCoupling {

 private:
   static SfcCoupling *DefaultSfcCoupling;

   static std::map<std::string, std::unique_ptr<SfcCoupling>> AllSfcCoupling;

   // Construct a new local coupling object
   SfcCoupling(const std::string &Name_, const HorzMesh *Mesh,
               const std::map<std::string, int> &ImportIdx,
               const std::map<std::string, int> &ExportIdx,
               TimeStepper *Stepper, const TimeInterval &CouplingTimeStep,
               const CouplingLayout &Layout);

   // Forbid copy and move construction
   SfcCoupling(const SfcCoupling &) = delete;
   SfcCoupling(SfcCoupling &&)      = delete;

   // Number of ocn timesteps acccumulated over the coupling interval
   I4 NAccumSteps;

   CouplingLayout Layout; ///< Coupling layout (MCT or MOAB)

   // Map of import/export variable names to index in the raw data arrays
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;

 public:
   std::string Name;

   I4 NCellsOwned; ///< Number of cells owned by this task

   // Coupling Variable containers
   CplToOcnFields CplToOcn; ///< Coupler to Ocean (x2o)
   OcnToCplFields OcnToCpl; ///< Ocean to Coupler (o2x)

   Alarm CouplingAlarm; ///< Alarm for coupling interval

   // Methods

   /// Create a new surface coupling by calling the constructor and put it
   /// in the AllSfcCoupling map
   static SfcCoupling *create(const std::string &Name, const HorzMesh *Mesh,
                              const std::map<std::string, int> &ImportIdx,
                              const std::map<std::string, int> &ExportIdx,
                              TimeStepper *Stepper,
                              const TimeInterval &CouplingTimeStep,
                              const CouplingLayout &CouplingLayout);

   /// Initlaize SfcCoupling
   static int init(const CouplingInitParams &CouplingInitParams);

   /// Destructor - deallocates all memory and deletes an SfcCoupling
   ~SfcCoupling();

   /// Dealllocates arrays
   static void clear();

   /// Remove surface coupling object by name
   static void erase(const std::string InName);

   static SfcCoupling *getDefault();

   static SfcCoupling *get(const std::string name);
};

} // end namespace OMEGA
#endif // defined OMEGA_SURFACECOUPLING_H
