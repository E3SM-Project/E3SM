#ifndef OMEGA_SURFACECOUPLING_H
#define OMEGA_SURFACECOUPLING_H
//===-- ocn/SfcCouling.h - surface coupling ----------------*- C++ -*-===//
//
/// \file
/// \brief Contains the coupling variables exchanged with the coupler
///
/// The SfcCouling class contains the variables exchanged with the coupler
/// for a sub-domain of the global horizontal mesh.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Forcing.h"
#include "HorzMesh.h"
#include "OceanState.h"
#include "TimeMgr.h"
#include "TimeStepper.h"

#include <string>

namespace OMEGA {

// Welford's online algorithm for computing a running average
KOKKOS_INLINE_FUNCTION Real updateAverage(const Real OldAvg,
                                          const Real NewValue,
                                          const I4 NAccumSteps) {
   return OldAvg + (NewValue - OldAvg) / (NAccumSteps + 1);
}

enum class CouplingLayout { MCT, MOAB };

// Parameters needed to initialize a SfcCoupling object. The information
// needed to initialize these parameters is provided by the coupler.
struct CouplingInitParams {
   int NImportFields;
   int NExportFields;
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
   HostArray1DReal SfcStressZonal; ///< Foxx_taux  [N m^-2]
   HostArray1DReal SfcStressMerid; ///< Foxx_tauy  [N m^-2]

   CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh);
};

// o2x: Ocean to Coupler
class OcnToCplFields {
 public:
   ///< So_t    [K], in-situ approx (potential temp at P=0)
   HostArray1DReal AvgSfcTemperatureH;

   /// TODO: Export practical salinity (unitless) to coupler
   ///< So_s    [g kg^-1], absolute salinity
   HostArray1DReal AvgSfcSalinityH;

   ///< So_u    [m s^-1]
   HostArray1DReal AvgSfcVelocityZonalH;

   ///< So_v    [m s^-1]
   HostArray1DReal AvgSfcVelocityMeridH;

   ///< So_ssh [m]
   /// instantaneous field, so no device mirror is needed
   HostArray1DReal InstSshCellH;

   // Accumulate one ocean timestep's contribution to the running averages
   void updateAverages(const OceanState *State, const Array3DReal &TracerArray,
                       I4 NAccumSteps, I4 NCellsOwned);

   // Copy device arrays into their host mirrros and do unit conversion.
   void copyToHost();

   // Reset all fields to 0
   void resetFields();

   OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh);

 private:
   // Device arrays are private to prevent any code from mirroring to host
   // except through copyToHost(), so that the different unit bewteen the
   // device and host mirrors of temperature and salinity are not exposed to
   // the rest of the code.
   Array1DReal AvgSfcTemperature; // [C], conservative temperature
   Array1DReal AvgSfcSalinity;    // [g kg^-1], absolute salinity
   Array1DReal AvgSfcVelocityZonal;
   Array1DReal AvgSfcVelocityMerid;

   // Scratch buffer for the in-situ Kelvin conversion in copyToHost()
   Array1DReal InSituTempScratch; // [K], in-situ approx (potential temp at P=0)
};

/// A class for interfacing with the coupler

/// The SfcCoupling class provides a container for the variables exchanged
/// to (o2x) and from (x2o) the coupler.
class SfcCoupling {

 private:
   static SfcCoupling *DefaultSfcCoupling;

   static std::map<std::string, std::unique_ptr<SfcCoupling>> AllSfcCoupling;

   // Construct a new local coupling object
   SfcCoupling(const std::string &Name_, const HorzMesh *Mesh,
               const int NImportFields_, const int NExportFields_,
               const std::map<std::string, int> &ImportIdx,
               const std::map<std::string, int> &ExportIdx,
               TimeStepper *Stepper, const TimeInterval &CouplingTimeStep,
               const CouplingLayout &Layout);

   // Forbid copy and move construction
   SfcCoupling(const SfcCoupling &) = delete;
   SfcCoupling(SfcCoupling &&)      = delete;

   // Create subview that only include the owned cells
   template <class View> auto ownedSubView(const View &V) const {
      return Kokkos::subview(V, std::make_pair(0, NCellsOwned));
   }

   // Number of ocn timesteps acccumulated over the coupling interval
   I4 NAccumSteps;

   CouplingLayout Layout; ///< Coupling layout (MCT or MOAB)

   // Map of import/export variable names to index in the raw data arrays
   std::map<std::string, int> ImportIdx;
   std::map<std::string, int> ExportIdx;

 public:
   std::string Name;

   I4 NCellsOwned; ///< Number of cells owned by this task

   // The values below will be larger than InportIdx.size() and ExportIdx.size()
   // because omega does not ingest all cpl fields (e.g. BGC, landice), yet...
   I4 NImportFields; ///< Num of fields in the x2o pointer array
   I4 NExportFields; ///< Num of fields in the o2x pointer array

   // Coupling Variable containers
   CplToOcnFields CplToOcn; ///< Coupler to Ocean (x2o)
   OcnToCplFields OcnToCpl; ///< Ocean to Coupler (o2x)

   Alarm CouplingAlarm; ///< Alarm for coupling interval

   /// View of Coupler to Ocean (x2o) raw data
   Kokkos::View<const Real **, Kokkos::LayoutStride, Kokkos::HostSpace,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
       CplToOcnView;

   /// View of Ocean to Coupler (o2x) raw data
   Kokkos::View<Real **, Kokkos::LayoutStride, Kokkos::HostSpace,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
       OcnToCplView;

   // Methods

   /// Create a new surface coupling by calling the constructor and put it
   /// in the AllSfcCoupling map
   static SfcCoupling *create(const std::string &Name, const HorzMesh *Mesh,
                              const int NImportFields, const int NExportFields,
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

   /// Get the default surface coupling object
   static SfcCoupling *getDefault();

   /// Get a surface coupling object by name
   static SfcCoupling *get(const std::string name);

   /// Getter for the number of ocean timestpes accumulated over the coupling
   /// interval
   I4 getNAccumSteps() const;

   /// Create views of the coupling data arrays
   void attachData(const Real *CplToOcnData, Real *OcnToCplData);

   /// Import data from the unmanaged view of x2o pointer into OcnToCpl object
   void importFromCoupler();

   /// Export data from OcnToCpl object into the unmanaged view of o2x pointer
   void exportToCoupler();

   /// Apply the imported data to the Forcing object
   void applyImportFields(Forcing *Forcing);

   /// Update the export fields
   void updateExportFields(const OceanState *State,
                           const Array3DReal &TracerArray);
};

} // end namespace OMEGA
#endif // defined OMEGA_SURFACECOUPLING_H
