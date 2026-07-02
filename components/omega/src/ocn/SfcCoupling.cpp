#include "SfcCoupling.h"
#include "Logging.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Tracers.h"
#include "VertCoord.h"

namespace OMEGA {

// create the static class member
SfcCoupling *SfcCoupling::DefaultSfcCoupling = nullptr;
std::map<std::string, std::unique_ptr<SfcCoupling>> SfcCoupling::AllSfcCoupling;

// Initalize the surface coupling. Assumes the default HorzMesh and
// TimeStepper have been initialized
int SfcCoupling::init(const CouplingInitParams &CouplingInitParams) {

   int Err = 0; // default successful return code

   // Retrieve the default horizontal mesh and timestepper
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   auto *DefTimeStepper  = TimeStepper::getDefault();

   TimeInterval OcnTimeStep = DefTimeStepper->getTimeStep();
   TimeInterval CplTimeStep = CouplingInitParams.CouplingTimeStep;

   R8 CplTimeStepSeconds, OcnTimeStepSeconds;
   CplTimeStep.get(CplTimeStepSeconds, TimeUnits::Seconds);
   OcnTimeStep.get(OcnTimeStepSeconds, TimeUnits::Seconds);

   if (CplTimeStepSeconds < OcnTimeStepSeconds) {
      LOG_ERROR("Coupling interval is: {} (seconds)", CplTimeStepSeconds);
      LOG_ERROR("Ocean timestep is:    {} (seconds)", OcnTimeStepSeconds);
      ABORT_ERROR(
          "The ocean timestep cannot be longer than the coupling interval.");
   }

   if (std::fmod(CplTimeStepSeconds, OcnTimeStepSeconds) > 1e-10) {
      LOG_ERROR("Coupling interval is: {} (seconds)", CplTimeStepSeconds);
      LOG_ERROR("Ocean timestep is:    {} (seconds)", OcnTimeStepSeconds);
      ABORT_ERROR(
          "Coupling interval must be evenly divisible by the ocean timestep.");
   }

   // Create the default surface coupling object and set pointer to it
   SfcCoupling::DefaultSfcCoupling = SfcCoupling::create(
       "Default", DefHorzMesh, CouplingInitParams.NImportFields,
       CouplingInitParams.NImportFields, CouplingInitParams.ImportIdx,
       CouplingInitParams.ExportIdx, DefTimeStepper, CplTimeStep,
       CouplingInitParams.Layout);

   return Err;
}

// Construct a new surface coupling object
SfcCoupling::SfcCoupling(const std::string &Name_, const HorzMesh *Mesh,
                         const int NImportFields_, const int NExportFields_,
                         const std::map<std::string, int> &ImportIdx,
                         const std::map<std::string, int> &ExportIdx,
                         TimeStepper *Stepper,
                         const TimeInterval &CouplingTimeStep,
                         const CouplingLayout &Layout)
    : Name(Name_), NImportFields(NImportFields_), NExportFields(NExportFields_),
      ImportIdx(ImportIdx), ExportIdx(ExportIdx), CplToOcn(Name_, Mesh),
      OcnToCpl(Name_, Mesh), Layout(Layout) {

   // Retrieve mesh cell count
   NCellsOwned = Mesh->NCellsOwned;

   NAccumSteps = 0;

   // Allocate variables on stack for creating the CouplingAlarm
   std::string AlarmName = "CouplingAlarm";
   Clock *StepperClock   = Stepper->getClock();
   TimeInstant StartTime = Stepper->getStartTime();

   // Avoid alarm name collisions on the shared clock for non-default instances
   if (Name_ != "Default")
      AlarmName += Name_;

   CouplingAlarm = Alarm(AlarmName, CouplingTimeStep, StartTime);
   StepperClock->attachAlarm(&CouplingAlarm);
}

// Create a new surface coupling object by calling the constructor and storing
// it in the AllSfcCoupling map
SfcCoupling *SfcCoupling::create(
    const std::string &Name, const HorzMesh *Mesh, const int NImportFields,
    const int NExportFields, const std::map<std::string, int> &ImportIdx,
    const std::map<std::string, int> &ExportIdx, TimeStepper *Stepper,
    const TimeInterval &CouplingTimeStep, const CouplingLayout &Layout) {

   // Check to see if a surface coupling of the same name already exists
   if (AllSfcCoupling.find(Name) != AllSfcCoupling.end()) {
      LOG_ERROR("Attempted to create a SfcCoupling with name {}, but a "
                "SfcCoupling with that name already exists",
                Name);
      return nullptr;
   }

   // create a new surface coupling on the heap and store it in the map of
   // unique_ptrs, which will manage its lifetime
   auto *NewSfcCoupling =
       new SfcCoupling(Name, Mesh, NImportFields, NExportFields, ImportIdx,
                       ExportIdx, Stepper, CouplingTimeStep, Layout);
   AllSfcCoupling.emplace(Name, NewSfcCoupling);

   return NewSfcCoupling;
} // end SfcCoupling create

// Get the default surface coupling object
SfcCoupling *SfcCoupling::getDefault() {
   return SfcCoupling::DefaultSfcCoupling;
}

// Get a surface coupling object by name
SfcCoupling *SfcCoupling::get(const std::string Name) {
   // look for an instance of this name
   auto it = AllSfcCoupling.find(Name);

   // if found, return the pointer
   if (it != AllSfcCoupling.end()) {
      return it->second.get();

      // otherwise print error message and return nullptr
   } else {
      LOG_ERROR("SfcCoupling::get: Attempted to retrieve non-existent "
                "surface coupling object:");
      LOG_ERROR("{} has not been defined or has been removed", Name);
      return nullptr;
   }
}

// Destructor
SfcCoupling::~SfcCoupling() {}

// Remove surface coupling object by name
void SfcCoupling::erase(const std::string Name) { AllSfcCoupling.erase(Name); }

// Remove all surface coupling objects
void SfcCoupling::clear() {
   AllSfcCoupling.clear();
   DefaultSfcCoupling = nullptr; // prevent dangling pointer
}

// Getter for private member NAccumSteps
I4 SfcCoupling::getNAccumSteps() const { return NAccumSteps; }

// Create views of the raw coupling data arrays
void SfcCoupling::attachData(const Real *CplToOcnData, Real *OcnToCplData) {

   // Kokkos::LayoutStride index math uses a runtime stride value, rather than
   // a compile-time-optimized stride value. Can switch to ifdefs if this
   // becomes a performance bottleneck
   Kokkos::LayoutStride CplToOcnLayout, OcnToCplLayout;

   if (Layout == CouplingLayout::MCT) {
      /// MCT layout: (NCellsOwned, NImportFields) - field idx strides faster
      CplToOcnLayout =
          Kokkos::LayoutStride(NImportFields, 1, NCellsOwned, NImportFields);
      OcnToCplLayout =
          Kokkos::LayoutStride(NExportFields, 1, NCellsOwned, NExportFields);
   } else if (Layout == CouplingLayout::MOAB) {
      /// MOAB layout: (NImportFields, NCellsOwned) - cell idx strides faster
      CplToOcnLayout =
          Kokkos::LayoutStride(NImportFields, NCellsOwned, NCellsOwned, 1);
      OcnToCplLayout =
          Kokkos::LayoutStride(NExportFields, NCellsOwned, NCellsOwned, 1);
   } else {
      ABORT_ERROR("SfcCoupling::attachData: Unknown coupling layout");
   }

   CplToOcnView = decltype(CplToOcnView)(CplToOcnData, CplToOcnLayout);
   OcnToCplView = decltype(OcnToCplView)(OcnToCplData, OcnToCplLayout);
}

void SfcCoupling::importFromCoupler() {

   if (CplToOcnView.data() == nullptr) {
      ABORT_ERROR(
          "CplToOcnView is not attached to data. The SfcCoupling::attachData "
          "method must be called before importing data from the coupler.");
   }

   // Get import field indices for surface stress components
   int TauxIdx = ImportIdx.at("Foxx_taux");
   int TauyIdx = ImportIdx.at("Foxx_tauy");

   // Copy Kokkos view handles
   auto CplToOcnView_   = CplToOcnView;
   auto SfcStressZonal_ = CplToOcn.SfcStressZonal;
   auto SfcStressMerid_ = CplToOcn.SfcStressMerid;

   /// TODO: Shouldn't be making direct calls to Kokkos here.
   ///       How often is threading used? Becuase this will be a serial loop
   ///       unless threading is used. But this has to be run on the host.
   auto Policy = Kokkos::RangePolicy<HostExecSpace, Kokkos::IndexType<int>>(
       0, NCellsOwned);
   Kokkos::parallel_for("importFromCoupler", Policy, [=](int Idx) {
      SfcStressZonal_(Idx) = CplToOcnView_(TauxIdx, Idx);
      SfcStressMerid_(Idx) = CplToOcnView_(TauyIdx, Idx);
   });
}

void SfcCoupling::exportToCoupler() {

   if (OcnToCplView.data() == nullptr) {
      ABORT_ERROR(
          "OcnToCplView is not attached to data. The SfcCoupling::attachData "
          "method must be called before exporting data to the coupler.");
   }

   // Copy the OcnToCpl fields to their host mirrors
   OcnToCpl.copyToHost();

   int TempIdx  = ExportIdx.at("So_t");
   int SalinIdx = ExportIdx.at("So_s");
   int VelUIdx  = ExportIdx.at("So_u");
   int VelVIdx  = ExportIdx.at("So_v");
   int SshIdx   = ExportIdx.at("So_ssh");

   // Copy Kokkos view handles
   auto OcnToCplView_        = OcnToCplView;
   auto AvgSfcTemperature_   = OcnToCpl.AvgSfcTemperatureH;
   auto AvgSfcSalinity_      = OcnToCpl.AvgSfcSalinityH;
   auto AvgSfcVelocityZonal_ = OcnToCpl.AvgSfcVelocityZonalH;
   auto AvgSfcVelocityMerid_ = OcnToCpl.AvgSfcVelocityMeridH;
   auto InstSshCellH_        = OcnToCpl.InstSshCellH;

   /// TODO: Shouldn't be making direct calls to Kokkos here.
   auto Policy = Kokkos::RangePolicy<HostExecSpace, Kokkos::IndexType<int>>(
       0, NCellsOwned);
   Kokkos::parallel_for("exportToCoupler", Policy, [=](int Idx) {
      OcnToCplView_(TempIdx, Idx)  = AvgSfcTemperature_(Idx);
      OcnToCplView_(SalinIdx, Idx) = AvgSfcSalinity_(Idx);
      OcnToCplView_(VelUIdx, Idx)  = AvgSfcVelocityZonal_(Idx);
      OcnToCplView_(VelVIdx, Idx)  = AvgSfcVelocityMerid_(Idx);
      OcnToCplView_(SshIdx, Idx)   = InstSshCellH_(Idx);
   });

   OcnToCpl.resetFields(); // Reset fields to 0 for the next coupling interval
   NAccumSteps = 0;        // Reset step counter for the next coupling interval
}
void SfcCoupling::applyImportFields(Forcing *Forcing) {

   // Copy the SfcCoupling host arrays into the Forcing device arrays.
   // Copy is only done over the owned cells, since thats all the SfcCoupling
   // data is defined over. Forcing will be responsible for halo exchanges.
   deepCopy(ownedSubView(Forcing->SfcStressForcing.ZonalStressCell),
            CplToOcn.SfcStressZonal);
   deepCopy(ownedSubView(Forcing->SfcStressForcing.MeridStressCell),
            CplToOcn.SfcStressMerid);
};

void SfcCoupling::updateExportFields(const OceanState *State,
                                     const Array3DReal &TracerArray) {

   I4 TemperatureIdx, SalinityIdx;
   Tracers::getIndex(TemperatureIdx, "Temperature");
   Tracers::getIndex(SalinityIdx, "Salinity");

   auto Temperature =
       Kokkos::subview(TracerArray, TemperatureIdx, Kokkos::ALL, Kokkos::ALL);
   auto Salinity =
       Kokkos::subview(TracerArray, SalinityIdx, Kokkos::ALL, Kokkos::ALL);

   VertCoord *DefVertCoord = VertCoord::getDefault();
   OMEGA_SCOPE(LocMinLayerCell, DefVertCoord->MinLayerCell);

   OMEGA_SCOPE(LocNAccumSteps, NAccumSteps);
   OMEGA_SCOPE(LocAvgSfcSalinity, OcnToCpl.AvgSfcSalinity);
   OMEGA_SCOPE(LocAvgSfcTemp, OcnToCpl.AvgSfcTemperature);
   OMEGA_SCOPE(LocAvgSfcVelZonal, OcnToCpl.AvgSfcVelocityZonal);
   OMEGA_SCOPE(LocAvgSfcVelMerid, OcnToCpl.AvgSfcVelocityMerid);

   // TODO: Implement vector reconsturction for velocity field.
   constexpr Real ConstSfcVelocity = 1e-4;

   parallelFor(
       {NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          const int KSfc = LocMinLayerCell(ICell);

          // Update the averaged fields using Welford's online algorithm
          LocAvgSfcTemp(ICell) = updateAverage(
              LocAvgSfcTemp(ICell), Temperature(ICell, KSfc), LocNAccumSteps);

          LocAvgSfcSalinity(ICell) = updateAverage(
              LocAvgSfcSalinity(ICell), Salinity(ICell, KSfc), LocNAccumSteps);

          LocAvgSfcVelZonal(ICell) = updateAverage(
              LocAvgSfcVelZonal(ICell), ConstSfcVelocity, LocNAccumSteps);

          LocAvgSfcVelMerid(ICell) = updateAverage(
              LocAvgSfcVelMerid(ICell), ConstSfcVelocity, LocNAccumSteps);
       });

   NAccumSteps++;
}

CplToOcnFields::CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh)
    : SfcStressZonal("SfcStressZonal" + Suffix, Mesh->NCellsOwned),
      SfcStressMerid("SfcStressMeridional" + Suffix, Mesh->NCellsOwned) {}

OcnToCplFields::OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh)
    : AvgSfcTemperature("AvgSfcTemperature" + Suffix, Mesh->NCellsOwned),
      AvgSfcSalinity("AvgSfcSalinity" + Suffix, Mesh->NCellsOwned),
      AvgSfcVelocityZonal("AvgSfcVelocityZonal" + Suffix, Mesh->NCellsOwned),
      AvgSfcVelocityMerid("AvgSfcVelocityMeridional" + Suffix,
                          Mesh->NCellsOwned),
      InstSshCellH("InstSshCellH" + Suffix, Mesh->NCellsOwned) {

   // Kokkok views created with a label are zero-initialized by default.
   // We reset the feilds here anyway to be explicit about the fact that the
   // OcnToCpl fields need to being a coupling interval with all zeros.
   resetFields();

   AvgSfcTemperatureH   = createHostMirrorCopy(AvgSfcTemperature);
   AvgSfcSalinityH      = createHostMirrorCopy(AvgSfcSalinity);
   AvgSfcVelocityZonalH = createHostMirrorCopy(AvgSfcVelocityZonal);
   AvgSfcVelocityMeridH = createHostMirrorCopy(AvgSfcVelocityMerid);
}

void OcnToCplFields::copyToHost() {

   deepCopy(AvgSfcTemperatureH, AvgSfcTemperature);
   deepCopy(AvgSfcSalinityH, AvgSfcSalinity);
   deepCopy(AvgSfcVelocityZonalH, AvgSfcVelocityZonal);
   deepCopy(AvgSfcVelocityMeridH, AvgSfcVelocityMerid);

   // SSH is an instantaneous field, so we don't bother with a device mirror of
   // our own. Instead, copy from the VertCoord, which owns SSH, host array.
   VertCoord *DefVertCoord = VertCoord::getDefault();

   auto SSHCellOwned = Kokkos::subview(
       DefVertCoord->SshCell, std::make_pair(0, (int)InstSshCellH.extent(0)));

   deepCopy(InstSshCellH, SSHCellOwned);
}

// OcnToCpl fields need to being a coupling interval with all values set to 0.
void OcnToCplFields::resetFields() {
   deepCopy(AvgSfcTemperature, 0.0);
   deepCopy(AvgSfcSalinity, 0.0);
   deepCopy(AvgSfcVelocityZonal, 0.0);
   deepCopy(AvgSfcVelocityMerid, 0.0);
}
} // namespace OMEGA
