#include "SfcCoupling.h"
#include "Logging.h"

namespace OMEGA {

// create the static class member
SfcCoupling *SfcCoupling::DefaultSfcCoupling = nullptr;
std::map<std::string, std::unique_ptr<SfcCoupling>> SfcCoupling::AllSfcCoupling;

// Initalize the surface coupling. Assumes the ... have been initialized
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
   SfcCoupling::DefaultSfcCoupling =
       SfcCoupling::create("Default", DefHorzMesh, CouplingInitParams.ImportIdx,
                           CouplingInitParams.ExportIdx, DefTimeStepper,
                           CplTimeStep, CouplingInitParams.Layout);

   return Err;
}

// Construct a new surface coupling object
SfcCoupling::SfcCoupling(const std::string &Name_, const HorzMesh *Mesh,
                         const std::map<std::string, int> &ImportIdx,
                         const std::map<std::string, int> &ExportIdx,
                         TimeStepper *Stepper,
                         const TimeInterval &CouplingTimeStep,
                         const CouplingLayout &Layout)
    : Name(Name_), ImportIdx(ImportIdx), ExportIdx(ExportIdx),
      CplToOcn(Name_, Mesh), OcnToCpl(Name_, Mesh), Layout(Layout) {

   // Retrieve mesh cell count
   NCellsOwned = Mesh->NCellsOwned;
   // Retrieve import/export field counts
   NImportFields = ImportIdx.size();
   NExportFields = ExportIdx.size();

   NAccumSteps = 0;

   // Allocate variables on stack for creating the CouplingAlarm
   std::string AlarmName = "CouplingAlarm";
   Clock *StepperClock   = Stepper->getClock();
   TimeInstant StartTime = Stepper->getStartTime();

   // Create a CouplingAlarm associated with CouplingTimeStep
   if (Name_ != "Default")
      AlarmName += Name_;

   CouplingAlarm = Alarm(AlarmName, CouplingTimeStep, StartTime);
   StepperClock->attachAlarm(&CouplingAlarm);
}

// Create a new surface coupling object by calling the constructor and storing
// it in the AllSfcCoupling map
SfcCoupling *SfcCoupling::create(const std::string &Name, const HorzMesh *Mesh,
                                 const std::map<std::string, int> &ImportIdx,
                                 const std::map<std::string, int> &ExportIdx,
                                 TimeStepper *Stepper,
                                 const TimeInterval &CouplingTimeStep,
                                 const CouplingLayout &Layout) {

   // Check to see if a surface coupling of the same name already exists
   if (AllSfcCoupling.find(Name) != AllSfcCoupling.end()) {
      LOG_ERROR("Attempted to create a SfcCoupling with name {}, but a "
                "SfcCoupling with that name already exists",
                Name);
      return nullptr;
   }

   // create a new surface coupling on the heap and store it in the map of
   // unique_ptrs, which will manage its lifetime
   auto *NewSfcCoupling = new SfcCoupling(Name, Mesh, ImportIdx, ExportIdx,
                                          Stepper, CouplingTimeStep, Layout);
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
          Kokkos::LayoutStride(NExportFields, 1, NCellsOwned, NImportFields);
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

   //
   int TauxIdx = ImportIdx.at("Foxx_taux");
   int TauyIdx = ImportIdx.at("Foxx_tauy");

   // Copy Kokkos view handles
   auto CplToOcnView_        = CplToOcnView;
   auto SfcStressZonal_      = CplToOcn.SfcStressZonal;
   auto SfcStressMeridional_ = CplToOcn.SfcStressMeridional;

   /// TODO: Shouldn't be making direct calls to Kokkos here.
   ///       How often is threading used? Becuase this will be a serial loop
   ///       unless threading is used. But this has to be run on the host.
   auto Policy = Kokkos::RangePolicy<HostExecSpace, Kokkos::IndexType<int>>(
       0, NCellsOwned);
   Kokkos::parallel_for("importFromCoupler", Policy, [=](int Idx) {
      SfcStressZonal_(Idx)      = CplToOcnView_(TauxIdx, Idx);
      SfcStressMeridional_(Idx) = CplToOcnView_(TauyIdx, Idx);
   });
}

void SfcCoupling::applyImportFields(Forcing *Forcing) {

   // Copy the SfcCoupling host arrays into the Forcing device arrays.
   // Copy is only done over the owned cells, since thats all the SfcCoupling
   // data is defined over. Forcing will be responsible for halo exchanges.
   deepCopy(ownedSubView(Forcing->SfcStressForcing.ZonalStressCell),
            CplToOcn.SfcStressZonal);
   deepCopy(ownedSubView(Forcing->SfcStressForcing.MeridStressCell),
            CplToOcn.SfcStressMeridional);
};

CplToOcnFields::CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh)
    : SfcStressZonal("SfcStressZonal" + Suffix, Mesh->NCellsOwned),
      SfcStressMeridional("SfcStressMeridional" + Suffix, Mesh->NCellsOwned) {}

OcnToCplFields::OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh)
    : SfcTemperature("SfcTemperature" + Suffix, Mesh->NCellsOwned),
      SfcVelocityZonal("SfcVelocityZonal" + Suffix, Mesh->NCellsOwned),
      SfcVelocityMeridional("SfcVelocityMeridional" + Suffix,
                            Mesh->NCellsOwned) {}
} // namespace OMEGA
