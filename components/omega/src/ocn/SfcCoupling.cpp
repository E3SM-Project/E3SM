#include "SfcCoupling.h"
#include "AuxiliaryState.h"
#include "Eos.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzOperators.h"
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
   OMEGA_REQUIRE(DefHorzMesh,
                 "Null default HorzMesh pointer in SfcCoupling::init");
   auto *DefTimeStepper = TimeStepper::getDefault();
   OMEGA_REQUIRE(DefTimeStepper,
                 "Null default TimeStepper pointer in SfcCoupling::init");

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
       CouplingInitParams.NExportFields, CouplingInitParams.ImportIdxMap,
       CouplingInitParams.ExportIdxMap, DefTimeStepper, CplTimeStep,
       CouplingInitParams.Layout);

   return Err;
}

// Construct a new surface coupling object
SfcCoupling::SfcCoupling(const std::string &Name_, const HorzMesh *Mesh,
                         const int NImportFields_, const int NExportFields_,
                         const std::map<std::string, int> &ImportIdxMap,
                         const std::map<std::string, int> &ExportIdxMap,
                         TimeStepper *Stepper,
                         const TimeInterval &CouplingTimeStep,
                         const CouplingLayout &Layout)
    : Name(Name_), NImportFields(NImportFields_), NExportFields(NExportFields_),
      ImportIdxMap(ImportIdxMap), ExportIdxMap(ExportIdxMap),
      CplToOcn(Name_, Mesh), OcnToCpl(Name_, Mesh), Layout(Layout) {

   // Retrieve mesh cell/edge count
   NCellsOwned = Mesh->NCellsOwned;
   NEdgesAll   = Mesh->NEdgesAll;

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
    const int NExportFields, const std::map<std::string, int> &ImportIdxMap,
    const std::map<std::string, int> &ExportIdxMap, TimeStepper *Stepper,
    const TimeInterval &CouplingTimeStep, const CouplingLayout &Layout) {

   OMEGA_REQUIRE(Mesh,
                 "Null HorzMesh pointer in SfcCoupling::create with Name = {}",
                 Name);
   OMEGA_REQUIRE(
       Stepper,
       "Null TimeStepper pointer in SfcCoupling::create with Name = {}", Name);

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
       new SfcCoupling(Name, Mesh, NImportFields, NExportFields, ImportIdxMap,
                       ExportIdxMap, Stepper, CouplingTimeStep, Layout);
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

// Getter for private member CouplingAlarm
Alarm *SfcCoupling::getCouplingAlarm() { return &CouplingAlarm; }

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

   // Get import field indices
   int TauxIdx   = ImportIdxMap.at("Foxx_taux");
   int TauyIdx   = ImportIdxMap.at("Foxx_tauy");
   int SwnetIdx  = ImportIdxMap.at("Foxx_swnet");
   int SenIdx    = ImportIdxMap.at("Foxx_sen");
   int LatIdx    = ImportIdxMap.at("Foxx_lat");
   int LwupIdx   = ImportIdxMap.at("Foxx_lwup");
   int LwdnIdx   = ImportIdxMap.at("Faxa_lwdn");
   int SaltIdx   = ImportIdxMap.at("Fioi_salt");
   int MelthIdx  = ImportIdxMap.at("Fioi_melth");
   int MeltwIdx  = ImportIdxMap.at("Fioi_meltw");
   int SnowIdx   = ImportIdxMap.at("Faxa_snow");
   int RainIdx   = ImportIdxMap.at("Faxa_rain");
   int EvapIdx   = ImportIdxMap.at("Foxx_evap");
   int RoflIdx   = ImportIdxMap.at("Foxx_rofl");
   int RofiIdx   = ImportIdxMap.at("Foxx_rofi");
   int BPressIdx = ImportIdxMap.at("Si_bpress");
   int PslvIdx   = ImportIdxMap.at("Sa_pslv");

   // Copy Kokkos view handles
   auto CplToOcnView_         = CplToOcnView;
   auto SfcStressZonal_       = CplToOcn.SfcStressZonalH;
   auto SfcStressMerid_       = CplToOcn.SfcStressMeridH;
   auto SnowFlux_             = CplToOcn.SnowFluxH;
   auto RainFlux_             = CplToOcn.RainFluxH;
   auto EvaporationFlux_      = CplToOcn.EvaporationFluxH;
   auto SeaIceFreshWaterFlux_ = CplToOcn.SeaIceFreshWaterFluxH;
   auto IceRunoffFlux_        = CplToOcn.IceRunoffFluxH;
   auto RiverRunoffFlux_      = CplToOcn.RiverRunoffFluxH;
   auto LatentHeatFluxEvap_   = CplToOcn.LatentHeatFluxEvapH;
   auto SensibleHeatFlux_     = CplToOcn.SensibleHeatFluxH;
   auto LongWaveHeatFluxUp_   = CplToOcn.LongWaveHeatFluxUpH;
   auto LongWaveHeatFluxDown_ = CplToOcn.LongWaveHeatFluxDownH;
   auto SeaIceHeatFlux_       = CplToOcn.SeaIceHeatFluxH;
   auto ShortWaveHeatFlux_    = CplToOcn.ShortWaveHeatFluxH;
   auto SeaIceSaltFlux_       = CplToOcn.SeaIceSaltFluxH;
   auto SurfacePressure_      = CplToOcn.SurfacePressureH;

   /// TODO: Shouldn't be making direct calls to Kokkos here.
   ///       How often is threading used? Becuase this will be a serial loop
   ///       unless threading is used. But this has to be run on the host.
   auto Policy = Kokkos::RangePolicy<HostExecSpace, Kokkos::IndexType<int>>(
       0, NCellsOwned);
   Kokkos::parallel_for("importFromCoupler", Policy, [=](int Idx) {
      SfcStressZonal_(Idx)       = CplToOcnView_(TauxIdx, Idx);
      SfcStressMerid_(Idx)       = CplToOcnView_(TauyIdx, Idx);
      SnowFlux_(Idx)             = CplToOcnView_(SnowIdx, Idx);
      RainFlux_(Idx)             = CplToOcnView_(RainIdx, Idx);
      EvaporationFlux_(Idx)      = CplToOcnView_(EvapIdx, Idx);
      SeaIceFreshWaterFlux_(Idx) = CplToOcnView_(MeltwIdx, Idx);
      IceRunoffFlux_(Idx)        = CplToOcnView_(RofiIdx, Idx);
      RiverRunoffFlux_(Idx)      = CplToOcnView_(RoflIdx, Idx);
      LatentHeatFluxEvap_(Idx)   = CplToOcnView_(LatIdx, Idx);
      SensibleHeatFlux_(Idx)     = CplToOcnView_(SenIdx, Idx);
      LongWaveHeatFluxUp_(Idx)   = CplToOcnView_(LwupIdx, Idx);
      LongWaveHeatFluxDown_(Idx) = CplToOcnView_(LwdnIdx, Idx);
      SeaIceHeatFlux_(Idx)       = CplToOcnView_(MelthIdx, Idx);
      ShortWaveHeatFlux_(Idx)    = CplToOcnView_(SwnetIdx, Idx);
      SeaIceSaltFlux_(Idx)       = CplToOcnView_(SaltIdx, Idx);

      // Limit sea ice basal pressure (Pa) to 5 m of seawater equivalent
      const Real SeaIcePressure =
          Kokkos::min(CplToOcnView_(BPressIdx, Idx), MaxSeaIcePressure);
      // Compute the relative sea level pressure (Pa)
      const Real SeaLevelPressure = CplToOcnView_(PslvIdx, Idx) - AtmRefP;
      // Compute the relative surface pressure (Pa);
      SurfacePressure_(Idx) = SeaIcePressure + SeaLevelPressure;
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

   int TempIdx        = ExportIdxMap.at("So_t");
   int SalinIdx       = ExportIdxMap.at("So_s");
   int VelUIdx        = ExportIdxMap.at("So_u");
   int VelVIdx        = ExportIdxMap.at("So_v");
   int SshIdx         = ExportIdxMap.at("So_ssh");
   int DhdxIdx        = ExportIdxMap.at("So_dhdx");
   int DhdyIdx        = ExportIdxMap.at("So_dhdy");
   int ExtraEnergyIdx = ExportIdxMap.at("Faoo_h2otemp");

   // Copy Kokkos view handles
   auto OcnToCplView_        = OcnToCplView;
   auto AvgSfcTemperature_   = OcnToCpl.AvgSfcTemperatureH;
   auto AvgSfcSalinity_      = OcnToCpl.AvgSfcSalinityH;
   auto AvgSfcVelocityZonal_ = OcnToCpl.AvgSfcVelocityZonalH;
   auto AvgSfcVelocityMerid_ = OcnToCpl.AvgSfcVelocityMeridH;
   auto AvgSfcSshGradZonal_  = OcnToCpl.AvgSfcSshGradZonalH;
   auto AvgSfcSshGradMerid_  = OcnToCpl.AvgSfcSshGradMeridH;
   auto InstSshCellH_        = OcnToCpl.InstSshCellH;
   auto AvgExtraEnergy_      = OcnToCpl.AvgExtraEnergyH;

   // Initalize all o2x fields to 0.0 for next coupling interval
   deepCopy(OcnToCplView_, 0.0_Real);

   /// TODO: Shouldn't be making direct calls to Kokkos here.
   auto Policy = Kokkos::RangePolicy<HostExecSpace, Kokkos::IndexType<int>>(
       0, NCellsOwned);
   Kokkos::parallel_for("exportToCoupler", Policy, [=](int Idx) {
      OcnToCplView_(TempIdx, Idx)        = AvgSfcTemperature_(Idx);
      OcnToCplView_(SalinIdx, Idx)       = AvgSfcSalinity_(Idx);
      OcnToCplView_(VelUIdx, Idx)        = AvgSfcVelocityZonal_(Idx);
      OcnToCplView_(VelVIdx, Idx)        = AvgSfcVelocityMerid_(Idx);
      OcnToCplView_(DhdxIdx, Idx)        = AvgSfcSshGradZonal_(Idx);
      OcnToCplView_(DhdyIdx, Idx)        = AvgSfcSshGradMerid_(Idx);
      OcnToCplView_(SshIdx, Idx)         = InstSshCellH_(Idx);
      OcnToCplView_(ExtraEnergyIdx, Idx) = AvgExtraEnergy_(Idx);
   });

   OcnToCpl.resetFields(); // Reset fields to 0 for the next coupling interval
   NAccumSteps = 0;        // Reset step counter for the next coupling interval
}
void SfcCoupling::applyImportFields(Forcing *Forcing, VertCoord *VertCoord) {

   // Copy the SfcCoupling host arrays into the Forcing/VertCoord device arrays
   // Copy is only done over the owned cells, since thats all the SfcCoupling
   // data is defined over. The halos are exchanged below, once all the owned
   // cells have been filled.
   deepCopy(ownedSubView(Forcing->SfcStressForcing.ZonalStressCell),
            CplToOcn.SfcStressZonalH);
   deepCopy(ownedSubView(Forcing->SfcStressForcing.MeridStressCell),
            CplToOcn.SfcStressMeridH);

   deepCopy(ownedSubView(Forcing->TracerForcing.SnowFluxCell),
            CplToOcn.SnowFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.RainFluxCell),
            CplToOcn.RainFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.EvaporationFluxCell),
            CplToOcn.EvaporationFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.SeaIceFreshWaterFluxCell),
            CplToOcn.SeaIceFreshWaterFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.IceRunoffFluxCell),
            CplToOcn.IceRunoffFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.RiverRunoffFluxCell),
            CplToOcn.RiverRunoffFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.LatentHeatFluxEvapCell),
            CplToOcn.LatentHeatFluxEvapH);
   deepCopy(ownedSubView(Forcing->TracerForcing.SensibleHeatFluxCell),
            CplToOcn.SensibleHeatFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.LongWaveHeatFluxUpCell),
            CplToOcn.LongWaveHeatFluxUpH);
   deepCopy(ownedSubView(Forcing->TracerForcing.LongWaveHeatFluxDownCell),
            CplToOcn.LongWaveHeatFluxDownH);
   deepCopy(ownedSubView(Forcing->TracerForcing.SeaIceHeatFluxCell),
            CplToOcn.SeaIceHeatFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.ShortWaveHeatFluxCell),
            CplToOcn.ShortWaveHeatFluxH);
   deepCopy(ownedSubView(Forcing->TracerForcing.SeaIceSaltFluxCell),
            CplToOcn.SeaIceSaltFluxH);

   deepCopy(ownedSubView(VertCoord->SurfacePressure),
            CplToOcn.SurfacePressureH);

   // Neither Forcing nor VertCoord re-exchange these arrays during the time
   // step loop, so their halos are exchanged here, immediately after the
   // owned cells have been filled from the coupler.
   Halo *MeshHalo = Halo::getDefault();

   VertCoord->updateSurfacePressure(MeshHalo);

   I4 HaloErr = Forcing->exchangeHalo();
   if (HaloErr != 0) {
      ABORT_ERROR("Error updating forcing halos after coupler import");
   }
};

void SfcCoupling::updateExportFields(const OceanState *State,
                                     const Array3DReal &TracerArray) {

   // SshCell is only written by computeMomVertAux, which the time steppers
   // last call on a stage state, so it is stale with respect to the state
   // being exported. Recompute it here, immediately before it is used.
   // TODO: once SSH is coupled directly, instead of passing the averaged SSH
   //       gradient, this will only need to be done once per coupling interval
   AuxiliaryState::getDefault()->computeMomVertAux(State, TracerArray, 0);

   const auto &ExtraEnergy =
       Forcing::getDefault()->TracerForcing.ExtraEnergyCell;
   OcnToCpl.updateFields(State, TracerArray, ExtraEnergy, NAccumSteps,
                         NCellsOwned, NEdgesAll);

   NAccumSteps++;
}

CplToOcnFields::CplToOcnFields(const std::string &Suffix, const HorzMesh *Mesh)
    : SfcStressZonalH("SfcStressZonal" + Suffix, Mesh->NCellsOwned),
      SfcStressMeridH("SfcStressMeridional" + Suffix, Mesh->NCellsOwned),
      SnowFluxH("SnowFlux" + Suffix, Mesh->NCellsOwned),
      RainFluxH("RainFlux" + Suffix, Mesh->NCellsOwned),
      EvaporationFluxH("EvaporationFlux" + Suffix, Mesh->NCellsOwned),
      SeaIceFreshWaterFluxH("SeaIceFreshWaterFlux" + Suffix, Mesh->NCellsOwned),
      IceRunoffFluxH("IceRunoffFlux" + Suffix, Mesh->NCellsOwned),
      RiverRunoffFluxH("RiverRunoffFlux" + Suffix, Mesh->NCellsOwned),
      LatentHeatFluxEvapH("LatentHeatFluxEvap" + Suffix, Mesh->NCellsOwned),
      SensibleHeatFluxH("SensibleHeatFlux" + Suffix, Mesh->NCellsOwned),
      LongWaveHeatFluxUpH("LongWaveHeatFluxUp" + Suffix, Mesh->NCellsOwned),
      LongWaveHeatFluxDownH("LongWaveHeatFluxDown" + Suffix, Mesh->NCellsOwned),
      SeaIceHeatFluxH("SeaIceHeatFlux" + Suffix, Mesh->NCellsOwned),
      ShortWaveHeatFluxH("ShortWaveHeatFlux" + Suffix, Mesh->NCellsOwned),
      SeaIceSaltFluxH("SeaIceSaltFlux" + Suffix, Mesh->NCellsOwned),
      SurfacePressureH("SurfacePressure" + Suffix, Mesh->NCellsOwned) {}

OcnToCplFields::OcnToCplFields(const std::string &Suffix, const HorzMesh *Mesh)
    : AvgSfcTemperature("AvgSfcTemperature" + Suffix, Mesh->NCellsOwned),
      AvgSfcSalinity("AvgSfcSalinity" + Suffix, Mesh->NCellsOwned),
      AvgSfcNormalVelocity("AvgSfcNormalVelocity" + Suffix, Mesh->NEdgesSize),
      AvgSfcVelocityZonalH("AvgSfcVelocityZonal" + Suffix, Mesh->NCellsOwned),
      AvgSfcVelocityMeridH("AvgSfcVelocityMeridional" + Suffix,
                           Mesh->NCellsOwned),
      AvgSfcSshGrad("AvgSfcSshGrad" + Suffix, Mesh->NEdgesSize),
      AvgSfcSshGradZonalH("AvgSfcSshGradZonal" + Suffix, Mesh->NCellsOwned),
      AvgSfcSshGradMeridH("AvgSfcSshGradMeridional" + Suffix,
                          Mesh->NCellsOwned),
      AvgExtraEnergy("AvgExtraEnergy" + Suffix, Mesh->NCellsOwned),
      InstSshCellH("InstSshCellH" + Suffix, Mesh->NCellsOwned),
      InSituTempScratch("InSituTempScratch" + Suffix, Mesh->NCellsOwned),
      PracSalinityScratch("PracSalinityScratch" + Suffix, Mesh->NCellsOwned),
      ReconZonalScratch("ReconZonalScratch" + Suffix, Mesh->NCellsOwned),
      ReconMeridScratch("ReconMeridScratch" + Suffix, Mesh->NCellsOwned) {

   // Kokkok views created with a label are zero-initialized by default.
   // We reset the fields here anyway to be explicit about the fact that the
   // OcnToCpl fields need to begin a coupling interval with all zeros.
   resetFields();

   AvgSfcTemperatureH = createHostMirrorCopy(AvgSfcTemperature);
   AvgSfcSalinityH    = createHostMirrorCopy(AvgSfcSalinity);
   AvgExtraEnergyH    = createHostMirrorCopy(AvgExtraEnergy);
}

void OcnToCplFields::updateFields(const OceanState *State,
                                  const Array3DReal &TracerArray,
                                  const Array1DReal &ExtraEnergy,
                                  const I4 NAccumSteps, const I4 NCellsOwned,
                                  const I4 NEdgesAll) {

   I4 TemperatureIdx, SalinityIdx;
   Tracers::getIndex(TemperatureIdx, "Temperature");
   Tracers::getIndex(SalinityIdx, "Salinity");

   // get normal velocity at current time level
   auto NormalVel = State->getNormalVelocity(0);
   auto Temperature =
       Kokkos::subview(TracerArray, TemperatureIdx, Kokkos::ALL, Kokkos::ALL);
   auto Salinity =
       Kokkos::subview(TracerArray, SalinityIdx, Kokkos::ALL, Kokkos::ALL);

   HorzMesh *DefHorzMesh   = HorzMesh::getDefault();
   VertCoord *DefVertCoord = VertCoord::getDefault();

   OMEGA_SCOPE(LocMinLayerCell, DefVertCoord->MinLayerCell);
   OMEGA_SCOPE(LocAvgSfcSalinity, AvgSfcSalinity);
   OMEGA_SCOPE(LocAvgSfcTemp, AvgSfcTemperature);
   OMEGA_SCOPE(LocAvgExtraEnergy, AvgExtraEnergy);

   parallelFor(
       {NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          const int KSfc = LocMinLayerCell(ICell);

          // Update the averaged fields using Welford's online algorithm
          LocAvgSfcTemp(ICell) = updateAverage(
              LocAvgSfcTemp(ICell), Temperature(ICell, KSfc), NAccumSteps);

          LocAvgSfcSalinity(ICell) = updateAverage(
              LocAvgSfcSalinity(ICell), Salinity(ICell, KSfc), NAccumSteps);

          LocAvgExtraEnergy(ICell) = updateAverage(
              LocAvgExtraEnergy(ICell), ExtraEnergy(ICell), NAccumSteps);
       });

   OMEGA_SCOPE(LocCellsOnEdge, DefHorzMesh->CellsOnEdge);
   OMEGA_SCOPE(LocDcEdge, DefHorzMesh->DcEdge);
   OMEGA_SCOPE(LocSshCell, DefVertCoord->SshCell);
   OMEGA_SCOPE(LocSurfacePressure, DefVertCoord->SurfacePressure);
   OMEGA_SCOPE(LocMinLayerEdgeBot, DefVertCoord->MinLayerEdgeBot);
   OMEGA_SCOPE(LocMaxLayerEdgeTop, DefVertCoord->MaxLayerEdgeTop);
   OMEGA_SCOPE(LocAvgSfcSshGrad, AvgSfcSshGrad);
   OMEGA_SCOPE(LocAvgSfcNormalVel, AvgSfcNormalVelocity);

   parallelFor(
       {NEdgesAll}, KOKKOS_LAMBDA(int IEdge) {
          const int KMin = LocMinLayerEdgeBot(IEdge);
          const int KMax = LocMaxLayerEdgeTop(IEdge);

          // exclude outer halo edges whose layer range is invalid
          if (KMin <= KMax) {
             LocAvgSfcNormalVel(IEdge) =
                 updateAverage(LocAvgSfcNormalVel(IEdge),
                               NormalVel(IEdge, KMin), NAccumSteps);

             const int ICell0 = LocCellsOnEdge(IEdge, 0);
             const int ICell1 = LocCellsOnEdge(IEdge, 1);

             const Real SshCell0 = pressureAdjustedSsh(
                 LocSshCell(ICell0), LocSurfacePressure(ICell0));
             const Real SshCell1 = pressureAdjustedSsh(
                 LocSshCell(ICell1), LocSurfacePressure(ICell1));

             const Real SshGrad = (SshCell1 - SshCell0) / LocDcEdge(IEdge);

             LocAvgSfcSshGrad(IEdge) =
                 updateAverage(LocAvgSfcSshGrad(IEdge), SshGrad, NAccumSteps);
          }
       });
}

// Conversion from conservative [C] to in-situ [K] temperature is done on
// the device, as well as conversion from practical to absolute salinity.
// This allows the expensive unit conversions to be called on device, once
// per coupling interval. Conversion is written to scratch buffers, to
// keep units consistent in time. Special care is paid to guard external
// access to the device arrays, so that inconsistent units between the device
// and host mirrors are not exposed to the rest of the code.
void OcnToCplFields::copyToHost() {

   // Convert averaged Conservative Temp to in-situ (approx by potential
   // temp at P=0) Kelvin, once per coupling interval, on device. Written
   // into a scratch buffer so AvgSfcTemperature always stays deg C.
   // A local Teos10Eos is constructed here rather than reusing Eos's
   // instance: calcPtFromCt() needs no config-derived state (fixed
   // polynomial coefficients only), so this avoids exposing Eos internals.
   EosType LocEosChoice = Eos::getInstance()->EosChoice;
   Teos10Eos LocTeos10(VertCoord::getDefault());
   OMEGA_SCOPE(LocAvgSfcTemp, AvgSfcTemperature);
   OMEGA_SCOPE(LocAvgSfcSalinity, AvgSfcSalinity);
   OMEGA_SCOPE(LocInSituTemp, InSituTempScratch);
   OMEGA_SCOPE(LocPracSalinity, PracSalinityScratch);

   // TEOS-10 conversion is applied once per coupling interval to the averaged
   // conservative temperature and absolute salinity. Therefore this computes
   // PtFromCt(mean(Sa), mean(Ct)), not mean(PtFromCt(Sa, Ct)); these are not
   // generally equivalent because the conversion is nonlinear. A true time
   // average of the converted quantity would require converting each timestep
   // before accumulating it. The same consideration applies to any future
   // nonlinear absolute-to-practical salinity conversion.
   parallelFor(
       {(int)AvgSfcTemperature.extent(0)}, KOKKOS_LAMBDA(int Cell) {
          const Real Ct         = LocAvgSfcTemp(Cell);
          const Real Sa         = LocAvgSfcSalinity(Cell);
          const Real Pt         = LocEosChoice == EosType::Teos10Eos
                                      ? LocTeos10.calcPtFromCt(Sa, Ct)
                                      : Ct;
          LocInSituTemp(Cell)   = Pt + TkFrz; // C to K temperature conversion
          LocPracSalinity(Cell) = Sa / Psu2Gpkg; // abs to prac sal. conversion
       });

   deepCopy(AvgSfcTemperatureH, InSituTempScratch);
   deepCopy(AvgSfcSalinityH, PracSalinityScratch);
   deepCopy(AvgExtraEnergyH, AvgExtraEnergy);

   // Retrieve the default horizontal mesh
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();

   OMEGA_SCOPE(LocNormalVelocity, AvgSfcNormalVelocity);
   OMEGA_SCOPE(LocSshGradEdge, AvgSfcSshGrad);
   OMEGA_SCOPE(LocReconZonal, ReconZonalScratch);
   OMEGA_SCOPE(LocReconMerid, ReconMeridScratch);

   VectorReconOnCell ReconCell(DefHorzMesh);

   parallelFor(
       {DefHorzMesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          ReconCell(LocReconZonal, LocReconMerid, ICell, LocNormalVelocity);
       });

   deepCopy(AvgSfcVelocityZonalH, ReconZonalScratch);
   deepCopy(AvgSfcVelocityMeridH, ReconMeridScratch);

   // Reuse the same scratch arrays; VectorReconOnCell overwrites every cell
   parallelFor(
       {DefHorzMesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          ReconCell(LocReconZonal, LocReconMerid, ICell, LocSshGradEdge);
       });

   deepCopy(AvgSfcSshGradZonalH, ReconZonalScratch);
   deepCopy(AvgSfcSshGradMeridH, ReconMeridScratch);

   // SSH is an instantaneous field, so we don't bother with a device mirror of
   // our own. Instead, copy from the VertCoord, which owns SSH, host array.
   VertCoord *DefVertCoord = VertCoord::getDefault();

   auto SSHCellOwned = Kokkos::subview(
       DefVertCoord->SshCell, std::make_pair(0, (int)InstSshCellH.extent(0)));

   deepCopy(InstSshCellH, SSHCellOwned);
}

// OcnToCpl fields need to begin a coupling interval with all values set to 0.
void OcnToCplFields::resetFields() {
   deepCopy(AvgSfcTemperature, 0.0_Real);
   deepCopy(AvgSfcSalinity, 0.0_Real);
   deepCopy(AvgSfcNormalVelocity, 0.0_Real);
   deepCopy(AvgSfcSshGrad, 0.0_Real);
   deepCopy(AvgExtraEnergy, 0.0_Real);
}
} // namespace OMEGA
