//===-- OMEGA Tracers Implementation ----------------------------*- C++ -*-===//
//
/// \file
/// \brief OMEGA Tracers Implementation
///
//
//===-----------------------------------------------------------------------===/

#include "Tracers.h"
#include "Config.h"
#include "Decomp.h"
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "TimeStepper.h"

namespace OMEGA {

// Initialize static member variables
std::vector<Array3DReal> Tracers::TracerArrays;
std::vector<HostArray3DReal> Tracers::TracerArraysH;

std::map<std::string, std::pair<I4, I4>> Tracers::TracerGroups;
std::map<std::string, I4> Tracers::TracerIndexes;
std::map<I4, std::string> Tracers::TracerNames;
std::vector<std::string> Tracers::TracerDimNames = {"NCells", "NVertLayers"};

Halo *Tracers::MeshHalo = nullptr;

I4 Tracers::NCellsOwned  = 0;
I4 Tracers::NCellsAll    = 0;
I4 Tracers::NCellsSize   = 0;
I4 Tracers::NTimeLevels  = 0;
I4 Tracers::NVertLayers  = 0;
I4 Tracers::CurTimeIndex = 0;
I4 Tracers::NumTracers   = 0;

//---------------------------------------------------------------------------
// Initialization
//---------------------------------------------------------------------------
void Tracers::init() {

   Error Err; // error code

   // Retrieve mesh cell/edge/vertex totals from Decomp
   HorzMesh *DefHorzMesh   = HorzMesh::getDefault();
   VertCoord *DefVertCoord = VertCoord::getDefault();

   NCellsOwned = DefHorzMesh->NCellsOwned;
   NCellsAll   = DefHorzMesh->NCellsAll;
   NCellsSize  = DefHorzMesh->NCellsSize;
   NVertLayers = DefVertCoord->NVertLayers;

   MeshHalo = Halo::getDefault();

   auto *DefTimeStepper = TimeStepper::getDefault();
   if (!DefTimeStepper) {
      ABORT_ERROR("TimeStepper needs to be initialized before Tracers");
   }
   NTimeLevels = DefTimeStepper->getNTimeLevels();

   if (NTimeLevels < 2) {
      ABORT_ERROR("Tracers: the number of time level is lower than 2");
   }

   CurTimeIndex = 0;

   // load Tracers configs
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TracersConfig("Tracers");
   Err += OmegaConfig->get(TracersConfig);
   CHECK_ERROR_ABORT(Err, "Tracers: Tracers group not found in Config");

   NumTracers     = 0;
   I4 TracerIndex = 0;

   // get tracers group and tracer names
   for (auto It = TracersConfig.begin(); It != TracersConfig.end(); ++It) {

      I4 GroupStartIndex = TracerIndex;

      std::string GroupName;
      Err += Config::getName(It, GroupName);
      CHECK_ERROR_ABORT(Err, "Tracers: error retrieving tracer group name");

      std::vector<std::string> _TracerNames;
      Err += TracersConfig.get(GroupName, _TracerNames);
      CHECK_ERROR_ABORT(Err,
                        "Tracers: error retrieving tracer names for group {}",
                        GroupName);

      for (auto _TracerName : _TracerNames) {
         TracerIndexes[_TracerName] = TracerIndex;
         TracerIndex++;
      }

      TracerGroups[GroupName] =
          std::pair<I4, I4>(GroupStartIndex, TracerIndex - GroupStartIndex);
      FieldGroup::create(GroupName);
   }

   // total number of tracers
   NumTracers = TracerIndex;

   // A tracer dimension is sometimes needed for aux variables so we
   // define it here
   std::shared_ptr<Dimension> TracerDim =
       Dimension::create("NTracers", NumTracers);

   // Initialize tracers arrays for device and host
   TracerArrays.resize(NTimeLevels);
   TracerArraysH.resize(NTimeLevels);

   // Allocate tracers data array and assign to tracers arrays
   for (I4 TimeIndex = 0; TimeIndex < NTimeLevels; ++TimeIndex) {
      TracerArrays[TimeIndex] =
          Array3DReal("TracerTime" + std::to_string(TimeIndex), NumTracers,
                      NCellsSize, NVertLayers);
      TracerArraysH[TimeIndex] =
          HostArray3DReal("TracerHTime" + std::to_string(TimeIndex), NumTracers,
                          NCellsSize, NVertLayers);
   }

   // Define tracers
   defineAllTracers();

   // Check if all tracers defined in config file are loaded
   if (TracerIndexes.size() != TracerNames.size()) {
      ABORT_ERROR("Tracer: not all tracers defined in config file is loaded.");
   }

   // Add Fields to FieldGroup
   // Create an overall Tracer and Restart group
   auto AllTracerGrp = FieldGroup::create("Tracers");
   if (!FieldGroup::exists("Restart"))
      auto RestartGroup = FieldGroup::create("Restart");

   for (const auto &GroupPair : TracerGroups) {
      auto GroupName                   = GroupPair.first;
      std::string TracerFieldGroupName = GroupName;
      auto TracerFieldGroup            = FieldGroup::get(TracerFieldGroupName);

      std::vector<std::string> _TracerNames;
      Err += TracersConfig.get(GroupName, _TracerNames);
      CHECK_ERROR_ABORT(Err,
                        "Tracers: error retrieving tracer names for group {}",
                        GroupName);

      for (auto _TracerName : _TracerNames) {
         std::string TracerFieldName = _TracerName;

         // add tracer Field to field group
         TracerFieldGroup->addField(TracerFieldName);

         // Add tracer Field to all tracer group
         AllTracerGrp->addField(TracerFieldName);

         // Add tracer Field to restart group
         FieldGroup::addFieldToGroup(TracerFieldName, "Restart");

         // Associate Field with data
         I4 TracerIndex                     = TracerIndexes[_TracerName];
         std::shared_ptr<Field> TracerField = Field::get(TracerFieldName);

         // Create a 2D subview by fixing the first dimension (TracerIndex)
         Array2DReal TracerSubview = Kokkos::subview(
             TracerArrays[CurTimeIndex], TracerIndex, Kokkos::ALL, Kokkos::ALL);
         TracerField->attachData<Array2DReal>(TracerSubview);
      }
   }

   return;
}

//---------------------------------------------------------------------------
// Define tracers
//---------------------------------------------------------------------------
I4 Tracers::define(const std::string &Name, const std::string &Description,
                   const std::string &Units, const std::string &StdName,
                   const Real ValidMin, const Real ValidMax,
                   const Real FillValue, I4 &Index) {

   // Do nothing if this tracer is not selected
   if (TracerIndexes.find(Name) == TracerIndexes.end()) {
      return 0;
   }

   auto TracerIndex = TracerIndexes[Name];

   // Return error if tracer already exists
   if (TracerNames.find(TracerIndex) != TracerNames.end()) {
      LOG_ERROR("Tracers: Tracer '{}' already exists", Name);
      return -1;
   }

   // set tracer index to name mapping
   TracerNames[TracerIndex] = Name;

   if (&Index != &IndxInvalid)
      Index = TracerIndex;

   // create a tracer field
   std::string TracerFieldName = Name;
   auto TracerField = Field::create(TracerFieldName, Description, Units,
                                    StdName, ValidMin, ValidMax, FillValue,
                                    TracerDimNames.size(), TracerDimNames);
   if (!TracerField) {
      LOG_ERROR("Tracers: Tracer field '{}' is not created", TracerFieldName);
      return -2;
   }

   return 0;
}

//---------------------------------------------------------------------------
// Deallocate tracer arrays
//---------------------------------------------------------------------------
I4 Tracers::clear() {

   // Deallocate memory for tracer arrays
   TracerArrays.clear();
   TracerArraysH.clear();

   TracerGroups.clear();
   TracerIndexes.clear();
   TracerNames.clear();

   NumTracers  = 0;
   NCellsOwned = 0;
   NCellsAll   = 0;
   NCellsSize  = 0;
   NTimeLevels = 0;
   NVertLayers = 0;

   return 0;
}

//---------------------------------------------------------------------------
// Query tracers
//---------------------------------------------------------------------------

I4 Tracers::getNumTracers() { return NumTracers; }

I4 Tracers::getIndex(I4 &TracerIndex, const std::string &TracerName) {
   // Check if tracer exists
   if (TracerIndexes.find(TracerName) != TracerIndexes.end()) {
      TracerIndex = TracerIndexes[TracerName];
      return 0; // Success
   }

   LOG_ERROR("Tracers: Tracer index for '{}' is not found.", TracerName);
   return -1; // Tracer not found
}

I4 Tracers::getName(std::string &TracerName, const I4 TracerIndex) {
   if (TracerNames.find(TracerIndex) != TracerNames.end()) {
      TracerName = TracerNames[TracerIndex];
      return 0; // Success
   }

   LOG_ERROR("Tracers: Tracer name for index '{}' is not found.", TracerIndex);
   return -1; // Tracer index not found
}

//---------------------------------------------------------------------------
// get Tracer arrays
// TimeLevel == [0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
Array3DReal Tracers::getAll(const I4 TimeLevel) {
   const I4 TimeIndex = getTimeIndex(TimeLevel);
   return TracerArrays[TimeIndex];
}

Array2DReal Tracers::getByIndex(const I4 TimeLevel, const I4 TracerIndex) {

   // Check if tracer index is valid
   OMEGA_REQUIRE(TracerIndex >= 0 && TracerIndex < NumTracers,
                 "Tracers: Tracer index {} is out of range", TracerIndex);

   const I4 TimeIndex = getTimeIndex(TimeLevel);
   return Kokkos::subview(TracerArrays[TimeIndex], TracerIndex, Kokkos::ALL,
                          Kokkos::ALL);
}

Array2DReal Tracers::getByName(const I4 TimeLevel,
                               const std::string &TracerName) {

   // Check if tracer exists
   OMEGA_REQUIRE(TracerIndexes.find(TracerName) != TracerIndexes.end(),
                 "Tracers: Tracer '{}' does not exist", TracerName);

   return getByIndex(TimeLevel, TracerIndexes[TracerName]);
}

HostArray3DReal Tracers::getAllHost(const I4 TimeLevel) {
   const I4 TimeIndex = getTimeIndex(TimeLevel);
   return TracerArraysH[TimeIndex];
}

HostArray2DReal Tracers::getHostByIndex(const I4 TimeLevel,
                                        const I4 TracerIndex) {

   // Check if tracer index is valid
   OMEGA_REQUIRE(TracerIndex >= 0 && TracerIndex < NumTracers,
                 "Tracers: Tracer index {} is out of range", TracerIndex);

   const I4 TimeIndex = getTimeIndex(TimeLevel);
   return Kokkos::subview(TracerArraysH[TimeIndex], TracerIndex, Kokkos::ALL,
                          Kokkos::ALL);
}
HostArray2DReal Tracers::getHostByName(const I4 TimeLevel,
                                       const std::string &TracerName) {
   // Check if tracer exists
   OMEGA_REQUIRE(TracerIndexes.find(TracerName) != TracerIndexes.end(),
                 "Tracers: Tracer '{}' does not exist", TracerName);

   // Get the index of the tracer
   I4 TracerIndex = TracerIndexes[TracerName];

   return getHostByIndex(TimeLevel, TracerIndex);
}

//---------------------------------------------------------------------------
// get Fields
//---------------------------------------------------------------------------
std::shared_ptr<Field> Tracers::getFieldByName(const std::string &TracerName) {
   // Check if tracer exists
   if (TracerIndexes.find(TracerName) == TracerIndexes.end()) {
      LOG_ERROR("Tracers: Tracer '{}' does not exist", TracerName);
      return nullptr;
   }

   return Field::get(TracerName);
}

std::shared_ptr<Field> Tracers::getFieldByIndex(const I4 TracerIndex) {
   // Check if tracer index is valid
   if (TracerIndex < 0 || TracerIndex >= NumTracers) {
      LOG_ERROR("Tracers: Tracer index {} is out of range", TracerIndex);
      return nullptr;
   }

   return getFieldByName(TracerNames[TracerIndex]);
}

//---------------------------------------------------------------------------
// Query tracer group
//---------------------------------------------------------------------------
std::vector<std::string> Tracers::getGroupNames() {
   std::vector<std::string> GroupNames;

   for (const auto &GroupPair : TracerGroups) {
      GroupNames.push_back(GroupPair.first);
   }

   return GroupNames;
}

I4 Tracers::getGroupRange(std::pair<I4, I4> &GroupRange,
                          const std::string &GroupName) {
   auto it = TracerGroups.find(GroupName);
   if (it != TracerGroups.end()) {
      GroupRange = it->second;
      return 0;
   }

   return -1;
}

bool Tracers::isGroupMemberByIndex(const I4 TracerIndex,
                                   const std::string GroupName) {
   auto it = TracerGroups.find(GroupName);
   if (it != TracerGroups.end()) {
      I4 StartIndex  = it->second.first;
      I4 GroupLength = it->second.second;
      return TracerIndex >= StartIndex &&
             TracerIndex < StartIndex + GroupLength;
   }

   return false;
}

bool Tracers::isGroupMemberByName(const std::string &TracerName,
                                  const std::string &GroupName) {

   I4 TracerIndex;
   if (getIndex(TracerIndex, TracerName) != 0) {
      return false;
   }

   return isGroupMemberByIndex(TracerIndex, GroupName);
}

//---------------------------------------------------------------------------
// deep copy at a time level
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
void Tracers::copyToDevice(const I4 TimeLevel) {

   const I4 TimeIndex = getTimeIndex(TimeLevel);

   deepCopy(TracerArrays[TimeIndex], TracerArraysH[TimeIndex]);
}

void Tracers::copyToHost(const I4 TimeLevel) {

   const I4 TimeIndex = getTimeIndex(TimeLevel);

   deepCopy(TracerArraysH[TimeIndex], TracerArrays[TimeIndex]);
}

//---------------------------------------------------------------------------
// halo exchange
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
I4 Tracers::exchangeHalo(const I4 TimeLevel) {

   I4 Err             = 0;
   const I4 TimeIndex = getTimeIndex(TimeLevel);

   Err = MeshHalo->exchangeFullArrayHalo(TracerArrays[TimeIndex], OnCell);
   if (Err != 0)
      return -1;

   return 0;
}

//---------------------------------------------------------------------------
//  update time level
//---------------------------------------------------------------------------
void Tracers::updateTimeLevels() {

   if (NTimeLevels == 1)
      ABORT_ERROR("Tracers: can't update time levels for NTimeLevels == 1");

   // Exchange halo
   exchangeHalo(1);

   CurTimeIndex = (CurTimeIndex + 1) % NTimeLevels;

   // Update TracerField data associations
   for (const auto &TracerPair : TracerIndexes) {
      auto TracerFieldName = TracerPair.first;
      auto TracerIndex     = TracerPair.second;

      std::shared_ptr<Field> TracerField = Field::get(TracerFieldName);

      Array2DReal TracerSubview = Kokkos::subview(
          TracerArrays[CurTimeIndex], TracerIndex, Kokkos::ALL, Kokkos::ALL);
      TracerField->attachData<Array2DReal>(TracerSubview);
   }

   return;
}

//---------------------------------------------------------------------------
// get index for time level
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
I4 Tracers::getTimeIndex(const I4 TimeLevel) {
   OMEGA_REQUIRE(NTimeLevels <= 1 ||
                     !(TimeLevel > 1 || (TimeLevel + NTimeLevels) <= 1),
                 "Tracers: Time level {} is out of range for NTimeLevels {}",
                 TimeLevel, NTimeLevels);

   return (TimeLevel + CurTimeIndex + NTimeLevels) % NTimeLevels;
}

} // namespace OMEGA
