//===-- OMEGA Tracers Implementation -----------------------------*- C++
//-*-===//
//
/// \file
/// \brief OMEGA Tracers Implementation
///
//
//===-----------------------------------------------------------------------===/

#include "Tracers.h"
#include "Decomp.h"
#include "IO.h"
#include "Logging.h"
#include "TimeStepper.h"

#include <iostream>

namespace OMEGA {

// Initialize static member variables
std::vector<Array3DReal> Tracers::TracerArrays;
std::vector<HostArray3DReal> Tracers::TracerArraysH;

std::map<std::string, std::pair<I4, I4>> Tracers::TracerGroups;
std::map<std::string, I4> Tracers::TracerIndexes;
std::map<I4, std::string> Tracers::TracerNames;
std::vector<std::string> Tracers::TracerDimNames = {"NCells", "NVertLevels"};

Halo *Tracers::MeshHalo = nullptr;

I4 Tracers::NCellsOwned  = 0;
I4 Tracers::NCellsAll    = 0;
I4 Tracers::NCellsSize   = 0;
I4 Tracers::NTimeLevels  = 0;
I4 Tracers::NVertLevels  = 0;
I4 Tracers::CurTimeIndex = 0;
I4 Tracers::NumTracers   = 0;

//---------------------------------------------------------------------------
// Initialization
//---------------------------------------------------------------------------
I4 Tracers::init() {

   I4 Err = 0;

   LOG_INFO("Tracers::init() called");

   // Retrieve mesh cell/edge/vertex totals from Decomp
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();

   NCellsOwned = DefHorzMesh->NCellsOwned;
   NCellsAll   = DefHorzMesh->NCellsAll;
   NCellsSize  = DefHorzMesh->NCellsSize;
   NVertLevels = DefHorzMesh->NVertLevels;

   MeshHalo = Halo::getDefault();

   auto *DefTimeStepper = TimeStepper::getDefault();
   if (!DefTimeStepper) {
      LOG_ERROR("TimeStepper needs to be initialized before Tracers");
      return -1;
   }
   NTimeLevels = DefTimeStepper->getNTimeLevels();

   if (NTimeLevels < 2) {
      LOG_ERROR("Tracers: the number of time level is lower than 2");
      return -2;
   }

   CurTimeIndex = 0;

   // load Tracers configs
   Config *OmegaConfig = Config::getOmegaConfig();
   Config TracersConfig("Tracers");
   Err = OmegaConfig->get(TracersConfig);
   if (Err != 0) {
      LOG_ERROR("Tracers: Tracers group not found in Config");
      return -3;
   }

   NumTracers     = 0;
   I4 TracerIndex = 0;

   // get tracers group and tracer names
   for (auto It = TracersConfig.begin(); It != TracersConfig.end(); ++It) {

      I4 GroupStartIndex = TracerIndex;

      std::string GroupName;
      I4 GroupNameErr = OMEGA::Config::getName(It, GroupName);
      if (GroupNameErr != 0) {
         LOG_ERROR("Tracers: {} tracer group name not found in TracersConfig",
                   GroupName);
         return -4;
      }

      std::vector<std::string> _TracerNames;
      I4 TracerNamesErr = TracersConfig.get(GroupName, _TracerNames);
      if (TracerNamesErr != 0) {
         LOG_ERROR("Tracers: {} group tracers not found in TracersConfig",
                   GroupName);
         return -5;
      }

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
                      NCellsSize, NVertLevels);
      TracerArraysH[TimeIndex] =
          HostArray3DReal("TracerHTime" + std::to_string(TimeIndex), NumTracers,
                          NCellsSize, NVertLevels);
   }

   // Define tracers
   defineAllTracers();

   // Check if all tracers defined in config file are loaded
   if (TracerIndexes.size() != TracerNames.size()) {
      LOG_ERROR("Tracer: not all tracers defined in config file is loaded.");
      return -6;
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
      TracersConfig.get(GroupName, _TracerNames);

      for (auto _TracerName : _TracerNames) {
         std::string TracerFieldName = _TracerName;

         // add tracer Field to field group
         Err = TracerFieldGroup->addField(TracerFieldName);
         if (Err != 0) {
            LOG_ERROR("Error adding {} to field group {}", TracerFieldName,
                      TracerFieldGroupName);
            return -7;
         }

         // Add tracer Field to all tracer group
         Err = AllTracerGrp->addField(TracerFieldName);
         if (Err != 0) {
            LOG_ERROR("Error adding {} to All Tracer group", TracerFieldName);
            return 8;
         }

         // Add tracer Field to restart group
         Err = FieldGroup::addFieldToGroup(TracerFieldName, "Restart");
         if (Err != 0) {
            LOG_ERROR("Error adding {} to Restart group", TracerFieldName);
            return 8;
         }

         // Associate Field with data
         I4 TracerIndex                     = TracerIndexes[_TracerName];
         std::shared_ptr<Field> TracerField = Field::get(TracerFieldName);

         // Create a 2D subview by fixing the first dimension (TracerIndex)
         Array2DReal TracerSubview = Kokkos::subview(
             TracerArrays[CurTimeIndex], TracerIndex, Kokkos::ALL, Kokkos::ALL);
         Err = TracerField->attachData<Array2DReal>(TracerSubview);
         if (Err != 0) {
            LOG_ERROR("Error attaching data array to field {}",
                      TracerFieldName);
            return -8;
         }
      }
   }

   return 0;
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

   LOG_INFO("Tracers::clear() called");

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
   NVertLevels = 0;

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
I4 Tracers::getAll(Array3DReal &TracerArray, const I4 TimeLevel) {
   I4 Err = 0;
   I4 TimeIndex;

   Err         = getTimeIndex(TimeIndex, TimeLevel);
   TracerArray = TracerArrays[TimeIndex];

   return Err;
}

I4 Tracers::getByIndex(Array2DReal &TracerArray, const I4 TimeLevel,
                       const I4 TracerIndex) {

   // Check if tracer index is valid
   if (TracerIndex < 0 || TracerIndex >= NumTracers) {
      LOG_ERROR("Tracers: Tracer index {} is out of range", TracerIndex);
      return -2;
   }

   I4 Err = 0;
   I4 TimeIndex;

   Err         = getTimeIndex(TimeIndex, TimeLevel);
   TracerArray = Kokkos::subview(TracerArrays[TimeIndex], TracerIndex,
                                 Kokkos::ALL, Kokkos::ALL);
   return Err;
}

I4 Tracers::getByName(Array2DReal &TracerArray, const I4 TimeLevel,
                      const std::string &TracerName) {
   // Check if tracer exists
   if (TracerIndexes.find(TracerName) == TracerIndexes.end()) {
      LOG_ERROR("Tracers: Tracer '{}' does not exist", TracerName);
      return -1;
   }

   int Err = getByIndex(TracerArray, TimeLevel, TracerIndexes[TracerName]);
   if (Err != 0)
      return -2;

   return 0;
}

I4 Tracers::getAllHost(HostArray3DReal &TracerArrayH, const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err          = getTimeIndex(TimeIndex, TimeLevel);
   TracerArrayH = TracerArraysH[TimeIndex];

   return Err;
}

I4 Tracers::getHostByIndex(HostArray2DReal &TracerArrayH, const I4 TimeLevel,
                           const I4 TracerIndex) {

   // Check if tracer index is valid
   if (TracerIndex < 0 || TracerIndex >= NumTracers) {
      LOG_ERROR("Tracers: Tracer index {} is out of range", TracerIndex);
      return -2;
   }

   I4 Err = 0;
   I4 TimeIndex;

   Err          = getTimeIndex(TimeIndex, TimeLevel);
   TracerArrayH = Kokkos::subview(TracerArraysH[TimeIndex], TracerIndex,
                                  Kokkos::ALL, Kokkos::ALL);
   return 0;
}

I4 Tracers::getHostByName(HostArray2DReal &TracerArrayH, const I4 TimeLevel,
                          const std::string &TracerName) {
   // Check if tracer exists
   if (TracerIndexes.find(TracerName) == TracerIndexes.end()) {
      LOG_ERROR("Tracers: Tracer '{}' does not exist", TracerName);
      return -1;
   }

   // Get the index of the tracer
   I4 TracerIndex = TracerIndexes[TracerName];
   I4 Err         = getHostByIndex(TracerArrayH, TimeLevel, TracerIndex);
   if (Err != 0)
      return -2;

   return 0;
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
I4 Tracers::copyToDevice(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err = getTimeIndex(TimeIndex, TimeLevel);
   deepCopy(TracerArrays[TimeIndex], TracerArraysH[TimeIndex]);

   return Err;
}

I4 Tracers::copyToHost(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err = getTimeIndex(TimeIndex, TimeLevel);
   deepCopy(TracerArraysH[TimeIndex], TracerArrays[TimeIndex]);

   return 0;
}

//---------------------------------------------------------------------------
// halo exchange
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
I4 Tracers::exchangeHalo(const I4 TimeLevel) {

   I4 Err = 0;
   I4 TimeIndex;

   Err = getTimeIndex(TimeIndex, TimeLevel);
   Err = MeshHalo->exchangeFullArrayHalo(TracerArrays[TimeIndex], OnCell);
   if (Err != 0)
      return -1;

   return 0;
}

//---------------------------------------------------------------------------
//  update time level
//---------------------------------------------------------------------------
I4 Tracers::updateTimeLevels() {

   if (NTimeLevels == 1) {
      LOG_ERROR("Tracers: can't update time levels for NTimeLevels == 1");
      return -1;
   }

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
      I4 Err = TracerField->attachData<Array2DReal>(TracerSubview);
      if (Err != 0) {
         LOG_ERROR("Error attaching data array to field {}", TracerFieldName);
         return Err;
      }
   }

   return 0;
}

//---------------------------------------------------------------------------
// get index for time level
// TimeLevel == [1:new, 0:current, -1:previous, -2:two times ago, ...]
//---------------------------------------------------------------------------
I4 Tracers::getTimeIndex(I4 &TimeIndex, const I4 TimeLevel) {

   // Check if time level is valid
   if (NTimeLevels > 1 && (TimeLevel > 1 || (TimeLevel + NTimeLevels) <= 1)) {
      LOG_ERROR("Tracers: Time level {} is out of range for NTimeLevels {}",
                TimeLevel, NTimeLevels);
      return -1;
   }

   TimeIndex = (TimeLevel + CurTimeIndex + NTimeLevels) % NTimeLevels;

   return 0;
}

} // namespace OMEGA
