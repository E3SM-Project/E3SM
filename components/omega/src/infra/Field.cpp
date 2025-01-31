//===-- infra/Field.cpp - OMEGA field class implementation ------*- C++ -*-===//
//
// Implements Field classes and methods
//
// This file implements classes and methods for Fields. Fields are the
// data and metadata for all fields that may be a part of an IO stream or any
// other module where it might be important to keep the data and metadata
// together. Fields are defined/created withing the modules that own them. It
// is also possible to define field groups to provide a shortcut for fields
// that are commonly used together. Once a field is defined/created, a data
// array can be attached detached or swapped to hold the latest data values
// (esp for fields with multiple time levels).
//
//===----------------------------------------------------------------------===//

#include "Field.h"
#include "DataTypes.h"
#include "Dimension.h"
#include "IO.h"
#include "Logging.h"
#include <iostream>
#include <map>
#include <memory>
#include <set>

namespace OMEGA {

// Initialize static variables
std::map<std::string, std::shared_ptr<Field>> Field::AllFields;
std::map<std::string, std::shared_ptr<FieldGroup>> FieldGroup::AllGroups;

//------------------------------------------------------------------------------
// Initializes the fields for global code and simulation metadata
int Field::init(const Clock *ModelClock // [in] default model clock
) {

   int Err = 0;

   // Create two fields to hold global metadata associated with either the
   // code or the simulation/experiment
   std::shared_ptr<Field> CodeField = create(CodeMeta);
   std::shared_ptr<Field> SimField  = create(SimMeta);

   // Define an unlimited time dimension for many time-dependent fields
   // for CF-compliant output
   std::shared_ptr<Dimension> TimeDim =
       Dimension::create("time", IO::Unlimited);

   // Define a time field with required metadata for CF-compliant output
   // It is defined here as a scalar field but the time axis will be added
   // during IO
   TimeInstant StartTime    = ModelClock->getStartTime();
   std::string StartTimeStr = StartTime.getString(4, 0, " ");
   std::string UnitString   = "seconds since " + StartTimeStr;
   CalendarKind CalKind     = Calendar::getKind();
   std::string CalName      = CalendarCFName[CalKind];
   std::vector<std::string> DimNames; // empty dim names vector
   std::shared_ptr<Field> TimeField = create("time", "time", UnitString, "time",
                                             0.0, 1.e20, -9.99e30, 0, DimNames);
   TimeField->addMetadata("calendar", CalName);

   return Err;
}

//------------------------------------------------------------------------------
// Create/destroy/query fields
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Checks to see if a field has been defined
bool Field::exists(const std::string &FieldName // [in] name of field
) {
   return (AllFields.find(FieldName) != AllFields.end());
}

//------------------------------------------------------------------------------
// Creates a field with standard metadata. This is the preferred
// interface for most fields in Omega. It enforces a list of required
// metadata. Note that if input parameters do not exist
// (eg stdName) or do not make sense (eg min/max or fill) for a
// given field, empty or 0 entries can be provided. Actual field data is
// attached in a separate call and additional metadata can be added later.

std::shared_ptr<Field>
Field::create(const std::string &FieldName,   // [in] Name of variable/field
              const std::string &Description, // [in] long Name or description
              const std::string &Units,       // [in] units
              const std::string &StdName,     // [in] CF standard Name
              const std::any ValidMin,        // [in] min valid field value
              const std::any ValidMax,        // [in] max valid field value
              const std::any FillValue, // [in] scalar for undefined entries
              const int NumDims,        // [in] number of dimensions
              const std::vector<std::string> &Dimensions, // [in] dim names
              bool TimeDependent // [in] flag for time dependent field
) {

   // Check to make sure a field of that name has not already been defined
   if (exists(FieldName)) {
      LOG_ERROR("Attempt to create field {} failed."
                " Field with that name already exists.",
                FieldName);
      return nullptr;
   }

   // Create an empty Field
   auto ThisField = std::make_shared<Field>();

   // Add field name to the instance (also added as metadata below)
   ThisField->FldName = FieldName;

   // Create empty metadata map: (name, value) pairs
   ThisField->FieldMeta;

   // Add standard metadata. For some CF standard attributes, we
   // also duplicate the metadata under the CF standard attribute name
   ThisField->FieldMeta["Name"]          = FieldName;
   ThisField->FieldMeta["name"]          = FieldName;
   ThisField->FieldMeta["Description"]   = Description;
   ThisField->FieldMeta["long_name"]     = Description;
   ThisField->FieldMeta["Units"]         = Units;
   ThisField->FieldMeta["units"]         = Units;
   ThisField->FieldMeta["StdName"]       = StdName;
   ThisField->FieldMeta["standard_name"] = StdName;
   ThisField->FieldMeta["ValidMin"]      = ValidMin;
   ThisField->FieldMeta["valid_min"]     = ValidMin;
   ThisField->FieldMeta["ValidMax"]      = ValidMax;
   ThisField->FieldMeta["valid_max"]     = ValidMax;
   ThisField->FieldMeta["FillValue"]     = FillValue;
   ThisField->FieldMeta["_FillValue"]    = FillValue;

   // Set the time-dependent flag
   ThisField->TimeDependent = TimeDependent;

   // Number of dimensions for the field
   ThisField->NDims = NumDims;

   // Dimension names for retrieval of dimension info
   // These must be in the same index order as the stored data
   // Also determine whether this is a distributed field - true if any of
   // the dimensions are distributed.
   ThisField->Distributed = false;
   ThisField->DimNames;
   if (NumDims > 0) {
      ThisField->DimNames.resize(NumDims);
      for (int I = 0; I < NumDims; ++I) {
         ThisField->DimNames[I] = Dimensions[I];
         if (Dimension::isDistributedDim(Dimensions[I]))
            ThisField->Distributed = true;
      }
   }

   // Initialize data info to Unknown - these will be determined when the
   // field data is attached.
   ThisField->DataType = ArrayDataType::Unknown;
   ThisField->MemLoc   = ArrayMemLoc::Unknown;

   // Initialize Data pointer as null until data is actually attached.
   ThisField->DataArray = nullptr;

   // Add to list of fields and return
   AllFields[FieldName] = ThisField;
   return ThisField;
}

//------------------------------------------------------------------------------
// Creates an empty field with a given name. This interface should only
// be used for collections of metadata with no attached data array.
// Two standard fields called CodeMeta and SimMeta are defined on init
// to hold global metadata not attached to fields.
std::shared_ptr<Field>
Field::create(const std::string &FieldName // [in] Name of field
) {

   // Check to make sure a field of that name has not already been defined
   if (exists(FieldName)) {
      LOG_ERROR("Attempt to create field {} failed."
                " Field with that name already exists.",
                FieldName);
      return nullptr;
   }

   // Create an empty Field
   auto ThisField = std::make_shared<Field>();

   // Field name
   ThisField->FldName = FieldName;

   // Metadata (name, value) pairs for descriptive metadata
   ThisField->FieldMeta;

   /// Number of dimensions is 0 for this field
   ThisField->NDims = 0;

   // Dimension name vector is empty
   ThisField->DimNames;

   // Initialize to Unknown or null - no data is attached
   ThisField->DataType  = ArrayDataType::Unknown;
   ThisField->MemLoc    = ArrayMemLoc::Unknown;
   ThisField->DataArray = nullptr;

   // Add to list of fields and return
   AllFields[FieldName] = ThisField;
   return ThisField;
}

//------------------------------------------------------------------------------
// Removes a single IOField from the list of available fields
// That process also decrements the reference counters for the
// shared pointers and removes them if those counters reach 0.
int Field::destroy(const std::string &FieldName // [in] Name of Field to remove
) {
   int Err = 0;

   // Check that the group exists
   if (exists(FieldName)) {

      // Erase the group from the list of all groups
      AllFields.erase(FieldName);

      // Group does not exist, exit with error
   } else {

      LOG_ERROR("Failed to destroy the field {}: field not found", FieldName);
      Err = 1;
   }

   return Err;
}

//------------------------------------------------------------------------------
// Removes all Fields. This must be called before exiting environments
// This removes all fields from the map structure and also
// decrements the reference counter for the shared pointers,
// removing them if the count has reached 0.
void Field::clear() { AllFields.clear(); }

//------------------------------------------------------------------------------
// Retrieve a field pointer by name
std::shared_ptr<Field>
Field::get(const std::string &FieldName // [in] name of field
) {

   if (exists(FieldName)) {
      return AllFields[FieldName];
   } else {
      LOG_ERROR("Unable to retrieve Field {}. Field not found.", FieldName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------
// Get field name
std::string Field::getName() const { return FldName; }

//------------------------------------------------------------------------------
// Determine type of a given field from instance
ArrayDataType Field::getType() const { return DataType; }

//------------------------------------------------------------------------------
// Determine type of a given field by name
ArrayDataType
Field::getFieldType(const std::string &FieldName // [in] name of field
) {
   // Retrieve field with given name
   auto It = AllFields.find(FieldName);

   // If found, use member function to return type
   if (It != AllFields.end()) {
      return It->second->getType();
   } else { // Could not find field with that name
      LOG_ERROR("Unable to retrieve field type. Field {} has not been created.",
                FieldName);
      return ArrayDataType::Unknown;
   }
}

//------------------------------------------------------------------------------
// Determine memory location of data from instance
ArrayMemLoc Field::getMemoryLocation() const { return MemLoc; }

//------------------------------------------------------------------------------
// Determine memory location of data by field name
ArrayMemLoc
Field::getFieldMemoryLocation(const std::string &FieldName // [in] name of field
) {
   // Retrieve field with given name
   auto It = AllFields.find(FieldName);

   // If found, use member function to return location
   if (It != AllFields.end()) {
      return It->second->getMemoryLocation();
   } else { // Could not find field with that name
      LOG_ERROR("Unable to retrieve memory location."
                " Field {} has not been created.",
                FieldName);
      return ArrayMemLoc::Unknown;
   }
}

//------------------------------------------------------------------------------
// Query whether field is located on host from instance
bool Field::isOnHost() const {
   if (MemLoc == ArrayMemLoc::Host or MemLoc == ArrayMemLoc::Both) {
      return true;
   } else {
      return false;
   }
}

//------------------------------------------------------------------------------
// Query whether field is located on host by field name
bool Field::isFieldOnHost(const std::string &FieldName // [in] name of field
) {
   // Retrieve field with given name
   auto It = AllFields.find(FieldName);

   // If found, use member function to return host flag
   if (It != AllFields.end()) {
      return It->second->isOnHost();
   } else { // Could not find field with that name
      LOG_ERROR("Unable to determine whether field is on host."
                " Field {} has not been created.",
                FieldName);
      return false;
   }
}

//------------------------------------------------------------------------------
// Returns the number of dimensions for the field
int Field::getNumDims() const { return NDims; }

//------------------------------------------------------------------------------
// Determines whether the field is time dependent and requires the unlimited
// time dimension during IO
bool Field::isTimeDependent() const { return TimeDependent; }

//------------------------------------------------------------------------------
// Determinse whether a field is distributed across tasks or whether a copy
// is entirely local. This is needed to determine whether a parallel IO or
// a non-distributed IO will be used.
bool Field::isDistributed() const { return Distributed; }

//------------------------------------------------------------------------------
// Returns a vector of dimension names associated with each dimension
// of an array field. Returns an error code.
int Field::getDimNames(
    std::vector<std::string> &Dimensions // [out] list of dimensions
) const {
   int Err = 0;

   // Make sure vector is correct length
   if (NDims > 0) {
      Dimensions.resize(NDims);
   } else {
      LOG_ERROR("Unable to retrieve dimension names for Field {}. NDims < 1",
                FldName);
      Err = -1;
   }

   // Fill vector
   for (int I = 0; I < NDims; ++I) {
      Dimensions[I] = DimNames[I];
   }

   return Err;
}

//------------------------------------------------------------------------------
// Metadata functions
//------------------------------------------------------------------------------
// Retrieves all Metadata for a Field instance. Returns a
// shared pointer to the map structure with metadata.
std::shared_ptr<Metadata> Field::getAllMetadata() const {
   return std::make_shared<Metadata>(FieldMeta);
}

//------------------------------------------------------------------------------
// Retrieves all Metadata for a Field given the field name. Returns a
// shared pointer to the map structure with metadata.
std::shared_ptr<Metadata>
Field::getFieldMetadata(const std::string &FieldName // [in] name of field
) {
   // Retrieve field with given name
   auto It = AllFields.find(FieldName);

   // If found, call member function to return Metadata
   if (It != AllFields.end()) {
      return It->second->getAllMetadata();
   } else { // Field not found
      LOG_ERROR("Unable to retrieve Metadata for Field {}."
                "No Field with that name could be found.",
                FieldName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------
// Checks for the existence of a metadata entry with the given name
bool Field::hasMetadata(
    const std::string &MetaName // [in] Name of metadata entry
) const {
   return (FieldMeta.find(MetaName) != FieldMeta.end());
}

//------------------------------------------------------------------------------
// Adds a metadata entry with the (name,value) pair
int Field::addMetadata(const std::string &MetaName, // [in] Name of new metadata
                       const std::any Value         // [in] Value of metadata
) {
   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      LOG_ERROR("Failed to add the metadata {} to field {} because the field "
                "already has an entry with that name.",
                MetaName, FldName);

      RetVal = -1;

   } else {
      FieldMeta[MetaName] = Value;
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Adds multiple metadata entries with a list of (name,value) pairs
int Field::addMetadata(
    const std::initializer_list<std::pair<std::string, std::any>> &MetaPairs) {

   int RetVal = 0;

   // Loop through pairs and add each individually
   for (const auto &MetaPair : MetaPairs) {
      RetVal += addMetadata(MetaPair.first, MetaPair.second);
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Updates a metadata entry with a new value
int Field::updateMetadata(const std::string &MetaName, // [in] Name of metadata
                          const std::any Value         // [in] Value of metadata
) {
   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      FieldMeta[MetaName] = Value;

   } else {
      LOG_ERROR("Failed to update metadata {} for field {} because the field "
                "has no entry with that name.",
                MetaName, FldName);

      RetVal = -1;
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Removes a metadata entry with the given name
int Field::removeMetadata(
    const std::string &MetaName // [in] Name of metadata entry to remove
) {

   int RetVal = 0;

   // Make sure metadata exists
   if (hasMetadata(MetaName)) {
      if (FieldMeta.erase(MetaName) != 1) {
         LOG_ERROR("Unknown error removing metadata {} in field {} ", MetaName,
                   FldName);
         RetVal = -2;
      }

   } else {
      LOG_ERROR("Failed to remove metadata {} for the field {}: "
                "entry does not exist.",
                MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;

} // end removeMetadata

//----------------------------------------------------------------------------//
// Removes all defined Metadata
void Field::removeAllMetadata() { FieldMeta.clear(); }

//------------------------------------------------------------------------------
// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an I4 type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    I4 &Value                    // [out] I4 Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<I4>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an I8 type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    I8 &Value                    // [out] I8 Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<I8>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an R4 type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    R4 &Value                    // [out] R4 Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<R4>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an R8 type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    R8 &Value                    // [out] R8 Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<R8>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an bool type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    bool &Value                  // [out] bool Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<bool>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an string type
int Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    std::string &Value           // [out] string Value of metadata
) {

   int RetVal = 0;

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<std::string>(FieldMeta[MetaName]);
   } else {
      LOG_ERROR("Metadata {} does not exist for field {}.", MetaName, FldName);
      RetVal = -1;
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Attach/retrieve data
//------------------------------------------------------------------------------
// These are template functions so implementation in header files
// Field::attachData
// Field::attachFieldData
// Field::getData
// Field::getFieldData
//------------------------------------------------------------------------------
// Field Group functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Creates an empty field group with a given name
std::shared_ptr<FieldGroup>
FieldGroup::create(const std::string &GroupName // [in] Name of group to create
) {

   // Check to make sure group has not already been defined
   if (exists(GroupName)) {
      LOG_ERROR("Attempt to create FieldGroup {} failed: group already exists",
                GroupName);
      return nullptr;

      // Create an empty group with the given name
   } else {

      auto Group     = std::make_shared<FieldGroup>();
      Group->GrpName = GroupName;
      // add to list of groups
      AllGroups[GroupName] = Group;

      return Group;
   }

} // end create group

//------------------------------------------------------------------------------
// Removes a field group
int FieldGroup::destroy(
    const std::string &GroupName // [in] Name of group to destroy
) {

   int RetVal = 0;

   // Check that the group exists
   if (exists(GroupName)) {

      // Erase the group from the list of all groups
      if (AllGroups.erase(GroupName) != 1) {
         LOG_ERROR("Unknown error trying to erase field group {}.", GroupName);
         RetVal = -2;
      }

      // Group does not exist, exit with error
   } else {

      LOG_ERROR("Failed to destroy the field group {}: group not found",
                GroupName);
      RetVal = -1;
   }

   return RetVal;

} // end destroy group

//------------------------------------------------------------------------------
// Removes all defined field groups
void FieldGroup::clear() { AllGroups.clear(); }

//------------------------------------------------------------------------------
// Determines whether a group of a given name exists
bool FieldGroup::exists(const std::string &GroupName // [in] Name of group
) {
   return (FieldGroup::AllGroups.find(GroupName) !=
           FieldGroup::AllGroups.end());
}

//------------------------------------------------------------------------------
// Determines whether a field of a given name exists in the group instance
bool FieldGroup::hasField(const std::string &FieldName // [in] Name of field
) const {
   return (Fields.find(FieldName) != Fields.end());
}

//------------------------------------------------------------------------------
// Determines whether a field of a given name exists in the group
// based on group name.
bool FieldGroup::isFieldInGroup(
    const std::string &FieldName, // [in] Name of field
    const std::string &GroupName  // [in] Name of group
) {

   // Check that the group exists
   if (exists(GroupName)) {

      // Now check for field in group
      auto Group = AllGroups[GroupName];
      return Group->hasField(FieldName);

      // Group not found
   } else {
      LOG_ERROR("Attempt to search for field {} in group {} failed: "
                "Group does not exist",
                FieldName, GroupName);
      return false;
   }
}

//------------------------------------------------------------------------------
// Adds a field to the group instance
int FieldGroup::addField(
    const std::string &FieldName // [in] Name of field to add
) {

   int RetVal = 0;

   // Add field name to the list of fields. If this is a duplicate, the
   // insert function knows not to repeat an entry so nothing happens.
   Fields.insert(FieldName);

   return RetVal;
}

//------------------------------------------------------------------------------
// Adds a field to a group based on group name
int FieldGroup::addFieldToGroup(
    const std::string &FieldName, // [in] Name of field to add
    const std::string &GroupName  // [in] Name of group
) {

   int RetVal = 0;

   // Check that the group exists
   if (exists(GroupName)) {

      auto Group = AllGroups[GroupName];
      Group->addField(FieldName);

   } else { // group does not exist
      LOG_ERROR("Failed to add field {} to group {}: Group does not exist",
                FieldName, GroupName);
      RetVal = 1;
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Removes a field from a group instance
int FieldGroup::removeField(
    const std::string &FieldName // [in] Name of field to remove
) {

   int RetVal = 0;

   // Removing a field is simply erasing the field from list
   Fields.erase(FieldName);

   return RetVal;
}

//------------------------------------------------------------------------------
// Removes a field from a group with a given name
int FieldGroup::removeFieldFromGroup(
    const std::string &FieldName, // [in] Name of field to remove
    const std::string &GroupName  // [in] Name of group holding field
) {

   int RetVal = 0;

   // Check that the group exists
   if (exists(GroupName)) {

      auto Group = AllGroups[GroupName];
      Group->removeField(FieldName);

   } else { // group does not exist
      LOG_ERROR("Failed to remove field {} from group {}: Group does not exist",
                FieldName, GroupName);
      RetVal = 1;
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Retrieves a pointer to a field group.
std::shared_ptr<FieldGroup>
FieldGroup::get(const std::string &GroupName // [in] Name of group to retrieve
) {

   // Check to make sure group exists
   if (exists(GroupName)) {
      // Return the pointer from the list of groups
      return AllGroups[GroupName];

   } else { // group does not exist
      LOG_ERROR("Failed to retrieve FieldGroup {}: group does not exist.",
                GroupName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------
// Returns the list of fields in a group instance. This function returns
// a copy of the field list.
std::set<std::string> FieldGroup::getFieldList() const {

   std::set<std::string> List;
   for (auto It = Fields.begin(); It != Fields.end(); ++It) {
      List.insert(*It);
   }

   return List;
}

//------------------------------------------------------------------------------
// Returns the list of fields in a group based on group name
std::set<std::string> FieldGroup::getFieldListFromGroup(
    const std::string &GroupName // [in] Name of group holding fields
) {

   std::set<std::string> List; // default empty list

   // Check for valid group name
   if (FieldGroup::exists(GroupName)) {

      // Retrieve group
      std::shared_ptr<FieldGroup> ThisGroup = AllGroups[GroupName];

      // Return list from member function
      List = ThisGroup->getFieldList();

   } else { // Group not found

      LOG_ERROR("Error getting field list from group {}. Group does not exist",
                GroupName);
   }

   return List;
}

//------------------------------------------------------------------------------
// Retrieves a full field from a group instance
std::shared_ptr<Field> FieldGroup::getField(
    const std::string &FieldName // [in] Name of field to retrieve
) const {

   // Check to see if the group has field
   if (hasField(FieldName)) {
      return Field::get(FieldName);
   } else {
      LOG_ERROR("Cannot find Field {} in Group {}", FieldName, GrpName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------
// Retrieves a full field from a group with a given name
std::shared_ptr<Field> FieldGroup::getFieldFromGroup(
    const std::string &FieldName, // [in] Name of field to retrieve
    const std::string &GroupName  // [in] Name of group holding field
) {

   // Check to see of group exists
   if (exists(GroupName)) {
      // Retrieve group - getField function performs checks
      std::shared_ptr<FieldGroup> Group = AllGroups[GroupName];
      return Group->getField(FieldName);
   } else {
      LOG_ERROR("Unable to retrieve Field {} from Group {}. Group not found.",
                FieldName, GroupName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------

} // namespace OMEGA
//===----------------------------------------------------------------------===//
