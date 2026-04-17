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
#include "Error.h"
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
void Field::init(const Clock *ModelClock // [in] default model clock
) {

   // Create two fields to hold global metadata associated with either the
   // code or the simulation/experiment
   std::shared_ptr<Field> CodeField = create(CodeMeta);
   std::shared_ptr<Field> SimField  = create(SimMeta);

   // Define an unlimited time dimension for many time-dependent fields
   // for CF-compliant output
   auto TimeDim = Dimension::create("time", IO::Unlimited);

   // Define a time field with required metadata for CF-compliant output
   // It is defined here as a scalar field but the time axis will be added
   // during IO
   TimeInstant StartTime    = ModelClock->getStartTime();
   std::string StartTimeStr = StartTime.getString(4, 0, " ");
   std::string UnitString   = "seconds since " + StartTimeStr;
   CalendarKind CalKind     = Calendar::getKind();
   std::string CalName      = CalendarCFName[CalKind];
   std::vector<std::string> DimNamesTmp; // empty dim names vector
   std::shared_ptr<Field> TimeField =
       create("time", "time", UnitString, "time", 0.0, 1.e20, -9.99e30, 0,
              DimNamesTmp, true, true);
   TimeField->addMetadata("calendar", CalName);
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
              bool TimeDependent,  // [in] flag for time dependent field
              bool RetainPrecision // [in] flag to retain full prec in IO
) {

   // Check to make sure a field of that name has not already been defined
   if (exists(FieldName))
      ABORT_ERROR("Attempt to create field {} failed."
                  " Field with that name already exists.",
                  FieldName);

   // Create an empty Field
   auto ThisField = std::make_shared<Field>();

   // Add field name to the instance (also added as metadata below)
   ThisField->FldName = FieldName;

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

   // Set the flag to retain full precision in IO operations even
   // if reduced precision for the stream/file is requested
   ThisField->RetainPrecision = RetainPrecision;

   // Number of dimensions for the field
   ThisField->NDims = NumDims;

   // Dimension names for retrieval of dimension info
   // These must be in the same index order as the stored data
   // Also determine whether this is a distributed field - true if any of
   // the dimensions are distributed.
   ThisField->Distributed = false;
   if (NumDims > 0) {
      if (static_cast<int>(Dimensions.size()) != NumDims)
         ABORT_ERROR("Attempt to create field {} failed. Expected {} "
                     "dimension names but received {}.",
                     FieldName, NumDims, Dimensions.size());

      ThisField->DimNames.resize(NumDims);
      for (int I = 0; I < NumDims; ++I) {
         if (Dimensions[I].empty())
            ABORT_ERROR("Attempt to create field {} failed. Dimension {} "
                        "has an empty name.",
                        FieldName, I);

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
   if (exists(FieldName))
      ABORT_ERROR("Attempt to create field {} failed."
                  " Field with that name already exists.",
                  FieldName);

   // Create an empty Field
   auto ThisField = std::make_shared<Field>();

   // Field name
   ThisField->FldName = FieldName;

   // Number of dimensions is 0 for this field
   ThisField->NDims = 0;

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
void Field::destroy(const std::string &FieldName // [in] Name of Field to remove
) {

   // Check that the group exists
   if (exists(FieldName)) {

      // Erase the group from the list of all groups
      AllFields.erase(FieldName);

      // Group does not exist, emit a warning but continue
   } else {
      LOG_WARN("Failed to destroy the field {}: field not found", FieldName);
   }

   return;
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

   // check that the field exists
   if (!exists(FieldName))
      ABORT_ERROR("Unable to retrieve Field {}. Field not found.", FieldName);

   return AllFields[FieldName];
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

   // If not found, abort
   if (It == AllFields.end())
      ABORT_ERROR("Field getFieldType error: Field {} does not exist.",
                  FieldName);

   // Otherwise use member function to return type
   return It->second->getType();
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

   // If not found, abort
   if (It == AllFields.end())
      ABORT_ERROR("Field getFieldMemoryLocation: Field {} does not exist",
                  FieldName);

   // Use member function to return location
   return It->second->getMemoryLocation();
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

   // If field not found, abort with error
   if (It == AllFields.end())
      ABORT_ERROR("Unable to determine whether field is on host."
                  " Field {} does not exist.",
                  FieldName);

   // Use member function to return host flag
   return It->second->isOnHost();
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
// Determinse whether full precision should be retained in IO operations
// when reduced precision for the stream/file is requested
bool Field::retainPrecision() const { return RetainPrecision; }

//------------------------------------------------------------------------------
// Returns a vector of dimension names associated with each dimension
// of an array field. Returns an error code.
void Field::getDimNames(
    std::vector<std::string> &Dimensions // [out] list of dimensions
) const {

   // Make sure vector is correct length
   if (NDims > 0) {
      Dimensions.resize(NDims);
   } else {
      ABORT_ERROR("Unable to retrieve dimension names for Field {}. NDims < 1",
                  FldName);
   }

   // Fill vector
   for (int I = 0; I < NDims; ++I) {
      Dimensions[I] = DimNames[I];
   }
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
   if (It == AllFields.end())
      ABORT_ERROR("Unable to retrieve Metadata for Field {}."
                  " Field does not exist.",
                  FieldName);

   return It->second->getAllMetadata();
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
void Field::addMetadata(const std::string &MetaName, // [in] Name of metadata
                        const std::any Value         // [in] Value of metadata
) {

   // If metadata already exists, abort (should use update instead)
   if (hasMetadata(MetaName))
      ABORT_ERROR("Failed to add the metadata {} to field {}. The field "
                  "already has an entry with that name. Use update instead.",
                  MetaName, FldName);

   FieldMeta[MetaName] = Value;

   return;
}

//------------------------------------------------------------------------------
// Adds multiple metadata entries with a list of (name,value) pairs
void Field::addMetadata(
    const std::initializer_list<std::pair<std::string, std::any>> &MetaPairs) {

   // Loop through pairs and add each individually
   for (const auto &MetaPair : MetaPairs) {
      addMetadata(MetaPair.first, MetaPair.second);
   }

   return;
}

//------------------------------------------------------------------------------
// Updates a metadata entry with a new value
void Field::updateMetadata(const std::string &MetaName, // [in] Metadata name
                           const std::any Value         // [in] Metadata value
) {

   if (hasMetadata(MetaName)) {
      FieldMeta[MetaName] = Value;

   } else {
      ABORT_ERROR("Unable to update metadata {} for field {}. "
                  "Field has no entry with that name.",
                  MetaName, FldName);
   }

   return;
}

//------------------------------------------------------------------------------
// Removes a metadata entry with the given name
void Field::removeMetadata(
    const std::string &MetaName // [in] Name of metadata entry to remove
) {

   // Make sure metadata exists
   if (hasMetadata(MetaName)) {
      FieldMeta.erase(MetaName);

   } else {
      LOG_WARN("Failed to remove metadata {} for the field {}: "
               "entry does not exist.",
               MetaName, FldName);
   }

   return;

} // end removeMetadata

//----------------------------------------------------------------------------//
// Removes all defined Metadata
void Field::removeAllMetadata() { FieldMeta.clear(); }

//------------------------------------------------------------------------------
// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an I4 type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    I4 &Value                    // [out] I4 Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<I4>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an I8 type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    I8 &Value                    // [out] I8 Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<I8>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an R4 type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    R4 &Value                    // [out] R4 Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<R4>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an R8 type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    R8 &Value                    // [out] R8 Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<R8>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an bool type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    bool &Value                  // [out] bool Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<bool>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
}

// Retrieves the value of the metadata associated with a given name
// This specific version of the overloaded interface coerces the value
// to an string type
Error Field::getMetadata(
    const std::string &MetaName, // [in] Name of metadata to get
    std::string &Value           // [out] string Value of metadata
) {

   Error Err; // returned error code, default success

   if (hasMetadata(MetaName)) {
      Value = std::any_cast<std::string>(FieldMeta[MetaName]);
   } else {
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "Metadata {} does not exist for field {}.", MetaName,
                   FldName);
   }

   return Err;
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
   if (exists(GroupName))
      ABORT_ERROR("Attempt to create FieldGroup {} that already exists",
                  GroupName);

   // Create an empty group with the given name
   auto Group     = std::make_shared<FieldGroup>();
   Group->GrpName = GroupName;
   // add to list of groups
   AllGroups[GroupName] = Group;

   return Group;

} // end create group

//------------------------------------------------------------------------------
// Removes a field group
void FieldGroup::destroy(
    const std::string &GroupName // [in] Name of group to destroy
) {

   // Check that the group exists
   if (exists(GroupName)) {

      // Erase the group from the list of all groups
      AllGroups.erase(GroupName);

      // Group does not exist, abort
   } else {
      ABORT_ERROR("Failed to destroy the field group {}: group not found",
                  GroupName);
   }

   return;

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
   if (!exists(GroupName))
      ABORT_ERROR("Attempt to search for field {} in group {} failed: "
                  "Group does not exist",
                  FieldName, GroupName);

   // Now check for field in group
   auto Group = AllGroups[GroupName];
   return Group->hasField(FieldName);
}

//------------------------------------------------------------------------------
// Adds a field to the group instance
void FieldGroup::addField(
    const std::string &FieldName // [in] Name of field to add
) {

   // Add field name to the list of fields. If this is a duplicate, the
   // insert function knows not to repeat an entry so nothing happens.
   Fields.insert(FieldName);
}

//------------------------------------------------------------------------------
// Adds a field to a group based on group name
void FieldGroup::addFieldToGroup(
    const std::string &FieldName, // [in] Name of field to add
    const std::string &GroupName  // [in] Name of group
) {

   // Check that the group exists
   if (!exists(GroupName))
      ABORT_ERROR("Failed to add field {} to group {}: Group does not exist",
                  FieldName, GroupName);

   // Add field name to group's list
   auto Group = AllGroups[GroupName];
   Group->addField(FieldName);

   return;
}

//------------------------------------------------------------------------------
// Removes a field from a group instance
void FieldGroup::removeField(
    const std::string &FieldName // [in] Name of field to remove
) {

   // Removing a field is simply erasing the field from list
   Fields.erase(FieldName);
   return;
}

//------------------------------------------------------------------------------
// Removes a field from a group with a given name
void FieldGroup::removeFieldFromGroup(
    const std::string &FieldName, // [in] Name of field to remove
    const std::string &GroupName  // [in] Name of group holding field
) {

   // Check that the group exists
   if (!exists(GroupName))
      ABORT_ERROR("Failed to remove field {} from group {}: "
                  "Group does not exist",
                  FieldName, GroupName);

   auto Group = AllGroups[GroupName];
   Group->removeField(FieldName);

   return;
}

//------------------------------------------------------------------------------
// Retrieves a pointer to a field group.
std::shared_ptr<FieldGroup>
FieldGroup::get(const std::string &GroupName // [in] Name of group to retrieve
) {

   // Check to make sure group exists
   if (!exists(GroupName))
      ABORT_ERROR("Failed to retrieve FieldGroup {}: group does not exist.",
                  GroupName);

   // Return the pointer from the list of groups
   return AllGroups[GroupName];
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
   if (!(FieldGroup::exists(GroupName)))
      ABORT_ERROR("Error getting field list from group {}. "
                  "Group does not exist",
                  GroupName);

   // Retrieve group
   std::shared_ptr<FieldGroup> ThisGroup = AllGroups[GroupName];

   // Return list from member function
   return ThisGroup->getFieldList();
}

//------------------------------------------------------------------------------
// Retrieves a full field from a group instance
std::shared_ptr<Field> FieldGroup::getField(
    const std::string &FieldName // [in] Name of field to retrieve
) const {

   // Check to see if the group has field
   if (!hasField(FieldName))
      ABORT_ERROR("Cannot find Field {} in Group {}", FieldName, GrpName);

   // Return the field
   return Field::get(FieldName);
}

//------------------------------------------------------------------------------
// Retrieves a full field from a group with a given name
std::shared_ptr<Field> FieldGroup::getFieldFromGroup(
    const std::string &FieldName, // [in] Name of field to retrieve
    const std::string &GroupName  // [in] Name of group holding field
) {

   // Check to see of group exists
   if (!exists(GroupName))
      ABORT_ERROR("Unable to retrieve Field {} from Group {}. Group not found.",
                  FieldName, GroupName);

   // Retrieve group - getField function performs checks
   std::shared_ptr<FieldGroup> Group = AllGroups[GroupName];
   return Group->getField(FieldName);
}

//------------------------------------------------------------------------------

} // namespace OMEGA
//===----------------------------------------------------------------------===//
