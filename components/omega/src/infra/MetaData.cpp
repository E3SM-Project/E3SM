//===-- infra/MetaData.cpp - Implements Omega MetaData functions --*- C++
//-*-===//
//
/// \file
/// \brief implements Omega metadata functions
///
/// This implements Omega metadata initialization.
//
//===----------------------------------------------------------------------===//

#include "MetaData.h"
#include "Logging.h"

namespace OMEGA {

std::map<std::string, std::shared_ptr<MetaData>> MetaData::AllFields;
std::map<std::string, std::shared_ptr<MetaDim>> MetaDim::AllDims;
std::map<std::string, std::shared_ptr<MetaGroup>> MetaGroup::AllGroups;

bool MetaData::has(const std::string Name /// Name of field
) {
   return (AllFields.find(Name) != AllFields.end());
}

/// Create an empty instance by FieldName, primarily for global
/// metadata not assigned to a variable/field. An exception
/// is occured if metadata of the same FieldName already exists
std::shared_ptr<MetaData> MetaData::create(const std::string Name) {

   if (has(Name)) {
      throw std::runtime_error("Failed to create a field instance because '" +
                               Name + "' already exists.");
   }

   auto Data = std::make_shared<MetaData>();

   AllFields[Name] = Data;

   return Data;
}

std::shared_ptr<MetaData> MetaData::create(
    const std::string Name,
    const std::initializer_list<std::pair<std::string, std::any>> &MetaPairs) {

   auto Data = create(Name);

   for (const auto &MetaPair : MetaPairs) {
      Data->addEntry(MetaPair.first, MetaPair.second);
   }

   return Data;
}

/// Create metadata for a variable/field using all
/// required metadata. Note that if input parameters don't exist
/// (eg stdName) or don't make sense (eg min/max or fill) for a
/// given field, empty or 0 entries can be provided. Similarly
/// for scalars, nDims can be 0 and an emtpy MetaDim vector can be
/// supplied.

std::shared_ptr<MetaData> ArrayMetaData::create(
    const std::string Name,        /// Name of var
    const std::string Description, /// Long Name or description
    const std::string Units,       /// Units
    const std::string StdName,     /// CF standard Name
    const std::any ValidMin,       /// Min valid field Value
    const std::any ValidMax,       /// Max valid field Value
    const std::any FillValue,      /// Scalar used for undefined entries
    const int NDims,               /// Number of dimensions
    const std::vector<std::shared_ptr<MetaDim>> Dimensions // MetaDim pointers
) {

   auto Data = MetaData::create(Name);

   Data->addEntry("Description", Description);
   Data->addEntry("Units", Units);
   Data->addEntry("StdName", StdName);
   Data->addEntry("ValidMin", ValidMin);
   Data->addEntry("ValidMax", ValidMax);
   Data->addEntry("FillValue", FillValue);
   Data->addEntry("NDims", NDims);
   Data->addEntry("Dimensions", Dimensions);

   return Data;
}

int MetaData::destroy(const std::string Name /// Name of field to remove
) {
   int RetVal = 0;

   if (!has(Name)) {
      LOG_ERROR("Failed to destroy the field '" + Name + "' because " +
                "the field name does not exist.");
      RetVal = -1;

   } else {
      if (AllFields.erase(Name) != 1) {
         LOG_ERROR("Field, '" + Name + "', is not correctly removed.");
         RetVal = -2;
      }
   }

   return RetVal;
}

std::shared_ptr<MetaData> MetaData::get(const std::string Name /// Name of field
) {
   if (!has(Name)) {
      throw std::runtime_error("Failed to retrieve a field instance because '" +
                               Name + "' does not exist.");
   }

   return AllFields[Name];
}

bool MetaData::hasEntry(const std::string Name // Name of metadata
) {
   return (MetaMap.find(Name) != MetaMap.end());
}

int MetaData::addEntry(const std::string Name, /// Name of new metadata to add
                       const std::any Value    /// Value of new metadata to add
) {
   int RetVal = 0;

   if (hasEntry(Name)) {
      LOG_ERROR("Failed to add the metadata '" + Name + "' because " +
                "the field already has the metadata.");

      RetVal = -1;

   } else {
      MetaMap[Name] = Value;
   }

   return RetVal;
}

int MetaData::removeEntry(const std::string Name /// Name of metadata to remove
) {

   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Failed to remove the metadata '" + Name + "' because " +
                "the metadata name does not exist.");
      RetVal = -1;

   } else {
      if (MetaMap.erase(Name) != 1) {
         LOG_ERROR("Metadata, '" + Name + "', is not correctly removed.");
         RetVal = -2;
      }
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       I4 &Value               /// I4 Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<I4>(MetaMap[Name]);
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       I8 &Value               /// I8 Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<I8>(MetaMap[Name]);
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       R4 &Value               /// R4 Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<R4>(MetaMap[Name]);
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       R8 &Value               /// R8 Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<R8>(MetaMap[Name]);
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       bool &Value             /// Bool Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<bool>(MetaMap[Name]);
   }

   return RetVal;
}

int MetaData::getEntry(const std::string Name, /// Name of metadata to get
                       std::string &Value      /// String Value of metadata
) {
   int RetVal = 0;

   if (!hasEntry(Name)) {
      LOG_ERROR("Metadata '" + Name + "' does not exist.");
      RetVal = -1;

   } else {
      Value = std::any_cast<std::string>(MetaMap[Name]);
   }

   return RetVal;
}

std::map<std::string, std::any> *MetaData::getAllEntries() { return &MetaMap; }

bool MetaDim::has(const std::string Name /// Name of dimension
) {
   return (AllDims.find(Name) != AllDims.end());
}

std::shared_ptr<MetaDim>
MetaDim::create(const std::string Name, // Name of dimension
                const I4 Length         // Length of dimension
) {
   if (has(Name)) {
      throw std::runtime_error("Failed to create a MetaDim instance because '" +
                               Name + "' already exists.");
   }

   auto Dim    = std::make_shared<MetaDim>();
   Dim->Length = Length;

   AllDims[Name] = Dim;

   return Dim;
}

int MetaDim::destroy(const std::string Name // Name of dimension
) {
   int RetVal = 0;

   if (!has(Name)) {
      LOG_ERROR("Failed to destroy the dimension '" + Name + "' because " +
                "the dimension name does not exist.");
      RetVal = -1;

   } else {
      if (AllDims.erase(Name) != 1) {
         LOG_ERROR("Dimension, '" + Name + "', was not successfully removed.");
         RetVal = -1;
      }
   }

   return RetVal;
}

std::shared_ptr<MetaDim>
MetaDim::get(const std::string Name /// Name of dimension
) {
   if (!has(Name)) {
      throw std::runtime_error(
          "Failed to retrieve a MetaDim instance because '" + Name +
          "' does not exist.");
   }

   return AllDims[Name];
}

int MetaDim::getLength(I4 &Length // Length of dimension
) {
   int RetVal = 0;

   Length = this->Length;

   return RetVal;
}

bool MetaGroup::has(const std::string Name /// Name of group
) {
   return (AllGroups.find(Name) != AllGroups.end());
}

std::shared_ptr<MetaGroup>
MetaGroup::create(const std::string Name // Name of group
) {
   if (has(Name)) {
      return get(Name);

   } else {
      auto Group = std::make_shared<MetaGroup>();

      AllGroups[Name] = Group;

      return Group;
   }
}

int MetaGroup::destroy(const std::string Name // Name of group
) {
   int RetVal = 0;

   if (!has(Name)) {
      LOG_ERROR("Failed to destroy the group '" + Name + "' because " +
                "the group name does not exist.");
      RetVal = -1;

   } else {
      if (AllGroups.erase(Name) != 1) {
         LOG_ERROR("Group, '" + Name + "', is not correctly removed.");
         RetVal = -2;
      }
   }

   return RetVal;
}

std::shared_ptr<MetaGroup>
MetaGroup::get(const std::string Name /// Name of dimension
) {
   if (!has(Name)) {
      throw std::runtime_error(
          "Failed to retrieve a MetaGroup instance because '" + Name +
          "' does not exist.");
   }

   return AllGroups[Name];
}

bool MetaGroup::hasField(const std::string FieldName /// Name of field
) {
   return (Fields.find(FieldName) != Fields.end());
}

int MetaGroup::addField(const std::string FieldName // Name of field to add
) {
   int RetVal = 0;

   if (hasField(FieldName)) {
      LOG_ERROR("Failed to add the field '" + FieldName + "' because " +
                "the group already has the field.");
      RetVal = -1; // The field already exists in the group.

   } else if (MetaData::has(FieldName)) {
      Fields[FieldName] = MetaData::get(FieldName);

   } else {
      LOG_ERROR("Failed to add the field '" + FieldName + "' because " +
                "the field does not exist.");
      RetVal = -2; // The field name does not exist.
   }

   return RetVal;
}

std::shared_ptr<MetaData>
MetaGroup::getField(const std::string FieldName /// Name of field to add
) {
   if (!hasField(FieldName)) {
      throw std::runtime_error("Failed to get a MetaData instance because "
                               "the group does not have '" +
                               FieldName + "'.");
   }

   return Fields[FieldName];
}

int MetaGroup::removeField(
    const std::string FieldName // Name of field to remove
) {
   int RetVal = 0;

   if (!hasField(FieldName)) {
      LOG_ERROR("Failed to get remove MetaData instance because "
                "the group does not have '" +
                FieldName + "'.");
      RetVal = -1;

   } else {
      if (Fields.erase(FieldName) != 1) {
         LOG_ERROR("Field, '" + FieldName + "', is not correctly removed.");
         RetVal = -2;
      }
   }

   return RetVal;
}

std::map<std::string, std::shared_ptr<MetaData>> *MetaGroup::getAllFields() {
   return &Fields;
}

} // Namespace OMEGA
