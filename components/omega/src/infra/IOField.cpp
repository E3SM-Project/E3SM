//===-- infra/IOField.cpp - IO field class implementation -------*- C++ -*-===//
//
// \file
// \brief Implements the IO Field class and methods
//
// This file implements the IOField class for OMEGA. IOFields are the
// data and metadata for all fields participating in IO at any point. Each
// module defines which fields are available with the field metadata. Modules
// are also responsible for assigning a correct data pointer for fields that
// they own and for updating that pointer if the location changes (eg if a
// time index changes).
//
//===----------------------------------------------------------------------===//

#include "IOField.h"
#include "DataTypes.h"
#include "Logging.h"
#include "MetaData.h"
#include <map>
#include <memory>

namespace OMEGA {

// Initialize static variables
std::map<std::string, std::shared_ptr<IOField>> IOField::AvailableFields;

//------------------------------------------------------------------------------
// Defines a field so it is available for IO. The caller must
// provide a name that matches a previously defined MetaData field.
// This specific interface defines only the metadata and the actual data
// array must be explicitly attached later.

int IOField::define(
    const std::string &FieldName ///< [in] Name of field to be defined
) {

   int Err = 0; // initialize return code

   // Perform some error checks.
   // If MetaData not defined, exit with an error
   if (!MetaData::has(FieldName)) {
      Err = -1;
      LOG_ERROR("IOField: error defining {} - MetaData does not exist",
                FieldName);
      return Err;
   }

   // If a field of the same name already exists, exit with an error
   if (IOField::isDefined(FieldName)) { // entry found
      Err = -2;
      LOG_ERROR("IOField: error defining {}. Field already exists", FieldName);
      return Err;
   }

   // Create an empty IOField
   auto ThisField = std::make_shared<IOField>();

   // Retrieve and attach the MetaData
   ThisField->MetadataPtr = MetaData::get(FieldName);

   // Add entry to available fields
   IOField::AvailableFields[FieldName] = ThisField;

   return Err;
}

//------------------------------------------------------------------------------
// Checks to see if a field exists with a given name
bool IOField::isDefined(const std::string &FieldName // [in] Name of field
) {
   auto it = IOField::AvailableFields.find(FieldName);
   if (it != IOField::AvailableFields.end()) { // entry found
      return true;
   } else {
      return false;
   }
}

//------------------------------------------------------------------------------
// Retrieves IOField MetaData by name. Returns a pointer to attached metadata.

std::shared_ptr<MetaData>
IOField::getMetaData(const std::string &FieldName ///< [in] name of IOField
) {

   if (IOField::isDefined(FieldName)) { // entry found

      auto ThisField = IOField::AvailableFields[FieldName];
      return ThisField->MetadataPtr;

   } else { // field not found

      LOG_ERROR("IOField: Attempted to get metadata failed, {} does not exist",
                FieldName);
      return nullptr;
   }
}

//------------------------------------------------------------------------------
// Removes a single IOField from the list of available fields

void IOField::erase(const std::string &Name /// Name of IOField to remove
) {
   IOField::AvailableFields.erase(Name);
}

//------------------------------------------------------------------------------
// Removes all IOFields. This must be called before exiting environments

void IOField::clear() { IOField::AvailableFields.clear(); }

} // end namespace OMEGA
//===----------------------------------------------------------------------===//
