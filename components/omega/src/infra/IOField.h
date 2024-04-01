#ifndef OMEGA_IOFIELD_H
#define OMEGA_IOFIELD_H
//===-- infra/IOField.h - IO field class ------------------------*- C++ -*-===//
//
/// \file
/// \brief Defines IO Field class and methods
///
/// This header defines classes and methods for IO Fields. IOFields are the
/// data and metadata for all fields participating in IO at any point. Each
/// module defines which fields are available with the field metadata. Modules
/// are also responsible for assigning a correct data pointer for fields that
/// they own and for updating that pointer if the location changes (eg if a
/// time index changes).
///
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Logging.h"
#include "MetaData.h"
#include <map>
#include <memory>

namespace OMEGA {

class IOField {

 private:
   /// Store and maintain all defined fields
   static std::map<std::string, std::shared_ptr<IOField>> AvailableFields;

   /// Metadata associated with this field
   std::shared_ptr<MetaData> MetadataPtr;

   /// Data assigned to this field. We use a void pointer to manage the
   /// many different Array types. These are cast to the appropriate type
   /// on retrieval.
   std::shared_ptr<void> Data;

 public:
   //---------------------------------------------------------------------------
   /// Checks to see if a field is defined
   static bool isDefined(const std::string &FieldName ///< [in] name of field
   );

   //---------------------------------------------------------------------------
   /// Defines a field so it is available for IO. The name of the
   /// field must match a previously defined MetaData field.
   /// This routine simply attaches the previously defined metadata
   /// and add the field to the list of available fields. The actual data
   /// array is attached in a separate call.
   static int define(const std::string &FieldName ///< [in] Name of field
   );

   //---------------------------------------------------------------------------
   /// Retrieves IOField MetaData by name. Returns a shared pointer to
   /// the attached metadata.
   static std::shared_ptr<MetaData>
   getMetaData(const std::string &FieldName ///< [in] name of IOField
   );

   // Template functions must have implementation in header files
   //---------------------------------------------------------------------------
   /// Attaches an array of data to an existing IOField. If a data array
   /// needs to be updated, calling the attach routine with new data
   /// will replaced the previously attached data array. This function is
   /// templated based on the array data type so a template argument with
   /// the proper array data type (eg <Array2DR4>) must be supplied.
   template <typename T>
   static int attachData(const std::string &FieldName, ///< [in] Name of IOField
                         const T &DataArray ///< [in] Array with data to attach
   ) {

      int Err = 0; // initialize return code

      // Check to make sure field exists
      if (isDefined(FieldName)) { // entry found
         // Retrieve the IO field
         auto ThisField = AvailableFields[FieldName];
         // Attach the data array
         ThisField->Data = std::make_shared<T>(DataArray);

      } else { // field has not yet been defined
         Err = -1;
         LOG_ERROR("IOField: error attaching data to {}. Field not defined",
                   FieldName);
      }
      return Err;
   };

   //---------------------------------------------------------------------------
   /// Retrieves IOField data array given the field name. Because all data
   /// arrays in OMEGA are Kokkos arrays, this is a shallow copy of the
   /// attached data array. This is a templated function on the supported
   /// OMEGA array types so a template argument with the proper type must
   /// also be supplied.
   template <typename T>
   static T getData(const std::string &FieldName ///< [in] name of IOField
   ) {

      // Check to see if field is defined
      if (!isDefined(FieldName)) { // no entry found return error
         LOG_ERROR("IOField: Attempted to get data failed, {} does not exist",
                   FieldName);
         // return an empty data object
         T Data;
         return Data;
      }

      // Retrieve the field by name
      auto ThisField = AvailableFields[FieldName];

      // Check to make sure data is attached
      if (ThisField->Data != nullptr) { // data is attached

         // Convert the data pointer to the appropriate type and
         // return data array dereferenced from the pointer
         return *(std::static_pointer_cast<T>(ThisField->Data));

      } else { // data has not yet been attached

         LOG_ERROR("IOField: Attempted to get data failed from field {}. "
                   "Data array has not been attached",
                   FieldName);
         // return an empty data object
         T Data;
         return Data;
      }
   };

   //---------------------------------------------------------------------------
   /// Removes a single IOField from the list of available fields
   /// That process also decrements the reference counters for the
   /// shared pointers and removes them if those counters reach 0.
   static void erase(const std::string &FieldName /// Name of IOField to remove
   );

   //---------------------------------------------------------------------------
   /// Removes all IOFields. This must be called before exiting environments
   /// This removes all fields from the map structure and also
   /// decrements the reference counter for the shared pointers,
   /// removing them if the count has reached 0.
   static void clear();

}; // end class IOField

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_IOFIELD_H
