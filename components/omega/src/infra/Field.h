#ifndef OMEGA_FIELD_H
#define OMEGA_FIELD_H
//===-- infra/Field.h - OMEGA field class -----------------------*- C++ -*-===//
//
/// \file
/// \brief Defines Field class and methods
///
/// This header defines classes and methods for Fields. Fields are the
/// data and metadata for all fields that may be a part of an IO stream or any
/// other module where it might be important to keep the data and metadata
/// together. Fields are defined/created withing the modules that own them. It
/// is also possible to define field groups to provide a shortcut for fields
/// that are commonly used together. Once a field is defined/created, a data
/// array can be attached detached or swapped to hold the latest data values
/// (esp for fields with multiple time levels).
///
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Dimension.h"
#include "Logging.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include <any>
#include <map>
#include <memory>
#include <set>

namespace OMEGA {

// Define the metadata data type
using Metadata = std::map<std::string, std::any>;

/// Field name to use for global code metadata
static const std::string CodeMeta{"code"}; ///< name for code metadata
/// Field name to use for global simulation metadata
static const std::string SimMeta{"simulation"}; ///< name for sim metatdata

//------------------------------------------------------------------------------
/// The Field class manages all metadata and attached data for OMEGA fields
/// used in IO and any other situation where both data and metadata are needed
class Field {

 private:
   /// Store and maintain all defined fields
   static std::map<std::string, std::shared_ptr<Field>> AllFields;

   /// Field name
   std::string FldName;

   /// Metadata (name, value) pairs for descriptive metadata
   Metadata FieldMeta;

   /// Number of dimensions for the field (0 if scalar field or global metadata)
   int NDims;

   /// Dimension names for retrieval of dimension info
   /// These must be in the same index order as the stored data
   std::vector<std::string> DimNames;

   /// Data type for field data
   ArrayDataType DataType;

   /// Location of data
   ArrayMemLoc MemLoc;

   /// Flag for whether this is a time-dependent field that needs the
   /// Unlimited time dimension added during IO
   bool TimeDependent;

   /// Flag for whether this is a field that is distributed across tasks
   /// or whether it is entirely local
   bool Distributed;

   /// Data attached to this field. This will be a pointer to the Kokkos
   /// array holding the data. We use a void pointer to manage all the
   /// various types and cast to the appropriate type when needed.
   std::shared_ptr<void> DataArray;

 public:
   //---------------------------------------------------------------------------
   // Initialization
   //---------------------------------------------------------------------------
   /// Initializes the fields for global code and simulation metadata
   /// It also initializes the unlimited time dimension needed by most fields
   static int init(const Clock *ModelClock ///< [in] the default model clock
   );

   //---------------------------------------------------------------------------
   // Create/destroy/query fields
   //---------------------------------------------------------------------------
   /// Checks to see if a field has been defined
   static bool exists(const std::string &FieldName ///< [in] name of field
   );

   //---------------------------------------------------------------------------
   /// Creates a field with standard metadata. This is the preferred
   /// interface for most fields in Omega. It enforces a list of required
   /// metadata. Note that if input parameters don't exist
   /// (eg stdName) or don't make sense (eg min/max or fill) for a
   /// given field, empty or 0 entries can be provided. Actual field data is
   /// attached in a separate call and additional metadata can be added later.
   static std::shared_ptr<Field>
   create(const std::string &FieldName,   ///< [in] Name of variable/field
          const std::string &Description, ///< [in] long Name or description
          const std::string &Units,       ///< [in] units
          const std::string &StdName,     ///< [in] CF standard Name
          const std::any ValidMin,        ///< [in] min valid field value
          const std::any ValidMax,        ///< [in] max valid field value
          const std::any FillValue,       ///< [in] scalar for undefined entries
          const int NumDims,              ///< [in] number of dimensions
          const std::vector<std::string> &Dimensions, ///< [in] dim names
          const bool TimeDependent = true ///< [in] opt flag for unlim time
   );

   //---------------------------------------------------------------------------
   /// Creates an empty field with a given name. This interface should only
   /// be used for collections of metadata with no attached data array.
   /// Two standard fields called CodeMeta and SimMeta are defined on init
   /// to hold global metadata not attached to fields.
   static std::shared_ptr<Field>
   create(const std::string &FieldName ///< [in] Name of field
   );

   //---------------------------------------------------------------------------
   /// Removes a single Field from the list of available fields
   /// That process also decrements the reference counters for the
   /// shared pointers and removes them if those counters reach 0.
   static int destroy(const std::string &FieldName /// Name of Field to remove
   );

   //---------------------------------------------------------------------------
   /// Removes all Fields. This must be called before exiting environments
   /// This removes all fields from the map structure and also
   /// decrements the reference counter for the shared pointers,
   /// removing them if the count has reached 0.
   static void clear();

   //---------------------------------------------------------------------------
   // Retrieval functions
   //---------------------------------------------------------------------------

   //---------------------------------------------------------------------------
   /// Retrieve pointer to a full field by name
   static std::shared_ptr<Field>
   get(const std::string &FieldName ///< [in] name of field to retrieve
   );

   //---------------------------------------------------------------------------
   /// Retrieve field name
   std::string getName() const;

   //---------------------------------------------------------------------------
   // Retrieve data type of field
   /// Determine type of a given field from instance
   ArrayDataType getType() const;

   /// Determine type of a given field by name
   static ArrayDataType
   getFieldType(const std::string &FieldName ///< [in] name of field
   );

   //---------------------------------------------------------------------------
   // Query location of field data
   /// Determine location of data from instance
   ArrayMemLoc getMemoryLocation() const;

   /// Determine location of data by field name
   static ArrayMemLoc
   getFieldMemoryLocation(const std::string &FieldName ///< [in] name of field
   );

   /// Query whether field is located on host from instance
   bool isOnHost() const;

   /// Query whether field is located on host by field name
   static bool
   isFieldOnHost(const std::string &FieldName ///< [in] name of field
   );

   //---------------------------------------------------------------------------
   // Query for dimension info
   //---------------------------------------------------------------------------
   /// Returns the number of dimensions for the field
   int getNumDims() const;

   /// Returns a vector of dimension names associated with each dimension
   /// of an array field. Returns an error code.
   int getDimNames(
       std::vector<std::string> &Dimensions ///< [out] list of dimensions
   ) const;

   //---------------------------------------------------------------------------
   // Query for other properties
   /// Determine whether this is a time-dependent field that requires the
   /// unlimited time dimension for IO
   bool isTimeDependent() const;

   /// Determine whether this is a distributed field or whether it is entirely
   /// local. This is needed to determine whether IO uses parallel read/write
   /// or an undistributed read/write.
   bool isDistributed() const;

   //---------------------------------------------------------------------------
   // Metadata functions
   //---------------------------------------------------------------------------
   /// Retrieves all Metadata for a Field instance. Returns a
   /// shared pointer to the map structure with metadata.
   std::shared_ptr<Metadata> getAllMetadata() const;

   /// Retrieves all Metadata for a Field given the field name. Returns a
   /// shared pointer to the map structure with metadata.
   static std::shared_ptr<Metadata>
   getFieldMetadata(const std::string &FieldName ///< [in] name of field
   );

   /// Checks for the existence of a metadata entry with the given name
   bool hasMetadata(const std::string &MetaName ///< [in] Name of metadata entry
   ) const;

   /// Adds a metadata entry with the (name,value) pair
   int addMetadata(const std::string &MetaName, ///< [in] Name of new metadata
                   const std::any Value         ///< [in] Value of new metadata
   );

   /// Adds multiple metadata entries with a list of (name,value) pairs
   int addMetadata(const std::initializer_list<std::pair<std::string, std::any>>
                       &MetaPairs);

   /// Updates a metadata entry with a new value
   int updateMetadata(const std::string &MetaName, ///< [in] Name of metadata
                      const std::any Value         ///< [in] Value of metadata
   );

   /// Removes a metadata entry with the given name
   int
   removeMetadata(const std::string &MetaName ///< [in] Name of entry to remove
   );

   /// Removes all metadata entries from a field
   void removeAllMetadata();

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to an I4 type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   I4 &Value                ///< [out] I4 Value of metadata
   );

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to an I8 type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   I8 &Value                ///< [out] I8 Value of metadata
   );

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to an R4 type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   R4 &Value                ///< [out]  R4 Value of metadata
   );

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to an R8 type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   R8 &Value                ///< [out] R8 Value of metadata
   );

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to a bool type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   bool &Value              ///< [out] bool Value of metadata
   );

   /// Retrieves the value of the metadata associated with a given name
   /// This specific version of the overloaded interface coerces the value
   /// to a string type
   int getMetadata(const std::string &Name, ///< [in] Name of metadata to get
                   std::string &Value       ///< [out] string Value of metadata
   );

   //---------------------------------------------------------------------------
   // Attach/detach data
   // Template functions must have implementation in header files
   //---------------------------------------------------------------------------
   /// Attaches an array of data to an existing Field instance. If a data array
   /// needs to be updated, calling the attach routine with new data
   /// will replace the previously attached data array. This function is
   /// templated based on the array data type so a template argument with
   /// the proper array data type (eg <Array2DR4>) must be supplied.
   template <typename T>
   int attachData(const T &InDataArray ///< [in] Array with data to attach
   ) {
      static_assert(isKokkosArray<T>,
                    "attachData requires Kokkos array as input");
      int Err = 0; // initialize return code

      // Attach the data array - this is a shallow copy
      DataArray = std::make_shared<T>(InDataArray);

      // Determine type and location
      DataType = Impl::checkArrayType<T>();
      MemLoc   = Impl::findArrayMemLoc<T>();

      return Err;
   };

   //---------------------------------------------------------------------------
   /// Attaches an array of data to an existing Field by name. If a data array
   /// needs to be updated, calling the attach routine with new data
   /// will replace the previously attached data array. This function is
   /// templated based on the array data type so a template argument with
   /// the proper array data type (eg <Array2DR4>) must be supplied.
   template <typename T>
   static int
   attachFieldData(const std::string &FieldName, ///< [in] Name of Field
                   const T &InDataArray ///< [in] Array with data to attach
   ) {
      static_assert(isKokkosArray<T>,
                    "attachFieldData requires Kokkos array as input");

      int Err = 0; // initialize return code

      // Check to make sure field exists
      if (exists(FieldName)) { // entry found
         // Retrieve the field
         auto ThisField = AllFields[FieldName];
         // Attach the data array
         Err = ThisField->attachData<T>(InDataArray);

      } else { // field has not yet been defined
         Err = -1;
         LOG_ERROR("Field: error attaching data to {}. Field not defined",
                   FieldName);
      }
      return Err;
   };

   //---------------------------------------------------------------------------
   /// Retrieves Field data from a Field instance. Because all data
   /// arrays in OMEGA are Kokkos arrays, this is a shallow copy of the
   /// attached data array. This is a templated function on the supported
   /// OMEGA array types so a template argument with the proper type must
   /// also be supplied.
   template <typename T> T getDataArray() {
      static_assert(
          isKokkosArray<T>,
          "getDataArray requires Kokkos array as its template argument");

      // Check to make sure data is attached
      if (DataArray != nullptr) { // data is attached

         // Convert the data pointer to the appropriate type and
         // return data array dereferenced from the pointer
         return *(std::static_pointer_cast<T>(DataArray));

      } else { // data has not yet been attached

         LOG_ERROR("Field: Attempt to get data failed from field {}. "
                   "Data array has not been attached",
                   FldName);
         // return an empty data object
         T Data;
         return Data;
      }
   };

   //---------------------------------------------------------------------------
   /// Retrieves Field data array given the field name. Because all data
   /// arrays in OMEGA are Kokkos arrays, this is a shallow copy of the
   /// attached data array. This is a templated function on the supported
   /// OMEGA array types so a template argument with the proper type must
   /// also be supplied.
   template <typename T>
   static T
   getFieldDataArray(const std::string &FieldName ///< [in] name of Field
   ) {

      T Data; // Set up empty array

      // Check to see if field is defined
      if (!exists(FieldName)) { // no entry found return error
         LOG_ERROR("Field: Attempted to get data failed, {} does not exist",
                   FieldName);
         return Data;
      }

      // Retrieve the field by name
      auto ThisField = AllFields[FieldName];

      // Retrieve the data
      Data = ThisField->getDataArray<T>();
      return Data;
   };

   //---------------------------------------------------------------------------
   // Field Group is a friend class so it can access field list
   friend class FieldGroup;

}; // end class Field

//------------------------------------------------------------------------------
/// The FieldGroup class allows fields to be grouped together and referred
/// to by a single name. This helps to reduce the length and complexity of
/// input lists and contents lists.
class FieldGroup {

 private:
   /// Store and maintain all defined groups in map container
   static std::map<std::string, std::shared_ptr<FieldGroup>> AllGroups;

   /// Name of group
   std::string GrpName;

   /// List of fields in group - only names are needed
   std::set<std::string> Fields;

 public:
   /// Creates an empty field group with a given name
   static std::shared_ptr<FieldGroup>
   create(const std::string &GroupName ///< [in] Name of group to create
   );

   /// Removes a field group
   static int
   destroy(const std::string &GroupName ///< [in] Name of group to destroy
   );

   /// Removes all defined field groups
   static void clear();

   /// Determines whether a group of a given name exists
   static bool exists(const std::string &GroupName ///< [in] Name of group
   );

   /// Determines whether a field of a given name exists in the group instance
   bool hasField(const std::string &FieldName ///< [in] Name of field
   ) const;

   /// Determines whether a field of a given name exists in the group
   /// based on group name.
   static bool
   isFieldInGroup(const std::string &FieldName, ///< [in] Name of field
                  const std::string &GroupName  ///< [in] Name of group
   );

   /// Adds a field to the group instance
   int addField(const std::string &FieldName ///< [in] Name of field to add
   );

   /// Adds a field to a group based on group name
   static int
   addFieldToGroup(const std::string &FieldName, ///< [in] Name of field to add
                   const std::string &GroupName  ///< [in] Name of group
   );

   /// Removes a field from a group instance
   int removeField(const std::string &FieldName ///< [in] Name of field removed
   );

   /// Removes a field from a group with a given name
   static int removeFieldFromGroup(
       const std::string &FieldName, ///< [in] Name of field to remove
       const std::string &GroupName  ///< [in] Name of group holding field
   );

   /// Retrieves a pointer to a field group
   static std::shared_ptr<FieldGroup>
   get(const std::string &GroupName ///< [in] Name of group to retrieve
   );

   /// Returns the list of fields in a group instance
   std::set<std::string> getFieldList() const;

   /// Returns the list of fields in a group based on group name
   static std::set<std::string> getFieldListFromGroup(
       const std::string &GroupName ///< [in] Name of group holding fields
   );

   /// Retrieves a full field from a group instance
   std::shared_ptr<Field>
   getField(const std::string &FieldName ///< [in] Name of field to retrieve
   ) const;

   /// Retrieves a full field from a group with a given name
   static std::shared_ptr<Field> getFieldFromGroup(
       const std::string &FieldName, ///< [in] Name of field to retrieve
       const std::string &GroupName  ///< [in] Name of group holding field
   );

}; // end class FieldGroup

//------------------------------------------------------------------------------

} // namespace OMEGA
//===----------------------------------------------------------------------===//
#endif // OMEGA_FIELD_H
