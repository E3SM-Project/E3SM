#ifndef OMEGA_METADATA_H
#define OMEGA_METADATA_H
//===-- infra/MetaData.h - metadata classes --*- C++ -*-===//
//
/// \file
/// \brief Defines metadata classes
///
/// This header defines classes for metadata.
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include <any>
#include <map>

namespace OMEGA {

static const std::string CodeMeta{"code"};
static const std::string SimMeta{"simulation"};

class MetaDim {

 private:
   /// The length of the dimension. Use 0 for unlimited length.
   I4 Length;

   /// Store and maintain all defined dimensions
   static std::map<std::string, std::shared_ptr<MetaDim>> AllDims;

 public:
   static std::shared_ptr<MetaDim>
   create(const std::string Name, /// Name of dimension
          const I4 Length         /// length of dimension
   );

   static int destroy(const std::string Name /// Name of dimension
   );

   static std::shared_ptr<MetaDim>
   get(const std::string Name /// Name of dimension
   );

   static bool has(const std::string Name /// Name of dimension
   );

   int getLength(I4 &Length /// length of dimension
   );
};

class MetaData {
 protected:
   /// metadata stored in map
   std::map<std::string, std::any> MetaMap;

   /// Store and maintain all defined metadata in this vector
   static std::map<std::string, std::shared_ptr<MetaData>> AllFields;

 public:
   static std::shared_ptr<MetaData>
   create(const std::string Name /// Name of field
   );

   static std::shared_ptr<MetaData>
   create(const std::string Name,
          const std::initializer_list<std::pair<std::string, std::any>>
              &MetaPairs);

   static int destroy(const std::string Name /// Name of field to remove
   );

   static std::shared_ptr<MetaData> get(const std::string Name /// Name of field
   );

   static bool has(const std::string Name /// Name of field
   );

   bool hasEntry(const std::string MetaName /// Name of metadata
   );

   int addEntry(const std::string MetaName, /// Name of new metadata to add
                const std::any Value        /// Value of new metadata to add
   );

   int removeEntry(const std::string MetaName /// Name of metadata to remove
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                I4 &Value                   /// I4 Value of metadata
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                I8 &Value                   /// I8 Value of metadata
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                R4 &Value                   ///  R4 Value of metadata
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                R8 &Value                   /// R8 Value of metadata
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                bool &Value                 /// bool Value of metadata
   );

   int getEntry(const std::string MetaName, /// Name of metadata to get
                std::string &Value          /// string Value of metadata
   );

   /// returns the pointer to the metadata map
   std::map<std::string, std::any> *getAllEntries();
};

class ArrayMetaData : public MetaData {

 public:
   static std::shared_ptr<MetaData>
   create(const std::string Name,        /// Name of var
          const std::string Description, /// long Name or description
          const std::string Units,       /// units
          const std::string StdName,     /// CF standard Name
          const std::any ValidMin,       /// min valid field value
          const std::any ValidMax,       /// max valid field value
          const std::any FillValue,      /// scalar used for undefined entries
          const int NDims,               /// number of dimensions
          const std::vector<std::shared_ptr<MetaDim>> Dimensions // dim pointers
   );
};

class MetaGroup {

 private:
   /// map of fields in group
   std::map<std::string, std::shared_ptr<MetaData>> Fields;

   /// Store and maintain all defined groups
   static std::map<std::string, std::shared_ptr<MetaGroup>> AllGroups;

 public:
   static std::shared_ptr<MetaGroup>
   create(const std::string Name /// Name of group
   );

   static int destroy(const std::string Name /// Name of group
   );

   static std::shared_ptr<MetaGroup>
   get(const std::string Name /// Name of group
   );

   static bool has(const std::string Name /// Name of group
   );

   bool hasField(const std::string FieldName /// Name of field to add
   );

   int addField(const std::string FieldName /// Name of field to add
   );

   std::shared_ptr<MetaData>
   getField(const std::string FieldName /// Name of field to add
   );

   int removeField(const std::string FieldName /// Name of field to remove
   );

   /// returns the pointer to the metagroup map
   std::map<std::string, std::shared_ptr<MetaData>> *getAllFields();
};

} // namespace OMEGA

#endif // OMEGA_METADATA_H
