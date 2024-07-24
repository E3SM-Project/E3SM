(omega-user-metadata)=

## Metadata

### Introduction to Metadata in Omega
- **Purpose**: Enhances Omega and E3SM data file management for streamlined analysis, archiving, and provenance.

### Working with Metadata
- **Inclusion of Metadata**: Utilize by including `MetaData.h` in the source
  code.
  ```
  #include "MetaData.h"
  ```

### Creating and Destroying MetaData Instances
- **Creation Methods**:
  1. Using `create` with a name string.
     ```
     auto Data1 = MetaData::create("MyVar1");
     ```
  2. Creating with name-value pairs using `std::make_pair`.
     ```
     auto Data2 = MetaData::create(
                    "MyVar2",
                    {
                        std::make_pair("Meta1", 1),
                        std::make_pair("Meta2", true),
                        std::make_pair("Meta3", "MyValue3"),
                    });
     ```
  3. For array-type variables, specific key-value pairs are mandatory.
     ```
     auto Data3 = ArrayMetaData::create(
                    "MyVar3",
                    "Description", /// long Name or description
                    "Units",       /// units
                    "StdName",     /// CF standard Name
                    std::numeric_limits<int>::min(), /// min valid value
                    std::numeric_limits<int>::max(), /// max valid value
                    FILL_VALUE,    /// scalar used for undefined entries
                    1,             /// number of dimensions
                    Dimensions     /// dim pointers
                    );
     ```
- **Destruction**: Utilize `destroy` function with the metadata name.
  ```
  int ret = MetaData::destroy("MyVar3");
  ```

### Managing Metadata in Instances
- **Adding Metadata**: Use `addEntry` function.
  ```
  const R8 AValue = 2.0;
  int ret = Data1->addEntry("NewMeta", AValue);
  ```
  The following data types are allowed as metadata values: I4, I8, R4, R8, std::string, and bool.
- **Removing Metadata**: Employ `removeEntry` function.
  ```
  int ret = Data1->removeEntry("NewMeta");
  ```
- **Retrieving Metadata**: Use `getEntry` function.
  ```
  R8 R8Value;
  int ret = Data1->getEntry("NewMeta", R8Value);
  ```

### Handling MetaData Instances
- **Retrieval**: Use `get` static function.
  ```
  auto Data4 = MetaData::get("MyVar2");
  ```

### Working with MetaDim Instances
- **Creating and Destroying**: Similar to MetaData instances.
  ```
  const I4 DimLength = 3;
  auto Dim1 = MetaDim::create("MyDim1", DimLength);
  int ret   = MetaDim::destroy("MyDim1");
  ```
- **Retrieving Dimensions**: Use `getLength` function.
  ```
  I4 DimLength;
  Dim1->getLength(DimLength);
  ```

### Managing MetaGroup Instances
- **Purpose**: For managing lengthy model configuration lists.
- **Creation**: Similar to MetaData and MetaDim classes.
  ```
  auto Group1 = MetaGroup::create("MyGroup1");
  int ret     = MetaGroup::destroy("MyGroup1");
  ```
- **Field Management**: Add, retrieve, or remove fields using `addField`,
`getField`, and `removeField`.
  ```
  int ret       = Group1->addField(FieldName);
  auto Data1    = Group1->getField(FieldName):
  int ret       = Group1->removeField(FieldName):
  ```

### Utility Functions
- **`has` Static Method**: Checks existence of an instance.
  ```
  bool exists = MetaData::has("MyVar");
  ```
- **`get` Static Method**: Retrieves an instance.
  ```
  auto Instance = MetaData::get("MyVar");
  ```

Refer to the Metadata Design and Developer Document for comprehensive details.
