(omega-dev-field)=

## Fields and Metadata (Field)

Omega includes a Field class that registers available fields and associated
metadata for use in IO (or any other part of Omega that may require the field
metadata). The module owning the field will define the Field together
with some required metadata for compliance with the Climate and Forecast
[(CF) metadata conventios.](http://cfconventions.org/). Groups of Fields
can also be defined to provide a short cut to groups that are commonly used
together, like the model state and tracer groups. Two special Fields, CodeMeta
and SimMeta (with default names "code" and "simulation", respectively), are
defined on initialization to store global metadata or scalars needed for
provenance or for restarting. The arrays containing Field data can be attached,
retrieved or updated as needed, typically during IO. The Field class is not
meant to be used during computations - the native Omega array types without
metadata are more appropriate for that context.

To use the Field class, the Field header file must be included and as part of
the application initialization, there should be a call to the init method:
Fields initialized with the init method:
```c++
int Err = Field::init();
```
which primarily defines the CodeMeta and SimMeta fields for later use.
For array fields, the appropriate Dimensions must be defined. The default
dimensions (eg NCells, NEdges, NVertLevels) will be defined by the relevant
Mesh initialization and should be done before any Fields are defined.
See {ref}`omega-dev-dimension`.

Fields are created with standard metadata using
```c++
   std::shared_ptr<Field> MyField =
   Field::create(FieldName,   ///< [in] Name of variable/field (string)
                 Description, ///< [in] long Name or description (string)
                 Units,       ///< [in] units (string)
                 StdName,     ///< [in] CF standard Name (string)
                 ValidMin,    ///< [in] min valid field value (same type as
                 ValidMax,    ///< [in] max valid field value  field data)
                 FillValue,   ///< [in] scalar used for undefined entries
                 NumDims,     ///< [in] number of dimensions (int)
                 Dimensions,  ///< [in] dim names (vector of strings)
                 InTimeDependent ///< [in] (opt, default true) if time varying
   );
```
This interface enforces a list of required metadata. If a CF standard name does
not exist, an empty string can be provided. This is uncommon for most fields
since the CF conventions maintain a fairly complete list, but can be the case
for some intermediate calculations or unique analyses. If there is no
restriction on valid range, an appropriately large range should be provided for
the data type. Similarly, if a FillValue is not being used, a very unique
number should be supplied to prevent accidentally treating valid data as a
FillValue. The optional TimeDependent argument can be omitted and is assumed
to be true by default. Fields with this attribute will be output with the
unlimited time dimension added. Time should not be added explicitly in the
dimension list since it will be added during I/O. Fields that do not change
with time should include this argument with the value false so that the time
dimension is not added. Actual field data stored in an array is attached in a
separate call as described below. Scalar fields can be added by setting the
NumDims to zero (the Dimensions vector is ignored but an empty vector
must still be supplied in the argument list). Scalar data is attached using
a 1D array with size 1. Fields without a data array can be created with:
```c++
   std::shared_ptr<Field> MyField =
   Field::create(FieldName ///< [in] Name of field
   );
```
but this interface should generally be not be used and global metadata
should be added to the existing CodeMeta or SimMeta fields for this
purpose.

Additional metadata in the form of a name-value pair can be added using:
```c++
   int Err = MyField->addMetadata(
                         MetaName, // [in] Name of new metadata (string)
                         Value     // [in] Value of new metadata
   );
```
where Value can be any supported data type (I4, I8, R4, R8, bool, string).
Multiple pairs can be added in a single call using:
``` c++
   Err = SimField->addMetadata(
         {std::make_pair("Name1", Val1),
          std::make_pair("Name2", Val2),
          std::make_pair("Name3", Val3),
          std::make_pair("Name4", Val4)});
```
where any number of pairs can be included and the values can be any supported
type.

As mentioned above, the actual data array is attached in a separate call
using a templated form. If the pointer to the Field is available, use the
member function:
```c++
   int Err = MyField->attachData<ArrayType>(InDataArray);
   // for example:
   int Err = MyField->attachData<HostArray1DI4>(CellID);
   int Err = MyField->attachData<Array2DR8>(NormalVelocity);
```
where ArrayType is one of the supported array types (Array1DI4, etc. or
HostArray1DI4, etc.). If the Field pointer has not been retrieved, an interface
is provided using the field name:
```c++
   int Err = Field::attachFieldData<ArrayType>(
                    FieldName,   // [in] Name of Field (string)
                    InDataArray  // [in] Array with data to attach
                    );
```
Note that the data is assumed to reside in only one location so if a mirror
array exists (eg if replicated on host and device), a separate Field may be
needed. However, it is is better to define only one location and allow the
IO or other modules to determine whether a transfer of data or mirror is needed.
If the data resides in an array whose location does not change (ie the
pointer always points to a fixed location), the attach can be performed when
the Field and array have been created and any updates to the data will be
captured correctly. If the location of the data changes (eg the time
level changes and the pointer points to a different time slice), the data must
be updated by calling the attach routine to replace the pointer to the new
location. It is up to the developer to insert the appropriate call to reattach
the data. As mentioned previously, scalar data should be attached using the
appropriate 1D HostArray with a size of 1. The attach function primarily sets
the pointer to the data location but it also sets the data type of the variable
and its memory location using two enum classes:
```c++
enum class FieldType {Unknown, I4, I8, R4, R8};
enum class FieldMemLoc {Unknown, Device, Host, Both};
```
These provide a simple way to query these properties when needed in some cases
(eg IO) where the type or location are unknown and need to be determined before
the field data is retrieved. The "Both" memory location is used in CPU-only
configurations where the host and device are the same or in future cases where
both may share a memory space. These are determined based on the data type of
the attached data array.

Once a field is defined, various query functions are available. The existence
of a field with a given name can be determined with:
```c++
   if (Field::exists(FieldName)) { do stuff; }
```
Other queries have two forms. If the Field pointer has already been retrieved,
a member function can be used. If not, a query by field name can be used
instead. The data type and memory location can be determined using:
```c++
   FieldType MyType1 = MyField->getType(); // member function
   FieldType MyType2 = Field::getFieldType(FieldName); // name version

   FieldMemLoc MyLoc1 = MyField->getMemoryLocation(); // member function
   FieldMemLoc MyLoc2 = Field::getFieldMemoryLocation(FieldName); // name vers
```
Sometimes it is only necessary to query whether the data exists on the
host rather than retrieving the full memory location.
```c++
   if (MyField->isOnHost()) { do stuff; } // member function
   if (Field::isFieldOnHost(FieldName)) { do stuff }; // name version
```
When looping through all defined fields, it is handy to have a query for the
field name, so we provide:
```c++
   std::string FieldName = MyField->getName();
```
The dimension information can be retrieved using:
```c++
   int NDims = MyField->getNumDims();
   std::vector<std::string> MyDimNames(NDims);
   int Err = MyField->getDimNames(MyDimNames);
```
Once the dimension names have been retrieved, the Dimension class API can be
used to extract further dimension information. Two other field quantities
can be retrieved, but are used only by the IOStream capability:
```c++
   bool IsTimeDependent = MyField->isTimeDependent();
   bool IsDistributed   = MyField->isDistributed();
```
The first determines whether the unlimited time dimension should be added
during IO operations. The second determines whether any of the dimensions
are distributed across MPI tasks so that parallel IO is required.

The data and metadata stored in a field can be retrieved using several
functions.  To retrieve a pointer to the full Field, use:
```c++
   std::shared_ptr<Field> MyField = Field::get(MyFieldName);
```
With this pointer all the member functions above can be used.
The Metadata associated with a field can be retrieved individually using:
```c++
   int Err = MyField->getMetadata(MetadataName, MetaValue);
```
where the MetaValue can be a scalar of any supported data type (I4, I8, R4, R8,
bool, std::string). If the value of a metadata entry needs to be changed,
an update function is provided:
```c++
   int Err = MyField->updateMetadata(MetadataName, NewMetaValue);
```
The existence of a metadata entry can be determined with:
```c++
   bool MetaExists = MyField->hasMetadata(MetadataName);
```
The entire group of metadata is stored in a Metadata type which is simply an
alias for a ``std::map<std::string, std::any>``. This collection of metadata
can be retrieved using:
```c++
   std::shared_ptr<Metadata> ThisMeta = ThisField->getAllMetadata();
   // or
   std::shared_ptr<Metadata> ThisMeta = Field::getFieldMeta(FieldName);
```
This can be useful if we need to extract all the metadata at once (eg during
the IO write phase). However, because the value is stored as a ``std::any``,
it must be coerced to the proper data type using the ``std::any`` type query
and the ``std::any_cast`` coercion function.

To retrieve the field data arrays, there are a few methods available. If the
array forms are needed, there are templated retrievals by either a member
function or a by-name interface:
```c++
   HostArray1DI4 MyData1 = MyField->getDataArray<HostArray1DI4>();
   Array2DR8 MyData2 = Field::getFieldDataArray<Array2DR8>(FieldName);
```
where all of the array and host array types are supported. If the array type
is not known in advance, the field can be queried for both type and memory
location as described previously.

Metadata can be removed from a Field using:
```c++
   int Err = MyField->removeMetadata(MetaName);
   MyField->removeAllMetadata();
```
depending on whether a single metadata entry or all metadata entries need to
be deleted. Entire fields can be removed using:
```c++
   int Err = Field::destroy(FieldName);
```
and before exiting, all fields should be removed using:
```c++
   Field::clear();
```

As mentioned above, Fields can be assigned to groups to provide an easy way
to reference fields that commonly appear together, especially when listing
contents of fields in IO files. Internally, a field group is implemented as
a simple set of field names stored as a ``std::set<std::string>``. A group
is created by first creating an empty group with the desired name:
```c++
   std::shared_ptr<FieldGroup> MyGroup = FieldGroup::create(GroupName);
```
Fields can then be added either through a member function if the group
pointer is available, or by group name:
```c++
   int Err = MyGroup->addField(FieldName);
   int Err = FieldGroup::addFieldToGroup(FieldName, GroupName);
```
The latter is useful especially if the group was created elsewhere. If the
field has already been added to the group, no additional entries are created.
To determine whether a group exists or whether a field is in a group, several
forms of these queries are available:
```c++
   bool MyGroupExists = FieldGroup::exists(GroupName);
   bool FieldIsInGroup = MyGroup->hasField(FieldName);
   bool FieldIsInGroup = FieldGroup::isFieldInGroup(FieldName,GroupName);
```
In addition, the entire list of fields current assigned to a group can be
retrieved using:
```c++
   std::set<std::string> MyFieldList = MyGroup->getFieldList();
   // or
   std::set<std::string> MyFieldList =
                   FieldGroup::getFieldListFromGroup(GroupName);
```
While this list can be used with the Field interfaces above to retrieve the
Field, we also provide a shortcut to retrieve a full field pointer from a
FieldGroup:
```c++
   std::shared_ptr<Field> MyField = MyGroup->getField(FieldName);
   // or
   std::shared_ptr<Field> MyField =
         FieldGroup::getFieldFromGroup(FieldName, GroupName);
```
The FieldGroup pointer can also be retrieved:
```c++
   std::shared_ptr<FieldGroup> MyGroup = FieldGroup::get(GroupName);
```

A field can be removed from a field group using:
```c++
   int Err = FieldGroup::removeField(FieldName);
   // or
   int Err = removeFieldFromGroup(FieldName, GroupName);
```
The removal of a field from a group does not remove the field itself, it only
removes the field name from the list of fields assigned to the group.
The entire group can be removed with:
```c++
   int Err = FieldGroup::destroy(GroupName);
```
and the usual ``FieldGroup::clear();`` should be used to remove all field
groups before exiting.
