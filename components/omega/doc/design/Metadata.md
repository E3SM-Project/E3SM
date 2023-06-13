
(omega-design-metadata)=

# Metadata

## 1 Overview

All data files from OMEGA and E3SM must be in a self-describing format
to enable easier analysis, archiving and provenance. We describe here
requirements and design for managing metadata within the OMEGA code
itself to enable appropriate metadata in all model output. The capabilities
described here must interact with other parts of the I/O within OMEGA,
particularly the lower-level I/O interfaces that
perform the actual reading/writing to disk and the management of
files or streams and their contents.

## 2 Requirements

### 2.1 Requirement: Data types

Metadata can be in all supported OMEGA data types (I4, I8, R4, R8,
Boolean/logical, strings). Initially only scalar metadata
is required, but extending to vectors/arrays may be desirable.

### 2.2 Requirement: Global and variable metadata

Metadata refers not only to field/variable metadata but
also metadata associated with the model, code base, simulation
or data file. This capability must be able to track metadata
for these global/general metadata in addition to metadata
attached to particular fields or variables.

### 2.3 Requirement: Metadata conventions

Where possible, all metadata must conform to the climate/forecast (CF)
metadata conventions at https://cfconventions.org. 

### 2.4 Requirement: Required metadata

To comply with conventions and E3SM practices, there will be a minimum
set of required metadata and interfaces must enforce this minimum set.
This minimum will be defined later, but for variables, this would
typically include a name (short), units, long name/description,
standard CF name (if exists), `_FillValue`, valid min/max and dimensions.

### 2.5 Requirement: Dimensions

For array variables, a definition of each dimension used is required.
In an unstructured model like OMEGA, the dimension describes the
global extent of the index space for each dimension that is then paired
with mesh fields for the full description of OMEGA coordinate locations.
For time dimensions, we require support for both a fixed length time
dimension (e.g. to support multiple time levels in a restart) as well
as an unlimited dimension to support time series.

### 2.6 Requirement: Available fields

To support other I/O components, we require a means for maintaining
a list of fields available for I/O based on the requested model
configuration. This enables the I/O streams to define contents of
a given file as part of model configuration input by supplying a
list of variable names to be included.

### 2.7 Desired: Metadata groups

Because lists of contents in the model config file can get long,
it would be useful to define and maintain a metadata group that
is simply a shorthand to a set of common fields (e.g. `meshFields` or
`prognosticVars`). Note that these would only be for the purpose
of shortening lists in the config file or managing lists internally.
The full list of fields would still appear in output file metadata.

## 3 Algorithmic Formulation

No algorithms are used.

## 4 Design

The general philosophy is that modules will define the metadata for
each field they own and add it to a list of available fields. This
internally generated list ensures that the list of available fields
remains consistent with the code base and specific model configuration.
Note that this requires most initialization to occur before processing
output streams so that the requested contents can be checked against
the list of defined fields.

Internally, the metadata class will make use of the C++ 
`std::map` that provides a useful capability for name, value
pairs of each data type.

### 4.1 Data types and parameters

#### 4.1.1 Parameters 

To enable general metadata associated with the code,
simulation or file, we define standard names for these
to be used as the "field" name during metadata definition.
For code and simulation metadata, we will use:

```c++
   constexpr std::string codeMeta{"code"};
   constexpr std::string simMeta {"simulation"};
```

For file metadata (eg. time stamps, convention, etc.),
we will use the stream name as the field name for metadata
definition.

#### 4.1.2 Class/structs/data types

The metadata class holds the collection of metadata associated
with a variable or general category (code, sim, file).
This class implemented as a map (name, value pairs) using `std::any`
to accommodate all supported types. It also stores all the
defined metadata in a static vector so the metadata can be
retrieved from anywhere.

```c++
   class Metadata {
   private:

      /// metadata stored in map
      std::map<std::string, std::any> metaMap;

      /// for arrays or vectors, the rank of the array
      int nDims; 

      /// for arrays or vectors, the dimensions
      std::vector<std::shared_ptr<MetaDim>> dimensions;
      
      /// Store and maintain all defined metadata in this vector
      static std::vector<std::shared_ptr<Metadata>> defFields;


   public:
      // [Methods described below]
   };
```

Another class is needed to describe dimensions for arrays
and vectors. Each is a name, length pair. Like the metadata
above, all defined dimensions are stored and maintained as a
static vector. 

```c++
   class MetaDim{

      private:
         /// name, length pair for the dimension
         std::map<std::string,I4> dimension;

         /// Store and maintain all defined dimensions
         static std::vector<std::shared_ptr<MetaDim>> allDims;

      public:
         [methods described below]
   };
```
Finally, to support common groups of metadata, we define a simple
class to define a list of variable names that are members of a
group. Like the classes above, we store and manage the group
definitions in a static vector of defined groups.

```c++
   class MetaGroup{

      private:
         /// name of group
         std::string groupName;

         /// list of fields in group
         std::vector<std::string> fieldList;

         /// Store and maintain all defined groups
         static std::vector<std::shared_ptr<MetaGroup>> allGroups;

      public:
         [methods described below]
   };
```

### 4.2 Methods

#### 4.2.1 Metadata creation

There will be two methods for creating metadata. Because we
are managing all instances of metadata within the class, these
are static functions rather than constructors. Also note in
the comments below that empty or zero values can be provided
in cases where input arguments don't make sense.

```c++
    /// Create an empty instance by name, primarily for global
    /// metadata not assigned to a variable/field. An error
    /// code is returned that is non-zero if metadata of the same
    /// name already exists
    int Metadata::create(const std::string name /// [in] name to assign
                         );

    /// Create metadata for a variable/field using all
    /// required metadata. Note that if input parameters don't exist
    /// (eg stdName) or don't make sense (eg min/max or fill) for a
    /// given field, empty or 0 entries can be provided. Similarly
    /// for scalars, nDims can be 0 and an emtpy MetaDim vector can be
    /// supplied.
    int Metadata::create
         (const std::string name,        /// name of var
          const std::string description, /// long name or description
          const std::string units,       /// units
          const std::string stdName,     /// CF standard name
          const std::any    validMin,    /// min valid field value
          const std::any    validMax,    /// max valid field value
          const std::any    fillValue,   /// scalar used for undefined entries 
          const int         nDims,       /// number of dimensions
          const std::vector<std::shared_ptr<MetaDims>> // dim pointers
          );
```

A destructor will be provided to free metadata and remove from the
defined list, though we anticipate most metadata will need to be
persistent through a given simulation.

#### 4.2.2 Add/remove metadata

An interface for adding new metadata to a previously defined Metadata
instance will be provided. This can be used to add global metadata to
the code, simulation, file metadata, but also for variables if additional
metadata is desired.  

```c++
   int Metadata::Add(const std::string varName, /// name of field to add data to
                     const std::string name,    /// name of new metadata to add
                     const std::any    value    /// value of new metadata to add
                     );
```

where varName is the name of the already defined variable, and
the other arguments are the usual name/value pair for the metadata
to be added. The value is a `std::any` type to support all types within
within the metadata map. This function will return 0 if successful and an error
if either the varName has not been defined or if a metadata pair with
that name already exists.

For symmetry, we will supply a remove function, though the
use case is likely rare.

```c++
   int Metadata::Remove(
            const std::string varName, /// name of field
            const std::string name.    /// name of metadata to remove
            ); 
```

#### 4.2.3 Retrieve metadata

The most common use case will be retrieving metadata. For a single
metadata entry, there will be an explicit get function:

```c++
   myValue = Metadata::Get(varName, name);
```

where varName is the field defined (or generic name for global
metadata) and name is the name associated with the metadata
to be retrieved. The internal `std::any` representation of the
value will be cast into the type expected.

During I/O stages, we will need a capability for retrieving
all metadata for a defined field or global category. We will
supply a function to retrieve the map associated with a given
data type. The `std::map` iterators and functionality can then
be used to extract all of the metadata for a given type. For
example, to retrieve all metadata from a defined field:

```c++
   std::map<std::string, std::any> *varMeta = MetadataGetAll(varName);
   std::string metaName;
   std::any    metaVal;

   // use std::map functionality to loop through all metadata
   for (auto it=varMeta.begin(); it != varMeta.end(), ++it){
      metaName = it->first;  // retrieve name part of meta pair
      metaVal  = it->second; // retrieve value part of meta pair

      // do whatever needed with metadata
      // metaVal can be cast into the appropriate type using
      //   std::any_cast<type>(metaVal)
   }
```

#### 4.2.4 Create dimension

For the dimension class, there will be a create function that will
both construct the instance and add it to the list of defined
dimensions. It will return an integer error code that is non-zero
if the dimension is already defined.

```c++
   int MetaDim::create(std::string name, // name of dimension
                       I4 length         // length of dimension
                       );
```

A delete function will be supplied as well.
```c++
   int MetaDim::delete(std::string name // name of dimension
                       );
```
We use create/delete functions rather than constructors so
that the class can store and manage all instances.

#### 4.2.5 Retrieve dim length

The only other function for a dimension will be a retrieval for
the dimension length:

```c++
   I4 myLength = MetaDim::getLength("name");
```

#### 4.2.6 Existence inquiry

For both metadata and dimensions, a function will be provided
to inquire whether metadata or dimensions with a given name
have already been defined:

```c++
   bool Metadata::isDefined(const std::string name);
   bool MetaDim::isDefined (const std::string name); 
```

#### 4.2.7 Metadata groups

Common groups of fields can be defined as metadata groups (MetaGroup).
We supply create/delete functions for defining a group of a given
name. Group names should start with GRP so that when groups are
included as contents in input config files, we can distinguish
a group from a field. Create will create an empty group to which
field names can be added. The functions will return an error code.
Note that in the case of groups, if a group of the same name
already exists, no error is returned - the create does nothing.
This is to support cases where multiple modules may be adding their
fields to a group during initialization and these may be called
in an arbitrary order.

```c++
   int MetaGroup::create(std::string name /// name of group
                         );
   int MetaGroup::delete(std::string name /// name of group
                         );
```

The members of a group can then be added individually. A remove
function is provided though use cases are rare. As in create,
if a field has already been added, the add function will do nothing.

```c++
   int MetaGroup::add(std::string name /// name of field to add
                      );
   int MetaGroup::remove(std::string name /// name of field to remove
                         );
```

Finally, we will need to be able to extract the names of all
fields that have been added to a group:

```c++
   std::vector<std::string> fieldList = 
       MetaGroup::get(std::string groupName);
```
## 5 Verification and Testing

### 5.1 Test creation/retrieval

Create dimensions for a sample multi-dimensional field and then
create field metadata with metadata entries for all supported types.
Test by retrieving all meta data and dimension information and
ensuring it is identical.
  - tests requirement 2.1, 2.4, 2.5

### 5.2 Test global metadata

Add metadata to supported global metadata fields "code" and
"simulation" and test by retrieving the same and comparing.
  - tests requirement 2.2

### 5.3 Test inquiry functions

Test available field and available dimensions by inquiring if
above metadata and dimensions have been defined. Also test names
that don't exist to test failure modes.
  - tests requirement 2.6

### 5.4 Test metadata groups

Define one or two groups of metadata with both defined fields
and non-existent fields (to check error modes). Retrieve the
list of fields in each group and check against the expected
names.
  - tests requirement 2.7

### 5.5 Test CF compliance

To test CF compliance, we will test output files using a CF checking
code. This may be a part of Polaris testing rather than the unit test
framework.
  - tests requirement 2.3
