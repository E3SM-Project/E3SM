(omega-dev-metadata)=

# Metadata

## MetaData Code Structure

The MetaData feature in Ometa is constructed using four distinct C++ classes.
These classes are designed to facilitate the organization and retrieval of
metadata throughout the Omega source code.

### Class Descriptions

#### 1. MetaData Class
- **Purpose**: Serves as the central class, holding key-value pair metadata
               for variables and other entities such as Ometa code or
               the model itself.
- **Structure**: Contains a static member that maintains references to all
               instances of the class, allowing for universal access
               throughout the codebase.

#### 2. MetaDim Class
- **Purpose**: Provides a structured representation of dimension information,
               primarily for array-type entities.
- **Structure**: Defines dimension-related attributes and functionalities to
               aid in handling array data.

#### 3. MetaGroup Class
- **Purpose**: Offers a mechanism to aggregate certain MetaData instances for
               streamlined management and retrieval.
- **Structure**: Facilitates grouping and accessing grouped metadata,
               enhancing organization and accessibility.

#### 4. ArrayMetaData Class
- **Purpose**: Specializes the MetaData class to handle array-type variables
               specifically, ensuring that required metadata is consistently
               provided.
- **Structure**: Inherits from MetaData and introduces additional constraints
               and functionalities pertinent to array-type data.

## Function Naming Convention
The naming of functions within these classes adheres to the LLVM coding style
and is categorized based on the nature of the operation and the type of member
variables they interact with. Function names begin with a lowercase letter.

### Specific Conventions
- **Static Member Access**: Functions interacting with static members are
               typically short verbs. For instance, `get("Name");` is utilized
               to instantiate a MetaData class.
- **Non-Static Member Access**: Functions dealing with non-static members often
               conclude with a specific descriptor indicating the targeted
               information, e.g., `getEntry("Name", Var);` retrieves
               metadata from a MetaData instance.
- **Instance Management**: The terms `create/destroy` are designated for
               functions that instantiate or delete an instance, respectively.
- **Instance Retrieval**: The `get` function is employed to obtain an instance
               from the static member holding all instances.
- **Existence Check**: The `has("Name");` function returns a boolean value,
               indicating whether an instance named "Name" exists.

## Return Type and Return Code

The create and get static functions are designed to return a `shared_ptr` of
the class instance. This design choice simplifies the usage and enhances
code readability by allowing the use of the auto keyword, eliminating the
need for repeatedly specifying the `shared_ptr` type.

Other Functions typically returns an integer code indicating the operation
status representing a successful operation with no errors. Negative Values
correspond to specific error cases, each value indicating a distinct error
type.

## Inheriting MetaData class

The MetaData class is a versatile construct intended for attaching metadata
not only to model variables but also to other elements related to the model
and code. It serves as a foundational class for metadata management.

The ArrayMetaData class is a specialized version of the MetaData class,
focusing specifically on array-type model variables. It defines a singular
method, `create`, which is an overridden version of the one in the MetaData
class. The create function in ArrayMetaData is tailored to ensure that all
necessary metadata for array-type variables is provided upon creation.
This enforcement guarantees the integrity and completeness of metadata
for array-type model variables.

## Create/Destroy a MetaData Instance

To add or retrieve metadata for a variable or other entity such as code or
a model, users have to create a MetaData instance for the variable or entity.

There are three ways to create a MetaData Instance. The simplest way is to use
the create method with a metadata name string:

```
auto Data1 = MetaData::create("MyVar1");

```

The create method is a static function that returns a C++ `shared_ptr` object
pointing to the MetaData instance named "MyVar1".

If the specified name already exists, a runtime error is thrown.

Once it is created successfully, users can add metadata later, as explained
in the section below.

Another way to create a MetaData Instance is to use the create function with
`std::make_pair`, consisting of name-value pairs:

```
auto Data2 = MetaData::create(
                "MyVar2",
                {
                    std::make_pair("Meta1", 1),
                    std::make_pair("Meta2", true),
                    std::make_pair("Meta3", "MyValue3"),
                }
            );
```
This allows multiple key-value pairs to be added at once during creation.

In the case of an array-type variable, it is mandatory to add a certain set
of key-value pairs. Therefore, users should use the following method to create
array-type variables:

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

The `destroy` function, used with the metadata name, removes a MetaData
instance so that the instance can no longer be accessed.

```
int ret;
ret = MetaData::destroy("MyVar3");
```
The return value of the `destroy` function is zero on success, -1 when there
is no metadata by that name, and -2 when the destroy action fails for
other reasons.

## Add/Remove a Metadata to/from a MetaData Instance

To add metadata to a `MetaData` instance, the `addEntry` function is used:

```
    int ret;
    const R8 AValue = 2.0;
    ret = Data1->addEntry("NewMeta", AValue);
```

The following data types are allowed as metadata values:

- I4
- I8
- R4
- R8
- std::string
- bool

`addEntry` returns an integer value of zero upon success and returns -1
if the metadata name already exists.

To remove metadata from a MetaData instance, the `removeEntry` function
is used.

```
    int ret;
    ret = Data1->removeEntry("NewMeta");
```

`removeEntry` returns an integer value of zero upon success, -1 when there
is no metadata by that name, and -2 when the removal action fails for other
reasons.

## Retreive a Metadata from a MetaData Instance

Users can retrieve metadata using the `getEntry` function.

```
    R8 R8Value;
    ret = Data1->getEntry("NewMeta", R8Value);
```

The first argument of `getEntry` is the name of the metadata. The second
argument is passed by reference, allowing `getEntry` to place
the metadata's value into it.

Note that `getEntry` is overloaded with several methods, each having
different data types for the second argument. It is the user's responsibility
to match the metadata name with the correct data type of the value. Otherwise,
the function may throw a type-casting error exception.

## Retreive a MetaData Instance

To retrieve an instance of MetaData for a variable or another entity, such
as code, the `get` static function can be used.

```
    auto Data4 = MetaData::get("MyVar2");
```

The `get` method is a static function that returns a C++ `shared_ptr` object
in the same way as the create function.

## Create/Destroy a MetaDim Instance

In the above example for array-type variables, users need to add metadata for
the dimension information of the variable. In Omega, the `MetaDim` class is
used for representing dimensions.

The `create` function returns a C++ `shared_ptr` object to the MetaDim
instance.

```
const I4 DimLength{3};
auto Dim1 = MetaDim::create("MyDim1", DimLength);
int ret   = MetaDim::destroy("MyDim1");
```

In case the specified name already exists, a runtime error is thrown.

## Retreive Dimension from a MetaDim Instance

Users can retrieve dimension information using `getLength`.

```
    I4 DimLength;
    Dim1->getLength(DimLength);
```

The argument is passed by reference, allowing `getLength` to assign the length
value to it.

As of now, the return value is always zero.

## Create/Destroy a MetaGroup Instance

Because lists of contents in the model configuration file can become lengthy,
it is useful to define and maintain a metadata group. This group is essentially
a shorthand for a set of common fields (e.g., `meshFields` or `prognosticVars`).
Common groups of fields can be defined as metadata groups (MetaGroup). Please
see the metadata design document for more details.

Similarly to the MetaData and MetaDim classes explained above, the create
function returns a C++ `shared_ptr` object to the MetaGroup instance.

```
auto Group1 = MetaGroup::create("MyGroup1");
int ret     = MetaGroup::destroy("MyGroup1");
```

In case the specified name already exists, a runtime error is thrown.

## Add, Retrieve, or Remove a MetaData Instance from a MetaGroup Instance

```
int ret;
const std::string FieldName{"MyField"};

auto Data1  = MetaData::create(FieldName);

ret         = Group1->addField(FieldName);
auto Data2  = Group1->getField(FieldName):
ret         = Group1->removeField(FieldName):
```

To add, retrieve, or remove a field (a MetaData instance) from
a `MetaGroup shared_ptr`, the `addField`, `getField`, and `removeField`
methods can be used, respectively. `addField` and `removeField` return zero
upon success, -1 if the group already contains the field name, and -2 if
another error occurs. `getField` throws an error if no field name is
specified in the argument.

## Utilities

All the classes described here have the following static utility functions:

### `has` static method

This static method returns a boolean value indicating whether an instance
specified by the argument string exists.

```
bool ret1, ret2, ret3;
ret1    = MetaData::has("MyVar");
ret2    = MetaDim::has("MyDim");
ret3    = MetaGroup::has("MyGroup");
```

### `get` static method

This static method returns an instance of the class specified by
the argument string.

```
auto Data   = MetaData::get("MyVar");
auto Dim    = MetaDim::get("MyDim");
auto Group  = MetaGroup::get("MyGroup");
```
